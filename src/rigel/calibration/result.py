"""CalibrationResult — the calibrator's output schema.

Per-region deconvolved gDNA / RNA mass across the region's three nodes (contained
plus the two boundary sides), the per-region gDNA geometric supports (contained
``gdna_region_eff_len`` + per-side ``gdna_boundary_len``), and the two library scalars
(``gdna_density_global``, ``rna_sense_frac``). The calibrator iterates the simplex sweep
to convergence but carries no per-pass convergence diagnostics in this schema. ``__post_init__``
enforces the intrinsic invariants (shapes, finiteness, sign); mass conservation against the
raw fragment counts is checked by the calibrator / tests (it needs the substrate, which the
result does not carry).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CalibrationConfig


def _check_region_array(arr: np.ndarray, name: str, n_regions: int) -> None:
    """Validate a per-region float64 array: shape, dtype, finite, non-negative."""
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"CalibrationResult.{name} must be a numpy array.")
    if arr.dtype != np.float64:
        raise ValueError(f"CalibrationResult.{name} must be float64; got {arr.dtype}.")
    if arr.shape != (n_regions,):
        raise ValueError(
            f"CalibrationResult.{name} has shape {arr.shape}; expected ({n_regions},)."
        )
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"CalibrationResult.{name} contains non-finite values.")
    if np.any(arr < 0.0):
        raise ValueError(f"CalibrationResult.{name} must be non-negative.")


@dataclass(frozen=True, slots=True)
class RnaWarmStart:
    """Per-region SOLVED per-strand RNA densities (ρ ≈ true abundance θ) — the EM warm-start signal.

    Consumed by the pipeline-layer warm-start builder (which owns the transcript index; calibration does
    not) to compute a capture-corrected density BOTTLENECK per transcript: it maps a transcript to its
    nodes and takes the min over three node ROLES, ready-to-use (the only arithmetic left to the builder is
    that min plus the capture correction ε it reconstructs from the gDNA mass fields). gDNA is NOT
    re-exposed — a node's gDNA density is reconstructible from ``mass_gdna_contained / gdna_region_eff_len``
    and the pooled seam ``(mass_gdna_right[r] + mass_gdna_left[r+1]) / mean(gdna_boundary_len[r, r+1])``.

    * **CONTAINED** (region node, nascent+mature): ``rho_contained_s[r] = f_s · M_u / E_rna_contained`` —
      per strand, from the region node's solved belief.
    * **CROSSING** (interior seam between region ``r`` and ``r+1``, LEFT-region keyed; nascent):
      ``rho_crossing_s[r] = f_s · M_seam / E_seam`` — per strand, from the seam BOUNDARY node's belief;
      ``0`` on each reference's last region (no seam to its right).
    * **SPLICED** (mature) — per region SIDE, in its OWN arrays so it is NEVER summed into the crossing
      density (a splice site is also an exon↔intron boundary; summing would let a neighbour's mature mass
      inflate a nascent shadow's bottleneck → the gDNA→nascent siphon). Each junction is single-stranded
      (the splice motif fixes the strand), so at a normal site exactly one of the pos/neg arrays is nonzero;
      both are nonzero only where an overlapping ``+`` and ``−`` gene share the exact splice coordinate (a
      dense-AMBIG coincidence), which per-strand storage keeps separate — the antisense isoform is read on
      its OWN strand, not crushed. ``rho_spliced_s_{left,right}[r] = spliced_s / E_spl``. A junction
      ``(r_left, r_right)`` reads the donor at ``rho_spliced_s_right[r_left]`` and the acceptor at
      ``rho_spliced_s_left[r_right]`` — the same per-strand ``_pick`` the contained/crossing roles use.

    All arrays are ``float64`` of length ``R`` (region-keyed, mirroring the ``mass_*`` layout).
    """

    rho_contained_pos: np.ndarray  # float64[R] — f_pos·M_u / E_rna (contained region node)
    rho_contained_neg: np.ndarray
    rho_crossing_pos: np.ndarray  # float64[R] — seam(r,r+1) f_pos·M_seam / E_seam (left-keyed; 0 at terminals)
    rho_crossing_neg: np.ndarray
    rho_spliced_pos_left: np.ndarray  # float64[R] — mature +ρ on region r's LEFT side (acceptor)
    rho_spliced_neg_left: np.ndarray
    rho_spliced_pos_right: np.ndarray  # float64[R] — mature +ρ on region r's RIGHT side (donor)
    rho_spliced_neg_right: np.ndarray

    def __post_init__(self) -> None:
        n = self.rho_contained_pos.shape[0]
        for name in (
            "rho_contained_pos",
            "rho_contained_neg",
            "rho_crossing_pos",
            "rho_crossing_neg",
            "rho_spliced_pos_left",
            "rho_spliced_neg_left",
            "rho_spliced_pos_right",
            "rho_spliced_neg_right",
        ):
            _check_region_array(getattr(self, name), f"RnaWarmStart.{name}", n)


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Per-region deconvolved mass + geometric gDNA length + library scalars."""

    # --- deconvolved mass across the region's 3 nodes (float64[R]) ---
    mass_gdna_contained: np.ndarray
    mass_rna_contained: np.ndarray
    mass_gdna_left: np.ndarray  # right side of the left boundary
    mass_rna_left: np.ndarray
    mass_gdna_right: np.ndarray  # left side of the right boundary
    mass_rna_right: np.ndarray

    # Spliced RNA mass per region (Σ over the 3 nodes). Carried so ``assemble_priors`` can
    # WITHHOLD it from ``rna_prior_count``: a spliced fragment has no gDNA candidate in the EM
    # (gDNA does not splice), so it is guaranteed-RNA and assigned directly — counting it in the
    # prior would double-count it and unfairly inflate the RNA side of the gDNA-vs-RNA *unspliced*
    # split. ``mass_rna_*`` themselves stay spliced-inclusive, so the per-node conservation
    # ``mass_gdna + mass_rna = total node mass`` is preserved; only the prior subtracts this.
    mass_rna_spliced: np.ndarray  # float64[R]

    # --- directional boundary per-side density length 𝓔(L)=E[min(ℓ,L)] per region (geometry) ---
    # The mass a boundary-crossing fragment deposits on one flank of region r. Doubles as the
    # POOLED-SEAM effective support: a seam between regions r and r+1 has support
    # ½·(gdna_boundary_len[r] + gdna_boundary_len[r+1]) — the deposition-faithful divisor that makes
    # the pooled seam mass ρ·support under uniform gDNA (see assemble_priors / capture_eff_length).
    gdna_boundary_len: np.ndarray  # float64[R]

    # --- region-contained gDNA effective support E[max(0,L−ℓ)] per region (geometry) ---
    # The count of contained-fragment start positions = the effective sampling support of the
    # region's contained gDNA mass. Under uniform genomic gDNA deposition at density ρ, the
    # expected contained mass is EXACTLY ρ·gdna_region_eff_len (a fragment must FIT to be
    # contained), so dividing the contained mass by this recovers the true density ρ — the
    # bedrock "factor = 1 under uniform gDNA" invariant. This, NOT region_size_bp, is the
    # density-correct IPR divisor for the region node: region_size_bp understates the density of
    # short regions (it ignores the fit-inside constraint), manufacturing a spurious contraction
    # in an unenriched library. Consumed by assemble_priors / capture_eff_length.
    gdna_region_eff_len: np.ndarray  # float64[R]

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
    rna_sense_frac: float  # in [0, 1], RNA sense fraction used by the strand clue
    gdna_strand_overdispersion: float  # in [0, 1), fitted gDNA strand Beta-Binomial dispersion
    rna_strand_overdispersion: float  # in [0, 1), fitted RNA strand Beta-Binomial dispersion

    # --- provenance ---
    n_regions: int
    config: CalibrationConfig

    # --- EM warm-start signal (Phase A; additive, default None) ---
    # Per-region SOLVED per-strand RNA densities for the pipeline-layer warm-start builder (see
    # ``RnaWarmStart``). Not read by ``assemble_priors`` or any serialized output, so populating it leaves
    # the default calibration output byte-identical. ``None`` on results that predate / skip the warm start.
    rna_warm_start: "RnaWarmStart | None" = None

    def __post_init__(self) -> None:
        n = self.n_regions
        if n < 0:
            raise ValueError(f"CalibrationResult.n_regions must be >= 0; got {n}.")
        if self.rna_warm_start is not None:
            if not isinstance(self.rna_warm_start, RnaWarmStart):
                raise ValueError("CalibrationResult.rna_warm_start must be an RnaWarmStart or None.")
            if self.rna_warm_start.rho_contained_pos.shape != (n,):
                raise ValueError(
                    "CalibrationResult.rna_warm_start arrays must have length n_regions "
                    f"({n}); got {self.rna_warm_start.rho_contained_pos.shape}."
                )

        for name in (
            "mass_gdna_contained",
            "mass_rna_contained",
            "mass_gdna_left",
            "mass_rna_left",
            "mass_gdna_right",
            "mass_rna_right",
            "mass_rna_spliced",
            "gdna_boundary_len",
            "gdna_region_eff_len",
        ):
            _check_region_array(getattr(self, name), name, n)

        if not np.isfinite(self.gdna_density_global) or self.gdna_density_global < 0.0:
            raise ValueError(
                f"CalibrationResult.gdna_density_global must be finite and >= 0; got {self.gdna_density_global}."
            )
        sense_frac = float(self.rna_sense_frac)
        if not 0.0 <= sense_frac <= 1.0:
            raise ValueError(
                f"CalibrationResult.rna_sense_frac must be in [0, 1]; got {sense_frac}."
            )
        overdispersion = float(self.gdna_strand_overdispersion)
        if not 0.0 <= overdispersion < 1.0:
            raise ValueError(
                "CalibrationResult.gdna_strand_overdispersion must be in [0, 1); "
                f"got {overdispersion}."
            )
        rna_overdispersion = float(self.rna_strand_overdispersion)
        if not 0.0 <= rna_overdispersion < 1.0:
            raise ValueError(
                "CalibrationResult.rna_strand_overdispersion must be in [0, 1); "
                f"got {rna_overdispersion}."
            )


__all__ = ["CalibrationResult", "RnaWarmStart"]
