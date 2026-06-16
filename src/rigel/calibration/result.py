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

    def __post_init__(self) -> None:
        n = self.n_regions
        if n < 0:
            raise ValueError(f"CalibrationResult.n_regions must be >= 0; got {n}.")

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


__all__ = ["CalibrationResult"]
