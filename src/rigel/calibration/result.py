"""CalibrationResult: the calibrator's output schema (layer 7).

The splice graph has three kinds of object, so the result has three axes::

    regions                N             deconvolved contained count and its geometric support
    contiguous boundaries  E = N − refs  deconvolved crossing count and its geometric support
    sj boundaries          J             the jumping flux: certified RNA, never deconvolved

A contiguous boundary is a 0-bp object with one set of numbers, never a left/right pair, so no consumer pools faces
and halves them. ``__post_init__`` enforces the intrinsic invariants (per-axis shape, dtype, finiteness, sign); mass
conservation against the raw counts needs the substrate, which the result does not carry.

Gate: ``tests/calibration/test_result_schema.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..config import CONSTANTS


def _check_axis_array(arr: np.ndarray, name: str, n: int) -> None:
    """Validate one per-object array: shape against its own axis, dtype (float64 or an integer count), finite and
    non-negative.

    The shape check matters most: ``E = N − n_refs`` differs from ``N`` by only a few hundred genome-wide, so an array
    keyed to the wrong axis has a plausible length and would silently read the wrong object's number.
    """
    if not isinstance(arr, np.ndarray):
        raise ValueError(f"CalibrationResult.{name} must be a numpy array.")
    if arr.dtype != np.float64 and arr.dtype.kind not in ("i", "u"):
        raise ValueError(
            f"CalibrationResult.{name} must be float64 or an integer count; got {arr.dtype}."
        )
    if arr.shape != (n,):
        raise ValueError(f"CalibrationResult.{name} has shape {arr.shape}; expected ({n},).")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"CalibrationResult.{name} contains non-finite values.")
    if np.any(arr < 0.0):
        raise ValueError(f"CalibrationResult.{name} must be non-negative.")


def _check_unit_interval(value: float, name: str, *, open_upper: bool = False) -> None:
    v = float(value)
    upper_ok = v < 1.0 if open_upper else v <= 1.0
    if not (v >= 0.0 and upper_ok):
        bound = "[0, 1)" if open_upper else "[0, 1]"
        raise ValueError(f"CalibrationResult.{name} must be in {bound}; got {v}.")


#: The per-object arrays on each axis, each checked against its own axis length (:func:`_check_axis_array`).
_REGION_ARRAYS = (
    "count_gdna_region",
    "count_rna_region",
    "gdna_region_eff_len",
    "gdna_frac_region",
    "rna_pos_frac_region",
    "rna_neg_frac_region",
    "gdna_capture_efficiency_region",
)
_BOUNDARY_ARRAYS = (
    "count_gdna_boundary",
    "count_rna_boundary",
    "count_rna_spliced_boundary",
    "boundary_mass_per_crossing",
    "boundary_spliced_mass_per_crossing",
    "gdna_boundary_eff_len",
    "gdna_boundary_conserved_len",
    "gdna_frac_boundary",
    "rna_pos_frac_boundary",
    "rna_neg_frac_boundary",
    "gdna_capture_efficiency_boundary",
)
_SJ_ARRAYS = ("count_rna_sj", "sj_mass_per_crossing")
#: The three-way composition on both axes: fractions, so bounded by 1.
_FRACTIONS = (
    "gdna_frac_region",
    "rna_pos_frac_region",
    "rna_neg_frac_region",
    "gdna_frac_boundary",
    "rna_pos_frac_boundary",
    "rna_neg_frac_boundary",
)


@dataclass(frozen=True, slots=True)
class CalibrationResult:
    """Deconvolved gDNA / RNA mass on the region and contiguous-boundary axes, the sj flux, the two
    gDNA geometric supports, and the library scalars."""

    # --- the deconvolved mixture, per region (float64[n_regions]) ---
    #: ``chain_region_deconv``: the region's contained unspliced count split by the converged belief.
    count_gdna_region: np.ndarray
    count_rna_region: np.ndarray

    # --- the deconvolved mixture, per contiguous boundary (float64[n_boundaries]) ---
    #: ``chain_boundary_deconv``: the boundary's unspliced crossing count split by the converged belief.
    #: ``count_rna_boundary`` includes the boundary's certified-RNA crossings (gDNA cannot splice), so
    #: ``count_gdna_boundary + count_rna_boundary == unspliced + spliced`` per boundary.
    count_gdna_boundary: np.ndarray
    count_rna_boundary: np.ndarray

    #: float64[n_boundaries]: the ``boundary_spliced`` part of ``count_rna_boundary``, molecules that crossed this
    #: boundary contiguously having spliced elsewhere, carried apart so a consumer can separate certified from
    #: deconvolved RNA. A region has no twin: its contained population cannot hold a spliced molecule.
    count_rna_spliced_boundary: np.ndarray

    #: float64[n_boundaries]: the incidence-to-fragment conversion, ``mass / count`` off the accumulator's
    #: conserved-mass bank, i.e. the mean fragment mass one crossing at this boundary carries. It is geometry, not a
    #: deconvolved mass: a property of the partition and the fragment-length law, identical under any gDNA/RNA split,
    #: so an oracle arm never overrides it. ``assemble_priors`` multiplies the gDNA boundary mass by it, because a
    #: fragment deposits +1 on each of the ``max(K, 1)`` boundaries it crosses and the EM reads a fragment count. It
    #: is 1.0 where both flanking regions exceed every fragment length and falls toward the region spacing otherwise.
    boundary_mass_per_crossing: np.ndarray

    # --- the jumping population, per sj boundary (float64[n_sj]) ---
    #: An incidence count, not a mass: ``sj_count`` summed over the genome-strand columns. A fragment deposits +1 on
    #: every sj it uses, so this over-counts the conserved mass in proportion to how many sj a fragment spans;
    #: :attr:`sj_conserved_mass` is the mass. Never deconvolved, since an sj boundary is pure RNA. ``assemble_priors``
    #: does not read it; it is exported for QC and reporting, because this population dominates a donor boundary.
    count_rna_sj: np.ndarray

    #: float64[n_boundaries] and float64[n_sj]: the ``boundary_spliced`` and sj twins of
    #: ``boundary_mass_per_crossing`` (the sj one is ``sj_mass / sj_count``), geometry in the same sense. The sj one
    #: makes a conserved library fragment count computable: a spliced fragment whose blocks all lie inside one region
    #: crosses no boundary and is not contained, so only the sj bank holds it.
    boundary_spliced_mass_per_crossing: np.ndarray
    sj_mass_per_crossing: np.ndarray

    # --- the gDNA geometric supports: expected admissible start positions, per object ---
    #: float64[n_regions]: gDNA's contained placements over the region's observable pieces
    #: (`effective_length.observable_contained_eff_length`), ``E_f[(region_len − w + 1)+]`` with nothing masked. Under
    #: uniform gDNA at density ρ the expected contained mass is exactly ``ρ · gdna_region_eff_len`` (a fragment must
    #: fit to be contained), so this, not the region length, is the density-correct divisor.
    gdna_region_eff_len: np.ndarray
    #: float64[n_boundaries]: gDNA's crossing placements over the observable genome
    #: (`effective_length.observable_crossing_eff_length`). gDNA's template is the chromosome, so each side's reach
    #: runs to the nearest masked base and the divisor is ``mu_g − 1`` where nothing is masked.
    gdna_boundary_eff_len: np.ndarray
    #: float64[n_boundaries]: gDNA's conserved share at each boundary: of the starts whose fragments cross it, the
    #: part of each placement's unit the deposit rule gives this boundary (`effective_length.observable_cut_shares` on
    #: the gDNA pmf). With ``gdna_region_eff_len`` it partitions gDNA's starts over the objects exactly, a start
    #: counted once however many boundaries its fragment crosses, so a locus's gDNA length is its objects' shares at
    #: their efficiencies (`priors`). Geometry: identical under any split.
    gdna_boundary_conserved_len: np.ndarray

    # --- the three-way composition ψ solves, per object ---
    #: float64[n_regions] / float64[n_boundaries]: the solved ``(f_g, f_pos, f_neg)`` at each object, the gDNA share
    #: and the two RNA strand shares of its unspliced population, over
    #: ``{gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}``. ``count_gdna_*`` and ``count_rna_*`` are this answer with
    #: the two RNA strands summed. The three close by construction (ψ's ``compose``, ``native/psi_kernel.h``), a slot
    #: the solve does not reach keeps its signature-binary belief, and they are published as solved, never
    #: renormalised.
    gdna_frac_region: np.ndarray
    rna_pos_frac_region: np.ndarray
    rna_neg_frac_region: np.ndarray
    gdna_frac_boundary: np.ndarray
    rna_pos_frac_boundary: np.ndarray
    rna_neg_frac_boundary: np.ndarray

    # --- library scalars ---
    gdna_density_global: float  # >= 0, global gDNA density (mass/bp); 0 in a zero-gDNA library
    #: The unit the capture efficiencies are read in: the peak density of the last landscape's largest basin, the
    #: bulk, in gDNA fragments per bp (`capture_efficiency`). ``None`` when no landscape was fitted, and then every
    #: efficiency is exactly 1; otherwise positive and finite.
    gdna_bulk_density: float | None
    #: The landscape's saturated level against the bulk: the most massive enriched basin's peak over the bulk peak,
    #: ``> 1``, or exactly 1 when nothing lies above the bulk or no landscape was fitted. The ruler's junction price
    #: saturates at it (`capture_eff_length`).
    gdna_capture_saturation: float
    #: float64[n_regions]: each region's capture efficiency ``exp(E[min(x, s(x))] − x_b)``, its own gDNA density
    #: against the bulk with the landscape as a two-group prior and its contained count as the likelihood
    #: (`capture_efficiency.capture_efficiencies`). A relative yield: 1 at the bulk, above 1 under capture, and exactly
    #: 1 everywhere with no bulk density. The ruler and the locus prior price every contained share at it.
    gdna_capture_efficiency_region: np.ndarray
    #: float64[n_boundaries]: each boundary's capture efficiency, from its crossing count on its crossing support, in
    #: the same unit. The ruler and the locus prior price every conserved boundary share at it, and the ruler prices a
    #: junction from the efficiencies beside it (`capture_eff_length`).
    gdna_capture_efficiency_boundary: np.ndarray
    rna_sense_frac: float  # in [0, 1], RNA sense fraction used by the strand clue

    # --- provenance: the three axis lengths, independent of each other ---
    n_regions: int
    n_boundaries: int
    n_sj: int

    def __post_init__(self) -> None:
        for axis in ("n_regions", "n_boundaries", "n_sj"):
            if int(getattr(self, axis)) < 0:
                raise ValueError(
                    f"CalibrationResult.{axis} must be >= 0; got {getattr(self, axis)}."
                )
        for names, axis in (
            (_REGION_ARRAYS, self.n_regions),
            (_BOUNDARY_ARRAYS, self.n_boundaries),
            (_SJ_ARRAYS, self.n_sj),
        ):
            for name in names:
                _check_axis_array(getattr(self, name), name, axis)
        self._check_fractions()
        self._check_capture()
        if not np.isfinite(self.gdna_density_global) or self.gdna_density_global < 0.0:
            raise ValueError(
                "CalibrationResult.gdna_density_global must be finite and >= 0; "
                f"got {self.gdna_density_global}."
            )
        _check_unit_interval(self.rna_sense_frac, "rna_sense_frac")

    def _check_fractions(self) -> None:
        """Each composition component is a fraction, bounded by 1 up to rounding; closure holds by construction (ψ's
        ``compose``)."""
        ceiling = 1.0 + CONSTANTS.calibration.fraction_tolerance
        for name in _FRACTIONS:
            arr = np.asarray(getattr(self, name), dtype=np.float64)
            if np.any(arr > ceiling):
                i = int(np.flatnonzero(arr > ceiling)[0])
                raise ValueError(
                    f"CalibrationResult.{name} is a fraction and must not exceed 1; index {i} is "
                    f"{float(arr[i])!r}."
                )

    def _check_capture(self) -> None:
        """The capture reading: a relative yield against the bulk (finite, non-negative, no upper bound), exactly 1
        everywhere when there is no bulk density, and a saturation level of at least the bulk itself."""
        c = np.asarray(self.gdna_capture_efficiency_region, dtype=np.float64)
        cb = np.asarray(self.gdna_capture_efficiency_boundary, dtype=np.float64)
        if self.gdna_bulk_density is None and (
            np.any(c != 1.0) or np.any(cb != 1.0) or self.gdna_capture_saturation != 1.0
        ):
            raise ValueError(
                "CalibrationResult.gdna_capture_efficiency_region and _boundary, and "
                "gdna_capture_saturation, must be exactly 1 when there is no bulk density: with no landscape "
                "nothing is read against anything."
            )
        if self.gdna_bulk_density is not None and not (
            np.isfinite(self.gdna_bulk_density) and self.gdna_bulk_density > 0.0
        ):
            raise ValueError(
                "CalibrationResult.gdna_bulk_density must be None or finite and > 0; "
                f"got {self.gdna_bulk_density}."
            )
        if not (np.isfinite(self.gdna_capture_saturation) and self.gdna_capture_saturation >= 1.0):
            raise ValueError(
                "CalibrationResult.gdna_capture_saturation must be finite and >= 1 (the bulk itself); "
                f"got {self.gdna_capture_saturation}."
            )

    # ── the conserved sj mass: the third axis in fragment units ──

    @property
    def sj_conserved_mass(self) -> np.ndarray:
        """float64[n_sj]: the conserved fragment mass at each sj, ``count_rna_sj · sj_mass_per_crossing``.

        This recovers the accumulator's ``sj_mass`` bank up to float rounding. It must stay a plain product:
        ``mass_per_crossing`` is 1.0 where the count is zero, so a branch that fell back to the factor would publish
        phantom mass at every sj nothing crossed. It is a property rather than a stored field because an oracle arm
        swaps ``count_rna_sj`` with ``dataclasses.replace``, and a cached array would survive the swap.
        """
        return np.asarray(self.count_rna_sj, dtype=np.float64) * np.asarray(
            self.sj_mass_per_crossing, dtype=np.float64
        )

    # ── the library fragment counts, in fragment units ──

    @property
    def library_gdna_fragments(self) -> float:
        """gDNA fragments in the library: a conserved count, not an object-incidence sum.

        gDNA cannot splice, so it sits on two axes; containment is exclusive, so only the crossing term needs
        converting.
        """
        return float(
            np.asarray(self.count_gdna_region, dtype=np.float64).sum()
            + (
                np.asarray(self.count_gdna_boundary, dtype=np.float64)
                * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
        )

    @property
    def library_rna_fragments(self) -> float:
        """RNA fragments in the library over all three axes, each converted by its own population's mass per crossing.

        The raw incidences must never be summed instead: a fragment books ``max(K, 1)`` boundary crossings and one
        incidence per sj it uses, so the sum over-counts, and differently for the two components because gDNA cannot
        splice. Adding the spliced crossings and the sj flux is not double counting: the deposit rule splits a spliced
        fragment's bases between ``boundary_spliced_mass`` and ``sj_mass``. A property for the same reason as
        :attr:`sj_conserved_mass`.
        """
        unspliced_boundary = np.maximum(
            np.asarray(self.count_rna_boundary, dtype=np.float64)
            - np.asarray(self.count_rna_spliced_boundary, dtype=np.float64),
            0.0,
        )
        return float(
            np.asarray(self.count_rna_region, dtype=np.float64).sum()
            + (
                unspliced_boundary * np.asarray(self.boundary_mass_per_crossing, dtype=np.float64)
            ).sum()
            + (
                np.asarray(self.count_rna_spliced_boundary, dtype=np.float64)
                * np.asarray(self.boundary_spliced_mass_per_crossing, dtype=np.float64)
            ).sum()
            # the sj conversion has one home, `sj_conserved_mass`
            + self.sj_conserved_mass.sum()
        )


__all__ = ["CalibrationResult"]
