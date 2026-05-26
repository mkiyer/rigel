"""Per-region density observation table (v4 Phase 2).

`DensityObservation` is a transient, sorted, frozen table of the gDNA-
relevant counts and FL-aware opportunities consumed by the density
model. It is built once inside the orchestrator from
``(region_arrays, ledger, fl_models.gdna)`` and dropped at function
exit; it is **not** stored on `CalibrationResult`.

The plan (`docs/fineregions/density_model_impl_plan_v4.md` §6) defers
the implementation to v3 §6:

* Counts come from the ledger's unspliced totals (POS + NEG).
* ``contained_leff`` reuses ``l_eff_contained(end - start, gdna_fl)``.
* ``boundary_left_leff`` and ``boundary_right_leff`` reuse
  ``fractional_boundary_side_exposure(end - start, gdna_fl)``. The
  helper depends only on the receiving-region length, so both sides
  share the same per-region opportunity.
* ``anchor_intergenic`` / ``anchor_intron`` use the
  ``fractional_evidence`` predicates; ``is_anchor`` is their union.
* ``spliced_count`` is diagnostic-only for the v5 density model. The additive
    v6 four-state tensor consumes it separately as expression evidence.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..frag_length_model import FragmentLengthModel
from ._arrays import RegionArrays
from ._exposure import fractional_boundary_side_exposure, l_eff_contained
from .fractional_evidence import is_intergenic, is_intron_only
from .region_count_ledger import RegionCountLedger


__all__ = [
    "DensityObservation",
    "build_density_observation",
]


@dataclass(frozen=True, slots=True)
class DensityObservation:
    """Per-region counts + FL-aware opportunities (sorted-row layout)."""

    # Counts (float32[R])
    contained_count: np.ndarray
    boundary_left_count: np.ndarray
    boundary_right_count: np.ndarray
    boundary_count: np.ndarray
    observed_compatible_count: np.ndarray

    # Opportunities (float64[R])
    contained_leff: np.ndarray
    boundary_left_leff: np.ndarray
    boundary_right_leff: np.ndarray
    boundary_leff: np.ndarray

    # Classification (bool[R])
    anchor_intergenic: np.ndarray
    anchor_intron: np.ndarray
    is_anchor: np.ndarray

    # Diagnostics
    spliced_count: np.ndarray
    region_length: np.ndarray


def build_density_observation(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    gdna_fl: FragmentLengthModel,
) -> DensityObservation:
    """Build a `DensityObservation` from sorted region/ledger arrays."""
    contained = ledger.contained_unspliced_total()
    bl = ledger.boundary_left_unspliced_total()
    br = ledger.boundary_right_unspliced_total()
    boundary = bl + br

    spans_bp = (region_arrays.end - region_arrays.start).astype(np.int64, copy=False)
    contained_leff = l_eff_contained(spans_bp, gdna_fl)
    side_leff = fractional_boundary_side_exposure(spans_bp, gdna_fl)
    boundary_leff = side_leff + side_leff

    sig = region_arrays.signature
    a_inter = is_intergenic(sig)
    a_intron = is_intron_only(sig)
    is_anchor = a_inter | a_intron

    spliced = ledger.spliced_total()

    return DensityObservation(
        contained_count=contained.astype(np.float32, copy=False),
        boundary_left_count=bl.astype(np.float32, copy=False),
        boundary_right_count=br.astype(np.float32, copy=False),
        boundary_count=boundary.astype(np.float32, copy=False),
        observed_compatible_count=(contained + boundary).astype(np.float32, copy=False),
        contained_leff=np.asarray(contained_leff, dtype=np.float64),
        boundary_left_leff=np.asarray(side_leff, dtype=np.float64),
        boundary_right_leff=np.asarray(side_leff, dtype=np.float64),
        boundary_leff=np.asarray(boundary_leff, dtype=np.float64),
        anchor_intergenic=np.asarray(a_inter, dtype=bool),
        anchor_intron=np.asarray(a_intron, dtype=bool),
        is_anchor=np.asarray(is_anchor, dtype=bool),
        spliced_count=spliced.astype(np.float32, copy=False),
        region_length=spans_bp.astype(np.float64, copy=False),
    )
