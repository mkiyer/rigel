"""CalibrationSubstrate — the calibrator-facing view of the accumulator payload.

The substrate is the only object that knows the 4-channel encoding and the
payload topology. It reduces the raw payload into the per-region sufficient
statistics the inference consumes (doc 01 §1, doc 03 §2), exposing **three
parallel per-region views** — *contained*, *left*, *right* — so the joint
deconvolution applies the identical code to all three views.

Channel encoding (``ch = (spliced?2:0) + (primary?0:1)``)::

    ch0 = unspliced & genome+     ch2 = spliced & sense
    ch1 = unspliced & genome−     ch3 = spliced & antisense

The substrate exposes channels **raw** — it does not pre-orient:

* **Unspliced** counts are genome strand (``pos``/``neg``). They are oriented
  to transcript sense **downstream** (the deconvolution, using the region strand +
  the library mode). Strand-ambiguous (AMBIG) regions have no valid unspliced
  orientation and are handled by the count-clue density imputation (D7).
* **Spliced** counts are already transcript-relative ``sense``/``antisense``
  (the scanner oriented them at deposit via the splice motif), valid even in
  AMBIG regions.

``strand_class`` is carried for downstream orientation + AMBIG routing.

Boundary views use the **per-side** flux (and mass): for region ``r``, the
*left* view is the right side of ``r``'s left boundary, the *right* view is
the left side of ``r``'s right boundary (the D1 side-attribution; see
``00_implementation_plan.md`` §4 D1, and the per-side flux in PR 2.5).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays, region_boundary_indices


def _make_view(counts4: np.ndarray, mass4: np.ndarray) -> "SubstrateView":
    """Split a 4-channel (counts, mass) pair into a per-region view."""
    c = counts4.astype(np.int64, copy=False)
    m = mass4.astype(np.float64, copy=False)
    return SubstrateView(
        n_unspliced_pos=c[:, 0],
        n_unspliced_neg=c[:, 1],
        n_spliced_sense=c[:, 2],
        n_spliced_antisense=c[:, 3],
        mass_unspliced=m[:, 0] + m[:, 1],
        mass_spliced=m[:, 2] + m[:, 3],
    )


@dataclass(frozen=True, slots=True)
class SubstrateView:
    """One per-region view's sufficient statistics (all length R)."""

    n_unspliced_pos: np.ndarray  # int64[R] — genome '+' (orient downstream)
    n_unspliced_neg: np.ndarray  # int64[R] — genome '−'
    n_spliced_sense: np.ndarray  # int64[R] — motif-relative sense
    n_spliced_antisense: np.ndarray  # int64[R] — motif-relative antisense
    mass_unspliced: np.ndarray  # float64[R] — strand-agnostic magnitude
    mass_spliced: np.ndarray  # float64[R]

    @property
    def n_unspliced(self) -> np.ndarray:
        """Total unspliced count (strand-agnostic), int64[R]."""
        return self.n_unspliced_pos + self.n_unspliced_neg

    @property
    def n_spliced(self) -> np.ndarray:
        """Total spliced count (strand-agnostic), int64[R]."""
        return self.n_spliced_sense + self.n_spliced_antisense


@dataclass(frozen=True, slots=True)
class CalibrationSubstrate:
    """Per-region sufficient statistics for the three model views.

    All arrays are length R, aligned to ``RegionArrays`` order. The boundary
    views are already projected onto regions via the D1 side-attribution.
    """

    n_regions: int
    region_len: np.ndarray  # float64[R] — physical region length (bp)
    strand_class: np.ndarray  # int8[R]    — region transcript-strand class
    contained: SubstrateView
    left: SubstrateView
    right: SubstrateView

    @classmethod
    def from_payload(
        cls,
        payload: AccumulatorPayload,
        region_arrays: RegionArrays,
    ) -> "CalibrationSubstrate":
        cls._check_alignment(payload, region_arrays)

        # Contained view: integer counts; mass equals the count.
        contained = _make_view(payload.region_contained, payload.region_contained)

        # Boundary projection (D1) over PER-SIDE flux + mass. region r:
        #   left  view ← right side of its LEFT boundary  (r is that boundary's
        #               right region): flux_right / mass_right at left_boundary.
        #   right view ← left side of its RIGHT boundary (r is that boundary's
        #               left region):  flux_left  / mass_left  at right_boundary.
        left_b, right_b = region_boundary_indices(
            payload.ref_region_offsets, payload.ref_boundary_offsets
        )
        left = _make_view(payload.boundary_flux_right[left_b], payload.boundary_mass_right[left_b])
        right = _make_view(payload.boundary_flux_left[right_b], payload.boundary_mass_left[right_b])

        return cls(
            n_regions=region_arrays.n_regions,
            region_len=np.ascontiguousarray(region_arrays.region_size_bp, dtype=np.float64),
            strand_class=np.ascontiguousarray(region_arrays.strand_class, dtype=np.int8),
            contained=contained,
            left=left,
            right=right,
        )

    @staticmethod
    def _check_alignment(payload: AccumulatorPayload, region_arrays: RegionArrays) -> None:
        """Enforce that the region geometry addresses the payload 1:1 (doc 04 §4.1)."""
        if region_arrays.n_regions != payload.r_total:
            raise CalibrationSubstrateError(
                f"region geometry has {region_arrays.n_regions} regions but the payload "
                f"has {payload.r_total}; rebuild the index."
            )
        expected = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
        actual = np.asarray(payload.ref_region_offsets, dtype=np.int64)
        if not np.array_equal(expected, actual):
            raise CalibrationSubstrateError(
                "region geometry per-reference offsets do not match the payload "
                "ref_region_offsets; the accumulator region order differs from "
                "index.region_df. Rebuild the index."
            )


__all__ = ["SubstrateView", "CalibrationSubstrate"]
