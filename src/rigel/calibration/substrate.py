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
the per-side flux in PR 2.5).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..scan_payload import AccumulatorPayload
from .errors import CalibrationSubstrateError
from .region_arrays import RegionArrays, boundary_region_indices, region_boundary_indices


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


@dataclass(frozen=True, slots=True)
class BoundarySubstrate:
    """Per-**boundary** sufficient statistics — the boundary nodes of the bipartite graph (PLAN v6 §2).

    The calibration graph is a linear bipartite chain ``R ↔ B ↔ R`` of region and **boundary** nodes; this is
    the boundary-indexed view of the payload (the twin of :class:`CalibrationSubstrate`, which is
    region-indexed). A boundary node ``b`` is one genomic point with two sides:

    * ``left``  — the side lying **inside** ``left_region[b]``  (the payload's ``boundary_*_left[b]``);
    * ``right`` — the side lying **inside** ``right_region[b]`` (the payload's ``boundary_*_right[b]``).

    A contiguous unspliced crossing credits **both** sides (it straddles the point); a spliced intron-skip
    credits **one** side (the exon flank) — so the boundary node owns the one-sided, motif-stranded spliced
    mass (the mature-RNA floor) and the two-sided unspliced crossing mass (its ``{f₊,f₋,f_g}`` pie). The two
    sides are mass-accounting, not two beliefs: the solver pools them into one belief and applies the result to
    each side's mass. ``left_region`` / ``right_region`` are ``-1`` on the off-edge side of a reference
    terminal (a terminal boundary has only one flanking region).

    This is the **same information** as ``CalibrationSubstrate.left`` / ``.right`` (which project the boundary
    sides onto regions), re-keyed by boundary: for an internal boundary ``b`` with
    ``(left_region=lr, right_region=rr)``, ``BoundarySubstrate.left[b] == CalibrationSubstrate.right[lr]`` and
    ``BoundarySubstrate.right[b] == CalibrationSubstrate.left[rr]`` (the re-indexing identity, asserted in the
    tests). No payload reshape, no C++ change — the boundary objects already exist.
    """

    n_boundaries: int
    left_region: (
        np.ndarray
    )  # int64[B] — region to the boundary's left; -1 at a reference-start terminal
    right_region: (
        np.ndarray
    )  # int64[B] — region to the boundary's right; -1 at a reference-end terminal
    left: SubstrateView  # the boundary's LEFT side (inside left_region): boundary_*_left[b]
    right: SubstrateView  # the boundary's RIGHT side (inside right_region): boundary_*_right[b]
    #: int8[B] — the boundary's splice-junction genomic strand (``Strand``: POS=1 /
    #: NEG=2, 0 = no junction). The mature-RNA anchor's strand, observed from the
    #: motif at deposit (one junction per boundary). Routes the one-sided spliced
    #: mass to its exon flank's strand — correct even at AMBIG / exon↔exon seams.
    junction_strand: np.ndarray

    @classmethod
    def from_payload(cls, payload: AccumulatorPayload) -> "BoundarySubstrate":
        # Boundary sides: counts = integer flux, mass = fractional crossing mass (the same _make_view
        # contract as the region-contained view, where counts == mass).
        left = _make_view(payload.boundary_flux_left, payload.boundary_mass_left)
        right = _make_view(payload.boundary_flux_right, payload.boundary_mass_right)
        left_region, right_region = boundary_region_indices(
            payload.ref_region_offsets, payload.ref_boundary_offsets
        )
        bjs = payload.boundary_junction_strand
        if bjs is None:  # stale payload predating the field → no junction strand
            bjs = np.zeros(payload.b_obj_total, dtype=np.int8)
        return cls(
            n_boundaries=payload.b_obj_total,
            left_region=left_region,
            right_region=right_region,
            left=left,
            right=right,
            junction_strand=np.ascontiguousarray(bjs, dtype=np.int8),
        )


__all__ = ["SubstrateView", "CalibrationSubstrate", "BoundarySubstrate"]
