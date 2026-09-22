"""Sorted region geometry, and the region-to-boundary index mapping.

Two pieces of pure geometry the calibrator builds on:

* :class:`RegionArrays` — a per-reference-CSR view of the region table, sorted by ``(ref_id, start)``
  so each reference's rows are contiguous and ascending. It carries the structural columns plus the
  int8 transcript-strand class derived from each region's signature, which is the strand-model input.

* The region-to-contiguous-boundary index mapping — :func:`region_right_boundary` and
  :func:`boundary_region_indices`.

A contiguous boundary is the boundary BETWEEN two adjacent regions, so a reference with ``k`` regions
owns exactly ``k - 1`` of them and a boundary always has a region on both sides. There are no terminal
boundary slots, and the boundary axis is therefore derivable from ``ref_id`` alone — no stored offsets.

That derivation rests on one fact: boundary ids are assigned per reference in genomic order, in
reference order, which is exactly the order adjacent same-reference region pairs appear in a
``(ref_id, start)``-sorted region table. :func:`~rigel.calibration.region_chain.build_region_chain`
lays out the same numbering by walking the payload's CSR offsets, and
``test_boundary_numbering_matches_the_chain_built_from_the_payload_offsets`` pins the two against
each other — a second algorithm, not a second call to the first.

No tunable parameters: this module is index arithmetic only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from .signature import transcript_strand_class


__all__ = [
    "RegionArrays",
    "boundary_region_indices",
    "overlapping_region_runs",
    "region_right_boundary",
]


@dataclass(frozen=True, slots=True)
class RegionArrays:
    """Per-reference-CSR view of the region table.

    Rows are sorted by ``(ref_id, start)`` so each ref's rows are contiguous
    and ascending. ``ref_offsets`` is the per-ref CSR boundary array; the
    slice ``[ref_offsets[r]:ref_offsets[r + 1]]`` is one reference's regions, sorted
    by ``start`` (and by ``end`` — regions tile a reference).
    """

    ref_id: np.ndarray  # int32, (R,)
    start: np.ndarray  # int64, (R,)
    end: np.ndarray  # int64, (R,)
    signature: np.ndarray  # uint8, (R,)
    strand_class: np.ndarray  # int8,  (R,) — TS_NONE/TS_POS/TS_NEG/TS_AMBIG
    region_size_bp: np.ndarray  # float64, (R,) — end - start
    ref_offsets: np.ndarray  # int32, (n_refs + 1,)
    n_refs: int

    @property
    def n_regions(self) -> int:
        return int(self.start.shape[0])

    @classmethod
    def from_index(cls, index) -> "RegionArrays":
        """Build the geometry for the partition the scanner actually deposits into.

        This and :func:`~rigel.calibration.splice_graph.build_region_partition_arrays` are the two
        halves of one contract — the calibration geometry must address the payload the scanner
        produced — so they read the same frame (``index.regions_df``, the splice graph) through
        one accessor. Passing a frame by hand is how the two drift apart, and nothing downstream
        detects that except as a shape error far from its cause.
        """
        return cls.from_frame(index.regions_df, index.ref_name_to_id)

    @classmethod
    def from_frame(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
    ) -> "RegionArrays":
        """Build from a partition frame directly — for tests that construct a partition by hand.

        Production uses :meth:`from_index`, which cannot address a different partition than the
        scanner did."""
        if "signature" not in region_df.columns:
            raise ValueError(
                "RegionArrays.from_frame: region_df is missing the 'signature' "
                "column. Rebuild the index against the current schema."
            )
        n_refs = len(ref_name_to_id)
        n_regions = len(region_df)

        ref_id = region_df["ref_name"].map(ref_name_to_id).to_numpy()
        if np.any(pd.isna(ref_id)):
            unknown = region_df.loc[pd.isna(ref_id), "ref_name"].unique().tolist()
            raise ValueError(
                f"RegionArrays.from_frame: region_df references {sorted(unknown)} "
                f"which are not in ref_name_to_id. Rebuild the index."
            )
        ref_id = ref_id.astype(np.int32, copy=False)

        order = np.lexsort((region_df["start"].to_numpy(), ref_id))
        ref_id = ref_id[order]
        start = region_df["start"].to_numpy().astype(np.int64, copy=False)[order]
        end = region_df["end"].to_numpy().astype(np.int64, copy=False)[order]
        signature = region_df["signature"].to_numpy().astype(np.uint8, copy=False)[order]
        strand_class = transcript_strand_class(signature)

        counts = np.bincount(ref_id, minlength=n_refs).astype(np.int32, copy=False)
        ref_offsets = np.empty(n_refs + 1, dtype=np.int32)
        ref_offsets[0] = 0
        np.cumsum(counts, out=ref_offsets[1:])
        if int(ref_offsets[-1]) != n_regions:
            raise RuntimeError(  # pragma: no cover — invariant guard
                "RegionArrays.from_frame: ref_offsets sum mismatch."
            )

        return cls(
            ref_id=ref_id,
            start=start,
            end=end,
            signature=signature,
            strand_class=strand_class,
            region_size_bp=(end - start).astype(np.float64, copy=False),
            ref_offsets=ref_offsets,
            n_refs=n_refs,
        )


# The index schema carries no per-boundary annotation flags: the solver reads an sj's strand from the
# accumulator's splice motif instead.


# ---------------------------------------------------------------------------
# Region ↔ contiguous-boundary index mapping
# ---------------------------------------------------------------------------


def region_right_boundary(ref_id: np.ndarray) -> np.ndarray:
    """``int64[N]`` — the contiguous boundary right of each region, ``-1`` at a reference's last region.

    Regions ``r`` and ``r + 1`` share a boundary exactly when they are in the same reference, so the
    boundary axis is the run of adjacent same-reference pairs, numbered in region order. A reference
    with one region owns no boundary; an empty reference contributes nothing.

    ``ref_id`` must be grouped — all of a reference's regions contiguous — which
    :class:`RegionArrays` guarantees by sorting on ``(ref_id, start)``. Ungrouped input would
    manufacture boundaries that straddle references, so it is refused rather than tolerated.
    """
    ref = np.asarray(ref_id)
    n = int(ref.shape[0])
    out = np.full(n, -1, dtype=np.int64)
    if n < 2:
        return out
    same = ref[:-1] == ref[1:]
    # grouped ⇔ each reference's rows form ONE run ⇔ the number of runs equals the number of
    # distinct references. `np.unique` counts distinct values; `same` counts run breaks.
    if int((~same).sum()) + 1 != int(np.unique(ref).shape[0]):
        raise ValueError(
            "ref_id is not grouped: each reference's regions must be contiguous. Build the geometry "
            "with RegionArrays.from_index / from_frame, which sorts on (ref_id, start)."
        )
    out[:-1][same] = np.arange(int(same.sum()), dtype=np.int64)
    return out


def boundary_region_indices(ref_id: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(lo_region, hi_region)`` — the two regions contiguous boundary ``e`` lies between, each ``int64[E]``.

    The exact inverse of :func:`region_right_boundary`. ``hi_region == lo_region + 1`` always, and
    both are in the same reference by construction — there is no terminal case and no ``-1``.
    """
    right = region_right_boundary(ref_id)
    lo = np.flatnonzero(right >= 0).astype(np.int64)
    return lo, lo + 1


def overlapping_region_runs(
    ref_id: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    region_start: np.ndarray,
    region_end: np.ndarray,
    ref_offsets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """``(lo, hi)`` per interval — the run of partition regions ``[lo, hi)`` that ``[start, end)`` overlaps.

    The regions tile each reference in order, so an interval overlaps one contiguous run of them: the
    first region ending after ``start`` through the last starting before ``end``. ``hi <= lo`` means the
    interval overlaps nothing. An interval whose ``ref_id`` is negative (a reference the partition does
    not carry) returns ``lo = hi = 0``. One ``searchsorted`` per reference present, not per interval.
    """
    ref_id = np.asarray(ref_id, dtype=np.int64)
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    region_start = np.asarray(region_start, dtype=np.int64)
    region_end = np.asarray(region_end, dtype=np.int64)
    ref_offsets = np.asarray(ref_offsets, dtype=np.int64)
    lo = np.zeros(ref_id.shape[0], dtype=np.int64)
    hi = np.zeros(ref_id.shape[0], dtype=np.int64)
    for r in np.unique(ref_id[ref_id >= 0]):
        rows = np.flatnonzero(ref_id == r)
        lo0, hi0 = int(ref_offsets[r]), int(ref_offsets[r + 1])
        lo[rows] = lo0 + np.searchsorted(region_end[lo0:hi0], start[rows], side="right")
        hi[rows] = lo0 + np.searchsorted(region_start[lo0:hi0], end[rows], side="left")
    return lo, hi
