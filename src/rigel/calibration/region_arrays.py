"""Sorted region geometry, the region-to-boundary index mapping, and the observable genome over both.

Three pieces of pure geometry the calibrator builds on:

* :class:`RegionArrays` — a per-reference-CSR view of the region table, sorted by ``(ref_id, start)``
  so each reference's rows are contiguous and ascending. It carries the structural columns plus the
  int8 transcript-strand class derived from each region's signature, which is the strand-model input,
  and the index's unobservable mask in the same reference ids.

* The region-to-contiguous-boundary index mapping — :func:`region_right_boundary` and
  :func:`boundary_region_indices`.

* The observable genome over the partition (`rigel.observable`): each region cut at the mask into its
  observable pieces (:func:`observable_pieces`), each boundary's bases to the nearest masked base on either
  side (:func:`boundary_reaches`), and each reference's maximal observable runs (:func:`observable_stretches`).
  The opportunity layer reads them. With no masked base they are the regions, unbounded reaches and the
  references themselves.

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

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np
import pandas as pd

from .signature import transcript_strand_class


__all__ = [
    "RegionArrays",
    "boundary_reaches",
    "boundary_region_indices",
    "observable_length",
    "observable_pieces",
    "observable_stretches",
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
    #: int64[M, 3] — the index's unobservable mask in this frame's reference ids, ``(ref_id, start, end)`` runs, merged
    #: and ascending by reference then start (`rigel.observable`); empty for a partition built by hand.
    unobservable: np.ndarray = field(default_factory=lambda: np.zeros((0, 3), dtype=np.int64))

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
        return cls.from_frame(index.regions_df, index.ref_name_to_id, index.unobservable_df)

    @classmethod
    def from_frame(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
        unobservable: pd.DataFrame | None = None,
    ) -> "RegionArrays":
        """Build from a partition frame directly — for tests that construct a partition by hand.

        Production uses :meth:`from_index`, which cannot address a different partition than the
        scanner did and carries the index's unobservable mask. *unobservable* is that mask's frame
        (``ref, start, end``); ``None`` is a partition with nothing masked."""
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

        mask = np.zeros((0, 3), dtype=np.int64)
        if unobservable is not None and len(unobservable):
            mask_ref = unobservable["ref"].astype(str).map(ref_name_to_id)
            if mask_ref.isna().any():
                unknown = sorted(unobservable.loc[mask_ref.isna(), "ref"].astype(str).unique())
                raise ValueError(
                    f"RegionArrays.from_frame: the unobservable mask names {unknown}, which are not in "
                    "ref_name_to_id. Rebuild the index."
                )
            mask = np.stack(
                [
                    mask_ref.to_numpy(np.int64),
                    unobservable["start"].to_numpy(np.int64),
                    unobservable["end"].to_numpy(np.int64),
                ],
                axis=1,
            )
            mask = mask[np.lexsort((mask[:, 1], mask[:, 0]))]

        return cls(
            ref_id=ref_id,
            start=start,
            end=end,
            signature=signature,
            strand_class=strand_class,
            region_size_bp=(end - start).astype(np.float64, copy=False),
            ref_offsets=ref_offsets,
            n_refs=n_refs,
            unobservable=mask,
        )


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


# ---------------------------------------------------------------------------
# The observable genome over the partition
# ---------------------------------------------------------------------------


def _mask_runs(region_arrays: RegionArrays, ref: int) -> tuple[np.ndarray, np.ndarray]:
    """``(starts, ends)`` — reference *ref*'s masked runs, ascending (merged, so the ends ascend too)."""
    mask = region_arrays.unobservable
    lo, hi = np.searchsorted(mask[:, 0], [ref, ref + 1])
    return mask[lo:hi, 1], mask[lo:hi, 2]


def observable_pieces(region_arrays: RegionArrays) -> tuple[np.ndarray, np.ndarray]:
    """``(piece_region, piece_len)`` int64 — every region cut at the mask into its observable pieces, in region order.

    A region with no masked base is one piece, itself; one wholly masked has none. A piece is a maximal run of
    observable bases inside one region, so a fragment contained in a region and touching no masked base lies in
    exactly one of its pieces.
    """
    ra = region_arrays
    regions, lengths = [], []
    for r in range(ra.n_refs):
        lo, hi = int(ra.ref_offsets[r]), int(ra.ref_offsets[r + 1])
        if hi <= lo:
            continue
        ms, me = _mask_runs(ra, r)
        if ms.size == 0:
            regions.append(np.arange(lo, hi, dtype=np.int64))
            lengths.append(ra.end[lo:hi] - ra.start[lo:hi])
            continue
        bounds = np.concatenate([ra.start[lo:hi], ra.end[hi - 1 : hi]])
        cuts = np.unique(np.concatenate([bounds, ms, me]))
        cuts = cuts[(cuts >= bounds[0]) & (cuts <= bounds[-1])]
        a, b = cuts[:-1], cuts[1:]
        j = np.searchsorted(ms, a, side="right") - 1
        masked = (j >= 0) & (a < me[np.maximum(j, 0)])
        keep = ~masked
        regions.append(lo + np.searchsorted(bounds, a[keep], side="right") - 1)
        lengths.append(b[keep] - a[keep])
    if not regions:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64)
    return (
        np.concatenate(regions).astype(np.int64, copy=False),
        np.concatenate(lengths).astype(np.int64, copy=False),
    )


def observable_length(region_arrays: RegionArrays) -> np.ndarray:
    """``float64[R]`` — each region's observable bases: its length less its masked bases."""
    piece_region, piece_len = observable_pieces(region_arrays)
    return np.bincount(
        piece_region, weights=piece_len.astype(np.float64), minlength=region_arrays.n_regions
    )


def boundary_reaches(region_arrays: RegionArrays) -> tuple[np.ndarray, np.ndarray]:
    """``(reach_lo, reach_hi)`` float64[E] — the bases from each contiguous boundary to the nearest masked base on its
    left and on its right, within its reference: 0 where the base beside the boundary is masked, ``inf`` where no
    masked base lies on that side, so gDNA's template runs on as it does with nothing masked."""
    ra = region_arrays
    lo, _hi = boundary_region_indices(ra.ref_id)
    x = ra.end[lo]
    reach_lo = np.full(lo.shape[0], np.inf)
    reach_hi = np.full(lo.shape[0], np.inf)
    ref = ra.ref_id[lo]
    for r in np.unique(ref):
        ms, me = _mask_runs(ra, int(r))
        if ms.size == 0:
            continue
        rows = np.flatnonzero(ref == r)
        xr = x[rows]
        # the last run that starts left of the boundary: the boundary's left base is masked iff it reaches it
        j = np.searchsorted(ms, xr, side="left") - 1
        has = j >= 0
        reach_lo[rows[has]] = xr[has] - np.minimum(me[j[has]], xr[has])
        # the first run that ends right of the boundary: the boundary's right base is masked iff it starts at or before it
        k = np.searchsorted(me, xr, side="right")
        has = k < ms.size
        reach_hi[rows[has]] = np.maximum(ms[k[has]], xr[has]) - xr[has]
    return reach_lo, reach_hi


def observable_stretches(region_arrays: RegionArrays) -> np.ndarray:
    """``int64[S]`` — the lengths of every reference's maximal observable runs: each reference cut at the mask, gDNA's
    templates. With nothing masked, the references' lengths."""
    ra = region_arrays
    out = []
    for r in range(ra.n_refs):
        lo, hi = int(ra.ref_offsets[r]), int(ra.ref_offsets[r + 1])
        if hi <= lo:
            continue
        ref_start, ref_end = int(ra.start[lo]), int(ra.end[hi - 1])
        ms, me = _mask_runs(ra, r)
        # a stretch runs from the reference's start, or a masked run's end, to the next masked run's start, or the
        # reference's end; a masked run at either end leaves an empty one there
        starts = np.concatenate([[ref_start], me])
        ends = np.concatenate([ms, [ref_end]])
        lengths = ends - starts
        out.append(lengths[lengths > 0])
    return np.concatenate(out).astype(np.int64) if out else np.zeros(0, dtype=np.int64)
