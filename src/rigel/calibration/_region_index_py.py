"""rigel.calibration._region_index_py — Per-reference region overlap query.

Pure-Python overlap index built on top of :class:`RegionArrays`.  Two
``np.searchsorted`` calls per query, exploiting the M2 invariant that
regions are non-overlapping within a reference (so per-ref ``end``
arrays are sorted ascending too).

See ``docs/calibration/m6_implementation_plan.md`` §3.3.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._arrays import RegionArrays


__all__ = ["RegionIndexPy"]


@dataclass(frozen=True, slots=True)
class RegionIndexPy:
    arrays: RegionArrays

    def overlap(self, ref_id: int, qstart: int, qend: int) -> np.ndarray:
        """Return sorted-position indices of regions overlapping
        ``[qstart, qend)`` on ``ref_id``.

        Returned indices are in :class:`RegionArrays` sorted-position
        space, NOT original ``region_id`` space.  They index correctly
        into both :class:`RegionArrays` and (post-reorder)
        :class:`PayloadArrays`.

        Returns an empty ``int64[]`` if ``ref_id`` is out of range or
        the query interval is empty / non-overlapping.
        """
        a = self.arrays
        if ref_id < 0 or ref_id >= a.n_refs:
            return np.empty(0, dtype=np.int64)
        lo = int(a.ref_offsets[ref_id])
        hi = int(a.ref_offsets[ref_id + 1])
        if lo == hi or qend <= qstart:
            return np.empty(0, dtype=np.int64)

        starts = a.start[lo:hi]
        ends = a.end[lo:hi]
        # First region with start >= qend : exclusive upper bound.
        j_hi = lo + int(np.searchsorted(starts, qend, side="left"))
        # First region with end > qstart : inclusive lower bound.
        # Per-ref `end` is sorted because regions are non-overlapping.
        j_lo = lo + int(np.searchsorted(ends, qstart, side="right"))
        if j_lo >= j_hi:
            return np.empty(0, dtype=np.int64)
        return np.arange(j_lo, j_hi, dtype=np.int64)
