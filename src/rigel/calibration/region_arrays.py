"""rigel.calibration.region_arrays — sorted region geometry + boundary mapping.

Two pieces of pure geometry the calibrator builds on:

* :class:`RegionArrays` — a per-reference-CSR view of the region table,
  sorted by ``(ref_id, start)`` so each reference's rows are contiguous and
  ascending. Carries the structural columns plus the int8 transcript-strand
  class derived from each region's signature (the D2 strand-model input).

* The **boundary↔region index mapping** — for each reference with ``k``
  regions there are ``k + 1`` boundary slots (the two terminals plus the
  ``k - 1`` internal seams). These functions map a region to the two
  boundaries flanking it, and a boundary to the regions on either side. They
  are computed purely from the accumulator payload's topology offsets
  (``ref_region_offsets`` / ``ref_boundary_offsets``), so they hold whenever
  the region arrays are aligned 1:1 with the payload.

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
    "region_boundary_indices",
    "boundary_region_indices",
]


@dataclass(frozen=True, slots=True)
class RegionArrays:
    """Per-reference-CSR view of the region table.

    Rows are sorted by ``(ref_id, start)`` so each ref's rows are contiguous
    and ascending. ``ref_offsets`` is the per-ref CSR boundary array; the
    slice ``[ref_offsets[r]:ref_offsets[r + 1]]`` is one reference's regions,
    sorted by ``start`` (and by ``end`` — regions are non-overlapping within a
    reference). ``order`` is the permutation applied to ``region_df`` rows.
    """

    ref_id: np.ndarray  # int32, (R,)
    start: np.ndarray  # int64, (R,)
    end: np.ndarray  # int64, (R,)
    signature: np.ndarray  # uint8, (R,)
    strand_class: np.ndarray  # int8,  (R,) — TS_NONE/TS_POS/TS_NEG/TS_AMBIG
    region_size_bp: np.ndarray  # float64, (R,) — end - start
    ref_offsets: np.ndarray  # int32, (n_refs + 1,)
    order: np.ndarray  # int64, (R,)
    n_refs: int

    @property
    def n_regions(self) -> int:
        return int(self.start.shape[0])

    @classmethod
    def from_region_df(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
    ) -> "RegionArrays":
        if "signature" not in region_df.columns:
            raise ValueError(
                "RegionArrays.from_region_df: region_df is missing the 'signature' "
                "column. Rebuild the index against the calibration-v6 schema."
            )
        n_refs = len(ref_name_to_id)
        n_regions = len(region_df)

        ref_id = region_df["ref_name"].map(ref_name_to_id).to_numpy()
        if np.any(pd.isna(ref_id)):
            unknown = region_df.loc[pd.isna(ref_id), "ref_name"].unique().tolist()
            raise ValueError(
                f"RegionArrays.from_region_df: region_df references {sorted(unknown)} "
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
                "RegionArrays.from_region_df: ref_offsets sum mismatch."
            )

        return cls(
            ref_id=ref_id,
            start=start,
            end=end,
            signature=signature,
            strand_class=strand_class,
            region_size_bp=(end - start).astype(np.float64, copy=False),
            ref_offsets=ref_offsets,
            order=order,
            n_refs=n_refs,
        )


# ---------------------------------------------------------------------------
# Boundary ↔ region index mapping
# ---------------------------------------------------------------------------


def _as_offsets(arr: np.ndarray, name: str) -> np.ndarray:
    offsets = np.asarray(arr, dtype=np.int64)
    if offsets.ndim != 1 or offsets.shape[0] < 1:
        raise ValueError(f"{name} must be a 1-D offset array of length n_refs + 1.")
    return offsets


def region_boundary_indices(
    ref_region_offsets: np.ndarray,
    ref_boundary_offsets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map each region to the boundaries flanking it.

    For region ``r`` in reference ``f`` with local index ``i``::

        left_boundary[r]  = ref_boundary_offsets[f] + i
        right_boundary[r] = ref_boundary_offsets[f] + i + 1

    Returns ``(left_boundary, right_boundary)``, each int64 of length
    ``R = ref_region_offsets[-1]``. Every region has both indices defined
    (the first/last region's outer boundary is its reference terminal).
    """
    rro = _as_offsets(ref_region_offsets, "ref_region_offsets")
    rbo = _as_offsets(ref_boundary_offsets, "ref_boundary_offsets")
    if rro.shape != rbo.shape:
        raise ValueError("ref_region_offsets and ref_boundary_offsets must share length.")

    n_refs = rro.shape[0] - 1
    r_total = int(rro[-1])
    region_ref = np.repeat(np.arange(n_refs, dtype=np.int64), np.diff(rro))
    local = np.arange(r_total, dtype=np.int64) - rro[region_ref]
    left_boundary = rbo[region_ref] + local
    return left_boundary, left_boundary + 1


def boundary_region_indices(
    ref_region_offsets: np.ndarray,
    ref_boundary_offsets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map each boundary to the regions on either side.

    For boundary slot ``j`` (0-based within its reference of ``k`` regions),
    ``left_region`` is the region to its left and ``right_region`` the region
    to its right, with ``-1`` on the off-edge side of a reference terminal
    (``j == 0`` has no left region; ``j == k`` has no right region).

    Returns ``(left_region, right_region)``, each int64 of length
    ``B = ref_boundary_offsets[-1]``. This is the exact inverse of
    :func:`region_boundary_indices` on internal seams and the natural
    one-sided attribution at terminals.
    """
    rro = _as_offsets(ref_region_offsets, "ref_region_offsets")
    rbo = _as_offsets(ref_boundary_offsets, "ref_boundary_offsets")
    if rro.shape != rbo.shape:
        raise ValueError("ref_region_offsets and ref_boundary_offsets must share length.")

    n_refs = rbo.shape[0] - 1
    b_total = int(rbo[-1])
    boundary_ref = np.repeat(np.arange(n_refs, dtype=np.int64), np.diff(rbo))
    local_b = np.arange(b_total, dtype=np.int64) - rbo[boundary_ref]
    regions_in_ref = np.diff(rro)[boundary_ref]
    region_base = rro[boundary_ref]
    left_region = np.where(local_b >= 1, region_base + local_b - 1, -1)
    right_region = np.where(local_b <= regions_in_ref - 1, region_base + local_b, -1)
    return left_region, right_region
