"""rigel.calibration._arrays \u2014 sorted region + payload views.

Pre-flattens region columns and the fractional ``region_counts`` matrix
into per-reference CSR sorted order so the locus-prior assembly pass
(Phase 4) can index both array families with one sorted-position index.

``RegionArrays`` carries the structural columns plus the per-region
transcript-strand class. ``PayloadArrays`` is a thin sorted view of the
12-channel ``region_counts`` matrix plus the unspliced channel
projections used by hot consumers.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from .fractional_evidence import transcript_strand_class
from .scan_payload import CalibrationScanPayload
from .signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    SPLICE_UNSPLICED,
    channel_index,
)


__all__ = ["RegionArrays", "PayloadArrays"]


@dataclass(frozen=True, slots=True)
class RegionArrays:
    """Pre-flattened, per-reference-CSR view of the region table.

    Rows are sorted by ``(ref_id, start)`` so each ref's rows are
    contiguous and ascending. ``ref_offsets`` is the per-ref CSR
    boundary array; the slice ``[ref_offsets[r]:ref_offsets[r+1]]`` is
    sorted by ``start`` (and by ``end`` too \u2014 regions are
    non-overlapping within a reference).
    """

    ref_id: np.ndarray  # int32, (R,)
    start: np.ndarray  # int64, (R,)
    end: np.ndarray  # int64, (R,)
    signature: np.ndarray  # uint8, (R,) \u2014 fine-region signature
    ts_class: np.ndarray  # int8,  (R,) \u2014 TS_NONE/TS_POS/TS_NEG/TS_AMBIG
    region_size_bp: np.ndarray  # float64, (R,) - end - start in bp (PR 03)
    ref_offsets: np.ndarray  # int32, (n_refs + 1,)
    order: np.ndarray  # int64, (R,)
    n_refs: int

    @classmethod
    def from_region_df(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
    ) -> "RegionArrays":
        if "signature" not in region_df.columns:
            raise ValueError(
                "RegionArrays.from_region_df: region_df is missing the "
                "'signature' column. Rebuild the index against the "
                "fractional accumulator schema."
            )
        n_refs = len(ref_name_to_id)
        n_regions = len(region_df)

        ref_id = region_df["ref_name"].map(ref_name_to_id).to_numpy()
        if np.any(pd.isna(ref_id)):
            unknown = region_df.loc[pd.isna(ref_id), "ref_name"].unique().tolist()
            raise ValueError(
                f"RegionArrays.from_region_df: region_df references "
                f"{sorted(unknown)} which are not in ref_name_to_id. "
                f"Rebuild the index."
            )
        ref_id = ref_id.astype(np.int32, copy=False)

        order = np.lexsort((region_df["start"].to_numpy(), ref_id))
        ref_id = ref_id[order]
        start = region_df["start"].to_numpy().astype(np.int64, copy=False)[order]
        end = region_df["end"].to_numpy().astype(np.int64, copy=False)[order]
        signature = region_df["signature"].to_numpy().astype(np.uint8, copy=False)[order]
        ts_class = transcript_strand_class(signature)

        counts = np.bincount(ref_id, minlength=n_refs).astype(np.int32, copy=False)
        ref_offsets = np.empty(n_refs + 1, dtype=np.int32)
        ref_offsets[0] = 0
        np.cumsum(counts, out=ref_offsets[1:])
        if int(ref_offsets[-1]) != n_regions:
            raise RuntimeError(  # pragma: no cover \u2014 invariant guard
                "RegionArrays.from_region_df: ref_offsets sum mismatch."
            )

        return cls(
            ref_id=ref_id,
            start=start,
            end=end,
            signature=signature,
            ts_class=ts_class,
            region_size_bp=(end - start).astype(np.float64, copy=False),
            ref_offsets=ref_offsets,
            order=order,
            n_refs=n_refs,
        )


@dataclass(frozen=True, slots=True)
class PayloadArrays:
    """Sorted view of the fractional ``region_counts`` matrix.

    Rows are reordered by ``RegionArrays.order`` so a single sorted
    position index can address both array families. The unspliced
    per-channel projections are materialised once for hot consumers.
    """

    region_counts_sorted: np.ndarray  # float32[R, 12]

    # Per-region physical fragment support, sorted by RegionArrays.order.
    region_unspliced_support_sorted: np.ndarray  # uint64[R]
    region_spliced_support_sorted: np.ndarray  # uint64[R]

    contained_unspliced_pos: np.ndarray  # float32[R]
    contained_unspliced_neg: np.ndarray  # float32[R]
    boundary_left_unspliced_pos: np.ndarray  # float32[R]
    boundary_left_unspliced_neg: np.ndarray  # float32[R]
    boundary_right_unspliced_pos: np.ndarray  # float32[R]
    boundary_right_unspliced_neg: np.ndarray  # float32[R]

    @classmethod
    def from_payload(
        cls,
        payload: CalibrationScanPayload,
        region_arrays: RegionArrays,
    ) -> "PayloadArrays":
        order = region_arrays.order
        rc = np.ascontiguousarray(payload.region_counts[order, :])
        u_sup = np.ascontiguousarray(payload.region_unspliced_support[order])
        s_sup = np.ascontiguousarray(payload.region_spliced_support[order])

        c_u_p = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        c_u_n = channel_index(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
        bl_u_p = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        bl_u_n = channel_index(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)
        br_u_p = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)
        br_u_n = channel_index(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)

        contained_pos = rc[:, c_u_p]
        contained_neg = rc[:, c_u_n]
        bl_pos = rc[:, bl_u_p]
        bl_neg = rc[:, bl_u_n]
        br_pos = rc[:, br_u_p]
        br_neg = rc[:, br_u_n]

        return cls(
            region_counts_sorted=rc,
            region_unspliced_support_sorted=u_sup,
            region_spliced_support_sorted=s_sup,
            contained_unspliced_pos=contained_pos,
            contained_unspliced_neg=contained_neg,
            boundary_left_unspliced_pos=bl_pos,
            boundary_left_unspliced_neg=bl_neg,
            boundary_right_unspliced_pos=br_pos,
            boundary_right_unspliced_neg=br_neg,
        )
