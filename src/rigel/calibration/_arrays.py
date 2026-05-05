"""rigel.calibration._arrays — Pre-flattened views of region_df + payload.

The locus-prior assembly pass touches the region table once per locus.
Materializing ``region_df["start"].to_numpy()`` per call dominates
runtime on indexes with millions of regions, so we flatten the five
hot columns once at :func:`assemble_priors` entry and pass the
resulting frozen views around.

See ``docs/calibration/m6_implementation_plan.md`` §3.2.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
import pandas as pd

from .scan_payload import (
    MASK_INTERGENIC,
    MASK_INTRON,
    CalibrationScanPayload,
)


__all__ = ["RegionArrays", "PayloadArrays"]


@dataclass(frozen=True, slots=True)
class RegionArrays:
    """Pre-flattened, per-reference-CSR view of the region table.

    Rows are sorted by ``(ref_id, start)`` so each ref's rows are
    contiguous and ascending.  ``ref_offsets`` is the per-ref CSR
    boundary array; the slice ``[ref_offsets[r]:ref_offsets[r+1]]``
    is sorted by ``start`` (and, by the M2 invariant, by ``end`` too —
    regions are non-overlapping within a reference).

    The ``order`` attribute is the permutation taking the original
    ``region_df`` row order to this sorted order, so payload columns
    (which live in original ``region_id`` order) can be reordered
    consistently via :meth:`PayloadArrays.from_payload`.
    """

    ref_id: np.ndarray        # int32, (R,)
    start: np.ndarray         # int64, (R,)
    end: np.ndarray           # int64, (R,)
    type: np.ndarray          # uint8, (R,)
    bf_left: np.ndarray       # bool,  (R,)
    bf_right: np.ndarray      # bool,  (R,)
    ref_offsets: np.ndarray   # int32, (n_refs + 1,)
    order: np.ndarray         # int64, (R,)  — region_df_row_order[i] in sorted slot i
    n_refs: int

    @classmethod
    def from_region_df(
        cls,
        region_df: pd.DataFrame,
        ref_name_to_id: Mapping[str, int],
    ) -> "RegionArrays":
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

        # Sort by (ref_id, start) so each ref's rows are contiguous + sorted.
        order = np.lexsort((region_df["start"].to_numpy(), ref_id))
        ref_id = ref_id[order]
        start = region_df["start"].to_numpy().astype(np.int64, copy=False)[order]
        end = region_df["end"].to_numpy().astype(np.int64, copy=False)[order]
        type_ = region_df["type"].to_numpy().astype(np.uint8, copy=False)[order]
        bf_left = region_df["boundary_flux_left"].to_numpy().astype(bool, copy=False)[order]
        bf_right = region_df["boundary_flux_right"].to_numpy().astype(bool, copy=False)[order]

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
            type=type_,
            bf_left=bf_left,
            bf_right=bf_right,
            ref_offsets=ref_offsets,
            order=order,
            n_refs=n_refs,
        )


@dataclass(frozen=True, slots=True)
class PayloadArrays:
    """Pre-flattened per-region count / flux columns, sorted to match
    :class:`RegionArrays`.

    The four columns are reordered via ``RegionArrays.order`` so that
    a single sorted-position index from :class:`RegionIndexPy` can
    safely index into both array families.
    """

    intergenic_per_region: np.ndarray  # int64, (R,)
    intron_per_region: np.ndarray      # int64, (R,)
    u_left: np.ndarray                 # int64, (R,)
    u_right: np.ndarray                # int64, (R,)

    @classmethod
    def from_payload(
        cls,
        payload: CalibrationScanPayload,
        region_arrays: RegionArrays,
    ) -> "PayloadArrays":
        order = region_arrays.order
        return cls(
            intergenic_per_region=np.ascontiguousarray(
                payload.per_region_counts[:, MASK_INTERGENIC][order]
            ),
            intron_per_region=np.ascontiguousarray(
                payload.per_region_counts[:, MASK_INTRON][order]
            ),
            u_left=np.ascontiguousarray(payload.u_left[order]),
            u_right=np.ascontiguousarray(payload.u_right[order]),
        )
