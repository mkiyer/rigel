"""Region-boundary table derived from the fractional region ledger.

The scanner currently stores boundary-crossing mass on each region side:
``boundary_left`` and ``boundary_right``.  The v6 calibration model wants a
biologically direct view where a boundary between two adjacent regions is the
primary object.  ``BoundaryTable`` is that view.  It is derived from existing
arrays and does not require scanner changes.
"""

from __future__ import annotations

from dataclasses import dataclass

import math

import numpy as np

from ._arrays import RegionArrays
from .region_count_ledger import RegionCountLedger


__all__ = [
    "BoundaryTable",
    "build_boundary_table",
    "validate_boundary_table",
]


@dataclass(frozen=True, slots=True)
class BoundaryTable:
    """One-dimensional boundary table over sorted regions.

    For each reference with ``N`` regions, there are ``N + 1`` boundary slots.
    Boundary slot 0 is the terminal boundary before the first region, slot N is
    the terminal boundary after the last region, and internal slot ``i`` is the
    boundary between local regions ``i - 1`` and ``i``.
    """

    boundary_pos: np.ndarray  # int64[B]
    ref_id: np.ndarray  # int32[B]
    ref_region_offsets: np.ndarray  # int64[n_refs + 1]
    ref_boundary_offsets: np.ndarray  # int64[n_refs + 1]
    is_terminal: np.ndarray  # bool[B]

    left_region_unspliced_pos: np.ndarray  # float32[B]
    left_region_unspliced_neg: np.ndarray
    left_region_spliced_pos: np.ndarray
    left_region_spliced_neg: np.ndarray

    right_region_unspliced_pos: np.ndarray  # float32[B]
    right_region_unspliced_neg: np.ndarray
    right_region_spliced_pos: np.ndarray
    right_region_spliced_neg: np.ndarray

    left_region_boundary_leff: np.ndarray  # float64[B]
    right_region_boundary_leff: np.ndarray  # float64[B]

    @property
    def n_boundaries(self) -> int:
        return int(self.boundary_pos.shape[0])

    def left_region_index(self) -> np.ndarray:
        """Return the region to the left of each boundary, or -1 at terminals."""
        out = np.full(self.n_boundaries, -1, dtype=np.int64)
        for ref_idx in range(int(self.ref_region_offsets.shape[0] - 1)):
            region_start = int(self.ref_region_offsets[ref_idx])
            region_end = int(self.ref_region_offsets[ref_idx + 1])
            boundary_start = int(self.ref_boundary_offsets[ref_idx])
            if region_end - region_start <= 1:
                continue
            count = region_end - region_start
            internal = boundary_start + np.arange(1, count, dtype=np.int64)
            out[internal] = region_start + np.arange(count - 1, dtype=np.int64)
        return out

    def right_region_index(self) -> np.ndarray:
        """Return the region to the right of each boundary, or -1 at terminals."""
        out = np.full(self.n_boundaries, -1, dtype=np.int64)
        for ref_idx in range(int(self.ref_region_offsets.shape[0] - 1)):
            region_start = int(self.ref_region_offsets[ref_idx])
            region_end = int(self.ref_region_offsets[ref_idx + 1])
            boundary_start = int(self.ref_boundary_offsets[ref_idx])
            if region_end - region_start <= 1:
                continue
            count = region_end - region_start
            internal = boundary_start + np.arange(1, count, dtype=np.int64)
            out[internal] = region_start + np.arange(1, count, dtype=np.int64)
        return out


def _validate_inputs(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    side_leff: np.ndarray,
) -> np.ndarray:
    region_count = int(region_arrays.ref_id.shape[0])
    side = np.asarray(side_leff, dtype=np.float64)
    if side.shape != (region_count,):
        raise ValueError(
            "build_boundary_table: side_leff must have shape (n_regions,); "
            f"got {side.shape}, expected {(region_count,)}."
        )
    for name in (
        "boundary_left_unspliced_pos",
        "boundary_left_unspliced_neg",
        "boundary_right_unspliced_pos",
        "boundary_right_unspliced_neg",
        "boundary_left_spliced_pos",
        "boundary_left_spliced_neg",
        "boundary_right_spliced_pos",
        "boundary_right_spliced_neg",
    ):
        arr = getattr(ledger, name)
        if arr.shape != (region_count,):
            raise ValueError(
                f"build_boundary_table: ledger.{name} has shape {arr.shape}; "
                f"expected {(region_count,)}."
            )
    return side


def build_boundary_table(
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
    side_leff: np.ndarray,
) -> BoundaryTable:
    """Derive a boundary table from sorted regions and current side slots."""
    side = _validate_inputs(region_arrays, ledger, side_leff)
    region_count = int(region_arrays.ref_id.shape[0])
    ref_region_offsets = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    n_refs = int(region_arrays.n_refs)
    if ref_region_offsets.shape != (n_refs + 1,):
        raise ValueError(
            "build_boundary_table: region_arrays.ref_offsets has unexpected shape "
            f"{ref_region_offsets.shape}; expected {(n_refs + 1,)}."
        )

    regions_per_ref = np.diff(ref_region_offsets)
    boundaries_per_ref = regions_per_ref + 1
    ref_boundary_offsets = np.empty(n_refs + 1, dtype=np.int64)
    ref_boundary_offsets[0] = 0
    np.cumsum(boundaries_per_ref, out=ref_boundary_offsets[1:])
    boundary_count = int(ref_boundary_offsets[-1])
    if boundary_count != region_count + n_refs:
        raise RuntimeError("build_boundary_table: boundary count invariant failed.")

    boundary_pos = np.zeros(boundary_count, dtype=np.int64)
    boundary_ref_id = np.repeat(np.arange(n_refs, dtype=np.int32), boundaries_per_ref)
    is_terminal = np.zeros(boundary_count, dtype=bool)

    def zeros32() -> np.ndarray:
        return np.zeros(boundary_count, dtype=np.float32)

    left_region_unspliced_pos = zeros32()
    left_region_unspliced_neg = zeros32()
    left_region_spliced_pos = zeros32()
    left_region_spliced_neg = zeros32()
    right_region_unspliced_pos = zeros32()
    right_region_unspliced_neg = zeros32()
    right_region_spliced_pos = zeros32()
    right_region_spliced_neg = zeros32()
    left_region_boundary_leff = np.zeros(boundary_count, dtype=np.float64)
    right_region_boundary_leff = np.zeros(boundary_count, dtype=np.float64)

    region_indices = np.arange(region_count, dtype=np.int64)
    region_ref = np.asarray(region_arrays.ref_id, dtype=np.int64)
    local_region = region_indices - ref_region_offsets[region_ref]
    left_boundary = ref_boundary_offsets[region_ref] + local_region
    right_boundary = left_boundary + 1
    internal_left_side = local_region > 0
    internal_right_side = local_region < (regions_per_ref[region_ref] - 1)

    # A region's right side lands in the left-region fields of its right boundary.
    left_region_unspliced_pos[right_boundary[internal_right_side]] = (
        ledger.boundary_right_unspliced_pos[internal_right_side]
    )
    left_region_unspliced_neg[right_boundary[internal_right_side]] = (
        ledger.boundary_right_unspliced_neg[internal_right_side]
    )
    left_region_spliced_pos[right_boundary[internal_right_side]] = (
        ledger.boundary_right_spliced_pos[internal_right_side]
    )
    left_region_spliced_neg[right_boundary[internal_right_side]] = (
        ledger.boundary_right_spliced_neg[internal_right_side]
    )
    left_region_boundary_leff[right_boundary[internal_right_side]] = side[internal_right_side]

    # A region's left side lands in the right-region fields of its left boundary.
    right_region_unspliced_pos[left_boundary[internal_left_side]] = (
        ledger.boundary_left_unspliced_pos[internal_left_side]
    )
    right_region_unspliced_neg[left_boundary[internal_left_side]] = (
        ledger.boundary_left_unspliced_neg[internal_left_side]
    )
    right_region_spliced_pos[left_boundary[internal_left_side]] = (
        ledger.boundary_left_spliced_pos[internal_left_side]
    )
    right_region_spliced_neg[left_boundary[internal_left_side]] = (
        ledger.boundary_left_spliced_neg[internal_left_side]
    )
    right_region_boundary_leff[left_boundary[internal_left_side]] = side[internal_left_side]

    for ref_idx in range(n_refs):
        region_start = int(ref_region_offsets[ref_idx])
        region_end = int(ref_region_offsets[ref_idx + 1])
        boundary_start = int(ref_boundary_offsets[ref_idx])
        boundary_end = int(ref_boundary_offsets[ref_idx + 1])
        is_terminal[boundary_start] = True
        is_terminal[boundary_end - 1] = True
        if region_end <= region_start:
            continue
        boundary_pos[boundary_start] = int(region_arrays.start[region_start])
        boundary_pos[boundary_end - 1] = int(region_arrays.end[region_end - 1])
        if region_end - region_start > 1:
            left_ends = np.asarray(region_arrays.end[region_start : region_end - 1], dtype=np.int64)
            right_starts = np.asarray(
                region_arrays.start[region_start + 1 : region_end], dtype=np.int64
            )
            if not np.array_equal(left_ends, right_starts):
                raise ValueError(
                    "build_boundary_table: regions must be contiguous within each reference."
                )
            boundary_pos[boundary_start + 1 : boundary_end - 1] = left_ends

    return BoundaryTable(
        boundary_pos=boundary_pos,
        ref_id=boundary_ref_id,
        ref_region_offsets=ref_region_offsets,
        ref_boundary_offsets=ref_boundary_offsets,
        is_terminal=is_terminal,
        left_region_unspliced_pos=left_region_unspliced_pos,
        left_region_unspliced_neg=left_region_unspliced_neg,
        left_region_spliced_pos=left_region_spliced_pos,
        left_region_spliced_neg=left_region_spliced_neg,
        right_region_unspliced_pos=right_region_unspliced_pos,
        right_region_unspliced_neg=right_region_unspliced_neg,
        right_region_spliced_pos=right_region_spliced_pos,
        right_region_spliced_neg=right_region_spliced_neg,
        left_region_boundary_leff=left_region_boundary_leff,
        right_region_boundary_leff=right_region_boundary_leff,
    )


def _assert_allclose(name: str, observed: np.ndarray, expected: np.ndarray) -> None:
    if not np.allclose(observed, expected, rtol=0.0, atol=0.0):
        raise AssertionError(f"validate_boundary_table: {name} does not match ledger.")


def validate_boundary_table(
    boundaries: BoundaryTable,
    region_arrays: RegionArrays,
    ledger: RegionCountLedger,
) -> None:
    """Validate exact reconstruction against current region-side ledger slots."""
    region_count = int(region_arrays.ref_id.shape[0])
    if boundaries.ref_boundary_offsets[-1] != boundaries.n_boundaries:
        raise AssertionError("validate_boundary_table: ref_boundary_offsets end mismatch.")

    region_indices = np.arange(region_count, dtype=np.int64)
    region_ref = np.asarray(region_arrays.ref_id, dtype=np.int64)
    ref_region_offsets = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    local_region = region_indices - ref_region_offsets[region_ref]
    left_boundary = boundaries.ref_boundary_offsets[region_ref] + local_region
    right_boundary = left_boundary + 1

    _assert_allclose(
        "right unspliced pos",
        boundaries.left_region_unspliced_pos[right_boundary],
        ledger.boundary_right_unspliced_pos,
    )
    _assert_allclose(
        "right unspliced neg",
        boundaries.left_region_unspliced_neg[right_boundary],
        ledger.boundary_right_unspliced_neg,
    )
    _assert_allclose(
        "right spliced pos",
        boundaries.left_region_spliced_pos[right_boundary],
        ledger.boundary_right_spliced_pos,
    )
    _assert_allclose(
        "right spliced neg",
        boundaries.left_region_spliced_neg[right_boundary],
        ledger.boundary_right_spliced_neg,
    )
    _assert_allclose(
        "left unspliced pos",
        boundaries.right_region_unspliced_pos[left_boundary],
        ledger.boundary_left_unspliced_pos,
    )
    _assert_allclose(
        "left unspliced neg",
        boundaries.right_region_unspliced_neg[left_boundary],
        ledger.boundary_left_unspliced_neg,
    )
    _assert_allclose(
        "left spliced pos",
        boundaries.right_region_spliced_pos[left_boundary],
        ledger.boundary_left_spliced_pos,
    )
    _assert_allclose(
        "left spliced neg",
        boundaries.right_region_spliced_neg[left_boundary],
        ledger.boundary_left_spliced_neg,
    )

    for ref_idx in range(int(region_arrays.n_refs)):
        region_start = int(ref_region_offsets[ref_idx])
        region_end = int(ref_region_offsets[ref_idx + 1])
        boundary_start = int(boundaries.ref_boundary_offsets[ref_idx])
        boundary_end = int(boundaries.ref_boundary_offsets[ref_idx + 1])
        if not boundaries.is_terminal[boundary_start] or not boundaries.is_terminal[boundary_end - 1]:
            raise AssertionError("validate_boundary_table: terminal flags missing.")
        if region_end - region_start <= 1:
            continue
        expected_pos = region_arrays.end[region_start : region_end - 1]
        observed_pos = boundaries.boundary_pos[boundary_start + 1 : boundary_end - 1]
        if not np.array_equal(observed_pos, expected_pos):
            raise AssertionError("validate_boundary_table: internal boundary positions mismatch.")

    if not math.isclose(
        float(boundaries.left_region_unspliced_pos.sum(dtype=np.float64)),
        float(ledger.boundary_right_unspliced_pos.sum(dtype=np.float64)),
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise AssertionError("validate_boundary_table: right-side mass conservation failed.")