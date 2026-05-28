"""Region-level count totals for the density model (v4 Phase 1).

`RegionCountLedger` is a thin, frozen companion to `PayloadArrays` that
exposes the 12 per-channel slices of ``region_counts_sorted`` as named
views and provides freshly-summed totals on demand. It replaces the
``contained_unspliced_total`` / ``boundary_*_unspliced_total`` fields
that used to be precomputed on `PayloadArrays` (see v4 plan §5).

The ledger stores views only — `np.shares_memory(ledger.contained_unspliced_pos,
payload_arrays.region_counts_sorted)` is true. Totals are returned as
fresh ``float32[R]`` arrays; consumers are free to cache.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._arrays import PayloadArrays
from .signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
)


__all__ = [
    "RegionCountLedger",
    "build_region_count_ledger",
]


@dataclass(frozen=True, slots=True)
class RegionCountLedger:
    """Per-region channel slices + on-demand totals.

    Every array is a non-copying view into
    ``PayloadArrays.region_counts_sorted``. Totals are computed on each
    call; consumers that need to reuse them should cache locally.
    """

    contained_unspliced_pos: np.ndarray  # float32[R]
    contained_unspliced_neg: np.ndarray
    boundary_left_unspliced_pos: np.ndarray
    boundary_left_unspliced_neg: np.ndarray
    boundary_right_unspliced_pos: np.ndarray
    boundary_right_unspliced_neg: np.ndarray
    contained_spliced_pos: np.ndarray
    contained_spliced_neg: np.ndarray
    boundary_left_spliced_pos: np.ndarray
    boundary_left_spliced_neg: np.ndarray
    boundary_right_spliced_pos: np.ndarray
    boundary_right_spliced_neg: np.ndarray

    # Per-region per-compartment physical fragment support counts,
    # partitioned by splice class. Views over the sorted support
    # vectors on the underlying ``PayloadArrays``. uint32. These are
    # the discrete sample-size denominators consumed by the FMA
    # fusion engine.
    contained_unspliced_support: np.ndarray       # uint32[R]
    boundary_left_unspliced_support: np.ndarray   # uint32[R]
    boundary_right_unspliced_support: np.ndarray  # uint32[R]
    contained_spliced_support: np.ndarray         # uint32[R]
    boundary_left_spliced_support: np.ndarray     # uint32[R]
    boundary_right_spliced_support: np.ndarray    # uint32[R]

    # Per-region aggregate physical fragment support: distinct fragments
    # touching each region, partitioned by splice class. Views over the
    # sorted aggregate vectors on the underlying ``PayloadArrays``;
    # integer ESS for the EB exposure / density model. A fragment
    # fully spanning a region contributes one increment here but two
    # increments to the per-compartment counters above.
    unspliced_support: np.ndarray  # uint64[R]
    spliced_support: np.ndarray    # uint64[R]

    def unspliced_compartment_sum(self) -> np.ndarray:
        """Per-region sum of unspliced support across all three
        compartments. NOT equal to ``unspliced_support`` when any
        fragments fully span their region (those increment both the
        left and right boundary cells). Useful as the per-compartment
        sample-size aggregate for FMA-style consumers.
        """
        out = self.contained_unspliced_support.astype(np.uint64, copy=True)
        out += self.boundary_left_unspliced_support.astype(np.uint64, copy=False)
        out += self.boundary_right_unspliced_support.astype(np.uint64, copy=False)
        return out

    def spliced_compartment_sum(self) -> np.ndarray:
        """Per-region sum of spliced support across all three compartments.
        See ``unspliced_compartment_sum`` for the spanning-fragment caveat.
        """
        out = self.contained_spliced_support.astype(np.uint64, copy=True)
        out += self.boundary_left_spliced_support.astype(np.uint64, copy=False)
        out += self.boundary_right_spliced_support.astype(np.uint64, copy=False)
        return out

    # --- Unspliced totals --------------------------------------------------

    def contained_unspliced_total(self) -> np.ndarray:
        return self.contained_unspliced_pos + self.contained_unspliced_neg

    def boundary_left_unspliced_total(self) -> np.ndarray:
        return self.boundary_left_unspliced_pos + self.boundary_left_unspliced_neg

    def boundary_right_unspliced_total(self) -> np.ndarray:
        return self.boundary_right_unspliced_pos + self.boundary_right_unspliced_neg

    def boundary_unspliced_total(self) -> np.ndarray:
        return self.boundary_left_unspliced_total() + self.boundary_right_unspliced_total()

    # --- Spliced totals ----------------------------------------------------

    def contained_spliced_total(self) -> np.ndarray:
        return self.contained_spliced_pos + self.contained_spliced_neg

    def boundary_spliced_total(self) -> np.ndarray:
        return (
            self.boundary_left_spliced_pos
            + self.boundary_left_spliced_neg
            + self.boundary_right_spliced_pos
            + self.boundary_right_spliced_neg
        )

    def spliced_total(self) -> np.ndarray:
        return self.contained_spliced_total() + self.boundary_spliced_total()


def build_region_count_ledger(payload_arrays: PayloadArrays) -> RegionCountLedger:
    """Construct a `RegionCountLedger` over ``payload_arrays.region_counts_sorted``.

    The returned ledger holds views (no copies); the underlying buffer
    is the sorted region-counts matrix on ``payload_arrays``.
    """
    rc = payload_arrays.region_counts_sorted
    if rc.ndim != 2:
        raise ValueError(
            f"build_region_count_ledger: region_counts_sorted must be 2D, got shape {rc.shape}."
        )

    def col(compartment: int, splice: int, strand: int) -> np.ndarray:
        return rc[:, channel_index(compartment, splice, strand)]

    return RegionCountLedger(
        contained_unspliced_pos=col(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS),
        contained_unspliced_neg=col(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG),
        boundary_left_unspliced_pos=col(
            COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS
        ),
        boundary_left_unspliced_neg=col(
            COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG
        ),
        boundary_right_unspliced_pos=col(
            COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_POS
        ),
        boundary_right_unspliced_neg=col(
            COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG
        ),
        contained_spliced_pos=col(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS),
        contained_spliced_neg=col(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_NEG),
        boundary_left_spliced_pos=col(
            COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_POS
        ),
        boundary_left_spliced_neg=col(
            COMPARTMENT_BOUNDARY_LEFT, SPLICE_SPLICED, CHANNEL_STRAND_NEG
        ),
        boundary_right_spliced_pos=col(
            COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_POS
        ),
        boundary_right_spliced_neg=col(
            COMPARTMENT_BOUNDARY_RIGHT, SPLICE_SPLICED, CHANNEL_STRAND_NEG
        ),
        contained_unspliced_support=payload_arrays.region_contained_unspliced_support_sorted,
        boundary_left_unspliced_support=payload_arrays.region_boundary_left_unspliced_support_sorted,
        boundary_right_unspliced_support=payload_arrays.region_boundary_right_unspliced_support_sorted,
        contained_spliced_support=payload_arrays.region_contained_spliced_support_sorted,
        boundary_left_spliced_support=payload_arrays.region_boundary_left_spliced_support_sorted,
        boundary_right_spliced_support=payload_arrays.region_boundary_right_spliced_support_sorted,
        unspliced_support=payload_arrays.region_unspliced_support_sorted,
        spliced_support=payload_arrays.region_spliced_support_sorted,
    )
