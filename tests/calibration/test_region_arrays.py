"""RegionArrays geometry + the region↔contiguous-boundary index mapping."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.region_chain import BOUNDARY, build_region_chain
from rigel.calibration.region_arrays import (
    RegionArrays,
    boundary_region_indices,
    region_right_boundary,
)
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_POS,
    TS_NEG,
    TS_NONE,
    TS_POS,
)


def _scrambled_region_df():
    # Two refs, rows deliberately out of (ref_id, start) order.
    rows = [
        ("chr2", 150, 400, BIT_EXON_NEG),
        ("chr1", 100, 200, BIT_EXON_POS),
        ("chr2", 0, 150, 0),
        ("chr1", 0, 100, 0),
        ("chr1", 200, 300, BIT_INTRON_POS),
    ]
    return pd.DataFrame(
        {
            "node_id": np.arange(len(rows), dtype=np.int64),
            "ref_name": pd.array([r[0] for r in rows], dtype="string"),
            "start": np.array([r[1] for r in rows], dtype=np.int64),
            "end": np.array([r[2] for r in rows], dtype=np.int64),
            "length": np.array([r[2] - r[1] for r in rows], dtype=np.int64),
            "signature": np.array([r[3] for r in rows], dtype=np.uint8),
        }
    )


def test_from_region_df_csr_ordering():
    df = _scrambled_region_df()
    ra = RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})

    # Sorted by (ref_id, start).
    np.testing.assert_array_equal(ra.ref_id, [0, 0, 0, 1, 1])
    np.testing.assert_array_equal(ra.start, [0, 100, 200, 0, 150])
    np.testing.assert_array_equal(ra.end, [100, 200, 300, 150, 400])
    np.testing.assert_array_equal(ra.ref_offsets, [0, 3, 5])
    np.testing.assert_array_equal(ra.region_size_bp, [100, 100, 100, 150, 250])
    np.testing.assert_array_equal(ra.strand_class, [TS_NONE, TS_POS, TS_POS, TS_NONE, TS_NEG])
    assert ra.n_regions == 5
    assert ra.n_refs == 2


def test_from_region_df_requires_signature():
    df = _scrambled_region_df().drop(columns=["signature"])
    try:
        RegionArrays.from_frame(df, {"chr1": 0, "chr2": 1})
    except ValueError as exc:
        assert "signature" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError for missing signature column")


# ---------------------------------------------------------------------------
# The region ↔ contiguous-boundary mapping (S5.f) — the k+1 boundary axis is retired.
# ---------------------------------------------------------------------------
#
# ⭐ Topology deliberately MULTI-REFERENCE with three different shapes, because the whole class of
# defect here is "the first reference happens to work". ref0 = 3 regions / 2 boundaries, ref1 = 2 regions /
# 1 boundary, ref2 = 1 region / 0 boundaries (a legal reference that owns no line at all). E == N − n_refs
# only counts non-empty refs: 6 regions, 3 refs, 3 boundaries.
REF_ID = np.array([0, 0, 0, 1, 1, 2], dtype=np.int32)
REF_REGION_OFFSETS = np.array([0, 3, 5, 6], dtype=np.int64)
REF_BOUNDARY_OFFSETS = np.array([0, 2, 3, 3], dtype=np.int64)


def test_region_right_boundary_is_minus_one_at_every_reference_end():
    # region 2 is chr0's last, region 4 chr1's last, region 5 chr2's only — none owns a line to its right.
    np.testing.assert_array_equal(region_right_boundary(REF_ID), [0, 1, -1, 2, -1, -1])


def test_boundary_region_indices_are_the_adjacent_pair():
    lo, hi = boundary_region_indices(REF_ID)
    np.testing.assert_array_equal(lo, [0, 1, 3])
    np.testing.assert_array_equal(hi, [1, 2, 4])
    # An boundary NEVER straddles two references — that is the invariant the old k+1 axis could not state.
    assert np.all(REF_ID[lo] == REF_ID[hi])


def test_a_single_region_reference_owns_no_boundary():
    lo, _ = boundary_region_indices(np.array([0, 1, 2], dtype=np.int32))
    assert lo.size == 0
    np.testing.assert_array_equal(
        region_right_boundary(np.array([0, 1, 2], dtype=np.int32)), [-1, -1, -1]
    )


def test_the_two_directions_round_trip():
    lo, hi = boundary_region_indices(REF_ID)
    right = region_right_boundary(REF_ID)
    np.testing.assert_array_equal(right[lo], np.arange(lo.size))  # boundary → its left region → itself
    has_boundary = right >= 0
    np.testing.assert_array_equal(lo[right[has_boundary]], np.flatnonzero(has_boundary))


def test_boundary_numbering_matches_the_chain_built_from_the_payload_offsets():
    """⭐ The gate that matters: re-derive the SAME numbering by a DIFFERENT algorithm.

    ``boundary_region_indices`` counts adjacent same-reference region pairs; ``build_region_chain`` walks the
    payload's two CSR offset arrays and lays out ``N E N E … N`` slot by slot. They must agree, or the
    calibration result's per-boundary arrays are keyed to a different axis than the payload's — the exact
    class of defect that once dropped 476,719 of 476,732 real fragments
    while every golden test passed. A validator that called the builder's own helper would prove
    nothing here (trap 1).
    """
    chain = build_region_chain(REF_REGION_OFFSETS, REF_BOUNDARY_OFFSETS)
    is_boundary = np.asarray(chain.kind) == BOUNDARY
    obj = np.asarray(chain.obj_idx)
    left_slot = np.asarray(chain.left)[is_boundary]
    right_slot = np.asarray(chain.right)[is_boundary]

    # the chain's own answer: boundary obj_idx e sits between these two region obj_idx values
    chain_edge_id = obj[is_boundary]
    chain_lo = obj[left_slot]
    chain_hi = obj[right_slot]
    order = np.argsort(chain_edge_id)

    lo, hi = boundary_region_indices(REF_ID)
    np.testing.assert_array_equal(chain_edge_id[order], np.arange(lo.size))
    np.testing.assert_array_equal(chain_lo[order], lo)
    np.testing.assert_array_equal(chain_hi[order], hi)


def test_ref_id_must_be_grouped():
    """A scrambled ``ref_id`` cannot produce a valid boundary axis, and must say so rather than
    silently manufacturing boundaries that straddle references."""
    with pytest.raises(ValueError, match="grouped"):
        region_right_boundary(np.array([0, 1, 0], dtype=np.int32))
