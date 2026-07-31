"""RegionArrays geometry + the node↔contiguous-edge index mapping."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.node_chain import EDGE, build_node_chain
from rigel.calibration.region_arrays import (
    RegionArrays,
    edge_node_indices,
    node_right_edge,
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
            "region_id": np.arange(len(rows), dtype=np.int64),
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
# The node ↔ contiguous-edge mapping (S5.f) — the k+1 boundary axis is retired.
# ---------------------------------------------------------------------------
#
# ⭐ Topology deliberately MULTI-REFERENCE with three different shapes, because the whole class of
# defect here is "the first reference happens to work". ref0 = 3 nodes / 2 edges, ref1 = 2 nodes /
# 1 edge, ref2 = 1 node / 0 edges (a legal reference that owns no line at all). E == N − n_refs
# only counts non-empty refs: 6 nodes, 3 refs, 3 edges.
REF_ID = np.array([0, 0, 0, 1, 1, 2], dtype=np.int32)
REF_NODE_OFFSETS = np.array([0, 3, 5, 6], dtype=np.int64)
REF_EDGE_OFFSETS = np.array([0, 2, 3, 3], dtype=np.int64)


def test_node_right_edge_is_minus_one_at_every_reference_end():
    # node 2 is chr0's last, node 4 chr1's last, node 5 chr2's only — none owns a line to its right.
    np.testing.assert_array_equal(node_right_edge(REF_ID), [0, 1, -1, 2, -1, -1])


def test_edge_node_indices_are_the_adjacent_pair():
    lo, hi = edge_node_indices(REF_ID)
    np.testing.assert_array_equal(lo, [0, 1, 3])
    np.testing.assert_array_equal(hi, [1, 2, 4])
    # An edge NEVER straddles two references — that is the invariant the old k+1 axis could not state.
    assert np.all(REF_ID[lo] == REF_ID[hi])


def test_a_single_node_reference_owns_no_edge():
    lo, _ = edge_node_indices(np.array([0, 1, 2], dtype=np.int32))
    assert lo.size == 0
    np.testing.assert_array_equal(
        node_right_edge(np.array([0, 1, 2], dtype=np.int32)), [-1, -1, -1]
    )


def test_the_two_directions_round_trip():
    lo, hi = edge_node_indices(REF_ID)
    right = node_right_edge(REF_ID)
    np.testing.assert_array_equal(right[lo], np.arange(lo.size))  # edge → its left node → itself
    has_edge = right >= 0
    np.testing.assert_array_equal(lo[right[has_edge]], np.flatnonzero(has_edge))


def test_edge_numbering_matches_the_chain_built_from_the_payload_offsets():
    """⭐ The gate that matters: re-derive the SAME numbering by a DIFFERENT algorithm.

    ``edge_node_indices`` counts adjacent same-reference node pairs; ``build_node_chain`` walks the
    payload's two CSR offset arrays and lays out ``N E N E … N`` slot by slot. They must agree, or the
    calibration result's per-edge arrays are keyed to a different axis than the payload's — the exact
    class of defect that once dropped 476,719 of 476,732 real fragments (`CARRY_FORWARD.md` §3 trap 20)
    while every golden test passed. A validator that called the builder's own helper would prove
    nothing here (trap 1).
    """
    chain = build_node_chain(REF_NODE_OFFSETS, REF_EDGE_OFFSETS)
    is_edge = np.asarray(chain.kind) == EDGE
    obj = np.asarray(chain.obj_idx)
    left_slot = np.asarray(chain.left)[is_edge]
    right_slot = np.asarray(chain.right)[is_edge]

    # the chain's own answer: edge obj_idx e sits between these two node obj_idx values
    chain_edge_id = obj[is_edge]
    chain_lo = obj[left_slot]
    chain_hi = obj[right_slot]
    order = np.argsort(chain_edge_id)

    lo, hi = edge_node_indices(REF_ID)
    np.testing.assert_array_equal(chain_edge_id[order], np.arange(lo.size))
    np.testing.assert_array_equal(chain_lo[order], lo)
    np.testing.assert_array_equal(chain_hi[order], hi)


def test_ref_id_must_be_grouped():
    """A scrambled ``ref_id`` cannot produce a valid edge axis, and must say so rather than
    silently manufacturing edges that straddle references."""
    with pytest.raises(ValueError, match="grouped"):
        node_right_edge(np.array([0, 1, 0], dtype=np.int32))
