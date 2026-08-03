"""The node<->edge chain: `N E N E ... E N`, and there are NO terminal slots.

     (S5.d)

⛔ **What this replaces.** The old chain was `B R B R ... R B` — a reference with `k` regions had `k + 1`
boundary slots, the two outermost of which were reference terminals carrying no data and existing only
so every region had something on both sides. The new partition has no such object: a reference with `k`
nodes owns exactly `k - 1` interior lines, so the chain **starts and ends with a NODE**.

That is not a renaming. It is one fewer object per reference in one direction and two fewer in the other,
and every off-by-one in this file is a real defect the old shape could not express.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.node_chain import EDGE, NODE, build_node_chain


def _chain(nodes_per_ref):
    """Build from per-reference node counts, deriving the edge offsets the way the payload does."""
    node_offsets = np.concatenate([[0], np.cumsum(nodes_per_ref)]).astype(np.int64)
    edges_per_ref = np.maximum(np.asarray(nodes_per_ref) - 1, 0)
    edge_offsets = np.concatenate([[0], np.cumsum(edges_per_ref)]).astype(np.int64)
    return build_node_chain(node_offsets, edge_offsets)


def test_a_reference_alternates_NODE_EDGE_and_ENDS_ON_A_NODE():
    chain = _chain([4])
    assert chain.n_nodes_total == 4 and chain.n_edges_total == 3
    assert list(chain.kind) == [NODE, EDGE, NODE, EDGE, NODE, EDGE, NODE]
    assert list(chain.obj_idx) == [0, 0, 1, 1, 2, 2, 3]


def test_the_chain_length_is_2k_minus_1_not_2k_plus_1():
    """⭐ The whole shape change in one assertion. The old chain had `2k + 1` slots for `k` regions
    (`k` regions + `k + 1` boundaries); the new one has `2k - 1`."""
    for k in (1, 2, 5, 50):
        assert _chain([k]).n_slots == 2 * k - 1


def test_a_SINGLE_NODE_reference_has_no_edges_and_is_still_a_slot():
    """1 bp and single-node references are legal — 15,687 nodes of length 1 exist at human scale — and a
    reference with one node has zero interior lines, not one."""
    chain = _chain([1])
    assert chain.n_slots == 1
    assert list(chain.kind) == [NODE]
    assert chain.left[0] == -1 and chain.right[0] == -1


def test_a_reference_with_NO_nodes_contributes_NOTHING():
    chain = _chain([3, 0, 2])
    assert chain.n_slots == (2 * 3 - 1) + 0 + (2 * 2 - 1)
    assert chain.n_nodes_total == 5 and chain.n_edges_total == 3


def test_references_do_not_BLEED_into_each_other():
    """⛔ The last node of one reference must not be adjacent to the first of the next. A chain that
    wrapped would relay density across a chromosome boundary, which is exactly the class of error that
    survives every aggregate check."""
    chain = _chain([3, 2])
    last_of_first = 2 * 3 - 2  # the final NODE slot of reference 0
    assert chain.kind[last_of_first] == NODE
    assert chain.right[last_of_first] == -1
    assert chain.left[last_of_first + 1] == -1


def test_adjacency_is_symmetric_and_alternates_TYPE():
    chain = _chain([4, 3])
    for slot in range(chain.n_slots):
        right = int(chain.right[slot])
        if right >= 0:
            assert int(chain.left[right]) == slot, "left/right must be inverses"
            assert chain.kind[right] != chain.kind[slot], "the chain is bipartite"


def test_an_EDGE_always_has_a_node_on_BOTH_sides():
    """An edge is the line BETWEEN two nodes, so it can never be a terminal. This is the invariant the
    old shape could not state: its terminal boundaries had exactly one flank."""
    chain = _chain([5, 2, 1])
    for slot in np.flatnonzero(np.asarray(chain.kind) == EDGE):
        assert chain.left[slot] >= 0 and chain.right[slot] >= 0


def test_edge_e_sits_between_node_e_and_node_e_plus_one():
    """The endpoints are IMPLICIT — that is the design's word — so this arithmetic is the contract."""
    chain = _chain([4])
    for slot in np.flatnonzero(np.asarray(chain.kind) == EDGE):
        e = int(chain.obj_idx[slot])
        assert int(chain.obj_idx[chain.left[slot]]) == e
        assert int(chain.obj_idx[chain.right[slot]]) == e + 1


def test_obj_idx_is_a_bijection_onto_each_axis():
    """Every node and every edge appears exactly once; nothing is visited twice or skipped."""
    chain = _chain([4, 1, 3])
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx)
    assert sorted(idx[kind == NODE]) == list(range(chain.n_nodes_total))
    assert sorted(idx[kind == EDGE]) == list(range(chain.n_edges_total))


def test_an_INCONSISTENT_edge_count_is_REFUSED_with_the_arithmetic_named():
    """Both offset arrays come from ONE payload, so a mismatch is an accumulator inconsistency rather
    than a stale index — and the error must say so, or the reader rebuilds the index for nothing."""
    with pytest.raises(ValueError, match="k nodes has exactly k-1"):
        build_node_chain(np.array([0, 4], np.int64), np.array([0, 4], np.int64))


def test_the_old_BOUNDARY_vocabulary_is_GONE():
    from rigel.calibration import node_chain as mod

    assert not hasattr(mod, "BOUNDARY")
    assert not hasattr(mod, "REGION")
