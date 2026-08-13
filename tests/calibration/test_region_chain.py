"""The region<->edge chain: `N E N E ... E N`, and there are NO terminal slots.

     (S5.d)

⛔ **What this replaces.** The old chain was `B R B R ... R B` — a reference with `k` regions had `k + 1`
boundary slots, the two outermost of which were reference terminals carrying no data and existing only
so every region had something on both sides. The new partition has no such object: a reference with `k`
regions owns exactly `k - 1` interior lines, so the chain **starts and ends with a REGION**.

That is not a renaming. It is one fewer object per reference in one direction and two fewer in the other,
and every off-by-one in this file is a real defect the old shape could not express.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_chain import EDGE, REGION, build_region_chain


def _chain(regions_per_ref):
    """Build from per-reference region counts, deriving the edge offsets the way the payload does."""
    region_offsets = np.concatenate([[0], np.cumsum(regions_per_ref)]).astype(np.int64)
    edges_per_ref = np.maximum(np.asarray(regions_per_ref) - 1, 0)
    edge_offsets = np.concatenate([[0], np.cumsum(edges_per_ref)]).astype(np.int64)
    return build_region_chain(region_offsets, edge_offsets)


def test_a_reference_alternates_REGION_EDGE_and_ENDS_ON_A_REGION():
    chain = _chain([4])
    assert chain.n_regions_total == 4 and chain.n_edges_total == 3
    assert list(chain.kind) == [REGION, EDGE, REGION, EDGE, REGION, EDGE, REGION]
    assert list(chain.obj_idx) == [0, 0, 1, 1, 2, 2, 3]


def test_the_chain_length_is_2k_minus_1_not_2k_plus_1():
    """⭐ The whole shape change in one assertion. The old chain had `2k + 1` slots for `k` regions
    (`k` regions + `k + 1` boundaries); the new one has `2k - 1`."""
    for k in (1, 2, 5, 50):
        assert _chain([k]).n_slots == 2 * k - 1


def test_a_SINGLE_REGION_reference_has_no_edges_and_is_still_a_slot():
    """1 bp and single-region references are legal — 15,687 regions of length 1 exist at human scale — and a
    reference with one region has zero interior lines, not one."""
    chain = _chain([1])
    assert chain.n_slots == 1
    assert list(chain.kind) == [REGION]
    assert chain.left[0] == -1 and chain.right[0] == -1


def test_a_reference_with_NO_regions_contributes_NOTHING():
    chain = _chain([3, 0, 2])
    assert chain.n_slots == (2 * 3 - 1) + 0 + (2 * 2 - 1)
    assert chain.n_regions_total == 5 and chain.n_edges_total == 3


def test_references_do_not_BLEED_into_each_other():
    """⛔ The last region of one reference must not be adjacent to the first of the next. A chain that
    wrapped would relay density across a chromosome boundary, which is exactly the class of error that
    survives every aggregate check."""
    chain = _chain([3, 2])
    last_of_first = 2 * 3 - 2  # the final REGION slot of reference 0
    assert chain.kind[last_of_first] == REGION
    assert chain.right[last_of_first] == -1
    assert chain.left[last_of_first + 1] == -1


def test_adjacency_is_symmetric_and_alternates_TYPE():
    chain = _chain([4, 3])
    for slot in range(chain.n_slots):
        right = int(chain.right[slot])
        if right >= 0:
            assert int(chain.left[right]) == slot, "left/right must be inverses"
            assert chain.kind[right] != chain.kind[slot], "the chain is bipartite"


def test_an_EDGE_always_has_a_region_on_BOTH_sides():
    """An edge is the line BETWEEN two regions, so it can never be a terminal. This is the invariant the
    old shape could not state: its terminal boundaries had exactly one flank."""
    chain = _chain([5, 2, 1])
    for slot in np.flatnonzero(np.asarray(chain.kind) == EDGE):
        assert chain.left[slot] >= 0 and chain.right[slot] >= 0


def test_edge_e_sits_between_region_e_and_region_e_plus_one():
    """The endpoints are IMPLICIT — that is the design's word — so this arithmetic is the contract."""
    chain = _chain([4])
    for slot in np.flatnonzero(np.asarray(chain.kind) == EDGE):
        e = int(chain.obj_idx[slot])
        assert int(chain.obj_idx[chain.left[slot]]) == e
        assert int(chain.obj_idx[chain.right[slot]]) == e + 1


def test_obj_idx_is_a_bijection_onto_each_axis():
    """Every region and every edge appears exactly once; nothing is visited twice or skipped."""
    chain = _chain([4, 1, 3])
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx)
    assert sorted(idx[kind == REGION]) == list(range(chain.n_regions_total))
    assert sorted(idx[kind == EDGE]) == list(range(chain.n_edges_total))


def test_an_INCONSISTENT_edge_count_is_REFUSED_with_the_arithmetic_named():
    """Both offset arrays come from ONE payload, so a mismatch is an accumulator inconsistency rather
    than a stale index — and the error must say so, or the reader rebuilds the index for nothing."""
    with pytest.raises(ValueError, match="k regions has exactly k-1"):
        build_region_chain(np.array([0, 4], np.int64), np.array([0, 4], np.int64))


def test_the_predecessors_TERMINAL_SLOT_SHAPE_is_GONE():
    """⛔⛔ **RE-POINTED FROM THE NAMES TO THE SHAPE (2026-08-13), and that is not a weakening.**

    This gate used to assert ``region_chain`` had no ``REGION`` and no ``BOUNDARY`` attribute, because
    those were the PREDECESSOR's names. The 2026-08-12 vocabulary ruling **re-adopts exactly those
    words** — so the old assertion now bans the live constants and can never pass again.

    ⭐ **The words were never the defect; the SHAPE was.** The predecessor ran ``B R B R … R B`` with
    ``k + 1`` boundary slots per reference, the two outermost carrying no data and existing only so every
    region had an object on each side. What replaced it is ``R B R … B R``: ``2k − 1`` slots, starting
    and ending with a REGION, and a BOUNDARY that always has a region on both sides — an invariant the
    old shape could not state. That is what must not regress, and it is what this now pins.
    """
    chain = _chain([4])
    assert chain.n_slots == 2 * 4 - 1, "k regions must give 2k-1 slots, never k + (k+1)"
    assert chain.kind[0] == REGION and chain.kind[-1] == REGION, "the chain starts and ends with a REGION"
    # ⭐ the invariant the predecessor could not state: every BOUNDARY has a region on BOTH sides.
    # ⚠ The constant is still spelled EDGE until stage 3 of the rename; the SHAPE is what is pinned.
    for i in range(chain.n_slots):
        if chain.kind[i] == EDGE:
            assert chain.left[i] >= 0 and chain.right[i] >= 0
