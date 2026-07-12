"""The unified region↔boundary node chain (`calibration.node_chain`).

Locks: the genomic B-R-B-…-R-B interleave, the per-reference adjacency (terminals → -1), and consistency with
the existing region↔boundary index maps (a region's chain-neighbours ARE its flanking boundaries; a boundary's
chain-neighbours ARE its flanking regions).
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain
from rigel.calibration.region_arrays import boundary_region_indices, region_boundary_indices


def test_interleave_and_adjacency_two_refs():
    # ref0: 2 regions (r0,r1) → boundaries b0,b1,b2 ; ref1: 1 region (r2) → boundaries b3,b4
    rro = np.array([0, 2, 3])
    rbo = np.array([0, 3, 5])
    ch = build_node_chain(rro, rbo)

    assert ch.n_nodes == 3 + 5  # R + B
    assert ch.n_regions == 3 and ch.n_boundaries == 5
    # genomic order ref0: B0 R0 B1 R1 B2  | ref1: B3 R2 B4
    assert list(ch.kind) == [
        BOUNDARY,
        REGION,
        BOUNDARY,
        REGION,
        BOUNDARY,
        BOUNDARY,
        REGION,
        BOUNDARY,
    ]
    assert list(ch.ref_idx) == [0, 0, 1, 1, 2, 3, 2, 4]
    # adjacency (node ids); reference terminals are -1
    assert list(ch.left) == [-1, 0, 1, 2, 3, -1, 5, 6]
    assert list(ch.right) == [1, 2, 3, 4, -1, 6, 7, -1]
    assert list(ch.order) == list(range(8))


def test_terminals_are_sinks():
    rro = np.array([0, 2, 3])
    rbo = np.array([0, 3, 5])
    ch = build_node_chain(rro, rbo)
    # the first node of each ref has no left, the last no right
    assert ch.left[0] == -1 and ch.left[5] == -1
    assert ch.right[4] == -1 and ch.right[7] == -1
    # both ref terminals are BOUNDARY nodes (the B-R-B-…-B shape)
    assert ch.kind[0] == BOUNDARY and ch.kind[4] == BOUNDARY
    assert ch.kind[5] == BOUNDARY and ch.kind[7] == BOUNDARY


def test_consistent_with_index_maps():
    # A richer topology: ref0 k=3, ref1 k=1, ref2 k=2.
    rro = np.array([0, 3, 4, 6])
    rbo = np.array([0, 4, 6, 9])
    ch = build_node_chain(rro, rbo)
    lb, rb = region_boundary_indices(rro, rbo)  # per region → flanking boundary array-indices
    blr, brr = boundary_region_indices(
        rro, rbo
    )  # per boundary → flanking region array-indices (-1 terminal)

    for node in range(ch.n_nodes):
        L, R = int(ch.left[node]), int(ch.right[node])
        if ch.kind[node] == REGION:
            r = int(ch.ref_idx[node])
            # a region's chain neighbours are BOUNDARY nodes == its flanking boundaries
            assert ch.kind[L] == BOUNDARY and ch.ref_idx[L] == lb[r]
            assert ch.kind[R] == BOUNDARY and ch.ref_idx[R] == rb[r]
        else:
            b = int(ch.ref_idx[node])
            # a boundary's chain neighbours are REGION nodes == its flanking regions (or -1 at a terminal)
            if blr[b] >= 0:
                assert L >= 0 and ch.kind[L] == REGION and ch.ref_idx[L] == blr[b]
            else:
                assert L == -1
            if brr[b] >= 0:
                assert R >= 0 and ch.kind[R] == REGION and ch.ref_idx[R] == brr[b]
            else:
                assert R == -1


def test_total_node_count_invariant():
    # n_nodes == R + B == 2R + n_refs (k+1 boundaries per k regions)
    rro = np.array([0, 3, 4, 6])
    rbo = np.array([0, 4, 6, 9])
    ch = build_node_chain(rro, rbo)
    n_refs = rro.shape[0] - 1
    assert ch.n_nodes == 2 * int(rro[-1]) + n_refs
