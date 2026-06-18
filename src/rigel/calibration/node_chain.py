"""The unified region↔boundary node chain (the bipartite belief-propagation graph).

`docs/calibration/bp_sweep_rebuild_plan.md` §1/§4. The calibration graph is a LINEAR bipartite chain of
**region** and **boundary** nodes, interleaved in genomic order. For a reference with ``k`` regions there are
``k + 1`` boundary slots, so its node sequence is::

    B0  R0  B1  R1  ...  R(k-1)  Bk        (boundary, region, boundary, region, ..., region, boundary)

This module builds that ordering + the left/right adjacency from the payload topology offsets alone (pure
index arithmetic — no C++, no payload reshape; the twin of :func:`region_arrays.boundary_region_indices`). The
sweep (P2) traverses ``order`` L→R then R→L; ``left``/``right`` give each node its single adjacent neighbour of
the OTHER type (``-1`` at a reference terminal — a propagation sink).

A node is addressed by ``(kind, idx)``: ``kind`` is :data:`REGION` or :data:`BOUNDARY`; ``idx`` indexes the
region arrays (``CalibrationSubstrate`` / ``RegionArrays``) or the boundary arrays (``BoundarySubstrate``)
respectively. This keeps the per-node sufficient statistics in their existing region-/boundary-keyed arrays —
the chain only sequences and links them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["REGION", "BOUNDARY", "NodeChain", "build_node_chain"]

REGION = 0
BOUNDARY = 1


@dataclass(frozen=True, slots=True)
class NodeChain:
    """The genomic-ordered region∪boundary chain + adjacency. All arrays length ``n_nodes = R + B``.

    ``order`` is the chain visiting sequence (genomic, per reference, B-R-B-…-R-B). ``kind``/``ref_idx`` give a
    node's type and its index into the region or boundary arrays. ``left``/``right`` are the adjacent nodes
    (positions in ``order``-space, i.e. node ids), ``-1`` at a reference terminal.
    """

    n_nodes: int
    kind: np.ndarray  # int8[n_nodes] — REGION or BOUNDARY
    ref_idx: np.ndarray  # int64[n_nodes] — index into the region (kind=REGION) or boundary (kind=BOUNDARY) arrays
    order: np.ndarray  # int64[n_nodes] — node ids in genomic visiting order
    left: np.ndarray  # int64[n_nodes] — left neighbour node id (-1 = reference start terminal)
    right: np.ndarray  # int64[n_nodes] — right neighbour node id (-1 = reference end terminal)

    @property
    def n_regions(self) -> int:
        return int(np.count_nonzero(self.kind == REGION))

    @property
    def n_boundaries(self) -> int:
        return int(np.count_nonzero(self.kind == BOUNDARY))


def build_node_chain(ref_region_offsets: np.ndarray, ref_boundary_offsets: np.ndarray) -> NodeChain:
    """Build the unified chain from the per-reference region/boundary CSR offsets.

    ``ref_region_offsets`` / ``ref_boundary_offsets`` are the payload's per-reference offset arrays (length
    ``n_refs + 1``); reference ``f`` owns regions ``[rro[f], rro[f+1])`` and boundaries ``[rbo[f], rbo[f+1])``
    with ``(rbo[f+1]-rbo[f]) == (rro[f+1]-rro[f]) + 1`` (k+1 boundaries for k regions). Node ids are assigned
    in genomic visiting order; a region/boundary keeps its own array index in ``ref_idx``.
    """
    rro = np.asarray(ref_region_offsets, dtype=np.int64)
    rbo = np.asarray(ref_boundary_offsets, dtype=np.int64)
    if rro.shape != rbo.shape:
        raise ValueError("ref_region_offsets and ref_boundary_offsets must share length (n_refs + 1).")
    n_refs = rro.shape[0] - 1
    r_total = int(rro[-1])
    b_total = int(rbo[-1])
    n_nodes = r_total + b_total

    kind = np.empty(n_nodes, dtype=np.int8)
    ref_idx = np.empty(n_nodes, dtype=np.int64)
    left = np.full(n_nodes, -1, dtype=np.int64)
    right = np.full(n_nodes, -1, dtype=np.int64)

    nid = 0  # next node id, assigned in genomic visiting order
    for f in range(n_refs):
        r0, r1 = int(rro[f]), int(rro[f + 1])
        b0, b1 = int(rbo[f]), int(rbo[f + 1])
        k = r1 - r0
        if (b1 - b0) != k + 1:
            raise ValueError(
                f"reference {f}: {b1 - b0} boundaries for {k} regions (expected {k + 1}); rebuild the index."
            )
        ref_start = nid
        # interleave B0 R0 B1 R1 ... R(k-1) Bk
        for i in range(k):
            kind[nid] = BOUNDARY
            ref_idx[nid] = b0 + i
            nid += 1
            kind[nid] = REGION
            ref_idx[nid] = r0 + i
            nid += 1
        kind[nid] = BOUNDARY  # trailing boundary Bk
        ref_idx[nid] = b1 - 1
        nid += 1
        # link consecutive nodes within this reference (terminals stay -1)
        ref_nodes = np.arange(ref_start, nid, dtype=np.int64)
        left[ref_nodes[1:]] = ref_nodes[:-1]
        right[ref_nodes[:-1]] = ref_nodes[1:]

    order = np.arange(n_nodes, dtype=np.int64)  # node ids already assigned in genomic order
    return NodeChain(n_nodes=n_nodes, kind=kind, ref_idx=ref_idx, order=order, left=left, right=right)
