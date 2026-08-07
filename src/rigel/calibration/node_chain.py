"""rigel.calibration.node_chain — the node<->edge chain the belief-propagation sweep traverses.

       Gate: ``tests/calibration/test_node_chain.py``

The calibration graph is a **linear bipartite chain of NODE and EDGE slots, interleaved in genomic
order**. A reference with ``k`` nodes owns exactly ``k − 1`` interior lines, so its slot sequence is::

    N0  E0  N1  E1  ...  E(k-2)  N(k-1)          2k − 1 slots

⭐ **It starts and ends with a NODE, and there are no terminal slots.** That is the whole shape change
from the predecessor, which ran ``B R B R … R B`` with ``k + 1`` boundary slots per reference — the two
outermost carrying no data and existing only so every region had an object on each side. A contiguous
edge is the line BETWEEN two adjacent nodes; there is no such line before the first or after the last, so
the object does not exist rather than existing empty. **An edge therefore always has a node on both
sides**, an invariant the old shape could not state.

Edge endpoints are **implicit**: edge ``i`` lies between node ``i`` and node ``i + 1``. Nothing stores
them, and this module is the one place that arithmetic lives.

A slot is addressed by ``(kind, obj_idx)``: ``kind`` is :data:`NODE` or :data:`EDGE`, and ``obj_idx``
indexes the node axis or the contiguous-edge axis respectively. That keeps every per-object statistic in
its own payload-shaped array — the chain only sequences and links them.

⚠ **Junction edges are NOT chain slots.** The graph is a DAG but not a polytree: every junction edge
closes an undirected loop, so a junction must be a FACTOR on its endpoint nodes and never a message
channel (— never break a cycle by dropping a junction edge,
that re-isolates the exon the edge exists for).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["EDGE", "NODE", "NodeChain", "NodeDeconv", "build_node_chain"]

NODE = 0
EDGE = 1


@dataclass(frozen=True, slots=True)
class NodeChain:
    """The genomic-ordered node∪edge chain and its adjacency. All arrays have length ``n_slots``.

    Slot ids are assigned in genomic visiting order, so ``order`` would be ``arange`` and is not stored.
    ``left``/``right`` give each slot its single adjacent slot of the OTHER kind, ``-1`` at a reference
    terminal — which is a propagation sink, and is now always a NODE.
    """

    kind: np.ndarray  # int8[n_slots] — NODE or EDGE
    obj_idx: np.ndarray  # int64[n_slots] — index into the node axis, or the contiguous-edge axis
    left: np.ndarray  # int64[n_slots] — adjacent slot id, -1 at a reference start
    right: np.ndarray  # int64[n_slots] — adjacent slot id, -1 at a reference end
    n_nodes_total: int
    n_edges_total: int

    @property
    def n_slots(self) -> int:
        return int(self.kind.shape[0])

    @property
    def is_node(self) -> np.ndarray:
        return self.kind == NODE

    @property
    def is_edge(self) -> np.ndarray:
        return self.kind == EDGE


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# ⭐⭐ ONE SLOT'S DECONVOLUTION RESULT — vocabulary, and it lives here because THREE LAYERS need it.
# It was defined in the STRAND family (layer 4) and imported by `node_geometry` and `simplex_logodds`
# (layer 3) and `sweep` (layer 6), so three layers reached UPWARD for a type. `module_census.py` is what
# made that visible, and the repair is the one a layering violation always asks for: the TYPE belongs at
# the bottom, not the code that happened to define it first. ⚠ It is not a strand concept — the pie
# `(f_pos, f_neg, f_g)` is the tool's central datum, and a slot is what carries it.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class NodeDeconv:
    """Per-node deconvolution result. TWO disjoint uses, hence the optional halves:

    * the per-node SOLVE (`simplex_logodds._solve_nodes_logodds_all`) returns the **composition** —
      ``*_frac`` + ``*_frac_var`` — and no mass (a node's mass is a per-FACE quantity; the solve is
      face-invariant, so a single ``*_mass`` here would be meaningless);
    * the chain PROJECTION (`sweep.chain_node_deconv` / `chain_edge_deconv`) returns the
      **mass** the downstream `CalibrationResult` consumes, and no precision.
    """

    gdna_frac: (
        np.ndarray
    )  # float64[K] — the node's gDNA composition (face-invariant; mass = frac·M_face)
    # per-strand RNA fractions of the UNSPLICED mass (posterior means; f_pos+f_neg+gdna_frac = 1), populated
    # by the simplex sweep for the per-strand RNA imputation (the bipartite R↔B↔R chain).
    rna_pos_frac: "np.ndarray | None" = None  # float64[K] — f_pos
    rna_neg_frac: "np.ndarray | None" = None  # float64[K] — f_neg
    # per-component posterior variances in LOG-FRACTION space — `Var(log f_c)`, NOT `Var(f_c)`. They are
    # grid moments of `log f_c` over the λ lattice (`simplex_logodds._solve_nodes_logodds`), because the
    # message currency is a log-density and the send precision `1/(Var(log f_c) + 1/n + σ²_transfer)` is
    # log-space throughout. ⚠ They are therefore NOT bounded by ¼ and routinely exceed it — a consumer that
    # needs the LINEAR `Var(f_c)` must convert (delta method: `Var(f_c) ≈ f_c²·Var(log f_c)`, as
    # `sweep.solve_chain` does when it builds `_var_fg` for `composition_logvar`). Set by the per-node
    # solve, consumed when a node emits a message. None on the chain region/boundary projections (precision
    # is a chain-node property, not needed by the downstream EM prior).
    # the PROJECTION's consumed output (calibrate/derive read ONLY these); None on the per-node solve.
    gdna_mass: "np.ndarray | None" = None  # float64[K]
    rna_mass: "np.ndarray | None" = None  # float64[K]  (= (1−gdna_frac)·M_unspliced + spliced mass)
    gdna_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_g)
    rna_pos_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_pos)
    rna_neg_frac_var: "np.ndarray | None" = None  # float64[K] — Var(log f_neg)


def build_node_chain(ref_node_offsets: np.ndarray, ref_edge_offsets: np.ndarray) -> NodeChain:
    """Build the chain from the payload's two per-reference CSR offset arrays.

    Reference ``f`` owns nodes ``[rno[f], rno[f+1])`` and contiguous edges ``[reo[f], reo[f+1])``, with
    ``edges == max(nodes − 1, 0)``. A reference with no nodes contributes nothing at all, which is legal.

    ⚠ Both arrays come from ONE accumulator payload, so a mismatch between them is an accumulator /
    payload inconsistency and **not** a stale index — rebuilding will not fix it, and the error says so.
    """
    node_offsets = np.asarray(ref_node_offsets, dtype=np.int64)
    edge_offsets = np.asarray(ref_edge_offsets, dtype=np.int64)
    if node_offsets.shape != edge_offsets.shape:
        raise ValueError(
            f"ref_node_offsets has shape {node_offsets.shape} and ref_edge_offsets "
            f"{edge_offsets.shape}; both are per-reference CSR arrays of length n_refs + 1."
        )
    n_refs = node_offsets.shape[0] - 1
    nodes_per_ref = np.diff(node_offsets)
    edges_per_ref = np.diff(edge_offsets)
    expected = np.maximum(nodes_per_ref - 1, 0)
    if not np.array_equal(edges_per_ref, expected):
        bad = int(np.argmax(edges_per_ref != expected))
        raise ValueError(
            f"reference {bad}: the payload reports {int(edges_per_ref[bad])} contiguous edges for "
            f"{int(nodes_per_ref[bad])} nodes, but a reference with k nodes has exactly k-1 interior "
            f"lines (expected {int(expected[bad])}). There are no terminal edge slots: an edge is the "
            f"line BETWEEN two adjacent nodes. Both offset arrays come from ONE accumulator payload, so "
            f"this is an accumulator/payload inconsistency, not a stale index — rebuilding will not fix it."
        )

    n_slots = int(nodes_per_ref.sum() + edges_per_ref.sum())
    kind = np.empty(n_slots, dtype=np.int8)
    obj_idx = np.empty(n_slots, dtype=np.int64)
    left = np.full(n_slots, -1, dtype=np.int64)
    right = np.full(n_slots, -1, dtype=np.int64)

    slot = 0
    for f in range(n_refs):
        k = int(nodes_per_ref[f])
        if k == 0:
            continue
        node_base, edge_base = int(node_offsets[f]), int(edge_offsets[f])
        first_slot = slot
        # N0 E0 N1 E1 ... E(k-2) N(k-1): a node, then the edge to its right, except after the last node
        for i in range(k):
            kind[slot] = NODE
            obj_idx[slot] = node_base + i
            slot += 1
            if i < k - 1:
                kind[slot] = EDGE
                obj_idx[slot] = edge_base + i
                slot += 1
        # link consecutive slots WITHIN this reference; the two terminals keep -1, so a sweep cannot
        # relay across a reference boundary
        ref_slots = np.arange(first_slot, slot, dtype=np.int64)
        left[ref_slots[1:]] = ref_slots[:-1]
        right[ref_slots[:-1]] = ref_slots[1:]

    return NodeChain(
        kind=kind,
        obj_idx=obj_idx,
        left=left,
        right=right,
        n_nodes_total=int(node_offsets[-1]),
        n_edges_total=int(edge_offsets[-1]),
    )
