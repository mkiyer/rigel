"""rigel.scan_payload — the BAM scan's accumulator tally, as a typed Python object.

    Spec: ``tests/native/_accumulator_reference.py``   ·   Design: ``docs/ACCUMULATOR_DESIGN.md`` §5.2

The C++ scanner returns ``result["calibration"]`` as a dict of **flat** ndarrays
(``BamScanner::build_result``); this module reshapes the two-column banks and validates the schema.
Everything downstream reads only this object.

⭐ **THE FIELD NAMES ARE THE SPECIFICATION'S ``Tally`` FIELD NAMES, CHARACTER FOR CHARACTER**, and so are
the C++ payload keys. One quantity, one name, in all three places — so there is no mapping table anywhere
that could drift. `tests/test_accumulator_payload.py` checks the schema against ``Tally`` itself rather
than against a written-out list, for the same reason.

THE AXES — three of them, off by one from each other per reference::

    cuts    0        100       200       600        c = 4 cuts on this reference
    nodes   [  n0  ][   n1   ][   n2   ]            c - 1 = 3 nodes
    lines            line 1    line 2               c - 2 = 2 contiguous edges

A reference contributing ``c`` cuts owns ``c − 1`` nodes and ``c − 2`` interior lines; one contributing
none owns neither, which is legal. Junction edges are their own axis, sliced by ``ref_sj_offsets``; the
flat slot order is the per-reference banks concatenated in reference order, which is what lets a
junction-edge id simply BE its slot.

WHAT THE NUMBERS MEAN. Every object stores an integer **count** and a fixed-point **density**
(``round(2^32 / placements)``), with ``placements = L`` at a node and ``L − 1`` at a 0-bp line. Both,
always: the count carries the statistical power, the density carries the level, and no fractional quantity
ever enters a likelihood. The trailing ``2`` on every bank is the **genome strand** — ``Strand.POS`` then
``Strand.NEG``, without exception. Sense/antisense is transcript-relative, derived by the consumer from the
junction's own strand, and never stored.

⚠ **OWNERSHIP: this object holds VIEWS, and it is the keep-alive.** ``np.ascontiguousarray(x, dtype=D)`` is
a **no-op** when the array already has dtype ``D``, so nothing here copies — the buffers belong to
capsules owned by the C++ side. Someone "adding a cast for safety" would silently double peak memory
against a 1.04 M-node partition. Do not.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Any

import numpy as np


#: The two columns of every bank ARE the two genome strands. ⚠ `Strand` has four values — OR semantics
#: make POS|NEG == AMBIGUOUS, and NONE means no strand — so only two of them name a column, and a fragment
#: carrying neither is rejected by the accumulator rather than filed under one.
N_STRAND_COLUMNS = 2

#: Five fragment-length pools, each pure by construction (design §8). The order is
#: `rigel::accumulator::FragmentPool` and `_accumulator_reference.FragmentPool`.
N_FRAGMENT_POOLS = 5


@dataclass(frozen=True, slots=True)
class ScanQC:
    """The denominators design §10.3 requires the accumulator to emit.

    Not optional and not derivable afterwards: **every conservation statement downstream must be able to
    name what it excluded.** Typed rather than a dict so that a misspelled denominator fails at this
    boundary instead of silently reading as zero somewhere far away.

    The field names are the specification's own ``Tally.qc`` keys.
    """

    deposited: int  # fragments that reached an object; == sum(node_start_count)
    dropped_too_long: int  # L above --max-fragment-length
    dropped_empty: int  # no path left after clipping to the reference
    dropped_strand_undefined: int  # align_strand named no column, so there was none to credit
    dropped_ambiguous_path: int  # the candidates imply >1 intron set, so L is undetermined (§9.1)
    unannotated_introns: int  # observed introns with no annotated junction
    contradictory_sj_strand: int  # the mates' motif tags disagreed; no splice trusted
    sj_implicit_fragments: int  # SPLICE_IMPLICIT: deposited, but barred from the pure-RNA pool
    introns_absorbed: int  # overlapping or abutting introns merged away

    @classmethod
    def from_dict(cls, qc: dict[str, Any]) -> "ScanQC":
        expected = {field.name for field in dataclasses.fields(cls)}
        missing = expected - set(qc)
        if missing:
            raise ValueError(
                f"the scan's qc block is missing {sorted(missing)}. Every one of these is a reported "
                f"denominator, so a missing key is a statement that cannot name what it excluded."
            )
        return cls(**{name: int(qc[name]) for name in expected})


@dataclass(frozen=True, slots=True)
class AccumulatorPayload:
    """One BAM scan's tally. Views over C++-owned buffers; this object is the keep-alive."""

    # -- the partition, echoed back so a consumer can locate every object without reloading the index --
    cut_positions: np.ndarray  # int64[n_cuts] — flat, reference-major, ascending within a reference
    ref_cut_offsets: np.ndarray  # int64[n_refs + 1] — CSR over cut_positions
    ref_node_offsets: np.ndarray  # int64[n_refs + 1]
    ref_edge_offsets: np.ndarray  # int64[n_refs + 1] — contiguous edges
    ref_sj_offsets: np.ndarray  # int64[n_refs + 1] — junction edges

    # -- nodes: two disjoint populations, each two genome-strand columns --
    node_contained_count: np.ndarray  # uint32[n_nodes, 2] — the whole path lies inside the node
    node_contained_density: np.ndarray  # uint64[n_nodes, 2]
    node_spanning_count: np.ndarray  # uint32[n_nodes, 2] — one segment covers the node whole
    node_spanning_density: np.ndarray  # uint64[n_nodes, 2]
    node_start_count: np.ndarray  # uint32[n_nodes] — THE invariant; sums to qc.deposited

    # -- contiguous edges: the 0-bp line between two adjacent nodes --
    edge_unspliced_count: np.ndarray  # uint32[n_edges, 2] — the mixture being deconvolved
    edge_unspliced_density: np.ndarray  # uint64[n_edges, 2]
    edge_spliced_count: np.ndarray  # uint32[n_edges, 2] — certified RNA: gDNA cannot be spliced
    edge_spliced_density: np.ndarray  # uint64[n_edges, 2]

    # -- junction edges: one exact donor->acceptor jump. Pure RNA by construction --
    sj_count: np.ndarray  # uint32[n_sj, 2]
    sj_density: np.ndarray  # uint64[n_sj, 2]

    # -- the fragment-length pools, binned at L, once per fragment --
    pool_lengths: np.ndarray  # int64[5, max_length + 1]

    qc: ScanQC
    max_length: int  # the fragment-length limit applied to L, and the pool-histogram width
    n_refs: int

    #: ⭐ Provenance, and it must cover **nodes AND edges**. The payload is edge-keyed by construction —
    #: its junction axis is meaningless against a different junction CSR — and `index.partition_hash`
    #: deliberately covers `nodes.feather` only. A 2026-07-29 flag fix rewrote every `edges.feather` while
    #: leaving every `nodes.feather` byte-identical, so a nodes-only key would have verified CLEAN against
    #: a stale cache. `None` when the scanner was driven without an index to hash.
    graph_hash: str | None = None

    # -- derived, never stored ------------------------------------------------------------------------

    @property
    def n_strand_columns(self) -> int:
        return N_STRAND_COLUMNS

    @property
    def n_nodes(self) -> int:
        return int(self.node_start_count.shape[0])

    @property
    def n_edges(self) -> int:
        return int(self.edge_unspliced_count.shape[0])

    @property
    def n_sj(self) -> int:
        return int(self.sj_count.shape[0])

    @classmethod
    def from_scan_result(
        cls, scan_result: dict[str, Any], graph_hash: str | None = None
    ) -> "AccumulatorPayload":
        cal = scan_result.get("calibration")
        if cal is None:
            raise ValueError(
                "scan_result['calibration'] is None; BamScanner.set_regions was not called"
            )

        n_refs = int(cal["n_refs"])
        n_cols = int(cal["n_strand_columns"])
        if n_cols != N_STRAND_COLUMNS:
            raise ValueError(
                f"the scan reports n_strand_columns={n_cols}, expected {N_STRAND_COLUMNS}. The two "
                f"columns are the two genome strands; any other count means the C++ and this schema "
                f"disagree about the channel axis."
            )
        max_length = int(cal["max_length"])

        offsets = {
            name: _offsets(cal, name, n_refs)
            for name in (
                "ref_cut_offsets",
                "ref_node_offsets",
                "ref_edge_offsets",
                "ref_sj_offsets",
            )
        }
        cut_positions = np.ascontiguousarray(cal["cut_positions"], dtype=np.int64)
        if int(offsets["ref_cut_offsets"][-1]) != cut_positions.shape[0]:
            raise ValueError(
                f"ref_cut_offsets ends at {int(offsets['ref_cut_offsets'][-1])} but cut_positions has "
                f"{cut_positions.shape[0]} entries"
            )

        # ⚠ Re-derived from the cut axis rather than trusted. An offset array of the right LENGTH can
        # still be inconsistent, and every consumer slices by these — a per-reference offset that drifts
        # is the defect class that once dropped 476,719 of 476,732 fragments while every golden passed.
        per_ref_cuts = np.diff(offsets["ref_cut_offsets"])
        for name, per_ref in (
            ("ref_node_offsets", np.maximum(per_ref_cuts - 1, 0)),
            ("ref_edge_offsets", np.maximum(per_ref_cuts - 2, 0)),
        ):
            expected = np.zeros(n_refs + 1, np.int64)
            np.cumsum(np.where(per_ref_cuts > 0, per_ref, 0), out=expected[1:])
            if not np.array_equal(offsets[name], expected):
                bad = int(np.argmax(offsets[name] != expected))
                raise ValueError(
                    f"{name}[{bad}] is {int(offsets[name][bad])} but the cut axis implies "
                    f"{int(expected[bad])}. A reference contributing c cuts owns c-1 nodes and c-2 "
                    f"interior lines, and none at all below two cuts."
                )

        n_nodes = int(offsets["ref_node_offsets"][-1])
        n_edges = int(offsets["ref_edge_offsets"][-1])
        n_sj = int(offsets["ref_sj_offsets"][-1])

        banks: dict[str, np.ndarray] = {}
        for name, rows, dtype in (
            ("node_contained_count", n_nodes, np.uint32),
            ("node_contained_density", n_nodes, np.uint64),
            ("node_spanning_count", n_nodes, np.uint32),
            ("node_spanning_density", n_nodes, np.uint64),
            ("edge_unspliced_count", n_edges, np.uint32),
            ("edge_unspliced_density", n_edges, np.uint64),
            ("edge_spliced_count", n_edges, np.uint32),
            ("edge_spliced_density", n_edges, np.uint64),
            ("sj_count", n_sj, np.uint32),
            ("sj_density", n_sj, np.uint64),
        ):
            banks[name] = _bank(cal, name, rows, dtype)

        node_start_count = np.ascontiguousarray(cal["node_start_count"], dtype=np.uint32)
        if node_start_count.shape != (n_nodes,):
            raise ValueError(
                f"node_start_count has shape {node_start_count.shape}, expected ({n_nodes},)"
            )
        pool_lengths = np.ascontiguousarray(cal["pool_lengths"], dtype=np.int64)
        if pool_lengths.size != N_FRAGMENT_POOLS * (max_length + 1):
            raise ValueError(
                f"pool_lengths has {pool_lengths.size} entries, expected "
                f"{N_FRAGMENT_POOLS} x (max_length + 1) = {N_FRAGMENT_POOLS * (max_length + 1)}"
            )

        return cls(
            cut_positions=cut_positions,
            **offsets,
            **banks,
            node_start_count=node_start_count,
            pool_lengths=pool_lengths.reshape(N_FRAGMENT_POOLS, max_length + 1),
            qc=ScanQC.from_dict(cal["qc"]),
            max_length=max_length,
            n_refs=n_refs,
            graph_hash=graph_hash,
        )


def _offsets(cal: dict[str, Any], name: str, n_refs: int) -> np.ndarray:
    array = np.ascontiguousarray(cal[name], dtype=np.int64)
    if array.shape != (n_refs + 1,):
        raise ValueError(f"{name} has shape {array.shape}, expected ({n_refs + 1},)")
    return array


def _bank(cal: dict[str, Any], name: str, rows: int, dtype: type) -> np.ndarray:
    """One two-column bank, reshaped from the flat array the C++ emits.

    ⚠ The dtype is asserted rather than coerced. ``ascontiguousarray`` would happily widen a uint32 count
    to uint64, which compares equal by value — so a schema drift would pass every value check downstream
    while the payload silently stopped matching the specification.
    """
    array = np.ascontiguousarray(cal[name])
    if array.dtype != dtype:
        raise ValueError(
            f"{name} has dtype {array.dtype}, expected {np.dtype(dtype)}. Counts are uint32 and "
            f"densities uint64; a widened array compares equal by value and would hide the change."
        )
    if array.size != rows * N_STRAND_COLUMNS:
        raise ValueError(
            f"{name} has {array.size} entries, which does not divide into {rows} rows x "
            f"{N_STRAND_COLUMNS} strand columns"
        )
    return array.reshape(rows, N_STRAND_COLUMNS)


__all__ = [
    "N_FRAGMENT_POOLS",
    "N_STRAND_COLUMNS",
    "AccumulatorPayload",
    "ScanQC",
]
