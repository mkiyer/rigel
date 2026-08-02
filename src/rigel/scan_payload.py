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

WHAT THE NUMBERS MEAN. Every object stores **three integer sums** over the fragments that landed on it::

    count           Sum 1
    inv_length_sum  Sum round(2^32 / placements)     placements = L at a node, L − 1 at a 0-bp line
    length_sum      Sum L

⚠ **``inv_length_sum`` is NOT called ``density`` on purpose.** It is an exact, model-free density at an
edge — the opportunity ``L−1`` and the deposit ``1/(L−1)`` cancel identically — and it is *not* a density
at a node, where the opportunity is ``(node − L + 1)₊`` and nothing cancels. One word for two concepts is
the defect this naming avoids; see ``docs/NODE_DENSITY_DERIVATION.md``.

⭐ **``length_sum`` exists because the other two are blind to one real case.** At an edge the count row is
``(mu_g − 1, mu_r − 1)`` and the inv-length row is ``(1, 1)``, so the determinant is ``mu_g − mu_r``: when
gDNA and RNA share a mean fragment length the pair carries *zero* information about the split, at any
depth. ``length_sum`` is an independent tilt and removes that blind spot.

The trailing ``2`` on every bank is the **genome strand** — ``Strand.POS`` then ``Strand.NEG``, without
exception. Sense/antisense is transcript-relative, derived by the consumer from the junction's own
strand, and never stored.

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

#: The pool axis, named. ⚠ **These live here, with the schema, and not in the consumer that indexes by
#: them** — they are the accumulator's own enum in a third language, and a consumer holding a private
#: copy is how three files end up disagreeing about which row is which. A disagreement here would fit the
#: gDNA length model from the RNA pool and nothing downstream would look wrong;
#: `tests/calibration/test_fl.py` pins them against the executable specification's enum itself.
#:
#: Purity is the point (design §8): the two DNA_* contained pools are ~99 % gDNA on real data, and
#: RNA_SPLICED used an ANNOTATED junction with the splice OBSERVED — gDNA cannot be spliced. ⭐ The two
#: *_EXON "splash" pools are the only ON-TARGET gDNA population, so they are named rather than folded
#: into the gDNA model: on-target gDNA runs ~42 bp shorter than off-target (§8.2), and the shipped model
#: read a gDNA mean of 146.05 against the pure intergenic pool's 88.0 precisely by pooling them in.
#: There is deliberately NO pool for an exonic contained fragment or a multi-line crossing — those are
#: gDNA/RNA mixtures, and an impure pool is worse than a missing one.
POOL_DNA_INTERGENIC = 0  # contained in an intergenic node — pure gDNA
POOL_DNA_INTRONIC = 1  # contained in an intronic node — pure gDNA
POOL_DNA_INTRON_EXON = 2  # crossing one line, flanks {intron, exon} — on-target gDNA
POOL_DNA_INTERGENIC_EXON = 3  # crossing one line, {intergenic, exon} — on-target gDNA
POOL_RNA_SPLICED = 4  # used an annotated junction, splice OBSERVED — pure RNA


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
    #: ⭐ >1 surviving hypothesis, so the fragment's gap is undetermined. ⚠ NOT dropped — it is held
    #: WHOLE for the second pass, and the identity is `deposited + deferred + dropped_* == offered`.
    deferred_undetermined_gap: int
    unannotated_introns: int  # observed introns with no annotated junction
    contradictory_sj_strand: int  # the mates' motif tags disagreed; no splice trusted
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
class GapCensus:
    """⭐ The umbrella census: how each fragment whose gap needed resolving was resolved.

    Every fragment for which the enumeration produced at least one **spliced** hypothesis is counted here —
    its ``L`` depends on whether a gap intron is cut. The subclasses are exhaustive and mutually exclusive,
    and the three ``gap_deferred_*`` are exactly ``qc.deferred_undetermined_gap``.

    ⛔ **Its own axis, and NOT a splice type.** The umbrella cuts ACROSS the splice census: a certified-RNA
    ``SPLICED_ANNOT`` fragment with an intron in its mate gap needs resolving exactly as much as an
    ``UNSPLICED`` one does, so putting these on ``splice_type`` would need two labels per fragment and would
    break C2's property that the splice census sums to the library.

    ⛔ **There is no ``gap_resolved_unspliced``, and that is not an omission.** A spliced hypothesis cuts
    bases the unspliced one keeps, so ``L_spliced <= L_unspliced`` always, and the one arbitration filter is
    ``L <= max_fragment_length``. If the unspliced path survives the filter then every spliced path does
    too — so the unspliced path can never be the sole survivor, and the class it would name is unreachable.
    The ordering is pinned by ``test_gap_hypothesis_arbitration.py``.

    The field names are the specification's own ``Tally.gap_resolution`` keys.
    """

    #: One hypothesis survived, and it necessarily cuts something: the gap intron is real and ``L``
    #: excludes it. ⚠ Classifies the ARBITRATION, not the deposit — such a fragment can still be rejected
    #: afterwards as ``TOO_LONG``, which is a different question with its own counter.
    gap_resolved_spliced: int
    #: ⛔ The unspliced path against exactly one spliced path. The open question is **RNA or gDNA** — one
    #: bit, and it is the composition question calibration exists to answer.
    gap_deferred_rna_or_gdna: int
    #: ⛔ Two or more spliced paths and no unspliced one: gDNA cannot be spliced, so the molecule is
    #: certified RNA and the open question is purely **which structure**.
    gap_deferred_which_introns: int
    #: ⛔ Both questions at once.
    gap_deferred_both: int

    @classmethod
    def zeros(cls) -> "GapCensus":
        """No fragment had a gap to resolve. ⚠ The ONE spelling of that, so a hand-built payload cannot
        invent a second — and it is a real state: a library with no annotated intron in any mate gap."""
        return cls(**{field.name: 0 for field in dataclasses.fields(cls)})

    @classmethod
    def from_dict(cls, census: dict[str, Any]) -> "GapCensus":
        expected = {field.name for field in dataclasses.fields(cls)}
        missing = expected - set(census)
        if missing:
            raise ValueError(
                f"the scan's gap_resolution block is missing {sorted(missing)}. The subclasses are "
                f"exhaustive by construction, so a missing one is a partition that does not close."
            )
        return cls(**{name: int(census[name]) for name in expected})

    @property
    def deferred(self) -> int:
        """The three ``gap_deferred_*`` summed — which must equal ``qc.deferred_undetermined_gap``."""
        return (
            self.gap_deferred_rna_or_gdna + self.gap_deferred_which_introns + self.gap_deferred_both
        )


#: The deferred bank's arrays, in the order :meth:`Tally.deferred_arrays` emits them. ⚠ Read off the
#: dataclass below rather than written out twice; the tuple exists because the validation needs an order.
_DEFERRED_RECORD_FIELDS = ("ref", "start", "end", "align_strand", "sj_strand")


@dataclass(frozen=True, slots=True)
class DeferredFragments:
    """⭐ **THE SIDE BUFFER.** Fragments whose gap has more than one surviving explanation, held WHOLE.

    A fragment's unsequenced mate gap may hold no intron, one, or several, and *which* cannot be observed —
    the bases are not there. Deciding needs a fragment-length distribution that does not exist until the
    first pass is over, so the first pass does not guess: it enumerates, and when more than one hypothesis
    survives it holds the fragment here. `docs/PLAN_TWO_PASS.md` §5.

    ⭐ **The FRAGMENT is stored, never its consequences.** Object ids are large, derived, and would have to
    be kept consistent with the partition; the fragment is small and replays exactly. The second pass
    re-enters ``Accumulator::deposit`` with the chosen hypothesis, so there is one tally path and
    byte-identity with the specification is preserved for free.

    Two nested variable-length levels — fragments hold hypotheses, hypotheses hold introns — so there are
    two offset arrays at each level. Offsets are cumulative and always start at 0, so an empty bank is
    ``[0]`` and never ``[]``, and ``n_fragments`` is ``len(offsets) - 1``.

    ⛔ **THIS IS THE ONE BANK WHOSE ORDER IS OBSERVABLE.** Every other is a sum of integers and integer
    addition is associative, so a per-worker merge is exact whatever order the chunks arrived in. This is a
    list, so the C++ export **sorts it on the record's own content** before it crosses the ABI — otherwise
    the same BAM would give a different byte sequence at 1, 2, 4 and 8 workers with identical contents. Two
    records that tie on that key are identical records, so no tie-break is needed or possible.

    ⚠ Every array is ``int64``, including the two strand columns, which are ``int32`` everywhere else in the
    scanner. One dtype for one bank: the parity gate compares dtypes, and a narrowing conversion at the
    boundary compares equal by value.
    """

    ref: np.ndarray  # int64[n] — which reference; the drain replays onto THAT cut axis
    start: np.ndarray  # int64[n] — the CLIPPED extent, because that is what the drain must replay
    end: np.ndarray  # int64[n]
    align_strand: np.ndarray  # int64[n]
    sj_strand: np.ndarray  # int64[n] — the OBSERVED motif strand, if any
    observed_intron_offsets: np.ndarray  # int64[n + 1] — in PAIRS, into observed_introns
    observed_introns: (
        np.ndarray
    )  # int64[2 * n_observed] — flat (start, end); cut under EVERY hypothesis
    hypothesis_offsets: np.ndarray  # int64[n + 1] — into the per-hypothesis arrays below
    hypothesis_sj_strand: (
        np.ndarray
    )  # int64[n_hypotheses] — INFERRED; an observed motif always wins
    hypothesis_intron_offsets: np.ndarray  # int64[n_hypotheses + 1] — in PAIRS
    hypothesis_introns: np.ndarray  # int64[2 * n_implied] — the IMPLIED introns; empty ⇒ unspliced
    hypothesis_t_offsets: np.ndarray  # int64[n_hypotheses + 1]
    hypothesis_t: np.ndarray  # int64[n_supporting] — which transcripts imply this path

    @property
    def n_fragments(self) -> int:
        return int(self.hypothesis_offsets.shape[0]) - 1

    @property
    def n_hypotheses(self) -> int:
        return int(self.hypothesis_sj_strand.shape[0])

    @classmethod
    def empty(cls) -> "DeferredFragments":
        """Nothing was deferred. ⚠ Offsets are cumulative and start at 0, so every offset array is ``[0]``
        and never ``[]`` — spelled once here so a hand-built payload cannot get that boundary wrong."""
        return cls.from_dict(
            {
                field.name: np.asarray(
                    [0] if field.name.endswith("_offsets") else [], dtype=np.int64
                )
                for field in dataclasses.fields(cls)
            }
        )

    @classmethod
    def from_dict(cls, deferred: dict[str, Any]) -> "DeferredFragments":
        """Validate and adopt the flat arrays the C++ emits.

        ⚠ The two nested CSRs are re-derived rather than trusted, exactly as ``ref_node_offsets`` is. An
        offset array of the right LENGTH can still be inconsistent, and the second pass indexes every one
        of these — a truncated bank would score a fragment against another fragment's hypotheses, which is
        a wrong answer that looks entirely plausible.
        """
        names = [field.name for field in dataclasses.fields(cls)]
        missing = set(names) - set(deferred)
        if missing:
            raise ValueError(
                f"the scan's deferred bank is missing {sorted(missing)}. Every array is part of one "
                f"record, so a missing one is a fragment the second pass cannot replay."
            )
        arrays: dict[str, np.ndarray] = {}
        for name in names:
            array = np.ascontiguousarray(deferred[name])
            if array.dtype != np.int64:
                raise ValueError(
                    f"deferred[{name!r}] has dtype {array.dtype}, expected int64. One bank, one dtype — "
                    f"a narrowed or widened array compares equal by value and would hide the change."
                )
            if array.ndim != 1:
                raise ValueError(
                    f"deferred[{name!r}] has shape {array.shape}, expected one dimension"
                )
            arrays[name] = array

        n = int(arrays["hypothesis_offsets"].shape[0]) - 1
        if n < 0:
            raise ValueError(
                "deferred['hypothesis_offsets'] is empty; offsets are cumulative and start at 0, so an "
                "empty bank is [0] and never []"
            )
        for name in _DEFERRED_RECORD_FIELDS:
            if arrays[name].shape != (n,):
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but the offsets describe "
                    f"{n} fragments"
                )
        for offsets_name, values_name, stride in (
            ("observed_intron_offsets", "observed_introns", 2),
            ("hypothesis_offsets", "hypothesis_sj_strand", 1),
            ("hypothesis_intron_offsets", "hypothesis_introns", 2),
            ("hypothesis_t_offsets", "hypothesis_t", 1),
        ):
            offsets = arrays[offsets_name]
            if offsets.shape[0] == 0 or int(offsets[0]) != 0:
                raise ValueError(
                    f"deferred[{offsets_name!r}] must start at 0; offsets are cumulative and an empty "
                    f"level is [0]"
                )
            if np.any(np.diff(offsets) < 0):
                bad = int(np.argmax(np.diff(offsets) < 0))
                raise ValueError(
                    f"deferred[{offsets_name!r}] decreases at {bad}; a CSR offset array cannot go "
                    f"backwards"
                )
            if int(offsets[-1]) * stride != arrays[values_name].shape[0]:
                raise ValueError(
                    f"deferred[{offsets_name!r}] ends at {int(offsets[-1])} but "
                    f"{values_name} has {arrays[values_name].shape[0]} entries "
                    f"({int(offsets[-1])} x {stride} expected)"
                )
        n_hypotheses = int(arrays["hypothesis_offsets"][-1])
        for name in ("hypothesis_intron_offsets", "hypothesis_t_offsets"):
            if arrays[name].shape[0] != n_hypotheses + 1:
                raise ValueError(
                    f"deferred[{name!r}] has {arrays[name].shape[0]} entries but there are "
                    f"{n_hypotheses} hypotheses, so it must have {n_hypotheses + 1}"
                )
        # ⭐ A fragment is deferred BECAUSE two or more hypotheses survived, so a record carrying fewer
        # than two is a bank that lost them — and the second pass would then "choose" from a set of one and
        # deposit an answer nothing supported.
        runs = np.diff(arrays["hypothesis_offsets"])
        if n and int(runs.min()) < 2:
            bad = int(np.argmin(runs))
            raise ValueError(
                f"deferred fragment {bad} carries {int(runs[bad])} hypotheses. A fragment is deferred "
                f"only when two or more survived, so every record must hold at least two."
            )
        return cls(**arrays)

    def observed_introns_of(self, i: int) -> np.ndarray:
        """Fragment ``i``'s observed introns as an ``[k, 2]`` view. Cut under **every** hypothesis."""
        lo, hi = int(self.observed_intron_offsets[i]), int(self.observed_intron_offsets[i + 1])
        return self.observed_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)

    def hypothesis_introns_of(self, h: int) -> np.ndarray:
        """Hypothesis ``h``'s IMPLIED introns as an ``[k, 2]`` view. Empty ⇒ the unspliced hypothesis."""
        lo, hi = int(self.hypothesis_intron_offsets[h]), int(self.hypothesis_intron_offsets[h + 1])
        return self.hypothesis_introns[2 * lo : 2 * hi].reshape(hi - lo, 2)

    def supporting_t_of(self, h: int) -> np.ndarray:
        """Which candidate transcripts imply hypothesis ``h``. Empty for the unspliced one."""
        lo, hi = int(self.hypothesis_t_offsets[h]), int(self.hypothesis_t_offsets[h + 1])
        return self.hypothesis_t[lo:hi]


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
    node_contained_inv_length_sum: np.ndarray  # uint64[n_nodes, 2]
    node_contained_length_sum: np.ndarray  # uint64[n_nodes, 2] — Sum L
    node_spanning_count: np.ndarray  # uint32[n_nodes, 2] — one segment covers the node whole
    node_spanning_inv_length_sum: np.ndarray  # uint64[n_nodes, 2]
    node_spanning_length_sum: np.ndarray  # uint64[n_nodes, 2] — Sum L
    node_start_count: np.ndarray  # uint32[n_nodes] — THE invariant; sums to qc.deposited

    # -- contiguous edges: the 0-bp line between two adjacent nodes --
    edge_unspliced_count: np.ndarray  # uint32[n_edges, 2] — the mixture being deconvolved
    edge_unspliced_inv_length_sum: np.ndarray  # uint64[n_edges, 2]
    edge_unspliced_length_sum: np.ndarray  # uint64[n_edges, 2] — Sum L
    edge_spliced_count: np.ndarray  # uint32[n_edges, 2] — certified RNA: gDNA cannot be spliced
    edge_spliced_inv_length_sum: np.ndarray  # uint64[n_edges, 2]
    edge_spliced_length_sum: np.ndarray  # uint64[n_edges, 2] — Sum L

    # -- junction edges: one exact donor->acceptor jump. Pure RNA by construction --
    sj_count: np.ndarray  # uint32[n_sj, 2]
    sj_inv_length_sum: np.ndarray  # uint64[n_sj, 2]
    sj_length_sum: np.ndarray  # uint64[n_sj, 2] — Sum L

    # -- the fragment-length pools, binned at L, once per fragment --
    pool_lengths: np.ndarray  # int64[5, max_length + 1]

    #: uint32[max_length + 1] — ⭐ **C1: EVERY deposited fragment, binned at its own L, no purity
    #: condition.** The five pools above are deliberately CONDITIONED (an impure pool is worse than a
    #: missing one), so none of them is an unconditional anchor — which is why the empirical-Bayes
    #: shrinkage in `calibration.fl` took its anchor from the SCANNER, which measures length by two
    #: other rules over another population (`docs/FRAGMENT_LENGTH_AUDIT.md`). This is that anchor, in
    #: the accumulator's own frame.
    #: ⚠ "Unconditional GIVEN DEPOSIT": it excludes what the accumulator rejects (too long, ambiguous
    #: path, strand-undefined, empty), each counted in ``qc``. That is exactly the population the pools
    #: are drawn from, which is what makes it the right anchor rather than merely a convenient one.
    deposited_lengths: np.ndarray

    #: ⭐ **The side buffer** — fragments whose gap has more than one surviving explanation, held WHOLE for
    #: the second pass. ⚠ NOT a loss: ``deposited + deferred + dropped_* == offered``, and this bank holds
    #: the fragments ``qc.deferred_undetermined_gap`` counts.
    deferred: DeferredFragments
    #: ⭐ How each gap was resolved — its own axis, cutting across the splice census.
    gap_resolution: GapCensus

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
            ("node_contained_inv_length_sum", n_nodes, np.uint64),
            ("node_contained_length_sum", n_nodes, np.uint64),
            ("node_spanning_count", n_nodes, np.uint32),
            ("node_spanning_inv_length_sum", n_nodes, np.uint64),
            ("node_spanning_length_sum", n_nodes, np.uint64),
            ("edge_unspliced_count", n_edges, np.uint32),
            ("edge_unspliced_inv_length_sum", n_edges, np.uint64),
            ("edge_unspliced_length_sum", n_edges, np.uint64),
            ("edge_spliced_count", n_edges, np.uint32),
            ("edge_spliced_inv_length_sum", n_edges, np.uint64),
            ("edge_spliced_length_sum", n_edges, np.uint64),
            ("sj_count", n_sj, np.uint32),
            ("sj_inv_length_sum", n_sj, np.uint64),
            ("sj_length_sum", n_sj, np.uint64),
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

        deposited_lengths = np.ascontiguousarray(cal["deposited_lengths"], dtype=np.uint32)
        if deposited_lengths.shape != (max_length + 1,):
            raise ValueError(
                f"deposited_lengths has shape {deposited_lengths.shape}, expected ({max_length + 1},)"
            )
        # ⭐ THE C1 INVARIANT, checked at the door. Same externally-checkable form as
        # ``sum(node_start_count) == deposited`` and a DIFFERENT statement: that one says every fragment
        # was located in space, this one that every fragment was binned by length. A histogram that is
        # the anchor for every FL model in the tool must not be allowed in one fragment short.
        n_binned = int(deposited_lengths.sum())
        if n_binned != int(cal["qc"]["deposited"]):
            raise ValueError(
                f"deposited_lengths sums to {n_binned} but {int(cal['qc']['deposited'])} fragments were "
                "deposited; the unconditional length histogram must bin every one of them exactly once."
            )

        qc = ScanQC.from_dict(cal["qc"])
        deferred = DeferredFragments.from_dict(cal["deferred"])
        # ⭐ THE CONSERVATION HALF, refused at the door. `qc.deferred_undetermined_gap` says how many
        # fragments were held; this bank is supposed to BE them. The identity
        # `deposited + deferred + dropped_* == offered` is worth nothing if the deferred term is a number
        # with no fragments behind it — and a cache can be truncated or partially written, which is
        # precisely how a bank would arrive short of the counter that describes it.
        if deferred.n_fragments != qc.deferred_undetermined_gap:
            raise ValueError(
                f"the deferred bank holds {deferred.n_fragments} fragments but "
                f"qc.deferred_undetermined_gap is {qc.deferred_undetermined_gap}; the counter and the "
                f"fragments it counts must be the same population, or the second pass silently drains a "
                f"different one."
            )
        gap_resolution = GapCensus.from_dict(cal["gap_resolution"])
        if gap_resolution.deferred != qc.deferred_undetermined_gap:
            raise ValueError(
                f"the gap census's three deferred_* sum to {gap_resolution.deferred} but "
                f"qc.deferred_undetermined_gap is {qc.deferred_undetermined_gap}; the subclasses are "
                f"exhaustive by construction, so this is a partition that does not close."
            )
        return cls(
            cut_positions=cut_positions,
            **offsets,
            **banks,
            node_start_count=node_start_count,
            pool_lengths=pool_lengths.reshape(N_FRAGMENT_POOLS, max_length + 1),
            deposited_lengths=deposited_lengths,
            deferred=deferred,
            gap_resolution=gap_resolution,
            qc=qc,
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
    "POOL_DNA_INTERGENIC",
    "POOL_DNA_INTERGENIC_EXON",
    "POOL_DNA_INTRONIC",
    "POOL_DNA_INTRON_EXON",
    "POOL_RNA_SPLICED",
    "AccumulatorPayload",
    "DeferredFragments",
    "GapCensus",
    "ScanQC",
]
