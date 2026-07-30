"""rigel accumulator — the executable SPECIFICATION.

The accumulator is the per-fragment tally built during the single-pass BAM scan; everything downstream
reads only its output. This module is the reference implementation, and the native accumulator is
required to reproduce it **byte for byte**. Where this file and a document disagree, this file wins.

    Design: ``docs/ACCUMULATOR_DESIGN.md``   ·   Plan: ``docs/IMPLEMENTATION_PLAN.md``
    Spec matrix: ``test_accumulator_spec.py``

THE MODEL
    The genome is a graph. **Nodes** are half-open intervals tiling each reference, numbered in genomic
    order. The 0-bp **lines** between adjacent nodes are **contiguous edges**; a **junction edge** is a
    directed donor→acceptor link taken from the annotation. A **fragment is a path** — its aligned
    blocks, joined across mate gaps and broken by introns::

        ref   |······ n0 ······|·· n1 ··|········ n2 ········|
        cuts  0               100      200                  600
        lines                  1        2
        path        [====== block ======]                       crosses line 1
        path        [= block =]~~intron~~[==== block ====]       crosses nothing; uses a junction

    Nodes count fragments **contained** (they fit inside) and **spanning** (they cover the node whole);
    edges count fragments **crossing**. Every object stores THREE sums over the fragments that landed
    on it: ``count`` = Sum 1, ``inv_length_sum`` = Sum 1/placements (fixed point), ``length_sum`` = Sum L.

WHAT EACH OBJECT'S NUMBERS MEAN
    With ``placements`` the number of admissible start positions — ``L`` at a node, ``L − 1`` at a 0-bp
    line — and a component of start-density ``rho`` and length distribution ``f``::

        E[count]          =  rho * E_f[placements]
        E[inv_length_sum]  =  rho * E_f[placements * (1/placements)]  =  rho   <- at an EDGE, exactly
        E[length_sum]      =  rho * E_f[placements * L]

    The opportunity factor cancels identically at an edge, for **any** length distribution, which is why
    no divisor and no length model appear there. It does **not** cancel at a node, where the opportunity
    is ``max(0, node − L + 1)``; there the term is a better-conditioned second moment and nothing more.

    ⚠ **This is why the channel is called ``inv_length_sum`` and NOT ``density``.** It IS a density at an
    edge, exactly; at a node it is not one. Naming it ``density`` put one word on two concepts and would
    mislead every consumer that reads a node. The three names are three sums and are honest everywhere.

    ⭐ **``length_sum`` is what separates two components that share a mean length.** At an edge the count
    row is ``(mu_g - 1, mu_r - 1)`` and the inv-length row is ``(1, 1)``, so the pair's determinant is
    ``mu_g - mu_r`` and it carries ZERO information about the gDNA/RNA split when the two means agree --
    at any depth. ``length_sum`` is a second, independent tilt and removes that blind spot.
    See ``docs/NODE_DENSITY_DERIVATION.md`` and ``scripts/design/observable_efficiency.py``.

⚠ NO PARTITIONING. Every crossed edge receives the FULL weight. The chance that a length-``L`` fragment
crosses a given line is proportional to ``L`` and the deposit is ``1/L``, so the two cancel and every
fragment length contributes equally to each edge. Dividing by the number of edges crossed destroys that
cancellation and makes the answer depend on node spacing — measured up to **3.6× low**. An edge count is
not a share of a conserved total; there is no total. Density is intensive, not extensive.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field

import numpy as np

from rigel.types import Strand


__all__ = [
    "INV_LENGTH_SCALE",
    "N_STRAND_COLUMNS",
    "STRAND_COLUMNS",
    "Accumulator",
    "DepositOutcome",
    "FragmentPool",
    "Partition",
    "Tally",
    "inv_length_quantum",
]

#: ⭐ THE STRAND CONVENTION, and it is the same one everywhere in this file.
#:
#: Every count and density array is ``[n, N_STRAND_COLUMNS]``, and the two columns **are** the two genome
#: strands — ``Strand.POS`` and ``Strand.NEG`` — following the index's own ``_pos``/``_neg`` column naming
#: (``reach_lo_pos``, ``reach_lo_neg``, …). Nothing is stored transcript-relative.
#:
#: **Sense / antisense is derived, never stored.** A junction edge carries its own genomic strand, so a
#: consumer that wants transcript-relative counts computes ``sense = (fragment strand == junction
#: strand)``. Storing it instead would put two conventions into one schema — which is exactly the defect
#: this replaced.
#:
#: ⚠ ``Strand`` has **four** values, not two: OR semantics make ``POS | NEG == AMBIGUOUS``, and ``NONE``
#: means no strand at all. Only POS and NEG name a column, so the mapping is a dict and a fragment whose
#: strand is neither raises ``KeyError`` here rather than being quietly filed under one of them. The
#: earlier version of this file mapped "anything that is not POS" to the minus column, which credited such
#: a fragment to the *wrong strand* instead of rejecting it.
#:
#: ⭐ **The column is indexed by ``align_strand`` — where the read ALIGNED — never by ``sj_strand``.**
#: The two are independent: ``align_strand`` is where the read sat on the genome, ``sj_strand`` is which way an
#: intron was spliced (from its ``GT..AG`` motif). Only a spliced read has the second, and it never selects
#: a column. Mixing them is what put one array into two conventions.
STRAND_COLUMNS: dict[int, int] = {Strand.POS: 0, Strand.NEG: 1}
N_STRAND_COLUMNS = len(STRAND_COLUMNS)


#: Densities accumulate as ``round(INV_LENGTH_SCALE / placements)`` in uint64. Integer addition is
#: associative, so a per-worker merge is bit-identical at any thread count — which float accumulation is
#: not, and that nondeterminism propagates to a ~2.6 % difference in the calibration output. The scale is
#: 2^32 because it keeps the quantisation error (≤ 6.9e-8 relative over ``L`` in [40, 1000]) below
#: float32's own epsilon while leaving ample headroom under the uint64 ceiling at realistic depth.
INV_LENGTH_SCALE = 1 << 32


def inv_length_quantum(placements: int) -> int:
    """``round(INV_LENGTH_SCALE / placements)``, rounding halves AWAY FROM ZERO.

    ⚠ The rounding mode is part of the contract: byte-identity across platforms and languages is
    undefined without it, and Python's built-in ``round`` is banker's rounding, which differs at ties.
    """
    if placements <= 0:
        raise ValueError(f"placements must be positive, got {placements}")
    return (2 * INV_LENGTH_SCALE + placements) // (2 * placements)


class DepositOutcome(enum.Enum):
    """Why a fragment did or did not deposit. Every rejection is counted, never silent."""

    DEPOSITED = "deposited"
    TOO_LONG = "dropped_too_long"
    EMPTY = "dropped_empty"
    #: The fragment has no single genome strand — NONE, or AMBIGUOUS (``POS | NEG``, which
    #: ``build_fragment`` produces when the mates agree in reference orientation). The column axis IS the
    #: genome strand, so there is no column to credit. Required by design §10.3, which lists
    #: strand-undefined fragments among the denominators the accumulator must emit.
    STRAND_UNDEFINED = "dropped_strand_undefined"
    #: The fragment's candidate transcripts imply **different intron sets**, so its path — and therefore
    #: ``L``, both quanta, the pool bin and the set of lines it crosses — is not determined. It deposits
    #: on nothing and waits for the second pass, which has the fragment length and the strand to
    #: discriminate with. Design §9.1.
    AMBIGUOUS_PATH = "dropped_ambiguous_path"


class FragmentPool(enum.IntEnum):
    """Fragment-length pools, chosen so that each is **pure by construction**.

    Purity is what removes the circularity: the length models are fitted from populations known to be
    one component, so nothing is ever estimated from the fragments it will later explain.

    There is deliberately **no pool** for an exonic contained fragment or a multi-line crossing — those
    are gDNA/RNA mixtures, and an impure pool is worse than a missing one.
    """

    DNA_INTERGENIC = 0  # contained in an intergenic node
    DNA_INTRONIC = 1  # contained in an intronic node
    DNA_INTRON_EXON = 2  # crossing exactly one line, flanks {intron, exon}
    DNA_INTERGENIC_EXON = 3  # crossing exactly one line, flanks {intergenic, exon}
    RNA_SPLICED = 4  # using an annotated junction, splice OBSERVED


#: Coarse node types, as ``signature.coarse_type_array`` emits them.
_TYPE_INTERGENIC, _TYPE_INTRON, _TYPE_EXON = 0, 1, 2

#: Flank-type pair -> splash pool. The key is sorted, so it is order-insensitive.
_SPLASH_POOL = {
    (_TYPE_INTRON, _TYPE_EXON): FragmentPool.DNA_INTRON_EXON,
    (_TYPE_INTERGENIC, _TYPE_EXON): FragmentPool.DNA_INTERGENIC_EXON,
}

#: Contained node type -> pure gDNA pool. An exonic node is a mixture and is absent by design.
_CONTAINED_POOL = {
    _TYPE_INTERGENIC: FragmentPool.DNA_INTERGENIC,
    _TYPE_INTRON: FragmentPool.DNA_INTRONIC,
}


# ---------------------------------------------------------------------------
# the partition: the cut axis, the node types, and the junction CSR
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Partition:
    """Everything the deposit addresses. Three axes, off by one from each other per reference.

    A reference contributing ``c`` cut positions owns ``c − 1`` nodes and ``c − 2`` contiguous edges
    (its interior lines). A reference with no nodes contributes no cuts at all::

        cuts    0        100       200       600        c = 4
        nodes   [  n0  ][   n1   ][   n2   ]            c - 1 = 3
        lines            line 1    line 2               c - 2 = 2

    Junctions are a CSR keyed by the **donor cut index**, which is the index the deposit has already
    computed while locating the lines its path crosses. It is cheap because every annotated intron has
    both endpoints as cuts, so "is this intron annotated?" reduces to "are both endpoints cuts, and is
    the pair registered?" — and if the start is not a cut, the table is never consulted.
    """

    cut_positions: np.ndarray  # int64[n_cuts] — flat, reference-major, ascending within a reference
    ref_cut_offsets: np.ndarray  # int64[n_refs + 1] — CSR over cut_positions
    node_types: np.ndarray  # uint8[n_nodes] — 0 intergenic / 1 intron / 2 exon
    ref_node_offsets: np.ndarray  # int64[n_refs + 1]
    ref_edge_offsets: np.ndarray  # int64[n_refs + 1]
    sj_offsets: np.ndarray  # int64[n_cuts + 1] — CSR over the donor cut index
    sj_acceptor_cut: np.ndarray  # int64[n_sj] — flat cut index of the intron's high end
    sj_strand: np.ndarray  # int8[n_sj]

    @property
    def n_refs(self) -> int:
        return int(self.ref_cut_offsets.shape[0]) - 1

    @property
    def n_nodes(self) -> int:
        return int(self.ref_node_offsets[-1])

    @property
    def n_edges(self) -> int:
        return int(self.ref_edge_offsets[-1])

    @property
    def n_sj(self) -> int:
        return int(self.sj_acceptor_cut.shape[0])

    @classmethod
    def from_cuts(cls, cuts_per_ref, node_types=None, junctions=()) -> "Partition":
        """Build from per-reference cut lists.

        ``junctions`` are ``(ref, intron_start, intron_end, sj_strand)``; both endpoints must be cuts on
        that reference. Junction ids are assigned by sorting on ``(donor cut, acceptor cut, sj_strand)``,
        so they are a deterministic function of the partition alone.
        """
        cuts_per_ref = [np.asarray(c, dtype=np.int64) for c in cuts_per_ref]
        n_refs = len(cuts_per_ref)
        for r, c in enumerate(cuts_per_ref):
            if c.size and (c.size < 2 or np.any(np.diff(c) <= 0)):
                raise ValueError(f"reference {r}: cuts must strictly increase, got {c.tolist()}")

        cut_offsets = np.zeros(n_refs + 1, np.int64)
        node_offsets = np.zeros(n_refs + 1, np.int64)
        edge_offsets = np.zeros(n_refs + 1, np.int64)
        for r, c in enumerate(cuts_per_ref):
            cut_offsets[r + 1] = cut_offsets[r] + c.size
            node_offsets[r + 1] = node_offsets[r] + max(c.size - 1, 0)
            edge_offsets[r + 1] = edge_offsets[r] + max(c.size - 2, 0)

        cut_positions = np.concatenate(cuts_per_ref) if n_refs else np.zeros(0, np.int64)
        if node_types is None:
            types = np.zeros(int(node_offsets[-1]), np.uint8)
        else:
            types = np.concatenate([np.asarray(t, np.uint8) for t in node_types])
        if types.shape[0] != int(node_offsets[-1]):
            raise ValueError(
                f"node_types has {types.shape[0]} entries; the partition has "
                f"{int(node_offsets[-1])} nodes"
            )

        donors, acceptors, sj_strands = [], [], []
        for ref, intron_start, intron_end, sj_strand in junctions:
            donor = _exact_cut(cut_positions, cut_offsets, ref, intron_start)
            acceptor = _exact_cut(cut_positions, cut_offsets, ref, intron_end)
            if donor < 0 or acceptor < 0:
                raise ValueError(
                    f"junction [{intron_start}, {intron_end}) on reference {ref} has an endpoint that "
                    f"is not a cut. Every annotated intron endpoint is a cut by construction, so this "
                    f"is a partition/annotation mismatch, not an unannotated junction."
                )
            donors.append(donor)
            acceptors.append(acceptor)
            sj_strands.append(int(sj_strand))
        donor_cut = np.asarray(donors, np.int64)
        acceptor_cut = np.asarray(acceptors, np.int64)
        sj_strand = np.asarray(sj_strands, np.int8)
        order = np.lexsort((sj_strand, acceptor_cut, donor_cut))
        donor_cut, acceptor_cut, sj_strand = donor_cut[order], acceptor_cut[order], sj_strand[order]

        n_cuts = int(cut_offsets[-1])
        sj_offsets = np.zeros(n_cuts + 1, np.int64)
        np.cumsum(np.bincount(donor_cut, minlength=n_cuts), out=sj_offsets[1:])
        return cls(
            cut_positions=cut_positions,
            ref_cut_offsets=cut_offsets,
            node_types=types,
            ref_node_offsets=node_offsets,
            ref_edge_offsets=edge_offsets,
            sj_offsets=sj_offsets,
            sj_acceptor_cut=acceptor_cut,
            sj_strand=sj_strand,
        )


def _exact_cut(cut_positions, cut_offsets, ref: int, position: int) -> int:
    """The flat cut index of ``position`` on ``ref``, or -1 if it is not a cut there."""
    first, last = int(cut_offsets[ref]), int(cut_offsets[ref + 1])
    if last <= first:
        return -1
    k = first + int(np.searchsorted(cut_positions[first:last], position))
    return k if k < last and int(cut_positions[k]) == position else -1


def _normalise_introns(introns, start: int, end: int) -> tuple[list[tuple[int, int]], int]:
    """The introns as a sorted, DISJOINT set clipped to ``[start, end)``, and how many were absorbed.

    ⚠ The caller is contracted to hand over a clean set, but a real BAM does not always allow it: the
    splice-motif strand tag is read once per RECORD (``XS`` for STAR, ``ts`` for minimap2 — the scanner
    auto-detects), so a pair whose mates disagree about an intron's acceptor produces two overlapping
    introns for one molecule. Measured on MO_3021: one read group in 875,670
    (``chr8:138206290-138206943``, three mutually overlapping introns).

    Normalising here — rather than trusting the caller — is what lets ``L`` be *defined* as the total
    length of the path's segments instead of computed by a second, independent formula. Two formulas for
    one quantity is how the same number ends up different in two places; with an overlap they disagreed
    by up to 1.5×, and a wide overlap drove ``L`` negative so a good fragment was silently discarded.
    """
    absorbed = 0
    merged: list[tuple[int, int]] = []
    for raw_start, raw_end in sorted((int(s), int(e)) for s, e in introns):
        intron_start, intron_end = max(raw_start, start), min(raw_end, end)
        if intron_end <= intron_start:  # zero-length, or entirely outside the fragment
            absorbed += 1
            continue
        # Overlapping OR ABUTTING introns are both malformed and both merge. Abutting means a
        # zero-length exon between them, which is physically impossible — there is no molecular
        # difference between a transcript with one and without it — so a single molecule can never
        # legitimately use two abutting introns. The index cannot produce the pair either: a
        # zero-length exon is dropped when the exon arrays are built, which fuses its two flanking
        # introns into one. So an abutting pair in a fragment is an alignment artifact.
        if merged and intron_start <= merged[-1][1]:
            previous_start, previous_end = merged[-1]
            merged[-1] = (previous_start, max(previous_end, intron_end))
            absorbed += 1
            continue
        merged.append((intron_start, intron_end))
    return merged, absorbed


def _segments(start: int, end: int, introns) -> list[tuple[int, int]]:
    """The path's contiguous genomic segments: ``[start, end)`` with the introns cut out.

    ``introns`` must already be sorted and disjoint (:func:`_normalise_introns`).
    """
    out, cursor = [], start
    for intron_start, intron_end in introns:
        if intron_start > cursor:
            out.append((cursor, intron_start))
        cursor = intron_end
    if end > cursor:
        out.append((cursor, end))
    return out


# ---------------------------------------------------------------------------
# the tally
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class Tally:
    """The accumulator's output. THREE arrays per population, all integer sums over its fragments.

    ``count`` = ``Sum 1`` carries the statistical power (a Beta-Binomial needs an integer);
    ``inv_length_sum`` = ``Sum 1/placements`` carries the level, and is an exact model-free density at an
    edge but NOT at a node; ``length_sum`` = ``Sum L`` is the second length tilt, and is what lets two
    components with the same mean fragment length be told apart at all.

    All three are integers, so **no fractional quantity ever enters a likelihood** and a per-worker merge
    is bit-identical at any thread count.
    """

    node_contained_count: np.ndarray  # uint32[n_nodes, 2]
    node_contained_inv_length_sum: np.ndarray  # uint64[n_nodes, 2]
    node_contained_length_sum: np.ndarray  # uint64[n_nodes, 2] — Sum L, the second length tilt
    node_spanning_count: np.ndarray  # uint32[n_nodes, 2]
    node_spanning_inv_length_sum: np.ndarray  # uint64[n_nodes, 2]
    node_spanning_length_sum: np.ndarray  # uint64[n_nodes, 2] — Sum L, the second length tilt
    node_start_count: np.ndarray  # uint32[n_nodes] — one per accepted fragment; THE invariant
    edge_unspliced_count: np.ndarray  # uint32[n_edges, 2]
    edge_unspliced_inv_length_sum: np.ndarray  # uint64[n_edges, 2]
    edge_unspliced_length_sum: np.ndarray  # uint64[n_edges, 2] — Sum L, the second length tilt
    edge_spliced_count: np.ndarray  # uint32[n_edges, 2] — certified RNA: gDNA cannot be spliced
    edge_spliced_inv_length_sum: np.ndarray  # uint64[n_edges, 2]
    edge_spliced_length_sum: np.ndarray  # uint64[n_edges, 2] — Sum L, the second length tilt
    edge_spliced_length_sum: np.ndarray  # uint64[n_edges, 2] — Sum L, the second length tilt
    sj_count: np.ndarray  # uint32[n_sj, 2]
    sj_inv_length_sum: np.ndarray  # uint64[n_sj, 2]
    sj_length_sum: np.ndarray  # uint64[n_sj, 2] — Sum L, the second length tilt
    sj_length_sum: np.ndarray  # uint64[n_sj, 2] — Sum L, the second length tilt
    pool_lengths: np.ndarray  # int64[5, max_fragment_length + 1] — binned at L, once per fragment
    qc: dict[str, int] = field(default_factory=dict)

    @classmethod
    def zeros(cls, n_nodes: int, n_edges: int, n_sj: int, max_length: int) -> "Tally":
        def counts(rows):
            return np.zeros((rows, N_STRAND_COLUMNS), np.uint32)

        def inv_length(rows):
            return np.zeros((rows, N_STRAND_COLUMNS), np.uint64)

        return cls(
            node_contained_count=counts(n_nodes),
            node_contained_inv_length_sum=inv_length(n_nodes),
            node_contained_length_sum=inv_length(n_nodes),
            node_spanning_count=counts(n_nodes),
            node_spanning_inv_length_sum=inv_length(n_nodes),
            node_spanning_length_sum=inv_length(n_nodes),
            node_start_count=np.zeros(n_nodes, np.uint32),
            edge_unspliced_count=counts(n_edges),
            edge_unspliced_inv_length_sum=inv_length(n_edges),
            edge_unspliced_length_sum=inv_length(n_edges),
            edge_spliced_count=counts(n_edges),
            edge_spliced_inv_length_sum=inv_length(n_edges),
            edge_spliced_length_sum=inv_length(n_edges),
            sj_count=counts(n_sj),
            sj_inv_length_sum=inv_length(n_sj),
            sj_length_sum=inv_length(n_sj),
            pool_lengths=np.zeros((len(FragmentPool), max_length + 1), np.int64),
            qc={outcome.value: 0 for outcome in DepositOutcome}
            | {
                "unannotated_introns": 0,
                "contradictory_sj_strand": 0,
                "sj_implicit_fragments": 0,
                "introns_absorbed": 0,
            },
        )


# ---------------------------------------------------------------------------
# the deposit
# ---------------------------------------------------------------------------


class Accumulator:
    """Deposit fragments into a :class:`Tally` against a :class:`Partition`.

    One public method. A fragment arrives as its genomic extent on ONE reference plus the introns inside
    it; mates are expected to be merged into that extent already, because the unsequenced gap between
    them is part of the molecule and must count toward ``L``.
    """

    __slots__ = ("partition", "max_fragment_length", "tally")

    def __init__(self, partition: Partition, max_fragment_length: int = 1000) -> None:
        self.partition = partition
        self.max_fragment_length = int(max_fragment_length)
        self.tally = Tally.zeros(
            partition.n_nodes, partition.n_edges, partition.n_sj, self.max_fragment_length
        )

    @property
    def n_edges(self) -> int:
        return self.partition.n_edges

    def deposit(
        self,
        ref: int,
        start: int,
        end: int,
        introns=(),
        align_strand: int = Strand.POS,
        sj_strand: int = Strand.NONE,
        sj_implicit: bool = False,
        path_ambiguous: bool = False,
    ) -> DepositOutcome:
        """Deposit one fragment; return why it did or did not land.

        ``[start, end)`` is the full genomic extent — leftmost block start to rightmost block end, mate
        gap included. ``introns`` are the excised gaps inside it as ``(start, end)`` pairs.
        ``sj_implicit`` marks a fragment whose splice was **implicit** rather than observed — the
        ``SPLICE_IMPLICIT`` class, where an annotated intron sits inside the unsequenced mate gap. It
        deposits normally but is barred from the pure-RNA length pool, because its splice is a product of
        the very model that pool is used to fit.

        ⚠ ``sj_implicit`` says the splice was *not sequenced*; it does **not** say the path is in doubt.
        A fragment whose candidate transcripts imply **different intron sets** has no determined ``L`` at
        all, and that is ``path_ambiguous`` — a rejection, not a flag on a deposit. The caller decides it,
        because only the caller has the candidate-transcript list; the accumulator's job is to make the
        loss COUNTED rather than silent, which is why it is an outcome here and not a `return` in the
        scanner.

        ⭐ **TWO STRANDS, AND THEY ARE INDEPENDENT.** Every read has the first; only a splice has the second.

        ``align_strand``
            The genomic strand the read **aligned** to, ``+`` or ``−``. Every read has one. It selects the
            array column and nothing else.
        ``sj_strand``
            A splice junction's strand, read from the aligner's **genomic-motif** tag (``XS`` for STAR,
            ``ts`` for minimap2 — auto-detected). A ``GT..AG`` intron is on ``+``, its reverse complement
            ``CT..AC`` on ``−``. Only spliced reads carry one, and it is used *only* to resolve an intron
            against the annotation (see :meth:`_sj_edge_id`).

        ⚠ The two are **unrelated**: an aligned strand says where the read sat, a splice strand says which
        way the intron was spliced, and neither constrains the other. Comparing them yields *sense* versus
        *antisense*, which is a derived quantity a consumer may compute and this accumulator never stores.
        ⭐ The shipped code collapses that comparison into one bool named ``primary``
        (``bam_scanner.cpp:1493``), which is how a dUTP library ended up with 0.6 % of spliced fragments in
        the column labelled *sense*.

        ⚠ ``sj_strand`` here is **observed**, per fragment. ``Partition.sj_strand`` is **annotated**, per
        junction edge. One quantity, two sources; :meth:`_sj_edge_id` compares them.
        """
        # ⚠ Checked FIRST, and before any geometry: the strand is a property of the fragment alone, and a
        # fragment with no single genome strand has no column in any bank. The scanner's gate at
        # `bam_scanner.cpp:1474-1480` did this before the old accumulator ever saw the fragment; doing it
        # here is what lets the loss be COUNTED instead of vanishing.
        column = STRAND_COLUMNS.get(align_strand)
        if column is None:
            return self._reject(DepositOutcome.STRAND_UNDEFINED)

        # ⚠ SECOND, and the order against the strand is part of the contract, because every fragment must
        # count exactly ONCE and a fragment can be both. It is filed under the strand, because that is
        # which denominator stays honest: `dropped_ambiguous_path` sizes the population the SECOND PASS
        # CAN RECOVER, and a fragment with no genome strand is not recoverable — the second pass resolves
        # which transcript, not which strand the read aligned to.
        if path_ambiguous:
            return self._reject(DepositOutcome.AMBIGUOUS_PATH)

        p = self.partition
        first_cut, last_cut = int(p.ref_cut_offsets[ref]), int(p.ref_cut_offsets[ref + 1])
        if last_cut - first_cut < 2:
            return self._reject(DepositOutcome.EMPTY)
        cuts = p.cut_positions[first_cut:last_cut]

        # clip to the reference; L is the CLIPPED length, so the placement count stays consistent
        start, end = max(int(start), int(cuts[0])), min(int(end), int(cuts[-1]))
        if end <= start:
            return self._reject(DepositOutcome.EMPTY)

        # ⚠ ONE definition of L: the total length of the path's segments. Deriving it any other way
        # invites two formulas for one quantity — see _normalise_introns.
        introns, absorbed = _normalise_introns(introns, start, end)
        segments = _segments(start, end, introns)
        length = sum(b - a for a, b in segments)
        if length <= 0:
            return self._reject(DepositOutcome.EMPTY)
        if length > self.max_fragment_length:
            return self._reject(DepositOutcome.TOO_LONG)
        self.tally.qc["introns_absorbed"] += absorbed

        # ── which annotated junctions does this path use? ─────────────────────────────────────────
        # ⚠ A contradictory motif strand disqualifies the whole fragment's splices, so it is checked once
        # here rather than per intron: `sj_strand` is the OR of the per-record strand-tag values, so
        # AMBIGUOUS means the mates disagreed about the same molecule. That is contradictory EVIDENCE, not
        # missing evidence, and it is counted on its own denominator — folding it into
        # `unannotated_introns` would poison the one metric whose job is measuring annotation coverage.
        if sj_strand == Strand.AMBIGUOUS:
            sj_ids: list[int] = []
            self.tally.qc["contradictory_sj_strand"] += 1
        else:
            sj_ids = [
                jid
                for jid in (
                    self._sj_edge_id(ref, intron_start, intron_end, sj_strand)
                    for intron_start, intron_end in introns
                )
                if jid >= 0
            ]
            self.tally.qc["unannotated_introns"] += len(introns) - len(sj_ids)
        spliced = bool(sj_ids)

        # ⚠ The path's own first and last COVERED base, not the fragment's extent. A leading or trailing
        # intron means the molecule does not begin at `start`: with introns [(150,480)] over [150,500) the
        # path is [480,500) and lives in a different node entirely. Using the extent would attribute the
        # fragment — and the start-count invariant — to a node it never touches.
        first_base, last_base = segments[0][0], segments[-1][1] - 1

        t = self.tally
        node_base, edge_base = int(p.ref_node_offsets[ref]), int(p.ref_edge_offsets[ref])
        t.node_start_count[node_base + self._local_node(cuts, first_base)] += 1
        t.qc[DepositOutcome.DEPOSITED.value] += 1
        if sj_implicit:
            t.qc["sj_implicit_fragments"] += 1

        # ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────
        # A line is crossed iff it lies strictly inside a segment, so per segment the crossed lines are a
        # contiguous index range and no container is needed. A node is SPANNED iff ONE segment crosses
        # both of its lines — not merely "both lines crossed", which would count a node the fragment
        # JUMPS OVER, whose two lines are touched by the two flanking segments from opposite sides.
        edge_count, edge_inv_length, edge_length = (
            (t.edge_spliced_count, t.edge_spliced_inv_length_sum, t.edge_spliced_length_sum)
            if spliced
            else (
                t.edge_unspliced_count,
                t.edge_unspliced_inv_length_sum,
                t.edge_unspliced_length_sum,
            )
        )
        quantum_edge = inv_length_quantum(length - 1) if length >= 2 else 0
        quantum_node = inv_length_quantum(length)
        n_crossed, sole_line = 0, -1
        for seg_start, seg_end in segments:
            first = int(np.searchsorted(cuts, seg_start, side="right"))
            last = int(np.searchsorted(cuts, seg_end, side="left"))
            for line in range(first, last):
                edge_count[edge_base + line - 1, column] += 1
                edge_inv_length[edge_base + line - 1, column] += quantum_edge
                edge_length[edge_base + line - 1, column] += length
            for line in range(first, last - 1):  # the node between two consecutive crossed lines
                t.node_spanning_count[node_base + line, column] += 1
                t.node_spanning_inv_length_sum[node_base + line, column] += quantum_node
                t.node_spanning_length_sum[node_base + line, column] += length
            if last > first:
                sole_line = first if (n_crossed == 0 and last - first == 1) else -1
                n_crossed += last - first

        for jid in sj_ids:
            t.sj_count[jid, column] += 1
            t.sj_inv_length_sum[jid, column] += quantum_edge
            t.sj_length_sum[jid, column] += length

        # ── contained: the WHOLE path lies inside ONE node ────────────────────────────────────────
        # ⚠ Not merely "crossed no line". An unannotated intron can swallow every line between two
        # blocks, leaving a fragment that crosses nothing yet straddles two nodes — crediting that as
        # contained would place its whole length in a node it only partly overlaps. Such a fragment is
        # neither contained nor crossing: it deposits on no object but is still counted (start_count),
        # so the loss is visible rather than silent.
        first_node = self._local_node(cuts, first_base)
        contained_node = -1
        if not sj_ids and first_node == self._local_node(cuts, last_base):
            contained_node = node_base + first_node
            t.node_contained_count[contained_node, column] += 1
            t.node_contained_inv_length_sum[contained_node, column] += quantum_node
            t.node_contained_length_sum[contained_node, column] += length

        pool = self._pool(spliced, sj_implicit, contained_node, sole_line, node_base)
        if pool is not None:
            t.pool_lengths[pool, length] += 1
        return DepositOutcome.DEPOSITED

    # -- helpers ----------------------------------------------------------------------------------

    def _reject(self, outcome: DepositOutcome) -> DepositOutcome:
        self.tally.qc[outcome.value] += 1
        return outcome

    @staticmethod
    def _local_node(cuts: np.ndarray, position: int) -> int:
        """Index within this reference of the node containing ``position``."""
        return min(max(int(np.searchsorted(cuts, position, side="right")) - 1, 0), cuts.size - 2)

    def _sj_edge_id(self, ref: int, intron_start: int, intron_end: int, sj_strand: int) -> int:
        """The junction-edge id for this intron, or -1 if it is not annotated.

        One to three iterations at human scale: 70.4 % of cuts are not a donor at all, and over those
        that are, the mean fan-out is 1.31.
        """
        p = self.partition
        donor = _exact_cut(p.cut_positions, p.ref_cut_offsets, ref, intron_start)
        if donor < 0:
            return -1
        acceptor = _exact_cut(p.cut_positions, p.ref_cut_offsets, ref, intron_end)
        if acceptor < 0:
            return -1
        for k in range(int(p.sj_offsets[donor]), int(p.sj_offsets[donor + 1])):
            if int(p.sj_acceptor_cut[k]) != acceptor:
                continue
            # ⚠ The strand filter applies only when the motif strand is DEFINITE. AMBIGUOUS means "no
            # strand information", not "a strand that matches nothing" — it used to be treated as the
            # latter, which silently demoted an annotated spliced fragment to an unspliced deposit,
            # credited `unannotated_introns` (poisoning the one metric that measures annotation coverage),
            # and dropped it from the pure RNA pool the length model is fitted from.
            #
            # The filter exists only to separate two junctions sharing a coordinate pair. ⭐ Measured: 0 of
            # 404,168 human junction coordinates are annotated on both strands, so it can only ever lose a
            # match, never disambiguate one. The resolver already reads a non-definite strand as "either"
            # (`sj_lookup_into` falls back to the union of POS and NEG), so this agrees with it.
            if sj_strand in STRAND_COLUMNS and int(p.sj_strand[k]) != int(sj_strand):
                continue
            return k
        return -1

    def _pool(self, spliced, sj_implicit, contained_node, sole_line, node_base):
        """The one length pool this fragment belongs to, or ``None``.

        Priority, so that every pool stays pure: an OBSERVED splice is unambiguously RNA; a contained
        fragment is typed by its node; a single-line crossing is a "splash" read typed by its two flank
        types. Anything else — an exonic contained fragment, a multi-line crossing — is a mixture and
        enters nothing.
        """
        types = self.partition.node_types
        if spliced:
            return None if sj_implicit else FragmentPool.RNA_SPLICED
        if contained_node >= 0:
            return _CONTAINED_POOL.get(int(types[contained_node]))
        if sole_line >= 0:
            flanks = (int(types[node_base + sole_line - 1]), int(types[node_base + sole_line]))
            return _SPLASH_POOL.get(tuple(sorted(flanks)))
        return None
