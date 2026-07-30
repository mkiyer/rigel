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
    edges count fragments **crossing**. Every object stores an integer count and a fixed-point density.

WHAT EACH OBJECT'S NUMBERS MEAN
    With ``placements`` the number of admissible start positions — ``L`` at a node, ``L − 1`` at a 0-bp
    line — and a component of start-density ``rho`` and length distribution ``f``::

        E[count]    =  rho * E_f[placements]
        E[density]  =  rho * E_f[placements * (1 / placements)]  =  rho     <- at an EDGE, exactly

    The opportunity factor cancels identically at an edge, for **any** length distribution, which is why
    no divisor and no length model appear there. It does **not** cancel at a node, where the opportunity
    is ``max(0, node − L + 1)``; there the density term is a better-conditioned second moment and
    nothing more.

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
    "CHANNEL_MINUS",
    "CHANNEL_PLUS",
    "DENSITY_SCALE",
    "N_STRANDS",
    "Accumulator",
    "DepositOutcome",
    "FragmentPool",
    "Partition",
    "Tally",
    "density_quantum",
    "genome_channel",
]

#: ⭐ THE STRAND CONVENTION, and it is the same one everywhere in this file.
#:
#: Every count and density array has exactly two channels, indexed by the fragment's own **genome**
#: strand: ``CHANNEL_PLUS`` for the ``+`` strand of the reference, ``CHANNEL_MINUS`` for ``−``. Nothing
#: is stored transcript-relative.
#:
#: **Sense / antisense is derived, never stored.** A junction edge carries its own genomic strand, so a
#: consumer that wants transcript-relative counts computes ``sense = (fragment strand == junction
#: strand)``. Storing it instead would put two conventions into one schema — which is exactly the defect
#: this replaced.
CHANNEL_PLUS, CHANNEL_MINUS = 0, 1
N_STRANDS = 2


def genome_channel(strand: int) -> int:
    """The channel index for a fragment aligned to genome strand ``strand``."""
    return CHANNEL_PLUS if strand == Strand.POS else CHANNEL_MINUS


#: Densities accumulate as ``round(DENSITY_SCALE / placements)`` in uint64. Integer addition is
#: associative, so a per-worker merge is bit-identical at any thread count — which float accumulation is
#: not, and that nondeterminism propagates to a ~2.6 % difference in the calibration output. The scale is
#: 2^32 because it keeps the quantisation error (≤ 6.9e-8 relative over ``L`` in [40, 1000]) below
#: float32's own epsilon while leaving ample headroom under the uint64 ceiling at realistic depth.
DENSITY_SCALE = 1 << 32


def density_quantum(placements: int) -> int:
    """``round(DENSITY_SCALE / placements)``, rounding halves AWAY FROM ZERO.

    ⚠ The rounding mode is part of the contract: byte-identity across platforms and languages is
    undefined without it, and Python's built-in ``round`` is banker's rounding, which differs at ties.
    """
    if placements <= 0:
        raise ValueError(f"placements must be positive, got {placements}")
    return (2 * DENSITY_SCALE + placements) // (2 * placements)


class DepositOutcome(enum.Enum):
    """Why a fragment did or did not deposit. Every rejection is counted, never silent."""

    DEPOSITED = "deposited"
    TOO_LONG = "dropped_too_long"
    EMPTY = "dropped_empty"


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
    junction_offsets: np.ndarray  # int64[n_cuts + 1] — CSR over the donor cut index
    junction_acceptor_cut: (
        np.ndarray
    )  # int64[n_junctions] — flat cut index of the intron's high end
    junction_strand: np.ndarray  # int8[n_junctions]

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
    def n_junctions(self) -> int:
        return int(self.junction_acceptor_cut.shape[0])

    @classmethod
    def from_cuts(cls, cuts_per_ref, node_types=None, junctions=()) -> "Partition":
        """Build from per-reference cut lists.

        ``junctions`` are ``(ref, intron_start, intron_end, strand)``; both endpoints must be cuts on
        that reference. Junction ids are assigned by sorting on ``(donor cut, acceptor cut, strand)``,
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

        donors, acceptors, strands = [], [], []
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
            strands.append(int(sj_strand))
        donor_cut = np.asarray(donors, np.int64)
        acceptor_cut = np.asarray(acceptors, np.int64)
        strand = np.asarray(strands, np.int8)
        order = np.lexsort((strand, acceptor_cut, donor_cut))
        donor_cut, acceptor_cut, strand = donor_cut[order], acceptor_cut[order], strand[order]

        n_cuts = int(cut_offsets[-1])
        junction_offsets = np.zeros(n_cuts + 1, np.int64)
        np.cumsum(np.bincount(donor_cut, minlength=n_cuts), out=junction_offsets[1:])
        return cls(
            cut_positions=cut_positions,
            ref_cut_offsets=cut_offsets,
            node_types=types,
            ref_node_offsets=node_offsets,
            ref_edge_offsets=edge_offsets,
            junction_offsets=junction_offsets,
            junction_acceptor_cut=acceptor_cut,
            junction_strand=strand,
        )


def _exact_cut(cut_positions, cut_offsets, ref: int, position: int) -> int:
    """The flat cut index of ``position`` on ``ref``, or -1 if it is not a cut there."""
    lo, hi = int(cut_offsets[ref]), int(cut_offsets[ref + 1])
    if hi <= lo:
        return -1
    k = lo + int(np.searchsorted(cut_positions[lo:hi], position))
    return k if k < hi and int(cut_positions[k]) == position else -1


def _normalise_introns(introns, lo: int, hi: int) -> tuple[list[tuple[int, int]], int]:
    """The introns as a sorted, DISJOINT set clipped to ``[lo, hi)``, and how many were absorbed.

    ⚠ The caller is contracted to hand over a clean set, but a real BAM does not always allow it: the
    scanner reads the ``XS`` tag once per RECORD, so a pair whose mates disagree about an intron's
    acceptor produces two overlapping introns for one molecule. Measured on MO_3021: one read group in
    875,670 (``chr8:138206290-138206943``, three mutually overlapping introns).

    Normalising here — rather than trusting the caller — is what lets ``L`` be *defined* as the total
    length of the path's segments instead of computed by a second, independent formula. Two formulas for
    one quantity is how the same number ends up different in two places; with an overlap they disagreed
    by up to 1.5×, and a wide overlap drove ``L`` negative so a good fragment was silently discarded.
    """
    absorbed = 0
    merged: list[tuple[int, int]] = []
    for start, end in sorted((int(s), int(e)) for s, e in introns):
        start, end = max(start, lo), min(end, hi)
        if end <= start:  # zero-length, or entirely outside the fragment
            absorbed += 1
            continue
        # Overlapping OR ABUTTING introns are both malformed and both merge. Abutting means a
        # zero-length exon between them, which is physically impossible — there is no molecular
        # difference between a transcript with one and without it — so a single molecule can never
        # legitimately use two abutting introns. The index cannot produce the pair either: a
        # zero-length exon is dropped when the exon arrays are built, which fuses its two flanking
        # introns into one. So an abutting pair in a fragment is an alignment artifact.
        if merged and start <= merged[-1][1]:
            previous_start, previous_end = merged[-1]
            merged[-1] = (previous_start, max(previous_end, end))
            absorbed += 1
            continue
        merged.append((start, end))
    return merged, absorbed


def _segments(lo: int, hi: int, introns) -> list[tuple[int, int]]:
    """The path's contiguous genomic segments: ``[lo, hi)`` with the introns cut out.

    ``introns`` must already be sorted and disjoint (:func:`_normalise_introns`).
    """
    out, cursor = [], lo
    for intron_start, intron_end in introns:
        if intron_start > cursor:
            out.append((cursor, intron_start))
        cursor = intron_end
    if hi > cursor:
        out.append((cursor, hi))
    return out


# ---------------------------------------------------------------------------
# the tally
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class Tally:
    """The accumulator's output. Two arrays per population: an integer count and a fixed-point density.

    Both, always — the count carries the statistical power (a Beta-Binomial needs an integer), the
    density carries the level, and **no fractional quantity ever enters a likelihood**.
    """

    node_contained_count: np.ndarray  # uint32[n_nodes, 2]
    node_contained_density: np.ndarray  # uint64[n_nodes, 2]
    node_spanning_count: np.ndarray  # uint32[n_nodes, 2]
    node_spanning_density: np.ndarray  # uint64[n_nodes, 2]
    node_start_count: np.ndarray  # uint32[n_nodes] — one per accepted fragment; THE invariant
    edge_unspliced_count: np.ndarray  # uint32[n_edges, 2]
    edge_unspliced_density: np.ndarray  # uint64[n_edges, 2]
    edge_spliced_count: np.ndarray  # uint32[n_edges, 2] — certified RNA: gDNA cannot be spliced
    edge_spliced_density: np.ndarray  # uint64[n_edges, 2]
    junction_count: np.ndarray  # uint32[n_junctions, 2]
    junction_density: np.ndarray  # uint64[n_junctions, 2]
    pool_lengths: np.ndarray  # int64[5, max_fragment_length + 1] — binned at L, once per fragment
    qc: dict[str, int] = field(default_factory=dict)

    @classmethod
    def zeros(cls, n_nodes: int, n_edges: int, n_junctions: int, max_length: int) -> "Tally":
        def counts(rows):
            return np.zeros((rows, N_STRANDS), np.uint32)

        def density(rows):
            return np.zeros((rows, N_STRANDS), np.uint64)

        return cls(
            node_contained_count=counts(n_nodes),
            node_contained_density=density(n_nodes),
            node_spanning_count=counts(n_nodes),
            node_spanning_density=density(n_nodes),
            node_start_count=np.zeros(n_nodes, np.uint32),
            edge_unspliced_count=counts(n_edges),
            edge_unspliced_density=density(n_edges),
            edge_spliced_count=counts(n_edges),
            edge_spliced_density=density(n_edges),
            junction_count=counts(n_junctions),
            junction_density=density(n_junctions),
            pool_lengths=np.zeros((len(FragmentPool), max_length + 1), np.int64),
            qc={outcome.value: 0 for outcome in DepositOutcome}
            | {"unannotated_introns": 0, "inferred_intron_fragments": 0, "introns_absorbed": 0},
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
            partition.n_nodes, partition.n_edges, partition.n_junctions, self.max_fragment_length
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
        motif_strand: int = Strand.NONE,
        introns_inferred: bool = False,
    ) -> DepositOutcome:
        """Deposit one fragment; return why it did or did not land.

        ``[start, end)`` is the full genomic extent — leftmost block start to rightmost block end, mate
        gap included. ``introns`` are the excised gaps inside it as ``(start, end)`` pairs.
        ``align_strand`` is the fragment's alignment strand and ``motif_strand`` its splice-motif strand;
        together they select the channel. ``introns_inferred`` marks a fragment whose splice was inferred
        rather than observed: it deposits normally but is barred from the pure-RNA length pool, because
        its splice is a product of the very model that pool is used to fit.
        """
        p = self.partition
        first_cut, last_cut = int(p.ref_cut_offsets[ref]), int(p.ref_cut_offsets[ref + 1])
        if last_cut - first_cut < 2:
            return self._reject(DepositOutcome.EMPTY)
        cuts = p.cut_positions[first_cut:last_cut]

        # clip to the reference; L is the CLIPPED length, so the placement count stays consistent
        lo, hi = max(int(start), int(cuts[0])), min(int(end), int(cuts[-1]))
        if hi <= lo:
            return self._reject(DepositOutcome.EMPTY)

        # ⚠ ONE definition of L: the total length of the path's segments. Deriving it any other way
        # invites two formulas for one quantity — see _normalise_introns.
        introns, absorbed = _normalise_introns(introns, lo, hi)
        segments = _segments(lo, hi, introns)
        length = sum(b - a for a, b in segments)
        if length <= 0:
            return self._reject(DepositOutcome.EMPTY)
        if length > self.max_fragment_length:
            return self._reject(DepositOutcome.TOO_LONG)
        self.tally.qc["introns_absorbed"] += absorbed

        # ── which annotated junctions does this path use? ─────────────────────────────────────────
        junction_ids = [
            jid
            for jid in (
                self._junction_of(ref, intron_start, intron_end, motif_strand)
                for intron_start, intron_end in introns
            )
            if jid >= 0
        ]
        self.tally.qc["unannotated_introns"] += len(introns) - len(junction_ids)
        spliced = bool(junction_ids)

        # ⭐ ONE strand convention, everywhere: the channel is the fragment's own GENOME strand.
        # Sense/antisense is a DERIVED quantity, not a stored one — a junction edge already carries its
        # own genomic strand, so a consumer computes `sense = (fragment strand == junction strand)` and
        # nothing is lost. Storing some banks by genome strand and others by sense is what put one array
        # into two conventions: measured on real cfRNA, 65-69 % of node_spanning deposits come from
        # spliced fragments and 40-44 % of those landed in the opposite column from the unspliced
        # fragments beside them.
        channel = genome_channel(align_strand)

        # ⚠ The path's own first and last COVERED base, not the fragment's extent. A leading or trailing
        # intron means the molecule does not begin at `lo`: with introns [(150,480)] over [150,500) the
        # path is [480,500) and lives in a different node entirely. Using the extent would attribute the
        # fragment — and the start-count invariant — to a node it never touches.
        first_base, last_base = segments[0][0], segments[-1][1] - 1

        t = self.tally
        node_base, edge_base = int(p.ref_node_offsets[ref]), int(p.ref_edge_offsets[ref])
        t.node_start_count[node_base + self._local_node(cuts, first_base)] += 1
        t.qc[DepositOutcome.DEPOSITED.value] += 1
        if introns_inferred:
            t.qc["inferred_intron_fragments"] += 1

        # ── crossings, per contiguous SEGMENT of the path ─────────────────────────────────────────
        # A line is crossed iff it lies strictly inside a segment, so per segment the crossed lines are a
        # contiguous index range and no container is needed. A node is SPANNED iff ONE segment crosses
        # both of its lines — not merely "both lines crossed", which would count a node the fragment
        # JUMPS OVER, whose two lines are touched by the two flanking segments from opposite sides.
        edge_count, edge_density = (
            (t.edge_spliced_count, t.edge_spliced_density)
            if spliced
            else (t.edge_unspliced_count, t.edge_unspliced_density)
        )
        quantum_edge = density_quantum(length - 1) if length >= 2 else 0
        quantum_node = density_quantum(length)
        n_crossed, sole_line = 0, -1
        for seg_start, seg_end in segments:
            first = int(np.searchsorted(cuts, seg_start, side="right"))
            last = int(np.searchsorted(cuts, seg_end, side="left"))
            for line in range(first, last):
                edge_count[edge_base + line - 1, channel] += 1
                edge_density[edge_base + line - 1, channel] += quantum_edge
            for line in range(first, last - 1):  # the node between two consecutive crossed lines
                t.node_spanning_count[node_base + line, channel] += 1
                t.node_spanning_density[node_base + line, channel] += quantum_node
            if last > first:
                sole_line = first if (n_crossed == 0 and last - first == 1) else -1
                n_crossed += last - first

        for jid in junction_ids:
            t.junction_count[jid, channel] += 1
            t.junction_density[jid, channel] += quantum_edge

        # ── contained: the WHOLE path lies inside ONE node ────────────────────────────────────────
        # ⚠ Not merely "crossed no line". An unannotated intron can swallow every line between two
        # blocks, leaving a fragment that crosses nothing yet straddles two nodes — crediting that as
        # contained would place its whole length in a node it only partly overlaps. Such a fragment is
        # neither contained nor crossing: it deposits on no object but is still counted (start_count),
        # so the loss is visible rather than silent.
        first_node = self._local_node(cuts, first_base)
        contained_node = -1
        if not junction_ids and first_node == self._local_node(cuts, last_base):
            contained_node = node_base + first_node
            t.node_contained_count[contained_node, channel] += 1
            t.node_contained_density[contained_node, channel] += quantum_node

        pool = self._pool(spliced, introns_inferred, contained_node, sole_line, node_base)
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

    def _junction_of(self, ref: int, intron_start: int, intron_end: int, motif_strand: int) -> int:
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
        for k in range(int(p.junction_offsets[donor]), int(p.junction_offsets[donor + 1])):
            if int(p.junction_acceptor_cut[k]) != acceptor:
                continue
            if motif_strand != Strand.NONE and int(p.junction_strand[k]) != int(motif_strand):
                continue
            return k
        return -1

    def _pool(self, spliced, inferred, contained_node, sole_line, node_base):
        """The one length pool this fragment belongs to, or ``None``.

        Priority, so that every pool stays pure: an OBSERVED splice is unambiguously RNA; a contained
        fragment is typed by its node; a single-line crossing is a "splash" read typed by its two flank
        types. Anything else — an exonic contained fragment, a multi-line crossing — is a mixture and
        enters nothing.
        """
        types = self.partition.node_types
        if spliced:
            return None if inferred else FragmentPool.RNA_SPLICED
        if contained_node >= 0:
            return _CONTAINED_POOL.get(int(types[contained_node]))
        if sole_line >= 0:
            flanks = (int(types[node_base + sole_line - 1]), int(types[node_base + sole_line]))
            return _SPLASH_POOL.get(tuple(sorted(flanks)))
        return None
