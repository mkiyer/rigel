"""THE ACCUMULATOR SPEC — the matrix the reference and the native build are both gated on.

    Design: ``docs/ACCUMULATOR_DESIGN.md``   ·   Plan: ``docs/IMPLEMENTATION_PLAN.md`` §3.6, §10.4

``_accumulator_reference.py`` is the executable specification; the native accumulator is required to
reproduce it byte for byte. This module is what "correct" means for both.

THE RULE UNDER TEST, in five lines. The genome is a graph: NODES are half-open intervals tiling each
reference, and the 0-bp LINES between adjacent nodes are CONTIGUOUS edges. A JUNCTION edge is a directed
donor→acceptor link from the annotation. A fragment is a PATH — its aligned blocks joined across mate
gaps and broken by introns — of length ``L = span − Σ intron``. Nodes count fragments CONTAINED and
SPANNING; edges count fragments CROSSING. Every object stores an integer count and a fixed-point density
``round(2^32 / placements)``, with ``placements = L`` at a node and ``L − 1`` at an edge.

⚠ **No partitioning.** Every crossed edge receives the FULL weight. The chance that a length-``L``
fragment crosses a given line is proportional to ``L`` and the deposit is ``1/L``, so the two cancel and
every fragment length contributes equally to each edge. Dividing by the number of edges crossed destroys
that cancellation and makes the answer depend on node spacing — measured up to **3.6× low**.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.types import Strand

from ._accumulator_reference import (
    DENSITY_SCALE,
    STRAND_COLUMNS,
    Accumulator,
    DepositOutcome,
    FragmentPool,
    Partition,
    density_quantum,
)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

#: chr1 cuts   0    100   200   201   400   900   1000
#: nodes        n0    n1    n2*   n3    n4    n5      (* n2 is 1 bp: [200,201))
#: lines           1     2     3     4     5          (local cut index)
CHR1_CUTS = [0, 100, 200, 201, 400, 900, 1000]
CHR2_CUTS = [0, 500, 1000]

#: coarse node type: 0 intergenic, 1 intron, 2 exon
CHR1_TYPES = [0, 2, 2, 1, 2, 0]
CHR2_TYPES = [0, 2]

#: an annotated intron whose endpoints are cuts 3 and 5, so it SWALLOWS the line at cut 4
JUNCTION = (0, 201, 900, Strand.POS)


def _partition(junctions=()):
    return Partition.from_cuts(
        [CHR1_CUTS, CHR2_CUTS], node_types=[CHR1_TYPES, CHR2_TYPES], junctions=junctions
    )


def _acc(junctions=(), **kw):
    return Accumulator(_partition(junctions), **kw)


def _edge(ref, line):
    """Global contiguous-edge id of the line at local cut index ``line``."""
    return (0 if ref == 0 else len(CHR1_CUTS) - 2) + line - 1


def _node(ref, local):
    return (0 if ref == 0 else len(CHR1_CUTS) - 1) + local


# ---------------------------------------------------------------------------
# the fixed-point density
# ---------------------------------------------------------------------------


def test_density_quantum_is_exact_and_rounds_half_away_from_zero():
    """The rounding mode is part of the contract — byte-identity is undefined without it, and Python's
    own ``round`` is banker's rounding, which differs at ties."""
    assert density_quantum(1) == DENSITY_SCALE
    assert density_quantum(2) == DENSITY_SCALE // 2
    assert density_quantum(512) * 512 == DENSITY_SCALE
    assert density_quantum(3) == (2 * DENSITY_SCALE + 3) // 6
    with pytest.raises(ValueError):
        density_quantum(0)


def test_one_fragment_recovers_its_own_reciprocal_length():
    acc = _acc()
    acc.deposit(0, 120, 320)
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 2), 0]) == 1
    assert int(t.edge_unspliced_density[_edge(0, 2), 0]) == density_quantum(200 - 1)


# ---------------------------------------------------------------------------
# contained / crossing / spanning
# ---------------------------------------------------------------------------


def test_a_contained_fragment_touches_ONE_node_and_no_edge():
    acc = _acc()
    assert acc.deposit(0, 220, 380) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 3), 0]) == 1
    assert int(t.node_contained_density[_node(0, 3), 0]) == density_quantum(160)
    assert t.node_contained_count.sum() == 1
    assert t.edge_unspliced_count.sum() == 0
    assert t.node_spanning_count.sum() == 0


def test_a_fragment_ENDING_AT_a_line_does_not_cross_it():
    """A line is 0 bp wide: crossing needs a base on each side."""
    acc = _acc()
    acc.deposit(0, 120, 200)
    assert acc.tally.edge_unspliced_count.sum() == 0
    assert int(acc.tally.node_contained_count[_node(0, 1), 0]) == 1


def test_a_fragment_STARTING_AT_a_line_does_not_cross_it():
    acc = _acc()
    acc.deposit(0, 201, 390)
    assert acc.tally.edge_unspliced_count.sum() == 0
    assert int(acc.tally.node_contained_count[_node(0, 3), 0]) == 1


def test_a_fragment_crossing_four_nodes_credits_exactly_THREE_edges_at_FULL_weight():
    acc = _acc()
    acc.deposit(0, 150, 500)  # touches n1 n2 n3 n4 -> lines at 200, 201, 400
    t = acc.tally
    assert [int(t.edge_unspliced_count[_edge(0, j), 0]) for j in (1, 2, 3, 4, 5)] == [0, 1, 1, 1, 0]
    quantum = density_quantum(350 - 1)
    assert all(int(t.edge_unspliced_density[_edge(0, j), 0]) == quantum for j in (2, 3, 4))
    assert t.node_contained_count.sum() == 0


def test_a_fragment_spanning_a_1bp_node_credits_BOTH_lines_and_the_spanning_slot():
    """1 bp nodes are legal — 15,687 of them at human scale — and nothing may assume length > 1."""
    acc = _acc()
    acc.deposit(0, 150, 300)
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 2), 0]) == 1
    assert int(t.edge_unspliced_count[_edge(0, 3), 0]) == 1
    assert int(t.node_spanning_count[_node(0, 2), 0]) == 1
    assert int(t.node_spanning_density[_node(0, 2), 0]) == density_quantum(150)
    assert int(t.node_contained_count[_node(0, 2), 0]) == 0


def test_contained_and_spanning_are_disjoint_and_L_equal_to_node_plus_one_is_in_NEITHER():
    """``contained`` needs ``L <= node``; ``spanning`` needs ``L >= node + 2``. The length between them
    crosses exactly one line and belongs to neither population."""
    acc = _acc()
    acc.deposit(0, 200, 202)  # n2 is 1 bp, L = 2
    t = acc.tally
    assert t.node_contained_count.sum() == 0
    assert t.node_spanning_count.sum() == 0
    assert int(t.edge_unspliced_count[_edge(0, 3), 0]) == 1


# ---------------------------------------------------------------------------
# splicing
# ---------------------------------------------------------------------------


def test_a_spliced_jump_deposits_NOTHING_on_the_lines_it_splices_over():
    """⭐ The defect this design removes. The intron [201,900) swallows the line at 400; the old rule,
    which asked only "does another slice follow?", could not tell that from a contiguous crossing."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    length = (950 - 150) - (900 - 201)
    assert int(t.sj_count[0, 0]) == 1
    assert int(t.sj_density[0, 0]) == density_quantum(length - 1)
    assert int(t.edge_unspliced_count[_edge(0, 4), 0]) == 0, "the swallowed line at 400"
    assert int(t.edge_spliced_count[_edge(0, 4), 0]) == 0
    assert int(t.node_spanning_count[_node(0, 4), 0]) == 0, "a node jumped OVER is not spanned"


def test_a_spliced_fragments_own_BLOCK_crossings_go_in_the_SPLICED_bank():
    """A spliced fragment is certified RNA — gDNA cannot be spliced — so a line its block genuinely
    crosses is the cleanest RNA marker available at a seam."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.edge_spliced_count[_edge(0, 2), 0]) == 1, "block [150,201) crosses the line at 200"
    assert t.edge_unspliced_count.sum() == 0


def test_a_node_is_SPANNED_only_when_ONE_segment_crosses_BOTH_its_lines():
    """⚠ Not "both lines crossed". Here an unannotated intron sits strictly inside n4 = [400,900), so
    segment 1 crosses the line at 400 and segment 2 crosses the line at 900 — both of n4's lines — yet
    the molecule never covers [500,600) and does not span the node."""
    acc = _acc()
    acc.deposit(0, 380, 950, introns=[(500, 600)])
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 4), 0]) == 1, "line 400, from segment 1"
    assert int(t.edge_unspliced_count[_edge(0, 5), 0]) == 1, "line 900, from segment 2"
    assert int(t.node_spanning_count[_node(0, 4), 0]) == 0, "but n4 is NOT spanned"


def test_an_UNANNOTATED_intron_credits_no_junction_and_nothing_across_the_gap():
    """Owner ruling: unannotated junctions are disproportionately artifactual, so they deposit on the
    UNSPLICED channel and compete with gDNA rather than being certified RNA."""
    acc = _acc(junctions=[JUNCTION])
    # [200,400) is NOT annotated; it swallows the line at 201
    acc.deposit(0, 50, 500, introns=[(200, 400)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.sj_count.sum() == 0
    assert t.edge_spliced_count.sum() == 0, "not certified RNA — it competes with gDNA"
    assert int(t.edge_unspliced_count[_edge(0, 1), 0]) == 1, (
        "block [50,200) crosses the line at 100"
    )
    assert int(t.edge_unspliced_count[_edge(0, 3), 0]) == 0, "the swallowed line at 201"
    assert t.qc["unannotated_introns"] == 1


def test_a_fragment_straddling_two_nodes_without_crossing_a_line_is_NOT_contained():
    """⚠ An unannotated intron can swallow every line between two blocks. The fragment then crosses
    nothing, yet it straddles two nodes — crediting it as *contained* would put its whole length in a
    node it only partly overlaps. It deposits on no object, and the start count is what keeps that
    visible rather than silent."""
    acc = _acc()
    acc.deposit(0, 120, 500, introns=[(200, 400)])  # blocks land in n1 and n4, crossing no line
    t = acc.tally
    assert t.edge_unspliced_count.sum() == 0
    assert t.node_contained_count.sum() == 0
    assert t.node_spanning_count.sum() == 0
    assert int(t.node_start_count.sum()) == 1, "still counted"


def test_an_unannotated_intron_inside_one_node_is_a_contained_unspliced_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, introns=[(300, 340)])
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 3), 0]) == 1
    assert int(t.node_contained_density[_node(0, 3), 0]) == density_quantum(180 - 40)
    assert t.qc["unannotated_introns"] == 1


def test_opposite_strand_junctions_at_the_same_coordinates_are_DISTINCT_edges():
    """Biologically impossible — splice motifs are not palindromic — so only a synthetic stress test can
    reach it, which is exactly why one exists."""
    acc = _acc(junctions=[(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)])
    acc.deposit(0, 150, 950, introns=[(201, 900)], strand=Strand.NEG, sj_strand=Strand.NEG)
    t = acc.tally
    assert t.sj_count.sum() == 1
    assert int(t.sj_count[1, STRAND_COLUMNS[Strand.NEG]]) == 1, "the NEG edge (id 1), genome minus"


@pytest.mark.parametrize("order", [("POS_first", 1), ("NEG_first", -1)])
def test_a_junction_id_is_a_function_of_the_PARTITION_not_of_argument_order(order):
    """⛔ The junction-edge id IS the rank in the sort, so the sort must be total — otherwise the id
    depends on the order the caller happened to list the junctions in, and the same graph gets two
    labellings.

    ⚠ This is the ONLY test that pins ``strand`` as part of the sort key. ``np.lexsort`` is stable, so a
    key of ``(acceptor, donor)`` alone gives the right answer for any input whose ties already arrive in
    the right order — which is every other test, and both real indexes. Reversing the argument order is
    what makes the missing key observable, and a strand-coincident pair is the only tie there is.
    """
    _, direction = order
    junctions = [(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)][::direction]
    part = Partition.from_cuts([CHR1_CUTS], node_types=[CHR1_TYPES], junctions=junctions)
    assert [int(s) for s in part.sj_strand] == [int(Strand.POS), int(Strand.NEG)], (
        "POS must sort to slot 0 whichever order it was passed in"
    )


def test_a_fragment_using_TWO_junctions_credits_BOTH():
    """Owner ruling: each edge owns its own expectation, and the strand model is fitted from a separate
    scan output, so crediting every junction distorts nothing.

    ⚠ The two introns must be separated by a real exon. Abutting introns imply a zero-length exon and are
    malformed (see ``test_ABUTTING_introns_are_MALFORMED_and_merge``), so this needs its own partition
    with room for an exon between them."""
    part = Partition.from_cuts(
        [[0, 100, 200, 300, 400, 500, 600]],
        node_types=[[0, 2, 1, 2, 1, 2]],
        junctions=[(0, 100, 200, Strand.POS), (0, 300, 400, Strand.POS)],
    )
    acc = Accumulator(part, max_fragment_length=10_000)
    acc.deposit(0, 50, 550, introns=[(100, 200), (300, 400)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.sj_count[1, STRAND_COLUMNS[Strand.POS]]) == 1
    assert t.qc["introns_absorbed"] == 0, "a real exon separates them; nothing is malformed"


# ---------------------------------------------------------------------------
# strand channels
# ---------------------------------------------------------------------------


def test_the_unspliced_bank_is_indexed_by_GENOME_strand():
    acc = _acc()
    acc.deposit(0, 150, 300, strand=Strand.NEG)
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 2), 1]) == 1
    assert int(t.edge_unspliced_count[_edge(0, 2), 0]) == 0


def test_EVERY_bank_including_the_junctions_is_indexed_by_GENOME_strand():
    """⭐ One convention throughout. Sense/antisense is DERIVED, never stored: the junction edge carries
    its own genomic strand, so a consumer computes ``sense = (fragment strand == junction strand)``.

    Here a genome-minus fragment splices across a ``+`` junction. Under a sense convention it would land
    in the antisense column; under the genome convention it lands in the minus column. Those happen to be
    the same index, so the discriminating case is the next test."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, introns=[(201, 900)], strand=Strand.NEG, sj_strand=Strand.POS)
    assert int(acc.tally.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1


def test_a_SENSE_fragment_on_the_minus_strand_is_still_booked_as_MINUS():
    """The discriminating case: sense-to-motif would say column 0, genome strand says column 1."""
    acc = _acc(junctions=[(0, 201, 900, Strand.NEG)])
    acc.deposit(0, 150, 950, introns=[(201, 900)], strand=Strand.NEG, sj_strand=Strand.NEG)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.edge_spliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1


# ---------------------------------------------------------------------------
# bounded influence, clipping, references
# ---------------------------------------------------------------------------


def test_a_fragment_over_the_length_limit_deposits_NOTHING_and_is_COUNTED():
    """⭐ Bounded influence. Unbounded, 1,000 read groups own 99.8 % of all edge crossings on a real
    library; with the limit on ``L`` they own 4.16 %. A silent drop would hide that."""
    acc = _acc(max_fragment_length=200)
    assert acc.deposit(0, 100, 500) is DepositOutcome.TOO_LONG
    t = acc.tally
    assert t.edge_unspliced_count.sum() == 0
    assert t.node_start_count.sum() == 0
    assert t.qc["dropped_too_long"] == 1


def test_the_limit_applies_to_L_and_NOT_to_the_SPAN():
    """⚠ A 300 bp molecule across a 10 kb intron has a 10 kb span. Limiting the span discards every
    spliced fragment — 37.96 % of read groups measured, against 5.45 % when the limit is on ``L``."""
    acc = _acc(junctions=[JUNCTION], max_fragment_length=200)
    out = acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.POS)
    assert out is DepositOutcome.DEPOSITED, "span 800, L = 101"
    assert int(acc.tally.sj_count[0, 0]) == 1


def test_a_fragment_is_clipped_to_its_reference_and_L_is_the_clipped_length():
    acc = _acc()
    acc.deposit(0, 950, 1200)  # chr1 ends at 1000
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 5), 0]) == 1
    assert int(t.node_contained_density[_node(0, 5), 0]) == density_quantum(50)


def test_a_single_node_reference_has_no_edges_and_still_accepts_a_fragment():
    acc = Accumulator(Partition.from_cuts([[0, 1000]], node_types=[[0]]))
    assert acc.n_edges == 0
    acc.deposit(0, 100, 300)
    assert int(acc.tally.node_contained_count[0, 0]) == 1


def test_the_per_reference_offsets_do_not_bleed():
    """chr1's fragment crosses its lines 2 and 3; chr2's crosses chr2's line 1. Nothing lands on a
    reference it did not come from — the failure mode that once dropped 476,719 of 476,732 fragments."""
    acc = _acc()
    acc.deposit(0, 150, 300)  # crosses the lines at 200 AND 201
    acc.deposit(1, 400, 700)  # crosses chr2's line at 500
    t = acc.tally
    assert [int(t.edge_unspliced_count[e, 0]) for e in range(acc.n_edges)] == [0, 1, 1, 0, 0, 1]
    assert int(t.edge_unspliced_count.sum()) == 3
    assert int(t.node_start_count[_node(0, 1)]) == 1
    assert int(t.node_start_count[_node(1, 0)]) == 1


# ---------------------------------------------------------------------------
# the invariant that can actually fire
# ---------------------------------------------------------------------------


def test_every_accepted_fragment_increments_exactly_ONE_start_count():
    """⚠ The crossing and contained totals are tautologies — they can only be evaluated by re-running
    the deposit. This one is checkable against a number the scanner knows independently."""
    acc = _acc(junctions=[JUNCTION])
    fragments = [(120, 320, ()), (220, 380, ()), (150, 950, [(201, 900)]), (950, 1200, ())]
    accepted = sum(
        acc.deposit(0, s, e, introns=i, sj_strand=Strand.POS) is DepositOutcome.DEPOSITED
        for s, e, i in fragments
    )
    assert accepted == 4
    assert int(acc.tally.node_start_count.sum()) == 4
    assert acc.tally.qc["deposited"] == 4


# ---------------------------------------------------------------------------
# the five fragment-length pools
# ---------------------------------------------------------------------------


def test_each_pool_is_reached_only_by_its_own_structural_class():
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 10, 90)  # contained in n0 — intergenic
    acc.deposit(0, 210, 390)  # contained in n3 — intronic
    acc.deposit(0, 380, 420)  # crosses the line at 400 only — flanks intron|exon
    acc.deposit(0, 950, 990)  # contained in n5 — intergenic
    acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.POS)  # annotated junction
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTERGENIC].sum()) == 2
    assert int(p[FragmentPool.DNA_INTRONIC].sum()) == 1
    assert int(p[FragmentPool.DNA_INTRON_EXON].sum()) == 1
    assert int(p[FragmentPool.RNA_SPLICED].sum()) == 1
    assert int(p[FragmentPool.DNA_INTERGENIC_EXON].sum()) == 0


def test_a_pool_is_binned_at_L_and_only_ONCE_per_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, introns=[(300, 340)])  # L = 180 - 40
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTRONIC, 140]) == 1
    assert int(p.sum()) == 1


def test_an_INFERRED_intron_is_kept_OUT_of_the_pure_RNA_pool():
    """Its splice is a model inference, not an observation, so certifying it as RNA would make the pool
    depend on the very length model it is used to fit."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.POS, introns_inferred=True)
    t = acc.tally
    assert int(t.sj_count[0, 0]) == 1, "it still deposits"
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 0
    assert t.qc["inferred_intron_fragments"] == 1


def test_a_multi_line_crossing_enters_NO_pool():
    """A splash read straddles ONE probe edge. A fragment crossing several lines has no single
    structural class, and an impure pool is worse than a missing one."""
    acc = _acc()
    acc.deposit(0, 150, 500)
    assert int(acc.tally.pool_lengths.sum()) == 0


# ---------------------------------------------------------------------------
# the two arithmetic gates — these are what make the design true
# ---------------------------------------------------------------------------


def _corpus(rng, n, ref_len):
    lengths = np.clip(rng.normal(200.0, 60.0, n).round().astype(np.int64), 30, 600)
    starts = rng.integers(0, ref_len - lengths)
    return starts, starts + lengths, lengths


def _uniform_accumulator(node_bp, ref_len):
    cuts = list(range(0, ref_len + 1, node_bp))
    return Accumulator(
        Partition.from_cuts([cuts], node_types=[[0] * (len(cuts) - 1)]), max_fragment_length=10_000
    )


@pytest.mark.parametrize("node_bp", [50, 200, 1000])
def test_the_crossing_DENSITY_recovers_the_true_density_with_NO_length_model(node_bp):
    """⭐ ``E[Σ 1/(L−1)] = ρ`` exactly, for ANY fragment-length distribution. This is the identity the
    whole design rests on and the reason no divisor and no length model appear at an edge.

    It must hold at every node spacing. Partitioning the weight by the number of edges crossed breaks
    it — measured 0.28× at 50 bp nodes, 0.54× at 100 bp, 0.91× at 200 bp — so this test is also what
    forbids partitioning.
    """
    ref_len, rho = 200_000, 0.05
    acc = _uniform_accumulator(node_bp, ref_len)
    rng = np.random.default_rng(7)
    starts, ends, _ = _corpus(rng, int(rho * ref_len), ref_len)
    for s, e in zip(starts, ends):
        acc.deposit(0, int(s), int(e))
    interior = slice(5, acc.n_edges - 5)
    estimate = (
        acc.tally.edge_unspliced_density[interior, :].sum() / DENSITY_SCALE / (acc.n_edges - 10)
    )
    assert 0.98 <= estimate / rho <= 1.02, f"{estimate / rho:.4f} at {node_bp} bp nodes"


def test_the_crossing_COUNT_recovers_density_times_mean_length():
    """The companion identity ``E[count] = ρ·(E[L] − 1)``. Together with the line above, this is the 2×2
    that separates gDNA from RNA by fragment length alone."""
    ref_len, rho = 200_000, 0.05
    acc = _uniform_accumulator(200, ref_len)
    rng = np.random.default_rng(11)
    starts, ends, lengths = _corpus(rng, int(rho * ref_len), ref_len)
    for s, e in zip(starts, ends):
        acc.deposit(0, int(s), int(e))
    interior = slice(5, acc.n_edges - 5)
    per_edge = acc.tally.edge_unspliced_count[interior, :].sum() / (acc.n_edges - 10)
    expected = rho * (lengths.mean() - 1.0)
    assert 0.98 <= per_edge / expected <= 1.02, f"{per_edge / expected:.4f}"


def test_the_deposit_is_independent_of_the_ORDER_fragments_arrive_in():
    """Integer addition is associative, which is what makes the per-worker merge bit-identical at any
    thread count — the property the float32 mass channels could not offer."""
    rng = np.random.default_rng(3)
    starts, ends, _ = _corpus(rng, 400, 900)
    order = rng.permutation(len(starts))

    def run(idx):
        acc = _acc()
        for k in idx:
            acc.deposit(0, int(starts[k]), int(ends[k]))
        return acc.tally

    a, b = run(range(len(starts))), run(order)
    for field in (
        "node_contained_count",
        "node_contained_density",
        "node_spanning_count",
        "node_spanning_density",
        "node_start_count",
        "edge_unspliced_count",
        "edge_unspliced_density",
        "pool_lengths",
    ):
        assert np.array_equal(getattr(a, field), getattr(b, field)), field


# ---------------------------------------------------------------------------
# malformed intron lists, and the ONE definition of L
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "introns,expected_length,expected_crossings,expected_absorbed",
    [
        ([(210, 260), (240, 300)], 360, 4, 1),  # overlapping
        ([(200, 400), (250, 300)], 250, 1, 1),  # nested
        ([(150, 480), (160, 470)], 120, 1, 1),  # wide overlap: the naive formula went NEGATIVE
        ([(210, 260), (210, 260)], 400, 4, 1),  # duplicated, as two disagreeing mates produce
        ([(300, 300)], 450, 4, 1),  # zero length
        ([(210, 260), (300, 340)], 360, 4, 0),  # well formed: nothing absorbed
    ],
)
def test_L_is_the_total_of_the_path_segments_even_when_the_intron_list_is_malformed(
    introns, expected_length, expected_crossings, expected_absorbed
):
    """⚠ ONE definition of ``L``: the total length of the path's segments.

    Computing it separately as ``span − Σ(intron lengths)`` is a *second* formula for the same quantity,
    and the two disagree the moment introns overlap — by up to 1.5×, and on a wide overlap the second one
    goes NEGATIVE, so a good fragment is silently discarded. A real BAM produces this: the scanner reads
    the ``XS`` tag once per record, so a pair whose mates disagree about an acceptor yields overlapping
    introns for one molecule (measured: 1 read group in 875,670 on MO_3021).

    It matters far beyond that rate, because a C++ author builds the segments first and will naturally
    sum them — so byte-identity, the only gate this design has, would come down to which of two
    contradictory rules the file happened to carry.
    """
    acc = _acc(max_fragment_length=10_000)
    assert acc.deposit(0, 50, 500, introns=introns) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert t.qc["introns_absorbed"] == expected_absorbed
    crossings = int(t.edge_unspliced_count.sum())
    assert crossings == expected_crossings
    assert int(t.edge_unspliced_density.sum()) == crossings * density_quantum(expected_length - 1)


def test_the_path_STARTS_where_its_first_covered_base_is_not_where_the_extent_begins():
    """A leading intron means the molecule does not begin at ``lo``. Attributing it to the node
    containing ``lo`` would credit the start-count invariant — and possibly the contained deposit — to a
    node the fragment never touches."""
    acc = _acc(max_fragment_length=10_000)
    acc.deposit(0, 150, 500, introns=[(150, 480)])  # the path is only [480,500), inside n4
    t = acc.tally
    assert int(t.node_start_count[_node(0, 4)]) == 1, "n4, where the path actually starts"
    assert int(t.node_start_count[_node(0, 1)]) == 0, "not n1, where the extent begins"
    assert int(t.node_contained_count[_node(0, 4), 0]) == 1
    assert int(t.node_contained_density[_node(0, 4), 0]) == density_quantum(20)


def test_a_duplicated_intron_credits_its_junction_ONCE():
    """Two mates reporting the same intron is one splice event, not two."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, introns=[(201, 900), (201, 900)], sj_strand=Strand.POS)
    assert int(acc.tally.sj_count[0, 0]) == 1
    assert acc.tally.qc["introns_absorbed"] == 1


def test_ABUTTING_introns_are_MALFORMED_and_merge():
    """⚠ Two introns sharing an endpoint imply a **zero-length exon** between them, which is physically
    impossible — a transcript with one is molecularly identical to a transcript without it. So a single
    molecule can never legitimately use both, and the pair is an alignment artifact.

    The index cannot produce it either: a zero-length exon is dropped when the exon arrays are built,
    which fuses its two flanking introns into one. Merged here, and counted."""
    acc = _acc(junctions=[(0, 201, 400, Strand.POS), (0, 400, 900, Strand.POS)])
    acc.deposit(0, 150, 950, introns=[(201, 400), (400, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.qc["introns_absorbed"] == 1
    assert t.sj_count.sum() == 0, "the merged span 201->900 is not an annotated junction"


def test_a_wide_overlap_no_longer_discards_a_good_fragment():
    """The naive formula gave L = −290 here and filed the fragment as ``dropped_empty`` — invisible to
    the start-count invariant, because a rejected fragment never reaches it."""
    acc = _acc(max_fragment_length=10_000)
    assert acc.deposit(0, 150, 500, introns=[(150, 480), (160, 470)]) is DepositOutcome.DEPOSITED
    assert acc.tally.qc["dropped_empty"] == 0
    assert int(acc.tally.node_start_count.sum()) == 1


# ---------------------------------------------------------------------------
# the node banks carry ONE strand convention
# ---------------------------------------------------------------------------

#: a junction far enough right that the first block still spans node n2 = [200,201)
SPAN_JUNCTION_POS = (0, 400, 900, Strand.POS)
SPAN_JUNCTION_NEG = (0, 400, 900, Strand.NEG)


def test_a_spliced_and_an_unspliced_fragment_of_the_SAME_genome_strand_share_a_column():
    """⚠ One array, one convention.

    A spliced fragment cannot be *contained* — both endpoints of an annotated intron are cuts, so it
    always crosses its junction edge — but its blocks routinely SPAN a node whole. Measured on real
    cfRNA, **65–69 % of all node_spanning deposits come from spliced fragments**. Indexing those by
    sense-relative-to-motif while the unspliced ones beside them use genome strand would put one array
    into two conventions, and 40–44 % of the spliced deposits would land in the opposite column from
    their unspliced neighbours.
    """
    acc = _acc(junctions=[SPAN_JUNCTION_POS])
    acc.deposit(0, 150, 300, strand=Strand.NEG)  # unspliced, genome minus
    acc.deposit(  # spliced, genome minus, ANTISENSE to its + junction
        0, 150, 950, introns=[(400, 900)], strand=Strand.NEG, sj_strand=Strand.POS
    )
    t = acc.tally
    assert int(t.node_spanning_count[_node(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 2, (
        "both genome minus"
    )
    assert int(t.node_spanning_count[_node(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1, "the junction bank too"


def test_a_spliced_SENSE_fragment_books_node_AND_junction_by_GENOME_strand():
    """The discriminating case: sense-to-motif would say column 0 for both; genome strand says 1."""
    acc = _acc(junctions=[SPAN_JUNCTION_NEG])
    acc.deposit(0, 150, 950, introns=[(400, 900)], strand=Strand.NEG, sj_strand=Strand.NEG)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1
    assert int(t.node_spanning_count[_node(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1


# ---------------------------------------------------------------------------
# a strand that is not a strand
# ---------------------------------------------------------------------------
#
# ⛔ Both cases below were found by an adversarial read of this spec before any C++ followed it, and both
# were confirmed by execution. They exist because `deposit` inherited no equivalent of the scanner's gate
# at ``bam_scanner.cpp:1474-1480``, which used to reject an undefined strand before the old accumulator
# ever saw it.


@pytest.mark.parametrize("undefined", [Strand.NONE, Strand.AMBIGUOUS])
def test_an_UNDEFINED_strand_is_REJECTED_not_silently_booked_as_MINUS(undefined):
    """⛔ The channel IS the genome strand, so a fragment without one has no channel.

    ``genome_channel`` is ``STRAND_COLUMNS[Strand.POS] if strand == POS else STRAND_COLUMNS[Strand.NEG]``, so every strand that is
    not POS — including NONE and AMBIGUOUS — used to land in the MINUS column. That is not a dropped
    fragment, it is a fragment **credited to the wrong strand**: exactly the class of error the one-strand
    convention exists to delete.

    ⚠ AMBIGUOUS is reachable in production, not hypothetical: ``build_fragment`` keys blocks by
    ``(ref, strand)``, so mates agreeing in reference orientation give ``strand = POS|NEG``.
    The design already requires the count (§10.3 lists strand-undefined fragments as a QC denominator the
    accumulator must emit); it was simply missing.
    """
    acc = _acc()
    assert acc.deposit(0, 150, 300, strand=undefined) is DepositOutcome.STRAND_UNDEFINED
    t = acc.tally
    assert t.qc["dropped_strand_undefined"] == 1
    assert int(t.node_contained_count.sum()) == 0, "must not be booked into either column"
    assert int(t.node_start_count.sum()) == 0, "a rejected fragment never reaches the invariant"


# ── the splice-junction motif strand is a THREE-WAY distinction, not a two-way one ─────────────────
#
# ``sj_strand`` is the OR of the per-record ``XS``/``ts`` tag values, so it carries three
# genuinely different states and each needs its own rule. The original spec had only two, treating
# AMBIGUOUS as though it were a definite strand:
#
#   NONE       no strand tag in the BAM at all  ->  no information; match on coordinates alone
#   POS / NEG  one definite observed strand      ->  must agree with the junction edge's own strand
#   AMBIGUOUS  the two mates' tags DISAGREE      ->  contradictory evidence; trust no splice
#
# ⚠ NONE must stay permissive. Aligners differ — STAR writes ``XS``, minimap2 writes ``ts``, and some
# write neither — so on an untagged BAM every spliced fragment arrives with NONE. Requiring a strand
# there would delete the entire spliced-RNA signal for that aligner.


def test_a_MISSING_sj_strand_MATCHES_on_coordinates_alone():
    """⛔ The case that makes untagged aligners work at all, so it is pinned before the two below.

    An aligner that writes neither ``XS`` nor ``ts`` gives every spliced fragment
    ``sj_strand = NONE``. If the junction lookup demanded a strand, that BAM would lose 100 % of its
    annotated junctions — and the loss would look like a stale annotation, not a convention bug.
    """
    acc = _acc(junctions=[JUNCTION])  # (0, 201, 900, POS)
    acc.deposit(0, 150, 950, introns=[(201, 900)], sj_strand=Strand.NONE)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 1
    assert t.qc["unannotated_introns"] == 0
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 1, "certified RNA"


def test_an_AMBIGUOUS_sj_strand_is_CONTRADICTORY_and_credits_NO_junction():
    """⛔ AMBIGUOUS is contradictory evidence, not missing evidence — and it is neither of the two things
    the original rule could express.

    ``sj_strand`` is the OR of a per-RECORD tag, so AMBIGUOUS (``POS | NEG``) means **the two mates
    disagreed about the same molecule**. That is a data-quality signal of the same family as mates agreeing
    in reference orientation, so the splice must not be trusted: no junction is credited and the fragment
    deposits on the unspliced channel, which is the safe direction the design already takes for
    unannotated junctions.

    ⚠ It must NOT be counted as an unannotated intron. That counter's whole purpose is measuring annotation
    coverage, so feeding it alignment disagreements makes the metric report a stale annotation whenever the
    aligner is inconsistent. It gets its own denominator.

    ⚠ Reachable today with no change of ours: ``collect_implicit_splice_introns`` stamps each PE gap's
    intron with the first matching candidate transcript's strand and the caller ORs them, so a two-gap
    fragment matching opposite-strand transcripts arrives here as AMBIGUOUS.
    """
    acc = _acc(junctions=[JUNCTION])  # (0, 201, 900, POS)
    outcome = acc.deposit(
        0, 150, 950, introns=[(201, 900)], strand=Strand.POS, sj_strand=Strand.AMBIGUOUS
    )
    t = acc.tally
    assert outcome is DepositOutcome.DEPOSITED, "the fragment is real; only its splice is untrusted"
    assert int(t.sj_count.sum()) == 0
    assert t.qc["contradictory_sj_strand"] == 1
    assert t.qc["unannotated_introns"] == 0, "the annotation-coverage metric must not be poisoned"
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 0, "not certified RNA"
    assert int(t.edge_unspliced_count.sum()) > 0, "it competes with gDNA, the safe direction"


def test_a_DEFINITE_but_WRONG_sj_strand_still_misses():
    """The third arm, so the rule above cannot be over-applied: a definite strand that disagrees with the
    junction edge's own strand is a real disagreement, and it IS an unannotated intron — that coordinate
    pair is not annotated on the strand it was observed on."""
    acc = _acc(junctions=[JUNCTION])  # (0, 201, 900, POS)
    acc.deposit(0, 150, 950, introns=[(201, 900)], strand=Strand.NEG, sj_strand=Strand.NEG)
    t = acc.tally
    assert int(t.sj_count.sum()) == 0
    assert t.qc["unannotated_introns"] == 1
    assert t.qc["contradictory_sj_strand"] == 0
    assert int(t.node_spanning_count[_node(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
