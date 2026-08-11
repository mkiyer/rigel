"""THE ACCUMULATOR SPEC — the matrix the reference and the native build are both gated on.

       §10.4

``_accumulator_reference.py`` is the executable specification; the native accumulator is required to
reproduce it byte for byte. This module is what "correct" means for both.

THE RULE UNDER TEST, in five lines. The genome is a graph: NODES are half-open intervals tiling each
reference, and the 0-bp LINES between adjacent nodes are CONTIGUOUS edges. A JUNCTION edge is a directed
donor→acceptor link from the annotation. A fragment is a PATH — its aligned blocks joined across mate
gaps and broken by introns — of length ``L = span − Σ intron``. Nodes count fragments CONTAINED; edges
count fragments CROSSING. Each population stores only the channels something READS — an integer count, a
fixed-point ``round(2^32 / placements)`` with ``placements = L`` at a node and ``L − 1`` at an edge, a
``Σ L``, and on the contiguous edges the CONSERVED MASS, which sums to one per fragment.

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
    INV_LENGTH_SCALE,
    STRAND_COLUMNS,
    Accumulator,
    DepositOutcome,
    FragmentPool,
    GapHypothesis,
    Partition,
    inv_length_quantum,
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



def _contained_quantum(ref, local, length):
    """The deposit a CONTAINED length-``length`` fragment makes on node ``(ref, local)``.

    ⭐⭐ **The deposit is ``1/OPPORTUNITY``, not ``1/length``** — a length-`w` fragment inside a node of
    length `ell` had `ell − w + 1` admissible start positions, so `1/(ell − w + 1)` cancels the
    opportunity identically and the channel is a DENSITY for any length distribution
    (`test_fragment_length_proof.test_the_node_deposit_is_the_RECIPROCAL_OPPORTUNITY_...`). ⚠ Derived from
    the fixture's own cuts rather than written as a number, so an assertion states the RULE.
    """
    cuts = CHR1_CUTS if ref == 0 else CHR2_CUTS
    return inv_length_quantum(cuts[local + 1] - cuts[local] - length + 1)


def _node(ref, local):
    return (0 if ref == 0 else len(CHR1_CUTS) - 1) + local


# ---------------------------------------------------------------------------
# the fixed-point density
# ---------------------------------------------------------------------------


def test_inv_length_quantum_is_exact_and_rounds_half_away_from_zero():
    """The rounding mode is part of the contract — byte-identity is undefined without it, and Python's
    own ``round`` is banker's rounding, which differs at ties."""
    assert inv_length_quantum(1) == INV_LENGTH_SCALE
    assert inv_length_quantum(2) == INV_LENGTH_SCALE // 2
    assert inv_length_quantum(512) * 512 == INV_LENGTH_SCALE
    assert inv_length_quantum(3) == (2 * INV_LENGTH_SCALE + 3) // 6
    with pytest.raises(ValueError):
        inv_length_quantum(0)


def test_one_fragment_recovers_its_own_reciprocal_length():
    acc = _acc()
    acc.deposit(0, 120, 320)
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 2), 0]) == 1
    assert int(t.edge_unspliced_inv_length_sum[_edge(0, 2)]) == inv_length_quantum(200 - 1)


# ---------------------------------------------------------------------------
# contained / crossing
# ---------------------------------------------------------------------------


def test_a_contained_fragment_touches_ONE_node_and_no_edge():
    acc = _acc()
    assert acc.deposit(0, 220, 380) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 3), 0]) == 1
    assert int(t.node_contained_inv_opportunity_sum[_node(0, 3)]) == _contained_quantum(0, 3, 160)
    assert t.node_contained_count.sum() == 1
    assert t.edge_unspliced_count.sum() == 0


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
    quantum = inv_length_quantum(350 - 1)
    assert all(int(t.edge_unspliced_inv_length_sum[_edge(0, j)]) == quantum for j in (2, 3, 4))
    assert t.node_contained_count.sum() == 0


def test_a_fragment_covering_a_1bp_node_credits_BOTH_lines_and_conserves_its_mass():
    """1 bp nodes are legal — 15,687 of them at human scale — and nothing may assume length > 1.

    ⭐ The 1 bp node is the sharpest case for the CONSERVED MASS: its slice is a single base shared
    between two bounding lines, so the fragment's three slices are 50 / 1 / 99 bases of a 150 bp
    molecule and must still sum to exactly one fragment.
    ⚠ This test used to assert a ``node_spanning`` deposit here. That bank was removed on evidence, and
    the mass is what now makes the case observable — a strictly stronger statement, since a count of 1
    survives any error in *how much* of the fragment was attributed."""
    acc = _acc()
    acc.deposit(0, 150, 300)
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 2), 0]) == 1
    assert int(t.edge_unspliced_count[_edge(0, 3), 0]) == 1
    assert int(t.node_contained_count[_node(0, 2), 0]) == 0
    mass = int(t.edge_unspliced_mass.sum())
    assert abs(mass - INV_LENGTH_SCALE) <= 2, (
        f"the fragment deposited {mass / INV_LENGTH_SCALE} fragments, not 1"
    )


def test_a_fragment_LONGER_than_a_node_is_not_contained_and_crosses_exactly_one_line():
    """``contained`` needs ``L <= node``. One base longer and the fragment is a CROSSING instead, and it
    touches the node axis nowhere at all — that axis has exactly one population now."""
    acc = _acc()
    acc.deposit(0, 200, 202)  # n2 is 1 bp, L = 2
    t = acc.tally
    assert t.node_contained_count.sum() == 0
    assert int(t.edge_unspliced_count[_edge(0, 3), 0]) == 1


# ---------------------------------------------------------------------------
# splicing
# ---------------------------------------------------------------------------


def test_a_spliced_jump_deposits_NOTHING_on_the_lines_it_splices_over():
    """⭐ The defect this design removes. The intron [201,900) swallows the line at 400; the old rule,
    which asked only "does another slice follow?", could not tell that from a contiguous crossing."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    length = (950 - 150) - (900 - 201)
    assert int(t.sj_count[0, 0]) == 1
    assert int(t.sj_inv_length_sum[0]) == inv_length_quantum(length - 1)
    assert int(t.edge_unspliced_count[_edge(0, 4), 0]) == 0, "the swallowed line at 400"
    assert int(t.edge_spliced_count[_edge(0, 4), 0]) == 0


def test_a_spliced_fragments_own_BLOCK_crossings_go_in_the_SPLICED_bank():
    """A spliced fragment is certified RNA — gDNA cannot be spliced — so a line its block genuinely
    crosses is the cleanest RNA marker available at a seam."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.edge_spliced_count[_edge(0, 2), 0]) == 1, "block [150,201) crosses the line at 200"
    assert t.edge_unspliced_count.sum() == 0


def test_a_MULTI_SEGMENT_unspliced_fragment_conserves_its_mass_across_BOTH_segments():
    """An unannotated intron sits strictly inside n4 = [400,900), so segment 1 crosses the line at 400
    and segment 2 crosses the line at 900 — two segments, each touching a line.

    ⭐ Every segment touches a line, so the conservation law applies unchanged and the mass sums to
    exactly one fragment. ⚠ This is the case the law's "deposited, unspliced, annotated" wording exists
    for: had either segment touched NO line, its bases would have had nowhere conserved to go and the
    total would be a PARTIAL."""
    acc = _acc()
    acc.deposit(0, 380, 950, observed_introns=[(500, 600)])
    t = acc.tally
    assert int(t.edge_unspliced_count[_edge(0, 4), 0]) == 1, "line 400, from segment 1"
    assert int(t.edge_unspliced_count[_edge(0, 5), 0]) == 1, "line 900, from segment 2"
    mass = int(t.edge_unspliced_mass.sum())
    assert abs(mass - INV_LENGTH_SCALE) <= 2, (
        f"a two-segment fragment deposited {mass / INV_LENGTH_SCALE} fragments, not 1"
    )


def test_an_UNANNOTATED_intron_credits_no_junction_and_nothing_across_the_gap():
    """Owner ruling: unannotated junctions are disproportionately artifactual, so they deposit on the
    UNSPLICED channel and compete with gDNA rather than being certified RNA."""
    acc = _acc(junctions=[JUNCTION])
    # [200,400) is NOT annotated; it swallows the line at 201
    acc.deposit(0, 50, 500, observed_introns=[(200, 400)], sj_strand=Strand.POS)
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
    acc.deposit(
        0, 120, 500, observed_introns=[(200, 400)]
    )  # blocks land in n1 and n4, crossing no line
    t = acc.tally
    assert t.edge_unspliced_count.sum() == 0
    assert t.node_contained_count.sum() == 0
    assert int(t.node_start_count.sum()) == 1, "still counted"


def test_an_unannotated_intron_inside_one_node_is_a_contained_unspliced_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, observed_introns=[(300, 340)])
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 3), 0]) == 1
    assert int(t.node_contained_inv_opportunity_sum[_node(0, 3)]) == _contained_quantum(0, 3, 180 - 40)
    assert t.qc["unannotated_introns"] == 1


def test_opposite_strand_junctions_at_the_same_coordinates_are_DISTINCT_edges():
    """Biologically impossible — splice motifs are not palindromic — so only a synthetic stress test can
    reach it, which is exactly why one exists."""
    acc = _acc(junctions=[(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
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
    acc.deposit(0, 50, 550, observed_introns=[(100, 200), (300, 400)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.sj_count[1, STRAND_COLUMNS[Strand.POS]]) == 1
    assert t.qc["introns_absorbed"] == 0, "a real exon separates them; nothing is malformed"


# ---------------------------------------------------------------------------
# strand channels
# ---------------------------------------------------------------------------


def test_the_unspliced_bank_is_indexed_by_GENOME_strand():
    acc = _acc()
    acc.deposit(0, 150, 300, align_strand=Strand.NEG)
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
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    assert int(acc.tally.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1


def test_the_junction_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION():
    """⛔⛔ **DO NOT COLLAPSE ``sj_count`` TO ONE COLUMN.** A junction is stranded by its genomic splicing
    MOTIF, so the strand of the *fragments* on it looks redundant — and every consumer today sums the two
    columns. It is retained anyway, on an owner ruling (2026-08-08), and this test is why.

    ⭐ **THE MECHANISM.** Aligners emit false-positive ``N`` CIGAR ops from plain genomic DNA.
    ``rigel.splice_blacklist`` catches the ones the sister tool ``alignable`` has already enumerated by
    coordinate — an a-priori list, and far from complete. A second, EMPIRICAL detector exists in a
    stranded library: real junctions inherit the library's global strand specificity, while alignment
    artifacts deposit onto BOTH strands and deviate from it. That test is only possible if the
    per-junction split by ALIGNED strand survives into the payload.

    ⚠ Unstranded data cannot use it — with κ = ½ there is no expectation to deviate from. That is a
    property of the detector, not a reason to drop the column.

    ⭐⭐ **AND THE DISCRIMINATING INFORMATION LIVES ONLY IN THE SPLIT**, which is the whole ruling: the
    clean junction and the artifactual one below carry the SAME total, so a collapsed bank cannot tell
    them apart at all. A tidying pass that removed the column — the same "store a channel where a named
    consumer reads it" principle that correctly removed six other banks — would delete this detector
    before it was built.
    """
    # a CLEAN junction: every fragment on one aligned strand, as a stranded library produces
    clean = _acc(junctions=[JUNCTION])
    for _ in range(20):
        clean.deposit(
            0, 150, 950, observed_introns=[(201, 900)],
            align_strand=Strand.NEG, sj_strand=Strand.POS,
        )
    # an ARTIFACTUAL junction: the aligner emitted it from both strands
    artifact = _acc(junctions=[JUNCTION])
    for i in range(20):
        artifact.deposit(
            0, 150, 950, observed_introns=[(201, 900)],
            align_strand=Strand.NEG if i % 2 else Strand.POS, sj_strand=Strand.POS,
        )

    clean_row = clean.tally.sj_count[0]
    artifact_row = artifact.tally.sj_count[0]
    assert list(clean_row) == [0, 20], "a clean junction sits entirely in one aligned-strand column"
    assert list(artifact_row) == [10, 10], "an artifact splits across both"

    # ⛔ THE RULING, as an assertion: collapsing the columns destroys the difference.
    assert int(clean_row.sum()) == int(artifact_row.sum()) == 20, (
        "the two junctions carry the same TOTAL, so a one-column sj_count cannot distinguish them — "
        "which is exactly why the strand split is retained"
    )


def test_a_SENSE_fragment_on_the_minus_strand_is_still_booked_as_MINUS():
    """The discriminating case: sense-to-motif would say column 0, genome strand says column 1."""
    acc = _acc(junctions=[(0, 201, 900, Strand.NEG)])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
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
    out = acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    assert out is DepositOutcome.DEPOSITED, "span 800, L = 101"
    assert int(acc.tally.sj_count[0, 0]) == 1


def test_a_fragment_is_clipped_to_its_reference_and_L_is_the_clipped_length():
    acc = _acc()
    acc.deposit(0, 950, 1200)  # chr1 ends at 1000
    t = acc.tally
    assert int(t.node_contained_count[_node(0, 5), 0]) == 1
    assert int(t.node_contained_inv_opportunity_sum[_node(0, 5)]) == _contained_quantum(0, 5, 50)


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
        acc.deposit(0, s, e, observed_introns=i, sj_strand=Strand.POS) is DepositOutcome.DEPOSITED
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
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS
    )  # annotated junction
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTERGENIC].sum()) == 2
    assert int(p[FragmentPool.DNA_INTRONIC].sum()) == 1
    assert int(p[FragmentPool.DNA_INTRON_EXON].sum()) == 1
    assert int(p[FragmentPool.RNA_SPLICED].sum()) == 1
    assert int(p[FragmentPool.DNA_INTERGENIC_EXON].sum()) == 0


def test_a_pool_is_binned_at_L_and_only_ONCE_per_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, observed_introns=[(300, 340)])  # L = 180 - 40
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTRONIC, 140]) == 1
    assert int(p.sum()) == 1


# ⛔ `test_an_IMPLICIT_splice_is_kept_OUT_of_the_pure_RNA_pool` was DELETED —
#
# It asserted the pool bar that `sj_implicit` existed to apply: a splice inferred rather than observed
# made the fragment's length "a product of the very model the pool is used to fit". The criterion is now
# DETERMINACY, not provenance — a fragment reaches the pool only when exactly ONE hypothesis survived, so
# its `L` is not in doubt at all, however it was arrived at.
#
# ⚠ The old criterion was MEASURED before it was deleted, because the two disagree and the disagreement
# is large and in the damaging direction: on the chr22 pilot the pool reads +0.67 % mean / +2.40 % sd
# against truth under determinacy and −9.58 % / −22.46 % under provenance. Barring inferred lengths
# preferentially bars fragments whose mates sit far apart — a purity filter on a length pool is a length
# filter. The inverse is now asserted by
# `test_ONE_SURVIVING_HYPOTHESIS_DEPOSITS_even_though_its_splice_was_never_sequenced`.


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
        acc.tally.edge_unspliced_inv_length_sum[interior].sum()
        / INV_LENGTH_SCALE
        / (acc.n_edges - 10)
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
        "node_contained_inv_opportunity_sum",
        "node_start_count",
        "edge_unspliced_count",
        "edge_unspliced_inv_length_sum",
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
    assert acc.deposit(0, 50, 500, observed_introns=introns) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert t.qc["introns_absorbed"] == expected_absorbed
    crossings = int(t.edge_unspliced_count.sum())
    assert crossings == expected_crossings
    assert int(t.edge_unspliced_inv_length_sum.sum()) == crossings * inv_length_quantum(
        expected_length - 1
    )


def test_the_path_STARTS_where_its_first_covered_base_is_not_where_the_extent_begins():
    """A leading intron means the molecule does not begin at ``lo``. Attributing it to the node
    containing ``lo`` would credit the start-count invariant — and possibly the contained deposit — to a
    node the fragment never touches."""
    acc = _acc(max_fragment_length=10_000)
    acc.deposit(0, 150, 500, observed_introns=[(150, 480)])  # the path is only [480,500), inside n4
    t = acc.tally
    assert int(t.node_start_count[_node(0, 4)]) == 1, "n4, where the path actually starts"
    assert int(t.node_start_count[_node(0, 1)]) == 0, "not n1, where the extent begins"
    assert int(t.node_contained_count[_node(0, 4), 0]) == 1
    assert int(t.node_contained_inv_opportunity_sum[_node(0, 4)]) == _contained_quantum(0, 4, 20)


def test_a_duplicated_intron_credits_its_junction_ONCE():
    """Two mates reporting the same intron is one splice event, not two."""
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900), (201, 900)], sj_strand=Strand.POS)
    assert int(acc.tally.sj_count[0, 0]) == 1
    assert acc.tally.qc["introns_absorbed"] == 1


def test_ABUTTING_introns_are_MALFORMED_and_merge():
    """⚠ Two introns sharing an endpoint imply a **zero-length exon** between them, which is physically
    impossible — a transcript with one is molecularly identical to a transcript without it. So a single
    molecule can never legitimately use both, and the pair is an alignment artifact.

    The index cannot produce it either: a zero-length exon is dropped when the exon arrays are built,
    which fuses its two flanking introns into one. Merged here, and counted."""
    acc = _acc(junctions=[(0, 201, 400, Strand.POS), (0, 400, 900, Strand.POS)])
    acc.deposit(0, 150, 950, observed_introns=[(201, 400), (400, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.qc["introns_absorbed"] == 1
    assert t.sj_count.sum() == 0, "the merged span 201->900 is not an annotated junction"


def test_a_wide_overlap_no_longer_discards_a_good_fragment():
    """The naive formula gave L = −290 here and filed the fragment as ``dropped_empty`` — invisible to
    the start-count invariant, because a rejected fragment never reaches it."""
    acc = _acc(max_fragment_length=10_000)
    assert (
        acc.deposit(0, 150, 500, observed_introns=[(150, 480), (160, 470)])
        is DepositOutcome.DEPOSITED
    )
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
    cfRNA, **65–69 % of all node_spanning deposits came from spliced fragments**. Indexing those by
    sense-relative-to-motif while the unspliced ones beside them use genome strand would put one array
    into two conventions, and 40–44 % of the spliced deposits would land in the opposite column from
    their unspliced neighbours.

    ⭐⭐ **RE-HOMED, NOT DELETED.** This test used to ride on ``node_spanning``, which was removed —
    and removing it took away the only NODE-axis population a spliced fragment can reach, since a
    spliced fragment can never be *contained* (both endpoints of an annotated intron are cuts). The
    claim is about the CONVENTION, not about that bank, so it now rides on the two banks a spliced
    fragment does reach: ``edge_spliced_count`` beside ``edge_unspliced_count`` at the SAME line, and
    ``sj_count``. ⛔ Deleting it with its old vehicle would have retired the only gate on a rule this
    codebase has already broken once.
    """
    acc = _acc(junctions=[SPAN_JUNCTION_POS])
    acc.deposit(0, 150, 300, align_strand=Strand.NEG)  # unspliced, genome minus
    acc.deposit(  # spliced, genome minus, ANTISENSE to its + junction
        0, 150, 950, observed_introns=[(400, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    t = acc.tally
    # ⭐ ONE line, TWO banks, one column: the unspliced fragment and the spliced one both cross the
    # line at 200, and both are genome-minus. A sense-relative convention would split them.
    assert int(t.edge_unspliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, "genome minus"
    assert int(t.edge_spliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, (
        "the spliced one is ANTISENSE to its + junction, and still books genome minus"
    )
    assert int(t.edge_unspliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.edge_spliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1, "the junction bank too"


def test_a_spliced_SENSE_fragment_books_EDGE_AND_junction_by_GENOME_strand():
    """The discriminating case: sense-to-motif would say column 0 for both; genome strand says 1."""
    acc = _acc(junctions=[SPAN_JUNCTION_NEG])
    acc.deposit(
        0, 150, 950, observed_introns=[(400, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1
    assert int(t.edge_spliced_count[_edge(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1


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
    assert acc.deposit(0, 150, 300, align_strand=undefined) is DepositOutcome.STRAND_UNDEFINED
    t = acc.tally
    assert t.qc["dropped_strand_undefined"] == 1
    assert int(t.node_contained_count.sum()) == 0, "must not be booked into either column"
    assert int(t.node_start_count.sum()) == 0, "a rejected fragment never reaches the invariant"


# ── an AMBIGUOUS PATH deposits nothing, and that is not the same thing as an IMPLICIT splice ────────
#
# Owner ruling, 2026-07-29 (design §9.1). A `SPLICE_IMPLICIT` fragment overlaps an annotated intron and
# matches in every other way, so it DOES deposit — the only thing missing is the sequenced motif, and
# `sj_implicit` records that. But when several candidate transcripts imply DIFFERENT INTRON SETS, the
# implied set fixes `L`, both quanta, the pool bin, the segment list and therefore which lines are
# crossed. There is no partial answer: it cannot deposit spliced (which junction is the unknown), and it
# cannot deposit unspliced either, because `L` involves an intron and does not fit the length
# distribution unless one candidate intron is cut out — the very choice in doubt. Forcing a choice is
# choosing an `L` at random. So it deposits NOTHING and waits for the second pass, which has the
# fragment length AND the strand to discriminate with.


def test_TWO_SURVIVING_HYPOTHESES_deposit_on_NOTHING_and_are_BUFFERED_WHOLE():
    """⛔ The whole point: an undetermined path is not a partial deposit, and not a loss either.

    The fragment used here would otherwise deposit richly — it crosses lines, uses an annotated junction
    and lands in a length pool — so a leak into any one bank is visible. And it must be RETAINED: this is
    the population the second pass drains, so a silent drop would understate what that pass owes.
    """
    acc = _acc(junctions=[JUNCTION])
    outcome = acc.deposit(
        0,
        150,
        950,
        sj_strand=Strand.POS,
        hypotheses=(GapHypothesis(((201, 900),), supporting_t_inds=(7,)), GapHypothesis()),
    )
    assert outcome is DepositOutcome.DEFERRED
    t = acc.tally
    assert t.qc["deferred_undetermined_gap"] == 1
    assert t.qc["deposited"] == 0
    for name in (
        "node_contained_count",
        "node_start_count",
        "edge_unspliced_count",
        "edge_spliced_count",
        "sj_count",
        "pool_lengths",
    ):
        assert int(getattr(t, name).sum()) == 0, f"{name} must be untouched"

    # ⭐ RETAINED WHOLE, and with every hypothesis — the second pass cannot choose between answers it
    # was not given, and it needs the supporting transcripts to weight them by abundance.
    assert len(t.deferred) == 1
    held = t.deferred[0]
    assert (held.ref, held.start, held.end) == (0, 150, 950)
    assert {path.introns for path in held.hypotheses} == {((201, 900),), ()}
    assert held.hypotheses[0].supporting_t_inds == (7,)
    # The genomic path competes against exactly one spliced path: the open question is RNA or gDNA.
    assert t.gap_resolution["gap_deferred_rna_or_gdna"] == 1


def test_ONE_SURVIVING_HYPOTHESIS_DEPOSITS_even_though_its_splice_was_never_sequenced():
    """⭐ The discriminating pair for the test above: **determinacy**, not provenance.

    Both fragments have a splice that was never sequenced. This one has only ONE explanation for its gap,
    so its ``L`` is not in doubt and it deposits — including into the pure-RNA length pool, which the
    deleted ``sj_implicit`` flag used to bar it from. Measured on the chr22 pilot, that bar cost the pool
    **−9.58 % mean / −22.46 % sd** against truth where determinacy reads **+0.67 % / +2.40 %**.
    """
    acc = _acc(junctions=[JUNCTION])
    outcome = acc.deposit(
        0, 150, 950, sj_strand=Strand.POS, hypotheses=(GapHypothesis(((201, 900),)),)
    )
    assert outcome is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.sj_count[0, 0]) == 1
    assert int(t.node_start_count.sum()) == 1
    assert t.qc["deferred_undetermined_gap"] == 0
    assert not t.deferred
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 1
    assert t.gap_resolution["gap_resolved_spliced"] == 1


def test_a_strand_undefined_AMBIGUOUS_fragment_is_counted_as_STRAND_UNDEFINED():
    """⚠ The precedence is part of the contract, because every fragment must count exactly ONCE.

    A fragment can be both. It is filed under the strand, and the reason is which denominator stays
    honest: the deferred queue sizes the population the **second pass can recover**, and a fragment with no
    genome strand is not recoverable — that pass resolves *which path*, not *which strand the read
    aligned to*. Buffering it would promise a recovery that cannot happen, and would put a fragment with
    no column into a queue whose drain needs one.
    """
    acc = _acc()
    outcome = acc.deposit(
        0,
        150,
        300,
        align_strand=Strand.AMBIGUOUS,
        hypotheses=(GapHypothesis(((200, 250),)), GapHypothesis()),
    )
    assert outcome is DepositOutcome.STRAND_UNDEFINED
    t = acc.tally
    assert t.qc["dropped_strand_undefined"] == 1
    assert t.qc["deferred_undetermined_gap"] == 0
    assert not t.deferred, "a fragment with no column must not enter the deferred queue"


def test_a_hypothesis_LONGER_THAN_THE_LIMIT_is_ruled_out_and_the_rest_stands():
    """⭐ `max_fragment_length` is not a new rule — it is the one that already makes TOO_LONG a rejection.

    Short-read chemistry does not sequence molecules past the limit, so a hypothesis implying a longer
    ``L`` is not a molecule this library contains. ⭐ Applied to the GENOMIC hypothesis it is exactly the
    owner's rule "a fragment whose genomic span exceeds the limit must be RNA": the genomic path's ``L``
    **is** that span. Here the span is 800 over a limit of 500, so the genomic path dies and the single
    spliced path — ``L`` = 800 − 699 = 101 — stands alone and deposits.
    """
    acc = _acc(max_fragment_length=500, junctions=[JUNCTION])
    outcome = acc.deposit(
        0,
        150,
        950,
        sj_strand=Strand.POS,
        hypotheses=(GapHypothesis(((201, 900),)), GapHypothesis()),
    )
    assert outcome is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.deposited_lengths.sum()) == 1
    assert int(np.nonzero(t.deposited_lengths)[0][0]) == 101
    assert t.gap_resolution["gap_resolved_spliced"] == 1
    assert not t.deferred


def test_when_the_limit_would_rule_out_EVERY_hypothesis_the_survivors_stand():
    """⚠ The escape clause, and it hands the fragment to the ordinary TOO_LONG rejection.

    Filtering to nothing would mean "this molecule cannot exist", which is not a conclusion the filter is
    entitled to draw — it is a statement about what the chemistry sequences, and the fragment is in the
    BAM. So the survivors stand and ``TOO_LONG`` counts it, exactly as it did before any of this.
    """
    acc = _acc(max_fragment_length=50)
    outcome = acc.deposit(0, 150, 950, hypotheses=(GapHypothesis(),))
    assert outcome is DepositOutcome.TOO_LONG
    assert acc.tally.qc["dropped_too_long"] == 1
    assert not acc.tally.deferred


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
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.NONE)
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
        0,
        150,
        950,
        observed_introns=[(201, 900)],
        align_strand=Strand.POS,
        sj_strand=Strand.AMBIGUOUS,
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
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count.sum()) == 0
    assert t.qc["unannotated_introns"] == 1
    assert t.qc["contradictory_sj_strand"] == 0
    assert int(t.edge_spliced_count.sum()) == 0, "not certified RNA — the strand disagreed"


# ---------------------------------------------------------------------------
# length_sum — the second length tilt
# ---------------------------------------------------------------------------


def test_length_sum_records_L_on_EVERY_object_the_fragment_touches():
    """One fragment, one ``L``, deposited whole on each object it lands on.

    ``length_sum`` is ``Sum L`` over that object's own fragments, so a fragment crossing three lines
    contributes its FULL ``L`` to each of them — the same no-partitioning rule the count and the
    reciprocal sum already follow. A share would re-create the partition dependence the redesign deletes.
    """
    acc = _acc()
    acc.deposit(0, 150, 500)  # L = 350, crosses lines 2, 3, 4 and spans nodes n2, n3
    t = acc.tally
    assert [int(t.edge_unspliced_length_sum[_edge(0, j)]) for j in (2, 3, 4)] == [350, 350, 350]
    assert int(t.node_contained_length_sum.sum()) == 0


def test_length_sum_is_L_and_NOT_the_genomic_span():
    """The molecule's length, so a cut intron does not count — the same ``L`` the pools bin at.

    Binning at the covered/genomic length instead is a defect this project has already paid for once:
    it collapsed the gDNA length pool to a spike at twice the read length
    trap 8), and here it would put a number in ``length_sum`` that no fragment-length model can explain.

    ⚠ **The intron here is UNANNOTATED on purpose.** ``length_sum`` survives only on the two banks the
    deconvolution consumes, and a SPLICED fragment reaches neither — the certified-RNA banks carry no
    length moments, because nothing deconvolves a fragment already known to be RNA. An unannotated
    intron cuts exactly the same bases out of ``L`` while routing to ``edge_unspliced``, so the claim is
    unchanged and observable.
    """
    acc = _acc()  # no junction table: (201, 900) is not annotated
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.POS)
    t = acc.tally
    length = (201 - 150) + (950 - 900)  # 101; the genomic span is 800
    # the only line either block crosses is at position 200 — the block [900, 950) crosses none, and
    # the intron is jumped, never crossed
    assert int(t.edge_unspliced_length_sum[_edge(0, 2)]) == length
    assert int(t.edge_unspliced_length_sum.sum()) == length


def test_BOTH_genome_strands_land_in_the_ONE_length_moment_slot():
    """⭐⭐ **The claim the strand collapse makes, stated directly.** ``count`` keeps two genome-strand
    columns; the length moments keep one, because which strand a read aligned to says nothing about
    whether the molecule was gDNA or RNA.

    ⛔ So a plus fragment and a minus fragment at the SAME object must both reach that single slot. The
    failure mode is a collapse that quietly keeps only one strand's deposits — and until this test
    existed, that defect was caught **only** by the C++ parity gate, which means a matching bug on both
    sides would have passed everything (``TRAPS: perturb-every-gate``).
    """
    acc = _acc()
    acc.deposit(0, 220, 380, align_strand=Strand.POS)
    acc.deposit(0, 220, 380, align_strand=Strand.NEG)
    t = acc.tally
    node = _node(0, 3)
    # the count still separates them, one per column
    assert int(t.node_contained_count[node, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.node_contained_count[node, STRAND_COLUMNS[Strand.NEG]]) == 1
    # ...and the moments pool them into the single slot
    assert int(t.node_contained_inv_opportunity_sum[node]) == 2 * _contained_quantum(0, 3, 160)
    assert int(t.node_contained_length_sum[node]) == 2 * 160

    edge_acc = _acc()
    edge_acc.deposit(0, 120, 320, align_strand=Strand.POS)
    edge_acc.deposit(0, 120, 320, align_strand=Strand.NEG)
    e = edge_acc.tally
    line = _edge(0, 2)
    assert int(e.edge_unspliced_count[line].sum()) == 2
    assert int(e.edge_unspliced_inv_length_sum[line]) == 2 * inv_length_quantum(200 - 1)
    assert int(e.edge_unspliced_length_sum[line]) == 2 * 200


def test_length_sum_and_count_SHARE_a_support():
    """``length_sum`` is zero exactly where ``count`` is zero, on every bank.

    ⚠ ``inv_length_sum`` does NOT have this property — an ``L == 1`` path on a junction books a count
    against a reciprocal sum of 0, the schema's one co-support violation — so a consumer must not assume
    the two behave alike. ``length_sum`` cannot do that: ``L >= 1`` always.
    """
    acc = _acc(junctions=[JUNCTION])
    acc.deposit(0, 150, 500)
    acc.deposit(0, 220, 380)
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    acc.deposit(0, 500, 501)  # L == 1: the co-support edge case
    t = acc.tally
    # ⚠ ``count`` is per genome strand and ``length_sum`` is not, so the support is compared on the
    # count's STRAND-SUMMED total — which is the population the single-column moment describes.
    for count, length_sum in (
        (t.node_contained_count, t.node_contained_length_sum),
        (t.edge_unspliced_count, t.edge_unspliced_length_sum),
    ):
        total = count.sum(axis=1)
        assert np.array_equal(total > 0, length_sum > 0)
        assert int(length_sum[total > 0].min(initial=1)) >= 1


def test_length_sum_over_count_is_the_MEAN_fragment_length_at_that_object():
    """The whole point of the channel: it makes an object's own mean fragment length observable.

    Two components are told apart by their length distributions, and this ratio is the first moment of
    the mixture actually present at this object — no fitted model anywhere in it.
    """
    acc = _acc()
    lengths = (120, 200, 260)
    for length in lengths:  # each crosses the line at 400 and nothing else
        acc.deposit(0, 400 - length // 2, 400 - length // 2 + length)
    t = acc.tally
    edge = _edge(0, 4)
    assert int(t.edge_unspliced_count[edge, 0]) == len(lengths)
    assert int(t.edge_unspliced_length_sum[edge]) == sum(lengths)
    mean = t.edge_unspliced_length_sum[edge] / t.edge_unspliced_count[edge, 0]
    assert mean == pytest.approx(sum(lengths) / len(lengths))


def test_length_sum_SEPARATES_two_populations_that_count_AND_inv_length_sum_CANNOT():
    """⭐ THE reason this channel exists, as an exact arithmetic identity rather than a simulation.

    Two fragment sets are constructed with the SAME count and a BYTE-IDENTICAL ``inv_length_sum``
    (``1/2 + 1/3 + 1/6 == 1/2 + 1/4 + 1/4``, and the fixed-point quanta happen to sum exactly), yet
    different total length. The shipped ``(count, inv_length_sum)`` pair is provably blind to the
    difference; ``length_sum`` sees it.

    This is the per-object shadow of the structural result: at an edge the count row is
    ``(mu_g − 1, mu_r − 1)`` and the reciprocal row is ``(1, 1)``, so the determinant is ``mu_g − mu_r``
    and the pair carries ZERO information about the gDNA/RNA split when the two means agree.
    """
    edge = _edge(0, 4)  # the line at 400, flanked by two wide nodes

    def deposit_lengths(lengths):
        acc = _acc()
        for length in lengths:
            acc.deposit(0, 399, 399 + length)  # every one crosses the line at 400, and only it
        t = acc.tally
        return (
            int(t.edge_unspliced_count[edge, 0]),
            int(t.edge_unspliced_inv_length_sum[edge]),
            int(t.edge_unspliced_length_sum[edge]),
        )

    count_a, inv_a, len_a = deposit_lengths((3, 4, 7))  # placements 2, 3, 6
    count_b, inv_b, len_b = deposit_lengths((3, 5, 5))  # placements 2, 4, 4

    assert count_a == count_b == 3, "same number of fragments"
    assert inv_a == inv_b == INV_LENGTH_SCALE, "the OLD pair cannot tell these apart, exactly"
    assert (len_a, len_b) == (14, 13), "the new channel can"


def test_the_density_FIELD_NAME_is_gone_everywhere():
    """The rename is complete, so no consumer can reach a half-migrated schema.

    ``inv_length_sum`` is an exact density at an edge and is NOT one at a node; keeping the old name
    would put one word on two concepts, which is.
    """
    t = _acc().tally
    stale = [name for name in t.__slots__ if name.endswith("_density")]
    assert stale == []
    assert {"node_contained_inv_opportunity_sum", "node_contained_length_sum"} <= set(t.__slots__)
