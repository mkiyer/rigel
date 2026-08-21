"""THE ACCUMULATOR SPEC — the matrix the reference and the native build are both gated on.

       §10.4

``_accumulator_reference.py`` is the executable specification; the native accumulator is required to
reproduce it byte for byte. This module is what "correct" means for both.

THE RULE UNDER TEST, in five boundaries. The genome is a graph: REGIONS are half-open intervals tiling each
reference, and the 0-bp BOUNDARIES between adjacent regions are CONTIGUOUS boundaries. A SJ boundary is a directed
donor→acceptor link from the annotation. A fragment is a PATH — its aligned blocks joined across mate
gaps and broken by introns — of length ``L = span − Σ intron``. Regions count fragments CONTAINED; boundaries
count fragments CROSSING; and every path books its FIRST/LAST covered base and the regions it STRICTLY
spans (the start/end/span banks, 2026-08-21). Each population stores only the channels something READS —
an integer count, a float64 reciprocal-opportunity sum ``Σ 1/A(w)`` (no fixed point anywhere, and no
``Σ L`` — deleted 2026-08-13), and on the contiguous boundaries the CONSERVED MASS, which sums to one
per fragment.

⚠ **No partitioning.** Every crossed boundary receives the FULL weight. The chance that a length-``L``
fragment crosses a given boundary is proportional to ``L − 1`` and the deposit is ``1/(L − 1)``, so the two
cancel and every fragment length contributes equally to each boundary. Dividing by the number of boundaries
crossed destroys that cancellation and makes the answer depend on region spacing — measured up to **3.6× low**.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.types import Strand

from ._accumulator_reference import (
    STRAND_COLUMNS,
    Accumulator,
    DepositOutcome,
    FragmentPool,
    GapHypothesis,
    Partition,
)

#: ⭐ ONE CONVENTION: fractions are float64. A sum of `n` round-to-nearest additions differs from the
#: real-arithmetic answer by at most `n` ulp. ⛔ DERIVED from the machine, never fitted.
EPS = float(np.finfo(np.float64).eps)


def close(got: float, want: float, deposits: int) -> bool:
    """`got == want` to within the representation error of `deposits` additions."""
    return abs(got - want) <= max(abs(want), 1.0) * deposits * EPS



# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

#: chr1 region_bounds   0    100   200   201   400   900   1000
#: regions        n0    n1    n2*   n3    n4    n5      (* n2 is 1 bp: [200,201))
#: boundaries           1     2     3     4     5          (local region_bound index)
CHR1_REGION_BOUNDS = [0, 100, 200, 201, 400, 900, 1000]
CHR2_REGION_BOUNDS = [0, 500, 1000]

#: coarse region type: 0 intergenic, 1 intron, 2 exon
CHR1_TYPES = [0, 2, 2, 1, 2, 0]
CHR2_TYPES = [0, 2]

#: an annotated intron whose endpoints are region_bounds 3 and 5, so it SWALLOWS the boundary at region_bound 4
SJ = (0, 201, 900, Strand.POS)


def _partition(sj=()):
    return Partition.from_region_bounds(
        [CHR1_REGION_BOUNDS, CHR2_REGION_BOUNDS], region_types=[CHR1_TYPES, CHR2_TYPES], sj=sj
    )


def _acc(sj=(), **kw):
    return Accumulator(_partition(sj), **kw)


def _boundary(ref, boundary):
    """Global contiguous-boundary id of the boundary at local region_bound index ``boundary``."""
    return (0 if ref == 0 else len(CHR1_REGION_BOUNDS) - 2) + boundary - 1



def _contained_quantum(ref, local, length):
    """The deposit a CONTAINED length-``length`` fragment makes on region ``(ref, local)``.

    ⭐⭐ **The deposit is ``1/OPPORTUNITY``, not ``1/length``** — a length-`w` fragment inside a region of
    length `ell` had `ell − w + 1` admissible start positions, so `1/(ell − w + 1)` cancels the
    opportunity identically and the channel is a DENSITY for any length distribution
    (`test_fragment_length_proof.test_the_region_deposit_is_the_RECIPROCAL_OPPORTUNITY_...`). ⚠ Derived from
    the fixture's own region_bounds rather than written as a number, so an assertion states the RULE.
    """
    region_bounds = CHR1_REGION_BOUNDS if ref == 0 else CHR2_REGION_BOUNDS
    return 1.0 / (region_bounds[local + 1] - region_bounds[local] - length + 1)


def _region(ref, local):
    return (0 if ref == 0 else len(CHR1_REGION_BOUNDS) - 1) + local


# ---------------------------------------------------------------------------
# the reciprocal-opportunity density
# ---------------------------------------------------------------------------
#
# ⚠ `test_inv_length_quantum_is_exact_and_rounds_half_away_from_zero` lived here and is GONE with the
# fixed point (owner, 2026-08-10: one numeric convention). It asserted a rounding mode that no longer
# exists. ⭐ Nothing is lost: it pinned the REPRESENTATION, while the tests below pin the THEOREM, and
# float64 satisfies the theorem 1e5-7e5x more tightly than the grid it replaced.


def test_one_fragment_recovers_its_own_reciprocal_length():
    acc = _acc()
    acc.deposit(0, 120, 320)
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 1
    assert close(float(t.boundary_unspliced_inv_length_sum[_boundary(0, 2)]), 1.0 / (200 - 1), 1)


# ---------------------------------------------------------------------------
# contained / crossing
# ---------------------------------------------------------------------------


def test_a_contained_fragment_touches_ONE_region_and_no_boundary():
    acc = _acc()
    assert acc.deposit(0, 220, 380) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.region_contained_count[_region(0, 3), 0]) == 1
    assert close(float(t.region_contained_inv_opportunity_sum[_region(0, 3)]), _contained_quantum(0, 3, 160), 1)
    assert t.region_contained_count.sum() == 1
    assert t.boundary_unspliced_count.sum() == 0


def test_a_fragment_ENDING_AT_a_boundary_does_not_cross_it():
    """A boundary is 0 bp wide: crossing needs a base on each side."""
    acc = _acc()
    acc.deposit(0, 120, 200)
    assert acc.tally.boundary_unspliced_count.sum() == 0
    assert int(acc.tally.region_contained_count[_region(0, 1), 0]) == 1


def test_a_fragment_STARTING_AT_a_boundary_does_not_cross_it():
    acc = _acc()
    acc.deposit(0, 201, 390)
    assert acc.tally.boundary_unspliced_count.sum() == 0
    assert int(acc.tally.region_contained_count[_region(0, 3), 0]) == 1


def test_a_fragment_crossing_four_regions_credits_exactly_THREE_boundaries_at_FULL_weight():
    acc = _acc()
    acc.deposit(0, 150, 500)  # touches n1 n2 n3 n4 -> boundaries at 200, 201, 400
    t = acc.tally
    assert [int(t.boundary_unspliced_count[_boundary(0, j), 0]) for j in (1, 2, 3, 4, 5)] == [0, 1, 1, 1, 0]
    quantum = 1.0 / (350 - 1)
    assert all(close(float(t.boundary_unspliced_inv_length_sum[_boundary(0, j)]), quantum, 1) for j in (2, 3, 4))
    assert t.region_contained_count.sum() == 0


def test_a_fragment_covering_a_1bp_region_credits_BOTH_boundaries_and_conserves_its_mass():
    """1 bp regions are legal — 15,687 of them at human scale — and nothing may assume length > 1.

    ⭐ The 1 bp region is the sharpest case for the CONSERVED MASS: its slice is a single base shared
    between two bounding boundaries, so the fragment's three slices are 50 / 1 / 99 bases of a 150 bp
    molecule and must still sum to exactly one fragment.
    ⚠ This test used to assert a ``region_spanning`` deposit here. That bank was removed on evidence, and
    the mass is what now makes the case observable — a strictly stronger statement, since a count of 1
    survives any error in *how much* of the fragment was attributed."""
    acc = _acc()
    acc.deposit(0, 150, 300)
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 1
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 1
    assert int(t.region_contained_count[_region(0, 2), 0]) == 0
    mass = int(t.boundary_unspliced_mass.sum())
    assert abs(mass - 1.0) <= 8 * EPS, (
        f"the fragment deposited {mass} fragments, not 1"
    )


def test_a_fragment_LONGER_than_a_region_is_not_contained_and_crosses_exactly_one_boundary():
    """``contained`` needs ``L <= region``. One base longer and the fragment is a CROSSING instead, and it
    touches the region axis nowhere at all — that axis has exactly one population now."""
    acc = _acc()
    acc.deposit(0, 200, 202)  # n2 is 1 bp, L = 2
    t = acc.tally
    assert t.region_contained_count.sum() == 0
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 1


# ---------------------------------------------------------------------------
# splicing
# ---------------------------------------------------------------------------


def test_a_spliced_jump_deposits_NOTHING_on_the_boundaries_it_splices_over():
    """⭐ The defect this design removes. The intron [201,900) swallows the boundary at 400; the old rule,
    which asked only "does another slice follow?", could not tell that from a contiguous crossing."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    length = (950 - 150) - (900 - 201)
    assert int(t.sj_count[0, 0]) == 1
    assert close(float(t.sj_inv_length_sum[0]), 1.0 / (length - 1), 1)
    assert int(t.boundary_unspliced_count[_boundary(0, 4), 0]) == 0, "the swallowed boundary at 400"
    assert int(t.boundary_spliced_count[_boundary(0, 4), 0]) == 0


def test_a_spliced_fragments_own_BLOCK_crossings_go_in_the_SPLICED_bank():
    """A spliced fragment is certified RNA — gDNA cannot be spliced — so a boundary its block genuinely
    crosses is the cleanest RNA marker available at a boundary."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.boundary_spliced_count[_boundary(0, 2), 0]) == 1, "block [150,201) crosses the boundary at 200"
    assert t.boundary_unspliced_count.sum() == 0


def test_a_MULTI_SEGMENT_unspliced_fragment_conserves_its_mass_across_BOTH_segments():
    """An unannotated intron sits strictly inside n4 = [400,900), so segment 1 crosses the boundary at 400
    and segment 2 crosses the boundary at 900 — two segments, each touching a boundary.

    ⭐ Every segment touches a boundary, so the conservation law applies unchanged and the mass sums to
    exactly one fragment. ⚠ This is the case the law's "deposited, unspliced, annotated" wording exists
    for: had either segment touched NO boundary, its bases would have had nowhere conserved to go and the
    total would be a PARTIAL."""
    acc = _acc()
    acc.deposit(0, 380, 950, observed_introns=[(500, 600)])
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 4), 0]) == 1, "boundary 400, from segment 1"
    assert int(t.boundary_unspliced_count[_boundary(0, 5), 0]) == 1, "boundary 900, from segment 2"
    mass = int(t.boundary_unspliced_mass.sum())
    assert abs(mass - 1.0) <= 8 * EPS, (
        f"a two-segment fragment deposited {mass} fragments, not 1"
    )


def test_an_UNANNOTATED_intron_credits_no_sj_and_nothing_across_the_gap():
    """Owner ruling: unannotated sj are disproportionately artifactual, so they deposit on the
    UNSPLICED channel and compete with gDNA rather than being certified RNA."""
    acc = _acc(sj=[SJ])
    # [200,400) is NOT annotated; it swallows the boundary at 201
    acc.deposit(0, 50, 500, observed_introns=[(200, 400)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.sj_count.sum() == 0
    assert t.boundary_spliced_count.sum() == 0, "not certified RNA — it competes with gDNA"
    assert int(t.boundary_unspliced_count[_boundary(0, 1), 0]) == 1, (
        "block [50,200) crosses the boundary at 100"
    )
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 0, "the swallowed boundary at 201"
    assert t.qc["unannotated_introns"] == 1


def test_a_fragment_straddling_two_regions_without_crossing_a_boundary_is_NOT_contained():
    """⚠ An unannotated intron can swallow every boundary between two blocks. The fragment then crosses
    nothing, yet it straddles two regions — crediting it as *contained* would put its whole length in a
    region it only partly overlaps. It deposits on no object, and the start count is what keeps that
    visible rather than silent."""
    acc = _acc()
    acc.deposit(
        0, 120, 500, observed_introns=[(200, 400)]
    )  # blocks land in n1 and n4, crossing no boundary
    t = acc.tally
    assert t.boundary_unspliced_count.sum() == 0
    assert t.region_contained_count.sum() == 0
    assert int(t.region_start_count.sum()) == 1, "still counted"


def test_an_unannotated_intron_inside_one_region_is_a_contained_unspliced_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, observed_introns=[(300, 340)])
    t = acc.tally
    assert int(t.region_contained_count[_region(0, 3), 0]) == 1
    assert close(float(t.region_contained_inv_opportunity_sum[_region(0, 3)]), _contained_quantum(0, 3, 180 - 40), 1)
    assert t.qc["unannotated_introns"] == 1


def test_opposite_strand_sj_at_the_same_coordinates_are_DISTINCT_boundaries():
    """Biologically impossible — splice motifs are not palindromic — so only a synthetic stress test can
    reach it, which is exactly why one exists."""
    acc = _acc(sj=[(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert t.sj_count.sum() == 1
    assert int(t.sj_count[1, STRAND_COLUMNS[Strand.NEG]]) == 1, "the NEG boundary (id 1), genome minus"


@pytest.mark.parametrize("order", [("POS_first", 1), ("NEG_first", -1)])
def test_a_sj_id_is_a_function_of_the_PARTITION_not_of_argument_order(order):
    """⛔ The sj-boundary id IS the rank in the sort, so the sort must be total — otherwise the id
    depends on the order the caller happened to list the sj in, and the same graph gets two
    labellings.

    ⚠ This is the ONLY test that pins ``strand`` as part of the sort key. ``np.lexsort`` is stable, so a
    key of ``(acceptor, donor)`` alone gives the right answer for any input whose ties already arrive in
    the right order — which is every other test, and both real indexes. Reversing the argument order is
    what makes the missing key observable, and a strand-coincident pair is the only tie there is.
    """
    _, direction = order
    sj = [(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)][::direction]
    part = Partition.from_region_bounds([CHR1_REGION_BOUNDS], region_types=[CHR1_TYPES], sj=sj)
    assert [int(s) for s in part.sj_strand] == [int(Strand.POS), int(Strand.NEG)], (
        "POS must sort to slot 0 whichever order it was passed in"
    )


def test_a_fragment_using_TWO_sj_credits_BOTH():
    """Owner ruling: each boundary owns its own expectation, and the strand model is fitted from a separate
    scan output, so crediting every sj distorts nothing.

    ⚠ The two introns must be separated by a real exon. Abutting introns imply a zero-length exon and are
    malformed (see ``test_ABUTTING_introns_are_MALFORMED_and_merge``), so this needs its own partition
    with room for an exon between them."""
    part = Partition.from_region_bounds(
        [[0, 100, 200, 300, 400, 500, 600]],
        region_types=[[0, 2, 1, 2, 1, 2]],
        sj=[(0, 100, 200, Strand.POS), (0, 300, 400, Strand.POS)],
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
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 1]) == 1
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 0


def test_EVERY_bank_including_the_sj_is_indexed_by_GENOME_strand():
    """⭐ One convention throughout. Sense/antisense is DERIVED, never stored: the sj boundary carries
    its own genomic strand, so a consumer computes ``sense = (fragment strand == sj strand)``.

    Here a genome-minus fragment splices across a ``+`` sj. Under a sense convention it would land
    in the antisense column; under the genome convention it lands in the minus column. Those happen to be
    the same index, so the discriminating case is the next test."""
    acc = _acc(sj=[SJ])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    assert int(acc.tally.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1


def test_the_sj_STRAND_SPLIT_IS_RETAINED_FOR_ALIGNER_ARTIFACT_DETECTION():
    """⛔⛔ **DO NOT COLLAPSE ``sj_count`` TO ONE COLUMN.** A sj is stranded by its genomic splicing
    MOTIF, so the strand of the *fragments* on it looks redundant — and every consumer today sums the two
    columns. It is retained anyway, on an owner ruling (2026-08-08), and this test is why.

    ⭐ **THE MECHANISM.** Aligners emit false-positive ``N`` CIGAR ops from plain genomic DNA.
    ``rigel.splice_blacklist`` catches the ones the sister tool ``alignable`` has already enumerated by
    coordinate — an a-priori list, and far from complete. A second, EMPIRICAL detector exists in a
    stranded library: real sj inherit the library's global strand specificity, while alignment
    artifacts deposit onto BOTH strands and deviate from it. That test is only possible if the
    per-sj split by ALIGNED strand survives into the payload.

    ⚠ Unstranded data cannot use it — with κ = ½ there is no expectation to deviate from. That is a
    property of the detector, not a reason to drop the column.

    ⭐⭐ **AND THE DISCRIMINATING INFORMATION LIVES ONLY IN THE SPLIT**, which is the whole ruling: the
    clean sj and the artifactual one below carry the SAME total, so a collapsed bank cannot tell
    them apart at all. A tidying pass that removed the column — the same "store a channel where a named
    consumer reads it" principle that correctly removed six other banks — would delete this detector
    before it was built.
    """
    # a CLEAN sj: every fragment on one aligned strand, as a stranded library produces
    clean = _acc(sj=[SJ])
    for _ in range(20):
        clean.deposit(
            0, 150, 950, observed_introns=[(201, 900)],
            align_strand=Strand.NEG, sj_strand=Strand.POS,
        )
    # an ARTIFACTUAL sj: the aligner emitted it from both strands
    artifact = _acc(sj=[SJ])
    for i in range(20):
        artifact.deposit(
            0, 150, 950, observed_introns=[(201, 900)],
            align_strand=Strand.NEG if i % 2 else Strand.POS, sj_strand=Strand.POS,
        )

    clean_row = clean.tally.sj_count[0]
    artifact_row = artifact.tally.sj_count[0]
    assert list(clean_row) == [0, 20], "a clean sj sits entirely in one aligned-strand column"
    assert list(artifact_row) == [10, 10], "an artifact splits across both"

    # ⛔ THE RULING, as an assertion: collapsing the columns destroys the difference.
    assert int(clean_row.sum()) == int(artifact_row.sum()) == 20, (
        "the two sj carry the same TOTAL, so a one-column sj_count cannot distinguish them — "
        "which is exactly why the strand split is retained"
    )


def test_the_sj_MASS_KEEPS_THE_SAME_SPLIT_and_it_is_the_ONLY_mass_that_does():
    """⭐⭐⭐ **THE MASS TWIN OF THE TEST ABOVE, AND THE REVERSED RULING MADE FALSIFIABLE.**
    ``accumulator.h`` ruled a mass bank is ONE value because *"nothing reads a mass per strand"*. The
    premise changed and the ruling was reversed on this axis alone (owner, 2026-08-12): the empirical
    artifact detector above needs the split, and a COUNT cannot separate a sj used by many short
    fragments from one used by few long ones — only a mass can. It is also what makes the filter
    single-pass instead of tally-filter-re-accumulate.

    ⛔ The gate is that the mass lands in the SAME column the count does. If the two used different
    conventions, ``mass[c]/count[c]`` would be a ratio across two different populations and the detector
    would be reading noise.
    ⚠ The ruling still STANDS for ``boundary_unspliced_mass`` and ``boundary_spliced_mass``, which have no such
    consumer — ``one-thing-varied``, checked by the shape gate in ``test_accumulator_payload``.
    """
    pos, neg = STRAND_COLUMNS[Strand.POS], STRAND_COLUMNS[Strand.NEG]

    def one(align):
        acc = _acc(sj=[SJ])
        acc.deposit(0, 150, 950, observed_introns=[(201, 900)], align_strand=align,
                    sj_strand=Strand.POS)
        return acc.tally

    t_pos, t_neg = one(Strand.POS), one(Strand.NEG)
    assert t_pos.sj_mass.ndim == 2, "the sj mass carries a strand"

    # ⭐ each deposit reaches its OWN column and only its own — the split is real, not decorative
    assert t_pos.sj_mass[0, pos] > 0.0 and t_pos.sj_mass[0, neg] == 0.0
    assert t_neg.sj_mass[0, neg] > 0.0 and t_neg.sj_mass[0, pos] == 0.0
    # ⛔ ...and it is the column the COUNT used, not merely "some column"
    assert int(t_pos.sj_count[0, pos]) == 1 and int(t_neg.sj_count[0, neg]) == 1
    # ⭐ the two agree, because the DEPOSIT RULE has no strand in it: the split is per strand, the rule
    # is not. A column-dependent share would be a different defect and this is what would catch it.
    assert close(float(t_pos.sj_mass[0, pos]), float(t_neg.sj_mass[0, neg]), 12)
    # ⭐⭐ and summing the columns returns what the one-column bank held — the property that left every
    # consumer below `substrate` unchanged by the schema move.
    assert close(float(t_pos.sj_mass.sum()), float(t_pos.sj_mass[0, pos]), 12)


def test_a_SENSE_fragment_on_the_minus_strand_is_still_booked_as_MINUS():
    """The discriminating case: sense-to-motif would say column 0, genome strand says column 1."""
    acc = _acc(sj=[(0, 201, 900, Strand.NEG)])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1


# ---------------------------------------------------------------------------
# bounded influence, clipping, references
# ---------------------------------------------------------------------------


def test_a_fragment_over_the_length_limit_deposits_NOTHING_and_is_COUNTED():
    """⭐ Bounded influence. Unbounded, 1,000 read groups own 99.8 % of all boundary crossings on a real
    library; with the limit on ``L`` they own 4.16 %. A silent drop would hide that."""
    acc = _acc(max_fragment_length=200)
    assert acc.deposit(0, 100, 500) is DepositOutcome.TOO_LONG
    t = acc.tally
    assert t.boundary_unspliced_count.sum() == 0
    assert t.region_start_count.sum() == 0
    assert t.qc["dropped_too_long"] == 1


def test_the_limit_applies_to_L_and_NOT_to_the_SPAN():
    """⚠ A 300 bp molecule across a 10 kb intron has a 10 kb span. Limiting the span discards every
    spliced fragment — 37.96 % of read groups measured, against 5.45 % when the limit is on ``L``."""
    acc = _acc(sj=[SJ], max_fragment_length=200)
    out = acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    assert out is DepositOutcome.DEPOSITED, "span 800, L = 101"
    assert int(acc.tally.sj_count[0, 0]) == 1


def test_a_fragment_is_clipped_to_its_reference_and_L_is_the_clipped_length():
    acc = _acc()
    acc.deposit(0, 950, 1200)  # chr1 ends at 1000
    t = acc.tally
    assert int(t.region_contained_count[_region(0, 5), 0]) == 1
    assert close(float(t.region_contained_inv_opportunity_sum[_region(0, 5)]), _contained_quantum(0, 5, 50), 1)


def test_a_single_region_reference_has_no_boundaries_and_still_accepts_a_fragment():
    acc = Accumulator(Partition.from_region_bounds([[0, 1000]], region_types=[[0]]))
    assert acc.n_boundaries == 0
    acc.deposit(0, 100, 300)
    assert int(acc.tally.region_contained_count[0, 0]) == 1


def test_the_per_reference_offsets_do_not_bleed():
    """chr1's fragment crosses its boundaries 2 and 3; chr2's crosses chr2's boundary 1. Nothing lands on a
    reference it did not come from — the failure mode that once dropped 476,719 of 476,732 fragments."""
    acc = _acc()
    acc.deposit(0, 150, 300)  # crosses the boundaries at 200 AND 201
    acc.deposit(1, 400, 700)  # crosses chr2's boundary at 500
    t = acc.tally
    assert [int(t.boundary_unspliced_count[e, 0]) for e in range(acc.n_boundaries)] == [0, 1, 1, 0, 0, 1]
    assert int(t.boundary_unspliced_count.sum()) == 3
    assert int(t.region_start_count[_region(0, 1)].sum()) == 1
    assert int(t.region_start_count[_region(1, 0)].sum()) == 1


# ---------------------------------------------------------------------------
# the invariant that can actually fire
# ---------------------------------------------------------------------------


def test_every_accepted_fragment_increments_exactly_ONE_start_count():
    """⚠ The crossing and contained totals are tautologies — they can only be evaluated by re-running
    the deposit. This one is checkable against a number the scanner knows independently."""
    acc = _acc(sj=[SJ])
    fragments = [(120, 320, ()), (220, 380, ()), (150, 950, [(201, 900)]), (950, 1200, ())]
    accepted = sum(
        acc.deposit(0, s, e, observed_introns=i, sj_strand=Strand.POS) is DepositOutcome.DEPOSITED
        for s, e, i in fragments
    )
    assert accepted == 4
    assert int(acc.tally.region_start_count.sum()) == 4
    assert acc.tally.qc["deposited"] == 4


# ---------------------------------------------------------------------------
# the five fragment-length pools
# ---------------------------------------------------------------------------


def test_each_pool_is_reached_only_by_its_own_structural_class():
    acc = _acc(sj=[SJ])
    acc.deposit(0, 10, 90)  # contained in n0 — intergenic
    acc.deposit(0, 210, 390)  # contained in n3 — intronic
    acc.deposit(0, 380, 420)  # crosses the boundary at 400 only — flanks intron|exon
    acc.deposit(0, 950, 990)  # contained in n5 — intergenic
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS
    )  # annotated sj
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


def test_a_multi_boundary_crossing_enters_NO_pool():
    """A splash read straddles ONE probe boundary. A fragment crossing several boundaries has no single
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


def _uniform_accumulator(region_bp, ref_len):
    region_bounds = list(range(0, ref_len + 1, region_bp))
    return Accumulator(
        Partition.from_region_bounds([region_bounds], region_types=[[0] * (len(region_bounds) - 1)]), max_fragment_length=10_000
    )


@pytest.mark.parametrize("region_bp", [50, 200, 1000])
def test_the_crossing_DENSITY_recovers_the_true_density_with_NO_length_model(region_bp):
    """⭐ ``E[Σ 1/(L−1)] = ρ`` exactly, for ANY fragment-length distribution. This is the identity the
    whole design rests on and the reason no divisor and no length model appear at an boundary.

    It must hold at every region spacing. Partitioning the weight by the number of boundaries crossed breaks
    it — measured 0.28× at 50 bp regions, 0.54× at 100 bp, 0.91× at 200 bp — so this test is also what
    forbids partitioning.
    """
    ref_len, rho = 200_000, 0.05
    acc = _uniform_accumulator(region_bp, ref_len)
    rng = np.random.default_rng(7)
    starts, ends, _ = _corpus(rng, int(rho * ref_len), ref_len)
    for s, e in zip(starts, ends):
        acc.deposit(0, int(s), int(e))
    interior = slice(5, acc.n_boundaries - 5)
    estimate = (
        acc.tally.boundary_unspliced_inv_length_sum[interior].sum() / (acc.n_boundaries - 10)
    )
    assert 0.98 <= estimate / rho <= 1.02, f"{estimate / rho:.4f} at {region_bp} bp regions"


def test_the_crossing_COUNT_recovers_density_times_mean_length():
    """The companion identity ``E[count] = ρ·(E[L] − 1)``. Together with the boundary above, this is the 2×2
    that separates gDNA from RNA by fragment length alone."""
    ref_len, rho = 200_000, 0.05
    acc = _uniform_accumulator(200, ref_len)
    rng = np.random.default_rng(11)
    starts, ends, lengths = _corpus(rng, int(rho * ref_len), ref_len)
    for s, e in zip(starts, ends):
        acc.deposit(0, int(s), int(e))
    interior = slice(5, acc.n_boundaries - 5)
    per_boundary = acc.tally.boundary_unspliced_count[interior, :].sum() / (acc.n_boundaries - 10)
    expected = rho * (lengths.mean() - 1.0)
    assert 0.98 <= per_boundary / expected <= 1.02, f"{per_boundary / expected:.4f}"


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
        "region_contained_count",
        "region_contained_inv_opportunity_sum",
        "region_start_count",
        "boundary_unspliced_count",
        "boundary_unspliced_inv_length_sum",
        "pool_lengths",
    ):
        got, want = getattr(a, field), getattr(b, field)
        # ⛔ Exact for the INTEGER banks; for the float64 fractions the deposits are re-associated by
        # the shuffle, and float addition is not associative. The tolerance is the representation's,
        # derived from the deposit count, never fitted.
        if getattr(got, "dtype", None) == np.float64:
            assert np.allclose(got, want, rtol=len(starts) * EPS, atol=0.0), field
        else:
            assert np.array_equal(got, want), field


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
    crossings = int(t.boundary_unspliced_count.sum())
    assert crossings == expected_crossings
    assert close(
        float(t.boundary_unspliced_inv_length_sum.sum()), crossings / (expected_length - 1), crossings
    )


def test_the_path_STARTS_where_its_first_covered_base_is_not_where_the_extent_begins():
    """A leading intron means the molecule does not begin at ``lo``. Attributing it to the region
    containing ``lo`` would credit the start-count invariant — and possibly the contained deposit — to a
    region the fragment never touches."""
    acc = _acc(max_fragment_length=10_000)
    acc.deposit(0, 150, 500, observed_introns=[(150, 480)])  # the path is only [480,500), inside n4
    t = acc.tally
    assert int(t.region_start_count[_region(0, 4)].sum()) == 1, "n4, where the path actually starts"
    assert int(t.region_start_count[_region(0, 1)].sum()) == 0, "not n1, where the extent begins"
    assert int(t.region_contained_count[_region(0, 4), 0]) == 1
    assert close(float(t.region_contained_inv_opportunity_sum[_region(0, 4)]), _contained_quantum(0, 4, 20), 1)


def test_a_duplicated_intron_credits_its_sj_ONCE():
    """Two mates reporting the same intron is one splice event, not two."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900), (201, 900)], sj_strand=Strand.POS)
    assert int(acc.tally.sj_count[0, 0]) == 1
    assert acc.tally.qc["introns_absorbed"] == 1


def test_ABUTTING_introns_are_MALFORMED_and_merge():
    """⚠ Two introns sharing an endpoint imply a **zero-length exon** between them, which is physically
    impossible — a transcript with one is molecularly identical to a transcript without it. So a single
    molecule can never legitimately use both, and the pair is an alignment artifact.

    The index cannot produce it either: a zero-length exon is dropped when the exon arrays are built,
    which fuses its two flanking introns into one. Merged here, and counted."""
    acc = _acc(sj=[(0, 201, 400, Strand.POS), (0, 400, 900, Strand.POS)])
    acc.deposit(0, 150, 950, observed_introns=[(201, 400), (400, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.qc["introns_absorbed"] == 1
    assert t.sj_count.sum() == 0, "the merged span 201->900 is not an annotated sj"


def test_a_wide_overlap_no_longer_discards_a_good_fragment():
    """The naive formula gave L = −290 here and filed the fragment as ``dropped_empty`` — invisible to
    the start-count invariant, because a rejected fragment never reaches it."""
    acc = _acc(max_fragment_length=10_000)
    assert (
        acc.deposit(0, 150, 500, observed_introns=[(150, 480), (160, 470)])
        is DepositOutcome.DEPOSITED
    )
    assert acc.tally.qc["dropped_empty"] == 0
    assert int(acc.tally.region_start_count.sum()) == 1


# ---------------------------------------------------------------------------
# the region banks carry ONE strand convention
# ---------------------------------------------------------------------------

#: a sj far enough right that the first block still spans region n2 = [200,201)
SPAN_SJ_POS = (0, 400, 900, Strand.POS)
SPAN_SJ_NEG = (0, 400, 900, Strand.NEG)


def test_a_spliced_and_an_unspliced_fragment_of_the_SAME_genome_strand_share_a_column():
    """⚠ One array, one convention.

    A spliced fragment cannot be *contained* — both endpoints of an annotated intron are region_bounds, so it
    always crosses its sj boundary — but its blocks routinely SPAN a region whole. Measured on real
    cfRNA, **65–69 % of all region_spanning deposits came from spliced fragments**. Indexing those by
    sense-relative-to-motif while the unspliced ones beside them use genome strand would put one array
    into two conventions, and 40–44 % of the spliced deposits would land in the opposite column from
    their unspliced neighbours.

    ⭐⭐ **RE-HOMED, NOT DELETED.** This test used to ride on ``region_spanning``, which was removed —
    and removing it took away the only REGION-axis population a spliced fragment can reach, since a
    spliced fragment can never be *contained* (both endpoints of an annotated intron are region_bounds). The
    claim is about the CONVENTION, not about that bank, so it now rides on the two banks a spliced
    fragment does reach: ``boundary_spliced_count`` beside ``boundary_unspliced_count`` at the SAME boundary, and
    ``sj_count``. ⛔ Deleting it with its old vehicle would have retired the only gate on a rule this
    codebase has already broken once.
    """
    acc = _acc(sj=[SPAN_SJ_POS])
    acc.deposit(0, 150, 300, align_strand=Strand.NEG)  # unspliced, genome minus
    acc.deposit(  # spliced, genome minus, ANTISENSE to its + sj
        0, 150, 950, observed_introns=[(400, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    t = acc.tally
    # ⭐ ONE boundary, TWO banks, one column: the unspliced fragment and the spliced one both cross the
    # boundary at 200, and both are genome-minus. A sense-relative convention would split them.
    assert int(t.boundary_unspliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, "genome minus"
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, (
        "the spliced one is ANTISENSE to its + sj, and still books genome minus"
    )
    assert int(t.boundary_unspliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1, "the sj bank too"


def test_a_spliced_SENSE_fragment_books_BOUNDARY_AND_sj_by_GENOME_strand():
    """The discriminating case: sense-to-motif would say column 0 for both; genome strand says 1."""
    acc = _acc(sj=[SPAN_SJ_NEG])
    acc.deposit(
        0, 150, 950, observed_introns=[(400, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1


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
    assert int(t.region_contained_count.sum()) == 0, "must not be booked into either column"
    assert int(t.region_start_count.sum()) == 0, "a rejected fragment never reaches the invariant"


# ── an AMBIGUOUS PATH deposits nothing, and that is not the same thing as an IMPLICIT splice ────────
#
# Owner ruling, 2026-07-29 (design §9.1). A `SPLICE_IMPLICIT` fragment overlaps an annotated intron and
# matches in every other way, so it DOES deposit — the only thing missing is the sequenced motif, and
# `sj_implicit` records that. But when several candidate transcripts imply DIFFERENT INTRON SETS, the
# implied set fixes `L`, both quanta, the pool bin, the segment list and therefore which boundaries are
# crossed. There is no partial answer: it cannot deposit spliced (which sj is the unknown), and it
# cannot deposit unspliced either, because `L` involves an intron and does not fit the length
# distribution unless one candidate intron is region_bound out — the very choice in doubt. Forcing a choice is
# choosing an `L` at random. So it deposits NOTHING and waits for the second pass, which has the
# fragment length AND the strand to discriminate with.


def test_TWO_SURVIVING_HYPOTHESES_deposit_on_NOTHING_and_are_BUFFERED_WHOLE():
    """⛔ The whole point: an undetermined path is not a partial deposit, and not a loss either.

    The fragment used here would otherwise deposit richly — it crosses boundaries, uses an annotated sj
    and lands in a length pool — so a leak into any one bank is visible. And it must be RETAINED: this is
    the population the second pass drains, so a silent drop would understate what that pass owes.
    """
    acc = _acc(sj=[SJ])
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
        "region_contained_count",
        "region_start_count",
        "boundary_unspliced_count",
        "boundary_spliced_count",
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
    acc = _acc(sj=[SJ])
    outcome = acc.deposit(
        0, 150, 950, sj_strand=Strand.POS, hypotheses=(GapHypothesis(((201, 900),)),)
    )
    assert outcome is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.sj_count[0, 0]) == 1
    assert int(t.region_start_count.sum()) == 1
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
    acc = _acc(max_fragment_length=500, sj=[SJ])
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
#   POS / NEG  one definite observed strand      ->  must agree with the sj boundary's own strand
#   AMBIGUOUS  the two mates' tags DISAGREE      ->  contradictory evidence; trust no splice
#
# ⚠ NONE must stay permissive. Aligners differ — STAR writes ``XS``, minimap2 writes ``ts``, and some
# write neither — so on an untagged BAM every spliced fragment arrives with NONE. Requiring a strand
# there would delete the entire spliced-RNA signal for that aligner.


def test_a_MISSING_sj_strand_MATCHES_on_coordinates_alone():
    """⛔ The case that makes untagged aligners work at all, so it is pinned before the two below.

    An aligner that writes neither ``XS`` nor ``ts`` gives every spliced fragment
    ``sj_strand = NONE``. If the sj lookup demanded a strand, that BAM would lose 100 % of its
    annotated sj — and the loss would look like a stale annotation, not a convention bug.
    """
    acc = _acc(sj=[SJ])  # (0, 201, 900, POS)
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.NONE)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 1
    assert t.qc["unannotated_introns"] == 0
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 1, "certified RNA"


def test_an_AMBIGUOUS_sj_strand_is_CONTRADICTORY_and_credits_NO_sj():
    """⛔ AMBIGUOUS is contradictory evidence, not missing evidence — and it is neither of the two things
    the original rule could express.

    ``sj_strand`` is the OR of a per-RECORD tag, so AMBIGUOUS (``POS | NEG``) means **the two mates
    disagreed about the same molecule**. That is a data-quality signal of the same family as mates agreeing
    in reference orientation, so the splice must not be trusted: no sj is credited and the fragment
    deposits on the unspliced channel, which is the safe direction the design already takes for
    unannotated sj.

    ⚠ It must NOT be counted as an unannotated intron. That counter's whole purpose is measuring annotation
    coverage, so feeding it alignment disagreements makes the metric report a stale annotation whenever the
    aligner is inconsistent. It gets its own denominator.

    ⚠ Reachable today with no change of ours: ``collect_implicit_splice_introns`` stamps each PE gap's
    intron with the first matching candidate transcript's strand and the caller ORs them, so a two-gap
    fragment matching opposite-strand transcripts arrives here as AMBIGUOUS.
    """
    acc = _acc(sj=[SJ])  # (0, 201, 900, POS)
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
    assert int(t.boundary_unspliced_count.sum()) > 0, "it competes with gDNA, the safe direction"


def test_a_DEFINITE_but_WRONG_sj_strand_still_misses():
    """The third arm, so the rule above cannot be over-applied: a definite strand that disagrees with the
    sj boundary's own strand is a real disagreement, and it IS an unannotated intron — that coordinate
    pair is not annotated on the strand it was observed on."""
    acc = _acc(sj=[SJ])  # (0, 201, 900, POS)
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count.sum()) == 0
    assert t.qc["unannotated_introns"] == 1
    assert t.qc["contradictory_sj_strand"] == 0
    assert int(t.boundary_spliced_count.sum()) == 0, "not certified RNA — the strand disagreed"


# ---------------------------------------------------------------------------
# the strand collapse — and the ONE channel that does not collapse
# ---------------------------------------------------------------------------


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
    region = _region(0, 3)
    # the count still separates them, one per column
    assert int(t.region_contained_count[region, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.region_contained_count[region, STRAND_COLUMNS[Strand.NEG]]) == 1
    # ...and the moments pool them into the single slot
    assert close(float(t.region_contained_inv_opportunity_sum[region]), 2 * _contained_quantum(0, 3, 160), 2)

    boundary_acc = _acc()
    boundary_acc.deposit(0, 120, 320, align_strand=Strand.POS)
    boundary_acc.deposit(0, 120, 320, align_strand=Strand.NEG)
    e = boundary_acc.tally
    boundary = _boundary(0, 2)
    assert int(e.boundary_unspliced_count[boundary].sum()) == 2
    assert close(float(e.boundary_unspliced_inv_length_sum[boundary]), 2.0 / (200 - 1), 2)


def test_the_density_FIELD_NAME_is_gone_everywhere():
    """The rename is complete, so no consumer can reach a half-migrated schema.

    ``inv_length_sum`` is an exact density at an boundary and is NOT one at a region; keeping the old name
    would put one word on two concepts, which is.
    """
    t = _acc().tally
    stale = [name for name in t.__slots__ if name.endswith("_density")]
    assert stale == []
    assert {"region_contained_inv_opportunity_sum"} <= set(t.__slots__)


# ---------------------------------------------------------------------------
# the START / END / SPAN region banks (2026-08-21) — the total-abundance carriers
# ---------------------------------------------------------------------------
# S and E have opportunity ℓ for EVERY fragment length (the walls are at template ends and are the
# CONSUMER's problem, side-selected there); V is a pmf functional per component BY DESIGN. A
# contained fragment is START ∧ END in one region and never SPAN. The ledger now closes twice over.


def test_every_deposited_fragment_has_exactly_one_START_and_one_END():
    """⭐ THE LEDGER, TWICE OVER: ΣS == ΣE == deposited — summed over BOTH strand columns, across a
    mixed population (contained, crossing, spliced) on two references."""
    acc = _acc(sj=[SJ])
    assert acc.deposit(0, 120, 180) is DepositOutcome.DEPOSITED  # contained in r0
    assert acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS) \
        is DepositOutcome.DEPOSITED  # spliced
    assert acc.deposit(0, 150, 260, align_strand=Strand.NEG) is DepositOutcome.DEPOSITED  # crossing
    assert acc.deposit(1, 20, 60) is DepositOutcome.DEPOSITED  # second reference
    t = acc.tally
    deposited = t.qc[DepositOutcome.DEPOSITED.value]
    assert int(t.region_start_count.sum()) == deposited
    assert int(t.region_end_count.sum()) == deposited


def test_START_and_END_are_the_PATHS_covered_bases_never_the_extent():
    """⚠ A leading intron means the molecule does not begin at ``start`` — the covered-base rule the
    start bank already keeps must hold for the END bank symmetrically (trailing intron)."""
    acc = _acc(sj=[SJ])
    # extent [150, 950) but the intron (201, 900) makes the path [150,201)+[900,950):
    # START in the region of base 150, END in the region of base 949 — never the extent's midpoints.
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    rb = np.asarray(CHR1_REGION_BOUNDS)
    r_start = int(np.searchsorted(rb, 150, side="right")) - 1
    r_end = int(np.searchsorted(rb, 949, side="right")) - 1
    assert t.region_start_count[r_start].sum() == 1
    assert t.region_end_count[r_end].sum() == 1
    assert t.region_end_count[r_start].sum() == 0, "END booked at the extent's start region"


def test_a_CONTAINED_fragment_increments_START_and_END_in_its_own_region_and_never_SPAN():
    acc = _acc()
    acc.deposit(0, 120, 180)  # wholly inside r0 = [100, 200) per CHR1_REGION_BOUNDS
    t = acc.tally
    r = int(np.argmax(t.region_contained_count.sum(axis=1)))
    assert t.region_start_count[r].sum() == 1
    assert t.region_end_count[r].sum() == 1
    assert t.region_span_count.sum() == 0, "a contained fragment spans nothing"
    assert (t.region_contained_count.sum(axis=1) <= np.minimum(
        t.region_start_count.sum(axis=1), t.region_end_count.sum(axis=1))).all()


def test_SPAN_counts_regions_STRICTLY_covered_and_a_JUMPED_region_gets_NOTHING():
    """⭐ span = every base covered by ONE segment, neither path endpoint inside. A region the path
    JUMPS over (intron exactly = the region) is not covered and must read 0 on every bank."""
    bounds = [0, 100, 200, 300, 600]
    part = Partition.from_region_bounds([bounds], sj=[(0, 200, 300, Strand.POS)])
    acc = Accumulator(part, max_fragment_length=10_000)
    # unspliced [50, 350): starts r0, ends r3, strictly covers r1 [100,200) and r2 [200,300)
    assert acc.deposit(0, 50, 350) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert t.region_span_count[1].sum() == 1 and t.region_span_count[2].sum() == 1
    assert t.region_span_count[0].sum() == 0 and t.region_span_count[3].sum() == 0
    # spliced [150, 350) with intron exactly r2 = [200, 300): r2 is JUMPED — no span, no start, no end
    acc2 = Accumulator(part, max_fragment_length=10_000)
    assert acc2.deposit(0, 150, 350, observed_introns=[(200, 300)], sj_strand=Strand.POS) \
        is DepositOutcome.DEPOSITED
    t2 = acc2.tally
    assert t2.region_span_count[2].sum() == 0, "a jumped region was counted as spanned"
    assert t2.region_start_count[1].sum() == 1 and t2.region_end_count[3].sum() == 1


def test_the_START_and_END_banks_carry_the_ALIGN_STRAND_column():
    """⛔ perturb-able: a NEG deposit must land column 1 on BOTH banks — a bank that sums the strands
    or fixes column 0 cannot pass."""
    acc = _acc()
    acc.deposit(0, 120, 180, align_strand=Strand.NEG)
    t = acc.tally
    assert t.region_start_count[:, STRAND_COLUMNS[Strand.NEG]].sum() == 1
    assert t.region_start_count[:, STRAND_COLUMNS[Strand.POS]].sum() == 0
    assert t.region_end_count[:, STRAND_COLUMNS[Strand.NEG]].sum() == 1
    assert t.region_end_count[:, STRAND_COLUMNS[Strand.POS]].sum() == 0


def test_SPAN_matches_the_strict_span_opportunity_ABSOLUTELY_on_a_uniform_field():
    """⭐ one fragment at every admissible start (rho = 1 per length): the interior region r1 of
    length ℓ = 100 must read V = (w − ℓ − 1)₊ EXACTLY — 0 at w = ℓ+1, w−101 above."""
    bounds = [0, 1000, 1100, 2100]
    part = Partition.from_region_bounds([bounds])
    ell = 100
    for w, want in ((ell, 0), (ell + 1, 0), (ell + 2, 1), (ell + 50, 49), (300, 199)):
        acc = Accumulator(part, max_fragment_length=10_000)
        for s in range(0, 2100 - w + 1):
            acc.deposit(0, s, s + w)
        got = int(acc.tally.region_span_count[1].sum())
        assert got == want, f"w={w}: V_r1 = {got}, want (w-ell-1)+ = {want}"
