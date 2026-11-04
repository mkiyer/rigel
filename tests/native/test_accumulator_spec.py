"""The accumulator's deposit rule, case by case: the matrix the reference and the native build are both held to.

``_accumulator_reference.py`` is the executable specification, and the native accumulator must reproduce it byte for
byte; this module defines "correct" for both. The genome is a graph whose regions are half-open intervals tiling each
reference, whose contiguous boundaries are the 0-bp gaps between adjacent regions, and whose sj boundaries are directed
donor→acceptor links from the annotation. A fragment is a path (its aligned blocks joined across mate gaps and broken by
introns) of length ``L``, the total of its segments. Regions count fragments contained; boundaries count fragments
crossing; every path books its first and last covered base and the regions it strictly spans. Each population stores
only the channels something reads: an integer count, a float64 reciprocal-opportunity sum ``Σ 1/A(w)``, and on the
contiguous boundaries the conserved mass, which sums to one per fragment.

Every crossed boundary receives the full weight. The chance that a length-``L`` fragment crosses a given boundary is
proportional to ``L − 1`` while the deposit is ``1/(L − 1)``, so the two cancel and every fragment length contributes
equally; dividing by the number of boundaries crossed would break that and make the answer depend on region spacing.
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

#: Fractions are float64, so a sum of `n` round-to-nearest additions differs from the exact answer by at most `n` ulp.
EPS = float(np.finfo(np.float64).eps)


def close(got: float, want: float, deposits: int) -> bool:
    """`got == want` to within the representation error of `deposits` additions."""
    return abs(got - want) <= max(abs(want), 1.0) * deposits * EPS


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

#: chr1 region_bounds   0    100   200   201   400   900   1000
#: regions        n0    n1    n2*   n3    n4    n5      (* n2 is 1 bp: [200,201))
#: boundaries           1     2     3     4     5          (local bound index)
CHR1_REGION_BOUNDS = [0, 100, 200, 201, 400, 900, 1000]
CHR2_REGION_BOUNDS = [0, 500, 1000]

#: coarse region type: 0 intergenic, 1 intron, 2 exon
CHR1_TYPES = [0, 2, 2, 1, 2, 0]
CHR2_TYPES = [0, 2]

#: an annotated intron whose endpoints are bounds 3 and 5, so it swallows the boundary at bound 4
SJ = (0, 201, 900, Strand.POS)


def _partition(sj=()):
    return Partition.from_region_bounds(
        [CHR1_REGION_BOUNDS, CHR2_REGION_BOUNDS], region_types=[CHR1_TYPES, CHR2_TYPES], sj=sj
    )


def _acc(sj=(), **kw):
    return Accumulator(_partition(sj), **kw)


def _boundary(ref, boundary):
    """Global contiguous-boundary id of the boundary at local bound index ``boundary``."""
    return (0 if ref == 0 else len(CHR1_REGION_BOUNDS) - 2) + boundary - 1


def _region(ref, local):
    return (0 if ref == 0 else len(CHR1_REGION_BOUNDS) - 1) + local


# ---------------------------------------------------------------------------
# the reciprocal-opportunity density
# ---------------------------------------------------------------------------
#
# The tests below check the identity, not the representation: a fraction is a plain float64.


def test_one_fragment_recovers_its_own_reciprocal_length():
    acc = _acc()
    acc.deposit(0, 120, 320)
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 1
    assert close(float(t.boundary_unspliced_inv_length_sum[_boundary(0, 2)]), 1.0 / (200 - 1), 1)


# ---------------------------------------------------------------------------
# contained / crossing
# ---------------------------------------------------------------------------


def test_a_contained_fragment_touches_one_region_and_no_boundary():
    acc = _acc()
    assert acc.deposit(0, 220, 380) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(t.region_contained_count[_region(0, 3), 0]) == 1
    assert t.region_contained_count.sum() == 1
    assert t.boundary_unspliced_count.sum() == 0


def test_a_fragment_ending_at_a_boundary_does_not_cross_it():
    """A boundary is 0 bp wide: crossing needs a base on each side."""
    acc = _acc()
    acc.deposit(0, 120, 200)
    assert acc.tally.boundary_unspliced_count.sum() == 0
    assert int(acc.tally.region_contained_count[_region(0, 1), 0]) == 1


def test_a_fragment_starting_at_a_boundary_does_not_cross_it():
    acc = _acc()
    acc.deposit(0, 201, 390)
    assert acc.tally.boundary_unspliced_count.sum() == 0
    assert int(acc.tally.region_contained_count[_region(0, 3), 0]) == 1


def test_a_fragment_crossing_four_regions_credits_three_boundaries_at_full_weight():
    acc = _acc()
    acc.deposit(0, 150, 500)  # touches n1 n2 n3 n4 -> boundaries at 200, 201, 400
    t = acc.tally
    assert [int(t.boundary_unspliced_count[_boundary(0, j), 0]) for j in (1, 2, 3, 4, 5)] == [
        0,
        1,
        1,
        1,
        0,
    ]
    quantum = 1.0 / (350 - 1)
    assert all(
        close(float(t.boundary_unspliced_inv_length_sum[_boundary(0, j)]), quantum, 1)
        for j in (2, 3, 4)
    )
    assert t.region_contained_count.sum() == 0


def test_a_fragment_covering_a_1bp_region_credits_both_boundaries_and_conserves_mass():
    """1 bp regions are legal (there are thousands at human scale), and nothing may assume a region is longer than one
    base. It is the sharpest case for the conserved mass: the 1 bp slice is shared between two bounding boundaries, so
    the three slices are 50 / 1 / 99 bases of a 150 bp molecule and must still sum to exactly one fragment. A count
    would survive any error in how much of the fragment was attributed; the mass does not.
    """
    acc = _acc()
    acc.deposit(0, 150, 300)
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 1
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 1
    assert int(t.region_contained_count[_region(0, 2), 0]) == 0
    mass = int(t.boundary_unspliced_mass.sum())
    assert abs(mass - 1.0) <= 8 * EPS, f"the fragment deposited {mass} fragments, not 1"


def test_a_fragment_longer_than_a_region_is_not_contained_and_crosses_one_boundary():
    """``contained`` needs ``L <= region``. One base longer and the fragment is a crossing instead, touching the region
    axis nowhere."""
    acc = _acc()
    acc.deposit(0, 200, 202)  # n2 is 1 bp, L = 2
    t = acc.tally
    assert t.region_contained_count.sum() == 0
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 1


# ---------------------------------------------------------------------------
# splicing
# ---------------------------------------------------------------------------


def test_a_spliced_jump_deposits_nothing_on_the_boundaries_it_splices_over():
    """The intron [201,900) swallows the boundary at 400; a rule that asked only "does another slice follow?" could not
    tell that from a contiguous crossing."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    length = (950 - 150) - (900 - 201)
    assert int(t.sj_count[0, 0]) == 1
    assert close(float(t.sj_inv_length_sum[0]), 1.0 / (length - 1), 1)
    assert int(t.boundary_unspliced_count[_boundary(0, 4), 0]) == 0, "the swallowed boundary at 400"
    assert int(t.boundary_spliced_count[_boundary(0, 4), 0]) == 0


def test_a_spliced_fragments_own_block_crossings_go_in_the_spliced_bank():
    """A spliced fragment is certified RNA (gDNA cannot be spliced), so a boundary its block crosses is the cleanest RNA
    marker available at a boundary."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert int(t.boundary_spliced_count[_boundary(0, 2), 0]) == 1, (
        "block [150,201) crosses the boundary at 200"
    )
    assert t.boundary_unspliced_count.sum() == 0


def test_a_multi_segment_unspliced_fragment_conserves_its_mass_across_both_segments():
    """An unannotated intron sits strictly inside n4 = [400,900), so segment 1 crosses the boundary at 400 and segment 2
    the boundary at 900. Every segment touches a boundary, so the conservation law applies unchanged and the mass sums
    to exactly one fragment; a segment touching no boundary would leave its bases nowhere conserved to go.
    """
    acc = _acc()
    acc.deposit(0, 380, 950, observed_introns=[(500, 600)])
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 4), 0]) == 1, "boundary 400, from segment 1"
    assert int(t.boundary_unspliced_count[_boundary(0, 5), 0]) == 1, "boundary 900, from segment 2"
    mass = int(t.boundary_unspliced_mass.sum())
    assert abs(mass - 1.0) <= 8 * EPS, f"a two-segment fragment deposited {mass} fragments, not 1"


def test_an_unannotated_intron_credits_no_sj_and_nothing_across_the_gap():
    """Unannotated sj are disproportionately artifactual, so they deposit on the unspliced channel and compete with gDNA
    rather than being certified RNA."""
    acc = _acc(sj=[SJ])
    # [200,400) is not annotated; it swallows the boundary at 201
    acc.deposit(0, 50, 500, observed_introns=[(200, 400)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.sj_count.sum() == 0
    assert t.boundary_spliced_count.sum() == 0, "not certified RNA — it competes with gDNA"
    assert int(t.boundary_unspliced_count[_boundary(0, 1), 0]) == 1, (
        "block [50,200) crosses the boundary at 100"
    )
    assert int(t.boundary_unspliced_count[_boundary(0, 3), 0]) == 0, "the swallowed boundary at 201"


def test_a_fragment_straddling_two_regions_without_a_crossing_is_not_contained():
    """An unannotated intron can swallow every boundary between two blocks. The fragment then crosses nothing yet
    straddles two regions, so crediting it as contained would put its whole length in a region it only partly overlaps.
    It deposits on no object, and the start count keeps that visible.
    """
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
    assert int(t.deposited_lengths[180 - 40]) == 1


def test_opposite_strand_sj_at_the_same_coordinates_are_distinct_boundaries():
    """Biologically impossible (splice motifs are not palindromic), so only a synthetic stress test can reach it."""
    acc = _acc(sj=[(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert t.sj_count.sum() == 1
    assert int(t.sj_count[1, STRAND_COLUMNS[Strand.NEG]]) == 1, (
        "the NEG boundary (id 1), genome minus"
    )


@pytest.mark.parametrize("order", [("POS_first", 1), ("NEG_first", -1)])
def test_a_sj_id_is_a_function_of_the_partition_not_of_argument_order(order):
    """The sj-boundary id is the rank in the sort, so the sort must be total; otherwise the id depends on the order the
    caller listed the sj in.

    This is the only test that pins ``strand`` as part of the sort key. ``np.lexsort`` is stable, so ``(acceptor,
    donor)`` alone is right whenever ties already arrive in order, as in every other test and both real indexes.
    Reversing the argument order exposes a missing key, and a strand-coincident pair is the only possible tie.
    """
    _, direction = order
    sj = [(0, 201, 900, Strand.POS), (0, 201, 900, Strand.NEG)][::direction]
    part = Partition.from_region_bounds([CHR1_REGION_BOUNDS], region_types=[CHR1_TYPES], sj=sj)
    assert [int(s) for s in part.sj_strand] == [int(Strand.POS), int(Strand.NEG)], (
        "POS must sort to slot 0 whichever order it was passed in"
    )


def test_a_fragment_using_two_sj_credits_both():
    """Each boundary owns its own expectation, and the strand model is fitted from a separate scan output, so crediting
    every sj a path uses distorts nothing. The two introns are separated by a real exon (abutting introns are
    malformed), so this needs its own partition.
    """
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


# ---------------------------------------------------------------------------
# strand channels
# ---------------------------------------------------------------------------


def test_the_unspliced_bank_is_indexed_by_genome_strand():
    acc = _acc()
    acc.deposit(0, 150, 300, align_strand=Strand.NEG)
    t = acc.tally
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 1]) == 1
    assert int(t.boundary_unspliced_count[_boundary(0, 2), 0]) == 0


def test_every_bank_including_the_sj_is_indexed_by_genome_strand():
    """One convention throughout. Sense/antisense is derived, never stored: the sj boundary carries its own genomic
    strand, so a consumer computes ``sense = (fragment strand == sj strand)``. Here a genome-minus fragment splices
    across a ``+`` sj; a sense convention would put it in the antisense column and the genome convention in the minus
    column, which are the same index, so the next test is the discriminating one.
    """
    acc = _acc(sj=[SJ])
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    assert int(acc.tally.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1


def test_the_sj_strand_split_is_retained_for_aligner_artifact_detection():
    """``sj_count`` keeps both aligned-strand columns. A sj is stranded by its splicing motif, so the strand of the
    fragments on it looks redundant, and every consumer sums the two columns; the split is kept for aligner-artifact
    detection.

    Aligners emit false-positive ``N`` CIGAR ops from plain genomic DNA. ``rigel.splice_blacklist`` catches those
    already listed by coordinate, which is far from complete. In a stranded library real sj inherit the library's strand
    specificity, while alignment artifacts deposit onto both strands and deviate from it, so the per-sj split by aligned
    strand is an empirical detector. Unstranded data cannot use it (with κ = ½ there is nothing to deviate from), but
    that is a property of the detector. The clean and the artifactual sj below carry the same total, so only the split
    tells them apart.
    """
    # a clean sj: every fragment on one aligned strand, as a stranded library produces
    clean = _acc(sj=[SJ])
    for _ in range(20):
        clean.deposit(
            0,
            150,
            950,
            observed_introns=[(201, 900)],
            align_strand=Strand.NEG,
            sj_strand=Strand.POS,
        )
    # an artifactual sj: the aligner emitted it from both strands
    artifact = _acc(sj=[SJ])
    for i in range(20):
        artifact.deposit(
            0,
            150,
            950,
            observed_introns=[(201, 900)],
            align_strand=Strand.NEG if i % 2 else Strand.POS,
            sj_strand=Strand.POS,
        )

    clean_row = clean.tally.sj_count[0]
    artifact_row = artifact.tally.sj_count[0]
    assert list(clean_row) == [0, 20], "a clean sj sits entirely in one aligned-strand column"
    assert list(artifact_row) == [10, 10], "an artifact splits across both"

    # Collapsing the columns destroys the difference.
    assert int(clean_row.sum()) == int(artifact_row.sum()) == 20, (
        "the two sj carry the same TOTAL, so a one-column sj_count cannot distinguish them — "
        "which is exactly why the strand split is retained"
    )


def test_the_sj_mass_keeps_the_same_split_and_it_is_the_only_mass_that_does():
    """The mass twin of the test above. Mass banks carry one value wherever nothing reads a mass per strand; the sj mass
    is the exception, because the artifact detector needs the split and only a mass can separate a sj used by many short
    fragments from one used by few long ones, in a single pass. The mass must land in the same column as the count, or
    ``mass[c]/count[c]`` would divide two different populations. ``boundary_unspliced_mass`` and
    ``boundary_spliced_mass`` have no such consumer and stay one-column.
    """
    pos, neg = STRAND_COLUMNS[Strand.POS], STRAND_COLUMNS[Strand.NEG]

    def one(align):
        acc = _acc(sj=[SJ])
        acc.deposit(
            0, 150, 950, observed_introns=[(201, 900)], align_strand=align, sj_strand=Strand.POS
        )
        return acc.tally

    t_pos, t_neg = one(Strand.POS), one(Strand.NEG)
    assert t_pos.sj_mass.ndim == 2, "the sj mass carries a strand"

    # each deposit reaches its own column and only its own
    assert t_pos.sj_mass[0, pos] > 0.0 and t_pos.sj_mass[0, neg] == 0.0
    assert t_neg.sj_mass[0, neg] > 0.0 and t_neg.sj_mass[0, pos] == 0.0
    # and it is the column the count used
    assert int(t_pos.sj_count[0, pos]) == 1 and int(t_neg.sj_count[0, neg]) == 1
    # the two agree because the deposit rule has no strand in it; a column-dependent share would show here.
    assert close(float(t_pos.sj_mass[0, pos]), float(t_neg.sj_mass[0, neg]), 12)
    # summing the columns returns what a one-column bank would hold, so every consumer below `substrate` can ignore the
    # split.
    assert close(float(t_pos.sj_mass.sum()), float(t_pos.sj_mass[0, pos]), 12)


def test_a_sense_fragment_on_the_minus_strand_is_booked_as_minus():
    """Sense-to-motif would say column 0; genome strand says column 1."""
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


def test_a_fragment_over_the_length_limit_deposits_nothing_and_is_counted():
    """Bounded influence: without the limit a handful of enormous read groups own almost every boundary crossing in a
    real library. The rejection is counted, so it is not hidden.
    """
    acc = _acc(max_fragment_length=200)
    assert acc.deposit(0, 100, 500) is DepositOutcome.TOO_LONG
    t = acc.tally
    assert t.boundary_unspliced_count.sum() == 0
    assert t.region_start_count.sum() == 0
    assert t.qc["dropped_too_long"] == 1


def test_the_limit_applies_to_the_length_not_the_span():
    """A 300 bp molecule across a 10 kb intron has a 10 kb span, so limiting the span would discard a large share of the
    spliced population while limiting ``L`` discards a small one.
    """
    acc = _acc(sj=[SJ], max_fragment_length=200)
    out = acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
    assert out is DepositOutcome.DEPOSITED, "span 800, L = 101"
    assert int(acc.tally.sj_count[0, 0]) == 1


def test_a_fragment_is_clipped_to_its_reference_and_its_length_is_the_clipped_one():
    acc = _acc()
    acc.deposit(0, 950, 1200)  # chr1 ends at 1000
    t = acc.tally
    assert int(t.region_contained_count[_region(0, 5), 0]) == 1
    assert int(t.deposited_lengths[50]) == 1


def test_a_single_region_reference_has_no_boundaries_and_accepts_a_fragment():
    acc = Accumulator(Partition.from_region_bounds([[0, 1000]], region_types=[[0]]))
    assert acc.n_boundaries == 0
    acc.deposit(0, 100, 300)
    assert int(acc.tally.region_contained_count[0, 0]) == 1


def test_the_per_reference_offsets_do_not_bleed():
    """chr1's fragment crosses its boundaries 2 and 3; chr2's crosses chr2's boundary 1. Nothing lands on a reference it
    did not come from.
    """
    acc = _acc()
    acc.deposit(0, 150, 300)  # crosses the boundaries at 200 and 201
    acc.deposit(1, 400, 700)  # crosses chr2's boundary at 500
    t = acc.tally
    assert [int(t.boundary_unspliced_count[e, 0]) for e in range(acc.n_boundaries)] == [
        0,
        1,
        1,
        0,
        0,
        1,
    ]
    assert int(t.boundary_unspliced_count.sum()) == 3
    assert int(t.region_start_count[_region(0, 1)].sum()) == 1
    assert int(t.region_start_count[_region(1, 0)].sum()) == 1


# ---------------------------------------------------------------------------
# the start-count invariant, checkable against the scanner
# ---------------------------------------------------------------------------


def test_every_accepted_fragment_increments_exactly_one_start_count():
    """The crossing and contained totals can only be checked by re-running the deposit; the start count is checkable
    against a number the scanner knows independently.
    """
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
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)  # annotated sj
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTERGENIC].sum()) == 2
    assert int(p[FragmentPool.DNA_INTRONIC].sum()) == 1
    assert int(p[FragmentPool.DNA_INTRON_EXON].sum()) == 1
    assert int(p[FragmentPool.RNA_SPLICED].sum()) == 1
    assert int(p[FragmentPool.DNA_INTERGENIC_EXON].sum()) == 0


def test_a_pool_is_binned_at_the_length_and_once_per_fragment():
    acc = _acc()
    acc.deposit(0, 210, 390, observed_introns=[(300, 340)])  # L = 180 - 40
    p = acc.tally.pool_lengths
    assert int(p[FragmentPool.DNA_INTRONIC, 140]) == 1
    assert int(p.sum()) == 1


# Pool membership is determinacy, not provenance: a fragment reaches a pool only when exactly one hypothesis survived,
# so its `L` is not in doubt however it was derived. Barring a splice that was inferred rather than sequenced would
# preferentially bar fragments whose mates sit far apart, filtering the pool by length.


def test_a_multi_boundary_crossing_enters_no_pool():
    """A splash read straddles one probe boundary. A fragment crossing several boundaries has no single structural
    class, and an impure pool is worse than a missing one."""
    acc = _acc()
    acc.deposit(0, 150, 500)
    assert int(acc.tally.pool_lengths.sum()) == 0


# ---------------------------------------------------------------------------
# the two arithmetic identities the design rests on
# ---------------------------------------------------------------------------


def _corpus(rng, n, ref_len):
    lengths = np.clip(rng.normal(200.0, 60.0, n).round().astype(np.int64), 30, 600)
    starts = rng.integers(0, ref_len - lengths)
    return starts, starts + lengths, lengths


def _uniform_accumulator(region_bp, ref_len):
    region_bounds = list(range(0, ref_len + 1, region_bp))
    return Accumulator(
        Partition.from_region_bounds(
            [region_bounds], region_types=[[0] * (len(region_bounds) - 1)]
        ),
        max_fragment_length=10_000,
    )


@pytest.mark.parametrize("region_bp", [50, 200, 1000])
def test_the_crossing_density_recovers_the_true_density_with_no_length_model(region_bp):
    """``E[Σ 1/(L−1)] = ρ`` exactly, for any fragment-length distribution: the identity that lets a boundary need no
    divisor and no length model. It must hold at every region spacing, so the test sweeps the spacing; dividing the
    weight by the number of boundaries crossed would shrink the answer as regions get closer.
    """
    ref_len, rho = 200_000, 0.05
    acc = _uniform_accumulator(region_bp, ref_len)
    rng = np.random.default_rng(7)
    starts, ends, _ = _corpus(rng, int(rho * ref_len), ref_len)
    for s, e in zip(starts, ends):
        acc.deposit(0, int(s), int(e))
    interior = slice(5, acc.n_boundaries - 5)
    estimate = acc.tally.boundary_unspliced_inv_length_sum[interior].sum() / (acc.n_boundaries - 10)
    assert 0.98 <= estimate / rho <= 1.02, f"{estimate / rho:.4f} at {region_bp} bp regions"


def test_the_crossing_count_recovers_density_times_mean_length():
    """The companion identity ``E[count] = ρ·(E[L] − 1)``. With the density above it forms the 2×2 that separates two
    components by fragment length alone."""
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


def test_the_deposit_is_independent_of_the_order_fragments_arrive_in():
    """Integer addition is associative, so the per-worker merge of the count banks is bit-identical at any thread count."""
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
        "region_start_count",
        "boundary_unspliced_count",
        "boundary_unspliced_inv_length_sum",
        "pool_lengths",
    ):
        got, want = getattr(a, field), getattr(b, field)
        # Exact for the integer banks; the float64 fractions are re-associated by the shuffle, so their tolerance is the
        # representation's, derived from the deposit count.
        if getattr(got, "dtype", None) == np.float64:
            assert np.allclose(got, want, rtol=len(starts) * EPS, atol=0.0), field
        else:
            assert np.array_equal(got, want), field


# ---------------------------------------------------------------------------
# malformed intron lists, and the one definition of L
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "introns,expected_length,expected_crossings",
    [
        ([(210, 260), (240, 300)], 360, 4),  # overlapping
        ([(200, 400), (250, 300)], 250, 1),  # nested
        ([(150, 480), (160, 470)], 120, 1),  # wide overlap: `span − Σ intron` would go negative
        ([(210, 260), (210, 260)], 400, 4),  # duplicated, as two disagreeing mates produce
        ([(300, 300)], 450, 4),  # zero length
        ([(210, 260), (300, 340)], 360, 4),  # well formed: nothing merged
    ],
)
def test_the_length_is_the_path_segments_total_even_with_malformed_introns(
    introns, expected_length, expected_crossings
):
    """One definition of ``L``: the total length of the path's segments.

    A second formula, ``span − Σ(intron lengths)``, disagrees as soon as introns overlap, and on a wide overlap goes
    negative, discarding a good fragment. Real BAMs produce such input: the scanner reads ``XS`` once per record, so a
    pair whose mates disagree about an acceptor yields overlapping introns for one molecule. The native code sums its
    segments too, so byte-identity depends on both sides using this one rule.
    """
    acc = _acc(max_fragment_length=10_000)
    assert acc.deposit(0, 50, 500, observed_introns=introns) is DepositOutcome.DEPOSITED
    t = acc.tally
    crossings = int(t.boundary_unspliced_count.sum())
    assert crossings == expected_crossings
    assert close(
        float(t.boundary_unspliced_inv_length_sum.sum()),
        crossings / (expected_length - 1),
        crossings,
    )


def test_the_path_starts_where_its_first_covered_base_is_not_where_the_extent_begins():
    """A leading intron means the molecule does not begin at ``lo``. Attributing it to the region containing ``lo``
    would credit the start count, and possibly the contained deposit, to a region the fragment never touches.
    """
    acc = _acc(max_fragment_length=10_000)
    acc.deposit(0, 150, 500, observed_introns=[(150, 480)])  # the path is only [480,500), inside n4
    t = acc.tally
    assert int(t.region_start_count[_region(0, 4)].sum()) == 1, "n4, where the path actually starts"
    assert int(t.region_start_count[_region(0, 1)].sum()) == 0, "not n1, where the extent begins"
    assert int(t.region_contained_count[_region(0, 4), 0]) == 1
    assert int(t.deposited_lengths[20]) == 1


def test_a_duplicated_intron_credits_its_sj_once():
    """Two mates reporting the same intron is one splice event, not two."""
    acc = _acc(sj=[SJ])
    acc.deposit(0, 150, 950, observed_introns=[(201, 900), (201, 900)], sj_strand=Strand.POS)
    assert int(acc.tally.sj_count[0, 0]) == 1


def test_abutting_introns_are_malformed_and_merge():
    """Two introns sharing an endpoint imply a zero-length exon between them, which is physically impossible (a
    transcript with one is molecularly identical to one without), so the pair is an alignment artifact. The index cannot
    produce it either: building the exon arrays drops a zero-length exon and fuses its flanking introns. Merged here.
    """
    acc = _acc(sj=[(0, 201, 400, Strand.POS), (0, 400, 900, Strand.POS)])
    acc.deposit(0, 150, 950, observed_introns=[(201, 400), (400, 900)], sj_strand=Strand.POS)
    t = acc.tally
    assert t.sj_count.sum() == 0, "the merged span 201->900 is not an annotated sj"


def test_a_wide_overlap_keeps_a_good_fragment():
    """The ``span − Σ intron`` formula would give a negative ``L`` here and file the fragment as ``dropped_empty``,
    invisible to the start-count invariant because a rejected fragment never reaches it. The path-segment total keeps
    it.
    """
    acc = _acc(max_fragment_length=10_000)
    assert (
        acc.deposit(0, 150, 500, observed_introns=[(150, 480), (160, 470)])
        is DepositOutcome.DEPOSITED
    )
    assert acc.tally.qc["dropped_empty"] == 0
    assert int(acc.tally.region_start_count.sum()) == 1


# ---------------------------------------------------------------------------
# the region banks carry one strand convention
# ---------------------------------------------------------------------------

#: a sj far enough right that the first block still spans region n2 = [200,201)
SPAN_SJ_POS = (0, 400, 900, Strand.POS)
SPAN_SJ_NEG = (0, 400, 900, Strand.NEG)


def test_a_spliced_and_an_unspliced_fragment_of_one_genome_strand_share_a_column():
    """One array, one convention.

    A spliced fragment cannot be contained (both endpoints of an annotated intron are region bounds, so it always
    crosses its sj boundary), but it shares boundary and sj banks with the unspliced fragments beside it. Indexing
    spliced deposits by sense-relative-to-motif while unspliced ones use genome strand would put one array in two
    conventions. The claim is about the convention, so it is asserted on both banks a spliced fragment reaches:
    ``boundary_spliced_count`` beside ``boundary_unspliced_count`` at the same boundary, and ``sj_count``.
    """
    acc = _acc(sj=[SPAN_SJ_POS])
    acc.deposit(0, 150, 300, align_strand=Strand.NEG)  # unspliced, genome minus
    acc.deposit(  # spliced, genome minus, antisense to its + sj
        0, 150, 950, observed_introns=[(400, 900)], align_strand=Strand.NEG, sj_strand=Strand.POS
    )
    t = acc.tally
    # One boundary, two banks, one column: the unspliced and the spliced fragment both cross the boundary at 200 and
    # both are genome-minus. A sense-relative convention would split them.
    assert int(t.boundary_unspliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, (
        "genome minus"
    )
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.NEG]]) == 1, (
        "the spliced one is ANTISENSE to its + sj, and still books genome minus"
    )
    assert int(t.boundary_unspliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.boundary_spliced_count[_boundary(0, 2), STRAND_COLUMNS[Strand.POS]]) == 0
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.NEG]]) == 1, "the sj bank too"


def test_a_spliced_sense_fragment_books_boundary_and_sj_by_genome_strand():
    """Sense-to-motif would say column 0 for both; genome strand says 1."""
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
# `deposit` has no equivalent of the scanner's own strand gate, so a fragment with no single genome strand is rejected
# here, where the rejection is counted, rather than upstream where the loss is invisible.


@pytest.mark.parametrize("undefined", [Strand.NONE, Strand.AMBIGUOUS])
def test_an_undefined_strand_is_rejected_not_silently_booked_as_minus(undefined):
    """The channel is the genome strand, so a fragment without one has no channel. Written as "the plus column if POS
    else the minus column", every strand that is not POS (NONE and AMBIGUOUS too) would land in the minus column,
    crediting the fragment to the wrong strand. The rejection is counted on its own QC denominator instead. AMBIGUOUS is
    reachable in production: ``build_fragment`` keys blocks by ``(ref, strand)``, so mates agreeing in reference
    orientation give ``strand = POS|NEG``.
    """
    acc = _acc()
    assert acc.deposit(0, 150, 300, align_strand=undefined) is DepositOutcome.STRAND_UNDEFINED
    t = acc.tally
    assert t.qc["dropped_strand_undefined"] == 1
    assert int(t.region_contained_count.sum()) == 0, "must not be booked into either column"
    assert int(t.region_start_count.sum()) == 0, "a rejected fragment never reaches the invariant"


# ── an ambiguous path deposits nothing, which is not the same thing as an implicit splice ────────
#
# A fragment whose splice was never sequenced but whose gap has one explanation still deposits: only the motif is
# missing, and the path is not in doubt. When several candidate transcripts imply different intron sets, the implied set
# fixes `L`, both quanta, the pool bin, the segment list and so the boundaries crossed. There is no partial answer: it
# cannot deposit spliced (which sj is unknown) nor unspliced (its `L` involves the intron in doubt). Forcing a choice
# would choose an `L` at random, so it deposits nothing and waits for the second pass, which has the fragment length and
# the strand to discriminate with.


def test_two_surviving_hypotheses_deposit_on_nothing_and_are_buffered_whole():
    """An undetermined path makes no partial deposit and is not lost. The fragment would otherwise deposit richly (it
    crosses boundaries, uses an annotated sj and lands in a length pool), so a leak into any bank shows. It must be
    retained, since this is the population the second pass drains.
    """
    acc = _acc(sj=[SJ])
    outcome = acc.deposit(
        0,
        150,
        950,
        sj_strand=Strand.POS,
        hypotheses=(GapHypothesis(((201, 900),)), GapHypothesis()),
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

    # Retained whole, with every hypothesis: the second pass cannot choose between answers it was not given.
    assert len(t.deferred) == 1
    held = t.deferred[0].placements[0]
    assert (held.ref, held.start, held.end) == (0, 150, 950)
    assert {path.introns for path in held.hypotheses} == {((201, 900),), ()}
    # The genomic path competes against exactly one spliced path: the open question is RNA or gDNA.
    assert t.gap_resolution["gap_deferred_rna_or_gdna"] == 1


def test_one_surviving_hypothesis_deposits_though_its_splice_was_never_sequenced():
    """The pair for the test above: determinacy, not provenance. Both fragments have a splice that was never sequenced;
    this one has only one explanation for its gap, so its ``L`` is not in doubt and it deposits, including into the
    pure-RNA length pool. Barring it on provenance would exclude fragments whose mates sit far apart and bias the pool's
    length distribution.
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


def test_a_strand_undefined_ambiguous_fragment_is_counted_as_strand_undefined():
    """Every fragment counts exactly once, so the precedence is part of the contract. A fragment that is both is filed
    under the strand: the deferred queue sizes the population the second pass can recover, and that pass resolves the
    path, not the strand the read aligned to. Holding the fragment would promise a recovery that cannot happen and put a
    fragment with no column into a queue whose drain needs one.
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


def test_a_hypothesis_longer_than_the_limit_is_ruled_out_and_the_rest_stands():
    """``max_fragment_length`` is the rule that already makes TOO_LONG a rejection. Short-read chemistry does not
    sequence molecules past the limit, so a hypothesis implying a longer ``L`` is not a molecule this library contains.
    Applied to the genomic hypothesis it is the rule "a fragment whose genomic span exceeds the limit must be RNA",
    since the genomic path's ``L`` is that span. Here the span is 800 over a limit of 500, so the genomic path dies and
    the single spliced path (``L`` = 800 − 699 = 101) deposits.
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


def test_when_the_limit_would_rule_out_every_hypothesis_the_survivors_stand():
    """Filtering to nothing would claim "this molecule cannot exist", which the filter cannot conclude: it describes
    what the chemistry sequences, and the fragment is in the BAM. So the survivors stand and the ordinary ``TOO_LONG``
    rejection counts the fragment.
    """
    acc = _acc(max_fragment_length=50)
    outcome = acc.deposit(0, 150, 950, hypotheses=(GapHypothesis(),))
    assert outcome is DepositOutcome.TOO_LONG
    assert acc.tally.qc["dropped_too_long"] == 1
    assert not acc.tally.deferred


# ── the splice-junction motif strand is a three-way distinction ───────────────────────────────────
#
# ``sj_strand`` is the OR of the per-record ``XS``/``ts`` tag values, so it carries three different
# states, each with its own rule:
#
#   NONE       no strand tag in the BAM at all  ->  no information; match on coordinates alone
#   POS / NEG  one definite observed strand      ->  must agree with the sj boundary's own strand
#   AMBIGUOUS  the two mates' tags disagree      ->  contradictory evidence; trust no splice
#
# NONE must stay permissive. STAR writes ``XS``, minimap2 writes ``ts``, and some aligners write
# neither, so on an untagged BAM every spliced fragment arrives with NONE; requiring a strand there
# would delete the whole spliced-RNA signal for that aligner.


def test_a_missing_sj_strand_matches_on_coordinates_alone():
    """An aligner that writes neither ``XS`` nor ``ts`` gives every spliced fragment ``sj_strand = NONE``. If the sj
    lookup demanded a strand, that BAM would lose all its annotated sj, and the loss would look like a stale annotation.
    """
    acc = _acc(sj=[SJ])  # (0, 201, 900, POS)
    acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.NONE)
    t = acc.tally
    assert int(t.sj_count[0, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 1, "certified RNA"


def test_an_ambiguous_sj_strand_is_contradictory_and_credits_no_sj():
    """AMBIGUOUS is contradictory evidence, not missing evidence.

    ``sj_strand`` is the OR of a per-record tag, so AMBIGUOUS (``POS | NEG``) means the two mates disagreed about the
    same molecule, a data-quality signal like mates agreeing in reference orientation. The splice is not trusted: no sj
    is credited and the fragment deposits on the unspliced channel, the safe direction already taken for unannotated sj.
    It is reachable in practice: ``collect_implicit_splice_introns`` stamps each PE gap's intron with the first matching
    candidate transcript's strand and the caller ORs them, so a two-gap fragment matching opposite-strand transcripts
    arrives as AMBIGUOUS.
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
    assert int(t.pool_lengths[FragmentPool.RNA_SPLICED].sum()) == 0, "not certified RNA"
    assert int(t.boundary_unspliced_count.sum()) > 0, "it competes with gDNA, the safe direction"


def test_a_definite_but_wrong_sj_strand_misses():
    """A definite strand that disagrees with the sj boundary's own strand is a real disagreement: that coordinate pair
    is not annotated on the observed strand, so it is an unannotated intron."""
    acc = _acc(sj=[SJ])  # (0, 201, 900, POS)
    acc.deposit(
        0, 150, 950, observed_introns=[(201, 900)], align_strand=Strand.NEG, sj_strand=Strand.NEG
    )
    t = acc.tally
    assert int(t.sj_count.sum()) == 0
    assert int(t.boundary_spliced_count.sum()) == 0, "not certified RNA — the strand disagreed"


# ---------------------------------------------------------------------------
# the strand collapse, and the one channel that does not collapse
# ---------------------------------------------------------------------------


def test_both_genome_strands_land_in_the_one_length_moment_slot():
    """``count`` keeps two genome-strand columns; the length moments keep one, because the strand a read aligned to says
    nothing about whether the molecule was gDNA or RNA. So a plus and a minus fragment at the same object must both
    reach the single slot. A collapse that kept only one strand's deposits would otherwise be caught only by C++ parity,
    which a matching bug on both sides passes.
    """
    acc = _acc()
    acc.deposit(0, 220, 380, align_strand=Strand.POS)
    acc.deposit(0, 220, 380, align_strand=Strand.NEG)
    t = acc.tally
    region = _region(0, 3)
    # the count still separates them, one per column
    assert int(t.region_contained_count[region, STRAND_COLUMNS[Strand.POS]]) == 1
    assert int(t.region_contained_count[region, STRAND_COLUMNS[Strand.NEG]]) == 1

    # and the moments pool them into the single slot
    boundary_acc = _acc()
    boundary_acc.deposit(0, 120, 320, align_strand=Strand.POS)
    boundary_acc.deposit(0, 120, 320, align_strand=Strand.NEG)
    e = boundary_acc.tally
    boundary = _boundary(0, 2)
    assert int(e.boundary_unspliced_count[boundary].sum()) == 2
    assert close(float(e.boundary_unspliced_inv_length_sum[boundary]), 2.0 / (200 - 1), 2)


# ---------------------------------------------------------------------------
# the START region bank: the ledger
# ---------------------------------------------------------------------------
# START has opportunity ℓ for every fragment length (the walls are at template ends). A contained
# fragment books its START in its own region.


def test_every_deposited_fragment_has_exactly_one_start():
    """The ledger: ΣS == deposited, summed over both strand columns, across a mixed population (contained, crossing,
    spliced) on two references."""
    acc = _acc(sj=[SJ])
    assert acc.deposit(0, 120, 180) is DepositOutcome.DEPOSITED  # contained in r0
    assert (
        acc.deposit(0, 150, 950, observed_introns=[(201, 900)], sj_strand=Strand.POS)
        is DepositOutcome.DEPOSITED
    )  # spliced
    assert acc.deposit(0, 150, 260, align_strand=Strand.NEG) is DepositOutcome.DEPOSITED  # crossing
    assert acc.deposit(1, 20, 60) is DepositOutcome.DEPOSITED  # second reference
    t = acc.tally
    deposited = t.qc[DepositOutcome.DEPOSITED.value]
    assert int(t.region_start_count.sum()) == deposited


def test_a_contained_fragment_increments_start_in_its_own_region():
    acc = _acc()
    acc.deposit(0, 120, 180)  # wholly inside r0 = [100, 200) per CHR1_REGION_BOUNDS
    t = acc.tally
    r = int(np.argmax(t.region_contained_count.sum(axis=1)))
    assert t.region_start_count[r].sum() == 1
    assert (t.region_contained_count.sum(axis=1) <= t.region_start_count.sum(axis=1)).all()


def test_the_start_bank_carries_the_align_strand_column():
    """A NEG deposit must land in column 1, so a bank that sums the strands or fixes column 0 fails."""
    acc = _acc()
    acc.deposit(0, 120, 180, align_strand=Strand.NEG)
    t = acc.tally
    assert t.region_start_count[:, STRAND_COLUMNS[Strand.NEG]].sum() == 1
    assert t.region_start_count[:, STRAND_COLUMNS[Strand.POS]].sum() == 0
