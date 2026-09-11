"""The implicit splice — a fragment whose annotated intron lies inside the unsequenced mate gap, so
the splice motif was never read. The first block gates the per-intron ``SPLICED_IMPLICIT``
discriminant, exercising the C++ ``FragmentResolver`` through the standard Python entry point
``rigel.resolution.resolve_fragment`` and the ``make_fragment`` helper over an acceptance matrix of
transcript geometries; one comprehensive GTF fixture built once per session supplies every geometry
the matrix needs, and each test sets the resolver's ``splicing_anchor_tolerance`` explicitly so the
cases do not depend on execution order. It also gates that gap introns are searched whatever the
splice type. The second block gates what the scanner does with such a fragment on a real scan: a
determined path deposits, and an undetermined one is held in the deferred bank with every hypothesis
on it.
"""

from __future__ import annotations

import textwrap

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType
from rigel.types import GenomicInterval, Strand

from _index_builder import build_test_index
from _resolution_reference import make_fragment, resolve_fragment


# Geometry of the fixture GTF below (0-based half-open, after the GTF parse):
#   g1 (+):
#     t_short   : exons (99,200),(299,400)            intron (200,299)
#     t_three   : exons (99,200),(299,400),(499,600)  introns (200,299),(400,499)
#   g2 (+):
#     t_micro   : exons (699,800),(804,900)           micro-intron (800,804)
#   g3 (+):
#     t_long    : exons (999,1100),(51099,52000)      long intron (1100,51099)
#   g4 (+):
#     t_no_intron   : exons (52999,53100)             single-exon
#   g5 (+):
#     t_one_exon: exons (59999,60500)                 single-exon, spans whole locus
#     t_split   : exons (59999,60100),(60399,60500)   intron (60100,60399)
#   g6 (+): the MIXED fragment — an observed splice AND a gap intron
#     t_mixed   : exons (62000,62200),(62400,62600),(62800,63000)
#                 introns (62200,62400) and (62600,62800)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

IMPLICIT_GTF = textwrap.dedent("""\
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_short"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_short"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t500\t600\t.\t+\t.\tgene_id "g1"; transcript_id "t_three"; gene_name "G1"; gene_type "protein_coding";
    chr1\ttest\texon\t700\t800\t.\t+\t.\tgene_id "g2"; transcript_id "t_micro"; gene_name "G2"; gene_type "protein_coding";
    chr1\ttest\texon\t805\t900\t.\t+\t.\tgene_id "g2"; transcript_id "t_micro"; gene_name "G2"; gene_type "protein_coding";
    chr1\ttest\texon\t1000\t1100\t.\t+\t.\tgene_id "g3"; transcript_id "t_long"; gene_name "G3"; gene_type "protein_coding";
    chr1\ttest\texon\t51100\t52000\t.\t+\t.\tgene_id "g3"; transcript_id "t_long"; gene_name "G3"; gene_type "protein_coding";
    chr1\ttest\texon\t53000\t53100\t.\t+\t.\tgene_id "g4"; transcript_id "t_no_intron"; gene_name "G4"; gene_type "protein_coding";
    chr1\ttest\texon\t60000\t60500\t.\t+\t.\tgene_id "g5"; transcript_id "t_one_exon"; gene_name "G5"; gene_type "protein_coding";
    chr1\ttest\texon\t60000\t60100\t.\t+\t.\tgene_id "g5"; transcript_id "t_split"; gene_name "G5"; gene_type "protein_coding";
    chr1\ttest\texon\t60400\t60500\t.\t+\t.\tgene_id "g5"; transcript_id "t_split"; gene_name "G5"; gene_type "protein_coding";
    chr1\ttest\texon\t62001\t62200\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed"; gene_name "G6"; gene_type "protein_coding";
    chr1\ttest\texon\t62401\t62600\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed"; gene_name "G6"; gene_type "protein_coding";
    chr1\ttest\texon\t62801\t63000\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed"; gene_name "G6"; gene_type "protein_coding";
""")


@pytest.fixture(scope="session")
def implicit_index(tmp_path_factory):
    """Build the comprehensive implicit-splice fixture index."""
    return build_test_index(tmp_path_factory, IMPLICIT_GTF, genome_size=70000, name="implicit_idx")


@pytest.fixture
def with_tolerance(implicit_index):
    """Yield (index, set_K) so tests choose K explicitly and we restore at teardown."""
    resolver = implicit_index.resolver

    def _set(k: int) -> None:
        resolver.set_splicing_anchor_tolerance(int(k))

    yield implicit_index, _set
    # Restore default K=0 between tests so ordering does not bleed state.
    resolver.set_splicing_anchor_tolerance(0)


def _exon(start: int, end: int) -> GenomicInterval:
    return GenomicInterval("chr1", start, end, Strand.POS)


def _resolve(index, exons, introns=()):
    frag = make_fragment(exons=exons, introns=introns)
    return resolve_fragment(frag, index)


# ---------------------------------------------------------------------------
# Acceptance matrix
# ---------------------------------------------------------------------------


class TestImplicitSpliceDiscriminant:
    """14-case acceptance matrix for the per-intron implicit-splice predicate."""

    # (a) intron strictly inside PE gap, K=0 → SPLICED_IMPLICIT
    def test_a_intron_strictly_in_gap_k0(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(0)
        # Fragment: blocks (149,200),(299,350); gap = (200,299) == t_short intron.
        result = _resolve(index, exons=(_exon(149, 200), _exon(299, 350)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

    # (b) PE gap does not overlap any candidate intron → UNSPLICED
    def test_b_intron_not_in_gap(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(0)
        # Two blocks both inside exon 1 of t_short / t_three (PE gap inside an exon).
        # blocks (120,150),(160,190); gap = (150,160) — entirely inside exon 1.
        # No candidate intron overlaps the gap.
        result = _resolve(index, exons=(_exon(120, 150), _exon(160, 190)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (c) intron protrudes 2 bp left of gap, K=3 → SPLICED_IMPLICIT
    def test_c_protrusion_within_K(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # t_short intron is (200,299). gap (202,299): left protrusion = 2.
        # blocks (149,202),(299,350).
        result = _resolve(index, exons=(_exon(149, 202), _exon(299, 350)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

    # (d) intron protrudes 5 bp left of gap, K=3 → UNSPLICED
    def test_d_protrusion_exceeds_K(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # gap (205,299): left protrusion = 5 bp > K=3.
        # blocks (149,205),(299,350).
        result = _resolve(index, exons=(_exon(149, 205), _exon(299, 350)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (e) 200 bp slice of a 50 kb intron lies in gap, K=3 → UNSPLICED
    def test_e_slice_of_long_intron(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # t_long intron is (1100, 51099). Fragment block 1 is in exon 1 (1050,1100),
        # block 2 falls inside the intron at (10000,10100); both intron boundaries
        # lie far outside the (1100,10000) gap.
        # block 2 lands intronic → empty t_set for that block; first block is
        # annotated → not chimeric per "one annotated + one intergenic" rule.
        result = _resolve(index, exons=(_exon(1050, 1100), _exon(10000, 10100)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (f) 4 bp microintron entirely in gap, K=3 → SPLICED_IMPLICIT
    def test_f_microintron_in_gap(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # t_micro micro-intron (800,804). gap (800,804); blocks (750,800),(804,850).
        result = _resolve(index, exons=(_exon(750, 800), _exon(804, 850)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

    # (g) 4 bp microintron entirely inside aligned block → UNSPLICED
    def test_g_microintron_in_block(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # blocks (750,810),(820,850); gap = (810,820). Microintron (800,804) lies
        # entirely in block 1, NOT in the gap.
        result = _resolve(index, exons=(_exon(750, 810), _exon(820, 850)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (h) two-intron transcript, only one intron contained → SPLICED_IMPLICIT
    def test_h_any_intron_satisfies(self, with_tolerance):
        """One of a transcript's introns in the gap is enough; the others need not be anywhere near it.

        Block 1 sits wholly inside exon 2, so intron 1 falls outside the fragment entirely and only
        intron 2 is in play, which is what "any intron satisfies" means. Running block 1 contiguously
        across intron 1 instead would make this a different case — a transcript the read CONTRADICTS,
        which is no longer a candidate for anything — and that case has its own test,
        ``test_h2_a_CONTRADICTED_transcript_is_not_a_candidate``.
        """
        index, set_k = with_tolerance
        set_k(3)
        # t_three has introns (200,299) and (400,499). Block 1 is inside exon 2, block 2 inside exon 3.
        # blocks (320,400),(499,550); gap = (400,499) == intron 2; intron 1 is outside the fragment.
        result = _resolve(index, exons=(_exon(320, 400), _exon(499, 550)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

    # (h2) a transcript the READ CONTRADICTS is not a candidate → UNSPLICED
    def test_h2_a_CONTRADICTED_transcript_is_not_a_candidate(self, with_tolerance):
        """`TRAPS: a-purity-filter-is-a-length-filter`, at the resolver level.

        A candidate is any transcript with a compatible combination of introns and exons, and
        ``cr.t_inds`` is NOT that set: it comes from ``merge_sets``, which falls back to a UNION when the
        intersection is empty, so a transcript the reads contradict can be in it.

        Block 1 runs contiguously from 149 to 400, straight across ``t_three``'s intron (200,299) — the
        read has 99 sequenced bases where ``t_three`` has no exon, so the molecule cannot be
        ``t_three``. Without the compatibility predicate the enumeration would take ``t_three``'s OTHER
        intron (400,499) as a gap hypothesis, taking 99 bp out of ``L`` on the authority of a transcript
        the read disproves, and report the fragment as an implicit splice.

        What is left is the honest answer: no annotated transcript explains this fragment as spliced, so
        the only hypothesis is the unspliced one — the molecule is gDNA, or nascent RNA that retains
        intron 1, and the gap is real template. ``t_three`` is still in ``t_inds``; the predicate is what
        stops it being believed.

        A partially spliced UNANNOTATED path — intron 2 removed, intron 1 retained — is physically real
        and is deliberately not a hypothesis, because unannotated sj are out of scope. That is a scope
        decision rather than an oversight, and it is why the answer here is UNSPLICED rather than a third
        possibility.
        """
        index, set_k = with_tolerance
        set_k(3)
        result = _resolve(index, exons=(_exon(149, 400), _exon(499, 550)))
        assert result is not None
        assert 1 in list(result.t_inds), (
            "t_three must still be in the overlap-derived candidate set, or this test proves nothing "
            "about the compatibility predicate — it would just be a fragment with no candidates"
        )
        assert result.n_gap_hypotheses == 1, (
            "the only hypothesis may be the unspliced one; a second means t_three's intron (400,499) was "
            "taken as a gap path on the authority of a transcript the read contradicts"
        )
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (i) single-block fragment → UNSPLICED (gate exons.size() >= 2 trips)
    def test_i_single_block(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # Single block over exon 1 only.
        result = _resolve(index, exons=(_exon(120, 180),))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (j) CIGAR-spliced annotated fragment → SPLICED_ANNOT (gate not entered)
    def test_j_cigar_spliced_annot_unaffected(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # Annotated SJ for t_short / t_three: (200, 299).
        result = _resolve(
            index,
            exons=(_exon(149, 200), _exon(299, 350)),
            introns=(GenomicInterval("chr1", 200, 299, Strand.POS),),
        )
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT)

    # (k) CIGAR-spliced unannotated fragment → SPLICED_UNANNOT (gate not entered)
    def test_k_cigar_spliced_unannot_unaffected(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # Use an unannotated SJ inside exon 1 region; matches no annotated sj.
        result = _resolve(
            index,
            exons=(_exon(120, 150), _exon(170, 195)),
            introns=(GenomicInterval("chr1", 150, 170, Strand.POS),),
        )
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_UNANNOT)

    # (l) two candidates, only one has contained intron → SPLICED_IMPLICIT
    def test_l_any_candidate_semantics(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # g5 locus: t_one_exon (single span 60000..60500) and t_split with intron
        # (60100, 60399). blocks (60050,60100),(60399,60450); gap = (60100,60399).
        # Only t_split has an intron in the gap.
        result = _resolve(index, exons=(_exon(60050, 60100), _exon(60399, 60450)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

    # (m) nearby but disjoint intron within K bp of gap → UNSPLICED
    def test_m_disjoint_intron_within_K(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # t_short intron (200,299). gap (302,350): intron_end 299 < gap_start 302.
        # No positive overlap (despite 3 bp distance ≤ K). Predicate requires
        # positive overlap before applying slack.
        # blocks (149,302),(350,399). Both blocks are in exon 1+2 footprints.
        result = _resolve(index, exons=(_exon(149, 302), _exon(350, 399)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)

    # (n) chimeric multi-block fragment → not promoted to SPLICED_IMPLICIT
    def test_n_chimera_gate_preserved(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        # Block 1 in t_long exon 1 (1050..1100); block 2 in t_no_intron (53050..53100).
        # Disjoint t_sets → chimera. Gap (1100, 53050) would otherwise contain
        # t_long intron (1100, 51099) but the chimera_type != CHIMERA_NONE
        # gate prevents implicit promotion. Result remains UNSPLICED.
        result = _resolve(index, exons=(_exon(1050, 1100), _exon(53050, 53100)))
        assert result is not None
        assert result.splice_type == int(SpliceType.UNSPLICED)


class TestGapIntronsAreSearchedWhateverTheSpliceType:
    """Gap-intron detection runs on EVERY fragment, not only unspliced ones.

    Looking for an annotated intron in the unsequenced mate gap only once the resolver has called the
    fragment ``SPLICE_UNSPLICED`` leaves that intron inside ``L`` for every fragment that also carries an
    observed CIGAR-N splice, which inflates the length tally's upper tail.

    The condition is on the SPLICE TYPE, never on the block count: an unspliced paired-end fragment
    already has two blocks and a mate gap, so the population at risk is SPLICED fragments that also have
    a gap intron, which are long by construction because they span two or more introns.

    Geometry (g6/t_mixed, 0-based half-open):

        exons     [62000,62200)   [62400,62600)   [62800,63000)
        introns             [62200,62400)   [62600,62800)
        fragment  ---=====|~~~~~~~~~|=====|.........|=====---
                  block1   OBSERVED  block2  mate gap  block3
                           CIGAR-N            holds [62600,62800)
    """

    OBSERVED = (62200, 62400)  #: sequenced as CIGAR-N; the detector must NOT re-derive it
    IMPLIED = (
        62600,
        62800,
    )  #: never sequenced; lies inside the mate gap and must be region_bound from L

    #: block1 · block2 · block3, with [62200,62400) crossed by an observed CIGAR-N splice and
    #: [62500,62900) an unsequenced mate gap.
    MIXED_BLOCKS = ((62100, 62200), (62400, 62500), (62900, 63000))

    #: The same molecule with the observed splice landing 2 bp INSIDE the annotated donor, so the
    #: observed intron is UNANNOTATED and the annotated one sits within the ±K anchor tolerance of it.
    NEAR_MISS_BLOCKS = ((62100, 62202), (62400, 62500), (62900, 63000))
    NEAR_MISS_OBSERVED = (62202, 62400)

    @staticmethod
    def _emitted(result):
        """Every hypothesis as a list of ``(start, end)`` pairs. ``[]`` is the unspliced one."""
        return [[(int(i[1]), int(i[2])) for i in path] for path in result.gap_hypotheses]

    def _resolve_mixed(self, index, blocks, observed):
        return _resolve(
            index,
            exons=tuple(_exon(*b) for b in blocks),
            introns=(GenomicInterval("chr1", observed[0], observed[1], Strand.POS),),
        )

    # U1 — the headline: a SPLICED fragment's gap intron is found.
    def test_U1_a_spliced_fragment_has_the_intron_in_its_mate_gap_found(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.MIXED_BLOCKS, self.OBSERVED)
        assert result is not None
        assert self._emitted(result) == [[self.IMPLIED]], (
            "the annotated intron inside the unsequenced mate gap was not found on a fragment that "
            "also carries an observed CIGAR-N splice — "
        )
        assert result.n_gap_hypotheses == 1

    # U2 — the observed CIGAR-N gap is NOT re-derived.
    def test_U2_the_observed_cigar_n_intron_is_not_re_derived(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.MIXED_BLOCKS, self.OBSERVED)
        assert result is not None
        assert all(self.OBSERVED not in path for path in self._emitted(result)), (
            "the gap finder walks consecutive aligned blocks, so a CIGAR-N intron is also a 'hole'. "
            "It must be dropped by EXACT (start, end) equality against the observed introns"
        )

    # U3 — the near-match trap: a different, nearby annotated intron must not be substituted.
    def test_U3_a_near_match_to_the_observed_gap_is_not_substituted(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.NEAR_MISS_BLOCKS, self.NEAR_MISS_OBSERVED)
        assert result is not None
        # The genuine mate-gap intron is still found — without this the assertion below would pass on a
        # detector that never ran at all.
        #
        # TWO hypotheses, and the empty one is not an accident: this fragment's observed splice is
        # UNANNOTATED, so it does NOT certify the molecule as RNA (an unannotated CIGAR-N may be a
        # misalignment, which is why `FragmentPool.RNA_SPLICED` requires an ANNOTATED sj). The
        # unspliced — genomic — hypothesis therefore stays live and the accumulator will defer.
        assert self._emitted(result) == [[self.IMPLIED], []], (
            "the observed gap [62202,62400) is within the K=3 anchor tolerance of the ANNOTATED intron "
            "[62200,62400), so dropping only exact matches is what stops a DIFFERENT intron being "
            "substituted for it. Two overlapping introns then normalise into one wider one and L comes "
            "out too SHORT — "
        )

    # U4 — the classification does not move. This work is about L, not about labelling.
    def test_U4_splice_type_does_not_move_when_a_gap_intron_is_found(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.MIXED_BLOCKS, self.OBSERVED)
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_ANNOT), (
            "detection is unconditional but the SPLICE_IMPLICIT PROMOTION stays unspliced-only: "
            "splice_type feeds scoring, the buffer, the strand training and the report's census, and "
            "re-labelling would silently move mass between reported categories"
        )
        near = self._resolve_mixed(index, self.NEAR_MISS_BLOCKS, self.NEAR_MISS_OBSERVED)
        assert near is not None
        assert near.splice_type == int(SpliceType.SPLICED_UNANNOT)


# ── The scanner: an implicit splice deposits when its path is determined, defers when it is not ─
#
# A ``SPLICE_IMPLICIT`` fragment has an annotated intron inside its unsequenced mate gap. The splice
# motif was never read, so ``sj_strand`` comes from the transcript that implied it — which is only
# legitimate when the candidates agree about WHICH intron is there. If they do not, ``L`` is
# undetermined and nothing about the fragment can be tallied. Both arms are pinned here on scenarios
# small enough to reason about, because on real cfRNA libraries the rule defers nearly every implicit
# fragment: that is either a true property of real candidate sets or a mistake in the unanimity test,
# and only the small scenarios let the real-data number be read as evidence rather than guessed at.


SEED = 4242

#: Exons far enough apart, and reads short enough, that a fragment spanning the sj puts BOTH mates
#: wholly inside an exon and the intron wholly inside the gap — which is what SPLICE_IMPLICIT means.
SIM = ReadSimConfig(
    frag_mean=260,
    frag_std=20,
    frag_min=220,
    frag_max=320,
    read_length=100,
    strand_specificity=1.0,
    seed=SEED,
)


def _scan(scenario, transcripts):
    scenario.add_gene("g1", "+", transcripts)
    result = scenario.build_oracle(n_fragments=1500, sim_config=SIM)
    _, _, _, payload = scan_and_buffer(
        str(result.bam_path), result.index, BamScanConfig(sj_strand_tag="auto")
    )
    return payload


@pytest.fixture
def scenario(tmp_path):
    sc = Scenario("implicit", genome_length=8000, seed=SEED, work_dir=tmp_path / "implicit")
    yield sc
    sc.cleanup()


def test_ONE_candidate_transcript_DEPOSITS_with_the_strand_inferred_from_it(scenario):
    """The unanimous case. One isoform, so the implied intron is not in doubt and the fragment tallies.

    If this ever goes to zero deposits, the unanimity test has become unsatisfiable and the whole
    implicit population is silently deferred, which on real data looks like a census rather than a bug.
    """
    payload = _scan(
        scenario, [{"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100}]
    )
    # Non-vacuity, and it is read off the UMBRELLA CENSUS rather than off a splice label. The census
    # counts every fragment whose gap needed resolving, which is exactly this scenario's population,
    # while `splice_type` is the scanner's record of what it SAW and is a different axis. If this ever
    # goes to zero the scenario stopped producing gap introns and everything below passes vacuously.
    assert payload.gap_resolution.gap_resolved_spliced > 0, (
        "no fragment had an intron in its unsequenced gap at all — the scenario stopped producing them, "
        "so this test would pass vacuously"
    )
    assert payload.qc.deferred_undetermined_gap == 0, (
        "a single candidate transcript cannot disagree with itself, so nothing here may be deferred"
    )
    assert payload.deferred.n_fragments == 0
    # It deposited: its sj was credited, and it is barred from the pure-RNA length pool.
    assert int(payload.sj_count.sum()) > 0


def test_TWO_candidates_implying_DIFFERENT_introns_are_DEFERRED(scenario):
    """The ambiguous case: two isoforms whose introns start at different places inside the same gap.

    Their implied ``L`` differs by 100 bp, so there is no answer to deposit — picking either is picking a
    fragment length at random.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1100), (3000, 3200)], "abundance": 100},
        ],
    )
    assert payload.qc.deferred_undetermined_gap > 0, (
        "two isoforms imply introns [1200,3000) and [1100,3000) in the same gap — a 100 bp difference in L "
        "— so those fragments have no determined path and must not deposit"
    )
    # And they are HELD, not dropped: the bank carries as many fragments as the counter claims, with
    # both paths on each, which is what makes them recoverable in the second pass.
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap
    assert int(np.diff(payload.deferred.hypothesis_offsets).min()) >= 2


def test_A_RETAINED_INTRON_ISOFORM_ALONE_IS_ENOUGH_TO_DEFER(scenario):
    """The case that makes the rule strict on real data.

    ``t2`` covers the whole locus as one exon, so for it the mate gap holds no intron at all: the
    fragment is unspliced with an ``L`` that includes the gap. That is a different hypothesis from
    ``t1``'s spliced ``L``, so "implies nothing here" has to count as a distinct answer — and GENCODE is
    full of retained-intron isoforms.

    The locus is TIGHT on purpose — a 300 bp intron, so the unspliced hypothesis's ``L`` is about
    560 bp and stays under ``max_fragment_length``. Spread the exons 1800 bp apart instead and the
    unspliced hypothesis is ruled out **by length** and the fragment deposits; that is the next test, and
    keeping the two apart is what stops either passing for the other's reason.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (1500, 1700)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1700)], "abundance": 100},
        ],
    )
    assert payload.qc.deferred_undetermined_gap > 0, (
        "one candidate implies an intron in the gap and the other implies none, which are two different "
        "fragment lengths for one molecule"
    )
    # The subclass names the question, and here it is the composition one: the unspliced path against
    # one spliced path is "RNA or gDNA", one bit. `t2` covering the locus as a single exon is what puts the
    # unspliced path in the set, and GENCODE is full of such retained-intron isoforms.
    assert payload.gap_resolution.gap_deferred_rna_or_gdna > 0
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap
    # Every held record carries BOTH paths, one of them empty — the unspliced hypothesis needs no flag,
    # because cutting nothing IS the statement that the gap is real template.
    empty = [
        h
        for h in range(payload.deferred.n_hypotheses)
        if payload.deferred.hypothesis_introns_of(h).size == 0
    ]
    assert len(empty) == payload.deferred.n_fragments, (
        "exactly one hypothesis per held fragment must be the empty (unspliced) path"
    )


def test_A_SPAN_OVER_THE_LIMIT_RULES_OUT_the_retained_intron_hypothesis(scenario):
    """The other side of the same rule, on a real scan — and it is not a separate rule.

    "If the genomic span exceeds ``max_fragment_length``, assume it is RNA" is the ordinary hypothesis
    filter applied to the unspliced path, whose ``L`` IS that span. Spread ``t1``'s exons 1800 bp apart
    and a sj-spanning fragment's unspliced ``L`` is ~2100 bp against a limit of 1000, so the
    retained-intron explanation is deleted and the spliced path stands alone and deposits.

    The filter is therefore not purely a cost gate: it changes CLASSIFICATION, so the same annotation
    defers in the test above and resolves here, depending only on how far apart the exons sit. Measuring
    that from both sides is the difference between a documented consequence and an assumption.
    """
    payload = _scan(
        scenario,
        [
            {"t_id": "t1", "exons": [(1000, 1200), (3000, 3200)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 3200)], "abundance": 100},
        ],
    )
    assert payload.gap_resolution.gap_resolved_spliced > 0, (
        "no fragment had a gap intron at all, so the span rule is not being exercised — this would pass "
        "vacuously on a scenario that stopped producing mate gaps"
    )
    assert payload.qc.deferred_undetermined_gap == 0, (
        "the unspliced hypothesis's L IS the genomic span, ~2100 bp here against a 1000 bp limit, so it "
        "must be filtered and the fragment's path determined"
    )
    assert payload.deferred.n_fragments == 0
    assert payload.qc.dropped_too_long == 0, (
        "and the SPLICED path is well under the limit, so nothing may be dropped as too long — a fragment "
        "rejected here would mean the filter deleted the wrong hypothesis"
    )
