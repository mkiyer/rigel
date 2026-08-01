"""Tests for the per-intron SPLICED_IMPLICIT discriminant.

Implements the 14-case acceptance matrix from
the implicit-splice bug fix (2026-05-09), section 4.

Strategy
--------
We exercise the C++ ``FragmentResolver`` directly through the standard
Python entry point ``rigel.resolution.resolve_fragment`` and the
``make_fragment`` helper.  A single comprehensive GTF fixture (built
once per session) provides every transcript geometry needed for the
matrix; each test sets the resolver's ``splicing_anchor_tolerance``
explicitly so tests are independent of execution order.

Geometry layout (0-based half-open, after GTF parse)
---------------------------------------------------
  g1 (+):
    t_short   : exons (99,200),(299,400)            intron (200,299)
    t_three   : exons (99,200),(299,400),(499,600)  introns (200,299),(400,499)
  g2 (+):
    t_micro   : exons (699,800),(804,900)           micro-intron (800,804)
  g3 (+):
    t_long    : exons (999,1100),(51099,52000)      long intron (1100,51099)
  g4 (+):
    t_no_intron   : exons (52999,53100)             single-exon
  g5 (+):
    t_one_exon: exons (59999,60500)                 single-exon, spans whole locus
    t_split   : exons (59999,60100),(60399,60500)   intron (60100,60399)
  g6 (+):  ⭐ SPEC_GAP_INTRONS.md — the MIXED fragment: an observed splice AND a gap intron
    t_mixed   : exons (62000,62200),(62400,62600),(62800,63000)
                introns (62200,62400) and (62600,62800)
"""

from __future__ import annotations

import textwrap

import pytest
from conftest import build_test_index

from rigel.types import GenomicInterval, Strand
from rigel.splice import SpliceType
from _resolution_reference import make_fragment, resolve_fragment


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
        index, set_k = with_tolerance
        set_k(3)
        # t_three has introns (200,299) and (400,499). Fragment is contiguous
        # over exons 1+2 in block 1, then jumps via PE gap over intron 2 to
        # exon 3 in block 2.
        # blocks (149,400),(499,550); gap = (400,499) == intron 2.
        result = _resolve(index, exons=(_exon(149, 400), _exon(499, 550)))
        assert result is not None
        assert result.splice_type == int(SpliceType.SPLICED_IMPLICIT)

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
        # Use an unannotated SJ inside exon 1 region; matches no annotated junction.
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
    """⭐ `docs/SPEC_GAP_INTRONS.md` — detection runs on EVERY fragment, not only unspliced ones.

    An annotated intron sitting in a fragment's **unsequenced mate gap** used to be looked for only when
    the resolver had already called the fragment ``SPLICE_UNSPLICED``, so a fragment carrying an observed
    CIGAR-N splice kept that intron inside ``L``. Measured on the chr22 pilot: the library's longest
    molecule is 713 bp and the tally reported 0.97 % of its mass at ≥ 700 bp.

    ⛔ **The gate is on the SPLICE TYPE, never on the block count.** An unspliced paired-end fragment
    already has two blocks and a mate gap and always worked; the missed population is *spliced* fragments
    that ALSO have a gap intron, which are long by construction because they span two or more introns.

    Geometry (g6/t_mixed, 0-based half-open):

        exons     [62000,62200)   [62400,62600)   [62800,63000)
        introns             [62200,62400)   [62600,62800)
        fragment  ---=====|~~~~~~~~~|=====|.........|=====---
                  block1   OBSERVED  block2  mate gap  block3
                           CIGAR-N            holds [62600,62800)
    """

    OBSERVED = (62200, 62400)  #: sequenced as CIGAR-N; the detector must NOT re-derive it
    IMPLIED = (62600, 62800)  #: never sequenced; lies inside the mate gap and must be cut from L

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
            "also carries an observed CIGAR-N splice — SPEC_GAP_INTRONS.md §0"
        )
        assert result.n_gap_hypotheses == 1

    # U2 — ⛔ the observed CIGAR-N gap is NOT re-derived (SPEC §1).
    def test_U2_the_observed_cigar_n_intron_is_not_re_derived(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.MIXED_BLOCKS, self.OBSERVED)
        assert result is not None
        assert all(self.OBSERVED not in path for path in self._emitted(result)), (
            "the gap finder walks consecutive aligned blocks, so a CIGAR-N intron is also a 'hole'. "
            "It must be dropped by EXACT (start, end) equality against the observed introns"
        )

    # U3 — ⛔ THE NEAR-MATCH TRAP. A different, nearby annotated intron must not be substituted.
    def test_U3_a_near_match_to_the_observed_gap_is_not_substituted(self, with_tolerance):
        index, set_k = with_tolerance
        set_k(3)
        result = self._resolve_mixed(index, self.NEAR_MISS_BLOCKS, self.NEAR_MISS_OBSERVED)
        assert result is not None
        # ⭐ The genuine mate-gap intron is still found — without this the assertion below would pass
        # on a detector that never ran at all.
        #
        # ⚠ TWO hypotheses, and the empty one is not an accident: this fragment's observed splice is
        # UNANNOTATED, so it does NOT certify the molecule as RNA (an unannotated CIGAR-N may be a
        # misalignment, which is why `FragmentPool.RNA_SPLICED` requires an ANNOTATED junction). The
        # unspliced — genomic — hypothesis therefore stays live and the accumulator will defer.
        assert self._emitted(result) == [[self.IMPLIED], []], (
            "the observed gap [62202,62400) is within the K=3 anchor tolerance of the ANNOTATED intron "
            "[62200,62400), so dropping only exact matches is what stops a DIFFERENT intron being "
            "substituted for it. Two overlapping introns then normalise into one wider one and L comes "
            "out too SHORT — SPEC_GAP_INTRONS.md §1"
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
