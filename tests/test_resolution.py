"""Tests for hulkrna.resolution — set merging, fragment length, ResolvedFragment, chimera."""

import pytest

from hulkrna.types import ChimeraType, MergeCriteria, Strand, GenomicInterval, IntervalType
from hulkrna.categories import SpliceType
from hulkrna.resolution import (
    merge_sets_with_criteria,
    compute_genomic_frag_length,
    ResolvedFragment,
    _detect_intrachromosomal_chimera,
    resolve_fragment,
    filter_by_overlap,
)
from hulkrna.fragment import Fragment


# =====================================================================
# merge_sets_with_criteria
# =====================================================================


class TestMergeSetsWithCriteria:
    """Progressive relaxation: intersection → nonempty intersection → union."""

    def test_empty_input(self):
        result = merge_sets_with_criteria([])
        assert result.is_empty
        assert result.criteria == MergeCriteria.EMPTY

    def test_single_set(self):
        result = merge_sets_with_criteria(
            [frozenset({1, 2})]
        )
        assert result.t_inds == frozenset({1, 2})
        assert result.criteria == MergeCriteria.INTERSECTION

    def test_intersection_succeeds(self):
        """When all sets share a common element, intersection wins."""
        result = merge_sets_with_criteria(
            [frozenset({1, 2, 3}), frozenset({2, 3, 4})],
        )
        assert result.t_inds == frozenset({2, 3})
        assert result.criteria == MergeCriteria.INTERSECTION

    def test_intersection_nonempty_when_empty_sets_cause_full_intersection_to_fail(self):
        """An empty set in the list collapses intersection to empty.
        Filtering to non-empty sets should recover the result."""
        result = merge_sets_with_criteria(
            [frozenset({1, 2}), frozenset()],
        )
        assert result.t_inds == frozenset({1, 2})
        assert result.criteria == MergeCriteria.INTERSECTION_NONEMPTY

    def test_union_fallback(self):
        """When non-empty sets are disjoint, fall through to union."""
        result = merge_sets_with_criteria(
            [frozenset({1}), frozenset({2})],
        )
        assert result.t_inds == frozenset({1, 2})
        assert result.criteria == MergeCriteria.UNION

    def test_all_empty_sets(self):
        """All empty sets → intersection is empty, nonempty filter finds nothing → union empty."""
        result = merge_sets_with_criteria(
            [frozenset(), frozenset()],
        )
        # All empty: intersection empty, nonempty list empty, union empty
        assert result.t_inds == frozenset()
        assert result.is_empty

    def test_t_only_intersection(self):
        """Result uses INTERSECTION when sets have overlap."""
        result = merge_sets_with_criteria(
            [frozenset({1, 2}), frozenset({2, 3})],
        )
        assert result.t_inds == frozenset({2})
        assert result.criteria == MergeCriteria.INTERSECTION


# =====================================================================
# compute_genomic_frag_length
# =====================================================================


class TestFragmentLength:
    def test_empty_fragment(self):
        frag = Fragment(exons=(), introns=())
        assert compute_genomic_frag_length(frag) == -1

    def test_single_exon(self):
        frag = Fragment(
            exons=(GenomicInterval("chr1", 100, 300, Strand.POS),),
            introns=(),
        )
        assert compute_genomic_frag_length(frag) == 200

    def test_two_exons_no_intron(self):
        """Gap between exons without annotated intron → footprint includes gap."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(),
        )
        # footprint = 400 - 100 = 300, no introns
        assert compute_genomic_frag_length(frag) == 300

    def test_two_exons_with_intron(self):
        """Observed intron subtracted from footprint."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 300, Strand.POS),),
        )
        # footprint = 400 - 100 = 300, intron = 100
        assert compute_genomic_frag_length(frag) == 200

    def test_overlapping_exons(self):
        """Overlapping exon blocks are treated naively (footprint only)."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 250, Strand.POS),
                GenomicInterval("chr1", 200, 400, Strand.POS),
            ),
            introns=(),
        )
        # footprint = 400 - 100 = 300
        assert compute_genomic_frag_length(frag) == 300


# =====================================================================
# ResolvedFragment properties
# =====================================================================


class TestResolvedFragment:
    def _make(self, **kwargs):
        defaults = dict(
            t_inds=frozenset({0}),
            n_genes=1,
            splice_type=SpliceType.UNSPLICED,
            exon_strand=Strand.POS,
            sj_strand=Strand.NONE,
            frag_lengths={0: 250},
            merge_criteria=MergeCriteria.INTERSECTION,
            num_hits=1,
            genomic_footprint=250,
            genomic_start=1000,
        )
        defaults.update(kwargs)
        return ResolvedFragment(**defaults)

    def test_is_unique_gene_single(self):
        rf = self._make(n_genes=1)
        assert rf.is_unique_gene is True

    def test_is_unique_gene_multi(self):
        rf = self._make(n_genes=2)
        assert rf.is_unique_gene is False

    def test_is_ambiguous_multi_gene(self):
        rf = self._make(n_genes=2, num_hits=1)
        assert (rf.n_genes > 1 or rf.num_hits > 1) is True

    def test_is_ambiguous_multimapped(self):
        rf = self._make(n_genes=1, num_hits=3)
        assert (rf.n_genes > 1 or rf.num_hits > 1) is True

    def test_not_ambiguous(self):
        rf = self._make(n_genes=1, num_hits=1)
        assert (rf.n_genes > 1 or rf.num_hits > 1) is False

    def test_has_annotated_sj(self):
        rf = self._make(splice_type=SpliceType.SPLICED_ANNOT)
        assert rf.splice_type == SpliceType.SPLICED_ANNOT

    def test_no_annotated_sj(self):
        for cat in (SpliceType.UNSPLICED, SpliceType.SPLICED_UNANNOT):
            rf = self._make(splice_type=cat)
            assert rf.splice_type != SpliceType.SPLICED_ANNOT

    def test_is_strand_qualified_all_criteria(self):
        rf = self._make(
            splice_type=SpliceType.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is True

    def test_not_strand_qualified_wrong_category(self):
        rf = self._make(
            splice_type=SpliceType.UNSPLICED,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_multi_gene(self):
        rf = self._make(
            splice_type=SpliceType.SPLICED_ANNOT,
            n_genes=2,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_ambiguous_exon(self):
        rf = self._make(
            splice_type=SpliceType.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.AMBIGUOUS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_none_sj(self):
        rf = self._make(
            splice_type=SpliceType.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NONE,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_chimeric(self):
        rf = self._make(
            splice_type=SpliceType.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
            chimera_type=ChimeraType.CIS_STRAND_SAME,
        )
        assert rf.is_strand_qualified is False

    def test_is_chimeric_none(self):
        rf = self._make()
        assert rf.is_chimeric is False

    def test_is_chimeric_trans(self):
        rf = self._make(chimera_type=ChimeraType.TRANS)
        assert rf.is_chimeric is True

    def test_is_chimeric_cis_strand_same(self):
        rf = self._make(chimera_type=ChimeraType.CIS_STRAND_SAME)
        assert rf.is_chimeric is True

    def test_is_chimeric_cis_strand_diff(self):
        rf = self._make(chimera_type=ChimeraType.CIS_STRAND_DIFF)
        assert rf.is_chimeric is True


# =====================================================================
# _detect_intrachromosomal_chimera
# =====================================================================


class TestDetectIntrachromosomalChimera:
    """Transcript-set disjointness detection for intrachromosomal chimeras."""

    def test_single_block_none(self):
        """One annotated block cannot be chimeric."""
        blocks = (GenomicInterval("chr1", 100, 200, Strand.POS),)
        t_sets = [frozenset({0, 1})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_empty_t_sets_none(self):
        """All empty transcript sets → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset(), frozenset()]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_one_empty_one_annotated(self):
        """One annotated, one intergenic → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset()]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_shared_transcript_none(self):
        """Two blocks sharing a transcript → connected → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({1, 2})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_disjoint_same_strand_chimera(self):
        """Two blocks with disjoint transcript sets, same strand → CIS_STRAND_SAME."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({2, 3})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        assert chimera_gap == 800  # 1000 - 200

    def test_disjoint_diff_strand_chimera(self):
        """Two blocks with disjoint sets, different strands → CIS_STRAND_DIFF."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.NEG),
        )
        t_sets = [frozenset({0}), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 300  # 500 - 200

    def test_overlapping_disjoint_blocks(self):
        """Disjoint transcript sets with overlapping genomic coordinates."""
        blocks = (
            GenomicInterval("chr1", 100, 500, Strand.POS),
            GenomicInterval("chr1", 300, 700, Strand.NEG),
        )
        t_sets = [frozenset({0}), frozenset({1})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 0  # overlapping

    def test_three_blocks_two_components(self):
        """Three blocks: two connected (shared tx), one disjoint."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
            GenomicInterval("chr1", 2000, 2100, Strand.NEG),
        )
        # Block 0 and 1 share transcript 1 → connected
        # Block 2 is disjoint from both → chimeric
        t_sets = [frozenset({0, 1}), frozenset({1, 2}), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        # Gap: min(2000-200, 2000-400) = 1600
        assert chimera_gap == 1600

    def test_three_disjoint_blocks(self):
        """Three blocks, all pairwise disjoint → three components → chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset({1}), frozenset({2})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        # Gap: min(500-200, 1000-200, 1000-600) = 300
        assert chimera_gap == 300

    def test_connected_via_transitive_overlap(self):
        """A-B overlap, B-C overlap, but A-C are disjoint → all connected."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.POS),
        )
        # A={0,1}, B={1,2}, C={2,3} → A-B share 1, B-C share 2 → all connected
        t_sets = [frozenset({0, 1}), frozenset({1, 2}), frozenset({2, 3})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_adjacent_blocks_zero_gap(self):
        """Disjoint transcript sets with blocks that are adjacent (gap=0)."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 200, 300, Strand.POS),
        )
        t_sets = [frozenset({0}), frozenset({1})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_SAME
        assert chimera_gap == 0  # adjacent, no gap

    def test_mixed_annotated_intergenic_blocks(self):
        """Annotated blocks are disjoint, intergenic block between them."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.NONE),  # intergenic
            GenomicInterval("chr1", 1000, 1100, Strand.NEG),
        )
        # Block 1 is intergenic (empty set), blocks 0 and 2 are disjoint
        t_sets = [frozenset({0}), frozenset(), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.CIS_STRAND_DIFF
        assert chimera_gap == 800  # 1000 - 200


# =====================================================================
# resolve_fragment
# =====================================================================


class _MockResolveIndex:
    def __init__(self):
        self.sj_map = {
            ("chr1", 200, 300, Strand.POS): frozenset({0})
        }
        self.t_to_g_arr = [0, 0, 0]

    def query(self, exon_block):
        key = (exon_block.ref, exon_block.start, exon_block.end)
        if key == ("chr1", 100, 200):
            return [
                (100, 200, IntervalType.EXON, frozenset({0, 1})),
                (100, 200, IntervalType.TRANSCRIPT, frozenset({0, 1})),
            ]
        if key == ("chr1", 300, 400):
            return [
                (300, 400, IntervalType.EXON, frozenset({1, 2})),
                (300, 400, IntervalType.TRANSCRIPT, frozenset({1, 2})),
            ]
        return []

    def query_gap_sjs(self, _ref, _start, _end):
        return []


class TestResolveFragment:
    def test_annotated_sj_prioritizes_sj_supported_transcripts(self):
        idx = _MockResolveIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 300, Strand.POS),),
        )

        result = resolve_fragment(frag, idx)

        assert result is not None
        assert result.splice_type == SpliceType.SPLICED_ANNOT
        assert result.t_inds == frozenset({0})

    def test_annotated_sj_matches_when_intron_strand_unknown(self):
        idx = _MockResolveIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 300, Strand.NONE),),
        )

        result = resolve_fragment(frag, idx)

        assert result is not None
        assert result.splice_type == SpliceType.SPLICED_ANNOT
        assert result.t_inds == frozenset({0})


# =====================================================================
# compute_overlap_profile
# =====================================================================


class _OverlapMockIndex:
    """Mock index with configurable exon and transcript-span intervals.

    Models a scenario with 3 transcripts:
      t0: exon 100-200, transcript 100-400 (span), exon 350-400
      t1: exon 100-200   (partial overlap with short fragments)
      t2: exon 350-500   (partial overlap with fragments near right end)

    All belong to gene 0.
    """

    def __init__(self):
        self.sj_map = {}
        self.t_to_g_arr = [0, 0, 0]
        # Collapsed intervals: (start, end, itype, t_set)
        self._intervals = [
            (100, 200, IntervalType.EXON, frozenset({0, 1})),
            (100, 400, IntervalType.TRANSCRIPT, frozenset({0})),
            (100, 200, IntervalType.TRANSCRIPT, frozenset({1})),
            (350, 400, IntervalType.EXON, frozenset({0})),
            (350, 500, IntervalType.EXON, frozenset({2})),
            (350, 500, IntervalType.TRANSCRIPT, frozenset({2})),
        ]

    def query(self, exon_block):
        results = []
        for h_start, h_end, itype, t_set in self._intervals:
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((h_start, h_end, itype, t_set))
        return results

    def query_gap_sjs(self, _ref, _start, _end):
        return []


class TestOverlapProfileViaResolve:
    """Test per-candidate overlap profiles produced by resolve_fragment."""

    def test_full_exon_overlap_single_block(self):
        """Fragment fully inside a reference exon → exon_bp=read_length, intron_bp=0."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 150, 200, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 50
        # t0 exon 100-200: covers 150-200 → 50bp exon
        assert result.overlap_bp[0] == (50, 0, 0)
        # t1 exon 100-200: covers 150-200 → 50bp exon
        assert result.overlap_bp[1] == (50, 0, 0)

    def test_intronic_overlap(self):
        """Fragment in intronic region → exon_bp=0, intron_bp=read_length."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 220, 320, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        # Only t0 has transcript span 100-400 covering 220-320 → intron_bp
        # transcript_bp = 100, exon_bp = 0 → intron_bp = 100
        assert result.overlap_bp[0] == (0, 100, 100)

    def test_mixed_exon_intron_overlap(self):
        """Fragment spanning exon-transcript boundary → partial exon + intron bp."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 180, 220, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 40
        # t0: exon 100-200 → 180-200 = 20bp; transcript 100-400 → 180-220 = 40bp
        # intron_bp = 40 - 20 = 20
        assert result.overlap_bp[0] == (20, 20, 20)
        # t1: exon 100-200 → 180-200 = 20bp; transcript 100-200 → 180-200 = 20bp
        # intron_bp = 20 - 20 = 0
        assert result.overlap_bp[1] == (20, 0, 0)

    def test_multi_block_overlap_sums(self):
        """Two exon blocks accumulate overlap per-transcript."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 350, 400, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 350, Strand.POS),),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        # t0: exon(150-200)=50bp + exon(350-400)=50bp = 100bp exon
        # transcript(150-200)=50bp + transcript(350-400)=50bp = 100bp tx
        # intron_bp = 100 - 100 = 0
        assert result.overlap_bp[0] == (100, 0, 0)

    def test_unambig_intron_subtracts_cross_transcript_exons(self):
        """Fragment in t0's intronic region that overlaps t1's exon → unambig_intron_bp=0.

        Models the core nRNA phantom scenario:
        - t0: exon 100-200, transcript 100-500, exon 400-500  (2-exon)
        - t1: exon 100-500  (single-exon, same locus)
        Fragment: 250-350 (block in t0's intronic region, but also in t1's exon)
        → t0 should have intron_bp=100, unambig_intron_bp=0
        → t1 should have exon_bp=100, unambig_intron_bp=0
        """
        class _CrossTxIndex:
            sj_map = {}
            t_to_g_arr = [0, 0]

            def query(self, exon_block):
                _intervals = [
                    (100, 200, IntervalType.EXON, frozenset({0})),
                    (100, 500, IntervalType.TRANSCRIPT, frozenset({0})),
                    (400, 500, IntervalType.EXON, frozenset({0})),
                    (100, 500, IntervalType.EXON, frozenset({1})),
                    (100, 500, IntervalType.TRANSCRIPT, frozenset({1})),
                ]
                results = []
                for h_start, h_end, itype, t_set in _intervals:
                    if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                        results.append((h_start, h_end, itype, t_set))
                return results

            def query_gap_sjs(self, _ref, _start, _end):
                return []

        idx = _CrossTxIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 250, 350, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        # t0: transcript 100-500 covers 250-350 → tx_bp=100, exon_bp=0
        # intron_bp = 100 - 0 = 100
        # But t1's exon 100-500 covers 250-350 → global exon union = [(250,350)]
        # Subtracting global exons from t0's tx (250,350) → 0 unambig
        assert result.overlap_bp[0] == (0, 100, 0)
        # t1: exon 100-500 covers 250-350 → 100bp exon, tx_bp=100
        # intron_bp = 100 - 100 = 0
        assert result.overlap_bp[1] == (100, 0, 0)

    def test_unambig_intron_partial_exon_overlap(self):
        """Partial overlap: some intronic bp are exonic for another transcript.

        - t0: exon 100-200, transcript 100-500, exon 400-500
        - t1: exon 280-350, transcript 280-350
        Fragment: 250-350
        → t0 intron_bp=100 (tx 100 - exon 0), but only 250-280 (30bp) is unambiguous
          (280-350 is exonic for t1)
        """
        class _PartialOverlapIndex:
            sj_map = {}
            t_to_g_arr = [0, 0]

            def query(self, exon_block):
                _intervals = [
                    (100, 200, IntervalType.EXON, frozenset({0})),
                    (100, 500, IntervalType.TRANSCRIPT, frozenset({0})),
                    (400, 500, IntervalType.EXON, frozenset({0})),
                    (280, 350, IntervalType.EXON, frozenset({1})),
                    (280, 350, IntervalType.TRANSCRIPT, frozenset({1})),
                ]
                results = []
                for h_start, h_end, itype, t_set in _intervals:
                    if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                        results.append((h_start, h_end, itype, t_set))
                return results

            def query_gap_sjs(self, _ref, _start, _end):
                return []

        idx = _PartialOverlapIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 250, 350, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.read_length == 100
        # t0: transcript 100-500 clipped to 250-350 → 100bp tx, exon_bp=0
        # intron_bp = 100 - 0 = 100
        # global exon union: t1 exon 280-350 clipped to 280-350 → [(280,350)]
        # unambig: (250,350) - [(280,350)] → 250-280 = 30bp
        assert result.overlap_bp[0] == (0, 100, 30)
        # t1: exon 280-350 clipped to 280-350 → 70bp exon, tx_bp=70
        # intron_bp = 70 - 70 = 0
        assert result.overlap_bp[1] == (70, 0, 0)


# =====================================================================
# filter_by_overlap
# =====================================================================


class TestFilterByOverlap:
    def test_filters_low_overlap_candidates(self):
        t_inds = frozenset({0, 1, 2})
        # BP counts with read_length=100: exon fracs are 100/100, 50/100, 30/100
        profiles = {0: (100, 0, 0), 1: (50, 30, 0), 2: (30, 60, 0)}
        # With min_frac_of_best=0.9 → threshold = 0.9 on exon_frac
        # Only t0 survives (1.0 >= 0.9)
        result = filter_by_overlap(t_inds, profiles, 100)
        assert result == frozenset({0})

    def test_keeps_close_candidates(self):
        t_inds = frozenset({0, 1, 2})
        profiles = {0: (100, 0, 0), 1: (95, 0, 0), 2: (30, 0, 0)}
        result = filter_by_overlap(t_inds, profiles, 100)
        assert result == frozenset({0, 1})

    def test_custom_threshold(self):
        t_inds = frozenset({0, 1, 2})
        profiles = {0: (100, 0, 0), 1: (60, 20, 0), 2: (30, 20, 0)}
        result = filter_by_overlap(t_inds, profiles, 100, min_frac_of_best=0.5)
        assert result == frozenset({0, 1})

    def test_no_fracs_returns_original(self):
        t_inds = frozenset({0, 1})
        result = filter_by_overlap(t_inds, {}, 100)
        assert result == t_inds

    def test_never_returns_empty(self):
        """Even with extreme threshold, should return original rather than empty."""
        t_inds = frozenset({0, 1})
        profiles = {0: (1, 50, 0), 1: (1, 50, 0)}
        # All candidates are equal and above zero → they should all survive
        result = filter_by_overlap(t_inds, profiles, 100, min_frac_of_best=2.0)
        # threshold = 0.01 * 2.0 = 0.02, both are below → would be empty
        # safety: returns original
        assert result == t_inds

    def test_single_candidate_unchanged(self):
        t_inds = frozenset({5})
        profiles = {5: (30, 50, 0)}
        result = filter_by_overlap(t_inds, profiles, 100)
        assert result == frozenset({5})


# =====================================================================
# Overlap filtering integration in resolve_fragment
# =====================================================================


class TestOverlapFiltering:
    """Verify overlap filtering is skipped for UNSPLICED fragments."""

    def test_unspliced_fragment_skips_filter(self):
        """An unspliced fragment spanning an exon-intron boundary.

        Fragment 150-250 overlaps:
          - t0: exon(100-200)=50bp + intron(200-350)=50bp
          - t1: exon(100-200)=50bp

        filter_by_overlap is skipped for UNSPLICED fragments, so both
        candidates are retained regardless of exon overlap differences.
        The EM's soft overlap penalty handles weighting instead.
        """
        idx = _OverlapMockIndex()
        # Single-exon frag 150-250: t0 covers fully, t1 only exon portion
        frag = Fragment(
            exons=(GenomicInterval("chr1", 150, 250, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.UNSPLICED
        # Both t0 and t1 have 50bp exon overlap → equal exon_frac → both kept
        assert 0 in result.t_inds
        assert 1 in result.t_inds

    def test_spliced_fragment_sj_plus_overlap(self):
        """Spliced fragments with annotated SJs: SJ narrows first, overlap applies too."""
        idx = _MockResolveIndex()
        # Standard spliced fragment — splice junction resolves to t0 only
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 300, Strand.POS),),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.SPLICED_ANNOT
        # SJ resolution gives {0} — single candidate, no overlap filtering needed
        assert result.t_inds == frozenset({0})

    def test_equal_overlap_keeps_all(self):
        """When all candidates have equal overlap, all are retained."""
        idx = _OverlapMockIndex()
        # Fragment 100-200 overlaps t0 (100/100) and t1 (100/100)
        frag = Fragment(
            exons=(GenomicInterval("chr1", 100, 200, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.UNSPLICED
        # Both t0 and t1 have equal overlap fractions → both kept
        assert 0 in result.t_inds
        assert 1 in result.t_inds


# =====================================================================
# Fragment-length discrimination (3-exon gene, 2 isoforms)
# =====================================================================


class _ThreeExonMockIndex:
    """Mock index for fragment-length discrimination test.

    Models the MINI_GTF 3-exon/2-isoform gene:
      t0: exons at 99-200, 299-400, 499-600  (3 exons)
      t1: exons at 99-200, 499-600           (2 exons, skips middle)
    Both belong to gene 0.

    Splice junctions:
      t0: (chr1, 200, 299, POS) and (chr1, 400, 499, POS)
      t1: (chr1, 200, 499, POS)
    """

    def __init__(self):
        # SJ map: key=(ref, start, end, strand) → frozenset of t_indices
        self.sj_map = {
            ("chr1", 200, 299, Strand.POS): frozenset({0}),
            ("chr1", 400, 499, Strand.POS): frozenset({0}),
            ("chr1", 200, 499, Strand.POS): frozenset({1}),
        }
        self.t_to_g_arr = [0, 0]
        # Collapsed intervals: (h_start, h_end, itype, t_set)
        self._intervals = [
            # Exon intervals
            (99, 200, IntervalType.EXON, frozenset({0, 1})),
            (299, 400, IntervalType.EXON, frozenset({0})),
            (499, 600, IntervalType.EXON, frozenset({0, 1})),
            # Transcript spans
            (99, 600, IntervalType.TRANSCRIPT, frozenset({0})),
            (99, 600, IntervalType.TRANSCRIPT, frozenset({1})),
        ]
        # SJ annotations for gap queries: (t_idx, strand, sj_start, sj_end)
        self._sj_entries = [
            (0, Strand.POS, 200, 299),
            (0, Strand.POS, 400, 499),
            (1, Strand.POS, 200, 499),
        ]

    def query(self, exon_block):
        results = []
        for h_start, h_end, itype, t_set in self._intervals:
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((h_start, h_end, itype, t_set))
        return results

    def query_gap_sjs(self, ref, start, end):
        """Return SJ entries fully contained in the gap."""
        results = []
        for t_idx, strand, sj_start, sj_end in self._sj_entries:
            if ref == "chr1" and sj_start >= start and sj_end <= end:
                results.append((t_idx, strand, sj_start, sj_end))
        return results


class TestFragLengthDiscrimination:
    """Verify fragment length helps distinguish isoforms with different splicing.

    Scenario: 3-exon gene with t0 (all 3) and t1 (exons 1+3, skipping #2).
    A fragment from t1 that spans exon1→exon3 will have a different
    genomic distance than the same position on t0.
    """

    def test_spliced_read_t0_specific_sj(self):
        """A read with t0-specific SJ (200-299) resolves to t0 only."""
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 299, Strand.POS),),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.SPLICED_ANNOT
        assert result.t_inds == frozenset({0})

    def test_spliced_read_t1_specific_sj(self):
        """A read with t1-specific SJ (200-499 skip) resolves to t1 only."""
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 499, 550, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 499, Strand.POS),),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.SPLICED_ANNOT
        assert result.t_inds == frozenset({1})

    def test_unspliced_read_in_shared_exon(self):
        """Unspliced read fully in exon1 (shared) → both t0 and t1."""
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.UNSPLICED
        # Both transcripts share exon1 fully → equal overlap → both kept
        assert result.t_inds == frozenset({0, 1})

    def test_unspliced_read_in_t0_only_exon(self):
        """Unspliced read in exon2 (t0-only at 299-400).

        With collapsed index, t1's TRANSCRIPT span (99-600) also covers
        this region, giving t1 an intronic overlap.  Both transcripts are
        candidates; downstream scoring discriminates via overlap profiles.
        """
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 320, 380, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.splice_type == SpliceType.UNSPLICED
        # t0: exon overlap, t1: intronic overlap via TRANSCRIPT span
        assert result.t_inds == frozenset({0, 1})

    def test_frag_length_differs_between_isoforms(self):
        """Fragment spanning exon1→exon3 has different fragment length per isoform.

        For t0 (3 exons): reads at 150-200, 499-550 with gap 200-499
          → gap = 299bp, gap contains introns 200-299 (99bp) and 400-499 (99bp)
          → insert = (200-150) + (550-499) + (299 - 99 - 99) = 50 + 51 + 101 = 202
          Actually, insert = genomic_span - sum(gap_introns) = 400 - 198 = 202

        For t1 (2 exons):  gap contains intron 200-499 (299bp)
          → insert = 400 - 299 = 101

        These differ → compute_frag_lengths returns differing values when they disagree.
        """
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 499, 550, Strand.POS),
            ),
            introns=(GenomicInterval("chr1", 200, 499, Strand.POS),),
        )
        # This has annotated SJ 200-499 → resolves to t1 specifically
        result = resolve_fragment(frag, idx)
        assert result is not None
        # SJ 200-499 matches t1 specifically
        assert result.t_inds == frozenset({1})
