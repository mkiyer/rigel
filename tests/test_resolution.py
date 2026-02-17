"""Tests for hulkrna.resolution — set merging, insert size, ResolvedFragment, chimera."""

import pytest

from hulkrna.types import ChimeraType, MergeCriteria, MergeResult, Strand, GenomicInterval, IntervalType
from hulkrna.categories import CountCategory
from hulkrna.resolution import (
    merge_sets_with_criteria,
    fragment_insert_size,
    ResolvedFragment,
    _detect_intrachromosomal_chimera,
    resolve_fragment,
    compute_exon_overlap,
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
# fragment_insert_size
# =====================================================================


class TestFragmentInsertSize:
    def test_empty_fragment(self):
        frag = Fragment(exons=(), introns=())
        assert fragment_insert_size(frag) == -1

    def test_single_exon(self):
        frag = Fragment(
            exons=(GenomicInterval("chr1", 100, 300, Strand.POS),),
            introns=(),
        )
        assert fragment_insert_size(frag) == 200

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
        assert fragment_insert_size(frag) == 300

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
        assert fragment_insert_size(frag) == 200

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
        assert fragment_insert_size(frag) == 300


# =====================================================================
# ResolvedFragment properties
# =====================================================================


class TestResolvedFragment:
    def _make(self, **kwargs):
        defaults = dict(
            t_inds=frozenset({0}),
            n_genes=1,
            count_cat=CountCategory.UNSPLICED,
            exon_strand=Strand.POS,
            sj_strand=Strand.NONE,
            insert_size=250,
            merge_criteria=MergeCriteria.INTERSECTION,
            num_hits=1,
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
        rf = self._make(count_cat=CountCategory.SPLICED_ANNOT)
        assert rf.count_cat == CountCategory.SPLICED_ANNOT

    def test_no_annotated_sj(self):
        for cat in (CountCategory.INTRON, CountCategory.UNSPLICED, CountCategory.SPLICED_UNANNOT):
            rf = self._make(count_cat=cat)
            assert rf.count_cat != CountCategory.SPLICED_ANNOT

    def test_is_strand_qualified_all_criteria(self):
        rf = self._make(
            count_cat=CountCategory.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is True

    def test_not_strand_qualified_wrong_category(self):
        rf = self._make(
            count_cat=CountCategory.UNSPLICED,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_multi_gene(self):
        rf = self._make(
            count_cat=CountCategory.SPLICED_ANNOT,
            n_genes=2,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_ambiguous_exon(self):
        rf = self._make(
            count_cat=CountCategory.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.AMBIGUOUS,
            sj_strand=Strand.NEG,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_none_sj(self):
        rf = self._make(
            count_cat=CountCategory.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NONE,
        )
        assert rf.is_strand_qualified is False

    def test_not_strand_qualified_chimeric(self):
        rf = self._make(
            count_cat=CountCategory.SPLICED_ANNOT,
            n_genes=1,
            exon_strand=Strand.POS,
            sj_strand=Strand.NEG,
            chimera_type=ChimeraType.STRAND_SAME,
        )
        assert rf.is_strand_qualified is False

    def test_is_chimeric_not_chimeric(self):
        rf = self._make()
        assert rf.is_chimeric is False

    def test_is_chimeric_interchromosomal(self):
        rf = self._make(chimera_type=ChimeraType.INTERCHROMOSOMAL)
        assert rf.is_chimeric is True

    def test_is_chimeric_strand_same(self):
        rf = self._make(chimera_type=ChimeraType.STRAND_SAME)
        assert rf.is_chimeric is True

    def test_is_chimeric_strand_diff(self):
        rf = self._make(chimera_type=ChimeraType.STRAND_DIFF)
        assert rf.is_chimeric is True


# =====================================================================
# _detect_intrachromosomal_chimera
# =====================================================================


class TestDetectIntrachromosomalChimera:
    """Transcript-set disjointness detection for intrachromosomal chimeras."""

    def test_single_block_not_chimeric(self):
        """One annotated block cannot be chimeric."""
        blocks = (GenomicInterval("chr1", 100, 200, Strand.POS),)
        t_sets = [frozenset({0, 1})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_empty_t_sets_not_chimeric(self):
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

    def test_shared_transcript_not_chimeric(self):
        """Two blocks sharing a transcript → connected → not chimeric."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 300, 400, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({1, 2})]
        assert _detect_intrachromosomal_chimera(blocks, t_sets) is None

    def test_disjoint_same_strand_chimera(self):
        """Two blocks with disjoint transcript sets, same strand → STRAND_SAME."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 1000, 1100, Strand.POS),
        )
        t_sets = [frozenset({0, 1}), frozenset({2, 3})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.STRAND_SAME
        assert chimera_gap == 800  # 1000 - 200

    def test_disjoint_diff_strand_chimera(self):
        """Two blocks with disjoint sets, different strands → STRAND_DIFF."""
        blocks = (
            GenomicInterval("chr1", 100, 200, Strand.POS),
            GenomicInterval("chr1", 500, 600, Strand.NEG),
        )
        t_sets = [frozenset({0}), frozenset({5})]
        result = _detect_intrachromosomal_chimera(blocks, t_sets)
        assert result is not None
        chimera_type, chimera_gap = result
        assert chimera_type == ChimeraType.STRAND_DIFF
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
        assert chimera_type == ChimeraType.STRAND_DIFF
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
        assert chimera_type == ChimeraType.STRAND_DIFF
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
        assert chimera_type == ChimeraType.STRAND_SAME
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
        assert chimera_type == ChimeraType.STRAND_SAME
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
        assert chimera_type == ChimeraType.STRAND_DIFF
        assert chimera_gap == 800  # 1000 - 200


# =====================================================================
# resolve_fragment
# =====================================================================


class _MockResolveIndex:
    def __init__(self):
        self.sj_map = {
            ("chr1", 200, 300, Strand.POS): (frozenset({0}), frozenset({0}))
        }
        self.t_to_g_arr = [0, 0, 0]

    def query_exon(self, exon_block):
        key = (exon_block.ref, exon_block.start, exon_block.end)
        if key == ("chr1", 100, 200):
            return [
                (0, 0, IntervalType.EXON),
                (1, 0, IntervalType.EXON),
            ]
        if key == ("chr1", 300, 400):
            return [
                (1, 0, IntervalType.EXON),
                (2, 0, IntervalType.EXON),
            ]
        return []

    def query_exon_with_coords(self, exon_block):
        key = (exon_block.ref, exon_block.start, exon_block.end)
        if key == ("chr1", 100, 200):
            return [
                (0, 0, IntervalType.EXON, 100, 200),
                (1, 0, IntervalType.EXON, 100, 200),
            ]
        if key == ("chr1", 300, 400):
            return [
                (1, 0, IntervalType.EXON, 300, 400),
                (2, 0, IntervalType.EXON, 300, 400),
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
        assert result.count_cat == CountCategory.SPLICED_ANNOT
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
        assert result.count_cat == CountCategory.SPLICED_ANNOT
        assert result.t_inds == frozenset({0})


# =====================================================================
# compute_exon_overlap
# =====================================================================


class _OverlapMockIndex:
    """Mock index with configurable exon intervals for overlap tests.

    Models a scenario with 3 transcripts:
      t0: exon 100-400 (one big exon)
      t1: exon 100-200   (partial overlap with short fragments)
      t2: exon 350-500   (partial overlap with fragments near right end)

    All belong to gene 0.
    """

    def __init__(self):
        self.sj_map = {}
        self.t_to_g_arr = [0, 0, 0]
        # Reference intervals: list of (t_idx, g_idx, itype, start, end)
        self._intervals = [
            (0, 0, IntervalType.EXON, 100, 400),
            (1, 0, IntervalType.EXON, 100, 200),
            (2, 0, IntervalType.EXON, 350, 500),
        ]

    def query_exon(self, exon_block):
        results = []
        for t_idx, g_idx, itype, h_start, h_end in self._intervals:
            # Check overlap
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((t_idx, g_idx, itype))
        return results

    def query_exon_with_coords(self, exon_block):
        results = []
        for t_idx, g_idx, itype, h_start, h_end in self._intervals:
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((t_idx, g_idx, itype, h_start, h_end))
        return results

    def query_gap_sjs(self, _ref, _start, _end):
        return []


class TestComputeExonOverlap:
    """Test per-candidate exon overlap fraction computation."""

    def test_full_overlap_single_block(self):
        """Fragment fully inside a reference exon → overlap = 1.0."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 150, 250, Strand.POS),),
            introns=(),
        )
        fracs = compute_exon_overlap(frag, idx)
        # t0 (100-400) fully covers 150-250 → 100/100 = 1.0
        assert fracs[0] == pytest.approx(1.0)
        # t1 (100-200) covers 150-200 → 50/100 = 0.5
        assert fracs[1] == pytest.approx(0.5)
        # t2 (350-500) no overlap with 150-250
        assert 2 not in fracs

    def test_partial_overlap_at_edge(self):
        """Fragment overlapping an exon edge → fractional overlap."""
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 190, 210, Strand.POS),),
            introns=(),
        )
        fracs = compute_exon_overlap(frag, idx)
        # t0 (100-400) fully covers 190-210 → 20/20 = 1.0
        assert fracs[0] == pytest.approx(1.0)
        # t1 (100-200) covers 190-200 → 10/20 = 0.5
        assert fracs[1] == pytest.approx(0.5)

    def test_empty_fragment(self):
        idx = _OverlapMockIndex()
        frag = Fragment(exons=(), introns=())
        fracs = compute_exon_overlap(frag, idx)
        assert fracs == {}

    def test_no_overlapping_intervals(self):
        idx = _OverlapMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 700, 800, Strand.POS),),
            introns=(),
        )
        fracs = compute_exon_overlap(frag, idx)
        assert fracs == {}

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
        fracs = compute_exon_overlap(frag, idx)
        # t0 (100-400): block1 150-200 → 50bp, block2 350-400 → 50bp → 100/100 = 1.0
        assert fracs[0] == pytest.approx(1.0)
        # t1 (100-200): block1 150-200 → 50bp, block2 no overlap → 50/100 = 0.5
        assert fracs[1] == pytest.approx(0.5)
        # t2 (350-500): block1 no overlap, block2 350-400 → 50bp → 50/100 = 0.5
        assert fracs[2] == pytest.approx(0.5)


# =====================================================================
# filter_by_overlap
# =====================================================================


class TestFilterByOverlap:
    def test_filters_low_overlap_candidates(self):
        t_inds = frozenset({0, 1, 2})
        fracs = {0: 1.0, 1: 0.5, 2: 0.3}
        # With min_frac_of_best=0.9 → threshold = 0.9
        # Only t0 survives (1.0 >= 0.9)
        result = filter_by_overlap(t_inds, fracs)
        assert result == frozenset({0})

    def test_keeps_close_candidates(self):
        t_inds = frozenset({0, 1, 2})
        fracs = {0: 1.0, 1: 0.95, 2: 0.3}
        result = filter_by_overlap(t_inds, fracs)
        assert result == frozenset({0, 1})

    def test_custom_threshold(self):
        t_inds = frozenset({0, 1, 2})
        fracs = {0: 1.0, 1: 0.6, 2: 0.3}
        result = filter_by_overlap(t_inds, fracs, min_frac_of_best=0.5)
        assert result == frozenset({0, 1})

    def test_no_fracs_returns_original(self):
        t_inds = frozenset({0, 1})
        result = filter_by_overlap(t_inds, {})
        assert result == t_inds

    def test_never_returns_empty(self):
        """Even with extreme threshold, should return original rather than empty."""
        t_inds = frozenset({0, 1})
        fracs = {0: 0.01, 1: 0.01}
        # All candidates are equal and above zero → they should all survive
        result = filter_by_overlap(t_inds, fracs, min_frac_of_best=2.0)
        # threshold = 0.01 * 2.0 = 0.02, both are below → would be empty
        # safety: returns original
        assert result == t_inds

    def test_single_candidate_unchanged(self):
        t_inds = frozenset({5})
        fracs = {5: 0.3}
        result = filter_by_overlap(t_inds, fracs)
        assert result == frozenset({5})


# =====================================================================
# Overlap filtering integration in resolve_fragment
# =====================================================================


class TestOverlapFiltering:
    """Verify overlap filtering is applied for UNSPLICED fragments only."""

    def test_unspliced_fragment_filtered_by_overlap(self):
        """An unspliced fragment with unequal overlap should be filtered."""
        idx = _OverlapMockIndex()
        # Single-exon frag 150-250: t0 covers fully, t1 only half
        frag = Fragment(
            exons=(GenomicInterval("chr1", 150, 250, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.count_cat == CountCategory.UNSPLICED
        # t1 has only 50% overlap, well below 90% of best (t0 at 100%)
        # → t1 should be filtered out
        assert 0 in result.t_inds
        assert 1 not in result.t_inds

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
        assert result.count_cat == CountCategory.SPLICED_ANNOT
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
        assert result.count_cat == CountCategory.UNSPLICED
        # Both t0 and t1 have equal overlap fractions → both kept
        assert 0 in result.t_inds
        assert 1 in result.t_inds


# =====================================================================
# Insert-size discrimination (3-exon gene, 2 isoforms)
# =====================================================================


class _ThreeExonMockIndex:
    """Mock index for insert-size discrimination test.

    Models the MINI_GTF 3-exon/2-isoform gene:
      t0: exons at 99-200, 299-400, 499-600  (3 exons)
      t1: exons at 99-200, 499-600           (2 exons, skips middle)
    Both belong to gene 0.

    Splice junctions:
      t0: (chr1, 200, 299, POS) and (chr1, 400, 499, POS)
      t1: (chr1, 200, 499, POS)
    """

    def __init__(self):
        # SJ map: key=(ref, start, end, strand) → (t_set, g_set)
        self.sj_map = {
            ("chr1", 200, 299, Strand.POS): (frozenset({0}), frozenset({0})),
            ("chr1", 400, 499, Strand.POS): (frozenset({0}), frozenset({0})),
            ("chr1", 200, 499, Strand.POS): (frozenset({1}), frozenset({0})),
        }
        self.t_to_g_arr = [0, 0]
        self._exons = [
            (0, 0, IntervalType.EXON, 99, 200),
            (0, 0, IntervalType.EXON, 299, 400),
            (0, 0, IntervalType.EXON, 499, 600),
            (1, 0, IntervalType.EXON, 99, 200),
            (1, 0, IntervalType.EXON, 499, 600),
        ]
        # Introns for gap SJ queries
        self._introns = [
            (0, 0, IntervalType.INTRON, 200, 299),
            (0, 0, IntervalType.INTRON, 400, 499),
            (1, 0, IntervalType.INTRON, 200, 499),
        ]

    def query_exon(self, exon_block):
        results = []
        for t_idx, g_idx, itype, h_start, h_end in self._exons:
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((t_idx, g_idx, itype))
        return results

    def query_exon_with_coords(self, exon_block):
        results = []
        for t_idx, g_idx, itype, h_start, h_end in self._exons:
            if exon_block.ref == "chr1" and exon_block.start < h_end and exon_block.end > h_start:
                results.append((t_idx, g_idx, itype, h_start, h_end))
        return results

    def query_gap_sjs(self, ref, start, end):
        """Return introns fully contained in the gap."""
        results = []
        for t_idx, g_idx, itype, h_start, h_end in self._introns:
            if ref == "chr1" and h_start >= start and h_end <= end:
                results.append((t_idx, g_idx, itype, h_start, h_end))
        return results


class TestInsertSizeDiscrimination:
    """Verify insert size helps distinguish isoforms with different splicing.

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
        assert result.count_cat == CountCategory.SPLICED_ANNOT
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
        assert result.count_cat == CountCategory.SPLICED_ANNOT
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
        assert result.count_cat == CountCategory.UNSPLICED
        # Both transcripts share exon1 fully → equal overlap → both kept
        assert result.t_inds == frozenset({0, 1})

    def test_unspliced_read_in_t0_only_exon(self):
        """Unspliced read in exon2 (t0-only at 299-400) → t0 only."""
        idx = _ThreeExonMockIndex()
        frag = Fragment(
            exons=(GenomicInterval("chr1", 320, 380, Strand.POS),),
            introns=(),
        )
        result = resolve_fragment(frag, idx)
        assert result is not None
        assert result.count_cat == CountCategory.UNSPLICED
        assert result.t_inds == frozenset({0})

    def test_insert_size_differs_between_isoforms(self):
        """Fragment spanning exon1→exon3 has different insert size per isoform.

        For t0 (3 exons): reads at 150-200, 499-550 with gap 200-499
          → gap = 299bp, gap contains introns 200-299 (99bp) and 400-499 (99bp)
          → insert = (200-150) + (550-499) + (299 - 99 - 99) = 50 + 51 + 101 = 202
          Actually, insert = genomic_span - sum(gap_introns) = 400 - 198 = 202

        For t1 (2 exons):  gap contains intron 200-499 (299bp)
          → insert = 400 - 299 = 101

        These differ → compute_insert_size returns -1 when they disagree.
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
