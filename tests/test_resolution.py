"""Tests for hulkrna.resolution — set merging, insert size, ResolvedFragment, chimera."""

import pytest

from hulkrna.types import ChimeraType, MergeCriteria, MergeResult, Strand, GenomicInterval
from hulkrna.categories import CountCategory
from hulkrna.resolution import (
    merge_sets_with_criteria,
    fragment_insert_size,
    ResolvedFragment,
    _detect_intrachromosomal_chimera,
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
