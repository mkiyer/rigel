"""Tests for hulkrna.count — fragment resolution, model training, and counting."""

import pytest
from unittest.mock import MagicMock

from hulkrna.core import (
    EMPTY_MERGE,
    CountCategory,
    CountStrand,
    CountType,
    GenomicInterval,
    IntervalType,
    MergeCriteria,
    Strand,
)
from hulkrna.fragment import Fragment
from hulkrna.index import HulkIndex
from hulkrna.insert_model import InsertSizeModels
from hulkrna.count import (
    _compute_insert_size,
    _fragment_insert_size,
    _resolve_fragment,
    pass1_learn,
    pass2_count,
    ReadCounter,
)
from hulkrna.strand_model import StrandModel


# ---------------------------------------------------------------------------
# Helpers: mock pysam reads
# ---------------------------------------------------------------------------

_REF_IDS = {"chr1": 0, "chr2": 1}


def _mock_read(
    query_name: str,
    reference_name: str,
    reference_start: int,
    cigartuples: list[tuple[int, int]],
    is_reverse: bool = False,
    is_read1: bool = True,
    is_read2: bool = False,
    is_paired: bool = True,
    is_duplicate: bool = False,
    is_proper_pair: bool = True,
    is_unmapped: bool = False,
    is_qcfail: bool = False,
    is_secondary: bool = False,
    is_supplementary: bool = False,
    mate_is_unmapped: bool = False,
    nh: int = 1,
    xs: str = "+",
    reference_id: int | None = None,
    next_reference_id: int = 0,
    next_reference_start: int = 0,
) -> MagicMock:
    """Build a minimal mock pysam.AlignedSegment."""
    read = MagicMock()
    read.query_name = query_name
    read.reference_name = reference_name
    read.reference_start = reference_start
    read.cigartuples = cigartuples
    read.is_reverse = is_reverse
    read.is_read1 = is_read1
    read.is_read2 = is_read2
    read.is_paired = is_paired
    read.is_duplicate = is_duplicate
    read.is_proper_pair = is_proper_pair
    read.is_unmapped = is_unmapped
    read.is_qcfail = is_qcfail
    read.is_secondary = is_secondary
    read.is_supplementary = is_supplementary
    read.mate_is_unmapped = mate_is_unmapped

    if reference_id is None:
        reference_id = _REF_IDS.get(reference_name, 0)
    read.reference_id = reference_id
    read.next_reference_id = next_reference_id
    read.next_reference_start = next_reference_start

    def _get_tag(tag):
        if tag == "NH":
            return nh
        elif tag == "XS":
            return xs
        raise KeyError(tag)

    read.get_tag = _get_tag
    read.has_tag = lambda tag: tag in ("NH", "XS")
    return read


# CIGAR ops
M = 0   # BAM_CMATCH
N = 3   # BAM_CREF_SKIP


def _make_pair(qname, ref, r1_start, r2_start, r1_len=50, r2_len=50,
               r2_reverse=True, xs="+"):
    """Create a simple unspliced read pair with reciprocal mate pointers."""
    ref_id = _REF_IDS.get(ref, 0)
    r1 = _mock_read(
        query_name=qname, reference_name=ref,
        reference_start=r1_start,
        cigartuples=[(M, r1_len)],
        is_reverse=False, is_read1=True, is_read2=False,
        reference_id=ref_id,
        next_reference_id=ref_id,
        next_reference_start=r2_start,
        xs=xs,
    )
    r2 = _mock_read(
        query_name=qname, reference_name=ref,
        reference_start=r2_start,
        cigartuples=[(M, r2_len)],
        is_reverse=r2_reverse, is_read1=False, is_read2=True,
        reference_id=ref_id,
        next_reference_id=ref_id,
        next_reference_start=r1_start,
        xs=xs,
    )
    return r1, r2


def _make_spliced_pair(qname, ref, r1_start, sj_start, sj_end,
                       r2_start, exon1_len=50, exon2_len=50, r2_len=50,
                       xs="+"):
    """Create a spliced read pair (r1 has N-op, r2 is unspliced)."""
    ref_id = _REF_IDS.get(ref, 0)
    intron_len = sj_end - sj_start
    r1 = _mock_read(
        query_name=qname, reference_name=ref,
        reference_start=r1_start,
        cigartuples=[(M, exon1_len), (N, intron_len), (M, exon2_len)],
        is_reverse=False, is_read1=True, is_read2=False,
        reference_id=ref_id,
        next_reference_id=ref_id,
        next_reference_start=r2_start,
        xs=xs,
    )
    r2 = _mock_read(
        query_name=qname, reference_name=ref,
        reference_start=r2_start,
        cigartuples=[(M, r2_len)],
        is_reverse=True, is_read1=False, is_read2=True,
        reference_id=ref_id,
        next_reference_id=ref_id,
        next_reference_start=r1_start,
        xs=xs,
    )
    return r1, r2


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def mini_index(mini_fasta_file, mini_gtf_file, tmp_path):
    """Build and load a HulkIndex from the mini fixtures."""
    index_dir = tmp_path / "index"
    HulkIndex.build(
        fasta_file=mini_fasta_file,
        gtf_file=mini_gtf_file,
        output_dir=index_dir,
        write_tsv=False,
    )
    return HulkIndex.load(index_dir)


# ---------------------------------------------------------------------------
# Tests: _resolve_fragment
# ---------------------------------------------------------------------------

class TestResolveFragment:
    """Test _resolve_fragment with hand-crafted Fragments."""

    def test_single_exon_in_gene(self, mini_index):
        """A single exon block fully inside gene A exon 1 (99-200).

        Should resolve to UNSPLICED with compatible transcripts.
        """
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.count_cat == CountCategory.UNSPLICED
        assert result.exon_strand == Strand.POS
        # Gene A (g_index=0) should be in results
        assert 0 in result.merge_result.g_inds

    def test_intergenic_region(self, mini_index):
        """Exon block in intergenic region (chr2 has no genes) → None."""
        frag = Fragment(
            exons=(GenomicInterval("chr2", 100, 200, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is None

    def test_spliced_fragment_with_sj_match(self, mini_index):
        """Fragment with an intron matching annotated SJ → SPLICED_ANNOT.

        t1 has intron (200, 299) in 0-based half-open coordinates.
        """
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ),
        )
        result = _resolve_fragment(frag, mini_index)

        assert result is not None
        assert result.count_cat == CountCategory.SPLICED_ANNOT
        assert result.sj_strand == Strand.POS
        # t1 (t_index=0) should be in the compatible set
        assert 0 in result.merge_result.t_inds

    def test_spliced_fragment_no_sj_match(self, mini_index):
        """Fragment with intron that doesn't match → SPLICED_UNANNOT."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 180, Strand.POS),
                GenomicInterval("chr1", 250, 280, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 180, 250, Strand.POS),
            ),
        )
        result = _resolve_fragment(frag, mini_index)

        assert result is not None
        assert result.count_cat == CountCategory.SPLICED_UNANNOT

    def test_empty_fragment(self, mini_index):
        """A fragment with no exons produces None."""
        frag = Fragment()
        result = _resolve_fragment(frag, mini_index)
        assert result is None

    def test_intronic_fragment(self, mini_index):
        """Fragment fully in intronic region (200-299 for t1) → INTRON.

        Block 230-270 overlaps t1's intron (200,299).
        """
        frag = Fragment(
            exons=(GenomicInterval("chr1", 230, 270, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)

        assert result is not None
        assert result.count_cat == CountCategory.INTRON

    def test_unique_gene_resolution(self, mini_index):
        """Two exon blocks in gene A → unique gene merge."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 160, Strand.POS),
                GenomicInterval("chr1", 520, 560, Strand.POS),
            ),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)

        assert result is not None
        assert result.merge_result.is_unique_gene
        # All transcripts should be from gene A (g_index=0)
        assert result.merge_result.g_inds == frozenset({0})

    def test_frag_strand_bitwise_or(self, mini_index):
        """exon_strand should be bitwise OR of exon block strands."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.exon_strand == Strand.POS

    def test_sj_narrows_to_specific_transcript(self, mini_index):
        """SJ match narrows compatible set.

        t1 (200,299) and t2 share exon (99,200) and (499,600).
        But only t1 has intron (200,299).
        A fragment with SJ (200,299) should narrow to t1.
        """
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        # t1 (index 0) has this intron, t2 (index 1) does not
        assert 0 in result.merge_result.t_inds
        assert 1 not in result.merge_result.t_inds


# ---------------------------------------------------------------------------
# Tests: _compute_insert_size
# ---------------------------------------------------------------------------

class TestInsertSize:
    """Test insert size computation."""

    def test_single_exon_block(self, mini_index):
        """Single exon block: insert = block length (no gaps)."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        size = _compute_insert_size(frag, frozenset({0, 1}), mini_index)
        assert size == 60

    def test_spliced_fragment(self, mini_index):
        """Spliced fragment: footprint minus observed intron.

        Exon blocks: (150,200) and (299,350).
        Observed intron: (200,299) = 99bp.
        Footprint: 350-150 = 200.
        Insert size: 200 - 99 = 101.
        """
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ),
        )
        size = _compute_insert_size(frag, frozenset({0}), mini_index)
        assert size == 101

    def test_paired_end_gap_with_ref_intron(self, mini_index):
        """Paired-end with gap containing reference intron.

        Fragment: exon blocks (120,180) and (320,380) — no observed introns.
        Footprint: 380 - 120 = 260.
        Gap: (180, 320).
        t1 intron (200,299) = 99bp fully contained in gap.
        insert_t1: 260 - 99 = 161.
        """
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 180, Strand.POS),
                GenomicInterval("chr1", 320, 380, Strand.POS),
            ),
            introns=(),
        )
        size = _compute_insert_size(frag, frozenset({0}), mini_index)
        assert size == 161

    def test_different_transcripts_ambiguous(self, mini_index):
        """Different transcripts yield different insert sizes → -1.

        Fragment: exon blocks (120,180) and (520,560).
        Footprint: 560 - 120 = 440.
        Gap: (180, 520).
        t1: introns (200,299)=99bp + (400,499)=99bp → 440 - 198 = 242.
        t2: intron (200,499)=299bp → 440 - 299 = 141.
        Ambiguous → returns -1.
        """
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 180, Strand.POS),
                GenomicInterval("chr1", 520, 560, Strand.POS),
            ),
            introns=(),
        )
        size = _compute_insert_size(frag, frozenset({0, 1}), mini_index)
        assert size == -1

    def test_single_transcript_gap_correction(self, mini_index):
        """Single compatible transcript with gap → unambiguous."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 180, Strand.POS),
                GenomicInterval("chr1", 520, 560, Strand.POS),
            ),
            introns=(),
        )
        # t1 only: 2 introns in gap → 440 - 198 = 242
        size = _compute_insert_size(frag, frozenset({0}), mini_index)
        assert size == 242
        # t2 only: 1 intron in gap → 440 - 299 = 141
        size = _compute_insert_size(frag, frozenset({1}), mini_index)
        assert size == 141

    def test_empty_compatible_set(self, mini_index):
        """Empty compatible set → -1."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        size = _compute_insert_size(frag, frozenset(), mini_index)
        assert size == -1

    def test_empty_fragment(self, mini_index):
        """Empty fragment → -1."""
        frag = Fragment()
        size = _compute_insert_size(frag, frozenset({0}), mini_index)
        assert size == -1

    def test_no_gap_all_same(self, mini_index):
        """No gaps → all compatible transcripts agree on insert size."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 160, Strand.POS),
                GenomicInterval("chr1", 140, 180, Strand.POS),
            ),
            introns=(),
        )
        size = _compute_insert_size(frag, frozenset({0, 1}), mini_index)
        assert size == 60

    def test_fragment_insert_size_simple(self, mini_index):
        """_fragment_insert_size: footprint minus observed introns."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 100, 200, Strand.POS),
                GenomicInterval("chr1", 300, 400, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 300, Strand.POS),
            ),
        )
        assert _fragment_insert_size(frag) == 200  # 300 - 100 = 200

    def test_fragment_insert_size_empty(self, mini_index):
        """_fragment_insert_size: empty fragment → -1."""
        frag = Fragment()
        assert _fragment_insert_size(frag) == -1


# ---------------------------------------------------------------------------
# Tests: ReadCounter
# ---------------------------------------------------------------------------

class TestReadCounter:
    """Test ReadCounter assignment logic."""

    def test_init_shapes(self):
        """Arrays have correct shapes."""
        counter = ReadCounter(100, 10)
        assert counter.t_counts.shape == (100, 12)
        assert counter.g_counts.shape == (10, 12)
        assert counter.t_counts.sum() == 0.0
        assert counter.g_counts.sum() == 0.0

    def test_unique_gene_sense(self, mini_index):
        """Unique gene with same-strand → SENSE count."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.merge_result.is_unique_gene

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index)

        # Gene A (g_index=0) is on + strand, fragment is POS → SENSE
        g_idx = next(iter(result.merge_result.g_inds))
        # count_type = cat * 3 + strand
        count_type = int(result.count_cat) * 3 + int(CountStrand.SENSE)
        assert counter.g_counts[g_idx, count_type] > 0.0

    def test_unique_gene_antisense(self, mini_index):
        """Unique gene with opposite-strand → ANTISENSE count."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.NEG),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.merge_result.is_unique_gene

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index)

        g_idx = next(iter(result.merge_result.g_inds))
        count_type = int(result.count_cat) * 3 + int(CountStrand.ANTISENSE)
        assert counter.g_counts[g_idx, count_type] > 0.0

    def test_ambiguous_strand(self, mini_index):
        """Fragment with AMBIGUOUS exon_strand → AMBIGUOUS count."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.AMBIGUOUS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index)

        g_idx = next(iter(result.merge_result.g_inds))
        count_type = int(result.count_cat) * 3 + int(CountStrand.AMBIGUOUS)
        assert counter.g_counts[g_idx, count_type] > 0.0

    def test_fractional_transcript_assignment(self, mini_index):
        """Multiple compatible transcripts → each gets 1/N."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        n_t = len(result.merge_result.t_inds)
        assert n_t > 1  # gene A has multiple transcripts

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index)

        # Each transcript should get 1/n_t
        for t_idx in result.merge_result.t_inds:
            assert counter.t_counts[t_idx].sum() == pytest.approx(1.0 / n_t)

    def test_gene_gets_full_count(self, mini_index):
        """Unique gene → gene gets 1.0 total count."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None
        assert result.merge_result.is_unique_gene

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index)

        g_idx = next(iter(result.merge_result.g_inds))
        assert counter.g_counts[g_idx].sum() == pytest.approx(1.0)

    def test_multimap_weighting(self, mini_index):
        """num_hits > 1 → each alignment pair weighted by 1/num_hits."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        result = _resolve_fragment(frag, mini_index)
        assert result is not None

        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        counter.assign(result, mini_index, num_hits=2)

        g_idx = next(iter(result.merge_result.g_inds))
        assert counter.g_counts[g_idx].sum() == pytest.approx(0.5)

    def test_get_dataframes(self, mini_index):
        """get_t_counts_df and get_g_counts_df return proper DataFrames."""
        counter = ReadCounter(mini_index.num_transcripts, mini_index.num_genes)
        t_df = counter.get_t_counts_df()
        g_df = counter.get_g_counts_df()
        assert list(t_df.columns) == list(CountType.columns())
        assert list(g_df.columns) == list(CountType.columns())
        assert len(t_df) == mini_index.num_transcripts
        assert len(g_df) == mini_index.num_genes


# ---------------------------------------------------------------------------
# Tests: pass1_learn (integration)
# ---------------------------------------------------------------------------

class TestPass1Learn:
    """Integration tests for pass1_learn using mock BAM iterators."""

    def test_returns_tuple(self, mini_index):
        """pass1_learn returns (stats, strand_model, insert_models)."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        result = pass1_learn([r1, r2], mini_index)
        assert isinstance(result, tuple)
        assert len(result) == 3
        stats, strand_model, insert_models = result
        assert isinstance(stats, dict)
        assert isinstance(strand_model, StrandModel)
        assert isinstance(insert_models, InsertSizeModels)

    def test_single_pair(self, mini_index):
        """One read pair → n_fragments=1."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        stats, _, _ = pass1_learn([r1, r2], mini_index)
        assert stats["unique"] == 1
        assert stats["n_fragments"] == 1

    def test_intergenic_counted(self, mini_index):
        """Intergenic fragment → n_intergenic stat."""
        r1, r2 = _make_pair("frag_ig", "chr2", 100, 200)
        stats, _, _ = pass1_learn([r1, r2], mini_index)
        assert stats["n_intergenic"] == 1

    def test_chimeric_skipped(self, mini_index):
        """Chimeric fragment (multi-ref) → n_chimeric stat."""
        ref_id_1 = _REF_IDS["chr1"]
        ref_id_2 = _REF_IDS["chr2"]
        r1 = _mock_read(
            query_name="chimera", reference_name="chr1",
            reference_start=120, cigartuples=[(M, 50)],
            is_reverse=False, is_read1=True, is_read2=False,
            reference_id=ref_id_1,
            next_reference_id=ref_id_2,
            next_reference_start=100,
        )
        r2 = _mock_read(
            query_name="chimera", reference_name="chr2",
            reference_start=100, cigartuples=[(M, 50)],
            is_reverse=True, is_read1=False, is_read2=True,
            reference_id=ref_id_2,
            next_reference_id=ref_id_1,
            next_reference_start=120,
        )
        stats, _, _ = pass1_learn([r1, r2], mini_index)
        assert stats["n_chimeric"] == 1

    def test_empty_bam(self, mini_index):
        """Empty BAM → zero stats."""
        stats, strand_model, insert_models = pass1_learn(iter([]), mini_index)
        assert stats["n_fragments"] == 0
        assert strand_model.n_observations == 0
        assert insert_models.n_observations == 0

    def test_duplicate_skipped(self, mini_index):
        """Duplicate reads should be skipped by default."""
        ref_id = _REF_IDS["chr1"]
        r1 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=120, cigartuples=[(M, 50)],
            is_reverse=False, is_read1=True, is_read2=False,
            is_duplicate=True,
            reference_id=ref_id, next_reference_id=ref_id,
            next_reference_start=150,
        )
        r2 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=150, cigartuples=[(M, 50)],
            is_reverse=True, is_read1=False, is_read2=True,
            is_duplicate=True,
            reference_id=ref_id, next_reference_id=ref_id,
            next_reference_start=120,
        )
        stats, _, _ = pass1_learn([r1, r2], mini_index, skip_duplicates=True)
        assert stats["n_fragments"] == 0 or stats.get("duplicate", 0) > 0

    def test_stats_keys_present(self, mini_index):
        """Stats dict should have all resolution counter keys."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        stats, _, _ = pass1_learn([r1, r2], mini_index)

        expected_keys = [
            "n_fragments", "n_chimeric", "n_intergenic",
            "n_with_exon", "n_with_intron_fallback",
            "n_with_annotated_sj", "n_with_unannotated_sj",
            "n_unique_gene", "n_multi_gene",
            "n_strand_skipped_no_sj", "n_strand_skipped_multi_gene",
            "n_strand_skipped_ambiguous",
            "n_insert_unambiguous", "n_insert_ambiguous",
            "n_insert_intergenic",
        ]
        for key in expected_keys:
            assert key in stats, f"Missing stats key: {key}"


# ---------------------------------------------------------------------------
# Tests: Strand model training
# ---------------------------------------------------------------------------

class TestStrandModelTraining:
    """Test that strand model is trained during pass1_learn."""

    def test_spliced_unique_gene_trains_model(self, mini_index):
        """Spliced fragment with annotated SJ + unique gene → strand obs.

        t1 intron (200,299) on + strand. Fragment with SJ matching t1.
        """
        r1, r2 = _make_spliced_pair(
            "spliced1", "chr1",
            r1_start=150, sj_start=200, sj_end=299,
            r2_start=320, exon1_len=50, exon2_len=50, r2_len=50,
        )
        stats, strand_model, _ = pass1_learn([r1, r2], mini_index)
        assert strand_model.n_observations > 0

    def test_unspliced_does_not_train_strand(self, mini_index):
        """Unspliced fragment → no strand observations."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        stats, strand_model, _ = pass1_learn([r1, r2], mini_index)
        assert strand_model.n_observations == 0
        assert stats["n_strand_skipped_no_sj"] == 1


# ---------------------------------------------------------------------------
# Tests: Insert size model training
# ---------------------------------------------------------------------------

class TestInsertModelTraining:
    """Test that insert size model is trained during pass1_learn."""

    def test_unambiguous_insert_trains_model(self, mini_index):
        """Single exon block → all compatible transcripts same insert → trains."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        stats, _, insert_models = pass1_learn([r1, r2], mini_index)
        assert insert_models.n_observations > 0
        assert stats["n_insert_unambiguous"] > 0


# ---------------------------------------------------------------------------
# Tests: pass2_count (integration)
# ---------------------------------------------------------------------------

class TestPass2Count:
    """Integration tests for pass2_count."""

    def test_returns_tuple(self, mini_index):
        """pass2_count returns (stats, ReadCounter)."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        result = pass2_count(
            [r1, r2], mini_index, strand_model, insert_models
        )
        assert isinstance(result, tuple)
        assert len(result) == 2
        stats, counter = result
        assert isinstance(stats, dict)
        assert isinstance(counter, ReadCounter)

    def test_single_pair_counted(self, mini_index):
        """One read pair → counts assigned."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        stats, counter = pass2_count(
            [r1, r2], mini_index, strand_model, insert_models
        )
        assert stats["n_fragments"] == 1
        # Gene A should have non-zero count
        assert counter.g_counts[0].sum() > 0.0

    def test_intergenic_not_counted(self, mini_index):
        """Intergenic fragment → no counts assigned."""
        r1, r2 = _make_pair("frag_ig", "chr2", 100, 200)
        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        stats, counter = pass2_count(
            [r1, r2], mini_index, strand_model, insert_models
        )
        assert stats["n_intergenic"] == 1
        assert counter.g_counts.sum() == 0.0

    def test_empty_bam(self, mini_index):
        """Empty BAM → zero counts."""
        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        stats, counter = pass2_count(
            iter([]), mini_index, strand_model, insert_models
        )
        assert stats["n_fragments"] == 0
        assert counter.g_counts.sum() == 0.0
        assert counter.t_counts.sum() == 0.0

    def test_multiple_fragments(self, mini_index):
        """Multiple fragments → counts accumulate."""
        pairs = []
        for i in range(5):
            r1, r2 = _make_pair(f"frag{i}", "chr1", 120, 150)
            pairs.extend([r1, r2])

        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        stats, counter = pass2_count(
            pairs, mini_index, strand_model, insert_models
        )
        assert stats["n_fragments"] == 5
        # gene A should have total ~5.0 counts
        assert counter.g_counts[0].sum() == pytest.approx(5.0)

    def test_pass2_stats_keys(self, mini_index):
        """Pass 2 stats dict should have resolution keys."""
        r1, r2 = _make_pair("frag1", "chr1", 120, 150)
        strand_model = StrandModel()
        insert_models = InsertSizeModels()
        stats, _ = pass2_count(
            [r1, r2], mini_index, strand_model, insert_models
        )
        expected_keys = [
            "n_fragments", "n_chimeric", "n_intergenic",
            "n_with_exon", "n_with_intron_fallback",
            "n_unique_gene", "n_multi_gene",
        ]
        for key in expected_keys:
            assert key in stats, f"Missing stats key: {key}"


# ---------------------------------------------------------------------------
# Tests: two-pass consistency
# ---------------------------------------------------------------------------

class TestTwoPassConsistency:
    """Verify that pass1 and pass2 agree on resolution stats."""

    def test_resolution_stats_match(self, mini_index):
        """Pass 1 and Pass 2 should produce identical resolution stats."""
        pairs = []
        for i in range(5):
            r1, r2 = _make_pair(f"frag{i}", "chr1", 120, 150)
            pairs.extend([r1, r2])

        # Run pass 1
        pass1_stats, strand_model, insert_models = pass1_learn(
            list(pairs), mini_index
        )

        # Run pass 2 with same input
        pass2_stats, counter = pass2_count(
            list(pairs), mini_index, strand_model, insert_models
        )

        # Resolution stats should match
        shared_keys = [
            "n_fragments", "n_chimeric", "n_intergenic",
            "n_with_exon", "n_with_intron_fallback",
            "n_unique_gene", "n_multi_gene",
        ]
        for key in shared_keys:
            assert pass1_stats[key] == pass2_stats[key], (
                f"Mismatch on {key}: pass1={pass1_stats[key]}, "
                f"pass2={pass2_stats[key]}"
            )
