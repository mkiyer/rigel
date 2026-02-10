"""Tests for hulkrna.index — reference index building."""

import pandas as pd
import pytest

from hulkrna.core import GenomicInterval, IntervalType, Strand
from hulkrna.index import (
    REF_LENGTHS_FEATHER,
    REF_LENGTHS_TSV,
    TRANSCRIPTS_FEATHER,
    TRANSCRIPTS_TSV,
    INTERVALS_FEATHER,
    INTERVALS_TSV,
    SJ_FEATHER,
    SJ_TSV,
    HulkIndex,
    build_genomic_intervals,
    build_splice_junctions,
    load_reference_lengths,
    read_transcripts,
    transcripts_to_dataframe,
)


class TestLoadReferenceLengths:
    def test_basic(self, mini_fasta_file):
        ref_lengths = load_reference_lengths(mini_fasta_file)
        assert ref_lengths == {"chr1": 2000, "chr2": 500}

    def test_order_preserved(self, mini_fasta_file):
        ref_lengths = load_reference_lengths(mini_fasta_file)
        assert list(ref_lengths.keys()) == ["chr1", "chr2"]


class TestReadTranscripts:
    def test_count(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        assert len(transcripts) == 3

    def test_t_index_sequential(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        assert [t.t_index for t in transcripts] == [0, 1, 2]

    def test_g_index_assigned(self, mini_gtf_file):
        """g1 and g2 get distinct g_index values."""
        transcripts = read_transcripts(mini_gtf_file)
        g_indices = {t.g_id: t.g_index for t in transcripts}
        # Two distinct genes
        assert len(set(g_indices.values())) == 2
        # t1 and t2 share g_index (same gene g1)
        t1 = next(t for t in transcripts if t.t_id == "t1")
        t2 = next(t for t in transcripts if t.t_id == "t2")
        assert t1.g_index == t2.g_index

    def test_sorted_by_position(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        keys = [(t.ref, t.start, t.end, t.strand) for t in transcripts]
        assert keys == sorted(keys)


class TestTranscriptsToDataframe:
    def test_shape(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        df = transcripts_to_dataframe(transcripts)
        assert len(df) == 3
        assert "t_id" in df.columns
        assert "g_id" in df.columns
        assert "ref" in df.columns

    def test_values(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        df = transcripts_to_dataframe(transcripts)
        assert list(df["t_id"]) == [t.t_id for t in transcripts]


class TestBuildSpliceJunctions:
    def test_count(self, mini_gtf_file):
        """t1 has 2 introns, t2 has 1, t3 has 1 → total 4."""
        transcripts = read_transcripts(mini_gtf_file)
        sj_df = build_splice_junctions(transcripts)
        assert len(sj_df) == 4

    def test_columns(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        sj_df = build_splice_junctions(transcripts)
        expected = [
            "ref", "start", "end", "strand",
            "interval_type", "t_index", "g_index",
        ]
        assert list(sj_df.columns) == expected

    def test_coordinates(self, mini_gtf_file):
        """Check intron coords for t1: (200,299) and (400,499)."""
        transcripts = read_transcripts(mini_gtf_file)
        sj_df = build_splice_junctions(transcripts)
        t1 = next(t for t in transcripts if t.t_id == "t1")
        t1_sjs = sj_df[sj_df["t_index"] == t1.t_index]
        starts = sorted(t1_sjs["start"].tolist())
        ends = sorted(t1_sjs["end"].tolist())
        assert starts == [200, 400]
        assert ends == [299, 499]

    def test_strand_preserved(self, mini_gtf_file):
        transcripts = read_transcripts(mini_gtf_file)
        sj_df = build_splice_junctions(transcripts)
        t3 = next(t for t in transcripts if t.t_id == "t3")
        t3_sjs = sj_df[sj_df["t_index"] == t3.t_index]
        assert all(t3_sjs["strand"] == Strand.NEG)


class TestBuildGenomicIntervals:
    def test_covers_genome(self, mini_gtf_file, mini_fasta_file):
        """Intervals should tile the genome without gaps."""
        ref_lengths = load_reference_lengths(mini_fasta_file)
        transcripts = read_transcripts(mini_gtf_file)
        iv_df = build_genomic_intervals(transcripts, ref_lengths)

        for ref, length in ref_lengths.items():
            ref_ivs = iv_df[iv_df["ref"] == ref].sort_values("start")
            # First interval starts at 0
            assert ref_ivs.iloc[0]["start"] == 0, f"{ref}: first start != 0"
            # Last interval ends at ref length
            assert ref_ivs.iloc[-1]["end"] == length, f"{ref}: last end != {length}"

    def test_intergenic_on_empty_chr(self, mini_gtf_file, mini_fasta_file):
        """chr2 has no genes — should be one intergenic interval [0, 500)."""
        ref_lengths = load_reference_lengths(mini_fasta_file)
        transcripts = read_transcripts(mini_gtf_file)
        iv_df = build_genomic_intervals(transcripts, ref_lengths)

        chr2 = iv_df[iv_df["ref"] == "chr2"]
        assert len(chr2) == 1
        row = chr2.iloc[0]
        assert row["start"] == 0
        assert row["end"] == 500
        assert row["interval_type"] == IntervalType.INTERGENIC

    def test_has_exon_intervals(self, mini_gtf_file, mini_fasta_file):
        ref_lengths = load_reference_lengths(mini_fasta_file)
        transcripts = read_transcripts(mini_gtf_file)
        iv_df = build_genomic_intervals(transcripts, ref_lengths)

        exons = iv_df[iv_df["interval_type"] == IntervalType.EXON]
        assert len(exons) > 0

    def test_has_intron_intervals(self, mini_gtf_file, mini_fasta_file):
        ref_lengths = load_reference_lengths(mini_fasta_file)
        transcripts = read_transcripts(mini_gtf_file)
        iv_df = build_genomic_intervals(transcripts, ref_lengths)

        introns = iv_df[iv_df["interval_type"] == IntervalType.INTRON]
        assert len(introns) > 0

    def test_has_intergenic_intervals(self, mini_gtf_file, mini_fasta_file):
        ref_lengths = load_reference_lengths(mini_fasta_file)
        transcripts = read_transcripts(mini_gtf_file)
        iv_df = build_genomic_intervals(transcripts, ref_lengths)

        intergenic = iv_df[iv_df["interval_type"] == IntervalType.INTERGENIC]
        # At minimum: before first gene, between g1 and g2 cluster, after g2, chr2
        assert len(intergenic) >= 2

    def test_missing_reference_raises(self, mini_gtf_file):
        """If a transcript references a chromosome not in FASTA, raise."""
        transcripts = read_transcripts(mini_gtf_file)
        fake_refs = {"chrX": 1000}  # chr1 missing!
        with pytest.raises(ValueError, match="not found in the FASTA"):
            build_genomic_intervals(transcripts, fake_refs)


class TestBuildIndex:
    """End-to-end test of the HulkIndex.build orchestrator."""

    def test_output_files_created(self, mini_gtf_file, mini_fasta_file, tmp_path):
        out_dir = tmp_path / "index_output"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)

        for fname in [
            REF_LENGTHS_FEATHER, TRANSCRIPTS_FEATHER,
            INTERVALS_FEATHER, SJ_FEATHER,
        ]:
            assert (out_dir / fname).exists(), f"{fname} missing"

    def test_tsv_mirrors(self, mini_gtf_file, mini_fasta_file, tmp_path):
        out_dir = tmp_path / "index_output"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir, write_tsv=True)

        for fname in [
            REF_LENGTHS_TSV, TRANSCRIPTS_TSV,
            INTERVALS_TSV, SJ_TSV,
        ]:
            assert (out_dir / fname).exists(), f"{fname} missing"

    def test_no_tsv_when_disabled(self, mini_gtf_file, mini_fasta_file, tmp_path):
        out_dir = tmp_path / "index_output"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir, write_tsv=False)

        for fname in [
            REF_LENGTHS_TSV, TRANSCRIPTS_TSV,
            INTERVALS_TSV, SJ_TSV,
        ]:
            assert not (out_dir / fname).exists(), f"{fname} should not exist"

    def test_feather_readable(self, mini_gtf_file, mini_fasta_file, tmp_path):
        """All Feather outputs should be readable as DataFrames."""
        out_dir = tmp_path / "index_output"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)

        ref_df = pd.read_feather(out_dir / REF_LENGTHS_FEATHER)
        assert len(ref_df) == 2
        assert list(ref_df.columns) == ["ref", "length"]

        t_df = pd.read_feather(out_dir / TRANSCRIPTS_FEATHER)
        assert len(t_df) == 3

        sj_df = pd.read_feather(out_dir / SJ_FEATHER)
        assert len(sj_df) == 4

        iv_df = pd.read_feather(out_dir / INTERVALS_FEATHER)
        assert len(iv_df) > 0

    def test_creates_output_dir(self, mini_gtf_file, mini_fasta_file, tmp_path):
        """output_dir should be created if it doesn't exist."""
        out_dir = tmp_path / "deeply" / "nested" / "dir"
        assert not out_dir.exists()
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)
        assert out_dir.exists()


class TestHulkIndexLoad:
    """Test HulkIndex.load() round-trip: build → load → inspect."""

    @pytest.fixture
    def hulk_index(self, mini_gtf_file, mini_fasta_file, tmp_path):
        """Build and load a HulkIndex from the mini fixtures."""
        out_dir = tmp_path / "idx"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)
        return HulkIndex.load(out_dir)

    def test_num_transcripts(self, hulk_index):
        assert hulk_index.num_transcripts == 3

    def test_num_genes(self, hulk_index):
        assert hulk_index.num_genes == 2

    def test_t_df_populated(self, hulk_index):
        assert hulk_index.t_df is not None
        assert "t_id" in hulk_index.t_df.columns
        assert list(hulk_index.t_df["t_id"]) == ["t1", "t2", "t3"]

    def test_g_df_populated(self, hulk_index):
        assert hulk_index.g_df is not None
        assert "g_id" in hulk_index.g_df.columns
        assert set(hulk_index.g_df["g_id"]) == {"g1", "g2"}

    def test_gene_table_aggregation(self, hulk_index):
        """Gene g1 has 2 transcripts, g2 has 1."""
        g1 = hulk_index.g_df[hulk_index.g_df["g_id"] == "g1"].iloc[0]
        g2 = hulk_index.g_df[hulk_index.g_df["g_id"] == "g2"].iloc[0]
        assert g1["num_transcripts"] == 2
        assert g2["num_transcripts"] == 1

    def test_lookup_arrays(self, hulk_index):
        assert hulk_index.t_to_g_arr is not None
        assert len(hulk_index.t_to_g_arr) == 3
        assert hulk_index.t_to_strand_arr is not None
        assert hulk_index.g_to_strand_arr is not None

    def test_cr_index_built(self, hulk_index):
        """Unified cgranges index should be populated."""
        assert hulk_index.cr is not None

    def test_sj_map_built(self, hulk_index):
        assert hulk_index.sj_map is not None
        assert len(hulk_index.sj_map) > 0

    def test_sj_cr_built(self, hulk_index):
        """SJ cgranges index should be populated."""
        assert hulk_index.sj_cr is not None

    def test_iv_lookup_arrays(self, hulk_index):
        """Interval lookup arrays should be populated with matching lengths."""
        assert hulk_index._iv_t_index is not None
        assert hulk_index._iv_g_index is not None
        assert hulk_index._iv_type is not None
        n = len(hulk_index._iv_t_index)
        assert len(hulk_index._iv_g_index) == n
        assert len(hulk_index._iv_type) == n

    def test_sj_lookup_arrays(self, hulk_index):
        """SJ lookup arrays should be populated with matching lengths."""
        assert hulk_index._sj_t_index is not None
        assert hulk_index._sj_g_index is not None
        assert hulk_index._sj_strand is not None
        n = len(hulk_index._sj_t_index)
        assert len(hulk_index._sj_g_index) == n
        assert len(hulk_index._sj_strand) == n


# ---------------------------------------------------------------------------
# query_exon() — unified interval overlap queries
# ---------------------------------------------------------------------------
# Mini-GTF reference (0-based half-open coordinates):
#
#   chr1 layout:
#     0..99     intergenic
#     99..200   exon  (t1,t2 → g1, + strand)
#     200..299  intron (t1 → g1)         / intergenic gap for t2
#     299..400  exon  (t1 → g1)
#     400..499  intron (t1 → g1)         / intergenic gap for t2
#     499..600  exon  (t1,t2 → g1)
#     600..999  intergenic
#     999..1100 exon  (t3 → g2, - strand)
#     1100..1199 intron (t3 → g2)
#     1199..1300 exon  (t3 → g2)
#     1300..2000 intergenic
#
#   chr2 layout:
#     0..500    intergenic (no genes)
#
#   SJs (introns):
#     t1: (200,299) + strand, (400,499) + strand
#     t2: (200,499) + strand
#     t3: (1100,1199) - strand

class TestQueryExon:
    """Test HulkIndex.query_exon() against the mini fixtures."""

    @pytest.fixture
    def hulk_index(self, mini_gtf_file, mini_fasta_file, tmp_path):
        out_dir = tmp_path / "idx"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)
        return HulkIndex.load(out_dir)

    def test_returns_tuples(self, hulk_index):
        """Results should be tuples (t_index, g_index, interval_type)."""
        hits = hulk_index.query_exon(GenomicInterval("chr2", 100, 200))
        assert len(hits) >= 1
        assert isinstance(hits[0], tuple)
        assert len(hits[0]) == 3

    def test_guaranteed_hit_on_valid_query(self, hulk_index):
        """Any query within a tiled reference must return ≥1 hit."""
        # Pure intergenic
        hits = hulk_index.query_exon(GenomicInterval("chr2", 50, 150))
        assert len(hits) >= 1
        # Pure exonic
        hits = hulk_index.query_exon(GenomicInterval("chr1", 120, 180))
        assert len(hits) >= 1

    def test_pure_intergenic(self, hulk_index):
        """Query in chr2 (no genes) should return only INTERGENIC hits."""
        hits = hulk_index.query_exon(GenomicInterval("chr2", 100, 200))
        # tuple format: (t_index, g_index, interval_type)
        assert all(h[2] == IntervalType.INTERGENIC for h in hits)
        assert all(h[0] == -1 for h in hits)
        assert all(h[1] == -1 for h in hits)

    def test_pure_exonic(self, hulk_index):
        """Query fully within t1/t2's shared exon (99,200).

        The exon appears once per transcript (t1 and t2 both have it),
        so we expect 2 EXON hits.
        """
        hits = hulk_index.query_exon(GenomicInterval("chr1", 120, 180))
        # tuple format: (t_index, g_index, interval_type)
        exon_hits = [h for h in hits if h[2] == IntervalType.EXON]
        # t1 and t2 both have exon (99,200) → 2 hits
        assert len(exon_hits) == 2
        # No intergenic hits
        ig_hits = [h for h in hits if h[2] == IntervalType.INTERGENIC]
        assert len(ig_hits) == 0

    def test_exon_intergenic_boundary(self, hulk_index):
        """Query spanning exon (99,200) and intergenic before it (0,99).

        A 150bp query from (50, 200) overlaps:
          intergenic (50..99)
          exon (99..200) × 2 transcripts (t1 and t2)
        """
        hits = hulk_index.query_exon(GenomicInterval("chr1", 50, 200))
        exon_hits = [h for h in hits if h[2] == IntervalType.EXON]
        ig_hits = [h for h in hits if h[2] == IntervalType.INTERGENIC]
        # 2 transcripts
        assert len(exon_hits) == 2
        # At least 1 intergenic hit
        assert len(ig_hits) >= 1

    def test_exon_intron_boundary(self, hulk_index):
        """Query spanning exon (99,200) and intron (200,299).

        Query (150, 250):
          exon (150..200) × 2 transcripts (t1, t2)
          intron from t1 (200..250)
          intron from t2 (200..250) — t2 intron is (200,499)
        """
        hits = hulk_index.query_exon(GenomicInterval("chr1", 150, 250))
        exon_hits = [h for h in hits if h[2] == IntervalType.EXON]
        intron_hits = [h for h in hits if h[2] == IntervalType.INTRON]
        # 2 exon hits (t1 and t2)
        assert len(exon_hits) == 2
        # 2 intron hits (t1: 200-299, t2: 200-499)
        assert len(intron_hits) == 2

    def test_always_has_hits(self, hulk_index):
        """Any query in a tiled reference always returns ≥1 hit."""
        queries = [
            GenomicInterval("chr1", 0, 100),
            GenomicInterval("chr1", 50, 200),
            GenomicInterval("chr1", 150, 350),
            GenomicInterval("chr1", 500, 1200),
            GenomicInterval("chr2", 0, 500),
        ]
        for q in queries:
            hits = hulk_index.query_exon(q)
            assert len(hits) >= 1, (
                f"Query ({q.ref}, {q.start}, {q.end}): no hits"
            )

    def test_mostly_intergenic_read(self, hulk_index):
        """A read with overlap in both intergenic and exonic regions.

        Query (10, 110): intergenic (10..99) + exon (99..110) × 2.
        """
        hits = hulk_index.query_exon(GenomicInterval("chr1", 10, 110))
        ig_hits = [h for h in hits if h[2] == IntervalType.INTERGENIC]
        exon_hits = [h for h in hits if h[2] == IntervalType.EXON]
        assert len(ig_hits) >= 1
        # 2 transcripts (t1, t2)
        assert len(exon_hits) == 2

    def test_gene2_exon(self, hulk_index):
        """Query within g2/t3's first exon (999,1100), neg strand."""
        hits = hulk_index.query_exon(GenomicInterval("chr1", 1020, 1080))
        exon_hits = [h for h in hits if h[2] == IntervalType.EXON]
        assert len(exon_hits) >= 1
        # All exon hits should belong to g2 (g_index = 1)
        g2_idx = hulk_index.g_df[
            hulk_index.g_df["g_id"] == "g2"
        ].iloc[0]["g_index"]
        for h in exon_hits:
            assert h[1] == g2_idx  # h[1] is g_index


# ---------------------------------------------------------------------------
# query_gap_sjs() — splice-junction containment queries
# ---------------------------------------------------------------------------

class TestQueryGapSJs:
    """Test HulkIndex.query_gap_sjs() for insert-size computation."""

    @pytest.fixture
    def hulk_index(self, mini_gtf_file, mini_fasta_file, tmp_path):
        out_dir = tmp_path / "idx"
        HulkIndex.build(mini_fasta_file, mini_gtf_file, out_dir)
        return HulkIndex.load(out_dir)

    def test_no_sjs_in_intergenic_gap(self, hulk_index):
        """Query in intergenic space should return no SJs."""
        hits = hulk_index.query_gap_sjs("chr2", 0, 500)
        assert len(hits) == 0

    def test_single_sj_contained(self, hulk_index):
        """Gap (150, 350) should contain t1's first intron (200, 299)."""
        hits = hulk_index.query_gap_sjs("chr1", 150, 350)
        assert len(hits) >= 1
        # At least one SJ should be (200, 299)
        # tuple format: (t_index, g_index, strand, sj_start, sj_end)
        coords = [(h[3], h[4]) for h in hits]
        assert (200, 299) in coords

    def test_returns_tuples_sj(self, hulk_index):
        """Results should be tuples (t_index, g_index, strand, sj_start, sj_end)."""
        hits = hulk_index.query_gap_sjs("chr1", 150, 350)
        for h in hits:
            assert isinstance(h, tuple)
            assert len(h) == 5

    def test_multiple_sjs_contained(self, hulk_index):
        """Gap (150, 550) should contain both t1 introns (200,299) and (400,499)
        plus t2's intron (200,499)."""
        hits = hulk_index.query_gap_sjs("chr1", 150, 550)
        coords = {(h[3], h[4]) for h in hits}
        assert (200, 299) in coords  # t1 intron 1
        assert (400, 499) in coords  # t1 intron 2
        assert (200, 499) in coords  # t2 intron

    def test_partial_overlap_excluded(self, hulk_index):
        """Gap (250, 350) should NOT contain (200, 299) since 200 < 250."""
        hits = hulk_index.query_gap_sjs("chr1", 250, 350)
        coords = [(h[3], h[4]) for h in hits]
        assert (200, 299) not in coords

    def test_exact_boundary_included(self, hulk_index):
        """Gap exactly matching an SJ's boundaries should include it."""
        hits = hulk_index.query_gap_sjs("chr1", 200, 299)
        coords = [(h[3], h[4]) for h in hits]
        assert (200, 299) in coords

    def test_sj_strand_preserved(self, hulk_index):
        """SJ hits should carry the correct XS-tag strand."""
        # t3's intron (1100, 1199) is on neg strand
        hits = hulk_index.query_gap_sjs("chr1", 1050, 1250)
        # tuple format: (t_index, g_index, strand, sj_start, sj_end)
        neg_hits = [h for h in hits if h[2] == Strand.NEG]
        assert len(neg_hits) >= 1
        assert any(h[3] == 1100 and h[4] == 1199 for h in neg_hits)

    def test_no_match_in_empty_region(self, hulk_index):
        """Region with no annotated SJs should return empty list."""
        hits = hulk_index.query_gap_sjs("chr1", 600, 999)
        assert len(hits) == 0

    def test_insert_size_subtraction(self, hulk_index):
        """Demonstrate the insert-size use case.

        Read 1: chr1:100-200  Read 2: chr1:500-600  (proper pair in g1)
        Gap: (200, 500)
        Genomic distance: 600 - 100 = 500
        Contained SJs in gap: (200,299)=99bp, (400,499)=99bp, (200,499)=299bp
        The correct SJ set depends on transcript assignment, but we should
        find all candidates.
        """
        hits = hulk_index.query_gap_sjs("chr1", 200, 500)
        # Should find t1's (200,299) and (400,499), plus t2's (200,499)
        coords = {(h[3], h[4]) for h in hits}
        assert (200, 299) in coords
        assert (400, 499) in coords
        assert (200, 499) in coords
