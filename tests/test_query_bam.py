"""Tests for hulkrna.query_bam — BAM scanning and index querying (Pass 1)."""

import pandas as pd
import pyarrow as pa
import pytest
from unittest.mock import MagicMock

from hulkrna.core import GenomicInterval, IntervalType, Strand
from hulkrna.fragment import Fragment
from hulkrna.index import HulkIndex
from hulkrna.query_bam import (
    HIT_COLUMNS,
    HIT_SCHEMA,
    _query_fragment,
    query_bam,
    read_hits,
)


# ---------------------------------------------------------------------------
# Helpers: mock pysam reads
# ---------------------------------------------------------------------------

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

    def _get_tag(tag):
        if tag == "NH":
            return nh
        elif tag == "XS":
            return xs
        raise KeyError(tag)

    read.get_tag = _get_tag
    return read


# CIGAR ops
M = 0   # BAM_CMATCH
N = 3   # BAM_CREF_SKIP


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
# Tests: _query_fragment
# ---------------------------------------------------------------------------

class TestQueryFragment:
    """Test _query_fragment with hand-crafted Fragments."""

    def test_single_exon_in_gene(self, mini_index):
        """A single exon block fully inside gene A exon 1 (99-200)."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        hits = _query_fragment(frag, 0, mini_index)
        assert len(hits) > 0
        # All hits should be for frag_id=0
        assert all(h[0] == 0 for h in hits)
        # Should have EXON-type hits
        exon_hits = [h for h in hits if h[5] == IntervalType.EXON]
        assert len(exon_hits) > 0
        # Should include gene A (g_index=0)
        g_indices = {h[8] for h in exon_hits}
        assert 0 in g_indices

    def test_intergenic_region(self, mini_index):
        """Exon block in intergenic region (chr2 has no genes)."""
        frag = Fragment(
            exons=(GenomicInterval("chr2", 100, 200, Strand.POS),),
            introns=(),
        )
        hits = _query_fragment(frag, 5, mini_index)
        assert len(hits) == 1
        # Should be intergenic
        assert hits[0][5] == IntervalType.INTERGENIC
        # frag_id should be 5
        assert hits[0][0] == 5
        # t_index and g_index should be -1 for intergenic
        assert hits[0][7] == -1
        assert hits[0][8] == -1

    def test_spliced_fragment_with_sj_match(self, mini_index):
        """Fragment with an intron that matches an annotated SJ.

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
        hits = _query_fragment(frag, 10, mini_index)

        # Should have exon overlap hits AND SJ hits
        exon_hits = [h for h in hits if h[5] == IntervalType.EXON]
        sj_hits = [h for h in hits if h[5] == IntervalType.SJ]
        assert len(exon_hits) > 0
        assert len(sj_hits) > 0

        # SJ hits should have overlap=0
        for sj in sj_hits:
            assert sj[6] == 0  # overlap column

        # SJ hits should have t_index for t1 (t_index=0)
        sj_t_indices = {h[7] for h in sj_hits}
        assert 0 in sj_t_indices

    def test_spliced_fragment_no_sj_match(self, mini_index):
        """Fragment with an intron that doesn't match any annotated SJ."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 180, Strand.POS),
                GenomicInterval("chr1", 250, 280, Strand.POS),
            ),
            introns=(
                # Novel splice junction (not in the annotation)
                GenomicInterval("chr1", 180, 250, Strand.POS),
            ),
        )
        hits = _query_fragment(frag, 0, mini_index)

        # Should have exon overlap hits but NO SJ hits
        sj_hits = [h for h in hits if h[5] == IntervalType.SJ]
        assert len(sj_hits) == 0

        # Should still have exon overlap hits
        exon_hits = [h for h in hits if h[5] != IntervalType.SJ]
        assert len(exon_hits) > 0

    def test_empty_fragment(self, mini_index):
        """A fragment with no exons or introns produces no hits."""
        frag = Fragment()
        hits = _query_fragment(frag, 0, mini_index)
        assert len(hits) == 0

    def test_chimeric_fragment_multi_ref(self, mini_index):
        """A chimeric fragment with exons on different chromosomes."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 120, 180, Strand.POS),
                GenomicInterval("chr2", 100, 200, Strand.POS),
            ),
            introns=(),
        )
        hits = _query_fragment(frag, 0, mini_index)

        # Hits on chr1 (gene A)
        chr1_hits = [h for h in hits if h[1] == "chr1"]
        assert len(chr1_hits) > 0

        # Hits on chr2 (intergenic)
        chr2_hits = [h for h in hits if h[1] == "chr2"]
        assert len(chr2_hits) == 1
        assert chr2_hits[0][5] == IntervalType.INTERGENIC

    def test_hit_tuple_columns(self, mini_index):
        """Verify hit tuple has exactly 9 elements matching HIT_COLUMNS."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        hits = _query_fragment(frag, 42, mini_index)
        assert len(hits) > 0
        for hit in hits:
            assert len(hit) == len(HIT_COLUMNS)
            # frag_id
            assert hit[0] == 42
            # ref
            assert isinstance(hit[1], str)
            # start, end
            assert isinstance(hit[2], int)
            assert isinstance(hit[3], int)

    def test_overlap_values(self, mini_index):
        """Overlap should equal query length when fully inside an interval."""
        # Exon block 120-180 is fully inside the EXON interval 99-200
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        hits = _query_fragment(frag, 0, mini_index)
        # For the EXON hits, overlap should be 60 (180-120)
        exon_hits = [h for h in hits if h[5] == IntervalType.EXON]
        overlaps = [h[6] for h in exon_hits]
        assert all(o > 0 for o in overlaps)
        # At least one hit should have overlap == 60
        assert 60 in overlaps

    def test_exon_block_spanning_intron(self, mini_index):
        """An exon block that spans both exon and intron regions of the index.

        Block 190-210 spans exon (99-200) and intron (200-299).
        """
        frag = Fragment(
            exons=(GenomicInterval("chr1", 190, 210, Strand.POS),),
            introns=(),
        )
        hits = _query_fragment(frag, 0, mini_index)

        # Should have both EXON and INTRON hits
        types = {h[5] for h in hits}
        assert IntervalType.EXON in types
        assert IntervalType.INTRON in types


# ---------------------------------------------------------------------------
# Tests: query_bam (integration with parse_bam_file)
# ---------------------------------------------------------------------------

class TestQueryBam:
    """Integration tests for query_bam using mock BAM iterators."""

    def _make_pair(self, qname, ref, r1_start, r2_start, r1_len=50, r2_len=50,
                   r2_reverse=True):
        """Create a simple unspliced read pair."""
        r1 = _mock_read(
            query_name=qname, reference_name=ref,
            reference_start=r1_start,
            cigartuples=[(M, r1_len)],
            is_reverse=False, is_read1=True, is_read2=False,
        )
        r2 = _mock_read(
            query_name=qname, reference_name=ref,
            reference_start=r2_start,
            cigartuples=[(M, r2_len)],
            is_reverse=r2_reverse, is_read1=False, is_read2=True,
        )
        return r1, r2

    def test_single_pair(self, mini_index, tmp_path):
        """One read pair → hits with frag_id=0."""
        r1, r2 = self._make_pair("frag1", "chr1", 120, 150)
        out = tmp_path / "hits.feather"

        stats = query_bam([r1, r2], mini_index, out)
        df = read_hits(out)

        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == HIT_COLUMNS
        assert len(df) > 0
        assert (df["frag_id"] == 0).all()
        assert stats["unique"] == 2
        assert stats["n_fragments"] == 1

    def test_multiple_pairs(self, mini_index, tmp_path):
        """Two read pairs → frag_id 0 and 1."""
        r1a, r2a = self._make_pair("fragA", "chr1", 120, 150)
        r1b, r2b = self._make_pair("fragB", "chr2", 100, 200)
        out = tmp_path / "hits.feather"

        stats = query_bam([r1a, r2a, r1b, r2b], mini_index, out)
        df = read_hits(out)

        assert set(df["frag_id"].unique()) == {0, 1}
        assert stats["n_fragments"] == 2
        # fragA on chr1 (genic), fragB on chr2 (intergenic)
        chr2_rows = df[df["ref"] == "chr2"]
        assert (chr2_rows["interval_type"] == IntervalType.INTERGENIC).all()

    def test_empty_bam(self, mini_index, tmp_path):
        """Empty BAM → empty file with correct schema."""
        out = tmp_path / "hits.feather"

        stats = query_bam(iter([]), mini_index, out)
        df = read_hits(out)

        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == HIT_COLUMNS
        assert len(df) == 0
        assert stats["n_fragments"] == 0
        assert stats["n_hits"] == 0

    def test_duplicate_skipped(self, mini_index, tmp_path):
        """Duplicate reads should be skipped by default."""
        r1 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=120, cigartuples=[(M, 50)],
            is_reverse=False, is_read1=True, is_read2=False,
            is_duplicate=True,
        )
        r2 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=150, cigartuples=[(M, 50)],
            is_reverse=True, is_read1=False, is_read2=True,
            is_duplicate=True,
        )
        out = tmp_path / "hits.feather"

        stats = query_bam([r1, r2], mini_index, out, skip_duplicates=True)
        df = read_hits(out)
        assert len(df) == 0

    def test_duplicate_kept(self, mini_index, tmp_path):
        """Duplicates retained when skip_duplicates=False."""
        r1 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=120, cigartuples=[(M, 50)],
            is_reverse=False, is_read1=True, is_read2=False,
            is_duplicate=True,
        )
        r2 = _mock_read(
            query_name="dup", reference_name="chr1",
            reference_start=150, cigartuples=[(M, 50)],
            is_reverse=True, is_read1=False, is_read2=True,
            is_duplicate=True,
        )
        out = tmp_path / "hits.feather"

        stats = query_bam([r1, r2], mini_index, out, skip_duplicates=False)
        df = read_hits(out)
        assert len(df) > 0

    def test_dtypes(self, mini_index, tmp_path):
        """Verify output dtypes match HIT_SCHEMA."""
        r1, r2 = self._make_pair("frag1", "chr1", 120, 150)
        out = tmp_path / "hits.feather"
        query_bam([r1, r2], mini_index, out)
        table = pa.ipc.open_file(pa.memory_map(str(out), "r")).read_all()
        assert table.schema.equals(HIT_SCHEMA)

    def test_chunked_write(self, mini_index, tmp_path):
        """Verify chunked writing produces correct output with tiny chunk_size."""
        pairs = []
        for i in range(10):
            r1, r2 = self._make_pair(f"frag{i}", "chr1", 120, 150)
            pairs.extend([r1, r2])
        out = tmp_path / "hits.feather"

        # Force flushing every 3 rows to exercise chunking
        stats = query_bam(pairs, mini_index, out, chunk_size=3)
        df = read_hits(out)

        assert stats["n_fragments"] == 10
        assert len(df) == stats["n_hits"]
        assert set(df["frag_id"].unique()) == set(range(10))

    def test_tsv_mirror(self, mini_index, tmp_path):
        """write_tsv=True produces a .tsv alongside the .feather."""
        r1, r2 = self._make_pair("frag1", "chr1", 120, 150)
        out = tmp_path / "hits.feather"
        query_bam([r1, r2], mini_index, out, write_tsv=True)
        assert out.exists()
        assert out.with_suffix(".tsv").exists()


# ---------------------------------------------------------------------------
# Tests: write_hits
# ---------------------------------------------------------------------------

class TestReadHits:
    """Test read_hits round-trip."""

    def test_read_hits_columns(self, tmp_path, mini_index):
        """read_hits returns DataFrame with correct columns."""
        frag = Fragment(
            exons=(GenomicInterval("chr1", 120, 180, Strand.POS),),
            introns=(),
        )
        # Write via query_bam with a single mock pair
        r1 = _mock_read(
            query_name="frag1", reference_name="chr1",
            reference_start=120, cigartuples=[(M, 60)],
            is_reverse=False, is_read1=True, is_read2=False,
        )
        r2 = _mock_read(
            query_name="frag1", reference_name="chr1",
            reference_start=140, cigartuples=[(M, 40)],
            is_reverse=True, is_read1=False, is_read2=True,
        )
        out = tmp_path / "hits.feather"
        query_bam([r1, r2], mini_index, out)

        df = read_hits(out)
        assert list(df.columns) == HIT_COLUMNS
        assert len(df) > 0

    def test_read_empty_file(self, tmp_path, mini_index):
        """read_hits on an empty-BAM output returns empty DataFrame."""
        out = tmp_path / "empty.feather"
        query_bam(iter([]), mini_index, out)
        df = read_hits(out)
        assert len(df) == 0
        assert list(df.columns) == HIT_COLUMNS


# ---------------------------------------------------------------------------
# Tests: SJ exact match through _query_fragment
# ---------------------------------------------------------------------------

class TestSJExactMatch:
    """Test SJ hit generation in _query_fragment."""

    def test_sj_overlap_is_zero(self, mini_index):
        """SJ hits should always have overlap=0."""
        # t1 intron: (200, 299) in 0-based half-open
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ),
        )
        hits = _query_fragment(frag, 0, mini_index)
        sj_hits = [h for h in hits if h[5] == IntervalType.SJ]
        assert len(sj_hits) > 0
        for sj in sj_hits:
            assert sj[6] == 0  # overlap must be 0

    def test_sj_hit_has_genomic_coords(self, mini_index):
        """SJ hits should carry the intron's genomic start/end."""
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.POS),
            ),
        )
        hits = _query_fragment(frag, 0, mini_index)
        sj_hits = [h for h in hits if h[5] == IntervalType.SJ]
        for sj in sj_hits:
            assert sj[1] == "chr1"   # ref
            assert sj[2] == 200      # start
            assert sj[3] == 299      # end

    def test_wrong_strand_sj_no_match(self, mini_index):
        """SJ with wrong strand should not match."""
        # t1 intron is on + strand; query with - strand
        frag = Fragment(
            exons=(
                GenomicInterval("chr1", 150, 200, Strand.POS),
                GenomicInterval("chr1", 299, 350, Strand.POS),
            ),
            introns=(
                GenomicInterval("chr1", 200, 299, Strand.NEG),
            ),
        )
        hits = _query_fragment(frag, 0, mini_index)
        sj_hits = [h for h in hits if h[5] == IntervalType.SJ]
        assert len(sj_hits) == 0
