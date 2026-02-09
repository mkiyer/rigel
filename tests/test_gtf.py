"""Tests for hulkrna.gtf — GTF parsing."""

import pytest

from hulkrna.gtf import GTF


# ---------------------------------------------------------------------------
# GTF.from_str — single-line parsing
# ---------------------------------------------------------------------------

class TestGTFFromStr:
    """Test parsing individual GTF lines."""

    EXON_LINE = (
        'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\t'
        'gene_id "G1"; transcript_id "T1"; gene_name "MyGene"; '
        'gene_type "protein_coding"; tag "basic"; tag "MANE_Select";'
    )

    def test_basic_fields(self):
        row = GTF.from_str(self.EXON_LINE)
        assert row.seqname == "chr1"
        assert row.source == "HAVANA"
        assert row.feature == "exon"
        assert row.strand == "+"
        assert row.score is None

    def test_coordinate_conversion(self):
        """GTF 1-based inclusive → 0-based half-open."""
        row = GTF.from_str(self.EXON_LINE)
        # GTF start=100 → 0-based = 99
        assert row.start == 99
        # GTF end=200 → stays 200 (half-open)
        assert row.end == 200

    def test_attributes(self):
        row = GTF.from_str(self.EXON_LINE)
        assert row.attrs["gene_id"] == "G1"
        assert row.attrs["transcript_id"] == "T1"
        assert row.attrs["gene_name"] == "MyGene"
        assert row.attrs["gene_type"] == "protein_coding"

    def test_tags(self):
        row = GTF.from_str(self.EXON_LINE)
        assert "basic" in row.tags
        assert "MANE_Select" in row.tags

    def test_score_parsed(self):
        line = 'chr1\ttest\texon\t10\t20\t0.5\t+\t.\tgene_id "G1"; transcript_id "T1";'
        row = GTF.from_str(line)
        assert row.score == 0.5

    def test_malformed_line_raises(self):
        with pytest.raises(ValueError, match="Malformed GTF"):
            GTF.from_str("chr1\tbad_line")


# ---------------------------------------------------------------------------
# GTF.parse — multi-line iteration
# ---------------------------------------------------------------------------

class TestGTFParse:
    def test_skips_comments_and_blanks(self):
        lines = [
            "# this is a comment",
            "",
            'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        ]
        results = list(GTF.parse(lines))
        assert len(results) == 1
        assert results[0].attrs["gene_id"] == "G1"

    def test_multiple_lines(self):
        lines = [
            'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
            'chr1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        ]
        results = list(GTF.parse(lines))
        assert len(results) == 2


# ---------------------------------------------------------------------------
# GTF.parse_file
# ---------------------------------------------------------------------------

class TestGTFParseFile:
    def test_parse_file(self, mini_gtf_file):
        results = list(GTF.parse_file(mini_gtf_file))
        # 7 exon lines in MINI_GTF
        assert len(results) == 7

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            list(GTF.parse_file(tmp_path / "nonexistent.gtf"))


# ---------------------------------------------------------------------------
# GTF.__str__ roundtrip
# ---------------------------------------------------------------------------

class TestGTFStr:
    def test_roundtrip_coordinates(self):
        """Parsing then stringifying should preserve 1-based coordinates."""
        original = 'chr1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";'
        row = GTF.from_str(original)
        output = str(row)
        # Re-parse the output
        row2 = GTF.from_str(output)
        assert row2.start == row.start
        assert row2.end == row.end
        assert row2.attrs["gene_id"] == "G1"
