"""Tests for hulkrna.transcript — Transcript parsing from GTF."""

import pytest

from hulkrna.core import Interval, Strand
from hulkrna.transcript import Transcript


class TestTranscriptReadGTF:
    """Test reading transcripts from a GTF file."""

    def test_read_gtf_count(self, mini_gtf_file):
        """Three transcripts (t1, t2, t3) should be parsed."""
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        assert len(transcripts) == 3

    def test_transcript_ids(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t_ids = {t.t_id for t in transcripts}
        assert t_ids == {"t1", "t2", "t3"}

    def test_gene_ids(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        g_ids = {t.g_id for t in transcripts}
        assert g_ids == {"g1", "g2"}

    def test_coordinate_conversion(self, mini_gtf_file):
        """GTF 1-based inclusive → 0-based half-open."""
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        # GTF exon 100-200 → 0-based (99, 200)
        assert t1.exons[0] == Interval(99, 200)
        # GTF exon 300-400 → 0-based (299, 400)
        assert t1.exons[1] == Interval(299, 400)

    def test_exon_count(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        t2 = next(t for t in transcripts if t.t_id == "t2")
        t3 = next(t for t in transcripts if t.t_id == "t3")
        assert len(t1.exons) == 3
        assert len(t2.exons) == 2
        assert len(t3.exons) == 2

    def test_strand(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        t3 = next(t for t in transcripts if t.t_id == "t3")
        assert t1.strand == Strand.POS
        assert t3.strand == Strand.NEG

    def test_gene_metadata(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        assert t1.g_name == "GeneA"
        assert t1.g_type == "protein_coding"
        assert t1.is_basic is True

    def test_exons_sorted(self, mini_gtf_file):
        """Exons within each transcript should be sorted by start."""
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        for t in transcripts:
            for i in range(len(t.exons) - 1):
                assert t.exons[i].start <= t.exons[i + 1].start

    def test_transcripts_sorted(self, mini_gtf_file):
        """Transcripts should be sorted by (ref, start, end, strand)."""
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        keys = [(t.ref, t.start, t.end, t.strand) for t in transcripts]
        assert keys == sorted(keys)

    def test_length_computed(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        # 3 exons: (99,200)=101bp + (299,400)=101bp + (499,600)=101bp = 303bp
        assert t1.length == 303


class TestTranscriptIntrons:
    """Test intron (splice junction) extraction."""

    def test_introns_from_three_exon_transcript(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t1 = next(t for t in transcripts if t.t_id == "t1")
        introns = list(t1.introns())
        # Between exon1(99,200) and exon2(299,400): intron (200, 299)
        # Between exon2(299,400) and exon3(499,600): intron (400, 499)
        assert introns == [(200, 299), (400, 499)]

    def test_introns_from_two_exon_transcript(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        t2 = next(t for t in transcripts if t.t_id == "t2")
        introns = list(t2.introns())
        # Between exon1(99,200) and exon2(499,600): intron (200, 499)
        assert introns == [(200, 499)]


class TestTranscriptToDict:
    def test_to_dict_keys(self, mini_gtf_file):
        transcripts = Transcript.read_gtf(str(mini_gtf_file))
        d = transcripts[0].to_dict()
        expected_keys = {
            "ref", "start", "end", "strand", "length",
            "t_id", "g_id", "t_index", "g_index",
            "g_name", "g_type",
            "is_basic", "is_mane", "is_ccds",
            "seq", "abundance",
        }
        assert set(d.keys()) == expected_keys
