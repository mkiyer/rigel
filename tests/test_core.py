"""Tests for hulkrna.core — shared types."""

import pytest

from hulkrna.core import Strand, Interval, GenomicInterval, IntervalType, RefInterval


# ---------------------------------------------------------------------------
# Strand
# ---------------------------------------------------------------------------

class TestStrand:
    """Tests for the Strand IntEnum."""

    def test_values(self):
        assert Strand.NONE == 0
        assert Strand.POS == 1
        assert Strand.NEG == 2
        assert Strand.AMBIGUOUS == 3

    def test_bitwise_or_produces_ambiguous(self):
        assert Strand.POS | Strand.NEG == Strand.AMBIGUOUS

    def test_bitwise_or_accumulation(self):
        combined = Strand.NONE
        combined |= Strand.POS
        assert combined == Strand.POS
        combined |= Strand.NEG
        assert combined == Strand.AMBIGUOUS

    def test_from_str_valid(self):
        assert Strand.from_str(".") == Strand.NONE
        assert Strand.from_str("+") == Strand.POS
        assert Strand.from_str("-") == Strand.NEG
        assert Strand.from_str("?") == Strand.AMBIGUOUS

    def test_from_str_invalid(self):
        with pytest.raises(ValueError, match="Invalid strand string"):
            Strand.from_str("X")

    def test_to_str_roundtrip(self):
        for s in (Strand.NONE, Strand.POS, Strand.NEG, Strand.AMBIGUOUS):
            assert Strand.from_str(s.to_str()) == s

    def test_from_is_reverse(self):
        assert Strand.from_is_reverse(False) == Strand.POS
        assert Strand.from_is_reverse(True) == Strand.NEG

    def test_opposite(self):
        assert Strand.POS.opposite() == Strand.NEG
        assert Strand.NEG.opposite() == Strand.POS
        assert Strand.NONE.opposite() == Strand.NONE
        assert Strand.AMBIGUOUS.opposite() == Strand.AMBIGUOUS


# ---------------------------------------------------------------------------
# Interval (formerly Exon)
# ---------------------------------------------------------------------------

class TestInterval:
    def test_creation(self):
        e = Interval(100, 200)
        assert e.start == 100
        assert e.end == 200

    def test_sorting(self):
        intervals = [Interval(500, 600), Interval(100, 200), Interval(300, 400)]
        assert sorted(intervals) == [
            Interval(100, 200), Interval(300, 400), Interval(500, 600)
        ]


# ---------------------------------------------------------------------------
# GenomicInterval
# ---------------------------------------------------------------------------

class TestGenomicInterval:
    def test_creation(self):
        gi = GenomicInterval("chr1", 100, 200, Strand.POS)
        assert gi.ref == "chr1"
        assert gi.start == 100
        assert gi.end == 200
        assert gi.strand == Strand.POS

    def test_default_strand(self):
        gi = GenomicInterval("chr1", 0, 1000)
        assert gi.strand == Strand.NONE


# ---------------------------------------------------------------------------
# IntervalType
# ---------------------------------------------------------------------------

class TestIntervalType:
    def test_values(self):
        assert IntervalType.EXON == 0
        assert IntervalType.INTRON == 1
        assert IntervalType.INTERGENIC == 2
        assert IntervalType.SJ == 3


# ---------------------------------------------------------------------------
# RefInterval (formerly Interval)
# ---------------------------------------------------------------------------

class TestRefInterval:
    def test_defaults(self):
        iv = RefInterval("chr1", 0, 1000)
        assert iv.strand == Strand.NONE
        assert iv.interval_type == IntervalType.INTERGENIC
        assert iv.t_index == -1
        assert iv.g_index == -1

    def test_exon_interval(self):
        iv = RefInterval("chr1", 100, 200, Strand.POS, IntervalType.EXON, 0, 0)
        assert iv.interval_type == IntervalType.EXON
        assert iv.t_index == 0
        assert iv.g_index == 0

    def test_splice_junction(self):
        """RefInterval can represent a splice junction via IntervalType.SJ."""
        sj = RefInterval("chr1", 200, 300, Strand.POS, IntervalType.SJ, 0, 0)
        assert sj.ref == "chr1"
        assert sj.start == 200
        assert sj.end == 300
        assert sj.strand == Strand.POS
        assert sj.interval_type == IntervalType.SJ
