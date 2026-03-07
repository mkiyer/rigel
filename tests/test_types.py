"""Tests for rigel.types — foundational data types."""

import pytest

from rigel.types import (
    GenomicInterval,
    Interval,
    IntervalType,
    MergeOutcome,
    AnnotatedInterval,
    Strand,
)


# =====================================================================
# Strand
# =====================================================================


class TestStrand:
    """Tests for the Strand IntEnum."""

    # -- values -----------------------------------------------------------

    def test_values(self):
        assert int(Strand.NONE) == 0
        assert int(Strand.POS) == 1
        assert int(Strand.NEG) == 2
        assert int(Strand.AMBIGUOUS) == 3

    def test_bitwise_or_pos_neg_gives_ambiguous(self):
        assert Strand.POS | Strand.NEG == Strand.AMBIGUOUS

    def test_bitwise_or_with_none(self):
        assert Strand.NONE | Strand.POS == Strand.POS
        assert Strand.NONE | Strand.NEG == Strand.NEG
        assert Strand.NONE | Strand.NONE == Strand.NONE

    def test_bitwise_or_idempotent(self):
        assert Strand.POS | Strand.POS == Strand.POS
        assert Strand.NEG | Strand.NEG == Strand.NEG

    # -- from_str / to_str ------------------------------------------------

    @pytest.mark.parametrize(
        "char, expected",
        [(".", Strand.NONE), ("+", Strand.POS), ("-", Strand.NEG), ("?", Strand.AMBIGUOUS)],
    )
    def test_from_str(self, char, expected):
        assert Strand.from_str(char) == expected

    def test_from_str_invalid_raises(self):
        with pytest.raises(ValueError, match="Invalid strand"):
            Strand.from_str("x")

    @pytest.mark.parametrize(
        "strand, expected",
        [(Strand.NONE, "."), (Strand.POS, "+"), (Strand.NEG, "-"), (Strand.AMBIGUOUS, "?")],
    )
    def test_to_str(self, strand, expected):
        assert strand.to_str() == expected

    def test_from_str_to_str_roundtrip(self):
        for ch in (".", "+", "-", "?"):
            assert Strand.from_str(ch).to_str() == ch

    # -- from_is_reverse ---------------------------------------------------

    def test_from_is_reverse_false_gives_pos(self):
        assert Strand.from_is_reverse(False) == Strand.POS

    def test_from_is_reverse_true_gives_neg(self):
        assert Strand.from_is_reverse(True) == Strand.NEG

    # -- opposite ----------------------------------------------------------

    def test_opposite_pos(self):
        assert Strand.POS.opposite() == Strand.NEG

    def test_opposite_neg(self):
        assert Strand.NEG.opposite() == Strand.POS

    def test_opposite_none_unchanged(self):
        assert Strand.NONE.opposite() == Strand.NONE

    def test_opposite_ambiguous_unchanged(self):
        assert Strand.AMBIGUOUS.opposite() == Strand.AMBIGUOUS

    def test_opposite_double_reversal(self):
        assert Strand.POS.opposite().opposite() == Strand.POS
        assert Strand.NEG.opposite().opposite() == Strand.NEG


# =====================================================================
# Interval
# =====================================================================


class TestInterval:
    def test_creation(self):
        iv = Interval(100, 200)
        assert iv.start == 100
        assert iv.end == 200

    def test_tuple_unpacking(self):
        start, end = Interval(10, 20)
        assert start == 10
        assert end == 20

    def test_comparison(self):
        assert Interval(10, 20) < Interval(10, 30)
        assert Interval(10, 20) < Interval(20, 30)
        assert Interval(10, 20) == Interval(10, 20)


# =====================================================================
# GenomicInterval
# =====================================================================


class TestGenomicInterval:
    def test_default_strand_is_none(self):
        gi = GenomicInterval("chr1", 100, 200)
        assert gi.strand == Strand.NONE

    def test_with_strand(self):
        gi = GenomicInterval("chr1", 100, 200, Strand.POS)
        assert gi.ref == "chr1"
        assert gi.start == 100
        assert gi.end == 200
        assert gi.strand == Strand.POS

    def test_tuple_unpacking(self):
        ref, start, end, strand = GenomicInterval("chr1", 0, 50, Strand.NEG)
        assert ref == "chr1"
        assert start == 0
        assert end == 50
        assert strand == Strand.NEG


# =====================================================================
# IntervalType
# =====================================================================


class TestIntervalType:
    def test_values(self):
        assert int(IntervalType.EXON) == 0
        assert int(IntervalType.TRANSCRIPT) == 1
        assert int(IntervalType.INTERGENIC) == 2
        assert int(IntervalType.SJ) == 3
        assert int(IntervalType.SJ_UNANNOT) == 4


# =====================================================================
# AnnotatedInterval
# =====================================================================


class TestAnnotatedInterval:
    def test_defaults(self):
        ri = AnnotatedInterval("chr1", 100, 200)
        assert ri.strand == Strand.NONE
        assert ri.interval_type == IntervalType.INTERGENIC
        assert ri.t_index == -1

    def test_with_all_fields(self):
        ri = AnnotatedInterval(
            "chr1", 100, 200,
            strand=Strand.POS,
            interval_type=IntervalType.EXON,
            t_index=3,
        )
        assert ri.ref == "chr1"
        assert ri.start == 100
        assert ri.end == 200
        assert ri.strand == Strand.POS
        assert ri.interval_type == IntervalType.EXON
        assert ri.t_index == 3


# =====================================================================
# MergeOutcome
# =====================================================================


class TestMergeOutcome:
    def test_values(self):
        assert int(MergeOutcome.INTERSECTION) == 0
        assert int(MergeOutcome.INTERSECTION_NONEMPTY) == 1
        assert int(MergeOutcome.UNION) == 2
        assert int(MergeOutcome.EMPTY) == 3


class TestCppConstantParity:
    """Verify C++ scoring constants match Python-computed values."""

    def test_log_half(self):
        import math
        from rigel._scoring_impl import LOG_HALF
        assert LOG_HALF == math.log(0.5)

    def test_tail_decay_lp(self):
        import math
        from rigel._scoring_impl import TAIL_DECAY_LP
        assert TAIL_DECAY_LP == math.log(0.99)


