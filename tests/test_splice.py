"""Tests for hulkrna.splice — splice classification enums."""

import pytest

from hulkrna.splice import (
    ANTISENSE_COLS,
    SpliceType,
    SpliceStrandCol,
    NUM_SPLICE_TYPES,
    NUM_SPLICE_STRAND_COLS,
    SPLICED_COLS,
)


# =====================================================================
# SpliceType
# =====================================================================


class TestSpliceType:
    def test_values(self):
        assert int(SpliceType.UNSPLICED) == 0
        assert int(SpliceType.SPLICED_UNANNOT) == 1
        assert int(SpliceType.SPLICED_ANNOT) == 2

    def test_length(self):
        assert len(SpliceType) == 3
        assert NUM_SPLICE_TYPES == 3


# =====================================================================
# SpliceStrandCol — 3 categories × 2 strands = 6 columns
# =====================================================================


class TestSpliceStrandCol:
    def test_num_splice_strand_cols(self):
        assert NUM_SPLICE_STRAND_COLS == 6

    def test_enum_values(self):
        assert int(SpliceStrandCol.UNSPLICED_SENSE) == 0
        assert int(SpliceStrandCol.UNSPLICED_ANTISENSE) == 1
        assert int(SpliceStrandCol.SPLICED_UNANNOT_SENSE) == 2
        assert int(SpliceStrandCol.SPLICED_UNANNOT_ANTISENSE) == 3
        assert int(SpliceStrandCol.SPLICED_ANNOT_SENSE) == 4
        assert int(SpliceStrandCol.SPLICED_ANNOT_ANTISENSE) == 5

    def test_from_category(self):
        """from_category(category, is_antisense) == category * 2 + is_antisense."""
        for cat in SpliceType:
            for anti in (False, True):
                expected = int(cat) * 2 + int(anti)
                assert SpliceStrandCol.from_category(cat, anti) == expected

    def test_is_antisense_property(self):
        for c in SpliceStrandCol:
            assert c.is_antisense == bool(c.value % 2)

    def test_category_property(self):
        for c in SpliceStrandCol:
            assert c.category == SpliceType(c.value // 2)

    def test_roundtrip(self):
        """from_category(col.category, col.is_antisense) == col."""
        for c in SpliceStrandCol:
            assert SpliceStrandCol.from_category(c.category, c.is_antisense) == c


# =====================================================================
# Pre-computed column subsets
# =====================================================================


class TestColumnSubsets:
    def test_antisense_cols(self):
        assert ANTISENSE_COLS == (1, 3, 5)

    def test_spliced_cols(self):
        assert SPLICED_COLS == (2, 3, 4, 5)

    def test_antisense_cols_match_enum(self):
        expected = tuple(c.value for c in SpliceStrandCol if c.is_antisense)
        assert ANTISENSE_COLS == expected

    def test_spliced_cols_match_enum(self):
        expected = tuple(
            c.value for c in SpliceStrandCol
            if c.category in (SpliceType.SPLICED_ANNOT, SpliceType.SPLICED_UNANNOT)
        )
        assert SPLICED_COLS == expected
