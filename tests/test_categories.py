"""Tests for hulkrna.categories — count classification enums."""

import pytest

from hulkrna.categories import (
    ANTISENSE_COLS,
    CountCategory,
    CountCol,
    NUM_CATEGORIES,
    NUM_COUNT_COLS,
    SPLICED_COLS,
)


# =====================================================================
# CountCategory
# =====================================================================


class TestCountCategory:
    def test_values(self):
        assert int(CountCategory.INTRON) == 0
        assert int(CountCategory.UNSPLICED) == 1
        assert int(CountCategory.SPLICED_UNANNOT) == 2
        assert int(CountCategory.SPLICED_ANNOT) == 3

    def test_length(self):
        assert len(CountCategory) == 4
        assert NUM_CATEGORIES == 4


# =====================================================================
# CountCol — 4 categories × 2 strands = 8 columns
# =====================================================================


class TestCountCol:
    def test_num_count_cols(self):
        assert NUM_COUNT_COLS == 8

    def test_enum_values(self):
        assert int(CountCol.INTRON_SENSE) == 0
        assert int(CountCol.INTRON_ANTISENSE) == 1
        assert int(CountCol.UNSPLICED_SENSE) == 2
        assert int(CountCol.UNSPLICED_ANTISENSE) == 3
        assert int(CountCol.SPLICED_UNANNOT_SENSE) == 4
        assert int(CountCol.SPLICED_UNANNOT_ANTISENSE) == 5
        assert int(CountCol.SPLICED_ANNOT_SENSE) == 6
        assert int(CountCol.SPLICED_ANNOT_ANTISENSE) == 7

    def test_from_category(self):
        """from_category(category, is_antisense) == category * 2 + is_antisense."""
        for cat in CountCategory:
            for anti in (False, True):
                expected = int(cat) * 2 + int(anti)
                assert CountCol.from_category(cat, anti) == expected

    def test_is_antisense_property(self):
        for c in CountCol:
            assert c.is_antisense == bool(c.value % 2)

    def test_category_property(self):
        for c in CountCol:
            assert c.category == CountCategory(c.value // 2)

    def test_roundtrip(self):
        """from_category(col.category, col.is_antisense) == col."""
        for c in CountCol:
            assert CountCol.from_category(c.category, c.is_antisense) == c


# =====================================================================
# Pre-computed column subsets
# =====================================================================


class TestColumnSubsets:
    def test_antisense_cols(self):
        assert ANTISENSE_COLS == (1, 3, 5, 7)

    def test_spliced_cols(self):
        assert SPLICED_COLS == (4, 5, 6, 7)

    def test_antisense_cols_match_enum(self):
        expected = tuple(c.value for c in CountCol if c.is_antisense)
        assert ANTISENSE_COLS == expected

    def test_spliced_cols_match_enum(self):
        expected = tuple(
            c.value for c in CountCol
            if c.category in (CountCategory.SPLICED_ANNOT, CountCategory.SPLICED_UNANNOT)
        )
        assert SPLICED_COLS == expected
