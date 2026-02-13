"""Tests for hulkrna.categories — count classification enums."""

import pytest

from hulkrna.categories import (
    CountCategory,
    CountStrand,
    CountType,
    NUM_COUNT_TYPES,
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


# =====================================================================
# CountStrand
# =====================================================================


class TestCountStrand:
    def test_values(self):
        assert int(CountStrand.AMBIGUOUS) == 0
        assert int(CountStrand.SENSE) == 1
        assert int(CountStrand.ANTISENSE) == 2

    def test_length(self):
        assert len(CountStrand) == 3


# =====================================================================
# CountType
# =====================================================================


class TestCountType:
    def test_num_count_types(self):
        assert NUM_COUNT_TYPES == 12
        assert len(CountType) == 12

    def test_encoding_formula(self):
        """Each CountType value should equal category * 3 + strand."""
        for cat in CountCategory:
            for strand in CountStrand:
                expected = int(cat) * 3 + int(strand)
                name = f"{cat.name}_{strand.name}"
                ct = CountType[name]
                assert int(ct) == expected, f"{name}: {int(ct)} != {expected}"

    def test_to_column_name(self):
        assert CountType.INTRON_SENSE.to_column_name() == "intron_sense"
        assert CountType.SPLICED_ANNOT_AMBIGUOUS.to_column_name() == "spliced_annot_ambiguous"

    def test_columns_yields_all(self):
        cols = list(CountType.columns())
        assert len(cols) == 12
        assert cols[0] == "intron_ambiguous"
        assert cols[-1] == "spliced_annot_antisense"

    def test_columns_all_lowercase(self):
        for col in CountType.columns():
            assert col == col.lower()
