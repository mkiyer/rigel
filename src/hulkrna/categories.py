"""
hulkrna.categories — Count classification enums for the counting pipeline.

Internal count arrays use 8 columns: 4 categories × 2 strands
(sense/antisense).  Strand is an internal signal used for gDNA
estimation (Iceberg model) and is not exposed in user-facing output.
"""

from enum import IntEnum


class CountCategory(IntEnum):
    """Category of a fragment's overlap with gene annotations."""
    INTRON = 0
    UNSPLICED = 1
    SPLICED_UNANNOT = 2
    SPLICED_ANNOT = 3


NUM_CATEGORIES = len(CountCategory)


class CountCol(IntEnum):
    """Internal column index for count arrays.

    Layout: 4 categories × 2 strands (sense/antisense) = 8 columns.
    Even indices are sense, odd indices are antisense.
    """
    INTRON_SENSE = 0
    INTRON_ANTISENSE = 1
    UNSPLICED_SENSE = 2
    UNSPLICED_ANTISENSE = 3
    SPLICED_UNANNOT_SENSE = 4
    SPLICED_UNANNOT_ANTISENSE = 5
    SPLICED_ANNOT_SENSE = 6
    SPLICED_ANNOT_ANTISENSE = 7

    @classmethod
    def from_category(cls, category: int, is_antisense: bool) -> "CountCol":
        """Look up column from category and strand."""
        return cls(int(category) * 2 + int(is_antisense))

    @property
    def is_antisense(self) -> bool:
        return bool(self.value % 2)

    @property
    def category(self) -> CountCategory:
        return CountCategory(self.value // 2)


NUM_COUNT_COLS = len(CountCol)

# Pre-computed column subsets (tuples for immutability and indexing).
ANTISENSE_COLS = tuple(c.value for c in CountCol if c.is_antisense)
SPLICED_COLS = tuple(
    c.value for c in CountCol
    if c.category in (CountCategory.SPLICED_ANNOT, CountCategory.SPLICED_UNANNOT)
)
