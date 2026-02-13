"""
hulkrna.categories — Count classification enums for the counting pipeline.

These enums encode how a fragment overlaps gene annotations (category),
the strand relationship between the fragment and the gene (strand),
and the flat cross-product used as array column indices (type).
"""

from enum import IntEnum


class CountCategory(IntEnum):
    """Category of a fragment's overlap with gene annotations."""
    INTRON = 0
    UNSPLICED = 1
    SPLICED_UNANNOT = 2
    SPLICED_ANNOT = 3


class CountStrand(IntEnum):
    """Strand classification for counting."""
    AMBIGUOUS = 0
    SENSE = 1
    ANTISENSE = 2


class CountType(IntEnum):
    """Flat encoding of CountCategory × CountStrand (4 × 3 = 12 values).

    ``index = category * 3 + strand_idx``
    """
    INTRON_AMBIGUOUS = 0
    INTRON_SENSE = 1
    INTRON_ANTISENSE = 2
    UNSPLICED_AMBIGUOUS = 3
    UNSPLICED_SENSE = 4
    UNSPLICED_ANTISENSE = 5
    SPLICED_UNANNOT_AMBIGUOUS = 6
    SPLICED_UNANNOT_SENSE = 7
    SPLICED_UNANNOT_ANTISENSE = 8
    SPLICED_ANNOT_AMBIGUOUS = 9
    SPLICED_ANNOT_SENSE = 10
    SPLICED_ANNOT_ANTISENSE = 11

    def to_column_name(self) -> str:
        """Return lower-cased name, e.g. ``'intron_sense'``."""
        return self.name.lower()

    @classmethod
    def columns(cls):
        """Yield all column names in enum order."""
        yield from (e.to_column_name() for e in cls)


NUM_COUNT_TYPES = len(CountType)


class AssignmentSource(IntEnum):
    """Provenance of a count assignment.

    ``UNIQUE`` — deterministic assignment from unambiguous fragments
    (1 gene, 1 transcript, NH = 1).

    ``EM`` — stochastic assignment from EM-resolved ambiguous fragments
    (isoform-ambiguous, gene-ambiguous, and multimapped, resolved
    simultaneously via Expectation-Maximization).
    """

    UNIQUE = 0
    EM = 1


NUM_ASSIGNMENT_SOURCES = len(AssignmentSource)
