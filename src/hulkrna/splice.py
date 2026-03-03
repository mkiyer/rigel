"""
hulkrna.splice — Splice classification enums for the counting pipeline.

Internal count arrays use 6 columns: 3 splice types × 2 strands
(sense/antisense).  Strand is an internal signal used for gDNA
estimation and is not exposed in user-facing output.

Fragment splice types are based on splice junction status only:
- UNSPLICED: no splice junctions detected
- SPLICED_UNANNOT: splice junctions present but not in the reference
- SPLICED_ANNOT: splice junctions exactly matching the reference

Intronic vs exonic overlap is captured separately in the per-candidate
overlap profile (n_exon_bp, n_intron_bp), not in the splice type.
"""

from enum import IntEnum


class SpliceType(IntEnum):
    """Category of a fragment based on splice junction status."""
    UNSPLICED = 0
    SPLICED_UNANNOT = 1
    SPLICED_ANNOT = 2


NUM_SPLICE_TYPES = len(SpliceType)

#: Pre-computed int constants for hot-path comparisons (avoid enum overhead).
SPLICE_UNSPLICED: int = int(SpliceType.UNSPLICED)         # 0
SPLICE_UNANNOT: int = int(SpliceType.SPLICED_UNANNOT)     # 1
SPLICE_ANNOT: int = int(SpliceType.SPLICED_ANNOT)         # 2


class SpliceStrandCol(IntEnum):
    """Internal column index for count arrays.

    Layout: 3 categories × 2 strands (sense/antisense) = 6 columns.
    Even indices are sense, odd indices are antisense.
    """
    UNSPLICED_SENSE = 0
    UNSPLICED_ANTISENSE = 1
    SPLICED_UNANNOT_SENSE = 2
    SPLICED_UNANNOT_ANTISENSE = 3
    SPLICED_ANNOT_SENSE = 4
    SPLICED_ANNOT_ANTISENSE = 5

    @classmethod
    def from_category(cls, category: int, is_antisense: bool) -> "SpliceStrandCol":
        """Look up column from category and strand."""
        return cls(int(category) * 2 + int(is_antisense))

    @property
    def is_antisense(self) -> bool:
        return bool(self.value % 2)

    @property
    def category(self) -> SpliceType:
        return SpliceType(self.value // 2)


NUM_SPLICE_STRAND_COLS = len(SpliceStrandCol)

# Pre-computed column subsets (tuples for immutability and indexing).
ANTISENSE_COLS = tuple(c.value for c in SpliceStrandCol if c.is_antisense)
SPLICED_COLS = tuple(
    c.value for c in SpliceStrandCol
    if c.category in (SpliceType.SPLICED_ANNOT, SpliceType.SPLICED_UNANNOT)
)
