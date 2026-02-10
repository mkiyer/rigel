"""
hulkrna.core — Shared data types used across multiple hulkrna modules.

This module defines the foundational types for genomic coordinates, strand
orientation, interval representation, and set-merging primitives used
throughout the pipeline.

Coordinate convention: all coordinates are 0-based, half-open (BED style).
"""

from dataclasses import dataclass
from enum import IntEnum
from typing import NamedTuple


# ---------------------------------------------------------------------------
# Strand
# ---------------------------------------------------------------------------

class Strand(IntEnum):
    """Genomic strand with bitwise OR semantics.

    Values are designed so that ``POS | NEG == AMBIGUOUS``, enabling
    efficient accumulation when a fragment spans both strands::

        combined = Strand.NONE
        combined |= Strand.POS   # → POS
        combined |= Strand.NEG   # → AMBIGUOUS
    """
    NONE = 0
    POS = 1
    NEG = 2
    AMBIGUOUS = 3  # POS | NEG

    # -- string conversion ---------------------------------------------------

    @classmethod
    def from_str(cls, s: str) -> "Strand":
        """Convert a single-character strand string to a Strand value.

        Accepted values: ``'.'``, ``'+'``, ``'-'``, ``'?'``.
        """
        _map = {".": 0, "+": 1, "-": 2, "?": 3}
        try:
            return cls(_map[s])
        except KeyError:
            raise ValueError(f"Invalid strand string: {s!r}") from None

    def to_str(self) -> str:
        """Return the single-character string representation."""
        return (".", "+", "-", "?")[self.value]

    # -- pysam helpers --------------------------------------------------------

    @classmethod
    def from_is_reverse(cls, is_reverse: bool) -> "Strand":
        """Convert a pysam ``is_reverse`` flag to a Strand."""
        return cls.NEG if is_reverse else cls.POS

    # -- arithmetic -----------------------------------------------------------
    def opposite(self) -> "Strand":
        """
        Return the opposite strand (POS <-> NEG). 
        NONE and AMBIGUOUS are returned unchanged.
        """
        # Swap bit 0 and bit 1
        return Strand(((self & 1) << 1) | (self >> 1))


# ---------------------------------------------------------------------------
# Interval
# ---------------------------------------------------------------------------

class Interval(NamedTuple):
    """A simple 0-based half-open interval (start, end).

    Used by ``Transcript`` to store exon coordinates. For intervals
    positioned on a specific chromosome/strand, use ``GenomicInterval``.
    """
    start: int
    end: int


# ---------------------------------------------------------------------------
# GenomicInterval
# ---------------------------------------------------------------------------

class GenomicInterval(NamedTuple):
    """An interval located on a specific chromosome and strand.

    Used by ``Fragment`` for aligned exons and splice junctions (introns),
    and anywhere a positioned interval is needed without annotation metadata.
    """
    ref: str
    start: int
    end: int
    strand: int = Strand.NONE


# ---------------------------------------------------------------------------
# IntervalType
# ---------------------------------------------------------------------------

class IntervalType(IntEnum):
    """Classification of a genomic interval relative to gene annotations.

    ``SJ`` is an annotated splice junction that exactly matches a known
    intron in the transcript reference.  ``SJ_UNANNOT`` is a splice
    junction observed in the CIGAR (N-operation) but not matching any
    annotated intron — these are recorded with ``t_index = -1`` and
    ``g_index = -1``.
    """
    EXON = 0
    INTRON = 1
    INTERGENIC = 2
    SJ = 3
    SJ_UNANNOT = 4


# ---------------------------------------------------------------------------
# RefInterval
# ---------------------------------------------------------------------------

class RefInterval(NamedTuple):
    """A reference-annotated genomic interval with index metadata.

    Used for both tiling intervals (EXON/INTRON/INTERGENIC in the cgranges
    overlap index) and splice junctions (SJ in the exact-match lookup).
    Replaces the former separate ``Interval`` and ``SpliceJunction`` types.
    """
    ref: str
    start: int
    end: int
    strand: int = Strand.NONE
    interval_type: int = IntervalType.INTERGENIC
    t_index: int = -1
    g_index: int = -1


# ---------------------------------------------------------------------------
# Merge criteria and result types
# ---------------------------------------------------------------------------

class MergeCriteria(IntEnum):
    """Which relaxation level succeeded in progressive set merging.

    During fragment resolution, transcript/gene index sets from
    individual exon blocks or splice junctions are merged via
    progressive relaxation:

    0. INTERSECTION — intersection of *all* sets (most specific)
    1. INTERSECTION_NONEMPTY — intersection of non-empty sets only
    2. UNION — union of all sets (most sensitive)
    3. EMPTY — no sets to merge (no hits of this type)
    """
    INTERSECTION = 0
    INTERSECTION_NONEMPTY = 1
    UNION = 2
    EMPTY = 3


@dataclass(frozen=True, slots=True)
class MergeResult:
    """Result of progressive set merging with criteria tracking.

    Wraps the merged transcript/gene index frozensets together with
    which relaxation level produced the result.
    """
    t_inds: frozenset[int]
    g_inds: frozenset[int]
    criteria: MergeCriteria

    @property
    def is_unique_gene(self) -> bool:
        """True if exactly one gene index in the merged result."""
        return len(self.g_inds) == 1

    @property
    def is_unique_transcript(self) -> bool:
        """True if exactly one transcript index in the merged result."""
        return len(self.t_inds) == 1

    @property
    def is_empty(self) -> bool:
        """True if the merge produced no transcript or gene indices."""
        return not self.t_inds and not self.g_inds


#: Sentinel for "no hits at all" — avoids creating new objects.
EMPTY_MERGE = MergeResult(frozenset(), frozenset(), MergeCriteria.EMPTY)


def merge_sets_with_criteria(
    t_sets: list[frozenset[int]],
    g_sets: list[frozenset[int]],
) -> MergeResult:
    """Merge transcript/gene index sets with progressive relaxation.

    Strategy
    --------
    1. **Intersection** of all sets (most specific).
    2. **Intersection of non-empty sets** — skip empty sets that would
       collapse everything to the empty set.
    3. **Union** of all sets (most sensitive).

    Parameters
    ----------
    t_sets : list of frozenset[int]
        Per-segment transcript index sets.
    g_sets : list of frozenset[int]
        Per-segment gene index sets.

    Returns
    -------
    MergeResult
    """
    if not t_sets and not g_sets:
        return EMPTY_MERGE

    # 1. Intersection of all sets (most specific)
    t_inds = frozenset.intersection(*t_sets) if t_sets else frozenset()
    g_inds = frozenset.intersection(*g_sets) if g_sets else frozenset()
    if t_inds or g_inds:
        return MergeResult(t_inds, g_inds, MergeCriteria.INTERSECTION)

    # 2. Intersection of non-empty sets
    t_nonempty = [s for s in t_sets if s]
    g_nonempty = [s for s in g_sets if s]
    t_inds = (
        frozenset.intersection(*t_nonempty) if t_nonempty else frozenset()
    )
    g_inds = (
        frozenset.intersection(*g_nonempty) if g_nonempty else frozenset()
    )
    if t_inds or g_inds:
        return MergeResult(t_inds, g_inds, MergeCriteria.INTERSECTION_NONEMPTY)

    # 3. Union of all sets (most sensitive)
    t_inds = frozenset.union(*t_sets) if t_sets else frozenset()
    g_inds = frozenset.union(*g_sets) if g_sets else frozenset()
    return MergeResult(t_inds, g_inds, MergeCriteria.UNION)


# ---------------------------------------------------------------------------
# Count enums
# ---------------------------------------------------------------------------

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
