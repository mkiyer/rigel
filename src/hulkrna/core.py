"""
hulkrna.core — Shared data types used across multiple hulkrna modules.

This module defines the foundational types for genomic coordinates, strand
orientation, and interval representation. Types that are specific to a single
command (e.g., counting-only types) live in their own modules instead.

Coordinate convention: all coordinates are 0-based, half-open (BED style).
"""

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
    """Classification of a genomic interval relative to gene annotations."""
    EXON = 0
    INTRON = 1
    INTERGENIC = 2
    SJ = 3


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
