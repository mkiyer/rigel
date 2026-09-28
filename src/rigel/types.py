"""
rigel.types — Foundational data types used across the rigel pipeline.

This module defines the core types for genomic coordinates, strand
orientation and interval representation.

Coordinate convention: all coordinates are 0-based, half-open (BED style).
"""

from enum import IntEnum
from typing import NamedTuple


# ---------------------------------------------------------------------------
# Strand
# ---------------------------------------------------------------------------

# Pre-computed mapping for Strand.from_str() — lives at module level
# because IntEnum's metaclass interprets class-level dicts as members.
_STRAND_STR_MAP: dict[str, int] = {".": 0, "+": 1, "-": 2, "?": 3}


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
        try:
            return cls(_STRAND_STR_MAP[s])
        except KeyError:
            raise ValueError(f"Invalid strand string: {s!r}") from None

    def to_str(self) -> str:
        """Return the single-character string representation."""
        return (".", "+", "-", "?")[self.value]


# ---------------------------------------------------------------------------
# Interval
# ---------------------------------------------------------------------------


class Interval(NamedTuple):
    """A simple 0-based half-open interval (start, end).

    Used by ``Transcript`` to store exon coordinates. For intervals
    positioned on a specific reference/strand, use ``GenomicInterval``.
    """

    start: int
    end: int


# ---------------------------------------------------------------------------
# GenomicInterval
# ---------------------------------------------------------------------------


class GenomicInterval(NamedTuple):
    """An interval located on a specific reference and strand.

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

    ``EXON`` is one exon of one transcript.  ``TRANSCRIPT`` marks the
    full transcript span ``[start, end)`` — intron overlap is derived as
    ``transcript_bp - exon_bp``.  ``INTERGENIC`` covers the reference outside
    every transcript.

    ``SJ`` is an annotated splice junction that exactly matches a known
    intron in the transcript reference.
    """

    EXON = 0
    TRANSCRIPT = 1
    INTERGENIC = 2
    SJ = 3


# ---------------------------------------------------------------------------
# AnnotatedInterval
# ---------------------------------------------------------------------------


class AnnotatedInterval(NamedTuple):
    """A reference-annotated genomic interval with index metadata.

    Used for both tiling intervals (EXON/TRANSCRIPT/INTERGENIC in the
    cgranges overlap index) and splice junctions (SJ in the exact-match
    lookup).  Gene index is derived from the transcript table at load
    time — not stored per-interval.
    """

    ref: str
    start: int
    end: int
    strand: int = Strand.NONE
    interval_type: int = IntervalType.INTERGENIC
    t_index: int = -1
