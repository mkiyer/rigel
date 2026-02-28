"""
hulkrna.types — Foundational data types used across the hulkrna pipeline.

This module defines the core types for genomic coordinates, strand
orientation, interval representation, and set-merging primitives.

Coordinate convention: all coordinates are 0-based, half-open (BED style).
"""

from dataclasses import dataclass
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

    ``EXON`` marks an individual exon boundary.  ``TRANSCRIPT`` marks the
    full transcript span ``[start, end)`` — intron overlap is derived as
    ``transcript_bp - exon_bp``.  ``INTERGENIC`` fills gaps between genes.

    ``SJ`` is an annotated splice junction that exactly matches a known
    intron in the transcript reference.  ``SJ_UNANNOT`` is a splice
    junction observed in the CIGAR (N-operation) but not matching any
    annotated intron — these are recorded with ``t_index = -1``.
    """
    EXON = 0
    TRANSCRIPT = 1
    INTERGENIC = 2
    SJ = 3
    SJ_UNANNOT = 4


# ---------------------------------------------------------------------------
# RefInterval
# ---------------------------------------------------------------------------

class RefInterval(NamedTuple):
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

    Wraps the merged transcript index frozenset together with
    which relaxation level produced the result.  Gene indices
    are derived from transcript indices via ``t_to_g_arr`` when
    needed, rather than being stored redundantly.
    """
    t_inds: frozenset[int]
    criteria: MergeCriteria

    @property
    def is_unique_transcript(self) -> bool:
        """True if exactly one transcript index in the merged result."""
        return len(self.t_inds) == 1

    @property
    def is_empty(self) -> bool:
        """True if the merge produced no transcript indices."""
        return not self.t_inds


#: Sentinel for "no hits at all" — avoids creating new objects.
EMPTY_MERGE = MergeResult(frozenset(), MergeCriteria.EMPTY)


# ---------------------------------------------------------------------------
# Chimera classification
# ---------------------------------------------------------------------------

class ChimeraType(IntEnum):
    """Classification of chimeric fragments.

    A chimeric fragment has exon blocks that map to disjoint transcript
    sets, indicating fusion or read-through events.

    Values
    ------
    NONE : int
        All exon blocks map to a connected set of transcripts.
    TRANS : int
        Exon blocks span multiple reference sequences (interchromosomal).
        Suggestive of trans-splicing or gene fusions.
    CIS_STRAND_SAME : int
        Intrachromosomal chimera where both disjoint exon-block
        clusters align to the same strand.  Suggestive of
        transcriptional read-through between adjacent genes.
    CIS_STRAND_DIFF : int
        Intrachromosomal chimera where the disjoint exon-block
        clusters align to different strands.  Suggestive of genomic
        rearrangement or trans-splicing.
    """
    NONE = 0
    TRANS = 1
    CIS_STRAND_SAME = 2
    CIS_STRAND_DIFF = 3
