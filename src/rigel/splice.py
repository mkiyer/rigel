"""
rigel.splice — Splice classification enums for the quantification pipeline.

Internal count arrays use 8 columns: 4 splice types × 2 strands
(sense/antisense).  The output sums the two strand columns; strand is
not exposed in user-facing output.

Fragment splice types are based on splice junction status only:
- UNSPLICED: no splice junctions detected
- SPLICED_ANNOT: splice junctions, every one held by one annotated transcript. There is no other
  spliced type: the scan withholds a hit whose junctions no one annotated transcript holds (the
  annotated-path rule), and the resolver refuses one
- SPLICED_IMPLICIT: the mate gap spans an annotated intron of a candidate
- SPLICE_ARTIFACT: a CIGAR junction the splice blacklist rejected

Intronic vs exonic overlap is recorded separately, per candidate (the
resolver's ``exon_bp`` / ``intron_bp``), not in the splice type.
"""

from enum import IntEnum


class SpliceType(IntEnum):
    """Category of a fragment based on splice junction status."""

    UNSPLICED = 0
    SPLICED_ANNOT = 1
    # the two derived classes — neither is a CIGAR splice:
    SPLICED_IMPLICIT = 2  # PE gap spans an annotated intron of any candidate
    SPLICE_ARTIFACT = 3  # CIGAR sj was rejected by the SJ blacklist


def census_field(stype: SpliceType) -> str:
    """The :class:`~rigel.stats.PipelineStats` field holding the scanner's count of ``stype``.

    The scanner classifies, so the scanner counts. The per-fragment splice breakdown is scanner QC —
    it has no algorithmic consumer, only the report — so it lives where it is generated and is passed
    through nothing to get there. Reading it off a fragment-length category model instead would count
    only the fragments that contributed a LENGTH observation, a population nobody states.

    There is no name table, deliberately. The C++ keys these counters off
    ``splice_type_label`` (``bam_scanner.cpp``), whose strings are exactly these member names
    lower-cased. Writing the correspondence as a mapping would create a second place for it to be
    wrong; deriving it means a category added to :class:`SpliceType` and forgotten in either
    language is caught by ``test_every_splice_type_is_censused`` rather than silently reading zero
    through ``dict.get(key, 0)``.
    """
    return f"n_census_{stype.name.lower()}"


#: ``SpliceType.SPLICED_ANNOT`` as a plain int.
SPLICE_ANNOT: int = int(SpliceType.SPLICED_ANNOT)  # 1


class SpliceStrandCol(IntEnum):
    """Internal column index for count arrays.

    Layout: 4 splice categories × 2 strands (sense/antisense) = 8 columns.

    Column formula::

        col = splice_type * 2 + int(is_antisense)

    where ``splice_type`` is a :class:`SpliceType` value (0..3).

    Even indices (0, 2, 4, 6) are sense; odd indices (1, 3, 5, 7) are
    antisense.
    """

    UNSPLICED_SENSE = 0
    UNSPLICED_ANTISENSE = 1
    SPLICED_ANNOT_SENSE = 2
    SPLICED_ANNOT_ANTISENSE = 3
    SPLICED_IMPLICIT_SENSE = 4
    SPLICED_IMPLICIT_ANTISENSE = 5
    SPLICE_ARTIFACT_SENSE = 6
    SPLICE_ARTIFACT_ANTISENSE = 7

    @property
    def category(self) -> SpliceType:
        return SpliceType(self.value // 2)


NUM_SPLICE_STRAND_COLS = len(SpliceStrandCol)

# Pre-computed column subset (a tuple for immutability and indexing).
SPLICED_COLS = tuple(c.value for c in SpliceStrandCol if c.category == SpliceType.SPLICED_ANNOT)
