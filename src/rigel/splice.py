"""rigel.splice — the splice classification enums.

A fragment's splice type depends on its junctions alone:

- UNSPLICED: no splice junction;
- SPLICED_ANNOT: junctions, every one held by one annotated transcript. It is the only spliced type, because the scan
  withholds a hit whose junctions no one annotated transcript holds (the annotated-path rule);
- SPLICED_IMPLICIT: the mate gap spans an annotated intron of a candidate;
- SPLICE_ARTIFACT: a CIGAR junction the splice blacklist rejected.

Exonic and intronic overlap are recorded per candidate (the resolver's ``exon_bp`` / ``intron_bp``), not here. Count
arrays have 8 columns, 4 splice types × 2 strands (sense, antisense); the output sums the two strands.
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

    The scanner classifies, so the scanner counts; the breakdown is QC for the report. The name is derived, not
    tabled, so there is no second correspondence to keep in step: the C++ keys these counters off
    ``splice_type_label`` (``bam_scanner.cpp``), whose strings are these member names lower-cased.
    """
    return f"n_census_{stype.name.lower()}"


class SpliceStrandCol(IntEnum):
    """A count array's column: ``col = splice_type * 2 + int(is_antisense)``, so even columns are sense and odd
    ones antisense."""

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

# the spliced-annotated columns
SPLICED_COLS = tuple(c.value for c in SpliceStrandCol if c.category == SpliceType.SPLICED_ANNOT)
