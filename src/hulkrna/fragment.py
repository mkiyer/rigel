"""
hulkrna.fragment — Fragment data structures for paired-end RNA-seq reads.

A Fragment consolidates the alignment information from a paired-end read pair
(read1 + read2) into a unified set of aligned exon blocks, splice junctions
(introns), and metadata needed for downstream Bayesian counting.

Aligned regions and splice junctions are stored as ``GenomicInterval``
objects (ref, start, end, strand). Exon block strands come from the read
genomic alignment; splice junction strands come from the XS tag.

R2 strand flip convention
-------------------------
Before merging, read 2's genomic strand is always flipped via
``Strand.opposite()``. In a proper FR pair this makes both reads
report the same strand (the transcript strand); in chimeric
same-strand pairs the flip produces AMBIGUOUS (POS | NEG).
"""

import collections
from dataclasses import dataclass

from .core import GenomicInterval, Strand
from .bam import parse_read


@dataclass(slots=True)
class Fragment:
    """Consolidated alignment information from a paired-end read pair.

    Attributes
    ----------
    exons : tuple[GenomicInterval, ...]
        Merged exon blocks (contiguous aligned regions on the reference).
        Strand is the read's genomic alignment strand (after r2 flip).
    introns : tuple[GenomicInterval, ...]
        Observed splice junctions from CIGAR N operations.
        Strand is from the XS tag (not the read alignment strand).
    insert_size : int | None
        Insert size (outer distance on the reference). None until
        computed from the reference index (requires knowing the
        reference coordinates of both reads).
    """

    exons: tuple[GenomicInterval, ...] = ()
    introns: tuple[GenomicInterval, ...] = ()
    insert_size: int | None = None

    # -- construction ---------------------------------------------------------

    @classmethod
    def from_reads(cls, r1, r2) -> "Fragment":
        """Construct a Fragment from paired-end reads.

        Read 2's genomic strand is flipped via ``Strand.opposite()``
        before merging so that proper FR pairs yield a single concordant
        strand and chimeric same-strand pairs yield AMBIGUOUS.

        Parameters
        ----------
        r1 : pysam.AlignedSegment or None
            Read 1 of the pair.
        r2 : pysam.AlignedSegment or None
            Read 2 of the pair.

        Returns
        -------
        Fragment
        """
        refs = set()
        exon_dict: dict[tuple[str, Strand], list[tuple[int, int]]] = (
            collections.defaultdict(list)
        )
        introns: set[GenomicInterval] = set()

        for i, read in enumerate((r1, r2)):
            if read is None:
                continue

            # Parse the read's CIGAR to get exon blocks and splice junctions
            rref, rstrand, rexons, rsjs = parse_read(read)

            # R2 strand flip: always flip read 2's genomic strand
            if i == 1:
                rstrand = rstrand.opposite()

            refs.add(rref)

            # Accumulate exon blocks keyed by (ref, strand) for merging
            exon_dict[(rref, rstrand)].extend(rexons)

            # Accumulate splice junctions (strand from XS tag, not alignment)
            for start, end, sj_strand in rsjs:
                introns.add(GenomicInterval(rref, start, end, sj_strand))

        # Merge overlapping/adjacent exon blocks within each (ref, strand) group
        merged_exons: list[GenomicInterval] = []
        for (eref, estrand), intervals in exon_dict.items():
            if not intervals:
                continue
            intervals = sorted(intervals)
            cur_start, cur_end = intervals[0]
            for start, end in intervals[1:]:
                if start <= cur_end:
                    cur_end = max(cur_end, end)
                else:
                    merged_exons.append(
                        GenomicInterval(eref, cur_start, cur_end, estrand)
                    )
                    cur_start, cur_end = start, end
            merged_exons.append(
                GenomicInterval(eref, cur_start, cur_end, estrand)
            )

        return cls(
            exons=tuple(merged_exons),
            introns=tuple(introns),
        )
