"""
hulkrna.bam — BAM read-pair iterator and CIGAR parsing.

Requires a **name-sorted** (``samtools sort -n``) or **collated**
(``samtools collate``) BAM file so that all alignments for a read
name arrive consecutively.  This eliminates the unbounded-memory
pair buffer required for coordinate-sorted input.

Read-name groups are processed atomically: all alignments for a
query name are collected, filtered, paired by mate position, and
yielded as a list of ``(r1, r2)`` tuples.  For unique mappers
(NH=1) the list has one element; for multimappers (NH>1, when
``include_multimap=True``) the list may have several.
"""

import logging
from itertools import groupby

import pysam

from .core import Strand

logger = logging.getLogger(__name__)

# N: intron (reference skip)
# M, D, =, X: advance reference within exon
CIGAR_OPS_ADVANCE_REF = {
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
    pysam.CDEL,
}


def parse_read(read):
    """Parse a single alignment's CIGAR into exon blocks and splice junctions.

    Returns
    -------
    tuple
        ``(ref, ref_strand, exons, sjs)`` where *exons* is a list of
        ``(start, end)`` tuples and *sjs* is a list of
        ``(start, end, sj_strand)`` tuples.
    """
    exons = []
    sjs = []
    ref = read.reference_name
    start = read.reference_start
    pos = start
    ref_strand = Strand.from_is_reverse(read.is_reverse)

    for (cig_op, cig_len) in read.cigartuples:
        # handle reference skip (intron)
        if cig_op == pysam.CREF_SKIP:
            sj_strand = Strand.from_str(read.get_tag('XS'))
            if pos > start:
                # add the exon from start to pos
                exons.append((start, pos))
            # add the splice junction from pos to (pos + cig_len)
            sjs.append((pos, pos + cig_len, sj_strand))
            start = pos + cig_len
            pos = start
        elif cig_op in CIGAR_OPS_ADVANCE_REF:
            # advance reference position
            pos += cig_len

    # add the last exon if it exists
    if pos > start:
        exons.append((start, pos))
    return ref, ref_strand, exons, sjs


def _pair_reads(r1_reads, r2_reads):
    """Pair read-1 and read-2 alignments within a group by mate position.

    Matching uses ``(reference_id, reference_start)`` of each read
    against the mate position encoded in ``(next_reference_id,
    next_reference_start)`` of its partner.

    Reads whose mate is unmapped (``mate_is_unmapped``) are yielded
    as singletons: ``(read, None)`` or ``(None, read)``.

    Parameters
    ----------
    r1_reads : list[pysam.AlignedSegment]
        Read-1 alignments (filtered, primary only).
    r2_reads : list[pysam.AlignedSegment]
        Read-2 alignments (filtered, primary only).

    Returns
    -------
    list[tuple]
        Each element is ``(r1, r2)``, ``(r1, None)``, or ``(None, r2)``.
    """
    pairs = []

    # Separate mate-unmapped singletons from pairable reads
    paired_r1 = []
    for r1 in r1_reads:
        if r1.mate_is_unmapped:
            pairs.append((r1, None))
        else:
            paired_r1.append(r1)

    paired_r2 = []
    for r2 in r2_reads:
        if r2.mate_is_unmapped:
            pairs.append((None, r2))
        else:
            paired_r2.append(r2)

    # Index r2 reads by their genomic position for O(n) matching
    r2_by_pos: dict[tuple[int, int], list[int]] = {}
    for i, r2 in enumerate(paired_r2):
        key = (r2.reference_id, r2.reference_start)
        r2_by_pos.setdefault(key, []).append(i)

    used_r2: set[int] = set()
    for r1 in paired_r1:
        mate_key = (r1.next_reference_id, r1.next_reference_start)
        candidates = r2_by_pos.get(mate_key, [])
        matched = False
        for idx in candidates:
            if idx in used_r2:
                continue
            r2 = paired_r2[idx]
            # Verify reciprocal mate pointers
            if (r2.next_reference_id == r1.reference_id
                    and r2.next_reference_start == r1.reference_start):
                pairs.append((r1, r2))
                used_r2.add(idx)
                matched = True
                break
        if not matched:
            # r1 with no matching r2 → singleton
            pairs.append((r1, None))

    # Remaining unmatched r2 → singletons
    for i, r2 in enumerate(paired_r2):
        if i not in used_r2:
            pairs.append((None, r2))

    return pairs


def parse_bam_file(
    bam_iter,
    stats,
    *,
    skip_duplicates=True,
    include_multimap=False,
):
    """Parse a name-sorted / collated BAM, yielding read-pair groups.

    All alignments for a read name (query name) must arrive
    consecutively — guaranteed by ``samtools sort -n`` or
    ``samtools collate``.

    Within each read-name group, r1 and r2 alignments are paired
    by matching ``(next_reference_id, next_reference_start)`` to the
    mate's ``(reference_id, reference_start)``.  Unpaired reads and
    reads whose mate is unmapped are yielded as singletons.

    Parameters
    ----------
    bam_iter : iterator
        Iterator over ``pysam.AlignedSegment`` objects.
    stats : dict
        Mutable dict populated with filtering statistics.
    skip_duplicates : bool
        If True (default), discard reads marked as PCR/optical
        duplicates.
    include_multimap : bool
        If True, yield multimapping read-name groups (NH > 1).
        If False (default), skip them (but still count in stats).

    Yields
    ------
    list[tuple[AlignedSegment | None, AlignedSegment | None]]
        A list of ``(r1, r2)`` pairs for one read name.  For unique
        mappers the list has one element; for multimappers it may
        have several.  Singletons appear as ``(r1, None)`` or
        ``(None, r2)``.
    """
    stat_keys = [
        'total', 'qc_fail', 'unmapped', 'secondary',
        'supplementary', 'duplicate',
        'n_read_names', 'unique', 'multimapping',
        'proper_pair', 'improper_pair', 'mate_unmapped', 'unpaired',
    ]
    for key in stat_keys:
        stats.setdefault(key, 0)

    for _qname, group_iter in groupby(
        bam_iter, key=lambda r: r.query_name
    ):
        r1_reads: list = []
        r2_reads: list = []
        nh = 1

        for read in group_iter:
            stats['total'] += 1

            # Skip unusable alignments
            if read.is_qcfail:
                stats['qc_fail'] += 1
                continue
            if read.is_secondary:
                stats['secondary'] += 1
                continue
            if read.is_supplementary:
                stats['supplementary'] += 1
                continue
            if read.is_unmapped:
                stats['unmapped'] += 1
                continue

            # Duplicate handling
            if read.is_duplicate:
                stats['duplicate'] += 1
                if skip_duplicates:
                    continue

            # Enforce paired-end
            assert read.is_paired, "Input BAM must be paired-end"

            # Read NH from first usable alignment in the group
            if read.has_tag('NH'):
                nh = read.get_tag('NH')

            if read.is_read1:
                r1_reads.append(read)
            elif read.is_read2:
                r2_reads.append(read)

        # Skip empty groups (everything filtered)
        if not r1_reads and not r2_reads:
            continue

        stats['n_read_names'] += 1

        # Multimap handling
        if nh > 1:
            stats['multimapping'] += 1
            if not include_multimap:
                continue
        else:
            stats['unique'] += 1

        # Pair r1 and r2 alignments by mate position
        pairs = _pair_reads(r1_reads, r2_reads)

        # Track pair-level stats
        for r1, r2 in pairs:
            if r1 is None or r2 is None:
                stats['mate_unmapped'] += 1
            elif r1.is_proper_pair:
                stats['proper_pair'] += 1
            else:
                stats['improper_pair'] += 1

        yield pairs
