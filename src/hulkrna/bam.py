"""hulkrna.bam — BAM read-pair iterator and CIGAR parsing."""

import pysam

from .core import Strand

# N: intron (reference skip)
# M, D, =, X: advance reference within exon
CIGAR_OPS_ADVANCE_REF = {
    pysam.CMATCH, 
    pysam.CEQUAL, 
    pysam.CDIFF, 
    pysam.CDEL
}

def parse_read(read):
    '''modified from pysam get_blocks to return alignment blocks'''
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


def parse_bam_file(
    bam_iter,
    stats,
    skip_duplicates=True,
    multimap_bamfh=None
):
    '''
    Generator function to parse a coordinate-sorted BAM file
    and group reads by query name.

    Args:
        bam_iter (iterator): An iterator over pysam AlignedSegment objects.
        stats (dict): A dictionary to populate with filtering statistics.
        skip_duplicates (bool): if true, skip reads marked as duplicates.
        multimap_bamfh (pysam.AlignmentFile|None): BAM handle to write
            multimapping reads (NH > 1). If None, multimappers are discarded.
    Yields:
        tuple: A tuple containing the (read1, read2) AlignedSegment objects
               for each pair, or (read, None)/(None, read) for singletons.
    '''
    # a buffer to hold reads waiting for their mate.
    pair_buf = {}

    # init stats
    stat_keys = [
        'total', 'qc_fail', 'unmapped', 'secondary',
        'supplementary', 'multimapping', 'unique', 
        'duplicate', 'improper_pair', 'proper_pair',
        'mate_unmapped', 'unpaired'
    ]
    for key in stat_keys:
        stats.setdefault(key, 0)

    for read in bam_iter:
        stats['total'] += 1

        # skip unusable reads
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

        # handle multimapping reads
        nh = read.get_tag('NH')
        if nh > 1:
            stats['multimapping'] += 1
            if multimap_bamfh is not None:
                multimap_bamfh.write(read)
            continue

        # optionally skip duplicate reads
        stats['duplicate'] += int(read.is_duplicate)
        if skip_duplicates and read.is_duplicate:
            continue
        stats['unique'] += 1

        # enforce paired-end
        assert read.is_paired, "Input BAM is not paired-end"
        assert read.is_read1 or read.is_read2, "Read is not paired correctly"
        
        # Handle reads with unmapped mates (yield as singletons)
        if read.mate_is_unmapped:
            stats['mate_unmapped'] += 1
            yield (read, None) if read.is_read1 else (None, read)
            continue

        # track proper pair flag
        if read.is_proper_pair:
            stats['proper_pair'] += 1
        else:
            stats['improper_pair'] += 1

        # pair reads using buffer
        k = read.query_name
        if k in pair_buf:
            # retrieve mate from buffer
            mate = pair_buf.pop(k)
            # yield the pair, ensuring read1 is always first in the tuple.
            yield (read, mate) if read.is_read1 else (mate, read)
        else:
            # store read in buffer and wait for its mate.
            pair_buf[k] = read
 

