

def parse_bam_file(bamfh, stats, skip_duplicates=True):
    '''
    Generator to yield properly paired reads from a coordinate-sorted BAM file.

    This function iterates through a BAM file and yields pairs of reads that
    are marked as being part of a "proper pair". It uses a buffer to hold
    one read of a pair until its mate is found.

    Args:
        bamfh (pysam.AlignmentFile): An opened pysam AlignmentFile object.
        stats (dict): A dictionary to populate with filtering statistics.
        skip_duplicates (bool): if true, skip reads marked as duplicates.
    Yields:
        tuple: A tuple containing the (read1, read2) AlignedSegment objects
               for each proper pair.
    '''
    # a buffer to hold reads waiting for their mate.
    pair_buf = {}

    # init stats
    stat_keys = [
        'total', 'qc_fail', 'unmapped', 'secondary', 
        'supplementary', 'multimapping', 'unique', 
        'duplicate', 'improper_pair', 'proper_pair'
    ]
    for key in stat_keys:
        stats.setdefault(key, 0)

    for read in bamfh.fetch(until_eof=True):
        # track read stats
        stats['total'] += 1
        stats['qc_fail'] += int(read.is_qcfail)
        stats['unmapped'] += int(read.is_unmapped)
        stats['secondary'] += int(read.is_secondary)
        stats['supplementary'] += int(read.is_supplementary)

        # skip unmapped, low quality, or secondary/supplementary reads
        # TODO: secondary/supplementary reads may be usable, but we skip them for now
        if (read.is_unmapped or read.is_qcfail or read.is_secondary or read.is_supplementary):
            continue

        # defer multimapping reads
        nh = read.get_tag('NH')
        if nh > 1:
            stats['multimapping'] += 1
            continue

        # optionally skip duplicate reads
        stats['unique'] += 1
        stats['duplicate'] += int(read.is_duplicate)
        if skip_duplicates and read.is_duplicate:
            continue

        # filter for proper pairs (paired, both mapped, correct orientation)
        # TODO: may be able to use single-ended reads or non-proper pairs, skip for now
        if not read.is_proper_pair:
            stats['improper_pair'] += 1
            continue
        stats['proper_pair'] += 1

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
 
#
# older version that yields singletons as well
#
# def parse_bam_file(bamfh, skip_duplicates=True):
#     '''Generator to yield usable reads from a BAM file'''
#     pair_buf = {}  # query_name: AlignedSegment waiting for mate
#     for read in bamfh.fetch(until_eof=True):
#         # filter reads that are not usable
#         # TODO: secondary/supplementary reads may be usable, but we skip them for now
#         if read.is_unmapped or read.is_qcfail or read.is_secondary or read.is_supplementary:
#             continue
#         if read.is_duplicate and skip_duplicates:
#             continue
#         # defer multimapping reads
#         nh = read.get_tag('NH')
#         if nh > 1:
#             continue
#         # here we have uniquely mapped reads
#         k = read.query_name
#         if read.is_paired:
#             assert read.is_read1 or read.is_read2, "Read is not paired correctly"
#             if read.mate_is_unmapped:
#                 # return singleton read
#                 yield (read, None) if read.is_read1 else (None, read)
#             elif k in pair_buf:
#                 # retrieve mate from buffer
#                 mate = pair_buf.pop(k)
#                 yield (read, mate) if read.is_read1 else (mate, read)
#             else:
#                 # store read in buffer and wait for mate
#                 pair_buf[k] = read
#         else:
#             # unpaired read, yield as singleton
#             yield (read, None)

