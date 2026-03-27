#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import collections
from enum import Enum, IntEnum

import numpy as np
import pandas as pd
import pysam

# import from local modules
from hulkrna.txref import Strand
from hulkrna.bam import parse_bam_file

# setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


CIGAR_OPS_ADVANCE_REF = {
    pysam.CMATCH,
    pysam.CEQUAL,
    pysam.CDIFF,
    pysam.CDEL
}


def infer_strand_from_sj(read, sj_strand_map):
    '''infer strand from annotated splice junctions in the read'''

    ref = read.reference_name
    start = read.reference_start
    pos = start
    strands = set()

    for (cig_op, cig_len) in read.cigartuples:
        # handle reference skip (intron)
        if cig_op == pysam.CREF_SKIP:
            start, end = pos, pos + cig_len
            annot_strands = sj_strand_map.get((ref, start, end), set())
            strands.update(annot_strands)
            start = pos + cig_len
            pos = start
        elif cig_op in CIGAR_OPS_ADVANCE_REF:  
            # advance reference position
            pos += cig_len

    if not strands:
        return Strand.NONE
    elif len(strands) > 1:
        # multiple strands found, return ambiguous
        return Strand.AMBIGUOUS
    return strands.pop()


def hulkrna_strand(index_dir,
                   input_bam_file, 
                   output_file, 
                   threads=1,
                   tmp_dir=None):    
    # build splice junction index
    logging.debug(f'Reading splice junctions')
    df = pd.read_feather(os.path.join(index_dir, 'sj.feather'))
    logging.debug(f'Building splice junction index')
    # map splice junctions to strand
    sj_strand_map = df.groupby(['ref', 'start', 'end'])['strand'].agg(set).to_dict()
    logging.debug(f'Splice junctions found: {len(sj_strand_map)}')

    # open input bam file
    if input_bam_file == '-':
        input_bam_file = sys.stdin
    bamfh = pysam.AlignmentFile(input_bam_file, "rb", threads=threads)

    # parse bam file and return paired reads
    logging.info(f'[START] Parsing BAM file {input_bam_file}')
    i = 0
    strand_counts = collections.defaultdict(int)
    for r1, r2 in parse_bam_file(bamfh):
        sj_strand = Strand.NONE
        r1_strand = Strand.NONE
        r2_strand = Strand.NONE
        if r1:
            r1_strand = Strand.from_is_reverse(r1.is_reverse)
            sj_strand = sj_strand | infer_strand_from_sj(r1, sj_strand_map)
        if r2:
            r2_strand = Strand.from_is_reverse(r2.is_reverse)
            sj_strand = sj_strand | infer_strand_from_sj(r2, sj_strand_map)

        if Strand.NONE < sj_strand < Strand.AMBIGUOUS:
            key = (r1_strand == sj_strand, r2_strand == sj_strand)
            strand_counts[key] += 1

        i += 1
        if i % 1000000 == 0:
            logging.info(f'Processed {i} reads, current strand counts: {dict(strand_counts)}')
    logging.info(f'[END] Parsed {i} pairs from BAM file {input_bam_file}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-@', '--threads', type=int, default=1, help='threads for compression')
    parser.add_argument('--tmpdir', type=str, default=None, help='temporary directory')
    parser.add_argument('--index', dest='index_dir', required=True, 
                        help='directory with index files')
    parser.add_argument('-i', '--input', dest='bam_file', required=True, 
                        help='input bam file ("-" for stdin)')
    parser.add_argument('-o', '--output', dest='output_file', required=True, 
                        help='output file')
    args = parser.parse_args()
    # check args
    if args.bam_file != '-' and not os.path.exists(args.bam_file):
        parser.error(f'input {args.bam_file} not found')
    return(
        hulkrna_strand(
            args.index_dir,
            args.bam_file, 
            args.output_file, 
            threads=args.threads,
            tmp_dir=args.tmpdir
        )
    )


if __name__ == "__main__":
    sys.exit(main())