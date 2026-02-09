#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import json
from pathlib import Path
from enum import Enum, IntEnum
from typing import NamedTuple, Set, List, Tuple

import numpy as np
import pandas as pd
import pysam
import cgranges


# import local modules
SCRIPT_PATH = Path(__file__).resolve().parent
HULKRNA_LIB_PATH = SCRIPT_PATH / '..' / 'lib'
if HULKRNA_LIB_PATH.exists() and str(HULKRNA_LIB_PATH) not in sys.path:
    sys.path.insert(0, str(HULKRNA_LIB_PATH))

# import from local modules
from hulkrna.txref import IntervalType, Strand, TxIndex, _merge_sets_with_relaxation
from hulkrna.count import CountCategory, CountStrand, CountType
from hulkrna.bam import parse_bam_file

# setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


TRANSCRIPT_COUNTS_TSV_GZ_FILE = 'transcript_counts.tsv.gz'
TRANSCRIPT_COUNTS_FEATHER_FILE = 'transcript_counts.feather'
GENE_COUNTS_TSV_GZ_FILE = 'gene_counts.tsv.gz'
GENE_COUNTS_FEATHER_FILE = 'gene_counts.feather'
SUMMARY_JSON_FILE = 'summary.json'
BAM_STATS_JSON_FILE = 'bam_stats.json'


class Fragment(NamedTuple):
    '''read fragment information'''
    refs: Set[str]
    strand: Strand
    exons: Set[Tuple[int, int]]
    sjs: Set[Tuple[int, int, int]]
    is_duplicate: bool
    size: int = 0  # size of the fragment, not used in counting


class StrandProtocol(Enum):
    UNSTRANDED = 'unstranded'
    FR = 'fr'  # read1 sense, read2 antisense
    RF = 'rf'  # read1 antisense, read2 sense
    FF = 'ff'  # both reads sense
    RR = 'rr'  # both reads antisense

    # mapping reads to strand using strand protocol
    __READ_TO_STRAND = {
        'unstranded': [(Strand.NONE, Strand.NONE), (Strand.NONE, Strand.NONE)], # unstranded
        'fr': [(Strand.POS, Strand.NEG), (Strand.NEG, Strand.POS)],  # read1 sense, read2 antisense
        'rf': [(Strand.NEG, Strand.POS), (Strand.POS, Strand.NEG)],  # read1 antisense, read2 sense
        'ff': [(Strand.POS, Strand.NEG), (Strand.POS, Strand.NEG)],  # both reads sense
        'rr': [(Strand.NEG, Strand.POS), (Strand.NEG, Strand.POS)],  # both reads antisense
    }

    def get_map(self):
        return self.__READ_TO_STRAND[self.value]

    @staticmethod
    def infer(read, strand_protocol_map):
        '''
        Infer the strand of a read based on protocol and read flags.
        Args:
            read: pysam.AlignedSegment
            strand_protocol: list of tuples (see STRAND_PROTOCOLS above)
        '''
        if read.is_paired:
            assert read.is_read1 != read.is_read2, "Read is not paired correctly"
            readnum = 0 if read.is_read1 else 1
        else:
            readnum = 0
        return strand_protocol_map[readnum][read.is_reverse]


# N: intron (reference skip)
# M, D, =, X: advance reference within exon
CIGAR_OPS_ADVANCE_REF = {
    pysam.CMATCH, 
    pysam.CEQUAL, 
    pysam.CDIFF, 
    pysam.CDEL
}


def parse_read(read, strand_protocol_map):
    '''modified from pysam get_blocks to return strand, exons, splice junctions'''
    exons = []
    sjs = []
    ref = read.reference_name
    start = read.reference_start
    pos = start
    read_strand = StrandProtocol.infer(read, strand_protocol_map)

    for (cig_op, cig_len) in read.cigartuples:
        # handle reference skip (intron)
        if cig_op == pysam.CREF_SKIP:
            assert read.has_tag('XS'), "Read must have XS tag for splice junctions"
            sj_strand = Strand.from_str(read.get_tag('XS'))
            # sj_strand = read.get_tag('XS')
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
    return ref, read_strand, exons, sjs


def merge_fragment_reads(r1, r2, strand_protocol_map):
    '''Merge alignment information from paired reads r1 and r2'''
    refs = set()
    exons = set()
    sjs = set()
    strand = Strand.NONE
    duplicate = False

    for read in (r1, r2):
        if read is None:
            continue
        # parse read and update sets
        rref, rstrand, rexons, rsjs = parse_read(read, strand_protocol_map)
        refs.add(rref)
        exons.update(rexons)
        sjs.update(rsjs)
        strand = strand | rstrand # bitwise OR
        duplicate = duplicate or read.is_duplicate
    return Fragment(refs, strand, exons, sjs, duplicate)


class ReadCounter:

    NUM_STRAND_TYPES = len(CountStrand)

    def __init__(self, num_transcripts: int, num_genes: int):
        # count attributes
        self.ambiguous_strand = 0
        self.duplicates = 0
        self.chimeras = 0
        self.intergenic = 0
        self.no_feature = 0
        self.spliced_no_feature = 0
        # use numpy arrays for counting
        self.t_counts_arr = np.zeros((num_transcripts, len(CountType)), dtype=np.float32)
        self.g_counts_arr = np.zeros((num_genes, len(CountType)), dtype=np.float32)

    @staticmethod
    def _increment(counts_arr, inds, count_type):
        counts_arr[inds, count_type] += (1.0 / len(inds))

    def _update_feature_counts(self, inds, ref_strand_arr, counts_arr, count_category, fragment_strand):
        """A general helper to determine strandedness and update counts for a feature (gene or transcript)."""
        if len(inds) == 1:
            # Unique feature match
            ref_strand = ref_strand_arr[inds[0]]
            strand_type = CountStrand.SENSE if (fragment_strand == ref_strand or fragment_strand == Strand.NONE) else CountStrand.ANTISENSE
        else:
            # Ambiguous feature match
            strand_type = CountStrand.AMBIGUOUS
        count_type = CountType(ReadCounter.NUM_STRAND_TYPES * count_category + strand_type).value
        self._increment(counts_arr, inds, count_type)

    def _assign_counts(self, t_inds, g_inds, count_category, strand, txindex):
        """Assigns counts by calling the generalized helper for transcripts and genes."""
        if t_inds:
            self._update_feature_counts(
                inds=t_inds,
                ref_strand_arr=txindex.t_to_strand_arr,
                counts_arr=self.t_counts_arr,
                count_category=count_category,
                fragment_strand=strand
            )
        if g_inds:
            self._update_feature_counts(
                inds=g_inds,
                ref_strand_arr=txindex.g_to_strand_arr,
                counts_arr=self.g_counts_arr,
                count_category=count_category,
                fragment_strand=strand
            )

    def get_t_counts_df(self):
        return pd.DataFrame(self.t_counts_arr, columns=CountType.columns())
    
    def get_g_counts_df(self):
        return pd.DataFrame(self.g_counts_arr, columns=CountType.columns())

    def summary(self):
        d = {
            "ambiguous_strand": self.ambiguous_strand,
            "duplicates": self.duplicates,
            "chimeras": self.chimeras,
            "intergenic": self.intergenic,
            "no_feature": self.no_feature,
            "spliced_no_feature": self.spliced_no_feature,
        }
        # Create pandas Series before calling .to_dict()
        g_sum = pd.Series(self.g_counts_arr.sum(axis=0).round(), index=CountType.columns()).to_dict()
        t_sum = pd.Series(self.t_counts_arr.sum(axis=0).round(), index=CountType.columns()).to_dict()
        d.update({
            "g_counts": g_sum,
            "t_counts": t_sum,
        })
        return d

    def summary_by_type(self):
        """Summarizes total counts by read type for genes and transcripts."""
        # Create pandas DataFrames before manipulating the index
        g_sum = pd.DataFrame([self.g_counts_arr.sum(axis=0)], columns=CountType.columns(), index=['g_counts'])
        t_sum = pd.DataFrame([self.t_counts_arr.sum(axis=0)], columns=CountType.columns(), index=['t_counts'])
        # Concatenate summaries
        summary = pd.concat([g_sum, t_sum], axis=0)
        return summary


def profile_strand_protocol(input_bam_file, stop_after=0, threads=1):
    # open input bam file
    if input_bam_file == '-':
        input_bam_file = sys.stdin

    strand_counts = {
        (True, True): 0,  # read1 sense, read2 sense
        (True, False): 0, # read1 sense, read2 antisense
        (False, True): 0, # read1 antisense, read2 sense
        (False, False): 0 # read1 antisense, read2 antisense
    }

    with pysam.AlignmentFile(input_bam_file, "rb", threads=threads) as bamfh:
        # parse bam file and return properly paired reads
        i = 0
        stats = {}
        for r1, r2 in parse_bam_file(bamfh, stats):
            # only consider reads with XS tag 
            # (requires aligner that sets this flag)
            # (spliced reads considered RNA reads, unspliced could be DNA)
            sj_strand = Strand.NONE
            for read in (r1, r2):
                if read and read.has_tag('XS'):
                    sj_strand |= Strand.from_str(read.get_tag('XS'))
            
            if Strand.NONE < sj_strand < Strand.AMBIGUOUS:
                r1_strand = Strand.from_is_reverse(r1.is_reverse) if r1 else Strand.NONE
                r2_strand = Strand.from_is_reverse(r2.is_reverse) if r2 else Strand.NONE
                key = (r1_strand == sj_strand, r2_strand == sj_strand)
                strand_counts[key] += 1

            i += 1
            if (stop_after > 0) and (i >= stop_after):
                break
            if i % 1000000 == 0:
                logging.info(f'Processed {i} reads, current strand counts: {dict(strand_counts)}')
    return strand_counts


def infer_strand_protocol(strand_counts, 
                          min_reads=1000,
                          threshold_frac=0.95):
    _STRAND_PROTOCOL_ALIASES = {
        (True, False): 'fr',  # read1 sense, read2 antisense
        (False, True): 'rf',  # read1 antisense, read2 sense
        (True, True): 'ff',   # both reads sense
        (False, False): 'rr'  # both reads antisense    
    }

    # sort counts by value
    sorted_strand_counts = sorted(strand_counts.items(), key=lambda kv: kv[1], reverse=True)
    # compute relative fraction of each strand protocol
    tot_counts = sum(count for _, count in sorted_strand_counts)
    sorted_strand_count_frac = [(k, v / tot_counts) for k, v in sorted_strand_counts]

    if tot_counts < min_reads:
        logging.warning(f'Insufficient reads ({tot_counts} < {min_reads}) to infer strand protocol')
        return 'unstranded', None
    # get the most common protocol
    strand_protocol = sorted_strand_counts[0][0]
    frac = sorted_strand_count_frac[0][1]
    if frac < threshold_frac:
        logging.warning(f'Most common strand protocol {strand_protocol} has fraction {frac:.2f}, '
                        f'which is below the threshold {threshold_frac}. '
                        f'Using unstranded protocol instead.')
        return 'unstranded', frac
    return _STRAND_PROTOCOL_ALIASES[strand_protocol], frac


def hulkrna_count(input_bam_file, index_dir, output_dir,
                  strand_protocol='infer',
                  ignore_duplicates=True,
                  threads=1,
                  feather_compression='lz4'):
    # setup feather options
    feather_kwargs = {'compression': feather_compression}
    
    # load index files
    logging.info(f'Loading index files from {index_dir}')
    txindex = TxIndex.load(index_dir)
    rcounter = ReadCounter(txindex.num_transcripts, txindex.num_genes)

    # infer strand protocol
    if strand_protocol == 'infer':
        logging.info(f'[START] Inferring strand protocol from {input_bam_file}')
        strand_counts = profile_strand_protocol(input_bam_file, stop_after=0, threads=threads)
        logging.info(f'[DONE] Strand counts: {dict(strand_counts)}')
        # TODO: setting min_reads=1 and threshold_frac=0 to always return a stranded protocol
        # but in the future we may want to make these parameters configurable
        strand_protocol, strand_frac = infer_strand_protocol(strand_counts, min_reads=1, threshold_frac=0)
        logging.info(f'Inferred strand protocol: {strand_protocol} with fraction {strand_frac:.2f}')

    # convert to StrandProtocol enum
    logging.info(f'Using strand protocol: {strand_protocol}')
    strand_protocol = StrandProtocol(strand_protocol)
    strand_protocol_map = strand_protocol.get_map()

    # open output bam file
    # obamfh = pysam.AlignmentFile(output_bam_file, "wb", template=bamfh)

    # track bam stats
    bam_stats = {}
    bam_stats['strand_protocol'] = strand_protocol.value
    bam_stats['strand_frac'] = strand_frac

    # open input bam file
    if input_bam_file == '-':
        input_bam_file = sys.stdin
    with pysam.AlignmentFile(input_bam_file, "rb", threads=threads) as bamfh:
        # parse bam file and return paired reads
        logging.info(f'[START] Parsing BAM file {input_bam_file}')
        i = 0
        for r1, r2 in parse_bam_file(bamfh, bam_stats):
            i += 1
            if i % 100000 == 0:
                logging.debug(f'Processed {i} fragments '
                            f'(no_feature={rcounter.no_feature}) '
                            f'(spliced_no_feature={rcounter.spliced_no_feature}) '
                            f'(intergenic={rcounter.intergenic}) '
                            f'(duplicates={rcounter.duplicates}) '
                            f'(chimeras={rcounter.chimeras}) '
                            f'(ambiguous_strand={rcounter.ambiguous_strand})')
                s = rcounter.summary_by_type()
                logging.debug(f'Counts summary:\n{s}')

            # parse reads into references, strand, exons, sjs
            fragment = merge_fragment_reads(r1, r2, strand_protocol_map)
            # refs, strand, exons, sjs, is_duplicate

            # duplicates
            if fragment.is_duplicate and ignore_duplicates:
                rcounter.duplicates += 1
                continue
            if len(fragment.refs) > 1:
                # reads align to multiple references
                rcounter.chimeras += 1
                continue
            if fragment.strand == Strand.AMBIGUOUS:
                # fragment maps to multiple strands suggestive of a chimera
                # or chromosomal rearrangement
                rcounter.ambiguous_strand += 1
                continue

            # single reference
            ref = fragment.refs.pop()

            # search for overlapping exons from transcripts and genes
            t_inds, g_inds = txindex.search_intervals(ref, fragment.exons, IntervalType.EXON)

            # check if fragment overlaps exon
            if t_inds or g_inds:
                # check for annotated splice junctions and exons
                if fragment.sjs:
                    # check if splice junctions match annotated splice junctions
                    sj_t_inds, sj_g_inds = txindex.search_splice_junctions(ref, fragment.sjs)
                    if sj_t_inds or sj_g_inds:
                        # merge exon and splice junction matches
                        t_inds, g_inds = _merge_sets_with_relaxation(
                            t_sets = [t_inds, sj_t_inds],
                            g_sets = [g_inds, sj_g_inds]
                        )
                        count_cat = CountCategory.SPLICED_ANNOT
                    else:
                        # no annotated splice junctions found (but fragment overlaps exons)
                        count_cat = CountCategory.SPLICED_UNANNOT
                else:
                    count_cat = CountCategory.UNSPLICED
                # assign counts to matching transcripts and genes
                if t_inds or g_inds:
                    rcounter._assign_counts(list(t_inds), list(g_inds), count_cat, fragment.strand, txindex)
                else:
                    # should not reach here
                    # TODO: remove this else block
                    assert False, "Annotated splice junctions found but no overlapping exons"
                    rcounter.spliced_no_feature += 1
            # fragment does not overlap exons
            else:
                # TODO: unannotated splice junctions that do not overlap exons may be included here
                # so intronic/intergenic reads include spliced and unspliced reads
                # TODO: may want to separate spliced/unspliced intronic/intergenic reads in the future
                # search for overlap with introns
                t_inds, g_inds = txindex.search_intervals(ref, fragment.exons, IntervalType.INTRON)
                if t_inds or g_inds:
                    rcounter._assign_counts(list(t_inds), list(g_inds), CountCategory.INTRON, fragment.strand, txindex)
                else:
                    # no overlap with introns, check intergenic regions
                    # TODO: should not need to check this, but doing it for safety/debugging
                    if txindex.search_intergenic(ref, fragment.exons):
                        rcounter.intergenic += 1
                    else:
                        # TODO: should never reach here
                        rcounter.no_feature += 1

        # final logging
        logging.debug(f'Processed {i} fragments '
                    f'(no_feature={rcounter.no_feature}) '
                    f'(spliced_no_feature={rcounter.spliced_no_feature}) '
                    f'(intergenic={rcounter.intergenic}) '
                    f'(duplicates={rcounter.duplicates}) '
                    f'(chimeras={rcounter.chimeras}) '
                    f'(ambiguous={rcounter.ambiguous_strand})')

        # create output directory
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        # write bam stats
        logging.info(f'[START] Writing bam stats')
        bam_stats_file = output_dir / BAM_STATS_JSON_FILE
        with open(bam_stats_file, 'w') as f:
            json.dump(bam_stats, f, indent=2)
        logging.info(f'[DONE] Wrote bam stats to {bam_stats_file}')

        # write read count summary to output directory
        logging.info(f'[START] Writing read count summary')
        summary_file = output_dir / SUMMARY_JSON_FILE
        with open(summary_file, 'w') as f:
            json.dump(rcounter.summary(), f, indent=2)
        logging.info(f'[DONE] Wrote read count summary to {summary_file}')

        # write transcript counts to output directory
        logging.info(f'[START] Writing transcripts')
        df = pd.concat([txindex.t_df, rcounter.get_t_counts_df()], axis=1)
        df.to_csv(output_dir / TRANSCRIPT_COUNTS_TSV_GZ_FILE, 
                  sep='\t', index=False, header=True)
        logging.info(f'[DONE] Wrote transcripts to TSV file')
        df.to_feather(output_dir / TRANSCRIPT_COUNTS_FEATHER_FILE, **feather_kwargs)
        logging.info(f'[DONE] Wrote transcripts to feather file')

        # write gene counts to output directory
        logging.info(f'[START] Writing genes')
        df = pd.concat([txindex.g_df, rcounter.get_g_counts_df()], axis=1)
        columns_to_output = [x for x in df.columns if x not in ['t_index', 't_id']]
        df = df[columns_to_output]
        df.to_csv(output_dir / GENE_COUNTS_TSV_GZ_FILE, 
                  sep='\t', index=False, header=True)
        logging.info(f'[DONE] Wrote genes to TSV file')
        df.to_feather(output_dir / GENE_COUNTS_FEATHER_FILE, **feather_kwargs)
        logging.info(f'[DONE] Wrote genes to feather file')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-@', '--threads', type=int, default=1, 
                        help='threads to use (default: 1)')
    parser.add_argument('--feather-compression', dest='feather_compression', default='lz4',
                        choices=['lz4', 'zstd', 'uncompressed'],
                        help='compression for feather files (default: lz4)')   
    parser.add_argument('--strand-protocol', dest='strand_protocol', 
                        type=str,
                        choices=['infer'] + [x.value for x in StrandProtocol],
                        default='infer',
                        help='strand protocol to use (default: infer)')
    parser.add_argument('--nodup', dest='ignore_duplicates', default=False, 
                        action='store_true', help='do not count duplicates')
    parser.add_argument('--index', dest='index_dir', required=True, 
                        help='directory with index files')
    parser.add_argument('-i', '--input', dest='bam_file', required=True, 
                        help='input bam file ("-" for stdin)')
    parser.add_argument('-o', '--output', dest='output_dir', required=True, 
                        help='output directory')
    args = parser.parse_args()
    # check args
    if not os.path.exists(args.index_dir):
        parser.error(f'index directory {args.index_dir} not found')
    if args.bam_file != '-' and not os.path.exists(args.bam_file):
        parser.error(f'input {args.bam_file} not found')

    return(hulkrna_count(
        input_bam_file = args.bam_file, 
        index_dir = args.index_dir, 
        output_dir = args.output_dir,
        strand_protocol = args.strand_protocol,
        ignore_duplicates = args.ignore_duplicates, 
        threads = args.threads,
        feather_compression = args.feather_compression
    ))


def snakemake_main():
    from snakemake.script import snakemake
    # Use the snakemake object to get outputs, params, and threads
    # Run the main build function with inputs from the Snakemake rule
    return (hulkrna_count(
        input_bam_file = snakemake.input.bam,
        index_dir = snakemake.input.index_dir,
        output_dir = snakemake.params.output_dir,
        strand_protocol = snakemake.params.strand_protocol,
        ignore_duplicates = snakemake.params.ignore_duplicates,
        threads = snakemake.threads,
        feather_compression = snakemake.params.feather_compression
    ))


def is_snakemake_run():
    """Check if the 'snakemake' global variable exists."""
    try:
        snakemake
        return True
    except NameError:
        return False

if __name__ == "__main__":
    if is_snakemake_run():
        sys.exit(snakemake_main())
    else:
        sys.exit(main())
