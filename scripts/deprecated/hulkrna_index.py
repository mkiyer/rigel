#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import collections
from pathlib import Path

import pandas as pd
import pysam

# setup logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# load hulkrna library path
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
HULKRNA_LIB_PATH = os.path.join(SCRIPT_PATH, '..', 'lib')
if os.path.exists(HULKRNA_LIB_PATH) and HULKRNA_LIB_PATH not in sys.path:
    sys.path.insert(0, HULKRNA_LIB_PATH)

# import from local modules
from hulkrna.txref import Interval, IntervalType, SJ
from hulkrna.transcript import Transcript


REF_LENGTHS_FEATHER_FILE = 'ref_lengths.feather'
REF_LENGTHS_TSV_FILE = 'ref_lengths.tsv'
TRANSCRIPTS_TSV_FILE = 'transcripts.tsv'
TRANSCRIPTS_FEATHER_FILE = 'transcripts.feather'
GENES_TSV_FILE = 'genes.tsv'
GENES_FEATHER_FILE = 'genes.feather'
INTERVALS_TSV_FILE = 'intervals.tsv'
INTERVALS_FEATHER_FILE = 'intervals.feather'
SJ_TSV_FILE = 'sj.tsv'
SJ_FEATHER_FILE = 'sj.feather'


def build_sj_map(transcripts):
    '''Build splice junction list from transcripts'''   
    df = [
        SJ(t.ref, start, end, t.strand, t.t_index, t.g_index)
        for t in transcripts
        for start, end in t.introns()
    ]
    # sort by ref, start, end, strand
    df.sort(key=lambda x: (x.ref, x.start, x.end, x.strand))
    df = pd.DataFrame(df, columns=SJ._fields)
    return df


def gen_exon_intervals(t):
    # Emit exons
    for e in t.exons:
        yield Interval(t.ref, e.start, e.end, t.strand, IntervalType.EXON, t.t_index, t.g_index)

def gen_intron_intervals(t):
    # Emit introns (between every pair of exons)
    for e1, e2 in zip(t.exons[:-1], t.exons[1:]):
        yield Interval(t.ref, e1.end, e2.start, t.strand, IntervalType.INTRON, t.t_index, t.g_index)

def gen_transcript_intervals(t):
    yield from gen_exon_intervals(t)
    yield from gen_intron_intervals(t)



def gen_genomic_intervals(transcripts, ref_lengths):
    '''
    Generate exon, intron, and intergenic intervals from sorted transcripts.
    transcripts: list of Transcript objects sorted by (ref, start, end)
    ref_lengths: dict of reference lengths {ref_name: length}
    '''    
    # group transcripts by reference
    ref_transcripts = collections.defaultdict(list)
    for t in transcripts:
        assert t.ref in ref_lengths, f"Transcript {t.t_id} has unknown reference {t.ref}"
        ref_transcripts[t.ref].append(t)

    # generate intervals for each reference
    for ref, ref_length in ref_lengths.items():
        t_list = ref_transcripts.get(ref, [])
        if not t_list:
            # no transcripts in ref, yield full interval as intergenic
            yield Interval(ref, 0, ref_length)
            continue
        # process sorted transcript intervals
        end = 0
        cluster = []
        for t in t_list:
            if t.start > end:
                # yield intergenic between end and t.start
                if end < t.start:
                    yield Interval(ref, end, t.start)
                # emit previous cluster intervals
                for tc in cluster:
                    yield from gen_transcript_intervals(tc)
                cluster = []
            end = max(end, t.end)
            cluster.append(t)
        # emit final cluster intervals
        for tc in cluster:
            yield from gen_transcript_intervals(tc)
        # yield intergenic to end
        if end < ref_length:
            yield Interval(ref, end, ref_length)
            

def build_interval_index(transcripts, ref_lengths):
    '''
    Find exon, intron, and intergenic intervals from transcripts
    transcripts: list of Transcript objects sorted by (ref, start, end)
    '''
    # generate intervals
    intervals = list(gen_genomic_intervals(transcripts, ref_lengths))
    intervals.sort(key=lambda x: (x.ref, x.start, x.end, x.strand))
    df = pd.DataFrame(intervals, columns=Interval._fields)
    return df


def build_gene_map(t_df):
    '''Build a gene dataframe from transcripts dataframe'''

    def _first_ensure_unique(x):
        return x.unique()[0] if len(x.unique()) == 1 else ValueError("Non-unique value found")

    agg_dict = {
        'ref': _first_ensure_unique,
        'start': 'min',
        'end': 'max',
        'strand': _first_ensure_unique,
        't_id': set,
        't_index': set,
        'g_id': _first_ensure_unique,
        'g_name': _first_ensure_unique,
        'g_type': _first_ensure_unique
    }
    g_df = (
        t_df
        .groupby('g_index')
        .agg(agg_dict)
        .reset_index()
    )
    return g_df


def hulkrna_index(fasta_file, gtf_file, output_dir, 
                  feather_compression='lz4'):
    # setup feather options
    feather_kwargs = {'compression': feather_compression}
    # setup output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    #
    # references
    #
    # load references from fasta file
    log_msg = f'Loading references from FASTA file {fasta_file}'
    logging.info(f'[START] {log_msg}')
    fastafh = pysam.FastaFile(fasta_file)
    ref_lengths = {}
    for i in range(fastafh.nreferences):
        ref_name = fastafh.references[i]
        ref_length = fastafh.lengths[i]
        ref_lengths[ref_name] = ref_length
    fastafh.close()
    logging.info(f'[DONE] Read {len(ref_lengths)} references')
    # write references to pandas dataframe
    logging.info(f'[START] Writing references')
    df = pd.DataFrame(list(ref_lengths.items()), columns=['ref', 'length'])
    df.to_feather(output_dir / REF_LENGTHS_FEATHER_FILE, **feather_kwargs)
    df.to_csv(output_dir / REF_LENGTHS_TSV_FILE, sep='\t', index=False, header=True)
    logging.info(f'[DONE] Wrote references')

    #
    # transcripts
    #
    # read transcripts
    logging.info(f'[START] Reading transcripts from gtf {gtf_file}')
    transcripts = Transcript.read_gtf(gtf_file)
    logging.info(f'[DONE] Read {len(transcripts)} transcripts from GTF')
    # assign gene and transcript index 
    g_id_to_index = {}
    for ti, t in enumerate(transcripts):
        t.t_index = ti
        if t.g_id not in g_id_to_index:
            g_id_to_index[t.g_id] = len(g_id_to_index)
        t.g_index = g_id_to_index[t.g_id]
    # write transcripts
    logging.info(f'[START] Writing transcripts')
    t_df = pd.DataFrame(t.to_dict() for t in transcripts)
    t_df.to_feather(output_dir / TRANSCRIPTS_FEATHER_FILE, **feather_kwargs)
    logging.info(f'[DONE] Wrote transcripts to feather file')
    t_df.to_csv(output_dir / TRANSCRIPTS_TSV_FILE, 
                sep='\t', index=False, header=True)
    logging.info(f'[DONE] Wrote transcripts to TSV file')

    #
    # genes
    #
    logging.info(f'[START] Building gene map')
    g_df = build_gene_map(t_df)
    logging.info(f'[DONE] Found {len(g_df)} genes')
    # write genes
    logging.info(f'[START] Writing genes')
    g_df.to_feather(output_dir / GENES_FEATHER_FILE, **feather_kwargs)
    logging.info(f'[DONE] Wrote transcripts to feather file')
    g_df.to_csv(output_dir / GENES_TSV_FILE, 
                sep='\t', index=False, header=True)
    logging.info(f'[DONE] Wrote transcripts to TSV file')

    #
    # splice junctions
    #
    logging.info(f'[START] Building splice junction mappings')
    df = build_sj_map(transcripts)
    logging.info(f'[DONE] Found {len(df)} splice junctions')
    # write splice junctions to output file
    logging.info(f'[START] Writing splice junctions to output files')
    df.to_csv(output_dir / SJ_TSV_FILE, 
              sep='\t', index=False, header=True)
    logging.info(f'[DONE] Wrote splice junctions to tsv')
    df.to_feather(output_dir / SJ_FEATHER_FILE, **feather_kwargs)
    logging.info(f'[DONE] Wrote splice junctions to feather')

    #
    # genomic intervals
    #
    logging.info(f'[START] Finding genomic intervals')
    df = build_interval_index(transcripts, ref_lengths)
    logging.info(f'[DONE] Found {len(df)} genomic intervals')
    # write intervals to output file
    logging.info(f'[START] Writing intervals')
    df.to_csv(output_dir / INTERVALS_TSV_FILE, 
              sep='\t', index=False, header=True)
    logging.info(f'[DONE] Wrote intervals to tsv file')
    df.to_feather(output_dir / INTERVALS_FEATHER_FILE, **feather_kwargs)
    logging.info(f'[DONE] Wrote binary intervals to feather file')
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', dest='fasta_file', required=True,
                        help='input fasta file')
    parser.add_argument('--gtf', dest='gtf_file', required=True, 
                        help='gtf file with gene annotations')
    parser.add_argument('-o', '--output-dir', dest='output_dir', required=True, 
                        help='output directory for index files')
    parser.add_argument('--feather-compression', dest='feather_compression', default='lz4',
                        choices=['lz4', 'zstd', 'uncompressed'],
                        help='compression for feather files (default: lz4)')
    args = parser.parse_args()

    # check args
    if not os.path.exists(args.fasta_file):
        parser.error(f'fasta file {args.fasta_file} not found')
    if not os.path.exists(args.gtf_file):
        parser.error(f'gtf file {args.gtf_file} not found')
    # run hulkrna index
    return(hulkrna_index(args.fasta_file, args.gtf_file, args.output_dir))


def snakemake_main():
    from snakemake.script import snakemake
    # Use the snakemake object to get outputs, params, and threads
    # Run the main build function with inputs from the Snakemake rule
    return (hulkrna_index(
        fasta_file=snakemake.input.genome_fasta,
        gtf_file=snakemake.input.genes_gtf,
        output_dir=snakemake.output.index_dir
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
