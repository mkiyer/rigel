#!/usr/bin/env python3
"""
hulkrna_count3 - Bayesian RNA-seq read counter

Two-pass approach:
    Pass 1: Learn strand protocol + collect initial abundance estimates
    Pass 2: Bayesian assignment using learned strand specificity + abundance priors

Unlike featureCounts/htseq-count which discard reads overlapping genes on
opposite strands, this tool uses Bayesian inference to probabilistically
assign reads based on:
    - Library strand specificity (learned from spliced reads)
    - Relative gene abundance (learned from Pass 1)
    - Splice junction concordance
"""

import os
import sys
import argparse
import logging
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

# import local modules
SCRIPT_PATH = Path(__file__).resolve().parent
HULKRNA_LIB_PATH = SCRIPT_PATH / '..' / 'lib'
if HULKRNA_LIB_PATH.exists() and str(HULKRNA_LIB_PATH) not in sys.path:
    sys.path.insert(0, str(HULKRNA_LIB_PATH))

from hulkrna.txref import TxIndex
from hulkrna.bam import parse_bam_file
from hulkrna.fragment import Fragment
from hulkrna.strand_model import StrandModel, ProtocolType
from hulkrna.bayes_count import BayesCounter

# setup logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# Output file names
TRANSCRIPT_COUNTS_TSV_GZ = 'transcript_counts.tsv.gz'
TRANSCRIPT_COUNTS_FEATHER = 'transcript_counts.feather'
GENE_COUNTS_TSV_GZ = 'gene_counts.tsv.gz'
GENE_COUNTS_FEATHER = 'gene_counts.feather'
SUMMARY_JSON = 'summary.json'
BAM_STATS_JSON = 'bam_stats.json'
STRAND_MODEL_JSON = 'strand_model.json'
MULTIMAP_BAM = 'multimap.bam'


def run_pass1(input_bam_file, txindex, ignore_duplicates, threads):
    """
    Pass 1: Learn strand protocol and collect initial abundance estimates.

    Iterates through the BAM file once to:
    1. Profile the strand protocol from spliced reads (StrandModel)
    2. Assign reads uniformly to overlapping genes to get initial abundances

    Parameters
    ----------
    input_bam_file : str
        Path to input BAM file.
    txindex : TxIndex
        Loaded transcript/gene index.
    ignore_duplicates : bool
        Whether to skip duplicate-flagged reads.
    threads : int
        BAM decompression threads.

    Returns
    -------
    strand_result : StrandModelResult
        Learned strand model with p_sense posterior.
    g_abundance : np.ndarray
        Gene abundance vector (total raw counts per gene).
    bam_stats : dict
        BAM parsing statistics.
    """
    logger.info('[PASS 1] Starting: learn strand + initial abundances')

    # Initialize strand model
    strand_model = StrandModel()

    # Initialize a BayesCounter with uniform prior for initial counting
    # We'll use a dummy strand result (unstranded) for Pass 1
    from hulkrna.strand_model import StrandModelResult
    dummy_strand = StrandModelResult(
        protocol=ProtocolType.UNSTRANDED,
        alpha=1.0,
        beta=1.0,
    )
    pass1_counter = BayesCounter(txindex, dummy_strand)

    # Open BAM file
    bam_input = sys.stdin if input_bam_file == '-' else input_bam_file
    bam_stats = {}

    with pysam.AlignmentFile(bam_input, "rb", threads=threads) as bamfh:
        i = 0
        for r1, r2 in parse_bam_file(bamfh, bam_stats, 
                                      skip_duplicates=ignore_duplicates):
            # Feed to strand model (learns from spliced reads)
            strand_model.observe(r1, r2)

            # Build fragment and count
            fragment = Fragment.from_reads(r1, r2)
            pass1_counter.count_fragment(fragment, ignore_duplicates)

            i += 1
            if i % 500000 == 0:
                logger.debug(f'Pass 1: {i} fragments processed')

    logger.info(f'[PASS 1] Done: {i} fragments processed')

    # Finalize strand model
    strand_result = strand_model.finalize(min_spliced_reads=100)
    logger.info(f'[PASS 1] Strand model: {strand_result.to_dict()}')

    # Extract gene abundance from Pass 1
    g_abundance = pass1_counter.get_gene_abundance_vector()
    n_genes_expressed = np.sum(g_abundance > 0)
    logger.info(f'[PASS 1] Genes with >0 counts: {n_genes_expressed} / {len(g_abundance)}')

    return strand_result, g_abundance, bam_stats


def run_pass2(input_bam_file, txindex, strand_result, g_abundance,
              ignore_duplicates, threads, output_dir, multimap_bam=True):
    """
    Pass 2: Bayesian assignment using learned strand + abundance priors.

    Parameters
    ----------
    input_bam_file : str
        Path to input BAM file.
    txindex : TxIndex
        Loaded transcript/gene index.
    strand_result : StrandModelResult
        Learned strand model from Pass 1.
    g_abundance : np.ndarray
        Gene abundance prior from Pass 1.
    ignore_duplicates : bool
        Whether to skip duplicate-flagged reads.
    threads : int
        BAM decompression threads.
    output_dir : Path
        Output directory.
    multimap_bam : bool
        Whether to write multimapping reads to a separate BAM file.

    Returns
    -------
    counter : BayesCounter
        Counter with final Bayesian counts.
    bam_stats : dict
        BAM parsing statistics for Pass 2.
    """
    logger.info('[PASS 2] Starting: Bayesian assignment')

    # Initialize Bayesian counter with learned strand + abundance
    counter = BayesCounter(
        txindex, strand_result,
        abundance_prior=g_abundance,
        abundance_pseudocount=1.0,
    )

    bam_input = sys.stdin if input_bam_file == '-' else input_bam_file
    bam_stats = {}

    multimap_bam_path = output_dir / MULTIMAP_BAM if multimap_bam else None

    with pysam.AlignmentFile(bam_input, "rb", threads=threads) as bamfh:
        mm_bamfh = None
        if multimap_bam_path:
            mm_bamfh = pysam.AlignmentFile(
                str(multimap_bam_path), "wb", template=bamfh
            )

        try:
            i = 0
            for r1, r2 in parse_bam_file(
                bamfh, bam_stats,
                skip_duplicates=ignore_duplicates,
                multimap_bamfh=mm_bamfh,
            ):
                fragment = Fragment.from_reads(r1, r2)
                counter.count_fragment(fragment, ignore_duplicates)

                i += 1
                if i % 500000 == 0:
                    logger.debug(
                        f'Pass 2: {i} fragments '
                        f'(exonic={counter.stats.n_exonic} '
                        f'intronic={counter.stats.n_intronic} '
                        f'intergenic={counter.stats.n_intergenic} '
                        f'multi_gene={counter.stats.n_multi_gene} '
                        f'opp_strand={counter.stats.n_opposite_strand_genes})'
                    )
                    logger.debug(f'Counts summary:\n{counter.summary_by_type()}')
        finally:
            if mm_bamfh is not None:
                mm_bamfh.close()

    logger.info(f'[PASS 2] Done: {i} fragments processed')
    logger.info(f'[PASS 2] Stats: {counter.stats.to_dict()}')

    return counter, bam_stats


def write_outputs(output_dir, txindex, counter, strand_result,
                  bam_stats_pass1, bam_stats_pass2, feather_compression):
    """Write all output files."""
    output_dir = Path(output_dir)
    feather_kwargs = {'compression': feather_compression}

    # Write strand model
    logger.info('Writing strand model')
    strand_file = output_dir / STRAND_MODEL_JSON
    with open(strand_file, 'w') as f:
        json.dump(strand_result.to_dict(), f, indent=2)

    # Write BAM stats (combine pass1 and pass2)
    logger.info('Writing BAM stats')
    combined_bam_stats = {
        'pass1': bam_stats_pass1,
        'pass2': bam_stats_pass2,
        'strand_model': strand_result.to_dict(),
    }
    with open(output_dir / BAM_STATS_JSON, 'w') as f:
        json.dump(combined_bam_stats, f, indent=2)

    # Write summary
    logger.info('Writing count summary')
    with open(output_dir / SUMMARY_JSON, 'w') as f:
        json.dump(counter.summary(), f, indent=2)

    # Write transcript counts
    logger.info('Writing transcript counts')
    df = pd.concat([txindex.t_df, counter.get_t_counts_df()], axis=1)
    df.to_csv(output_dir / TRANSCRIPT_COUNTS_TSV_GZ,
              sep='\t', index=False, header=True)
    df.to_feather(output_dir / TRANSCRIPT_COUNTS_FEATHER, **feather_kwargs)

    # Write gene counts
    logger.info('Writing gene counts')
    df = pd.concat([txindex.g_df, counter.get_g_counts_df()], axis=1)
    columns_to_output = [x for x in df.columns if x not in ['t_index', 't_id']]
    df = df[columns_to_output]
    df.to_csv(output_dir / GENE_COUNTS_TSV_GZ,
              sep='\t', index=False, header=True)
    df.to_feather(output_dir / GENE_COUNTS_FEATHER, **feather_kwargs)

    logger.info('All outputs written')


def hulkrna_count(
    input_bam_file,
    index_dir,
    output_dir,
    ignore_duplicates=True,
    threads=1,
    feather_compression='lz4',
):
    """
    Main entry point for Bayesian RNA-seq counting.

    Two-pass approach:
        Pass 1: Learn strand protocol + initial gene abundances
        Pass 2: Bayesian assignment with strand + abundance priors
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Load index
    logger.info(f'Loading index from {index_dir}')
    txindex = TxIndex.load(index_dir)
    logger.info(f'Index loaded: {txindex.num_transcripts} transcripts, '
                f'{txindex.num_genes} genes')

    # ---- PASS 1 ----
    strand_result, g_abundance, bam_stats_pass1 = run_pass1(
        input_bam_file, txindex, ignore_duplicates, threads
    )

    # ---- PASS 2 ----
    counter, bam_stats_pass2 = run_pass2(
        input_bam_file, txindex, strand_result, g_abundance,
        ignore_duplicates, threads, output_dir
    )

    # ---- OUTPUTS ----
    write_outputs(
        output_dir, txindex, counter, strand_result,
        bam_stats_pass1, bam_stats_pass2, feather_compression
    )

    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Bayesian RNA-seq read counter with strand-aware '
                    'probabilistic gene assignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        '-@', '--threads', type=int, default=1,
        help='threads for BAM decompression',
    )
    parser.add_argument(
        '--feather-compression', dest='feather_compression', default='lz4',
        choices=['lz4', 'zstd', 'uncompressed'],
        help='compression for feather output files',
    )
    parser.add_argument(
        '--nodup', dest='ignore_duplicates', default=False,
        action='store_true',
        help='do not count duplicate-flagged reads',
    )
    parser.add_argument(
        '--index', dest='index_dir', required=True,
        help='directory with hulkrna index files',
    )
    parser.add_argument(
        '-i', '--input', dest='bam_file', required=True,
        help='input BAM file ("-" for stdin)',
    )
    parser.add_argument(
        '-o', '--output', dest='output_dir', required=True,
        help='output directory',
    )
    args = parser.parse_args()

    # Validate arguments
    if not os.path.exists(args.index_dir):
        parser.error(f'index directory {args.index_dir} not found')
    if args.bam_file != '-' and not os.path.exists(args.bam_file):
        parser.error(f'input {args.bam_file} not found')

    return hulkrna_count(
        input_bam_file=args.bam_file,
        index_dir=args.index_dir,
        output_dir=args.output_dir,
        ignore_duplicates=args.ignore_duplicates,
        threads=args.threads,
        feather_compression=args.feather_compression,
    )


def snakemake_main():
    from snakemake.script import snakemake
    return hulkrna_count(
        input_bam_file=snakemake.input.bam,
        index_dir=snakemake.input.index_dir,
        output_dir=snakemake.params.output_dir,
        ignore_duplicates=snakemake.params.ignore_duplicates,
        threads=snakemake.threads,
        feather_compression=snakemake.params.feather_compression,
    )


def is_snakemake_run():
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
