"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna count   — Two-pass Bayesian read counting
    hulkrna gather  — (not yet implemented)
    hulkrna pileup  — (not yet implemented)
"""

import logging
import sys
import argparse
from pathlib import Path


def get_version() -> str:
    try:
        from . import __version__
        return f"hulkrna {__version__}"
    except ImportError:
        return "hulkrna (unknown version)"


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def index_command(args: argparse.Namespace) -> int:
    """Run the ``hulkrna index`` subcommand."""
    from .index import HulkIndex

    fasta = Path(args.fasta_file)
    gtf = Path(args.gtf_file)

    if not fasta.exists():
        sys.exit(f"Error: FASTA file not found: {fasta}")
    if not gtf.exists():
        sys.exit(f"Error: GTF file not found: {gtf}")

    HulkIndex.build(
        fasta_file=fasta,
        gtf_file=gtf,
        output_dir=args.output_dir,
        feather_compression=args.feather_compression,
        write_tsv=not args.no_tsv,
    )
    return 0


def count_command(args: argparse.Namespace) -> int:
    """Run the ``hulkrna count`` subcommand.

    Pass 1: Scan BAM, resolve fragments, train strand/insert models.
    Pass 2: Re-scan BAM, assign fractional counts using trained models.
    """
    import json
    import pysam
    from .index import HulkIndex
    from .count import pass1_learn, pass2_count

    bam_path = Path(args.bam_file)
    index_dir = Path(args.index_dir)
    output_dir = Path(args.output_dir)

    if not bam_path.exists():
        sys.exit(f"Error: BAM file not found: {bam_path}")
    if not index_dir.is_dir():
        sys.exit(f"Error: Index directory not found: {index_dir}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define output file paths
    strand_json = output_dir / "strand_model.json"
    insert_json = output_dir / "insert_size_models.json"
    t_counts_path = output_dir / "transcript_counts.feather"
    g_counts_path = output_dir / "gene_counts.feather"
    stats_path = output_dir / "stats.json"

    # Load reference index
    logging.info(f"[START] Loading index from {index_dir}")
    index = HulkIndex.load(index_dir)
    logging.info(f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
                 f"{index.num_genes} genes")

    # -- Pass 1: BAM scan + resolution + model training --------------------
    # Note: Pass 1 always excludes multimappers for model training
    bamfh = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        pass1_stats, strand_model, insert_models = pass1_learn(
            bamfh.fetch(until_eof=True),
            index,
            skip_duplicates=not args.keep_duplicates,
        )
    finally:
        bamfh.close()

    # Log pass 1 stats
    for key, val in sorted(pass1_stats.items()):
        logging.info(f"  {key}: {val:,}")

    # Write models to JSON
    strand_model.write_json(strand_json)
    insert_models.write_json(insert_json)

    # -- Pass 2: BAM scan + count assignment --------------------------------
    bamfh = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        pass2_stats, counter = pass2_count(
            bamfh.fetch(until_eof=True),
            index,
            strand_model,
            insert_models,
            skip_duplicates=not args.keep_duplicates,
            include_multimap=args.include_multimap,
        )
    finally:
        bamfh.close()

    # Write count tables
    counter.get_t_counts_df().to_feather(str(t_counts_path))
    counter.get_g_counts_df().to_feather(str(g_counts_path))
    logging.info(f"[DONE] Wrote {t_counts_path} and {g_counts_path}")

    # Write TSV mirrors if requested
    if not args.no_tsv:
        counter.get_t_counts_df().to_csv(
            t_counts_path.with_suffix(".tsv"), sep="\t", index=False
        )
        counter.get_g_counts_df().to_csv(
            g_counts_path.with_suffix(".tsv"), sep="\t", index=False
        )

    # Write combined stats
    all_stats = {"pass1": pass1_stats, "pass2": pass2_stats}
    with open(stats_path, "w") as f:
        json.dump(all_stats, f, indent=2)
    logging.info(f"  Stats written to {stats_path}")

    return 0


def gather_command(args: argparse.Namespace) -> int:
    """Placeholder for ``hulkrna gather``."""
    sys.exit("Error: 'hulkrna gather' is not yet implemented.")


def pileup_command(args: argparse.Namespace) -> int:
    """Placeholder for ``hulkrna pileup``."""
    sys.exit("Error: 'hulkrna pileup' is not yet implemented.")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="hulkrna",
        description="hulkrna: Bayesian RNA-seq read counting tool",
    )
    parser.add_argument(
        "--version", action="version", version=get_version(),
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- INDEX ---------------------------------------------------------------
    idx = subparsers.add_parser(
        "index", help="Build reference index from FASTA + GTF",
    )
    idx.add_argument(
        "--fasta", dest="fasta_file", required=True,
        help="Genome FASTA file (must be indexed with samtools faidx)",
    )
    idx.add_argument(
        "--gtf", dest="gtf_file", required=True,
        help="Gene annotation GTF file (GENCODE recommended)",
    )
    idx.add_argument(
        "-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for index files",
    )
    idx.add_argument(
        "--feather-compression", dest="feather_compression", default="lz4",
        choices=["lz4", "zstd", "uncompressed"],
        help="Compression for Feather files (default: lz4)",
    )
    idx.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV mirror files",
    )
    idx.set_defaults(func=index_command)

    # --- COUNT ---------------------------------------------------------------
    cnt = subparsers.add_parser(
        "count", help="Count reads: model training (Pass 1) + count assignment (Pass 2)",
    )
    cnt.add_argument(
        "--bam", dest="bam_file", required=True,
        help="Name-sorted or collated BAM file (STAR-aligned, with XS/NH tags)",
    )
    cnt.add_argument(
        "--index", dest="index_dir", required=True,
        help="Directory containing hulkrna index files",
    )
    cnt.add_argument(
        "-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for counts and models",
    )
    cnt.add_argument(
        "--keep-duplicates", dest="keep_duplicates",
        action="store_true", default=False,
        help="Keep reads marked as PCR/optical duplicates (default: discard)",
    )
    cnt.add_argument(
        "--include-multimap", dest="include_multimap",
        action="store_true", default=False,
        help="Include multimapping reads (NH > 1) in output (default: discard)",
    )
    cnt.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV count files",
    )
    cnt.set_defaults(func=count_command)

    # --- GATHER (stub) -------------------------------------------------------
    gth = subparsers.add_parser(
        "gather", help="Aggregate per-sample counts (not yet implemented)",
    )
    gth.set_defaults(func=gather_command)

    # --- PILEUP (stub) -------------------------------------------------------
    pup = subparsers.add_parser(
        "pileup", help="Generate coverage tracks (not yet implemented)",
    )
    pup.set_defaults(func=pileup_command)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """CLI entry point registered as ``hulkrna`` in pyproject.toml."""
    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    parser = build_parser()
    args = parser.parse_args()
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()