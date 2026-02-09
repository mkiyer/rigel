"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna count   — (not yet implemented)
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
    """Run the ``hulkrna count`` subcommand (Pass 1: BAM scan + index query)."""
    import pysam
    from .index import HulkIndex
    from .query_bam import query_bam

    bam_path = Path(args.bam_file)
    index_dir = Path(args.index_dir)
    output = Path(args.output)

    if not bam_path.exists():
        sys.exit(f"Error: BAM file not found: {bam_path}")
    if not index_dir.is_dir():
        sys.exit(f"Error: Index directory not found: {index_dir}")

    # Load reference index
    logging.info(f"[START] Loading index from {index_dir}")
    index = HulkIndex.load(index_dir)
    logging.info(f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
                 f"{index.num_genes} genes")

    # Open BAM and optional multimap output
    multimap_bamfh = None
    bamfh = pysam.AlignmentFile(str(bam_path), "rb")
    try:
        if args.multimap_bam:
            multimap_bamfh = pysam.AlignmentFile(
                args.multimap_bam, "wb", template=bamfh
            )

        # Run Pass 1: scan BAM, query index, stream hits to Feather
        stats = query_bam(
            bamfh.fetch(until_eof=True),
            index,
            output,
            skip_duplicates=not args.keep_duplicates,
            multimap_bamfh=multimap_bamfh,
            compression=args.feather_compression,
            write_tsv=not args.no_tsv,
        )

        # Log stats
        for key, val in sorted(stats.items()):
            logging.info(f"  {key}: {val:,}")

    finally:
        bamfh.close()
        if multimap_bamfh is not None:
            multimap_bamfh.close()

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
        "count", help="Scan BAM and query reference index (Pass 1)",
    )
    cnt.add_argument(
        "--bam", dest="bam_file", required=True,
        help="Coordinate-sorted BAM file (STAR-aligned, with XS/NH tags)",
    )
    cnt.add_argument(
        "--index", dest="index_dir", required=True,
        help="Directory containing hulkrna index files",
    )
    cnt.add_argument(
        "-o", "--output", dest="output", required=True,
        help="Output path for intermediate hit Feather file",
    )
    cnt.add_argument(
        "--keep-duplicates", dest="keep_duplicates",
        action="store_true", default=False,
        help="Keep reads marked as PCR/optical duplicates (default: discard)",
    )
    cnt.add_argument(
        "--multimap-bam", dest="multimap_bam", default=None,
        help="Write multimapping reads (NH > 1) to this BAM file",
    )
    cnt.add_argument(
        "--feather-compression", dest="feather_compression", default="lz4",
        choices=["lz4", "zstd", "uncompressed"],
        help="Compression for Feather output (default: lz4)",
    )
    cnt.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV mirror file",
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