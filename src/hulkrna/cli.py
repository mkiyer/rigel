"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna count   — Single-pass Bayesian read counting
    hulkrna sim     — Generate synthetic test scenarios
    hulkrna gather  — (not yet implemented)
    hulkrna pileup  — (not yet implemented)
"""

import logging
import sys
import time
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

    Single-pass pipeline: scan BAM, resolve fragments, train
    strand/insert models, then assign counts via unified EM.
    """
    import json
    from .index import HulkIndex
    from .pipeline import run_pipeline

    bam_path = Path(args.bam_file)
    index_dir = Path(args.index_dir)
    output_dir = Path(args.output_dir)

    if not bam_path.exists():
        sys.exit(f"Error: BAM file not found: {bam_path}")
    if not index_dir.is_dir():
        sys.exit(f"Error: Index directory not found: {index_dir}")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Generate seed from time if not provided
    if args.seed is None:
        seed = int(time.time())
        logging.info(f"No seed provided, using timestamp: {seed}")
    else:
        seed = args.seed
        logging.info(f"Using provided seed: {seed}")

    alpha = args.alpha
    logging.info(f"Using Dirichlet pseudocount alpha: {alpha}")

    # Define output file paths
    strand_json = output_dir / "strand_model.json"
    insert_json = output_dir / "insert_size_models.json"
    t_counts_path = output_dir / "transcript_counts.feather"
    g_counts_path = output_dir / "gene_counts.feather"
    t_sparse_path = output_dir / "transcript_counts_sparse.feather"
    g_sparse_path = output_dir / "gene_counts_sparse.feather"
    stats_path = output_dir / "stats.json"
    config_path = output_dir / "config.json"

    # Load reference index
    logging.info(f"[START] Loading index from {index_dir}")
    index = HulkIndex.load(index_dir)
    logging.info(f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
                 f"{index.num_genes} genes")

    # Resolve sj_strand_tag from CLI list to str | tuple
    sj_tag_list = args.sj_strand_tag
    if len(sj_tag_list) == 1:
        sj_strand_tag = sj_tag_list[0]              # "auto", "XS", or "ts"
    else:
        sj_strand_tag = tuple(sj_tag_list)           # ("XS", "ts")

    # Write config file with all parameters
    config = {
        "command": "hulkrna count",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "parameters": {
            "bam_file": str(bam_path.resolve()),
            "index_dir": str(index_dir.resolve()),
            "output_dir": str(output_dir.resolve()),
            "seed": seed,
            "alpha": alpha,
            "em_iterations": args.em_iterations,
            "skip_duplicates": not args.keep_duplicates,
            "include_multimap": args.include_multimap,
            "sj_strand_tag": sj_strand_tag if isinstance(sj_strand_tag, str) else list(sj_strand_tag),
        },
    }
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    logging.info(f"[CONFIG] Written to {config_path}")

    # -- Run pipeline --
    result = run_pipeline(
        bam_path,
        index,
        skip_duplicates=not args.keep_duplicates,
        include_multimap=args.include_multimap,
        sj_strand_tag=sj_strand_tag,
        seed=seed,
        alpha=alpha,        em_iterations=args.em_iterations,    )

    # Log stats
    for key, val in sorted(result.stats.to_dict().items()):
        if isinstance(val, (int, float)):
            logging.info(f"  {key}: {val:,}")

    # Write models to JSON
    result.strand_models.write_json(strand_json)
    result.insert_models.write_json(insert_json)

    # Write count tables
    # Write count tables (wide format)
    result.counter.get_t_counts_df().to_feather(str(t_counts_path))
    result.counter.get_g_counts_df(index.t_to_g_arr).to_feather(str(g_counts_path))
    logging.info(f"[DONE] Wrote {t_counts_path} and {g_counts_path}")

    # Write sparse count tables (long format with provenance)
    result.counter.get_sparse_t_counts_df().to_feather(str(t_sparse_path))
    result.counter.get_sparse_g_counts_df(index.t_to_g_arr).to_feather(str(g_sparse_path))
    logging.info(f"[DONE] Wrote {t_sparse_path} and {g_sparse_path}")

    # Write TSV mirrors if requested
    if not args.no_tsv:
        result.counter.get_t_counts_df().to_csv(
            t_counts_path.with_suffix(".tsv"), sep="\t", index=False
        )
        result.counter.get_g_counts_df(index.t_to_g_arr).to_csv(
            g_counts_path.with_suffix(".tsv"), sep="\t", index=False
        )
        result.counter.get_sparse_t_counts_df().to_csv(
            t_sparse_path.with_suffix(".tsv"), sep="\t", index=False
        )
        result.counter.get_sparse_g_counts_df(index.t_to_g_arr).to_csv(
            g_sparse_path.with_suffix(".tsv"), sep="\t", index=False
        )

    # Write stats
    with open(stats_path, "w") as f:
        json.dump(result.stats.to_dict(), f, indent=2)
    logging.info(f"  Stats written to {stats_path}")

    return 0


def sim_command(args: argparse.Namespace) -> int:
    """Run the ``hulkrna sim`` subcommand."""
    import yaml
    from .sim import Scenario, SimConfig

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load scenario config from YAML
    with open(args.config) as f:
        cfg = yaml.safe_load(f)

    genome_length = cfg.get("genome_length", args.genome_length)
    seed = cfg.get("seed", args.seed)
    ref_name = cfg.get("ref_name", "chr1")

    sim_config = SimConfig(
        frag_mean=cfg.get("frag_mean", 250),
        frag_std=cfg.get("frag_std", 50),
        frag_min=cfg.get("frag_min", 50),
        frag_max=cfg.get("frag_max", 1000),
        read_length=cfg.get("read_length", 150),
        error_rate=cfg.get("error_rate", 0.0),
        seed=seed,
    )

    sc = Scenario(
        name=cfg.get("name", "scenario"),
        genome_length=genome_length,
        seed=seed,
        work_dir=output_dir,
        ref_name=ref_name,
    )

    # Add genes from config
    for gene in cfg.get("genes", []):
        sc.add_gene(
            gene_id=gene["gene_id"],
            strand=gene["strand"],
            transcripts=gene["transcripts"],
            gene_name=gene.get("gene_name"),
            gene_type=gene.get("gene_type", "protein_coding"),
        )

    result = sc.build(
        n_fragments=cfg.get("n_fragments", args.num_reads),
        sim_config=sim_config,
    )

    logging.info(f"Scenario artifacts written to {output_dir}")
    logging.info(f"  FASTA: {result.fasta_path}")
    logging.info(f"  GTF:   {result.gtf_path}")
    logging.info(f"  BAM:   {result.bam_path}")
    logging.info(f"  Index: {result.index_dir}")
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
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Enable verbose (DEBUG-level) logging",
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
        help="Name-sorted or collated BAM file (with NH tag)",
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
        "--sj-strand-tag", dest="sj_strand_tag",
        nargs="+", default=["auto"],
        help="BAM tag(s) for splice-junction strand. 'auto' (default) "
             "detects the tag from the BAM. Use 'XS' for STAR, 'ts' for "
             "minimap2, or list multiple to check in order (e.g. XS ts).",
    )
    cnt.add_argument(
        "--seed", dest="seed", type=int, default=None,
        help="Random seed for reproducibility (default: use current timestamp)",
    )
    cnt.add_argument(
        "--alpha", dest="alpha", type=float, default=1.0,
        help="Dirichlet pseudocount for prior smoothing (default: 1.0)",
    )
    cnt.add_argument(
        "--em-iterations", dest="em_iterations", type=int, default=10,
        help="Maximum EM iterations for ambiguous fragment resolution "
             "(default: 10). Set to 0 for unique-only priors.",
    )
    cnt.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV count files",
    )
    cnt.set_defaults(func=count_command)

    # --- SIM -----------------------------------------------------------------
    sim = subparsers.add_parser(
        "sim", help="Generate synthetic test scenarios",
    )
    sim.add_argument(
        "--config", dest="config", required=True,
        help="YAML configuration file defining the scenario",
    )
    sim.add_argument(
        "-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for scenario artifacts",
    )
    sim.add_argument(
        "--genome-length", dest="genome_length", type=int, default=5000,
        help="Genome length in bases (default: 5000, overridden by YAML)",
    )
    sim.add_argument(
        "--seed", dest="seed", type=int, default=42,
        help="Random seed (default: 42, overridden by YAML)",
    )
    sim.add_argument(
        "--num-reads", dest="num_reads", type=int, default=1000,
        help="Number of fragments to simulate (default: 1000, overridden by YAML)",
    )
    sim.set_defaults(func=sim_command)

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
    parser = build_parser()
    args = parser.parse_args()

    log_level = logging.DEBUG if getattr(args, 'verbose', False) else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()