"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna count   — Single-pass Bayesian fragment abundance estimation
    hulkrna sim     — Generate synthetic test scenarios
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
        gtf_parse_mode=args.gtf_parse_mode,
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
    from .config import EMConfig, PipelineConfig, ScanConfig, ScoringConfig
    from .scoring import (
        overhang_alpha_to_log_penalty,
        GDNA_SPLICE_PENALTIES,
        SPLICE_UNANNOT,
    )

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

    alpha = args.em_prior_alpha
    gamma = args.em_prior_gamma
    logging.info(
        f"Using EM prior: alpha={alpha}, gamma={gamma}"
    )

    # Define output file paths
    strand_json = output_dir / "strand_model.json"
    fl_json = output_dir / "frag_length_models.json"
    counts_path = output_dir / "counts.feather"
    gene_counts_path = output_dir / "gene_counts.feather"
    detail_path = output_dir / "counts_detail.feather"
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
            "em_prior_alpha": alpha,
            "em_prior_gamma": gamma,
            "em_iterations": args.em_iterations,
            "em_convergence_delta": args.em_convergence_delta,
            "gdna_splice_penalty_unannot": args.gdna_splice_penalty_unannot,
            "confidence_threshold": args.confidence_threshold,
            "skip_duplicates": not args.keep_duplicates,
            "include_multimap": args.include_multimap,
            "sj_strand_tag": sj_strand_tag if isinstance(sj_strand_tag, str) else list(sj_strand_tag),
            "annotated_bam": args.annotated_bam,
        },
    }
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    logging.info(f"[CONFIG] Written to {config_path}")

    # -- Build pipeline config --
    overhang_log_penalty = overhang_alpha_to_log_penalty(args.overhang_alpha)
    mismatch_log_penalty = overhang_alpha_to_log_penalty(args.mismatch_alpha)

    gdna_splice_penalties = dict(GDNA_SPLICE_PENALTIES)
    gdna_splice_penalties[SPLICE_UNANNOT] = args.gdna_splice_penalty_unannot

    pipeline_config = PipelineConfig(
        em=EMConfig(
            seed=seed,
            prior_alpha=alpha,
            prior_gamma=gamma,
            iterations=args.em_iterations,
            convergence_delta=args.em_convergence_delta,
            confidence_threshold=args.confidence_threshold,
        ),
        scan=ScanConfig(
            skip_duplicates=not args.keep_duplicates,
            include_multimap=args.include_multimap,
            sj_strand_tag=sj_strand_tag,
            overlap_min_frac=args.overlap_min_frac,
        ),
        scoring=ScoringConfig(
            overhang_log_penalty=overhang_log_penalty,
            mismatch_log_penalty=mismatch_log_penalty,
            gdna_splice_penalties=gdna_splice_penalties,
        ),
        annotated_bam_path=getattr(args, 'annotated_bam', None),
    )

    # -- Run pipeline --
    result = run_pipeline(bam_path, index, config=pipeline_config)

    # Log stats
    for key, val in sorted(result.stats.to_dict().items()):
        if isinstance(val, (int, float)):
            logging.info(f"  {key}: {val:,}")

    # Write models to JSON
    result.strand_models.write_json(strand_json)
    result.frag_length_models.write_json(fl_json)

    # Write count tables
    counter = result.estimator
    counts_df = counter.get_counts_df(index)
    gene_counts_df = counter.get_gene_counts_df(index)
    detail_df = counter.get_detail_df(index)

    counts_df.to_feather(str(counts_path))
    gene_counts_df.to_feather(str(gene_counts_path))
    detail_df.to_feather(str(detail_path))
    logging.info(f"[DONE] Wrote {counts_path}, {gene_counts_path}, "
                 f"{detail_path}")

    # Write TSV mirrors if requested
    if not args.no_tsv:
        counts_df.to_csv(
            counts_path.with_suffix(".tsv"), sep="\t", index=False
        )
        gene_counts_df.to_csv(
            gene_counts_path.with_suffix(".tsv"), sep="\t", index=False
        )
        detail_df.to_csv(
            detail_path.with_suffix(".tsv"), sep="\t", index=False
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


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="hulkrna",
        description="hulkrna: Bayesian RNA-seq fragment abundance estimation",
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
    idx.add_argument(
        "--gtf-parse-mode",
        dest="gtf_parse_mode",
        default="strict",
        choices=["strict", "warn-skip"],
        help=(
            "GTF parse mode: 'strict' (default) fails on malformed lines; "
            "'warn-skip' logs warnings and skips malformed lines"
        ),
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
        help="Include multimapping reads in output (default: discard). "
             "Detected via NH tag (STAR) or secondary BAM flag (minimap2).",
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
        "--em-prior-alpha", dest="em_prior_alpha", type=float, default=0.01,
        help="Flat Dirichlet pseudocount per eligible EM component "
             "(default: 0.01).  Provides baseline regularisation "
             "independent of coverage geometry.",
    )
    cnt.add_argument(
        "--em-prior-gamma", dest="em_prior_gamma", type=float, default=1.0,
        help="OVR prior scale factor (gamma).  Default 1.0.  "
             "Controls how strongly coverage-weighted One Virtual Read "
             "priors influence the EM.",
    )
    cnt.add_argument(
        "--em-iterations", dest="em_iterations", type=int, default=1000,
        help="Maximum EM iterations for ambiguous fragment resolution "
             "(default: 1000). Set to 0 for unique-only counting.",
    )
    cnt.add_argument(
        "--em-convergence-delta", dest="em_convergence_delta",
        type=float, default=1e-6,
        help="(Advanced) Convergence threshold for EM parameter updates "
             "(default: 1e-6). Raising to 1e-5 provides a very modest "
             "speedup with negligible accuracy impact.",
    )
    from .scoring import DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
    cnt.add_argument(
        "--gdna-splice-penalty-unannot",
        dest="gdna_splice_penalty_unannot",
        type=float,
        default=DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
        help=(
            "gDNA splice penalty for SPLICED_UNANNOT fragments "
            f"(default: {DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT:g}). "
            "Higher values make gDNA more competitive for unannotated "
            "spliced reads."
        ),
    )
    cnt.add_argument(
        "--confidence-threshold", dest="confidence_threshold",
        type=float, default=0.95,
        help="Minimum RNA-normalized posterior for an EM assignment "
             "to be counted as high-confidence (default: 0.95). "
             "Affects the count_high_conf column.",
    )
    cnt.add_argument(
        "--overlap-min-frac", dest="overlap_min_frac",
        type=float, default=0.99,
        help="Min fraction of best exon overlap to retain a candidate "
             "transcript during resolution (default: 0.99 = within 1%% of best). "
             "Lower values (e.g. 0.9) keep candidates within 10%% of best.",
    )
    from .scoring import DEFAULT_OVERHANG_ALPHA, DEFAULT_MISMATCH_ALPHA
    cnt.add_argument(
        "--overhang-alpha", dest="overhang_alpha",
        type=float, default=DEFAULT_OVERHANG_ALPHA,
        help=(
            "Per-base overhang penalty alpha in [0, 1]. "
            "Each base outside the target boundary multiplies the "
            "probability by alpha. "
            f"(default: {DEFAULT_OVERHANG_ALPHA}). "
            "0 = hard binary gate, 1 = off (no penalty)."
        ),
    )
    cnt.add_argument(
        "--mismatch-alpha", dest="mismatch_alpha",
        type=float, default=DEFAULT_MISMATCH_ALPHA,
        help=(
            "Per-mismatch (NM tag) penalty alpha in [0, 1]. "
            "Each edit-distance mismatch multiplies the "
            "probability by alpha. "
            f"(default: {DEFAULT_MISMATCH_ALPHA}). "
            "0 = hard binary gate, 1 = off (no penalty)."
        ),
    )
    cnt.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV count files",
    )
    cnt.add_argument(
        "--annotated-bam", dest="annotated_bam", default=None,
        help="Write an annotated BAM with per-fragment assignment tags "
             "(ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS) to this path. "
             "Enables read-level introspection of pipeline decisions. "
             "Requires a second BAM pass (modest runtime overhead).",
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