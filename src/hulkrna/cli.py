"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna quant   — Single-pass Bayesian fragment abundance estimation
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
    from .index import TranscriptIndex

    fasta = Path(args.fasta_file)
    gtf = Path(args.gtf_file)

    if not fasta.exists():
        sys.exit(f"Error: FASTA file not found: {fasta}")
    if not gtf.exists():
        sys.exit(f"Error: GTF file not found: {gtf}")

    TranscriptIndex.build(
        fasta_file=fasta,
        gtf_file=gtf,
        output_dir=args.output_dir,
        feather_compression=args.feather_compression,
        write_tsv=not args.no_tsv,
        gtf_parse_mode=args.gtf_parse_mode,
    )
    return 0


def quant_command(args: argparse.Namespace) -> int:
    """Run the ``hulkrna quant`` subcommand.

    Single-pass pipeline: scan BAM, resolve fragments, train
    strand/insert models, then quantify via unified EM.
    """
    import json
    from .index import TranscriptIndex
    from .pipeline import run_pipeline
    from .config import EMConfig, PipelineConfig, BamScanConfig, FragmentScoringConfig
    from .scoring import (
        overhang_alpha_to_log_penalty,
        GDNA_SPLICE_PENALTIES,
    )
    from .splice import SPLICE_UNANNOT

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
    quant_path = output_dir / "quant.feather"
    gene_quant_path = output_dir / "gene_quant.feather"
    detail_path = output_dir / "quant_detail.feather"
    stats_path = output_dir / "stats.json"
    config_path = output_dir / "config.json"

    # Load reference index
    logging.info(f"[START] Loading index from {index_dir}")
    index = TranscriptIndex.load(index_dir)
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
        "command": "hulkrna quant",
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
        scan=BamScanConfig(
            skip_duplicates=not args.keep_duplicates,
            include_multimap=args.include_multimap,
            sj_strand_tag=sj_strand_tag,
        ),
        scoring=FragmentScoringConfig(
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

    # Write quant tables
    estimator = result.estimator
    quant_df = estimator.get_counts_df(index)
    gene_quant_df = estimator.get_gene_counts_df(index)
    detail_df = estimator.get_detail_df(index)

    quant_df.to_feather(str(quant_path))
    gene_quant_df.to_feather(str(gene_quant_path))
    detail_df.to_feather(str(detail_path))
    logging.info(f"[DONE] Wrote {quant_path}, {gene_quant_path}, "
                 f"{detail_path}")

    # Write TSV mirrors if requested
    if not args.no_tsv:
        quant_df.to_csv(
            quant_path.with_suffix(".tsv"), sep="\t", index=False
        )
        gene_quant_df.to_csv(
            gene_quant_path.with_suffix(".tsv"), sep="\t", index=False
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
        **{k: cfg[k] for k in (
            "frag_mean", "frag_std", "frag_min", "frag_max",
            "read_length", "error_rate",
        ) if k in cfg},
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

    # --- QUANT ---------------------------------------------------------------
    quant_parser = subparsers.add_parser(
        "quant", help="Quantify transcripts: model training (Pass 1) + abundance estimation (Pass 2)",
    )
    quant_parser.add_argument(
        "--bam", dest="bam_file", required=True,
        help="Name-sorted or collated BAM file (with NH tag)",
    )
    quant_parser.add_argument(
        "--index", dest="index_dir", required=True,
        help="Directory containing hulkrna index files",
    )
    quant_parser.add_argument(
        "-o", "--output-dir", dest="output_dir", required=True,
        help="Output directory for quantification results and models",
    )
    quant_parser.add_argument(
        "--keep-duplicates", dest="keep_duplicates",
        action="store_true", default=False,
        help="Keep reads marked as PCR/optical duplicates (default: discard)",
    )
    quant_parser.add_argument(
        "--include-multimap", dest="include_multimap",
        action="store_true", default=False,
        help="Include multimapping reads in output (default: discard). "
             "Detected via NH tag (STAR) or secondary BAM flag (minimap2).",
    )
    quant_parser.add_argument(
        "--sj-strand-tag", dest="sj_strand_tag",
        nargs="+", default=["auto"],
        help="BAM tag(s) for splice-junction strand. 'auto' (default) "
             "detects the tag from the BAM. Use 'XS' for STAR, 'ts' for "
             "minimap2, or list multiple to check in order (e.g. XS ts).",
    )
    quant_parser.add_argument(
        "--seed", dest="seed", type=int, default=None,
        help="Random seed for reproducibility (default: use current timestamp)",
    )
    quant_parser.add_argument(
        "--em-prior-alpha", dest="em_prior_alpha", type=float, default=0.01,
        help="Flat Dirichlet pseudocount per eligible EM component "
             "(default: 0.01).  Provides baseline regularisation "
             "independent of coverage geometry.",
    )
    quant_parser.add_argument(
        "--em-prior-gamma", dest="em_prior_gamma", type=float, default=1.0,
        help="OVR prior scale factor (gamma).  Default 1.0.  "
             "Controls how strongly coverage-weighted One Virtual Read "
             "priors influence the EM.",
    )
    quant_parser.add_argument(
        "--em-iterations", dest="em_iterations", type=int, default=1000,
        help="Maximum EM iterations for ambiguous fragment resolution "
             "(default: 1000). Set to 0 for unique-only quantification.",
    )
    quant_parser.add_argument(
        "--em-convergence-delta", dest="em_convergence_delta",
        type=float, default=1e-6,
        help="(Advanced) Convergence threshold for EM parameter updates "
             "(default: 1e-6). Raising to 1e-5 provides a very modest "
             "speedup with negligible accuracy impact.",
    )
    from .scoring import DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
    quant_parser.add_argument(
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
    quant_parser.add_argument(
        "--confidence-threshold", dest="confidence_threshold",
        type=float, default=0.95,
        help="Minimum RNA-normalized posterior for an EM assignment "
             "to be considered high-confidence (default: 0.95). "
             "Affects the count_high_conf column.",
    )
    from .scoring import DEFAULT_OVERHANG_ALPHA, DEFAULT_MISMATCH_ALPHA
    quant_parser.add_argument(
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
    quant_parser.add_argument(
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
    quant_parser.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=False,
        help="Skip writing human-readable TSV quant files",
    )
    quant_parser.add_argument(
        "--annotated-bam", dest="annotated_bam", default=None,
        help="Write an annotated BAM with per-fragment assignment tags "
             "(ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS) to this path. "
             "Enables read-level introspection of pipeline decisions. "
             "Requires a second BAM pass (modest runtime overhead).",
    )
    quant_parser.set_defaults(func=quant_command)

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