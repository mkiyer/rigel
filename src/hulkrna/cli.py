"""
hulkrna.cli — Unified command-line interface.

Entry point: ``hulkrna`` (registered in pyproject.toml).

Subcommands:
    hulkrna index   — Build reference index from FASTA + GTF
    hulkrna quant   — Single-pass Bayesian fragment abundance estimation
    hulkrna sim     — Generate synthetic test scenarios
"""

import argparse
import logging
import sys
import time
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
    import yaml
    from .index import TranscriptIndex
    from .pipeline import run_pipeline
    from .config import EMConfig, PipelineConfig, BamScanConfig, FragmentScoringConfig
    from .scoring import (
        overhang_alpha_to_log_penalty,
        GDNA_SPLICE_PENALTIES,
        DEFAULT_OVERHANG_ALPHA,
        DEFAULT_MISMATCH_ALPHA,
        DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    )
    from .splice import SPLICE_UNANNOT

    # -- Resolve parameters: CLI > YAML > hardcoded defaults --
    _DEFAULTS = {
        "seed": None,
        "em_prior_alpha": 0.01,
        "em_prior_gamma": 1.0,
        "em_iterations": 1000,
        "em_convergence_delta": 1e-6,
        "prune_threshold": 0.1,
        "confidence_threshold": 0.95,
        "include_multimap": True,
        "keep_duplicates": False,
        "no_tsv": False,
        "sj_strand_tag": ["auto"],
        "annotated_bam": None,
        "overhang_alpha": DEFAULT_OVERHANG_ALPHA,
        "mismatch_alpha": DEFAULT_MISMATCH_ALPHA,
        "gdna_splice_penalty_unannot": DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
        "tss_window": 200,
        "nrna_frac_kappa_global": None,
        "nrna_frac_kappa_locus": None,
        "nrna_frac_kappa_tss": None,
        "nrna_frac_mom_min_evidence_global": 50.0,
        "nrna_frac_mom_min_evidence_locus": 30.0,
        "nrna_frac_mom_min_evidence_tss": 20.0,
        "nrna_frac_kappa_min": 2.0,
        "nrna_frac_kappa_max": 200.0,
        "nrna_frac_kappa_fallback": 5.0,
        "nrna_frac_kappa_min_obs": 20,
        "gdna_kappa_chrom": None,
        "gdna_kappa_locus": None,
        "gdna_mom_min_evidence_chrom": 50.0,
        "gdna_mom_min_evidence_locus": 30.0,
        "strand_prior_kappa": 2.0,
    }
    _resolve_quant_args(args, _DEFAULTS)

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

    # Normalise sj_strand_tag from YAML (may be str)
    if isinstance(args.sj_strand_tag, str):
        args.sj_strand_tag = [args.sj_strand_tag]

    # Define output file paths
    quant_path = output_dir / "quant.feather"
    gene_quant_path = output_dir / "gene_quant.feather"
    loci_path = output_dir / "loci.feather"
    detail_path = output_dir / "quant_detail.feather"
    summary_path = output_dir / "summary.json"
    config_path = output_dir / "config.json"

    # Load reference index
    logging.info(f"[START] Loading index from {index_dir}")
    index = TranscriptIndex.load(index_dir)
    logging.info(f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
                 f"{index.num_genes} genes")

    # Resolve sj_strand_tag from list to str | tuple
    sj_tag_list = args.sj_strand_tag
    if len(sj_tag_list) == 1:
        sj_strand_tag = sj_tag_list[0]              # "auto", "XS", or "ts"
    else:
        sj_strand_tag = tuple(sj_tag_list)           # ("XS", "ts")

    # Write config file with all parameters
    config = {
        "command": "hulkrna quant",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "source_config": str(Path(args.config).resolve()) if args.config else None,
        "parameters": {
            "bam_file": str(bam_path.resolve()),
            "index_dir": str(index_dir.resolve()),
            "output_dir": str(output_dir.resolve()),
            "seed": seed,
            "em_prior_alpha": alpha,
            "em_prior_gamma": gamma,
            "em_iterations": args.em_iterations,
            "em_convergence_delta": args.em_convergence_delta,
            "prune_threshold": args.prune_threshold,
            "gdna_splice_penalty_unannot": args.gdna_splice_penalty_unannot,
            "confidence_threshold": args.confidence_threshold,
            "skip_duplicates": not args.keep_duplicates,
            "include_multimap": args.include_multimap,
            "sj_strand_tag": sj_strand_tag if isinstance(sj_strand_tag, str) else list(sj_strand_tag),
            "annotated_bam": args.annotated_bam,
            "tss_window": args.tss_window,
            "nrna_frac_kappa_global": args.nrna_frac_kappa_global,
            "nrna_frac_kappa_locus": args.nrna_frac_kappa_locus,
            "nrna_frac_kappa_tss": args.nrna_frac_kappa_tss,
            "nrna_frac_mom_min_evidence_global": args.nrna_frac_mom_min_evidence_global,
            "nrna_frac_mom_min_evidence_locus": args.nrna_frac_mom_min_evidence_locus,
            "nrna_frac_mom_min_evidence_tss": args.nrna_frac_mom_min_evidence_tss,
            "nrna_frac_kappa_min": args.nrna_frac_kappa_min,
            "nrna_frac_kappa_max": args.nrna_frac_kappa_max,
            "nrna_frac_kappa_fallback": args.nrna_frac_kappa_fallback,
            "nrna_frac_kappa_min_obs": args.nrna_frac_kappa_min_obs,
            "gdna_kappa_chrom": args.gdna_kappa_chrom,
            "gdna_kappa_locus": args.gdna_kappa_locus,
            "gdna_mom_min_evidence_chrom": args.gdna_mom_min_evidence_chrom,
            "gdna_mom_min_evidence_locus": args.gdna_mom_min_evidence_locus,
            "strand_prior_kappa": args.strand_prior_kappa,
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
            prune_threshold=args.prune_threshold,
            confidence_threshold=args.confidence_threshold,
            tss_window=args.tss_window,
            nrna_frac_kappa_global=args.nrna_frac_kappa_global,
            nrna_frac_kappa_locus=args.nrna_frac_kappa_locus,
            nrna_frac_kappa_tss=args.nrna_frac_kappa_tss,
            nrna_frac_mom_min_evidence_global=args.nrna_frac_mom_min_evidence_global,
            nrna_frac_mom_min_evidence_locus=args.nrna_frac_mom_min_evidence_locus,
            nrna_frac_mom_min_evidence_tss=args.nrna_frac_mom_min_evidence_tss,
            nrna_frac_kappa_min=args.nrna_frac_kappa_min,
            nrna_frac_kappa_max=args.nrna_frac_kappa_max,
            nrna_frac_kappa_fallback=args.nrna_frac_kappa_fallback,
            nrna_frac_kappa_min_obs=args.nrna_frac_kappa_min_obs,
            gdna_kappa_chrom=args.gdna_kappa_chrom,
            gdna_kappa_locus=args.gdna_kappa_locus,
            gdna_mom_min_evidence_chrom=args.gdna_mom_min_evidence_chrom,
            gdna_mom_min_evidence_locus=args.gdna_mom_min_evidence_locus,
        ),
        scan=BamScanConfig(
            skip_duplicates=not args.keep_duplicates,
            include_multimap=args.include_multimap,
            sj_strand_tag=sj_strand_tag,
            strand_prior_kappa=args.strand_prior_kappa,
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

    # Write quant tables
    estimator = result.estimator
    quant_df = estimator.get_counts_df(index)
    gene_quant_df = estimator.get_gene_counts_df(index)
    loci_df = estimator.get_loci_df()
    detail_df = estimator.get_detail_df(index)

    quant_df.to_feather(str(quant_path))
    gene_quant_df.to_feather(str(gene_quant_path))
    loci_df.to_feather(str(loci_path))
    detail_df.to_feather(str(detail_path))
    logging.info(f"[DONE] Wrote {quant_path}, {gene_quant_path}, "
                 f"{loci_path}, {detail_path}")

    # Write TSV mirrors if requested
    if not args.no_tsv:
        quant_df.to_csv(
            quant_path.with_suffix(".tsv"), sep="\t", index=False
        )
        gene_quant_df.to_csv(
            gene_quant_path.with_suffix(".tsv"), sep="\t", index=False
        )
        loci_df.to_csv(
            loci_path.with_suffix(".tsv"), sep="\t", index=False
        )
        detail_df.to_csv(
            detail_path.with_suffix(".tsv"), sep="\t", index=False
        )

    # Build and write summary.json
    stats = result.stats
    sm = result.strand_models
    flm = result.frag_length_models
    total_mrna = float(quant_df["mrna"].sum())
    total_nrna = float(quant_df["nrna"].sum())
    total_gdna = float(estimator.gdna_total)
    total_rna = total_mrna + total_nrna
    total_all = total_rna + total_gdna + stats.n_intergenic

    summary = {
        "hulkrna_version": "0.1.0",
        "command": "hulkrna quant",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input": {
            "bam_file": str(bam_path.resolve()),
            "index_dir": str(index_dir.resolve()),
        },
        "library": {
            "protocol": (
                "R1-sense" if sm.read1_sense else "R1-antisense"
            ),
            "strand_specificity": round(sm.strand_specificity, 6),
            "p_r1_sense": round(sm.p_r1_sense, 6),
            "read1_sense": bool(sm.read1_sense),
            "frag_length_mean": round(flm.global_model.mean, 1),
            "frag_length_median": round(flm.global_model.median, 1),
            "frag_length_std": round(flm.global_model.std, 1),
            "frag_length_mode": int(flm.global_model.mode),
        },
        "alignment": {
            "total_reads": stats.total,
            "mapped_reads": stats.total - stats.unmapped,
            "unique_reads": stats.unique,
            "multimapping_reads": stats.multimapping,
            "proper_pairs": stats.proper_pair,
            "duplicate_reads": stats.duplicate,
            "qc_fail_reads": stats.qc_fail,
        },
        "fragments": {
            "total": stats.n_fragments,
            "genic": stats.n_fragments - stats.n_intergenic - stats.n_chimeric,
            "intergenic": stats.n_intergenic,
            "chimeric": stats.n_chimeric,
        },
        "quantification": {
            "n_transcripts": index.num_transcripts,
            "n_genes": index.num_genes,
            "n_loci": len(estimator.locus_results),
            "n_unambig_assigned": stats.deterministic_unambig_units,
            "n_em_assigned": (
                stats.em_routed_unambig_units
                + stats.em_routed_ambig_same_strand_units
                + stats.em_routed_ambig_opp_strand_units
                + stats.em_routed_multimapper_units
            ),
            "mrna_total": round(total_mrna, 2),
            "nrna_total": round(total_nrna, 2),
            "gdna_total": round(total_gdna, 2),
            "intergenic_total": stats.n_intergenic,
            "mrna_fraction": round(total_mrna / total_all, 6) if total_all > 0 else 0.0,
            "nrna_fraction": round(total_nrna / total_all, 6) if total_all > 0 else 0.0,
            "gdna_fraction": round(total_gdna / total_all, 6) if total_all > 0 else 0.0,
        },
        "strand_models": sm.to_dict(),
        "frag_length_models": flm.to_dict(),
        "pipeline_stats": stats.to_dict(),
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info(f"[DONE] Summary written to {summary_path}")

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
# YAML config / CLI merge helper
# ---------------------------------------------------------------------------

def _resolve_quant_args(
    args: argparse.Namespace,
    defaults: dict,
) -> None:
    """Merge CLI arguments, YAML config file, and hardcoded defaults.

    Priority (highest → lowest):
        1. Explicit CLI flags
        2. Values from the ``--config`` YAML file
        3. Hardcoded *defaults*

    Detection of "explicitly set on CLI" relies on all overridable
    argparse defaults being ``None``.  A non-None ``getattr(args, key)``
    means the user provided it on the command line.

    The *args* namespace is modified **in place**.
    """
    import yaml

    yaml_config: dict = {}
    config_path = getattr(args, "config", None)
    if config_path:
        with open(config_path) as fh:
            raw = yaml.safe_load(fh) or {}
        # Normalise: allow hyphens, convert to underscores
        yaml_config = {k.replace("-", "_"): v for k, v in raw.items()}
        unknown = set(yaml_config) - set(defaults)
        if unknown:
            logging.warning(
                "Unknown config keys ignored: %s", sorted(unknown),
            )

    for key, default_val in defaults.items():
        cli_val = getattr(args, key, None)
        if cli_val is not None:
            continue  # CLI wins
        yaml_val = yaml_config.get(key)
        if yaml_val is not None:
            setattr(args, key, yaml_val)
        else:
            setattr(args, key, default_val)


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
        "quant",
        help="Quantify transcripts: model training (Pass 1) + abundance estimation (Pass 2)",
    )

    # Required I/O
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

    # Optional config file
    quant_parser.add_argument(
        "--config", dest="config", default=None,
        help="YAML configuration file.  Any key matching a CLI option "
             "(using underscores) may be set here; explicit CLI flags "
             "override the file.",
    )

    # Boolean flags — use BooleanOptionalAction so both --flag and
    # --no-flag forms exist, and default=None lets us detect "not set
    # on CLI" for proper YAML override.
    quant_parser.add_argument(
        "--include-multimap", "--no-include-multimap",
        dest="include_multimap",
        action=argparse.BooleanOptionalAction, default=None,
        help="Include multimapping reads (default: yes).  "
             "Detected via NH tag (STAR) or secondary BAM flag (minimap2).",
    )
    quant_parser.add_argument(
        "--keep-duplicates", "--no-keep-duplicates",
        dest="keep_duplicates",
        action=argparse.BooleanOptionalAction, default=None,
        help="Keep reads marked as PCR/optical duplicates (default: no).",
    )
    quant_parser.add_argument(
        "--no-tsv", dest="no_tsv", action="store_true", default=None,
        help="Skip writing human-readable TSV quant files.",
    )

    # Common options (all default=None for YAML override)
    quant_parser.add_argument(
        "--sj-strand-tag", dest="sj_strand_tag",
        nargs="+", default=None,
        help="BAM tag(s) for splice-junction strand. 'auto' (default) "
             "detects the tag from the BAM. Use 'XS' for STAR, 'ts' for "
             "minimap2, or list multiple to check in order (e.g. XS ts).",
    )
    quant_parser.add_argument(
        "--seed", dest="seed", type=int, default=None,
        help="Random seed for reproducibility (default: use current timestamp).",
    )
    quant_parser.add_argument(
        "--em-prior-alpha", dest="em_prior_alpha", type=float, default=None,
        help="Flat Dirichlet pseudocount per eligible EM component "
             "(default: 0.01).",
    )
    quant_parser.add_argument(
        "--em-prior-gamma", dest="em_prior_gamma", type=float, default=None,
        help="OVR prior scale factor (gamma).  Default 1.0.  "
             "Controls how strongly coverage-weighted One Virtual Read "
             "priors influence the EM.",
    )
    quant_parser.add_argument(
        "--em-iterations", dest="em_iterations", type=int, default=None,
        help="Maximum EM iterations for ambiguous fragment resolution "
             "(default: 1000). Set to 0 for unambig-only quantification.",
    )
    quant_parser.add_argument(
        "--confidence-threshold", dest="confidence_threshold",
        type=float, default=None,
        help="Minimum RNA-normalized posterior for high-confidence "
             "assignment (default: 0.95).",
    )
    quant_parser.add_argument(
        "--annotated-bam", dest="annotated_bam", default=None,
        help="Write an annotated BAM with per-fragment assignment tags "
             "(ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS) to this path.  "
             "Requires a second BAM pass (modest runtime overhead).",
    )

    # -- Advanced parameters --------------------------------------------------
    adv = quant_parser.add_argument_group("advanced options")
    adv.add_argument(
        "--prune-threshold", dest="prune_threshold",
        type=float, default=None,
        help="Post-EM pruning evidence-ratio threshold (default: 0.1). "
             "Components with zero unambig evidence and evidence ratio "
             "(data_count / alpha) below this value are zeroed out and "
             "the EM is re-run to redistribute mass.  Set to -1 to "
             "disable pruning entirely.",
    )
    adv.add_argument(
        "--em-convergence-delta", dest="em_convergence_delta",
        type=float, default=None,
        help="Convergence threshold for EM parameter updates "
             "(default: 1e-6).",
    )
    adv.add_argument(
        "--gdna-splice-penalty-unannot",
        dest="gdna_splice_penalty_unannot",
        type=float, default=None,
        help="gDNA splice penalty for SPLICED_UNANNOT fragments "
             "(default: 0.01).",
    )
    adv.add_argument(
        "--overhang-alpha", dest="overhang_alpha",
        type=float, default=None,
        help="Per-base overhang penalty alpha in [0,1] (default: 0.01). "
             "0 = hard gate, 1 = no penalty.",
    )
    adv.add_argument(
        "--mismatch-alpha", dest="mismatch_alpha",
        type=float, default=None,
        help="Per-mismatch (NM tag) penalty alpha in [0,1] (default: 0.1). "
             "0 = hard gate, 1 = no penalty.",
    )
    adv.add_argument(
        "--tss-window", dest="tss_window",
        type=int, default=None,
        help="Fuzzy TSS grouping window in bp (default: 200). "
             "Transcripts whose 5' ends lie within this distance are "
             "grouped for the nrna_frac prior hierarchy.",
    )
    adv.add_argument(
        "--nrna-frac-kappa-global", dest="nrna_frac_kappa_global",
        type=float, default=None,
        help="Shrinkage κ pulling locus-strand nrna_frac toward the global "
             "prior.  Default: auto-estimate via Method of Moments.",
    )
    adv.add_argument(
        "--nrna-frac-kappa-locus", dest="nrna_frac_kappa_locus",
        type=float, default=None,
        help="Shrinkage κ pulling TSS-group nrna_frac toward the "
             "locus-strand estimate.  Default: auto-estimate.",
    )
    adv.add_argument(
        "--nrna-frac-kappa-tss", dest="nrna_frac_kappa_tss",
        type=float, default=None,
        help="Shrinkage κ pulling transcript nrna_frac toward the "
             "TSS-group estimate.  Default: auto-estimate.",
    )
    adv.add_argument(
        "--nrna-frac-mom-min-evidence-global",
        dest="nrna_frac_mom_min_evidence_global",
        type=float, default=None,
        help="Min fragment evidence for global MoM κ estimation "
             "(default: 50).",
    )
    adv.add_argument(
        "--nrna-frac-mom-min-evidence-locus",
        dest="nrna_frac_mom_min_evidence_locus",
        type=float, default=None,
        help="Min fragment evidence for locus MoM κ estimation "
             "(default: 30).",
    )
    adv.add_argument(
        "--nrna-frac-mom-min-evidence-tss",
        dest="nrna_frac_mom_min_evidence_tss",
        type=float, default=None,
        help="Min fragment evidence for TSS-level MoM κ estimation "
             "(default: 20).",
    )
    adv.add_argument(
        "--nrna-frac-kappa-min", dest="nrna_frac_kappa_min",
        type=float, default=None,
        help="Lower clamp for MoM-estimated κ (default: 2.0).",
    )
    adv.add_argument(
        "--nrna-frac-kappa-max", dest="nrna_frac_kappa_max",
        type=float, default=None,
        help="Upper clamp for MoM-estimated κ (default: 200.0).",
    )
    adv.add_argument(
        "--nrna-frac-kappa-fallback", dest="nrna_frac_kappa_fallback",
        type=float, default=None,
        help="Fallback κ when too few features for MoM (default: 5.0).",
    )
    adv.add_argument(
        "--nrna-frac-kappa-min-obs", dest="nrna_frac_kappa_min_obs",
        type=int, default=None,
        help="Minimum features required for MoM κ estimation; "
             "fewer triggers fallback (default: 20).",
    )
    adv.add_argument(
        "--gdna-kappa-chrom", dest="gdna_kappa_chrom",
        type=float, default=None,
        help="Shrinkage κ pulling chromosome gDNA rate toward the "
             "global estimate.  Default: auto-estimate via MoM.",
    )
    adv.add_argument(
        "--gdna-kappa-locus", dest="gdna_kappa_locus",
        type=float, default=None,
        help="Shrinkage κ pulling locus gDNA rate toward the "
             "chromosome estimate.  Default: auto-estimate via MoM.",
    )
    adv.add_argument(
        "--gdna-mom-min-evidence-chrom",
        dest="gdna_mom_min_evidence_chrom",
        type=float, default=None,
        help="Min fragment evidence for chromosome gDNA MoM κ "
             "estimation (default: 50).",
    )
    adv.add_argument(
        "--gdna-mom-min-evidence-locus",
        dest="gdna_mom_min_evidence_locus",
        type=float, default=None,
        help="Min fragment evidence for locus gDNA MoM κ "
             "estimation (default: 30).",
    )
    adv.add_argument(
        "--strand-prior-kappa", dest="strand_prior_kappa",
        type=float, default=None,
        help="Strand model prior pseudocount κ. Beta prior is "
             "Beta(κ/2, κ/2), shrinking toward 0.5 (max entropy). "
             "Default: 2.0 (uniform Beta(1,1) prior).",
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