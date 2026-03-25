"""
rigel.cli — Unified command-line interface.

Entry point: ``rigel`` (registered in pyproject.toml).

Subcommands:
    rigel index   — Build reference index from FASTA + GTF
    rigel quant   — Single-pass Bayesian fragment abundance estimation
    rigel sim     — Generate synthetic test scenarios
    rigel export  — Convert feather outputs to TSV or Parquet
"""

import argparse
import logging
import sys
import time
from pathlib import Path


def get_version() -> str:
    try:
        from . import __version__
        return f"rigel {__version__}"
    except ImportError:
        return "rigel (unknown version)"


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------

def index_command(args: argparse.Namespace) -> int:
    """Run the ``rigel index`` subcommand."""
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
        nrna_tolerance=args.nrna_tolerance,
    )
    return 0


def quant_command(args: argparse.Namespace) -> int:
    """Run the ``rigel quant`` subcommand.

    Single-pass pipeline: scan BAM, resolve fragments, train
    strand/insert models, then quantify via unified EM.
    """
    from .index import TranscriptIndex
    from .pipeline import run_pipeline

    # -- Resolve parameters: CLI > YAML > hardcoded defaults --
    _resolve_quant_args(args, _build_quant_defaults())

    # -- Validate required I/O args (may come from CLI or YAML config) --
    for attr, flag in [
        ("bam_file", "--bam"),
        ("index_dir", "--index"),
        ("output_dir", "-o/--output-dir"),
    ]:
        if getattr(args, attr, None) is None:
            sys.exit(f"Error: {flag} is required (via CLI or --config YAML)")

    bam_path = Path(args.bam_file)
    index_dir = Path(args.index_dir)
    output_dir = Path(args.output_dir)

    if not bam_path.exists():
        sys.exit(f"Error: BAM file not found: {bam_path}")
    if not index_dir.is_dir():
        sys.exit(f"Error: Index directory not found: {index_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # -- Resolve seed --
    seed = _resolve_seed(args)

    # -- Normalise sj_strand_tag (YAML may deliver a bare string) --
    if isinstance(args.sj_strand_tag, str):
        args.sj_strand_tag = [args.sj_strand_tag]
    sj_tag_list = args.sj_strand_tag
    sj_strand_tag = sj_tag_list[0] if len(sj_tag_list) == 1 else tuple(sj_tag_list)

    # -- Persist run config (now part of summary.json, written at end) --

    # -- Load reference index --
    logging.info(f"[START] Loading index from {index_dir}")
    index = TranscriptIndex.load(index_dir)
    logging.info(f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
                 f"{index.num_genes} genes")

    # -- Build pipeline config + run --
    pipeline_config = _build_pipeline_config(args, seed, sj_strand_tag)
    result = run_pipeline(bam_path, index, config=pipeline_config)

    # -- Write outputs --
    _write_quant_outputs(result, index, output_dir, args)

    return 0


# ---------------------------------------------------------------------------
# quant_command helpers
# ---------------------------------------------------------------------------

def _resolve_seed(args: argparse.Namespace) -> int:
    """Return an explicit seed or generate one from the current timestamp."""
    if args.seed is None:
        seed = int(time.time())
        logging.info(f"No seed provided, using timestamp: {seed}")
    else:
        seed = args.seed
        logging.info(f"Using provided seed: {seed}")
    return seed


def _build_pipeline_config(
    args: argparse.Namespace,
    seed: int,
    sj_strand_tag: str | tuple[str, ...],
) -> "PipelineConfig":
    """Translate resolved CLI args into a ``PipelineConfig``.

    Field mapping is driven by ``_PARAM_SPECS`` — see the declarative
    registry below ``sim_command``.
    """
    from .config import EMConfig, PipelineConfig, BamScanConfig, FragmentScoringConfig

    em_kw: dict = {}
    scan_kw: dict = {}
    scoring_kw: dict = {}
    top_kw: dict = {}
    _section = {"em": em_kw, "scan": scan_kw, "scoring": scoring_kw}

    for spec in _PARAM_SPECS:
        cli_val = getattr(args, spec.cli_dest)
        config_val = _cli_to_config(cli_val, spec.transform)
        if "." in spec.config_path:
            section, field_name = spec.config_path.split(".", 1)
            _section[section][field_name] = config_val
        else:
            top_kw[spec.config_path] = config_val

    # Override with pre-resolved values (seed and sj_strand_tag are
    # normalised in quant_command before reaching here).
    em_kw["seed"] = seed
    scan_kw["sj_strand_tag"] = sj_strand_tag

    logging.info(
        f"Using EM prior: pseudocount={em_kw.get('prior_pseudocount', 1.0)}"
    )

    return PipelineConfig(
        em=EMConfig(**em_kw),
        scan=BamScanConfig(**scan_kw),
        scoring=FragmentScoringConfig(**scoring_kw),
        **top_kw,
    )


def _write_quant_outputs(result, index, output_dir: Path, args) -> None:
    """Write quantification tables and summary.json."""
    import json

    quant_path = output_dir / "quant.feather"
    gene_quant_path = output_dir / "gene_quant.feather"
    nrna_quant_path = output_dir / "nrna_quant.feather"
    loci_path = output_dir / "loci.feather"
    summary_path = output_dir / "summary.json"
    config_yaml_path = output_dir / "config.yaml"

    # Log stats
    for key, val in sorted(result.stats.to_dict().items()):
        if isinstance(val, (int, float)):
            logging.info(f"  {key}: {val:,}")

    # Write quant tables (Feather v2 + ZSTD compression)
    estimator = result.estimator
    quant_df = estimator.get_counts_df(index)
    gene_quant_df = estimator.get_gene_counts_df(index)
    nrna_quant_df = estimator.get_nrna_counts_df(index)
    loci_df = estimator.get_loci_df(index)

    feather_kw = {"compression": "zstd"}
    quant_df.to_feather(str(quant_path), **feather_kw)
    gene_quant_df.to_feather(str(gene_quant_path), **feather_kw)
    nrna_quant_df.to_feather(str(nrna_quant_path), **feather_kw)
    loci_df.to_feather(str(loci_path), **feather_kw)
    logging.info(f"[DONE] Wrote {quant_path}, {gene_quant_path}, "
                 f"{nrna_quant_path}, {loci_path}")

    # Write TSV mirrors if requested
    if getattr(args, "tsv", False):
        for df, path in [
            (quant_df, quant_path),
            (gene_quant_df, gene_quant_path),
            (nrna_quant_df, nrna_quant_path),
            (loci_df, loci_path),
        ]:
            df.to_csv(str(path.with_suffix(".tsv")), sep="\t", index=False)

    # Build and write summary.json
    stats = result.stats
    sm = result.strand_models
    flm = result.frag_length_models
    total_mrna = float(quant_df["mrna"].sum())
    total_nrna = float(quant_df["nrna"].sum())
    total_gdna = float(estimator.gdna_em_count)
    total_rna = total_mrna + total_nrna
    total_all = total_rna + total_gdna + stats.n_intergenic

    from . import __version__

    # Strand model summary (exonic_spliced only)
    sm_primary = sm.exonic_spliced
    ci_lo, ci_hi = sm_primary.posterior_95ci()

    # Fragment length: full histograms + summary statistics per category
    fl_dict = flm.to_dict()

    # Calibration section
    cal_dict = None
    if result.calibration is not None:
        cal_dict = result.calibration.to_summary_dict()

    # Command section — record CLI arguments
    cmd_params: dict = {
        "bam_file": str(Path(args.bam_file).resolve()),
        "index_dir": str(Path(args.index_dir).resolve()),
        "output_dir": str(output_dir.resolve()),
    }
    seen = set()
    for spec in _PARAM_SPECS:
        dest = spec.cli_dest
        if dest in seen:
            continue
        seen.add(dest)
        val = getattr(args, dest, None)
        if dest == "sj_strand_tag" and not isinstance(val, str):
            val = list(val) if val else val
        cmd_params[dest] = val

    summary = {
        "rigel_version": __version__,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "command": {
            "subcommand": "quant",
            "arguments": cmd_params,
            "config_file": (
                str(Path(args.config).resolve()) if getattr(args, "config", None) else None
            ),
        },
        "configuration": result.pipeline_config.to_dict(),
        "input": {
            "bam_file": str(Path(args.bam_file).resolve()),
            "index_dir": str(Path(args.index_dir).resolve()),
        },
        "alignment_stats": {
            "total_reads": stats.total,
            "mapped_reads": stats.total - stats.unmapped,
            "unique_reads": stats.unique,
            "multimapping_reads": stats.multimapping,
            "proper_pairs": stats.proper_pair,
            "duplicate_reads": stats.duplicate,
            "qc_fail_reads": stats.qc_fail,
        },
        "fragment_stats": {
            "total": stats.n_fragments,
            "genic": stats.n_fragments - stats.n_intergenic - stats.n_chimeric,
            "intergenic": stats.n_intergenic,
            "chimeric": stats.n_chimeric,
            "chimeric_trans": stats.n_chimeric_trans,
            "chimeric_cis_same": stats.n_chimeric_cis_strand_same,
            "chimeric_cis_diff": stats.n_chimeric_cis_strand_diff,
            "with_annotated_sj": stats.n_with_annotated_sj,
            "with_unannotated_sj": stats.n_with_unannotated_sj,
        },
        "strand_model": {
            "protocol": "R1-sense" if sm.read1_sense else "R1-antisense",
            "strand_specificity": round(sm.strand_specificity, 6),
            "p_r1_sense": round(sm.p_r1_sense, 6),
            "read1_sense": bool(sm.read1_sense),
            "n_training_fragments": sm.n_observations,
            "posterior_variance": round(sm_primary.posterior_variance(), 8),
            "ci_95": [round(ci_lo, 6), round(ci_hi, 6)],
        },
        "calibration": cal_dict,
        "fragment_length": fl_dict,
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
    }
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    logging.info(f"[DONE] Summary written to {summary_path}")

    # Write config.yaml — reproducible run configuration
    _write_config_yaml(config_yaml_path, args)


def _write_config_yaml(config_yaml_path: Path, args: argparse.Namespace) -> None:
    """Write a config.yaml that can reproduce the current run.

    The file uses the same key format accepted by ``--config``, so a user
    can rerun with::

        rigel quant --config results/config.yaml
    """
    import yaml

    from . import __version__

    data: dict = {}

    # I/O paths (resolved to absolute)
    for key in ("bam_file", "index_dir", "output_dir"):
        val = getattr(args, key, None)
        if val is not None:
            data[key] = str(Path(val).resolve())

    # All tunable parameters from the declarative registry
    seen: set[str] = set()
    for spec in _PARAM_SPECS:
        dest = spec.cli_dest
        if dest in seen:
            continue
        seen.add(dest)
        val = getattr(args, dest, None)
        # Normalise types for clean YAML serialisation
        if isinstance(val, Path):
            val = str(val)
        elif isinstance(val, tuple):
            val = list(val)
        data[dest] = val

    # Extra flags not in _PARAM_SPECS
    data["tsv"] = bool(getattr(args, "tsv", False))

    header = (
        f"# Rigel run configuration\n"
        f"# Generated by {__version__} on {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
        f"#\n"
        f"# Rerun with:\n"
        f"#   rigel quant --config {config_yaml_path.name}\n"
        f"# CLI arguments override values in this file.\n"
    )
    with open(config_yaml_path, "w") as f:
        f.write(header)
        yaml.dump(
            data, f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )
    logging.info(f"[DONE] Config written to {config_yaml_path}")


def sim_command(args: argparse.Namespace) -> int:
    """Run the ``rigel sim`` subcommand."""
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


def export_command(args: argparse.Namespace) -> int:
    """Run the ``rigel export`` subcommand — convert feather outputs to TSV or Parquet."""
    import pandas as pd

    output_dir = Path(args.output_dir)
    fmt = args.format
    feather_files = sorted(output_dir.glob("*.feather"))
    if not feather_files:
        logging.error(f"No .feather files found in {output_dir}")
        return 1

    for fpath in feather_files:
        df = pd.read_feather(str(fpath))
        if fmt == "tsv":
            out = fpath.with_suffix(".tsv")
            df.to_csv(str(out), sep="\t", index=False)
        elif fmt == "parquet":
            out = fpath.with_suffix(".parquet")
            df.to_parquet(str(out), compression="zstd", index=False)
        else:
            logging.error(f"Unknown format: {fmt}")
            return 1
        logging.info(f"  {fpath.name} -> {out.name}")

    logging.info(f"[DONE] Exported {len(feather_files)} file(s) to {fmt}")
    return 0


# ---------------------------------------------------------------------------
# Declarative CLI ↔ config registry
# ---------------------------------------------------------------------------
#
# _PARAM_SPECS is the single mapping from argparse dest names to config
# dataclass fields.  Adding a new tuneable parameter requires:
#   1. Add the dataclass field in config.py
#   2. Add the argparse argument below in build_parser()
#   3. Add a _ParamSpec entry here
# Everything else (_build_quant_defaults, _build_pipeline_config,
# _write_run_config, _resolve_quant_args) is driven by this registry.

from dataclasses import dataclass as _dataclass
import math as _math


@_dataclass(frozen=True)
class _ParamSpec:
    """Maps a CLI argparse dest to a config dataclass field."""
    cli_dest: str               # argparse dest, e.g. "prior_pseudocount"
    config_path: str            # dotted config path, e.g. "em.prior_pseudocount"
    transform: str = "direct"   # "direct" | "invert_bool" | "log_penalty"
    #                             | "gdna_splice" | "path_or_none" | "sj_tag"


_PARAM_SPECS: tuple[_ParamSpec, ...] = (
    # -- EMConfig: direct mappings --
    _ParamSpec("seed", "em.seed"),
    _ParamSpec("prior_pseudocount", "em.prior_pseudocount"),
    _ParamSpec("em_iterations", "em.iterations"),
    _ParamSpec("em_convergence_delta", "em.convergence_delta"),
    _ParamSpec("assignment_mode", "em.assignment_mode"),
    _ParamSpec("assignment_min_posterior", "em.assignment_min_posterior"),
    _ParamSpec("em_mode", "em.mode"),
    # -- BamScanConfig: direct --
    _ParamSpec("include_multimap", "scan.include_multimap"),
    # -- BamScanConfig: transformed --
    _ParamSpec("keep_duplicates", "scan.skip_duplicates", "invert_bool"),
    _ParamSpec("sj_strand_tag", "scan.sj_strand_tag", "sj_tag"),
    _ParamSpec("tmpdir", "scan.spill_dir", "path_or_none"),
    # -- FragmentScoringConfig: transformed --
    _ParamSpec("overhang_alpha", "scoring.overhang_log_penalty", "log_penalty"),
    _ParamSpec("mismatch_alpha", "scoring.mismatch_log_penalty", "log_penalty"),
    _ParamSpec("gdna_splice_penalty_unannot", "scoring.gdna_splice_penalties", "gdna_splice"),
    _ParamSpec("pruning_min_posterior", "scoring.pruning_min_posterior"),
    # -- Fan-out: threads → both EM and scan --
    _ParamSpec("threads", "em.n_threads"),
    _ParamSpec("threads", "scan.n_scan_threads"),
    # -- Top-level PipelineConfig --
    _ParamSpec("annotated_bam", "annotated_bam_path"),
)


def _resolve_config_path(cfg: object, path: str) -> object:
    """Follow a dotted path like ``'em.prior_pseudocount'`` on a config object."""
    obj = cfg
    for attr in path.split("."):
        obj = getattr(obj, attr)
    return obj


def _config_to_cli(val: object, transform: str) -> object:
    """Transform a config-space default to CLI-space."""
    if transform == "direct":
        return val
    if transform == "invert_bool":
        return not val
    if transform == "sj_tag":
        return [val] if isinstance(val, str) else list(val)
    if transform == "log_penalty":
        return _math.exp(val) if val > float("-inf") else 0.0
    if transform == "gdna_splice":
        from .scoring import DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
        if val is None:
            return DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
        from .splice import SPLICE_UNANNOT
        return val.get(SPLICE_UNANNOT, DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT)
    if transform == "path_or_none":
        return str(val) if val else None
    raise ValueError(f"Unknown transform: {transform!r}")


def _cli_to_config(val: object, transform: str) -> object:
    """Transform a CLI-space value to config-space."""
    if transform == "direct":
        return val
    if transform == "invert_bool":
        return not val
    if transform == "sj_tag":
        # Normalised by quant_command before _build_pipeline_config;
        # the overwritten value is used, but keep this for completeness.
        if isinstance(val, str):
            return val
        return val[0] if len(val) == 1 else tuple(val)
    if transform == "log_penalty":
        from .scoring import overhang_alpha_to_log_penalty
        return overhang_alpha_to_log_penalty(val)
    if transform == "gdna_splice":
        from .scoring import GDNA_SPLICE_PENALTIES
        from .splice import SPLICE_UNANNOT
        penalties = dict(GDNA_SPLICE_PENALTIES)
        penalties[SPLICE_UNANNOT] = val
        return penalties
    if transform == "path_or_none":
        return Path(val) if val else None
    raise ValueError(f"Unknown transform: {transform!r}")


# ---------------------------------------------------------------------------
# YAML config / CLI merge helpers
# ---------------------------------------------------------------------------


def _build_quant_defaults() -> dict:
    """Build quant subcommand defaults from config dataclasses.

    Single source of truth: all default values derived from the frozen
    dataclasses in ``config.py`` via the ``_PARAM_SPECS`` registry.

    Scoring fields (overhang_alpha, mismatch_alpha, gdna_splice_penalty)
    use exact CLI-space constants from ``scoring.py`` to avoid float
    noise from log/exp round-tripping.
    """
    from .config import PipelineConfig
    from .scoring import (
        DEFAULT_OVERHANG_ALPHA,
        DEFAULT_MISMATCH_ALPHA,
        DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    )
    cfg = PipelineConfig()

    # Scoring fields have different representations in CLI (alpha) vs
    # config (log-penalty / dict).  Use exact constants to avoid float
    # noise from math.exp(math.log(alpha)).
    _scoring_defaults = {
        "overhang_alpha": DEFAULT_OVERHANG_ALPHA,
        "mismatch_alpha": DEFAULT_MISMATCH_ALPHA,
        "gdna_splice_penalty_unannot": DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    }

    defaults: dict = {}
    for spec in _PARAM_SPECS:
        if spec.cli_dest in defaults:
            continue  # skip fan-out duplicates (e.g. threads)
        if spec.cli_dest in _scoring_defaults:
            defaults[spec.cli_dest] = _scoring_defaults[spec.cli_dest]
        else:
            val = _resolve_config_path(cfg, spec.config_path)
            defaults[spec.cli_dest] = _config_to_cli(val, spec.transform)
    return defaults


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
        # I/O keys and extra flags are valid in config YAML
        valid_keys = set(defaults) | {"bam_file", "index_dir", "output_dir", "tsv"}
        unknown = set(yaml_config) - valid_keys
        if unknown:
            logging.warning(
                "Unknown config keys ignored: %s", sorted(unknown),
            )

    # Populate I/O args and extra flags from YAML when not set on CLI
    for io_key in ("bam_file", "index_dir", "output_dir", "tsv"):
        cli_val = getattr(args, io_key, None)
        if cli_val is not None:
            continue
        yaml_val = yaml_config.get(io_key)
        if yaml_val is not None:
            setattr(args, io_key, yaml_val)

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
        prog="rigel",
        description="rigel: Bayesian RNA-seq fragment abundance estimation",
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
    idx.add_argument(
        "--nrna-tolerance",
        dest="nrna_tolerance",
        type=int,
        default=20,
        help=(
            "Max distance (bp) for clustering transcript start/end sites "
            "when building synthetic nascent RNA transcripts (default: 20)"
        ),
    )
    idx.set_defaults(func=index_command)

    # --- QUANT ---------------------------------------------------------------
    quant_parser = subparsers.add_parser(
        "quant",
        help="Single-pass BAM scan and Bayesian transcript abundance estimation",
    )

    # I/O (required via CLI or config YAML)
    quant_parser.add_argument(
        "--bam", dest="bam_file", default=None,
        help="Name-sorted BAM file. Minimap2 and STAR aligners are "
             "supported. The NH tag is used when present; otherwise "
             "multimappers are detected from secondary BAM flags.",
    )
    quant_parser.add_argument(
        "--index", dest="index_dir", default=None,
        help="Directory containing rigel index files (from 'rigel index').",
    )
    quant_parser.add_argument(
        "-o", "--output-dir", dest="output_dir", default=None,
        help="Output directory for quantification results.",
    )

    # Optional config file
    quant_parser.add_argument(
        "--config", dest="config", default=None,
        help="YAML configuration file. Accepts the same keys as CLI "
             "options (with underscores). CLI flags override config file "
             "values. A config.yaml is written to the output directory "
             "after each run for easy reproducibility.",
    )

    # Boolean flags — use BooleanOptionalAction so both --flag and
    # --no-flag forms exist, and default=None lets us detect "not set
    # on CLI" for proper YAML override.
    quant_parser.add_argument(
        "--include-multimap",
        dest="include_multimap",
        action=argparse.BooleanOptionalAction, default=None,
        help="Include multimapping reads (default: yes). "
             "Detected via NH tag or secondary BAM flag.",
    )
    quant_parser.add_argument(
        "--keep-duplicates",
        dest="keep_duplicates",
        action=argparse.BooleanOptionalAction, default=None,
        help="Keep reads marked as PCR/optical duplicates (default: no).",
    )
    quant_parser.add_argument(
        "--tsv", dest="tsv", action="store_true", default=False,
        help="Also write TSV (.tsv) mirrors of quant tables.",
    )

    # Common options (all default=None for YAML override)
    quant_parser.add_argument(
        "--sj-strand-tag", dest="sj_strand_tag",
        nargs="+", default=None,
        help="BAM tag(s) for splice-junction strand (default: auto). "
             "'auto' detects the tag from the first 10,000 reads. "
             "Use 'XS' for STAR, 'ts' for minimap2, or list multiple "
             "tags to check in order (e.g. XS ts).",
    )
    quant_parser.add_argument(
        "--seed", dest="seed", type=int, default=None,
        help="Random seed for reproducibility (default: use current timestamp).",
    )
    quant_parser.add_argument(
        "--prior-pseudocount", dest="prior_pseudocount",
        type=float, default=None,
        help="Total OVR prior budget C in virtual fragments (default: 1.0). "
             "Distributed as gamma*C to gDNA and (1-gamma)*C coverage-weighted "
             "across RNA components.",
    )
    quant_parser.add_argument(
        "--em-iterations", dest="em_iterations", type=int, default=None,
        help="Maximum EM iterations (default: 1000). "
             "Set to 0 for unambiguous-only quantification.",
    )
    quant_parser.add_argument(
        "--assignment-mode", dest="assignment_mode",
        choices=["sample", "fractional", "map"], default=None,
        help="Post-EM fragment assignment mode (default: sample). "
             "'sample' draws from the posterior distribution. "
             "'fractional' preserves EM posterior weights (traditional). "
             "'map' assigns each fragment to its highest-posterior component.",
    )
    quant_parser.add_argument(
        "--annotated-bam", dest="annotated_bam", default=None,
        help="Write an annotated BAM with per-fragment assignment tags "
             "(ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS) to this path. "
             "Requires a second pass over the BAM.",
    )
    quant_parser.add_argument(
        "--tmpdir", default=None,
        help="Directory for temporary buffer spill files when memory "
             "limits are exceeded (default: system temp directory).",
    )

    # -- Advanced parameters --------------------------------------------------
    adv = quant_parser.add_argument_group("advanced options")
    adv.add_argument(
        "--assignment-min-posterior", dest="assignment_min_posterior",
        type=float, default=None,
        help="Minimum posterior for a component to be eligible for "
             "discrete assignment (map/sample modes). Default: 0.01.",
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
        "--pruning-min-posterior", dest="pruning_min_posterior",
        type=float, default=None,
        help="Minimum posterior threshold for candidate pruning "
             "(default: 1e-4). Lower values keep more candidates "
             "(conservative). Set to 0 to disable pruning entirely.",
    )
    adv.add_argument(
        "--em-mode", dest="em_mode",
        choices=["vbem", "map"], default=None,
        help="EM algorithm variant (default: vbem). "
             "'vbem' uses Variational Bayes EM with digamma-based soft "
             "updates. 'map' uses MAP-EM with hard max(0, n+a-1) updates.",
    )
    quant_parser.add_argument(
        "--threads", dest="threads", type=int, default=None,
        help="Number of threads for BAM scanning and locus EM "
             "(default: 0 = all available cores). Set to 1 for "
             "sequential execution.",
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

    # --- EXPORT --------------------------------------------------------------
    export = subparsers.add_parser(
        "export", help="Convert feather output files to TSV or Parquet",
    )
    export.add_argument(
        "output_dir",
        help="Directory containing .feather output files from rigel quant",
    )
    export.add_argument(
        "-f", "--format", dest="format", default="tsv",
        choices=["tsv", "parquet"],
        help="Output format (default: tsv)",
    )
    export.set_defaults(func=export_command)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    """CLI entry point registered as ``rigel`` in pyproject.toml."""
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