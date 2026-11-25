"""The ``rigel`` command line (the entry point registered in pyproject.toml).

Subcommands:
    rigel index   — Build reference index from FASTA + GTF
    rigel quant   — Single-pass BAM scan and Bayesian transcript abundance estimation
    rigel sim     — Generate synthetic test scenarios
    rigel export  — Convert feather outputs to TSV or Parquet
    rigel report  — Build a self-contained HTML QC report from a quant output directory

``rigel quant``'s tunable options map onto :class:`rigel.config.PipelineConfig` through one table, ``PARAM_SPECS``:
the parser's defaults, the keys a ``--config`` YAML may set, the run's configuration and the ``config.yaml`` a run
writes are all read from it.
"""

from __future__ import annotations

import argparse
import logging
import math
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, get_args

from . import outputs

if TYPE_CHECKING:
    from .config import PipelineConfig


def get_version() -> str:
    from . import __version__

    return f"rigel {__version__}"


# ---------------------------------------------------------------------------
# Subcommand handlers
# ---------------------------------------------------------------------------


def index_command(args: argparse.Namespace) -> int:
    """Run the ``rigel index`` subcommand: the arguments resolved (CLI flags over the ``--config`` YAML over
    :class:`~rigel.config.IndexConfig`'s defaults), checked, then :meth:`TranscriptIndex.build`."""
    from .index import TranscriptIndex

    try:
        resolve_index_args(args)
    except (ValueError, OSError) as exc:
        sys.exit(f"Error: {exc}")
    for key in ("fasta", "gtf", "output_dir"):
        if getattr(args, key) is None:
            sys.exit(f"Error: --{key.replace('_', '-')} is required (via CLI or --config YAML)")
    fasta = Path(args.fasta)
    gtf = Path(args.gtf)

    if not fasta.exists():
        sys.exit(f"Error: FASTA file not found: {fasta}")
    if not gtf.exists():
        sys.exit(f"Error: GTF file not found: {gtf}")

    if args.no_mappability and args.alignable_zarr is not None:
        sys.exit("Error: --no-mappability and --alignable-zarr are mutually exclusive.")
    if not args.no_mappability and args.alignable_zarr is None:
        sys.exit(
            "Error: --alignable-zarr PATH is required (use --no-mappability to opt out "
            "explicitly, e.g. on a synthetic genome)."
        )

    try:
        TranscriptIndex.build(
            fasta_file=fasta,
            gtf_file=gtf,
            output_dir=args.output_dir,
            feather_compression=args.feather_compression,
            write_tsv=not args.no_tsv,
            gtf_parse_mode=args.gtf_parse_mode,
            collapse_duplicate_transcripts=args.collapse_duplicate_transcripts,
            end_tolerance=args.end_tolerance,
            alignable_zarr_path=(Path(args.alignable_zarr) if args.alignable_zarr else None),
        )
    # a malformed annotation, an unreadable store, references that disagree
    except (ValueError, OSError) as exc:
        sys.exit(f"Error: {exc}")
    return 0


def quant_command(args: argparse.Namespace) -> int:
    """Run the ``rigel quant`` subcommand.

    Resolves the arguments (CLI flags over the ``--config`` YAML over the config defaults), loads the index, runs
    :func:`rigel.pipeline.run_pipeline` — one BAM scan, calibration, then the per-locus EM — and writes the output
    tables, ``summary.json`` and ``config.yaml``; with ``--annotated-bam``, a second pass over the BAM then writes the
    annotated BAM, checked against the scan. Every setting is checked before the index loads, and every refusal is an
    error, not a traceback.
    """
    from .index import IndexLoadError, TranscriptIndex
    from .pipeline import BamInputError, run_pipeline

    try:
        resolve_quant_args(args)
        bam_path, index_dir, output_dir = _quant_paths(args)
        output_dir.mkdir(parents=True, exist_ok=True)
        # YAML may deliver a bare string; the CLI always delivers a list
        if isinstance(args.sj_strand_tag, str):
            args.sj_strand_tag = [args.sj_strand_tag]
        pipeline_config = build_pipeline_config(args)
    # a TypeError is a setting of the wrong type, as YAML's string `1e3` for an integer; an OSError a --config this
    # run cannot read or an output directory it cannot make
    except (ValueError, TypeError, OSError) as exc:
        sys.exit(f"Error: {exc}")
    _check_annotated_bam_target(pipeline_config.annotated_bam_path, bam_path, index_dir, output_dir)

    logging.info(f"[START] Loading index from {index_dir}")
    try:
        index = TranscriptIndex.load(index_dir)
    except IndexLoadError as exc:
        sys.exit(f"Error: {exc}")
    logging.info(
        f"[DONE] Loaded index: {index.num_transcripts} transcripts, "
        f"{index.num_annotated_genes} genes"
    )

    # The tables are written before the optional annotated-BAM pass, so a failure there loses nothing
    try:
        run_pipeline(
            bam_path,
            index,
            config=pipeline_config,
            before_annotated_bam=lambda result: write_quant_outputs(
                result, index, output_dir, args
            ),
        )
    except BamInputError as exc:
        sys.exit(f"Error: {exc}")
    return 0


def _quant_paths(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    """The BAM, the index directory and the output directory: each required on the command line or in the
    ``--config`` YAML, and the first two must exist."""
    for attr, flag in (
        ("bam_file", "--bam"),
        ("index_dir", "--index"),
        ("output_dir", "-o/--output-dir"),
    ):
        if getattr(args, attr) is None:
            raise ValueError(f"{flag} is required (via CLI or --config YAML)")
    bam_path, index_dir = Path(args.bam_file), Path(args.index_dir)
    if not bam_path.exists():
        raise ValueError(f"BAM file not found: {bam_path}")
    if not index_dir.is_dir():
        raise ValueError(f"Index directory not found: {index_dir}")
    return bam_path, index_dir, Path(args.output_dir)


def _check_annotated_bam_target(
    annotated, bam_path: Path, index_dir: Path, output_dir: Path
) -> None:
    """Refuse, before the run, an annotated BAM the second pass could not write, or one that names a file this run
    reads or writes: the input BAM, a file of the index, or a product of the output directory (:mod:`rigel.outputs`),
    by path or through a link. That pass runs after the quantification, so the failure would otherwise come only at
    the end, and over any of those files it would destroy it."""
    from .pipeline import BamInputError, check_annotated_bam_path

    if annotated is None:
        return
    try:
        check_annotated_bam_path(bam_path, annotated)
    except BamInputError as exc:
        sys.exit(f"Error: {exc}")
    target = Path(annotated)
    for path in (*outputs.managed_paths(output_dir), *sorted(Path(index_dir).iterdir())):
        if target.resolve() == path.resolve() or (
            target.exists() and path.exists() and target.samefile(path)
        ):
            sys.exit(
                f"Error: --annotated-bam {annotated} is {path}, a file this run reads or writes: name another path"
            )
    target_dir = Path(annotated).resolve().parent
    if not target_dir.is_dir() or not os.access(target_dir, os.W_OK):
        sys.exit(
            f"Error: --annotated-bam {annotated}: {target_dir} is not a directory this run can write to"
        )


def build_pipeline_config(args: argparse.Namespace) -> PipelineConfig:
    """The resolved arguments as a ``PipelineConfig``, field by field through ``PARAM_SPECS``."""
    from .config import (
        BamScanConfig,
        CalibrationConfig,
        EMConfig,
        FragmentScoringConfig,
        PipelineConfig,
    )

    sections: dict[str, dict] = {"em": {}, "scan": {}, "scoring": {}, "calibration": {}}
    top: dict = {}
    for spec in PARAM_SPECS:
        value = _cli_to_config(getattr(args, spec.cli_dest), spec.transform)
        section, _, field_name = spec.config_path.rpartition(".")
        (sections[section] if section else top)[field_name] = value
    return PipelineConfig(
        em=EMConfig(**sections["em"]),
        scan=BamScanConfig(**sections["scan"]),
        scoring=FragmentScoringConfig(**sections["scoring"]),
        calibration=CalibrationConfig(**sections["calibration"]),
        **top,
    )


# ---------------------------------------------------------------------------
# rigel quant's outputs: the tables, summary.json and config.yaml
# ---------------------------------------------------------------------------


def write_quant_outputs(result, index, output_dir: Path, args: argparse.Namespace) -> None:
    """``rigel quant``'s outputs (:func:`rigel.quant_writer.write_quant_outputs`), with the command and ``config.yaml``
    as this run's arguments resolved them."""
    from . import quant_writer

    command = {
        "subcommand": "quant",
        "arguments": _command_arguments(args, output_dir),
        "config_file": str(Path(args.config).resolve()) if args.config else None,
    }
    quant_writer.write_quant_outputs(
        result,
        index,
        output_dir,
        tsv=args.tsv,
        command=command,
        write_config=lambda path: write_config_yaml(path, args),
    )


def _command_arguments(args: argparse.Namespace, output_dir: Path) -> dict:
    """The command's arguments: the three paths resolved, then every tunable option as the run resolved it."""
    arguments: dict = {
        "bam_file": str(Path(args.bam_file).resolve()),
        "index_dir": str(Path(args.index_dir).resolve()),
        "output_dir": str(output_dir.resolve()),
    }
    for dest in _OPTION_DESTS:
        val = getattr(args, dest)
        if dest == "sj_strand_tag" and not isinstance(val, str):
            val = list(val) if val else val
        arguments[dest] = val
    return arguments


def write_config_yaml(config_yaml_path: Path, args: argparse.Namespace) -> None:
    """Write a config.yaml that reproduces the run.

    The file uses the keys ``--config`` accepts, every path absolute, so a user can rerun it from any directory::

        rigel quant --config results/config.yaml

    The junction-strand tag is the one the run read with, as the list the CLI takes.
    """
    import yaml

    from . import __version__

    data: dict = {
        key: str(Path(getattr(args, key)).resolve())
        for key in ("bam_file", "index_dir", "output_dir")
        if getattr(args, key) is not None
    }
    quant: dict = {}
    for dest in _OPTION_DESTS:
        val = getattr(args, dest)
        if dest in ("annotated_bam", "tmpdir") and val is not None:
            val = str(Path(val).resolve())
        elif dest == "sj_strand_tag":
            tags = _cli_to_config(val, "sj_tag")
            val = list(tags) if isinstance(tags, tuple) else tags.split(",")
        quant[dest] = val
    data["quant"] = quant
    data["tsv"] = bool(args.tsv)

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
        yaml.dump(data, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    logging.info(f"[DONE] Config written to {config_yaml_path}")


# ---------------------------------------------------------------------------
# rigel sim, export and report
# ---------------------------------------------------------------------------

#: The read-model keys a ``rigel sim`` scenario passes to ``ReadSimConfig``.
_SIM_READ_KEYS = ("frag_mean", "frag_std", "frag_min", "frag_max", "read_length", "error_rate")
#: Every top-level key ``sim_command`` reads.
_SIM_KEYS = frozenset(
    {"name", "ref_name", "genome_length", "seed", "n_fragments", "genes", *_SIM_READ_KEYS}
)
#: Every key a gene entry reads, and the ones it must have.
_SIM_GENE_KEYS = frozenset({"gene_id", "strand", "transcripts", "gene_name", "gene_type"})
_SIM_GENE_REQUIRED = frozenset({"gene_id", "strand", "transcripts"})


def _read_sim_scenario(path) -> dict:
    """A ``rigel sim`` scenario, checked before anything is written: every key one this command reads (a misspelled
    one would run on its default), at least one gene, and every gene's and transcript's keys."""
    import yaml

    from .sim.annotation import check_transcript_keys

    with open(path) as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError(f"{path} must hold a mapping of scenario keys")
    unknown = sorted(set(cfg) - _SIM_KEYS)
    if unknown:
        raise ValueError(f"unknown scenario key(s): {unknown}; known: {sorted(_SIM_KEYS)}")
    genes = cfg.get("genes") or []
    if not genes:
        raise ValueError("a scenario needs at least one gene under `genes:`")
    for gene in genes:
        if not isinstance(gene, dict):
            raise ValueError(f"a gene entry must be a mapping, not {gene!r}")
        unknown = sorted(set(gene) - _SIM_GENE_KEYS)
        missing = sorted(_SIM_GENE_REQUIRED - set(gene))
        if unknown or missing:
            raise ValueError(
                f"gene {gene.get('gene_id')!r}: "
                + (f"unknown key(s) {unknown}" if unknown else f"missing {missing}")
                + f"; a gene takes {sorted(_SIM_GENE_KEYS)}"
            )
        for tdef in gene["transcripts"] or []:
            check_transcript_keys(gene["gene_id"], tdef)
    cfg["genes"] = genes
    return cfg


def sim_command(args: argparse.Namespace) -> int:
    """Run the ``rigel sim`` subcommand; a scenario's YAML value wins over the command-line flag of its name."""
    from .config import CONSTANTS
    from .sim import ReadSimConfig, Scenario

    try:
        cfg = _read_sim_scenario(args.config)
    except ValueError as exc:
        sys.exit(f"Error: {exc}")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    seed = cfg.get("seed", args.seed)
    sim_config = ReadSimConfig(**{k: cfg[k] for k in _SIM_READ_KEYS if k in cfg}, seed=seed)
    sc = Scenario(
        name=cfg.get("name", "scenario"),
        genome_length=cfg.get("genome_length", args.genome_length),
        seed=seed,
        work_dir=output_dir,
        ref_name=cfg.get("ref_name", CONSTANTS.simulator.ref_name),
    )
    for gene in cfg["genes"]:
        # a gene's optional keys are passed only when given, so the simulator's own defaults stand
        sc.add_gene(
            gene_id=gene["gene_id"],
            strand=gene["strand"],
            transcripts=gene["transcripts"],
            **{key: gene[key] for key in ("gene_name", "gene_type") if key in gene},
        )
    result = sc.build(n_fragments=cfg.get("n_fragments", args.num_reads), sim_config=sim_config)

    logging.info(f"Scenario artifacts written to {output_dir}")
    logging.info(f"  FASTA: {result.fasta_path}")
    logging.info(f"  GTF:   {result.gtf_path}")
    logging.info(f"  BAM:   {result.bam_path}")
    logging.info(f"  Index: {result.index_dir}")
    return 0


def export_command(args: argparse.Namespace) -> int:
    """Run the ``rigel export`` subcommand — convert a complete run's tables to TSV or Parquet. Only the run's own
    products are converted (:mod:`rigel.outputs`), each published whole; the next ``rigel quant`` in the directory
    removes them with its other products, so an export never outlives the run it was made from."""
    import pandas as pd

    output_dir = Path(args.output_dir)
    fmt = args.format
    if not outputs.is_complete(output_dir):
        logging.error(
            f"{output_dir} holds no complete run (no {outputs.SUMMARY}): nothing to export"
        )
        return 1
    feather_files = [
        p for p in (outputs.table_path(output_dir, n) for n in outputs.TABLES) if p.is_file()
    ]

    for fpath in feather_files:
        df = pd.read_feather(str(fpath))
        if fmt == "tsv":
            out = fpath.with_suffix(".tsv")
            outputs.publish(out, lambda p: df.to_csv(str(p), sep="\t", index=False))
        else:
            out = fpath.with_suffix(".parquet")
            outputs.publish(
                out,
                lambda p: df.to_parquet(str(p), compression=outputs.TABLE_COMPRESSION, index=False),
            )
        logging.info(f"  {fpath.name} -> {out.name}")

    logging.info(f"[DONE] Exported {len(feather_files)} file(s) to {fmt}")
    return 0


def report_command(args: argparse.Namespace) -> int:
    """Run the ``rigel report`` subcommand — build an HTML QC report from a run."""
    from .report import build_report
    from .report.substrate import SubstrateError

    try:
        out = build_report(args.output_dir, out_path=args.output, title=args.title)
    except SubstrateError as exc:
        logging.error(f"{exc}")
        return 1
    logging.info(f"[DONE] Report written to {out}")
    return 0


# ---------------------------------------------------------------------------
# The quant options ↔ PipelineConfig registry
# ---------------------------------------------------------------------------
#
# A new tunable needs its config field (config.py), its option (the _add_quant_*_options builders) and a ParamSpec
# row below; its default, the config the run reads, summary.json's arguments and config.yaml all follow from the row.


@dataclass(frozen=True)
class ParamSpec:
    """Maps a CLI argparse dest to a config dataclass field."""

    cli_dest: str  # argparse dest, e.g. "em_iterations"
    config_path: str  # dotted config path, e.g. "em.iterations"
    transform: str = (
        "direct"  # "direct" | "invert_bool" | "sj_tag" | "path_or_none" | "gib_to_bytes"
    )


PARAM_SPECS: tuple[ParamSpec, ...] = (
    # -- EMConfig: direct mappings --
    ParamSpec("seed", "em.seed"),
    ParamSpec("em_iterations", "em.iterations"),
    ParamSpec("em_convergence_delta", "em.convergence_delta"),
    ParamSpec("assignment_mode", "em.assignment_mode"),
    ParamSpec("em_mode", "em.mode"),
    # -- BamScanConfig: direct --
    ParamSpec("include_multimap", "scan.include_multimap"),
    ParamSpec("intron_overhang_tolerance", "scan.intron_overhang_tolerance"),
    ParamSpec("scan_bgzf_threads", "scan.bgzf_threads"),
    ParamSpec("scan_fragments_per_chunk", "scan.fragments_per_chunk"),
    ParamSpec("scan_read_name_batch_size", "scan.read_name_batch_size"),
    # -- BamScanConfig: transformed --
    ParamSpec("keep_duplicates", "scan.skip_duplicates", "invert_bool"),
    ParamSpec("sj_strand_tag", "scan.sj_strand_tag", "sj_tag"),
    ParamSpec("tmpdir", "scan.spill_dir", "path_or_none"),
    # -- FragmentScoringConfig --
    ParamSpec("overhang_alpha", "scoring.overhang_alpha"),
    ParamSpec("pruning_min_posterior", "scoring.pruning_min_posterior"),
    # -- CalibrationConfig: advanced --
    ParamSpec("calib_refit_iters", "calibration.calib_refit_iters"),
    ParamSpec("sweep_block_slots", "calibration.sweep_block_slots"),
    ParamSpec("intron_departure_prior", "calibration.intron_departure_prior"),
    # -- Fan-out: total threads → the EM's, the scan's and calibration's budgets --
    ParamSpec("threads", "em.n_threads"),
    ParamSpec("threads", "scan.total_threads"),
    ParamSpec("threads", "calibration.n_threads"),
    # -- BamScanConfig: buffer sizing --
    ParamSpec("scan_buffer_size", "scan.buffer_size_bytes", "gib_to_bytes"),
    # -- Top-level PipelineConfig --
    ParamSpec("annotated_bam", "annotated_bam_path"),
    ParamSpec("emit_locus_stats", "emit_locus_stats"),
)

#: Each tunable option once, in table order: ``threads`` fans out to three config fields but is one option.
_OPTION_DESTS: tuple[str, ...] = tuple(dict.fromkeys(spec.cli_dest for spec in PARAM_SPECS))

#: The keys a ``--config`` YAML may set besides the tunable options, each with the value it takes when neither the
#: command line nor the YAML sets it.
_CONFIG_FILE_KEYS = {"bam_file": None, "index_dir": None, "output_dir": None, "tsv": False}


def _config_to_cli(val: object, transform: str) -> object:
    """Transform a config-space default to CLI-space."""
    if transform == "direct":
        return val
    if transform == "invert_bool":
        return not val
    if transform == "sj_tag":
        return [val] if isinstance(val, str) else list(val)
    if transform == "path_or_none":
        return str(val) if val else None
    if transform == "gib_to_bytes":
        return val / (1024**3)  # bytes → GiB for CLI display
    raise ValueError(f"Unknown transform: {transform!r}")


def _cli_to_config(val: object, transform: str) -> object:
    """Transform a CLI-space value to config-space. A value a transform would read as another kind is refused, since a
    YAML value reaches here without the parser's types: ``not "false"`` is a valid ``False``, and a GiB count that is
    a string or a bool would multiply."""
    if transform == "direct":
        return val
    if transform == "invert_bool":
        if not isinstance(val, bool):
            raise TypeError(f"a switch must be true or false; got {val!r}")
        return not val
    if transform == "sj_tag":
        if isinstance(val, str):
            return val
        return val[0] if len(val) == 1 else tuple(val)
    if transform == "path_or_none":
        return Path(val) if val else None
    if transform == "gib_to_bytes":
        if isinstance(val, bool) or not isinstance(val, (int, float)) or not math.isfinite(val):
            raise TypeError(f"a size in GiB must be a finite number; got {val!r}")
        size = int(val * 1024**3)  # GiB → bytes
        if val < 0 or (val > 0 and size == 0):
            # 0 is a setting of its own (no spilling); a size that rounds to it would switch spilling off unasked
            raise ValueError(
                f"a size in GiB must be 0 (no spilling) or at least one byte; got {val!r}"
            )
        return size
    raise ValueError(f"Unknown transform: {transform!r}")


def build_quant_defaults() -> dict:
    """Every tunable option's default, read off ``PipelineConfig()`` through ``PARAM_SPECS``, so the config
    dataclasses are the one source of every default; a fanned-out option takes its first field's."""
    from .config import PipelineConfig

    cfg = PipelineConfig()
    defaults: dict = {}
    for spec in PARAM_SPECS:
        if spec.cli_dest in defaults:
            continue
        val = cfg
        for attr in spec.config_path.split("."):
            val = getattr(val, attr)
        defaults[spec.cli_dest] = _config_to_cli(val, spec.transform)
    return defaults


# ---------------------------------------------------------------------------
# The --config YAML: CLI > YAML > defaults
# ---------------------------------------------------------------------------


def resolve_quant_args(args: argparse.Namespace) -> None:
    """Fill every option in place: an explicit CLI flag wins, then the ``--config`` YAML, then
    :func:`build_quant_defaults`.

    Every overridable argparse default is ``None``, so a non-``None`` value on *args* is one the user gave on the
    command line.
    """
    defaults = build_quant_defaults()
    yaml_config = _read_config_file(args.config, defaults) if args.config else {}
    for key, default_val in {**_CONFIG_FILE_KEYS, **defaults}.items():
        if getattr(args, key) is not None:
            continue  # CLI wins
        yaml_val = yaml_config.get(key)
        setattr(args, key, default_val if yaml_val is None else yaml_val)


def build_index_defaults() -> dict:
    """Every ``rigel index`` option's default: the paths none, the switches off, the settings
    :class:`~rigel.config.IndexConfig`'s, so the config dataclass is the one source of every default."""
    from .config import IndexConfig

    cfg = IndexConfig()
    return {
        "fasta": None,
        "gtf": None,
        "output_dir": None,
        "alignable_zarr": None,
        "no_mappability": False,
        "no_tsv": not cfg.write_tsv,
        "feather_compression": cfg.feather_compression,
        "gtf_parse_mode": cfg.gtf_parse_mode,
        "collapse_duplicate_transcripts": cfg.collapse_duplicate_transcripts,
        "end_tolerance": cfg.end_tolerance,
    }


def resolve_index_args(args: argparse.Namespace) -> None:
    """Fill every ``rigel index`` option in place: an explicit CLI flag wins, then the ``--config`` YAML (its
    ``index:`` block or its top level), then :func:`build_index_defaults`."""
    defaults = build_index_defaults()
    yaml_config = _read_config_file(args.config, defaults, block="index") if args.config else {}
    for key, default_val in defaults.items():
        if getattr(args, key) is not None:
            continue  # CLI wins
        yaml_val = yaml_config.get(key)
        setattr(args, key, default_val if yaml_val is None else yaml_val)


def _read_config_file(config_path, defaults: dict, block: str = "quant") -> dict:
    """A ``--config`` YAML's keys and values for one command: hyphens read as underscores, its *block*'s keys
    (``quant:`` or ``index:``) under the top-level ones, the other command's block left alone, and a key that is
    neither an option nor, for ``quant``, one of ``_CONFIG_FILE_KEYS`` refused."""
    import yaml

    with open(config_path) as fh:
        text = fh.read()
    try:
        raw = yaml.safe_load(text) or {}
    except yaml.YAMLError as exc:
        raise ValueError(f"{config_path} is not valid YAML: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(
            f"{config_path} must hold a mapping of keys to values, not a {type(raw).__name__}"
        )
    raw_normalized = {str(k).replace("-", "_"): v for k, v in raw.items()}
    own = raw_normalized.pop(block, None)
    raw_normalized.pop("index" if block == "quant" else "quant", None)
    yaml_config: dict = {}
    if own is not None:
        if not isinstance(own, dict):
            raise ValueError(f"{block} config block must be a mapping.")
        yaml_config.update({str(k).replace("-", "_"): v for k, v in own.items()})
    yaml_config.update(raw_normalized)
    valid_keys = set(defaults) | (set(_CONFIG_FILE_KEYS) if block == "quant" else set())
    unknown = set(yaml_config) - valid_keys
    if unknown:
        raise ValueError(_unknown_config_keys(config_path, text, unknown, valid_keys))
    return yaml_config


def _unknown_config_keys(config_path, text: str, unknown: set, valid_keys: set) -> str:
    """The refusal of a ``--config`` file's unknown keys: those keys, the ones this version knows, and, for a
    ``config.yaml`` another rigel wrote (whose removed keys are the usual cause), the version that wrote it."""
    import re

    from . import __version__

    message = (
        f"unknown config key(s) in {config_path}: {sorted(unknown)}. "
        f"rigel {__version__} knows: {', '.join(sorted(valid_keys))}."
    )
    written = re.search(r"^# Generated by (\S+) on ", text, flags=re.MULTILINE)
    if written is not None and written.group(1) != __version__:
        message += (
            f" The file was written by rigel {written.group(1)}, and a config.yaml reruns an analysis exactly"
            " only on the version that wrote it: rerun with that version, or delete the keys this one does not know."
        )
    return message


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Construct the top-level argument parser with all subcommands."""
    parser = argparse.ArgumentParser(
        prog="rigel",
        description="rigel: Bayesian RNA-seq fragment abundance estimation",
    )
    parser.add_argument("--version", action="version", version=get_version())
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose (DEBUG-level) logging"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_index_parser(subparsers)
    _add_quant_parser(subparsers)
    _add_sim_parser(subparsers)
    _add_export_parser(subparsers)
    _add_report_parser(subparsers)
    return parser


def _add_index_parser(subparsers) -> None:
    """``rigel index``'s options. Every overridable default is ``None``, so ``resolve_index_args`` can tell a flag the
    user gave from one the YAML or :class:`~rigel.config.IndexConfig` supplies; every default a help string names is
    read from the config, never restated here."""
    from .config import FeatherCompression, GtfParseMode, IndexConfig

    idx = subparsers.add_parser("index", help="Build reference index from FASTA + GTF")
    idx.add_argument("--fasta", help="Genome FASTA file (must be indexed with samtools faidx)")
    idx.add_argument("--gtf", help="Gene annotation GTF file (GENCODE recommended)")
    idx.add_argument("-o", "--output-dir", help="Output directory for index files")
    idx.add_argument(
        "--config",
        help=(
            "YAML configuration file: the same keys as these options, with underscores, at the top level or in an "
            "'index:' block (a 'quant:' block is left to 'rigel quant', so one file can hold both). CLI flags "
            "override it."
        ),
    )
    idx.add_argument(
        "--feather-compression",
        choices=list(get_args(FeatherCompression)),
        help=f"Compression for Feather files (default: {IndexConfig.feather_compression})",
    )
    idx.add_argument(
        "--no-tsv",
        action="store_true",
        default=None,
        help="Skip writing human-readable TSV mirror files",
    )
    idx.add_argument(
        "--gtf-parse-mode",
        choices=list(get_args(GtfParseMode)),
        help=(
            "GTF parse mode: 'strict' fails on a malformed record; "
            f"'warn-skip' logs a warning and skips it (default: {IndexConfig.gtf_parse_mode})"
        ),
    )
    idx.add_argument(
        "--collapse-duplicate-transcripts",
        action="store_true",
        default=None,
        help=(
            "Collapse transcripts sharing identical exon coordinates (same ref, "
            "strand, and every exon start/end) instead of failing. Keeps the "
            "lexicographically-smallest transcript ID per group and drops the rest. "
            "Such transcripts are mathematically unidentifiable in quantification, "
            "so collapsing is loss-free. Default: fail with a report of the "
            "offending groups. Useful for GENCODE, which contains a handful of "
            "byte-identical transcript annotations."
        ),
    )
    idx.add_argument(
        "--end-tolerance",
        type=int,
        help=(
            "Multi-exon transcripts whose start and end sites each lie within this many bp of another's share one "
            "synthetic span, the nascent RNA they model, transitively; a single-exon transcript within it of a span "
            f"is that span, being its own nascent RNA (default: {IndexConfig.end_tolerance})"
        ),
    )
    idx.add_argument(
        "--alignable-zarr",
        help=(
            "Path to an alignable Zarr store built for the same genome+aligner. "
            "Provides the unobservable genome: the bases no read pair can be placed "
            "back at, at the store's longest read length. "
            "Required unless --no-mappability is set explicitly."
        ),
    )
    idx.add_argument(
        "--no-mappability",
        action="store_true",
        default=None,
        help=(
            "Explicitly opt out of the alignable store: the unobservable genome is the "
            "FASTA's N runs alone. Use this for synthetic "
            "genomes or any setting where running "
            "'alignable' is unnecessary. Mutually exclusive with --alignable-zarr."
        ),
    )
    idx.set_defaults(func=index_command)


def _add_quant_parser(subparsers) -> None:
    """``rigel quant``'s options, by group. Every overridable default is ``None``, so ``resolve_quant_args`` can tell
    a flag the user gave from one the YAML or the config supplies; every default a help string names is read from the
    config, never restated here."""
    quant = subparsers.add_parser(
        "quant",
        help="Single-pass BAM scan and Bayesian transcript abundance estimation",
    )
    defaults = build_quant_defaults()
    _add_quant_io_options(quant.add_argument_group("input / output"))
    _add_quant_alignment_options(quant.add_argument_group("alignment options"), defaults)
    _add_quant_model_options(quant.add_argument_group("model parameters"), defaults)
    _add_quant_performance_options(quant.add_argument_group("performance"), defaults)
    _add_quant_advanced_options(quant.add_argument_group("advanced options"), defaults)
    quant.set_defaults(func=quant_command)


def _yes_no(flag: bool) -> str:
    return "yes" if flag else "no"


def _add_quant_io_options(grp) -> None:
    grp.add_argument(
        "--bam",
        dest="bam_file",
        help="Name-sorted or collated paired-end BAM file. STAR and minimap2 "
        "are supported. A read pair's multiplicity is its usable alignments; "
        "the NH tag is read only to find a truncated output.",
    )
    grp.add_argument(
        "--index",
        dest="index_dir",
        help="Directory containing rigel index files (from 'rigel index').",
    )
    grp.add_argument("-o", "--output-dir", help="Output directory for quantification results.")
    grp.add_argument(
        "--config",
        help="YAML configuration file. Accepts the same keys as CLI "
        "options (with underscores). CLI flags override config file "
        "values. A config.yaml is written to the output directory "
        "after each run for easy reproducibility.",
    )
    grp.add_argument(
        "--annotated-bam",
        help="Write an annotated BAM with per-fragment assignment tags "
        "(ZT, ZG, ZR, ZI, ZJ, ZF, ZW, ZC, ZH, ZN, ZS, ZL) to this path. "
        "Requires a second pass over the BAM.",
    )
    grp.add_argument(
        "--tsv",
        action=argparse.BooleanOptionalAction,
        help=f"Also write TSV (.tsv) mirrors of quant tables (default: {_yes_no(_CONFIG_FILE_KEYS['tsv'])}).",
    )


def _add_quant_alignment_options(grp, qd: dict) -> None:
    grp.add_argument(
        "--include-multimap",
        action=argparse.BooleanOptionalAction,
        help="Include read pairs with more than one usable alignment "
        f"(default: {_yes_no(qd['include_multimap'])}). Off, each is excluded, counted in summary.json and "
        "written to the annotated BAM as dropped.",
    )
    grp.add_argument(
        "--keep-duplicates",
        action=argparse.BooleanOptionalAction,
        help=f"Keep reads marked as PCR/optical duplicates (default: {_yes_no(qd['keep_duplicates'])}).",
    )
    grp.add_argument(
        "--sj-strand-tag",
        nargs="+",
        help=f"BAM tag(s) for splice-junction strand (default: {' '.join(qd['sj_strand_tag'])}). "
        "'auto' reads XS and minimap2's ts on every spliced record; when both are present and disagree, the "
        "junction's strand is ambiguous. Name 'XS' (HISAT2, or STAR with --outSAMstrandField intronMotif) or "
        "'ts' (minimap2) to read only that tag, or list two to check in order (e.g. XS ts); 'none' reads no tag.",
    )


def _add_quant_model_options(grp, qd: dict) -> None:
    from .config import AssignmentMode, EMMode

    grp.add_argument(
        "--seed",
        type=int,
        help=f"Seed of the sampled assignment's draw (default: {qd['seed']}); another seed is another draw from "
        "the same posterior.",
    )
    grp.add_argument(
        "--em-iterations",
        type=int,
        help=f"EM budget in EM steps, at least 1 (default: {qd['em_iterations']}, a measured early stop). The EM runs "
        "accelerated (SQUAREM) steps of three EM updates each, up to a third of this budget and at least one, until "
        "it converges.",
    )
    grp.add_argument(
        "--em-mode",
        choices=list(get_args(EMMode)),
        help=f"EM algorithm variant (default: {qd['em_mode']}). "
        "'vbem' uses Variational Bayes EM with digamma-based soft "
        "updates. 'map' uses MAP-EM: each step adds the calibration's gDNA and RNA "
        "pseudocounts to the components' fragment counts and normalises.",
    )
    grp.add_argument(
        "--assignment-mode",
        choices=list(get_args(AssignmentMode)),
        help=f"Post-EM fragment assignment (default: {qd['assignment_mode']}). 'sample' gives every fragment to one "
        "transcript: each transcript's fractional count is rounded within its locus, then each fragment "
        "is drawn from its own posterior toward the transcripts still short of their count. "
        "'fractional' keeps each fragment's posterior weights.",
    )


def _add_quant_performance_options(grp, qd: dict) -> None:
    from .config import CONSTANTS

    grp.add_argument(
        "--threads",
        type=int,
        help=f"Total thread budget for Rigel (default: {qd['threads']} = all available "
        "cores). During BAM scan, this budget is split between scan "
        "workers and --scan-bgzf-threads; the calibration solve and the locus EM use the same budget "
        "because the stages run serially.",
    )
    grp.add_argument(
        "--scan-bgzf-threads",
        type=int,
        help="BGZF decompression threads reserved from --threads during "
        "BAM scan (default: derived from the budget, one per "
        f"{CONSTANTS.resources.scan_workers_per_bgzf_thread} threads). "
        "Set to 0 to disable htslib threaded decompression.",
    )
    grp.add_argument(
        "--scan-buffer-size",
        type=float,
        help="Memory in GiB for the scanned fragments held for quantification before the oldest chunks "
        f"spill to disk (default: {qd['scan_buffer_size']:g}); 0 disables spilling, holding every chunk in "
        "memory. It bounds that buffer only, not the chunks the scan workers hold in flight. Increase if you have "
        "ample RAM to reduce spill I/O.",
    )
    grp.add_argument(
        "--tmpdir",
        help="Directory for temporary buffer spill files when memory "
        "limits are exceeded (default: system temp directory).",
    )
    grp.add_argument(
        "--scan-fragments-per-chunk",
        type=int,
        help="Buffered fragments per scan chunk before the native scanner "
        f"hands data to the Python fragment buffer (default: {qd['scan_fragments_per_chunk']}).",
    )
    grp.add_argument(
        "--scan-read-name-batch-size",
        type=int,
        help="Advanced: read-name groups per native scanner input queue "
        f"item (default: {qd['scan_read_name_batch_size']}). Lower values are useful for queue-boundary "
        "debugging; larger values may reduce queue overhead at the cost "
        "of slightly more in-flight memory.",
    )


def _add_quant_advanced_options(grp, qd: dict) -> None:
    grp.add_argument(
        "--em-convergence-delta",
        type=float,
        help=f"Convergence threshold for EM parameter updates (default: {qd['em_convergence_delta']:g}).",
    )
    grp.add_argument(
        "--calib-refit-iters",
        type=int,
        help=f"Calibration prior-bootstrap iterations (default {qd['calib_refit_iters']}). Each one re-fits the population gDNA "
        "landscape on the current solve, resets the belief, and re-solves; later iterations change the "
        "answer little, and each costs one extra full sweep. 0 = the prior-free pass-0 alone, which fits no gDNA landscape and so applies no "
        "capture correction. Advanced calibration knob.",
    )
    grp.add_argument(
        "--sweep-block-slots",
        type=int,
        help="Calibration working set: the chain is solved one locus block at a time, the pieces between "
        f"intergenic regions merged up to this many slots per block (default {qd['sweep_block_slots']}; the whole "
        "chain as one block is set only from Python, CalibrationConfig(sweep_block_slots=None)). "
        "Performance only — the answer is the same for every value; it sizes the solver's per-thread "
        "arena, so smaller blocks use less memory per sweep. Advanced "
        "calibration knob.",
    )
    grp.add_argument(
        "--intron-departure-prior",
        type=float,
        help="The intron abundance model's departure prior, strictly between 0 and 1 "
        f"(default {qd['intron_departure_prior']:g}): the probability that a region's gDNA is not at the intergenic "
        "background (a probe's reach, a copy-number change, unannotated transcription). It caps how far the background "
        "can argue against an intron's own strand evidence, and keeps a region the background cannot explain from "
        "setting it. A stranded library measures its introns' own value from their strands, with this as its prior; "
        "an unstranded library uses this value, and answers are flat across a wide range around the default. Raising "
        "it toward 1 turns the intron abundance model off, for libraries whose intergenic regions are not a gDNA "
        "reference. Advanced calibration knob.",
    )
    grp.add_argument(
        "--overhang-alpha",
        type=float,
        help=f"Per-base overhang penalty alpha in [0,1] (default: {qd['overhang_alpha']:g}). "
        "0 = hard gate, 1 = no penalty.",
    )
    grp.add_argument(
        "--pruning-min-posterior",
        type=float,
        help="Minimum posterior threshold for candidate pruning, in [0,1), against the best "
        "candidate at the same alignment; a multimapper's alignments are pruned each on its own "
        f"(default: {qd['pruning_min_posterior']:g}). Lower values keep more candidates "
        "(conservative). Set to 0 to disable pruning entirely.",
    )
    grp.add_argument(
        "--intron-overhang-tolerance",
        type=int,
        help="Intron-overhang tolerance K (bp): how many bases an aligned read end may lie inside "
        "an annotated intron, at the edge of a gap the reads did not sequence or at a fragment's "
        "outer end, as when an aligner runs a read a base or two across a splice site. A transcript "
        "is contradicted only by an intron holding more than K aligned bases, and an intron counts "
        "as lying in an unsequenced gap when it extends past neither end of the gap by more than "
        "K bp. Along a witness transcript it decides whether two mates are one molecule; it also "
        "shapes the candidates' lengths and the gap hypotheses calibration arbitrates. Junctions "
        "are matched exactly. "
        f"Default: {qd['intron_overhang_tolerance']}.",
    )
    grp.add_argument(
        "--emit-locus-stats",
        action=argparse.BooleanOptionalAction,
        help="Write per-locus EM convergence profiling data to "
        f"locus_stats.feather in the output directory (default: {_yes_no(qd['emit_locus_stats'])}). Includes "
        "iteration counts, timing, and equivalence class statistics "
        "for every locus. Useful for debugging convergence.",
    )


def _add_sim_parser(subparsers) -> None:
    from .config import CONSTANTS

    defaults = CONSTANTS.simulator
    sim = subparsers.add_parser("sim", help="Generate synthetic test scenarios")
    sim.add_argument(
        "--config", required=True, help="YAML configuration file defining the scenario"
    )
    sim.add_argument(
        "-o", "--output-dir", required=True, help="Output directory for scenario artifacts"
    )
    sim.add_argument(
        "--genome-length",
        type=int,
        default=defaults.scenario_genome_length,
        help=f"Genome length in bases (default: {defaults.scenario_genome_length}, overridden by YAML)",
    )
    sim.add_argument(
        "--seed",
        type=int,
        default=defaults.seed,
        help=f"Random seed (default: {defaults.seed}, overridden by YAML)",
    )
    sim.add_argument(
        "--num-reads",
        type=int,
        default=defaults.scenario_fragments,
        help=f"Number of fragments to simulate (default: {defaults.scenario_fragments}, overridden by YAML)",
    )
    sim.set_defaults(func=sim_command)


def _add_export_parser(subparsers) -> None:
    export = subparsers.add_parser("export", help="Convert feather output files to TSV or Parquet")
    export.add_argument(
        "output_dir", help="Directory containing .feather output files from rigel quant"
    )
    export.add_argument(
        "-f",
        "--format",
        default="tsv",
        choices=["tsv", "parquet"],
        help="Output format (default: tsv)",
    )
    export.set_defaults(func=export_command)


def _add_report_parser(subparsers) -> None:
    report = subparsers.add_parser(
        "report",
        help="Build a self-contained HTML QC report from a rigel quant output directory",
    )
    report.add_argument("output_dir", help="A rigel quant --output-dir (must contain summary.json)")
    report.add_argument(
        "-o", "--output", help="Output HTML path (default: <output_dir>/report.html)"
    )
    report.add_argument("--title", help="Report title (default: 'Rigel QC · <sample>')")
    report.set_defaults(func=report_command)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main() -> None:
    """CLI entry point registered as ``rigel`` in pyproject.toml."""
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
