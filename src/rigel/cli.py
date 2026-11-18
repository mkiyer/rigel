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
import os
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .config import PipelineConfig


def get_version() -> str:
    from . import __version__

    return f"rigel {__version__}"


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
            nascent_tolerance=args.nascent_tolerance,
            alignable_zarr_path=(Path(args.alignable_zarr) if args.alignable_zarr else None),
            splice_blacklist_min_count=args.splice_blacklist_min_count,
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
        resolve_quant_args(args, build_quant_defaults())
        bam_path, index_dir, output_dir = _quant_paths(args)
        output_dir.mkdir(parents=True, exist_ok=True)
        # YAML may deliver a bare string; the CLI always delivers a list
        if isinstance(args.sj_strand_tag, str):
            args.sj_strand_tag = [args.sj_strand_tag]
        pipeline_config = build_pipeline_config(args)
    # a TypeError is a setting of the wrong type, as YAML's string `1e3` for an integer
    except (ValueError, TypeError) as exc:
        sys.exit(f"Error: {exc}")
    _check_annotated_bam_target(pipeline_config.annotated_bam_path, bam_path)

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


def _check_annotated_bam_target(annotated, bam_path: Path) -> None:
    """Refuse, before the run, an annotated BAM the second pass could not write, or one that is the input BAM: that
    pass runs after the quantification, so the failure would otherwise come only at the end, and over the input it
    would destroy it."""
    from .pipeline import BamInputError, check_annotated_bam_path

    if annotated is None:
        return
    try:
        check_annotated_bam_path(bam_path, annotated)
    except BamInputError as exc:
        sys.exit(f"Error: {exc}")
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

#: summary.json schema version, which a consumer should read before parsing. The fragment-length
#: histograms live in ``fragment_lengths.feather``, not the JSON; the JSON carries the splice and
#: strand-contamination diagnostics; ``gdna_em_effective_length`` summarises the loci table's column of that name;
#: ``length_tables`` is the EM's two length tables and the decision behind them.
SUMMARY_SCHEMA_VERSION = 8

#: summary.json's floats (fractions, strand rates, quantiles) are rounded to this many decimals.
SUMMARY_DECIMALS = 6

#: summary.json's fragment counts are rounded to this many decimals: a hundredth of a fragment.
SUMMARY_COUNT_DECIMALS = 2

#: summary.json's lengths in bp (a length law's mean, median and standard deviation) are rounded to this many
#: decimals: a hundredth of a base.
SUMMARY_LENGTH_DECIMALS = 2

#: summary.json's strand posterior variance keeps this many significant digits: it is p(1 − p)/n, of order 1e-8 at a
#: million fragments, where fixed decimals print 0, and at most 0.25, so 8 significant digits are at least 8 decimals.
SUMMARY_VARIANCE_DIGITS = 8

#: The tables every run writes, each as ``<name>.feather`` with a TSV mirror under ``--tsv``, in the order written.
_QUANT_TABLES = ("quant", "gene_quant", "nascent_quant", "loci")

#: ``fragment_lengths.feather``'s pure-pool categories, in ``FragmentPool`` order. Named here, not derived from the
#: enum, because they are an output contract: a rename in the enum must not rename a reported category.
_POOL_REPORT_NAMES = (
    "pool_gdna_intergenic",
    "pool_gdna_intronic",
    "pool_gdna_intron_exon",
    "pool_gdna_intergenic_exon",
    "pool_rna_spliced",
)

#: What a store-built index records of the store its unobservable mask was read from (`rigel.observable`).
_MASK_STORE_KEYS = (
    "aligner",
    "aligner_version",
    "read_length_bins",
    "frag_len_mean",
    "frag_len_sd",
    "error_rate",
    "references_from_fasta",
)


def write_quant_outputs(result, index, output_dir: Path, args: argparse.Namespace) -> None:
    """Write everything ``rigel quant`` produces: the four quantification tables, ``locus_stats.feather`` when the
    run collected it, ``summary.json``, ``fragment_lengths.feather``, the calibration track and ``config.yaml``."""
    import json

    estimator = result.estimator
    _remove_stale_outputs(output_dir, args.tsv, bool(estimator.locus_stats))
    for key, val in sorted(result.stats.to_dict().items()):
        if isinstance(val, (int, float)):
            logging.info(f"  {key}: {val:,}")

    tables = _write_tables(estimator, index, output_dir, args.tsv)
    fl_summary, fl_histogram_df = fragment_length_report(result.fl_models)

    summary_path = output_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(_summary(result, index, output_dir, args, tables, fl_summary), f, indent=2)
    logging.info(f"[DONE] Summary written to {summary_path}")

    # The raw fragment-length histograms, long form: one (category, length, count) row per 1-bp bin
    fragment_lengths_path = output_dir / "fragment_lengths.feather"
    _write_table(fl_histogram_df, fragment_lengths_path, args.tsv)
    logging.info(f"[DONE] Wrote {fragment_lengths_path} ({len(fl_histogram_df)} rows)")

    _write_calibration_track(result.calibration_track, output_dir)
    write_config_yaml(output_dir / "config.yaml", args)


def _remove_stale_outputs(output_dir: Path, tsv: bool, locus_stats: bool) -> None:
    """A reused output directory keeps no product this run does not write: a quant.tsv from an earlier --tsv run
    beside this run's quant.feather would disagree with it silently."""
    stale = []
    if not tsv:
        stale += [output_dir / f"{name}.tsv" for name in (*_QUANT_TABLES, "fragment_lengths")]
    if not locus_stats:
        stale.append(output_dir / "locus_stats.feather")
    for path in stale:
        path.unlink(missing_ok=True)


def _write_table(df, path: Path, tsv: bool) -> None:
    """One output table: Feather v2 with zstd, and its TSV mirror under ``--tsv``."""
    df.to_feather(str(path), compression="zstd")
    if tsv:
        df.to_csv(str(path.with_suffix(".tsv")), sep="\t", index=False)


def _write_tables(estimator, index, output_dir: Path, tsv: bool) -> dict:
    """Write the four quantification tables, and the per-locus EM profile when the run collected one; return the
    tables by name."""
    import pandas as pd

    tables = {
        "quant": estimator.get_counts_df(index),
        "gene_quant": estimator.get_gene_counts_df(index),
        "nascent_quant": estimator.get_nascent_quant_df(index),
        "loci": estimator.get_loci_df(index),
    }
    paths = [output_dir / f"{name}.feather" for name in _QUANT_TABLES]
    for name, path in zip(_QUANT_TABLES, paths):
        _write_table(tables[name], path, tsv)
    logging.info(f"[DONE] Wrote {', '.join(str(path) for path in paths)}")

    if estimator.locus_stats:
        locus_stats_path = output_dir / "locus_stats.feather"
        locus_stats_df = pd.DataFrame(estimator.locus_stats)
        locus_stats_df.to_feather(str(locus_stats_path), compression="zstd")
        logging.info(f"[DONE] Wrote {locus_stats_path} ({len(locus_stats_df)} loci)")
    return tables


def _write_calibration_track(track, output_dir: Path) -> None:
    """The genome-wide gDNA track, calibration's per-region gDNA level: a Feather (the report's input) and a bedGraph
    (for a genome browser)."""
    if track is None or len(track) == 0:
        return
    from .calibration.track import write_bedgraph

    track_path = output_dir / "calibration_track.feather"
    bedgraph_path = output_dir / "calibration_track.bedgraph"
    track.to_feather(str(track_path), compression="zstd")
    write_bedgraph(track, bedgraph_path)
    logging.info(f"[DONE] Wrote {track_path} + {bedgraph_path.name} ({len(track)} regions)")


def _summary(result, index, output_dir: Path, args, tables: dict, fl_summary: dict) -> dict:
    """``summary.json``: the run's command and configuration, then what the scan saw, the strand and calibration
    answers, the length laws and the quantification's pools."""
    from . import __version__

    stats = result.stats
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "rigel_version": __version__,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "command": {
            "subcommand": "quant",
            "arguments": _command_arguments(args, output_dir),
            "config_file": str(Path(args.config).resolve()) if args.config else None,
        },
        "configuration": result.pipeline_config.to_dict(),
        "input": {
            "bam_file": str(Path(args.bam_file).resolve()),
            "index_dir": str(Path(args.index_dir).resolve()),
        },
        # the observable genome every opportunity of the run was read over: the index's mask and where it came from
        "observable_genome": observable_genome_summary(index),
        "alignment_stats": _alignment_summary(stats),
        "read_pairs": _read_pair_summary(stats),
        "fragment_stats": _fragment_summary(stats, index),
        "strand_model": _strand_model_summary(result.strand_models),
        "calibration": _calibration_summary(result.calibration, result.junction_cap),
        "gdna_em_effective_length": _series_summary(tables["loci"], "gdna_em_effective_length"),
        "fragment_length": fl_summary,
        "length_tables": _length_tables_summary(result.length_tables),
        "quantification": _quantification_summary(
            stats, result.estimator, index, float(tables["quant"]["count"].sum())
        ),
    }


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


def observable_genome_summary(index) -> dict:
    """The observable genome the run read (`rigel.observable`): the index's record of what its unobservable mask is —
    ``rule``, and the store's ``read_length`` bin it reads (``None`` without a store) — the observable and unobservable
    bases over the genome and ``by_reference``, and the ``store`` the mask was read from (``None`` without one, when the
    mask is the reference's unknown bases)."""
    from .index import load_manifest

    record = (load_manifest(index.index_dir) or {}).get("observable") or {}
    df = index.unobservable_df
    masked = {ref: 0 for ref in index.ref_lengths}
    if df is not None and len(df):
        for ref, n in (df["end"] - df["start"]).groupby(df["ref"].astype(str)).sum().items():
            masked[str(ref)] += int(n)
    by_reference = {
        ref: {"observable_bases": int(n) - masked[ref], "unobservable_bases": masked[ref]}
        for ref, n in index.ref_lengths.items()
    }
    return {
        "rule": record.get("rule"),
        "read_length": record.get("read_length"),
        "observable_bases": sum(v["observable_bases"] for v in by_reference.values()),
        "unobservable_bases": sum(v["unobservable_bases"] for v in by_reference.values()),
        "store": {k: record[k] for k in _MASK_STORE_KEYS} if "aligner" in record else None,
        "by_reference": by_reference,
    }


def _alignment_summary(stats) -> dict:
    """The records: each flag counter counts every record carrying its flag (``proper_pairs`` the primary records
    flagged as a proper pair), and every record has one status."""
    return {
        "total_reads": stats.total,
        "mapped_reads": stats.total - stats.unmapped,
        "unmapped_reads": stats.unmapped,
        "secondary_reads": stats.secondary,
        "supplementary_reads": stats.supplementary,
        "proper_pairs": stats.proper_pair,
        "duplicate_reads": stats.duplicate,
        "qc_fail_reads": stats.qc_fail,
        # One status per record, summing to ``total_reads``: filtered (unmapped, QC-failed or a skipped duplicate),
        # used (a usable alignment's), unused (in a placed read pair, in no usable alignment), or unplaced (in a read
        # pair withheld or excluded by the multimap option)
        "records": {
            "filtered": stats.filtered,
            "used": stats.n_records_used,
            "unused": stats.n_records_unused,
            "unplaced": stats.n_records_unplaced,
        },
        # Each spliced record's junction-strand tag, read on that record: given a strand, given none, or ambiguous
        # because its XS and ts disagree; ``invalid`` counts the tag values that were neither '+' nor '-', each read
        # as absent.
        "sj_strand_tag": {
            "tagged": stats.n_sj_tag_tagged,
            "untagged": stats.n_sj_tag_untagged,
            "conflicting": stats.n_sj_tag_conflicting,
            "invalid": stats.n_sj_tag_invalid,
        },
        # NM and NH on each record: ``missing`` without the tag, ``invalid`` holding it as anything but an integer,
        # each read as absent (NM 0, NH 1); ``conflicting`` counts the alignments whose records disagree on NM, each
        # read with the larger
        "nm_tag": {
            "missing": stats.n_nm_tag_missing,
            "invalid": stats.n_nm_tag_invalid,
            "conflicting": stats.n_nm_conflicting,
        },
        "nh_tag": {"missing": stats.n_nh_tag_missing, "invalid": stats.n_nh_tag_invalid},
    }


def _read_pair_summary(stats) -> dict:
    """One outcome per read pair (a run of records sharing a name), summing to ``total``: excluded (every record
    filtered), withheld with no usable alignment by reason (an unpaired one also by the furthest check any of its
    alignments reached), excluded by the multimap option, or placed, unique or multimapping."""
    return {
        "total": stats.n_read_pairs,
        "excluded": stats.n_excluded,
        "withheld": {
            "unpaired": stats.n_withheld_unpaired,
            "unannotated_splice": stats.n_withheld_unannotated_splice,
            "off_index": stats.n_withheld_off_index,
            "chimeric": stats.n_withheld_chimeric,
        },
        "unpaired_by_reason": {
            "missing_mate": stats.n_withheld_unpaired_missing_mate,
            "reference": stats.n_withheld_unpaired_reference,
            "strand": stats.n_withheld_unpaired_strand,
            "outward": stats.n_withheld_unpaired_outward,
            "length": stats.n_withheld_unpaired_length,
        },
        "excluded_by_option": stats.n_excluded_by_option,
        "placed": {
            "unique": stats.n_placed_unique,
            "multimapping": stats.n_placed_multimapping,
        },
    }


def _fragment_summary(stats, index) -> dict:
    """The fragments: genic and intergenic, the hits rejected by reason, the splice census and the mismatch rate."""
    rate = stats.mismatch_rate
    return {
        "total": stats.n_fragments,
        "genic": stats.n_fragments - stats.n_intergenic,
        "intergenic": stats.n_intergenic,
        "with_annotated_sj": stats.n_with_annotated_sj,
        # Templates whose NH exceeds the alignments reported for them (a truncated output)
        "nh_above_reported": stats.n_nh_above_reported,
        "hits_rejected": {
            "missing_mate": stats.n_hits_rejected_missing_mate,
            "reference": stats.n_hits_rejected_reference,
            "strand": stats.n_hits_rejected_strand,
            "outward": stats.n_hits_rejected_outward,
            "length": stats.n_hits_rejected_length,
            "contradictory": stats.n_hits_rejected_contradictory,
            "novel_junction": stats.n_hits_rejected_novel_junction,
            "unsupported": stats.n_hits_rejected_unsupported,
            "off_index": stats.n_hits_rejected_off_index,
            "chimeric": stats.n_hits_rejected_chimeric,
        },
        "splice": _splice_summary(stats, index),
        # The library's mismatch rate, read on its unique unspliced hits (their NM over the read bases they align):
        # the price of a multimapper's hit's mismatches and soft-clipped bases; null with no such hit
        "mismatch_rate": {
            "rate": None if rate is None else round(rate, SUMMARY_DECIMALS),
            "mismatches": stats.n_unique_unspliced_mismatches,
            "bases": stats.n_unique_unspliced_bases,
        },
    }


def _splice_summary(stats, index) -> dict:
    """The splice-type breakdown: the scanner's per-fragment census, which counts every fragment it offers the
    accumulator (and closes against the accumulator's qc block), plus the blacklist tally. It reports the
    ``spliced_implicit`` and ``splice_artifact`` classes the ``with_*_sj`` counters omit. The blacklist's size tells
    "detection off" (no blacklist in the index: 0) from "detection on, 0 artifacts found"."""
    from .splice import SpliceType, census_field

    blacklist_size = int(index.sj_blacklist_size)
    census = {
        "unspliced": SpliceType.UNSPLICED,
        "spliced_annotated": SpliceType.SPLICED_ANNOT,
        "spliced_implicit": SpliceType.SPLICED_IMPLICIT,
        "splice_artifact": SpliceType.SPLICE_ARTIFACT,
    }
    return {
        **{label: int(getattr(stats, census_field(stype))) for label, stype in census.items()},
        "sj_blacklisted": int(stats.n_sj_blacklisted),
        "sj_blacklist_size": blacklist_size,
        "sj_blacklist_loaded": bool(blacklist_size),
    }


def _strand_model_summary(sm) -> dict:
    """The strand model: the protocol and its specificity from the spliced reads, the junction table behind it, the
    junction fit (:func:`strand_fit_summary`) and the all-exonic diagnostic."""
    spliced = sm.exonic_spliced
    ci_lo, ci_hi = spliced.posterior_95ci()
    return {
        "protocol": "R1-sense" if sm.read1_sense else "R1-antisense",
        "strand_specificity": round(sm.strand_specificity, SUMMARY_DECIMALS),
        "p_r1_sense": round(sm.p_r1_sense, SUMMARY_DECIMALS),
        "read1_sense": bool(sm.read1_sense),
        "n_training_fragments": sm.n_observations,
        "posterior_variance": _significant(spliced.posterior_variance(), SUMMARY_VARIANCE_DIGITS),
        "ci_95": [round(ci_lo, SUMMARY_DECIMALS), round(ci_hi, SUMMARY_DECIMALS)],
        # The per-sj strand table this 2×2 is the marginal of, and the population the junction fit reads κ from.
        "sj": sm.sj_table.to_dict(),
        **strand_fit_summary(sm),
        # The all-exonic (RNA + gDNA) diagnostic model, never used for scoring. A positive ``contamination_gap`` means
        # unstranded gDNA pulls the all-exonic estimate toward 0.5 relative to the spliced-only model: a sign of gDNA
        # contamination.
        "diagnostics": {
            "exonic_all_specificity": round(sm.exonic.strand_specificity, SUMMARY_DECIMALS),
            "exonic_all_p_r1_sense": round(sm.exonic.p_r1_sense, SUMMARY_DECIMALS),
            "exonic_all_n_fragments": sm.exonic.n_observations,
            "contamination_gap": round(
                spliced.strand_specificity - sm.exonic.strand_specificity, SUMMARY_DECIMALS
            ),
        },
    }


def _significant(x: float, digits: int) -> float:
    """``x`` rounded to ``digits`` significant digits."""
    return float(f"{x:.{digits}g}")


def strand_fit_summary(strand_models) -> dict:
    """``summary.json``'s view of the junction fit behind κ.

    κ, the spliced reads' sense fraction, pooled over every junction and over the junctions the three-class fit
    credits to genuine RNA; each class's share of the junctions and of their fragments (shallow splice artifacts
    are most of the junctions and few of the fragments); whether the strand channel is live, the protocol decision
    calibration's strand channel reads, made here by the same two functions on the same 2×2; and the strand
    overdispersion, 0 by policy for gDNA and RNA alike, so the channel is binomial. Where the channel is dead the
    classes are not identified and calibration reads the pooled counts (``StrandModel.genuine_n_same``), so the
    genuine κ is the pooled one and there are no class shares."""
    from .calibration.strand_balance import fit_strand_balance, strand_discriminability
    from .junction_fit import junction_classes

    table = strand_models.sj_table
    n_sense, depth = float(table.n_sense.sum()), float(table.depth.sum())
    pooled = round((n_sense + 1.0) / (depth + 2.0), SUMMARY_DECIMALS)
    names = ("genuine", "artifact", "reversed")

    def shares(values):
        return (
            None
            if values is None
            else {k: round(v, SUMMARY_DECIMALS) for k, v in zip(names, values)}
        )

    balance = fit_strand_balance(strand_models)
    live = bool(strand_discriminability(balance.rna_sense_frac, balance.n_observations) > 0.0)
    if live:
        classes = junction_classes(table.n_sense, table.depth)
        genuine = round(classes.kappa, SUMMARY_DECIMALS)
        by_junction, by_fragment = shares(classes.by_junction), shares(classes.by_fragment)
    else:
        genuine, by_junction, by_fragment = pooled, None, None
    return {
        "kappa": {"pooled": pooled, "genuine": genuine},
        "junction_classes": {"by_junction": by_junction, "by_fragment": by_fragment},
        "strand_live": live,
        "strand_overdispersion": {"value": 0.0, "model": "binomial"},
    }


def _calibration_summary(cal, junction_cap: float) -> dict | None:
    """Calibration's library scalars and the junction cap fitted from them, never a per-region dict."""
    if cal is None:
        return None
    return {
        # full precision: a density is a few fragments per Mb, which a fixed count of decimals rounds to 0
        "gdna_density_global": float(cal.gdna_density_global),
        # calibration's capture answer, which the ruler reads: the gDNA bulk density the capture efficiencies are
        # read against (None = no gDNA landscape, nothing contracted), the weight the reading is raised to (the share
        # of the reads whose own gDNA would show; 1 keeps the reading, 0 claims no capture) and the landscape's
        # saturated level against that bulk at that weight (1 = no level above it)
        "gdna_bulk_density": (
            None if cal.gdna_bulk_density is None else float(cal.gdna_bulk_density)
        ),
        "gdna_capture_weight": round(float(cal.gdna_capture_weight), SUMMARY_DECIMALS),
        "gdna_capture_saturation": float(cal.gdna_capture_saturation),
        # the cap on a junction's capture price in the EM's lengths: the saturated level, raised as far as RNA's own
        # spliced fragments support (equal to it where they raise nothing)
        "junction_capture_cap": float(junction_cap),
        "fraction_rna_sense": round(float(cal.rna_sense_frac), SUMMARY_DECIMALS),
        "n_regions": int(cal.n_regions),
        "n_boundaries": int(cal.n_boundaries),
        "n_sj": int(cal.n_sj),
    }


def _series_summary(df, column: str) -> dict[str, float]:
    """A table column's min, median, 95th percentile and max; zeros when the table is empty or lacks it."""
    if column not in df.columns or df.empty:
        return {"min": 0.0, "median": 0.0, "p95": 0.0, "max": 0.0}
    values = df[column]
    return {
        "min": round(float(values.min()), SUMMARY_DECIMALS),
        "median": round(float(values.median()), SUMMARY_DECIMALS),
        "p95": round(float(values.quantile(0.95)), SUMMARY_DECIMALS),
        "max": round(float(values.max()), SUMMARY_DECIMALS),
    }


def _length_tables_summary(tables) -> dict | None:
    """The EM's two length tables and the decision behind them (`calibration.shared_length_law`): ``weight`` 0 is one
    shared law for both origins, 1 gDNA's own law brought into RNA's captured frame, continuous between. ``kappa`` is
    RNA's captured fraction, ``captured`` the capture mixture's fraction, ``contrast`` and ``bound`` the Hellinger
    distances the weight compares, ``noise`` the squared distance sampling noise alone adds; the three means are the
    two tables' and gDNA's uncaptured law's, in bp."""
    if tables is None:
        return None
    import numpy as np

    d = tables.decision
    L = np.arange(tables.rna.size, dtype=np.float64)
    return {
        "weight": round(float(d.weight), SUMMARY_DECIMALS),
        "kappa": round(float(d.kappa), SUMMARY_DECIMALS),
        "captured": round(float(d.captured), SUMMARY_DECIMALS),
        "contrast": round(float(d.contrast), SUMMARY_DECIMALS),
        "bound": round(float(d.bound), SUMMARY_DECIMALS),
        "noise": round(float(d.noise), SUMMARY_DECIMALS),
        "rna_table_mean": round(float(tables.rna @ L), SUMMARY_LENGTH_DECIMALS),
        "gdna_table_mean": round(float(tables.gdna @ L), SUMMARY_LENGTH_DECIMALS),
        "gdna_uncaptured_mean": round(float(tables.anchor @ L), SUMMARY_LENGTH_DECIMALS),
    }


def fragment_length_report(fl_models):
    """Split the fragment-length models into a lean summary + a tidy histogram table.

    Returns ``(summary_by_category, histogram_df)`` where

    * ``summary_by_category`` is ``{category: {n_observations, mean, std, median, mode, max_size, overflow_count,
      overflow_fraction}}`` (:func:`length_summary`) — small, human-readable, destined for ``summary.json`` (no raw
      per-bin arrays).
    * ``histogram_df`` is a long ``(category, length, count)`` DataFrame — the raw 1-bp histograms
      (:func:`length_histogram`), destined for ``fragment_lengths.feather``, which keeps them out of ``summary.json``.

    Every category measures one quantity, the accumulator's ``L``, over a stated population:

    ``global``
        the unconditional anchor: every deposited fragment, no purity condition.
    ``gdna`` / ``rna``
        the two modelled pools, the empirical views of what scoring and calibration actually use.
    the five ``pool_*`` rows
        the five pools unaggregated (none is pure). ``gdna`` is the fitted uniform-frame gDNA
        histogram and ``rna`` the de-tilted ``pool_rna_spliced``; the two ``*_exon`` crossing rows are
        the on-target gDNA, shown on their own so the off-target / on-target comparison is an output.

    There are no per-splice-type histograms: the scanner measures those by other rules than ``L``. The splice
    counts come from the scanner's census instead.
    """
    import pandas as pd

    models = {
        "global": fl_models.global_model(),
        "gdna": fl_models.gdna_model(),
        "rna": fl_models.rna_model(),
        **{
            name: fl_models._empirical(fl_models.pool_counts[i])
            for i, name in enumerate(_POOL_REPORT_NAMES)
        },
    }
    cats: list = []
    lengths: list = []
    counts: list = []
    for name, model in models.items():
        bins, values = length_histogram(model)
        cats += [name] * len(bins)
        lengths += bins.tolist()
        counts += values.tolist()
    histogram_df = pd.DataFrame(
        {
            "category": pd.Series(cats, dtype="category"),
            "length": pd.Series(lengths, dtype="int32"),
            "count": pd.Series(counts, dtype="float64"),
        }
    )
    return {name: length_summary(model) for name, model in models.items()}, histogram_df


def length_summary(model) -> dict:
    """A length law's statistics as ``summary.json`` reports them. ``mean``, ``std`` and ``mode`` describe the
    in-range lengths; the right-censored tail (lengths ``>= max_size``) is ``overflow_count`` and
    ``overflow_fraction``."""
    total = model.total_weight
    overflow = model.overflow_count
    return {
        "n_observations": int(model.n_observations),
        "mean": float(round(model.mean, SUMMARY_LENGTH_DECIMALS)),
        "std": float(round(model.std, SUMMARY_LENGTH_DECIMALS)),
        "median": float(round(model.median, SUMMARY_LENGTH_DECIMALS)),
        "mode": int(model.mode),
        "max_size": int(model.max_size),
        "overflow_count": float(round(overflow, SUMMARY_COUNT_DECIMALS)),
        "overflow_fraction": float(round(overflow / total, SUMMARY_DECIMALS)) if total > 0 else 0.0,
    }


def length_histogram(model):
    """``(lengths, counts)``: a length law's 1-bp histogram, trimmed to its first and last non-empty bins."""
    import numpy as np

    nonzero = np.flatnonzero(model.counts)
    if nonzero.size == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.float64)
    lo, hi = int(nonzero[0]), int(nonzero[-1]) + 1
    return np.arange(lo, hi), model.counts[lo:hi]


def _quantification_summary(stats, estimator, index, count_annotated: float) -> dict:
    """The quantification: what was counted and assigned, then the pools (:func:`quantification_pools`)."""
    is_synthetic = index.t_df["is_synthetic"].to_numpy(dtype=bool)
    n_nascent = int(is_synthetic.sum())
    return {
        # the annotated transcripts, single-exon included; ``n_nascent`` counts the synthetic spans
        "n_transcripts": int(is_synthetic.size) - n_nascent,
        "n_nascent": n_nascent,
        "n_genes": index.num_annotated_genes,
        "n_loci": len(estimator.locus_results),
        "n_unambig_assigned": stats.deterministic_unambig_units,
        "n_em_assigned": (
            stats.em_routed_unambig_units
            + stats.em_routed_ambig_same_strand_units
            + stats.em_routed_ambig_opp_strand_units
            + stats.em_routed_multimapper_units
        ),
        # dropped by a hard scoring gate (an alpha of 0): every candidate failed it, so no pool holds them
        "n_gated_out": stats.n_gated_out,
        **quantification_pools(
            stats,
            count_annotated,
            float(estimator.synthetic_count),
            float(estimator.gdna_em_count),
            float(estimator.gdna_intergenic_em_count),
        ),
    }


def quantification_pools(
    stats,
    count_annotated: float,
    count_nascent: float,
    count_gdna_em: float,
    count_gdna_intergenic_em: float,
) -> dict:
    """``summary.json``'s pools and their fractions of every fragment counted.

    An intergenic fragment reaches no locus. It is unspliced (a spliced hit the scan keeps has a transcript holding
    its junctions), so gDNA by construction: ``count_gdna`` is the EM's genic gDNA (``count_gdna_em``) plus the
    intergenic fragments; leaving them out would halve the observed contamination of a library with substantial
    gDNA. The annotated-transcript (``fraction_annotated``), synthetic-span (``fraction_nascent``) and gDNA fractions
    sum to one; ``fraction_gdna_em`` is the EM's genic part alone, and ``fraction_intergenic`` the intergenic part."""
    # The templates that reach no locus (``PipelineStats.n_intergenic``), and the gDNA the EM put at a multimapper's
    # intergenic place: neither is in a locus.
    intergenic = stats.n_intergenic + count_gdna_intergenic_em
    total_all = count_annotated + count_nascent + count_gdna_em + intergenic
    count_gdna = count_gdna_em + intergenic

    def share(x: float) -> float:
        return round(x / total_all, SUMMARY_DECIMALS) if total_all > 0 else 0.0

    return {
        "count_annotated": round(count_annotated, SUMMARY_COUNT_DECIMALS),
        "count_nascent": round(count_nascent, SUMMARY_COUNT_DECIMALS),
        "count_gdna": round(count_gdna, SUMMARY_COUNT_DECIMALS),
        "count_gdna_em": round(count_gdna_em, SUMMARY_COUNT_DECIMALS),
        "fraction_annotated": share(count_annotated),
        "fraction_nascent": share(count_nascent),
        "fraction_gdna": share(count_gdna),
        "fraction_gdna_em": share(count_gdna_em),
        "fraction_intergenic": share(intergenic),
    }


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
        ref_name=cfg.get("ref_name", "chr1"),
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
        else:
            out = fpath.with_suffix(".parquet")
            df.to_parquet(str(out), compression="zstd", index=False)
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
        "direct"  # "direct" | "invert_bool" | "sj_tag" | "path_or_none" | "gb_to_bytes"
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
    # -- Fan-out: total threads → the EM's, the scan's and calibration's budgets --
    ParamSpec("threads", "em.n_threads"),
    ParamSpec("threads", "scan.total_threads"),
    ParamSpec("threads", "calibration.n_threads"),
    # -- BamScanConfig: buffer sizing --
    ParamSpec("scan_buffer_size", "scan.buffer_size_bytes", "gb_to_bytes"),
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
    if transform == "gb_to_bytes":
        return val / (1024**3)  # bytes → GiB for CLI display
    raise ValueError(f"Unknown transform: {transform!r}")


def _cli_to_config(val: object, transform: str) -> object:
    """Transform a CLI-space value to config-space."""
    if transform == "direct":
        return val
    if transform == "invert_bool":
        return not val
    if transform == "sj_tag":
        if isinstance(val, str):
            return val
        return val[0] if len(val) == 1 else tuple(val)
    if transform == "path_or_none":
        return Path(val) if val else None
    if transform == "gb_to_bytes":
        return int(val * 1024**3)  # GiB → bytes
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


def resolve_quant_args(args: argparse.Namespace, defaults: dict) -> None:
    """Fill every option in place: an explicit CLI flag wins, then the ``--config`` YAML, then *defaults*.

    Every overridable argparse default is ``None``, so a non-``None`` value on *args* is one the user gave on the
    command line.
    """
    yaml_config = _read_config_file(args.config, defaults) if args.config else {}
    for key, default_val in {**_CONFIG_FILE_KEYS, **defaults}.items():
        if getattr(args, key) is not None:
            continue  # CLI wins
        yaml_val = yaml_config.get(key)
        setattr(args, key, default_val if yaml_val is None else yaml_val)


def _read_config_file(config_path, defaults: dict) -> dict:
    """A ``--config`` YAML's keys and values: hyphens read as underscores, a ``quant:`` block's keys under the
    top-level ones, and a key that is neither an option nor one of ``_CONFIG_FILE_KEYS`` refused."""
    import yaml

    with open(config_path) as fh:
        text = fh.read()
    raw = yaml.safe_load(text) or {}
    if not isinstance(raw, dict):
        raise ValueError(
            f"{config_path} must hold a mapping of keys to values, not a {type(raw).__name__}"
        )
    raw_normalized = {k.replace("-", "_"): v for k, v in raw.items()}
    quant_block = raw_normalized.pop("quant", None)
    yaml_config: dict = {}
    if quant_block is not None:
        if not isinstance(quant_block, dict):
            raise ValueError("quant config block must be a mapping.")
        yaml_config.update({str(k).replace("-", "_"): v for k, v in quant_block.items()})
    yaml_config.update(raw_normalized)
    valid_keys = set(defaults) | set(_CONFIG_FILE_KEYS)
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
    from .config import IndexConfig

    idx = subparsers.add_parser("index", help="Build reference index from FASTA + GTF")
    idx.add_argument(
        "--fasta",
        dest="fasta_file",
        required=True,
        help="Genome FASTA file (must be indexed with samtools faidx)",
    )
    idx.add_argument(
        "--gtf",
        dest="gtf_file",
        required=True,
        help="Gene annotation GTF file (GENCODE recommended)",
    )
    idx.add_argument("-o", "--output-dir", required=True, help="Output directory for index files")
    idx.add_argument(
        "--feather-compression",
        default=IndexConfig.feather_compression,
        choices=["lz4", "zstd", "uncompressed"],
        help="Compression for Feather files (default: %(default)s)",
    )
    idx.add_argument(
        "--no-tsv",
        action="store_true",
        default=not IndexConfig.write_tsv,
        help="Skip writing human-readable TSV mirror files",
    )
    idx.add_argument(
        "--gtf-parse-mode",
        default=IndexConfig.gtf_parse_mode,
        choices=["strict", "warn-skip"],
        help=(
            "GTF parse mode: 'strict' fails on a malformed record; "
            "'warn-skip' logs a warning and skips it (default: %(default)s)"
        ),
    )
    idx.add_argument(
        "--collapse-duplicate-transcripts",
        action="store_true",
        default=IndexConfig.collapse_duplicate_transcripts,
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
        "--nascent-tolerance",
        type=int,
        default=IndexConfig.nascent_tolerance,
        help=(
            "Max distance (bp) for clustering multi-exon transcripts' start and end sites "
            "into the synthetic spans that model their nascent RNA (default: %(default)s)"
        ),
    )
    idx.add_argument(
        "--alignable-zarr",
        help=(
            "Path to an alignable Zarr store built for the same genome+aligner. "
            "Provides the splice-junction artifact blacklist (for filtering "
            "spurious sj at BAM-scan time) and the unobservable genome: the bases "
            "no read pair can be placed back at, at the store's longest read length. "
            "Required unless --no-mappability is set explicitly."
        ),
    )
    idx.add_argument(
        "--no-mappability",
        action="store_true",
        help=(
            "Explicitly opt out of the alignable store: no splice blacklist, and the "
            "unobservable genome is the FASTA's N runs alone. Use this for synthetic "
            "genomes or any setting where running "
            "'alignable' is unnecessary. Mutually exclusive with --alignable-zarr."
        ),
    )
    idx.add_argument(
        "--splice-blacklist-min-count",
        type=int,
        default=IndexConfig.splice_blacklist_min_count,
        help=(
            "(Advanced) Minimum unique-fragment count per (chrom, intron, "
            "read_length) for a sj to enter the splice-artifact blacklist "
            "(default: %(default)s, the alignable tool's own threshold). Lower values "
            "(e.g. 1) admit more singleton artifacts; higher values keep only "
            "the most reproducible. Ignored when --no-mappability is set."
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
        "(ZT, ZG, ZR, ZI, ZJ, ZF, ZW, ZC, ZH, ZN, ZS, ZL, ZB) to this path. "
        "Requires a second pass over the BAM.",
    )
    grp.add_argument(
        "--tsv",
        action=argparse.BooleanOptionalAction,
        help="Also write TSV (.tsv) mirrors of quant tables (default: no).",
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
        "'ts' (minimap2) to read only that tag, or list two to check in order (e.g. XS ts).",
    )


def _add_quant_model_options(grp, qd: dict) -> None:
    grp.add_argument(
        "--seed",
        type=int,
        help=f"Seed of the sampled assignment's draw (default: {qd['seed']}); another seed is another draw from "
        "the same posterior.",
    )
    grp.add_argument(
        "--em-iterations",
        type=int,
        help=f"EM iteration budget (default: {qd['em_iterations']}). The EM runs accelerated (SQUAREM) steps of "
        "three EM updates each, up to a third of this budget, until it converges. At least one "
        "step always runs, so 0 does not skip the EM: every fragment is still assigned.",
    )
    grp.add_argument(
        "--em-mode",
        choices=["vbem", "map"],
        help=f"EM algorithm variant (default: {qd['em_mode']}). "
        "'vbem' uses Variational Bayes EM with digamma-based soft "
        "updates. 'map' uses MAP-EM: each step adds the calibration's gDNA and RNA "
        "pseudocounts to the components' fragment counts and normalises.",
    )
    grp.add_argument(
        "--assignment-mode",
        choices=["sample", "fractional"],
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
        "workers and --scan-bgzf-threads; locus EM uses the same budget "
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
        help="Maximum scan buffer size in GiB before chunks are spilled "
        f"to disk (default: {qd['scan_buffer_size']:g}). Increase if you have ample RAM to reduce "
        "spill I/O.",
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
        f"intergenic regions merged up to this many slots per block (default {qd['sweep_block_slots']}; 'none' is not "
        "accepted here — the whole chain is a config choice). "
        "Performance only — the answer is the same for every value; it sizes the solver's per-thread "
        "arena, so smaller blocks use less memory per sweep. Advanced "
        "calibration knob.",
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
        help="Minimum posterior threshold for candidate pruning, in [0,1) "
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
        default=defaults.scenario_seed,
        help=f"Random seed (default: {defaults.scenario_seed}, overridden by YAML)",
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
