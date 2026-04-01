#!/usr/bin/env python3
"""Comprehensive benchmark analysis for rigel simulation studies.

Assumes pre-computed outputs exist in a unified benchmark directory:

    <benchmark_dir>/
        manifest.json
        truth_abundances_nrna_none.tsv
        truth_abundances_nrna_rand.tsv
        sample_sheet.tsv
        rigel_index/
        salmon_index/
        runs/<condition>/
            sim_oracle.bam
            sim_R1.fq.gz, sim_R2.fq.gz
            rigel_star/          (or rigel_oracle/, rigel_<config>/, ...)
                quant.feather
                gene_quant.feather
                nrna_quant.feather
                summary.json
                loci.feather
                config.yaml
            salmon/
                quant.sf.gz
                quant.genes.sf.gz

Compares tool outputs against ground truth and produces detailed reports
at transcript, gene, and pool levels with breakdowns by condition
parameters (gDNA fraction, strand specificity, nRNA contamination).

Usage
-----
::

    conda activate rigel
    python scripts/benchmark_analysis.py \\
        --benchmark-dir /path/to/hulkrna_benchmarks \\
        -o results/benchmark_report

Output
------
::

    <outdir>/
        transcript_metrics.csv      # per-tool, per-condition transcript metrics
        gene_metrics.csv            # per-tool, per-condition gene metrics
        pool_summary.csv            # pool-level (mRNA/nRNA/gDNA) truth vs pred
        per_transcript_detail.parquet  # full per-transcript truth vs pred (all tools, all conditions)
        condition_summary.csv       # condition metadata + fragment counts
        report.md                   # human-readable summary
        figures/                    # (optional) if matplotlib available
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════

BENCHMARK_DIR_ENV = "RIGEL_BENCHMARK_DIR"
DEFAULT_BENCHMARK_DIR = (
    "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks"
)


@dataclass
class ToolOutput:
    """Descriptor for a tool's output location and format."""

    name: str
    subdir: str
    # transcript-level
    tx_file: str = ""
    tx_id_col: str = ""
    tx_count_col: str = ""
    tx_tpm_col: str = ""
    # gene-level
    gene_file: str = ""
    gene_id_col: str = ""
    gene_count_col: str = ""
    gene_tpm_col: str = ""
    # format
    file_format: str = "feather"  # "feather", "tsv", "tsv.gz"
    # summary / metadata
    summary_file: str = ""


# Built-in tool output definitions
RIGEL_TOOL = ToolOutput(
    name="rigel",
    subdir="rigel_star",
    tx_file="quant.feather",
    tx_id_col="transcript_id",
    tx_count_col="mrna",
    tx_tpm_col="tpm",
    gene_file="gene_quant.feather",
    gene_id_col="gene_id",
    gene_count_col="mrna",
    gene_tpm_col="tpm",
    file_format="feather",
    summary_file="summary.json",
)

SALMON_TOOL = ToolOutput(
    name="salmon",
    subdir="salmon",
    tx_file="quant.sf.gz",
    tx_id_col="Name",
    tx_count_col="NumReads",
    tx_tpm_col="TPM",
    gene_file="quant.genes.sf.gz",
    gene_id_col="Name",
    gene_count_col="NumReads",
    gene_tpm_col="TPM",
    file_format="tsv.gz",
)


def _parse_condition_name(name: str) -> dict:
    """Parse condition name → metadata dict.

    Format: gdna_{none|low|high}_ss_{0.50|0.90|1.00}_nrna_{none|rand}
    """
    parts = name.split("_")
    meta = {"condition": name}
    try:
        gdna_idx = parts.index("gdna") + 1
        meta["gdna_label"] = parts[gdna_idx]
    except (ValueError, IndexError):
        meta["gdna_label"] = "unknown"
    try:
        ss_idx = parts.index("ss") + 1
        meta["strand_specificity"] = float(parts[ss_idx])
    except (ValueError, IndexError):
        meta["strand_specificity"] = 1.0
    try:
        nrna_idx = parts.index("nrna") + 1
        meta["nrna_label"] = parts[nrna_idx]
    except (ValueError, IndexError):
        meta["nrna_label"] = "none"
    return meta


# ═══════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════


def load_truth(benchmark_dir: Path, nrna_label: str) -> pd.DataFrame:
    """Load truth abundances for a given nRNA condition."""
    fname = f"truth_abundances_nrna_{nrna_label}.tsv"
    path = benchmark_dir / fname
    if not path.exists():
        raise FileNotFoundError(f"Truth file not found: {path}")
    df = pd.read_csv(path, sep="\t")
    return df


def load_manifest(benchmark_dir: Path) -> dict:
    """Load manifest.json."""
    path = benchmark_dir / "manifest.json"
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with open(path) as f:
        return json.load(f)


def discover_conditions(benchmark_dir: Path) -> list[str]:
    """Discover available condition directories."""
    runs_dir = benchmark_dir / "runs"
    if not runs_dir.exists():
        raise FileNotFoundError(f"Runs directory not found: {runs_dir}")
    conditions = sorted([
        d.name for d in runs_dir.iterdir()
        if d.is_dir() and d.name.startswith("gdna_")
    ])
    return conditions


def discover_tools(benchmark_dir: Path, condition: str) -> list[ToolOutput]:
    """Discover available tool outputs for a condition."""
    cond_dir = benchmark_dir / "runs" / condition
    tools = []

    # Check for rigel outputs (any rigel_* subdirectory)
    for d in sorted(cond_dir.iterdir()):
        if d.is_dir() and d.name.startswith("rigel_"):
            variant = d.name  # e.g. rigel_star, rigel_oracle, rigel_custom
            tool = ToolOutput(
                name=variant,
                subdir=variant,
                tx_file="quant.feather",
                tx_id_col="transcript_id",
                tx_count_col="mrna",
                tx_tpm_col="tpm",
                gene_file="gene_quant.feather",
                gene_id_col="gene_id",
                gene_count_col="mrna",
                gene_tpm_col="tpm",
                file_format="feather",
                summary_file="summary.json",
            )
            # Verify the quant file exists
            if (d / tool.tx_file).exists():
                tools.append(tool)

    # Check for salmon
    salmon_dir = cond_dir / "salmon"
    if salmon_dir.exists() and (salmon_dir / "quant.sf.gz").exists():
        tools.append(SALMON_TOOL)

    # Check for kallisto (future)
    kallisto_dir = cond_dir / "kallisto"
    if kallisto_dir.exists():
        abundance = kallisto_dir / "abundance.tsv"
        if abundance.exists():
            tools.append(ToolOutput(
                name="kallisto",
                subdir="kallisto",
                tx_file="abundance.tsv",
                tx_id_col="target_id",
                tx_count_col="est_counts",
                tx_tpm_col="tpm",
                file_format="tsv",
            ))

    return tools


def load_tool_tx_output(
    benchmark_dir: Path, condition: str, tool: ToolOutput
) -> pd.DataFrame | None:
    """Load transcript-level output for a tool."""
    path = benchmark_dir / "runs" / condition / tool.subdir / tool.tx_file
    if not path.exists():
        return None
    if tool.file_format == "feather":
        return pd.read_feather(path)
    elif tool.file_format == "tsv.gz":
        return pd.read_csv(path, sep="\t")
    elif tool.file_format == "tsv":
        return pd.read_csv(path, sep="\t")
    return None


def load_tool_gene_output(
    benchmark_dir: Path, condition: str, tool: ToolOutput
) -> pd.DataFrame | None:
    """Load gene-level output for a tool."""
    if not tool.gene_file:
        return None
    path = benchmark_dir / "runs" / condition / tool.subdir / tool.gene_file
    if not path.exists():
        return None
    if tool.file_format == "feather":
        return pd.read_feather(path)
    elif tool.file_format == "tsv.gz":
        return pd.read_csv(path, sep="\t")
    elif tool.file_format == "tsv":
        return pd.read_csv(path, sep="\t")
    return None


def load_rigel_summary(
    benchmark_dir: Path, condition: str, tool: ToolOutput
) -> dict | None:
    """Load rigel summary.json for pool-level stats."""
    if not tool.summary_file:
        return None
    path = benchmark_dir / "runs" / condition / tool.subdir / tool.summary_file
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def parse_truth_from_fastq(r1_path: Path) -> tuple[dict[str, int], dict[str, int], int]:
    """Parse ground-truth fragment counts from FASTQ read names.

    Returns (mrna_counts, nrna_counts, n_gdna).
    """
    mrna_counts: Counter[str] = Counter()
    nrna_counts: Counter[str] = Counter()
    n_gdna = 0

    opener = gzip.open if str(r1_path).endswith(".gz") else open
    with opener(r1_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            qname = line[1:].strip()
            if qname.endswith("/1"):
                qname = qname[:-2]
            t_id = qname.split(":")[0]
            if t_id.startswith("gdna"):
                n_gdna += 1
            elif t_id.startswith("nrna_"):
                nrna_counts[t_id[5:]] += 1
            else:
                mrna_counts[t_id] += 1

    return dict(mrna_counts), dict(nrna_counts), n_gdna


# ═══════════════════════════════════════════════════════════════════
# Metrics computation
# ═══════════════════════════════════════════════════════════════════


def _spearman_r(x: np.ndarray, y: np.ndarray) -> float:
    """Spearman rank correlation (no scipy dependency)."""
    n = len(x)
    if n < 2:
        return np.nan
    rx = np.empty(n)
    ry = np.empty(n)
    for src, dst in ((x, rx), (y, ry)):
        order = np.argsort(src)
        i = 0
        while i < n:
            j = i + 1
            while j < n and src[order[j]] == src[order[i]]:
                j += 1
            avg_rank = (i + j - 1) / 2.0
            for k in range(i, j):
                dst[order[k]] = avg_rank
            i = j
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def compute_metrics(
    truth: np.ndarray,
    predicted: np.ndarray,
    *,
    pseudocount: float = 0.0,
) -> dict:
    """Compute comprehensive accuracy metrics between truth and predicted arrays.

    Both arrays must be aligned (same length, same transcript ordering).
    """
    n = len(truth)
    assert len(predicted) == n

    # Basic error
    residual = predicted - truth
    abs_err = np.abs(residual)

    total_truth = float(truth.sum())
    total_pred = float(predicted.sum())

    mae = float(abs_err.mean()) if n > 0 else np.nan
    rmse = float(np.sqrt(np.mean(residual**2))) if n > 0 else np.nan
    total_abs_error = float(abs_err.sum())
    median_abs_error = float(np.median(abs_err)) if n > 0 else np.nan

    # Relative error (for nonzero truth)
    mask_pos = truth > 0
    n_expressed_truth = int(mask_pos.sum())
    if n_expressed_truth > 0:
        rel_err = abs_err[mask_pos] / truth[mask_pos]
        mean_rel_error = float(rel_err.mean())
        median_rel_error = float(np.median(rel_err))
    else:
        mean_rel_error = np.nan
        median_rel_error = np.nan

    # MAPE (mean absolute percentage error) on expressed transcripts
    mape = mean_rel_error * 100 if not np.isnan(mean_rel_error) else np.nan

    # Correlation
    if n > 1 and truth.std() > 0 and predicted.std() > 0:
        pearson_r = float(np.corrcoef(truth, predicted)[0, 1])
        spearman_r = _spearman_r(truth, predicted)
    else:
        pearson_r = np.nan
        spearman_r = np.nan

    # Log-space correlation (for TPM-like quantities)
    log_pearson = np.nan
    log_spearman = np.nan
    if pseudocount > 0:
        lt = np.log2(truth + pseudocount)
        lp = np.log2(predicted + pseudocount)
        if lt.std() > 0 and lp.std() > 0:
            log_pearson = float(np.corrcoef(lt, lp)[0, 1])
            log_spearman = _spearman_r(lt, lp)

    # Classification-style metrics (expressed/not expressed)
    pred_pos = predicted > 0
    tp = int((mask_pos & pred_pos).sum())
    fp = int((~mask_pos & pred_pos).sum())
    fn = int((mask_pos & ~pred_pos).sum())
    tn = int((~mask_pos & ~pred_pos).sum())
    precision = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    recall = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else np.nan

    # Coefficient of variation of RMSE
    cv_rmse = rmse / total_truth * n if total_truth > 0 and not np.isnan(rmse) else np.nan

    # Weighted absolute relative error
    if total_truth > 0 and n_expressed_truth > 0:
        ware = float((abs_err[mask_pos] / truth[mask_pos] * truth[mask_pos]).sum() / total_truth)
    else:
        ware = np.nan

    return {
        "n_transcripts": n,
        "n_expressed_truth": n_expressed_truth,
        "n_expressed_pred": int(pred_pos.sum()),
        "total_truth": total_truth,
        "total_predicted": total_pred,
        "total_abs_error": total_abs_error,
        "mae": mae,
        "median_abs_error": median_abs_error,
        "rmse": rmse,
        "mean_rel_error": mean_rel_error,
        "median_rel_error": median_rel_error,
        "mape": mape,
        "ware": ware,
        "pearson_r": pearson_r,
        "spearman_r": spearman_r,
        "log2_pearson_r": log_pearson,
        "log2_spearman_r": log_spearman,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
    }


def compute_stratified_metrics(
    truth: np.ndarray,
    predicted: np.ndarray,
    *,
    expression_bins: list[tuple[str, float, float]] | None = None,
) -> dict[str, dict]:
    """Compute metrics stratified by expression level.

    Returns dict mapping bin_name → metrics dict.
    """
    if expression_bins is None:
        expression_bins = [
            ("zero", 0.0, 0.0),
            ("low_1-10", 1.0, 10.0),
            ("mid_10-100", 10.0, 100.0),
            ("high_100-1000", 100.0, 1000.0),
            ("very_high_1000+", 1000.0, np.inf),
        ]

    results = {}
    for bin_name, lo, hi in expression_bins:
        if lo == 0.0 and hi == 0.0:
            mask = truth == 0
        elif hi == np.inf:
            mask = truth >= lo
        else:
            mask = (truth >= lo) & (truth < hi)

        if mask.sum() == 0:
            continue

        results[bin_name] = compute_metrics(
            truth[mask], predicted[mask], pseudocount=1.0
        )

    return results


# ═══════════════════════════════════════════════════════════════════
# Pool-level analysis (mRNA / nRNA / gDNA)
# ═══════════════════════════════════════════════════════════════════


def compute_pool_metrics(
    benchmark_dir: Path,
    condition: str,
    tool: ToolOutput,
    truth_df: pd.DataFrame,
    manifest: dict,
    fastq_truth: tuple[dict[str, int], dict[str, int], int] | None = None,
) -> dict | None:
    """Compute pool-level (mRNA/nRNA/gDNA) metrics for rigel outputs."""
    summary = load_rigel_summary(benchmark_dir, condition, tool)
    if summary is None:
        return None

    meta = _parse_condition_name(condition)
    quant = summary.get("quantification", {})
    frag_stats = summary.get("fragment_stats", {})

    # From rigel summary
    mrna_pred = quant.get("mrna_total", 0.0)
    nrna_pred = quant.get("nrna_total", 0.0)
    gdna_pred = quant.get("gdna_total", 0.0)
    intergenic_pred = quant.get("intergenic_total", 0.0)

    # Ground truth from manifest
    sim_params = manifest.get("simulation", {})
    n_rna_frags = sim_params.get("n_rna_fragments", 0)

    # Find this condition in manifest conditions list
    cond_meta = {}
    for c in manifest.get("conditions", []):
        if c.get("name") == condition:
            cond_meta = c
            break

    n_rna_truth = cond_meta.get("n_rna", 0)
    n_gdna_truth = cond_meta.get("n_gdna", 0)

    # mRNA / nRNA truth from truth_df abundances, scaled to fragment counts
    mrna_truth_total = float(truth_df["mrna_abundance"].sum())
    nrna_truth_total = float(truth_df["nrna_abundance"].sum())
    total_rna_truth = mrna_truth_total + nrna_truth_total

    # Scale: truth abundances are fractional (sum to ~1M);
    # actual fragment counts = abundance * (n_rna / total_abundance)
    mrna_frag_truth = 0
    nrna_frag_truth = 0
    gdna_frag_truth = 0
    if fastq_truth is not None:
        mrna_frag_truth = sum(fastq_truth[0].values())
        nrna_frag_truth = sum(fastq_truth[1].values())
        gdna_frag_truth = fastq_truth[2]
    elif total_rna_truth > 0 and n_rna_truth > 0:
        scale = n_rna_truth / total_rna_truth
        mrna_frag_truth = int(round(mrna_truth_total * scale))
        nrna_frag_truth = int(round(nrna_truth_total * scale))
        gdna_frag_truth = n_gdna_truth
    else:
        mrna_frag_truth = n_rna_truth
        gdna_frag_truth = n_gdna_truth

    total_frags = frag_stats.get("total", 0)

    result = {
        "condition": condition,
        "tool": tool.name,
        "gdna_label": meta["gdna_label"],
        "strand_specificity": meta["strand_specificity"],
        "nrna_label": meta["nrna_label"],
        # Fragment-level truth
        "mrna_frag_truth": mrna_frag_truth,
        "nrna_frag_truth": nrna_frag_truth,
        "gdna_frag_truth": gdna_frag_truth,
        "total_frag_truth": mrna_frag_truth + nrna_frag_truth + gdna_frag_truth,
        # Rigel predictions
        "mrna_pred": mrna_pred,
        "nrna_pred": nrna_pred,
        "gdna_pred": gdna_pred,
        "intergenic_pred": intergenic_pred,
        "total_frag_assigned": total_frags,
        # Errors
        "mrna_error": mrna_pred - mrna_frag_truth if mrna_frag_truth > 0 else np.nan,
        "nrna_error": nrna_pred - nrna_frag_truth,
        "gdna_error": gdna_pred - gdna_frag_truth if gdna_frag_truth > 0 else np.nan,
        "mrna_rel_error": (
            (mrna_pred - mrna_frag_truth) / mrna_frag_truth
            if mrna_frag_truth > 0
            else np.nan
        ),
        "nrna_rel_error": (
            (nrna_pred - nrna_frag_truth) / nrna_frag_truth
            if nrna_frag_truth > 0
            else np.nan
        ),
        "gdna_rel_error": (
            (gdna_pred - gdna_frag_truth) / gdna_frag_truth
            if gdna_frag_truth > 0
            else np.nan
        ),
        # Calibration
        "strand_specificity_estimated": summary.get("strand_model", {}).get(
            "strand_specificity", np.nan
        ),
        "gdna_density_estimated": summary.get("calibration", {}).get(
            "gdna_density_global", np.nan
        ),
    }
    return result


# ═══════════════════════════════════════════════════════════════════
# Rigel-specific detailed analysis
# ═══════════════════════════════════════════════════════════════════


def rigel_detailed_analysis(
    benchmark_dir: Path,
    condition: str,
    tool: ToolOutput,
) -> dict | None:
    """Extract rigel-specific details: calibration, strand model, etc."""
    summary = load_rigel_summary(benchmark_dir, condition, tool)
    if summary is None:
        return None

    strand = summary.get("strand_model", {})
    calib = summary.get("calibration", {})
    frag_len = summary.get("fragment_length", {})
    quant = summary.get("quantification", {})

    return {
        "condition": condition,
        "tool": tool.name,
        # Strand model
        "strand_protocol": strand.get("protocol", ""),
        "strand_specificity": strand.get("strand_specificity", np.nan),
        "p_r1_sense": strand.get("p_r1_sense", np.nan),
        "strand_n_training": strand.get("n_training_fragments", 0),
        "strand_posterior_var": strand.get("posterior_variance", np.nan),
        # Calibration
        "gdna_density_global": calib.get("gdna_density_global", np.nan),
        "gdna_mixing_prop": calib.get("mixing_proportion", np.nan),
        "kappa_strand": calib.get("kappa_strand", np.nan),
        "gdna_fl_mean": calib.get("gdna_fl_mean", np.nan),
        "gdna_fl_mode": calib.get("gdna_fl_mode", np.nan),
        "n_eligible_regions": calib.get("n_eligible_regions", 0),
        "calibration_converged": calib.get("converged", False),
        # Fragment length
        "fl_global_mean": frag_len.get("global", {}).get("mean", np.nan),
        "fl_global_mode": frag_len.get("global", {}).get("mode", np.nan),
        "fl_rna_mean": frag_len.get("rna", {}).get("mean", np.nan),
        "fl_intergenic_mean": frag_len.get("intergenic", {}).get("mean", np.nan),
        # Quantification summary
        "n_transcripts": quant.get("n_transcripts", 0),
        "n_genes": quant.get("n_genes", 0),
        "n_loci": quant.get("n_loci", 0),
        "n_unambig": quant.get("n_unambig_assigned", 0),
        "n_em": quant.get("n_em_assigned", 0),
        "mrna_fraction": quant.get("mrna_fraction", np.nan),
    }


# ═══════════════════════════════════════════════════════════════════
# Main analysis pipeline
# ═══════════════════════════════════════════════════════════════════


def run_analysis(
    benchmark_dir: Path,
    outdir: Path,
    *,
    parse_fastq_truth: bool = False,
    log2_pseudocount: float = 1.0,
) -> None:
    """Run full benchmark analysis.

    Parameters
    ----------
    benchmark_dir : Path
        Root of the unified benchmark directory.
    outdir : Path
        Directory for output reports.
    parse_fastq_truth : bool
        If True, parse FASTQ read names for per-fragment ground truth
        (slow but precise). If False, use truth TSV only.
    log2_pseudocount : float
        Pseudocount for log-space correlation metrics.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    logger.info("Benchmark directory: %s", benchmark_dir)
    logger.info("Output directory: %s", outdir)

    # Load manifest
    manifest = load_manifest(benchmark_dir)
    sim_params = manifest.get("simulation", {})
    logger.info(
        "Simulation: %d RNA fragments, seed=%d",
        sim_params.get("n_rna_fragments", 0),
        sim_params.get("sim_seed", 0),
    )

    # Load truth tables
    truth_cache: dict[str, pd.DataFrame] = {}
    for nrna_label in ("none", "rand"):
        try:
            truth_cache[nrna_label] = load_truth(benchmark_dir, nrna_label)
            logger.info(
                "Loaded truth_%s: %d transcripts",
                nrna_label,
                len(truth_cache[nrna_label]),
            )
        except FileNotFoundError:
            logger.warning("Truth file for nrna_%s not found", nrna_label)

    # Discover conditions
    conditions = discover_conditions(benchmark_dir)
    logger.info("Found %d conditions: %s", len(conditions), conditions)

    # ── Collect results ──────────────────────────────────────────
    all_tx_metrics: list[dict] = []
    all_gene_metrics: list[dict] = []
    all_pool_metrics: list[dict] = []
    all_rigel_details: list[dict] = []
    all_per_tx_rows: list[pd.DataFrame] = []
    all_stratified: list[dict] = []
    condition_summaries: list[dict] = []

    # Precompute tx→gene mapping from truth (needed for salmon gene-level)
    tx2gene: dict[str, str] = {}
    for nrna_label, tdf in truth_cache.items():
        for row in tdf[["transcript_id", "gene_id"]].itertuples(index=False):
            tx2gene[row.transcript_id] = row.gene_id

    for ci, condition in enumerate(conditions):
        meta = _parse_condition_name(condition)
        nrna_label = meta["nrna_label"]
        truth_df = truth_cache.get(nrna_label)
        if truth_df is None:
            logger.warning("No truth for %s (nrna=%s), skipping", condition, nrna_label)
            continue

        print(
            f"[{ci + 1}/{len(conditions)}] {condition} "
            f"(gdna={meta['gdna_label']}, ss={meta['strand_specificity']}, "
            f"nrna={meta['nrna_label']})",
            flush=True,
        )

        # Parse FASTQ truth if requested
        fastq_truth = None
        if parse_fastq_truth:
            r1_path = benchmark_dir / "runs" / condition / "sim_R1.fq.gz"
            if r1_path.exists():
                logger.info("  Parsing FASTQ truth from %s...", r1_path)
                fastq_truth = parse_truth_from_fastq(r1_path)
                logger.info(
                    "  FASTQ truth: mRNA=%d, nRNA=%d, gDNA=%d",
                    sum(fastq_truth[0].values()),
                    sum(fastq_truth[1].values()),
                    fastq_truth[2],
                )

        # Condition metadata from manifest
        cond_manifest = {}
        for c in manifest.get("conditions", []):
            if c.get("name") == condition:
                cond_manifest = c
                break

        n_rna_total = cond_manifest.get("n_rna", 0)
        n_gdna_total = cond_manifest.get("n_gdna", 0)

        # Compute truth scaling: truth abundances are fractional (sum to 1M)
        # but simulation generates n_rna_fragments actual fragments.
        # Scale truth to fragment counts.
        truth_mrna_sum = float(truth_df["mrna_abundance"].sum())
        truth_nrna_sum = float(truth_df["nrna_abundance"].sum())
        truth_total_sum = truth_mrna_sum + truth_nrna_sum
        if truth_total_sum > 0 and n_rna_total > 0:
            truth_scale = n_rna_total / truth_total_sum
        else:
            truth_scale = 1.0

        # Scaled fragment counts
        mrna_frag_truth = int(round(truth_mrna_sum * truth_scale))
        nrna_frag_truth = int(round(truth_nrna_sum * truth_scale))

        cond_summary = {
            **meta,
            "n_rna": n_rna_total,
            "n_gdna": n_gdna_total,
            "gdna_rate": cond_manifest.get("gdna_rate", 0.0),
            "mrna_frag_truth": mrna_frag_truth,
            "nrna_frag_truth": nrna_frag_truth,
            "truth_scale": truth_scale,
        }
        if fastq_truth is not None:
            cond_summary["fastq_mrna"] = sum(fastq_truth[0].values())
            cond_summary["fastq_nrna"] = sum(fastq_truth[1].values())
            cond_summary["fastq_gdna"] = fastq_truth[2]
        condition_summaries.append(cond_summary)

        # Discover tools
        tools = discover_tools(benchmark_dir, condition)
        logger.info("  Tools: %s", [t.name for t in tools])

        for tool in tools:
            print(f"  {tool.name}...", end="", flush=True)

            # ── Transcript-level ─────────────────────────────────
            tx_df = load_tool_tx_output(benchmark_dir, condition, tool)
            if tx_df is not None:
                # Build aligned arrays: truth vs predicted
                # For rigel: merge on transcript_id
                # For salmon: merge on Name
                tx_id_col = tool.tx_id_col
                tx_count_col = tool.tx_count_col

                # Merge truth with predictions
                merged = truth_df.merge(
                    tx_df[[tx_id_col, tx_count_col]].rename(
                        columns={tx_id_col: "transcript_id", tx_count_col: "predicted"}
                    ),
                    on="transcript_id",
                    how="left",
                )
                merged["predicted"] = merged["predicted"].fillna(0.0)

                # Scale truth abundances to fragment counts
                # Truth mrna_abundance values are fractional (sum ~ 753K-1M);
                # actual fragment counts = truth * truth_scale
                truth_arr = merged["mrna_abundance"].values.astype(np.float64) * truth_scale
                pred_arr = merged["predicted"].values.astype(np.float64)

                metrics = compute_metrics(truth_arr, pred_arr, pseudocount=log2_pseudocount)
                metrics["condition"] = condition
                metrics["tool"] = tool.name
                metrics.update(meta)
                all_tx_metrics.append(metrics)

                # Stratified metrics
                strat = compute_stratified_metrics(truth_arr, pred_arr)
                for bin_name, bin_metrics in strat.items():
                    bin_metrics["condition"] = condition
                    bin_metrics["tool"] = tool.name
                    bin_metrics["expression_bin"] = bin_name
                    bin_metrics.update(meta)
                    all_stratified.append(bin_metrics)

                # Per-transcript detail (keep concise: only expressed or predicted)
                merged["truth_scaled"] = merged["mrna_abundance"] * truth_scale
                detail_mask = (merged["truth_scaled"] > 0) | (merged["predicted"] > 0)
                if detail_mask.any():
                    detail = merged.loc[detail_mask, [
                        "transcript_id", "gene_id", "gene_name",
                        "mrna_abundance", "nrna_abundance",
                        "truth_scaled", "predicted",
                    ]].copy()
                    detail["condition"] = condition
                    detail["tool"] = tool.name
                    detail["residual"] = detail["predicted"] - detail["truth_scaled"]
                    detail["abs_error"] = detail["residual"].abs()
                    detail["rel_error"] = np.where(
                        detail["truth_scaled"] > 0,
                        detail["abs_error"] / detail["truth_scaled"],
                        np.nan,
                    )
                    all_per_tx_rows.append(detail)

            # ── Gene-level ───────────────────────────────────────
            if "gene_id" in truth_df.columns:
                # Aggregate truth to gene level (scaled)
                gene_truth = (
                    truth_df.groupby("gene_id")["mrna_abundance"]
                    .sum()
                    .reset_index()
                )
                gene_truth["mrna_scaled"] = gene_truth["mrna_abundance"] * truth_scale

                # For rigel: use gene_quant.feather directly
                # For salmon: quant.genes.sf.gz has transcript IDs as gene names
                #   (broken tx2gene mapping), so aggregate tx-level counts instead
                gene_pred: dict[str, float] | None = None

                if tool.name.startswith("rigel"):
                    gene_df = load_tool_gene_output(benchmark_dir, condition, tool)
                    if gene_df is not None:
                        gene_pred = dict(zip(
                            gene_df[tool.gene_id_col],
                            gene_df[tool.gene_count_col],
                        ))
                elif tx_df is not None:
                    # Aggregate transcript-level counts to genes via tx2gene
                    gene_pred = {}
                    for _, row in tx_df.iterrows():
                        tid = row[tool.tx_id_col]
                        gid = tx2gene.get(tid)
                        if gid is not None:
                            gene_pred[gid] = gene_pred.get(gid, 0.0) + float(row[tool.tx_count_col])

                if gene_pred is not None:
                    gene_truth["predicted"] = gene_truth["gene_id"].map(
                        lambda g: gene_pred.get(g, 0.0)
                    )
                    g_truth = gene_truth["mrna_scaled"].values.astype(np.float64)
                    g_pred = gene_truth["predicted"].values.astype(np.float64)

                    g_metrics = compute_metrics(g_truth, g_pred, pseudocount=log2_pseudocount)
                    g_metrics["condition"] = condition
                    g_metrics["tool"] = tool.name
                    g_metrics.update(meta)
                    all_gene_metrics.append(g_metrics)

            # ── Pool-level (rigel only) ──────────────────────────
            if tool.name.startswith("rigel"):
                pool = compute_pool_metrics(
                    benchmark_dir, condition, tool, truth_df, manifest, fastq_truth
                )
                if pool:
                    all_pool_metrics.append(pool)

                detail = rigel_detailed_analysis(benchmark_dir, condition, tool)
                if detail:
                    all_rigel_details.append(detail)

            print(" done", flush=True)

    # ── Write outputs ────────────────────────────────────────────
    print("\nWriting outputs...", flush=True)

    tx_metrics_df = pd.DataFrame(all_tx_metrics)
    tx_metrics_df.to_csv(outdir / "transcript_metrics.csv", index=False)
    logger.info("Wrote transcript_metrics.csv (%d rows)", len(tx_metrics_df))

    gene_metrics_df = pd.DataFrame(all_gene_metrics)
    gene_metrics_df.to_csv(outdir / "gene_metrics.csv", index=False)
    logger.info("Wrote gene_metrics.csv (%d rows)", len(gene_metrics_df))

    if all_pool_metrics:
        pool_df = pd.DataFrame(all_pool_metrics)
        pool_df.to_csv(outdir / "pool_summary.csv", index=False)
        logger.info("Wrote pool_summary.csv (%d rows)", len(pool_df))

    if all_rigel_details:
        detail_df = pd.DataFrame(all_rigel_details)
        detail_df.to_csv(outdir / "rigel_details.csv", index=False)
        logger.info("Wrote rigel_details.csv (%d rows)", len(detail_df))

    if all_stratified:
        strat_df = pd.DataFrame(all_stratified)
        strat_df.to_csv(outdir / "stratified_metrics.csv", index=False)
        logger.info("Wrote stratified_metrics.csv (%d rows)", len(strat_df))

    if all_per_tx_rows:
        per_tx_df = pd.concat(all_per_tx_rows, ignore_index=True)
        per_tx_df.to_parquet(outdir / "per_transcript_detail.parquet", index=False)
        logger.info(
            "Wrote per_transcript_detail.parquet (%d rows)", len(per_tx_df)
        )

    cond_df = pd.DataFrame(condition_summaries)
    cond_df.to_csv(outdir / "condition_summary.csv", index=False)

    # ── Generate markdown report ─────────────────────────────────
    write_report(
        outdir / "report.md",
        tx_metrics_df,
        gene_metrics_df,
        pd.DataFrame(all_pool_metrics) if all_pool_metrics else None,
        pd.DataFrame(all_stratified) if all_stratified else None,
        pd.DataFrame(all_rigel_details) if all_rigel_details else None,
        manifest,
    )

    print(f"Analysis complete. Results in: {outdir}", flush=True)


# ═══════════════════════════════════════════════════════════════════
# Report generation
# ═══════════════════════════════════════════════════════════════════


def write_report(
    path: Path,
    tx_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    pool_df: pd.DataFrame | None,
    strat_df: pd.DataFrame | None,
    rigel_df: pd.DataFrame | None,
    manifest: dict,
) -> None:
    """Write comprehensive markdown report."""
    with open(path, "w") as f:
        f.write("# Benchmark Analysis Report\n\n")
        f.write(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        # Simulation info
        sim = manifest.get("simulation", {})
        f.write("## Simulation Parameters\n\n")
        f.write(f"- RNA fragments: {sim.get('n_rna_fragments', 'N/A'):,}\n")
        f.write(f"- Fragment length: {sim.get('frag_mean', 'N/A')} ± {sim.get('frag_std', 'N/A')}\n")
        f.write(f"- Read length: {sim.get('read_length', 'N/A')}\n")
        f.write(f"- Error rate: {sim.get('error_rate', 'N/A')}\n")
        f.write(f"- Seed: {sim.get('sim_seed', 'N/A')}\n\n")

        # Overview by tool
        if not tx_df.empty:
            f.write("## Transcript-Level Summary by Tool\n\n")
            tools = tx_df["tool"].unique()
            for tool in sorted(tools):
                tdf = tx_df[tx_df["tool"] == tool]
                f.write(f"### {tool}\n\n")
                f.write(f"- Conditions evaluated: {len(tdf)}\n")
                f.write(f"- Mean Pearson R: {tdf['pearson_r'].mean():.4f} "
                        f"(range: {tdf['pearson_r'].min():.4f}–{tdf['pearson_r'].max():.4f})\n")
                f.write(f"- Mean Spearman R: {tdf['spearman_r'].mean():.4f}\n")
                f.write(f"- Mean MAPE: {tdf['mape'].mean():.1f}%\n")
                f.write(f"- Mean RMSE: {tdf['rmse'].mean():.2f}\n")
                f.write(f"- Mean MAE: {tdf['mae'].mean():.2f}\n")
                f.write(f"- Mean WARE: {tdf['ware'].mean():.4f}\n")
                f.write("\n")

        # Transcript-level table: all conditions × tools
        if not tx_df.empty:
            f.write("## Transcript-Level Metrics (All Conditions)\n\n")
            cols = [
                "condition", "tool", "gdna_label", "strand_specificity",
                "nrna_label", "pearson_r", "spearman_r", "rmse", "mae",
                "mape", "ware", "precision", "recall", "f1",
            ]
            avail = [c for c in cols if c in tx_df.columns]
            _write_md_table(f, tx_df[avail])

        # Gene-level table
        if not gene_df.empty:
            f.write("\n## Gene-Level Metrics (All Conditions)\n\n")
            cols = [
                "condition", "tool", "gdna_label", "strand_specificity",
                "nrna_label", "pearson_r", "spearman_r", "rmse", "mae",
                "mape", "ware",
            ]
            avail = [c for c in cols if c in gene_df.columns]
            _write_md_table(f, gene_df[avail])

        # Pool-level
        if pool_df is not None and not pool_df.empty:
            f.write("\n## Pool-Level Summary (Rigel)\n\n")
            cols = [
                "condition", "tool",
                "mrna_frag_truth", "mrna_pred", "mrna_rel_error",
                "nrna_frag_truth", "nrna_pred", "nrna_rel_error",
                "gdna_frag_truth", "gdna_pred", "gdna_rel_error",
            ]
            avail = [c for c in cols if c in pool_df.columns]
            _write_md_table(f, pool_df[avail])

        # Stratified metrics
        if strat_df is not None and not strat_df.empty:
            f.write("\n## Stratified Metrics by Expression Level\n\n")
            # Pivot: for each tool, show metrics by expression bin across conditions
            for tool in sorted(strat_df["tool"].unique()):
                f.write(f"### {tool}\n\n")
                tdf = strat_df[strat_df["tool"] == tool]
                pivot = tdf.groupby("expression_bin").agg(
                    mean_pearson=("pearson_r", "mean"),
                    mean_spearman=("spearman_r", "mean"),
                    mean_mape=("mape", "mean"),
                    mean_mae=("mae", "mean"),
                    mean_ware=("ware", "mean"),
                    n_conditions=("condition", "nunique"),
                ).reset_index()
                _write_md_table(f, pivot)
                f.write("\n")

        # Rigel calibration details
        if rigel_df is not None and not rigel_df.empty:
            f.write("\n## Rigel Calibration & Model Details\n\n")
            cols = [
                "condition", "tool",
                "strand_specificity", "strand_protocol",
                "gdna_density_global", "gdna_mixing_prop",
                "kappa_strand", "calibration_converged",
                "fl_global_mean", "fl_rna_mean", "fl_intergenic_mean",
                "n_loci", "n_unambig", "n_em", "mrna_fraction",
            ]
            avail = [c for c in cols if c in rigel_df.columns]
            _write_md_table(f, rigel_df[avail])

        # Cross-tool comparison
        if not tx_df.empty and len(tx_df["tool"].unique()) > 1:
            f.write("\n## Cross-Tool Comparison\n\n")
            f.write("### Transcript-Level Pearson R by Condition\n\n")
            pivot = tx_df.pivot_table(
                index="condition", columns="tool", values="pearson_r"
            ).reset_index()
            _write_md_table(f, pivot)

            f.write("\n### Transcript-Level MAPE by Condition\n\n")
            pivot = tx_df.pivot_table(
                index="condition", columns="tool", values="mape"
            ).reset_index()
            _write_md_table(f, pivot)

            f.write("\n### Transcript-Level WARE by Condition\n\n")
            pivot = tx_df.pivot_table(
                index="condition", columns="tool", values="ware"
            ).reset_index()
            _write_md_table(f, pivot)

        # Sensitivity analysis: metrics vs condition params
        if not tx_df.empty:
            f.write("\n## Sensitivity Analysis\n\n")

            for tool in sorted(tx_df["tool"].unique()):
                tdf = tx_df[tx_df["tool"] == tool]
                f.write(f"### {tool}\n\n")

                # By gDNA level
                if "gdna_label" in tdf.columns:
                    f.write("#### By gDNA Contamination Level\n\n")
                    gdna_agg = tdf.groupby("gdna_label").agg(
                        mean_pearson=("pearson_r", "mean"),
                        mean_mape=("mape", "mean"),
                        mean_ware=("ware", "mean"),
                        mean_rmse=("rmse", "mean"),
                        n=("condition", "count"),
                    ).reset_index()
                    _write_md_table(f, gdna_agg)
                    f.write("\n")

                # By strand specificity
                if "strand_specificity" in tdf.columns:
                    f.write("#### By Strand Specificity\n\n")
                    ss_agg = tdf.groupby("strand_specificity").agg(
                        mean_pearson=("pearson_r", "mean"),
                        mean_mape=("mape", "mean"),
                        mean_ware=("ware", "mean"),
                        mean_rmse=("rmse", "mean"),
                        n=("condition", "count"),
                    ).reset_index()
                    _write_md_table(f, ss_agg)
                    f.write("\n")

                # By nRNA contamination
                if "nrna_label" in tdf.columns:
                    f.write("#### By nRNA Contamination\n\n")
                    nrna_agg = tdf.groupby("nrna_label").agg(
                        mean_pearson=("pearson_r", "mean"),
                        mean_mape=("mape", "mean"),
                        mean_ware=("ware", "mean"),
                        mean_rmse=("rmse", "mean"),
                        n=("condition", "count"),
                    ).reset_index()
                    _write_md_table(f, nrna_agg)
                    f.write("\n")

    logger.info("Wrote report: %s", path)


def _write_md_table(f, df: pd.DataFrame) -> None:
    """Write a pandas DataFrame as a markdown table."""
    if df.empty:
        f.write("(no data)\n")
        return

    # Format numeric columns
    formatted = df.copy()
    for col in formatted.columns:
        if formatted[col].dtype in (np.float64, np.float32, float):
            formatted[col] = formatted[col].map(
                lambda x: f"{x:.4f}" if pd.notna(x) and abs(x) < 100
                else (f"{x:.1f}" if pd.notna(x) else "—")
            )

    headers = list(formatted.columns)
    f.write("| " + " | ".join(str(h) for h in headers) + " |\n")
    f.write("| " + " | ".join("---" for _ in headers) + " |\n")
    for _, row in formatted.iterrows():
        f.write("| " + " | ".join(str(v) for v in row) + " |\n")


# ═══════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Analyze benchmark results from the unified benchmark directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--benchmark-dir", "-b",
        default=DEFAULT_BENCHMARK_DIR,
        help="Root benchmark directory (default: %(default)s)",
    )
    p.add_argument(
        "--output", "-o",
        default="benchmark_report",
        help="Output directory for reports (default: %(default)s)",
    )
    p.add_argument(
        "--parse-fastq-truth",
        action="store_true",
        help="Parse FASTQ read names for exact per-fragment ground truth "
             "(slow; ~10 min for 10M reads per condition)",
    )
    p.add_argument(
        "--log2-pseudocount",
        type=float,
        default=1.0,
        help="Pseudocount for log2-space correlation metrics (default: 1.0)",
    )
    p.add_argument(
        "--conditions",
        nargs="*",
        help="Subset of conditions to analyze (default: all)",
    )
    p.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose logging",
    )
    return p


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%H:%M:%S",
    )

    benchmark_dir = Path(args.benchmark_dir)
    if not benchmark_dir.exists():
        print(f"Error: benchmark directory not found: {benchmark_dir}", file=sys.stderr)
        return 1

    outdir = Path(args.output)

    print(f"Benchmark analysis", flush=True)
    print(f"  Benchmark dir: {benchmark_dir}", flush=True)
    print(f"  Output:        {outdir}", flush=True)
    print(f"  FASTQ truth:   {args.parse_fastq_truth}", flush=True)

    t0 = time.monotonic()
    run_analysis(
        benchmark_dir,
        outdir,
        parse_fastq_truth=args.parse_fastq_truth,
        log2_pseudocount=args.log2_pseudocount,
    )
    elapsed = time.monotonic() - t0
    print(f"\nCompleted in {elapsed:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
