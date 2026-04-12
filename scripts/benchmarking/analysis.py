"""Analysis: compare tool outputs against ground truth.

Refactored from ``scripts/benchmark/benchmark_analysis.py``.
"""
from __future__ import annotations

import gzip
import json
import logging
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import BenchmarkConfig

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════
# Tool output descriptors
# ═══════════════════════════════════════════════════════════════════


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
    file_format: str = "feather"  # "feather" | "tsv" | "tsv.gz"
    # summary / metadata
    summary_file: str = ""


def _rigel_tool(subdir: str, name: str | None = None) -> ToolOutput:
    """Create a ToolOutput descriptor for a rigel output directory."""
    return ToolOutput(
        name=name or subdir,
        subdir=subdir,
        tx_file="quant.feather",
        tx_id_col="transcript_id",
        tx_count_col="count",
        tx_tpm_col="tpm",
        gene_file="gene_quant.feather",
        gene_id_col="gene_id",
        gene_count_col="count",
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


# ═══════════════════════════════════════════════════════════════════
# Condition parsing
# ═══════════════════════════════════════════════════════════════════


def parse_condition(name: str) -> dict:
    """Parse condition name → metadata dict.

    Format: ``gdna_{none|low|high}_ss_{0.50|0.90|1.00}_nrna_{none|rand}``
    """
    parts = name.split("_")
    meta = {"condition": name}
    try:
        meta["gdna_label"] = parts[parts.index("gdna") + 1]
    except (ValueError, IndexError):
        meta["gdna_label"] = "unknown"
    try:
        meta["strand_specificity"] = float(parts[parts.index("ss") + 1])
    except (ValueError, IndexError):
        meta["strand_specificity"] = 1.0
    try:
        meta["nrna_label"] = parts[parts.index("nrna") + 1]
    except (ValueError, IndexError):
        meta["nrna_label"] = "none"
    return meta


# ═══════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════


def load_truth(benchmark_dir: Path, nrna_label: str) -> pd.DataFrame:
    """Load truth abundances for a given nRNA condition."""
    path = benchmark_dir / f"truth_abundances_nrna_{nrna_label}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Truth file not found: {path}")
    return pd.read_csv(path, sep="\t")


def load_truth_for_condition(
    benchmark_dir: Path, condition: str, nrna_label: str,
) -> pd.DataFrame:
    """Load truth: prefer per-condition file, fall back to nrna-label file."""
    # Per-condition truth (e.g. truth_abundances_pristine.tsv)
    per_cond = benchmark_dir / f"truth_abundances_{condition}.tsv"
    if per_cond.exists():
        return pd.read_csv(per_cond, sep="\t")
    # Legacy format (truth_abundances_nrna_none.tsv)
    return load_truth(benchmark_dir, nrna_label)


def discover_tools(cond_dir: Path) -> list[ToolOutput]:
    """Discover available tool outputs for a condition directory."""
    tools = []

    # rigel/<config>/ directories (new layout)
    rigel_dir = cond_dir / "rigel"
    if rigel_dir.exists() and rigel_dir.is_dir():
        for d in sorted(rigel_dir.iterdir()):
            if d.is_dir() and (d / "quant.feather").exists():
                tools.append(_rigel_tool(f"rigel/{d.name}", name=f"rigel/{d.name}"))

    # rigel_* directories (legacy layout — rigel_star, rigel_oracle, etc.)
    for d in sorted(cond_dir.iterdir()):
        if d.is_dir() and d.name.startswith("rigel_"):
            if (d / "quant.feather").exists():
                tools.append(_rigel_tool(d.name))

    # salmon
    salmon_dir = cond_dir / "salmon"
    if salmon_dir.exists() and (salmon_dir / "quant.sf.gz").exists():
        tools.append(SALMON_TOOL)

    # kallisto
    kallisto_dir = cond_dir / "kallisto"
    if kallisto_dir.exists() and (kallisto_dir / "abundance.tsv").exists():
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


def _read_table(path: Path, fmt: str) -> pd.DataFrame | None:
    """Read a table file, returning None if not found."""
    if not path.exists():
        return None
    if fmt == "feather":
        return pd.read_feather(path)
    elif fmt in ("tsv", "tsv.gz"):
        return pd.read_csv(path, sep="\t")
    return None


def load_tool_tx(cond_dir: Path, tool: ToolOutput) -> pd.DataFrame | None:
    return _read_table(cond_dir / tool.subdir / tool.tx_file, tool.file_format)


def load_tool_gene(cond_dir: Path, tool: ToolOutput) -> pd.DataFrame | None:
    if not tool.gene_file:
        return None
    return _read_table(cond_dir / tool.subdir / tool.gene_file, tool.file_format)


def load_summary(cond_dir: Path, tool: ToolOutput) -> dict | None:
    if not tool.summary_file:
        return None
    path = cond_dir / tool.subdir / tool.summary_file
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def parse_truth_from_fastq(
    r1_path: Path,
) -> tuple[dict[str, int], dict[str, int], int]:
    """Parse ground-truth fragment counts from FASTQ read names.

    Returns ``(mrna_counts, nrna_counts, n_gdna)``.
    """
    mrna: Counter[str] = Counter()
    nrna: Counter[str] = Counter()
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
                nrna[t_id[5:]] += 1
            else:
                mrna[t_id] += 1
    return dict(mrna), dict(nrna), n_gdna


# ═══════════════════════════════════════════════════════════════════
# Metrics
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
    """Compute accuracy metrics between aligned truth and predicted arrays."""
    n = len(truth)
    assert len(predicted) == n

    residual = predicted - truth
    abs_err = np.abs(residual)
    total_truth = float(truth.sum())
    total_pred = float(predicted.sum())
    mae = float(abs_err.mean()) if n > 0 else np.nan
    rmse = float(np.sqrt(np.mean(residual**2))) if n > 0 else np.nan
    total_abs_error = float(abs_err.sum())
    median_abs_error = float(np.median(abs_err)) if n > 0 else np.nan

    mask_pos = truth > 0
    n_expressed_truth = int(mask_pos.sum())
    if n_expressed_truth > 0:
        rel_err = abs_err[mask_pos] / truth[mask_pos]
        mean_rel_error = float(rel_err.mean())
        median_rel_error = float(np.median(rel_err))
    else:
        mean_rel_error = np.nan
        median_rel_error = np.nan

    mape = mean_rel_error * 100 if not np.isnan(mean_rel_error) else np.nan

    if n > 1 and truth.std() > 0 and predicted.std() > 0:
        pearson_r = float(np.corrcoef(truth, predicted)[0, 1])
        spearman_r = _spearman_r(truth, predicted)
    else:
        pearson_r = np.nan
        spearman_r = np.nan

    log_pearson = np.nan
    log_spearman = np.nan
    if pseudocount > 0:
        lt = np.log2(truth + pseudocount)
        lp = np.log2(predicted + pseudocount)
        if lt.std() > 0 and lp.std() > 0:
            log_pearson = float(np.corrcoef(lt, lp)[0, 1])
            log_spearman = _spearman_r(lt, lp)

    pred_pos = predicted > 0
    tp = int((mask_pos & pred_pos).sum())
    fp = int((~mask_pos & pred_pos).sum())
    fn = int((mask_pos & ~pred_pos).sum())
    tn = int((~mask_pos & ~pred_pos).sum())
    precision = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    recall = tp / (tp + fn) if (tp + fn) > 0 else np.nan
    f1 = (
        2 * precision * recall / (precision + recall)
        if (precision + recall) > 0
        else np.nan
    )

    if total_truth > 0 and n_expressed_truth > 0:
        ware = float(
            (abs_err[mask_pos] / truth[mask_pos] * truth[mask_pos]).sum()
            / total_truth
        )
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
    """Compute metrics stratified by expression level."""
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
        results[bin_name] = compute_metrics(truth[mask], predicted[mask], pseudocount=1.0)
    return results


# ═══════════════════════════════════════════════════════════════════
# Pool-level analysis (mRNA / nRNA / gDNA)
# ═══════════════════════════════════════════════════════════════════


def compute_pool_metrics(
    cond_dir: Path,
    condition: str,
    tool: ToolOutput,
    truth_df: pd.DataFrame,
    manifest: dict,
    fastq_truth: tuple[dict[str, int], dict[str, int], int] | None = None,
) -> dict | None:
    """Compute pool-level (mRNA/nRNA/gDNA) metrics for rigel outputs."""
    summary = load_summary(cond_dir, tool)
    if summary is None:
        return None

    meta = parse_condition(condition)
    quant = summary.get("quantification", {})

    mrna_pred = quant.get("mrna_total", 0.0)
    nrna_pred = quant.get("nrna_total", 0.0)
    gdna_em_pred = quant.get("gdna_total", 0.0)
    intergenic_pred = quant.get("intergenic_total", 0.0)
    # Intergenic fragments are gDNA that don't overlap any transcript
    # annotation and bypass the EM.  Combine with EM-assigned gDNA for
    # the total gDNA prediction.
    gdna_pred = gdna_em_pred + intergenic_pred

    # Find condition metadata in manifest
    cond_meta = {}
    for c in manifest.get("conditions", []):
        if c.get("name") == condition:
            cond_meta = c
            break

    n_rna_truth = cond_meta.get("n_rna", 0)
    n_gdna_truth = cond_meta.get("n_gdna", 0)

    mrna_truth_total = float(truth_df["mrna_abundance"].sum())
    nrna_truth_total = float(truth_df["nrna_abundance"].sum())
    total_rna_truth = mrna_truth_total + nrna_truth_total

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
        nrna_frag_truth = 0
        gdna_frag_truth = n_gdna_truth

    return {
        "condition": condition,
        "tool": tool.name,
        "gdna_label": meta["gdna_label"],
        "strand_specificity": meta["strand_specificity"],
        "nrna_label": meta["nrna_label"],
        # Truth
        "mrna_frag_truth": mrna_frag_truth,
        "nrna_frag_truth": nrna_frag_truth,
        "gdna_frag_truth": gdna_frag_truth,
        "total_frag_truth": mrna_frag_truth + nrna_frag_truth + gdna_frag_truth,
        # Predictions
        "mrna_pred": mrna_pred,
        "nrna_pred": nrna_pred,
        "gdna_pred": gdna_pred,
        "gdna_em_pred": gdna_em_pred,
        "intergenic_pred": intergenic_pred,
        # Errors
        "mrna_error": mrna_pred - mrna_frag_truth if mrna_frag_truth > 0 else np.nan,
        "nrna_error": nrna_pred - nrna_frag_truth,
        "gdna_error": gdna_pred - gdna_frag_truth if gdna_frag_truth > 0 else np.nan,
        "mrna_rel_error": (
            (mrna_pred - mrna_frag_truth) / mrna_frag_truth
            if mrna_frag_truth > 0 else np.nan
        ),
        "nrna_rel_error": (
            (nrna_pred - nrna_frag_truth) / nrna_frag_truth
            if nrna_frag_truth > 0 else np.nan
        ),
        "gdna_rel_error": (
            (gdna_pred - gdna_frag_truth) / gdna_frag_truth
            if gdna_frag_truth > 0 else np.nan
        ),
        # Calibration
        "strand_specificity_estimated": summary.get("strand_model", {}).get(
            "strand_specificity", np.nan
        ),
        "gdna_density_estimated": summary.get("calibration", {}).get(
            "gdna_density_global", np.nan
        ),
    }


def rigel_detailed_analysis(
    cond_dir: Path,
    condition: str,
    tool: ToolOutput,
) -> dict | None:
    """Extract rigel-specific details: calibration, strand model, etc."""
    summary = load_summary(cond_dir, tool)
    if summary is None:
        return None

    strand = summary.get("strand_model", {})
    calib = summary.get("calibration", {})
    frag_len = summary.get("fragment_length", {})
    quant = summary.get("quantification", {})

    return {
        "condition": condition,
        "tool": tool.name,
        "strand_protocol": strand.get("protocol", ""),
        "strand_specificity": strand.get("strand_specificity", np.nan),
        "p_r1_sense": strand.get("p_r1_sense", np.nan),
        "strand_n_training": strand.get("n_training_fragments", 0),
        "strand_posterior_var": strand.get("posterior_variance", np.nan),
        "gdna_density_global": calib.get("gdna_density_global", np.nan),
        "gdna_mixing_prop": calib.get("mixing_proportion", np.nan),
        "kappa_strand": calib.get("kappa_strand", np.nan),
        "gdna_fl_mean": calib.get("gdna_fl_mean", np.nan),
        "gdna_fl_mode": calib.get("gdna_fl_mode", np.nan),
        "n_eligible_regions": calib.get("n_eligible_regions", 0),
        "calibration_converged": calib.get("converged", False),
        "fl_global_mean": frag_len.get("global", {}).get("mean", np.nan),
        "fl_global_mode": frag_len.get("global", {}).get("mode", np.nan),
        "fl_rna_mean": frag_len.get("rna", {}).get("mean", np.nan),
        "fl_intergenic_mean": frag_len.get("intergenic", {}).get("mean", np.nan),
        "n_transcripts": quant.get("n_transcripts", 0),
        "n_genes": quant.get("n_genes", 0),
        "n_loci": quant.get("n_loci", 0),
        "n_unambig": quant.get("n_unambig_assigned", 0),
        "n_em": quant.get("n_em_assigned", 0),
        "mrna_fraction": quant.get("mrna_fraction", np.nan),
    }


# ═══════════════════════════════════════════════════════════════════
# Report generation
# ═══════════════════════════════════════════════════════════════════


def _write_md_table(f, df: pd.DataFrame) -> None:
    """Write a DataFrame as a markdown table."""
    if df.empty:
        f.write("(no data)\n")
        return
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

        sim = manifest.get("simulation", {})
        f.write("## Simulation Parameters\n\n")
        f.write(f"- RNA fragments: {sim.get('n_rna_fragments', 'N/A'):,}\n")
        f.write(f"- Fragment length: {sim.get('frag_mean', 'N/A')} "
                f"± {sim.get('frag_std', 'N/A')}\n")
        f.write(f"- Read length: {sim.get('read_length', 'N/A')}\n")
        f.write(f"- Error rate: {sim.get('error_rate', 'N/A')}\n")
        f.write(f"- Seed: {sim.get('sim_seed', 'N/A')}\n\n")

        # Tool summaries
        if not tx_df.empty:
            f.write("## Transcript-Level Summary by Tool\n\n")
            for tool in sorted(tx_df["tool"].unique()):
                tdf = tx_df[tx_df["tool"] == tool]
                f.write(f"### {tool}\n\n")
                f.write(f"- Conditions evaluated: {len(tdf)}\n")
                f.write(f"- Mean Spearman R: {tdf['spearman_r'].mean():.4f}\n")
                f.write(f"- Mean Pearson R: {tdf['pearson_r'].mean():.4f} "
                        f"(range: {tdf['pearson_r'].min():.4f}"
                        f"–{tdf['pearson_r'].max():.4f})\n")
                if "log2_spearman_r" in tdf.columns and tdf["log2_spearman_r"].notna().any():
                    f.write(f"- Mean log₂ Spearman R: "
                            f"{tdf['log2_spearman_r'].mean():.4f}\n")
                    f.write(f"- Mean log₂ Pearson R: "
                            f"{tdf['log2_pearson_r'].mean():.4f}\n")
                f.write(f"- Mean MAPE: {tdf['mape'].mean():.1f}%\n")
                f.write(f"- Mean RMSE: {tdf['rmse'].mean():.2f}\n")
                f.write(f"- Mean MAE: {tdf['mae'].mean():.2f}\n")
                f.write(f"- Mean WARE: {tdf['ware'].mean():.4f}\n")
                # Count-level metrics if available
                if "count_spearman_r" in tdf.columns and tdf["count_spearman_r"].notna().any():
                    f.write(f"- Count Spearman R: "
                            f"{tdf['count_spearman_r'].mean():.4f}\n")
                    f.write(f"- Count MAE: "
                            f"{tdf['count_mae'].mean():.1f}\n")
                    f.write(f"- Count RMSE: "
                            f"{tdf['count_rmse'].mean():.1f}\n")
                f.write("\n")

        # Full tables
        if not tx_df.empty:
            f.write("## Transcript-Level Metrics (All Conditions)\n\n")
            cols = [
                "condition", "tool", "gdna_label", "strand_specificity",
                "nrna_label", "spearman_r", "pearson_r",
                "log2_spearman_r", "log2_pearson_r",
                "rmse", "mae", "mape", "ware",
                "precision", "recall", "f1",
            ]
            # Add count columns if present
            count_cols = [
                "count_spearman_r", "count_pearson_r",
                "count_mae", "count_rmse",
            ]
            cols.extend(c for c in count_cols if c in tx_df.columns)
            _write_md_table(f, tx_df[[c for c in cols if c in tx_df.columns]])

        if not gene_df.empty:
            f.write("\n## Gene-Level Metrics (All Conditions)\n\n")
            cols = [
                "condition", "tool", "gdna_label", "strand_specificity",
                "nrna_label", "spearman_r", "pearson_r",
                "log2_spearman_r", "log2_pearson_r",
                "rmse", "mae", "mape", "ware",
            ]
            _write_md_table(f, gene_df[[c for c in cols if c in gene_df.columns]])

        if pool_df is not None and not pool_df.empty:
            f.write("\n## Pool-Level Summary (Rigel)\n\n")
            cols = [
                "condition", "tool",
                "mrna_frag_truth", "mrna_pred", "mrna_rel_error",
                "nrna_frag_truth", "nrna_pred", "nrna_rel_error",
                "gdna_frag_truth", "gdna_pred", "gdna_rel_error",
            ]
            _write_md_table(f, pool_df[[c for c in cols if c in pool_df.columns]])

        if strat_df is not None and not strat_df.empty:
            f.write("\n## Stratified Metrics by Expression Level\n\n")
            for tool in sorted(strat_df["tool"].unique()):
                f.write(f"### {tool}\n\n")
                tdf = strat_df[strat_df["tool"] == tool]
                pivot = tdf.groupby("expression_bin").agg(
                    mean_spearman=("spearman_r", "mean"),
                    mean_pearson=("pearson_r", "mean"),
                    mean_mape=("mape", "mean"),
                    mean_mae=("mae", "mean"),
                    mean_ware=("ware", "mean"),
                    n_conditions=("condition", "nunique"),
                ).reset_index()
                _write_md_table(f, pivot)
                f.write("\n")

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
            _write_md_table(f, rigel_df[[c for c in cols if c in rigel_df.columns]])

        # Cross-tool comparison
        if not tx_df.empty and len(tx_df["tool"].unique()) > 1:
            f.write("\n## Cross-Tool Comparison\n\n")
            for metric, label in [
                ("spearman_r", "Spearman R"),
                ("pearson_r", "Pearson R"),
                ("mape", "MAPE"),
                ("ware", "WARE"),
            ]:
                f.write(f"### Transcript-Level {label} by Condition\n\n")
                pivot = tx_df.pivot_table(
                    index="condition", columns="tool", values=metric
                ).reset_index()
                _write_md_table(f, pivot)
                f.write("\n")

        # Sensitivity analysis
        if not tx_df.empty:
            f.write("\n## Sensitivity Analysis\n\n")
            for tool in sorted(tx_df["tool"].unique()):
                tdf = tx_df[tx_df["tool"] == tool]
                f.write(f"### {tool}\n\n")
                for dim_col, dim_name in [
                    ("gdna_label", "gDNA Contamination Level"),
                    ("strand_specificity", "Strand Specificity"),
                    ("nrna_label", "nRNA Contamination"),
                ]:
                    if dim_col not in tdf.columns:
                        continue
                    f.write(f"#### By {dim_name}\n\n")
                    agg = tdf.groupby(dim_col).agg(
                        mean_pearson=("pearson_r", "mean"),
                        mean_mape=("mape", "mean"),
                        mean_ware=("ware", "mean"),
                        mean_rmse=("rmse", "mean"),
                        n=("condition", "count"),
                    ).reset_index()
                    _write_md_table(f, agg)
                    f.write("\n")

    logger.info("Wrote report: %s", path)


# ═══════════════════════════════════════════════════════════════════
# Main analysis pipeline
# ═══════════════════════════════════════════════════════════════════


def run_analysis(
    cfg: BenchmarkConfig,
    *,
    outdir: Path | None = None,
    condition_filter: list[str] | None = None,
) -> None:
    """Run full benchmark analysis.

    Parameters
    ----------
    cfg : BenchmarkConfig
        Benchmark configuration.
    outdir : Path or None
        Output directory (overrides config).
    condition_filter : list[str] or None
        Subset of conditions (overrides config).
    """
    outdir = outdir or Path(cfg.analysis.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    benchmark_dir = cfg.benchmark_dir

    logger.info("Benchmark directory: %s", benchmark_dir)
    logger.info("Output directory: %s", outdir)

    # Load manifest
    manifest = cfg.load_manifest()
    sim_params = manifest.get("simulation", {})
    logger.info(
        "Simulation: %d RNA fragments, seed=%d",
        sim_params.get("n_rna_fragments", 0),
        sim_params.get("sim_seed", 0),
    )

    # Load truth tables (legacy nrna-label format)
    truth_cache: dict[str, pd.DataFrame] = {}
    for nrna_label in ("none", "rand"):
        try:
            truth_cache[nrna_label] = load_truth(benchmark_dir, nrna_label)
            logger.info(
                "Loaded truth_%s: %d transcripts",
                nrna_label, len(truth_cache[nrna_label]),
            )
        except FileNotFoundError:
            logger.debug("Truth file for nrna_%s not found", nrna_label)

    # Discover conditions
    conditions = condition_filter or cfg.get_conditions()
    logger.info("Found %d conditions: %s", len(conditions), conditions)

    # Per-condition truth cache (loaded on demand below)
    cond_truth_cache: dict[str, pd.DataFrame] = {}

    # Precompute tx→gene mapping from truth
    tx2gene: dict[str, str] = {}
    for tdf in truth_cache.values():
        for row in tdf[["transcript_id", "gene_id"]].itertuples(index=False):
            tx2gene[row.transcript_id] = row.gene_id

    # Collectors
    all_tx_metrics: list[dict] = []
    all_gene_metrics: list[dict] = []
    all_pool_metrics: list[dict] = []
    all_rigel_details: list[dict] = []
    all_per_tx_rows: list[pd.DataFrame] = []
    all_stratified: list[dict] = []
    condition_summaries: list[dict] = []

    log2_pc = cfg.analysis.log2_pseudocount
    parse_fq = cfg.analysis.parse_fastq_truth

    for ci, condition in enumerate(conditions):
        meta = parse_condition(condition)
        nrna_label = meta["nrna_label"]

        # Try per-condition truth, then fall back to nrna-label truth
        if condition not in cond_truth_cache:
            try:
                cond_truth_cache[condition] = load_truth_for_condition(
                    benchmark_dir, condition, nrna_label,
                )
            except FileNotFoundError:
                pass
        truth_df = cond_truth_cache.get(condition)
        if truth_df is None:
            truth_df = truth_cache.get(nrna_label)
        if truth_df is None:
            logger.warning("No truth for %s (nrna=%s), skipping", condition, nrna_label)
            continue

        # Ensure tx2gene has entries from per-condition truth
        if condition in cond_truth_cache:
            for row in truth_df[["transcript_id", "gene_id"]].itertuples(index=False):
                tx2gene.setdefault(row.transcript_id, row.gene_id)

        print(
            f"[{ci + 1}/{len(conditions)}] {condition} "
            f"(gdna={meta['gdna_label']}, ss={meta['strand_specificity']}, "
            f"nrna={meta['nrna_label']})",
            flush=True,
        )

        cond_dir = cfg.condition_dir(condition)

        # FASTQ truth
        fastq_truth = None
        if parse_fq:
            r1 = cond_dir / "sim_R1.fq.gz"
            if r1.exists():
                logger.info("  Parsing FASTQ truth from %s...", r1)
                fastq_truth = parse_truth_from_fastq(r1)

        # Condition manifest metadata
        cond_manifest = {}
        for c in manifest.get("conditions", []):
            if c.get("name") == condition:
                cond_manifest = c
                break

        n_rna_total = cond_manifest.get("n_rna", 0)
        n_gdna_total = cond_manifest.get("n_gdna", 0)

        # Truth scaling
        truth_mrna_sum = float(truth_df["mrna_abundance"].sum())
        truth_nrna_sum = float(truth_df["nrna_abundance"].sum())
        truth_total_sum = truth_mrna_sum + truth_nrna_sum
        if truth_total_sum > 0 and n_rna_total > 0:
            truth_scale = n_rna_total / truth_total_sum
        else:
            truth_scale = 1.0

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

        # Discover tools for this condition
        tools = discover_tools(cond_dir)
        logger.info("  Tools: %s", [t.name for t in tools])

        for tool in tools:
            print(f"  {tool.name}...", end="", flush=True)

            # ── Transcript-level (TPM-to-TPM) ────────────────────
            tx_df = load_tool_tx(cond_dir, tool)
            if tx_df is not None:
                # Use TPM columns for comparison: truth mrna_abundance
                # is already in TPM units (sum ≈ 1,000,000).
                pred_tpm_col = tool.tx_tpm_col or tool.tx_count_col
                pred_count_col = tool.tx_count_col

                # Build merge columns
                merge_cols = [tool.tx_id_col, pred_tpm_col]
                rename_map = {
                    tool.tx_id_col: "transcript_id",
                    pred_tpm_col: "predicted",
                }
                # Also grab counts if available and different from TPM col
                has_counts = (
                    pred_count_col
                    and pred_count_col != pred_tpm_col
                    and pred_count_col in tx_df.columns
                )
                if has_counts:
                    merge_cols.append(pred_count_col)
                    rename_map[pred_count_col] = "pred_count"

                merged = truth_df.merge(
                    tx_df[merge_cols].rename(columns=rename_map),
                    on="transcript_id",
                    how="left",
                )
                merged["predicted"] = merged["predicted"].fillna(0.0)
                if has_counts:
                    merged["pred_count"] = merged["pred_count"].fillna(0.0)
                elif pred_count_col and pred_count_col in tx_df.columns:
                    # count col is same as tpm col (e.g. salmon NumReads)
                    pass

                # Renormalize truth mRNA TPM to sum to 1M for fair comparison.
                # In nrna_rand conditions, truth mrna_abundance sums to <1M
                # (because some TPM budget goes to nRNA), but tool TPM always
                # sums to ~1M. Without renormalization, every transcript would
                # appear systematically overestimated.
                truth_raw = merged["mrna_abundance"].values.astype(np.float64)
                truth_sum = truth_raw.sum()
                if truth_sum > 0:
                    truth_arr = truth_raw * (1e6 / truth_sum)
                else:
                    truth_arr = truth_raw
                pred_arr = merged["predicted"].values.astype(np.float64)

                metrics = compute_metrics(truth_arr, pred_arr, pseudocount=log2_pc)
                metrics["condition"] = condition
                metrics["tool"] = tool.name
                metrics.update(meta)
                all_tx_metrics.append(metrics)

                # Count-based metrics (if tool provides counts)
                count_metrics = None
                if has_counts:
                    pred_counts = merged["pred_count"].values.astype(np.float64)
                    # Compute truth counts: scale mrna_abundance by truth_scale
                    truth_counts = truth_raw * truth_scale
                    count_metrics = compute_metrics(truth_counts, pred_counts)
                    # Prefix count metrics
                    for k, v in count_metrics.items():
                        metrics[f"count_{k}"] = v

                # Stratified
                strat = compute_stratified_metrics(truth_arr, pred_arr)
                for bin_name, bin_metrics in strat.items():
                    bin_metrics["condition"] = condition
                    bin_metrics["tool"] = tool.name
                    bin_metrics["expression_bin"] = bin_name
                    bin_metrics.update(meta)
                    all_stratified.append(bin_metrics)

                # Per-transcript detail (TPM, renormalized)
                truth_norm_factor = (1e6 / truth_sum) if truth_sum > 0 else 1.0
                merged["truth_tpm"] = merged["mrna_abundance"] * truth_norm_factor
                detail_mask = (merged["truth_tpm"] > 0) | (merged["predicted"] > 0)
                if detail_mask.any():
                    detail_cols = [
                        "transcript_id", "gene_id", "gene_name",
                        "mrna_abundance", "nrna_abundance",
                        "truth_tpm", "predicted",
                    ]
                    detail = merged.loc[detail_mask, detail_cols].copy()
                    detail["condition"] = condition
                    detail["tool"] = tool.name
                    detail["residual"] = detail["predicted"] - detail["truth_tpm"]
                    detail["abs_error"] = detail["residual"].abs()
                    detail["rel_error"] = np.where(
                        detail["truth_tpm"] > 0,
                        detail["abs_error"] / detail["truth_tpm"],
                        np.nan,
                    )
                    # Add count columns if available
                    if has_counts:
                        detail["pred_count"] = merged.loc[
                            detail_mask, "pred_count"
                        ].values
                        detail["truth_count"] = (
                            merged.loc[detail_mask, "mrna_abundance"].values
                            * truth_scale
                        )
                        detail["count_abs_error"] = np.abs(
                            detail["pred_count"] - detail["truth_count"]
                        )
                    all_per_tx_rows.append(detail)

            # ── Gene-level (TPM-to-TPM) ────────────────────────
            # Always aggregate transcript-level TPM to gene level for
            # consistency. Rigel's gene_quant.feather uses a different
            # effective_length which produces gene TPM ≠ sum(tx TPM).
            if "gene_id" in truth_df.columns and tx_df is not None:
                gene_truth = (
                    truth_df.groupby("gene_id")["mrna_abundance"]
                    .sum()
                    .reset_index()
                )

                tx_pred_col = tool.tx_tpm_col or tool.tx_count_col
                gene_pred: dict[str, float] = {}
                for _, row in tx_df.iterrows():
                    tid = row[tool.tx_id_col]
                    gid = tx2gene.get(tid)
                    if gid is not None:
                        gene_pred[gid] = (
                            gene_pred.get(gid, 0.0) + float(row[tx_pred_col])
                        )

                if gene_pred:
                    gene_truth["predicted"] = gene_truth["gene_id"].map(
                        lambda g: gene_pred.get(g, 0.0)
                    )
                    # Renormalize gene-level truth TPM to 1M (same as tx-level)
                    g_raw = gene_truth["mrna_abundance"].values.astype(np.float64)
                    g_sum = g_raw.sum()
                    g_truth = g_raw * (1e6 / g_sum) if g_sum > 0 else g_raw
                    g_pred = gene_truth["predicted"].values.astype(np.float64)

                    g_metrics = compute_metrics(g_truth, g_pred, pseudocount=log2_pc)
                    g_metrics["condition"] = condition
                    g_metrics["tool"] = tool.name
                    g_metrics.update(meta)
                    all_gene_metrics.append(g_metrics)

            # ── Pool-level (rigel only) ──────────────────────────
            if tool.name.startswith("rigel"):
                pool = compute_pool_metrics(
                    cond_dir, condition, tool, truth_df, manifest, fastq_truth
                )
                if pool:
                    all_pool_metrics.append(pool)

                detail = rigel_detailed_analysis(cond_dir, condition, tool)
                if detail:
                    all_rigel_details.append(detail)

            print(" done", flush=True)

    # ── Write outputs ────────────────────────────────────────────
    print("\nWriting outputs...", flush=True)

    tx_metrics_df = pd.DataFrame(all_tx_metrics)
    tx_metrics_df.to_csv(outdir / "transcript_metrics.csv", index=False)

    gene_metrics_df = pd.DataFrame(all_gene_metrics)
    gene_metrics_df.to_csv(outdir / "gene_metrics.csv", index=False)

    if all_pool_metrics:
        pd.DataFrame(all_pool_metrics).to_csv(outdir / "pool_summary.csv", index=False)

    if all_rigel_details:
        pd.DataFrame(all_rigel_details).to_csv(outdir / "rigel_details.csv", index=False)

    if all_stratified:
        pd.DataFrame(all_stratified).to_csv(outdir / "stratified_metrics.csv", index=False)

    if all_per_tx_rows:
        per_tx = pd.concat(all_per_tx_rows, ignore_index=True)
        per_tx.to_parquet(outdir / "per_transcript_detail.parquet", index=False)

    pd.DataFrame(condition_summaries).to_csv(outdir / "condition_summary.csv", index=False)

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
