#!/usr/bin/env python3
"""Publication-quality benchmark report generator.

Generates matplotlib figures and an integrated HTML report from
benchmark analysis outputs (transcript_metrics.csv, gene_metrics.csv,
per_transcript_detail.parquet, etc.).

Usage:
    python -m scripts.benchmarking report -c config.yaml -o results/report
    python -m scripts.benchmarking report -c config.yaml -o results/report --figs-only
"""
from __future__ import annotations

import json
import logging
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# ── Palette ────────────────────────────────────────────────────────
# Colorblind-safe palette (Wong, 2011) extended for tools.
TOOL_PALETTE = {
    "rigel/oracle_vbem": "#000000",
    "rigel/oracle_map": "#444444",
    "rigel/star_vbem": "#0072B2",
    "rigel/star_map": "#56B4E9",
    "rigel/vbem": "#009E73",
    "rigel/map": "#CC79A7",
    "kallisto": "#E69F00",
    "salmon": "#D55E00",
}

# Display names for tools
TOOL_LABELS = {
    "rigel/oracle_vbem": "Rigel Oracle (VBEM)",
    "rigel/oracle_map": "Rigel Oracle (MAP)",
    "rigel/star_vbem": "Rigel+STAR (VBEM)",
    "rigel/star_map": "Rigel+STAR (MAP)",
    "rigel/vbem": "Rigel+mm2 (VBEM)",
    "rigel/map": "Rigel+mm2 (MAP)",
    "kallisto": "Kallisto",
    "salmon": "Salmon",
}


def _tool_color(tool: str) -> str:
    return TOOL_PALETTE.get(tool, "#888888")


def _tool_label(tool: str) -> str:
    return TOOL_LABELS.get(tool, tool)


def _setup_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.titleweight": "bold",
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linewidth": 0.5,
    })


# ═══════════════════════════════════════════════════════════════════
# Data loading
# ═══════════════════════════════════════════════════════════════════


def load_analysis(analysis_dir: Path) -> dict:
    """Load all analysis outputs from a results directory."""
    data = {}
    for name, fname, loader in [
        ("tx_metrics", "transcript_metrics.csv", pd.read_csv),
        ("gene_metrics", "gene_metrics.csv", pd.read_csv),
        ("stratified", "stratified_metrics.csv", pd.read_csv),
        ("pool", "pool_summary.csv", pd.read_csv),
        ("rigel_details", "rigel_details.csv", pd.read_csv),
        ("condition_summary", "condition_summary.csv", pd.read_csv),
    ]:
        path = analysis_dir / fname
        if path.exists():
            data[name] = loader(path)
        else:
            data[name] = pd.DataFrame()

    pq_path = analysis_dir / "per_transcript_detail.parquet"
    if pq_path.exists():
        data["per_tx"] = pd.read_parquet(pq_path)
    else:
        data["per_tx"] = pd.DataFrame()

    return data


def _get_tools(data: dict) -> list[str]:
    """Get ordered list of tools present in the data."""
    if data["tx_metrics"].empty:
        return []
    tools = data["tx_metrics"]["tool"].unique().tolist()
    # Order: oracle first, then star, then mm2, then external
    priority = {
        "rigel/oracle_vbem": 0, "rigel/oracle_map": 1,
        "rigel/star_vbem": 2, "rigel/star_map": 3,
        "rigel/vbem": 4, "rigel/map": 5,
        "kallisto": 6, "salmon": 7,
    }
    return sorted(tools, key=lambda t: priority.get(t, 99))


# ═══════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════


def fig_scatter_grid(data: dict, out_dir: Path, condition: str | None = None):
    """Transcript-level scatter: truth vs predicted log₂(TPM+1)."""
    per_tx = data["per_tx"]
    if per_tx.empty:
        return
    if condition:
        per_tx = per_tx[per_tx.condition == condition]

    tools = _get_tools(data)
    tools_in_data = [t for t in tools if t in per_tx.tool.unique()]
    n = len(tools_in_data)
    if n == 0:
        return

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    pc = 1.0
    for ax, tool in zip(axes, tools_in_data):
        sub = per_tx[per_tx.tool == tool]
        truth_log2 = np.log2(sub["truth_tpm"].values + pc)
        pred_log2 = np.log2(sub["predicted"].values + pc)

        # Spearman R
        from scipy.stats import spearmanr
        sr, _ = spearmanr(truth_log2, pred_log2)
        pr = np.corrcoef(truth_log2, pred_log2)[0, 1]

        hb = ax.hexbin(
            truth_log2, pred_log2,
            gridsize=80, cmap="YlOrRd", mincnt=1,
            norm=mcolors.LogNorm(vmin=1, vmax=5000),
            linewidths=0.1, edgecolors="face",
        )
        lims = [-0.5, 14]
        ax.plot(lims, lims, "k-", lw=1, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Truth log₂(TPM + 1)")
        ax.set_ylabel("Predicted log₂(TPM + 1)")
        ax.set_title(f"{_tool_label(tool)}\nρ = {sr:.4f}  (r = {pr:.4f})")
        ax.set_aspect("equal")

    fig.suptitle(
        "Transcript-Level Accuracy — log₂(TPM + 1)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig_scatter_grid.png")
    plt.close(fig)
    logger.info("  scatter grid")


def fig_metrics_comparison(data: dict, out_dir: Path, level: str = "tx"):
    """Grouped bar chart comparing key metrics across tools."""
    df = data["tx_metrics"] if level == "tx" else data["gene_metrics"]
    if df.empty:
        return

    tools = _get_tools(data)
    tools = [t for t in tools if t in df.tool.unique()]
    if not tools:
        return

    metrics_spec = [
        ("spearman_r", "Spearman ρ", True),
        ("pearson_r", "Pearson r", True),
        ("ware", "WARE", False),
        ("mae", "MAE (TPM)", False),
        ("mape", "MAPE (%)", False),
    ]
    if level == "tx":
        metrics_spec.append(("f1", "F1 Score", True))

    n_metrics = len(metrics_spec)
    fig, axes = plt.subplots(1, n_metrics, figsize=(4.5 * n_metrics, 5))

    for ax, (col, label, higher_better) in zip(axes, metrics_spec):
        if col not in df.columns:
            ax.set_visible(False)
            continue
        vals = []
        for t in tools:
            v = df[df.tool == t][col].mean()
            vals.append(v)

        bars = ax.bar(
            range(len(tools)), vals,
            color=[_tool_color(t) for t in tools],
            edgecolor="white", linewidth=0.5,
        )
        ax.set_xticks(range(len(tools)))
        ax.set_xticklabels([_tool_label(t) for t in tools], rotation=35, ha="right")
        ax.set_ylabel(label)
        ax.set_title(label)

        # Value annotations
        for bar, v in zip(bars, vals):
            fmt = f"{v:.4f}" if abs(v) < 10 else f"{v:.1f}"
            if col == "mape":
                fmt = f"{v:.1f}%"
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (bar.get_height() * 0.02 + 0.001),
                fmt, ha="center", va="bottom", fontsize=8, fontweight="bold",
            )

        # Auto y-limits
        if higher_better and all(v > 0.5 for v in vals if not np.isnan(v)):
            ax.set_ylim(min(vals) * 0.95, 1.02)

    level_label = "Transcript" if level == "tx" else "Gene"
    fig.suptitle(
        f"{level_label}-Level Metrics Comparison",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / f"fig_{level}_metrics.png")
    plt.close(fig)
    logger.info("  %s metrics comparison", level)


def fig_stratified_bars(data: dict, out_dir: Path):
    """Expression-stratified Spearman R and WARE."""
    strat = data["stratified"]
    if strat.empty:
        return

    tools = _get_tools(data)
    tools = [t for t in tools if t in strat.tool.unique()]

    bin_order = ["low_1-10", "mid_10-100", "high_100-1000", "very_high_1000+"]
    bin_labels = ["1–10", "10–100", "100–1K", "≥1K"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    metrics = [
        ("spearman_r", "Spearman ρ"),
        ("pearson_r", "Pearson r"),
        ("ware", "WARE"),
    ]

    for ax, (metric, ylabel) in zip(axes, metrics):
        x = np.arange(len(bin_order))
        w = 0.8 / len(tools)
        offsets = np.linspace(-0.4 + w / 2, 0.4 - w / 2, len(tools))

        for i, tool in enumerate(tools):
            tool_data = strat[strat.tool == tool]
            vals = []
            for b in bin_order:
                row = tool_data[tool_data.expression_bin == b]
                vals.append(float(row[metric].iloc[0]) if len(row) else 0)
            ax.bar(
                x + offsets[i], vals, w,
                label=_tool_label(tool), color=_tool_color(tool),
                edgecolor="white", linewidth=0.3,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(bin_labels)
        ax.set_xlabel("Truth Expression Bin (TPM)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{ylabel} by Expression Level")
        ax.legend(fontsize=7, ncol=2)

    fig.suptitle(
        "Expression-Stratified Performance",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig_stratified.png")
    plt.close(fig)
    logger.info("  stratified bars")


def fig_pool_level(data: dict, out_dir: Path):
    """Pool-level fragment assignment for Rigel tools."""
    pool = data["pool"]
    if pool.empty:
        return

    rigel_tools = [t for t in pool.tool.unique() if t.startswith("rigel")]
    if not rigel_tools:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Stacked bar: fragment decomposition
    ax = axes[0]
    bar_labels = ["Truth"] + [_tool_label(t) for t in rigel_tools]
    n_bars = len(bar_labels)

    # Build arrays
    mrna = [float(pool.iloc[0].get("mrna_frag_truth", 0))]
    nrna = [float(pool.iloc[0].get("nrna_frag_truth", 0))]
    gdna = [float(pool.iloc[0].get("gdna_frag_truth", 0))]
    intergenic = [0.0]

    for tool in rigel_tools:
        row = pool[pool.tool == tool].iloc[0]
        mrna.append(float(row["mrna_pred"]))
        nrna.append(float(row["nrna_pred"]))
        gdna.append(float(row.get("gdna_em_pred", 0)))
        intergenic.append(float(row.get("intergenic_pred", 0)))

    mrna = np.array(mrna) / 1e6
    nrna = np.array(nrna) / 1e6
    gdna = np.array(gdna) / 1e6
    intergenic = np.array(intergenic) / 1e6

    x = np.arange(n_bars)
    w = 0.6
    ax.bar(x, mrna, w, label="mRNA", color="#0072B2")
    ax.bar(x, nrna, w, bottom=mrna, label="nRNA", color="#E69F00")
    ax.bar(x, gdna, w, bottom=mrna + nrna, label="gDNA", color="#D55E00")
    ax.bar(x, intergenic, w, bottom=mrna + nrna + gdna,
           label="Intergenic", color="#999999")
    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels, rotation=20, ha="right")
    ax.set_ylabel("Fragments (millions)")
    ax.set_title("Fragment Pool Decomposition")
    ax.legend(loc="upper right")

    # Relative error bars
    ax = axes[1]
    cats = ["mRNA", "gDNA"]
    tool_data = {}
    for tool in rigel_tools:
        row = pool[pool.tool == tool].iloc[0]
        mrna_re = float(row.get("mrna_rel_error", 0)) * 100
        gdna_re = float(row.get("gdna_rel_error", 0)) * 100
        tool_data[tool] = [mrna_re, gdna_re]

    x = np.arange(len(cats))
    w2 = 0.8 / len(rigel_tools)
    offsets = np.linspace(-0.4 + w2 / 2, 0.4 - w2 / 2, len(rigel_tools))
    for i, tool in enumerate(rigel_tools):
        bars = ax.bar(
            x + offsets[i], tool_data[tool], w2,
            label=_tool_label(tool), color=_tool_color(tool),
            edgecolor="white", linewidth=0.3,
        )
        for bar, v in zip(bars, tool_data[tool]):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + (-1.5 if v < 0 else 0.3),
                f"{v:.1f}%", ha="center", va="top" if v < 0 else "bottom",
                fontsize=8, fontweight="bold",
            )
    ax.set_xticks(x)
    ax.set_xticklabels(cats)
    ax.set_ylabel("Relative Error (%)")
    ax.set_title("Pool-Level Estimation Error")
    ax.legend()
    ax.axhline(y=0, color="black", lw=0.5)

    fig.suptitle("Pool-Level Fragment Assignment", fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "fig_pool_level.png")
    plt.close(fig)
    logger.info("  pool level")


def fig_precision_recall(data: dict, out_dir: Path):
    """Precision vs recall scatter."""
    df = data["tx_metrics"]
    if df.empty or "precision" not in df.columns:
        return

    tools = _get_tools(data)
    tools = [t for t in tools if t in df.tool.unique()]

    fig, ax = plt.subplots(figsize=(8, 6))

    for tool in tools:
        row = df[df.tool == tool].iloc[0]
        if pd.isna(row.get("precision")) or pd.isna(row.get("recall")):
            continue
        ax.scatter(
            row["recall"], row["precision"],
            s=200, color=_tool_color(tool), zorder=5,
            edgecolors="black", linewidth=0.8,
        )
        ax.annotate(
            f"  {_tool_label(tool)}\n  F1={row['f1']:.3f}",
            (row["recall"], row["precision"]),
            fontsize=9, va="center",
        )

    # F1 iso-lines
    for f1_val in [0.7, 0.75, 0.8, 0.85, 0.9]:
        r_range = np.linspace(0.01, 1.0, 200)
        p_line = (f1_val * r_range) / (2 * r_range - f1_val)
        valid = (p_line > 0) & (p_line <= 1)
        ax.plot(r_range[valid], p_line[valid], "k-", alpha=0.1, lw=0.8)

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title("Transcript Detection — Precision vs Recall\n(threshold: TPM > 0)")
    ax.set_xlim(0.5, 1.0)
    ax.set_ylim(0.5, 1.0)
    plt.tight_layout()
    fig.savefig(out_dir / "fig_precision_recall.png")
    plt.close(fig)
    logger.info("  precision-recall")


def fig_gene_scatter(data: dict, out_dir: Path, condition: str | None = None):
    """Gene-level scatter: truth vs predicted log₂(TPM+1)."""
    per_tx = data["per_tx"]
    if per_tx.empty:
        return

    if condition:
        per_tx = per_tx[per_tx.condition == condition]

    tools = _get_tools(data)
    tools_in_data = [t for t in tools if t in per_tx.tool.unique()]
    n = len(tools_in_data)
    if n == 0:
        return

    # Aggregate to gene level
    gene_dfs = {}
    for tool in tools_in_data:
        sub = per_tx[per_tx.tool == tool]
        gene = sub.groupby("gene_id").agg(
            truth_tpm=("truth_tpm", "sum"),
            predicted=("predicted", "sum"),
        ).reset_index()
        gene_dfs[tool] = gene

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    pc = 1.0
    for ax, tool in zip(axes, tools_in_data):
        gene = gene_dfs[tool]
        truth_log2 = np.log2(gene["truth_tpm"].values + pc)
        pred_log2 = np.log2(gene["predicted"].values + pc)

        from scipy.stats import spearmanr
        sr, _ = spearmanr(truth_log2, pred_log2)
        pr = np.corrcoef(truth_log2, pred_log2)[0, 1]

        hb = ax.hexbin(
            truth_log2, pred_log2,
            gridsize=60, cmap="YlGn", mincnt=1,
            norm=mcolors.LogNorm(vmin=1, vmax=2000),
            linewidths=0.1, edgecolors="face",
        )
        lims = [-0.5, 16]
        ax.plot(lims, lims, "k-", lw=1, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Truth log₂(TPM + 1)")
        ax.set_ylabel("Predicted log₂(TPM + 1)")
        ax.set_title(f"{_tool_label(tool)}\nρ = {sr:.4f}  (r = {pr:.4f})")
        ax.set_aspect("equal")

    fig.suptitle(
        "Gene-Level Accuracy — log₂(TPM + 1)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig_gene_scatter.png")
    plt.close(fig)
    logger.info("  gene scatter grid")


def fig_count_scatter(data: dict, out_dir: Path, condition: str | None = None):
    """Count-level scatter: truth vs predicted counts (log₁₀)."""
    per_tx = data["per_tx"]
    if per_tx.empty or "pred_count" not in per_tx.columns:
        return

    if condition:
        per_tx = per_tx[per_tx.condition == condition]

    tools = _get_tools(data)
    tools_in_data = [
        t for t in tools
        if t in per_tx.tool.unique()
        and per_tx.loc[per_tx.tool == t, "pred_count"].notna().any()
    ]
    n = len(tools_in_data)
    if n == 0:
        return

    fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
    if n == 1:
        axes = [axes]

    for ax, tool in zip(axes, tools_in_data):
        sub = per_tx[(per_tx.tool == tool) & per_tx.pred_count.notna()]
        truth_log = np.log10(sub["truth_count"].values + 1)
        pred_log = np.log10(sub["pred_count"].values + 1)

        from scipy.stats import spearmanr
        sr, _ = spearmanr(truth_log, pred_log)

        hb = ax.hexbin(
            truth_log, pred_log,
            gridsize=60, cmap="YlOrBr", mincnt=1,
            norm=mcolors.LogNorm(vmin=1, vmax=5000),
            linewidths=0.1, edgecolors="face",
        )
        lims = [-0.2, 6]
        ax.plot(lims, lims, "k-", lw=1, alpha=0.5)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Truth log₁₀(count + 1)")
        ax.set_ylabel("Predicted log₁₀(count + 1)")
        ax.set_title(f"{_tool_label(tool)}\nρ = {sr:.4f}")
        ax.set_aspect("equal")

    fig.suptitle(
        "Transcript-Level Count Accuracy — log₁₀(count + 1)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig_count_scatter.png")
    plt.close(fig)
    logger.info("  count scatter grid")


# ═══════════════════════════════════════════════════════════════════
# HTML Report
# ═══════════════════════════════════════════════════════════════════


def _fmt_metric(v, fmt=".4f"):
    if pd.isna(v):
        return "—"
    return f"{v:{fmt}}"


def _build_metrics_table_html(
    df: pd.DataFrame,
    tools: list[str],
    metrics: list[tuple[str, str, str]],  # (col, label, fmt)
) -> str:
    """Build an HTML metrics table."""
    rows = []
    for col, label, fmt in metrics:
        if col not in df.columns:
            continue
        cells = f"<td><b>{label}</b></td>"
        vals = []
        for t in tools:
            v = df[df.tool == t][col].mean() if t in df.tool.values else np.nan
            vals.append(v)
        # Find best value
        valid = [(i, v) for i, v in enumerate(vals) if not np.isnan(v)]
        if valid:
            if col in ("spearman_r", "pearson_r", "f1", "precision", "recall",
                        "log2_spearman_r", "log2_pearson_r"):
                best_idx = max(valid, key=lambda x: x[1])[0]
            else:
                best_idx = min(valid, key=lambda x: x[1])[0]
        else:
            best_idx = -1

        for i, v in enumerate(vals):
            bold = " font-weight:bold; color:#0072B2;" if i == best_idx else ""
            cells += f'<td style="text-align:right;{bold}">{_fmt_metric(v, fmt)}</td>'
        rows.append(f"<tr>{cells}</tr>")
    return "\n".join(rows)


def generate_html_report(data: dict, out_dir: Path, title: str = "Benchmark Report"):
    """Generate integrated HTML report with embedded figures and tables."""
    tools = _get_tools(data)
    tx = data["tx_metrics"]
    gene = data["gene_metrics"]
    pool = data["pool"]
    strat = data["stratified"]
    condition_summary = data.get("condition_summary", pd.DataFrame())

    # Figure paths (relative to report)
    figs = [f for f in sorted(out_dir.glob("fig_*.png"))]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{title}</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
       max-width: 1200px; margin: 0 auto; padding: 20px 40px; color: #333; line-height: 1.6; }}
h1 {{ font-size: 28px; margin: 30px 0 10px; border-bottom: 2px solid #0072B2; padding-bottom: 8px; }}
h2 {{ font-size: 22px; margin: 25px 0 10px; color: #0072B2; }}
h3 {{ font-size: 16px; margin: 20px 0 8px; }}
table {{ border-collapse: collapse; margin: 15px 0; width: 100%; font-size: 13px; }}
th {{ background: #f0f4f8; padding: 8px 12px; text-align: right; border-bottom: 2px solid #ccc; font-weight: 600; }}
th:first-child {{ text-align: left; }}
td {{ padding: 6px 12px; border-bottom: 1px solid #eee; text-align: right; }}
td:first-child {{ text-align: left; font-weight: 500; }}
tr:hover {{ background: #f8fbfd; }}
img {{ max-width: 100%; margin: 15px 0; border: 1px solid #eee; border-radius: 4px; }}
.figure-caption {{ font-size: 12px; color: #666; margin: -10px 0 20px; font-style: italic; }}
.metric-highlight {{ background: #e8f4fd; padding: 2px 6px; border-radius: 3px; font-weight: bold; }}
.summary-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 15px; margin: 15px 0; }}
.summary-card {{ background: #f8fbfd; border: 1px solid #dee; border-radius: 6px; padding: 15px; }}
.summary-card h4 {{ margin: 0 0 8px; color: #0072B2; font-size: 14px; }}
.summary-card .value {{ font-size: 24px; font-weight: bold; }}
.summary-card .label {{ font-size: 12px; color: #666; }}
.timestamp {{ color: #999; font-size: 12px; margin-bottom: 20px; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="timestamp">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
"""

    # ── Simulation parameters ──
    if not condition_summary.empty:
        html += "<h2>Simulation Parameters</h2>\n"
        row = condition_summary.iloc[0]
        html += '<div class="summary-grid">\n'
        for label, key in [("Conditions", None), ("gDNA", "gdna_label"),
                           ("Strand Specificity", "strand_specificity"),
                           ("nRNA", "nrna_label")]:
            if key and key in row.index:
                val = row[key]
            elif key is None:
                val = len(condition_summary)
            else:
                continue
            html += f'<div class="summary-card"><h4>{label}</h4><div class="value">{val}</div></div>\n'
        html += '</div>\n'

    # ── Summary cards ──
    if not tx.empty:
        html += "<h2>Summary</h2>\n"
        html += '<div class="summary-grid">\n'
        for tool in tools:
            t_row = tx[tx.tool == tool]
            if t_row.empty:
                continue
            sr = t_row["spearman_r"].mean()
            pr = t_row["pearson_r"].mean()
            ware = t_row["ware"].mean()
            html += f"""<div class="summary-card">
<h4>{_tool_label(tool)}</h4>
<div class="value">ρ = {sr:.4f}</div>
<div class="label">Spearman | r = {pr:.4f} Pearson | WARE = {ware:.4f}</div>
</div>\n"""
        html += '</div>\n'

    # ── Transcript-level metrics table ──
    if not tx.empty:
        html += "<h2>Transcript-Level Metrics</h2>\n"
        html += "<table>\n<thead><tr><th>Metric</th>"
        for t in tools:
            html += f"<th>{_tool_label(t)}</th>"
        html += "</tr></thead>\n<tbody>\n"

        tx_metrics_spec = [
            ("spearman_r", "Spearman ρ", ".4f"),
            ("pearson_r", "Pearson r", ".4f"),
            ("log2_spearman_r", "log₂ Spearman ρ", ".4f"),
            ("log2_pearson_r", "log₂ Pearson r", ".4f"),
            ("ware", "WARE", ".4f"),
            ("mae", "MAE (TPM)", ".3f"),
            ("rmse", "RMSE (TPM)", ".2f"),
            ("mape", "MAPE (%)", ".1f"),
            ("f1", "Detection F1", ".3f"),
            ("precision", "Precision", ".3f"),
            ("recall", "Recall", ".3f"),
        ]
        # Count metrics
        if "count_spearman_r" in tx.columns:
            tx_metrics_spec.extend([
                ("count_spearman_r", "Count Spearman ρ", ".4f"),
                ("count_pearson_r", "Count Pearson r", ".4f"),
                ("count_mae", "Count MAE", ".1f"),
                ("count_rmse", "Count RMSE", ".1f"),
            ])

        html += _build_metrics_table_html(tx, tools, tx_metrics_spec)
        html += "\n</tbody></table>\n"

    # ── Gene-level metrics table ──
    if not gene.empty:
        html += "<h2>Gene-Level Metrics</h2>\n"
        html += "<table>\n<thead><tr><th>Metric</th>"
        for t in tools:
            html += f"<th>{_tool_label(t)}</th>"
        html += "</tr></thead>\n<tbody>\n"

        gene_metrics_spec = [
            ("spearman_r", "Spearman ρ", ".4f"),
            ("pearson_r", "Pearson r", ".4f"),
            ("ware", "WARE", ".4f"),
            ("mae", "MAE (TPM)", ".3f"),
            ("rmse", "RMSE (TPM)", ".2f"),
            ("mape", "MAPE (%)", ".1f"),
        ]
        html += _build_metrics_table_html(gene, tools, gene_metrics_spec)
        html += "\n</tbody></table>\n"

    # ── Pool-level table ──
    if not pool.empty:
        html += "<h2>Pool-Level Fragment Assignment (Rigel)</h2>\n"
        html += "<table>\n<thead><tr>"
        html += "<th>Tool</th><th>mRNA pred</th><th>mRNA err</th>"
        html += "<th>nRNA pred</th><th>gDNA pred</th><th>gDNA err</th>"
        html += "</tr></thead>\n<tbody>\n"
        for _, row in pool.iterrows():
            mrna_err = row.get("mrna_rel_error", np.nan)
            gdna_err = row.get("gdna_rel_error", np.nan)
            html += f"""<tr>
<td>{_tool_label(row['tool'])}</td>
<td>{row['mrna_pred']:,.0f}</td>
<td>{_fmt_metric(mrna_err * 100 if pd.notna(mrna_err) else np.nan, '.2f')}%</td>
<td>{row['nrna_pred']:,.0f}</td>
<td>{row['gdna_pred']:,.0f}</td>
<td>{_fmt_metric(gdna_err * 100 if pd.notna(gdna_err) else np.nan, '.1f')}%</td>
</tr>\n"""
        html += "</tbody></table>\n"

    # ── Stratified metrics table ──
    if not strat.empty:
        html += "<h2>Expression-Stratified Metrics</h2>\n"
        bin_order = ["low_1-10", "mid_10-100", "high_100-1000", "very_high_1000+"]
        bin_labels_map = {
            "low_1-10": "1–10 TPM",
            "mid_10-100": "10–100 TPM",
            "high_100-1000": "100–1K TPM",
            "very_high_1000+": "≥1K TPM",
        }

        for metric, label in [("spearman_r", "Spearman ρ"), ("ware", "WARE")]:
            html += f"<h3>{label} by Expression Level</h3>\n"
            html += "<table>\n<thead><tr><th>Expression</th>"
            for t in tools:
                html += f"<th>{_tool_label(t)}</th>"
            html += "</tr></thead>\n<tbody>\n"

            for b in bin_order:
                html += f"<tr><td>{bin_labels_map.get(b, b)}</td>"
                vals = []
                for t in tools:
                    row = strat[(strat.tool == t) & (strat.expression_bin == b)]
                    v = float(row[metric].iloc[0]) if len(row) else np.nan
                    vals.append(v)
                # Best value
                valid = [(i, v) for i, v in enumerate(vals) if not np.isnan(v)]
                if valid:
                    if metric in ("spearman_r", "pearson_r"):
                        best_idx = max(valid, key=lambda x: x[1])[0]
                    else:
                        best_idx = min(valid, key=lambda x: x[1])[0]
                else:
                    best_idx = -1

                for i, v in enumerate(vals):
                    bold = " font-weight:bold; color:#0072B2;" if i == best_idx else ""
                    html += f'<td style="text-align:right;{bold}">{_fmt_metric(v)}</td>'
                html += "</tr>\n"
            html += "</tbody></table>\n"

    # ── Figures ──
    if figs:
        html += "<h2>Figures</h2>\n"
        fig_titles = {
            "fig_scatter_grid.png": "Transcript-Level Scatter (log₂ TPM + 1)",
            "fig_tx_metrics.png": "Transcript-Level Metrics Comparison",
            "fig_gene_metrics.png": "Gene-Level Metrics Comparison",
            "fig_gene_scatter.png": "Gene-Level Scatter (log₂ TPM + 1)",
            "fig_stratified.png": "Expression-Stratified Performance",
            "fig_pool_level.png": "Pool-Level Fragment Assignment",
            "fig_precision_recall.png": "Precision vs Recall",
            "fig_count_scatter.png": "Transcript-Level Count Accuracy",
        }
        for fig_path in figs:
            cap = fig_titles.get(fig_path.name, fig_path.stem)
            html += f'<h3>{cap}</h3>\n'
            html += f'<img src="{fig_path.name}" alt="{cap}">\n'

    html += "\n</body></html>"

    report_path = out_dir / "report.html"
    with open(report_path, "w") as f:
        f.write(html)
    logger.info("Wrote HTML report: %s", report_path)


# ═══════════════════════════════════════════════════════════════════
# Main entry point
# ═══════════════════════════════════════════════════════════════════


def generate_report(
    analysis_dir: Path,
    out_dir: Path | None = None,
    *,
    title: str = "Benchmark Report",
    figs_only: bool = False,
) -> None:
    """Generate publication-quality figures and HTML report.

    Parameters
    ----------
    analysis_dir : Path
        Directory containing analysis outputs (transcript_metrics.csv, etc.)
    out_dir : Path or None
        Output directory for figures + HTML. Defaults to analysis_dir.
    title : str
        Report title.
    figs_only : bool
        If True, generate only figures (no HTML report).
    """
    out_dir = out_dir or analysis_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    _setup_style()

    logger.info("Loading analysis data from %s", analysis_dir)
    data = load_analysis(analysis_dir)
    if data["tx_metrics"].empty:
        logger.error("No transcript metrics found in %s", analysis_dir)
        return

    conditions = data["tx_metrics"]["condition"].unique()
    condition = conditions[0] if len(conditions) == 1 else None

    logger.info("Generating figures...")
    fig_scatter_grid(data, out_dir, condition)
    fig_metrics_comparison(data, out_dir, "tx")
    fig_metrics_comparison(data, out_dir, "gene")
    fig_gene_scatter(data, out_dir, condition)
    fig_stratified_bars(data, out_dir)
    fig_pool_level(data, out_dir)
    fig_precision_recall(data, out_dir)
    fig_count_scatter(data, out_dir, condition)

    if not figs_only:
        logger.info("Generating HTML report...")
        generate_html_report(data, out_dir, title=title)

    logger.info("Report generation complete: %s", out_dir)
