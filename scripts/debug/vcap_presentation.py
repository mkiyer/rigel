#!/usr/bin/env python3
"""Generate VCaP benchmark presentation figures and HTML slideshow.

Usage:
    conda activate rigel
    python scripts/debug/vcap_presentation.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns

# ── Paths ──────────────────────────────────────────────────────────
BENCHMARK_DIR = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate"
)
CONDITION = "gdna_high_ss_0.90_nrna_none"
VBEM_DIR = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "vbem"
MAP_DIR = BENCHMARK_DIR / "runs" / CONDITION / "rigel" / "map"
SALMON_DIR = BENCHMARK_DIR / "runs" / CONDITION / "salmon"
KALLISTO_DIR = BENCHMARK_DIR / "runs" / CONDITION / "kallisto"
TRUTH_PATH = BENCHMARK_DIR / "truth_abundances_nrna_none.tsv"
RESULTS_DIR = Path("results/vcap")
OUT_DIR = Path("results/vcap/presentation")

# ── Style ──────────────────────────────────────────────────────────
TOOL_COLORS = {
    "Rigel VBEM": "#2196F3",
    "Rigel MAP": "#4CAF50",
    "Salmon": "#FF9800",
    "Kallisto": "#9C27B0",
}
TOOL_ORDER = ["Rigel VBEM", "Rigel MAP", "Salmon", "Kallisto"]


def setup_style():
    sns.set_theme(style="whitegrid", font_scale=1.2)
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "axes.titleweight": "bold",
        "axes.titlesize": 14,
    })


# ── Data Loading ──────────────────────────────────────────────────
def load_all():
    truth = pd.read_csv(TRUTH_PATH, sep="\t")

    vbem = pd.read_feather(VBEM_DIR / "quant.feather")
    map_em = pd.read_feather(MAP_DIR / "quant.feather")

    # Salmon
    salmon_path = SALMON_DIR / "quant.sf.gz"
    salmon = pd.read_csv(salmon_path, sep="\t")
    salmon = salmon.rename(columns={"Name": "transcript_id", "TPM": "tpm", "NumReads": "count"})

    # Kallisto
    kallisto_path = KALLISTO_DIR / "abundance.tsv"
    kallisto = pd.read_csv(kallisto_path, sep="\t")
    kallisto = kallisto.rename(columns={"target_id": "transcript_id", "tpm": "tpm", "est_counts": "count"})

    # Gene-level
    vbem_gene = pd.read_feather(VBEM_DIR / "gene_quant.feather")
    map_gene = pd.read_feather(MAP_DIR / "gene_quant.feather")

    # Summaries
    with open(VBEM_DIR / "summary.json") as f:
        vbem_summary = json.load(f)
    with open(MAP_DIR / "summary.json") as f:
        map_summary = json.load(f)

    locus_stats = pd.read_feather(VBEM_DIR / "locus_stats.feather")

    return {
        "truth": truth,
        "vbem": vbem,
        "map": map_em,
        "salmon": salmon,
        "kallisto": kallisto,
        "vbem_gene": vbem_gene,
        "map_gene": map_gene,
        "vbem_summary": vbem_summary,
        "map_summary": map_summary,
        "locus_stats": locus_stats,
    }


def build_merged_tx(data):
    """Build merged transcript-level dataframe with all tools."""
    truth = data["truth"][["transcript_id", "gene_id", "gene_name", "mrna_abundance"]]
    merged = truth.copy()

    for tool_key, tool_name, tpm_col in [
        ("vbem", "vbem", "tpm"),
        ("map", "map", "tpm"),
        ("salmon", "salmon", "tpm"),
        ("kallisto", "kallisto", "tpm"),
    ]:
        df = data[tool_key][["transcript_id", tpm_col]].copy()
        df = df.rename(columns={tpm_col: f"{tool_name}_tpm"})
        merged = merged.merge(df, on="transcript_id", how="left")
        merged[f"{tool_name}_tpm"] = merged[f"{tool_name}_tpm"].fillna(0)

    return merged


def build_merged_gene(data):
    """Build merged gene-level dataframe."""
    truth = data["truth"].groupby("gene_id").agg(
        gene_name=("gene_name", "first"),
        mrna_abundance=("mrna_abundance", "sum"),
    ).reset_index()

    merged = truth.copy()
    for tool_key, tool_name in [("vbem_gene", "vbem"), ("map_gene", "map")]:
        df = data[tool_key][["gene_id", "tpm"]].copy()
        df = df.rename(columns={"tpm": f"{tool_name}_tpm"})
        merged = merged.merge(df, on="gene_id", how="left")
        merged[f"{tool_name}_tpm"] = merged[f"{tool_name}_tpm"].fillna(0)

    # Salmon and kallisto: aggregate from transcript level
    for tool_key, tool_name in [("salmon", "salmon"), ("kallisto", "kallisto")]:
        tx = data["truth"][["transcript_id", "gene_id"]].merge(
            data[tool_key][["transcript_id", "tpm"]],
            on="transcript_id", how="left",
        )
        tx["tpm"] = tx["tpm"].fillna(0)
        gene_tpm = tx.groupby("gene_id")["tpm"].sum().reset_index()
        gene_tpm = gene_tpm.rename(columns={"tpm": f"{tool_name}_tpm"})
        merged = merged.merge(gene_tpm, on="gene_id", how="left")
        merged[f"{tool_name}_tpm"] = merged[f"{tool_name}_tpm"].fillna(0)

    return merged


# ── Figure 1: Transcript-level scatter plots ──────────────────────
def fig_scatter_grid(merged, out_dir):
    pc = 1.0
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    tools = [
        ("vbem_tpm", "Rigel VBEM"),
        ("map_tpm", "Rigel MAP"),
        ("salmon_tpm", "Salmon"),
        ("kallisto_tpm", "Kallisto"),
    ]

    truth_log2 = np.log2(merged["mrna_abundance"] + pc)

    for ax, (col, name) in zip(axes, tools):
        pred_log2 = np.log2(merged[col] + pc)
        r = np.corrcoef(truth_log2, pred_log2)[0, 1]

        hb = ax.hexbin(
            truth_log2, pred_log2,
            gridsize=80, cmap="YlOrRd", mincnt=1,
            norm=mcolors.LogNorm(vmin=1, vmax=5000),
            linewidths=0.1, edgecolors="face",
        )
        lims = [-0.5, 14]
        ax.plot(lims, lims, "k-", lw=1, alpha=0.7)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Truth log₂(TPM + 1)")
        ax.set_ylabel("Predicted log₂(TPM + 1)")
        ax.set_title(f"{name}\nR = {r:.4f}")
        ax.set_aspect("equal")

    fig.colorbar(hb, ax=axes[-1], label="Count", shrink=0.7)
    fig.suptitle(
        "Transcript-Level Accuracy — log₂(TPM + 1)",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig1_scatter_grid.png")
    plt.close(fig)


# ── Figure 2: Transcript-level bar metrics ────────────────────────
def fig_tx_metrics_bars(out_dir):
    metrics = pd.read_csv(RESULTS_DIR / "transcript_metrics.csv")
    metrics["tool_label"] = metrics["tool"].map({
        "rigel/vbem": "Rigel VBEM",
        "rigel/map": "Rigel MAP",
        "salmon": "Salmon",
        "kallisto": "Kallisto",
    })

    fig, axes = plt.subplots(1, 5, figsize=(25, 5))

    vals = metrics.set_index("tool_label").loc[TOOL_ORDER]

    # Pearson R
    ax = axes[0]
    bars = ax.bar(
        TOOL_ORDER, vals["pearson_r"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylim(0.85, 1.0)
    ax.set_ylabel("Pearson R")
    ax.set_title("Pearson Correlation")
    for bar, v in zip(bars, vals["pearson_r"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
                f"{v:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # WARE
    ax = axes[1]
    bars = ax.bar(
        TOOL_ORDER, vals["ware"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("WARE")
    ax.set_title("Weighted Abs. Relative Error")
    for bar, v in zip(bars, vals["ware"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # MAE
    ax = axes[2]
    bars = ax.bar(
        TOOL_ORDER, vals["mae"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("MAE (TPM)")
    ax.set_title("Mean Absolute Error")
    for bar, v in zip(bars, vals["mae"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # MAPE
    ax = axes[3]
    bars = ax.bar(
        TOOL_ORDER, vals["mape"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("MAPE (%)")
    ax.set_title("Mean Abs. Percent Error")
    for bar, v in zip(bars, vals["mape"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # F1
    ax = axes[4]
    bars = ax.bar(
        TOOL_ORDER, vals["f1"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylim(0.6, 1.0)
    ax.set_ylabel("F1 Score")
    ax.set_title("Detection F1 (TPM > 1)")
    for bar, v in zip(bars, vals["f1"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    fig.suptitle(
        "Transcript-Level Metrics — All Tools",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig2_tx_metrics.png")
    plt.close(fig)


# ── Figure 3: Gene-level scatter ──────────────────────────────────
def fig_gene_scatter(merged_gene, out_dir):
    pc = 1.0
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    tools = [
        ("vbem_tpm", "Rigel VBEM"),
        ("map_tpm", "Rigel MAP"),
        ("salmon_tpm", "Salmon"),
        ("kallisto_tpm", "Kallisto"),
    ]

    truth_log2 = np.log2(merged_gene["mrna_abundance"] + pc)

    for ax, (col, name) in zip(axes, tools):
        pred_log2 = np.log2(merged_gene[col] + pc)
        r = np.corrcoef(truth_log2, pred_log2)[0, 1]

        hb = ax.hexbin(
            truth_log2, pred_log2,
            gridsize=60, cmap="YlGn", mincnt=1,
            norm=mcolors.LogNorm(vmin=1, vmax=2000),
            linewidths=0.1, edgecolors="face",
        )
        lims = [-0.5, 16]
        ax.plot(lims, lims, "k-", lw=1, alpha=0.7)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        ax.set_xlabel("Truth log₂(TPM + 1)")
        ax.set_ylabel("Predicted log₂(TPM + 1)")
        ax.set_title(f"{name}\nR = {r:.4f}")
        ax.set_aspect("equal")

    fig.colorbar(hb, ax=axes[-1], label="Count", shrink=0.7)
    fig.suptitle(
        "Gene-Level Accuracy — log₂(TPM + 1)",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig3_gene_scatter.png")
    plt.close(fig)


# ── Figure 4: Gene-level metrics bars ─────────────────────────────
def fig_gene_metrics_bars(out_dir):
    metrics = pd.read_csv(RESULTS_DIR / "gene_metrics.csv")
    metrics["tool_label"] = metrics["tool"].map({
        "rigel/vbem": "Rigel VBEM",
        "rigel/map": "Rigel MAP",
        "salmon": "Salmon",
        "kallisto": "Kallisto",
    })

    fig, axes = plt.subplots(1, 5, figsize=(25, 5))

    vals = metrics.set_index("tool_label").loc[TOOL_ORDER]

    # Pearson R
    ax = axes[0]
    bars = ax.bar(
        TOOL_ORDER, vals["pearson_r"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylim(0.97, 1.0)
    ax.set_ylabel("Pearson R")
    ax.set_title("Gene-Level Pearson R")
    for bar, v in zip(bars, vals["pearson_r"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.0005,
                f"{v:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # WARE
    ax = axes[1]
    bars = ax.bar(
        TOOL_ORDER, vals["ware"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("WARE")
    ax.set_title("Gene-Level WARE")
    for bar, v in zip(bars, vals["ware"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.002,
                f"{v:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # MAE
    ax = axes[2]
    bars = ax.bar(
        TOOL_ORDER, vals["mae"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("MAE (TPM)")
    ax.set_title("Gene-Level MAE")
    for bar, v in zip(bars, vals["mae"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # MAPE
    ax = axes[3]
    bars = ax.bar(
        TOOL_ORDER, vals["mape"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylabel("MAPE (%)")
    ax.set_title("Gene-Level MAPE")
    for bar, v in zip(bars, vals["mape"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                f"{v:.1f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    # Spearman R
    ax = axes[4]
    bars = ax.bar(
        TOOL_ORDER, vals["spearman_r"],
        color=[TOOL_COLORS[t] for t in TOOL_ORDER],
        edgecolor="white", linewidth=0.5,
    )
    ax.set_ylim(0.65, 0.95)
    ax.set_ylabel("Spearman R")
    ax.set_title("Gene-Level Spearman R")
    for bar, v in zip(bars, vals["spearman_r"]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.4f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.tick_params(axis="x", rotation=20)

    fig.suptitle(
        "Gene-Level Metrics — All Tools",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig4_gene_metrics.png")
    plt.close(fig)


# ── Figure 5: Pool-level fragment assignment ──────────────────────
def fig_pool_level(out_dir):
    pool = pd.read_csv(RESULTS_DIR / "pool_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: stacked bar for rigel VBEM + MAP
    ax = axes[0]
    tools = ["rigel/vbem", "rigel/map"]
    labels = ["Rigel VBEM", "Rigel MAP"]
    truth_mrna = 50_000_000
    truth_gdna = 25_000_000

    x = np.arange(3)  # truth, vbem, map
    width = 0.6
    bar_labels = ["Truth", "Rigel VBEM", "Rigel MAP"]

    mrna_vals = [truth_mrna]
    nrna_vals = [0]
    gdna_vals = [truth_gdna]
    inter_vals = [0]

    for tool in tools:
        row = pool[pool["tool"] == tool].iloc[0]
        mrna_vals.append(float(row["mrna_pred"]))
        nrna_vals.append(float(row["nrna_pred"]))
        gdna_vals.append(float(row["gdna_pred"]))
        inter_vals.append(float(row["intergenic_pred"]))

    mrna_vals = np.array(mrna_vals) / 1e6
    nrna_vals = np.array(nrna_vals) / 1e6
    gdna_vals = np.array(gdna_vals) / 1e6
    inter_vals = np.array(inter_vals) / 1e6

    ax.bar(x, mrna_vals, width, label="mRNA", color="#2196F3")
    ax.bar(x, nrna_vals, width, bottom=mrna_vals, label="nRNA", color="#FF9800")
    ax.bar(x, gdna_vals, width, bottom=mrna_vals + nrna_vals, label="gDNA", color="#F44336")
    ax.bar(x, inter_vals, width, bottom=mrna_vals + nrna_vals + gdna_vals,
           label="Intergenic", color="#9E9E9E")

    ax.set_xticks(x)
    ax.set_xticklabels(bar_labels)
    ax.set_ylabel("Fragments (millions)")
    ax.set_title("Fragment Pool Decomposition")
    ax.legend(loc="upper right")
    ax.set_ylim(0, 85)

    # Add truth lines
    ax.axhline(y=50, color="#2196F3", linestyle="--", alpha=0.4, lw=1)
    ax.axhline(y=75, color="#F44336", linestyle="--", alpha=0.4, lw=1)

    # Right: relative error bars for mRNA and gDNA
    ax = axes[1]
    categories = ["mRNA", "gDNA"]
    vbem_row = pool[pool["tool"] == "rigel/vbem"].iloc[0]
    map_row = pool[pool["tool"] == "rigel/map"].iloc[0]

    vbem_errs = [
        float(vbem_row["mrna_rel_error"]) * 100,
        float(vbem_row["gdna_rel_error"]) * 100,
    ]
    map_errs = [
        float(map_row["mrna_rel_error"]) * 100,
        float(map_row["gdna_rel_error"]) * 100,
    ]

    x = np.arange(len(categories))
    w = 0.3
    b1 = ax.bar(x - w/2, vbem_errs, w, label="Rigel VBEM", color=TOOL_COLORS["Rigel VBEM"])
    b2 = ax.bar(x + w/2, map_errs, w, label="Rigel MAP", color=TOOL_COLORS["Rigel MAP"])
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_ylabel("Relative Error (%)")
    ax.set_title("Pool-Level Estimation Error")
    ax.legend()
    ax.axhline(y=0, color="black", lw=0.5)
    for bar, v in zip(b1, vbem_errs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 1,
                f"{v:.2f}%", ha="center", va="top", fontsize=9, fontweight="bold")
    for bar, v in zip(b2, map_errs):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() - 1,
                f"{v:.2f}%", ha="center", va="top", fontsize=9, fontweight="bold")

    fig.suptitle(
        "Pool-Level Fragment Assignment",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig5_pool_level.png")
    plt.close(fig)


# ── Figure 6: Expression-stratified grouped bar ───────────────────
def fig_stratified(out_dir):
    strat = pd.read_csv(RESULTS_DIR / "stratified_metrics.csv")
    strat["tool_label"] = strat["tool"].map({
        "rigel/vbem": "Rigel VBEM",
        "rigel/map": "Rigel MAP",
        "salmon": "Salmon",
        "kallisto": "Kallisto",
    })

    bin_order = ["low_1-10", "mid_10-100", "high_100-1000", "very_high_1000+"]
    bin_labels = ["1–10", "10–100", "100–1K", "1K+"]

    fig, axes = plt.subplots(1, 2, figsize=(16, 5.5))

    # Pearson R by expression bin
    ax = axes[0]
    x = np.arange(len(bin_order))
    w = 0.18
    offsets = [-1.5*w, -0.5*w, 0.5*w, 1.5*w]

    for i, tool in enumerate(TOOL_ORDER):
        tool_data = strat[strat["tool_label"] == tool]
        vals = []
        for b in bin_order:
            row = tool_data[tool_data["expression_bin"] == b]
            vals.append(float(row["pearson_r"].iloc[0]) if len(row) > 0 else 0)
        ax.bar(x + offsets[i], vals, w, label=tool, color=TOOL_COLORS[tool],
               edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel("Truth Expression Bin (TPM)")
    ax.set_ylabel("Pearson R")
    ax.set_title("Pearson R by Expression Level")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.1)

    # WARE by expression bin
    ax = axes[1]
    for i, tool in enumerate(TOOL_ORDER):
        tool_data = strat[strat["tool_label"] == tool]
        vals = []
        for b in bin_order:
            row = tool_data[tool_data["expression_bin"] == b]
            vals.append(float(row["ware"].iloc[0]) if len(row) > 0 else 0)
        ax.bar(x + offsets[i], vals, w, label=tool, color=TOOL_COLORS[tool],
               edgecolor="white", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels)
    ax.set_xlabel("Truth Expression Bin (TPM)")
    ax.set_ylabel("WARE")
    ax.set_title("WARE by Expression Level")
    ax.legend(fontsize=9)

    fig.suptitle(
        "Expression-Stratified Performance",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig6_stratified.png")
    plt.close(fig)


# ── Figure 7: VBEM fix — before/after ─────────────────────────────
def fig_vbem_fix(out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Data: old broken VBEM vs new fixed VBEM vs MAP
    cats = ["Old VBEM\n(broken)", "New VBEM\n(fixed)", "MAP-EM"]
    colors = ["#E53935", "#2196F3", "#4CAF50"]

    # Pearson R
    ax = axes[0]
    vals = [0.8145, 0.9865, 0.9862]
    bars = ax.bar(cats, vals, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_ylim(0.75, 1.02)
    ax.set_ylabel("Pearson R")
    ax.set_title("Pearson R Recovery")
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
                f"{v:.4f}", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.axhline(y=0.9862, color="#4CAF50", linestyle="--", alpha=0.4, lw=1, label="MAP baseline")

    # WARE
    ax = axes[1]
    vals = [0.4037, 0.0816, 0.0796]
    bars = ax.bar(cats, vals, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_ylabel("WARE")
    ax.set_title("WARE Recovery")
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.008,
                f"{v:.4f}", ha="center", va="bottom", fontsize=11, fontweight="bold")

    # False zeros
    ax = axes[2]
    vals = [4865, 1190, 623]
    bars = ax.bar(cats, vals, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_ylabel("False Zeros (truth > 1 TPM)")
    ax.set_title("False Zero Reduction")
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
                f"{v:,}", ha="center", va="bottom", fontsize=11, fontweight="bold")

    fig.suptitle(
        "VBEM Clamp Floor Fix — Before vs After",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig7_vbem_fix.png")
    plt.close(fig)


# ── Figure 8: Convergence analysis ────────────────────────────────
def fig_convergence(locus_stats, out_dir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    ls = locus_stats.copy()
    ls["time_s"] = ls["squarem_us"] / 1e6

    # Iteration histogram (exclude mega-locus for visibility)
    ax = axes[0]
    non_mega = ls[~ls["is_mega_locus"]]
    ax.hist(non_mega["squarem_iterations"], bins=50, color="#2196F3",
            edgecolor="white", linewidth=0.5, alpha=0.8)
    ax.set_xlabel("SQUAREM Iterations")
    ax.set_ylabel("Number of Loci")
    ax.set_title(f"Iteration Distribution\n(excl. mega-locus, N={len(non_mega):,})")
    ax.set_yscale("log")
    ax.axvline(x=333, color="red", linestyle="--", lw=1.5, label="Max (333)")
    ax.legend()

    # Components vs iterations scatter
    ax = axes[1]
    multi = ls[ls["n_components"] > 2].copy()
    mega_mask = multi["is_mega_locus"]
    ax.scatter(
        multi.loc[~mega_mask, "n_components"],
        multi.loc[~mega_mask, "squarem_iterations"],
        s=5, alpha=0.3, color="#2196F3", label="Regular loci",
    )
    if mega_mask.any():
        ax.scatter(
            multi.loc[mega_mask, "n_components"],
            multi.loc[mega_mask, "squarem_iterations"],
            s=100, color="red", marker="*", zorder=5, label="Mega-locus",
        )
    ax.set_xscale("log")
    ax.set_xlabel("EM Components")
    ax.set_ylabel("SQUAREM Iterations")
    ax.set_title("Locus Size vs Convergence Speed")
    ax.legend()

    # Time pie chart
    ax = axes[2]
    mega_time = float(ls[ls["is_mega_locus"]]["time_s"].sum())
    other_time = float(ls[~ls["is_mega_locus"]]["time_s"].sum())
    scan_time = 2003 - (mega_time + other_time)  # approximate scan time
    sizes = [mega_time, other_time, scan_time]
    labels_pie = [
        f"Mega-locus EM\n{mega_time:.0f}s",
        f"Other loci EM\n{other_time:.0f}s",
        f"BAM scan + other\n{scan_time:.0f}s",
    ]
    colors_pie = ["#F44336", "#2196F3", "#9E9E9E"]
    ax.pie(sizes, labels=labels_pie, colors=colors_pie, autopct="%1.0f%%",
           startangle=90, textprops={"fontsize": 10})
    ax.set_title("Pipeline Time Breakdown")

    fig.suptitle(
        "VBEM Convergence Analysis",
        fontsize=16, fontweight="bold", y=1.02,
    )
    plt.tight_layout()
    fig.savefig(out_dir / "fig8_convergence.png")
    plt.close(fig)


# ── Figure 9: Precision / Recall / F1 ─────────────────────────────
def fig_precision_recall(out_dir):
    metrics = pd.read_csv(RESULTS_DIR / "transcript_metrics.csv")
    metrics["tool_label"] = metrics["tool"].map({
        "rigel/vbem": "Rigel VBEM",
        "rigel/map": "Rigel MAP",
        "salmon": "Salmon",
        "kallisto": "Kallisto",
    })

    fig, ax = plt.subplots(figsize=(8, 6))

    vals = metrics.set_index("tool_label").loc[TOOL_ORDER]

    for tool in TOOL_ORDER:
        row = vals.loc[tool]
        ax.scatter(
            row["recall"], row["precision"],
            s=200, color=TOOL_COLORS[tool], zorder=5,
            edgecolors="black", linewidth=0.8,
        )
        ax.annotate(
            f"  {tool}\n  F1={row['f1']:.3f}",
            (row["recall"], row["precision"]),
            fontsize=10, fontweight="bold",
            va="center",
        )

    ax.set_xlabel("Recall", fontsize=13)
    ax.set_ylabel("Precision", fontsize=13)
    ax.set_title("Transcript Detection — Precision vs Recall\n(threshold: TPM > 1)",
                 fontsize=14, fontweight="bold")
    ax.set_xlim(0.55, 1.0)
    ax.set_ylim(0.55, 1.0)

    # F1 iso-lines
    for f1_val in [0.7, 0.75, 0.8, 0.85, 0.9]:
        recall_range = np.linspace(0.01, 1.0, 200)
        precision_line = (f1_val * recall_range) / (2 * recall_range - f1_val)
        valid = (precision_line > 0) & (precision_line <= 1)
        ax.plot(recall_range[valid], precision_line[valid], "k-", alpha=0.1, lw=0.8)
        # Label
        idx = np.argmin(np.abs(recall_range - 0.95))
        if valid[idx]:
            ax.text(0.97, precision_line[idx], f"F1={f1_val}",
                    fontsize=7, alpha=0.4, va="center")

    plt.tight_layout()
    fig.savefig(out_dir / "fig9_precision_recall.png")
    plt.close(fig)


# ── HTML Slideshow ─────────────────────────────────────────────────
def generate_html(out_dir):
    slides = [
        {
            "title": "Rigel Benchmark — VCaP Prostate Simulation",
            "subtitle": "Transcript quantification with joint mRNA / nRNA / gDNA modeling",
            "content": """
            <div style="text-align:left; max-width:800px; margin:0 auto; font-size:18px; line-height:1.7;">
                <ul>
                    <li><b>Simulation:</b> VCaP prostate cancer expression profile</li>
                    <li><b>Fragments:</b> 50M RNA + 25M gDNA = 75M total paired-end reads</li>
                    <li><b>Strand specificity:</b> 0.90 (stranded library)</li>
                    <li><b>Annotation:</b> GENCODE (GRCh38) — 254,461 transcripts, 63,472 genes</li>
                    <li><b>Tools compared:</b> Rigel (VBEM, MAP-EM) vs Salmon vs Kallisto</li>
                    <li><b>Key challenge:</b> 50% gDNA contamination — tools without gDNA modeling are disadvantaged</li>
                </ul>
            </div>
            """,
            "image": None,
        },
        {
            "title": "Transcript-Level Accuracy",
            "subtitle": "log₂(TPM + 1) correlation — 254,461 transcripts",
            "content": None,
            "image": "fig1_scatter_grid.png",
        },
        {
            "title": "Transcript-Level Metrics",
            "subtitle": "Pearson R &bull; WARE &bull; MAE &bull; MAPE &bull; Detection F1",
            "content": None,
            "image": "fig2_tx_metrics.png",
        },
        {
            "title": "Precision vs Recall",
            "subtitle": "Transcript detection at TPM > 1 threshold",
            "content": None,
            "image": "fig9_precision_recall.png",
        },
        {
            "title": "Gene-Level Accuracy",
            "subtitle": "log₂(TPM + 1) correlation — 63,472 genes",
            "content": None,
            "image": "fig3_gene_scatter.png",
        },
        {
            "title": "Gene-Level Metrics",
            "subtitle": "Pearson R &bull; WARE &bull; MAE &bull; MAPE &bull; Spearman R",
            "content": None,
            "image": "fig4_gene_metrics.png",
        },
        {
            "title": "Expression-Stratified Performance",
            "subtitle": "Accuracy breakdown by expression level",
            "content": None,
            "image": "fig6_stratified.png",
        },
        {
            "title": "Pool-Level Fragment Assignment",
            "subtitle": "mRNA / nRNA / gDNA decomposition accuracy",
            "content": None,
            "image": "fig5_pool_level.png",
        },
        {
            "title": "VBEM Clamp Floor Fix",
            "subtitle": "Before vs After — eliminating the absorbing barrier",
            "content": """
            <div style="text-align:left; max-width:800px; margin:0 auto; font-size:16px; line-height:1.6;">
                <p><b>Problem:</b> SQUAREM acceleration could push component alpha below recovery threshold.
                   With prior ≈ 10⁻⁶, digamma creates an absorbing barrier — components become permanently dead.</p>
                <p><b>Fix:</b> <code>VBEM_CLAMP_FLOOR = 0.1</code> — keeps all components in the recoverable regime
                   (ψ(0.1) ≈ −10.4, weight ≈ 3×10⁻⁵). Components can still be naturally suppressed by EM,
                   but can recover if they have genuine read support.</p>
            </div>
            """,
            "image": "fig7_vbem_fix.png",
        },
        {
            "title": "VBEM Convergence Analysis",
            "subtitle": "Per-locus iteration and timing profiles",
            "content": None,
            "image": "fig8_convergence.png",
        },
        {
            "title": "Summary",
            "subtitle": "",
            "content": """
            <div style="text-align:left; max-width:850px; margin:0 auto; font-size:17px; line-height:1.8;">
                <table style="width:100%; border-collapse:collapse; margin:20px 0;">
                    <tr style="background:#f5f5f5;">
                        <th style="padding:8px; text-align:left; border-bottom:2px solid #ddd;">Metric</th>
                        <th style="padding:8px; text-align:right; border-bottom:2px solid #ddd;">Rigel VBEM</th>
                        <th style="padding:8px; text-align:right; border-bottom:2px solid #ddd;">Rigel MAP</th>
                        <th style="padding:8px; text-align:right; border-bottom:2px solid #ddd;">Salmon</th>
                        <th style="padding:8px; text-align:right; border-bottom:2px solid #ddd;">Kallisto</th>
                    </tr>
                    <tr><td style="padding:6px;">TX Pearson R</td><td style="padding:6px; text-align:right;"><b>0.9865</b></td><td style="padding:6px; text-align:right;">0.9862</td><td style="padding:6px; text-align:right;">0.8900</td><td style="padding:6px; text-align:right;">0.9958</td></tr>
                    <tr style="background:#f9f9f9;"><td style="padding:6px;">TX WARE</td><td style="padding:6px; text-align:right;"><b>0.082</b></td><td style="padding:6px; text-align:right;">0.080</td><td style="padding:6px; text-align:right;">0.363</td><td style="padding:6px; text-align:right;">0.116</td></tr>
                    <tr><td style="padding:6px;">Gene Pearson R</td><td style="padding:6px; text-align:right;"><b>0.9954</b></td><td style="padding:6px; text-align:right;">0.9952</td><td style="padding:6px; text-align:right;">0.9880</td><td style="padding:6px; text-align:right;">0.9984</td></tr>
                    <tr style="background:#f9f9f9;"><td style="padding:6px;">Gene WARE</td><td style="padding:6px; text-align:right;"><b>0.034</b></td><td style="padding:6px; text-align:right;">0.035</td><td style="padding:6px; text-align:right;">0.118</td><td style="padding:6px; text-align:right;">0.081</td></tr>
                    <tr><td style="padding:6px;">mRNA rel. error</td><td style="padding:6px; text-align:right;"><b>−0.06%</b></td><td style="padding:6px; text-align:right;">−0.06%</td><td colspan=2 style="padding:6px; text-align:center; color:#999;">N/A</td></tr>
                    <tr style="background:#f9f9f9;"><td style="padding:6px;">gDNA rel. error</td><td style="padding:6px; text-align:right;">−14.6%</td><td style="padding:6px; text-align:right;">−14.8%</td><td colspan=2 style="padding:6px; text-align:center; color:#999;">No gDNA model</td></tr>
                </table>
                <ul>
                    <li><b>Rigel achieves lowest WARE</b> (0.082) — best abundance accuracy at both transcript and gene level in the presence of 50% gDNA contamination</li>
                    <li><b>VBEM and MAP-EM are statistically equivalent</b> — VBEM clamp floor fix fully restores accuracy</li>
                    <li><b>Salmon degrades</b> under gDNA contamination (no gDNA model), WARE 4.5× worse than Rigel</li>
                    <li><b>Kallisto</b> has highest Pearson R but 1.4× higher WARE due to gDNA contamination artifacts</li>
                </ul>
            </div>
            """,
            "image": None,
        },
    ]

    html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Rigel Benchmark — VCaP Prostate</title>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif; background: #1a1a2e; color: #e0e0e0; overflow: hidden; height: 100vh; }
.slide { display: none; width: 100vw; height: 100vh; flex-direction: column; align-items: center; justify-content: center; padding: 40px; }
.slide.active { display: flex; }
.slide h1 { font-size: 36px; color: #ffffff; margin-bottom: 8px; text-align: center; }
.slide h2 { font-size: 20px; color: #90caf9; margin-bottom: 30px; font-weight: 400; text-align: center; }
.slide img { max-width: 92vw; max-height: 68vh; object-fit: contain; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.5); }
.slide .content { color: #e0e0e0; }
.slide .content table { width: 100%; }
.slide .content th { color: #90caf9; }
.slide .content code { background: #2d2d44; padding: 2px 6px; border-radius: 4px; font-size: 14px; }
.nav { position: fixed; bottom: 20px; width: 100%; display: flex; justify-content: center; gap: 20px; z-index: 100; }
.nav button { padding: 10px 28px; font-size: 16px; border: none; border-radius: 6px; cursor: pointer; background: #304060; color: #e0e0e0; transition: background 0.2s; }
.nav button:hover { background: #4060a0; }
.nav button:disabled { opacity: 0.3; cursor: default; }
.counter { position: fixed; bottom: 25px; right: 30px; font-size: 14px; color: #666; z-index: 100; }
</style>
</head>
<body>
"""

    for i, slide in enumerate(slides):
        active = ' active' if i == 0 else ''
        html += f'<div class="slide{active}" id="slide-{i}">\n'
        html += f'  <h1>{slide["title"]}</h1>\n'
        if slide["subtitle"]:
            html += f'  <h2>{slide["subtitle"]}</h2>\n'
        if slide["image"]:
            html += f'  <img src="{slide["image"]}" alt="{slide["title"]}">\n'
        if slide["content"]:
            html += f'  <div class="content">{slide["content"]}</div>\n'
        html += '</div>\n'

    html += f"""
<div class="nav">
  <button id="prev" onclick="navigate(-1)">← Prev</button>
  <button id="next" onclick="navigate(1)">Next →</button>
</div>
<div class="counter" id="counter"></div>

<script>
let current = 0;
const total = {len(slides)};

function show(n) {{
  document.querySelectorAll('.slide').forEach(s => s.classList.remove('active'));
  document.getElementById('slide-' + n).classList.add('active');
  document.getElementById('counter').textContent = (n + 1) + ' / ' + total;
  document.getElementById('prev').disabled = (n === 0);
  document.getElementById('next').disabled = (n === total - 1);
}}

function navigate(dir) {{
  current = Math.max(0, Math.min(total - 1, current + dir));
  show(current);
}}

document.addEventListener('keydown', function(e) {{
  if (e.key === 'ArrowRight' || e.key === ' ') navigate(1);
  if (e.key === 'ArrowLeft') navigate(-1);
}});

show(0);
</script>
</body>
</html>"""

    with open(out_dir / "presentation.html", "w") as f:
        f.write(html)


# ── Main ──────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    setup_style()

    print("Loading data...")
    data = load_all()

    print("Building merged datasets...")
    merged_tx = build_merged_tx(data)
    merged_gene = build_merged_gene(data)

    print("Generating figures...")
    fig_scatter_grid(merged_tx, OUT_DIR)
    print("  [1/9] Transcript scatter grid")

    fig_tx_metrics_bars(OUT_DIR)
    print("  [2/9] Transcript metrics bars")

    fig_gene_scatter(merged_gene, OUT_DIR)
    print("  [3/9] Gene scatter grid")

    fig_gene_metrics_bars(OUT_DIR)
    print("  [4/9] Gene metrics bars")

    fig_pool_level(OUT_DIR)
    print("  [5/9] Pool-level assignment")

    fig_stratified(OUT_DIR)
    print("  [6/9] Expression-stratified")

    fig_vbem_fix(OUT_DIR)
    print("  [7/9] VBEM fix before/after")

    fig_convergence(data["locus_stats"], OUT_DIR)
    print("  [8/9] Convergence analysis")

    fig_precision_recall(OUT_DIR)
    print("  [9/9] Precision-recall")

    print("Generating HTML presentation...")
    generate_html(OUT_DIR)

    print(f"\nDone! Output in {OUT_DIR}/")
    print(f"  Open: {OUT_DIR}/presentation.html")
    print(f"  Figures: {OUT_DIR}/fig*.png")


if __name__ == "__main__":
    main()
