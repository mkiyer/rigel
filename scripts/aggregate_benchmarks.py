#!/usr/bin/env python3
"""Aggregate multi-seed combinatorial benchmark results into unified reports.

Reads summary.json from each seed_*/  subdirectory (produced by
benchmark_region_competition.py), computes per-condition and overall
statistics (mean, median, std across seeds), and writes:

- aggregate_summary.json — structured aggregate results
- aggregate_summary.csv  — tidy long-form CSV
- aggregate_summary.md   — human-readable Markdown report
- aggregate_per_tx.csv   — per-transcript detail across seeds & conditions

Usage:
    PYTHONPATH=src conda run -n hulkrna python scripts/aggregate_benchmarks.py \\
        --input-dir /path/to/bench_output \\
        --output-dir /path/to/bench_output
"""
from __future__ import annotations

import argparse
import json
import logging
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

TRANSCRIPT_TOOLS = ("hulkrna", "hulkrna_mm", "salmon", "kallisto")
GENE_TOOLS = ("hulkrna", "hulkrna_mm", "salmon", "kallisto", "htseq")
METRICS = ("mae", "rmse", "pearson", "spearman")


# ── Loading ──────────────────────────────────────────────────────────


def load_seed_results(input_dir: Path) -> list[tuple[int, list[dict]]]:
    """Load summary.json from each seed_* subdirectory."""
    results = []
    for seed_dir in sorted(input_dir.glob("seed_*")):
        sj = seed_dir / "summary.json"
        if not sj.exists():
            logger.warning("No summary.json in %s", seed_dir)
            continue
        seed = int(seed_dir.name.split("_")[-1])
        data = json.loads(sj.read_text())
        if data:
            results.append((seed, data))
            logger.info("Loaded seed %d: %d condition results", seed, len(data))
    return results


def load_per_tx_results(input_dir: Path) -> pd.DataFrame:
    """Load per_transcript_counts.csv from all seed/region/condition dirs."""
    parts: list[pd.DataFrame] = []
    for seed_dir in sorted(input_dir.glob("seed_*")):
        seed = int(seed_dir.name.split("_")[-1])
        for ptx in seed_dir.rglob("per_transcript_counts.csv"):
            # Path: seed_X/region/condition/per_transcript_counts.csv
            condition = ptx.parent.name
            region = ptx.parent.parent.name
            df = pd.read_csv(ptx)
            df["seed"] = seed
            df["region"] = region
            df["condition"] = condition
            # Parse condition dir name: gdna_<label>_ss_<value>
            cond_parts = condition.split("_ss_")
            if len(cond_parts) == 2:
                df["gdna_label"] = cond_parts[0].replace("gdna_", "")
                df["strand_specificity"] = float(cond_parts[1])
            parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def load_per_gene_results(input_dir: Path) -> pd.DataFrame:
    """Load per_gene_counts.csv from all seed/region/condition dirs."""
    parts: list[pd.DataFrame] = []
    for seed_dir in sorted(input_dir.glob("seed_*")):
        seed = int(seed_dir.name.split("_")[-1])
        for pgene in seed_dir.rglob("per_gene_counts.csv"):
            condition = pgene.parent.name
            region = pgene.parent.parent.name
            df = pd.read_csv(pgene)
            df["seed"] = seed
            df["region"] = region
            df["condition"] = condition
            cond_parts = condition.split("_ss_")
            if len(cond_parts) == 2:
                df["gdna_label"] = cond_parts[0].replace("gdna_", "")
                df["strand_specificity"] = float(cond_parts[1])
            parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


# ── Aggregation ──────────────────────────────────────────────────────


def aggregate(
    seed_results: list[tuple[int, list[dict]]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build tidy DataFrames for transcript-level and gene-level metrics.

    Returns (tx_metrics_df, gene_metrics_df), each with columns:
        seed, region, gdna_label, gdna_abundance, strand_specificity,
        tool, mae, rmse, pearson, spearman, elapsed_sec
    """
    tx_rows: list[dict] = []
    gene_rows: list[dict] = []

    for seed, conditions in seed_results:
        for r in conditions:
            base = {
                "seed": seed,
                "region": r["region"],
                "n_transcripts": r["n_transcripts"],
                "n_genes": r["n_genes"],
                "n_fragments": r["n_fragments"],
                "n_gdna_actual": r["n_gdna_actual"],
                "gdna_label": r["gdna_label"],
                "gdna_abundance": r["gdna_abundance"],
                "strand_specificity": r["strand_specificity"],
            }
            for tool, tm in r.get("transcript_metrics", {}).items():
                tx_rows.append({
                    **base,
                    "tool": tool,
                    "mae": tm["mean_abs_error"],
                    "rmse": tm["rmse"],
                    "pearson": tm["pearson"],
                    "spearman": tm["spearman"],
                    "elapsed_sec": tm["elapsed_sec"],
                    "total_abs_error": tm["total_abs_error"],
                })
            for tool, tm in r.get("gene_metrics", {}).items():
                gene_rows.append({
                    **base,
                    "tool": tool,
                    "mae": tm["mean_abs_error"],
                    "rmse": tm["rmse"],
                    "pearson": tm["pearson"],
                    "spearman": tm["spearman"],
                    "elapsed_sec": tm["elapsed_sec"],
                    "total_abs_error": tm["total_abs_error"],
                })

    tx_df = pd.DataFrame(tx_rows) if tx_rows else pd.DataFrame()
    gene_df = pd.DataFrame(gene_rows) if gene_rows else pd.DataFrame()
    return tx_df, gene_df


# ── Report writing ───────────────────────────────────────────────────


def _fmt(val: float, decimals: int = 3) -> str:
    if math.isnan(val):
        return "—"
    return f"{val:.{decimals}f}"


def _agg_table(
    df: pd.DataFrame,
    group_cols: list[str],
    tools: tuple[str, ...],
    metric: str = "mae",
    decimals: int = 3,
) -> list[str]:
    """Build a Markdown table of mean metric values grouped by group_cols.

    Columns: *group_cols | tool1 | tool2 | ...
    """
    if df.empty:
        return ["(no data)"]

    header_parts = [col.replace("_", " ").title() for col in group_cols] + list(tools)
    header = "| " + " | ".join(header_parts) + " |"
    sep_parts = ["---"] * len(group_cols) + ["---:"] * len(tools)
    sep = "| " + " | ".join(sep_parts) + " |"
    lines = [header, sep]

    # Get unique group combos
    groups = df[group_cols].drop_duplicates().sort_values(group_cols)
    for _, group_row in groups.iterrows():
        mask = pd.Series(True, index=df.index)
        for col in group_cols:
            mask &= df[col] == group_row[col]
        sub = df[mask]

        row_parts = [str(group_row[col]) for col in group_cols]
        for tool in tools:
            tsub = sub[sub["tool"] == tool]
            if len(tsub) == 0:
                row_parts.append("—")
            else:
                val = tsub[metric].mean()
                row_parts.append(_fmt(val, decimals))
        lines.append("| " + " | ".join(row_parts) + " |")

    return lines


def write_aggregate_report(
    tx_df: pd.DataFrame,
    gene_df: pd.DataFrame,
    per_tx: pd.DataFrame,
    output_dir: Path,
    n_seeds: int,
    include_htseq: bool,
) -> None:
    """Write aggregate JSON, CSV, and Markdown reports."""
    output_dir.mkdir(parents=True, exist_ok=True)

    tx_tools = TRANSCRIPT_TOOLS
    gene_tools = GENE_TOOLS if include_htseq else TRANSCRIPT_TOOLS

    # ── Aggregate JSON ──────────────────────────────────────────────

    agg_json: dict = {"n_seeds": n_seeds}

    if not tx_df.empty:
        overall_tx = (
            tx_df.groupby("tool")[list(METRICS)]
            .agg(["mean", "std"])
            .to_dict()
        )
        agg_json["transcript_level_overall"] = {
            tool: {
                m: {
                    "mean": float(tx_df[tx_df["tool"] == tool][m].mean()),
                    "std": float(tx_df[tx_df["tool"] == tool][m].std()),
                }
                for m in METRICS
            }
            for tool in tx_tools
            if tool in tx_df["tool"].values
        }

    if not gene_df.empty:
        agg_json["gene_level_overall"] = {
            tool: {
                m: {
                    "mean": float(gene_df[gene_df["tool"] == tool][m].mean()),
                    "std": float(gene_df[gene_df["tool"] == tool][m].std()),
                }
                for m in METRICS
            }
            for tool in gene_tools
            if tool in gene_df["tool"].values
        }

    with open(output_dir / "aggregate_summary.json", "w") as f:
        json.dump(agg_json, f, indent=2)

    # ── Aggregate CSV ───────────────────────────────────────────────

    if not tx_df.empty:
        tx_df.to_csv(output_dir / "aggregate_transcript_metrics.csv", index=False)
    if not gene_df.empty:
        gene_df.to_csv(output_dir / "aggregate_gene_metrics.csv", index=False)
    if not per_tx.empty:
        per_tx.to_csv(output_dir / "aggregate_per_tx.csv", index=False)

    # ── Aggregate Markdown ──────────────────────────────────────────

    lines = [
        f"# Aggregate Benchmark Summary ({n_seeds} seeds)",
        "",
    ]

    # Overall transcript-level
    lines.extend([
        "## Overall Transcript-Level Metrics (mean across all regions, conditions, seeds)",
        "",
    ])
    if not tx_df.empty:
        overall = tx_df.groupby("tool")[list(METRICS)].mean().sort_values("mae")
        lines.extend([
            "| Tool | MAE | RMSE | Pearson | Spearman |",
            "| --- | ---: | ---: | ---: | ---: |",
        ])
        for tool, row in overall.iterrows():
            if tool in tx_tools:
                lines.append(
                    f"| {tool} | {_fmt(row['mae'])} | {_fmt(row['rmse'])} | "
                    f"{_fmt(row['pearson'], 4)} | {_fmt(row['spearman'], 4)} |"
                )
    lines.append("")

    # Overall gene-level
    lines.extend([
        "## Overall Gene-Level Metrics (mean across all regions, conditions, seeds)",
        "",
    ])
    if not gene_df.empty:
        overall_g = gene_df.groupby("tool")[list(METRICS)].mean().sort_values("mae")
        lines.extend([
            "| Tool | MAE | RMSE | Pearson | Spearman |",
            "| --- | ---: | ---: | ---: | ---: |",
        ])
        for tool, row in overall_g.iterrows():
            if tool in gene_tools:
                lines.append(
                    f"| {tool} | {_fmt(row['mae'])} | {_fmt(row['rmse'])} | "
                    f"{_fmt(row['pearson'], 4)} | {_fmt(row['spearman'], 4)} |"
                )
    lines.append("")

    # Transcript-level MAE by gDNA level
    lines.extend([
        "## Transcript-Level MAE by gDNA Level",
        "",
    ])
    if not tx_df.empty:
        lines.extend(_agg_table(tx_df, ["gdna_label"], tx_tools, "mae"))
    lines.append("")

    # Transcript-level MAE by strand specificity
    lines.extend([
        "## Transcript-Level MAE by Strand Specificity",
        "",
    ])
    if not tx_df.empty:
        lines.extend(_agg_table(tx_df, ["strand_specificity"], tx_tools, "mae"))
    lines.append("")

    # Transcript-level MAE by condition (gDNA × strand)
    lines.extend([
        "## Transcript-Level MAE by Condition (gDNA × Strand Specificity)",
        "",
    ])
    if not tx_df.empty:
        lines.extend(
            _agg_table(tx_df, ["gdna_label", "strand_specificity"], tx_tools, "mae")
        )
    lines.append("")

    # Per-region breakdown (mean across seeds)
    lines.extend([
        "## Per-Region Transcript-Level MAE (mean across seeds & conditions)",
        "",
    ])
    if not tx_df.empty:
        lines.extend(_agg_table(tx_df, ["region"], tx_tools, "mae"))
    lines.append("")

    # Correlation by gDNA level
    lines.extend([
        "## Transcript-Level Pearson by gDNA Level",
        "",
    ])
    if not tx_df.empty:
        lines.extend(
            _agg_table(tx_df, ["gdna_label"], tx_tools, "pearson", decimals=4)
        )
    lines.append("")

    # Gene-level MAE by gDNA level (includes htseq if available)
    lines.extend([
        "## Gene-Level MAE by gDNA Level",
        "",
    ])
    if not gene_df.empty:
        lines.extend(_agg_table(gene_df, ["gdna_label"], gene_tools, "mae"))
    lines.append("")

    # ── Per-transcript diagnostics ──────────────────────────────────

    if not per_tx.empty:
        for tool in tx_tools:
            if tool in per_tx.columns:
                per_tx[f"ae_{tool}"] = (per_tx[tool] - per_tx["truth"]).abs()

        # Error by abundance quartile
        if "abundance" in per_tx.columns:
            per_tx["ab_q"] = pd.qcut(
                per_tx["abundance"].rank(method="first"),
                4,
                labels=["Q1_low", "Q2", "Q3", "Q4_high"],
            )
            lines.extend([
                "## Mean Absolute Error by Abundance Quartile",
                "",
            ])
            ae_cols = {tool: f"ae_{tool}" for tool in tx_tools if f"ae_{tool}" in per_tx.columns}
            header = "| Quartile | " + " | ".join(ae_cols.keys()) + " |"
            sep = "| --- | " + " | ".join(["---:"] * len(ae_cols)) + " |"
            lines.extend([header, sep])
            for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
                sub = per_tx[per_tx["ab_q"] == q]
                if len(sub) == 0:
                    continue
                vals = [_fmt(sub[col].mean()) for col in ae_cols.values()]
                lines.append(f"| {q} | " + " | ".join(vals) + " |")
            lines.append("")

        # Dropout rate
        lines.extend([
            "## Dropout Rate (truth > 0 but predicted ≤ 0)",
            "",
            "| Tool | Dropout Rate |",
            "| --- | ---: |",
        ])
        for tool in tx_tools:
            if tool in per_tx.columns:
                rate = float(
                    ((per_tx["truth"] > 0) & (per_tx[tool] <= 0)).mean()
                )
                lines.append(f"| {tool} | {_fmt(rate, 4)} |")
        lines.append("")

        # Worst transcripts (highest mean AE across all conditions/seeds)
        lines.extend([
            "## Top 20 Worst Transcripts (highest mean abs error, hulkrna)",
            "",
        ])
        if "ae_hulkrna" in per_tx.columns:
            worst = (
                per_tx.groupby("transcript_id")["ae_hulkrna"]
                .mean()
                .sort_values(ascending=False)
                .head(20)
            )
            lines.extend([
                "| Transcript | Mean AE (hulkrna) | Mean Truth | Mean Abundance |",
                "| --- | ---: | ---: | ---: |",
            ])
            for tid, ae in worst.items():
                sub = per_tx[per_tx["transcript_id"] == tid]
                mean_truth = sub["truth"].mean()
                mean_ab = sub["abundance"].mean() if "abundance" in sub.columns else 0
                lines.append(
                    f"| {tid} | {_fmt(ae)} | {_fmt(mean_truth)} | {_fmt(mean_ab)} |"
                )
            lines.append("")

    # Write
    with open(output_dir / "aggregate_summary.md", "w") as f:
        f.write("\n".join(lines) + "\n")

    logger.info("Wrote aggregate reports to %s", output_dir)


# ── CLI ──────────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Aggregate multi-seed benchmark results"
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing seed_* subdirectories",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Output directory for aggregate reports",
    )
    args = parser.parse_args()

    seed_results = load_seed_results(args.input_dir)
    if not seed_results:
        logger.error("No seed results found in %s", args.input_dir)
        return 1

    n_seeds = len(seed_results)
    logger.info("Loaded %d seeds", n_seeds)

    # Build tidy DataFrames
    tx_df, gene_df = aggregate(seed_results)
    logger.info(
        "Aggregated: %d transcript-level rows, %d gene-level rows",
        len(tx_df),
        len(gene_df),
    )

    # Detect whether htseq was included
    include_htseq = "htseq" in gene_df["tool"].values if not gene_df.empty else False

    # Load per-transcript data
    per_tx = load_per_tx_results(args.input_dir)
    logger.info("Loaded %d per-transcript rows", len(per_tx))

    write_aggregate_report(
        tx_df,
        gene_df,
        per_tx,
        args.output_dir,
        n_seeds=n_seeds,
        include_htseq=include_htseq,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
