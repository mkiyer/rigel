#!/usr/bin/env python3
"""
Aggregate multi-seed benchmark results into a unified report.

Reads summary.json from each seed_*/  subdirectory, computes per-region
and overall statistics (mean, median, std across seeds), and writes:

- aggregate_summary.md  — human-readable report
- aggregate_summary.csv — machine-readable results
- aggregate_per_tx.csv  — per-transcript detail across seeds

Usage:
    PYTHONPATH=src conda run -n hulkrna python scripts/aggregate_benchmarks.py \
        --input-dir /path/to/bench_output \
        --output-dir /path/to/bench_output
"""

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
            logger.info("Loaded seed %d: %d regions", seed, len(data))
    return results


def load_per_tx_results(input_dir: Path) -> pd.DataFrame:
    """Load per_transcript_counts.csv from each seed/region subdirectory."""
    parts = []
    for seed_dir in sorted(input_dir.glob("seed_*")):
        seed = int(seed_dir.name.split("_")[-1])
        for ptx in seed_dir.glob("*/per_transcript_counts.csv"):
            region = ptx.parent.name
            df = pd.read_csv(ptx)
            df["seed"] = seed
            df["region"] = region
            parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def aggregate(seed_results: list[tuple[int, list[dict]]]) -> dict:
    """Compute aggregate statistics across seeds."""
    # Collect per-region per-tool metrics across seeds
    region_tool_metrics = defaultdict(lambda: defaultdict(list))

    for seed, regions in seed_results:
        for r in regions:
            region = r["region"]
            for tool_name in ("hulkrna", "salmon", "kallisto"):
                tm = r[tool_name]
                region_tool_metrics[region][tool_name].append({
                    "seed": seed,
                    "mae": tm["mean_abs_error"],
                    "rmse": tm["rmse"],
                    "pearson": tm["pearson"],
                    "spearman": tm["spearman"],
                    "total_abs_error": tm["total_abs_error"],
                    "total_truth": tm["total_truth"],
                    "total_observed": tm["total_observed"],
                })

    return dict(region_tool_metrics)


def write_aggregate_report(
    agg: dict,
    per_tx: pd.DataFrame,
    output_dir: Path,
    n_seeds: int,
) -> None:
    """Write the aggregate markdown report."""
    output_dir.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Aggregate Benchmark Report",
        "",
        f"**Seeds**: {n_seeds}",
        f"**Configuration**: zero gDNA, strand-specificity=1.0, abundance-mode=random",
        "",
    ]

    # === Overall summary across all regions ===
    tool_all = defaultdict(list)  # tool -> list of MAEs across all region-seeds
    for region in sorted(agg.keys()):
        for tool_name in ("hulkrna", "salmon", "kallisto"):
            metrics = agg[region][tool_name]
            for m in metrics:
                tool_all[tool_name].append(m)

    lines.extend([
        "## Overall Summary (all regions, all seeds)",
        "",
        "| Tool | Mean MAE | Median MAE | Std MAE | Mean RMSE | Mean Pearson | Mean Spearman |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for tool in ("hulkrna", "salmon", "kallisto"):
        vals = tool_all[tool]
        maes = [v["mae"] for v in vals]
        rmses = [v["rmse"] for v in vals]
        pears = [v["pearson"] for v in vals]
        spears = [v["spearman"] for v in vals]
        lines.append(
            f"| {tool} | {np.mean(maes):.3f} | {np.median(maes):.3f} | "
            f"{np.std(maes):.3f} | {np.mean(rmses):.3f} | "
            f"{np.mean(pears):.5f} | {np.mean(spears):.5f} |"
        )

    # === Ratios ===
    lines.extend(["", "### hulkrna / salmon ratio (MAE)", ""])
    for region in sorted(agg.keys()):
        h_maes = [m["mae"] for m in agg[region]["hulkrna"]]
        s_maes = [m["mae"] for m in agg[region]["salmon"]]
        if s_maes and np.mean(s_maes) > 0:
            ratio = np.mean(h_maes) / np.mean(s_maes)
            lines.append(f"- **{region}**: {ratio:.2f}x")
    h_all = [m["mae"] for m in tool_all["hulkrna"]]
    s_all = [m["mae"] for m in tool_all["salmon"]]
    if np.mean(s_all) > 0:
        lines.append(f"- **Overall**: {np.mean(h_all) / np.mean(s_all):.2f}x")

    # === Per-region breakdown ===
    lines.extend([
        "",
        "## Per-Region Results (mean across seeds)",
        "",
        "| Region | hulkrna MAE | salmon MAE | kallisto MAE | h/s ratio | hulkrna Pearson | salmon Pearson |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for region in sorted(agg.keys()):
        h = agg[region]["hulkrna"]
        s = agg[region]["salmon"]
        k = agg[region]["kallisto"]
        h_mae = np.mean([m["mae"] for m in h])
        s_mae = np.mean([m["mae"] for m in s])
        k_mae = np.mean([m["mae"] for m in k])
        h_pear = np.mean([m["pearson"] for m in h])
        s_pear = np.mean([m["pearson"] for m in s])
        ratio = h_mae / s_mae if s_mae > 0 else float("inf")
        lines.append(
            f"| {region} | {h_mae:.3f} | {s_mae:.3f} | {k_mae:.3f} | "
            f"{ratio:.2f}x | {h_pear:.4f} | {s_pear:.4f} |"
        )

    # === Per-transcript diagnosis ===
    if not per_tx.empty and "truth" in per_tx.columns:
        lines.extend([
            "",
            "## Diagnostic Breakdown",
            "",
        ])

        # Dropout analysis
        for tool in ("hulkrna", "salmon", "kallisto"):
            if tool not in per_tx.columns:
                continue
            expressed = per_tx[per_tx["truth"] > 0]
            if len(expressed) == 0:
                continue
            dropout = float((expressed[tool] <= 0).mean())
            lines.append(f"- **{tool}** dropout rate (truth>0, predicted<=0): {dropout:.4f}")

        # Error by abundance quartile
        expressed = per_tx[per_tx["truth"] > 0].copy()
        if len(expressed) > 4:
            expressed["ab_q"] = pd.qcut(
                expressed["truth"].rank(method="first"),
                4,
                labels=["Q1_low", "Q2", "Q3", "Q4_high"],
            )
            lines.extend([
                "",
                "### Mean absolute error by truth-count quartile",
                "",
                "| Quartile | n | hulkrna | salmon | kallisto |",
                "| --- | ---: | ---: | ---: | ---: |",
            ])
            for tool in ("hulkrna", "salmon", "kallisto"):
                if tool not in expressed.columns:
                    continue
                expressed[f"ae_{tool}"] = (expressed[tool] - expressed["truth"]).abs()
            for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
                sub = expressed[expressed["ab_q"] == q]
                if len(sub) == 0:
                    continue
                row = f"| {q} | {len(sub)} |"
                for tool in ("hulkrna", "salmon", "kallisto"):
                    col = f"ae_{tool}"
                    if col in sub.columns:
                        row += f" {sub[col].mean():.3f} |"
                    else:
                        row += " — |"
                lines.append(row)

        # Top error transcripts
        if "hulkrna" in per_tx.columns:
            per_tx_mean = (
                per_tx.groupby(["region", "transcript_id"])
                .agg(
                    truth=("truth", "mean"),
                    hulkrna=("hulkrna", "mean"),
                    salmon=("salmon", "mean"),
                    kallisto=("kallisto", "mean"),
                )
                .reset_index()
            )
            per_tx_mean["hulk_ae"] = (per_tx_mean["hulkrna"] - per_tx_mean["truth"]).abs()
            per_tx_mean["salm_ae"] = (per_tx_mean["salmon"] - per_tx_mean["truth"]).abs()
            worst = per_tx_mean.nlargest(20, "hulk_ae")

            lines.extend([
                "",
                "### Top 20 worst hulkrna transcripts (by mean absolute error)",
                "",
                "| Region | Transcript | Truth | hulkrna | salmon | kallisto | hulk_AE | salm_AE |",
                "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
            ])
            for _, row in worst.iterrows():
                lines.append(
                    f"| {row['region']} | {row['transcript_id']} | "
                    f"{row['truth']:.0f} | {row['hulkrna']:.1f} | "
                    f"{row['salmon']:.1f} | {row['kallisto']:.1f} | "
                    f"{row['hulk_ae']:.1f} | {row['salm_ae']:.1f} |"
                )

    md_path = output_dir / "aggregate_summary.md"
    with open(md_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    logger.info("Wrote aggregate report to %s", md_path)

    # === CSV output ===
    rows = []
    for region in sorted(agg.keys()):
        for tool_name in ("hulkrna", "salmon", "kallisto"):
            for m in agg[region][tool_name]:
                rows.append({
                    "region": region,
                    "tool": tool_name,
                    "seed": m["seed"],
                    "mae": m["mae"],
                    "rmse": m["rmse"],
                    "pearson": m["pearson"],
                    "spearman": m["spearman"],
                    "total_abs_error": m["total_abs_error"],
                })
    csv_path = output_dir / "aggregate_summary.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    logger.info("Wrote CSV to %s", csv_path)

    # === Per-transcript aggregate ===
    if not per_tx.empty:
        ptx_path = output_dir / "aggregate_per_tx.csv"
        per_tx.to_csv(ptx_path, index=False)
        logger.info("Wrote per-transcript detail to %s", ptx_path)


def main():
    parser = argparse.ArgumentParser(description="Aggregate multi-seed benchmark results")
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    args = parser.parse_args()

    if args.output_dir is None:
        args.output_dir = args.input_dir

    seed_results = load_seed_results(args.input_dir)
    if not seed_results:
        logger.error("No seed results found in %s", args.input_dir)
        return 1

    per_tx = load_per_tx_results(args.input_dir)
    agg = aggregate(seed_results)
    n_seeds = len(seed_results)
    write_aggregate_report(agg, per_tx, args.output_dir, n_seeds)

    # Print quick summary to stdout
    print(f"\n{'='*60}")
    print(f"Aggregated {n_seeds} seeds, {len(agg)} regions")
    print(f"{'='*60}")
    for tool in ("hulkrna", "salmon", "kallisto"):
        all_maes = []
        for region in agg:
            for m in agg[region][tool]:
                all_maes.append(m["mae"])
        if all_maes:
            print(f"  {tool:12s}: mean MAE={np.mean(all_maes):8.3f}  median={np.median(all_maes):8.3f}")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
