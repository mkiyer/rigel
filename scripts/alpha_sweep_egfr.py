#!/usr/bin/env python3
"""Sweep overhang alpha values on EGFR region and aggregate results.

Runs benchmark_region_competition.py for the EGFR region at each alpha
value, then aggregates transcript-level metrics across all alpha values.

Usage:
    PYTHONPATH=src conda run -n hulkrna python scripts/alpha_sweep_egfr.py
"""
from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Configuration ────────────────────────────────────────────────────

GENOME = Path("/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genome_controls.fasta.bgz")
GTF = Path("/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz")
OUTDIR_ROOT = Path("/Users/mkiyer/Downloads/hulkrna_runs/alpha_sweep_egfr")

ALPHA_VALUES = [0.0, 0.001, 0.01, 0.1, 0.5, 1.0]
ALPHA_LABELS = {
    0.0: "0.0 (binary)",
    0.001: "0.001",
    0.01: "0.01 (proposed)",
    0.1: "0.1",
    0.5: "0.5",
    1.0: "1.0 (off)",
}

EGFR_REGION = "chr7:55019017-55211628"

N_FRAGMENTS = 50000
SIM_SEED = 101
PIPELINE_SEED = 101
ABUNDANCE_SEED = 101
GDNA_LEVELS = "none,low,moderate,high"
STRAND_SPECIFICITIES = "0.95,0.99,1.0"
THREADS = 4


def run_benchmark(alpha: float, outdir: Path) -> bool:
    """Run benchmark_region_competition.py for one alpha value."""
    cmd = [
        sys.executable,
        "scripts/benchmark_region_competition.py",
        "--genome", str(GENOME),
        "--gtf", str(GTF),
        "--region", EGFR_REGION,
        "--outdir", str(outdir),
        "--n-fragments", str(N_FRAGMENTS),
        "--sim-seed", str(SIM_SEED),
        "--pipeline-seed", str(PIPELINE_SEED),
        "--abundance-seed", str(ABUNDANCE_SEED),
        "--gdna-levels", GDNA_LEVELS,
        "--strand-specificities", STRAND_SPECIFICITIES,
        "--abundance-mode", "random",
        "--abundance-min", "0.01",
        "--abundance-max", "100000",
        "--threads", str(THREADS),
        "--include-htseq", "--htseq-conda-env", "htseq",
        "--keep-going",
        "--verbose",
    ]
    if alpha is not None:
        cmd.extend(["--overhang-alpha", str(alpha)])

    print(f"\n{'='*72}")
    print(f"Running alpha={alpha} -> {outdir}")
    print(f"{'='*72}\n")
    env = os.environ.copy()
    env["PYTHONPATH"] = "src"
    result = subprocess.run(cmd, env=env)
    return result.returncode == 0


def load_summary(outdir: Path) -> list[dict] | None:
    """Load summary.json from a benchmark run."""
    sj = outdir / "summary.json"
    if not sj.exists():
        return None
    return json.loads(sj.read_text())


def aggregate_results(results: dict[float, list[dict]]) -> pd.DataFrame:
    """Build tidy DataFrame from all alpha results."""
    rows = []
    for alpha, conditions in results.items():
        for r in conditions:
            base = {
                "alpha": alpha,
                "alpha_label": ALPHA_LABELS.get(alpha, str(alpha)),
                "region": r["region"],
                "gdna_label": r["gdna_label"],
                "strand_specificity": r["strand_specificity"],
                "n_fragments": r["n_fragments"],
                "n_gdna_actual": r["n_gdna_actual"],
            }
            for tool, tm in r.get("transcript_metrics", {}).items():
                rows.append({
                    **base,
                    "tool": tool,
                    "mae": tm["mean_abs_error"],
                    "rmse": tm["rmse"],
                    "pearson": tm["pearson"],
                    "spearman": tm["spearman"],
                    "elapsed_sec": tm["elapsed_sec"],
                    "total_abs_error": tm["total_abs_error"],
                })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def load_per_tx(outdir: Path, alpha: float) -> pd.DataFrame:
    """Load per-transcript counts from a benchmark run.

    Directory layout: outdir/<region>/<gdna_label>_ss_<ss>/per_transcript_counts.csv
    """
    parts = []
    for ptx in outdir.rglob("per_transcript_counts.csv"):
        condition = ptx.parent.name          # e.g. "none_ss_0.95"
        region = ptx.parent.parent.name      # e.g. "EGFR"
        df = pd.read_csv(ptx)
        df["alpha"] = alpha
        df["condition"] = condition
        df["region"] = region
        cond_parts = condition.split("_ss_")
        if len(cond_parts) == 2:
            df["gdna_label"] = cond_parts[0]
            df["strand_specificity"] = float(cond_parts[1])
        parts.append(df)
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def fmt(val, decimals=3):
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return "—"
    return f"{val:.{decimals}f}"


def write_report(df: pd.DataFrame, per_tx: pd.DataFrame, out_path: Path):
    """Write the alpha sweep analysis report."""
    lines = [
        "# Exponential Overhang Penalty: Alpha Sweep on EGFR Region",
        "",
        "## Background",
        "",
        "The exponential overhang penalty replaces the old `overlap_exponent × log(overlap_frac)` model.",
        "For each base of a fragment outside the transcript's exon boundary, probability is multiplied",
        "by `alpha`. In log-space: `log_penalty = b_out × log(alpha)`.",
        "",
        "| Alpha | log(alpha)/base | 1-base penalty | 5-base penalty |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for a in ALPHA_VALUES:
        if a <= 0:
            lp = "−∞"
            p1 = "0"
            p5 = "0"
        elif a >= 1.0:
            lp = "0"
            p1 = "1.0"
            p5 = "1.0"
        else:
            lp = fmt(np.log(a), 3)
            p1 = fmt(a, 6)
            p5 = fmt(a**5, 6)
        lines.append(f"| {a} | {lp} | {p1} | {p5} |")

    lines.extend(["", "## Configuration", ""])
    lines.append(f"- Region: EGFR ({EGFR_REGION})")
    lines.append(f"- Fragments: {N_FRAGMENTS:,}")
    lines.append(f"- Seeds: sim={SIM_SEED}, pipeline={PIPELINE_SEED}, abundance={ABUNDANCE_SEED}")
    lines.append(f"- gDNA levels: {GDNA_LEVELS}")
    lines.append(f"- Strand specificities: {STRAND_SPECIFICITIES}")
    lines.append("")

    if df.empty:
        lines.append("*No results available.*")
        out_path.write_text("\n".join(lines) + "\n")
        return

    # ── Overall MAE by alpha ────────────────────────────────────────

    lines.extend([
        "## Overall Transcript-Level MAE by Alpha",
        "",
        "Mean across all conditions (gDNA levels × strand specificities).",
        "",
    ])

    # Pivot: tool as columns, alpha as rows
    tools = ["hulkrna", "hulkrna_mm", "salmon", "kallisto"]
    lines.append("| Alpha | " + " | ".join(tools) + " |")
    lines.append("| ---: | " + " | ".join(["---:"] * len(tools)) + " |")
    for alpha in ALPHA_VALUES:
        sub = df[df["alpha"] == alpha]
        if sub.empty:
            continue
        vals = []
        for tool in tools:
            tsub = sub[sub["tool"] == tool]
            vals.append(fmt(tsub["mae"].mean()) if len(tsub) else "—")
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")

    # ── RMSE by alpha ───────────────────────────────────────────────

    lines.extend([
        "## Overall Transcript-Level RMSE by Alpha",
        "",
    ])
    lines.append("| Alpha | " + " | ".join(tools) + " |")
    lines.append("| ---: | " + " | ".join(["---:"] * len(tools)) + " |")
    for alpha in ALPHA_VALUES:
        sub = df[df["alpha"] == alpha]
        if sub.empty:
            continue
        vals = []
        for tool in tools:
            tsub = sub[sub["tool"] == tool]
            vals.append(fmt(tsub["rmse"].mean()) if len(tsub) else "—")
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")

    # ── Pearson by alpha ────────────────────────────────────────────

    lines.extend([
        "## Overall Transcript-Level Pearson by Alpha",
        "",
    ])
    lines.append("| Alpha | " + " | ".join(tools) + " |")
    lines.append("| ---: | " + " | ".join(["---:"] * len(tools)) + " |")
    for alpha in ALPHA_VALUES:
        sub = df[df["alpha"] == alpha]
        if sub.empty:
            continue
        vals = []
        for tool in tools:
            tsub = sub[sub["tool"] == tool]
            vals.append(fmt(tsub["pearson"].mean(), 5) if len(tsub) else "—")
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")

    # ── MAE by gDNA level ───────────────────────────────────────────

    lines.extend([
        "## hulkrna MAE by Alpha × gDNA Level",
        "",
    ])
    gdna_levels = sorted(df["gdna_label"].unique(), key=lambda x: ["none", "low", "moderate", "high"].index(x) if x in ["none", "low", "moderate", "high"] else 99)
    lines.append("| Alpha | " + " | ".join(gdna_levels) + " |")
    lines.append("| ---: | " + " | ".join(["---:"] * len(gdna_levels)) + " |")
    for alpha in ALPHA_VALUES:
        sub = df[(df["alpha"] == alpha) & (df["tool"] == "hulkrna")]
        if sub.empty:
            continue
        vals = []
        for gl in gdna_levels:
            gsub = sub[sub["gdna_label"] == gl]
            vals.append(fmt(gsub["mae"].mean()) if len(gsub) else "—")
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")

    # ── MAE by strand specificity ───────────────────────────────────

    lines.extend([
        "## hulkrna MAE by Alpha × Strand Specificity",
        "",
    ])
    ss_values = sorted(df["strand_specificity"].unique())
    ss_labels = [f"SS={s:.2f}" for s in ss_values]
    lines.append("| Alpha | " + " | ".join(ss_labels) + " |")
    lines.append("| ---: | " + " | ".join(["---:"] * len(ss_values)) + " |")
    for alpha in ALPHA_VALUES:
        sub = df[(df["alpha"] == alpha) & (df["tool"] == "hulkrna")]
        if sub.empty:
            continue
        vals = []
        for ss in ss_values:
            ssub = sub[sub["strand_specificity"] == ss]
            vals.append(fmt(ssub["mae"].mean()) if len(ssub) else "—")
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | " + " | ".join(vals) + " |")
    lines.append("")

    # ── Per-isoform analysis ────────────────────────────────────────

    if not per_tx.empty and "hulkrna" in per_tx.columns:
        lines.extend([
            "## Per-Isoform Error (hulkrna): EGFR T1 vs T2 by Alpha",
            "",
            "T1 = ENST00000275493 (main transcript, longer)",
            "T2 = ENST00000450046 (shorter isoform)",
            "",
            "| Alpha | Condition | T1 Truth | T1 Est | T1 Bias | T2 Truth | T2 Est | T2 Bias |",
            "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
        ])
        t1_ids = per_tx["transcript_id"].str.startswith("ENST00000275493")
        t2_ids = per_tx["transcript_id"].str.startswith("ENST00000450046")

        for alpha in ALPHA_VALUES:
            asub = per_tx[per_tx["alpha"] == alpha]
            if asub.empty:
                continue
            # Aggregate across conditions for summary
            t1 = asub[t1_ids & (asub["alpha"] == alpha)]
            t2 = asub[t2_ids & (asub["alpha"] == alpha)]
            if not t1.empty and not t2.empty:
                label = ALPHA_LABELS.get(alpha, str(alpha))
                lines.append(
                    f"| {label} | (mean) | "
                    f"{fmt(t1['truth'].mean(), 0)} | "
                    f"{fmt(t1['hulkrna'].mean(), 1)} | "
                    f"{fmt(t1['hulkrna'].mean() - t1['truth'].mean(), 1)} | "
                    f"{fmt(t2['truth'].mean(), 0)} | "
                    f"{fmt(t2['hulkrna'].mean(), 1)} | "
                    f"{fmt(t2['hulkrna'].mean() - t2['truth'].mean(), 1)} |"
                )
        lines.append("")

        # Summary table: mean absolute T1/T2 bias by alpha
        lines.extend([
            "## Summary: Mean Absolute T1/T2 Bias by Alpha",
            "",
            "| Alpha | |T1 Bias| | |T2 Bias| | Combined |",
            "| ---: | ---: | ---: | ---: |",
        ])
        for alpha in ALPHA_VALUES:
            asub = per_tx[per_tx["alpha"] == alpha]
            if asub.empty:
                continue
            t1 = asub[t1_ids & (asub["alpha"] == alpha)]
            t2 = asub[t2_ids & (asub["alpha"] == alpha)]
            if not t1.empty and not t2.empty:
                b1 = abs(t1["hulkrna"].mean() - t1["truth"].mean())
                b2 = abs(t2["hulkrna"].mean() - t2["truth"].mean())
                label = ALPHA_LABELS.get(alpha, str(alpha))
                lines.append(f"| {label} | {fmt(b1, 1)} | {fmt(b2, 1)} | {fmt(b1 + b2, 1)} |")
        lines.append("")

    # ── Comparison with salmon ──────────────────────────────────────

    lines.extend([
        "## hulkrna vs Salmon: Overall MAE by Alpha",
        "",
        "| Alpha | hulkrna MAE | salmon MAE | Ratio | Better? |",
        "| ---: | ---: | ---: | ---: | --- |",
    ])
    for alpha in ALPHA_VALUES:
        h = df[(df["alpha"] == alpha) & (df["tool"] == "hulkrna")]
        s = df[(df["alpha"] == alpha) & (df["tool"] == "salmon")]
        if h.empty or s.empty:
            continue
        h_mae = h["mae"].mean()
        s_mae = s["mae"].mean()
        ratio = h_mae / s_mae if s_mae > 0 else float("inf")
        better = "hulkrna" if h_mae < s_mae else "salmon"
        label = ALPHA_LABELS.get(alpha, str(alpha))
        lines.append(f"| {label} | {fmt(h_mae)} | {fmt(s_mae)} | {fmt(ratio, 2)} | {better} |")
    lines.append("")

    out_path.write_text("\n".join(lines) + "\n")
    print(f"\nReport written to {out_path}")


def main():
    OUTDIR_ROOT.mkdir(parents=True, exist_ok=True)

    # Run benchmarks for each alpha
    all_results: dict[float, list[dict]] = {}
    all_per_tx: list[pd.DataFrame] = []

    for alpha in ALPHA_VALUES:
        alpha_dir = OUTDIR_ROOT / f"alpha_{alpha}"
        alpha_dir.mkdir(parents=True, exist_ok=True)

        ok = run_benchmark(alpha, alpha_dir)
        if not ok:
            print(f"WARNING: benchmark failed for alpha={alpha}")
            continue

        summary = load_summary(alpha_dir)
        if summary:
            all_results[alpha] = summary

        ptx = load_per_tx(alpha_dir, alpha)
        if not ptx.empty:
            all_per_tx.append(ptx)

    # Aggregate
    df = aggregate_results(all_results)
    per_tx = pd.concat(all_per_tx, ignore_index=True) if all_per_tx else pd.DataFrame()

    # Save raw data
    if not df.empty:
        df.to_csv(OUTDIR_ROOT / "alpha_sweep_metrics.csv", index=False)
    if not per_tx.empty:
        per_tx.to_csv(OUTDIR_ROOT / "alpha_sweep_per_tx.csv", index=False)

    # Write report
    report_path = Path("docs/alpha_sweep_egfr_report.md")
    write_report(df, per_tx, report_path)

    print(f"\nDone! Report at {report_path}")
    print(f"Raw data at {OUTDIR_ROOT}")


if __name__ == "__main__":
    main()
