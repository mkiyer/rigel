#!/usr/bin/env python3
"""Phase A analysis: pseudocount sweep for MAP-EM vs VBEM.

Reads per_transcript_counts.csv from every condition and produces a
compact table showing how each metric varies with em_pseudocount for
both MAP and VBEM.
"""
from __future__ import annotations

import itertools
import os
import sys
import re

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# ── Configuration ─────────────────────────────────────────────────
DATA_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_phase_a_pseudocount"
REGIONS = ["FGFR2", "BRCA1"]

GDNA_LABELS = ["none", "r20"]
NRNA_LABELS = ["none", "r30"]
SS_VALUES = [1.0]

EM_MODES = ["map", "vbem"]
PSEUDOCOUNTS = [0.01, 0.05, 0.1, 0.25, 0.5, 1.0]
PC_LABELS = ["001", "005", "010", "025", "050", "100"]

FP_THRESHOLD = 1.0


def tool_col(mode, pc_label):
    return f"hulkrna_{mode}_pc{pc_label}_oracle"


ALL_TOOLS = [tool_col(m, p) for m in EM_MODES for p in PC_LABELS]


def load_condition(region, gdna_label, nrna_label, ss):
    cond = f"gdna_{gdna_label}_nrna_{nrna_label}_ss_{ss:.2f}"
    csv_path = os.path.join(DATA_DIR, region, cond, "per_transcript_counts.csv")
    if not os.path.exists(csv_path):
        return None
    df = pd.read_csv(csv_path)
    df["region"] = region
    df["gdna_label"] = gdna_label
    df["nrna_label"] = nrna_label
    df["ss"] = ss
    df["condition"] = cond
    return df


def compute_metrics(df, tool):
    truth = df["truth"]
    pred = df[tool]
    err = pred - truth
    zero_truth = df[truth == 0]
    expressed = df[truth > 0]

    mae = np.abs(err).mean()
    rmse = np.sqrt((err ** 2).mean())
    r_sp = spearmanr(truth, pred)[0] if len(df) > 2 else np.nan
    r_pe = pearsonr(truth, pred)[0] if len(df) > 2 else np.nan

    fp_count = (zero_truth[tool] > FP_THRESHOLD).sum()
    fp_mass = zero_truth.loc[zero_truth[tool] > FP_THRESHOLD, tool].sum()
    dropout = (expressed[tool] < 0.5).sum() if len(expressed) > 0 else 0

    return {
        "mae": mae, "rmse": rmse, "spearman": r_sp, "pearson": r_pe,
        "fp_count": fp_count, "fp_mass": fp_mass, "dropout": dropout,
    }


def main():
    # ── Load all conditions ──────────────────────────────────────
    all_frames = []
    conditions = list(itertools.product(REGIONS, GDNA_LABELS, NRNA_LABELS, SS_VALUES))

    for region, gdna, nrna, ss in conditions:
        df = load_condition(region, gdna, nrna, ss)
        if df is not None:
            all_frames.append(df)

    if not all_frames:
        print("ERROR: No data found. Did the benchmark complete?")
        sys.exit(1)

    combined = pd.concat(all_frames, ignore_index=True)

    # Verify all tool columns exist
    available = [c for c in ALL_TOOLS if c in combined.columns]
    missing = [c for c in ALL_TOOLS if c not in combined.columns]
    if missing:
        print(f"WARNING: Missing columns: {missing}")
    print(f"Loaded {len(combined)} rows across {len(all_frames)} conditions")
    print(f"Available tools: {len(available)}/{len(ALL_TOOLS)}\n")

    # ── SECTION 1: Aggregate sweep table ─────────────────────────
    print("=" * 110)
    print("AGGREGATE PSEUDOCOUNT SWEEP (all regions × conditions pooled)")
    print("=" * 110)
    header = (
        f"{'Mode':5s} {'PC':>5s}  "
        f"{'MAE':>7s} {'RMSE':>8s} {'Spearman':>9s} {'Pearson':>8s} "
        f"{'FP':>4s} {'FP mass':>8s} {'Drop':>5s}"
    )
    print(header)
    print("-" * len(header))

    agg_rows = []
    for mode in EM_MODES:
        for pc, pcl in zip(PSEUDOCOUNTS, PC_LABELS):
            tool = tool_col(mode, pcl)
            if tool not in combined.columns:
                continue
            m = compute_metrics(combined, tool)
            print(
                f"{mode:5s} {pc:5.2f}  "
                f"{m['mae']:7.2f} {m['rmse']:8.2f} {m['spearman']:9.4f} {m['pearson']:8.4f} "
                f"{m['fp_count']:4d} {m['fp_mass']:8.1f} {m['dropout']:5d}"
            )
            agg_rows.append({"mode": mode, "pc": pc, **m})
        print()

    # ── SECTION 2: Per-condition breakdown ───────────────────────
    print("\n" + "=" * 110)
    print("PER-CONDITION BREAKDOWN")
    print("=" * 110)

    for region, gdna, nrna, ss in conditions:
        cond_label = f"{region} gdna={gdna} nrna={nrna}"
        df = combined[
            (combined["region"] == region)
            & (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty:
            continue

        print(f"\n  {cond_label}")
        print(f"  {'Mode':5s} {'PC':>5s}  {'MAE':>7s} {'Spear':>7s} {'FP':>4s} {'FP mass':>8s} {'Drop':>5s}")
        print(f"  {'-'*55}")

        for mode in EM_MODES:
            for pc, pcl in zip(PSEUDOCOUNTS, PC_LABELS):
                tool = tool_col(mode, pcl)
                if tool not in df.columns:
                    continue
                m = compute_metrics(df, tool)
                print(
                    f"  {mode:5s} {pc:5.2f}  "
                    f"{m['mae']:7.2f} {m['spearman']:7.4f} {m['fp_count']:4d} "
                    f"{m['fp_mass']:8.1f} {m['dropout']:5d}"
                )
            print()

    # ── SECTION 3: Sensitivity curves ────────────────────────────
    print("\n" + "=" * 110)
    print("SENSITIVITY TO PSEUDOCOUNT (aggregate)")
    print("=" * 110)

    for mode in EM_MODES:
        print(f"\n  {mode.upper()} mode:")
        print(f"  {'PC':>5s}  {'ΔMAE vs 0.5':>12s} {'ΔSpear':>8s} {'ΔFP':>5s} {'ΔFP mass':>10s} {'ΔDrop':>6s}")
        print(f"  {'-'*55}")

        # Reference: pc=0.5 for this mode
        ref_tool = tool_col(mode, "050")
        if ref_tool not in combined.columns:
            continue
        ref = compute_metrics(combined, ref_tool)

        for pc, pcl in zip(PSEUDOCOUNTS, PC_LABELS):
            tool = tool_col(mode, pcl)
            if tool not in combined.columns:
                continue
            m = compute_metrics(combined, tool)
            d_mae = m["mae"] - ref["mae"]
            d_sp = m["spearman"] - ref["spearman"]
            d_fp = m["fp_count"] - ref["fp_count"]
            d_fpm = m["fp_mass"] - ref["fp_mass"]
            d_drop = m["dropout"] - ref["dropout"]
            marker = " ◄" if pc == 0.5 else ""
            print(
                f"  {pc:5.2f}  {d_mae:+12.2f} {d_sp:+8.4f} {d_fp:+5d} "
                f"{d_fpm:+10.1f} {d_drop:+6d}{marker}"
            )

    # ── SECTION 4: Best config per metric ────────────────────────
    print("\n" + "=" * 110)
    print("BEST CONFIGURATION PER METRIC")
    print("=" * 110)

    agg_df = pd.DataFrame(agg_rows)
    for metric, ascending in [
        ("mae", True), ("spearman", False),
        ("fp_count", True), ("fp_mass", True), ("dropout", True),
    ]:
        best = agg_df.sort_values(metric, ascending=ascending).iloc[0]
        print(
            f"  Best {metric:10s}: {best['mode']:5s} pc={best['pc']:.2f}  "
            f"({metric}={best[metric]:.4f})"
        )

    # ── SECTION 5: Pareto frontier ───────────────────────────────
    print("\n" + "=" * 110)
    print("PARETO FRONTIER: FP count vs Dropout")
    print("=" * 110)
    print(f"  {'Mode':5s} {'PC':>5s}  {'FP':>4s} {'Drop':>5s} {'MAE':>7s} {'Spear':>7s}  Pareto?")
    print(f"  {'-'*60}")

    for _, row in agg_df.sort_values(["fp_count", "dropout"]).iterrows():
        # Simple Pareto check: not dominated on both FP and dropout
        dominated = agg_df[
            (agg_df["fp_count"] <= row["fp_count"])
            & (agg_df["dropout"] <= row["dropout"])
            & ((agg_df["fp_count"] < row["fp_count"]) | (agg_df["dropout"] < row["dropout"]))
        ]
        is_pareto = "✓" if len(dominated) == 0 else ""
        print(
            f"  {row['mode']:5s} {row['pc']:5.2f}  "
            f"{int(row['fp_count']):4d} {int(row['dropout']):5d} "
            f"{row['mae']:7.2f} {row['spearman']:7.4f}  {is_pareto}"
        )

    # ── SECTION 6: MAP vs VBEM matched-pseudocount comparison ────
    print("\n" + "=" * 110)
    print("MATCHED PSEUDOCOUNT: MAP vs VBEM at each level")
    print("=" * 110)
    print(f"  {'PC':>5s}  {'MAP MAE':>8s} {'VBEM MAE':>9s} {'MAP FP':>7s} {'VBEM FP':>8s} "
          f"{'MAP Drop':>9s} {'VBEM Drop':>10s} {'MAP Sp':>7s} {'VBEM Sp':>8s}")
    print(f"  {'-'*90}")

    for pc, pcl in zip(PSEUDOCOUNTS, PC_LABELS):
        map_tool = tool_col("map", pcl)
        vbem_tool = tool_col("vbem", pcl)
        if map_tool not in combined.columns or vbem_tool not in combined.columns:
            continue
        mm = compute_metrics(combined, map_tool)
        vm = compute_metrics(combined, vbem_tool)
        print(
            f"  {pc:5.2f}  {mm['mae']:8.2f} {vm['mae']:9.2f} "
            f"{mm['fp_count']:7d} {vm['fp_count']:8d} "
            f"{mm['dropout']:9d} {vm['dropout']:10d} "
            f"{mm['spearman']:7.4f} {vm['spearman']:8.4f}"
        )


if __name__ == "__main__":
    main()
