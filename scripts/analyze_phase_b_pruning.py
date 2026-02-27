#!/usr/bin/env python3
"""Phase B analysis: post-EM pruning sweep for MAP-EM and VBEM.

Reads per_transcript_counts.csv from every condition and produces a
compact table showing how FP, dropout, MAE, and Spearman vary with
pruning evidence-ratio threshold.
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# ── Configuration ─────────────────────────────────────────────────
DATA_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_phase_b_pruning"
REGIONS = ["FGFR2", "BRCA1"]

GDNA_LABELS = ["none", "r20"]
NRNA_LABELS = ["none", "r30"]
SS_VALUES = [1.0]

# Config names match the YAML hulkrna_configs keys
CONFIGS = [
    ("map_noprune",  "MAP  none"),
    ("map_prune010", "MAP  0.10"),
    ("map_prune020", "MAP  0.20"),
    ("map_prune030", "MAP  0.30"),
    ("map_prune050", "MAP  0.50"),
    ("vbem_noprune", "VBEM none"),
    ("vbem_prune030","VBEM 0.30"),
]

FP_THRESHOLD = 1.0


def tool_col(config_key):
    return f"hulkrna_{config_key}_oracle"


ALL_TOOLS = [tool_col(k) for k, _ in CONFIGS]


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

    available = [c for c in ALL_TOOLS if c in combined.columns]
    missing = [c for c in ALL_TOOLS if c not in combined.columns]
    if missing:
        print(f"WARNING: Missing columns: {missing}")
    print(f"Loaded {len(combined)} rows across {len(all_frames)} conditions")
    print(f"Available tools: {len(available)}/{len(ALL_TOOLS)}\n")

    # ── SECTION 1: Aggregate sweep table ─────────────────────────
    print("=" * 110)
    print("AGGREGATE PRUNING SWEEP (all regions × conditions pooled)")
    print("=" * 110)
    header = (
        f"{'Mode':5s} {'Prune':>6s}  "
        f"{'MAE':>7s} {'RMSE':>8s} {'Spearman':>9s} {'Pearson':>8s} "
        f"{'FP':>4s} {'FP mass':>8s} {'Drop':>5s}"
    )
    print(header)
    print("-" * len(header))

    agg_rows = []
    for config_key, label in CONFIGS:
        tool = tool_col(config_key)
        if tool not in combined.columns:
            continue
        m = compute_metrics(combined, tool)
        mode_str = label.split()[0]
        prune_str = label.split()[1]
        print(
            f"{mode_str:5s} {prune_str:>6s}  "
            f"{m['mae']:7.2f} {m['rmse']:8.2f} {m['spearman']:9.4f} {m['pearson']:8.4f} "
            f"{m['fp_count']:4d} {m['fp_mass']:8.1f} {m['dropout']:5d}"
        )
        agg_rows.append({"config": config_key, "label": label, **m})

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
        print(f"  {'Label':12s}  {'MAE':>7s} {'Spear':>7s} {'FP':>4s} {'FP mass':>8s} {'Drop':>5s}")
        print(f"  {'-'*55}")

        for config_key, label in CONFIGS:
            tool = tool_col(config_key)
            if tool not in df.columns:
                continue
            m = compute_metrics(df, tool)
            print(
                f"  {label:12s}  "
                f"{m['mae']:7.2f} {m['spearman']:7.4f} {m['fp_count']:4d} "
                f"{m['fp_mass']:8.1f} {m['dropout']:5d}"
            )
        print()

    # ── SECTION 3: Delta from MAP baseline ───────────────────────
    print("\n" + "=" * 110)
    print("DELTA FROM MAP BASELINE (no-prune)")
    print("=" * 110)
    ref_tool = tool_col("map_noprune")
    if ref_tool in combined.columns:
        ref = compute_metrics(combined, ref_tool)
        print(f"  {'Label':12s}  {'ΔMAE':>8s} {'ΔSpear':>8s} {'ΔFP':>5s} {'ΔFP mass':>10s} {'ΔDrop':>6s}")
        print(f"  {'-'*55}")
        for config_key, label in CONFIGS:
            tool = tool_col(config_key)
            if tool not in combined.columns:
                continue
            m = compute_metrics(combined, tool)
            marker = " ◄" if config_key == "map_noprune" else ""
            print(
                f"  {label:12s}  {m['mae']-ref['mae']:+8.2f} "
                f"{m['spearman']-ref['spearman']:+8.4f} "
                f"{m['fp_count']-ref['fp_count']:+5d} "
                f"{m['fp_mass']-ref['fp_mass']:+10.1f} "
                f"{m['dropout']-ref['dropout']:+6d}{marker}"
            )

    # ── SECTION 4: Pareto frontier ───────────────────────────────
    print("\n" + "=" * 110)
    print("PARETO FRONTIER: FP count vs Dropout")
    print("=" * 110)
    agg_df = pd.DataFrame(agg_rows)
    print(f"  {'Label':12s}  {'FP':>4s} {'Drop':>5s} {'MAE':>7s} {'Spear':>7s}  Pareto?")
    print(f"  {'-'*55}")

    for _, row in agg_df.sort_values(["fp_count", "dropout"]).iterrows():
        dominated = agg_df[
            (agg_df["fp_count"] <= row["fp_count"])
            & (agg_df["dropout"] <= row["dropout"])
            & ((agg_df["fp_count"] < row["fp_count"]) | (agg_df["dropout"] < row["dropout"]))
        ]
        is_pareto = "✓" if len(dominated) == 0 else ""
        print(
            f"  {row['label']:12s}  "
            f"{int(row['fp_count']):4d} {int(row['dropout']):5d} "
            f"{row['mae']:7.2f} {row['spearman']:7.4f}  {is_pareto}"
        )

    # ── SECTION 5: MAP+prune vs VBEM head-to-head ────────────────
    print("\n" + "=" * 110)
    print("HEAD-TO-HEAD: Best MAP+prune vs VBEM (per condition)")
    print("=" * 110)

    map_prune_configs = [k for k, _ in CONFIGS if k.startswith("map_prune")]
    vbem_ref = tool_col("vbem_noprune")

    if vbem_ref in combined.columns and map_prune_configs:
        print(f"\n  {'Condition':30s}  {'Best MAP+prune':14s} {'MAP FP':>7s} {'VBEM FP':>8s} "
              f"{'MAP MAE':>8s} {'VBEM MAE':>9s} {'MAP Drop':>9s} {'VBEM Drop':>10s}")
        print(f"  {'-'*110}")

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

            vm = compute_metrics(df, vbem_ref)

            # Find best MAP+prune by FP then MAE
            best_config = None
            best_m = None
            for ck in map_prune_configs:
                t = tool_col(ck)
                if t not in df.columns:
                    continue
                mm = compute_metrics(df, t)
                if best_m is None or mm["fp_count"] < best_m["fp_count"] or (
                    mm["fp_count"] == best_m["fp_count"] and mm["mae"] < best_m["mae"]
                ):
                    best_config = ck
                    best_m = mm

            if best_m is not None:
                print(
                    f"  {cond_label:30s}  {best_config:14s} "
                    f"{best_m['fp_count']:7d} {vm['fp_count']:8d} "
                    f"{best_m['mae']:8.2f} {vm['mae']:9.2f} "
                    f"{best_m['dropout']:9d} {vm['dropout']:10d}"
                )


if __name__ == "__main__":
    main()
