#!/usr/bin/env python3
"""Phase C analysis: 10-region validation of MAP pc=0.01 vs old default and VBEM.

Reads per_transcript_counts.csv from every region × condition and produces
tables comparing MAP pc=0.01 (new default), MAP pc=0.50 (old default),
and VBEM pc=0.25 across all 10 regions.
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# ── Configuration ─────────────────────────────────────────────────
DATA_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_phase_c_10region"
REGIONS = [
    "FGFR2", "EGFR", "FHB", "BRCA1", "HBB",
    "HOXA", "GAPDH", "ELANE", "BCR", "HES4",
]

GDNA_LABELS = ["none", "r20"]
NRNA_LABELS = ["none", "r30"]
SS_VALUES = [1.0]

CONFIGS = [
    ("map_pc001",  "MAP 0.01"),
    ("map_pc050",  "MAP 0.50"),
    ("vbem_pc025", "VBEM 0.25"),
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
        "n_transcripts": len(df), "n_expressed": len(expressed),
        "n_zero": len(zero_truth),
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
    print(f"Regions: {len(set(combined['region']))} | Tools: {len(available)}/{len(ALL_TOOLS)}\n")

    # ── SECTION 1: Grand aggregate ───────────────────────────────
    print("=" * 100)
    print("GRAND AGGREGATE (all 10 regions × 4 conditions pooled)")
    print("=" * 100)
    header = (
        f"{'Config':10s}  "
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
        print(
            f"{label:10s}  "
            f"{m['mae']:7.2f} {m['rmse']:8.2f} {m['spearman']:9.4f} {m['pearson']:8.4f} "
            f"{m['fp_count']:4d} {m['fp_mass']:8.1f} {m['dropout']:5d}"
        )
        agg_rows.append({"config": config_key, "label": label, **m})

    # Delta from old default
    agg_df = pd.DataFrame(agg_rows)
    old_ref = tool_col("map_pc050")
    if old_ref in combined.columns:
        ref = compute_metrics(combined, old_ref)
        print(f"\n  Delta from old default (MAP 0.50):")
        print(f"  {'Config':10s}  {'ΔMAE':>8s} {'ΔSpear':>8s} {'ΔFP':>5s} {'ΔFP mass':>10s} {'ΔDrop':>6s}")
        print(f"  {'-'*55}")
        for config_key, label in CONFIGS:
            tool = tool_col(config_key)
            if tool not in combined.columns:
                continue
            m = compute_metrics(combined, tool)
            print(
                f"  {label:10s}  {m['mae']-ref['mae']:+8.2f} "
                f"{m['spearman']-ref['spearman']:+8.4f} "
                f"{m['fp_count']-ref['fp_count']:+5d} "
                f"{m['fp_mass']-ref['fp_mass']:+10.1f} "
                f"{m['dropout']-ref['dropout']:+6d}"
            )

    # ── SECTION 2: Per-condition aggregate (all regions) ─────────
    print("\n" + "=" * 100)
    print("PER-CONDITION AGGREGATE (all 10 regions pooled per condition)")
    print("=" * 100)

    cond_combos = list(itertools.product(GDNA_LABELS, NRNA_LABELS, SS_VALUES))
    for gdna, nrna, ss in cond_combos:
        cond_label = f"gdna={gdna} nrna={nrna}"
        df = combined[
            (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty:
            continue
        print(f"\n  {cond_label} ({len(df)} transcripts across {len(set(df['region']))} regions)")
        print(f"  {'Config':10s}  {'MAE':>7s} {'Spear':>7s} {'FP':>4s} {'FP mass':>8s} {'Drop':>5s}")
        print(f"  {'-'*55}")
        for config_key, label in CONFIGS:
            tool = tool_col(config_key)
            if tool not in df.columns:
                continue
            m = compute_metrics(df, tool)
            print(
                f"  {label:10s}  "
                f"{m['mae']:7.2f} {m['spearman']:7.4f} {m['fp_count']:4d} "
                f"{m['fp_mass']:8.1f} {m['dropout']:5d}"
            )

    # ── SECTION 3: Per-region breakdown ──────────────────────────
    print("\n" + "=" * 100)
    print("PER-REGION RESULTS (all 4 conditions pooled per region)")
    print("=" * 100)
    header_r = (
        f"  {'Region':8s} {'Config':10s}  "
        f"{'MAE':>7s} {'Spear':>7s} {'FP':>4s} {'FP mass':>8s} {'Drop':>5s} "
        f"{'#tx':>4s} {'#expr':>5s} {'#zero':>5s}"
    )
    print(header_r)
    print("  " + "-" * (len(header_r) - 2))

    region_rows = []
    for region in REGIONS:
        df = combined[combined["region"] == region]
        if df.empty:
            continue
        for config_key, label in CONFIGS:
            tool = tool_col(config_key)
            if tool not in df.columns:
                continue
            m = compute_metrics(df, tool)
            print(
                f"  {region:8s} {label:10s}  "
                f"{m['mae']:7.2f} {m['spearman']:7.4f} {m['fp_count']:4d} "
                f"{m['fp_mass']:8.1f} {m['dropout']:5d} "
                f"{m['n_transcripts']:4d} {m['n_expressed']:5d} {m['n_zero']:5d}"
            )
            region_rows.append({"region": region, "config": config_key, "label": label, **m})
        print()

    # ── SECTION 4: Per-region wins/losses ────────────────────────
    print("\n" + "=" * 100)
    print("HEAD-TO-HEAD: MAP 0.01 vs MAP 0.50 per region (all conditions)")
    print("=" * 100)
    print(f"  {'Region':8s}  {'ΔFP':>5s} {'ΔDrop':>6s} {'ΔMAE':>8s} {'ΔSpear':>8s}  Winner")
    print(f"  {'-'*55}")

    new_tool = tool_col("map_pc001")
    old_tool = tool_col("map_pc050")
    wins = {"new": 0, "old": 0, "tie": 0}

    for region in REGIONS:
        df = combined[combined["region"] == region]
        if df.empty or new_tool not in df.columns or old_tool not in df.columns:
            continue
        m_new = compute_metrics(df, new_tool)
        m_old = compute_metrics(df, old_tool)
        d_fp = m_new["fp_count"] - m_old["fp_count"]
        d_drop = m_new["dropout"] - m_old["dropout"]
        d_mae = m_new["mae"] - m_old["mae"]
        d_sp = m_new["spearman"] - m_old["spearman"]

        # Winner by FP first, then MAE
        if d_fp < 0 or (d_fp == 0 and d_mae < -0.01):
            winner = "NEW ✓"
            wins["new"] += 1
        elif d_fp > 0 or (d_fp == 0 and d_mae > 0.01):
            winner = "OLD"
            wins["old"] += 1
        else:
            winner = "TIE"
            wins["tie"] += 1

        print(
            f"  {region:8s}  {d_fp:+5d} {d_drop:+6d} {d_mae:+8.2f} {d_sp:+8.4f}  {winner}"
        )

    print(f"\n  Score: NEW wins {wins['new']}, OLD wins {wins['old']}, TIES {wins['tie']}")

    # ── SECTION 5: Per-region: MAP 0.01 vs VBEM 0.25 ────────────
    print("\n" + "=" * 100)
    print("HEAD-TO-HEAD: MAP 0.01 vs VBEM 0.25 per region (all conditions)")
    print("=" * 100)
    print(f"  {'Region':8s}  {'ΔFP':>5s} {'ΔDrop':>6s} {'ΔMAE':>8s} {'ΔSpear':>8s}  Winner")
    print(f"  {'-'*55}")

    vbem_tool = tool_col("vbem_pc025")
    wins2 = {"map": 0, "vbem": 0, "tie": 0}

    for region in REGIONS:
        df = combined[combined["region"] == region]
        if df.empty or new_tool not in df.columns or vbem_tool not in df.columns:
            continue
        m_map = compute_metrics(df, new_tool)
        m_vbem = compute_metrics(df, vbem_tool)
        d_fp = m_map["fp_count"] - m_vbem["fp_count"]
        d_drop = m_map["dropout"] - m_vbem["dropout"]
        d_mae = m_map["mae"] - m_vbem["mae"]
        d_sp = m_map["spearman"] - m_vbem["spearman"]

        # Winner: fewer FPs first, then lower MAE
        if d_fp < 0 or (d_fp == 0 and d_mae < -0.01):
            winner = "MAP ✓"
            wins2["map"] += 1
        elif d_fp > 0 or (d_fp == 0 and d_mae > 0.01):
            winner = "VBEM"
            wins2["vbem"] += 1
        else:
            winner = "TIE"
            wins2["tie"] += 1

        print(
            f"  {region:8s}  {d_fp:+5d} {d_drop:+6d} {d_mae:+8.2f} {d_sp:+8.4f}  {winner}"
        )

    print(f"\n  Score: MAP wins {wins2['map']}, VBEM wins {wins2['vbem']}, TIES {wins2['tie']}")

    # ── SECTION 6: Summary ───────────────────────────────────────
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)
    if agg_rows:
        agg_df = pd.DataFrame(agg_rows)
        for _, row in agg_df.iterrows():
            print(
                f"  {row['label']:10s}  MAE={row['mae']:.2f}  Spear={row['spearman']:.4f}  "
                f"FP={int(row['fp_count'])}  Drop={int(row['dropout'])}"
            )


if __name__ == "__main__":
    main()
