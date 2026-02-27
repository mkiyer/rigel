#!/usr/bin/env python3
"""Comprehensive VBEM diagnostic: MAP vs VBEM across all conditions.

Reads per_transcript_counts.csv from every (region × gdna × nrna × ss)
condition, compares MAP vs VBEM on:
  - False positives (truth=0, predicted > threshold)
  - Negative control leakage (specific t_ctrl transcripts)
  - Dropout (truth>0, predicted<0.5)
  - MAE, RMSE, Pearson, Spearman
  - gDNA/nRNA separation accuracy
"""
from __future__ import annotations

import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# ── Configuration ─────────────────────────────────────────────────
DATA_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_vbem_comprehensive"
REGIONS = ["FGFR2", "BRCA1"]

GDNA_LABELS = ["none", "r20", "r100"]
NRNA_LABELS = ["none", "r30"]
SS_VALUES = [1.0, 0.9]

MAP_COL = "hulkrna_map_oracle"
VBEM_COL = "hulkrna_vbem_oracle"
TOOLS = [MAP_COL, VBEM_COL, "salmon", "kallisto"]

FP_THRESHOLD = 1.0  # count threshold for false positives


def load_condition(region, gdna_label, nrna_label, ss):
    """Load per_transcript_counts.csv for one condition."""
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
    """Compute standard metrics for a tool column."""
    truth = df["truth"]
    pred = df[tool]
    err = pred - truth
    zero_truth = df[truth == 0]
    expressed = df[truth > 0]

    n_total = len(df)
    mae = np.abs(err).mean()
    rmse = np.sqrt((err ** 2).mean())
    r_pe = pearsonr(truth, pred)[0] if n_total > 2 else np.nan
    r_sp = spearmanr(truth, pred)[0] if n_total > 2 else np.nan

    fp_count = (zero_truth[tool] > FP_THRESHOLD).sum()
    fp_mass = zero_truth.loc[zero_truth[tool] > FP_THRESHOLD, tool].sum()
    dropout = (expressed[tool] < 0.5).sum() if len(expressed) > 0 else 0
    dropout_mass = expressed.loc[expressed[tool] < 0.5, "truth"].sum() if len(expressed) > 0 else 0.0

    return {
        "n_total": n_total,
        "n_expressed": len(expressed),
        "n_zero_truth": len(zero_truth),
        "mae": mae,
        "rmse": rmse,
        "pearson": r_pe,
        "spearman": r_sp,
        "fp_count": fp_count,
        "fp_mass": fp_mass,
        "dropout": dropout,
        "dropout_mass": dropout_mass,
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
    print(f"Loaded {len(combined)} transcript rows across {len(all_frames)} conditions\n")

    # ── SECTION 1: Per-condition summary ─────────────────────────
    print("=" * 100)
    print("PER-CONDITION SUMMARY: MAP vs VBEM")
    print("=" * 100)
    header = (
        f"{'Condition':48s} {'Tool':30s} "
        f"{'MAE':>7s} {'Spearman':>9s} {'FP':>4s} {'FP mass':>8s} "
        f"{'Drop':>5s} {'Pearson':>8s}"
    )
    print(header)
    print("-" * len(header))

    summary_rows = []
    for region, gdna, nrna, ss in conditions:
        cond_label = f"{region} gdna={gdna} nrna={nrna} ss={ss:.2f}"
        cond_name = f"gdna_{gdna}_nrna_{nrna}_ss_{ss:.2f}"
        df = combined[
            (combined["region"] == region)
            & (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty:
            continue

        for tool in [MAP_COL, VBEM_COL]:
            if tool not in df.columns:
                continue
            m = compute_metrics(df, tool)
            short_tool = tool.replace("hulkrna_", "").replace("_oracle", "")
            print(
                f"  {cond_label:46s} {short_tool:30s} "
                f"{m['mae']:7.2f} {m['spearman']:9.4f} {m['fp_count']:4d} "
                f"{m['fp_mass']:8.1f} {m['dropout']:5d} {m['pearson']:8.4f}"
            )
            summary_rows.append({
                "region": region, "gdna": gdna, "nrna": nrna, "ss": ss,
                "tool": short_tool, **m,
            })
        print()

    # ── SECTION 2: VBEM improvement summary ──────────────────────
    print("\n" + "=" * 100)
    print("VBEM vs MAP: IMPROVEMENT SUMMARY")
    print("=" * 100)
    print(
        f"{'Condition':48s} "
        f"{'ΔMAE':>8s} {'ΔSpear':>8s} {'ΔFP':>5s} {'ΔFP mass':>10s} "
        f"{'ΔDrop':>6s}"
    )
    print("-" * 95)

    for region, gdna, nrna, ss in conditions:
        cond_label = f"{region} gdna={gdna} nrna={nrna} ss={ss:.2f}"
        df = combined[
            (combined["region"] == region)
            & (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty or MAP_COL not in df.columns or VBEM_COL not in df.columns:
            continue

        m_map = compute_metrics(df, MAP_COL)
        m_vbem = compute_metrics(df, VBEM_COL)

        d_mae = m_vbem["mae"] - m_map["mae"]
        d_sp = m_vbem["spearman"] - m_map["spearman"]
        d_fp = m_vbem["fp_count"] - m_map["fp_count"]
        d_fpm = m_vbem["fp_mass"] - m_map["fp_mass"]
        d_drop = m_vbem["dropout"] - m_map["dropout"]

        print(
            f"  {cond_label:46s} "
            f"{d_mae:+8.2f} {d_sp:+8.4f} {d_fp:+5d} {d_fpm:+10.1f} "
            f"{d_drop:+6d}"
        )

    # ── SECTION 3: Top false positives per condition ─────────────
    print("\n" + "=" * 100)
    print("TOP FALSE POSITIVES (truth=0, per condition, top 10)")
    print("=" * 100)

    for region, gdna, nrna, ss in conditions:
        cond_label = f"{region} gdna={gdna} nrna={nrna} ss={ss:.2f}"
        df = combined[
            (combined["region"] == region)
            & (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty:
            continue
        zero = df[df["truth"] == 0]
        if zero.empty:
            continue

        # Only print conditions where at least one tool has FPs
        has_fp = False
        for tool in TOOLS:
            if tool in zero.columns and (zero[tool] > FP_THRESHOLD).any():
                has_fp = True
                break
        if not has_fp:
            continue

        print(f"\n  {cond_label}")
        # Sort by max of MAP/VBEM
        zero = zero.copy()
        fp_cols = [c for c in [MAP_COL, VBEM_COL] if c in zero.columns]
        zero["_max_pred"] = zero[fp_cols].max(axis=1)
        top = zero.nlargest(10, "_max_pred")

        for _, row in top.iterrows():
            tid = str(row["transcript_id"])[:30]
            parts = [f"{tid:30s}"]
            for tool in TOOLS:
                if tool in row.index:
                    short = tool.replace("hulkrna_", "").replace("_oracle", "")
                    parts.append(f"{short}={row[tool]:7.1f}")
            print("    " + "  ".join(parts))

    # ── SECTION 4: Negative control / gDNA scenario analysis ─────
    print("\n" + "=" * 100)
    print("NEGATIVE CONTROL ANALYSIS: gDNA stress scenarios")
    print("=" * 100)
    print("(Focus: Do zero-truth transcripts stay near zero under gDNA contamination?)\n")

    for region in REGIONS:
        for nrna in NRNA_LABELS:
            for ss in SS_VALUES:
                # Compare across gDNA levels
                rows_per_gdna = {}
                for gdna in GDNA_LABELS:
                    df = combined[
                        (combined["region"] == region)
                        & (combined["gdna_label"] == gdna)
                        & (combined["nrna_label"] == nrna)
                        & (combined["ss"] == ss)
                    ]
                    if not df.empty:
                        rows_per_gdna[gdna] = df

                if len(rows_per_gdna) < 2:
                    continue

                group = f"{region} nrna={nrna} ss={ss:.2f}"
                print(f"  {group}:")
                print(
                    f"    {'gDNA':8s} {'Tool':10s} "
                    f"{'FP>1':>5s} {'FP>5':>5s} {'FP>25':>6s} "
                    f"{'max(pred)':>10s} {'FP mass':>8s}"
                )

                for gdna in GDNA_LABELS:
                    if gdna not in rows_per_gdna:
                        continue
                    df = rows_per_gdna[gdna]
                    zero = df[df["truth"] == 0]
                    if zero.empty:
                        continue

                    for tool in [MAP_COL, VBEM_COL]:
                        if tool not in zero.columns:
                            continue
                        short = tool.replace("hulkrna_", "").replace("_oracle", "")
                        fp1 = (zero[tool] > 1).sum()
                        fp5 = (zero[tool] > 5).sum()
                        fp25 = (zero[tool] > 25).sum()
                        max_pred = zero[tool].max()
                        fp_mass = zero.loc[zero[tool] > 1, tool].sum()
                        print(
                            f"    {gdna:8s} {short:10s} "
                            f"{fp1:5d} {fp5:5d} {fp25:6d} "
                            f"{max_pred:10.1f} {fp_mass:8.1f}"
                        )
                print()

    # ── SECTION 5: Expressed transcript accuracy under stress ────
    print("\n" + "=" * 100)
    print("EXPRESSED TRANSCRIPT ACCURACY UNDER STRESS")
    print("=" * 100)
    print("(Focus: Does VBEM hurt expressed transcripts under gDNA/nRNA?)\n")

    for region in REGIONS:
        for gdna in GDNA_LABELS:
            for nrna in NRNA_LABELS:
                for ss in SS_VALUES:
                    df = combined[
                        (combined["region"] == region)
                        & (combined["gdna_label"] == gdna)
                        & (combined["nrna_label"] == nrna)
                        & (combined["ss"] == ss)
                    ]
                    if df.empty:
                        continue
                    expressed = df[df["truth"] > 0]
                    if len(expressed) < 3:
                        continue

                    cond = f"{region} gdna={gdna} nrna={nrna} ss={ss:.2f}"
                    map_err = np.abs(expressed[MAP_COL] - expressed["truth"])
                    vbem_err = np.abs(expressed[VBEM_COL] - expressed["truth"])

                    improved = (vbem_err < map_err - 0.01).sum()
                    same = ((vbem_err >= map_err - 0.01) & (vbem_err <= map_err + 0.01)).sum()
                    worse = (vbem_err > map_err + 0.01).sum()

                    print(
                        f"  {cond:48s}  "
                        f"VBEM better={improved:3d}  same={same:3d}  worse={worse:3d}  "
                        f"MAP_MAE={map_err.mean():7.2f}  VBEM_MAE={vbem_err.mean():7.2f}"
                    )

    # ── SECTION 6: Aggregate across all conditions ───────────────
    print("\n" + "=" * 100)
    print("AGGREGATE ACROSS ALL CONDITIONS")
    print("=" * 100)

    agg_map = compute_metrics(combined, MAP_COL)
    agg_vbem = compute_metrics(combined, VBEM_COL)

    print(f"\n  {'Metric':20s} {'MAP':>12s} {'VBEM':>12s} {'Delta':>12s}")
    print(f"  {'-'*60}")
    for metric in ["mae", "rmse", "spearman", "pearson", "fp_count", "fp_mass", "dropout"]:
        m = agg_map[metric]
        v = agg_vbem[metric]
        d = v - m
        fmt = ".2f" if metric in ("mae", "rmse", "fp_mass") else (
            ".4f" if metric in ("spearman", "pearson") else "d"
        )
        print(f"  {metric:20s} {m:>12{fmt}} {v:>12{fmt}} {d:>+12{fmt}}")

    # ── SECTION 7: Decision matrix ───────────────────────────────
    print("\n" + "=" * 100)
    print("DECISION MATRIX: Should VBEM replace MAP as default?")
    print("=" * 100)

    # Count conditions where VBEM is better/worse
    n_better_fp = 0
    n_worse_fp = 0
    n_better_mae = 0
    n_worse_mae = 0
    n_more_dropout = 0
    n_worse_neg_ctrl = 0
    n_conditions = 0

    for region, gdna, nrna, ss in conditions:
        df = combined[
            (combined["region"] == region)
            & (combined["gdna_label"] == gdna)
            & (combined["nrna_label"] == nrna)
            & (combined["ss"] == ss)
        ]
        if df.empty:
            continue
        n_conditions += 1

        m_map = compute_metrics(df, MAP_COL)
        m_vbem = compute_metrics(df, VBEM_COL)

        if m_vbem["fp_count"] < m_map["fp_count"]:
            n_better_fp += 1
        elif m_vbem["fp_count"] > m_map["fp_count"]:
            n_worse_fp += 1

        if m_vbem["mae"] < m_map["mae"] - 0.1:
            n_better_mae += 1
        elif m_vbem["mae"] > m_map["mae"] + 0.1:
            n_worse_mae += 1

        if m_vbem["dropout"] > m_map["dropout"]:
            n_more_dropout += 1

        # Negative control check: worst FP under gDNA stress
        zero = df[df["truth"] == 0]
        if not zero.empty and gdna != "none":
            map_max = zero[MAP_COL].max()
            vbem_max = zero[VBEM_COL].max()
            if vbem_max > map_max * 1.5:  # VBEM 50% worse on worst FP
                n_worse_neg_ctrl += 1

    print(f"\n  Total conditions tested: {n_conditions}")
    print(f"  VBEM fewer FPs:    {n_better_fp:3d} / {n_conditions}")
    print(f"  VBEM more FPs:     {n_worse_fp:3d} / {n_conditions}")
    print(f"  VBEM better MAE:   {n_better_mae:3d} / {n_conditions}")
    print(f"  VBEM worse MAE:    {n_worse_mae:3d} / {n_conditions}")
    print(f"  VBEM more dropout: {n_more_dropout:3d} / {n_conditions}")
    print(f"  VBEM worse neg ctrl (gDNA stress): {n_worse_neg_ctrl}")

    verdict = "INCONCLUSIVE"
    if n_worse_neg_ctrl > 0:
        verdict = "REJECT — VBEM breaks negative controls under gDNA stress"
    elif n_worse_fp > 0:
        verdict = "CAUTION — VBEM increases FPs in some conditions"
    elif n_more_dropout > 0:
        verdict = "CAUTION — VBEM increases dropout in some conditions"
    elif n_better_fp > n_conditions // 2 and n_worse_mae == 0:
        verdict = "ACCEPT — VBEM consistently fewer FPs with no MAE regression"
    elif n_better_fp > 0 and n_worse_mae == 0:
        verdict = "LEAN ACCEPT — VBEM helps FPs in some conditions, no regression"
    print(f"\n  >>> VERDICT: {verdict}")
    print()


if __name__ == "__main__":
    main()
