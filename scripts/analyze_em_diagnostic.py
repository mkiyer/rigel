#!/usr/bin/env python3
"""Analyze EM diagnostic benchmark results: MAP vs MAP_PRUNE vs VBEM."""

import sys
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr

DATA_DIR = "/Users/mkiyer/Downloads/hulkrna_runs/bench_em_diagnostic"
REGIONS = ["FGFR2", "BRCA1"]
TOOLS = [
    "hulkrna_map_oracle",
    "hulkrna_map_prune_oracle",
    "hulkrna_vbem_oracle",
    "salmon",
    "kallisto",
]

frames = []
for region in REGIONS:
    csv = f"{DATA_DIR}/{region}/gdna_none_nrna_none_ss_1.00/per_transcript_counts.csv"
    df = pd.read_csv(csv)
    df["region"] = region
    frames.append(df)
all_df = pd.concat(frames, ignore_index=True)

print(f"Total transcripts: {len(all_df)}")
print(f"  With truth > 0: {(all_df['truth'] > 0).sum()}")
print(f"  With truth == 0: {(all_df['truth'] == 0).sum()}")
print()

# ── False Positive Analysis ──────────────────────────────────────
print("=" * 70)
print("FALSE POSITIVE ANALYSIS (truth == 0, predicted > threshold)")
print("=" * 70)
zero_truth = all_df[all_df["truth"] == 0].copy()
print(f"\nTranscripts with truth=0: {len(zero_truth)}")

for threshold in [0.5, 1.0, 2.0, 5.0]:
    print(f"\n  Threshold > {threshold}:")
    for tool in TOOLS:
        fp = (zero_truth[tool] > threshold).sum()
        total_fp_mass = zero_truth.loc[zero_truth[tool] > threshold, tool].sum()
        print(f"    {tool:35s}: {fp:3d} FP (total mass: {total_fp_mass:8.1f})")

# ── Per-region FP detail ─────────────────────────────────────────
for region in REGIONS:
    reg_zero = zero_truth[zero_truth["region"] == region]
    print(f"\n  {region} (truth=0 transcripts: {len(reg_zero)}):")
    for tool in TOOLS:
        fp = (reg_zero[tool] > 1.0).sum()
        total_mass = reg_zero.loc[reg_zero[tool] > 1.0, tool].sum()
        print(f"    {tool:35s}: {fp:3d} FP>1 (mass: {total_mass:8.1f})")

# ── Top FP transcripts comparison ────────────────────────────────
print("\n" + "=" * 70)
print("TOP 15 FALSE POSITIVES (sorted by MAP prediction)")
print("=" * 70)
zero_top = zero_truth.nlargest(15, "hulkrna_map_oracle")
for _, row in zero_top.iterrows():
    tid = row["transcript_id"][:25]
    print(
        f"  {tid:25s} {row['region']:6s}  "
        f"MAP={row['hulkrna_map_oracle']:7.1f}  "
        f"PRUNE={row['hulkrna_map_prune_oracle']:7.1f}  "
        f"VBEM={row['hulkrna_vbem_oracle']:7.1f}  "
        f"SAL={row['salmon']:7.1f}  "
        f"KAL={row['kallisto']:7.1f}"
    )

# ── Dropout Analysis ─────────────────────────────────────────────
print("\n" + "=" * 70)
print("DROPOUT ANALYSIS (truth > 0, predicted < 0.5)")
print("=" * 70)
expressed = all_df[all_df["truth"] > 0].copy()
print(f"\nExpressed transcripts: {len(expressed)}")

for tool in TOOLS:
    dropout = (expressed[tool] < 0.5).sum()
    dropout_mass = expressed.loc[expressed[tool] < 0.5, "truth"].sum()
    print(f"  {tool:35s}: {dropout:3d} dropouts (truth mass: {dropout_mass:8.1f})")

# ── Overcount Analysis ───────────────────────────────────────────
print("\n" + "=" * 70)
print("OVERCOUNT ANALYSIS (truth > 0, predicted > 2× truth)")
print("=" * 70)

for tool in TOOLS:
    oc = expressed[expressed[tool] > 2 * expressed["truth"]]
    print(f"  {tool:35s}: {len(oc):3d} overcounted (>2× truth)")

# ── Summary Statistics ───────────────────────────────────────────
print("\n" + "=" * 70)
print("AGGREGATE METRICS (all transcripts)")
print("=" * 70)
print(f"\n{'Tool':35s} {'MAE':>8s} {'RMSE':>8s} {'Pearson':>8s} {'Spearman':>8s} {'FP>1':>6s} {'Dropout':>8s}")
print("-" * 85)

for tool in TOOLS:
    err = all_df[tool] - all_df["truth"]
    mae = np.abs(err).mean()
    rmse = np.sqrt((err ** 2).mean())
    
    r_pe, _ = pearsonr(all_df["truth"], all_df[tool])
    r_sp, _ = spearmanr(all_df["truth"], all_df[tool])
    
    fp = (zero_truth[tool] > 1.0).sum()
    dropout = ((expressed["truth"] > 0) & (expressed[tool] < 0.5)).sum()
    
    print(f"  {tool:35s} {mae:8.2f} {rmse:8.2f} {r_pe:8.4f} {r_sp:8.4f} {fp:6d} {dropout:8d}")

# ── Per-region metrics ───────────────────────────────────────────
for region in REGIONS:
    reg = all_df[all_df["region"] == region]
    reg_zero = reg[reg["truth"] == 0]
    reg_exp = reg[reg["truth"] > 0]
    
    print(f"\n  {region} ({len(reg)} tx, {len(reg_exp)} expressed, {len(reg_zero)} zero-truth):")
    print(f"  {'Tool':35s} {'MAE':>8s} {'Spearman':>8s} {'FP>1':>6s} {'FP mass':>10s}")
    print(f"  {'-'*75}")
    
    for tool in TOOLS:
        err = reg[tool] - reg["truth"]
        mae = np.abs(err).mean()
        r_sp, _ = spearmanr(reg["truth"], reg[tool])
        fp = (reg_zero[tool] > 1.0).sum()
        fp_mass = reg_zero.loc[reg_zero[tool] > 1.0, tool].sum()
        print(f"  {tool:35s} {mae:8.2f} {r_sp:8.4f} {fp:6d} {fp_mass:10.1f}")

# ── VBEM vs MAP improvement on expressed transcripts ─────────────
print("\n" + "=" * 70)
print("VBEM vs MAP: PER-TRANSCRIPT IMPROVEMENTS (expressed only)")
print("=" * 70)

for region in REGIONS:
    reg_exp = all_df[(all_df["region"] == region) & (all_df["truth"] > 0)]
    map_err = np.abs(reg_exp["hulkrna_map_oracle"] - reg_exp["truth"])
    vbem_err = np.abs(reg_exp["hulkrna_vbem_oracle"] - reg_exp["truth"])
    
    improved = (vbem_err < map_err).sum()
    same = (vbem_err == map_err).sum()
    worse = (vbem_err > map_err).sum()
    
    print(f"\n  {region} (expressed: {len(reg_exp)}):")
    print(f"    VBEM better: {improved} transcripts")
    print(f"    Same:        {same} transcripts")
    print(f"    VBEM worse:  {worse} transcripts")
    print(f"    MAP MAE:     {map_err.mean():.2f}")
    print(f"    VBEM MAE:    {vbem_err.mean():.2f}")
    
    # Worst VBEM regressions
    reg_exp = reg_exp.copy()
    reg_exp["map_ae"] = map_err.values
    reg_exp["vbem_ae"] = vbem_err.values
    reg_exp["delta"] = (vbem_err - map_err).values
    
    worst = reg_exp.nlargest(5, "delta")
    if not worst.empty:
        print(f"    Worst VBEM regressions:")
        for _, row in worst.iterrows():
            tid = row["transcript_id"][:25]
            print(
                f"      {tid:25s} truth={row['truth']:.0f}  "
                f"MAP={row['hulkrna_map_oracle']:.1f}  "
                f"VBEM={row['hulkrna_vbem_oracle']:.1f}  "
                f"Δ_AE={row['delta']:+.1f}"
            )

# ── MAP_PRUNE investigation ──────────────────────────────────────
print("\n" + "=" * 70)
print("MAP_PRUNE INVESTIGATION: Why no improvement?")
print("=" * 70)
# Check if map and map_prune are truly identical
for region in REGIONS:
    reg = all_df[all_df["region"] == region]
    map_vals = reg["hulkrna_map_oracle"].values
    prune_vals = reg["hulkrna_map_prune_oracle"].values
    diff = np.abs(map_vals - prune_vals)
    n_diff = (diff > 0.01).sum()
    max_diff = diff.max()
    print(f"\n  {region}: {n_diff} transcripts differ by >0.01 (max diff: {max_diff:.6f})")
    if n_diff > 0:
        idx = np.argsort(-diff)[:10]
        for i in idx[:5]:
            row = reg.iloc[i]
            print(
                f"    {row['transcript_id'][:25]:25s} truth={row['truth']:.0f}  "
                f"MAP={map_vals[i]:.4f}  PRUNE={prune_vals[i]:.4f}  "
                f"diff={diff[i]:.6f}"
            )
    else:
        print("    → Pruning threshold 2×prior is too lenient;")
        print("      zombie transcripts accumulate alpha >> 2×prior from shared exons")
        # Show alpha values for top FP transcripts
        reg_zero = reg[reg["truth"] == 0]
        top_fp = reg_zero.nlargest(5, "hulkrna_map_oracle")
        for _, row in top_fp.iterrows():
            tid = row["transcript_id"][:25]
            pred = row["hulkrna_map_oracle"]
            print(
                f"    {tid:25s} truth=0 MAP={pred:.1f}"
                f"  (alpha ≫ 2×0.01 → not pruned)"
            )
