#!/usr/bin/env python
"""Compare sweep results before and after gDNA double-counting fix.

Usage:
    python scripts/analyze_sweep_comparison.py \
        /Users/mkiyer/Downloads/hulkrna_runs/sweep_full_grid.tsv \
        /Users/mkiyer/Downloads/hulkrna_runs/sweep_full_grid_v2.tsv
"""

import sys
import pandas as pd
import numpy as np

pd.set_option("display.max_columns", 40)
pd.set_option("display.width", 200)
pd.set_option("display.float_format", "{:.4f}".format)

v1_path = sys.argv[1]
v2_path = sys.argv[2]

v1 = pd.read_csv(v1_path, sep="\t")
v2 = pd.read_csv(v2_path, sep="\t")

# Exclude degenerate mrna=0 cases (nrna_abundance override only applies when
# transcript abundance > 0, so mrna=0 → nRNA and gDNA levels meaningless)
v1 = v1[v1["mrna"] > 0].copy()
v2 = v2[v2["mrna"] > 0].copy()

print(f"V1 rows: {len(v1)}, V2 rows: {len(v2)}")
print()

# Merge on grid keys for paired comparison
merged = v1.merge(v2, on=["gdna", "nrna", "mrna", "ss"], suffixes=("_v1", "_v2"))
print(f"Matched rows: {len(merged)}")
print()

# ── SECTION 1: Overall error distribution ──
print("=" * 70)
print("SECTION 1: Overall error distribution (|total_rna_rel_err|)")
print("=" * 70)
for label, col in [("V1 (before fix)", "total_rna_rel_err_v1"),
                    ("V2 (after fix)",  "total_rna_rel_err_v2")]:
    abs_err = merged[col].abs()
    pct_gt1 = (abs_err > 0.01).mean() * 100
    pct_gt5 = (abs_err > 0.05).mean() * 100
    pct_gt10 = (abs_err > 0.10).mean() * 100
    print(f"\n  {label}:")
    print(f"    Mean |err|      = {abs_err.mean():.4f}")
    print(f"    Median |err|    = {abs_err.median():.4f}")
    print(f"    Max |err|       = {abs_err.max():.4f}")
    print(f"    % > 1% error    = {pct_gt1:.1f}%")
    print(f"    % > 5% error    = {pct_gt5:.1f}%")
    print(f"    % > 10% error   = {pct_gt10:.1f}%")

# ── SECTION 2: Phantom gDNA (gdna=0 cases) ──
print()
print("=" * 70)
print("SECTION 2: Phantom gDNA (gdna=0 cases, false gDNA > 1 frag)")
print("=" * 70)
g0 = merged[merged["gdna"] == 0]
for ss_val in sorted(g0["ss"].unique()):
    subset = g0[g0["ss"] == ss_val]
    for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
        col = f"gdna_pipe{sfx}"
        n_phantom = (subset[col] > 1.0).sum()
        mean_phantom = subset[col].mean()
        max_phantom = subset[col].max()
        print(f"  SS={ss_val} {ver}: phantom_cases={n_phantom}/{len(subset)}  "
              f"mean={mean_phantom:.1f}  max={max_phantom:.1f}")

# ── SECTION 3: gDNA accuracy when gdna > 0 ──
print()
print("=" * 70)
print("SECTION 3: gDNA accuracy (gdna > 0) — relative error")
print("=" * 70)
g_pos = merged[merged["gdna"] > 0]
for ss_val in sorted(g_pos["ss"].unique()):
    subset = g_pos[g_pos["ss"] == ss_val]
    for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
        col = f"gdna_rel_err{sfx}"
        vals = subset[col]
        print(f"  SS={ss_val} {ver}: mean_rel_err={vals.mean():.4f}  "
              f"median={vals.median():.4f}  min={vals.min():.4f}  max={vals.max():.4f}")

# ── SECTION 4: gDNA over-absorption by gdna level at SS=1.0 ──
print()
print("=" * 70)
print("SECTION 4: gDNA over-absorption by level (SS=1.0)")
print("=" * 70)
ss10 = merged[merged["ss"] == 1.0]
for g_level in [20, 50, 100, 200, 500]:
    subset = ss10[ss10["gdna"] == g_level]
    if len(subset) == 0:
        continue
    v1_err = subset["gdna_rel_err_v1"].mean()
    v2_err = subset["gdna_rel_err_v2"].mean()
    print(f"  g={g_level:>3d}: V1 mean_rel_err={v1_err:+.4f}  V2 mean_rel_err={v2_err:+.4f}  "
          f"Δ={v2_err - v1_err:+.4f}")

# ── SECTION 5: nRNA accuracy ──
print()
print("=" * 70)
print("SECTION 5: nRNA accuracy (nrna > 0)")
print("=" * 70)
n_pos = merged[merged["nrna"] > 0]
for ss_val in sorted(n_pos["ss"].unique()):
    subset = n_pos[n_pos["ss"] == ss_val]
    for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
        col = f"nrna_rel_err{sfx}"
        vals = subset[col]
        print(f"  SS={ss_val} {ver}: mean_rel_err={vals.mean():.4f}  "
              f"median={vals.median():.4f}  min={vals.min():.4f}  max={vals.max():.4f}")

# ── SECTION 6: Total RNA at gdna=0 ──
print()
print("=" * 70)
print("SECTION 6: Total RNA accuracy at gdna=0 (mRNA+nRNA conservation)")
print("=" * 70)
g0_valid = merged[(merged["gdna"] == 0)]
for ss_val in sorted(g0_valid["ss"].unique()):
    subset = g0_valid[g0_valid["ss"] == ss_val]
    for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
        col = f"total_rna_rel_err{sfx}"
        vals = subset[col].abs()
        print(f"  SS={ss_val} {ver}: mean_abs_err={vals.mean():.4f}  max={vals.max():.4f}")

# ── SECTION 7: Negative control t2 ──
print()
print("=" * 70)
print("SECTION 7: Negative control t2 (should be 0)")
print("=" * 70)
for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
    col = f"t2_observed{sfx}"
    n_pos_t2 = (merged[col] > 0.5).sum()
    max_t2 = merged[col].max()
    mean_t2 = merged[col].mean()
    print(f"  {ver}: cases_with_t2>0.5 = {n_pos_t2}/{len(merged)}  "
          f"mean={mean_t2:.2f}  max={max_t2:.2f}")

# ── SECTION 8: Pure mRNA cases (gdna=0, nrna=0) ──
print()
print("=" * 70)
print("SECTION 8: Pure mRNA cases (gdna=0, nrna=0)")
print("=" * 70)
pure = merged[(merged["gdna"] == 0) & (merged["nrna"] == 0)]
for ss_val in sorted(pure["ss"].unique()):
    subset = pure[pure["ss"] == ss_val]
    for ver, sfx in [("V1", "_v1"), ("V2", "_v2")]:
        gdna_col = f"gdna_pipe{sfx}"
        rna_col = f"total_rna_rel_err{sfx}"
        t1_col = f"t1_rel_err{sfx}"
        print(f"  SS={ss_val} {ver}: false_gdna={subset[gdna_col].mean():.1f}  "
              f"total_rna_err={subset[rna_col].abs().mean():.4f}  "
              f"t1_rel_err={subset[t1_col].abs().mean():.4f}")

# ── SECTION 9: Worst cases comparison ──
print()
print("=" * 70)
print("SECTION 9: Top 10 worst |total_rna_rel_err| — V1 vs V2")
print("=" * 70)
merged["v1_abs"] = merged["total_rna_rel_err_v1"].abs()
merged["v2_abs"] = merged["total_rna_rel_err_v2"].abs()
worst_v1 = merged.nlargest(10, "v1_abs")[
    ["gdna", "nrna", "mrna", "ss", "total_rna_rel_err_v1", "total_rna_rel_err_v2",
     "gdna_pipe_v1", "gdna_pipe_v2"]]
print("\nWorst V1 cases:")
print(worst_v1.to_string(index=False))

worst_v2 = merged.nlargest(10, "v2_abs")[
    ["gdna", "nrna", "mrna", "ss", "total_rna_rel_err_v1", "total_rna_rel_err_v2",
     "gdna_pipe_v1", "gdna_pipe_v2"]]
print("\nWorst V2 cases:")
print(worst_v2.to_string(index=False))

# ── SECTION 10: Big improvement summary ──
print()
print("=" * 70)
print("SECTION 10: Improvement summary")
print("=" * 70)
merged["improvement"] = merged["v1_abs"] - merged["v2_abs"]
improved = (merged["improvement"] > 0.001).sum()
worsened = (merged["improvement"] < -0.001).sum()
unchanged = len(merged) - improved - worsened
print(f"  Cases improved (> 0.1% less total RNA error):   {improved} ({improved/len(merged)*100:.1f}%)")
print(f"  Cases worsened (> 0.1% more total RNA error):   {worsened} ({worsened/len(merged)*100:.1f}%)")
print(f"  Cases unchanged (within 0.1%):                  {unchanged} ({unchanged/len(merged)*100:.1f}%)")
print()
mean_v1 = merged["v1_abs"].mean()
mean_v2 = merged["v2_abs"].mean()
print(f"  Mean |total_rna_rel_err|: V1={mean_v1:.4f}  V2={mean_v2:.4f}  Δ={mean_v2-mean_v1:+.4f}")
print()

# Breakdown by SS
for ss_val in sorted(merged["ss"].unique()):
    ss_sub = merged[merged["ss"] == ss_val]
    m1 = ss_sub["v1_abs"].mean()
    m2 = ss_sub["v2_abs"].mean()
    print(f"  SS={ss_val}: V1={m1:.4f}  V2={m2:.4f}  Δ={m2-m1:+.4f}")
