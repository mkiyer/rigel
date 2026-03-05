#!/usr/bin/env python
"""Analyze the full sweep grid TSV and print error summary tables."""

import sys
import pandas as pd
import numpy as np

TSV = "/Users/mkiyer/Downloads/hulkrna_runs/sweep_full_grid.tsv"

df = pd.read_csv(TSV, sep="\t", comment="#")
# Replace inf with NaN for numeric ops
for c in df.columns:
    if df[c].dtype == object:
        df[c] = pd.to_numeric(df[c], errors="coerce")

print(f"Total rows: {len(df)}")
print()

# =================================================================
# 1. OVERALL ERROR DISTRIBUTION
# =================================================================
print("=" * 70)
print("1. OVERALL ERROR DISTRIBUTION")
print("=" * 70)

# Total RNA relative error (signed)
df["total_rna_signed_rel"] = np.where(
    df["total_rna_gt"] > 0,
    (df["total_rna_pipe"] - df["total_rna_gt"]) / df["total_rna_gt"],
    np.nan,
)

for metric, label in [
    ("total_rna_rel_err", "Total RNA (mRNA+nRNA) rel err"),
    ("mrna_rel_err", "mRNA signed rel err"),
    ("t1_rel_err", "t1 rel err"),
]:
    vals = df[metric].replace([np.inf, -np.inf], np.nan).dropna()
    print(f"\n{label} (N={len(vals)}):")
    print(f"  mean={vals.mean():.4f}  median={vals.median():.4f}")
    print(f"  min={vals.min():.4f}  max={vals.max():.4f}")
    print(f"  P5={vals.quantile(0.05):.4f}  P95={vals.quantile(0.95):.4f}")
    # How many > 5%, > 10%, > 20%?
    for thresh in [0.01, 0.05, 0.10, 0.20, 0.50]:
        n_above = (vals.abs() > thresh).sum()
        print(f"  |err| > {thresh:.0%}: {n_above} ({n_above/len(vals):.1%})")

# =================================================================
# 2. FALSE gDNA (phantom gDNA when gdna_gt=0)
# =================================================================
print("\n" + "=" * 70)
print("2. FALSE gDNA: cases with gdna=0 but gdna_pipe > 0")
print("=" * 70)

no_gdna = df[df["gdna"] == 0].copy()
no_gdna["false_gdna_frac"] = no_gdna["gdna_pipe"] / no_gdna["n_frags_pipeline"]
phantom = no_gdna[no_gdna["gdna_pipe"] > 1.0]
print(f"\n{len(phantom)} / {len(no_gdna)} cases have phantom gDNA > 1 frag")
if len(phantom) > 0:
    print("\nPhantom gDNA by SS:")
    for ss in sorted(phantom["ss"].unique()):
        sub = phantom[phantom["ss"] == ss]
        print(f"  SS={ss}: N={len(sub)}, "
              f"mean_false_gdna={sub['gdna_pipe'].mean():.0f}, "
              f"max_false_gdna={sub['gdna_pipe'].max():.0f}, "
              f"mean_frac={sub['false_gdna_frac'].mean():.3f}")

# =================================================================
# 3. ERROR HEATMAP: t1_rel_err by (gdna, nrna, mrna) at each SS
# =================================================================
print("\n" + "=" * 70)
print("3. t1 RELATIVE ERROR by SS level")
print("=" * 70)

for ss in sorted(df["ss"].unique()):
    sub = df[(df["ss"] == ss) & (df["mrna"] > 0)].copy()
    print(f"\n--- SS = {ss} ---")
    # Pivot: rows = gdna, cols = mrna, values = mean t1_rel_err across nrna levels
    pivot = sub.pivot_table(
        index="gdna", columns="mrna",
        values="t1_rel_err", aggfunc="mean",
    )
    print("Mean |t1_rel_err| by gDNA (rows) × mRNA (cols), averaged over nRNA:")
    print(pivot.to_string(float_format=lambda x: f"{x:.3f}"))

# =================================================================
# 4. TOTAL RNA ACCOUNTABILITY by SS
# =================================================================
print("\n" + "=" * 70)
print("4. TOTAL RNA (mRNA+nRNA) ACCOUNTABILITY by SS")
print("=" * 70)

has_rna = df[df["total_rna_gt"] > 0].copy()
has_rna["total_rna_signed"] = (has_rna["total_rna_pipe"] - has_rna["total_rna_gt"]) / has_rna["total_rna_gt"]

for ss in sorted(has_rna["ss"].unique()):
    sub = has_rna[has_rna["ss"] == ss]
    signed = sub["total_rna_signed"]
    print(f"\n  SS={ss} (N={len(sub)}): "
          f"mean={signed.mean():+.4f}  median={signed.median():+.4f}  "
          f"P5={signed.quantile(0.05):+.4f}  P95={signed.quantile(0.95):+.4f}")

# =================================================================
# 5. WORST CASES: top 20 by |total_rna_rel_err|
# =================================================================
print("\n" + "=" * 70)
print("5. WORST 30 CASES by |total RNA relative error|")
print("=" * 70)

has_rna_sorted = has_rna.copy()
has_rna_sorted["abs_total_rna_rel"] = has_rna_sorted["total_rna_rel_err"].abs()
worst = has_rna_sorted.nlargest(30, "abs_total_rna_rel")
cols_show = ["gdna", "nrna", "mrna", "ss",
             "mrna_gt", "mrna_pipe", "nrna_gt", "nrna_pipe",
             "gdna_gt", "gdna_pipe", "total_rna_gt", "total_rna_pipe",
             "total_rna_rel_err"]
print(worst[cols_show].to_string(index=False))

# =================================================================
# 6. NEGATIVE CONTROL: t2 false positives
# =================================================================
print("\n" + "=" * 70)
print("6. NEGATIVE CONTROL (t2 always abundance=0)")
print("=" * 70)

t2_pos = df[df["t2_observed"] > 0.5]
print(f"\n{len(t2_pos)} / {len(df)} cases have t2_observed > 0.5")
if len(t2_pos) > 0:
    print(f"  max t2_observed: {t2_pos['t2_observed'].max():.1f}")
    print(f"  mean t2_observed: {t2_pos['t2_observed'].mean():.1f}")
    # Show worst
    worst_t2 = t2_pos.nlargest(10, "t2_observed")
    print("\nWorst 10:")
    print(worst_t2[["gdna", "nrna", "mrna", "ss", "t2_observed"]].to_string(index=False))

# =================================================================
# 7. PURE mRNA ONLY (gdna=0, nrna=0): should be perfect
# =================================================================
print("\n" + "=" * 70)
print("7. PURE mRNA ONLY (gdna=0, nrna=0)")
print("=" * 70)

pure = df[(df["gdna"] == 0) & (df["nrna"] == 0) & (df["mrna"] > 0)].copy()
for ss in sorted(pure["ss"].unique()):
    sub = pure[pure["ss"] == ss]
    print(f"\n  SS={ss}: "
          f"mean_gdna_pipe={sub['gdna_pipe'].mean():.1f}  "
          f"t1_mean_rel_err={sub['t1_rel_err'].mean():.4f}  "
          f"max t1_abs_diff={sub['t1_abs_diff'].max():.0f}")

# =================================================================
# 8. gDNA OVER-ABSORPTION: gdna > 0 cases
# =================================================================
print("\n" + "=" * 70)
print("8. gDNA OVER-ABSORPTION: mean gdna_rel_err by gdna × ss")
print("=" * 70)

has_gdna = df[df["gdna"] > 0].copy()
pivot_gdna = has_gdna.pivot_table(
    index="gdna", columns="ss",
    values="gdna_rel_err", aggfunc="mean",
)
print(pivot_gdna.to_string(float_format=lambda x: f"{x:+.3f}"))

# =================================================================
# 9. nRNA ACCURACY: cases with nrna > 0 and mrna > 0
# =================================================================
print("\n" + "=" * 70)
print("9. nRNA ACCURACY: signed nrna_rel_err by nrna × ss (gdna=0, mrna>0)")
print("=" * 70)

has_nrna = df[(df["nrna"] > 0) & (df["gdna"] == 0) & (df["mrna"] > 0)].copy()
pivot_nrna = has_nrna.pivot_table(
    index="nrna", columns="ss",
    values="nrna_rel_err", aggfunc="mean",
)
print("Mean nrna_rel_err by nRNA abundance × SS:")
print(pivot_nrna.to_string(float_format=lambda x: f"{x:+.3f}"))

# =================================================================
# 10. CLEAN CASES (SS=1.0) ERROR SUMMARY
# =================================================================
print("\n" + "=" * 70)
print("10. CLEAN CASES: SS=1.0, any gDNA/nRNA/mRNA")
print("=" * 70)

clean = df[df["ss"] == 1.0].copy()
has_rna_clean = clean[clean["total_rna_gt"] > 0]
print(f"\nN={len(has_rna_clean)} cases with SS=1.0 and RNA present")
signed = (has_rna_clean["total_rna_pipe"] - has_rna_clean["total_rna_gt"]) / has_rna_clean["total_rna_gt"]
print(f"  Total RNA rel err: mean={signed.mean():+.4f}  "
      f"median={signed.median():+.4f}  "
      f"P5={signed.quantile(0.05):+.4f}  P95={signed.quantile(0.95):+.4f}  "
      f"max_abs={signed.abs().max():.4f}")
# Worst at SS=1.0
worst_clean = has_rna_clean.copy()
worst_clean["abs_tr"] = signed.abs()
worst10 = worst_clean.nlargest(10, "abs_tr")
print("\nWorst 10 at SS=1.0:")
print(worst10[["gdna", "nrna", "mrna", "mrna_gt", "mrna_pipe",
               "nrna_gt", "nrna_pipe", "gdna_gt", "gdna_pipe",
               "total_rna_gt", "total_rna_pipe"]].to_string(index=False))
