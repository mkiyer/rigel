#!/usr/bin/env python
"""Focused analysis of sweep grid — excludes degenerate mrna=0 cases."""

import pandas as pd
import numpy as np

TSV = "/Users/mkiyer/Downloads/hulkrna_runs/sweep_full_grid.tsv"

df = pd.read_csv(TSV, sep="\t", comment="#")
for c in df.columns:
    if df[c].dtype == object:
        df[c] = pd.to_numeric(df[c], errors="coerce")

# Exclude mrna=0 (degenerate: nrna_abundance override doesn't apply when
# t.abundance=0, so the simulator falls back to uniform → t2 gets
# garbage counts.  Not a pipeline bug, just a test setup artifact.)
df = df[df["mrna"] > 0].copy()
print(f"Rows with mrna>0: {len(df)}")

# Computed columns
df["total_rna_gt"] = df["mrna_gt"] + df["nrna_gt"]
df["total_rna_pipe"] = df["mrna_pipe"] + df["nrna_pipe"]
df["total_rna_signed_rel"] = np.where(
    df["total_rna_gt"] > 0,
    (df["total_rna_pipe"] - df["total_rna_gt"]) / df["total_rna_gt"],
    np.nan,
)
df["gdna_signed_rel"] = np.where(
    df["gdna_gt"] > 0,
    (df["gdna_pipe"] - df["gdna_gt"]) / df["gdna_gt"],
    np.nan,
)
df["nrna_signed_rel"] = np.where(
    df["nrna_gt"] > 0,
    (df["nrna_pipe"] - df["nrna_gt"]) / df["nrna_gt"],
    np.nan,
)

# ═══════════════════════════════════════════════════════════════════
# FINDING 1: PHANTOM gDNA (false positive gDNA when none exists)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 1: PHANTOM gDNA (gdna=0, nrna=0, mrna>0 — pure mRNA)")
print("═" * 70)

pure_mrna = df[(df["gdna"] == 0) & (df["nrna"] == 0)].copy()
print(f"\n{'SS':>6} {'N':>4} {'mean_false_gDNA':>16} {'as_%_of_total':>14} {'t1_mean_rel_err':>16}")
for ss in sorted(pure_mrna["ss"].unique()):
    sub = pure_mrna[pure_mrna["ss"] == ss]
    print(f"{ss:>6.2f} {len(sub):>4d} {sub['gdna_pipe'].mean():>16.1f} "
          f"{sub['gdna_pipe'].mean() / 10000:>14.3%} {sub['t1_rel_err'].mean():>16.4f}")

print("\nThis shows that at SS=0.9/0.95, ~1.5% of reads become phantom gDNA")
print("even when there is zero gDNA contamination. At SS=1.0 and SS=0.5 it's fine.")

# ═══════════════════════════════════════════════════════════════════
# FINDING 2: gDNA OVER-ABSORPTION (absorbs RNA as gDNA)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 2: gDNA OVER-ABSORPTION")
print("═" * 70)

has_gdna = df[df["gdna"] > 0].copy()
print("\nMean gDNA OVER-estimation (signed rel err) by gDNA × SS:")
pivot = has_gdna.pivot_table(index="gdna", columns="ss", values="gdna_signed_rel", aggfunc="mean")
print(pivot.to_string(float_format=lambda x: f"{x:+.1%}"))

print("\ngDNA is ALWAYS over-estimated. Lower gDNA → worse relative error.")
print("Even at SS=1.0, g=20 has +30% over-estimation, g=500 has +12%.")

# ═══════════════════════════════════════════════════════════════════
# FINDING 3: nRNA SYSTEMATIC UNDER-ESTIMATION
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 3: nRNA SYSTEMATIC UNDER-ESTIMATION (gdna=0, mrna>0)")
print("═" * 70)

has_nrna_no_gdna = df[(df["nrna"] > 0) & (df["gdna"] == 0)].copy()
print("\nMean nRNA signed rel err by nRNA × SS (no gDNA):")
pivot_n = has_nrna_no_gdna.pivot_table(index="nrna", columns="ss", values="nrna_signed_rel", aggfunc="mean")
print(pivot_n.to_string(float_format=lambda x: f"{x:+.1%}"))

print("\nnRNA is ALWAYS under-estimated even at SS=1.0.")
print("At SS=1.0, under-estimation is 17-24% — this is purely exonic nRNA leaking to mRNA.")

# Now also show total RNA at SS=1.0, gdna=0 — should be perfect
clean = has_nrna_no_gdna[has_nrna_no_gdna["ss"] == 1.0].copy()
print("\nBut total RNA (mRNA+nRNA) with SS=1.0, gdna=0:")
for nrna_level in sorted(clean["nrna"].unique()):
    sub = clean[clean["nrna"] == nrna_level]
    trel = sub["total_rna_signed_rel"]
    print(f"  nRNA={nrna_level}: mean_totalRNA_rel={trel.mean():+.4f}  max_abs={trel.abs().max():.4f}")

# ═══════════════════════════════════════════════════════════════════
# FINDING 4: TOTAL RNA BY SS AT gdna=0 (no gDNA confusion)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 4: TOTAL RNA ACCOUNTABILITY at gdna=0 by SS")
print("═" * 70)

no_gdna = df[df["gdna"] == 0].copy()
for ss in sorted(no_gdna["ss"].unique()):
    sub = no_gdna[no_gdna["ss"] == ss]
    trel = sub["total_rna_signed_rel"].dropna()
    print(f"\n  SS={ss} (N={len(trel)}): "
          f"mean={trel.mean():+.4f}  median={trel.median():+.4f}  "
          f"P5/P95=[{trel.quantile(0.05):+.4f}, {trel.quantile(0.95):+.4f}]  "
          f"max_abs={trel.abs().max():.4f}")

# ═══════════════════════════════════════════════════════════════════
# FINDING 5: DETAILED PHANTOM gDNA INVESTIGATION
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 5: PHANTOM gDNA by (gdna, nrna, SS) — gdna=0 only, mrna=100")
print("═" * 70)

phantom_detail = df[(df["gdna"] == 0) & (df["mrna"] == 100)].copy()
print(f"\n{'nrna':>6} {'ss':>6} {'false_gDNA':>10} {'nrna_gt':>8} {'nrna_pipe':>10} "
      f"{'mrna_gt':>8} {'mrna_pipe':>10} {'total_rna_gt':>13} {'total_rna_pipe':>15}")
for _, r in phantom_detail.sort_values(["nrna", "ss"]).iterrows():
    print(f"{int(r['nrna']):>6d} {r['ss']:>6.2f} {r['gdna_pipe']:>10.1f} "
          f"{int(r['nrna_gt']):>8d} {r['nrna_pipe']:>10.1f} "
          f"{int(r['mrna_gt']):>8d} {r['mrna_pipe']:>10.1f} "
          f"{int(r['total_rna_gt']):>13d} {r['total_rna_pipe']:>15.1f}")

# ═══════════════════════════════════════════════════════════════════
# FINDING 6: t2 NEGATIVE CONTROL (excluding mrna=0)
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 6: NEGATIVE CONTROL t2 (mrna>0)")
print("═" * 70)

t2_pos = df[df["t2_observed"] > 0.5]
print(f"\n{len(t2_pos)} / {len(df)} cases have t2_observed > 0.5")
if len(t2_pos) > 0:
    print(f"  max: {t2_pos['t2_observed'].max():.1f}")
    print(f"  mean: {t2_pos['t2_observed'].mean():.1f}")
    worst = t2_pos.nlargest(15, "t2_observed")
    print("\nWorst 15:")
    print(worst[["gdna", "nrna", "mrna", "ss", "t2_expected", "t2_observed",
                 "gdna_gt", "gdna_pipe"]].to_string(index=False))

# ═══════════════════════════════════════════════════════════════════
# FINDING 7: PERFECT CONDITIONS — SS=1.0, gdna=0 
# ═══════════════════════════════════════════════════════════════════
print("\n" + "═" * 70)
print("FINDING 7: IDEAL CONDITIONS (SS=1.0, gdna=0)")
print("═" * 70)

ideal = df[(df["ss"] == 1.0) & (df["gdna"] == 0)].copy()
print(f"\nN={len(ideal)} cases")
print(f"  max |t1_rel_err| = {ideal['t1_rel_err'].max():.4f}")
print(f"  mean |t1_rel_err| = {ideal['t1_rel_err'].mean():.4f}")
print(f"  max |total_rna_signed_rel| = {ideal['total_rna_signed_rel'].abs().max():.4f}")
print(f"  max phantom gDNA = {ideal['gdna_pipe'].max():.1f}")
print(f"  max t2 false positive = {ideal['t2_observed'].max():.1f}")
print("\nAll ideal cases:")
for _, r in ideal.sort_values(["mrna", "nrna"]).iterrows():
    if r["nrna_gt"] > 0:
        nrna_tag = f" nrna={int(r['nrna_gt'])}→{r['nrna_pipe']:.0f}"
    else:
        nrna_tag = ""
    print(f"  m{int(r['mrna']):>3d} n{int(r['nrna']):>3d}: "
          f"t1={int(r['t1_expected'])}→{r['t1_observed']:.0f} "
          f"(err={r['t1_rel_err']:.3f})"
          f"{nrna_tag}  "
          f"totalRNA={int(r['total_rna_gt'])}→{r['total_rna_pipe']:.0f} "
          f"gDNA={r['gdna_pipe']:.0f}")
