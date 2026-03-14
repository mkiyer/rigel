#!/usr/bin/env python
"""Comprehensive analysis of T+N+1 ablation sweep (v2 — simulator fix applied).

Sections:
  1. Dataset overview
  2. Accuracy by strand specificity
  3. Accuracy by gDNA fraction
  4. Cross-table: SS × gdna_fraction
  5. False positive analysis
  6. nRNA-gDNA confusion analysis
  7. Worst cases
  8. Fragment length recovery
  9. Weighted accuracy
  10. Edge cases: high nRNA / low mRNA and vice versa
"""
import numpy as np
import pandas as pd
from pathlib import Path

TSV = Path("sweep_results/tn1_ablation_v2/sweep_results.tsv")
SEP = "=" * 74


def main():
    df = pd.read_csv(TSV, sep="\t")
    print(f"Loaded {len(df)} runs, {df.shape[1]} columns")

    # Signed relative errors
    df["gdna_rel_signed"] = np.where(
        df.gdna_expected > 0,
        (df.gdna_observed - df.gdna_expected) / df.gdna_expected, np.nan)
    df["nrna_rel_signed"] = np.where(
        df.nrna_expected > 0,
        (df.nrna_observed - df.nrna_expected) / df.nrna_expected, np.nan)
    df["mrna_rel_signed"] = np.where(
        df.total_mrna_expected > 0,
        (df.total_mrna_observed - df.total_mrna_expected) / df.total_mrna_expected, np.nan)

    df["gdna_excess"] = df.gdna_observed - df.gdna_expected
    df["nrna_excess"] = df.nrna_observed - df.nrna_expected
    df["mrna_excess"] = df.total_mrna_observed - df.total_mrna_expected

    # ── Section 1: Overview ──
    print(f"\n{SEP}")
    print("1. DATASET OVERVIEW")
    print(SEP)
    print(f"  Runs: {len(df)}")
    print(f"  SS values: {sorted(df.strand_specificity.unique())}")
    print(f"  gdna_fraction values: {sorted(df.gdna_fraction.unique())}")
    print(f"  TA1 range: {df.TA1.min()} .. {df.TA1.max()}")
    print(f"  NTA1 range: {df.NTA1.min()} .. {df.NTA1.max()}")
    # Check for degenerate cases
    degen = df[(df.TA1 == 0) & (df.NTA1 == 0)]
    print(f"  TA1=0 & NTA1=0 runs (gDNA-only): {len(degen)}")
    print(f"  n_fragments_actual range: {df.n_fragments_actual.min()} .. {df.n_fragments_actual.max()}")

    # ── Section 2: Accuracy by SS ──
    print(f"\n{SEP}")
    print("2. ACCURACY BY STRAND SPECIFICITY")
    print(SEP)
    for ss in sorted(df.strand_specificity.unique()):
        sub = df[df.strand_specificity == ss]
        m = sub.mrna_rel_signed.dropna()
        n = sub.nrna_rel_signed.dropna()
        g = sub.gdna_rel_signed.dropna()
        print(f"\n  SS = {ss:.2f} ({len(sub)} runs)")
        if len(m) > 0:
            print(f"    mRNA: mean={m.mean():+.4f}  median={m.median():+.4f}  "
                  f"|mean|={m.abs().mean():.4f}  p10={m.quantile(0.1):+.4f}  p90={m.quantile(0.9):+.4f}")
        if len(n) > 0:
            print(f"    nRNA: mean={n.mean():+.4f}  median={n.median():+.4f}  "
                  f"|mean|={n.abs().mean():.4f}  p10={n.quantile(0.1):+.4f}  p90={n.quantile(0.9):+.4f}")
        if len(g) > 0:
            print(f"    gDNA: mean={g.mean():+.4f}  median={g.median():+.4f}  "
                  f"|mean|={g.abs().mean():.4f}  p10={g.quantile(0.1):+.4f}  p90={g.quantile(0.9):+.4f}")

    # ── Section 3: Accuracy by gDNA fraction ──
    print(f"\n{SEP}")
    print("3. ACCURACY BY gDNA FRACTION")
    print(SEP)
    for gf in sorted(df.gdna_fraction.unique()):
        sub = df[df.gdna_fraction == gf]
        m = sub.mrna_rel_signed.dropna()
        n = sub.nrna_rel_signed.dropna()
        g = sub.gdna_rel_signed.dropna()
        print(f"\n  gdna_fraction = {gf:.1f} ({len(sub)} runs)")
        if len(m) > 0:
            print(f"    mRNA: mean={m.mean():+.4f}  |mean|={m.abs().mean():.4f}")
        if len(n) > 0:
            print(f"    nRNA: mean={n.mean():+.4f}  |mean|={n.abs().mean():.4f}")
        if len(g) > 0:
            print(f"    gDNA: mean={g.mean():+.4f}  |mean|={g.abs().mean():.4f}")

    # ── Section 4: Cross-table SS × gDNA ──
    print(f"\n{SEP}")
    print("4. CROSS-TABLE: SS × gdna_fraction (mean signed relative error)")
    print(SEP)
    for component, col in [("mRNA", "mrna_rel_signed"),
                           ("nRNA", "nrna_rel_signed"),
                           ("gDNA", "gdna_rel_signed")]:
        print(f"\n  {component}:")
        gf_vals = sorted(df.gdna_fraction.unique())
        header = "  SS\\gf  " + "".join(f"{gf:>8.1f}" for gf in gf_vals)
        print(header)
        for ss in sorted(df.strand_specificity.unique()):
            row = f"  {ss:.2f}   "
            for gf in gf_vals:
                sub = df[(df.strand_specificity == ss) & (df.gdna_fraction == gf)]
                val = sub[col].dropna()
                if len(val) > 0:
                    row += f"{val.mean():+8.4f}"
                else:
                    row += "     N/A"
            print(row)

    # ── Section 5: False positives ──
    print(f"\n{SEP}")
    print("5. FALSE POSITIVE ANALYSIS")
    print(SEP)

    # A. nRNA FPs: NTA1=0, TA1>0
    nrna_fp_pool = df[(df.NTA1 == 0) & (df.TA1 > 0)]
    nrna_fp = nrna_fp_pool[nrna_fp_pool.nrna_observed > 1.0]
    print(f"\n  A. nRNA false positives (NTA1=0, TA1>0, nRNA_obs > 1):")
    print(f"     {len(nrna_fp)} / {len(nrna_fp_pool)} runs ({100*len(nrna_fp)/max(1,len(nrna_fp_pool)):.0f}%)")
    for ss in sorted(nrna_fp.strand_specificity.unique()):
        sub = nrna_fp[nrna_fp.strand_specificity == ss]
        print(f"     SS={ss:.2f}: {len(sub)} FPs, mean={sub.nrna_observed.mean():.1f}, "
              f"median={sub.nrna_observed.median():.1f}, max={sub.nrna_observed.max():.1f}")

    # B. gDNA FPs: gdna_frac=0
    has_rna = (df.TA1 > 0) | (df.NTA1 > 0)
    gdna_fp_pool = df[(df.gdna_fraction == 0) & has_rna]
    gdna_fp = gdna_fp_pool[gdna_fp_pool.gdna_observed > 1.0]
    print(f"\n  B. gDNA false positives (gdna_frac=0, has RNA, gDNA_obs > 1):")
    print(f"     {len(gdna_fp)} / {len(gdna_fp_pool)} runs ({100*len(gdna_fp)/max(1,len(gdna_fp_pool)):.0f}%)")
    for ss in sorted(gdna_fp.strand_specificity.unique()):
        sub = gdna_fp[gdna_fp.strand_specificity == ss]
        print(f"     SS={ss:.2f}: {len(sub)} FPs, mean gDNA_obs={sub.gdna_observed.mean():.1f}")

    # C. mRNA FPs: TA1=0 (TD=TE=0 always), NTA1>0 (simulator now gives 0 mRNA frags for TA1=0)
    mrna_fp_pool = df[(df.TA1 == 0) & (df.NTA1 > 0)]
    mrna_fp = mrna_fp_pool[mrna_fp_pool.total_mrna_observed > 1.0]
    print(f"\n  C. mRNA false positives (TA1=0, NTA1>0, mRNA_obs > 1):")
    print(f"     {len(mrna_fp)} / {len(mrna_fp_pool)} runs ({100*len(mrna_fp)/max(1,len(mrna_fp_pool)):.0f}%)")
    if len(mrna_fp) > 0:
        for ss in sorted(mrna_fp.strand_specificity.unique()):
            sub = mrna_fp[mrna_fp.strand_specificity == ss]
            print(f"     SS={ss:.2f}: {len(sub)} FPs, mean mRNA_obs={sub.total_mrna_observed.mean():.1f}")

    # ── Section 6: nRNA-gDNA confusion ──
    print(f"\n{SEP}")
    print("6. nRNA-gDNA CONFUSION ANALYSIS")
    print(SEP)
    both = df[(df.NTA1 > 0) & (df.gdna_fraction > 0) & (df.TA1 > 0)].copy()
    print(f"\n  Runs with mRNA + nRNA + gDNA all present: {len(both)}")
    for ss in sorted(both.strand_specificity.unique()):
        sub = both[both.strand_specificity == ss]
        corr_gn = sub[["gdna_excess", "nrna_excess"]].corr().iloc[0, 1]
        print(f"\n    SS={ss:.2f} ({len(sub)} runs):")
        print(f"      Mean gDNA excess: {sub.gdna_excess.mean():+.1f}")
        print(f"      Mean nRNA excess: {sub.nrna_excess.mean():+.1f}")
        print(f"      Mean mRNA excess: {sub.mrna_excess.mean():+.1f}")
        print(f"      Conservation:     {(sub.gdna_excess + sub.nrna_excess + sub.mrna_excess).mean():+.2f}")
        print(f"      Corr(gDNA,nRNA):  {corr_gn:.4f}")

    # ── Section 7: Worst cases ──
    print(f"\n{SEP}")
    print("7. WORST CASES")
    print(SEP)

    # Top gDNA overestimation
    has_gdna = df[df.gdna_expected > 0].copy()
    has_gdna_sorted = has_gdna.sort_values("gdna_rel_signed", ascending=False)
    print("\n  Top 5 gDNA overestimation:")
    for _, r in has_gdna_sorted.head(5).iterrows():
        print(f"    SS={r.strand_specificity:.2f} gf={r.gdna_fraction:.1f} "
              f"TA1={int(r.TA1)} NTA1={int(r.NTA1)} → "
              f"gDNA: exp={r.gdna_expected:.0f} obs={r.gdna_observed:.0f} "
              f"err={r.gdna_rel_signed:+.2%}")

    # Top nRNA underestimation
    has_nrna = df[df.nrna_expected > 0].copy()
    has_nrna_sorted = has_nrna.sort_values("nrna_rel_signed")
    print("\n  Top 5 nRNA underestimation:")
    for _, r in has_nrna_sorted.head(5).iterrows():
        print(f"    SS={r.strand_specificity:.2f} gf={r.gdna_fraction:.1f} "
              f"TA1={int(r.TA1)} NTA1={int(r.NTA1)} → "
              f"nRNA: exp={r.nrna_expected:.0f} obs={r.nrna_observed:.0f} "
              f"err={r.nrna_rel_signed:+.2%}")

    # Top nRNA FP (spurious nRNA when NTA1=0)
    nrna_fp_sorted = nrna_fp.sort_values("nrna_observed", ascending=False)
    print("\n  Top 5 spurious nRNA (NTA1=0):")
    for _, r in nrna_fp_sorted.head(5).iterrows():
        print(f"    SS={r.strand_specificity:.2f} gf={r.gdna_fraction:.1f} "
              f"TA1={int(r.TA1)} → nRNA_obs={r.nrna_observed:.0f}")

    # ── Section 8: Fragment length recovery ──
    print(f"\n{SEP}")
    print("8. FRAGMENT LENGTH RECOVERY")
    print(SEP)
    has_rna_fl = df[df.rna_fl_n_obs > 0]
    if len(has_rna_fl) > 0:
        rna_err = has_rna_fl.rna_fl_mean_err
        print(f"  RNA FL mean error: {rna_err.mean():+.2f} bp (std={rna_err.std():.2f})")
    has_gdna_fl = df[df.gdna_fl_n_obs > 0]
    if len(has_gdna_fl) > 0:
        gfl = has_gdna_fl.gdna_fl_mean_err
        print(f"  gDNA FL mean error: {gfl.mean():+.2f} bp (std={gfl.std():.2f})")

    # ── Section 9: Weighted accuracy ──
    print(f"\n{SEP}")
    print("9. WEIGHTED ACCURACY (by true count)")
    print(SEP)
    valid = df[(df.gdna_fraction > 0) & (df.NTA1 > 0) & (df.TA1 > 0)].copy()
    for ss in sorted(valid.strand_specificity.unique()):
        sub = valid[valid.strand_specificity == ss]
        mrna_w = (sub.total_mrna_observed.sum() - sub.total_mrna_expected.sum()) / sub.total_mrna_expected.sum()
        nrna_w = (sub.nrna_observed.sum() - sub.nrna_expected.sum()) / sub.nrna_expected.sum()
        gdna_w = (sub.gdna_observed.sum() - sub.gdna_expected.sum()) / sub.gdna_expected.sum()
        print(f"  SS={ss:.2f}: mRNA={mrna_w:+.4f}  nRNA={nrna_w:+.4f}  gDNA={gdna_w:+.4f}")

    # ── Section 10: Edge cases ──
    print(f"\n{SEP}")
    print("10. EDGE CASES")
    print(SEP)

    # A: NTA1 >> TA1 (high nRNA, low mRNA)
    print("\n  A. High nRNA / low mRNA (NTA1 >= 200, 0 < TA1 <= 20, gdna=0.5):")
    edge_a = df[(df.NTA1 >= 200) & (df.TA1 > 0) & (df.TA1 <= 20) & (df.gdna_fraction == 0.5)]
    for ss in sorted(edge_a.strand_specificity.unique()):
        sub = edge_a[edge_a.strand_specificity == ss]
        me = sub.mrna_rel_signed.mean()
        ne = sub.nrna_rel_signed.mean()
        print(f"    SS={ss:.2f}: mRNA err={me:+.4f}  nRNA err={ne:+.4f} ({len(sub)} runs)")

    # B: TA1 >> NTA1 (high mRNA, low nRNA)
    print("\n  B. High mRNA / low nRNA (TA1 >= 200, 0 < NTA1 <= 20, gdna=0.5):")
    edge_b = df[(df.TA1 >= 200) & (df.NTA1 > 0) & (df.NTA1 <= 20) & (df.gdna_fraction == 0.5)]
    for ss in sorted(edge_b.strand_specificity.unique()):
        sub = edge_b[edge_b.strand_specificity == ss]
        ne = sub.nrna_rel_signed.mean()
        ge = sub.gdna_rel_signed.mean()
        print(f"    SS={ss:.2f}: nRNA err={ne:+.4f}  gDNA err={ge:+.4f} ({len(sub)} runs)")

    # C: gDNA-only (TA1=0, NTA1=0, gdna_frac > 0)
    print("\n  C. gDNA-only runs (TA1=0, NTA1=0, gdna_frac > 0):")
    gdna_only = df[(df.TA1 == 0) & (df.NTA1 == 0) & (df.gdna_fraction > 0)]
    if len(gdna_only) > 0:
        for ss in sorted(gdna_only.strand_specificity.unique()):
            sub = gdna_only[gdna_only.strand_specificity == ss]
            ge = sub.gdna_rel_signed.dropna()
            mobs = sub.total_mrna_observed.mean()
            nobs = sub.nrna_observed.mean()
            gobs = sub.gdna_observed.mean()
            gexp = sub.gdna_expected.mean()
            print(f"    SS={ss:.2f}: gDNA err={ge.mean():+.4f}  "
                  f"spurious mRNA={mobs:.1f}  spurious nRNA={nobs:.1f}  "
                  f"({len(sub)} runs)")
    else:
        print("    No gDNA-only runs found")

    print(f"\n{SEP}")
    print("ANALYSIS COMPLETE")
    print(SEP)


if __name__ == "__main__":
    main()
