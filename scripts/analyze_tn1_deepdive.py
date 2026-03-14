#!/usr/bin/env python
"""Deep-dive analysis: false positives, SS=0.5 dynamics, identifiability."""
import numpy as np
import pandas as pd
from pathlib import Path

TSV = Path("sweep_results/tn1_ablation/sweep_results.tsv")

def main():
    df = pd.read_csv(TSV, sep="\t")

    # Signed errors
    df["gdna_rel_signed"] = np.where(
        df.gdna_expected > 0,
        (df.gdna_observed - df.gdna_expected) / df.gdna_expected, 0.0)
    df["nrna_rel_signed"] = np.where(
        df.nrna_expected > 0,
        (df.nrna_observed - df.nrna_expected) / df.nrna_expected, 0.0)
    df["mrna_rel_signed"] = np.where(
        df.total_mrna_expected > 0,
        (df.total_mrna_observed - df.total_mrna_expected) / df.total_mrna_expected, 0.0)

    # ── A: FALSE POSITIVE deep-dive ──
    print("=" * 74)
    print("A. FALSE POSITIVE DEEP-DIVE")
    print("=" * 74)

    # 1. gDNA false positives: gdna_fraction=0, nRNA present
    print("\n--- gDNA false positives when gdna_frac=0 ---")
    no_gdna = df[df.gdna_fraction == 0].copy()
    for ss in sorted(no_gdna.strand_specificity.unique()):
        sub = no_gdna[no_gdna.strand_specificity == ss]
        fp = sub[sub.gdna_observed > 1.0]
        print(f"\n  SS={ss:.2f}: {len(fp)}/{len(sub)} runs have spurious gDNA > 1")
        if len(fp) > 0:
            print(f"    mean spurious gDNA obs: {fp.gdna_observed.mean():.1f}")
            print(f"    These runs have NTA1: {sorted(fp.NTA1.unique().tolist())}")
            # Show pattern
            for _, r in fp.head(5).iterrows():
                print(f"    TA1={r.TA1:.0f} NTA1={r.NTA1:.0f} -> gDNA_obs={r.gdna_observed:.1f}  nRNA_obs={r.nrna_observed:.1f}  nRNA_exp={r.nrna_expected:.0f}")

    # 2. nRNA false positives: NTA1=0, mRNA present
    print("\n\n--- nRNA false positives when NTA1=0 ---")
    no_nrna = df[df.NTA1 == 0].copy()
    for ss in sorted(no_nrna.strand_specificity.unique()):
        sub = no_nrna[no_nrna.strand_specificity == ss]
        fp = sub[sub.nrna_observed > 1.0]
        print(f"\n  SS={ss:.2f}: {len(fp)}/{len(sub)} runs have spurious nRNA > 1")
        if len(fp) > 0:
            print(f"    mean spurious nRNA obs: {fp.nrna_observed.mean():.1f}")
            # Show by gdna_fraction
            for gf in sorted(fp.gdna_fraction.unique()):
                s = fp[fp.gdna_fraction == gf]
                print(f"    gdna_frac={gf:.1f}: {len(s)} FPs, mean nRNA_obs={s.nrna_observed.mean():.1f}")

    # 3. mRNA false positives: TA1=0 but observed > 0
    print("\n\n--- mRNA false positives when TA1=0 ---")
    no_mrna = df[df.TA1 == 0].copy()
    for ss in sorted(no_mrna.strand_specificity.unique()):
        sub = no_mrna[no_mrna.strand_specificity == ss]
        fp = sub[sub.TA1_observed > 1.0]
        print(f"\n  SS={ss:.2f}: {len(fp)}/{len(sub)} runs have spurious mRNA > 1")
        if len(fp) > 0:
            print(f"    mean spurious mRNA obs: {fp.TA1_observed.mean():.1f}")
            top = fp.nlargest(3, "TA1_observed")
            for _, r in top.iterrows():
                print(f"    NTA1={r.NTA1:.0f} gdna_f={r.gdna_fraction:.1f} -> "
                      f"TA1_obs={r.TA1_observed:.1f} nRNA_obs={r.nrna_observed:.1f} nRNA_exp={r.nrna_expected:.0f}")

    # ── B: SS=0.5 identifiability collapse ──
    print("\n\n" + "=" * 74)
    print("B. SS=0.5 IDENTIFIABILITY COLLAPSE")
    print("=" * 74)

    ss05 = df[df.strand_specificity == 0.5].copy()

    # Where does nRNA go at SS=0.5?
    print("\nSS=0.5: Where does nRNA mass go?")
    ss05_nrna = ss05[ss05.nrna_expected > 0].copy()
    ss05_nrna["nrna_to_gdna"] = ss05_nrna.gdna_observed - ss05_nrna.gdna_expected
    ss05_nrna["nrna_to_mrna"] = ss05_nrna.total_mrna_observed - ss05_nrna.total_mrna_expected
    ss05_nrna["nrna_missing"] = ss05_nrna.nrna_expected - ss05_nrna.nrna_observed

    print(f"  mean nRNA missing: {ss05_nrna.nrna_missing.mean():.1f}")
    print(f"  mean excess gDNA:  {ss05_nrna.nrna_to_gdna.mean():.1f}")
    print(f"  mean excess mRNA:  {ss05_nrna.nrna_to_mrna.mean():.1f}")
    print(f"  Conservation check (missing - excess_gdna - excess_mrna): "
          f"{(ss05_nrna.nrna_missing - ss05_nrna.nrna_to_gdna - ss05_nrna.nrna_to_mrna).mean():.1f}")

    # SS=0.5, no gDNA — pure nRNA->mRNA confusion
    ss05_no_g = ss05[(ss05.gdna_fraction == 0) & (ss05.nrna_expected > 0)]
    if len(ss05_no_g) > 0:
        print(f"\n  SS=0.5, no gDNA, nRNA present ({len(ss05_no_g)} runs):")
        print(f"    nRNA signed err: mean={ss05_no_g.nrna_rel_signed.mean():+.4f}")
        print(f"    mRNA signed err: mean={ss05_no_g.mrna_rel_signed.mean():+.4f}")
        print(f"    -> At SS=0.5 nRNA is indistinguishable from mRNA (no strand signal)")

    # ── C: NTA1=0 false positive root cause ──
    print("\n\n" + "=" * 74)
    print("C. NTA1=0 FALSE POSITIVE ROOT CAUSE: nRNA EMERGING FROM gDNA")
    print("=" * 74)

    # Filter: NTA1=0, gDNA present, nRNA false positive
    fp_data = df[(df.NTA1 == 0) & (df.gdna_fraction > 0)].copy()
    fp_data["nrna_excess"] = fp_data.nrna_observed  # all nRNA is spurious
    fp_data["gdna_deficit"] = fp_data.gdna_expected - fp_data.gdna_observed

    for ss in sorted(fp_data.strand_specificity.unique()):
        sub = fp_data[fp_data.strand_specificity == ss]
        sub_fp = sub[sub.nrna_observed > 1.0]
        if len(sub_fp) > 0:
            print(f"\n  SS={ss:.2f}: {len(sub_fp)} FP runs")
            print(f"    mean spurious nRNA: {sub_fp.nrna_excess.mean():.1f}")
            print(f"    mean gDNA deficit:  {sub_fp.gdna_deficit.mean():.1f}")
            print(f"    ratio (nRNA_excess / gDNA_deficit): {(sub_fp.nrna_excess / sub_fp.gdna_deficit.replace(0, np.nan)).mean():.3f}")

    # ── D: gDNA underestimation at SS >= 0.75 ──
    print("\n\n" + "=" * 74)
    print("D. gDNA UNDERESTIMATION AT SS >= 0.75")
    print("=" * 74)

    good_ss = df[(df.strand_specificity >= 0.75) & (df.gdna_expected > 0)].copy()
    good_ss["gdna_deficit_pct"] = (good_ss.gdna_expected - good_ss.gdna_observed) / good_ss.gdna_expected * 100

    for ss in sorted(good_ss.strand_specificity.unique()):
        sub = good_ss[good_ss.strand_specificity == ss]
        under = sub[sub.gdna_rel_signed < -0.1]
        print(f"\n  SS={ss:.2f}: {len(under)}/{len(sub)} runs have gDNA underestimation > 10%")
        if len(under) > 0:
            print(f"    mean deficit: {under.gdna_deficit_pct.mean():.1f}%")
            print(f"    median deficit: {under.gdna_deficit_pct.median():.1f}%")
            # Where does the missing gDNA go?
            under_nrna = under[under.nrna_expected > 0]
            if len(under_nrna) > 0:
                excess_nrna = (under_nrna.nrna_observed - under_nrna.nrna_expected).mean()
                deficit_gdna = (under_nrna.gdna_expected - under_nrna.gdna_observed).mean()
                print(f"    when nRNA present: mean gDNA deficit={deficit_gdna:.1f}, mean nRNA excess={excess_nrna:.1f}")

    # ── E: Comparison SS=1.0 vs SS=0.9: the best case ──
    print("\n\n" + "=" * 74)
    print("E. BEST-CASE PERFORMANCE: SS=0.9 and SS=1.0")
    print("=" * 74)

    for ss in [0.9, 1.0]:
        sub = df[(df.strand_specificity == ss)]
        sub_g = sub[(sub.gdna_expected > 0)]
        sub_n = sub[(sub.nrna_expected > 0)]
        sub_m = sub[(sub.total_mrna_expected > 0)]
        sub_both = sub[(sub.gdna_expected > 0) & (sub.nrna_expected > 0)]

        print(f"\nSS={ss:.1f}:")
        if len(sub_g) > 0:
            ge = sub_g.gdna_rel_signed
            print(f"  gDNA: mean={ge.mean():+.4f}  |mean|={ge.abs().mean():.4f}  "
                  f"good(<10%)={((ge.abs() < 0.1).mean()*100):.0f}%  max={ge.abs().max():.4f}")
        if len(sub_n) > 0:
            ne = sub_n.nrna_rel_signed
            print(f"  nRNA: mean={ne.mean():+.4f}  |mean|={ne.abs().mean():.4f}  "
                  f"good(<10%)={((ne.abs() < 0.1).mean()*100):.0f}%  max={ne.abs().max():.4f}")
        if len(sub_m) > 0:
            me = sub_m.mrna_rel_signed
            print(f"  mRNA: mean={me.mean():+.4f}  |mean|={me.abs().mean():.4f}  "
                  f"good(<10%)={((me.abs() < 0.1).mean()*100):.0f}%  max={me.abs().max():.4f}")

    # ── F: Balanced accuracy summary table ──
    print("\n\n" + "=" * 74)
    print("F. ACCURACY SUMMARY BY SS (gdna_fraction=0.3, TA1>0, NTA1>0)")
    print("=" * 74)

    sub = df[(df.gdna_fraction == 0.3) & (df.TA1 > 0) & (df.NTA1 > 0)]
    for ss in sorted(sub.strand_specificity.unique()):
        s = sub[sub.strand_specificity == ss]
        me = (s.total_mrna_observed - s.total_mrna_expected) / s.total_mrna_expected
        ne = (s.nrna_observed - s.nrna_expected) / s.nrna_expected
        ge = (s.gdna_observed - s.gdna_expected) / s.gdna_expected
        print(f"  SS={ss:.2f}:  mRNA={me.mean():+.3f} ({me.abs().mean():.3f})  "
              f"nRNA={ne.mean():+.3f} ({ne.abs().mean():.3f})  "
              f"gDNA={ge.mean():+.3f} ({ge.abs().mean():.3f})")

    # ── G: Conservation of mass check ──
    print("\n\n" + "=" * 74)
    print("G. CONSERVATION OF MASS CHECK")
    print("=" * 74)

    df["total_obs"] = df.total_mrna_observed + df.nrna_observed + df.gdna_observed
    df["total_exp"] = df.total_mrna_expected + df.nrna_expected + df.gdna_expected
    df["mass_diff"] = df.total_obs - df.total_exp
    df["mass_rel"] = df.mass_diff / df.total_exp

    print(f"\n  Total observed vs expected:")
    print(f"    mean mass_rel_err: {df.mass_rel.mean():+.4f}")
    print(f"    median mass_rel_err: {df.mass_rel.median():+.4f}")
    print(f"    p5/p95: {df.mass_rel.quantile(0.05):+.4f} / {df.mass_rel.quantile(0.95):+.4f}")

    # By SS
    for ss in sorted(df.strand_specificity.unique()):
        sub = df[df.strand_specificity == ss]
        print(f"    SS={ss:.2f}: mass_rel_err mean={sub.mass_rel.mean():+.4f}")

    # ── H: Fragment length model impact ──
    print("\n\n" + "=" * 74)
    print("H. FRAGMENT LENGTH MODEL RECOVERY")
    print("=" * 74)

    valid_fl = df[df.rna_fl_est_mean != ""].copy()
    if len(valid_fl) > 0:
        valid_fl["rna_fl_est_mean"] = pd.to_numeric(valid_fl["rna_fl_est_mean"])
        valid_fl["rna_fl_mean_err"] = pd.to_numeric(valid_fl["rna_fl_mean_err"])
        print(f"\n  RNA FL model (true mean=250):")
        print(f"    mean est: {valid_fl.rna_fl_est_mean.mean():.2f}")
        print(f"    mean err: {valid_fl.rna_fl_mean_err.mean():.2f}")

    valid_gdna_fl = df[(df.gdna_fl_est_mean != "") & (df.gdna_fraction > 0)].copy()
    if len(valid_gdna_fl) > 0:
        valid_gdna_fl["gdna_fl_est_mean"] = pd.to_numeric(valid_gdna_fl["gdna_fl_est_mean"])
        valid_gdna_fl["gdna_fl_mean_err"] = pd.to_numeric(valid_gdna_fl["gdna_fl_mean_err"])
        print(f"\n  gDNA FL model (true mean=250):")
        print(f"    mean est: {valid_gdna_fl.gdna_fl_est_mean.mean():.2f}")
        print(f"    mean err: {valid_gdna_fl.gdna_fl_mean_err.mean():.2f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
