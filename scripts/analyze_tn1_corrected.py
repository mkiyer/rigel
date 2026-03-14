#!/usr/bin/env python
"""Corrected false positive analysis and nRNA/gDNA confusion patterns."""
import numpy as np
import pandas as pd
from pathlib import Path

TSV = Path("sweep_results/tn1_ablation/sweep_results.tsv")

def main():
    df = pd.read_csv(TSV, sep="\t")

    # ── Degenerate case identification ──
    # When TA1=0, NTA1=0, TD=0, TE=0: the simulator distributes mRNA uniformly
    # These are NOT real "false positives" — they're simulator artifacts
    df["all_zero_rna"] = (df.TA1 == 0) & (df.NTA1 == 0)

    print("=" * 74)
    print("1. SIMULATOR DEGENERATE CASE CHECK")
    print("=" * 74)
    degen = df[df.all_zero_rna]
    print("When TA1=0 AND NTA1=0:")
    print("  The simulator in n_rna_fragments mode still generates 10K fragments")
    print("  distributed uniformly across transcripts. These are NOT zero-expression")
    print("  scenarios — they're 'equal expression' scenarios.")
    print("  => Excluding these from false positive analysis.")
    print("  ({} runs affected)".format(len(degen)))

    # ── REAL false positives ──
    print("\n" + "=" * 74)
    print("2. REAL FALSE POSITIVES")
    print("=" * 74)

    # A. nRNA false positives: NTA1=0, TA1>0 (real scenario)
    real_nrna_fp_pool = df[(df.NTA1 == 0) & (df.TA1 > 0)]
    fp = real_nrna_fp_pool[real_nrna_fp_pool.nrna_observed > 1.0]
    print("\nA. nRNA false positives (NTA1=0, TA1>0, nRNA_obs > 1):")
    print("   {} / {} runs".format(len(fp), len(real_nrna_fp_pool)))
    if len(fp) > 0:
        for ss in sorted(fp.strand_specificity.unique()):
            sub = fp[fp.strand_specificity == ss]
            print("   SS={:.2f}: {} FPs, mean={:.1f}, max={:.1f}".format(
                ss, len(sub), sub.nrna_observed.mean(), sub.nrna_observed.max()))
            for gf in sorted(sub.gdna_fraction.unique()):
                s = sub[sub.gdna_fraction == gf]
                print("     gdna_frac={:.1f}: {} runs, mean nRNA_obs={:.1f}".format(
                    gf, len(s), s.nrna_observed.mean()))

    # B. gDNA false positives: gdna_frac=0, any_rna=True
    real_gdna_fp_pool = df[(df.gdna_fraction == 0) & ~df.all_zero_rna]
    fp_g = real_gdna_fp_pool[real_gdna_fp_pool.gdna_observed > 1.0]
    print("\nB. gDNA false positives (gdna_frac=0, has RNA, gDNA_obs > 1):")
    print("   {} / {} runs".format(len(fp_g), len(real_gdna_fp_pool)))
    if len(fp_g) > 0:
        for ss in sorted(fp_g.strand_specificity.unique()):
            sub = fp_g[fp_g.strand_specificity == ss]
            print("   SS={:.2f}: {} FPs, mean gDNA_obs={:.1f}".format(
                ss, len(sub), sub.gdna_observed.mean()))

    # ── nRNA-gDNA confusion matrix ──
    print("\n" + "=" * 74)
    print("3. nRNA-gDNA CONFUSION: WHERE DOES MISASSIGNED MASS GO?")
    print("=" * 74)

    # Only consider runs where both are present and non-zero
    both = df[(df.NTA1 > 0) & (df.gdna_fraction > 0) & (df.TA1 > 0)].copy()
    both["gdna_over"] = both.gdna_observed - both.gdna_expected
    both["nrna_over"] = both.nrna_observed - both.nrna_expected
    both["mrna_over"] = both.total_mrna_observed - both.total_mrna_expected

    print("\nWith all three (mRNA, nRNA, gDNA) present ({} runs):".format(len(both)))

    for ss in sorted(both.strand_specificity.unique()):
        sub = both[both.strand_specificity == ss]
        print("\n  SS={:.2f}:".format(ss))
        print("    Mean gDNA excess: {:+.1f}".format(sub.gdna_over.mean()))
        print("    Mean nRNA excess: {:+.1f}".format(sub.nrna_over.mean()))
        print("    Mean mRNA excess: {:+.1f}".format(sub.mrna_over.mean()))
        print("    Sum (should = 0): {:+.1f}".format(
            (sub.gdna_over + sub.nrna_over + sub.mrna_over).mean()))

        # Correlation patterns
        corr_gn = sub[["gdna_over", "nrna_over"]].corr().iloc[0, 1]
        corr_gm = sub[["gdna_over", "mrna_over"]].corr().iloc[0, 1]
        corr_nm = sub[["nrna_over", "mrna_over"]].corr().iloc[0, 1]
        print("    Corr(gDNA, nRNA): {:.4f}".format(corr_gn))
        print("    Corr(gDNA, mRNA): {:.4f}".format(corr_gm))
        print("    Corr(nRNA, mRNA): {:.4f}".format(corr_nm))

    # ── Accuracy excluding problematic cases ──
    print("\n" + "=" * 74)
    print("4. ACCURACY EXCLUDING DEGENERATE CASES")
    print("=" * 74)

    clean = df[~df.all_zero_rna & (df.TA1 > 0)]
    clean["gdna_rel"] = np.where(
        clean.gdna_expected > 0,
        (clean.gdna_observed - clean.gdna_expected) / clean.gdna_expected, np.nan)
    clean["nrna_rel"] = np.where(
        clean.nrna_expected > 0,
        (clean.nrna_observed - clean.nrna_expected) / clean.nrna_expected, np.nan)
    clean["mrna_rel"] = np.where(
        clean.total_mrna_expected > 0,
        (clean.total_mrna_observed - clean.total_mrna_expected) / clean.total_mrna_expected, np.nan)

    for ss in sorted(clean.strand_specificity.unique()):
        sub = clean[clean.strand_specificity == ss]
        me = sub.mrna_rel.dropna()
        ne = sub.nrna_rel.dropna()
        ge = sub.gdna_rel.dropna()
        print("\nSS={:.2f} ({} runs):".format(ss, len(sub)))
        if len(me) > 0:
            print("  mRNA: mean={:+.4f}  |mean|={:.4f}  p10={:+.4f}  p90={:+.4f}".format(
                me.mean(), me.abs().mean(), me.quantile(0.1), me.quantile(0.9)))
        if len(ne) > 0:
            print("  nRNA: mean={:+.4f}  |mean|={:.4f}  p10={:+.4f}  p90={:+.4f}".format(
                ne.mean(), ne.abs().mean(), ne.quantile(0.1), ne.quantile(0.9)))
        if len(ge) > 0:
            print("  gDNA: mean={:+.4f}  |mean|={:.4f}  p10={:+.4f}  p90={:+.4f}".format(
                ge.mean(), ge.abs().mean(), ge.quantile(0.1), ge.quantile(0.9)))

    # ── Key pattern: high NTA1, low TA1 ──
    print("\n" + "=" * 74)
    print("5. HIGH nRNA / LOW mRNA PATTERN (NTA1 >> TA1)")
    print("=" * 74)

    high_nrna = df[(df.NTA1 >= 200) & (df.TA1 <= 20) & (df.TA1 > 0) &
                   (df.gdna_fraction == 0.5)].copy()
    if len(high_nrna) > 0:
        for ss in sorted(high_nrna.strand_specificity.unique()):
            sub = high_nrna[high_nrna.strand_specificity == ss]
            mrna_err = (sub.total_mrna_observed - sub.total_mrna_expected) / sub.total_mrna_expected
            nrna_err = np.where(sub.nrna_expected > 0,
                               (sub.nrna_observed - sub.nrna_expected) / sub.nrna_expected, 0)
            print("  SS={:.2f}: mRNA_err={:+.4f}  nRNA_err={:+.4f} ({} runs)".format(
                ss, mrna_err.mean(), nrna_err.mean(), len(sub)))

    # ── Key pattern: high TA1, low NTA1 ──
    print("\n" + "=" * 74)
    print("6. HIGH mRNA / LOW nRNA PATTERN (TA1 >> NTA1)")
    print("=" * 74)

    high_mrna = df[(df.TA1 >= 200) & (df.NTA1 <= 20) & (df.NTA1 > 0) &
                   (df.gdna_fraction == 0.5)].copy()
    if len(high_mrna) > 0:
        for ss in sorted(high_mrna.strand_specificity.unique()):
            sub = high_mrna[high_mrna.strand_specificity == ss]
            nrna_err = (sub.nrna_observed - sub.nrna_expected) / sub.nrna_expected
            gdna_err = (sub.gdna_observed - sub.gdna_expected) / sub.gdna_expected
            print("  SS={:.2f}: nRNA_err={:+.4f}  gDNA_err={:+.4f} ({} runs)".format(
                ss, nrna_err.mean(), gdna_err.mean(), len(sub)))

    # ── Revenue-weighted accuracy (weighted by true abundance) ──
    print("\n" + "=" * 74)
    print("7. WEIGHTED ACCURACY (by true count, excluding degenerate)")
    print("=" * 74)

    valid = df[~df.all_zero_rna & (df.gdna_fraction > 0) & (df.NTA1 > 0) & (df.TA1 > 0)].copy()
    for ss in sorted(valid.strand_specificity.unique()):
        sub = valid[valid.strand_specificity == ss]
        # Weighted errors
        total_mrna_exp = sub.total_mrna_expected.sum()
        total_mrna_obs = sub.total_mrna_observed.sum()
        total_nrna_exp = sub.nrna_expected.sum()
        total_nrna_obs = sub.nrna_observed.sum()
        total_gdna_exp = sub.gdna_expected.sum()
        total_gdna_obs = sub.gdna_observed.sum()

        mrna_w = (total_mrna_obs - total_mrna_exp) / total_mrna_exp
        nrna_w = (total_nrna_obs - total_nrna_exp) / total_nrna_exp
        gdna_w = (total_gdna_obs - total_gdna_exp) / total_gdna_exp

        print("  SS={:.2f}: mRNA={:+.4f}  nRNA={:+.4f}  gDNA={:+.4f}".format(
            ss, mrna_w, nrna_w, gdna_w))


if __name__ == "__main__":
    main()
