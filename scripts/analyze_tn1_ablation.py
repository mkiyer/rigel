#!/usr/bin/env python
"""Analyze T+N+1 ablation sweep results."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TSV = Path("sweep_results/tn1_ablation/sweep_results.tsv")

def main():
    df = pd.read_csv(TSV, sep="\t")
    print(f"Loaded {len(df)} runs, {df.shape[1]} columns")

    # Compute signed relative errors (direction matters)
    df["gdna_rel_err_signed"] = np.where(
        df.gdna_expected > 0,
        (df.gdna_observed - df.gdna_expected) / df.gdna_expected,
        0.0,
    )
    df["nrna_rel_err_signed"] = np.where(
        df.nrna_expected > 0,
        (df.nrna_observed - df.nrna_expected) / df.nrna_expected,
        0.0,
    )
    df["mrna_rel_err_signed"] = np.where(
        df.total_mrna_expected > 0,
        (df.total_mrna_observed - df.total_mrna_expected) / df.total_mrna_expected,
        0.0,
    )

    # ── Section 1: Overall statistics ──
    print("\n" + "=" * 74)
    print("SECTION 1: OVERALL ERROR STATISTICS")
    print("=" * 74)

    has_gdna = df[df.gdna_fraction > 0].copy()
    no_gdna = df[df.gdna_fraction == 0].copy()
    has_nrna = df[df.NTA1 > 0].copy()
    has_mrna = df[df.TA1 > 0].copy()

    for label, sub in [("All runs", df), ("No gDNA", no_gdna), ("With gDNA", has_gdna)]:
        print(f"\n--- {label} ({len(sub)} runs) ---")
        for col, name in [
            ("total_mrna_rel_err", "mRNA |rel_err|"),
            ("mrna_rel_err_signed", "mRNA signed_err"),
            ("nrna_rel_err", "|nRNA rel_err|"),
            ("nrna_rel_err_signed", "nRNA signed_err"),
        ]:
            vals = sub[col].replace([np.inf, -np.inf], np.nan).dropna()
            if len(vals) > 0:
                print(f"  {name:25s}: mean={vals.mean():+.4f}  med={vals.median():+.4f}  p5={vals.quantile(0.05):+.4f}  p95={vals.quantile(0.95):+.4f}")

    # gDNA-specific stats
    g = has_gdna[has_gdna.gdna_expected > 0]
    if len(g) > 0:
        print(f"\n--- gDNA estimation ({len(g)} runs with gDNA > 0) ---")
        gr = g.gdna_rel_err_signed
        print(f"  gDNA signed_err  : mean={gr.mean():+.4f}  med={gr.median():+.4f}  p5={gr.quantile(0.05):+.4f}  p95={gr.quantile(0.95):+.4f}")
        print(f"  Over  (>+10%): {(gr > 0.1).sum():4d} / {len(gr)} = {(gr > 0.1).mean():.1%}")
        print(f"  Under (<-10%): {(gr < -0.1).sum():4d} / {len(gr)} = {(gr < -0.1).mean():.1%}")
        print(f"  Good  (|e|<10%): {(gr.abs() < 0.1).sum():4d} / {len(gr)} = {(gr.abs() < 0.1).mean():.1%}")

    # ── Section 2: By strand specificity ──
    print("\n" + "=" * 74)
    print("SECTION 2: ERROR BY STRAND SPECIFICITY")
    print("=" * 74)

    for ss in sorted(df.strand_specificity.unique()):
        sub = df[df.strand_specificity == ss]
        sub_g = sub[(sub.gdna_fraction > 0) & (sub.gdna_expected > 0)]
        sub_n = sub[sub.nrna_expected > 0]
        sub_m = sub[sub.total_mrna_expected > 0]

        print(f"\nSS={ss:.2f} ({len(sub)} runs):")
        if len(sub_m) > 0:
            me = sub_m.mrna_rel_err_signed
            print(f"  mRNA  : mean={me.mean():+.4f}  med={me.median():+.4f}")
        if len(sub_n) > 0:
            ne = sub_n.nrna_rel_err_signed
            print(f"  nRNA  : mean={ne.mean():+.4f}  med={ne.median():+.4f}")
        if len(sub_g) > 0:
            ge = sub_g.gdna_rel_err_signed
            print(f"  gDNA  : mean={ge.mean():+.4f}  med={ge.median():+.4f}")

    # ── Section 3: By gDNA fraction ──
    print("\n" + "=" * 74)
    print("SECTION 3: ERROR BY gDNA FRACTION")
    print("=" * 74)

    for gf in sorted(df.gdna_fraction.unique()):
        sub = df[df.gdna_fraction == gf]
        sub_g = sub[sub.gdna_expected > 0]
        sub_n = sub[sub.nrna_expected > 0]
        sub_m = sub[sub.total_mrna_expected > 0]

        print(f"\ngdna_fraction={gf:.1f} ({len(sub)} runs):")
        if len(sub_m) > 0:
            me = sub_m.mrna_rel_err_signed
            print(f"  mRNA  : mean={me.mean():+.4f}  med={me.median():+.4f}")
        if len(sub_n) > 0:
            ne = sub_n.nrna_rel_err_signed
            print(f"  nRNA  : mean={ne.mean():+.4f}  med={ne.median():+.4f}")
        if len(sub_g) > 0:
            ge = sub_g.gdna_rel_err_signed
            print(f"  gDNA  : mean={ge.mean():+.4f}  med={ge.median():+.4f}")

    # ── Section 4: Cross-table: SS × gDNA fraction ──
    print("\n" + "=" * 74)
    print("SECTION 4: gDNA SIGNED ERROR — SS × gDNA_FRACTION CROSS TABLE")
    print("=" * 74)

    g_only = df[(df.gdna_expected > 0)].copy()
    if len(g_only) > 0:
        pt = g_only.pivot_table(
            values="gdna_rel_err_signed",
            index="strand_specificity",
            columns="gdna_fraction",
            aggfunc="mean",
        )
        print("\nMean gDNA signed relative error:")
        print(pt.to_string(float_format=lambda x: f"{x:+.3f}"))

    print("\n")
    n_only = df[(df.nrna_expected > 0)].copy()
    if len(n_only) > 0:
        pt = n_only.pivot_table(
            values="nrna_rel_err_signed",
            index="strand_specificity",
            columns="gdna_fraction",
            aggfunc="mean",
        )
        print("Mean nRNA signed relative error:")
        print(pt.to_string(float_format=lambda x: f"{x:+.3f}"))

    print("\n")
    m_only = df[(df.total_mrna_expected > 0)].copy()
    if len(m_only) > 0:
        pt = m_only.pivot_table(
            values="mrna_rel_err_signed",
            index="strand_specificity",
            columns="gdna_fraction",
            aggfunc="mean",
        )
        print("Mean mRNA signed relative error:")
        print(pt.to_string(float_format=lambda x: f"{x:+.3f}"))

    # ── Section 5: Worst cases ──
    print("\n" + "=" * 74)
    print("SECTION 5: WORST CASES")
    print("=" * 74)

    # Worst gDNA overestimation
    if len(g_only) > 0:
        worst_g = g_only.nlargest(10, "gdna_rel_err_signed")
        print("\nTop 10 gDNA OVERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'gDNA_exp':>9} {'gDNA_obs':>9} {'gDNA_err':>9}")
        for _, r in worst_g.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.gdna_expected:9.0f} {r.gdna_observed:9.1f} {r.gdna_rel_err_signed:+9.3f}")

        worst_g_under = g_only.nsmallest(10, "gdna_rel_err_signed")
        print("\nTop 10 gDNA UNDERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'gDNA_exp':>9} {'gDNA_obs':>9} {'gDNA_err':>9}")
        for _, r in worst_g_under.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.gdna_expected:9.0f} {r.gdna_observed:9.1f} {r.gdna_rel_err_signed:+9.3f}")

    # Worst nRNA (where nRNA expected > 0)
    if len(n_only) > 0:
        worst_n = n_only.nlargest(10, "nrna_rel_err_signed")
        print("\nTop 10 nRNA OVERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'nRNA_exp':>9} {'nRNA_obs':>9} {'nRNA_err':>9}")
        for _, r in worst_n.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.nrna_expected:9.0f} {r.nrna_observed:9.1f} {r.nrna_rel_err_signed:+9.3f}")

        worst_n_under = n_only.nsmallest(10, "nrna_rel_err_signed")
        print("\nTop 10 nRNA UNDERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'nRNA_exp':>9} {'nRNA_obs':>9} {'nRNA_err':>9}")
        for _, r in worst_n_under.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.nrna_expected:9.0f} {r.nrna_observed:9.1f} {r.nrna_rel_err_signed:+9.3f}")

    # Worst mRNA
    if len(m_only) > 0:
        worst_m = m_only.nlargest(10, "mrna_rel_err_signed")
        print("\nTop 10 mRNA OVERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'mRNA_exp':>9} {'mRNA_obs':>9} {'mRNA_err':>9}")
        for _, r in worst_m.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.total_mrna_expected:9.0f} {r.total_mrna_observed:9.1f} {r.mrna_rel_err_signed:+9.3f}")

        worst_m_under = m_only.nsmallest(10, "mrna_rel_err_signed")
        print("\nTop 10 mRNA UNDERESTIMATION:")
        print(f"{'SS':>5} {'gdna_f':>7} {'TA1':>5} {'NTA1':>5} {'mRNA_exp':>9} {'mRNA_obs':>9} {'mRNA_err':>9}")
        for _, r in worst_m_under.iterrows():
            print(f"{r.strand_specificity:5.2f} {r.gdna_fraction:7.1f} {r.TA1:5.0f} {r.NTA1:5.0f} "
                  f"{r.total_mrna_expected:9.0f} {r.total_mrna_observed:9.1f} {r.mrna_rel_err_signed:+9.3f}")

    # ── Section 6: nRNA vs gDNA confusion analysis ──
    print("\n" + "=" * 74)
    print("SECTION 6: nRNA vs gDNA CONFUSION (where both present)")
    print("=" * 74)

    both = df[(df.nrna_expected > 0) & (df.gdna_expected > 0)].copy()
    if len(both) > 0:
        # Check if gDNA overestimation correlates with nRNA underestimation
        both["gdna_over"] = both.gdna_observed - both.gdna_expected
        both["nrna_over"] = both.nrna_observed - both.nrna_expected
        corr = both[["gdna_over", "nrna_over"]].corr().iloc[0, 1]
        print(f"\nCorrelation(gDNA_overcount, nRNA_overcount) = {corr:.4f}")
        print("  (Negative = gDNA siphoning from nRNA, the old T+N+2 problem)")

        # By SS
        for ss in sorted(both.strand_specificity.unique()):
            sub = both[both.strand_specificity == ss]
            c = sub[["gdna_over", "nrna_over"]].corr().iloc[0, 1]
            print(f"  SS={ss:.2f}: corr={c:+.4f} ({len(sub)} runs)")

    # ── Section 7: Zero-expression scenarios (false positives) ──
    print("\n" + "=" * 74)
    print("SECTION 7: FALSE POSITIVES (expected=0, observed>0)")
    print("=" * 74)

    # nRNA false positives (NTA1=0 but nRNA observed > 0)
    nrna_fp = df[(df.NTA1 == 0) & (df.nrna_observed > 1.0)]
    print(f"\nnRNA false positives (NTA1=0, observed>1): {len(nrna_fp)} / {len(df[df.NTA1 == 0])}")
    if len(nrna_fp) > 0:
        print(f"  median spurious nRNA: {nrna_fp.nrna_observed.median():.1f}")
        print(f"  max spurious nRNA: {nrna_fp.nrna_observed.max():.1f}")
        top_fp = nrna_fp.nlargest(5, "nrna_observed")
        for _, r in top_fp.iterrows():
            print(f"    SS={r.strand_specificity:.2f} gdna_f={r.gdna_fraction:.1f} TA1={r.TA1:.0f} -> nRNA_obs={r.nrna_observed:.1f}")

    # gDNA false positives (gdna_fraction=0 but gDNA observed > 0)
    gdna_fp = df[(df.gdna_fraction == 0) & (df.gdna_observed > 1.0)]
    print(f"\ngDNA false positives (gdna_frac=0, observed>1): {len(gdna_fp)} / {len(df[df.gdna_fraction == 0])}")
    if len(gdna_fp) > 0:
        print(f"  median spurious gDNA: {gdna_fp.gdna_observed.median():.1f}")
        print(f"  max spurious gDNA: {gdna_fp.gdna_observed.max():.1f}")

    # mRNA false positives (TA1=0 but TA1 observed > 0)
    mrna_fp = df[(df.TA1 == 0) & (df.TA1_observed > 1.0)]
    print(f"\nmRNA false positives (TA1=0, observed>1): {len(mrna_fp)} / {len(df[df.TA1 == 0])}")
    if len(mrna_fp) > 0:
        print(f"  median spurious mRNA: {mrna_fp.TA1_observed.median():.1f}")
        print(f"  max spurious mRNA: {mrna_fp.TA1_observed.max():.1f}")

    # ── Section 8: Detailed NTA1 vs TA1 interaction ──
    print("\n" + "=" * 74)
    print("SECTION 8: nRNA ACCURACY BY (NTA1, TA1) — SS=1.0, gdna_frac=0.3")
    print("=" * 74)

    sub = df[(df.strand_specificity == 1.0) & (df.gdna_fraction == 0.3) & (df.NTA1 > 0)]
    if len(sub) > 0:
        pt = sub.pivot_table(
            values="nrna_rel_err_signed",
            index="NTA1",
            columns="TA1",
            aggfunc="mean",
        )
        print("\nnRNA signed relative error (rows=NTA1, cols=TA1):")
        print(pt.to_string(float_format=lambda x: f"{x:+.3f}"))

    sub2 = df[(df.strand_specificity == 1.0) & (df.gdna_fraction == 0.3) & (df.gdna_expected > 0)]
    if len(sub2) > 0:
        pt2 = sub2.pivot_table(
            values="gdna_rel_err_signed",
            index="NTA1",
            columns="TA1",
            aggfunc="mean",
        )
        print("\ngDNA signed relative error (rows=NTA1, cols=TA1):")
        print(pt2.to_string(float_format=lambda x: f"{x:+.3f}"))

    # ── Section 9: SS=0.5 deep dive (hardest case) ──
    print("\n" + "=" * 74)
    print("SECTION 9: SS=0.5 DEEP DIVE (non-strand-specific)")
    print("=" * 74)

    ss05 = df[df.strand_specificity == 0.5].copy()
    if len(ss05) > 0:
        # With gDNA
        ss05g = ss05[(ss05.gdna_expected > 0) & (ss05.nrna_expected > 0)]
        if len(ss05g) > 0:
            print(f"\nSS=0.5, both nRNA+gDNA present ({len(ss05g)} runs):")
            print(f"  gDNA signed err: mean={ss05g.gdna_rel_err_signed.mean():+.4f}, med={ss05g.gdna_rel_err_signed.median():+.4f}")
            print(f"  nRNA signed err: mean={ss05g.nrna_rel_err_signed.mean():+.4f}, med={ss05g.nrna_rel_err_signed.median():+.4f}")
            print(f"  mRNA signed err: mean={ss05g.mrna_rel_err_signed.mean():+.4f}, med={ss05g.mrna_rel_err_signed.median():+.4f}")

    # ── Section 10: High gDNA burden (gDNA_fraction >= 1.0) ──
    print("\n" + "=" * 74)
    print("SECTION 10: HIGH gDNA BURDEN (gdna_fraction >= 1.0)")
    print("=" * 74)

    high_g = df[df.gdna_fraction >= 1.0].copy()
    for gf in sorted(high_g.gdna_fraction.unique()):
        sub = high_g[high_g.gdna_fraction == gf]
        sub_g = sub[sub.gdna_expected > 0]
        sub_n = sub[sub.nrna_expected > 0]
        sub_m = sub[sub.total_mrna_expected > 0]
        print(f"\ngdna_fraction={gf:.1f}:")
        if len(sub_m) > 0:
            me = sub_m.mrna_rel_err_signed
            print(f"  mRNA: mean={me.mean():+.4f}  med={me.median():+.4f}  max_abs={me.abs().max():.4f}")
        if len(sub_n) > 0:
            ne = sub_n.nrna_rel_err_signed
            print(f"  nRNA: mean={ne.mean():+.4f}  med={ne.median():+.4f}  max_abs={ne.abs().max():.4f}")
        if len(sub_g) > 0:
            ge = sub_g.gdna_rel_err_signed
            print(f"  gDNA: mean={ge.mean():+.4f}  med={ge.median():+.4f}  max_abs={ge.abs().max():.4f}")

    # ── Section 11: Abundance-dependent bias ──
    print("\n" + "=" * 74)
    print("SECTION 11: ABUNDANCE-DEPENDENT BIAS (SS=1.0, gdna_frac=0.5)")
    print("=" * 74)

    sub = df[(df.strand_specificity == 1.0) & (df.gdna_fraction == 0.5)]
    if len(sub) > 0:
        for nta in sorted(sub.NTA1.unique()):
            s = sub[(sub.NTA1 == nta) & (sub.nrna_expected > 0)]
            if len(s) > 0:
                ne = s.nrna_rel_err_signed
                print(f"  NTA1={nta:4.0f}: nRNA_err={ne.mean():+.4f}  (n={len(s)})")

        print()
        for ta in sorted(sub.TA1.unique()):
            s = sub[(sub.TA1 == ta) & (sub.total_mrna_expected > 0)]
            if len(s) > 0:
                me = s.mrna_rel_err_signed
                print(f"  TA1={ta:4.0f}: mRNA_err={me.mean():+.4f}  (n={len(s)})")


if __name__ == "__main__":
    main()
