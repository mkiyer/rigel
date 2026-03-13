#!/usr/bin/env python3
"""Deep analysis of locus_config_simple baseline sweep results.

Grid: 4 SS × 6 gdna_frac × 9 TA1 × 9 NTA1 = 1944 runs
Scenario: 3 transcripts (TA1 multi-exon, TD single-exon, TE neg-ctrl)
          TD=0, TE=0, only TA1 and NTA1 vary
          RNA FL == gDNA FL (identical distributions)
          n_rna_fragments=10000, gdna as fraction on top
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

TSV = Path("/Users/mkiyer/Downloads/rigel_runs/locus_simple_ablations/baseline/sweep_results.tsv")


def load():
    df = pd.read_csv(TSV, sep="\t")
    # Derived columns
    df["gdna_rel_err"] = np.where(
        df["gdna_expected"] > 0,
        (df["gdna_observed"] - df["gdna_expected"]) / df["gdna_expected"],
        np.where(df["gdna_observed"] > 0, np.inf, 0.0),
    )
    df["gdna_abs_rel_err"] = np.abs(df["gdna_rel_err"])
    df["nrna_abs_rel_err"] = df["nrna_rel_err"].abs()
    df["mrna_abs_rel_err"] = df["total_mrna_rel_err"].abs()

    # TA1-specific signed error
    df["TA1_signed_err"] = np.where(
        df["TA1_expected"] > 0,
        (df["TA1_observed"] - df["TA1_expected"]) / df["TA1_expected"],
        np.where(df["TA1_observed"] > 0, np.inf, 0.0),
    )
    # nRNA signed error
    df["nrna_signed_err"] = np.where(
        df["nrna_expected"] > 0,
        (df["nrna_observed"] - df["nrna_expected"]) / df["nrna_expected"],
        np.where(df["nrna_observed"] > 0, np.inf, 0.0),
    )

    # Useful categories
    df["has_mrna"] = df["TA1"] > 0
    df["has_nrna"] = df["NTA1"] > 0
    df["has_gdna"] = df["gdna_fraction"] > 0
    df["scenario"] = df.apply(lambda r: (
        ("mRNA" if r["TA1"] > 0 else "") +
        ("+nRNA" if r["NTA1"] > 0 else "") +
        ("+gDNA" if r["gdna_fraction"] > 0 else "") or "empty"
    ), axis=1)

    return df


def sep(char="=", n=90):
    print(char * n)


def section(title):
    print(f"\n{'='*90}")
    print(f"  {title}")
    print(f"{'='*90}\n")


def main():
    df = load()
    n = len(df)

    # =====================================================================
    section("1. DATASET OVERVIEW")
    # =====================================================================
    print(f"Total runs: {n}")
    print(f"Grid: SS={sorted(df.strand_specificity.unique())} × "
          f"gdna_frac={sorted(df.gdna_fraction.unique())} × "
          f"TA1={sorted(df.TA1.unique())} × NTA1={sorted(df.NTA1.unique())}")
    print(f"TD={df.TD.unique()}, TE={df.TE.unique()} (always 0)")

    # =====================================================================
    section("2. ZERO-EXPRESSION BASELINE (TA1=0, NTA1=0)")
    # =====================================================================
    # When no RNA is expressed, ALL fragments are gDNA (or nothing)
    zero = df[(df["TA1"] == 0) & (df["NTA1"] == 0)]
    print(f"Runs with TA1=0, NTA1=0: {len(zero)}")

    # No-gDNA case
    zero_nogdna = zero[zero["gdna_fraction"] == 0]
    print(f"\n  gdna_fraction=0 ({len(zero_nogdna)} runs):")
    print(f"    n_fragments_actual: {zero_nogdna.n_fragments_actual.unique()}")
    print(f"    TA1_observed: {zero_nogdna.TA1_observed.unique()}")

    # With gDNA
    zero_gdna = zero[zero["gdna_fraction"] > 0]
    print(f"\n  gdna_fraction>0 ({len(zero_gdna)} runs):")
    for ss in sorted(zero_gdna.strand_specificity.unique()):
        sub = zero_gdna[zero_gdna.strand_specificity == ss]
        print(f"    SS={ss}:")
        for gf in sorted(sub.gdna_fraction.unique()):
            row = sub[sub.gdna_fraction == gf].iloc[0]
            gdna_exp = row["gdna_expected"]
            gdna_obs = row["gdna_observed"]
            nrna_obs = row["nrna_observed"]
            mrna_obs = row["total_mrna_observed"]
            err = (gdna_obs - gdna_exp) / gdna_exp * 100 if gdna_exp > 0 else 0
            print(f"      gf={gf}: gdna_exp={gdna_exp:.0f} gdna_obs={gdna_obs:.0f} "
                  f"err={err:+.1f}%  |  nrna_leak={nrna_obs:.0f}  mrna_leak={mrna_obs:.0f}")

    # =====================================================================
    section("3. PURE mRNA (NTA1=0, gdna_fraction=0)")
    # =====================================================================
    pure_mrna = df[(df["NTA1"] == 0) & (df["gdna_fraction"] == 0) & (df["TA1"] > 0)]
    print(f"Runs: {len(pure_mrna)}")
    for ss in sorted(pure_mrna.strand_specificity.unique()):
        sub = pure_mrna[pure_mrna.strand_specificity == ss]
        print(f"\n  SS={ss}: (TA1 values: {sorted(sub.TA1.unique())})")
        for _, row in sub.iterrows():
            ta1 = row["TA1"]
            obs = row["TA1_observed"]
            exp = row["TA1_expected"]
            nrna_obs = row["nrna_observed"]
            gdna_obs = row["gdna_observed"]
            err = (obs - exp) / exp * 100 if exp > 0 else 0
            print(f"    TA1={ta1:>5.0f}: expected={exp:>7.0f} observed={obs:>9.1f} "
                  f"err={err:>+6.1f}%  nrna_leak={nrna_obs:>6.1f}  gdna_leak={gdna_obs:>6.1f}")

    # =====================================================================
    section("4. PURE nRNA (TA1=0, NTA1>0, gdna_fraction=0)")
    # =====================================================================
    pure_nrna = df[(df["TA1"] == 0) & (df["NTA1"] > 0) & (df["gdna_fraction"] == 0)]
    print(f"Runs: {len(pure_nrna)}")
    for ss in sorted(pure_nrna.strand_specificity.unique()):
        sub = pure_nrna[pure_nrna.strand_specificity == ss]
        print(f"\n  SS={ss}:")
        for _, row in sub.iterrows():
            nta1 = row["NTA1"]
            nrna_exp = row["nrna_expected"]
            nrna_obs = row["nrna_observed"]
            mrna_obs = row["total_mrna_observed"]
            gdna_obs = row["gdna_observed"]
            err = (nrna_obs - nrna_exp) / nrna_exp * 100 if nrna_exp > 0 else 0
            print(f"    NTA1={nta1:>5.0f}: nrna_exp={nrna_exp:>7.0f} nrna_obs={nrna_obs:>9.1f} "
                  f"err={err:>+6.1f}%  mrna_leak={mrna_obs:>6.1f}  gdna_leak={gdna_obs:>6.1f}")

    # =====================================================================
    section("5. mRNA + gDNA (NTA1=0, gdna_fraction>0, TA1>0)")
    # =====================================================================
    mrna_gdna = df[(df["NTA1"] == 0) & (df["gdna_fraction"] > 0) & (df["TA1"] > 0)]
    print(f"Runs: {len(mrna_gdna)}")
    # Focus on high TA1 to see gDNA error clearly
    for ta1 in [100, 500, 1000]:
        sub = mrna_gdna[mrna_gdna.TA1 == ta1]
        if len(sub) == 0:
            continue
        print(f"\n  TA1={ta1}:")
        for ss in sorted(sub.strand_specificity.unique()):
            ssub = sub[sub.strand_specificity == ss]
            print(f"    SS={ss}:")
            for _, row in ssub.iterrows():
                gf = row["gdna_fraction"]
                gdna_exp = row["gdna_expected"]
                gdna_obs = row["gdna_observed"]
                mrna_exp = row["total_mrna_expected"]
                mrna_obs = row["total_mrna_observed"]
                nrna_obs = row["nrna_observed"]
                g_err = (gdna_obs - gdna_exp) / gdna_exp * 100 if gdna_exp > 0 else 0
                m_err = (mrna_obs - mrna_exp) / mrna_exp * 100 if mrna_exp > 0 else 0
                print(f"      gf={gf}: gdna {gdna_exp:.0f}→{gdna_obs:.0f} ({g_err:+.1f}%)  "
                      f"mrna {mrna_exp:.0f}→{mrna_obs:.1f} ({m_err:+.1f}%)  "
                      f"nrna_leak={nrna_obs:.1f}")

    # =====================================================================
    section("6. nRNA + gDNA (TA1=0, NTA1>0, gdna_fraction>0)")
    # =====================================================================
    nrna_gdna = df[(df["TA1"] == 0) & (df["NTA1"] > 0) & (df["gdna_fraction"] > 0)]
    print(f"Runs: {len(nrna_gdna)}")
    for nta1 in [100, 500, 1000]:
        sub = nrna_gdna[nrna_gdna.NTA1 == nta1]
        if len(sub) == 0:
            continue
        print(f"\n  NTA1={nta1}:")
        for ss in sorted(sub.strand_specificity.unique()):
            ssub = sub[sub.strand_specificity == ss]
            print(f"    SS={ss}:")
            for _, row in ssub.iterrows():
                gf = row["gdna_fraction"]
                gdna_exp = row["gdna_expected"]
                gdna_obs = row["gdna_observed"]
                nrna_exp = row["nrna_expected"]
                nrna_obs = row["nrna_observed"]
                mrna_obs = row["total_mrna_observed"]
                g_err = (gdna_obs - gdna_exp) / gdna_exp * 100 if gdna_exp > 0 else 0
                n_err = (nrna_obs - nrna_exp) / nrna_exp * 100 if nrna_exp > 0 else 0
                print(f"      gf={gf}: gdna {gdna_exp:.0f}→{gdna_obs:.0f} ({g_err:+.1f}%)  "
                      f"nrna {nrna_exp:.0f}→{nrna_obs:.1f} ({n_err:+.1f}%)  "
                      f"mrna_leak={mrna_obs:.1f}")

    # =====================================================================
    section("7. mRNA + nRNA + gDNA (all three present)")
    # =====================================================================
    all_three = df[(df["TA1"] > 0) & (df["NTA1"] > 0) & (df["gdna_fraction"] > 0)]
    print(f"Runs: {len(all_three)}")
    # Show TA1=500, NTA1=500 as representative
    for ta1 in [100, 500]:
        for nta1 in [100, 500]:
            sub = all_three[(all_three.TA1 == ta1) & (all_three.NTA1 == nta1)]
            if len(sub) == 0:
                continue
            print(f"\n  TA1={ta1}, NTA1={nta1}:")
            for ss in sorted(sub.strand_specificity.unique()):
                ssub = sub[sub.strand_specificity == ss]
                print(f"    SS={ss}:")
                for _, row in ssub.iterrows():
                    gf = row["gdna_fraction"]
                    gdna_exp = row["gdna_expected"]
                    gdna_obs = row["gdna_observed"]
                    mrna_exp = row["total_mrna_expected"]
                    mrna_obs = row["total_mrna_observed"]
                    nrna_exp = row["nrna_expected"]
                    nrna_obs = row["nrna_observed"]
                    g_err = (gdna_obs - gdna_exp) / gdna_exp * 100 if gdna_exp > 0 else 0
                    m_err = (mrna_obs - mrna_exp) / mrna_exp * 100 if mrna_exp > 0 else 0
                    n_err = (nrna_obs - nrna_exp) / nrna_exp * 100 if nrna_exp > 0 else 0
                    print(f"      gf={gf}: gdna {gdna_exp:.0f}→{gdna_obs:.0f} ({g_err:+.1f}%)  "
                          f"mrna {mrna_exp:.0f}→{mrna_obs:.1f} ({m_err:+.1f}%)  "
                          f"nrna {nrna_exp:.0f}→{nrna_obs:.1f} ({n_err:+.1f}%)")

    # =====================================================================
    section("8. AGGREGATE ERROR STATISTICS")
    # =====================================================================
    # Only runs with nonzero expectations
    has_g = df[df["gdna_expected"] > 0]
    has_n = df[df["nrna_expected"] > 0]
    has_m = df[df["total_mrna_expected"] > 0]

    print("  gDNA error (runs with gDNA > 0):")
    print(f"    N={len(has_g)}")
    print(f"    Mean |err|: {has_g.gdna_abs_rel_err.mean():.1%}")
    print(f"    Median |err|: {has_g.gdna_abs_rel_err.median():.1%}")
    print(f"    Mean signed err: {has_g.gdna_rel_err.mean():+.1%}")
    print(f"    P90 |err|: {has_g.gdna_abs_rel_err.quantile(0.9):.1%}")
    print(f"    P99 |err|: {has_g.gdna_abs_rel_err.quantile(0.99):.1%}")

    print(f"\n  nRNA error (runs with nRNA > 0):")
    print(f"    N={len(has_n)}")
    print(f"    Mean |err|: {has_n.nrna_abs_rel_err.mean():.1%}")
    print(f"    Median |err|: {has_n.nrna_abs_rel_err.median():.1%}")
    print(f"    Mean signed err: {has_n.nrna_signed_err.mean():+.1%}")

    print(f"\n  mRNA error (runs with mRNA > 0):")
    print(f"    N={len(has_m)}")
    print(f"    Mean |err|: {has_m.mrna_abs_rel_err.mean():.1%}")
    print(f"    Median |err|: {has_m.mrna_abs_rel_err.median():.1%}")
    print(f"    Mean signed err: {((has_m.total_mrna_observed - has_m.total_mrna_expected) / has_m.total_mrna_expected).mean():+.1%}")

    # =====================================================================
    section("9. gDNA ERROR BY SS AND GDNA_FRACTION")
    # =====================================================================
    print("  Mean gDNA signed error (%):")
    pivot = has_g.pivot_table(
        values="gdna_rel_err", index="strand_specificity",
        columns="gdna_fraction", aggfunc="mean"
    ) * 100
    print(pivot.to_string(float_format="%.1f"))

    print("\n\n  Mean gDNA |error| (%):")
    pivot2 = has_g.pivot_table(
        values="gdna_abs_rel_err", index="strand_specificity",
        columns="gdna_fraction", aggfunc="mean"
    ) * 100
    print(pivot2.to_string(float_format="%.1f"))

    # =====================================================================
    section("10. FALSE POSITIVE LEAKAGE")
    # =====================================================================
    # Cases where component should be zero but isn't
    no_mrna = df[df["total_mrna_expected"] == 0]
    no_nrna = df[df["nrna_expected"] == 0]
    no_gdna = df[df["gdna_expected"] == 0]

    print(f"  mRNA false positives (mrna_expected=0, N={len(no_mrna)}):")
    if len(no_mrna) > 0:
        fp_mrna = no_mrna[no_mrna["total_mrna_observed"] > 1]
        print(f"    Runs with mrna_observed > 1: {len(fp_mrna)} ({len(fp_mrna)/len(no_mrna)*100:.0f}%)")
        if len(fp_mrna) > 0:
            print(f"    Mean false mrna: {fp_mrna.total_mrna_observed.mean():.1f}")
            print(f"    Max false mrna: {fp_mrna.total_mrna_observed.max():.1f}")

    print(f"\n  nRNA false positives (nrna_expected=0, N={len(no_nrna)}):")
    if len(no_nrna) > 0:
        fp_nrna = no_nrna[no_nrna["nrna_observed"] > 1]
        print(f"    Runs with nrna_observed > 1: {len(fp_nrna)} ({len(fp_nrna)/len(no_nrna)*100:.0f}%)")
        if len(fp_nrna) > 0:
            print(f"    Mean false nrna: {fp_nrna.nrna_observed.mean():.1f}")
            print(f"    Max false nrna: {fp_nrna.nrna_observed.max():.1f}")

    print(f"\n  gDNA false positives (gdna_expected=0, N={len(no_gdna)}):")
    if len(no_gdna) > 0:
        fp_gdna = no_gdna[no_gdna["gdna_observed"] > 1]
        print(f"    Runs with gdna_observed > 1: {len(fp_gdna)} ({len(fp_gdna)/len(no_gdna)*100:.0f}%)")
        if len(fp_gdna) > 0:
            print(f"    Mean false gdna: {fp_gdna.gdna_observed.mean():.1f}")
            print(f"    Max false gdna: {fp_gdna.gdna_observed.max():.1f}")

    # =====================================================================
    section("11. WORST CASES — TOP 20 gDNA ERRORS")
    # =====================================================================
    worst = has_g.nlargest(20, "gdna_abs_rel_err")
    cols = ["strand_specificity", "gdna_fraction", "TA1", "NTA1",
            "gdna_expected", "gdna_observed", "gdna_rel_err",
            "nrna_expected", "nrna_observed", "total_mrna_expected", "total_mrna_observed"]
    for _, row in worst.iterrows():
        print(f"  SS={row.strand_specificity} gf={row.gdna_fraction} TA1={row.TA1:.0f} NTA1={row.NTA1:.0f}: "
              f"gdna {row.gdna_expected:.0f}→{row.gdna_observed:.0f} ({row.gdna_rel_err:+.1%})  "
              f"nrna {row.nrna_expected:.0f}→{row.nrna_observed:.1f}  "
              f"mrna {row.total_mrna_expected:.0f}→{row.total_mrna_observed:.1f}")

    # =====================================================================
    section("12. CONSERVATION CHECK — TOTAL FRAGMENTS")
    # =====================================================================
    df["total_assigned"] = df["total_mrna_observed"] + df["nrna_observed"] + df["gdna_observed"]
    df["total_expected_all"] = df["total_mrna_expected"] + df["nrna_expected"] + df["gdna_expected"]
    # Intergenic fragments are NOT assigned, so compare to (actual - intergenic - chimeric)
    df["assignable"] = df["n_fragments_actual"] - df["n_intergenic"] - df["n_chimeric"]
    df["conservation_err"] = (df["total_assigned"] - df["assignable"]) / df["assignable"]

    print(f"  Total assigned vs assignable fragments:")
    print(f"    Mean conservation error: {df.conservation_err.mean():.6f}")
    print(f"    Max |conservation error|: {df.conservation_err.abs().max():.6f}")
    print(f"    Runs with |error| > 0.01: {(df.conservation_err.abs() > 0.01).sum()}")

    # =====================================================================
    section("13. UNEXPRESSED TRANSCRIPT LEAKAGE (TD, TE)")
    # =====================================================================
    print(f"  TD (single-exon, unexpressed):")
    print(f"    Runs with TD_observed > 0.1: {(df.TD_observed > 0.1).sum()}")
    print(f"    Max TD_observed: {df.TD_observed.max():.2f}")
    print(f"\n  TE (2-exon neg strand, unexpressed):")
    print(f"    Runs with TE_observed > 0.1: {(df.TE_observed > 0.1).sum()}")
    print(f"    Max TE_observed: {df.TE_observed.max():.2f}")

    # =====================================================================
    section("14. SIMPLEST CASE: PURE gDNA ONLY (TA1=0, NTA1=0, gf>0)")
    # =====================================================================
    # This is the TRIVIAL case — no RNA at all, everything is gDNA
    pure_gdna = df[(df["TA1"] == 0) & (df["NTA1"] == 0) & (df["gdna_fraction"] > 0)]
    print(f"Runs: {len(pure_gdna)}")
    print("  All fragments should be classified as gDNA, nothing as mRNA/nRNA.")
    for ss in sorted(pure_gdna.strand_specificity.unique()):
        sub = pure_gdna[pure_gdna.strand_specificity == ss]
        for _, row in sub.iterrows():
            gf = row["gdna_fraction"]
            ge = row["gdna_expected"]
            go = row["gdna_observed"]
            no = row["nrna_observed"]
            mo = row["total_mrna_observed"]
            frags = row["n_fragments_actual"]
            inter = row["n_intergenic"]
            err = (go - ge) / ge * 100 if ge > 0 else 0
            print(f"  SS={ss} gf={gf}: {frags:.0f} frags ({inter:.0f} intergenic)  "
                  f"gdna {ge:.0f}→{go:.0f} ({err:+.1f}%)  "
                  f"mrna_fp={mo:.1f}  nrna_fp={no:.1f}")


if __name__ == "__main__":
    main()
