#!/usr/bin/env python
"""analyze_ablation.py — Analyze comprehensive ablation sweep results.

Usage:
    python scripts/analyze_ablation.py sweep_results/ablation/sweep_results.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def load_results(tsv_path: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    # Compute pool-level signed errors
    t_ids = [c.replace("_expected", "") for c in df.columns if c.endswith("_expected")
             and c not in ("nrna_expected", "gdna_expected",
                           "total_mrna_expected", "total_rna_expected")]
    df["pool_mrna_exp"] = sum(df[f"{t}_expected"] for t in t_ids)
    df["pool_mrna_obs"] = sum(df[f"{t}_observed"] for t in t_ids)
    df["pool_mrna_err"] = df["pool_mrna_obs"] - df["pool_mrna_exp"]
    df["pool_nrna_err"] = df["nrna_observed"] - df["nrna_expected"]
    df["pool_gdna_err"] = df["gdna_observed"] - df["gdna_expected"]

    # Relative errors (absolute)
    df["mrna_rel_err"] = np.where(
        df["pool_mrna_exp"] > 0,
        np.abs(df["pool_mrna_err"]) / df["pool_mrna_exp"], 0.0)
    df["nrna_rel_abs"] = np.where(
        df["nrna_expected"] > 0,
        np.abs(df["pool_nrna_err"]) / df["nrna_expected"], 0.0)
    df["gdna_rel_abs"] = np.where(
        df["gdna_expected"] > 0,
        np.abs(df["pool_gdna_err"]) / df["gdna_expected"], 0.0)
    # Signed relative
    df["nrna_rel_sgn"] = np.where(
        df["nrna_expected"] > 0,
        df["pool_nrna_err"] / df["nrna_expected"], 0.0)
    df["gdna_rel_sgn"] = np.where(
        df["gdna_expected"] > 0,
        df["pool_gdna_err"] / df["gdna_expected"], 0.0)

    df["ss"] = df["strand_specificity"]
    if "strand_symmetry_kappa" in df.columns:
        df["kappa"] = df["strand_symmetry_kappa"]
    else:
        df["kappa"] = 6.0  # default
    return df


def section(title: str):
    w = 80
    print(f"\n{'=' * w}")
    print(f"  {title}")
    print(f"{'=' * w}")


def report_overall(df: pd.DataFrame):
    section("OVERALL ACCURACY SUMMARY")
    n = len(df)
    has_nrna = df["nrna_expected"] > 0
    has_gdna = df["gdna_expected"] > 0

    print(f"Total runs: {n}")
    print(f"  With nRNA: {has_nrna.sum()}, With gDNA: {has_gdna.sum()}, "
          f"Both: {(has_nrna & has_gdna).sum()}, Clean: {(~has_nrna & ~has_gdna).sum()}")

    for label, mask, col in [
        ("mRNA", pd.Series(True, index=df.index), "mrna_rel_err"),
        ("nRNA", has_nrna, "nrna_rel_abs"),
        ("gDNA", has_gdna, "gdna_rel_abs"),
    ]:
        sub = df[mask]
        if len(sub) == 0:
            continue
        v = sub[col].values
        print(f"  {label} |rel err|: median={np.median(v):.4f}, "
              f"mean={np.mean(v):.4f}, p95={np.percentile(v, 95):.4f}, "
              f"max={np.max(v):.4f}")


def report_grid(df: pd.DataFrame):
    section("ACCURACY GRID: nRNA level × gDNA level")
    print("  (median mRNA relative error across all SS/kappa/pattern combos)")
    agg = df.groupby(["NTA", "gdna"]).agg(
        n=("mrna_rel_err", "count"),
        mrna_med=("mrna_rel_err", "median"),
        mrna_max=("mrna_rel_err", "max"),
        nrna_med=("nrna_rel_abs", "median"),
        gdna_med=("gdna_rel_abs", "median"),
    ).reset_index()

    print(f"{'NTA':>6} {'gDNA':>6} {'n':>4}  "
          f"{'mRNA_med':>9} {'mRNA_max':>9}  "
          f"{'nRNA_med':>9} {'gDNA_med':>9}")
    print("-" * 70)
    for _, r in agg.iterrows():
        nrna_s = f"{r['nrna_med']:>9.4f}" if r["NTA"] > 0 else f"{'---':>9}"
        gdna_s = f"{r['gdna_med']:>9.4f}" if r["gdna"] > 0 else f"{'---':>9}"
        print(f"{r['NTA']:>6.0f} {r['gdna']:>6.0f} {r['n']:>4.0f}  "
              f"{r['mrna_med']:>9.4f} {r['mrna_max']:>9.4f}  "
              f"{nrna_s} {gdna_s}")


def report_by_ss(df: pd.DataFrame):
    section("ACCURACY BY STRAND SPECIFICITY (contaminated runs only)")
    has_contam = (df["nrna_expected"] > 0) | (df["gdna_expected"] > 0)
    sub = df[has_contam]
    agg = sub.groupby("ss").agg(
        n=("mrna_rel_err", "count"),
        mrna_med=("mrna_rel_err", "median"),
        mrna_p95=("mrna_rel_err", lambda x: np.percentile(x, 95)),
        nrna_med=("nrna_rel_abs", "median"),
        gdna_med=("gdna_rel_abs", "median"),
    ).reset_index()

    print(f"{'SS':>6} {'n':>4}  {'mRNA_med':>9} {'mRNA_p95':>9}  "
          f"{'nRNA_med':>9} {'gDNA_med':>9}")
    print("-" * 60)
    for _, r in agg.iterrows():
        print(f"{r['ss']:>6.2f} {r['n']:>4.0f}  "
              f"{r['mrna_med']:>9.4f} {r['mrna_p95']:>9.4f}  "
              f"{r['nrna_med']:>9.4f} {r['gdna_med']:>9.4f}")


def report_by_kappa(df: pd.DataFrame):
    section("KAPPA COMPARISON (contaminated runs only)")
    has_contam = (df["nrna_expected"] > 0) | (df["gdna_expected"] > 0)
    sub = df[has_contam]
    agg = sub.groupby("kappa").agg(
        n=("mrna_rel_err", "count"),
        mrna_med=("mrna_rel_err", "median"),
        mrna_p95=("mrna_rel_err", lambda x: np.percentile(x, 95)),
        nrna_med=("nrna_rel_abs", "median"),
        gdna_med=("gdna_rel_abs", "median"),
    ).reset_index()

    print(f"{'kappa':>6} {'n':>4}  {'mRNA_med':>9} {'mRNA_p95':>9}  "
          f"{'nRNA_med':>9} {'gDNA_med':>9}")
    print("-" * 60)
    for _, r in agg.iterrows():
        print(f"{r['kappa']:>6.1f} {r['n']:>4.0f}  "
              f"{r['mrna_med']:>9.4f} {r['mrna_p95']:>9.4f}  "
              f"{r['nrna_med']:>9.4f} {r['gdna_med']:>9.4f}")


def report_kappa_x_gdna(df: pd.DataFrame):
    section("KAPPA × gDNA INTERACTION (runs with nRNA present)")
    sub = df[df["nrna_expected"] > 0]
    if len(sub) == 0:
        return
    agg = sub.groupby(["gdna", "kappa"]).agg(
        n=("mrna_rel_err", "count"),
        mrna_med=("mrna_rel_err", "median"),
        nrna_med=("nrna_rel_abs", "median"),
        gdna_med=("gdna_rel_abs", "median"),
        nrna_sgn=("nrna_rel_sgn", "median"),
    ).reset_index()

    print(f"{'gDNA':>6} {'kappa':>6} {'n':>4}  "
          f"{'mRNA_med':>9} {'nRNA_med':>9} {'gDNA_med':>9}  {'nRNA_dir':>9}")
    print("-" * 65)
    for _, r in agg.iterrows():
        gdna_s = f"{r['gdna_med']:>9.4f}" if r["gdna"] > 0 else f"{'---':>9}"
        print(f"{r['gdna']:>6.0f} {r['kappa']:>6.1f} {r['n']:>4.0f}  "
              f"{r['mrna_med']:>9.4f} {r['nrna_med']:>9.4f} {gdna_s}  "
              f"{r['nrna_sgn']:>+9.4f}")


def report_mass_flow(df: pd.DataFrame):
    section("MASS FLOW ANALYSIS (signed errors: where do fragments leak?)")
    has_nrna = df["nrna_expected"] > 0
    has_gdna = df["gdna_expected"] > 0

    for name, mask in [
        ("Clean (no contam)", ~has_nrna & ~has_gdna),
        ("nRNA only", has_nrna & ~has_gdna),
        ("gDNA only", ~has_nrna & has_gdna),
        ("nRNA + gDNA", has_nrna & has_gdna),
    ]:
        sub = df[mask]
        if len(sub) == 0:
            continue
        print(f"\n  {name} ({len(sub)} runs):")
        print(f"    mRNA signed err: med={sub['pool_mrna_err'].median():+.0f}, "
              f"mean={sub['pool_mrna_err'].mean():+.0f}, "
              f"max_abs={sub['pool_mrna_err'].abs().max():.0f}")
        if sub["nrna_expected"].sum() > 0:
            print(f"    nRNA signed err: med={sub['pool_nrna_err'].median():+.0f}, "
                  f"mean={sub['pool_nrna_err'].mean():+.0f}")
        if sub["gdna_expected"].sum() > 0:
            print(f"    gDNA signed err: med={sub['pool_gdna_err'].median():+.0f}, "
                  f"mean={sub['pool_gdna_err'].mean():+.0f}")
        # Conservation: all errors should sum to ~ 0
        total_err = sub["pool_mrna_err"] + sub["pool_nrna_err"] + sub["pool_gdna_err"]
        print(f"    Conservation residual: med={total_err.median():+.0f}, "
              f"max_abs={total_err.abs().max():.0f}")


def report_nrna_direction(df: pd.DataFrame):
    section("nRNA ESTIMATION DIRECTION (over vs under)")
    has_nrna = df["nrna_expected"] > 0
    sub = df[has_nrna].copy()
    n_over = (sub["pool_nrna_err"] > 0).sum()
    n_under = (sub["pool_nrna_err"] < 0).sum()
    n_exact = (sub["pool_nrna_err"] == 0).sum()
    print(f"  Over: {n_over}, Under: {n_under}, Exact: {n_exact}")

    print(f"\n  {'NTA':>6} {'gDNA':>6} {'n':>4}  {'nRNA_err_med':>12} {'mRNA_err_med':>12}  "
          f"{'%over':>6}")
    print("  " + "-" * 60)
    for nta in sorted(sub["NTA"].unique()):
        for gdna_lev in sorted(sub["gdna"].unique()):
            s = sub[(sub["NTA"] == nta) & (sub["gdna"] == gdna_lev)]
            if len(s) == 0:
                continue
            pct_over = (s["pool_nrna_err"] > 0).mean() * 100
            print(f"  {nta:>6.0f} {gdna_lev:>6.0f} {len(s):>4}  "
                  f"{s['pool_nrna_err'].median():>+12.0f} "
                  f"{s['pool_mrna_err'].median():>+12.0f}  "
                  f"{pct_over:>5.0f}%")


def report_pattern(df: pd.DataFrame):
    section("TRANSCRIPT PATTERN SENSITIVITY (contaminated only)")
    t_cols = [c for c in ["TA1", "TA2", "TA3", "TA4"] if c in df.columns]
    if not t_cols:
        print("  (no pattern columns)")
        return
    df = df.copy()
    df["pattern"] = df[t_cols].idxmax(axis=1)
    has_contam = (df["nrna_expected"] > 0) | (df["gdna_expected"] > 0)
    sub = df[has_contam]

    agg = sub.groupby("pattern").agg(
        n=("mrna_rel_err", "count"),
        mrna_med=("mrna_rel_err", "median"),
        mrna_max=("mrna_rel_err", "max"),
    ).reset_index()

    print(f"{'Pattern':>8} {'n':>4}  {'mRNA_med':>9} {'mRNA_max':>9}")
    print("-" * 35)
    for _, r in agg.iterrows():
        print(f"{r['pattern']:>8} {r['n']:>4.0f}  "
              f"{r['mrna_med']:>9.4f} {r['mrna_max']:>9.4f}")


def report_worst(df: pd.DataFrame, n: int = 20):
    section(f"TOP {n} WORST mRNA RELATIVE ERRORS")
    worst = df.nlargest(n, "mrna_rel_err")
    keys = ["mrna_rel_err", "NTA", "gdna", "ss", "kappa",
            "pool_mrna_exp", "pool_mrna_obs",
            "nrna_expected", "nrna_observed",
            "gdna_expected", "gdna_observed"]
    avail = [c for c in keys if c in worst.columns]
    pd.set_option("display.float_format", lambda x: f"{x:.2f}")
    print(worst[avail].to_string(index=False))
    pd.reset_option("display.float_format")


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <sweep_results.tsv>")
        sys.exit(1)

    df = load_results(sys.argv[1])

    report_overall(df)
    report_grid(df)
    report_by_ss(df)
    report_by_kappa(df)
    report_kappa_x_gdna(df)
    report_mass_flow(df)
    report_nrna_direction(df)
    report_pattern(df)
    report_worst(df)

    print(f"\n{'=' * 80}")
    print("  ANALYSIS COMPLETE")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
