#!/usr/bin/env python3
"""Analyze hyperparameter ablation sweep results.

Reads sweep_results.tsv from each HP ablation directory and produces
summary tables showing how gDNA, nRNA, and mRNA errors vary with each HP.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path("/Users/mkiyer/Downloads/rigel_runs/hp_ablation")

SWEEPS = {
    "strand_symmetry_kappa": {
        "values": [2.0, 3.0, 6.0, 12.0, 24.0, 50.0],
        "default": 6.0,
    },
    "prior_alpha": {
        "values": [0.001, 0.01, 0.1, 1.0, 5.0, 10.0],
        "default": 0.01,
    },
    "prune_threshold": {
        "values": [-1.0, 0.01, 0.1, 1.0, 5.0, 10.0],
        "default": 0.1,
    },
}

PATTERN_NAMES = {
    0: "pure_mRNA",
    1: "mixed_A_mod_nRNA",
    2: "convergent_heavy_nRNA",
    3: "everything_moderate",
}


def load_sweep(hp_name):
    path = BASE / hp_name / "sweep_results.tsv"
    df = pd.read_csv(path, sep="\t")

    # Compute gDNA relative error (signed): (obs - exp) / exp
    df["gdna_rel_err"] = (df["gdna_observed"] - df["gdna_expected"]) / df["gdna_expected"]
    df["gdna_abs_rel_err"] = df["gdna_rel_err"].abs()

    # nrna_rel_err already exists in the TSV
    df["nrna_abs_rel_err"] = df["nrna_rel_err"].abs()

    # total_mrna_rel_err already exists (signed)
    df["mrna_abs_rel_err"] = df["total_mrna_rel_err"].abs()

    # Assign pattern index based on row order: 96 rows = 4 patterns × 2 SS × 2 gdna × 6 HP
    # Sweep order: patterns (outermost) → SS → gdna_frac → HP values (innermost)
    n_hp = len(SWEEPS[hp_name]["values"])
    n_gdna = 2
    n_ss = 2
    n_pat = 4
    inner = n_hp * n_gdna * n_ss  # per pattern
    df["pattern_idx"] = np.arange(len(df)) // inner
    df["pattern"] = df["pattern_idx"].map(PATTERN_NAMES)

    return df


def print_separator(char="=", width=100):
    print(char * width)


def analyze_one(hp_name, info):
    df = load_sweep(hp_name)
    hp_col = hp_name
    default_val = info["default"]

    print_separator()
    print(f"  HYPERPARAMETER: {hp_name}  (default = {default_val})")
    print_separator()

    # ── 1. Overall summary by HP value ──
    print(f"\n{'─'*80}")
    print("1. OVERALL METRICS BY HP VALUE (mean across all patterns/SS/gdna)")
    print(f"{'─'*80}")
    agg = (
        df.groupby(hp_col)
        .agg(
            gdna_err_mean=("gdna_abs_rel_err", "mean"),
            gdna_err_med=("gdna_abs_rel_err", "median"),
            nrna_err_mean=("nrna_abs_rel_err", "mean"),
            nrna_err_med=("nrna_abs_rel_err", "median"),
            mrna_err_mean=("mrna_abs_rel_err", "mean"),
            mrna_err_med=("mrna_abs_rel_err", "median"),
        )
        .reset_index()
    )
    for _, row in agg.iterrows():
        marker = " <-- default" if row[hp_col] == default_val else ""
        print(
            f"  {hp_col}={row[hp_col]:<8.3f}  "
            f"gDNA |err|: mean={row['gdna_err_mean']:6.1%} med={row['gdna_err_med']:6.1%}  "
            f"nRNA |err|: mean={row['nrna_err_mean']:6.1%} med={row['nrna_err_med']:6.1%}  "
            f"mRNA |err|: mean={row['mrna_err_mean']:6.1%} med={row['mrna_err_med']:6.1%}"
            f"{marker}"
        )

    # ── 2. Breakdown by SS and gdna_fraction ──
    print(f"\n{'─'*80}")
    print("2. gDNA |error| BY HP VALUE × SS × gdna_fraction")
    print(f"{'─'*80}")
    pivot = (
        df.groupby([hp_col, "strand_specificity", "gdna_fraction"])["gdna_abs_rel_err"]
        .mean()
        .reset_index()
    )
    for ss in sorted(df["strand_specificity"].unique()):
        for gf in sorted(df["gdna_fraction"].unique()):
            sub = pivot[(pivot["strand_specificity"] == ss) & (pivot["gdna_fraction"] == gf)]
            vals_str = "  ".join(
                f"{row[hp_col]:.3g}:{row['gdna_abs_rel_err']:5.1%}"
                for _, row in sub.iterrows()
            )
            print(f"  SS={ss}  gDNA_frac={gf}:  {vals_str}")

    # ── 3. Breakdown by pattern ──
    print(f"\n{'─'*80}")
    print("3. gDNA |error| BY HP VALUE × PATTERN (mean over SS/gdna_frac)")
    print(f"{'─'*80}")
    for pidx in sorted(PATTERN_NAMES.keys()):
        pname = PATTERN_NAMES[pidx]
        sub = df[df["pattern_idx"] == pidx]
        agg_p = sub.groupby(hp_col)["gdna_abs_rel_err"].mean()
        vals_str = "  ".join(f"{v:.3g}:{agg_p[v]:5.1%}" for v in info["values"])
        print(f"  {pname:<28s}  {vals_str}")

    # ── 4. nRNA error by pattern ──
    print(f"\n{'─'*80}")
    print("4. nRNA |error| BY HP VALUE × PATTERN (mean over SS/gdna_frac)")
    print(f"{'─'*80}")
    for pidx in sorted(PATTERN_NAMES.keys()):
        pname = PATTERN_NAMES[pidx]
        sub = df[df["pattern_idx"] == pidx]
        agg_p = sub.groupby(hp_col)["nrna_abs_rel_err"].mean()
        vals_str = "  ".join(f"{v:.3g}:{agg_p[v]:5.1%}" for v in info["values"])
        print(f"  {pname:<28s}  {vals_str}")

    # ── 5. Sensitivity metric ──
    print(f"\n{'─'*80}")
    print("5. SENSITIVITY (range of mean |error| across HP values)")
    print(f"{'─'*80}")
    gdna_range = agg["gdna_err_mean"].max() - agg["gdna_err_mean"].min()
    nrna_range = agg["nrna_err_mean"].max() - agg["nrna_err_mean"].min()
    mrna_range = agg["mrna_err_mean"].max() - agg["mrna_err_mean"].min()
    best_gdna = agg.loc[agg["gdna_err_mean"].idxmin()]
    best_nrna = agg.loc[agg["nrna_err_mean"].idxmin()]
    print(f"  gDNA |err| range: {gdna_range:6.1%}  (best={best_gdna['gdna_err_mean']:5.1%} @ {hp_col}={best_gdna[hp_col]:.3g})")
    print(f"  nRNA |err| range: {nrna_range:6.1%}  (best={best_nrna['nrna_err_mean']:5.1%} @ {hp_col}={best_nrna[hp_col]:.3g})")
    print(f"  mRNA |err| range: {mrna_range:6.1%}")

    return agg


def main():
    results = {}
    for hp_name, info in SWEEPS.items():
        results[hp_name] = analyze_one(hp_name, info)
        print("\n")

    # ── Cross-HP comparison ──
    print_separator("=")
    print("  CROSS-HYPERPARAMETER SENSITIVITY SUMMARY")
    print_separator("=")
    print(f"\n  {'HP':<28s}  {'gDNA range':>12s}  {'nRNA range':>12s}  {'mRNA range':>12s}  {'Best gDNA val':>14s}  {'Best gDNA err':>14s}")
    print(f"  {'─'*28}  {'─'*12}  {'─'*12}  {'─'*12}  {'─'*14}  {'─'*14}")
    for hp_name, agg in results.items():
        hp_col = hp_name
        gdna_r = agg["gdna_err_mean"].max() - agg["gdna_err_mean"].min()
        nrna_r = agg["nrna_err_mean"].max() - agg["nrna_err_mean"].min()
        mrna_r = agg["mrna_err_mean"].max() - agg["mrna_err_mean"].min()
        best = agg.loc[agg["gdna_err_mean"].idxmin()]
        print(
            f"  {hp_name:<28s}  {gdna_r:>11.1%}  {nrna_r:>11.1%}  {mrna_r:>11.1%}"
            f"  {best[hp_col]:>14.3g}  {best['gdna_err_mean']:>13.1%}"
        )

    print(f"\n  Interpretation:")
    print(f"  - Small range → error is NOT sensitive to this HP (structural issue)")
    print(f"  - Large range → error IS tunable via this HP")
    print()


if __name__ == "__main__":
    main()
