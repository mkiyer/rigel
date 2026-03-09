#!/usr/bin/env python3
"""Analyze sweep results from kappa/pseudo grid search."""
import pandas as pd
import numpy as np
import sys
import os

def analyze_grid(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Compute composite score
    df['mrna_abs_diff'] = abs(df['total_mrna_observed'] - df['total_mrna_expected'])
    df['composite'] = df['mrna_abs_diff'] + df['gdna_abs_diff'] + df['nrna_abs_diff']
    
    pd.set_option('display.width', 200)
    pd.set_option('display.max_rows', 60)
    
    print("=== COMPOSITE SCORE (|mRNA_diff| + |gDNA_diff| + |nRNA_diff|, lower=better) ===")
    pivot = df.pivot_table(values='composite', 
                           index='strand_symmetry_kappa',
                           columns='strand_symmetry_pseudo',
                           aggfunc='first')
    print(pivot.round(0).to_string())
    print()
    
    print("=== mRNA RELATIVE ERROR % ===")
    pivot_mrna = df.pivot_table(values='total_mrna_rel_err',
                                index='strand_symmetry_kappa',
                                columns='strand_symmetry_pseudo',
                                aggfunc='first')
    print((pivot_mrna * 100).round(1).to_string())
    print()
    
    print("=== gDNA ABSOLUTE DIFF (from expected) ===")
    pivot_gdna = df.pivot_table(values='gdna_abs_diff',
                                index='strand_symmetry_kappa',
                                columns='strand_symmetry_pseudo',
                                aggfunc='first')
    print(pivot_gdna.round(0).to_string())
    print()
    
    print("=== nRNA ABSOLUTE DIFF (from expected) ===")
    pivot_nrna = df.pivot_table(values='nrna_abs_diff',
                                index='strand_symmetry_kappa',
                                columns='strand_symmetry_pseudo',
                                aggfunc='first')
    print(pivot_nrna.round(0).to_string())
    print()
    
    # Find best
    best_idx = df['composite'].idxmin()
    best = df.loc[best_idx]
    print(f"=== BEST CONFIG ===")
    print(f"  kappa={best['strand_symmetry_kappa']}, pseudo={best['strand_symmetry_pseudo']}")
    print(f"  composite={best['composite']:.0f}")
    print(f"  mRNA: exp={best['total_mrna_expected']:.0f}, obs={best['total_mrna_observed']:.0f}, rel_err={best['total_mrna_rel_err']*100:.1f}%")
    print(f"  gDNA: exp={best['gdna_expected']:.0f}, obs={best['gdna_observed']:.0f}, diff={best['gdna_abs_diff']:.0f}")
    print(f"  nRNA: exp={best['nrna_expected']:.0f}, obs={best['nrna_observed']:.0f}, diff={best['nrna_abs_diff']:.0f}")


def analyze_baseline(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['mrna_abs_diff'] = abs(df['total_mrna_observed'] - df['total_mrna_expected'])
    df['composite'] = df['mrna_abs_diff'] + df['gdna_abs_diff'] + df['nrna_abs_diff']
    
    pd.set_option('display.width', 200)
    
    # Only runs with gdna > 0
    gdna_df = df[df['gdna'] > 0].copy()
    
    print("=== BASELINE COMPARISON: kappa=2.0 (disabled) vs kappa=6.0 ===")
    print()
    
    for kappa in [2.0, 6.0]:
        k_df = gdna_df[gdna_df['strand_symmetry_kappa'] == kappa]
        print(f"--- kappa={kappa} {'(disabled)' if kappa == 2.0 else '(default)'} ---")
        print(f"  Mean mRNA rel error: {k_df['total_mrna_rel_err'].mean()*100:.1f}%")
        print(f"  Mean gDNA abs diff:  {k_df['gdna_abs_diff'].mean():.0f}")
        print(f"  Mean nRNA abs diff:  {k_df['nrna_abs_diff'].mean():.0f}")
        print(f"  Mean composite:      {k_df['composite'].mean():.0f}")
        print()
    
    # Paired comparison
    k2 = gdna_df[gdna_df['strand_symmetry_kappa'] == 2.0].sort_values(['TA1','TA2','TA3','TA4','NTA','gdna']).reset_index(drop=True)
    k6 = gdna_df[gdna_df['strand_symmetry_kappa'] == 6.0].sort_values(['TA1','TA2','TA3','TA4','NTA','gdna']).reset_index(drop=True)
    
    print("=== PER-SCENARIO COMPARISON ===")
    print(f"{'TA1':>4} {'TA2':>4} {'TA3':>4} {'TA4':>4} {'NTA':>4} | {'mRNA_err_k2':>12} {'mRNA_err_k6':>12} {'Δ':>8} | {'gDNA_diff_k2':>12} {'gDNA_diff_k6':>12} | {'nRNA_diff_k2':>12} {'nRNA_diff_k6':>12}")
    print("-" * 120)
    for i in range(len(k2)):
        r2, r6 = k2.iloc[i], k6.iloc[i]
        delta_mrna = r6['total_mrna_rel_err'] - r2['total_mrna_rel_err']
        better = "✓" if abs(r6['total_mrna_rel_err']) < abs(r2['total_mrna_rel_err']) else "✗"
        print(f"{r2['TA1']:>4.0f} {r2['TA2']:>4.0f} {r2['TA3']:>4.0f} {r2['TA4']:>4.0f} {r2['NTA']:>4.0f} | "
              f"{r2['total_mrna_rel_err']*100:>11.1f}% {r6['total_mrna_rel_err']*100:>11.1f}% {delta_mrna*100:>7.1f}% | "
              f"{r2['gdna_abs_diff']:>12.0f} {r6['gdna_abs_diff']:>12.0f} | "
              f"{r2['nrna_abs_diff']:>12.0f} {r6['nrna_abs_diff']:>12.0f} {better}")


def analyze_robustness(tsv_path):
    df = pd.read_csv(tsv_path, sep='\t')
    df['mrna_abs_diff'] = abs(df['total_mrna_observed'] - df['total_mrna_expected'])
    df['composite'] = df['mrna_abs_diff'] + df['gdna_abs_diff'] + df['nrna_abs_diff']
    
    pd.set_option('display.width', 200)
    
    gdna_df = df[df['gdna'] > 0].copy()
    kappas = sorted(gdna_df['strand_symmetry_kappa'].unique())
    
    # Summary by kappa across all conditions
    print("=== ROBUSTNESS ACROSS ALL gdna>0 CONDITIONS ===")
    for kappa in kappas:
        k_df = gdna_df[gdna_df['strand_symmetry_kappa'] == kappa]
        tag = "(disabled)" if kappa == 2.0 else "(default)" if kappa == 6.0 else ""
        print(f"\n--- kappa={kappa} {tag} ---")
        print(f"  Mean mRNA rel error: {k_df['total_mrna_rel_err'].mean()*100:.1f}%")
        print(f"  Max mRNA rel error:  {k_df['total_mrna_rel_err'].max()*100:.1f}%")
        print(f"  Mean gDNA abs diff:  {k_df['gdna_abs_diff'].mean():.0f}")
        print(f"  Mean nRNA abs diff:  {k_df['nrna_abs_diff'].mean():.0f}")
        print(f"  Mean composite:      {k_df['composite'].mean():.0f}")
    
    # Breakdown by SS
    print("\n\n=== BY STRAND SPECIFICITY ===")
    for ss in sorted(gdna_df['strand_specificity'].unique()):
        ss_df = gdna_df[gdna_df['strand_specificity'] == ss]
        print(f"\nSS={ss}:")
        for kappa in kappas:
            k_df = ss_df[ss_df['strand_symmetry_kappa'] == kappa]
            if len(k_df) > 0:
                print(f"  kappa={kappa}: mRNA_err={k_df['total_mrna_rel_err'].mean()*100:.1f}%, gdna_diff={k_df['gdna_abs_diff'].mean():.0f}, nrna_diff={k_df['nrna_abs_diff'].mean():.0f}, composite={k_df['composite'].mean():.0f}")
    
    # Breakdown by gDNA level
    print("\n\n=== BY gDNA LEVEL ===")
    for gdna_level in sorted(gdna_df['gdna'].unique()):
        gl_df = gdna_df[gdna_df['gdna'] == gdna_level]
        print(f"\ngDNA={gdna_level}:")
        for kappa in kappas:
            k_df = gl_df[gl_df['strand_symmetry_kappa'] == kappa]
            if len(k_df) > 0:
                print(f"  kappa={kappa}: mRNA_err={k_df['total_mrna_rel_err'].mean()*100:.1f}%, gdna_diff={k_df['gdna_abs_diff'].mean():.0f}, nrna_diff={k_df['nrna_abs_diff'].mean():.0f}, composite={k_df['composite'].mean():.0f}")
    
    # Breakdown by NTA
    print("\n\n=== BY NASCENT RNA ===")
    for nta in sorted(gdna_df['NTA'].unique()):
        nta_df = gdna_df[gdna_df['NTA'] == nta]
        print(f"\nNTA={nta}:")
        for kappa in kappas:
            k_df = nta_df[nta_df['strand_symmetry_kappa'] == kappa]
            if len(k_df) > 0:
                print(f"  kappa={kappa}: mRNA_err={k_df['total_mrna_rel_err'].mean()*100:.1f}%, gdna_diff={k_df['gdna_abs_diff'].mean():.0f}, nrna_diff={k_df['nrna_abs_diff'].mean():.0f}, composite={k_df['composite'].mean():.0f}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["grid", "baseline", "robustness"])
    parser.add_argument("--dir", default=None, help="Override sweep results directory")
    args = parser.parse_args()

    defaults = {
        "grid": "sweep_results/kappa_pseudo_grid",
        "baseline": "sweep_results/baseline_comparison",
        "robustness": "sweep_results/robustness",
    }
    base_dir = args.dir if args.dir else defaults[args.mode]
    tsv_path = os.path.join(base_dir, "sweep_results.tsv")

    if args.mode == "grid":
        analyze_grid(tsv_path)
    elif args.mode == "baseline":
        analyze_baseline(tsv_path)
    elif args.mode == "robustness":
        analyze_robustness(tsv_path)
