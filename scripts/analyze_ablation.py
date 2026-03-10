#!/usr/bin/env python3
"""Analyze extreme ablation study results."""
import pandas as pd
import numpy as np
import sys

def load_data():
    gdna = pd.read_csv('sweep_results/ablation_gdna/sweep_results.tsv', sep='\t')
    nrna = pd.read_csv('sweep_results/ablation_nrna/sweep_results.tsv', sep='\t')
    mrna = pd.read_csv('sweep_results/ablation_mrna/sweep_results.tsv', sep='\t')
    comb = pd.read_csv('sweep_results/ablation_combined/sweep_results.tsv', sep='\t')
    return gdna, nrna, mrna, comb

def analyze_gdna(gdna):
    print("=" * 95)
    print("ABLATION 1: gDNA EXTREME STRESS TEST")
    print("Fixed: TA1=100, TA4=3, NTA=100 | Vary: gDNA=[0..10000], SS=[0.9,1.0], kappa=[2,6]")
    print("=" * 95)
    
    for ss in [1.0, 0.9]:
        for k in [6.0, 2.0]:
            sub = gdna[(gdna.strand_specificity == ss) & (gdna.strand_symmetry_kappa == k)]
            label = f"SS={ss}, kappa={k:.0f}"
            if ss == 1.0 and k == 6.0:
                label += " (symmetric gDNA + penalty = HARDEST)"
            elif ss == 1.0 and k == 2.0:
                label += " (symmetric gDNA, no penalty)"
            elif ss == 0.9 and k == 6.0:
                label += " (asymmetric gDNA + penalty)"
            else:
                label += " (asymmetric gDNA, no penalty)"
            print(f"\n--- {label} ---")
            print(f"  {'gDNA':>8s} | {'mRNA_err':>9s} | {'gDNA_diff':>10s} | {'nRNA_diff':>10s} | {'TA1_err':>9s} | {'gDNA_exp':>10s} | {'gDNA_obs':>10s}")
            print(f"  {'-'*8} | {'-'*9} | {'-'*10} | {'-'*10} | {'-'*9} | {'-'*10} | {'-'*10}")
            for _, r in sub.sort_values('gdna').iterrows():
                print(f"  {int(r.gdna):>8d} | {r.total_mrna_rel_err:>8.1%} | "
                      f"{r.gdna_abs_diff:>10.0f} | {r.nrna_abs_diff:>10.0f} | "
                      f"{r.TA1_rel_err:>8.1%} | {r.gdna_expected:>10.0f} | {r.gdna_observed:>10.0f}")

def analyze_nrna(nrna):
    print("\n" + "=" * 95)
    print("ABLATION 2: nRNA EXTREME STRESS TEST")
    print("Fixed: TA1=100, TA4=3 | Vary: NTA=[0..10000], gDNA=[0,100], kappa=[2,6]")
    print("=" * 95)
    
    for g in [0, 100]:
        for k in [6.0, 2.0]:
            sub = nrna[(nrna.gdna == g) & (nrna.strand_symmetry_kappa == k)]
            label = f"gDNA={g}, kappa={k:.0f}"
            print(f"\n--- {label} ---")
            print(f"  {'NTA':>8s} | {'mRNA_err':>9s} | {'nRNA_diff':>10s} | {'nRNA_exp':>10s} | {'nRNA_obs':>10s} | {'gDNA_diff':>10s} | {'TA1_err':>9s}")
            print(f"  {'-'*8} | {'-'*9} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*9}")
            for _, r in sub.sort_values('NTA').iterrows():
                print(f"  {int(r.NTA):>8d} | {r.total_mrna_rel_err:>8.1%} | "
                      f"{r.nrna_abs_diff:>10.0f} | {r.nrna_expected:>10.0f} | {r.nrna_observed:>10.0f} | "
                      f"{r.gdna_abs_diff:>10.0f} | {r.TA1_rel_err:>8.1%}")

def analyze_mrna(mrna):
    print("\n" + "=" * 95)
    print("ABLATION 3: mRNA DYNAMIC RANGE")
    print("Fixed: TA4=3, NTA=100, gDNA=100, SS=1.0 | Vary: TA1=[1..10000], kappa=[2,6]")
    print("=" * 95)
    
    for k in [6.0, 2.0]:
        sub = mrna[mrna.strand_symmetry_kappa == k]
        print(f"\n--- kappa={k:.0f} ---")
        print(f"  {'TA1':>8s} | {'mRNA_err':>9s} | {'TA1_err':>9s} | {'TA1_exp':>10s} | {'TA1_obs':>10s} | "
              f"{'gDNA_diff':>10s} | {'nRNA_diff':>10s}")
        print(f"  {'-'*8} | {'-'*9} | {'-'*9} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10}")
        for _, r in sub.sort_values('TA1').iterrows():
            print(f"  {int(r.TA1):>8d} | {r.total_mrna_rel_err:>8.1%} | {r.TA1_rel_err:>8.1%} | "
                  f"{r.TA1_expected:>10.0f} | {r.TA1_observed:>10.0f} | "
                  f"{r.gdna_abs_diff:>10.0f} | {r.nrna_abs_diff:>10.0f}")

def analyze_combined(comb):
    print("\n" + "=" * 95)
    print("ABLATION 4: COMBINED EXTREMES (NTA x gDNA matrix)")
    print("Fixed: TA1=100, TA4=3, SS=1.0, kappa=6.0 | Vary: NTA=[0..10000], gDNA=[0..10000]")
    print("=" * 95)
    
    # Pivot tables
    nta_vals = sorted(comb.NTA.unique())
    gdna_vals = sorted(comb.gdna.unique())
    
    print("\n--- mRNA Total Relative Error (%) ---")
    print(f"{'NTA\\gDNA':>10s}", end="")
    for g in gdna_vals:
        print(f" | {int(g):>8d}", end="")
    print()
    for n in nta_vals:
        print(f"{int(n):>10d}", end="")
        for g in gdna_vals:
            row = comb[(comb.NTA == n) & (comb.gdna == g)]
            if len(row) == 1:
                v = row.total_mrna_rel_err.values[0]
                print(f" | {v:>7.1%}", end="")
            else:
                print(f" | {'N/A':>8s}", end="")
        print()
    
    print("\n--- gDNA Absolute Diff ---")
    print(f"{'NTA\\gDNA':>10s}", end="")
    for g in gdna_vals:
        print(f" | {int(g):>8d}", end="")
    print()
    for n in nta_vals:
        print(f"{int(n):>10d}", end="")
        for g in gdna_vals:
            row = comb[(comb.NTA == n) & (comb.gdna == g)]
            if len(row) == 1:
                v = row.gdna_abs_diff.values[0]
                print(f" | {v:>8.0f}", end="")
            else:
                print(f" | {'N/A':>8s}", end="")
        print()
    
    print("\n--- nRNA Absolute Diff ---")
    print(f"{'NTA\\gDNA':>10s}", end="")
    for g in gdna_vals:
        print(f" | {int(g):>8d}", end="")
    print()
    for n in nta_vals:
        print(f"{int(n):>10d}", end="")
        for g in gdna_vals:
            row = comb[(comb.NTA == n) & (comb.gdna == g)]
            if len(row) == 1:
                v = row.nrna_abs_diff.values[0]
                print(f" | {v:>8.0f}", end="")
            else:
                print(f" | {'N/A':>8s}", end="")
        print()
    
    print("\n--- TA1 Relative Error (%) ---")
    print(f"{'NTA\\gDNA':>10s}", end="")
    for g in gdna_vals:
        print(f" | {int(g):>8d}", end="")
    print()
    for n in nta_vals:
        print(f"{int(n):>10d}", end="")
        for g in gdna_vals:
            row = comb[(comb.NTA == n) & (comb.gdna == g)]
            if len(row) == 1:
                v = row.TA1_rel_err.values[0]
                print(f" | {v:>7.1%}", end="")
            else:
                print(f" | {'N/A':>8s}", end="")
        print()

def convergence_analysis(gdna, nrna, mrna, comb):
    print("\n" + "=" * 95)
    print("CONVERGENCE & ERROR SOURCE ANALYSIS")
    print("=" * 95)
    
    # Check if model converges by looking at total_rna_rel_err (mass conservation)
    all_data = pd.concat([
        gdna.assign(ablation='gDNA'),
        nrna.assign(ablation='nRNA'),
        mrna.assign(ablation='mRNA'),
        comb.assign(ablation='combined')
    ])
    
    print("\n--- Mass Conservation (total_rna_rel_err) by ablation ---")
    for abl in ['gDNA', 'nRNA', 'mRNA', 'combined']:
        sub = all_data[all_data.ablation == abl]
        print(f"  {abl:>10s}: mean={sub.total_rna_rel_err.mean():.4%}, "
              f"max={sub.total_rna_rel_err.max():.4%}, "
              f"worst_condition: gdna={int(sub.loc[sub.total_rna_rel_err.idxmax(), 'gdna'])}, "
              f"NTA={int(sub.loc[sub.total_rna_rel_err.idxmax(), 'NTA'])}")
    
    print("\n--- Where does mRNA error exceed 20%? ---")
    bad = all_data[all_data.total_mrna_rel_err > 0.20].sort_values('total_mrna_rel_err', ascending=False)
    if len(bad) > 0:
        print(f"  {len(bad)} conditions with mRNA error > 20%:")
        for _, r in bad.head(15).iterrows():
            print(f"    [{r.ablation:>8s}] TA1={int(r.TA1)}, NTA={int(r.NTA)}, gDNA={int(r.gdna)}, "
                  f"SS={r.strand_specificity}, k={r.strand_symmetry_kappa:.0f} => mRNA_err={r.total_mrna_rel_err:.1%}, "
                  f"gDNA_diff={r.gdna_abs_diff:.0f}, nRNA_diff={r.nrna_abs_diff:.0f}")
    
    print("\n--- Error decomposition: gDNA siphon vs nRNA bleed ---")
    print("  (positive gdna_obs - gdna_exp means gDNA OVER-estimated, siphoning FROM mRNA)")
    print("  (positive nrna_obs - nrna_exp means nRNA OVER-estimated, siphoning FROM mRNA)")
    
    # For extreme gDNA cases with SS=1.0, kappa=6
    extreme_g = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == 6.0) & (gdna.gdna >= 1000)]
    print("\n  Extreme gDNA (SS=1.0, k=6, gDNA>=1000):")
    for _, r in extreme_g.sort_values('gdna').iterrows():
        gdna_err_dir = "OVER" if r.gdna_observed > r.gdna_expected else "UNDER"
        nrna_err_dir = "OVER" if r.nrna_observed > r.nrna_expected else "UNDER"
        print(f"    gDNA={int(r.gdna):>6d}: gDNA {gdna_err_dir} by {abs(r.gdna_observed - r.gdna_expected):.0f}, "
              f"nRNA {nrna_err_dir} by {abs(r.nrna_observed - r.nrna_expected):.0f}, "
              f"mRNA_err={r.total_mrna_rel_err:.1%}")
    
    # For extreme nRNA cases with gDNA=100
    extreme_n = nrna[(nrna.gdna == 100) & (nrna.strand_symmetry_kappa == 6.0) & (nrna.NTA >= 1000)]
    print("\n  Extreme nRNA (gDNA=100, k=6, NTA>=1000):")
    for _, r in extreme_n.sort_values('NTA').iterrows():
        gdna_err_dir = "OVER" if r.gdna_observed > r.gdna_expected else "UNDER"
        nrna_err_dir = "OVER" if r.nrna_observed > r.nrna_expected else "UNDER"
        print(f"    NTA={int(r.NTA):>6d}: nRNA {nrna_err_dir} by {abs(r.nrna_observed - r.nrna_expected):.0f}, "
              f"gDNA {gdna_err_dir} by {abs(r.gdna_observed - r.gdna_expected):.0f}, "
              f"mRNA_err={r.total_mrna_rel_err:.1%}")
    
    # Ratio analysis: error as fraction of contamination level
    print("\n--- Error scaling: does error grow sub-linearly, linearly, or super-linearly? ---")
    g_k6_ss10 = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == 6.0) & (gdna.gdna > 0)]
    print("\n  gDNA scaling (SS=1.0, k=6):")
    prev_gdna = None
    prev_err = None
    for _, r in g_k6_ss10.sort_values('gdna').iterrows():
        ratio = r.gdna_abs_diff / r.gdna_expected if r.gdna_expected > 0 else float('inf')
        if prev_gdna is not None and prev_err is not None:
            gdna_mult = r.gdna / prev_gdna
            err_mult = r.total_mrna_rel_err / prev_err if prev_err > 0 else float('inf')
            print(f"    gDNA={int(r.gdna):>6d}: mRNA_err={r.total_mrna_rel_err:.1%}, "
                  f"gDNA_err_ratio={ratio:.2%}, gDNA_x{gdna_mult:.0f} => err_x{err_mult:.1f}")
        else:
            print(f"    gDNA={int(r.gdna):>6d}: mRNA_err={r.total_mrna_rel_err:.1%}, "
                  f"gDNA_err_ratio={ratio:.2%}")
        prev_gdna = r.gdna
        prev_err = r.total_mrna_rel_err

if __name__ == '__main__':
    gdna, nrna, mrna, comb = load_data()
    analyze_gdna(gdna)
    analyze_nrna(nrna)
    analyze_mrna(mrna)
    analyze_combined(comb)
    convergence_analysis(gdna, nrna, mrna, comb)
