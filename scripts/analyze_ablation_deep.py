#!/usr/bin/env python3
"""Deep-dive analysis: error flow, mass transfer, and root causes."""
import pandas as pd
import numpy as np

def load_data():
    gdna = pd.read_csv('sweep_results/ablation_gdna/sweep_results.tsv', sep='\t')
    nrna = pd.read_csv('sweep_results/ablation_nrna/sweep_results.tsv', sep='\t')
    mrna = pd.read_csv('sweep_results/ablation_mrna/sweep_results.tsv', sep='\t')
    comb = pd.read_csv('sweep_results/ablation_combined/sweep_results.tsv', sep='\t')
    return gdna, nrna, mrna, comb

def error_flow_analysis(gdna):
    """Where does misattributed mass flow?"""
    print("=" * 100)
    print("ERROR FLOW ANALYSIS: Where does misattributed mass go?")
    print("=" * 100)
    
    # For each gDNA level, compute: mRNA_delta, gDNA_delta, nRNA_delta
    # delta = observed - expected (positive = over-estimated, negative = under-estimated)
    for ss in [1.0, 0.9]:
        for k in [6.0, 2.0]:
            sub = gdna[(gdna.strand_specificity == ss) & (gdna.strand_symmetry_kappa == k)]
            label = f"SS={ss}, kappa={k:.0f}"
            print(f"\n--- {label} ---")
            print(f"  {'gDNA':>8s} | {'mRNA_delta':>12s} | {'gDNA_delta':>12s} | {'nRNA_delta':>12s} | {'SUM':>10s} | {'mRNA%':>8s}")
            print(f"  {'-'*8} | {'-'*12} | {'-'*12} | {'-'*12} | {'-'*10} | {'-'*8}")
            for _, r in sub.sort_values('gdna').iterrows():
                mrna_delta = r.total_mrna_observed - r.total_mrna_expected
                gdna_delta = r.gdna_observed - r.gdna_expected
                nrna_delta = r.nrna_observed - r.nrna_expected
                total = mrna_delta + gdna_delta + nrna_delta
                mrna_pct = r.total_mrna_rel_err
                print(f"  {int(r.gdna):>8d} | {mrna_delta:>+12.0f} | {gdna_delta:>+12.0f} | "
                      f"{nrna_delta:>+12.0f} | {total:>+10.0f} | {mrna_pct:>7.1%}")

def signal_to_noise_analysis(gdna, mrna):
    """How does mRNA signal-to-noise affect accuracy?"""
    print("\n" + "=" * 100)
    print("SIGNAL-TO-NOISE ANALYSIS")
    print("(mRNA fraction of total fragments vs. estimation accuracy)")
    print("=" * 100)
    
    # gDNA ablation: compute mRNA fraction
    g_k6_ss10 = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == 6.0)]
    print(f"\n--- gDNA sweep (SS=1.0, k=6): mRNA diluted by increasing gDNA ---")
    print(f"  {'gDNA':>8s} | {'mRNA_frac':>10s} | {'gDNA_frac':>10s} | {'nRNA_frac':>10s} | {'mRNA_err':>9s} | {'TA1_err':>9s}")
    print(f"  {'-'*8} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*9} | {'-'*9}")
    for _, r in g_k6_ss10.sort_values('gdna').iterrows():
        total_exp = r.total_mrna_expected + r.gdna_expected + r.nrna_expected
        mrna_frac = r.total_mrna_expected / total_exp if total_exp > 0 else 0
        gdna_frac = r.gdna_expected / total_exp if total_exp > 0 else 0
        nrna_frac = r.nrna_expected / total_exp if total_exp > 0 else 0
        print(f"  {int(r.gdna):>8d} | {mrna_frac:>9.1%} | {gdna_frac:>9.1%} | "
              f"{nrna_frac:>9.1%} | {r.total_mrna_rel_err:>8.1%} | {r.TA1_rel_err:>8.1%}")
    
    # mRNA ablation: compute mRNA fraction
    m_k6 = mrna[mrna.strand_symmetry_kappa == 6.0]
    print(f"\n--- mRNA sweep (k=6): mRNA grows, gDNA/nRNA fixed ---")
    print(f"  {'TA1':>8s} | {'mRNA_frac':>10s} | {'gDNA_frac':>10s} | {'nRNA_frac':>10s} | {'mRNA_err':>9s} | {'TA1_err':>9s}")
    print(f"  {'-'*8} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*9} | {'-'*9}")
    for _, r in m_k6.sort_values('TA1').iterrows():
        total_exp = r.total_mrna_expected + r.gdna_expected + r.nrna_expected
        mrna_frac = r.total_mrna_expected / total_exp if total_exp > 0 else 0
        gdna_frac = r.gdna_expected / total_exp if total_exp > 0 else 0
        nrna_frac = r.nrna_expected / total_exp if total_exp > 0 else 0
        print(f"  {int(r.TA1):>8d} | {mrna_frac:>9.1%} | {gdna_frac:>9.1%} | "
              f"{nrna_frac:>9.1%} | {r.total_mrna_rel_err:>8.1%} | {r.TA1_rel_err:>8.1%}")

def penalty_effectiveness(gdna):
    """Compare kappa=2 vs kappa=6 penalty impact."""
    print("\n" + "=" * 100)
    print("PENALTY EFFECTIVENESS: kappa=6 (targeted penalty) vs kappa=2 (no penalty)")
    print("=" * 100)
    
    for ss in [1.0, 0.9]:
        print(f"\n--- SS={ss} ---")
        print(f"  {'gDNA':>8s} | {'k2_mRNA%':>9s} | {'k6_mRNA%':>9s} | {'k2_gDNA_d':>10s} | {'k6_gDNA_d':>10s} | "
              f"{'k2_nRNA_d':>10s} | {'k6_nRNA_d':>10s} | {'penalty_helps':>14s}")
        print(f"  {'-'*8} | {'-'*9} | {'-'*9} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*10} | {'-'*14}")
        
        k2 = gdna[(gdna.strand_specificity == ss) & (gdna.strand_symmetry_kappa == 2.0)].set_index('gdna')
        k6 = gdna[(gdna.strand_specificity == ss) & (gdna.strand_symmetry_kappa == 6.0)].set_index('gdna')
        
        for g in sorted(k2.index):
            r2, r6 = k2.loc[g], k6.loc[g]
            helps = "YES" if r6.total_mrna_rel_err < r2.total_mrna_rel_err else "NO"
            if abs(r6.total_mrna_rel_err - r2.total_mrna_rel_err) < 0.005:
                helps = "NEUTRAL"
            print(f"  {int(g):>8d} | {r2.total_mrna_rel_err:>8.1%} | {r6.total_mrna_rel_err:>8.1%} | "
                  f"{r2.gdna_abs_diff:>10.0f} | {r6.gdna_abs_diff:>10.0f} | "
                  f"{r2.nrna_abs_diff:>10.0f} | {r6.nrna_abs_diff:>10.0f} | {helps:>14s}")

def nrna_gdna_interaction(comb):
    """Analyze how nRNA and gDNA interact at extremes."""
    print("\n" + "=" * 100)
    print("nRNA-gDNA INTERACTION ANALYSIS (combined extremes matrix)")
    print("=" * 100)
    
    nta_vals = sorted(comb.NTA.unique())
    gdna_vals = sorted(comb.gdna.unique())
    
    # Error direction matrix
    print("\n--- mRNA over/under-estimation direction ---")
    print(f"{'NTA\\gDNA':>10s}", end="")
    for g in gdna_vals:
        print(f" | {int(g):>10d}", end="")
    print()
    for n in nta_vals:
        print(f"{int(n):>10d}", end="")
        for g in gdna_vals:
            row = comb[(comb.NTA == n) & (comb.gdna == g)]
            if len(row) == 1:
                r = row.iloc[0]
                delta = r.total_mrna_observed - r.total_mrna_expected
                direction = "OVER" if delta > 0 else "UNDER"
                print(f" | {direction:>5s} {abs(delta):>4.0f}", end="")
            else:
                print(f" | {'N/A':>10s}", end="")
        print()
    
    # Where does mass go when gDNA is extreme?
    print("\n--- Mass transfer at extreme gDNA (gDNA=10000 row) ---")
    for _, r in comb[comb.gdna == 10000].sort_values('NTA').iterrows():
        mrna_d = r.total_mrna_observed - r.total_mrna_expected
        gdna_d = r.gdna_observed - r.gdna_expected
        nrna_d = r.nrna_observed - r.nrna_expected
        print(f"  NTA={int(r.NTA):>6d}: mRNA={mrna_d:>+8.0f}, gDNA={gdna_d:>+8.0f}, nRNA={nrna_d:>+8.0f}")
    
    print("\n--- Mass transfer at extreme nRNA (NTA=10000 column) ---")
    for _, r in comb[comb.NTA == 10000].sort_values('gdna').iterrows():
        mrna_d = r.total_mrna_observed - r.total_mrna_expected
        gdna_d = r.gdna_observed - r.gdna_expected
        nrna_d = r.nrna_observed - r.nrna_expected
        print(f"  gDNA={int(r.gdna):>6d}: mRNA={mrna_d:>+8.0f}, gDNA={gdna_d:>+8.0f}, nRNA={nrna_d:>+8.0f}")

def convergence_check(gdna, nrna, mrna, comb):
    """Check if the model actually converges at extreme values."""
    print("\n" + "=" * 100)
    print("CONVERGENCE CHECK: Mass conservation (total_rna_rel_err)")
    print("(Small = good convergence, Large = mass leak)")
    print("=" * 100)
    
    all_data = pd.concat([
        gdna.assign(ablation='gDNA'),
        nrna.assign(ablation='nRNA'),
        mrna.assign(ablation='mRNA'),
        comb.assign(ablation='combined')
    ])
    
    # Most concerning: where total_rna_rel_err is very large
    bad = all_data[all_data.total_rna_rel_err > 0.10].sort_values('total_rna_rel_err', ascending=False)
    print(f"\n  {len(bad)} conditions with >10% total RNA error (mass conservation failure):")
    for _, r in bad.head(20).iterrows():
        total_exp = r.total_rna_expected
        total_obs = r.total_rna_observed
        leak_pct = (total_obs - total_exp) / total_exp if total_exp > 0 else 0
        print(f"    [{r.ablation:>8s}] gDNA={int(r.gdna)}, NTA={int(r.NTA)}, SS={r.strand_specificity}, "
              f"k={r.strand_symmetry_kappa:.0f}: total_rna_err={r.total_rna_rel_err:.1%} "
              f"(expected={total_exp:.0f}, obs={total_obs:.0f}, n_frag={r.n_fragments_actual:.0f})")
    
    # Check n_intergenic and n_chimeric at extremes
    print(f"\n--- Intergenic & chimeric fragment counts at extremes ---")
    extreme = all_data[(all_data.gdna >= 5000) | (all_data.NTA >= 5000)]
    for _, r in extreme.sort_values('total_rna_rel_err', ascending=False).head(10).iterrows():
        print(f"    [{r.ablation:>8s}] gDNA={int(r.gdna)}, NTA={int(r.NTA)}: "
              f"n_intergenic={int(r.n_intergenic)}, n_chimeric={int(r.n_chimeric)}, "
              f"total_rna_err={r.total_rna_rel_err:.1%}")

if __name__ == '__main__':
    gdna, nrna, mrna, comb = load_data()
    error_flow_analysis(gdna)
    signal_to_noise_analysis(gdna, mrna)
    penalty_effectiveness(gdna)
    nrna_gdna_interaction(comb)
    convergence_check(gdna, nrna, mrna, comb)
