#!/usr/bin/env python3
"""Final synthesis analysis: root cause and ratios."""
import pandas as pd
import numpy as np

def load_data():
    gdna = pd.read_csv('sweep_results/ablation_gdna/sweep_results.tsv', sep='\t')
    nrna = pd.read_csv('sweep_results/ablation_nrna/sweep_results.tsv', sep='\t')
    mrna = pd.read_csv('sweep_results/ablation_mrna/sweep_results.tsv', sep='\t')
    comb = pd.read_csv('sweep_results/ablation_combined/sweep_results.tsv', sep='\t')
    return gdna, nrna, mrna, comb

def genic_fraction_analysis(gdna):
    """Key insight: how many gDNA fragments are *genic* (compete with mRNA in EM)?"""
    print("=" * 100)
    print("GENIC FRACTION ANALYSIS")
    print("At extreme gDNA, most fragments are intergenic (ignored by EM)")
    print("Only genic gDNA fragments compete with mRNA/nRNA in the EM")
    print("=" * 100)
    
    g = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == 6.0)]
    print(f"\n  {'gDNA':>7s} | {'n_frag':>8s} | {'n_intergenic':>13s} | {'n_genic':>8s} | {'genic%':>8s} | "
          f"{'RNA_exp':>8s} | {'gDNA_genic':>11s} | {'noise_ratio':>12s} | {'penalty_leak':>13s} | {'leak/RNA%':>10s}")
    for _, r in g.sort_values('gdna').iterrows():
        n_genic = r.n_fragments_actual - r.n_intergenic - r.n_chimeric
        genic_pct = n_genic / r.n_fragments_actual if r.n_fragments_actual > 0 else 0
        rna_exp = r.total_mrna_expected + r.nrna_expected
        gdna_genic = n_genic - rna_exp  # Approximate genic gDNA
        noise_ratio = gdna_genic / rna_exp if rna_exp > 0 else 0
        # Penalty leak = gDNA under-estimation (mass transferred to nRNA)
        penalty_leak = -(r.gdna_observed - r.gdna_expected)  # positive = mass removed from gDNA
        leak_rna_pct = penalty_leak / rna_exp if rna_exp > 0 else 0
        print(f"  {int(r.gdna):>7d} | {int(r.n_fragments_actual):>8d} | "
              f"{int(r.n_intergenic):>13d} | {n_genic:>8.0f} | {genic_pct:>7.1%} | "
              f"{rna_exp:>8.0f} | {gdna_genic:>11.0f} | {noise_ratio:>11.1f}:1 | "
              f"{penalty_leak:>13.0f} | {leak_rna_pct:>9.0%}")

def per_transcript_analysis(gdna):
    """Individual transcript accuracy at extreme gDNA levels."""
    print("\n" + "=" * 100)
    print("PER-TRANSCRIPT ANALYSIS: How does each transcript fare at extremes?")
    print("(SS=1.0, kappa=6)")
    print("=" * 100)
    
    g = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == 6.0)]
    print(f"\n  {'gDNA':>7s} | {'TA1_exp':>8s} | {'TA1_obs':>8s} | {'TA1_err':>8s} | "
          f"{'TA4_exp':>8s} | {'TA4_obs':>8s} | {'TA4_err':>8s} | {'mRNA_err':>9s}")
    for _, r in g.sort_values('gdna').iterrows():
        print(f"  {int(r.gdna):>7d} | {r.TA1_expected:>8.0f} | {r.TA1_observed:>8.0f} | {r.TA1_rel_err:>7.1%} | "
              f"{r.TA4_expected:>8.0f} | {r.TA4_observed:>8.0f} | {r.TA4_rel_err:>7.1%} | {r.total_mrna_rel_err:>8.1%}")

def kappa2_vs_kappa6_detailed(gdna):
    """Detailed comparison showing WHY k=2 beats k=6 at gDNA=10000."""
    print("\n" + "=" * 100)
    print("WHY kappa=2 BEATS kappa=6 AT EXTREME gDNA (SS=1.0)")
    print("=" * 100)
    
    for g_val in [5000, 10000]:
        print(f"\n--- gDNA={g_val} ---")
        for k in [2.0, 6.0]:
            r = gdna[(gdna.strand_specificity == 1.0) & (gdna.strand_symmetry_kappa == k) & (gdna.gdna == g_val)].iloc[0]
            mrna_d = r.total_mrna_observed - r.total_mrna_expected
            gdna_d = r.gdna_observed - r.gdna_expected
            nrna_d = r.nrna_observed - r.nrna_expected
            print(f"  k={k:.0f}: mRNA({r.total_mrna_expected:.0f}→{r.total_mrna_observed:.0f}, {mrna_d:+.0f}) | "
                  f"gDNA({r.gdna_expected:.0f}→{r.gdna_observed:.0f}, {gdna_d:+.0f}) | "
                  f"nRNA({r.nrna_expected:.0f}→{r.nrna_observed:.0f}, {nrna_d:+.0f}) | "
                  f"mRNA_err={r.total_mrna_rel_err:.1%}")

def nrna_protective_effect(nrna):
    """How does high nRNA affect the gDNA estimation?"""
    print("\n" + "=" * 100)
    print("nRNA PROTECTIVE EFFECT: Does more nRNA help or hurt mRNA accuracy?")
    print("(kappa=6, gDNA=100)")
    print("=" * 100)
    
    n = nrna[(nrna.gdna == 100) & (nrna.strand_symmetry_kappa == 6.0)]
    print(f"\n  {'NTA':>7s} | {'mRNA_frac':>10s} | {'nRNA_frac':>10s} | {'gDNA_frac':>10s} | "
          f"{'mRNA_err':>9s} | {'TA1_err':>9s} | {'nRNA_err':>9s}")
    for _, r in n.sort_values('NTA').iterrows():
        total_exp = r.total_mrna_expected + r.nrna_expected + r.gdna_expected
        mrna_f = r.total_mrna_expected / total_exp if total_exp > 0 else 0
        nrna_f = r.nrna_expected / total_exp if total_exp > 0 else 0
        gdna_f = r.gdna_expected / total_exp if total_exp > 0 else 0
        print(f"  {int(r.NTA):>7d} | {mrna_f:>9.1%} | {nrna_f:>9.1%} | {gdna_f:>9.1%} | "
              f"{r.total_mrna_rel_err:>8.1%} | {r.TA1_rel_err:>8.1%} | {r.nrna_rel_err:>8.1%}")

def combined_worst_cases(comb):
    """Identify the absolute worst conditions and understand why."""
    print("\n" + "=" * 100)
    print("COMBINED EXTREMES: WORST CASES ANALYZED")
    print("=" * 100)
    
    c = comb.copy()
    c['mrna_delta'] = c.total_mrna_observed - c.total_mrna_expected
    c['gdna_delta'] = c.gdna_observed - c.gdna_expected
    c['nrna_delta'] = c.nrna_observed - c.nrna_expected
    
    worst = c.sort_values('total_mrna_rel_err', ascending=False).head(10)
    print(f"\n  {'NTA':>6s} | {'gDNA':>6s} | {'mRNA_err':>9s} | {'mRNA_d':>8s} | {'gDNA_d':>8s} | "
          f"{'nRNA_d':>8s} | {'n_inter':>8s} | {'TA1_err':>8s}")
    for _, r in worst.iterrows():
        print(f"  {int(r.NTA):>6d} | {int(r.gdna):>6d} | {r.total_mrna_rel_err:>8.1%} | "
              f"{r.mrna_delta:>+8.0f} | {r.gdna_delta:>+8.0f} | {r.nrna_delta:>+8.0f} | "
              f"{int(r.n_intergenic):>8d} | {r.TA1_rel_err:>7.1%}")

def summary_statistics(gdna, nrna, mrna, comb):
    """Summary statistics across all ablations."""
    print("\n" + "=" * 100)
    print("SUMMARY STATISTICS")
    print("=" * 100)
    
    all_data = pd.concat([
        gdna.assign(ablation='gDNA'),
        nrna.assign(ablation='nRNA'),
        mrna.assign(ablation='mRNA'),
        comb.assign(ablation='combined')
    ])
    
    # Thresholds
    for thresh in [0.05, 0.10, 0.20, 0.50]:
        n_bad = (all_data.total_mrna_rel_err > thresh).sum()
        n_total = len(all_data)
        pct = n_bad / n_total
        print(f"  mRNA error > {thresh:.0%}: {n_bad}/{n_total} ({pct:.0%}) conditions")
    
    print(f"\n  Overall median mRNA error: {all_data.total_mrna_rel_err.median():.1%}")
    print(f"  Overall mean mRNA error: {all_data.total_mrna_rel_err.mean():.1%}")
    
    # By ablation
    for abl in ['gDNA', 'nRNA', 'mRNA', 'combined']:
        sub = all_data[all_data.ablation == abl]
        print(f"  [{abl:>8s}] median={sub.total_mrna_rel_err.median():.1%}, "
              f"mean={sub.total_mrna_rel_err.mean():.1%}, "
              f"max={sub.total_mrna_rel_err.max():.1%}")

if __name__ == '__main__':
    gdna, nrna, mrna, comb = load_data()
    genic_fraction_analysis(gdna)
    per_transcript_analysis(gdna)
    kappa2_vs_kappa6_detailed(gdna)
    nrna_protective_effect(nrna)
    combined_worst_cases(comb)
    summary_statistics(gdna, nrna, mrna, comb)
