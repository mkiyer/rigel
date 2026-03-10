#!/usr/bin/env python3
"""Compare v1 and v2 ablation results."""
import pandas as pd
import numpy as np

df = pd.read_csv('sweep_results/ablation_v2/sweep_results.tsv', sep='\t')
df['pool_mrna_exp'] = df['TA1_expected'] + df['TA2_expected'] + df['TA3_expected'] + df['TA4_expected']
df['pool_mrna_obs'] = df['TA1_observed'] + df['TA2_observed'] + df['TA3_observed'] + df['TA4_observed']
df['mrna_err'] = df['pool_mrna_obs'] - df['pool_mrna_exp']
df['mrna_rel_err'] = abs(df['mrna_err']) / df['pool_mrna_exp'].clip(lower=1)
df['nrna_err'] = df['nrna_observed'] - df['nrna_expected']

# 1. kappa=4 focused analysis at gDNA=1000
print('=== kappa=4 at gDNA=1000 (nRNA present) ===')
mask = (df['strand_symmetry_kappa'] == 4.0) & (df['gdna'] == 1000) & (df['NTA'] > 0)
sub = df[mask]
nrna_rel = abs(sub['nrna_err'] / sub['nrna_expected'].clip(lower=1))
print(f'  n={len(sub)}, mRNA_med_rel={sub["mrna_rel_err"].median():.4f}, nRNA_med_rel={nrna_rel.median():.4f}')
print(f'  nRNA direction: over={sum(sub["nrna_err"]>0)}, under={sum(sub["nrna_err"]<0)}')
print(f'  nRNA_err_med={sub["nrna_err"].median():.0f}')

# 2. Compare v1 vs v2 for the SAME kappa values (2 and 6)
df1 = pd.read_csv('sweep_results/ablation/sweep_results.tsv', sep='\t')
df1['pool_mrna_exp'] = df1['TA1_expected'] + df1['TA2_expected'] + df1['TA3_expected'] + df1['TA4_expected']
df1['pool_mrna_obs'] = df1['TA1_observed'] + df1['TA2_observed'] + df1['TA3_observed'] + df1['TA4_observed']
df1['mrna_rel_err'] = abs(df1['pool_mrna_obs'] - df1['pool_mrna_exp']) / df1['pool_mrna_exp'].clip(lower=1)

print()
print('=== V1 vs V2 matched comparison (kappa=2 and kappa=6 only) ===')
for k in [2.0, 6.0]:
    m1 = df1['strand_symmetry_kappa'] == k
    m2 = df['strand_symmetry_kappa'] == k
    contam1 = m1 & ((df1['NTA'] > 0) | (df1['gdna'] > 0))
    contam2 = m2 & ((df['NTA'] > 0) | (df['gdna'] > 0))
    print(f'  kappa={k}: V1 mRNA_med={df1.loc[contam1,"mrna_rel_err"].median():.4f}  V2 mRNA_med={df.loc[contam2,"mrna_rel_err"].median():.4f}')

# 3. Per-transcript accuracy at worst case (gDNA=1000, kappa=4)
print()
print('=== Per-transcript errors at gDNA=1000, kappa=4, NTA=100 ===')
mask = (df['strand_symmetry_kappa'] == 4.0) & (df['gdna'] == 1000) & (df['NTA'] == 100)
sub = df[mask]
for t in ['TA1', 'TA2', 'TA3', 'TA4']:
    rel = abs(sub[f'{t}_observed'] - sub[f'{t}_expected']) / sub[f'{t}_expected'].clip(lower=1)
    print(f'  {t}: med_rel_err={rel.median():.4f}, max_rel_err={rel.max():.4f}')

# 4. nRNA at low levels with gDNA
print()
print('=== nRNA=10 survival rates ===')
for g in [0, 10, 100, 1000]:
    for k in [2.0, 4.0, 6.0]:
        mask = (df['NTA'] == 10) & (df['gdna'] == g) & (df['strand_symmetry_kappa'] == k)
        sub = df[mask]
        zeroed = (sub['nrna_observed'] < 1.0).sum()
        print(f'  gDNA={g:4d} kappa={k}: {zeroed}/{len(sub)} runs nRNA zeroed, nRNA_med={sub["nrna_observed"].median():.1f} (exp={sub["nrna_expected"].median():.0f})')

# 5. gDNA accuracy improvement
print()
print('=== gDNA accuracy V1 vs V2 ===')
for g in [10, 100, 1000]:
    m1 = (df1['gdna'] == g)
    m2 = (df['gdna'] == g)
    gdna_rel1 = abs(df1.loc[m1, 'gdna_observed'] - df1.loc[m1, 'gdna_expected']) / df1.loc[m1, 'gdna_expected'].clip(lower=1)
    gdna_rel2 = abs(df.loc[m2, 'gdna_observed'] - df.loc[m2, 'gdna_expected']) / df.loc[m2, 'gdna_expected'].clip(lower=1)
    print(f'  gDNA={g:4d}: V1 gdna_err_med={gdna_rel1.median():.4f}  V2 gdna_err_med={gdna_rel2.median():.4f}')

# 6. Mass flow comparison V1 vs V2
print()
print('=== Mass flow V1 vs V2 (nRNA+gDNA scenarios) ===')
df1['nrna_err'] = df1['nrna_observed'] - df1['nrna_expected']
df1['mrna_err'] = (df1['TA1_observed'] + df1['TA2_observed'] + df1['TA3_observed'] + df1['TA4_observed']) - df1['pool_mrna_exp']

m1 = (df1['NTA'] > 0) & (df1['gdna'] > 0)
m2 = (df['NTA'] > 0) & (df['gdna'] > 0)
print(f'  V1: nRNA_err_med={df1.loc[m1,"nrna_err"].median():+.0f}, mRNA_err_med={df1.loc[m1,"mrna_err"].median():+.0f}')
print(f'  V2: nRNA_err_med={df.loc[m2,"nrna_err"].median():+.0f}, mRNA_err_med={df.loc[m2,"mrna_err"].median():+.0f}')

# By kappa
for k in [2.0, 6.0]:
    m1k = m1 & (df1['strand_symmetry_kappa'] == k)
    m2k = m2 & (df['strand_symmetry_kappa'] == k)
    print(f'  kappa={k} V1: nRNA_err_med={df1.loc[m1k,"nrna_err"].median():+.0f}, mRNA_err_med={df1.loc[m1k,"mrna_err"].median():+.0f}')
    print(f'  kappa={k} V2: nRNA_err_med={df.loc[m2k,"nrna_err"].median():+.0f}, mRNA_err_med={df.loc[m2k,"mrna_err"].median():+.0f}')
