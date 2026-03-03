#!/usr/bin/env python3
"""Analyze pristine and realistic benchmark results for the report."""
import pandas as pd
import numpy as np
import json

pristine_dir = '/Users/mkiyer/Downloads/hulkrna_runs/bench_realdata_pristine'
realistic_dir = '/Users/mkiyer/Downloads/hulkrna_runs/bench_realdata_realistic'

def load_csv(path):
    return pd.read_csv(path)

# Load CSVs
prist = load_csv(f'{pristine_dir}/summary.csv')
real = load_csv(f'{realistic_dir}/summary.csv')

print('=== PRISTINE BENCHMARK (gDNA=0, nRNA=0, SS=1.0) ===')
print(f'Rows: {len(prist)}, Columns: {list(prist.columns)}')
print()

# Transcript-level averages per tool
for level in ['transcript', 'gene']:
    level_col = f'{level}_mae' if f'{level}_mae' in prist.columns else None
    if level_col is None:
        # Try alternate column names
        mae_col = [c for c in prist.columns if 'mae' in c.lower() and level[0] in c.lower()]
        if mae_col:
            level_col = mae_col[0]
    
# Let's look at the column structure
print('Columns:', list(prist.columns))
print()
print('Sample rows:')
print(prist.head(3).to_string())
print()

# Group by level and tool to get averages
for level in prist['level'].unique():
    subset = prist[prist['level'] == level]
    print(f'\n=== PRISTINE: {level}-level averages ===')
    avg = subset.groupby('tool')[['mean_abs_error', 'rmse', 'pearson', 'spearman']].mean()
    print(avg.round(4).to_string())

print('\n\n' + '='*70)
print('=== REALISTIC BENCHMARK (gDNA=0.34, nRNA=0.10, SS=0.997) ===')
print(f'Rows: {len(real)}, Columns: {list(real.columns)}')

for level in real['level'].unique():
    subset = real[real['level'] == level]
    print(f'\n=== REALISTIC: {level}-level averages ===')
    avg = subset.groupby('tool')[['mean_abs_error', 'rmse', 'pearson', 'spearman']].mean()
    print(avg.round(4).to_string())

# Per-region comparison pristine vs realistic
print('\n\n=== PER-REGION TRANSCRIPT-LEVEL COMPARISON ===')
print(f'{"Region":12s} | {"Hulk_prist":>12s} {"Sal_prist":>12s} {"Kal_prist":>12s} | {"Hulk_real":>12s} {"Sal_real":>12s} {"Kal_real":>12s}')
print('-' * 100)

regions = prist[prist['level'] == 'transcript']['region'].unique()
for region in regions:
    prow = prist[(prist['level'] == 'transcript') & (prist['region'] == region)]
    rrow = real[(real['level'] == 'transcript') & (real['region'] == region)]
    
    hulk_p = prow[prow['tool'].str.contains('hulkrna')]['mean_abs_error'].values
    sal_p = prow[prow['tool'] == 'salmon']['mean_abs_error'].values
    kal_p = prow[prow['tool'] == 'kallisto']['mean_abs_error'].values
    
    hulk_r = rrow[rrow['tool'].str.contains('hulkrna')]['mean_abs_error'].values
    sal_r = rrow[rrow['tool'] == 'salmon']['mean_abs_error'].values
    kal_r = rrow[rrow['tool'] == 'kallisto']['mean_abs_error'].values
    
    hp = f'{hulk_p[0]:.2f}' if len(hulk_p) > 0 else 'N/A'
    sp = f'{sal_p[0]:.2f}' if len(sal_p) > 0 else 'N/A'
    kp = f'{kal_p[0]:.2f}' if len(kal_p) > 0 else 'N/A'
    hr = f'{hulk_r[0]:.2f}' if len(hulk_r) > 0 else 'N/A'
    sr = f'{sal_r[0]:.2f}' if len(sal_r) > 0 else 'N/A'
    kr = f'{kal_r[0]:.2f}' if len(kal_r) > 0 else 'N/A'
    
    print(f'{region:12s} | {hp:>12s} {sp:>12s} {kp:>12s} | {hr:>12s} {sr:>12s} {kr:>12s}')

# Gene-level comparison
print('\n\n=== PER-REGION GENE-LEVEL COMPARISON ===')
print(f'{"Region":12s} | {"Hulk_prist":>12s} {"Sal_prist":>12s} {"Kal_prist":>12s} {"HTS_prist":>12s} | {"Hulk_real":>12s} {"Sal_real":>12s} {"Kal_real":>12s} {"HTS_real":>12s}')
print('-' * 140)

for region in regions:
    prow = prist[(prist['level'] == 'gene') & (prist['region'] == region)]
    rrow = real[(real['level'] == 'gene') & (real['region'] == region)]
    
    get = lambda df, tool: df[df['tool'].str.contains(tool)]['mean_abs_error'].values
    
    hp = get(prow, 'hulkrna')
    sp = get(prow, 'salmon')
    kp = get(prow, 'kallisto')
    htp = get(prow, 'htseq')
    
    hr = get(rrow, 'hulkrna')
    sr = get(rrow, 'salmon')
    kr = get(rrow, 'kallisto')
    htr = get(rrow, 'htseq')
    
    fmt = lambda v: f'{v[0]:.1f}' if len(v) > 0 else 'N/A'
    
    print(f'{region:12s} | {fmt(hp):>12s} {fmt(sp):>12s} {fmt(kp):>12s} {fmt(htp):>12s} | {fmt(hr):>12s} {fmt(sr):>12s} {fmt(kr):>12s} {fmt(htr):>12s}')

# Compute improvement ratios (realistic)
print('\n\n=== HULKRNA ADVANTAGE RATIOS (realistic, transcript-level) ===')
real_tx = real[real['level'] == 'transcript']
for region in regions:
    sub = real_tx[real_tx['region'] == region]
    hulk_mae = sub[sub['tool'].str.contains('hulkrna')]['mean_abs_error'].values
    sal_mae = sub[sub['tool'] == 'salmon']['mean_abs_error'].values
    kal_mae = sub[sub['tool'] == 'kallisto']['mean_abs_error'].values
    
    if len(hulk_mae) > 0 and len(sal_mae) > 0 and len(kal_mae) > 0:
        sal_ratio = sal_mae[0] / max(hulk_mae[0], 0.01)
        kal_ratio = kal_mae[0] / max(hulk_mae[0], 0.01)
        print(f'{region:12s}  hulkrna MAE={hulk_mae[0]:>8.2f}  salmon/hulk={sal_ratio:>5.1f}x  kallisto/hulk={kal_ratio:>5.1f}x')

# Average across regions
print('\n=== AGGREGATE AVERAGES ===')
for label, df in [('pristine', prist), ('realistic', real)]:
    for level in ['transcript', 'gene']:
        sub = df[df['level'] == level]
        tools = sub['tool'].unique()
        print(f'\n{label} {level}-level:')
        for tool in sorted(tools):
            t_sub = sub[sub['tool'] == tool]
            mae = t_sub['mean_abs_error'].mean()
            rmse = t_sub['rmse'].mean()
            pearson = t_sub['pearson'].mean()
            spearman = t_sub['spearman'].mean()
            print(f'  {tool:25s}  MAE={mae:>8.2f}  RMSE={rmse:>8.2f}  Pearson={pearson:.4f}  Spearman={spearman:.4f}')
