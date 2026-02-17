#!/usr/bin/env python3
from pathlib import Path
import json
import statistics as st

root = Path('/Users/mkiyer/Downloads/hulkrna_runs/fgfr2_depth_runs_postfix')
rows = []
for seed_dir in sorted(root.glob('seed_*')):
    sj = seed_dir / 'summary.json'
    if not sj.exists():
        continue
    data = json.loads(sj.read_text())
    if not data:
        continue
    rec = data[0]
    seed = int(seed_dir.name.split('_')[-1])
    for tool in ('hulkrna','salmon','kallisto'):
        t = rec[tool]
        rows.append({
            'seed': seed,
            'tool': tool,
            'mae': float(t['mean_abs_error']),
            'pearson': float(t['pearson']),
            'spearman': float(t['spearman']),
        })

if not rows:
    print('no data')
    raise SystemExit(1)

for tool in ('hulkrna','salmon','kallisto'):
    vals = [r for r in rows if r['tool']==tool]
    maes = [r['mae'] for r in vals]
    pe = [r['pearson'] for r in vals]
    sp = [r['spearman'] for r in vals]
    print(f"{tool}: n={len(vals)} mean_mae={st.mean(maes):.3f} median_mae={st.median(maes):.3f} mean_pearson={st.mean(pe):.5f} mean_spearman={st.mean(sp):.5f}")

print('\nPer-seed hulkrna MAE:')
for r in sorted([x for x in rows if x['tool']=='hulkrna'], key=lambda z: z['seed']):
    print(f"seed={r['seed']} mae={r['mae']:.3f}")
