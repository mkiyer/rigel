#!/usr/bin/env python3
"""Compare fragment stats and per-transcript quant across all 3 versions:
  v1 (original): per-hit chimera counting + CROSS-PAIR
  v2 (no cross-pair): per-fragment chimera + singletons instead of cross-pair
  v3 (restored cross-pair): per-fragment chimera + CROSS-PAIR restored
"""
import json
import pyarrow.feather as feather

BASE = '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file'

versions = {
    'v1_original': f'{BASE}/rigel_minimap2',
    'v2_no_xpair': f'{BASE}/rigel_minimap2_v2',
    'v3_xpair_ok': f'{BASE}/rigel_minimap2_v3',
    'oracle': f'{BASE}/rigel_oracle',
}

# Fragment stats
print("=" * 80)
print("FRAGMENT STATS COMPARISON")
print("=" * 80)
metrics = ['total', 'genic', 'intergenic', 'chimeric', 'chimeric_trans',
           'chimeric_cis_same', 'chimeric_cis_diff']
print(f"{'metric':25s}", end="")
for name in versions:
    print(f"  {name:>14s}", end="")
print()
print("-" * 90)

stats_data = {}
for name, path in versions.items():
    with open(f'{path}/summary.json') as f:
        stats_data[name] = json.load(f)

for m in metrics:
    print(f"{m:25s}", end="")
    for name in versions:
        val = stats_data[name]['fragment_stats'].get(m, 'N/A')
        print(f"  {val:>14s}" if isinstance(val, str) else f"  {val:>14,d}", end="")
    print()

# Also show n_multimapper_alignments
print()
print("Multimapper alignments entering EM:")
for name in versions:
    q = stats_data[name]['quantification']
    print(f"  {name}: n_em_assigned={q['n_em_assigned']}, mrna_total={q['mrna_total']:.0f}")

# Per-transcript comparison
print()
print("=" * 80)
print("PER-TRANSCRIPT QUANT (GAPDH basic isoforms)")
print("=" * 80)

gapdh = ['ENST00000229239.10', 'ENST00000396859.5', 'ENST00000396861.8',
         'ENST00000492719.2', 'ENST00000556624.6', 'ENST00000619601.1']

quant = {}
for name, path in versions.items():
    df = feather.read_table(f'{path}/quant.feather').to_pandas()
    quant[name] = df

print(f"{'transcript':25s} {'oracle':>8s} {'v1_orig':>8s} {'v2_noXP':>8s} {'v3_XP':>8s} {'v3/oracle':>10s}")
print("-" * 75)
for tid in gapdh:
    vals = {}
    for name in versions:
        row = quant[name].loc[quant[name]['transcript_id'] == tid]
        vals[name] = float(row['mrna'].iloc[0]) if len(row) else 0
    ratio = vals['v3_xpair_ok'] / vals['oracle'] if vals['oracle'] > 0 else float('inf')
    print(f"{tid:25s} {vals['oracle']:8.0f} {vals['v1_original']:8.0f} "
          f"{vals['v2_no_xpair']:8.0f} {vals['v3_xpair_ok']:8.0f} {ratio:10.2f}x")

# Pseudogene leakage
print()
print("NON-GAPDH LEAKAGE (transcripts with >10 mrna counts):")
for name in ['v1_original', 'v2_no_xpair', 'v3_xpair_ok']:
    df = quant[name]
    leak = df.loc[(~df['transcript_id'].isin(gapdh)) & (df['mrna'] > 10)]
    total_leak = float(leak['mrna'].sum())
    print(f"  {name}: {len(leak)} transcripts, {total_leak:.0f} total leaked counts")
    for _, row in leak.nlargest(5, 'mrna').iterrows():
        print(f"    {row['transcript_id']:30s} gene={row['gene_name']:15s} mrna={row['mrna']:6.0f}")
