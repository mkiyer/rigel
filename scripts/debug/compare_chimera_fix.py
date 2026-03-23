#!/usr/bin/env python3
"""Compare per-transcript quant before and after chimera fix."""
import pyarrow.feather as feather

# Before fix
old = feather.read_table(
    '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file/rigel_minimap2/quant.feather'
).to_pandas()
# After fix
new = feather.read_table(
    '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file/rigel_minimap2_v2/quant.feather'
).to_pandas()
# Oracle
oracle = feather.read_table(
    '/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file/rigel_oracle/quant.feather'
).to_pandas()

gapdh = ['ENST00000229239.10', 'ENST00000396859.5', 'ENST00000396861.8',
         'ENST00000492719.2', 'ENST00000556624.6', 'ENST00000619601.1']

print(f"{'transcript':25s}  {'oracle':>8s}  {'old_mm2':>8s}  {'new_mm2':>8s}  {'delta':>8s}")
print("-" * 70)
for tid in gapdh:
    orow = oracle.loc[oracle['transcript_id'] == tid]
    oldrow = old.loc[old['transcript_id'] == tid]
    newrow = new.loc[new['transcript_id'] == tid]
    o_count = float(orow['mrna'].iloc[0]) if len(orow) else 0
    old_count = float(oldrow['mrna'].iloc[0]) if len(oldrow) else 0
    new_count = float(newrow['mrna'].iloc[0]) if len(newrow) else 0
    diff = new_count - old_count
    print(f"{tid:25s}  {o_count:8.0f}  {old_count:8.0f}  {new_count:8.0f}  {diff:+8.0f}")

print()
o_total = float(oracle['mrna'].sum())
old_total = float(old['mrna'].sum())
new_total = float(new['mrna'].sum())
print(f"Total mrna: oracle={o_total:.0f}  old_mm2={old_total:.0f}  new_mm2={new_total:.0f}")

# Also check pseudogene leakage
print("\n--- Non-GAPDH transcripts with mrna_count > 0 ---")
non_gapdh = new.loc[(~new['transcript_id'].isin(gapdh)) & (new['mrna'] > 10)]
non_gapdh_old = old.loc[(~old['transcript_id'].isin(gapdh)) & (old['mrna'] > 10)]
print(f"Old: {len(non_gapdh_old)} transcripts with >10 counts")
print(f"New: {len(non_gapdh)} transcripts with >10 counts")
for _, row in non_gapdh.nlargest(10, 'mrna').iterrows():
    old_match = old.loc[old['transcript_id'] == row['transcript_id']]
    old_c = float(old_match['mrna'].iloc[0]) if len(old_match) else 0
    print(f"  {row['transcript_id']:30s} g={row['gene_id']:20s} old={old_c:6.0f} new={row['mrna']:6.0f}")
