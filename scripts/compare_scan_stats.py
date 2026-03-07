#!/usr/bin/env python3
"""Compare scan-phase stats across multiple rigel output directories."""
import json
import sys

dirs = [
    '/Users/mkiyer/Downloads/rigel_runs/baseline_output',
    '/Users/mkiyer/Downloads/rigel_runs/phase1_output',
    '/Users/mkiyer/Downloads/rigel_runs/phase1_output_run2',
]
labels = ['baseline', 'phase1_r1', 'phase1_r2']

summaries = []
for d in dirs:
    with open(f'{d}/summary.json') as f:
        summaries.append(json.load(f))

# Scan-phase stats: these are deterministic and grouping-sensitive
scan_keys = [
    'total_reads', 'mapped_reads', 'unique_alignments',
    'n_fragments', 'n_proper_pairs', 'n_chimeric',
    'n_intergenic', 'n_intergenic_spliced', 'n_intergenic_unspliced',
    'n_with_exon', 'n_with_annotated_sj',
    'n_strand_trained', 'n_same_strand', 'n_ambig_strand',
    'n_frag_length_unambiguous', 'n_frag_length_ambiguous',
    'n_frag_length_intergenic',
]

print('Scan-phase stats (MUST be exact for grouping correctness):')
header = f'{"Key":<42s}'
for l in labels:
    header += f' {l:>14s}'
header += '  Match?'
print(header)
print('-' * 100)

all_scan_ok = True
for key in scan_keys:
    vals = []
    for s in summaries:
        v = s.get('alignment', {}).get(key,
            s.get('quantification', {}).get(key, 'N/A'))
        vals.append(v)
    match = vals[0] == vals[1] == vals[2]
    tag = 'YES' if match else 'NO <<<'
    if not match:
        all_scan_ok = False
    row = f'{key:<42s}'
    for v in vals:
        row += f' {str(v):>14s}'
    row += f'  {tag}'
    print(row)

# EM-phase stats: may vary due to threading non-determinism
em_keys = ['n_gdna_em', 'n_gdna_total', 'n_gdna_unambig']
print()
print('EM-phase stats (may vary due to threading):')
for key in em_keys:
    vals = []
    for s in summaries:
        v = s.get('quantification', {}).get(key, 'N/A')
        vals.append(v)
    match = vals[0] == vals[1] == vals[2]
    tag = 'YES' if match else 'varies'
    row = f'{key:<42s}'
    for v in vals:
        row += f' {str(v):>14s}'
    row += f'  {tag}'
    print(row)

# Library model stats (grouping-sensitive, deterministic)
print()
print('Library model stats (deterministic):')
for key in ['strand_specificity', 'p_r1_sense', 'frag_length_mean',
            'frag_length_median', 'frag_length_std', 'frag_length_mode']:
    vals = []
    for s in summaries:
        v = s.get('library', {}).get(key, 'N/A')
        vals.append(v)
    match = vals[0] == vals[1] == vals[2]
    tag = 'YES' if match else 'NO <<<'
    if not match:
        all_scan_ok = False
    row = f'{key:<42s}'
    for v in vals:
        row += f' {str(v):>14s}'
    row += f'  {tag}'
    print(row)

print()
if all_scan_ok:
    print('*** ALL SCAN-PHASE AND LIBRARY STATS MATCH ***')
    print('Qname grouping is preserved.')
else:
    print('*** SCAN-PHASE STATS MISMATCH — GROUPING MAY BE BROKEN ***')
    sys.exit(1)
