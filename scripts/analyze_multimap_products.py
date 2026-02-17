"""Analyze R1×R2 cross-product sizes in real benchmark BAMs."""
import glob
import sys
from collections import defaultdict
import pysam

bam_files = sorted(glob.glob('/Users/mkiyer/Downloads/hulkrna_runs/bench_zero_gdna_v2/*/align/reads.bam'))

for bam_path in bam_files:
    region = bam_path.split('/')[-3]
    bam = pysam.AlignmentFile(bam_path, 'rb')
    qnames = defaultdict(lambda: {'r1': 0, 'r2': 0, 'r1_refs': set(), 'r2_refs': set()})
    total = 0
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_supplementary:
            continue
        total += 1
        key = 'r1' if r.is_read1 else 'r2'
        qnames[r.query_name][key] += 1
        ref_key = 'r1_refs' if r.is_read1 else 'r2_refs'
        qnames[r.query_name][ref_key].add(r.reference_id)

    n_multimap = 0
    max_product = 0
    product_dist = defaultdict(int)
    asymmetric = 0
    for q, c in qnames.items():
        if c['r1'] > 1 or c['r2'] > 1:
            n_multimap += 1
            product = c['r1'] * c['r2']
            max_product = max(max_product, product)
            product_dist[product] += 1
            if c['r1'] != c['r2']:
                asymmetric += 1

    bam.close()
    print(f'\n{region}: {total} rec, {len(qnames)} reads, {n_multimap} multimap ({100*n_multimap/max(len(qnames),1):.1f}%), max_product={max_product}, asymm={asymmetric}')
    for p in sorted(product_dist)[:10]:
        cnt = product_dist[p]
        print(f'  product={p}: {cnt} reads ({100*cnt/max(n_multimap,1):.1f}%)')
    if len(product_dist) > 10:
        remaining = sum(v for k, v in product_dist.items() if k > sorted(product_dist.keys())[9])
        print(f'  ... {remaining} more with product>{sorted(product_dist.keys())[9]}')
