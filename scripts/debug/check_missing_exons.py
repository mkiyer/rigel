#!/usr/bin/env python3
"""Check if any transcripts are missing from _t_exon_intervals.

If some transcripts have 0 exons in the CSR, they'll get FL=0 in v4,
meaning they're excluded from FL scoring. This could change the EM
solution compared to v3 where they DID get FL values.
"""

import sys
sys.path.insert(0, "/Users/mkiyer/proj/rigel/src")

from rigel.index import TranscriptIndex
import numpy as np

idx_dir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v4/rigel_index"
print(f"Loading index from {idx_dir}...")
idx = TranscriptIndex.load(idx_dir)

n_t = len(idx.t_to_g_arr)
n_with_exons = sum(1 for t in range(n_t) if idx.get_exon_intervals(t) is not None)
n_without = n_t - n_with_exons

print(f"Total transcripts: {n_t}")
print(f"With exon intervals: {n_with_exons}")
print(f"Without exon intervals: {n_without}")

if n_without > 0:
    # Find which transcripts are missing
    missing = [t for t in range(n_t) if idx.get_exon_intervals(t) is None]
    print(f"\nFirst 20 missing transcript indices: {missing[:20]}")
    
    # Get their names from transcript_id array
    for t in missing[:20]:
        tid = idx.transcript_ids[t] if t < len(idx.transcript_ids) else "?"
        gid = idx.gene_ids[idx.t_to_g_arr[t]] if idx.t_to_g_arr[t] >= 0 else "?"
        print(f"  t_idx={t}  {tid}  gene={gid}")
