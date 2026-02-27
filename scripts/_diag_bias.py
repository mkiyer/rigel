#!/usr/bin/env python3
"""Trace why some mRNA candidates have ll=-703 after bias correction."""
import sys
import math
import numpy as np
from pathlib import Path

sys.path.insert(0, "src")

from hulkrna.index import HulkIndex
from hulkrna.config import PipelineConfig
import hulkrna.estimator as estimator_mod

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")

TARGET_GLOBAL_T = 290

captured = {}

# Capture BEFORE bias correction
original_run_locus_em = estimator_mod.AbundanceEstimator.run_locus_em

def capture_run_locus_em(self, locus_em, **kwargs):
    t_arr = locus_em.local_to_global_t
    if TARGET_GLOBAL_T in t_arr:
        n_t = locus_em.n_transcripts
        local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
        nrna_local = n_t + local_idx  # 311
        mrna_local = local_idx  # 118
        
        # Save log_liks BEFORE bias correction
        captured['log_liks_before'] = locus_em.log_liks.copy()
        captured['t_indices'] = locus_em.t_indices.copy()
        captured['tx_starts'] = locus_em.tx_starts.copy()
        captured['tx_ends'] = locus_em.tx_ends.copy()
        captured['offsets'] = locus_em.offsets.copy()
        captured['bias_profiles'] = locus_em.bias_profiles
        captured['nrna_local'] = nrna_local
        captured['mrna_local'] = mrna_local
        captured['n_t'] = n_t
        
    result = original_run_locus_em(self, locus_em, **kwargs)
    
    if TARGET_GLOBAL_T in t_arr:
        captured['log_liks_after'] = locus_em.log_liks.copy()
    
    return result

estimator_mod.AbundanceEstimator.run_locus_em = capture_run_locus_em

print("Running pipeline...")
from hulkrna.pipeline import run_pipeline
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=cfg)

# Analyze the captured data
ll_before = captured['log_liks_before']
ll_after = captured['log_liks_after']
t_indices = captured['t_indices']
tx_starts = captured['tx_starts']
tx_ends = captured['tx_ends']
offsets = captured['offsets']
profiles = captured['bias_profiles']
nrna_local = captured['nrna_local']
mrna_local = captured['mrna_local']
n_t = captured['n_t']

n_units = len(offsets) - 1

print(f"\nnrna_local={nrna_local}, mrna_local={mrna_local}")

# Find units where comp mrna_local has very negative ll_after
bad_mrna_count = 0
good_mrna_count = 0

for u in range(n_units):
    s = int(offsets[u])
    e = int(offsets[u + 1])
    if s == e:
        continue
    seg_t = t_indices[s:e]
    seg_ll_before = ll_before[s:e]
    seg_ll_after = ll_after[s:e]
    seg_tx_s = tx_starts[s:e]
    seg_tx_e = tx_ends[s:e]
    
    # Check if mrna_local and nrna_local are in this unit
    mrna_pos = -1
    nrna_pos = -1
    for j in range(len(seg_t)):
        if seg_t[j] == mrna_local:
            mrna_pos = j
        if seg_t[j] == nrna_local:
            nrna_pos = j
    
    if mrna_pos < 0 or nrna_pos < 0:
        continue
    
    mrna_ll_b = seg_ll_before[mrna_pos]
    mrna_ll_a = seg_ll_after[mrna_pos]
    nrna_ll_b = seg_ll_before[nrna_pos]
    nrna_ll_a = seg_ll_after[nrna_pos]
    
    if mrna_ll_a < -100:
        bad_mrna_count += 1
        if bad_mrna_count <= 5:
            # Compute bias correction manually
            mrna_fl = int(seg_tx_e[mrna_pos]) - int(seg_tx_s[mrna_pos])
            nrna_fl = int(seg_tx_e[nrna_pos]) - int(seg_tx_s[nrna_pos])
            
            mrna_profile = profiles[mrna_local]
            nrna_profile = profiles[nrna_local]
            
            mrna_eff = mrna_profile.effective_length(mrna_fl)
            nrna_eff = nrna_profile.effective_length(nrna_fl)
            
            mrna_log_eff = math.log(max(mrna_eff, 1e-300))
            nrna_log_eff = math.log(max(nrna_eff, 1e-300))
            
            print(f"\n  BAD unit #{bad_mrna_count} (unit {u}):")
            print(f"    mRNA: tx_start={seg_tx_s[mrna_pos]}, tx_end={seg_tx_e[mrna_pos]}, "
                  f"frag_len={mrna_fl}")
            print(f"    mRNA: profile_length={mrna_profile.length}, "
                  f"eff_len={mrna_eff}, log_eff={mrna_log_eff:.4f}")
            print(f"    mRNA: ll_before={mrna_ll_b:.4f}, correction=-{mrna_log_eff:.4f}, "
                  f"ll_after={mrna_ll_a:.4f}")
            print(f"    nRNA: tx_start={seg_tx_s[nrna_pos]}, tx_end={seg_tx_e[nrna_pos]}, "
                  f"frag_len={nrna_fl}")
            print(f"    nRNA: profile_length={nrna_profile.length}, "
                  f"eff_len={nrna_eff}, log_eff={nrna_log_eff:.4f}")
            print(f"    nRNA: ll_before={nrna_ll_b:.4f}, correction=-{nrna_log_eff:.4f}, "
                  f"ll_after={nrna_ll_a:.4f}")
    else:
        good_mrna_count += 1

print(f"\n\nUnits with GOOD mRNA ll (> -100): {good_mrna_count}")
print(f"Units with BAD mRNA ll (< -100): {bad_mrna_count}")
print(f"Total units with both comps: {good_mrna_count + bad_mrna_count}")
