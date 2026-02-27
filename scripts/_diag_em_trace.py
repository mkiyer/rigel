#!/usr/bin/env python3
"""Deep EM trace: dump alpha_dir / theta per iteration for dead nRNA comp."""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, "src")

from hulkrna.index import HulkIndex
from hulkrna.config import PipelineConfig
import hulkrna.estimator as estimator_mod

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")

TARGET_GLOBAL_T = 290  # ENST00000624314.1

# Monkey-patch _em_step to trace
original_em_step = estimator_mod._em_step
em_step_call_count = [0]
trace_comp = [None]  # (nrna_local, mrna_local) when set

def traced_em_step(theta, ec_data, log_eff_len, unique_totals, prior, em_totals):
    result = original_em_step(theta, ec_data, log_eff_len, unique_totals, prior, em_totals)
    if trace_comp[0] is not None:
        nrna_i, mrna_i = trace_comp[0]
        em_step_call_count[0] += 1
        n = em_step_call_count[0]
        if n <= 30 or n % 10 == 0:
            print(
                f"  _em_step #{n:3d}: "
                f"theta_in[nRNA]={theta[nrna_i]:.6e}  "
                f"em_totals[nRNA]={em_totals[nrna_i]:.6e}  "
                f"prior[nRNA]={prior[nrna_i]:.6e}  "
                f"unique[nRNA]={unique_totals[nrna_i]:.6e}  "
                f"theta_out[nRNA]={result[nrna_i]:.6e}  ||  "
                f"theta_in[mRNA]={theta[mrna_i]:.6e}  "
                f"em_totals[mRNA]={em_totals[mrna_i]:.6e}  "
                f"theta_out[mRNA]={result[mrna_i]:.6e}"
            )
    return result

estimator_mod._em_step = traced_em_step

# Monkey-patch run_locus_em to enable tracing for our target locus
original_run_locus_em = estimator_mod.AbundanceEstimator.run_locus_em

def instrumented_run_locus_em(self, locus_em, **kwargs):
    t_arr = locus_em.local_to_global_t
    n_t = locus_em.n_transcripts
    
    if TARGET_GLOBAL_T in t_arr:
        local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
        nrna_local = n_t + local_idx
        
        print(f"\n{'='*70}")
        print(f"LOCUS EM for ENST00000624314.1 — n_t={n_t}, local={local_idx}, nRNA={nrna_local}")
        print(f"  prior[nRNA]={locus_em.prior[nrna_local]:.6e}")
        print(f"  prior[mRNA]={locus_em.prior[local_idx]:.6e}")
        print(f"  unique_totals[nRNA]={locus_em.unique_totals[nrna_local]:.6e}")
        print(f"  unique_totals[mRNA]={locus_em.unique_totals[local_idx]:.6e}")
        print(f"  em_config.mode={self.em_config.mode}")
        print(f"  em_config.prune_threshold={self.em_config.prune_threshold}")
        print(f"  em_config.prior_alpha={self.em_config.prior_alpha}")
        print(f"  em_config.prior_gamma={self.em_config.prior_gamma}")
        
        # Count candidates
        n_nrna_cands = int((locus_em.t_indices == nrna_local).sum())
        n_mrna_cands = int((locus_em.t_indices == local_idx).sum())
        print(f"  # nRNA candidates in CSR={n_nrna_cands}")
        print(f"  # mRNA candidates in CSR={n_mrna_cands}")
        
        # Enable tracing
        trace_comp[0] = (nrna_local, local_idx)
        em_step_call_count[0] = 0
    
    theta, alpha = original_run_locus_em(self, locus_em, **kwargs)
    
    if TARGET_GLOBAL_T in t_arr:
        local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
        nrna_local = n_t + local_idx
        gdna_idx = 2 * n_t
        
        print(f"\n  FINAL theta[nRNA]={theta[nrna_local]:.6e}")
        print(f"  FINAL theta[mRNA]={theta[local_idx]:.6e}")
        print(f"  FINAL theta[gDNA]={theta[gdna_idx]:.6e}")
        
        # Check how many nRNA components have non-zero theta
        nrna_thetas = theta[n_t:2*n_t]
        print(f"  nRNA comps with theta > 1e-10: {int((nrna_thetas > 1e-10).sum())}")
        print(f"  nRNA theta sum: {nrna_thetas.sum():.6e}")
        top_nrna = np.argsort(nrna_thetas)[::-1][:5]
        for i, k in enumerate(top_nrna):
            print(f"    top nRNA [{i}]: local={k}, theta={nrna_thetas[k]:.6e}, "
                  f"global_t={t_arr[k]}")
        
        trace_comp[0] = None
        print(f"{'='*70}\n")
    
    return theta, alpha

estimator_mod.AbundanceEstimator.run_locus_em = instrumented_run_locus_em

print("Running pipeline with EM iteration tracing...")
from hulkrna.pipeline import run_pipeline
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=cfg)

print(f"\nFinal: nrna_em_counts[290] = {pipe.estimator.nrna_em_counts[290]:.1f}")
