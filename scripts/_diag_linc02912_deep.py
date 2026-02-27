#!/usr/bin/env python3
"""Deep diagnostic: instrument locus EM for ENST00000624314.1 (t_idx=290)."""
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

# Monkey-patch run_locus_em to instrument the locus containing our target
original_run_locus_em = estimator_mod.AbundanceEstimator.run_locus_em

def instrumented_run_locus_em(self, locus_em, **kwargs):
    t_arr = locus_em.local_to_global_t
    n_t = locus_em.n_transcripts
    
    if TARGET_GLOBAL_T in t_arr:
        local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
        nrna_local = n_t + local_idx
        gdna_idx = 2 * n_t
        
        print(f"\n{'='*60}")
        print(f"INSTRUMENTED LOCUS EM — contains ENST00000624314.1")
        print(f"  n_transcripts = {n_t}")
        print(f"  n_components = {locus_em.n_components}")
        print(f"  n_units = {len(locus_em.offsets) - 1}")
        print(f"  local_idx of ENST00000624314.1 = {local_idx}")
        print(f"  nRNA component index = {nrna_local}")
        print(f"  gdna component index = {gdna_idx}")
        
        # Check prior and unique_totals
        print(f"\n  prior[mRNA ENST000624314] = {locus_em.prior[local_idx]}")
        print(f"  prior[nRNA ENST000624314] = {locus_em.prior[nrna_local]}")
        print(f"  unique_totals[mRNA ENST000624314] = {locus_em.unique_totals[local_idx]}")
        print(f"  unique_totals[nRNA ENST000624314] = {locus_em.unique_totals[nrna_local]}")
        print(f"  nrna_init[local] = {locus_em.nrna_init[local_idx]}")
        
        # Count nRNA candidates for ENST00000624314.1
        n_nrna_cands = int((locus_em.t_indices == nrna_local).sum())
        n_mrna_cands = int((locus_em.t_indices == local_idx).sum())
        print(f"\n  # mRNA candidates for ENST000624314 in CSR = {n_mrna_cands}")
        print(f"  # nRNA candidates for ENST000624314 in CSR = {n_nrna_cands}")
        
        # Show prior distribution
        prior = locus_em.prior
        n_zero_prior = int((prior == 0).sum())
        n_nonzero_prior = int((prior > 0).sum())
        print(f"\n  prior: {n_nonzero_prior} non-zero, {n_zero_prior} zero")
        print(f"  prior[nRNA range] zeros: {int((prior[n_t:2*n_t] == 0).sum())} / {n_t}")
        print(f"  prior[nRNA range] nonzero: {int((prior[n_t:2*n_t] > 0).sum())} / {n_t}")
    
    # Call original
    theta, alpha = original_run_locus_em(self, locus_em, **kwargs)
    
    if TARGET_GLOBAL_T in t_arr:
        print(f"\n  CONVERGED theta:")
        print(f"    theta[mRNA ENST000624314] = {theta[local_idx]:.6e}")
        print(f"    theta[nRNA ENST000624314] = {theta[nrna_local]:.6e}")
        print(f"    theta[gDNA] = {theta[gdna_idx]:.6e}")
        
        # Show top theta values
        top = np.argsort(theta)[::-1][:15]
        print(f"\n  Top 15 theta:")
        from hulkrna.index import HulkIndex
        idx = HulkIndex.load(INDEX_DIR)
        for rank, ci in enumerate(top):
            if ci < n_t:
                pool = "mRNA"
                gt = t_arr[ci]
                tid = idx.t_df.loc[gt, "t_id"]
            elif ci < gdna_idx:
                pool = "nRNA"
                gt = t_arr[ci - n_t]
                tid = idx.t_df.loc[gt, "t_id"]
            else:
                pool = "gDNA"
                tid = "N/A"
                gt = -1
            print(f"    [{rank}] comp={ci} {pool} {tid} theta={theta[ci]:.6e}")
        
        # Show which component would win for fragments that are candidates
        # for nRNA ENST000624314
        print(f"\n  Assignment analysis for nRNA ENST000624314 candidates:")
        offsets = locus_em.offsets
        t_indices = locus_em.t_indices
        log_liks = locus_em.log_liks
        eff_len = locus_em.effective_lengths
        log_weights = np.log(theta + 1e-300) - np.log(eff_len)
        
        n_units = len(offsets) - 1
        nrna_assigned_count = 0
        for u in range(n_units):
            s = int(offsets[u])
            e = int(offsets[u + 1])
            if e <= s:
                continue
            seg_t = t_indices[s:e]
            if nrna_local not in seg_t:
                continue
            # This unit has an nRNA candidate for ENST000624314
            seg_ll = log_liks[s:e]
            seg_lw = log_weights[seg_t]
            seg_lp = seg_ll + seg_lw
            seg_lp -= seg_lp.max()
            seg_p = np.exp(seg_lp)
            seg_p /= seg_p.sum()
            
            best = int(np.argmax(seg_p))
            best_comp = seg_t[best]
            nrna_idx_in_seg = int(np.where(seg_t == nrna_local)[0][0])
            nrna_p = seg_p[nrna_idx_in_seg]
            
            if nrna_p > 0.001:
                nrna_assigned_count += 1
                if nrna_assigned_count <= 5:
                    print(f"    Unit: {len(seg_t)} candidates, nRNA_ENST posterior={nrna_p:.6f}")
                    for j in range(len(seg_t)):
                        ci = seg_t[j]
                        if ci < n_t:
                            pool = "mRNA"
                            gt = t_arr[ci]
                        elif ci < gdna_idx:
                            pool = "nRNA"
                            gt = t_arr[ci - n_t]
                        else:
                            pool = "gDNA"
                            gt = -1
                        tid = idx.t_df.loc[gt, "t_id"] if gt >= 0 else "N/A"
                        print(f"      comp={ci} {pool:4s} {tid:30s} ll={seg_ll[j]:8.2f} lw={seg_lw[j]:10.2f} post={seg_p[j]:.6f}")
        
        print(f"\n  Total units with nRNA ENST000624314 posterior > 0.1%: {nrna_assigned_count}")
        print(f"{'='*60}\n")
    
    return theta, alpha

estimator_mod.AbundanceEstimator.run_locus_em = instrumented_run_locus_em

print("Running pipeline with instrumented locus EM...")
from hulkrna.pipeline import run_pipeline
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=cfg)

print(f"\nFinal: nrna_em_counts[290] = {pipe.estimator.nrna_em_counts[290]:.1f}")
print(f"Final: t_counts[290] sum = {pipe.estimator.t_counts[290].sum():.1f}")
