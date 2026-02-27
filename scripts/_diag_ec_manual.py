#!/usr/bin/env python3
"""Directly compute what em_totals[311] should be in the first E-step."""
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, "src")

from hulkrna.index import HulkIndex
from hulkrna.config import PipelineConfig
import hulkrna.estimator as estimator_mod
import hulkrna.locus as locus_mod

INDEX_DIR = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/hulkrna_index")
BAM_PATH = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam")

TARGET_GLOBAL_T = 290

# Capture the locus_em and the estimator
captured = {}

original_run_locus_em = estimator_mod.AbundanceEstimator.run_locus_em

def capture_run_locus_em(self, locus_em, **kwargs):
    t_arr = locus_em.local_to_global_t
    if TARGET_GLOBAL_T in t_arr:
        captured['locus_em'] = locus_em
        captured['estimator'] = self
    return original_run_locus_em(self, locus_em, **kwargs)

estimator_mod.AbundanceEstimator.run_locus_em = capture_run_locus_em

print("Running pipeline to capture locus EM data...")
from hulkrna.pipeline import run_pipeline
cfg = PipelineConfig()
pipe = run_pipeline(BAM_PATH, HulkIndex.load(INDEX_DIR), config=cfg)

locus_em = captured['locus_em']
est = captured['estimator']

t_arr = locus_em.local_to_global_t
n_t = locus_em.n_transcripts
local_idx = int(np.where(t_arr == TARGET_GLOBAL_T)[0][0])
nrna_local = n_t + local_idx
mrna_local = local_idx
gdna_idx = 2 * n_t
n_total = locus_em.n_components

print(f"n_t={n_t}, local_idx={local_idx}, nrna={nrna_local}, mrna={mrna_local}")

# Reproduce the warm start to get theta
prior_orig = locus_em.prior.copy()
unique_totals = locus_em.unique_totals.copy()

# Apply bias correction (this mutates log_liks in place — already done)
# The captured locus_em already has bias-corrected log_liks

# Build equiv classes
ec_list = estimator_mod._build_equiv_classes(
    locus_em.offsets, locus_em.t_indices, locus_em.log_liks,
    locus_em.coverage_weights,
)
ec_data = []
ec_wt = []
for comp_idx, ll_matrix, wt_matrix in ec_list:
    scratch = np.empty_like(ll_matrix)
    ec_data.append((comp_idx, ll_matrix, scratch))
    ec_wt.append(wt_matrix)

# Warm start
eligible = prior_orig > 0.0
theta_init = unique_totals.copy()
coverage_totals = np.zeros(n_total, dtype=np.float64)
n_ambiguous = 0

for i, (comp_idx, ll_matrix, _scratch) in enumerate(ec_data):
    wt_matrix = ec_wt[i]
    n = ll_matrix.shape[0]
    n_ambiguous += n
    has_p = eligible[comp_idx]
    masked_wt = wt_matrix * has_p[np.newaxis, :]
    row_sums = masked_wt.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    shares = masked_wt / row_sums
    per_comp = shares.sum(axis=0)
    theta_init[comp_idx] += per_comp
    coverage_totals[comp_idx] += per_comp

alpha = est.em_config.prior_alpha
gamma = est.em_config.prior_gamma
if n_ambiguous > 0 and gamma > 0.0:
    ovr = gamma * coverage_totals / n_ambiguous
else:
    ovr = np.zeros(n_total, dtype=np.float64)
prior = np.where(eligible, alpha + ovr, 0.0)

theta = theta_init + prior
total = theta.sum()
if total > 0:
    theta /= total

print(f"\ntheta[mRNA]={theta[mrna_local]:.10e}")
print(f"theta[nRNA]={theta[nrna_local]:.10e}")
print(f"theta[nRNA] is zero? {theta[nrna_local] == 0.0}")
print(f"theta[nRNA] bits: {theta[nrna_local].hex()}")

# Now compute what the first E-step would give
_EM_LOG_EPSILON = 1e-300
eff_len = locus_em.effective_lengths
log_eff_len = np.log(eff_len)
log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len

print(f"\nlog_weight[mRNA]={log_weights[mrna_local]:.6f}")
print(f"log_weight[nRNA]={log_weights[nrna_local]:.6f}")

# Manually compute E-step for comp 311
em_totals_manual = np.zeros(n_total, dtype=np.float64)

total_nrna_posterior = 0.0
n_ecs_with_nrna = 0

for i, (comp_idx, ll_matrix, scratch) in enumerate(ec_data):
    col_list = comp_idx.tolist()
    if nrna_local not in col_list:
        continue
    n_ecs_with_nrna += 1
    nrna_col = col_list.index(nrna_local)
    mrna_col = col_list.index(mrna_local) if mrna_local in col_list else -1
    
    # Manual E-step for this EC
    lw = log_weights[comp_idx]
    s = ll_matrix + lw[np.newaxis, :]  # (n, k)
    max_r = s.max(axis=1, keepdims=True)
    s = s - max_r
    s = np.exp(s)
    row_sum = s.sum(axis=1, keepdims=True)
    bad = (row_sum == 0) | ~np.isfinite(row_sum)
    if bad.any():
        row_sum[bad] = 1.0
        s_normed = s / row_sum
        s_normed[bad.ravel()] = 0.0
    else:
        s_normed = s / row_sum
    
    nrna_posterior_sum = float(s_normed[:, nrna_col].sum())
    total_nrna_posterior += nrna_posterior_sum
    
    if n_ecs_with_nrna <= 3:
        n_rows = ll_matrix.shape[0]
        print(f"\n  EC #{n_ecs_with_nrna}: {n_rows} rows x {len(comp_idx)} cols")
        print(f"    nrna_col={nrna_col}, mrna_col={mrna_col}")
        print(f"    lw[nrna_col]={lw[nrna_col]:.4f}")
        if mrna_col >= 0:
            print(f"    lw[mrna_col]={lw[mrna_col]:.4f}")
        # Check first row
        row0_ll = ll_matrix[0]
        row0_scratch = row0_ll + lw
        row0_max = row0_scratch.max()
        row0_after = np.exp(row0_scratch - row0_max)
        row0_sum = row0_after.sum()
        row0_post = row0_after / row0_sum if row0_sum > 0 else row0_after
        mrna_ll = row0_ll[mrna_col] if mrna_col >= 0 else float('nan')
        mrna_sc = row0_scratch[mrna_col] if mrna_col >= 0 else float('nan')
        print(f"    Row 0: ll[nrna]={row0_ll[nrna_col]:.4f}, ll[mrna]={mrna_ll:.4f}")
        print(f"    Row 0: scratch[nrna]={row0_scratch[nrna_col]:.4f}, scratch[mrna]={mrna_sc:.4f}")
        print(f"    Row 0: max_r={row0_max:.4f}")
        mrna_post = row0_post[mrna_col] if mrna_col >= 0 else float('nan')
        print(f"    Row 0: post[nrna]={row0_post[nrna_col]:.10e}")
        if mrna_col >= 0:
            print(f"    Row 0: post[mrna]={mrna_post:.10e}")
        
        # Stats
        nrna_posteriors = s_normed[:, nrna_col]
        print(f"    nrna posterior: max={nrna_posteriors.max():.10e}, "
              f"mean={nrna_posteriors.mean():.10e}, sum={nrna_posterior_sum:.10e}")
        # Top 5 posterior values
        top5 = np.argsort(nrna_posteriors)[::-1][:5]
        for j in top5:
            mrna_ll_j = ll_matrix[j, mrna_col] if mrna_col >= 0 else float('nan')
            print(f"      row {j}: nrna_post={nrna_posteriors[j]:.10e}, "
                  f"ll_nrna={ll_matrix[j, nrna_col]:.4f}, "
                  f"ll_mrna={mrna_ll_j:.4f}")

print(f"\n\nTotal ECs with nRNA comp {nrna_local}: {n_ecs_with_nrna}")
print(f"Total manual em_totals[nRNA]: {total_nrna_posterior:.10e}")
print(f"Expected em_totals from actual _em_step: ~1.826")
