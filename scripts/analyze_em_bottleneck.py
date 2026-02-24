#!/usr/bin/env python
"""Analyze EM bottleneck and test optimizations."""
import sys
import time
from pathlib import Path
from collections import Counter

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from hulkrna.index import HulkIndex
from hulkrna.pipeline import (
    scan_and_buffer, _scan_and_build_em_data, _build_loci,
    _build_locus_em_data, _compute_eb_gdna_priors, _compute_nrna_init,
    _GDNA_SPLICE_PENALTIES, _SPLICE_UNANNOT,
    _DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT,
    overhang_alpha_to_log_penalty, DEFAULT_OVERHANG_ALPHA, DEFAULT_MISMATCH_ALPHA,
)
from hulkrna.estimator import AbundanceEstimator

bam = '/Users/mkiyer/Downloads/hulkrna_runs/bench_pristine_10_regions/chr17_43044295_43170245/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam'
idx_dir = '/Users/mkiyer/Downloads/hulkrna_runs/bench_pristine_10_regions/chr17_43044295_43170245/hulkrna_index/'

index = HulkIndex.load(idx_dir)
bamfh = pysam.AlignmentFile(bam, 'rb')
stats, sm, flm, buf = scan_and_buffer(
    bamfh.fetch(until_eof=True), index, sj_strand_tag='XS'
)
bamfh.close()
sm.finalize()
flm.finalize()

exl = index.t_df['length'].values.astype(np.float64)
mf = flm.global_model.mean
eff = flm.global_model.compute_all_transcript_eff_lens(exl.astype(np.int64))
gs = (index.g_df['end'].values - index.g_df['start'].values).astype(np.float64)
ts = (index.t_df['end'].values - index.t_df['start'].values).astype(np.float64)
isp = index.intron_span_per_gene()

gdna_sp = dict(_GDNA_SPLICE_PENALTIES)
gdna_sp[_SPLICE_UNANNOT] = _DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT
olp = overhang_alpha_to_log_penalty(DEFAULT_OVERHANG_ALPHA)
mlp = overhang_alpha_to_log_penalty(DEFAULT_MISMATCH_ALPHA)

counter = AbundanceEstimator(
    index.num_transcripts, index.num_genes, alpha=0.5,
    effective_lengths=eff, exonic_lengths=exl, t_to_g=index.t_to_g_arr,
    gene_spans=gs, mean_frag=mf, intronic_spans=isp, transcript_spans=ts,
)

em_data = _scan_and_build_em_data(
    buf, index, sm, flm, counter, stats, 100000,
    gdna_splice_penalties=gdna_sp, overhang_log_penalty=olp,
    mismatch_log_penalty=mlp,
)

nrna_init = _compute_nrna_init(
    counter.transcript_intronic_sense,
    counter.transcript_intronic_antisense,
    ts, exl, mf, sm,
)
counter.nrna_init = nrna_init

loci = _build_loci(em_data, index)
gdna_inits = _compute_eb_gdna_priors(loci, em_data, counter, index, sm)

# Analyze the mega-locus
print("=" * 60)
print("LOCUS ANALYSIS")
print("=" * 60)
for i, l in enumerate(loci):
    n_t = len(l.transcript_indices)
    n_u = len(l.unit_indices)
    offs = em_data.offsets
    n_cand = sum(int(offs[u + 1] - offs[u]) for u in l.unit_indices)
    print(f"Locus {i}: {n_t} transcripts, {n_u} units, {n_cand} candidates, "
          f"cand/unit={n_cand / max(n_u, 1):.1f}")
    patterns = set()
    for u in l.unit_indices:
        start = int(offs[u])
        end = int(offs[u + 1])
        t_set = tuple(sorted(em_data.t_indices[start:end].tolist()))
        patterns.add(t_set)
    print(f"  Unique candidate patterns (equiv classes): {len(patterns)}")

# Now instrument the EM
locus = loci[0]
locus_em = _build_locus_em_data(
    locus, em_data, counter, index, mf, gdna_init=gdna_inits[0],
)

print(f"\nMega-locus EM data:")
print(f"  n_components:      {locus_em.n_components}")
print(f"  n_transcripts:     {locus_em.n_transcripts}")
print(f"  offsets shape:     {locus_em.offsets.shape}")
print(f"  t_indices shape:   {locus_em.t_indices.shape}")
n_units_l = len(locus_em.offsets) - 1
print(f"  n_units:           {n_units_l}")
print(f"  candidates/unit:   {len(locus_em.t_indices) / n_units_l:.1f}")

# EM parameters
_EM_LOG_EPSILON = 1e-300
n_total = locus_em.n_components
offsets_l = locus_em.offsets
t_indices_l = locus_em.t_indices
log_liks_l = locus_em.log_liks
seg_lengths = np.diff(offsets_l)
eff_len = locus_em.effective_lengths
log_eff_len = np.log(eff_len)

# Warm-start theta
theta = locus_em.unique_totals + locus_em.prior
total = theta.sum()
theta = theta / total if total > 0 else theta

print("\n" + "=" * 60)
print("EM ITERATION BENCHMARK (10 iterations each)")
print("=" * 60)

# Method 1: Current (add.at + repeat)
theta1 = theta.copy()
t0 = time.perf_counter()
for _ in range(10):
    log_weights = np.log(theta1 + _EM_LOG_EPSILON) - log_eff_len
    log_posteriors = log_liks_l + log_weights[t_indices_l]
    seg_max = np.maximum.reduceat(log_posteriors, offsets_l[:-1])
    log_posteriors -= np.repeat(seg_max, seg_lengths)
    posteriors = np.exp(log_posteriors)
    seg_sum = np.add.reduceat(posteriors, offsets_l[:-1])
    posteriors /= np.repeat(seg_sum, seg_lengths)
    em_totals = np.zeros(n_total, dtype=np.float64)
    np.add.at(em_totals, t_indices_l, posteriors)
    theta1 = (locus_em.unique_totals + em_totals + locus_em.prior)
    total = theta1.sum()
    theta1 /= total
t1 = time.perf_counter()
print(f"  add.at + repeat:       {t1 - t0:.3f}s  ({(t1 - t0) / 10 * 1000:.1f} ms/iter)")

# Method 2: bincount + precomputed unit_idx
theta2 = theta.copy()
unit_idx = np.repeat(np.arange(n_units_l, dtype=np.intp), seg_lengths)
t0 = time.perf_counter()
for _ in range(10):
    log_weights = np.log(theta2 + _EM_LOG_EPSILON) - log_eff_len
    log_posteriors = log_liks_l + log_weights[t_indices_l]
    seg_max = np.maximum.reduceat(log_posteriors, offsets_l[:-1])
    log_posteriors -= seg_max[unit_idx]
    posteriors = np.exp(log_posteriors)
    seg_sum = np.add.reduceat(posteriors, offsets_l[:-1])
    posteriors /= seg_sum[unit_idx]
    em_totals = np.bincount(
        t_indices_l, weights=posteriors, minlength=n_total,
    ).astype(np.float64)
    theta2 = (locus_em.unique_totals + em_totals + locus_em.prior)
    total = theta2.sum()
    theta2 /= total
t1 = time.perf_counter()
print(f"  bincount + unit_idx:   {t1 - t0:.3f}s  ({(t1 - t0) / 10 * 1000:.1f} ms/iter)")
print(f"  theta match: {np.allclose(theta1, theta2, atol=1e-8)}")

# Convergence analysis: run full EM and track delta per iteration
print("\n" + "=" * 60)
print("EM CONVERGENCE PROFILE")
print("=" * 60)

theta_c = theta.copy()
deltas = []
for it in range(1000):
    log_weights = np.log(theta_c + _EM_LOG_EPSILON) - log_eff_len
    log_posteriors = log_liks_l + log_weights[t_indices_l]
    seg_max = np.maximum.reduceat(log_posteriors, offsets_l[:-1])
    log_posteriors -= seg_max[unit_idx]
    posteriors = np.exp(log_posteriors)
    seg_sum = np.add.reduceat(posteriors, offsets_l[:-1])
    posteriors /= seg_sum[unit_idx]
    em_totals = np.bincount(
        t_indices_l, weights=posteriors, minlength=n_total,
    ).astype(np.float64)
    theta_new = (locus_em.unique_totals + em_totals + locus_em.prior)
    total = theta_new.sum()
    theta_new /= total
    delta = np.abs(theta_new - theta_c).sum()
    deltas.append(delta)
    theta_c = theta_new
    if delta < 1e-6:
        break

print(f"  Converged after {len(deltas)} iterations")
print(f"  Final delta: {deltas[-1]:.2e}")
for milestone in [10, 50, 100, 200, 500, len(deltas) - 1]:
    if milestone < len(deltas):
        print(f"    iter {milestone + 1:>4d}: delta = {deltas[milestone]:.2e}")

# Equivalence class analysis
print("\n" + "=" * 60)
print("EQUIVALENCE CLASS ANALYSIS (local EM indices)")
print("=" * 60)

pattern_counts = Counter()
for u in range(n_units_l):
    start = int(offsets_l[u])
    end = int(offsets_l[u + 1])
    pattern = tuple(sorted(t_indices_l[start:end].tolist()))
    pattern_counts[pattern] += 1

print(f"Total EM units:     {n_units_l:,}")
print(f"Unique patterns:    {len(pattern_counts):,}")
print(f"Compression ratio:  {n_units_l / len(pattern_counts):.1f}x")
print(f"\nTop 15 patterns by frequency:")
for rank, (p, c) in enumerate(pattern_counts.most_common(15)):
    print(f"  #{rank + 1:>2d}: {c:>6,} units, {len(p):>3d} candidates  "
          f"indices={list(p)[:8]}{'...' if len(p) > 8 else ''}")

buf.cleanup()
