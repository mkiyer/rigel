#!/usr/bin/env python
"""Benchmark different EM iteration strategies."""
import sys
import time
from pathlib import Path

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
from hulkrna.estimator import AbundanceEstimator, _build_equiv_classes

_EM_LOG_EPSILON = 1e-300
_EM_CONVERGENCE_DELTA = 1e-6


def setup(bam, idx_dir):
    """Load data and build locus EM input."""
    index = HulkIndex.load(idx_dir)
    bamfh = pysam.AlignmentFile(bam, 'rb')
    stats, sm, flm, buf = scan_and_buffer(
        bamfh.fetch(until_eof=True), index, sj_strand_tag='XS',
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
        counter.transcript_intronic_antisense, ts, exl, mf, sm,
    )
    counter.nrna_init = nrna_init
    loci = _build_loci(em_data, index)
    gdna_inits = _compute_eb_gdna_priors(loci, em_data, counter, index, sm)
    locus_em = _build_locus_em_data(
        loci[0], em_data, counter, index, mf, gdna_init=gdna_inits[0],
    )
    buf.cleanup()
    return locus_em


def common_init(locus_em):
    """Warm-start theta (original CSR approach)."""
    n_total = locus_em.n_components
    prior = locus_em.prior.copy()
    unique_totals = locus_em.unique_totals.copy()
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks
    n_units = len(offsets) - 1
    seg_lengths = np.diff(offsets)
    eff_len = locus_em.effective_lengths
    log_eff_len = np.log(eff_len)

    theta_init = unique_totals.copy()
    if len(seg_lengths) > 0 and len(t_indices) > 0:
        cand_weights = 1.0 / np.repeat(
            seg_lengths.astype(np.float64), seg_lengths,
        )
        has_prior = prior[t_indices] > 0.0
        cand_weights *= has_prior
        np.add.at(theta_init, t_indices, cand_weights)

    theta = theta_init + prior
    total = theta.sum()
    if total > 0:
        theta /= total
    return theta, prior, unique_totals, log_eff_len


def approach_a_original(locus_em, n_iters):
    """Original CSR: repeat + add.at."""
    theta, prior, unique_totals, log_eff_len = common_init(locus_em)
    n_total = locus_em.n_components
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks
    seg_lengths = np.diff(offsets)

    t0 = time.perf_counter()
    em_totals = np.zeros(n_total, dtype=np.float64)
    for _ in range(n_iters):
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        log_posteriors = log_liks + log_weights[t_indices]
        seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
        log_posteriors -= np.repeat(seg_max, seg_lengths)
        posteriors = np.exp(log_posteriors)
        seg_sum = np.add.reduceat(posteriors, offsets[:-1])
        bad_seg = (seg_sum == 0) | ~np.isfinite(seg_sum)
        seg_sum[bad_seg] = 1.0
        posteriors /= np.repeat(seg_sum, seg_lengths)
        if bad_seg.any():
            bad_mask = np.repeat(bad_seg, seg_lengths)
            posteriors[bad_mask] = 0.0
        em_totals = np.zeros(n_total, dtype=np.float64)
        np.add.at(em_totals, t_indices, posteriors)
        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total
        delta = np.abs(theta_new - theta).sum()
        theta = theta_new
        if delta < _EM_CONVERGENCE_DELTA:
            break
    elapsed = time.perf_counter() - t0
    return theta, unique_totals + em_totals + prior, elapsed


def approach_b_equiv_class(locus_em, n_iters):
    """Equivalence-class: per-class dense matrix."""
    theta, prior, unique_totals, log_eff_len = common_init(locus_em)
    n_total = locus_em.n_components
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks

    ec_list = _build_equiv_classes(offsets, t_indices, log_liks)

    t0 = time.perf_counter()
    em_totals = np.zeros(n_total, dtype=np.float64)
    for _ in range(n_iters):
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        em_totals[:] = 0.0
        for comp_idx, ll_matrix in ec_list:
            lw = log_weights[comp_idx]
            log_post = ll_matrix + lw
            max_r = log_post.max(axis=1, keepdims=True)
            post = np.exp(log_post - max_r)
            row_sum = post.sum(axis=1, keepdims=True)
            bad = (row_sum == 0) | ~np.isfinite(row_sum)
            if bad.any():
                row_sum[bad] = 1.0
                post /= row_sum
                post[bad.ravel()] = 0.0
            else:
                post /= row_sum
            em_totals[comp_idx] += post.sum(axis=0)
        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total
        delta = np.abs(theta_new - theta).sum()
        theta = theta_new
        if delta < _EM_CONVERGENCE_DELTA:
            break
    elapsed = time.perf_counter() - t0
    return theta, unique_totals + em_totals + prior, elapsed


def approach_c_prealloc_csr(locus_em, n_iters):
    """Pre-allocated flat CSR: np.take + bincount."""
    theta, prior, unique_totals, log_eff_len = common_init(locus_em)
    n_total = locus_em.n_components
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks
    n_units = len(offsets) - 1
    seg_lengths = np.diff(offsets)

    unit_idx = np.repeat(np.arange(n_units, dtype=np.intp), seg_lengths)
    offsets_int = offsets[:-1].astype(np.intp)
    n_cand = len(t_indices)
    _buf = np.empty(n_cand, dtype=np.float64)
    _buf2 = np.empty(n_cand, dtype=np.float64)

    t0 = time.perf_counter()
    em_totals = np.zeros(n_total, dtype=np.float64)
    for _ in range(n_iters):
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        np.take(log_weights, t_indices, out=_buf)
        _buf += log_liks
        seg_max = np.maximum.reduceat(_buf, offsets_int)
        np.take(seg_max, unit_idx, out=_buf2)
        _buf -= _buf2
        np.exp(_buf, out=_buf)
        seg_sum = np.add.reduceat(_buf, offsets_int)
        bad_seg = (seg_sum == 0) | ~np.isfinite(seg_sum)
        seg_sum[bad_seg] = 1.0
        np.take(seg_sum, unit_idx, out=_buf2)
        _buf /= _buf2
        if bad_seg.any():
            _buf[bad_seg[unit_idx]] = 0.0
        em_totals = np.bincount(
            t_indices, weights=_buf, minlength=n_total,
        ).astype(np.float64)
        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total
        delta = np.abs(theta_new - theta).sum()
        theta = theta_new
        if delta < _EM_CONVERGENCE_DELTA:
            break
    elapsed = time.perf_counter() - t0
    return theta, unique_totals + em_totals + prior, elapsed


def approach_d_equiv_prealloc(locus_em, n_iters):
    """Equiv-class with pre-allocated scratch arrays per class."""
    theta, prior, unique_totals, log_eff_len = common_init(locus_em)
    n_total = locus_em.n_components
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks

    ec_list = _build_equiv_classes(offsets, t_indices, log_liks)
    # Pre-allocate scratch for each class
    ec_data = []
    for comp_idx, ll_matrix in ec_list:
        scratch = np.empty_like(ll_matrix)
        row_sum_buf = np.empty((ll_matrix.shape[0], 1), dtype=np.float64)
        ec_data.append((comp_idx, ll_matrix, scratch, row_sum_buf))

    t0 = time.perf_counter()
    em_totals = np.zeros(n_total, dtype=np.float64)
    for _ in range(n_iters):
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        em_totals[:] = 0.0
        for comp_idx, ll_matrix, scratch, rs_buf in ec_data:
            lw = log_weights[comp_idx]
            np.add(ll_matrix, lw, out=scratch)
            max_r = scratch.max(axis=1, keepdims=True)
            scratch -= max_r
            np.exp(scratch, out=scratch)
            np.sum(scratch, axis=1, keepdims=True, out=rs_buf)
            scratch /= rs_buf
            em_totals[comp_idx] += scratch.sum(axis=0)
        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total
        delta = np.abs(theta_new - theta).sum()
        theta = theta_new
        if delta < _EM_CONVERGENCE_DELTA:
            break
    elapsed = time.perf_counter() - t0
    return theta, unique_totals + em_totals + prior, elapsed


def approach_e_hybrid(locus_em, n_iters):
    """Hybrid: dense matrices for large classes, flat CSR for small ones."""
    theta, prior, unique_totals, log_eff_len = common_init(locus_em)
    n_total = locus_em.n_components
    offsets = locus_em.offsets
    t_indices = locus_em.t_indices
    log_liks = locus_em.log_liks
    n_units = len(offsets) - 1
    seg_lengths = np.diff(offsets)

    ec_list = _build_equiv_classes(offsets, t_indices, log_liks)

    # Split into large and small
    THRESHOLD = 100
    large_cls = []
    small_offsets_list = [0]
    small_t_list = []
    small_ll_list = []
    for comp_idx, ll_matrix in ec_list:
        n = ll_matrix.shape[0]
        if n >= THRESHOLD:
            scratch = np.empty_like(ll_matrix)
            large_cls.append((comp_idx, ll_matrix, scratch))
        else:
            k = len(comp_idx)
            for row_i in range(n):
                small_t_list.extend(comp_idx.tolist())
                small_ll_list.extend(ll_matrix[row_i].tolist())
                small_offsets_list.append(len(small_t_list))

    has_small = len(small_t_list) > 0
    if has_small:
        s_offsets = np.array(small_offsets_list, dtype=np.int64)
        s_t_idx = np.array(small_t_list, dtype=np.int32)
        s_ll = np.array(small_ll_list, dtype=np.float64)
        s_n = len(s_offsets) - 1
        s_seg_len = np.diff(s_offsets)
        s_unit_idx = np.repeat(np.arange(s_n, dtype=np.intp), s_seg_len)
        s_off_int = s_offsets[:-1].astype(np.intp)
        s_buf = np.empty(len(s_t_idx), dtype=np.float64)
        s_buf2 = np.empty(len(s_t_idx), dtype=np.float64)

    n_large = sum(m.shape[0] for _, m, _ in large_cls)
    n_small = sum(1 for _ in range(len(small_offsets_list) - 1)) if has_small else 0
    print(f"  Hybrid split: {len(large_cls)} large classes ({n_large} units), "
          f"CSR batch ({n_small} units)")

    t0 = time.perf_counter()
    em_totals = np.zeros(n_total, dtype=np.float64)
    for _ in range(n_iters):
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        em_totals[:] = 0.0

        # Large classes: dense matrix
        for comp_idx, ll_matrix, scratch in large_cls:
            lw = log_weights[comp_idx]
            np.add(ll_matrix, lw, out=scratch)
            max_r = scratch.max(axis=1, keepdims=True)
            scratch -= max_r
            np.exp(scratch, out=scratch)
            row_sum = scratch.sum(axis=1, keepdims=True)
            scratch /= row_sum
            em_totals[comp_idx] += scratch.sum(axis=0)

        # Small classes: flat CSR batch
        if has_small:
            np.take(log_weights, s_t_idx, out=s_buf)
            s_buf += s_ll
            seg_max = np.maximum.reduceat(s_buf, s_off_int)
            np.take(seg_max, s_unit_idx, out=s_buf2)
            s_buf -= s_buf2
            np.exp(s_buf, out=s_buf)
            seg_sum = np.add.reduceat(s_buf, s_off_int)
            bad = (seg_sum == 0) | ~np.isfinite(seg_sum)
            seg_sum[bad] = 1.0
            np.take(seg_sum, s_unit_idx, out=s_buf2)
            s_buf /= s_buf2
            if bad.any():
                s_buf[bad[s_unit_idx]] = 0.0
            np.add.at(em_totals, s_t_idx, s_buf)

        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total
        delta = np.abs(theta_new - theta).sum()
        theta = theta_new
        if delta < _EM_CONVERGENCE_DELTA:
            break
    elapsed = time.perf_counter() - t0
    return theta, unique_totals + em_totals + prior, elapsed


def main():
    bam = '/Users/mkiyer/Downloads/hulkrna_runs/bench_pristine_10_regions/chr17_43044295_43170245/gdna_none_nrna_none_ss_1.00/align_oracle/reads_namesort.bam'
    idx_dir = '/Users/mkiyer/Downloads/hulkrna_runs/bench_pristine_10_regions/chr17_43044295_43170245/hulkrna_index/'

    print("Loading data...")
    locus_em = setup(bam, idx_dir)
    n_units = len(locus_em.offsets) - 1
    n_cand = len(locus_em.t_indices)
    print(f"  {n_units:,} units, {n_cand:,} candidates, "
          f"{locus_em.n_components} components\n")

    N = 100  # iterations for timing

    print("=" * 60)
    print(f"Benchmarking {N} EM iterations each")
    print("=" * 60)

    approaches = [
        ("A: Original CSR (repeat+add.at)", approach_a_original),
        ("B: Equiv-class (current)", approach_b_equiv_class),
        ("C: Pre-alloc CSR (take+bincount)", approach_c_prealloc_csr),
        ("D: Equiv-class + pre-alloc scratch", approach_d_equiv_prealloc),
        ("E: Hybrid (dense large + CSR small)", approach_e_hybrid),
    ]

    results = {}
    for name, fn in approaches:
        print(f"\n{name}")
        theta, alpha, elapsed = fn(locus_em, N)
        ms_per_iter = elapsed / N * 1000
        print(f"  {elapsed:.3f}s total, {ms_per_iter:.2f} ms/iter")
        results[name] = (theta, alpha, elapsed)

    # Cross-check: all should produce the same theta
    print("\n" + "=" * 60)
    print("Correctness check")
    print("=" * 60)
    ref_theta = results[approaches[0][0]][0]
    for name, _ in approaches[1:]:
        other_theta = results[name][0]
        diff = np.abs(ref_theta - other_theta).max()
        match = "OK" if diff < 1e-6 else f"MISMATCH max_diff={diff:.2e}"
        print(f"  {name}: {match}")


if __name__ == "__main__":
    main()
