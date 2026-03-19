#!/usr/bin/env python3
"""Diagnose the nRNA siphon effect.

Runs a single simulation scenario with NTA1 >> TA1 and extracts
per-fragment scoring data to identify exactly why nRNA absorbs
mRNA fragments.

Usage:
    conda activate rigel
    python scripts/debug/diagnose_nrna_siphon.py
"""

import gc
import logging
import sys
import tempfile
from collections import Counter
from dataclasses import replace as _replace
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure rigel is importable
_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import (
    BamScanConfig,
    EMConfig,
    FragmentScoringConfig,
    PipelineConfig,
)
from rigel.sim import GDNAConfig, Scenario, SimConfig

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def build_scenario(mrna_ab=20, nrna_ab=500, n_rna=10000, ss=1.0, seed=42):
    """Build a scenario matching the benchmark geometry."""
    sc = Scenario(
        "siphon",
        genome_length=50000,
        seed=seed,
        work_dir=Path(tempfile.mkdtemp(prefix="rigel_siphon_")),
    )
    # Gene A: multi-exon (same as benchmark config)
    sc.add_gene("GA", "+", [
        {
            "t_id": "TA1",
            "exons": [(1000, 1020), (5000, 5500), (12000, 13000)],
            "abundance": mrna_ab,
            "nrna_abundance": nrna_ab,
        },
    ])
    # Gene D: single-exon (unexpressed control)
    sc.add_gene("GD", "+", [
        {
            "t_id": "TD",
            "exons": [(28000, 29000)],
            "abundance": 0,
            "nrna_abundance": 0,
        },
    ])
    # Gene E: 2-exon negative-strand (unexpressed control)
    sc.add_gene("GE", "-", [
        {
            "t_id": "TE",
            "exons": [(40000, 41000), (45000, 46000)],
            "abundance": 0,
            "nrna_abundance": 0,
        },
    ])
    sim_cfg = SimConfig(
        frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
        read_length=150, strand_specificity=ss, seed=seed,
    )
    result = sc.build_oracle(
        n_fragments=n_rna,
        sim_config=sim_cfg,
        nrna_abundance=0.0,        # use per-transcript values
        n_rna_fragments=n_rna,
        gdna_fraction=0.0,
    )
    return sc, result


def extract_scored_fragments(result, pipe_cfg):
    """Run scan + scoring to get ScoredFragments WITHOUT running EM."""
    import os
    from rigel.buffer import FragmentBuffer
    from rigel.estimator import AbundanceEstimator
    from rigel.frag_length_model import FragmentLengthModels
    from rigel.pipeline import (
        PipelineStats,
        scan_and_buffer,
        _setup_geometry_and_estimator,
    )
    from rigel.scan import FragmentRouter
    from rigel.scoring import FragmentScorer
    from rigel.strand_model import StrandModels

    index = result.index
    bam_path = str(result.bam_path)
    scan = pipe_cfg.scan

    # Detect sj_strand_tag auto
    if scan.sj_strand_tag == "auto":
        from rigel.native import detect_sj_strand_tag as _detect
        detected = _detect(bam_path)
        scan = _replace(scan, sj_strand_tag=detected)

    stats, strand_models, fl_models, buffer, region_counts, fl_table = scan_and_buffer(
        bam_path, index, scan,
    )

    # Finalize models
    strand_models.finalize()
    fl_models.build_scoring_models()
    fl_models.finalize()

    # Geometry + estimator
    geometry, estimator = _setup_geometry_and_estimator(index, fl_models, pipe_cfg.em)

    # Build scorer
    ctx = FragmentScorer.from_models(
        strand_models, fl_models, index, estimator,
        overhang_log_penalty=pipe_cfg.scoring.overhang_log_penalty,
        mismatch_log_penalty=pipe_cfg.scoring.mismatch_log_penalty,
        gdna_splice_penalties=pipe_cfg.scoring.gdna_splice_penalties,
    )

    # Build router and scan fragments
    router = FragmentRouter(
        ctx, estimator, stats, index, strand_models, annotations=None,
    )
    em_data = router.scan(buffer, log_every=1_000_000)

    return em_data, stats, strand_models, fl_models, estimator, geometry


def analyze_scoring(em_data, index, fl_models, estimator, geometry):
    """Analyze per-fragment mRNA vs nRNA scoring differences."""
    offsets = em_data.offsets
    ti = em_data.t_indices
    ll = em_data.log_liks
    cw = em_data.coverage_weights
    splice = em_data.splice_type
    n_units = em_data.n_units
    nrna_base = em_data.nrna_base_index
    n_t = nrna_base  # num_transcripts

    print(f"\n{'=' * 70}")
    print(f"SCORED FRAGMENTS ANALYSIS")
    print(f"{'=' * 70}")
    print(f"  Total EM units: {n_units}")
    print(f"  Total candidates: {em_data.n_candidates}")
    print(f"  nrna_base_index: {nrna_base}")
    print(f"  Deterministic unambig (from stats): skipped, only EM units shown")

    # Classify each unit
    both_count = 0
    mrna_only_count = 0
    nrna_only_count = 0

    records = []  # per-unit analysis

    for u in range(n_units):
        s, e = offsets[u], offsets[u + 1]
        u_ti = ti[s:e]
        u_ll = ll[s:e]
        u_cw = cw[s:e]

        m_mask = u_ti < n_t
        n_mask = u_ti >= n_t

        has_mrna = m_mask.any()
        has_nrna = n_mask.any()

        if has_mrna and has_nrna:
            both_count += 1
            m_best = u_ll[m_mask].max()
            n_best = u_ll[n_mask].max()
            m_cw = u_cw[m_mask][u_ll[m_mask].argmax()]
            n_cw = u_cw[n_mask][u_ll[n_mask].argmax()]
            records.append({
                "unit": u,
                "splice": int(splice[u]),
                "n_cands": int(e - s),
                "mrna_ll": m_best,
                "nrna_ll": n_best,
                "ll_diff": n_best - m_best,  # positive = nRNA favored
                "mrna_cw": m_cw,
                "nrna_cw": n_cw,
                "cw_diff": n_cw - m_cw,
            })
        elif has_mrna:
            mrna_only_count += 1
        elif has_nrna:
            nrna_only_count += 1

    print(f"\n--- Unit routing ---")
    print(f"  Both mRNA + nRNA candidates: {both_count}")
    print(f"  mRNA only:                   {mrna_only_count}")
    print(f"  nRNA only:                   {nrna_only_count}")

    if not records:
        print("  No fragments with both candidates — nothing to analyze.")
        return

    df = pd.DataFrame(records)

    # LL difference analysis
    print(f"\n--- Log-likelihood comparison (nRNA_ll - mRNA_ll) ---")
    print(f"  {len(df)} units with both mRNA and nRNA candidates:")
    print(f"  mean diff:   {df['ll_diff'].mean():+.6f}")
    print(f"  median diff: {df['ll_diff'].median():+.6f}")
    print(f"  std diff:    {df['ll_diff'].std():.6f}")
    print(f"  min diff:    {df['ll_diff'].min():+.6f}")
    print(f"  max diff:    {df['ll_diff'].max():+.6f}")

    eps = 1e-10
    n_equal = (df["ll_diff"].abs() < eps).sum()
    n_nrna_better = (df["ll_diff"] > eps).sum()
    n_mrna_better = (df["ll_diff"] < -eps).sum()
    print(f"\n  Exactly equal:   {n_equal:6d} ({n_equal/len(df)*100:.1f}%)")
    print(f"  nRNA favored:    {n_nrna_better:6d} ({n_nrna_better/len(df)*100:.1f}%)")
    print(f"  mRNA favored:    {n_mrna_better:6d} ({n_mrna_better/len(df)*100:.1f}%)")

    # Histogram
    bins = [-50, -10, -5, -2, -1, -0.5, -0.1, -0.01, 0.01, 0.1, 0.5, 1, 2, 5, 10, 50]
    hist, edges = np.histogram(df["ll_diff"], bins=bins)
    print(f"\n  Distribution of ll_diff (nRNA - mRNA):")
    mx = max(hist) if max(hist) > 0 else 1
    for i in range(len(hist)):
        bar = "#" * min(hist[i] * 60 // mx, 60)
        print(f"    [{edges[i]:+7.2f}, {edges[i+1]:+7.2f}): {hist[i]:6d} {bar}")

    # Coverage weight analysis
    print(f"\n--- Coverage weight comparison ---")
    print(f"  mRNA cw: mean={df['mrna_cw'].mean():.6f}, "
          f"std={df['mrna_cw'].std():.6f}")
    print(f"  nRNA cw: mean={df['nrna_cw'].mean():.6f}, "
          f"std={df['nrna_cw'].std():.6f}")
    print(f"  cw_diff (nRNA - mRNA): mean={df['cw_diff'].mean():+.6f}, "
          f"std={df['cw_diff'].std():.6f}")

    # Effective length analysis
    print(f"\n--- Effective length analysis ---")
    t_lengths = index.t_df["length"].values
    t_spans = (index.t_df["end"].values - index.t_df["start"].values)
    print(f"  TA1 exonic length: {t_lengths[0]}")
    print(f"  TA1 genomic span:  {t_spans[0]}")

    eff_lens = geometry.effective_lengths
    tx_spans = geometry.transcript_spans
    print(f"  TA1 effective length (mRNA): {eff_lens[0]:.1f}")

    if index.nrna_df is not None and len(index.nrna_df) > 0:
        nrna_spans = index.nrna_df["end"].values - index.nrna_df["start"].values
        print(f"  NTA1 nRNA span: {nrna_spans[0]}")

        fl_mean = fl_models.rna_model.mean
        eff_nrna = max(nrna_spans[0] - fl_mean + 1, 1)
        print(f"  NTA1 nRNA effective length: {eff_nrna:.1f}")
        print(f"  Ratio (nRNA_eff / mRNA_eff): {eff_nrna / eff_lens[0]:.2f}x")

    # Break down by actual LL diff bucket
    print(f"\n--- Fragment examples (sorted by ll_diff) ---")
    df_sorted = df.sort_values("ll_diff")
    print(f"\n  Most mRNA-favored (top 5):")
    for _, r in df_sorted.head(5).iterrows():
        print(f"    unit={int(r['unit']):5d} spl={int(r['splice'])} "
              f"mrna_ll={r['mrna_ll']:.4f} nrna_ll={r['nrna_ll']:.4f} "
              f"diff={r['ll_diff']:+.4f} "
              f"mrna_cw={r['mrna_cw']:.6f} nrna_cw={r['nrna_cw']:.6f}")

    print(f"\n  Most nRNA-favored (top 5):")
    for _, r in df_sorted.tail(5).iterrows():
        print(f"    unit={int(r['unit']):5d} spl={int(r['splice'])} "
              f"mrna_ll={r['mrna_ll']:.4f} nrna_ll={r['nrna_ll']:.4f} "
              f"diff={r['ll_diff']:+.4f} "
              f"mrna_cw={r['mrna_cw']:.6f} nrna_cw={r['nrna_cw']:.6f}")

    print(f"\n  Near-equal examples (|diff| < 0.01, first 5):")
    near_eq = df[df["ll_diff"].abs() < 0.01]
    for _, r in near_eq.head(5).iterrows():
        print(f"    unit={int(r['unit']):5d} spl={int(r['splice'])} "
              f"mrna_ll={r['mrna_ll']:.4f} nrna_ll={r['nrna_ll']:.4f} "
              f"diff={r['ll_diff']:+.6f} "
              f"mrna_cw={r['mrna_cw']:.6f} nrna_cw={r['nrna_cw']:.6f}")

    return df


def run_pipeline_and_compare(result, pipe_cfg):
    """Run full pipeline and compare with ground truth."""
    from rigel.pipeline import run_pipeline

    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)

    # Ground truth
    gt_mrna = result.ground_truth_from_bam()
    gt_nrna = result.ground_truth_nrna_count_from_bam()

    print(f"\n{'=' * 70}")
    print(f"PIPELINE RESULTS vs GROUND TRUTH")
    print(f"{'=' * 70}")

    t_df = pr.estimator.get_counts_df(result.index)
    print(f"\n  Transcript-level output:")
    print(f"  {'t_id':10s} {'mrna_obs':>10s} {'mrna_exp':>10s} "
          f"{'nrna_obs':>10s} {'nrna_exp':>10s}")
    for _, row in t_df.iterrows():
        tid = row["transcript_id"]
        mrna_exp = gt_mrna.get(tid, 0)
        print(f"  {tid:10s} {row['mrna']:10.1f} {mrna_exp:10d} "
              f"{row['nrna']:10.1f} {gt_nrna:10d}")

    mrna_obs = t_df["mrna"].sum()
    mrna_exp = sum(gt_mrna.values())
    nrna_obs = t_df["nrna"].sum()
    nrna_exp = gt_nrna

    print(f"\n  Summary:")
    print(f"    mRNA: expected={mrna_exp}, observed={mrna_obs:.1f}, "
          f"error={abs(mrna_obs - mrna_exp) / max(mrna_exp, 1):.3f}")
    print(f"    nRNA: expected={nrna_exp}, observed={nrna_obs:.1f}, "
          f"error={abs(nrna_obs - nrna_exp) / max(nrna_exp, 1):.3f}")
    siphon = mrna_exp - mrna_obs
    print(f"    Siphon: {siphon:+.1f} fragments moved from mRNA → nRNA")

    return pr


def analyze_post_bias(em_data, index):
    """Compute post-bias-correction LL differences and theoretical equilibrium."""
    offsets = em_data.offsets
    ti = em_data.t_indices
    ll = em_data.log_liks
    tx_s = em_data.tx_starts
    tx_e = em_data.tx_ends
    n_units = em_data.n_units
    nrna_base = em_data.nrna_base_index
    n_t = nrna_base

    # Profile lengths for bias correction
    t_lengths = index.t_df["length"].values  # exonic lengths
    t_starts = index.t_df["start"].values
    t_ends = index.t_df["end"].values

    nrna_starts = index.nrna_df["start"].values if index.nrna_df is not None else np.array([])
    nrna_ends = index.nrna_df["end"].values if index.nrna_df is not None else np.array([])

    # Build profile_lengths array: [mRNA_0, ..., mRNA_{n_t-1}, nRNA_0, ...]
    n_nrna = len(nrna_starts) if len(nrna_starts) > 0 else 0

    print(f"\n{'=' * 70}")
    print(f"POST-BIAS-CORRECTION ANALYSIS")
    print(f"{'=' * 70}")
    print(f"  n_transcripts={n_t}, n_nrna_spans={n_nrna}")

    for i in range(n_t):
        print(f"  mRNA[{i}]: exonic_len={t_lengths[i]}, "
              f"genomic=[{t_starts[i]},{t_ends[i]})")
    for i in range(n_nrna):
        span = nrna_ends[i] - nrna_starts[i]
        print(f"  nRNA[{i}]: span={span}, range=[{nrna_starts[i]},{nrna_ends[i]})")

    # Compute per-candidate bias correction and post-bias LLs
    both_records = []
    nrna_only_count = 0
    mrna_only_count = 0

    for u in range(n_units):
        s, e = offsets[u], offsets[u + 1]
        u_ti = ti[s:e]
        u_ll = ll[s:e].copy()
        u_tx_s = tx_s[s:e]
        u_tx_e = tx_e[s:e]

        m_mask = u_ti < n_t
        n_mask = u_ti >= n_t

        if not m_mask.any() or not n_mask.any():
            if m_mask.any():
                mrna_only_count += 1
            else:
                nrna_only_count += 1
            continue

        # Apply bias correction to copies
        corrected = np.zeros(len(u_ll))
        for j in range(len(u_ll)):
            frag_len = max(u_tx_e[j] - u_tx_s[j], 0)
            c_idx = u_ti[j]
            if c_idx < n_t:
                prof_len = t_lengths[c_idx]
            else:
                nrna_idx = c_idx - n_t
                prof_len = nrna_ends[nrna_idx] - nrna_starts[nrna_idx]
            eff_len = max(prof_len - frag_len + 1, 1)
            corrected[j] = u_ll[j] - np.log(eff_len)

        m_best_raw = u_ll[m_mask].max()
        n_best_raw = u_ll[n_mask].max()
        m_best_corr = corrected[m_mask].max()
        n_best_corr = corrected[n_mask].max()

        # Fragment length (from the mRNA candidate)
        m_idx = np.where(m_mask)[0][u_ll[m_mask].argmax()]
        frag_len = u_tx_e[m_idx] - u_tx_s[m_idx]

        both_records.append({
            "unit": u,
            "frag_len": int(frag_len),
            "raw_diff": n_best_raw - m_best_raw,
            "corr_diff": n_best_corr - m_best_corr,
            "m_raw": m_best_raw,
            "n_raw": n_best_raw,
            "m_corr": m_best_corr,
            "n_corr": n_best_corr,
        })

    df = pd.DataFrame(both_records)

    print(f"\n--- Post-bias-correction comparison ({len(df)} units) ---")
    print(f"  nRNA-only: {nrna_only_count}, mRNA-only: {mrna_only_count}")

    eps = 1e-6
    n_eq = (df["corr_diff"].abs() < eps).sum()
    n_nrna_fav = (df["corr_diff"] > eps).sum()
    n_mrna_fav = (df["corr_diff"] < -eps).sum()
    print(f"\n  After bias correction:")
    print(f"    Equal:         {n_eq:6d}")
    print(f"    nRNA favored:  {n_nrna_fav:6d}")
    print(f"    mRNA favored:  {n_mrna_fav:6d}")
    print(f"    corr_diff stats: mean={df['corr_diff'].mean():+.4f}, "
          f"median={df['corr_diff'].median():+.4f}, "
          f"std={df['corr_diff'].std():.4f}")

    # Bucket the corrected differences
    bins2 = [-2000, -100, -10, -5, -2, -1, 0, 1, 2, 5, 10, 100, 2000]
    hist2, edges2 = np.histogram(df["corr_diff"], bins=bins2)
    print(f"\n  Distribution of corr_diff (nRNA - mRNA, post-bias):")
    mx = max(hist2) if max(hist2) > 0 else 1
    for i in range(len(hist2)):
        bar = "#" * min(hist2[i] * 50 // mx, 50)
        print(f"    [{edges2[i]:+8.0f}, {edges2[i+1]:+8.0f}): {hist2[i]:6d} {bar}")

    # Fragment length distribution for "equal" raw fragments
    raw_equal = df[df["raw_diff"].abs() < 1e-10]
    if len(raw_equal) > 0:
        print(f"\n  Fragment lengths for {len(raw_equal)} raw-equal fragments:")
        print(f"    min={raw_equal['frag_len'].min()}, "
              f"max={raw_equal['frag_len'].max()}, "
              f"mean={raw_equal['frag_len'].mean():.1f}")
        print(f"    corr_diff for raw-equal: mean={raw_equal['corr_diff'].mean():+.4f}")

    # ===== Theoretical equilibrium =====
    print(f"\n--- Theoretical MAP-EM equilibrium (2-component model) ---")
    # For raw-equal fragments: compute average bias correction advantage
    if len(raw_equal) > 0:
        avg_m_corr = raw_equal["m_corr"].mean()
        avg_n_corr = raw_equal["n_corr"].mean()
        bias_adv = avg_m_corr - avg_n_corr
        print(f"  Avg mRNA corr LL:  {avg_m_corr:.4f}")
        print(f"  Avg nRNA corr LL:  {avg_n_corr:.4f}")
        print(f"  Avg mRNA bias advantage: {bias_adv:+.4f} nats")

    n_det = 5  # deterministic spliced

    # Build arrays for vectorized equilibrium computation
    m_corr_arr = df["m_corr"].values
    n_corr_arr = df["n_corr"].values

    # Initialize theta (warm-start: distribute coverage proportionally)
    theta_m = float(n_det + len(df) * 0.5)  # initial guess
    theta_n = float(nrna_only_count + len(df) * 0.5)

    print(f"\n  EM iteration trace (2-component MAP-EM):")
    print(f"  {'iter':>4s} {'theta_m':>10s} {'theta_n':>10s} "
          f"{'em_m':>10s} {'em_n':>10s} {'delta':>10s}")

    for iteration in range(500):
        log_tm = np.log(max(theta_m, 1e-300))
        log_tn = np.log(max(theta_n, 1e-300))

        # Vectorized E-step
        log_pm = m_corr_arr + log_tm
        log_pn = n_corr_arr + log_tn
        max_lp = np.maximum(log_pm, log_pn)
        pm = np.exp(log_pm - max_lp)
        pn = np.exp(log_pn - max_lp)
        total_p = pm + pn
        mask = total_p > 0
        post_m = np.where(mask, pm / total_p, 0.0)

        em_m = post_m.sum()
        em_n = (1.0 - post_m[mask]).sum()

        theta_m_new = n_det + em_m
        theta_n_new = nrna_only_count + em_n

        delta = abs(theta_m_new - theta_m)
        if iteration < 10 or iteration % 50 == 0 or delta < 0.001:
            print(f"  {iteration:4d} {theta_m:10.3f} {theta_n:10.3f} "
                  f"{em_m:10.3f} {em_n:10.3f} {delta:10.6f}")

        if delta < 0.001 and iteration > 5:
            print(f"  ** Converged at iteration {iteration}")
            break

        theta_m = theta_m_new
        theta_n = theta_n_new

    p_final = em_m / len(df) if len(df) > 0 else 0

    print(f"\n  Final equilibrium:")
    print(f"    theta_mRNA = {theta_m:.4f}")
    print(f"    theta_nRNA = {theta_n:.4f}")
    print(f"    EM→mRNA from shared = {em_m:.4f}")
    print(f"    EM→nRNA from shared = {em_n:.4f}")
    print(f"    Total mRNA = {n_det} (det) + {em_m:.2f} (EM) = {n_det + em_m:.2f}")
    print(f"    Avg mRNA posterior for shared = {p_final:.4f}")

    # Also test: what if nRNA-only fragments DON'T dilute theta?
    # (hypothetical: if nRNA-only fragments were in a separate component)
    print(f"\n--- Sensitivity: what if nRNA-only excluded from theta_n? ---")
    theta_m2 = float(n_det + len(df) * 0.5)
    theta_n2 = float(len(df) * 0.5)
    for iteration in range(500):
        log_tm = np.log(max(theta_m2, 1e-300))
        log_tn = np.log(max(theta_n2, 1e-300))
        log_pm = m_corr_arr + log_tm
        log_pn = n_corr_arr + log_tn
        max_lp = np.maximum(log_pm, log_pn)
        pm = np.exp(log_pm - max_lp)
        pn = np.exp(log_pn - max_lp)
        total_p = pm + pn
        mask = total_p > 0
        post_m = np.where(mask, pm / total_p, 0.0)
        em_m2 = post_m.sum()
        em_n2 = (1.0 - post_m[mask]).sum()
        theta_m2_new = n_det + em_m2
        theta_n2_new = em_n2  # NO nrna_only_count
        delta = abs(theta_m2_new - theta_m2)
        if delta < 0.001 and iteration > 5:
            break
        theta_m2 = theta_m2_new
        theta_n2 = theta_n2_new
    print(f"    Total mRNA (no dilution) = {n_det + em_m2:.2f}")
    print(f"    theta_m={theta_m2:.2f}, theta_n={theta_n2:.2f}")
    print(f"    Avg mRNA posterior = {em_m2/len(df):.4f}")

    return df


def main():
    # =====================================================================
    # Scenario parameters
    # =====================================================================
    MRNA_AB = 20
    NRNA_AB = 500
    N_RNA = 10000
    SS = 1.0
    SEED = 42

    print("=" * 70)
    print(f"nRNA SIPHON DIAGNOSTIC")
    print(f"Scenario: mRNA_ab={MRNA_AB}, nRNA_ab={NRNA_AB}, "
          f"n_rna={N_RNA}, ss={SS}")
    print("=" * 70)

    # Build scenario
    sc, result = build_scenario(MRNA_AB, NRNA_AB, N_RNA, SS, SEED)

    try:
        # ========= Run with VBEM =========
        pipe_cfg_vbem = PipelineConfig(
            em=EMConfig(seed=SEED, mode="vbem"),
            scan=BamScanConfig(sj_strand_tag="auto"),
            scoring=FragmentScoringConfig(),
        )
        print(f"\n{'#' * 70}")
        print(f"# MODE: VBEM")
        print(f"{'#' * 70}")
        pr_vbem = run_pipeline_and_compare(result, pipe_cfg_vbem)

        # ========= Run with MAP-EM =========
        pipe_cfg_map = PipelineConfig(
            em=EMConfig(seed=SEED, mode="map"),
            scan=BamScanConfig(sj_strand_tag="auto"),
            scoring=FragmentScoringConfig(),
        )
        print(f"\n{'#' * 70}")
        print(f"# MODE: MAP-EM")
        print(f"{'#' * 70}")
        pr_map = run_pipeline_and_compare(result, pipe_cfg_map)

        # ========= Extract scored fragments for deep analysis =========
        em_data, stats, strand_models, fl_models, estimator, geometry = \
            extract_scored_fragments(result, pipe_cfg_vbem)

        # Pre-bias scoring analysis
        df_raw = analyze_scoring(
            em_data, result.index, fl_models, estimator, geometry,
        )

        # Post-bias correction analysis + theoretical equilibrium
        df_bias = analyze_post_bias(em_data, result.index)

        # FL model analysis
        print(f"\n{'=' * 70}")
        print(f"FRAGMENT LENGTH MODEL")
        print(f"{'=' * 70}")
        fl_model = fl_models.rna_model
        print(f"  mean={fl_model.mean:.1f}, n_obs={fl_model.n_observations}")
        log_prob = fl_model._log_prob
        test_flens = [100, 150, 200, 250, 300, 400, 600, 1000]
        max_sz = fl_model.max_size
        print(f"  max_size={max_sz}")
        for fl in test_flens:
            if fl <= max_sz:
                lp = log_prob[fl]
            else:
                lp = float('-inf')
            print(f"    flen={fl:6d}: log_prob={lp:.4f}")

        # ========= Check calibration / gDNA info =========
        print(f"\n{'=' * 70}")
        print(f"gDNA CALIBRATION INFO")
        print(f"{'=' * 70}")
        est = pr_vbem.estimator
        print(f"  gdna_em_total (VBEM): {est._gdna_em_total:.4f}")
        est2 = pr_map.estimator
        print(f"  gdna_em_total (MAP):  {est2._gdna_em_total:.4f}")

        # Per-locus breakdown
        loci_df = est.get_loci_df()
        print(f"\n  Per-locus breakdown (VBEM):")
        print(loci_df.to_string(index=False))

        loci_df2 = est2.get_loci_df()
        print(f"\n  Per-locus breakdown (MAP):")
        print(loci_df2.to_string(index=False))

    finally:
        sc.cleanup()


if __name__ == "__main__":
    main()
