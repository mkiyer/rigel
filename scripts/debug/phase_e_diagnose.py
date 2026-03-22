#!/usr/bin/env python3
"""Diagnose the nRNA/gDNA pool shift in Phase C vs baseline.

Examines per-locus differences to understand where gDNA fragments
are being reassigned to nRNA transcripts.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────
BENCH_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")
OUTPUT_DIR = BENCH_DIR / "phase_e_results"
INDEX_DIR = OUTPUT_DIR / "rigel_index"

# Use the hardest condition: ss=1.00 (fully stranded) with high gDNA
COND = "gdna_high_ss_1.00_nrna_default"
BAM = BENCH_DIR / COND / "align_oracle" / "reads_namesort.bam"


def main():
    from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import run_pipeline

    print("Loading index...", flush=True)
    index = TranscriptIndex.load(str(INDEX_DIR))
    print(f"  {index.num_transcripts} transcripts, "
          f"{int(index.t_df['is_synthetic_nrna'].sum())} synthetic")

    print(f"\nRunning rigel on {COND}...", flush=True)
    cfg = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto", include_multimap=True),
        scoring=FragmentScoringConfig(),
    )
    t0 = time.monotonic()
    pipe = run_pipeline(BAM, index, config=cfg)
    elapsed = time.monotonic() - t0
    print(f"  Done in {elapsed:.1f}s")

    est = pipe.estimator
    stats = pipe.stats

    # Pool totals
    counts_df = est.get_counts_df(index)
    syn_mask = index.t_df["is_synthetic_nrna"].values
    t_total = est.t_counts.sum(axis=1)  # raw transcript-level totals

    mrna_total = float(t_total[~syn_mask].sum())
    nrna_total = float(t_total[syn_mask].sum())
    gdna_total = float(stats.n_gdna_total)
    print(f"\n  Pool totals:")
    print(f"    mRNA (non-synthetic): {mrna_total:.1f}")
    print(f"    nRNA (synthetic):     {nrna_total:.1f}")
    print(f"    gDNA:                 {gdna_total:.1f}")
    print(f"    Sum:                  {mrna_total + nrna_total + gdna_total:.1f}")

    # Loci summary
    loci_df = est.get_loci_df()
    print(f"\n  Loci: {len(loci_df)} total")
    print(f"    Loci columns: {loci_df.columns.tolist()}")

    # How many synthetic transcripts have >0 counts?
    syn_counts = t_total[syn_mask]
    n_active_syn = int((syn_counts > 0).sum())
    print(f"\n  Synthetic nRNA transcripts with >0 counts: "
          f"{n_active_syn} / {syn_mask.sum()}")
    print(f"    Total nRNA count: {syn_counts.sum():.1f}")
    print(f"    Mean count (active): {syn_counts[syn_counts > 0].mean():.2f}" if n_active_syn > 0 else "")
    print(f"    Max count: {syn_counts.max():.2f}")
    print(f"    Median count (active): {np.median(syn_counts[syn_counts > 0]):.2f}" if n_active_syn > 0 else "")

    # Top synthetic nRNA transcripts by count
    syn_indices = np.where(syn_mask)[0]
    top_idx = syn_indices[np.argsort(-syn_counts)][:20]
    print(f"\n  Top 20 synthetic nRNA transcripts by count:")
    for i, ti in enumerate(top_idx):
        tid = index.t_df.iloc[ti]["t_id"]
        gid = index.t_df.iloc[ti]["g_id"]
        gname = index.t_df.iloc[ti]["g_name"]
        cnt = t_total[ti]
        print(f"    {i+1:3d}. {tid:<50s} gene={gname:<15s} count={cnt:.1f}")

    # Fragment count conservation
    total_frags = stats.total
    assigned_total = mrna_total + nrna_total + gdna_total
    print(f"\n  Fragment conservation:")
    print(f"    Input fragments: {total_frags}")
    print(f"    Assigned:        {assigned_total:.1f}")
    print(f"    n_chimeric:      {stats.n_chimeric}")
    print(f"    n_gated_out:     {stats.n_gated_out}")

    # EM counts breakdown
    em_total = float(est.em_counts.sum())
    unambig_total = float(est.unambig_counts.sum())
    print(f"\n  EM breakdown:")
    print(f"    Unambiguous: {unambig_total:.1f}")
    print(f"    EM-assigned: {em_total:.1f}")
    print(f"    unambig nRNA: {float(est.unambig_counts.sum(axis=1)[syn_mask].sum()):.1f}")
    print(f"    EM nRNA:      {float(est.em_counts.sum(axis=1)[syn_mask].sum()):.1f}")


if __name__ == "__main__":
    main()
