#!/usr/bin/env python3
"""Diagnose per-transcript differences between baseline and Phase E results.

Focus: What kind of transcripts drive the MAE increase?
"""
from __future__ import annotations

import collections
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

BENCH_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")
SIM_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_salmon")
OUTPUT_DIR = BENCH_DIR / "phase_e_results"
GTF = "/Users/mkiyer/Downloads/rigel_runs/refs/human/genes_controls.gtf.gz"

COND = "gdna_high_ss_1.00_nrna_default"


def parse_truth(r1_path: Path) -> dict[str, int]:
    """Parse mRNA ground-truth from FASTQ."""
    counts: collections.Counter = collections.Counter()
    opener = gzip.open if str(r1_path).endswith(".gz") else open
    with opener(r1_path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            qname = line[1:].strip()
            if qname.endswith("/1"):
                qname = qname[:-2]
            t_id = qname.split(":")[0]
            if not t_id.startswith("gdna") and not t_id.startswith("nrna_"):
                counts[t_id] += 1
    return dict(counts)


def main():
    print(f"Condition: {COND}\n")

    # --- Load truth ---
    r1_path = SIM_DIR / COND / "sim_R1.fq.gz"
    print(f"Parsing truth from {r1_path.name}...", flush=True)
    truth = parse_truth(r1_path)
    print(f"  mRNA truth: {len(truth)} transcripts, {sum(truth.values())} frags")

    # --- Load baseline results ---
    base_csv = BENCH_DIR / COND / "per_transcript_counts_oracle.csv"
    if base_csv.exists():
        base_df = pd.read_csv(base_csv)
        print(f"\n  Baseline transcript counts: {len(base_df)} rows")
        print(f"  Columns: {base_df.columns.tolist()}")
        # The baseline column is 'rigel_default_oracle'
        base_df = base_df.rename(columns={"rigel_default_oracle": "mrna"})
        print(f"  mrna total: {base_df['mrna'].sum():.0f}")
    else:
        print(f"  Baseline CSV not found at {base_csv}")
        base_df = None

    # --- Load new results ---
    from rigel.index import TranscriptIndex

    index = TranscriptIndex.load(str(OUTPUT_DIR / "rigel_index"))
    
    from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
    from rigel.pipeline import run_pipeline

    bam_path = BENCH_DIR / COND / "align_oracle" / "reads_namesort.bam"
    cfg = PipelineConfig(
        em=EMConfig(seed=42),
        scan=BamScanConfig(sj_strand_tag="auto", include_multimap=True),
        scoring=FragmentScoringConfig(),
    )
    print(f"\n  Running rigel on {bam_path.name}...", flush=True)
    pipe = run_pipeline(bam_path, index, config=cfg)
    new_df = pipe.estimator.get_counts_df(index)
    
    # Filter out synthetics
    new_df_filt = new_df[~new_df["transcript_id"].str.startswith("RIGEL_NRNA_")].copy()
    print(f"  New transcript counts: {len(new_df_filt)} rows (after filtering synthetics)")
    
    # --- Compare new vs truth ---
    new_counts = dict(zip(new_df_filt["transcript_id"], new_df_filt["mrna"]))
    
    all_tids = sorted(set(truth) | set(new_counts))
    n = len(all_tids)
    truth_arr = np.array([truth.get(t, 0.0) for t in all_tids])
    new_arr = np.array([new_counts.get(t, 0.0) for t in all_tids])
    abs_err = np.abs(truth_arr - new_arr)
    
    # Top error transcripts
    err_df = pd.DataFrame({
        "t_id": all_tids,
        "truth": truth_arr,
        "new": new_arr,
        "abs_err": abs_err,
        "signed_err": new_arr - truth_arr,
    })
    err_df = err_df.sort_values("abs_err", ascending=False)
    
    print(f"\n  === TOP 30 ERROR TRANSCRIPTS (new vs truth) ===")
    print(f"  {'t_id':<25s} {'truth':>10s} {'new':>10s} {'abs_err':>10s} {'signed':>10s}")
    print("  " + "-" * 67)
    for _, row in err_df.head(30).iterrows():
        print(f"  {row['t_id']:<25s} {row['truth']:>10.1f} {row['new']:>10.1f} "
              f"{row['abs_err']:>10.1f} {row['signed_err']:>+10.1f}")
    
    # --- Error breakdown: positive vs negative ---
    over = err_df[err_df["signed_err"] > 0.5]
    under = err_df[err_df["signed_err"] < -0.5]
    print(f"\n  Over-estimated (new > truth): {len(over)} transcripts, total excess = {over['signed_err'].sum():.0f}")
    print(f"  Under-estimated (new < truth): {len(under)} transcripts, total deficit = {under['signed_err'].sum():.0f}")
    
    # --- Now compare baseline vs truth (if we have baseline) ---
    if base_df is not None and 'mrna' in base_df.columns:
        tid_col = 'transcript_id' if 'transcript_id' in base_df.columns else 't_id'
        base_counts = dict(zip(base_df[tid_col], base_df["mrna"]))
        base_arr = np.array([base_counts.get(t, 0.0) for t in all_tids])
        base_err = np.abs(truth_arr - base_arr)
        
        err_comp = pd.DataFrame({
            "t_id": all_tids,
            "truth": truth_arr,
            "base": base_arr,
            "new": new_arr,
            "base_err": base_err,
            "new_err": abs_err,
            "err_delta": abs_err - base_err,  # positive = regression
        })
        
        # Transcripts where error got WORSE
        worse = err_comp[err_comp["err_delta"] > 0.5].sort_values("err_delta", ascending=False)
        better = err_comp[err_comp["err_delta"] < -0.5].sort_values("err_delta")
        
        print(f"\n  === ERROR DELTA (new - baseline): ===")
        print(f"  Transcripts where error INCREASED: {len(worse)}, total regression = {worse['err_delta'].sum():.0f}")
        print(f"  Transcripts where error DECREASED: {len(better)}, total improvement = {better['err_delta'].sum():.0f}")
        
        print(f"\n  TOP 30 REGRESSIONS (error increased most):")
        print(f"  {'t_id':<25s} {'truth':>8s} {'base':>8s} {'new':>8s} {'b_err':>8s} {'n_err':>8s} {'delta':>8s}")
        print("  " + "-" * 85)
        for _, row in worse.head(30).iterrows():
            print(f"  {row['t_id']:<25s} {row['truth']:>8.1f} {row['base']:>8.1f} {row['new']:>8.1f} "
                  f"{row['base_err']:>8.1f} {row['new_err']:>8.1f} {row['err_delta']:>+8.1f}")
        
        print(f"\n  TOP 30 IMPROVEMENTS (error decreased most):")
        print(f"  {'t_id':<25s} {'truth':>8s} {'base':>8s} {'new':>8s} {'b_err':>8s} {'n_err':>8s} {'delta':>8s}")
        print("  " + "-" * 85)
        for _, row in better.head(30).iterrows():
            print(f"  {row['t_id']:<25s} {row['truth']:>8.1f} {row['base']:>8.1f} {row['new']:>8.1f} "
                  f"{row['base_err']:>8.1f} {row['new_err']:>8.1f} {row['err_delta']:>+8.1f}")
        
        # Statistics on the delta distribution
        print(f"\n  Error delta statistics:")
        print(f"    Mean delta: {err_comp['err_delta'].mean():.4f}")
        print(f"    Median delta: {err_comp['err_delta'].median():.4f}")
        print(f"    P25/P75: {err_comp['err_delta'].quantile(0.25):.4f} / {err_comp['err_delta'].quantile(0.75):.4f}")
        print(f"    Total abs error (baseline): {base_err.sum():.0f}")
        print(f"    Total abs error (new): {abs_err.sum():.0f}")
        print(f"    Net delta: {abs_err.sum() - base_err.sum():.0f}")
        
        # Check: is the increase driven by many small changes or a few large ones?
        n_regression = (err_comp["err_delta"] > 0.5).sum()
        n_improve = (err_comp["err_delta"] < -0.5).sum()
        n_same = len(err_comp) - n_regression - n_improve
        print(f"\n  Change distribution:")
        print(f"    Regression (>0.5): {n_regression} transcripts")
        print(f"    Improved (<-0.5): {n_improve} transcripts")
        print(f"    Same (±0.5): {n_same} transcripts")


if __name__ == "__main__":
    main()
