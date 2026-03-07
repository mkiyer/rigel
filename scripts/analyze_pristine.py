#!/usr/bin/env python3
"""Analyze pristine benchmark results."""
import pandas as pd
import numpy as np

def spearman_r(x, y):
    rx = np.empty_like(x)
    ry = np.empty_like(y)
    rx[np.argsort(x)] = np.arange(len(x), dtype=float)
    ry[np.argsort(y)] = np.arange(len(y), dtype=float)
    d = rx - ry
    n = len(x)
    return 1.0 - 6.0 * np.sum(d**2) / (n * (n**2 - 1))

def analyze(path, label):
    df = pd.read_csv(path)
    print(f"Columns: {list(df.columns)}")
    print(f"Shape: {df.shape}")
    print()

    tools = [c for c in df.columns if c not in
             ('transcript_id', 'gene_id', 'gene_name',
              'mrna_abundance', 'nrna_abundance', 'mrna_truth', 'nrna_truth')]
    print(f"Tools: {tools}")
    print()

    truth = df['mrna_truth'].astype(float)
    n_expressed = (truth > 0).sum()
    n_total = len(truth)
    print(f"=== {label} ===")
    print(f"Transcripts: {n_total}")
    print(f"Expressed (truth > 0): {n_expressed}")
    print(f"Total truth fragments: {truth.sum():.0f}")
    print()

    for tool in tools:
        pred = df[tool].astype(float)
        total_pred = pred.sum()

        # Filter to union of expressed
        mask = (truth > 0) | (pred > 0)
        t = truth[mask].values
        p = pred[mask].values
        n = int(mask.sum())

        mae = np.mean(np.abs(t - p))
        rmse = np.sqrt(np.mean((t - p) ** 2))

        pearson = np.corrcoef(t, p)[0, 1] if t.std() > 0 and p.std() > 0 else 0.0
        spearman = spearman_r(t, p) if len(t) > 2 else 0.0

        # MAPE on expressed only
        expressed = truth > 0
        t_exp = truth[expressed].values
        p_exp = pred[expressed].values
        mape = np.mean(np.abs(t_exp - p_exp) / np.maximum(t_exp, 1)) * 100

        # False positive / false negative counts
        fp = int(((truth == 0) & (pred > 0)).sum())
        fn = int(((truth > 0) & (pred == 0)).sum())

        print(f"{tool}:")
        print(f"  Total predicted:      {total_pred:.0f}")
        print(f"  N transcripts (union): {n}")
        print(f"  MAE:                  {mae:.2f}")
        print(f"  RMSE:                 {rmse:.2f}")
        print(f"  Pearson:              {pearson:.6f}")
        print(f"  Spearman:             {spearman:.6f}")
        print(f"  MAPE (expressed):     {mape:.2f}%")
        print(f"  False positives:      {fp}")
        print(f"  False negatives:      {fn}")
        print()

base = "/Users/mkiyer/Downloads/rigel_runs/benchmark_output/gdna_none_ss_1.00_nrna_none"

print("=" * 60)
analyze(f"{base}/per_transcript_counts_oracle.csv", "Pristine (Oracle BAM)")

import os
minimap2_path = f"{base}/per_transcript_counts_minimap2.csv"
if os.path.exists(minimap2_path):
    print("=" * 60)
    analyze(minimap2_path, "Pristine (Minimap2 BAM)")
