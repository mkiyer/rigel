#!/usr/bin/env python3
"""Re-parse kallisto output, regenerate per-transcript CSV, and re-analyze."""
import pandas as pd
import numpy as np
import csv
from pathlib import Path

base = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output/gdna_none_ss_1.00_nrna_none")

# ── 1. Re-parse kallisto with ID fix ────────────────────────────
kallisto_tsv = base / "kallisto" / "quant" / "abundance.tsv"
kdf = pd.read_csv(kallisto_tsv, sep="\t")
# Strip pipe-delimited GENCODE headers to just transcript ID
ids = kdf["target_id"].str.split("|").str[0]
kallisto_counts = dict(zip(ids, kdf["est_counts"]))
n_nonzero = sum(1 for v in kallisto_counts.values() if v > 0)
print(f"Kallisto: {len(kallisto_counts)} transcripts, {n_nonzero} non-zero")

# ── 2. Load existing CSV and replace kallisto column ─────────────
csv_path = base / "per_transcript_counts_oracle.csv"
df = pd.read_csv(csv_path)
print(f"Loaded CSV: {df.shape}")
print(f"Old kallisto non-zero: {(df['kallisto'] > 0).sum()}")

# Map kallisto counts by transcript_id
df['kallisto'] = df['transcript_id'].map(kallisto_counts).fillna(0.0)
print(f"New kallisto non-zero: {(df['kallisto'] > 0).sum()}")

# ── 3. Write updated CSV ────────────────────────────────────────
backup_path = csv_path.with_suffix('.csv.bak')
csv_path.rename(backup_path)
df.to_csv(csv_path, index=False)
print(f"Wrote updated CSV ({backup_path.name} backed up)")

# ── 4. Analysis ─────────────────────────────────────────────────
def spearman_r(x, y):
    rx = np.empty_like(x, dtype=float)
    ry = np.empty_like(y, dtype=float)
    rx[np.argsort(x)] = np.arange(len(x), dtype=float)
    ry[np.argsort(y)] = np.arange(len(y), dtype=float)
    d = rx - ry
    n = len(x)
    return 1.0 - 6.0 * np.sum(d**2) / (n * (n**2 - 1))

truth = df['mrna_truth'].astype(float).values
n_expressed = (truth > 0).sum()
print(f"\n{'='*70}")
print(f"Pristine Oracle: {len(truth)} transcripts, {n_expressed} expressed, {truth.sum():.0f} total")
print(f"{'='*70}")

tools = ['rigel_oracle', 'salmon', 'kallisto']

# Two passes: raw and floored (< 1 → 0)
for floor_label, do_floor in [("RAW", False), ("FLOORED (<1 → 0)", True)]:
    print(f"\n{'─'*70}")
    print(f"  {floor_label}")
    print(f"{'─'*70}")

    for tool in tools:
        pred = df[tool].astype(float).values.copy()
        if do_floor:
            pred[pred < 1.0] = 0.0

        total_pred = pred.sum()

        # Union mask
        mask = (truth > 0) | (pred > 0)
        t = truth[mask]
        p = pred[mask]
        n = int(mask.sum())

        mae = np.mean(np.abs(t - p))
        rmse = np.sqrt(np.mean((t - p) ** 2))
        pearson = np.corrcoef(t, p)[0, 1] if t.std() > 0 and p.std() > 0 else 0.0
        sp = spearman_r(t, p) if len(t) > 2 else 0.0

        fp = int(((truth == 0) & (pred > 0)).sum())
        fn = int(((truth > 0) & (pred == 0)).sum())

        # MAPE on expressed
        t_exp = truth[truth > 0]
        p_exp = pred[truth > 0]
        mape = np.mean(np.abs(t_exp - p_exp) / np.maximum(t_exp, 1)) * 100

        print(f"\n  {tool}:")
        print(f"    Total predicted:  {total_pred:>12,.0f}")
        print(f"    N (union):        {n:>12,}")
        print(f"    Pearson:          {pearson:>12.6f}")
        print(f"    Spearman:         {sp:>12.6f}")
        print(f"    MAE:              {mae:>12.2f}")
        print(f"    RMSE:             {rmse:>12.2f}")
        print(f"    MAPE (expressed): {mape:>11.2f}%")
        print(f"    False positives:  {fp:>12,}")
        print(f"    False negatives:  {fn:>12,}")
