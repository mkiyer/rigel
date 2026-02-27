#!/usr/bin/env python3
"""Analyze chr8 PVT1/MYC benchmark results."""
import pandas as pd
import numpy as np

# Load per-transcript counts
df = pd.read_csv(
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/per_transcript_counts.csv"
)

# Load abundance file to get negative control info
ab = pd.read_csv(
    "/Users/mkiyer/proj/hulkrna/scripts/benchmarking/chr8_pvt1_myc_abundance.tsv",
    sep="\t",
)
neg_ids = set(ab.loc[ab["is_negative_control"], "transcript_id"])

# Separate negative controls from expressed
df["is_neg"] = df["transcript_id"].isin(neg_ids)
neg = df[df["is_neg"]]
pos = df[~df["is_neg"]]

tools = ["hulkrna_oracle", "hulkrna_minimap2", "salmon", "kallisto"]

print(f"=== NEGATIVE CONTROLS (truth=0, n={len(neg)}) ===")
for t in tools:
    mean_pred = neg[t].mean()
    max_pred = neg[t].max()
    nonzero = (neg[t] > 0.5).sum()
    print(
        f"  {t:25s}  mean={mean_pred:.3f}  max={max_pred:.1f}"
        f"  false_positives(>0.5)={nonzero}"
    )

print()
print(f"=== EXPRESSED TRANSCRIPTS (truth>0, n={len(pos)}) ===")
for t in tools:
    errs = (pos[t] - pos["truth"]).abs()
    mae = errs.mean()
    rmse = np.sqrt((errs**2).mean())
    pearson = pos[["truth", t]].corr().iloc[0, 1]
    dropout = (pos[t] <= 0.5).sum()
    print(
        f"  {t:25s}  MAE={mae:.2f}  RMSE={rmse:.2f}"
        f"  r={pearson:.4f}  dropouts={dropout}"
    )

print()
print("=== TOP-10 WORST TRANSCRIPTS (by hulkrna_oracle absolute error) ===")
pos2 = pos.copy()
pos2["err_oracle"] = (pos2["hulkrna_oracle"] - pos2["truth"]).abs()
worst = pos2.nlargest(10, "err_oracle")
for _, r in worst.iterrows():
    print(
        f"  {r.transcript_id}  truth={r.truth:.1f}"
        f"  hulk_oracle={r.hulkrna_oracle:.1f}"
        f"  salmon={r.salmon:.1f}"
        f"  kallisto={r.kallisto:.1f}"
        f"  err={r.err_oracle:.1f}"
    )

print()
print("=== ABUNDANCE QUARTILE ANALYSIS (expressed only) ===")
pos2["quartile"] = pd.qcut(
    pos2["truth"], 4, labels=["Q1_low", "Q2", "Q3", "Q4_high"]
)
for q in ["Q1_low", "Q2", "Q3", "Q4_high"]:
    qdf = pos2[pos2["quartile"] == q]
    row = f"  {q:10s} n={len(qdf):3d}"
    for t in tools:
        mae_q = (qdf[t] - qdf["truth"]).abs().mean()
        row += f"  {t}={mae_q:.2f}"
    print(row)

# Gene-level negative control analysis
print()
print("=== PER-GENE NEGATIVE CONTROL FALSE POSITIVES (top 10 by hulkrna_oracle FP count) ===")
neg2 = neg.merge(
    ab[["transcript_id", "gene_name"]], on="transcript_id", how="left"
)
for gname, gdf in neg2.groupby("gene_name"):
    fp_oracle = (gdf["hulkrna_oracle"] > 0.5).sum()
    fp_mm2 = (gdf["hulkrna_minimap2"] > 0.5).sum()
    fp_salmon = (gdf["salmon"] > 0.5).sum()
    fp_kall = (gdf["kallisto"] > 0.5).sum()
    max_oracle = gdf["hulkrna_oracle"].max()
    if fp_oracle > 0 or fp_salmon > 0 or max_oracle > 1.0:
        print(
            f"  {gname:15s}  n_neg={len(gdf)}  "
            f"FP: hulk_orc={fp_oracle} hulk_mm2={fp_mm2} "
            f"salmon={fp_salmon} kall={fp_kall}  "
            f"max_hulk_orc={max_oracle:.1f}"
        )
