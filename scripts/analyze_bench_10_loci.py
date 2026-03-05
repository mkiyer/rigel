#!/usr/bin/env python3
"""Analyze benchmark_10_loci results and produce a diagnostic report."""

import pandas as pd
import numpy as np
import sys
import json
from pathlib import Path

CSV = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_10_loci/summary.csv")
JSON = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_10_loci/summary.json")
DIAG = Path("/Users/mkiyer/Downloads/hulkrna_runs/bench_10_loci/diagnostics.md")

df = pd.read_csv(CSV)

# Separate transcript/gene level metrics from pool-level
tx = df[df["level"] == "transcript"].copy()
gene = df[df["level"] == "gene"].copy()
pool = df[df["level"] == "pool"].copy()

print("=" * 80)
print("BENCHMARK 10-LOCI RESULTS ANALYSIS")
print("=" * 80)

print(f"\nTotal rows: {len(df)}")
print(f"Transcript rows: {len(tx)}, Gene rows: {len(gene)}, Pool rows: {len(pool)}")
print(f"Regions ({tx['region'].nunique()}): {sorted(tx['region'].unique())}")
print(f"Tools: {sorted(tx['tool'].unique())}")
print(f"Strand specs: {sorted(tx['strand_specificity'].unique())}")
print(f"gDNA labels: {sorted(tx['gdna_label'].unique())}")
print(f"nRNA labels: {sorted(tx['nrna_label'].unique())}")

# ── 1. TRANSCRIPT-LEVEL OVERALL ──────────────────────────────────────────
print("\n" + "=" * 80)
print("1. TRANSCRIPT-LEVEL: OVERALL TOOL PERFORMANCE")
print("=" * 80)
agg = (
    tx.groupby("tool")
    .agg(
        n=("pearson", "count"),
        pearson=("pearson", "mean"),
        spearman=("spearman", "mean"),
        rmse=("rmse", "mean"),
        mae=("mean_abs_error", "mean"),
    )
    .round(4)
)
print(agg.to_string())

# ── 2. GENE-LEVEL OVERALL ────────────────────────────────────────────────
print("\n" + "=" * 80)
print("2. GENE-LEVEL: OVERALL TOOL PERFORMANCE")
print("=" * 80)
agg_g = (
    gene.groupby("tool")
    .agg(
        n=("pearson", "count"),
        pearson=("pearson", "mean"),
        spearman=("spearman", "mean"),
        rmse=("rmse", "mean"),
        mae=("mean_abs_error", "mean"),
    )
    .round(4)
)
print(agg_g.to_string())

# ── 3. TRANSCRIPT-LEVEL BY CONDITION (gDNA × nRNA × strand) ─────────────
print("\n" + "=" * 80)
print("3. TRANSCRIPT-LEVEL: TOOL × gDNA × nRNA × STRAND")
print("=" * 80)
cond = (
    tx.groupby(["tool", "gdna_label", "nrna_label", "strand_specificity"])
    .agg(
        pearson=("pearson", "mean"),
        spearman=("spearman", "mean"),
        rmse=("rmse", "mean"),
    )
    .round(4)
)
print(cond.to_string())

# ── 4. TRANSCRIPT-LEVEL BY REGION ────────────────────────────────────────
print("\n" + "=" * 80)
print("4. TRANSCRIPT-LEVEL: TOOL × REGION (averaged over conditions)")
print("=" * 80)
region_tool = (
    tx.groupby(["region", "tool"])
    .agg(
        pearson=("pearson", "mean"),
        spearman=("spearman", "mean"),
        rmse=("rmse", "mean"),
    )
    .round(4)
)
print(region_tool.to_string())

# ── 5. GENE-LEVEL BY REGION ─────────────────────────────────────────────
print("\n" + "=" * 80)
print("5. GENE-LEVEL: TOOL × REGION (averaged over conditions)")
print("=" * 80)
gene_region_tool = (
    gene.groupby(["region", "tool"])
    .agg(
        pearson=("pearson", "mean"),
        spearman=("spearman", "mean"),
        rmse=("rmse", "mean"),
    )
    .round(4)
)
print(gene_region_tool.to_string())

# ── 6. POOL-LEVEL ANALYSIS (mRNA vs gDNA classification) ────────────────
print("\n" + "=" * 80)
print("6. POOL-LEVEL: SIGNED ERROR BY TOOL × POOL × gDNA CONDITION")
print("=" * 80)
pool_agg = (
    pool.groupby(["tool", "pool", "gdna_label", "nrna_label"])
    .agg(
        n=("signed_error", "count"),
        mean_signed_err=("signed_error", "mean"),
        mean_abs_err=("signed_error", lambda x: x.abs().mean()),
        truth_mean=("truth", "mean"),
        obs_mean=("observed", "mean"),
    )
    .round(2)
)
print(pool_agg.to_string())

# ── 7. WORST TRANSCRIPT-LEVEL CASES ─────────────────────────────────────
print("\n" + "=" * 80)
print("7. WORST 20 TRANSCRIPT-LEVEL SPEARMAN SCORES (hulkrna)")
print("=" * 80)
hulk_tx = tx[tx["tool"].str.startswith("hulkrna")].copy()
worst = hulk_tx.nsmallest(20, "spearman")[
    [
        "region",
        "tool",
        "gdna_label",
        "nrna_label",
        "strand_specificity",
        "spearman",
        "pearson",
        "rmse",
    ]
]
print(worst.to_string(index=False))

# ── 8. HULKRNA vs BEST COMPETITOR: HEAD-TO-HEAD ─────────────────────────
print("\n" + "=" * 80)
print("8. HEAD-TO-HEAD: hulkrna_oracle vs salmon vs kallisto (transcript)")
print("=" * 80)

# Define condition key
tx["cond_key"] = (
    tx["region"]
    + "|"
    + tx["gdna_label"]
    + "|"
    + tx["nrna_label"]
    + "|"
    + tx["strand_specificity"].astype(str)
)

# Pivot by tool
pvt = tx.pivot_table(
    index="cond_key", columns="tool", values="spearman", aggfunc="first"
)
tools = [c for c in pvt.columns if c in ("hulkrna_oracle", "hulkrna_minimap2", "salmon", "kallisto")]
pvt_sub = pvt[tools].dropna()

for t in tools:
    if t == "hulkrna_oracle":
        continue
    if "hulkrna_oracle" in pvt_sub.columns:
        wins = (pvt_sub["hulkrna_oracle"] > pvt_sub[t]).sum()
        losses = (pvt_sub["hulkrna_oracle"] < pvt_sub[t]).sum()
        ties = (pvt_sub["hulkrna_oracle"] == pvt_sub[t]).sum()
        print(f"hulkrna_oracle vs {t}: wins={wins}, losses={losses}, ties={ties}")

print()
for t in tools:
    if t == "hulkrna_minimap2":
        continue
    if "hulkrna_minimap2" in pvt_sub.columns:
        wins = (pvt_sub["hulkrna_minimap2"] > pvt_sub[t]).sum()
        losses = (pvt_sub["hulkrna_minimap2"] < pvt_sub[t]).sum()
        ties = (pvt_sub["hulkrna_minimap2"] == pvt_sub[t]).sum()
        print(f"hulkrna_minimap2 vs {t}: wins={wins}, losses={losses}, ties={ties}")

# ── 9. CONDITIONS WHERE HULKRNA LOSES TO SALMON OR KALLISTO ──────────────
print("\n" + "=" * 80)
print("9. CONDITIONS WHERE hulkrna_oracle LOSES TO salmon OR kallisto")
print("=" * 80)
for competitor in ["salmon", "kallisto"]:
    if competitor not in pvt_sub.columns or "hulkrna_oracle" not in pvt_sub.columns:
        continue
    losses_mask = pvt_sub["hulkrna_oracle"] < pvt_sub[competitor]
    if losses_mask.sum() == 0:
        print(f"  hulkrna_oracle never loses to {competitor}")
        continue
    loss_rows = pvt_sub[losses_mask][["hulkrna_oracle", competitor]].copy()
    loss_rows["gap"] = loss_rows[competitor] - loss_rows["hulkrna_oracle"]
    loss_rows = loss_rows.sort_values("gap", ascending=False)
    print(f"\nhulkrna_oracle loses to {competitor} ({losses_mask.sum()} conditions):")
    print(loss_rows.head(20).to_string())

# ── 10. POOL ERROR PATTERNS PER REGION ───────────────────────────────────
print("\n" + "=" * 80)
print("10. POOL-LEVEL: mRNA signed error per tool × region (gDNA=high only)")
print("=" * 80)
pool_hi = pool[(pool["gdna_label"] == "high") & (pool["pool"] == "mRNA")]
if len(pool_hi) > 0:
    pool_region = (
        pool_hi.groupby(["region", "tool"])
        .agg(
            mean_signed_err=("signed_error", "mean"),
        )
        .round(1)
    )
    print(pool_region.to_string())
else:
    print("No high-gDNA mRNA pool data found")

# ── 11. WORST CONDITIONS: hulkrna_minimap2 ──────────────────────────────
print("\n" + "=" * 80)
print("11. WORST 20 hulkrna_minimap2 TRANSCRIPT SPEARMAN SCORES")
print("=" * 80)
hulk_mm = tx[tx["tool"] == "hulkrna_minimap2"].copy()
worst_mm = hulk_mm.nsmallest(20, "spearman")[
    [
        "region",
        "gdna_label",
        "nrna_label",
        "strand_specificity",
        "spearman",
        "pearson",
        "rmse",
        "n_transcripts",
    ]
]
print(worst_mm.to_string(index=False))

# ── 12. REGION DIFFICULTY RANKING ────────────────────────────────────────
print("\n" + "=" * 80)
print("12. REGION DIFFICULTY (mean spearman across ALL tools)")
print("=" * 80)
diff = (
    tx.groupby("region")
    .agg(
        mean_spearman=("spearman", "mean"),
        mean_pearson=("pearson", "mean"),
        n_tx=("n_transcripts", "first"),
    )
    .sort_values("mean_spearman")
    .round(4)
)
print(diff.to_string())

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
