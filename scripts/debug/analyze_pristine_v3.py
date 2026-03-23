#!/usr/bin/env python3
"""Deep analysis of pristine benchmark v3 results.

Compares rigel_oracle, rigel_minimap2, salmon, kallisto at transcript level.
Identifies root causes of minimap2 gap and improvement opportunities.
"""
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ── Paths ────────────────────────────────────────────────────────
V2 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v2/gdna_none_ss_0.95_nrna_none")
V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
SIM = Path("/Users/mkiyer/Downloads/rigel_runs/sim_pristine")

# ── Load truth ───────────────────────────────────────────────────
truth_path = SIM / "truth_abundances_nrna_none.tsv"
truth_df = pd.read_csv(truth_path, sep="\t")
print(f"Truth: {len(truth_df)} transcripts, {truth_df.columns.tolist()}")
print(truth_df.head())

# Build truth dict — mrna_abundance is the fragment count column
truth = dict(zip(truth_df["transcript_id"], truth_df["mrna_abundance"]))
n_expressed = sum(1 for v in truth.values() if v > 0)
print(f"\nExpressed transcripts: {n_expressed} / {len(truth)}")

# ── Load tool results ────────────────────────────────────────────
def load_rigel_csv(path, tool_col=None):
    df = pd.read_csv(path)
    # Columns: transcript_id, gene_id, gene_name, mrna_abundance, nrna_abundance,
    #          mrna_truth, nrna_truth, kallisto, rigel_<aligner>, salmon
    if tool_col is None:
        # Auto-detect rigel column
        rigel_cols = [c for c in df.columns if c.startswith("rigel_")]
        tool_col = rigel_cols[0] if rigel_cols else "mrna_abundance"
    return df, tool_col

def load_salmon(path):
    df = pd.read_csv(path / "quant" / "quant.sf", sep="\t")
    return dict(zip(df["Name"], df["NumReads"]))

def load_kallisto(path):
    df = pd.read_csv(path / "quant" / "abundance.tsv", sep="\t")
    ids = df["target_id"].str.split("|").str[0]
    return dict(zip(ids, df["est_counts"]))

orc_df, orc_col = load_rigel_csv(V3 / "per_transcript_counts_oracle.csv")
mm3_df, mm3_col = load_rigel_csv(V3 / "per_transcript_counts_minimap2.csv")
mm2_df, mm2_col = load_rigel_csv(V2 / "per_transcript_counts_minimap2.csv")

print(f"Rigel CSV columns: {orc_df.columns.tolist()}")
print(f"Oracle col: {orc_col}, MM2 v3 col: {mm3_col}, MM2 v2 col: {mm2_col}")

# Build dicts from the per-transcript CSV (which has all tools)
oracle_v3 = dict(zip(orc_df["transcript_id"], orc_df[orc_col]))
mm2_v3 = dict(zip(mm3_df["transcript_id"], mm3_df[mm3_col]))
mm2_v2 = dict(zip(mm2_df["transcript_id"], mm2_df[mm2_col]))

# Salmon and kallisto are embedded in the oracle CSV
if "salmon" in orc_df.columns:
    salmon = dict(zip(orc_df["transcript_id"], orc_df["salmon"]))
else:
    salmon = load_salmon(V3)
if "kallisto" in orc_df.columns:
    kallisto = dict(zip(orc_df["transcript_id"], orc_df["kallisto"]))
else:
    kallisto = load_kallisto(V3)

print(f"\nTool transcript counts:")
for name, d in [("oracle_v3", oracle_v3), ("mm2_v3", mm2_v3), ("mm2_v2", mm2_v2),
                ("salmon", salmon), ("kallisto", kallisto)]:
    print(f"  {name}: {len(d)} transcripts, total={sum(d.values()):.0f}")

# ── Unified dataframe ────────────────────────────────────────────
# Truth comes from the benchmark CSV (mrna_truth column = fragment counts from FASTQ)
truth_from_csv = dict(zip(orc_df["transcript_id"], orc_df["mrna_truth"]))
all_tids = sorted(set(truth_from_csv) | set(oracle_v3) | set(mm2_v3) | set(salmon) | set(kallisto))
df = pd.DataFrame({
    "tid": all_tids,
    "truth": [truth_from_csv.get(t, 0) for t in all_tids],
    "oracle": [oracle_v3.get(t, 0) for t in all_tids],
    "mm2_v3": [mm2_v3.get(t, 0) for t in all_tids],
    "mm2_v2": [mm2_v2.get(t, 0) for t in all_tids],
    "salmon": [salmon.get(t, 0) for t in all_tids],
    "kallisto": [kallisto.get(t, 0) for t in all_tids],
})

# Compute errors
for tool in ["oracle", "mm2_v3", "mm2_v2", "salmon", "kallisto"]:
    df[f"{tool}_err"] = df[tool] - df["truth"]
    df[f"{tool}_abs_err"] = df[f"{tool}_err"].abs()
    df[f"{tool}_rel_err"] = df[f"{tool}_err"] / df["truth"].clip(lower=1)

print(f"\nUnified DF: {len(df)} transcripts")

# ══════════════════════════════════════════════════════════════════
# 1. TOP-LEVEL COMPARISON
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("1. TOP-LEVEL COMPARISON (v2 pre-fix vs v3 post-fix)")
print("="*80)

metrics = {}
for tool in ["oracle", "mm2_v3", "mm2_v2", "salmon", "kallisto"]:
    t_arr = df["truth"].values.astype(float)
    o_arr = df[tool].values.astype(float)
    abs_err = np.abs(t_arr - o_arr)
    mae = abs_err.mean()
    rmse = np.sqrt(np.mean((t_arr - o_arr)**2))
    total = o_arr.sum()
    total_abs = abs_err.sum()
    metrics[tool] = {"MAE": mae, "RMSE": rmse, "total": total, "total_abs_err": total_abs}

print(f"\n{'Tool':<15} {'MAE':>8} {'RMSE':>10} {'Total Predicted':>16} {'Total |Err|':>14}")
print("-"*65)
for tool in ["oracle", "mm2_v3", "mm2_v2", "salmon", "kallisto"]:
    m = metrics[tool]
    print(f"{tool:<15} {m['MAE']:>8.2f} {m['RMSE']:>10.2f} {m['total']:>16,.0f} {m['total_abs_err']:>14,.0f}")

print(f"\nv2→v3 improvement:")
print(f"  MAE: {metrics['mm2_v2']['MAE']:.2f} → {metrics['mm2_v3']['MAE']:.2f} "
      f"({(1 - metrics['mm2_v3']['MAE']/metrics['mm2_v2']['MAE'])*100:.1f}% reduction)")
print(f"  RMSE: {metrics['mm2_v2']['RMSE']:.2f} → {metrics['mm2_v3']['RMSE']:.2f} "
      f"({(1 - metrics['mm2_v3']['RMSE']/metrics['mm2_v2']['RMSE'])*100:.1f}% reduction)")

# ══════════════════════════════════════════════════════════════════
# 2. POOL-LEVEL ANALYSIS (nRNA siphon, gDNA leakage)
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("2. POOL-LEVEL ANALYSIS")
print("="*80)

# Load nRNA counts from rigel
# Use already-loaded dataframes
print(f"\n{'Metric':<30} {'Oracle v3':>12} {'MM2 v2':>12} {'MM2 v3':>12} {'Improvement':>12}")
print("-"*80)
nrna_col = "nrna_abundance"
oracle_nrna_total = orc_df[nrna_col].sum() if nrna_col in orc_df.columns else 0
mm2_v2_nrna_total = mm2_df[nrna_col].sum() if nrna_col in mm2_df.columns else 0
mm2_v3_nrna_total = mm3_df[nrna_col].sum() if nrna_col in mm3_df.columns else 0
print(f"{'nRNA siphon (should be 0)':<30} {oracle_nrna_total:>12,.0f} {mm2_v2_nrna_total:>12,.0f} {mm2_v3_nrna_total:>12,.0f} {(1-mm2_v3_nrna_total/max(mm2_v2_nrna_total,1))*100:>11.1f}%")

mrna_col = "mrna_abundance"
oracle_mrna_total = orc_df[mrna_col].sum()
mm2_v2_mrna_total = mm2_df[mrna_col].sum()
mm2_v3_mrna_total = mm3_df[mrna_col].sum()
print(f"{'mRNA recovered':<30} {oracle_mrna_total:>12,.0f} {mm2_v2_mrna_total:>12,.0f} {mm2_v3_mrna_total:>12,.0f} {mm2_v3_mrna_total-mm2_v2_mrna_total:>+12,.0f}")

# ══════════════════════════════════════════════════════════════════
# 3. ERROR DISTRIBUTION BY EXPRESSION LEVEL
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("3. ERROR DISTRIBUTION BY EXPRESSION LEVEL")
print("="*80)

# Bin transcripts by truth count
bins = [0, 1, 10, 100, 1000, 10000, float("inf")]
labels = ["0", "1-9", "10-99", "100-999", "1K-10K", ">10K"]
df["expr_bin"] = pd.cut(df["truth"], bins=bins, labels=labels, right=False)

print(f"\n{'Bin':<10} {'N_tx':>8}", end="")
for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    print(f"  {tool+' MAE':>14}", end="")
print()
print("-"*80)

for label in labels:
    subset = df[df["expr_bin"] == label]
    n = len(subset)
    print(f"{label:<10} {n:>8}", end="")
    for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
        mae = subset[f"{tool}_abs_err"].mean()
        print(f"  {mae:>14.2f}", end="")
    print()

# ══════════════════════════════════════════════════════════════════
# 4. WORST TRANSCRIPTS FOR MINIMAP2 v3
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("4. TOP 30 WORST TRANSCRIPTS FOR rigel_minimap2 v3")
print("="*80)

worst = df.nlargest(30, "mm2_v3_abs_err")
print(f"\n{'Transcript':<25} {'Truth':>8} {'MM2v3':>8} {'Oracle':>8} {'Salmon':>8} {'Kalli':>8} {'MM2 Err':>9} {'Orc Err':>9}")
print("-"*100)
for _, row in worst.iterrows():
    print(f"{row['tid']:<25} {row['truth']:>8.0f} {row['mm2_v3']:>8.0f} {row['oracle']:>8.0f} "
          f"{row['salmon']:>8.0f} {row['kallisto']:>8.0f} {row['mm2_v3_err']:>+9.0f} {row['oracle_err']:>+9.0f}")

# ══════════════════════════════════════════════════════════════════
# 5. IMPROVEMENT ANALYSIS: v2 → v3
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("5. IMPROVEMENT ANALYSIS: v2 → v3 (per-transcript)")
print("="*80)

df["improvement"] = df["mm2_v2_abs_err"] - df["mm2_v3_abs_err"]
improved = df[df["improvement"] > 0]
degraded = df[df["improvement"] < 0]
unchanged = df[df["improvement"] == 0]

print(f"\n  Transcripts improved:  {len(improved):>8} (total improvement: {improved['improvement'].sum():>+12,.0f})")
print(f"  Transcripts degraded:  {len(degraded):>8} (total degradation: {degraded['improvement'].sum():>+12,.0f})")
print(f"  Transcripts unchanged: {len(unchanged):>8}")
print(f"  Net improvement:                            {df['improvement'].sum():>+12,.0f}")

print(f"\nTop 20 most improved transcripts:")
most_improved = df.nlargest(20, "improvement")
print(f"{'Transcript':<25} {'Truth':>8} {'v2 Err':>9} {'v3 Err':>9} {'Δ':>9}")
print("-"*65)
for _, row in most_improved.iterrows():
    print(f"{row['tid']:<25} {row['truth']:>8.0f} {row['mm2_v2_abs_err']:>9.0f} {row['mm2_v3_abs_err']:>9.0f} {row['improvement']:>+9.0f}")

print(f"\nTop 20 most degraded transcripts:")
most_degraded = df.nsmallest(20, "improvement")
print(f"{'Transcript':<25} {'Truth':>8} {'v2 Err':>9} {'v3 Err':>9} {'Δ':>9}")
print("-"*65)
for _, row in most_degraded.iterrows():
    print(f"{row['tid']:<25} {row['truth']:>8.0f} {row['mm2_v2_abs_err']:>9.0f} {row['mm2_v3_abs_err']:>9.0f} {row['improvement']:>+9.0f}")

# ══════════════════════════════════════════════════════════════════
# 6. GAP ANALYSIS: MM2 v3 vs ORACLE
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("6. GAP ANALYSIS: rigel_minimap2 v3 vs oracle")
print("="*80)

# For each transcript, how much of the mm2 error is also present in oracle?
df["mm2_excess"] = df["mm2_v3_abs_err"] - df["oracle_abs_err"]
df["mm2_excess_positive"] = df["mm2_excess"].clip(lower=0)

total_mm2_err = df["mm2_v3_abs_err"].sum()
total_orc_err = df["oracle_abs_err"].sum()
total_excess = df["mm2_excess_positive"].sum()

print(f"\n  Total |error| in mm2 v3:     {total_mm2_err:>12,.0f}")
print(f"  Total |error| in oracle:     {total_orc_err:>12,.0f}")
print(f"  Excess mm2 error over oracle:{total_excess:>12,.0f} ({total_excess/total_mm2_err*100:.1f}% of mm2 error)")

# Classify: is the error direction same or opposite?
both_positive = ((df["mm2_v3_err"] > 0) & (df["oracle_err"] > 0)).sum()
both_negative = ((df["mm2_v3_err"] < 0) & (df["oracle_err"] < 0)).sum()
opposite_sign = ((df["mm2_v3_err"] * df["oracle_err"]) < 0).sum()
both_zero = ((df["mm2_v3_err"] == 0) & (df["oracle_err"] == 0)).sum()
print(f"\n  Error direction (mm2 vs oracle):")
print(f"    Same direction (both over):   {both_positive:>8}")
print(f"    Same direction (both under):  {both_negative:>8}")
print(f"    Opposite sign:                {opposite_sign:>8}")
print(f"    Both exact:                   {both_zero:>8}")

# ══════════════════════════════════════════════════════════════════
# 7. MINIMAP2 ERROR DECOMPOSITION
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("7. MINIMAP2 ERROR DECOMPOSITION")
print("="*80)

# Transcripts where mm2 has error > 10 but oracle doesn't
mm2_specific = df[(df["mm2_v3_abs_err"] > 10) & (df["oracle_abs_err"] <= 2)]
print(f"\nTranscripts with mm2 |error|>10 but oracle |error|≤2: {len(mm2_specific)}")

# Are these over or under estimated?
n_over = (mm2_specific["mm2_v3_err"] > 0).sum()
n_under = (mm2_specific["mm2_v3_err"] < 0).sum()
total_over = mm2_specific[mm2_specific["mm2_v3_err"] > 0]["mm2_v3_err"].sum()
total_under = mm2_specific[mm2_specific["mm2_v3_err"] < 0]["mm2_v3_err"].sum()
print(f"  Over-estimated:  {n_over:>6} transcripts, total excess: {total_over:>+10,.0f}")
print(f"  Under-estimated: {n_under:>6} transcripts, total deficit:{total_under:>+10,.0f}")

# What about the salmon/kallisto error at these same transcripts?
print(f"\n  At these same transcripts:")
print(f"    Salmon MAE:   {mm2_specific['salmon_abs_err'].mean():.2f}")
print(f"    Kallisto MAE: {mm2_specific['kallisto_abs_err'].mean():.2f}")
print(f"    MM2 v3 MAE:   {mm2_specific['mm2_v3_abs_err'].mean():.2f}")

# ══════════════════════════════════════════════════════════════════
# 8. SALMON vs RIGEL COMPARISON
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("8. HEAD-TO-HEAD: RIGEL ORACLE vs SALMON")
print("="*80)

# Transcripts where rigel oracle beats salmon
rigel_better = df[df["oracle_abs_err"] < df["salmon_abs_err"]]
salmon_better = df[df["salmon_abs_err"] < df["oracle_abs_err"]]
tied = df[df["oracle_abs_err"] == df["salmon_abs_err"]]
print(f"\n  Rigel oracle wins: {len(rigel_better):>8} transcripts")
print(f"  Salmon wins:       {len(salmon_better):>8} transcripts")
print(f"  Tied:              {len(tied):>8} transcripts")
print(f"\n  MAE where rigel wins: rigel={rigel_better['oracle_abs_err'].mean():.3f}, salmon={rigel_better['salmon_abs_err'].mean():.3f}")
if len(salmon_better) > 0:
    print(f"  MAE where salmon wins: rigel={salmon_better['oracle_abs_err'].mean():.3f}, salmon={salmon_better['salmon_abs_err'].mean():.3f}")

# Total error margin
print(f"\n  Total |error| — rigel oracle: {df['oracle_abs_err'].sum():>12,.0f}")
print(f"  Total |error| — salmon:       {df['salmon_abs_err'].sum():>12,.0f}")
print(f"  Rigel advantage:              {df['salmon_abs_err'].sum() - df['oracle_abs_err'].sum():>+12,.0f}")

# Gene-level comparison
# Need gene info — load from rigel output
gid_map = dict(zip(orc_df["transcript_id"], orc_df["gene_id"]))
df["gid"] = df["tid"].map(gid_map).fillna("unknown")

gene_df = df.groupby("gid").agg(
    truth=("truth", "sum"),
    oracle=("oracle", "sum"),
    mm2_v3=("mm2_v3", "sum"),
    salmon=("salmon", "sum"),
    kallisto=("kallisto", "sum"),
).reset_index()
for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    gene_df[f"{tool}_abs_err"] = (gene_df[tool] - gene_df["truth"]).abs()

print(f"\n  Gene-level MAE:")
for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    print(f"    {tool:<15} {gene_df[f'{tool}_abs_err'].mean():.3f}")

# ══════════════════════════════════════════════════════════════════
# 9. FRACTION OF TOTAL ERROR BY TOP TX
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("9. ERROR CONCENTRATION")
print("="*80)

for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    sorted_err = df[f"{tool}_abs_err"].sort_values(ascending=False)
    total = sorted_err.sum()
    top10 = sorted_err.head(10).sum()
    top50 = sorted_err.head(50).sum()
    top100 = sorted_err.head(100).sum()
    top500 = sorted_err.head(500).sum()
    print(f"\n  {tool}:")
    print(f"    Top 10 transcripts:  {top10:>10,.0f} ({top10/total*100:.1f}% of total error)")
    print(f"    Top 50 transcripts:  {top50:>10,.0f} ({top50/total*100:.1f}%)")
    print(f"    Top 100 transcripts: {top100:>10,.0f} ({top100/total*100:.1f}%)")
    print(f"    Top 500 transcripts: {top500:>10,.0f} ({top500/total*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════
# 10. QUANTILE ANALYSIS
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("10. ERROR QUANTILES (absolute error)")
print("="*80)

quantiles = [0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 1.0]
print(f"\n{'Quantile':<10}", end="")
for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    print(f"  {tool:>12}", end="")
print()
print("-"*60)
for q in quantiles:
    print(f"{q:<10.3f}", end="")
    for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
        val = df[f"{tool}_abs_err"].quantile(q)
        print(f"  {val:>12.1f}", end="")
    print()

# ══════════════════════════════════════════════════════════════════
# 11. SALMON "MISSING" COUNTS ANALYSIS
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("11. TOTAL COUNT RECONCILIATION")
print("="*80)

print(f"\n  Truth total:           {df['truth'].sum():>12,.0f}")
print(f"  Rigel oracle total:    {df['oracle'].sum():>12,.0f} (diff: {df['oracle'].sum() - df['truth'].sum():>+10,.0f})")
print(f"  Rigel mm2 v3 total:    {df['mm2_v3'].sum():>12,.0f} (diff: {df['mm2_v3'].sum() - df['truth'].sum():>+10,.0f})")
print(f"  Salmon total:          {df['salmon'].sum():>12,.0f} (diff: {df['salmon'].sum() - df['truth'].sum():>+10,.0f})")
print(f"  Kallisto total:        {df['kallisto'].sum():>12,.0f} (diff: {df['kallisto'].sum() - df['truth'].sum():>+10,.0f})")

# ══════════════════════════════════════════════════════════════════
# 12. SHARED HIGH-ERROR TRANSCRIPTS
# ══════════════════════════════════════════════════════════════════
print("\n" + "="*80)
print("12. SHARED vs TOOL-SPECIFIC HIGH-ERROR TRANSCRIPTS (|err|>50)")
print("="*80)

for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    df[f"{tool}_hi"] = df[f"{tool}_abs_err"] > 50

n_any = df[["oracle_hi", "mm2_v3_hi", "salmon_hi", "kallisto_hi"]].any(axis=1).sum()
n_all = df[["oracle_hi", "mm2_v3_hi", "salmon_hi", "kallisto_hi"]].all(axis=1).sum()
print(f"\n  Transcripts with |err|>50 in ANY tool:  {n_any}")
print(f"  Transcripts with |err|>50 in ALL tools: {n_all}")

for tool in ["oracle", "mm2_v3", "salmon", "kallisto"]:
    n = df[f"{tool}_hi"].sum()
    unique = (df[f"{tool}_hi"] & ~df[[f"{t}_hi" for t in ["oracle", "mm2_v3", "salmon", "kallisto"] if t != tool]].any(axis=1)).sum()
    print(f"  {tool}: {n} high-error tx, {unique} unique to this tool")

print("\n\n=== ANALYSIS COMPLETE ===")
