#!/usr/bin/env python3
"""Find transcripts where oracle+salmon/kallisto are accurate but rigel-minimap2 fails.

Goal: identify cases where the problem is specifically minimap2 alignment quality
(not a fundamental identifiability problem), so rigel can be improved.

Accuracy criterion for a tool:
  - err_ratio = |obs - truth| / max(truth, 1)
  - "accurate" = err_ratio <= ACCURATE_THRESHOLD
  - "inaccurate" = err_ratio > INACCURATE_THRESHOLD  (mm2 must be THIS bad)

Only consider transcripts with truth >= MIN_TRUTH_COUNT to avoid noisy low-count cases.
"""
import pandas as pd
import numpy as np

COND = "/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine/gdna_none_ss_0.95_nrna_none"

# Thresholds
MIN_TRUTH_COUNT = 50        # ignore very-low-expression transcripts
ACCURATE_THRESHOLD = 0.25   # oracle/salmon/kallisto must be within 25% of truth
MM2_INACCURATE_THRESHOLD = 0.5  # mm2 must deviate by at least 50%

# ── Load data ────────────────────────────────────────────────────
oracle_df = pd.read_csv(f"{COND}/per_transcript_counts_oracle.csv")
mm2_df    = pd.read_csv(f"{COND}/per_transcript_counts_minimap2.csv")

df = oracle_df.merge(mm2_df[["transcript_id", "rigel_minimap2"]], on="transcript_id", how="left")
df["rigel_minimap2"] = df["rigel_minimap2"].fillna(0)

# ── Compute error ratios ─────────────────────────────────────────
truth = df["mrna_truth"].values.astype(float)
denom = np.maximum(truth, 1.0)

df["err_oracle"]  = np.abs(df["rigel_oracle"]  - truth) / denom
df["err_kallisto"] = np.abs(df["kallisto"]       - truth) / denom
df["err_salmon"]  = np.abs(df["salmon"]         - truth) / denom
df["err_mm2"]     = np.abs(df["rigel_minimap2"] - truth) / denom

# Absolute error (counts)
df["abs_err_mm2"] = np.abs(df["rigel_minimap2"] - truth)

# ── Filter: exclude very-low-truth transcripts ──────────────────
df_hi = df[df["mrna_truth"] >= MIN_TRUTH_COUNT].copy()

# ── Apply truth criterion ─────────────────────────────────────────
oracle_good     = df_hi["err_oracle"]   <= ACCURATE_THRESHOLD
salmon_good     = df_hi["err_salmon"]   <= ACCURATE_THRESHOLD
kallisto_good   = df_hi["err_kallisto"] <= ACCURATE_THRESHOLD
mm2_bad         = df_hi["err_mm2"]      >= MM2_INACCURATE_THRESHOLD

# Main filter: oracle accurate AND (salmon OR kallisto) accurate AND mm2 bad
target = df_hi[oracle_good & (salmon_good | kallisto_good) & mm2_bad].copy()
target = target.sort_values("abs_err_mm2", ascending=False)

print(f"Filter criteria:")
print(f"  truth >= {MIN_TRUTH_COUNT}")
print(f"  oracle err_ratio <= {ACCURATE_THRESHOLD}  (within {int(100*ACCURATE_THRESHOLD)}% of truth)")
print(f"  salmon OR kallisto err_ratio <= {ACCURATE_THRESHOLD}")
print(f"  rigel_mm2 err_ratio >= {MM2_INACCURATE_THRESHOLD}  (off by at least {int(100*MM2_INACCURATE_THRESHOLD)}%)")
print(f"\nTotal transcripts with truth >= {MIN_TRUTH_COUNT}: {len(df_hi):,}")
print(f"Transcripts matching ALL criteria: {len(target):,}")

if len(target) == 0:
    print("\nNo transcripts found. Relaxing thresholds...")
    for mm2_thresh in [0.4, 0.3, 0.2]:
        t2 = df_hi[oracle_good & (salmon_good | kallisto_good) & (df_hi["err_mm2"] >= mm2_thresh)]
        print(f"  mm2 err >= {mm2_thresh}: {len(t2)} transcripts")

print("\n" + "=" * 100)
print(f"TOP 20 TRANSCRIPTS (oracle+salmon/kallisto accurate, mm2 bad)")
print("=" * 100)

cols = ["transcript_id", "gene_id", "gene_name", "mrna_truth",
        "rigel_oracle", "salmon", "kallisto", "rigel_minimap2",
        "err_oracle", "err_salmon", "err_kallisto", "err_mm2", "abs_err_mm2"]
print(f"\n  {'transcript_id':42s}  {'gene_name':18s}  {'truth':>7}  {'oracle':>7}  {'salmon':>7}  "
      f"{'kalli':>7}  {'mm2':>7}  {'err_ora':>7}  {'err_sal':>7}  {'err_mm2':>7}  abs_err")
print("  " + "-" * 115)
for _, r in target.head(20).iterrows():
    print(
        f"  {r.transcript_id:42s}  {str(r.gene_name):18s}  {r.mrna_truth:7.0f}  "
        f"{r.rigel_oracle:7.0f}  {r.salmon:7.0f}  {r.kallisto:7.0f}  "
        f"{r.rigel_minimap2:7.0f}  {r.err_oracle:7.3f}  {r.err_salmon:7.3f}  "
        f"{r.err_mm2:7.3f}  {r.abs_err_mm2:8.0f}"
    )

# Gene-level aggregation of the target set
print("\n" + "=" * 100)
print("TOP 10 GENES (by total abs_err_mm2 in target set)")
print("=" * 100)
gene_summary = target.groupby(["gene_id", "gene_name"]).agg(
    n_tx              = ("transcript_id", "count"),
    truth_total       = ("mrna_truth", "sum"),
    oracle_total      = ("rigel_oracle", "sum"),
    salmon_total      = ("salmon", "sum"),
    kallisto_total    = ("kallisto", "sum"),
    mm2_total         = ("rigel_minimap2", "sum"),
    abs_err_mm2_total = ("abs_err_mm2", "sum"),
).reset_index()
gene_summary["mm2_ratio"] = gene_summary["mm2_total"] / gene_summary["truth_total"].replace(0, np.nan)
gene_summary = gene_summary.sort_values("abs_err_mm2_total", ascending=False)

print(f"\n  {'gene_id':25s}  {'gene_name':18s}  ntx  {'truth':>8}  {'oracle':>8}  {'salmon':>8}  "
      f"{'kalli':>8}  {'mm2':>8}  {'mm2/truth':>9}  abs_err")
print("  " + "-" * 120)
for _, r in gene_summary.head(10).iterrows():
    ratio_str = f"{r.mm2_ratio:.3f}" if not np.isnan(r.mm2_ratio) else "    nan"
    print(
        f"  {r.gene_id:25s}  {str(r.gene_name):18s}  {r.n_tx:3d}  {r.truth_total:8.0f}  "
        f"{r.oracle_total:8.0f}  {r.salmon_total:8.0f}  {r.kallisto_total:8.0f}  "
        f"{r.mm2_total:8.0f}  {ratio_str:>9}  {r.abs_err_mm2_total:8.0f}"
    )

# Save top targets for deep investigation
top10_genes = gene_summary.head(10)["gene_id"].tolist()
print("\n" + "=" * 100)
print("DEEP DIVE: Per-transcript breakdown for top 10 genes")
print("=" * 100)
for gid in top10_genes:
    gene_rows = df[df["gene_id"].astype(str) == str(gid)].copy()
    gene_name = gene_rows["gene_name"].iloc[0]
    gene_rows["err_mm2_r"] = np.abs(gene_rows["rigel_minimap2"] - gene_rows["mrna_truth"]) / np.maximum(gene_rows["mrna_truth"], 1)
    gene_rows = gene_rows.sort_values("mrna_truth", ascending=False)
    print(f"\n  {gid}  ({gene_name})  [n_transcripts={len(gene_rows)}, total_truth={gene_rows.mrna_truth.sum():.0f}]")
    print(f"  {'transcript_id':42s}  {'truth':>7}  {'oracle':>7}  {'salmon':>7}  {'kalli':>7}  {'mm2':>7}  {'err_ora%':>8}  {'err_sal%':>8}  {'err_mm2%':>8}")
    for _, r in gene_rows.iterrows():
        denom = max(r.mrna_truth, 1)
        def pct(val):
            return f"{100*abs(val - r.mrna_truth)/denom:6.1f}%"
        print(
            f"    {r.transcript_id:40s}  {r.mrna_truth:7.0f}  "
            f"{r.rigel_oracle:7.0f}  {r.salmon:7.0f}  {r.kallisto:7.0f}  "
            f"{r.rigel_minimap2:7.0f}  {pct(r.rigel_oracle):>8}  {pct(r.salmon):>8}  "
            f"{pct(r.rigel_minimap2):>8}"
        )
