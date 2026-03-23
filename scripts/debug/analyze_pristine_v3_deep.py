#!/usr/bin/env python3
"""Deep root-cause analysis of minimap2-specific errors in pristine v3 benchmark."""
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

V3 = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/gdna_none_ss_0.95_nrna_none")
IDX = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v3/rigel_index")

# Load data
orc_df = pd.read_csv(V3 / "per_transcript_counts_oracle.csv")
mm3_df = pd.read_csv(V3 / "per_transcript_counts_minimap2.csv")

# Build unified df
df = pd.DataFrame({
    "tid": orc_df["transcript_id"],
    "gid": orc_df["gene_id"],
    "truth": orc_df["mrna_truth"],
    "oracle": orc_df["rigel_oracle"],
    "mm2": mm3_df["rigel_minimap2"],
    "salmon": orc_df["salmon"],
    "kallisto": orc_df["kallisto"],
})
df["orc_err"] = df["oracle"] - df["truth"]
df["mm2_err"] = df["mm2"] - df["truth"]
df["mm2_abs_err"] = df["mm2_err"].abs()
df["orc_abs_err"] = df["orc_err"].abs()
df["mm2_excess"] = df["mm2_abs_err"] - df["orc_abs_err"]

# Load transcript metadata
tx_df = pd.read_feather(IDX / "transcripts.feather")
# Columns: ref, start, end, strand, length, t_id, g_id, t_index, g_index, g_name,
#          g_type, is_basic, is_mane, is_ccds, is_synthetic_nrna, ...
print(f"Index transcripts: {len(tx_df)}, columns: {tx_df.columns.tolist()}")

# Need to map t_id to per_transcript CSV transcript_id
# The t_id in index doesn't have version; let's just use n_exons from intervals
# Actually, use what we can from index
tx_meta = {}
for _, row in tx_df.iterrows():
    tx_meta[row["t_id"]] = {
        "ref": row["ref"], "strand": row["strand"],
        "length": row["length"],
        "is_nrna": bool(row["is_synthetic_nrna"]),
    }

# Count exons from intervals
intv_df = pd.read_feather(IDX / "intervals.feather")
print(f"Intervals: {len(intv_df)}, columns: {intv_df.columns.tolist()}")
# Count exon intervals per transcript
if "itype" in intv_df.columns:
    exon_counts = intv_df[intv_df["itype"] == 0].groupby("t_index").size().to_dict()
else:
    exon_counts = {}

for _, row in tx_df.iterrows():
    tid = row["t_id"]
    if tid in tx_meta:
        tx_meta[tid]["n_exons"] = exon_counts.get(row["t_index"], 1)
        tx_meta[tid]["genomic_span"] = row["end"] - row["start"]

def get_meta(tid, key, default=None):
    return tx_meta.get(tid, {}).get(key, default)

df["ref"] = df["tid"].map(lambda t: get_meta(t, "ref", "?"))
df["strand"] = df["tid"].map(lambda t: get_meta(t, "strand", "?"))
df["n_exons"] = df["tid"].map(lambda t: get_meta(t, "n_exons", 0))
df["spliced_len"] = df["tid"].map(lambda t: get_meta(t, "length", 0))
df["genomic_span"] = df["tid"].map(lambda t: get_meta(t, "genomic_span", 0))
df["is_nrna"] = df["tid"].map(lambda t: get_meta(t, "is_nrna", False))

# Filter to non-nRNA annotated transcripts
real = df[~df["is_nrna"]].copy()
print(f"Annotated transcripts: {len(real)}, nRNA synthetics: {df['is_nrna'].sum()}")

# ── 1. MM2-specific errors: where mm2 is bad but oracle is fine ────
print("\n" + "="*80)
print("1. MM2-SPECIFIC HIGH-ERROR TRANSCRIPTS (mm2 |err|>50, oracle |err|<10)")
print("="*80)

mm2_specific = real[(real["mm2_abs_err"] > 50) & (real["orc_abs_err"] < 10)].copy()
print(f"\nCount: {len(mm2_specific)} transcripts")
print(f"Total mm2 excess error at these tx: {mm2_specific['mm2_abs_err'].sum():,.0f}")

over = mm2_specific[mm2_specific["mm2_err"] > 0]
under = mm2_specific[mm2_specific["mm2_err"] < 0]
print(f"\n  Over-predicted:  {len(over)} tx, total excess: +{over['mm2_err'].sum():,.0f}")
print(f"  Under-predicted: {len(under)} tx, total deficit: {under['mm2_err'].sum():,.0f}")

# ── 2. Are over-predictions concentrated at certain genes? ────
print("\n" + "="*80)
print("2. GENE-LEVEL PATTERN OF MM2-SPECIFIC OVER-PREDICTIONS")
print("="*80)

over_by_gene = over.groupby("gid").agg(
    n_tx=("tid", "count"),
    total_excess=("mm2_err", "sum"),
    total_truth=("truth", "sum"),
    total_mm2=("mm2", "sum"),
).sort_values("total_excess", ascending=False)
print(f"\nGenes with over-predicted transcripts: {len(over_by_gene)}")
print(f"\nTop 20 genes by mm2 excess:")
print(f"{'Gene':<30} {'N_tx':>5} {'Truth':>8} {'MM2':>8} {'Excess':>8}")
print("-"*62)
for gid, row in over_by_gene.head(20).iterrows():
    print(f"{gid:<30} {row['n_tx']:>5} {row['total_truth']:>8,.0f} {row['total_mm2']:>8,.0f} {row['total_excess']:>+8,.0f}")

# ── 3. n_exons and spliced_length vs error ────
print("\n" + "="*80)
print("3. MM2 ERROR vs TRANSCRIPT STRUCTURE")
print("="*80)

# For mm2-specific errors
print(f"\nMM2-specific over-predicted (n={len(over)}):")
print(f"  Median n_exons: {over['n_exons'].median():.0f} (all tx median: {real['n_exons'].median():.0f})")
print(f"  Median spliced_len: {over['spliced_len'].median():.0f} (all: {real['spliced_len'].median():.0f})")
print(f"  Median genomic_span: {over['genomic_span'].median():.0f} (all: {real['genomic_span'].median():.0f})")
print(f"  Single-exon fraction: {(over['n_exons']==1).mean():.2%} (all: {(real['n_exons']==1).mean():.2%})")

print(f"\nMM2-specific under-predicted (n={len(under)}):")
print(f"  Median n_exons: {under['n_exons'].median():.0f}")
print(f"  Median spliced_len: {under['spliced_len'].median():.0f}")
print(f"  Median genomic_span: {under['genomic_span'].median():.0f}")
print(f"  Single-exon fraction: {(under['n_exons']==1).mean():.2%}")

# Correlation of n_exons with excess error
from scipy.stats import spearmanr
expressed = real[real["truth"] > 0].copy()
r, p = spearmanr(expressed["n_exons"], expressed["mm2_excess"])
print(f"\n  Spearman(n_exons, mm2_excess): r={r:.4f}, p={p:.2e}")
r2, p2 = spearmanr(expressed["genomic_span"], expressed["mm2_excess"])
print(f"  Spearman(genomic_span, mm2_excess): r={r2:.4f}, p={p2:.2e}")

# ── 4. Isoform redistribution: genes where oracle is fine but mm2 shuffles ────
print("\n" + "="*80)
print("4. ISOFORM REDISTRIBUTION (gene OK, isoforms wrong)")
print("="*80)

gene_df = real.groupby("gid").agg(
    n_tx=("tid", "count"),
    truth=("truth", "sum"),
    oracle=("oracle", "sum"),
    mm2=("mm2", "sum"),
    salmon=("salmon", "sum"),
    orc_abs_err=("orc_abs_err", "sum"),
    mm2_abs_err=("mm2_abs_err", "sum"),
).reset_index()
gene_df["gene_orc_err"] = (gene_df["oracle"] - gene_df["truth"]).abs()
gene_df["gene_mm2_err"] = (gene_df["mm2"] - gene_df["truth"]).abs()

# Genes where gene-level is fine but transcript-level is bad
redistributed = gene_df[
    (gene_df["gene_mm2_err"] < 20) &  # gene level OK
    (gene_df["mm2_abs_err"] > 100) &  # transcript level bad
    (gene_df["n_tx"] > 1)  # multi-isoform
].sort_values("mm2_abs_err", ascending=False)

print(f"\nGenes with gene-level |err|<20 but tx-level sum|err|>100 (multi-isoform): {len(redistributed)}")
print(f"Total isoform redistribution error: {redistributed['mm2_abs_err'].sum():,.0f}")

print(f"\nTop 20 redistributed genes:")
print(f"{'Gene':<30} {'N_tx':>5} {'Truth':>8} {'Gene Err':>9} {'Tx |Err|':>9}")
print("-"*65)
for _, row in redistributed.head(20).iterrows():
    print(f"{row['gid']:<30} {row['n_tx']:>5} {row['truth']:>8,.0f} {row['gene_mm2_err']:>9.0f} {row['mm2_abs_err']:>9,.0f}")

total_iso_err = redistributed["mm2_abs_err"].sum()
total_mm2_err = real["mm2_abs_err"].sum()
print(f"\nIsoform redistribution: {total_iso_err:,.0f} / {total_mm2_err:,.0f} = {total_iso_err/total_mm2_err:.1%} of total mm2 error")

# ── 5. Chromosomal distribution ────
print("\n" + "="*80)
print("5. ERROR BY CHROMOSOME")
print("="*80)

chr_errs = real.groupby("ref").agg(
    n_tx=("tid", "count"),
    truth=("truth", "sum"),
    mm2_abs_err=("mm2_abs_err", "sum"),
    orc_abs_err=("orc_abs_err", "sum"),
    mm2_excess=("mm2_excess", "sum"),
).sort_values("mm2_excess", ascending=False)
print(f"\n{'Chrom':<15} {'N_tx':>6} {'Truth':>10} {'Orc |Err|':>10} {'MM2 |Err|':>10} {'MM2 Excess':>10}")
print("-"*65)
for ref, row in chr_errs.head(25).iterrows():
    print(f"{ref:<15} {row['n_tx']:>6} {row['truth']:>10,.0f} {row['orc_abs_err']:>10,.0f} {row['mm2_abs_err']:>10,.0f} {row['mm2_excess']:>+10,.0f}")

# ── 6. nRNA siphon details ────
print("\n" + "="*80)
print("6. nRNA SIPHON DETAIL")
print("="*80)

nrna_tx = df[df["is_nrna"]].copy()
nrna_with_counts = nrna_tx[nrna_tx["mm2"] > 0]
print(f"\nnRNA transcripts with mm2 count > 0: {len(nrna_with_counts)} / {len(nrna_tx)}")
print(f"Total mm2 counts in nRNA: {nrna_tx['mm2'].sum():,.0f}")
print(f"Total oracle counts in nRNA: {nrna_tx['oracle'].sum():,.0f}")

# ── 7. Fragments missing from mm2 ────
print("\n" + "="*80)
print("7. FRAGMENT ACCOUNTING (where do 336K missing fragments go?)")
print("="*80)

total_truth = real["truth"].sum()
total_mm2 = real["mm2"].sum()
total_orc = real["oracle"].sum()
total_mm2_nrna = nrna_tx["mm2"].sum()
total_orc_nrna = nrna_tx["oracle"].sum()

print(f"\n  Truth (mRNA):           {total_truth:>12,.0f}")
print(f"  Oracle mRNA:            {total_orc:>12,.0f}")
print(f"  Oracle nRNA:            {total_orc_nrna:>12,.0f}")
print(f"  Oracle gDNA:            {10000000 - total_orc - total_orc_nrna:>12,.0f}")
print(f"  MM2 mRNA:               {total_mm2:>12,.0f}")
print(f"  MM2 nRNA:               {total_mm2_nrna:>12,.0f}")
print(f"  MM2 gDNA:               {10000000 - total_mm2 - total_mm2_nrna:>12,.0f}")
print(f"  Salmon total:           {real['salmon'].sum():>12,.0f}")
print(f"  Kallisto total:         {real['kallisto'].sum():>12,.0f}")
print(f"\n  MM2 mRNA deficit:       {total_mm2 - total_truth:>+12,.0f}")
print(f"  MM2 nRNA siphon:        {total_mm2_nrna:>+12,.0f}")

# How many fragments does mm2 lose to unmapped/intergenic?
mm2_total = total_mm2 + total_mm2_nrna  # + gDNA from summary
print(f"\n  Fragments quantified by mm2 (mRNA+nRNA): {mm2_total:>12,.0f}")
print(f"  Fragments lost (unmapped/intergenic):     {total_truth - mm2_total:>12,.0f}")

# ── 8. Top genes where mm2 errors are concentrated ────
print("\n" + "="*80)
print("8. TOP 30 GENES BY MM2 |ERROR| (expressed genes only)")
print("="*80)

expressed_genes = gene_df[gene_df["truth"] > 0].sort_values("mm2_abs_err", ascending=False)
print(f"\n{'Gene':<30} {'N_tx':>5} {'Truth':>8} {'Oracle':>8} {'MM2':>8} {'Gene Err':>9} {'Tx |Err|':>9}")
print("-"*80)
for _, row in expressed_genes.head(30).iterrows():
    print(f"{row['gid']:<30} {row['n_tx']:>5} {row['truth']:>8,.0f} {row['oracle']:>8,.0f} {row['mm2']:>8,.0f} {row['gene_mm2_err']:>9.0f} {row['mm2_abs_err']:>9,.0f}")

# ── 9. Summary of where the remaining gap comes from ────
print("\n" + "="*80)
print("9. REMAINING GAP DECOMPOSITION SUMMARY")
print("="*80)

mm2_total_err = real["mm2_abs_err"].sum()
orc_total_err = real["orc_abs_err"].sum()

# Shared error (present in both oracle and mm2)
shared = real.copy()
shared["shared_err"] = shared[["mm2_abs_err", "orc_abs_err"]].min(axis=1)
shared_total = shared["shared_err"].sum()

# Isoform redistribution error (gene OK, tx bad)
iso_err = redistributed["mm2_abs_err"].sum()

# Pure mm2-specific (oracle < 10, mm2 > 50)
mm2_only_err = mm2_specific["mm2_abs_err"].sum()

print(f"\n  Total mm2 |error|:                {mm2_total_err:>12,.0f}")
print(f"  Total oracle |error|:             {orc_total_err:>12,.0f}")
print(f"  Shared error (min of both):       {shared_total:>12,.0f} ({shared_total/mm2_total_err:.1%})")
print(f"  MM2-specific error:               {mm2_total_err - shared_total:>12,.0f} ({(mm2_total_err-shared_total)/mm2_total_err:.1%})")
print(f"  Of which isoform redistribution:  {iso_err:>12,.0f}")
print(f"  Of which pure mm2 (orc<10):       {mm2_only_err:>12,.0f}")

print("\n=== ROOT CAUSE ANALYSIS COMPLETE ===")
