#!/usr/bin/env python3
"""Generate comprehensive statistics for the merged human genome region table."""

import numpy as np
import pandas as pd

df = pd.read_parquet("/tmp/human_regions_merged.parquet")

# Classify regions
df["is_intergenic"] = ~df["tx_pos"] & ~df["tx_neg"]
df["is_exonic_only_pos"] = df["exon_pos"] & ~df["exon_neg"] & ~df["tx_neg"]
df["is_exonic_only_neg"] = df["exon_neg"] & ~df["exon_pos"] & ~df["tx_pos"]
df["is_exonic_ambig"] = df["exon_pos"] & df["exon_neg"]
df["is_exonic_any"] = df["exon_pos"] | df["exon_neg"]
df["is_intronic_only_pos"] = df["tx_pos"] & ~df["exon_pos"] & ~df["tx_neg"]
df["is_intronic_only_neg"] = df["tx_neg"] & ~df["exon_neg"] & ~df["tx_pos"]
df["is_intronic_mixed"] = (df["tx_pos"] & df["tx_neg"]) & ~df["is_exonic_any"]


def classify(row):
    if row["is_intergenic"]:
        return "intergenic"
    if row["is_exonic_ambig"]:
        return "exon_ambig"
    if row["is_exonic_only_pos"]:
        return "exon_pos"
    if row["is_exonic_only_neg"]:
        return "exon_neg"
    if row["exon_pos"] and row["tx_neg"]:
        return "exon_pos+intronic_neg"
    if row["exon_neg"] and row["tx_pos"]:
        return "exon_neg+intronic_pos"
    if row["is_intronic_mixed"]:
        return "intronic_ambig"
    if row["is_intronic_only_pos"]:
        return "intronic_pos"
    if row["is_intronic_only_neg"]:
        return "intronic_neg"
    return "other"


df["region_type"] = df.apply(classify, axis=1)

total_bp = df["length"].sum()
n_regions = len(df)

print("=" * 80)
print("HUMAN GENOME REGION PARTITION STATISTICS (MERGED)")
print("=" * 80)
print(f"Total regions:     {n_regions:>12,}")
print(f"Total genome bp:   {total_bp:>12,}")
print(f"References:        {df['ref'].nunique():>12,}")
print()

# Overall size distribution
print("-" * 80)
print("OVERALL REGION SIZE DISTRIBUTION")
print("-" * 80)
pcts = [0, 1, 5, 10, 25, 50, 75, 90, 95, 99, 100]
vals = np.percentile(df["length"], pcts)
for p, v in zip(pcts, vals):
    print(f"  {p:>3}th percentile:  {v:>12,.0f} bp")
print(f"  Mean:             {df['length'].mean():>12,.0f} bp")
print(f"  Std:              {df['length'].std():>12,.0f} bp")
print()

# Region type breakdown
print("-" * 80)
print("REGION TYPE BREAKDOWN")
print("-" * 80)
grp = (
    df.groupby("region_type")
    .agg(
        count=("length", "size"),
        total_bp=("length", "sum"),
        mean_bp=("length", "mean"),
        median_bp=("length", "median"),
        min_bp=("length", "min"),
        max_bp=("length", "max"),
    )
    .sort_values("total_bp", ascending=False)
)
grp["pct_regions"] = 100.0 * grp["count"] / n_regions
grp["pct_genome"] = 100.0 * grp["total_bp"] / total_bp

for rt, row in grp.iterrows():
    print(
        f"  {rt:<25s}  n={row['count']:>8,}  ({row['pct_regions']:5.1f}%)  "
        f"bp={row['total_bp']:>14,}  ({row['pct_genome']:5.1f}%)  "
        f"med={row['median_bp']:>8,.0f}  mean={row['mean_bp']:>10,.0f}  "
        f"range=[{row['min_bp']:,}-{row['max_bp']:,}]"
    )
print()

# Small region analysis
print("-" * 80)
print("SMALL REGION ANALYSIS (cumulative counts and genome fraction)")
print("-" * 80)
thresholds = [10, 25, 50, 100, 200, 500, 1000, 2000, 5000]
for t in thresholds:
    mask = df["length"] < t
    n = mask.sum()
    bp = df.loc[mask, "length"].sum()
    print(
        f"  < {t:>5,} bp:  {n:>8,} regions ({100*n/n_regions:5.1f}%),"
        f"  {bp:>12,} bp ({100*bp/total_bp:.4f}%)"
    )
print()

# Small region breakdown by type
print("-" * 80)
print("SMALL REGIONS (<200 bp) BY TYPE")
print("-" * 80)
small = df[df["length"] < 200]
if len(small) > 0:
    sg = (
        small.groupby("region_type")
        .agg(
            count=("length", "size"),
            total_bp=("length", "sum"),
            median_bp=("length", "median"),
        )
        .sort_values("count", ascending=False)
    )
    for rt, row in sg.iterrows():
        print(
            f"  {rt:<25s}  n={row['count']:>8,}"
            f"  bp={row['total_bp']:>10,}  med={row['median_bp']:>6,.0f}"
        )
print()

# Size distribution per type
print("-" * 80)
print("SIZE PERCENTILES BY REGION TYPE (bp)")
print("-" * 80)
for rt in grp.index:
    sub = df.loc[df["region_type"] == rt, "length"]
    ps = np.percentile(sub, [5, 25, 50, 75, 95])
    print(
        f"  {rt:<25s}  5%={ps[0]:>8,.0f}  25%={ps[1]:>8,.0f}"
        f"  50%={ps[2]:>8,.0f}  75%={ps[3]:>8,.0f}  95%={ps[4]:>8,.0f}"
    )
print()

# Flag combination distribution
print("-" * 80)
print("FLAG COMBINATIONS (unique flag tuples)")
print("-" * 80)
flag_cols = ["exon_pos", "exon_neg", "tx_pos", "tx_neg"]
fc = (
    df.groupby(flag_cols)
    .agg(count=("length", "size"), total_bp=("length", "sum"))
    .sort_values("count", ascending=False)
)
for flags, row in fc.iterrows():
    ep, en, tp, tn = flags
    label = f"ep={int(ep)} en={int(en)} tp={int(tp)} tn={int(tn)}"
    print(f"  {label:<30s}  n={row['count']:>8,}  bp={row['total_bp']:>14,}")
print()

# Intergenic region sizes
print("-" * 80)
print("INTERGENIC REGIONS (no transcript overlap)")
print("-" * 80)
ig = df[df["is_intergenic"]]
print(f"  Count:  {len(ig):,}")
print(
    f"  Total bp: {ig['length'].sum():,}"
    f" ({100*ig['length'].sum()/total_bp:.1f}% of genome)"
)
ps = np.percentile(ig["length"], [5, 25, 50, 75, 95])
print(
    f"  Size percentiles: 5%={ps[0]:,.0f}  25%={ps[1]:,.0f}"
    f"  50%={ps[2]:,.0f}  75%={ps[3]:,.0f}  95%={ps[4]:,.0f}"
)
print()

# Exonic region sizes
print("-" * 80)
print("EXONIC REGIONS (any strand)")
print("-" * 80)
ex = df[df["is_exonic_any"]]
print(f"  Count:  {len(ex):,}")
print(
    f"  Total bp: {ex['length'].sum():,}"
    f" ({100*ex['length'].sum()/total_bp:.1f}% of genome)"
)
ps = np.percentile(ex["length"], [5, 25, 50, 75, 95])
print(
    f"  Size percentiles: 5%={ps[0]:,.0f}  25%={ps[1]:,.0f}"
    f"  50%={ps[2]:,.0f}  75%={ps[3]:,.0f}  95%={ps[4]:,.0f}"
)
print()

# Intronic region sizes (tx but not exon)
print("-" * 80)
print("INTRONIC REGIONS (tx but not exon, any strand)")
print("-" * 80)
intr = df[(df["tx_pos"] | df["tx_neg"]) & ~df["is_exonic_any"]]
print(f"  Count:  {len(intr):,}")
print(
    f"  Total bp: {intr['length'].sum():,}"
    f" ({100*intr['length'].sum()/total_bp:.1f}% of genome)"
)
ps = np.percentile(intr["length"], [5, 25, 50, 75, 95])
print(
    f"  Size percentiles: 5%={ps[0]:,.0f}  25%={ps[1]:,.0f}"
    f"  50%={ps[2]:,.0f}  75%={ps[3]:,.0f}  95%={ps[4]:,.0f}"
)
print()

# Pre-merge vs merged comparison
print("-" * 80)
print("MERGE IMPACT")
print("-" * 80)
pre_merge_count = 1_043_881  # from previous run
print(f"  Pre-merge regions:   {pre_merge_count:>10,}")
print(f"  Post-merge regions:  {n_regions:>10,}")
print(f"  Reduction:           {pre_merge_count - n_regions:>10,} ({100*(pre_merge_count - n_regions)/pre_merge_count:.1f}%)")
