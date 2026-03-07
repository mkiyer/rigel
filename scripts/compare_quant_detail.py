#!/usr/bin/env python3
"""Compare quant_detail.feather between two runs."""
import sys
import pyarrow.feather as feather
import pandas as pd

run1 = sys.argv[1]
run2 = sys.argv[2]

d1 = feather.read_table(f"{run1}/quant_detail.feather").to_pandas()
d2 = feather.read_table(f"{run2}/quant_detail.feather").to_pandas()

print(f"Run1 rows: {len(d1)}, Run2 rows: {len(d2)}")
print(f"Categories in run1: {sorted(d1['category'].unique())}")
print(f"Categories in run2: {sorted(d2['category'].unique())}")
print(f"Sources in run1: {sorted(d1['source'].unique())}")
print(f"Sources in run2: {sorted(d2['source'].unique())}")

# Convert to plain strings for concat
for col in ["transcript_id", "category", "source"]:
    d1[col] = d1[col].astype(str)
    d2[col] = d2[col].astype(str)

# Create keys
d1["key"] = d1["transcript_id"] + "|" + d1["category"] + "|" + d1["source"]
d2["key"] = d2["transcript_id"] + "|" + d2["category"] + "|" + d2["source"]

# Find rows only in one run
only_in_1 = set(d1["key"]) - set(d2["key"])
only_in_2 = set(d2["key"]) - set(d1["key"])
print(f"\nKeys only in run1: {len(only_in_1)}")
for k in sorted(only_in_1)[:20]:
    row = d1[d1["key"] == k].iloc[0]
    print(f"  {k}: count={row['count']}")

print(f"\nKeys only in run2: {len(only_in_2)}")
for k in sorted(only_in_2)[:20]:
    row = d2[d2["key"] == k].iloc[0]
    print(f"  {k}: count={row['count']}")

# Merge on key and compare counts
merged = d1.merge(d2, on="key", suffixes=("_1", "_2"))
diff = merged[merged["count_1"] != merged["count_2"]]
print(f"\nShared keys with different counts: {len(diff)}")
if len(diff) > 0:
    diff["abs_diff"] = (diff["count_1"] - diff["count_2"]).abs()
    diff["rel_diff"] = diff["abs_diff"] / diff[["count_1", "count_2"]].abs().max(axis=1).clip(lower=1e-30)
    print(f"  Max abs diff: {diff['abs_diff'].max():.6f}")
    print(f"  Max rel diff: {diff['rel_diff'].max():.6f}")
    print(f"  Mean abs diff: {diff['abs_diff'].mean():.6f}")
    
    # Show biggest diffs
    top = diff.nlargest(20, "abs_diff")
    for _, row in top.iterrows():
        print(f"  {row['key']}: {row['count_1']:.6f} vs {row['count_2']:.6f} (diff={row['abs_diff']:.6f})")

# Count by category
print("\n--- Counts by category ---")
for cat in sorted(d1["category"].unique()):
    c1 = d1[d1["category"] == cat]["count"].sum()
    c2 = d2[d2["category"] == cat]["count"].sum()
    print(f"  {cat}: {c1:.4f} vs {c2:.4f} (diff={c1-c2:.4f})")
