#!/usr/bin/env python3
"""Compare two benchmark summary.json files side-by-side."""
import json
import sys

v4_path = sys.argv[1]
v5_path = sys.argv[2]

with open(v4_path) as f:
    v4 = json.load(f)
with open(v5_path) as f:
    v5 = json.load(f)

print("=" * 90)
print("COMPARISON: v4 (baseline) vs v5 (Phase 3b)")
print("=" * 90)

for v4_cond, v5_cond in zip(v4, v5):
    aligner = v4_cond["aligner"]
    print(f"\n--- {v4_cond['dataset_name']} / {aligner} ---")

    for tool in v4_cond["transcript_metrics"]:
        m4 = v4_cond["transcript_metrics"][tool]
        m5 = v5_cond["transcript_metrics"].get(tool)
        if m5 is None:
            print(f"  {tool}: MISSING in v5!")
            continue

        print(f"  {tool} (transcript-level):")
        keys = [
            "mean_abs_error",
            "rmse",
            "pearson",
            "spearman",
            "total_observed",
            "total_abs_error",
            "elapsed_sec",
            "peak_rss_mb",
            "throughput_frags_per_sec",
        ]
        for key in keys:
            val4 = m4.get(key, "N/A")
            val5 = m5.get(key, "N/A")
            if isinstance(val4, (int, float)) and isinstance(val5, (int, float)):
                delta = val5 - val4
                pct = (delta / val4 * 100) if val4 != 0 else 0
                print(
                    f"    {key:30s}  v4={val4:>14.3f}  v5={val5:>14.3f}  "
                    f"delta={delta:>+12.3f} ({pct:>+6.1f}%)"
                )

    # Gene-level
    for tool in v4_cond.get("gene_metrics", {}):
        g4 = v4_cond["gene_metrics"][tool]
        g5 = v5_cond["gene_metrics"].get(tool)
        if g5 is None:
            continue
        print(f"  {tool} (gene-level):")
        for key in ["mean_abs_error", "rmse", "pearson", "spearman"]:
            val4 = g4.get(key, "N/A")
            val5 = g5.get(key, "N/A")
            if isinstance(val4, (int, float)) and isinstance(val5, (int, float)):
                delta = val5 - val4
                print(
                    f"    {key:30s}  v4={val4:>14.4f}  v5={val5:>14.4f}  "
                    f"delta={delta:>+12.4f}"
                )

    # Pool counts
    for tool in v4_cond.get("pool_counts", {}):
        p4 = v4_cond["pool_counts"][tool]
        p5 = v5_cond.get("pool_counts", {}).get(tool, {})
        print(f"  {tool} (pool counts):")
        for k in ["mature_rna", "nascent_rna", "genomic_dna"]:
            v4v = p4.get(k, 0)
            v5v = p5.get(k, 0)
            print(f"    {k:20s}  v4={v4v:>14.0f}  v5={v5v:>14.0f}  delta={v5v - v4v:>+10.0f}")
