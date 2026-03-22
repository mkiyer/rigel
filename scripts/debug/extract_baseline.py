#!/usr/bin/env python3
"""Extract baseline benchmark metrics for the 3 stress-test conditions."""
import json

SUMMARY = "/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune/summary.json"
TARGETS = [
    "gdna_high_ss_1.00_nrna_default",
    "gdna_high_ss_0.90_nrna_default",
    "gdna_high_ss_0.50_nrna_default",
]

with open(SUMMARY) as f:
    data = json.load(f)

for entry in data:
    if entry["dataset_name"] not in TARGETS:
        continue
    ds = entry["dataset_name"]
    al = entry["aligner"]
    print(f"=== {ds} | aligner={al} ===")
    print(f"  gdna_rate={entry['gdna_rate']}, ss={entry['strand_specificity']}")
    print(
        f"  n_rna={entry['n_rna_fragments']}, n_gdna={entry.get('n_gdna','?')}, "
        f"n_total={entry['n_total_fragments']}"
    )
    print(
        f"  n_mrna_truth={entry['n_mrna_truth']}, "
        f"  n_nrna_truth={entry['n_nrna_truth']}, "
        f"  n_gdna_truth={entry['n_gdna_truth']}"
    )
    for tool_name, m in entry["transcript_metrics"].items():
        if "rigel" in tool_name:
            print(
                f"  TRANSCRIPT {tool_name}: "
                f"MAE={m['mean_abs_error']:.4f}, "
                f"RMSE={m['rmse']:.4f}, "
                f"pearson={m['pearson']:.6f}, "
                f"spearman={m['spearman']:.4f}"
            )
            print(
                f"    total_truth={m['total_truth']:.1f}, "
                f"total_observed={m['total_observed']:.1f}, "
                f"total_abs_error={m['total_abs_error']:.1f}"
            )
    for tool_name, m in entry["gene_metrics"].items():
        if "rigel" in tool_name:
            print(
                f"  GENE      {tool_name}: "
                f"MAE={m['mean_abs_error']:.4f}, "
                f"RMSE={m['rmse']:.4f}, "
                f"pearson={m['pearson']:.6f}, "
                f"spearman={m['spearman']:.4f}"
            )
    if "pool_counts" in entry:
        for tool_name, pc in entry["pool_counts"].items():
            if "rigel" in tool_name:
                print(f"  POOLS     {tool_name}: {pc}")
    print()
