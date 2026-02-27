#!/usr/bin/env python3
"""Analyze wrong-pool fragments from per-read accuracy CSV."""
import csv
from collections import Counter

csv_path = (
    "/Users/mkiyer/Downloads/hulkrna_runs/bench_chr8_pvt1_myc/PVT1_MYC/"
    "gdna_none_nrna_none_ss_1.00/align_oracle/accuracy_default/per_read_accuracy.csv"
)

with open(csv_path) as f:
    reader = csv.DictReader(f)
    print("Columns:", reader.fieldnames)

    wrong_pool = Counter()
    assigned_pool = Counter()
    assigned_tid = Counter()
    assigned_gene = Counter()
    for r in reader:
        if r["verdict"] == "wrong_pool":
            wrong_pool[r["truth_tid"]] += 1
            assigned_pool[r["assigned_pool"]] += 1
            assigned_tid[r.get("assigned_tid", "N/A")] += 1
            assigned_gene[r.get("assigned_gid", "N/A")] += 1

    print(f"\nTotal wrong_pool fragments: {sum(wrong_pool.values())}")

    print(f"\nWrong pool by truth transcript (top 10):")
    for tid, count in wrong_pool.most_common(10):
        print(f"  {tid}: {count}")

    print(f"\nWrong pool assigned to which pool:")
    for pool, count in assigned_pool.most_common():
        print(f"  {pool}: {count}")

    print(f"\nWrong pool assigned to which transcript (top 10):")
    for tid, count in assigned_tid.most_common(10):
        print(f"  {tid}: {count}")

    print(f"\nWrong pool assigned to which gene (top 10):")
    for gid, count in assigned_gene.most_common(10):
        print(f"  {gid}: {count}")
