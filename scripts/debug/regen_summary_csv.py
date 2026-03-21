#!/usr/bin/env python3
"""Regenerate summary.csv from updated summary.json."""
import json
from pathlib import Path

import pandas as pd

OUTDIR = Path("~/Downloads/rigel_runs/benchmark_output_v6_prune").expanduser()

with open(OUTDIR / "summary.json") as f:
    data = json.load(f)

rows = []
for r in data:
    base = {
        "dataset": r["dataset_name"],
        "gdna_label": r.get("gdna_label", ""),
        "gdna_rate": r.get("gdna_rate", 0),
        "strand_specificity": r.get("strand_specificity", 0),
        "aligner": r["aligner"],
        "n_rna_fragments": r.get("n_rna_fragments", 0),
        "n_gdna": r.get("n_gdna_truth", 0),
        "alignment_elapsed_sec": r.get("alignment_elapsed", 0),
    }
    for tool, m in r.get("transcript_metrics", {}).items():
        row = {**base, "level": "transcript", "tool": tool}
        if isinstance(m, dict):
            row.update({k: v for k, v in m.items() if k != "tool"})
        rows.append(row)
    for tool, m in r.get("gene_metrics", {}).items():
        row = {**base, "level": "gene", "tool": tool}
        if isinstance(m, dict):
            row.update({k: v for k, v in m.items() if k != "tool"})
        rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUTDIR / "summary.csv", index=False)
print(f"Wrote summary.csv ({len(df)} rows)")
