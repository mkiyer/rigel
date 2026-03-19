#!/usr/bin/env python3
"""Analyze benchmark golden results."""
import pandas as pd
import os

golden = "scripts/benchmark/golden"
for name in sorted(os.listdir(golden)):
    tsv = os.path.join(golden, name, "sweep_results.tsv")
    if not os.path.exists(tsv):
        continue
    df = pd.read_csv(tsv, sep="\t")
    print(f"\n=== {name} ({len(df)} runs) ===")

    # Overall mRNA accuracy
    if "total_mrna_rel_err" in df.columns:
        err = df["total_mrna_rel_err"]
        print(f"  mRNA rel_err: mean={err.mean():.3f}, median={err.median():.3f}, max={err.max():.3f}")

    # nRNA accuracy
    if "nrna_rel_err" in df.columns:
        nrna_runs = df[df["nrna_expected"] > 0]
        if len(nrna_runs) > 0:
            nerr = nrna_runs["nrna_rel_err"]
            print(f"  nRNA rel_err (n={len(nrna_runs)}): mean={nerr.mean():.3f}, median={nerr.median():.3f}, max={nerr.max():.3f}")

    # gDNA accuracy
    if "gdna_abs_diff" in df.columns:
        gdna_runs = df[df["gdna_expected"] > 0]
        if len(gdna_runs) > 0:
            gabs = gdna_runs["gdna_abs_diff"].abs()
            gexp = gdna_runs["gdna_expected"]
            grel = gabs / gexp
            print(f"  gDNA rel_err (n={len(gdna_runs)}): mean={grel.mean():.3f}, max={grel.max():.3f}")

    # Worst cases
    if "TA1_rel_err" in df.columns:
        worst = df.nlargest(5, "total_mrna_rel_err")
        print("  Worst 5 by mRNA rel_err:")
        key_cols = [c for c in df.columns if c in ("strand_specificity", "gdna_fraction",
            "TA1", "NTA1", "prior_pseudocount", "mode", "prune_threshold",
            "strand_specificity", "overhang_log_penalty", "mismatch_log_penalty")]
        for _, row in worst.iterrows():
            cfg = ", ".join(f"{c}={row[c]}" for c in key_cols if pd.notna(row[c]))
            print(f"    {cfg} => mRNA_err={row['total_mrna_rel_err']:.3f}, "
                  f"nRNA_diff={row.get('nrna_abs_diff', 0):.0f}, "
                  f"gDNA_diff={row.get('gdna_abs_diff', 0):.0f}")
