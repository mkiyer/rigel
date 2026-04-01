#!/usr/bin/env python3
"""Analyze loci.feather from MAP vs VBEM benchmark runs."""
import pandas as pd

BDIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/runs"

for cond in ["gdna_none_ss_1.00_nrna_none", "gdna_high_ss_0.90_nrna_none"]:
    for mode in ["map", "vbem"]:
        path = f"{BDIR}/{cond}/rigel/{mode}/loci.feather"
        loci = pd.read_feather(path)
        print(f"=== {cond} / {mode} ===")
        print(f"  Columns: {list(loci.columns)}")
        print(f"  Shape: {loci.shape}")
        # Find iteration/convergence related columns
        iter_cols = [
            c
            for c in loci.columns
            if "iter" in c.lower()
            or "converg" in c.lower()
            or "time" in c.lower()
            or "squarem" in c.lower()
        ]
        if iter_cols:
            for col in iter_cols:
                s = loci[col]
                print(
                    f"  {col}: mean={s.mean():.3f}, max={s.max()}, sum={s.sum():.1f}"
                )
        # Top 5 by n_transcripts
        size_col = (
            "n_transcripts" if "n_transcripts" in loci.columns else loci.columns[0]
        )
        top = loci.nlargest(5, size_col)
        show_cols = [
            c
            for c in [
                "locus_id",
                "n_transcripts",
                "n_fragments",
                "n_ec",
                "n_components",
            ]
            + iter_cols
            if c in top.columns
        ]
        print(f"  Top 5 loci by {size_col}:")
        print(top[show_cols].to_string(index=False))
        print()
