"""How many candidate transcripts/nRNAs compete with gDNA per FP fragment?

For each hotspot window we already have per-fragment FP tables; tabulate the
ZN distribution (n_candidates) and the cumulative RNA posterior implied by:
   sum_k theta_k/eff_len_k vs theta_g/(2*gdna_eff_len)
using the dominant nRNA + median target eff_len.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DEEP = Path("/Users/mkiyer/proj/rigel/results/vcap_hotspot_deepdive_2026-05-17")

for tsv in sorted(DEEP.glob("chr*_fp.tsv")):
    df = pd.read_csv(tsv, sep="\t")
    if df.empty:
        continue
    name = tsv.stem.replace("_fp", "")
    print(f"\n=== {name} ===  ({len(df):,} FP fragments)")
    print("ZN (n_candidates) distribution:")
    print(df["zn"].describe(percentiles=[0.5, 0.9, 0.99]).to_string())
    print("ZN value_counts (top 10):")
    print(df["zn"].value_counts().head(10).to_string())
    # ZW distribution stratified by ZN
    print("\nmedian ZW by ZN-bin:")
    df["zn_bin"] = pd.cut(df["zn"], bins=[0, 2, 5, 10, 50, 200, 1000, 10000])
    print(df.groupby("zn_bin", observed=True)["zw"].agg(["count", "median", "mean"]).to_string())
