#!/usr/bin/env python
"""Compare two real-data rigel runs to quantify non-determinism."""
import pyarrow.feather as feather
import numpy as np

import sys
base = sys.argv[1] if len(sys.argv) > 1 else "/Users/mkiyer/Downloads/rigel_runs/baseline_output"
orig = sys.argv[2] if len(sys.argv) > 2 else "/Users/mkiyer/Downloads/rigel_runs/orig_run1"

for fname in ["quant.feather", "gene_quant.feather", "loci.feather"]:
    print(f"=== {fname} ===")
    df1 = feather.read_table(f"{base}/{fname}").to_pandas()
    df2 = feather.read_table(f"{orig}/{fname}").to_pandas()
    print(f"  rows: {len(df1)} vs {len(df2)}")

    for col in df1.columns:
        if df1[col].dtype in ["float64", "float32", "int64", "int32"]:
            a = df1[col].values.astype(float)
            b = df2[col].values.astype(float)
            if np.array_equal(a, b, equal_nan=True):
                continue
            mask = ~(np.isnan(a) & np.isnan(b))
            diff = np.abs(a[mask] - b[mask])
            maxd = np.max(diff)
            denom = np.maximum(np.abs(a[mask]), np.abs(b[mask]))
            denom = np.where(denom == 0, 1, denom)
            rel = np.max(diff / denom)
            ndiff = np.sum(diff > 0)
            print(f"  {col}: max_abs={maxd:.6e}, max_rel={rel:.6e}, n_diff={ndiff}/{len(a[mask])}")
        elif df1[col].dtype == "object":
            if not (df1[col] == df2[col]).all():
                n = (df1[col] != df2[col]).sum()
                print(f"  {col}: {n} string values differ")
