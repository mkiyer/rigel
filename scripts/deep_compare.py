#!/usr/bin/env python
"""Deep comparison of two rigel runs: check pre-EM and post-EM values separately."""
import sys
import pyarrow.feather as feather
import numpy as np

run1 = sys.argv[1]
run2 = sys.argv[2]

for fname in ["quant.feather", "quant_detail.feather", "loci.feather"]:
    print(f"\n{'='*60}")
    print(f"=== {fname} ===")
    print(f"{'='*60}")
    try:
        df1 = feather.read_table(f"{run1}/{fname}").to_pandas()
        df2 = feather.read_table(f"{run2}/{fname}").to_pandas()
    except Exception as e:
        print(f"  SKIP: {e}")
        continue

    print(f"  rows: {len(df1)} vs {len(df2)}")

    for col in df1.columns:
        if df1[col].dtype not in ["float64", "float32", "int64", "int32"]:
            continue
        a = df1[col].values.astype(float)
        b = df2[col].values.astype(float)
        both_nan = np.isnan(a) & np.isnan(b)
        if np.array_equal(a, b, equal_nan=True):
            print(f"  {col}: BIT-EXACT")
        else:
            mask = ~both_nan
            diff = np.abs(a[mask] - b[mask])
            maxd = np.max(diff)
            denom = np.maximum(np.abs(a[mask]), np.abs(b[mask]))
            denom = np.where(denom == 0, 1, denom)
            rel = np.max(diff / denom)
            ndiff = int(np.sum(diff > 0))
            print(f"  {col}: max_abs={maxd:.6e}, max_rel={rel:.6e}, n_diff={ndiff}/{len(a[mask])}")
            # Show worst cases
            worst_idx = np.argmax(diff)
            print(f"    worst at row {worst_idx}: {a[mask][worst_idx]:.18e} vs {b[mask][worst_idx]:.18e}")
