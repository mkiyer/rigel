#!/usr/bin/env python3
"""Compare two rigel quant output directories for numerical equivalence."""
import sys
import json
import pandas as pd
import numpy as np

def compare_feather(path1, path2, label, tol=1e-9):
    df1 = pd.read_feather(path1)
    df2 = pd.read_feather(path2)
    print(f"\n=== {label} ===")
    print(f"Shape: {df1.shape} vs {df2.shape}")
    assert df1.shape == df2.shape, "Shape mismatch!"

    all_ok = True
    for col in df1.columns:
        dtype = df1[col].dtype
        is_numeric = pd.api.types.is_numeric_dtype(dtype)
        if not is_numeric:
            eq = (df1[col] == df2[col]).all()
            print(f"  {col}: exact match = {eq}")
            if not eq:
                all_ok = False
        else:
            v1 = df1[col].to_numpy(dtype=float, na_value=0.0)
            v2 = df2[col].to_numpy(dtype=float, na_value=0.0)
            diff = np.abs(v1 - v2)
            maxdiff = diff.max()
            denom = np.maximum(np.abs(v1), 1e-30)
            reldiff = (diff / denom).max()
            corr = np.corrcoef(v1, v2)[0, 1] if v1.std() > 0 else 1.0
            ok = maxdiff < tol
            tag = "OK" if ok else "FAIL"
            print(f"  {col}: max_abs_diff={maxdiff:.2e}, "
                  f"max_rel_diff={reldiff:.2e}, corr={corr:.12f} [{tag}]")
            if not ok:
                all_ok = False
    return all_ok


def compare_json(path1, path2):
    with open(path1) as f:
        j1 = json.load(f)
    with open(path2) as f:
        j2 = json.load(f)

    print("\n=== summary.json ===")
    # Compare alignment stats (should be exact)
    for key in ["alignment", "library"]:
        if key in j1 and key in j2:
            for k, v1 in j1[key].items():
                v2 = j2[key].get(k)
                if isinstance(v1, float):
                    diff = abs(v1 - v2) if v2 is not None else float("inf")
                    ok = diff < 1e-6
                    print(f"  {key}.{k}: {v1} vs {v2} (diff={diff:.2e}) "
                          f"[{'OK' if ok else 'DIFF'}]")
                else:
                    eq = v1 == v2
                    if not eq:
                        print(f"  {key}.{k}: {v1} vs {v2} [DIFF]")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <baseline_dir> <test_dir>")
        sys.exit(1)

    base_dir = sys.argv[1]
    test_dir = sys.argv[2]

    ok1 = compare_feather(f"{base_dir}/quant.feather",
                          f"{test_dir}/quant.feather", "quant.feather")
    ok2 = compare_feather(f"{base_dir}/gene_quant.feather",
                          f"{test_dir}/gene_quant.feather", "gene_quant.feather")
    ok3 = compare_feather(f"{base_dir}/loci.feather",
                          f"{test_dir}/loci.feather", "loci.feather")

    compare_json(f"{base_dir}/summary.json", f"{test_dir}/summary.json")

    if ok1 and ok2 and ok3:
        print("\n*** ALL CHECKS PASSED ***")
    else:
        print("\n*** SOME CHECKS FAILED (see above) ***")
        sys.exit(1)


if __name__ == "__main__":
    main()
