#!/usr/bin/env python3
"""Compare two rigel quant outputs for regression testing.

Usage:
    python3 scripts/regression_check.py <baseline_dir> <test_dir> [--tol 1e-6]

Compares quant.feather, gene_quant.feather, and quant_detail.feather
between baseline and test directories. Reports max absolute and relative
differences for all numeric columns.
"""

import argparse
import sys

import numpy as np
import pyarrow.feather as pf


def compare_feather(baseline_path, test_path, label, atol=1e-6, rtol=1e-6):
    """Compare two feather files column-by-column."""
    try:
        df_base = pf.read_table(baseline_path).to_pandas()
        df_test = pf.read_table(test_path).to_pandas()
    except FileNotFoundError as e:
        print(f"  SKIP {label}: {e}")
        return True

    if df_base.shape != df_test.shape:
        print(f"  FAIL {label}: shape mismatch {df_base.shape} vs {df_test.shape}")
        return False

    all_pass = True
    numeric_cols = df_base.select_dtypes(include=[np.number]).columns
    for col in numeric_cols:
        a = df_base[col].values.astype(np.float64)
        b = df_test[col].values.astype(np.float64)
        # Handle NaNs: both NaN = match, one NaN = mismatch
        nan_a = np.isnan(a)
        nan_b = np.isnan(b)
        nan_mismatch = np.sum(nan_a != nan_b)
        if nan_mismatch > 0:
            print(f"  FAIL {label}.{col}: {nan_mismatch} NaN mismatches")
            all_pass = False
            continue
        # Compare non-NaN values
        mask = ~nan_a
        if mask.sum() == 0:
            continue
        a_valid = a[mask]
        b_valid = b[mask]
        abs_diff = np.abs(a_valid - b_valid)
        max_abs = np.max(abs_diff)
        denom = np.maximum(np.abs(a_valid), 1e-15)
        max_rel = np.max(abs_diff / denom)
        close = np.allclose(a_valid, b_valid, atol=atol, rtol=rtol)
        status = "PASS" if close else "FAIL"
        if not close:
            all_pass = False
            n_diff = np.sum(~np.isclose(a_valid, b_valid, atol=atol, rtol=rtol))
            print(f"  {status} {label}.{col}: max_abs={max_abs:.2e} max_rel={max_rel:.2e} n_diff={n_diff}/{len(a_valid)}")
        else:
            print(f"  {status} {label}.{col}: max_abs={max_abs:.2e} max_rel={max_rel:.2e}")

    # Check string columns match exactly
    str_cols = df_base.select_dtypes(include=["object", "string"]).columns
    for col in str_cols:
        if not (df_base[col] == df_test[col]).all():
            n_diff = (df_base[col] != df_test[col]).sum()
            print(f"  FAIL {label}.{col}: {n_diff} string mismatches")
            all_pass = False
        else:
            print(f"  PASS {label}.{col}: exact match")

    return all_pass


def main():
    parser = argparse.ArgumentParser(description="Rigel regression check")
    parser.add_argument("baseline", help="Baseline output directory")
    parser.add_argument("test", help="Test output directory")
    parser.add_argument("--atol", type=float, default=1e-6, help="Absolute tolerance")
    parser.add_argument("--rtol", type=float, default=1e-6, help="Relative tolerance")
    args = parser.parse_args()

    files = ["quant.feather", "gene_quant.feather", "quant_detail.feather"]
    all_pass = True

    for f in files:
        print(f"\n=== {f} ===")
        base_path = f"{args.baseline}/{f}"
        test_path = f"{args.test}/{f}"
        ok = compare_feather(base_path, test_path, f, atol=args.atol, rtol=args.rtol)
        if not ok:
            all_pass = False

    print()
    if all_pass:
        print("REGRESSION CHECK: ALL PASS")
        return 0
    else:
        print("REGRESSION CHECK: FAILURES DETECTED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
