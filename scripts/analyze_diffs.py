"""Analyze magnitude of differences between rigel runs."""
import pyarrow.feather as feather
import numpy as np

b = '/Users/mkiyer/Downloads/rigel_runs/baseline_output'
p1 = '/Users/mkiyer/Downloads/rigel_runs/phase1_output'
p2 = '/Users/mkiyer/Downloads/rigel_runs/phase1_output_run2'

pairs = [
    ('baseline vs phase1', b, p1),
    ('phase1 vs phase1_r2', p1, p2),
]

for name, d1, d2 in pairs:
    df1 = feather.read_table(f'{d1}/quant.feather').to_pandas()
    df2 = feather.read_table(f'{d2}/quant.feather').to_pandas()
    num_cols = [c for c in df1.columns if df1[c].dtype in [np.float64, np.float32, np.int64, np.int32]]
    for col in num_cols:
        v1 = df1[col].to_numpy(dtype=float)
        v2 = df2[col].to_numpy(dtype=float)
        diff = np.abs(v1 - v2)
        max_diff = diff.max()
        mean_diff = diff.mean()
        max_val = max(np.abs(v1).max(), np.abs(v2).max())
        n_nonzero = np.sum(diff > 0)
        n_total = len(v1)
        rel_max = max_diff / max_val if max_val > 0 else 0
        print(f'{name} | {col}: max_abs_diff={max_diff:.6g}, mean_abs_diff={mean_diff:.6g}, '
              f'max_val={max_val:.6g}, rel_max_diff={rel_max:.6g}, '
              f'n_differ={n_nonzero}/{n_total}')
    print()

# Also check gene-level
print("=== Gene-level ===")
for name, d1, d2 in pairs:
    df1 = feather.read_table(f'{d1}/gene_quant.feather').to_pandas()
    df2 = feather.read_table(f'{d2}/gene_quant.feather').to_pandas()
    num_cols = [c for c in df1.columns if df1[c].dtype in [np.float64, np.float32, np.int64, np.int32]]
    for col in num_cols:
        v1 = df1[col].to_numpy(dtype=float)
        v2 = df2[col].to_numpy(dtype=float)
        diff = np.abs(v1 - v2)
        max_diff = diff.max()
        max_val = max(np.abs(v1).max(), np.abs(v2).max())
        n_nonzero = np.sum(diff > 0)
        rel_max = max_diff / max_val if max_val > 0 else 0
        print(f'{name} | {col}: max_abs_diff={max_diff:.6g}, max_val={max_val:.6g}, '
              f'rel_max_diff={rel_max:.6g}, n_differ={n_nonzero}/{len(v1)}')
    print()

# Check locus-level
print("=== Locus-level ===")
for name, d1, d2 in pairs:
    df1 = feather.read_table(f'{d1}/loci.feather').to_pandas()
    df2 = feather.read_table(f'{d2}/loci.feather').to_pandas()
    for col in df1.select_dtypes(include=[np.number]).columns:
        v1 = df1[col].to_numpy(dtype=float)
        v2 = df2[col].to_numpy(dtype=float)
        diff = np.abs(v1 - v2)
        max_diff = diff.max()
        if max_diff > 0:
            max_val = max(np.abs(v1).max(), np.abs(v2).max())
            rel_max = max_diff / max_val if max_val > 0 else 0
            n_nonzero = np.sum(diff > 0)
            print(f'{name} | {col}: max_abs_diff={max_diff:.6g}, max_val={max_val:.6g}, '
                  f'rel_max_diff={rel_max:.6g}, n_differ={n_nonzero}/{len(v1)}')
    print()
