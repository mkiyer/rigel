#!/usr/bin/env python
"""Analyze calibration validation sweep results."""
import csv
import sys
from collections import defaultdict

import numpy as np

tsv_path = sys.argv[1] if len(sys.argv) > 1 else "sweep_results/cal_validation/sweep_results.tsv"

with open(tsv_path) as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

print(f"Total runs: {len(rows)}")

# Parse numeric fields
def num(v, default=np.nan):
    try:
        return float(v) if v not in ("", None) else default
    except (ValueError, TypeError):
        return default

# ── 1. Calibration kappa recovery ────────────────────────────
print("\n" + "="*80)
print("1. KAPPA RECOVERY (κ_true vs κ_est)")
print("="*80)

# Group by (gdna_fraction, gdna_strand_kappa)
kappa_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    kt = num(r["cal_kappa_true"])
    ke = num(r["cal_kappa_est"])
    kappa_groups[(gf, kt)].append(ke)

print(f"\n{'gF':>5} {'κ_true':>7} {'κ_est_mean':>10} {'κ_est_std':>9} {'n':>3}")
print("-"*40)
for (gf, kt), vals in sorted(kappa_groups.items()):
    vals = np.array(vals)
    valid = vals[~np.isnan(vals)]
    if len(valid) > 0:
        print(f"{gf:5.1f} {kt:7.0f} {valid.mean():10.1f} {valid.std():9.1f} {len(valid):3d}")
    else:
        print(f"{gf:5.1f} {kt:7.0f} {'N/A':>10} {'N/A':>9} {len(vals):3d}")

# ── 2. Calibration density estimation ────────────────────────
print("\n" + "="*80)
print("2. DENSITY ESTIMATION")
print("="*80)

den_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    ss = num(r["strand_specificity"])
    den = num(r["cal_density_est"])
    den_groups[(gf, ss)].append(den)

print(f"\n{'gF':>5} {'SS':>5} {'density_mean':>12} {'density_std':>12} {'n':>3}")
print("-"*45)
for (gf, ss), vals in sorted(den_groups.items()):
    vals = np.array(vals)
    valid = vals[~np.isnan(vals)]
    if len(valid) > 0:
        print(f"{gf:5.1f} {ss:5.2f} {valid.mean():12.4e} {valid.std():12.4e} {len(valid):3d}")

# ── 3. Fragment length recovery ──────────────────────────────
print("\n" + "="*80)
print("3. gDNA FRAGMENT LENGTH RECOVERY")
print("="*80)

fl_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    fm_true = num(r["cal_gdna_fl_true_mean"])
    fm_est = num(r["cal_gdna_fl_est_mean"])
    fl_groups[(gf, fm_true)].append(fm_est)

print(f"\n{'gF':>5} {'FL_true':>7} {'FL_est_mean':>11} {'FL_err_mean':>11} {'n':>3}")
print("-"*45)
for (gf, fmt), vals in sorted(fl_groups.items()):
    vals = np.array(vals)
    valid = vals[~np.isnan(vals)]
    if len(valid) > 0:
        err = valid - fmt
        print(f"{gf:5.1f} {fmt:7.0f} {valid.mean():11.1f} {err.mean():+11.1f} {len(valid):3d}")

# ── 4. Seed selection ────────────────────────────────────────
print("\n" + "="*80)
print("4. SEED REGION SELECTION")
print("="*80)

seed_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    ss = num(r["strand_specificity"])
    ns = num(r["cal_n_seed"])
    nr = num(r["cal_n_regions"])
    seed_groups[(gf, ss)].append((ns, nr))

print(f"\n{'gF':>5} {'SS':>5} {'seed_mean':>9} {'seed_std':>9} {'n_reg':>5}")
print("-"*40)
for (gf, ss), vals in sorted(seed_groups.items()):
    seeds = np.array([v[0] for v in vals])
    nregs = np.array([v[1] for v in vals])
    valid = seeds[~np.isnan(seeds)]
    vr = nregs[~np.isnan(nregs)]
    if len(valid) > 0:
        print(f"{gf:5.1f} {ss:5.2f} {valid.mean():9.1f} {valid.std():9.1f} {vr.mean():5.0f}")

# ── 5. gDNA quantification accuracy ─────────────────────────
print("\n" + "="*80)
print("5. gDNA QUANTIFICATION ACCURACY")
print("="*80)

gdna_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    ss = num(r["strand_specificity"])
    ge = num(r["gdna_expected"])
    go = num(r["gdna_observed"])
    gdna_groups[(gf, ss)].append((ge, go))

print(f"\n{'gF':>5} {'SS':>5} {'exp_mean':>9} {'obs_mean':>9} {'err_mean':>9} {'rel_err':>7}")
print("-"*50)
for (gf, ss), vals in sorted(gdna_groups.items()):
    exps = np.array([v[0] for v in vals])
    obss = np.array([v[1] for v in vals])
    diff = obss - exps
    rel = np.where(exps > 0, np.abs(diff) / exps, 0)
    print(f"{gf:5.1f} {ss:5.2f} {exps.mean():9.0f} {obss.mean():9.0f} {diff.mean():+9.0f} {rel.mean():7.1%}")

# ── 6. Convergence ───────────────────────────────────────────
print("\n" + "="*80)
print("6. CONVERGENCE STATISTICS")
print("="*80)

conv_true = sum(1 for r in rows if r.get("cal_converged") == "True")
conv_false = sum(1 for r in rows if r.get("cal_converged") == "False")
conv_na = len(rows) - conv_true - conv_false

iters = [num(r["cal_n_iterations"]) for r in rows]
iters = np.array([i for i in iters if not np.isnan(i)])

print(f"Converged: {conv_true}, Not converged: {conv_false}, N/A: {conv_na}")
if len(iters) > 0:
    print(f"Iterations: mean={iters.mean():.1f}, min={iters.min():.0f}, max={iters.max():.0f}")

# ── 7. Mean weight distribution ──────────────────────────────
print("\n" + "="*80)
print("7. MEAN WEIGHT by gdna_fraction × SS")
print("="*80)

mw_groups = defaultdict(list)
for r in rows:
    gf = num(r["gdna_fraction"])
    ss = num(r["strand_specificity"])
    mw = num(r["cal_mean_weight"])
    mw_groups[(gf, ss)].append(mw)

print(f"\n{'gF':>5} {'SS':>5} {'mw_mean':>8} {'mw_std':>8}")
print("-"*30)
for (gf, ss), vals in sorted(mw_groups.items()):
    vals = np.array(vals)
    valid = vals[~np.isnan(vals)]
    if len(valid) > 0:
        print(f"{gf:5.1f} {ss:5.2f} {valid.mean():8.3f} {valid.std():8.3f}")

# ── 8. Kappa recovery by SS (with gDNA present) ─────────────
print("\n" + "="*80)
print("8. KAPPA RECOVERY by SS (gDNA > 0 only)")
print("="*80)

for ss in [0.5, 0.75, 0.9, 1.0]:
    print(f"\n  SS={ss}:")
    kappa_ss = defaultdict(list)
    for r in rows:
        if num(r["gdna_fraction"]) == 0:
            continue
        if abs(num(r["strand_specificity"]) - ss) > 0.01:
            continue
        kt = num(r["cal_kappa_true"])
        ke = num(r["cal_kappa_est"])
        gf = num(r["gdna_fraction"])
        kappa_ss[(kt, gf)].append(ke)

    print(f"  {'κ_true':>7} {'gF':>5} {'κ_est':>8} {'err':>8}")
    print(f"  {'-'*30}")
    for (kt, gf), vals in sorted(kappa_ss.items()):
        vals = np.array(vals)
        valid = vals[~np.isnan(vals)]
        if len(valid) > 0:
            print(f"  {kt:7.0f} {gf:5.1f} {valid.mean():8.1f} {(valid.mean()-kt):+8.1f}")

print("\n" + "="*80)
print("ANALYSIS COMPLETE")
print("="*80)
