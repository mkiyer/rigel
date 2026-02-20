#!/usr/bin/env python
"""Analyze diagnostic results for gDNA/nRNA confusion patterns."""
import json
import numpy as np

with open("docs/two_exon_ctrl_diagnostic.json") as f:
    data = json.load(f)

# Sort for readable output
data_sorted = sorted(data, key=lambda x: (
    x["strand_specificity"], x["gdna_abundance"], x["nrna_fraction"]))

# === 1. gDNA + nRNA confusion (both > 0) ===
print("=== gDNA + nRNA confusion (both > 0) ===")
print(f"{'ss':>5} {'gdna':>4} {'nrna':>4} | "
      f"{'g_exp':>5} {'g_pip':>5} {'g_err':>6} | "
      f"{'n_exp':>5} {'n_pip':>5} {'n_err':>6} | "
      f"{'t1_err':>6} {'g_init':>6} {'n_init':>6}")
for r in data_sorted:
    if r["gdna_abundance"] > 0 and r["nrna_fraction"] > 0:
        print(f"{r['strand_specificity']:5.2f} {r['gdna_abundance']:4d} "
              f"{r['nrna_fraction']:4.1f} | "
              f"{r['gdna_expected']:5d} {r['gdna_pipeline']:5.0f} "
              f"{r['gdna_signed_err']:+6.0f} | "
              f"{r['nrna_expected']:5d} {r['nrna_pipeline']:5.0f} "
              f"{r['nrna_signed_err']:+6.0f} | "
              f"{r['t1_signed_err']:+6.1f} {r['gdna_init_total']:6.0f} "
              f"{r['nrna_init_total']:6.0f}")

# === 2. SS=0.5 nRNA collapse (nrna > 0, no gDNA) ===
print("\n=== SS=0.5 nRNA collapse (nrna > 0, gdna=0) ===")
for r in data_sorted:
    if r["gdna_abundance"] == 0 and r["nrna_fraction"] > 0 and r["strand_specificity"] == 0.5:
        print(f"  nrna_frac={r['nrna_fraction']:.1f}: "
              f"n_exp={r['nrna_expected']:3d} n_pip={r['nrna_pipeline']:5.0f} "
              f"gdna_pip={r['gdna_pipeline']:5.0f} "
              f"g_init={r['gdna_init_total']:.0f} n_init={r['nrna_init_total']:.0f} "
              f"intronic_s={r.get('intronic_sense', 0):.0f} intronic_a={r.get('intronic_antisense', 0):.0f}"
              f" sense_all={r.get('gene_sense_all', 0):.0f} anti_all={r.get('gene_antisense_all', 0):.0f}")

# === 3. Key question: are gDNA+nRNA errors antisymmetric? ===
print("\n=== Error antisymmetry check (gdna_err ≈ -nrna_err?) ===")
for r in data_sorted:
    if r["gdna_abundance"] > 0 and r["nrna_fraction"] > 0:
        g = r["gdna_signed_err"]
        n = r["nrna_signed_err"]
        residual = g + n  # if antisymmetric, this should be ~0
        if abs(residual) > 5:
            print(f"  ss={r['strand_specificity']:.2f} gdna={r['gdna_abundance']:3d} "
                  f"nrna={r['nrna_fraction']:.1f}: "
                  f"g_err={g:+.0f} n_err={n:+.0f} residual={residual:+.0f}")

# === 4. Worst t2 false positives ===
print("\n=== t2 false positives > 0 ===")
for r in data_sorted:
    if r["t2_observed"] > 0.5:
        print(f"  ss={r['strand_specificity']:.2f} gdna={r['gdna_abundance']:3d} "
              f"nrna={r['nrna_fraction']:.1f}: "
              f"t2_fp={r['t2_observed']:.0f}")

# === 5. Aggregate by SS for gDNA-only runs ===
print("\n=== gDNA accuracy by SS (nRNA=0, gDNA > 0) ===")
for ss in [0.5, 0.75, 0.9, 0.95, 1.0]:
    subset = [r for r in data if r["strand_specificity"] == ss
              and r["gdna_abundance"] > 0 and r["nrna_fraction"] == 0.0]
    if subset:
        errs = [r["gdna_signed_err"] for r in subset]
        print(f"  SS={ss:.2f}: gDNA errors = {[f'{e:+.0f}' for e in errs]}")

# === 6. nRNA accuracy at SS>=0.9 with no gDNA ===
print("\n=== nRNA accuracy at SS>=0.9, gDNA=0 ===")
for r in data_sorted:
    if r["gdna_abundance"] == 0 and r["nrna_fraction"] > 0 and r["strand_specificity"] >= 0.9:
        print(f"  ss={r['strand_specificity']:.2f} nrna_frac={r['nrna_fraction']:.1f}: "
              f"n_exp={r['nrna_expected']:3d} n_pip={r['nrna_pipeline']:5.0f} "
              f"n_err={r['nrna_signed_err']:+.0f} t1_err={r['t1_signed_err']:+.1f}")
