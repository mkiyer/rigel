#!/usr/bin/env python
"""Analyze diagnostic sweep results."""
import json
from pathlib import Path

base = Path("/Users/mkiyer/Downloads/rigel_runs")

def load(name):
    with open(base / name / "sweep_results.json") as f:
        return json.load(f)

print("=" * 80)
print("gDNA TITRATION (TA1=100, nRNA=0, SS=1.0)")
print("=" * 80)
for r in load("diag_gdna_titration"):
    print(
        f"gdna={r['gdna']:>3}  "
        f"TA1_exp={r['TA1_expected']:>5}  TA1_obs={r['TA1_observed']:>7.1f}  "
        f"TA1_rel={r['TA1_rel_err']*100:>5.1f}%  "
        f"gdna_exp={r['gdna_expected']:>5}  gdna_obs={r['gdna_observed']:>8.1f}  "
        f"gdna_diff={r['gdna_abs_diff']:>6.1f}"
    )

print()
print("=" * 80)
print("STRAND SPECIFICITY (TA1=100, nRNA=100, gDNA=100)")
print("=" * 80)
for r in load("diag_ss_sweep"):
    print(
        f"ss={r['strand_specificity']:.2f}  "
        f"TA1_exp={r['TA1_expected']:>4}  TA1_obs={r['TA1_observed']:>7.1f}  "
        f"TA1_rel={r['TA1_rel_err']*100:>5.1f}%  "
        f"nrna_exp={r['nrna_expected']:>5}  nrna_obs={r['nrna_observed']:>7.1f}  "
        f"nrna_diff={r['nrna_abs_diff']:>6.1f}  "
        f"gdna_exp={r['gdna_expected']:>5}  gdna_obs={r['gdna_observed']:>7.1f}  "
        f"gdna_diff={r['gdna_abs_diff']:>6.1f}"
    )

print()
print("=" * 80)
print("FRAGMENT COUNT (TA1=100, nRNA=100, gDNA=100, SS=1.0)")
print("=" * 80)
for r in load("diag_frag_count"):
    print(
        f"n_frag={int(r['n_fragments']):>5}  "
        f"TA1_exp={r['TA1_expected']:>5}  TA1_obs={r['TA1_observed']:>7.1f}  "
        f"TA1_rel={r['TA1_rel_err']*100:>5.1f}%  "
        f"nrna_diff={r['nrna_abs_diff']:>6.1f}  "
        f"gdna_diff={r['gdna_abs_diff']:>6.1f}  "
        f"mrna_rel={r['total_mrna_rel_err']*100:>5.1f}%"
    )

print()
print("=" * 80)
print("PRUNE THRESHOLD (TA1=100, TA4=3, nRNA=100, gDNA=100, SS=1.0)")
print("=" * 80)
for r in load("diag_prune_threshold"):
    pt = r.get("prune_threshold", "?")
    print(
        f"pt={pt}  "
        f"TA1_obs={r['TA1_observed']:>7.1f}  "
        f"TA4_exp={r['TA4_expected']:>3}  TA4_obs={r['TA4_observed']:>6.1f}  "
        f"TA4_rel={r['TA4_rel_err']*100:>5.1f}%  "
        f"gdna_diff={r['gdna_abs_diff']:>6.1f}  "
        f"mrna_rel={r['total_mrna_rel_err']*100:>5.1f}%"
    )
