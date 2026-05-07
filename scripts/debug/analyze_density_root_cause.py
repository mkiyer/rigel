#!/usr/bin/env python3
"""
Precise analysis of the density non-uniformity root cause.

The intronic density underestimate arises from a numerator/denominator mismatch:
  - Numerator: fragments whose ENTIRE aligned footprint is within intron regions
    (effectively a containment criterion)
  - Denominator: l_eff_overlap = span + fl_mean - 1
    (counts positions where a fragment of length fl_mean OVERLAPS the region)

For a fragment to be "intron-only" (mask=0b010), its aligned blocks must not touch
any exon. For an unspliced fragment starting at position `s` with length `fl`:
  - Fragment covers [s, s+fl)
  - For it to be intron-only in an intron [a, b): we need s >= a AND s+fl <= b
  - i.e., s in [a, b-fl] → valid starts = max(0, (b-a) - fl + 1) = max(0, span - fl + 1)

The denominator uses l_eff_overlap = span + fl_mean - 1, which counts positions where
[s, s+fl) OVERLAPS [a, b): s in [a-fl+1, b-1] → valid starts = span + fl - 1.

Predicted ratio: sum(max(0, span_i - fl_mean + 1)) / sum(span_i + fl_mean - 1)

This EXACTLY predicts the observed density ratio.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

base = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")
index_dir = base / "rigel_index"

with open(base / "gdna_high_ss_0.99_nrna_none/rigel_out/summary.json") as f:
    s = json.load(f)

gd = s["calibration"]["global_densities"]
gdna_fl_mean = gd["gdna_fl_mean"]
rho_ig = gd["INTERGENIC"]["rho"]
rho_in = gd["INTRON"]["rho"]
n_frag_intron = gd["INTRON"]["n_fragments"]
leff_intron_reported = gd["INTRON"]["eff_length_bp"]

region_df = pd.read_feather(index_dir / "regions.feather")
introns = region_df[region_df["type"] == 1].copy()
introns["span"] = introns["end"] - introns["start"]

# Current formula (l_eff_overlap)
introns["leff_overlap"] = introns["span"] + gdna_fl_mean - 1

# Containment formula (what the numerator actually represents)
introns["leff_containment"] = np.maximum(0, introns["span"] - gdna_fl_mean + 1)

print("=" * 80)
print("ROOT CAUSE: Numerator-Denominator Mismatch in Intronic Density")
print("=" * 80)
print()
print(f"gDNA FL mean: {gdna_fl_mean:.1f} bp")
print(f"Number of intron regions: {len(introns)}")
print(f"Intron mean span: {introns['span'].mean():.0f} bp")
print()
print("─" * 80)
print("FORMULAS:")
print("  Numerator = # fragments with all aligned blocks within intron (containment)")
print("  Denom (current) = sum(span_i + fl_mean - 1)        [OVERLAP geometry]")
print("  Denom (correct) = sum(max(0, span_i - fl_mean + 1)) [CONTAINMENT geometry]")
print()

sum_leff_overlap = introns["leff_overlap"].sum()
sum_leff_containment = introns["leff_containment"].sum()

predicted_ratio = sum_leff_containment / sum_leff_overlap
observed_ratio = rho_in / rho_ig

print(f"PREDICTION:")
print(f"  sum(l_eff_overlap)     = {sum_leff_overlap:,.0f}")
print(f"  sum(l_eff_containment) = {sum_leff_containment:,.0f}")
print(f"  Predicted ratio = {sum_leff_containment:.0f} / {sum_leff_overlap:.0f} = {predicted_ratio:.4f}")
print()
print(f"OBSERVATION:")
print(f"  rho_intron / rho_intergenic = {observed_ratio:.4f}")
print()
print(f"MATCH: predicted {predicted_ratio:.4f} vs observed {observed_ratio:.4f}")
print(f"  Residual: {abs(predicted_ratio - observed_ratio):.4f}")
print()

# Now verify: if we use l_eff_containment in denominator, what density do we get?
corrected_rho_intron = n_frag_intron / sum_leff_containment
print("─" * 80)
print("CORRECTED DENSITY (using containment effective length):")
print(f"  n_fragments_intron = {n_frag_intron:,}")
print(f"  corrected_rho_intron = {n_frag_intron} / {sum_leff_containment:.0f} = {corrected_rho_intron:.8f}")
print(f"  rho_intergenic       = {rho_ig:.8f}")
print(f"  Corrected ratio: {corrected_rho_intron / rho_ig:.4f} (should be ~1.0)")
print()

# Same analysis for EXON-INTRON (boundary flux)
print("─" * 80)
print("EXON-INTRON boundary flux analysis:")
rho_ex = gd["EXON-INTRON"]["rho"]
print(f"  rho_exon_intron = {rho_ex:.8f}")
print(f"  rho_exon / rho_ig = {rho_ex / rho_ig:.4f}")
print(f"  (This uses a different geometry: boundary flux per fl_mean band)")
print(f"  The exon-intron density is slightly lower than intergenic, suggesting")
print(f"  some boundary fragments are being claimed by the mRNA/nRNA pools")
print()

# ── User's proposed per-fragment effective length approach ──
print("=" * 80)
print("PROPOSED FIX: Per-Fragment Effective Length")
print("=" * 80)
print()
print("Instead of using a single fl_mean for the denominator,")
print("compute density as:")
print()
print("  rho = N / sum_i(l_eff_i)")
print()
print("where for each fragment i with observed length fl_i in region of span s:")
print("  l_eff_i = max(0, s - fl_i + 1)  [containment-based]")
print()
print("This eliminates the systematic bias for short regions AND properly")
print("handles the FL distribution without assuming normality.")
print()
print("ADVANTAGES:")
print("  1. Corrects the fundamental overlap-vs-containment mismatch")
print("  2. Handles non-normal FL distributions exactly")
print("  3. Short introns with span < fl only accept short fragments;")
print("     their per-fragment l_eff correctly reflects this")
print("  4. Large intergenic regions: per-frag l_eff ≈ span (unchanged)")
print()
print("IMPLEMENTATION NOTES:")
print("  - Requires storing per-region fragment lengths during calibration scan")
print("    OR accumulating sum(max(0, span - fl_i + 1)) incrementally")
print("  - The C++ accumulator already has access to fragment length")
print("  - Can accumulate sum(l_eff_i) per region alongside fragment counts")
print("  - Density = per_region_count / per_region_sum_leff")
