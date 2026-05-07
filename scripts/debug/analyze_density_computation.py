#!/usr/bin/env python3
"""Analyze intronic density computation and effective lengths."""

import json
import numpy as np
import pandas as pd
from pathlib import Path

base = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")
index_dir = base / "rigel_index"

# Load calibration summary
with open(base / "gdna_high_ss_0.99_nrna_none/rigel_out/summary.json") as f:
    s = json.load(f)

cal = s["calibration"]
gd = cal["global_densities"]
gdna_fl_mean = gd["gdna_fl_mean"]

print(f"gDNA FL mean: {gdna_fl_mean:.2f}")
print()
for rtype in ["INTERGENIC", "INTRON", "EXON-INTRON"]:
    d = gd[rtype]
    print(f"{rtype}:")
    print(f"  rho = {d['rho']:.8f}")
    print(f"  n_fragments = {d['n_fragments']:,}")
    print(f"  eff_length_bp = {d['eff_length_bp']:,.0f}")
    print(f"  n_regions_used = {d['n_regions_used']}")
    print(f"  kappa = {d['kappa']:.2f}")
    print()

rho_ig = gd["INTERGENIC"]["rho"]
rho_in = gd["INTRON"]["rho"]
rho_ex = gd["EXON-INTRON"]["rho"]
print(f"Ratios:")
print(f"  rho_intron / rho_intergenic = {rho_in / rho_ig:.4f}")
print(f"  rho_exon / rho_intergenic   = {rho_ex / rho_ig:.4f}")
print(f"  rho_exon / rho_intron       = {rho_ex / rho_in:.4f}")
print()

# ── Load region_df from index to understand effective lengths ──
region_df = pd.read_feather(index_dir / "regions.feather")
print(f"Total regions: {len(region_df)}")
print(f"Region type distribution:")
type_map = {0: "INTERGENIC", 1: "INTRON", 2: "EXON"}
for t in [0, 1, 2]:
    mask = region_df["type"] == t
    subset = region_df[mask]
    spans = subset["end"] - subset["start"]
    leffs = spans + gdna_fl_mean - 1
    print(f"  {type_map[t]:>12}: n={mask.sum():>5}, "
          f"total_span={spans.sum():>12,}, "
          f"total_leff={leffs.sum():>12,.0f}, "
          f"mean_span={spans.mean():>8.0f}, "
          f"min_span={spans.min():>6}, "
          f"max_span={spans.max():>8}")

print()
print("=" * 80)
print("DETAILED INTRONIC ANALYSIS")
print("=" * 80)
print()

# Focus on introns
introns = region_df[region_df["type"] == 1].copy()
introns["span"] = introns["end"] - introns["start"]
introns["leff_overlap"] = introns["span"] + gdna_fl_mean - 1  # l_eff_overlap formula

# The user's point: for a truly intronic fragment to NOT overlap any exon,
# it needs start position such that [start, start+fl) is fully within the intron.
# Valid start positions: 0 to (intron_span - fl).
# Number of valid starts = max(0, intron_span - fl + 1).
# But the l_eff_overlap formula gives span + fl_mean - 1, which is the NUMBER
# OF GENOMIC POSITIONS where a fragment of length fl_mean can have its
# midpoint/start and still overlap the region.
#
# Wait - l_eff_overlap = span + fl - 1. This is the number of genomic positions
# where a fragment of given length overlaps with the region.
# For CONTAINMENT (fragment fully inside intron): valid starts = span - fl + 1.
#
# The calibration uses l_eff_overlap because it counts fragments that OVERLAP
# the intron region. But what gets COUNTED as intronic? Only fragments fully
# CONTAINED within the intron! Fragments that overlap both intron and exon
# get classified differently (as exon-intron boundary or exonic).

print("The l_eff_overlap formula: leff = span + fl_mean - 1")
print("This counts fragments that OVERLAP the region from any side.")
print()
print("But INTRONIC fragments are classified as those CONTAINED within the intron.")
print("For containment: valid_starts = max(0, span - fl + 1)")
print()
print("This creates a MISMATCH: the denominator (leff = span + fl - 1) is TOO LARGE")
print("relative to the numerator (fragments fully contained in intron).")
print()

# Calculate the discrepancy
introns["leff_containment"] = np.maximum(0, introns["span"] - gdna_fl_mean + 1)
introns["ratio"] = introns["leff_containment"] / introns["leff_overlap"]

print(f"Intron statistics:")
print(f"  Mean span:               {introns['span'].mean():.0f} bp")
print(f"  Median span:             {introns['span'].median():.0f} bp")
print(f"  Mean l_eff_overlap:      {introns['leff_overlap'].mean():.0f}")
print(f"  Mean l_eff_containment:  {introns['leff_containment'].mean():.0f}")
print(f"  Mean ratio (contain/overlap): {introns['ratio'].mean():.4f}")
print(f"  n_introns where span < fl_mean: {(introns['span'] < gdna_fl_mean).sum()}")
print()

# The expected density ratio if using l_eff_overlap but counting contained frags:
# Observed: n_contained_fragments / leff_overlap
# True:     n_contained_fragments / leff_containment
# Ratio of observed/true = leff_containment / leff_overlap < 1
# This predicts rho_intron < rho_intergenic by this ratio
weighted_ratio = (introns["leff_containment"].sum() / introns["leff_overlap"].sum())
print(f"Predicted rho_intron / rho_intergenic:")
print(f"  = sum(leff_containment) / sum(leff_overlap)")
print(f"  = {introns['leff_containment'].sum():,.0f} / {introns['leff_overlap'].sum():,.0f}")
print(f"  = {weighted_ratio:.4f}")
print()
print(f"Actual observed: rho_intron / rho_intergenic = {rho_in / rho_ig:.4f}")
print()

# For intergenic regions (typically much larger), the ratio is ~1.0
igenic = region_df[region_df["type"] == 0].copy()
igenic["span"] = igenic["end"] - igenic["start"]
igenic["leff_overlap"] = igenic["span"] + gdna_fl_mean - 1
igenic["leff_containment"] = np.maximum(0, igenic["span"] - gdna_fl_mean + 1)
ig_ratio = igenic["leff_containment"].sum() / igenic["leff_overlap"].sum()
print(f"Intergenic containment/overlap ratio: {ig_ratio:.6f}")
print(f"  (close to 1.0 because intergenic regions are large)")
print()

# Show distribution of intron spans
print("Intron span distribution:")
bins = [0, 100, 200, 350, 500, 1000, 2000, 5000, 10000, 50000]
for i in range(len(bins) - 1):
    mask = (introns["span"] >= bins[i]) & (introns["span"] < bins[i+1])
    n = mask.sum()
    if n > 0:
        sub = introns[mask]
        print(f"  [{bins[i]:>5}, {bins[i+1]:>5}): n={n:>4}, "
              f"mean_leff_overlap={sub['leff_overlap'].mean():>8.0f}, "
              f"mean_leff_contain={sub['leff_containment'].mean():>8.0f}, "
              f"mean_ratio={sub['ratio'].mean():.3f}")

print()
print("=" * 80)
print("USER's SUGGESTION: Per-fragment effective length")
print("=" * 80)
print()
print("Currently: density = n_fragments / sum(leff_overlap)")
print("  where leff_overlap uses a SINGLE fl_mean for all regions.")
print()
print("Proposed: density = sum(1) / sum(leff_per_fragment)")
print("  where leff_per_fragment = max(0, region_span - actual_frag_len + 1)")
print("  for each fragment, using its ACTUAL observed length.")
print()
print("This correctly handles:")
print("  1. Non-normal FL distributions (heavy tails)")
print("  2. Short introns where only short fragments fit")
print("  3. Per-fragment weight = 1/leff gives proper density estimate")
print()
print("Expected improvement: eliminates the systematic underestimate")
print("of intronic density caused by using mean FL in the denominator")
print("when actual contained fragments are biased toward shorter FLs.")
