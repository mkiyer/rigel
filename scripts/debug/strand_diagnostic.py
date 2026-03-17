#!/usr/bin/env python
"""Quick strand channel diagnostic for unstranded data."""
import sys
import numpy as np
sys.path.insert(0, "scripts/debug")

from calibration_stress_test import ScenarioConfig, generate_scenario
from rigel.calibration import (
    compute_region_stats, compute_sense_fraction,
    build_strand_histogram, _compute_strand_llr_binomial_mixture,
)

cfg = ScenarioConfig(
    n_rna_regions=800, n_gdna_regions=200,
    strand_specificity=0.50, gdna_count_mean=30.0,
    seed=42, label="test",
)
data = generate_scenario(cfg)
stats = compute_region_stats(data["region_counts"], data["region_df"])
sf = compute_sense_fraction(stats)
truth = data["truth_labels"]
n = len(truth)

# Build seed gamma
gamma = np.full(n, 0.5)
gamma[stats["n_spliced"] > 0] = 0.0
gamma[truth == 1] = 1.0
gamma[stats["n_spliced"] > 0] = 0.0

is_gdna = truth == 1
soft = (stats["n_spliced"] == 0) & (stats["n_total"] > 0) & (stats["region_length"] > 0)

# --- Default: gDNA symmetric, RNA not symmetric ---
gdna_h = build_strand_histogram(sf, gamma, symmetric=True)
rna_h = build_strand_histogram(sf, 1.0 - gamma, symmetric=False)

n_bins = len(gdna_h[1])
mid = n_bins // 2
print("=== gDNA histogram (symmetric=True) ===")
print(f"  peak density: {gdna_h[1].max():.4f} at bin {np.argmax(gdna_h[1])}")
print(f"  center density: {gdna_h[1][mid]:.4f}")

print("=== RNA histogram (symmetric=False) ===")
print(f"  peak density: {rna_h[1].max():.4f} at bin {np.argmax(rna_h[1])}")
print(f"  center density: {rna_h[1][mid]:.4f}")

llr = _compute_strand_llr_binomial_mixture(stats, gdna_h, rna_h, n)
print(f"\nStrand LLR for gDNA: mean={llr[is_gdna & soft].mean():.3f}, "
      f"std={llr[is_gdna & soft].std():.3f}")

# --- Both symmetric ---
rna_h2 = build_strand_histogram(sf, 1.0 - gamma, symmetric=True)
llr2 = _compute_strand_llr_binomial_mixture(stats, gdna_h, rna_h2, n)
print(f"Strand LLR (both sym): mean={llr2[is_gdna & soft].mean():.3f}, "
      f"std={llr2[is_gdna & soft].std():.3f}")

# --- Key insight: histogram widths ---
print(f"\ngDNA sf: mean={sf[is_gdna].mean():.3f}, std={sf[is_gdna & np.isfinite(sf)].std():.3f}")
is_rna = truth == 0
print(f"RNA  sf: mean={sf[is_rna].mean():.3f}, std={sf[is_rna & np.isfinite(sf)].std():.3f}")
print(f"\ngDNA reads/region: ~{stats['n_total'][is_gdna].mean():.0f}")
print(f"RNA  reads/region: ~{stats['n_total'][is_rna].mean():.0f}")
