#!/usr/bin/env python
"""Trace pi_soft iterations through calibrate_gdna with both fixes."""
import sys
import numpy as np
sys.path.insert(0, "scripts/debug")

from calibration_stress_test import ScenarioConfig, generate_scenario, evaluate_scenario
from rigel.calibration import calibrate_gdna

cfg = ScenarioConfig(
    n_rna_regions=800, n_gdna_regions=200,
    strand_specificity=0.50, gdna_count_mean=30.0,
    gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42,
    label="low_gdna",
)
data = generate_scenario(cfg)
truth = data["truth_labels"]
is_gdna = truth == 1

result = calibrate_gdna(
    data["region_counts"], data["fl_table"], data["region_df"],
    0.50, diagnostics=True,
)

print(f"Final: pi={result.mixing_proportion:.6f}, "
      f"gDNA_gamma_mean={result.region_posteriors[is_gdna].mean():.4f}, "
      f"converged={result.converged}, n_iter={result.n_iterations}")

for i, h in enumerate(result.iteration_history):
    gamma = h["gamma"]
    gdna_gamma = gamma[is_gdna]
    print(f"  iter {i+1:2d}: pi={h['pi']:.4f}, pi_soft={h['pi_soft']:.4f}, "
          f"gDNA_gamma_mean={gdna_gamma.mean():.4f}, "
          f"gDNA_gamma>0.5: {(gdna_gamma > 0.5).sum()}/200, "
          f"delta_pi_soft={h['delta_pi']:.6f}")

ev = evaluate_scenario(result, truth, cfg, 0.0)
print(f"\nAUC={ev.auroc:.4f}, gF1={ev.gdna_f1:.4f}, rF1={ev.rna_f1:.4f}")
