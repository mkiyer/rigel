#!/usr/bin/env python
"""Quick check of scenario results with pi_soft + strand fixes."""
import sys
sys.path.insert(0, "scripts/debug")
from calibration_stress_test import ScenarioConfig, generate_scenario, evaluate_scenario
from rigel.calibration import calibrate_gdna

scenarios = [
    ("low_gdna SS=0.50 FL=identical", ScenarioConfig(
        n_rna_regions=800, n_gdna_regions=200, strand_specificity=0.50,
        gdna_count_mean=30.0, gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42)),
    ("med_gdna SS=0.50 FL=identical", ScenarioConfig(
        n_rna_regions=600, n_gdna_regions=400, strand_specificity=0.50,
        gdna_count_mean=60.0, gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42)),
    ("equal_gdna SS=0.50 (same density)", ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500, strand_specificity=0.50,
        gdna_count_mean=200.0, gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42)),
    ("equal_gdna_lower SS=0.50", ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500, strand_specificity=0.50,
        gdna_count_mean=50.0, gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42)),
    ("low_gdna SS=0.95", ScenarioConfig(
        n_rna_regions=800, n_gdna_regions=200, strand_specificity=0.95,
        gdna_count_mean=30.0, gdna_fl_mean=250.0, gdna_fl_std=50.0, seed=42)),
]

for label, cfg in scenarios:
    data = generate_scenario(cfg)
    result = calibrate_gdna(
        data["region_counts"], data["fl_table"],
        data["region_df"], cfg.strand_specificity,
    )
    ev = evaluate_scenario(result, data["truth_labels"], cfg, 0.0)
    print(
        f"{label:40s}  AUC={ev.auroc:.4f}  gF1={ev.gdna_f1:.4f}  "
        f"rF1={ev.rna_f1:.4f}  pi_err={ev.pi_error:+.4f}  "
        f"gdna_post={ev.mean_gdna_posterior:.4f}"
    )
