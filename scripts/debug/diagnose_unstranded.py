#!/usr/bin/env python
"""Diagnose posterior collapse for SS=0.50 (unstranded) scenarios.

Traces through calibration EM step-by-step to find where posteriors
collapse when only density should be discriminating.
"""
from __future__ import annotations

import numpy as np
import sys
sys.path.insert(0, "scripts/debug")

from calibration_stress_test import ScenarioConfig, generate_scenario, evaluate_scenario

from rigel.calibration import (
    calibrate_gdna,
    compute_region_stats,
    compute_log_density,
    compute_sense_fraction,
    build_density_histogram,
    build_strand_histogram,
    build_gdna_fl_model,
    _seed_initial_partition,
    _e_step,
    _m_step,
    _compute_density_llr,
    _compute_strand_llr_binomial_mixture,
    _DENSITY_BINS,
)


def trace_scenario(label: str, cfg: ScenarioConfig):
    """Full EM trace for a single scenario."""
    print(f"\n{'='*80}")
    print(f"SCENARIO: {label}")
    print(f"  SS={cfg.strand_specificity}, n_rna={cfg.n_rna_regions}, "
          f"n_gdna={cfg.n_gdna_regions}, true_pi={cfg.true_pi:.3f}")
    print(f"  rna_count_mean={cfg.rna_count_mean}, gdna_count_mean={cfg.gdna_count_mean}")
    print(f"  rna_fl_mean={cfg.rna_fl_mean}, gdna_fl_mean={cfg.gdna_fl_mean}")
    print(f"{'='*80}")

    data = generate_scenario(cfg)
    truth = data["truth_labels"]
    rc = data["region_counts"]
    rdf = data["region_df"]
    fl = data["fl_table"]

    # Phase 1: statistics
    stats = compute_region_stats(rc, rdf)
    n_regions = len(stats["n_total"])
    eligible = (stats["n_total"] > 0) & (stats["region_length"] > 0)
    sense_frac = compute_sense_fraction(stats)
    log_d, epsilon = compute_log_density(stats, eligible)

    is_gdna = truth == 1
    is_rna = truth == 0
    has_splice = stats["n_spliced"] > 0

    print(f"\n--- Data summary ---")
    print(f"  Total regions: {n_regions}, eligible: {eligible.sum()}")
    print(f"  RNA regions with spliced reads: {(is_rna & has_splice).sum()}/{is_rna.sum()}")
    print(f"  gDNA regions with spliced reads: {(is_gdna & has_splice).sum()}/{is_gdna.sum()}")
    print(f"  RNA unspliced-only: {(is_rna & ~has_splice & eligible).sum()}")

    print(f"\n--- Density (n/L) ---")
    rna_densities = stats["n_total"][is_rna] / stats["region_length"][is_rna]
    gdna_densities = stats["n_total"][is_gdna] / stats["region_length"][is_gdna]
    print(f"  RNA  density: mean={rna_densities.mean():.4f}, "
          f"median={np.median(rna_densities):.4f}, "
          f"p10={np.percentile(rna_densities, 10):.4f}, "
          f"p90={np.percentile(rna_densities, 90):.4f}")
    print(f"  gDNA density: mean={gdna_densities.mean():.4f}, "
          f"median={np.median(gdna_densities):.4f}, "
          f"p10={np.percentile(gdna_densities, 10):.4f}, "
          f"p90={np.percentile(gdna_densities, 90):.4f}")

    rna_log_d = log_d[is_rna]
    gdna_log_d = log_d[is_gdna]
    print(f"  RNA  log_d:   mean={rna_log_d.mean():.3f}, "
          f"std={rna_log_d.std():.3f}")
    print(f"  gDNA log_d:   mean={gdna_log_d.mean():.3f}, "
          f"std={gdna_log_d.std():.3f}")
    overlap = np.sum((gdna_log_d[:, None] >= rna_log_d[None, :]).any(axis=1))
    print(f"  gDNA regions with log_d >= min(RNA log_d): {overlap}/{is_gdna.sum()}")

    print(f"\n--- Strand ---")
    rna_sf = sense_frac[is_rna & np.isfinite(sense_frac)]
    gdna_sf = sense_frac[is_gdna & np.isfinite(sense_frac)]
    print(f"  RNA  sense_frac: mean={rna_sf.mean():.3f}, std={rna_sf.std():.3f}")
    print(f"  gDNA sense_frac: mean={gdna_sf.mean():.3f}, std={gdna_sf.std():.3f}")

    # Phase 2: bin edges
    log_d_eligible = log_d[eligible]
    hist_min = float(np.min(log_d_eligible)) - 0.5
    hist_max = float(np.max(log_d_eligible)) + 0.5
    density_bin_edges = np.linspace(hist_min, hist_max, _DENSITY_BINS + 1)
    print(f"\n--- Density histogram ---")
    print(f"  Bin range: [{hist_min:.2f}, {hist_max:.2f}], n_bins={_DENSITY_BINS}")

    # Phase 3: seed
    gamma, pi_init, init_diag = _seed_initial_partition(
        stats, log_d, sense_frac, cfg.strand_specificity, eligible,
    )
    print(f"\n--- Seed initialization ---")
    print(f"  n_expressed_seed={init_diag['n_expressed_seed']}, "
          f"n_gdna_seed={init_diag['n_gdna_seed']}, pi_init={pi_init:.4f}")
    print(f"  gDNA gamma after seed: mean={gamma[is_gdna].mean():.3f}, "
          f"median={np.median(gamma[is_gdna]):.3f}")
    print(f"  RNA  gamma after seed: mean={gamma[is_rna].mean():.3f}")
    # Count gDNA regions in seed partition
    print(f"  gDNA regions with gamma=1 (seeded): {(gamma[is_gdna] == 1.0).sum()}")
    print(f"  gDNA regions with gamma=0 (spliced): {(gamma[is_gdna] == 0.0).sum()}")
    print(f"  gDNA regions with gamma=0.5 (unknown): "
          f"{((gamma[is_gdna] > 0) & (gamma[is_gdna] < 1)).sum()}")

    # Build initial histograms
    gdna_density_hist = build_density_histogram(
        log_d[eligible], gamma[eligible] * stats["n_total"][eligible],
        density_bin_edges,
    )
    rna_density_hist = build_density_histogram(
        log_d[eligible],
        (1.0 - gamma)[eligible] * stats["n_total"][eligible],
        density_bin_edges,
    )
    gdna_strand_hist = build_strand_histogram(sense_frac, gamma, symmetric=True)
    rna_strand_hist = build_strand_histogram(sense_frac, 1.0 - gamma, symmetric=False)

    # FL models
    fl_region_ids = fl["region_id"].values.astype(np.intp)
    fl_frag_lens = fl["frag_len"].values.astype(np.intp)
    gdna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, gamma)
    rna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, 1.0 - gamma)

    # Show initial density histogram mass
    print(f"\n--- Initial density histograms ---")
    gdna_d = gdna_density_hist[1]
    rna_d = rna_density_hist[1]
    # Find peak bins
    gdna_peak = np.argmax(gdna_d)
    rna_peak = np.argmax(rna_d)
    bin_centers = 0.5 * (density_bin_edges[:-1] + density_bin_edges[1:])
    print(f"  gDNA peak: bin {gdna_peak} (center={bin_centers[gdna_peak]:.2f}), "
          f"density={gdna_d[gdna_peak]:.4f}")
    print(f"  RNA  peak: bin {rna_peak} (center={bin_centers[rna_peak]:.2f}), "
          f"density={rna_d[rna_peak]:.4f}")

    # Compute LLR from initial histograms
    llr_density = _compute_density_llr(
        log_d, gdna_density_hist, rna_density_hist, eligible, n_regions,
    )
    llr_strand = _compute_strand_llr_binomial_mixture(
        stats, gdna_strand_hist, rna_strand_hist, n_regions,
    )

    print(f"\n--- Initial LLR (before first E-step) ---")
    soft = eligible & ~has_splice
    if soft.any():
        print(f"  Soft regions: {soft.sum()} "
              f"(gDNA: {(soft & is_gdna).sum()}, RNA: {(soft & is_rna).sum()})")
        print(f"  Density LLR (gDNA regions): "
              f"mean={llr_density[is_gdna & soft].mean():.3f}, "
              f"std={llr_density[is_gdna & soft].std():.3f}, "
              f"min={llr_density[is_gdna & soft].min():.3f}, "
              f"max={llr_density[is_gdna & soft].max():.3f}")
        if (is_rna & soft).any():
            print(f"  Density LLR (RNA soft regions): "
                  f"mean={llr_density[is_rna & soft].mean():.3f}")
        print(f"  Strand LLR  (gDNA regions): "
              f"mean={llr_strand[is_gdna & soft].mean():.3f}, "
              f"std={llr_strand[is_gdna & soft].std():.3f}")
        if (is_rna & soft).any():
            print(f"  Strand LLR  (RNA soft regions): "
                  f"mean={llr_strand[is_rna & soft].mean():.3f}")

        # Show what combined log-odds would be
        pi_safe = np.clip(pi_init, 1e-12, 1 - 1e-12)
        log_prior = np.log(pi_safe / (1 - pi_safe))
        log_odds_gdna = log_prior + llr_density[is_gdna & soft] + llr_strand[is_gdna & soft]
        gamma_gdna = 1.0 / (1.0 + np.exp(-np.clip(log_odds_gdna, -500, 500)))
        print(f"\n  Log prior odds: {log_prior:.3f}")
        print(f"  gDNA log-odds: mean={log_odds_gdna.mean():.3f}, "
              f"min={log_odds_gdna.min():.3f}, max={log_odds_gdna.max():.3f}")
        print(f"  gDNA gamma from first E-step: mean={gamma_gdna.mean():.4f}, "
              f"min={gamma_gdna.min():.4f}, max={gamma_gdna.max():.4f}")
        print(f"  gDNA gamma > 0.5: {(gamma_gdna > 0.5).sum()}/{len(gamma_gdna)}")

    # Now run the full EM with diagnostics and trace each iteration
    print(f"\n--- EM iterations ---")
    pi = pi_init
    for iteration in range(20):
        gamma_new = _e_step(
            stats, pi, log_d, sense_frac, eligible,
            gdna_density_hist=gdna_density_hist,
            rna_density_hist=rna_density_hist,
            gdna_strand_hist=gdna_strand_hist,
            rna_strand_hist=rna_strand_hist,
            fl_region_ids=fl_region_ids,
            fl_frag_lens=fl_frag_lens,
            gdna_fl_model=gdna_fl,
            rna_fl_model=rna_fl,
        )

        (new_pi, lG, lE,
         gdna_density_hist, rna_density_hist,
         gdna_strand_hist, rna_strand_hist) = _m_step(
            stats, gamma_new, log_d, sense_frac, eligible, density_bin_edges,
        )
        gdna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, gamma_new)
        rna_fl = build_gdna_fl_model(fl_region_ids, fl_frag_lens, 1.0 - gamma_new)

        delta = abs(new_pi - pi)
        gdna_gamma_vals = gamma_new[is_gdna & soft] if soft.any() else gamma_new[is_gdna]
        rna_gamma_vals = gamma_new[is_rna]

        # Recompute LLR for diagnostics
        llr_d = _compute_density_llr(
            log_d, gdna_density_hist, rna_density_hist, eligible, n_regions,
        )
        llr_s = _compute_strand_llr_binomial_mixture(
            stats, gdna_strand_hist, rna_strand_hist, n_regions,
        )

        print(f"  iter {iteration+1:2d}: π={new_pi:.4f} (Δ={delta:.6f}), "
              f"gDNA_γ={gdna_gamma_vals.mean():.4f} "
              f"[{gdna_gamma_vals.min():.4f},{gdna_gamma_vals.max():.4f}], "
              f"RNA_γ={rna_gamma_vals.mean():.4f}, "
              f"dens_llr_gdna={llr_d[is_gdna & soft].mean():.2f} "
              f"strand_llr_gdna={llr_s[is_gdna & soft].mean():.2f}" if soft.any() else "")

        pi = new_pi
        gamma = gamma_new
        if delta < 1e-4:
            print(f"  Converged at iteration {iteration+1}")
            break

    # Final evaluation
    import time
    result = calibrate_gdna(rc, fl, rdf, cfg.strand_specificity)
    eval_result = evaluate_scenario(result, truth, cfg, 0.0)
    print(f"\n--- Final result ---")
    print(f"  AUC={eval_result.auroc:.4f}, gF1={eval_result.gdna_f1:.4f}, "
          f"rF1={eval_result.rna_f1:.4f}")
    print(f"  π={result.mixing_proportion:.4f}, π_err={eval_result.pi_error:.4f}")
    print(f"  mean_gdna_post={eval_result.mean_gdna_posterior:.4f}, "
          f"mean_rna_post={eval_result.mean_rna_posterior:.4f}")
    print(f"  separation={eval_result.separation:.4f}")
    return eval_result


def main():
    print("DIAGNOSING UNSTRANDED (SS=0.50) POSTERIOR COLLAPSE")
    print("=" * 80)

    # Scenario 1: low gDNA, identical FL — should be easy for density
    trace_scenario("low_gdna_fl_identical", ScenarioConfig(
        n_rna_regions=800, n_gdna_regions=200,
        strand_specificity=0.50,
        gdna_count_mean=30.0,
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        seed=42,
        label="diag_low_gdna",
    ))

    # Scenario 2: med gDNA, identical FL
    trace_scenario("med_gdna_fl_identical", ScenarioConfig(
        n_rna_regions=600, n_gdna_regions=400,
        strand_specificity=0.50,
        gdna_count_mean=60.0,
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        seed=42,
        label="diag_med_gdna",
    ))

    # Scenario 3: equal gDNA, identical FL — hardest case
    trace_scenario("equal_gdna_fl_identical", ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.50,
        gdna_count_mean=200.0,  # same as RNA
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        seed=42,
        label="diag_equal_gdna",
    ))

    # Scenario 4: equal gDNA but LOWER density — realistic additive model
    # In reality, expressed = gDNA + RNA, so gDNA-only should have LOWER density
    trace_scenario("equal_gdna_lower_density", ScenarioConfig(
        n_rna_regions=500, n_gdna_regions=500,
        strand_specificity=0.50,
        gdna_count_mean=50.0,  # gDNA contributes far less than RNA
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        seed=42,
        label="diag_equal_lower",
    ))

    # Scenario 5: high strand specificity for comparison
    trace_scenario("low_gdna_ss0.95", ScenarioConfig(
        n_rna_regions=800, n_gdna_regions=200,
        strand_specificity=0.95,
        gdna_count_mean=30.0,
        gdna_fl_mean=250.0, gdna_fl_std=50.0,
        seed=42,
        label="diag_ss0.95",
    ))


if __name__ == "__main__":
    main()
