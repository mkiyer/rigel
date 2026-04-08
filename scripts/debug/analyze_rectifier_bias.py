#!/usr/bin/env python3
"""Quantify rectifier bias in the strand formula.

Computes the per-region vs global strand estimates for pure RNA
at various SS levels using the actual diagnostic scenario data.
Also models finite-sample noise to understand how often the global
estimator produces exact zero vs nonzero.
"""
import numpy as np
from scipy import stats as sp_stats

np.set_printoptions(precision=4)


def per_region_estimate(n_anti_arr, n_unspliced_arr, SS):
    """Current code: per-region max(0,...), then sum."""
    denom = max(SS - 0.5, 1e-12)
    e_per = np.maximum(0.0, (n_anti_arr - n_unspliced_arr * (1 - SS)) / denom)
    return e_per.sum()


def global_estimate(n_anti_arr, n_unspliced_arr, SS):
    """Proposed: sum first, then max(0,...)."""
    denom = max(SS - 0.5, 1e-12)
    global_anti = n_anti_arr.sum()
    global_unspliced = n_unspliced_arr.sum()
    return max(0.0, (global_anti - global_unspliced * (1 - SS)) / denom)


def simulate_pure_rna(n_regions, frags_per_region, SS, n_trials=10000):
    """Simulate pure RNA and compare per-region vs global estimators."""
    per_region_vals = []
    global_vals = []

    for _ in range(n_trials):
        # Each region has frags_per_region unspliced fragments
        # Antisense fraction = (1 - SS) for pure RNA
        n_unspliced = np.full(n_regions, frags_per_region)
        n_anti = np.random.binomial(n_unspliced, 1 - SS)

        per_region_vals.append(per_region_estimate(n_anti, n_unspliced, SS))
        global_vals.append(global_estimate(n_anti, n_unspliced, SS))

    per_arr = np.array(per_region_vals)
    glob_arr = np.array(global_vals)

    return {
        "per_region_mean": per_arr.mean(),
        "per_region_zero_frac": (per_arr == 0).mean(),
        "global_mean": glob_arr.mean(),
        "global_zero_frac": (glob_arr == 0).mean(),
    }


print("=" * 70)
print("Part 1: Rectifier bias quantification (Monte Carlo)")
print("=" * 70)
print()

# Scenario: 2 regions with ~1000 fragments each (like the failing test)
for SS in [0.50, 0.65, 0.90, 0.95, 1.0]:
    result = simulate_pure_rna(n_regions=2, frags_per_region=1000, SS=SS, n_trials=50000)
    print(f"SS={SS:.2f} | 2 regions × 1000 frags:")
    print(f"  Per-region: mean E[gDNA]={result['per_region_mean']:.2f}, "
          f"P(exactly 0)={result['per_region_zero_frac']:.3f}")
    print(f"  Global:     mean E[gDNA]={result['global_mean']:.2f}, "
          f"P(exactly 0)={result['global_zero_frac']:.3f}")
    print()

# What about real data scale? 500 regions
print("-" * 70)
print("Real data scale: 500 regions × 1000 frags each")
print("-" * 70)
for SS in [0.65, 0.90, 0.95]:
    result = simulate_pure_rna(n_regions=500, frags_per_region=1000, SS=SS, n_trials=5000)
    print(f"SS={SS:.2f} | 500 regions × 1000 frags:")
    print(f"  Per-region: mean E[gDNA]={result['per_region_mean']:.1f}, "
          f"P(exactly 0)={result['per_region_zero_frac']:.3f}")
    print(f"  Global:     mean E[gDNA]={result['global_mean']:.1f}, "
          f"P(exactly 0)={result['global_zero_frac']:.3f}")
    print()


print("=" * 70)
print("Part 2: SS measurement noise effect")
print("=" * 70)
print()

# The measured SS differs from true SS due to sampling noise.
# If SS_measured > SS_true, the formula expects fewer antisense reads
# than observed, producing false gDNA signal even with global estimator.
print("Effect of SS measurement error on global estimator:")
print("True SS=0.90, 2000 total unspliced fragments")
print()

true_SS = 0.90
n_total = 2000
n_trials = 50000

for ss_error in [-0.03, -0.01, 0.0, 0.01, 0.03]:
    measured_SS = true_SS + ss_error
    vals = []
    for _ in range(n_trials):
        n_anti = np.random.binomial(n_total, 1 - true_SS)
        denom = max(measured_SS - 0.5, 1e-12)
        e_gdna = max(0.0, (n_anti - n_total * (1 - measured_SS)) / denom)
        vals.append(e_gdna)
    vals = np.array(vals)
    print(f"  SS_measured={measured_SS:.3f} (error={ss_error:+.3f}): "
          f"mean E[gDNA]={vals.mean():.1f}, "
          f"P(zero)={float((vals == 0).mean()):.3f}")


print()
print("=" * 70)
print("Part 3: Does C magnitude matter for degenerate EM?")
print("=" * 70)
print()

# Simulate simple 2-component EM with degenerate likelihoods
# to verify that final allocation = prior ratio regardless of C.
def em_degenerate(n_frags, alpha_1, alpha_2, n_iter=200):
    """Run MAP-EM with identical likelihoods for all components."""
    # All fragments have likelihood 1.0 for both components
    theta = np.array([0.5, 0.5])  # initial
    alpha = np.array([alpha_1, alpha_2])
    for _ in range(n_iter):
        # E-step: degenerate → responsibility = theta / sum(theta) = theta
        resp = theta.copy()
        # M-step: MAP
        numerator = n_frags * resp + alpha - 1.0
        numerator = np.maximum(numerator, 1e-12)
        theta = numerator / numerator.sum()
    return theta


print("2-component EM with degenerate likelihoods, γ=0.10:")
for C in [0.5, 1.0, 2.0, 5.0, 10.0, 100.0, 1000.0]:
    gamma = 0.10
    alpha_g, alpha_r = gamma * C, (1 - gamma) * C
    theta = em_degenerate(2000, alpha_g, alpha_r)
    print(f"  C={C:>7.1f}: α_gDNA={alpha_g:.2f}, α_RNA={alpha_r:.2f} → "
          f"θ_gDNA={theta[0]:.6f}")

print()
print("But what about partially degenerate? (gDNA LLR = -0.1 per fragment)")
print()


def em_partial_degenerate(n_frags, alpha_1, alpha_2, llr_per_frag=-0.1, n_iter=200):
    """MAP-EM where component 2 has slightly higher likelihood."""
    theta = np.array([0.5, 0.5])
    alpha = np.array([alpha_1, alpha_2])
    # log-likelihood ratio: component 1 (gDNA) is llr_per_frag worse than component 2 (RNA)
    lik_ratio = np.exp(llr_per_frag)
    for _ in range(n_iter):
        # E-step: p(c=1|x) ∝ theta[0] * lik_ratio, p(c=2|x) ∝ theta[1] * 1.0
        resp = np.array([theta[0] * lik_ratio, theta[1] * 1.0])
        resp /= resp.sum()
        # M-step
        numerator = n_frags * resp + alpha - 1.0
        numerator = np.maximum(numerator, 1e-12)
        theta = numerator / numerator.sum()
    return theta


for llr in [-0.01, -0.05, -0.1, -0.5]:
    print(f"  LLR={llr} per fragment (gDNA disfavored):")
    for C in [1.0, 2.0, 5.0, 100.0, 1000.0]:
        gamma = 0.10
        alpha_g, alpha_r = gamma * C, (1 - gamma) * C
        theta = em_partial_degenerate(2000, alpha_g, alpha_r, llr)
        print(f"    C={C:>7.1f}: θ_gDNA={theta[0]:.6f}")
    print()
