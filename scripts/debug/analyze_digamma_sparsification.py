#!/usr/bin/env python3
"""
VBEM Digamma Sparsification Analysis
=====================================

Key finding from Experiment 1: warm-start at γ works for large N but FAILS
for small N with tiny C, even though the initialization is correct.

Root cause: exp(ψ(n + α)) ≈ n + α - 0.5 for large n. When α < 0.5, the
VBEM systematically shrinks the component by ~(0.5 - α) per iteration.
Over hundreds of iterations, this drift drains small components to zero.

The threshold is α = 0.5 (the Jeffreys prior), below which VBEM sparsifies.

This script:
1. Quantifies the per-iteration drift as a function of α
2. Tests the "α_gDNA floor" solution: α_gDNA = max(0.5, γ*C)
3. Tests across the full range of N, γ, and n_components
"""
import numpy as np
from scipy.special import digamma


def vbem(N, alphas, theta_init, n_iter=500):
    """VBEM with degenerate likelihoods (all components identical)."""
    theta = theta_init.copy().astype(np.float64)
    for _ in range(n_iter):
        psi = digamma(theta + alphas)
        psi -= psi.max()
        r = np.exp(psi)
        r /= r.sum()
        theta_new = N * r
        if np.allclose(theta_new, theta, atol=1e-12):
            break
        theta = theta_new
    return theta


def vbem_partial_ll(N, alphas, theta_init, log_liks, n_iter=500):
    """VBEM with per-component log-likelihood offsets."""
    theta = theta_init.copy().astype(np.float64)
    for _ in range(n_iter):
        psi = digamma(theta + alphas)
        psi -= psi.max()
        r = np.exp(psi + log_liks)
        r /= r.sum()
        theta_new = N * r
        if np.allclose(theta_new, theta, atol=1e-12):
            break
        theta = theta_new
    return theta


# ==================================================================
print("=" * 70)
print("Part 1: Per-iteration drift as a function of α")
print("=" * 70)
print()
print("For a 2-component system with θ_1=200, θ_2=1800 (N=2000),")
print("compare exp(ψ(θ+α))/Σ vs θ/N to measure sparsification drift.")
print()
print(f"{'α':>8} {'exp(ψ(θ+α))':>14} {'θ ratio':>10} {'drift/iter':>12} "
      f"{'iters to collapse':>18}")
for alpha_g in [0.01, 0.05, 0.10, 0.20, 0.30, 0.50, 0.75, 1.00, 2.00]:
    theta_g = 200.0
    theta_r = 1800.0
    alpha_r = 10.0 - alpha_g  # keep total C=10 for comparison
    
    psi_g = digamma(theta_g + alpha_g)
    psi_r = digamma(theta_r + alpha_r)
    exp_g = np.exp(psi_g - max(psi_g, psi_r))
    exp_r = np.exp(psi_r - max(psi_g, psi_r))
    ratio = exp_g / (exp_g + exp_r)
    theta_new = 2000 * ratio
    drift = theta_new - theta_g
    iters = abs(theta_g / drift) if abs(drift) > 1e-10 else float('inf')
    print(f"{alpha_g:>8.2f} {theta_new:>14.4f} {theta_g/2000:>10.4f} "
          f"{drift:>+12.4f} {iters:>18.0f}")

print()
print("KEY: drift > 0 means component grows, < 0 means sparsification.")
print("     Neutral point is at α ≈ 0.5 (Jeffreys prior).")


# ==================================================================
print()
print("=" * 70)
print("Part 2: α_gDNA floor solution")
print("=" * 70)
print()
print("Design: α_gDNA = max(α_floor, γ*C), α_RNA = (1-γ)*C")
print("Combined with warm-start θ_init[gDNA] = γ*N")
print()

gamma = 0.10
alpha_floors = [0.0, 0.25, 0.50, 0.75, 1.0]
C_base = 2.0  # small base C

for N in [100, 500, 2000, 10000, 100000, 1000000]:
    print(f"N={N:>10,}, γ={gamma}, C_base={C_base}")
    theta_init = np.array([gamma * N, (1 - gamma) * N])
    
    for floor in alpha_floors:
        alpha_g = max(floor, gamma * C_base)
        alpha_r = (1 - gamma) * C_base
        alphas = np.array([alpha_g, alpha_r])
        theta = vbem(N, alphas, theta_init)
        frac = theta[0] / theta.sum()
        error = abs(frac - gamma) / gamma
        print(f"  floor={floor:.2f}: α_gDNA={alpha_g:.2f}, α_RNA={alpha_r:.2f}, "
              f"θ_gDNA={frac:.6f} (err={error:.1%})")
    print()


# ==================================================================
print()
print("=" * 70)
print("Part 3: Floor + warm-start across γ values")
print("=" * 70)
print()
print("Does the floor work for low γ (heavy RNA) and high γ (heavy gDNA)?")
print()

N = 2000
C_base = 2.0
alpha_floor = 0.50

for gamma in [0.001, 0.01, 0.05, 0.10, 0.20, 0.50, 0.80]:
    theta_init = np.array([gamma * N, (1 - gamma) * N])
    
    # With floor
    alpha_g = max(alpha_floor, gamma * C_base)
    alpha_r = max(alpha_floor, (1 - gamma) * C_base)  # floor for RNA too
    alphas = np.array([alpha_g, alpha_r])
    theta = vbem(N, alphas, theta_init)
    frac = theta[0] / theta.sum()
    error = abs(frac - gamma) / max(gamma, 1e-10)
    
    # Without floor (standard C)
    alpha_g_nf = gamma * C_base
    alpha_r_nf = (1 - gamma) * C_base
    alphas_nf = np.array([alpha_g_nf, alpha_r_nf])
    theta_nf = vbem(N, alphas_nf, theta_init)
    frac_nf = theta_nf[0] / theta_nf.sum()
    error_nf = abs(frac_nf - gamma) / max(gamma, 1e-10)
    
    print(f"γ={gamma:.3f}: WITH floor θ_gDNA={frac:.6f} (err={error:.1%}), "
          f"WITHOUT floor θ_gDNA={frac_nf:.6f} (err={error_nf:.1%})")


# ==================================================================
print()
print("=" * 70)
print("Part 4: Multi-component test (gDNA + 4 mRNA + 4 nRNA)")
print("=" * 70)
print()
print("Test that the floor approach prevents gDNA collapse")
print("while keeping RNA pseudocounts small (no siphoning)")
print()

gamma = 0.10
n_mrna = 4
n_nrna = 4
n_c = 1 + n_mrna + n_nrna
C_base = 2.0
alpha_floor = 0.50

for N in [100, 2000, 100000]:
    print(f"N={N:>10,}, γ={gamma}, C_base={C_base}, floor={alpha_floor}")
    
    # Warm-start: gDNA at γ*N, mRNA gets 65%, nRNA gets 35% of RNA
    theta_rna_total = (1 - gamma) * N
    theta_init = np.zeros(n_c)
    theta_init[0] = gamma * N
    theta_init[1:1+n_mrna] = 0.65 * theta_rna_total / n_mrna
    theta_init[1+n_mrna:] = 0.35 * theta_rna_total / n_nrna
    
    # Prior with floor on gDNA only
    alpha_g = max(alpha_floor, gamma * C_base)
    alpha_per_rna = (1 - gamma) * C_base / (n_mrna + n_nrna)
    alphas = np.zeros(n_c)
    alphas[0] = alpha_g
    alphas[1:] = alpha_per_rna
    
    theta = vbem(N, alphas, theta_init)
    total = theta.sum()
    gdna_frac = theta[0] / total
    mrna_frac = theta[1:1+n_mrna].sum() / total
    nrna_frac = theta[1+n_mrna:].sum() / total
    
    # Compare with no floor
    alphas_nf = np.zeros(n_c)
    alphas_nf[0] = gamma * C_base
    alphas_nf[1:] = alpha_per_rna
    theta_nf = vbem(N, alphas_nf, theta_init)
    total_nf = theta_nf.sum()
    
    print(f"  WITH floor: gDNA={gdna_frac:.4f} mRNA={mrna_frac:.4f} "
          f"nRNA={nrna_frac:.4f} (C_eff={alphas.sum():.2f})")
    print(f"  NO floor:   gDNA={theta_nf[0]/total_nf:.4f} "
          f"mRNA={theta_nf[1:1+n_mrna].sum()/total_nf:.4f} "
          f"nRNA={theta_nf[1+n_mrna:].sum()/total_nf:.4f} "
          f"(C_eff={alphas_nf.sum():.2f})")
    
    # Also compare with old κ×N approach
    kappa = 0.50
    C_old = kappa * N
    alphas_old = np.zeros(n_c)
    alphas_old[0] = gamma * C_old
    alpha_per_rna_old = (1 - gamma) * C_old / (n_mrna + n_nrna)
    alphas_old[1:] = alpha_per_rna_old
    theta_old = vbem(N, alphas_old, theta_init)
    total_old = theta_old.sum()
    
    print(f"  κ×N (old):  gDNA={theta_old[0]/total_old:.4f} "
          f"mRNA={theta_old[1:1+n_mrna].sum()/total_old:.4f} "
          f"nRNA={theta_old[1+n_mrna:].sum()/total_old:.4f} "
          f"(C_eff={alphas_old.sum():.1f})")
    print()


# ==================================================================
print()
print("=" * 70)
print("Part 5: Partially degenerate with floor (strand signal present)")
print("=" * 70)
print()
print("When strand signal gives gDNA a per-fragment LLR penalty,")
print("does the floor interfere? (It shouldn't — likelihood dominates.)")
print()

gamma = 0.10
N = 2000
alpha_floor = 0.50
C_base = 2.0

for SS in [0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90]:
    if SS == 0.50:
        avg_llr = 0.0
    else:
        llr_sense = np.log(0.5 / SS)
        llr_anti = np.log(0.5 / (1 - SS))
        avg_llr = SS * llr_sense + (1 - SS) * llr_anti
    
    # With floor
    alpha_g = max(alpha_floor, gamma * C_base)
    alpha_r = (1 - gamma) * C_base
    alphas = np.array([alpha_g, alpha_r])
    theta_init = np.array([gamma * N, (1 - gamma) * N])
    log_liks = np.array([avg_llr * N, 0.0])  # total LLR for gDNA component
    
    # For partially degenerate: use per-fragment LLR
    theta = vbem_partial_ll(N, alphas, theta_init, np.array([avg_llr, 0.0]))
    frac = theta[0] / theta.sum()
    
    # Pure RNA (γ=0, no gDNA in data)
    # With floor, but no gDNA in the data: likelihood kills gDNA?
    alpha_g_pure = max(alpha_floor, 0.0 * C_base)  # γ=0 → floor applies
    alphas_pure = np.array([alpha_g_pure, C_base])
    theta_init_pure = np.array([0.0 * N + 1.0, (1.0) * N - 1.0])  # minimal gDNA init
    theta_pure = vbem_partial_ll(N, alphas_pure, theta_init_pure,
                                  np.array([avg_llr, 0.0]))
    frac_pure = theta_pure[0] / theta_pure.sum()
    
    print(f"SS={SS:.2f} (LLR={avg_llr:+.4f}/frag): "
          f"θ_gDNA={frac:.6f} (γ=0.10 in data), "
          f"θ_gDNA_pure={frac_pure:.6f} (γ=0 in data)")


# ==================================================================
print()
print("=" * 70)
print("Part 6: The γ=0 case — does the floor create false gDNA?")
print("=" * 70)
print()
print("Critical test: when calibration says γ=0, the gDNA component")
print("should be disabled. With α_floor=0.5, does gDNA leak?")
print()

gamma_cal = 0.0
N = 2000

# Case A: floor applies, but we check if gDNA gate (α_gDNA=0 → disabled)
# should prevent this
print("Case A: Calibration γ=0, gDNA gate disables component (current code)")
print("  → gDNA prior set to 0.0, component disabled. Floor irrelevant.")
print()

# Case B: Calibration γ=0.001 (near-zero but not exactly 0)
gamma_cal = 0.001
print(f"Case B: Calibration γ={gamma_cal} (near-zero)")
alpha_g = max(0.50, gamma_cal * C_base)
alpha_r = (1 - gamma_cal) * C_base
print(f"  α_gDNA = max(0.5, {gamma_cal*C_base:.4f}) = {alpha_g:.2f}")
print(f"  α_RNA = {alpha_r:.4f}")
print(f"  Implied prior gDNA fraction: {alpha_g/(alpha_g+alpha_r):.4f} "
      f"(vs calibrated γ={gamma_cal:.4f})")
print()
print("  This biases gDNA upward from 0.1% to "
      f"{alpha_g/(alpha_g+alpha_r)*100:.1f}%. Is this acceptable?")
print()

# What fraction does VBEM converge to?
for SS in [0.50, 0.65, 0.90]:
    if SS == 0.50:
        avg_llr = 0.0
    else:
        llr_sense = np.log(0.5 / SS)
        llr_anti = np.log(0.5 / (1 - SS))
        avg_llr = SS * llr_sense + (1 - SS) * llr_anti
    
    alphas = np.array([alpha_g, alpha_r])
    theta_init = np.array([gamma_cal * N, (1 - gamma_cal) * N])
    theta = vbem_partial_ll(N, alphas, theta_init, np.array([avg_llr, 0.0]))
    frac = theta[0] / theta.sum()
    print(f"  SS={SS:.2f}: θ_gDNA = {frac:.6f} "
          f"(true γ={gamma_cal}, floor-amplified prior)")

print()
print()

# Case C: threshold approach — only apply floor when γ > some minimum
gamma_threshold = 0.01
print(f"Case C: Alternative — apply floor only when γ > {gamma_threshold}")
print("  When γ ≤ threshold: α_gDNA = γ*C (may sparsify, but γ is tiny)")
print("  When γ > threshold: α_gDNA = max(0.5, γ*C)")
print()

for gamma_test in [0.001, 0.005, 0.01, 0.05, 0.10, 0.30]:
    if gamma_test > gamma_threshold:
        alpha_g = max(0.50, gamma_test * C_base)
    else:
        alpha_g = gamma_test * C_base
    alpha_r = (1 - gamma_test) * C_base
    
    alphas = np.array([alpha_g, alpha_r])
    theta_init = np.array([gamma_test * N, (1 - gamma_test) * N])
    theta = vbem(N, alphas, theta_init)
    frac = theta[0] / theta.sum()
    print(f"  γ={gamma_test:.3f}: α_gDNA={alpha_g:.4f}, θ_gDNA={frac:.6f}")


# ==================================================================
print()
print("=" * 70)
print("SYNTHESIS: Proposed design")
print("=" * 70)
print("""
The VBEM has an inherent sparsification bias from the digamma function:
  exp(ψ(n + α)) ≈ n + α - 0.5
When α < 0.5, the effective VBEM weight is LESS than the count, causing
systematic downward drift of ~(0.5 - α) per iteration.

This means:
- Components with α < 0.5 get driven toward zero
- The drift is ABSOLUTE, not relative to N
- For small N: drift dominates → false sparsification
- For large N: drift negligible → warm-start works

PROPOSED DESIGN:
1. Global strand calibration (fix rectifier bias)
2. Warm-start: θ_gDNA_init = γ × N_locus
3. Prior: α_gDNA = max(0.5, γ * C_base) with small C_base (e.g., 2.0)
4. Prior: α_RNA_j proportional to coverage, summing to (1-γ)*C_base
5. Gate: if γ = 0.0, disable gDNA component entirely (existing behavior)

Properties:
- Not brittle: single parameter C_base ∈ [1, 5], insensitive to value
- Works for any N: warm-start for large N, floor for small N
- No siphoning: total RNA α ≈ 1.8, spread across n_rna components
- gDNA protected: α_gDNA = 0.5 prevents false sparsification
- Compatible with gDNA gate: when γ=0, gDNA disabled
""")
