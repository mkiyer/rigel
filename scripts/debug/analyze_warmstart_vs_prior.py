#!/usr/bin/env python3
"""Analyze the interaction between warm-start and prior in VBEM.

Key question: Is the warm-start (initialization) or the prior (C)
the operative mechanism for anchoring gDNA in the degenerate case?

If warm-start is sufficient, we can use a small fixed C and avoid
the entire κ×N scaling problem.
"""
import numpy as np
from scipy.special import digamma


def vbem(N, alphas, theta_init, n_iter=500, verbose=False):
    """Run VBEM with given initialization and priors.
    
    All components have identical likelihoods (fully degenerate).
    """
    n_c = len(alphas)
    theta = theta_init.copy().astype(np.float64)
    
    for it in range(n_iter):
        # E-step: responsibilities ∝ exp(digamma(theta + alpha))
        # (likelihoods are all 1.0)
        psi = digamma(theta + alphas)
        psi -= psi.max()  # numerical stability
        r = np.exp(psi)
        r /= r.sum()
        # M-step
        theta_new = N * r
        if np.allclose(theta_new, theta, atol=1e-10):
            if verbose:
                print(f"  Converged at iteration {it}")
            break
        theta = theta_new
    
    return theta / theta.sum()


def warm_start_gdna(N, gamma, n_rna_components):
    """Reproduce rigel's warm-start: θ_gDNA = (γ/(1-γ)) × Σ(θ_RNA)."""
    # Initially, RNA components split N equally
    theta_rna_total = N  # all starts as RNA
    theta_gdna = (gamma / max(1.0 - gamma, 1e-12)) * theta_rna_total
    theta_rna_each = theta_rna_total / n_rna_components
    theta = np.zeros(n_rna_components + 1)
    theta[0] = theta_gdna  # gDNA
    theta[1:] = theta_rna_each
    return theta


def coverage_start(N, n_components):
    """Start from coverage-proportional split (≈uniform for degenerate)."""
    return np.full(n_components, N / n_components)


print("=" * 70)
print("Experiment 1: Warm-start at γ with tiny C vs large C")
print("=" * 70)
print()
print("Degenerate likelihoods, γ=0.10, 2 components (gDNA + RNA)")
print()

gamma = 0.10
for N in [100, 2000, 100000, 1000000]:
    print(f"N = {N:>10,}")
    
    # Warm-start from γ
    theta_init = np.array([gamma * N, (1 - gamma) * N])
    
    for C in [0.5, 1.0, 2.0, 5.0, 10.0, 100.0]:
        alphas = np.array([gamma * C, (1 - gamma) * C])
        result = vbem(N, alphas, theta_init)
        print(f"  C={C:>6.1f}: θ_gDNA={result[0]:.6f} "
              f"(target={gamma:.2f}, init=γ)")
    
    # Also test coverage-proportional start (uniform)
    theta_init_unif = np.array([N / 2.0, N / 2.0])
    for C in [1.0, 10.0, 100.0]:
        alphas = np.array([gamma * C, (1 - gamma) * C])
        result = vbem(N, alphas, theta_init_unif)
        print(f"  C={C:>6.1f}: θ_gDNA={result[0]:.6f} "
              f"(target={gamma:.2f}, init=UNIFORM)")
    print()


print("=" * 70)
print("Experiment 2: Multi-component (gDNA + 4 mRNA + 4 nRNA = 9 components)")
print("=" * 70)
print()
print("Degenerate likelihoods, γ=0.10, 9 components")
print("nRNA gets 35% of RNA coverage (mimicking SS=0.65)")
print()

gamma = 0.10
n_mrna = 4
n_nrna = 4
n_c = 1 + n_mrna + n_nrna

for N in [2000, 100000]:
    print(f"N = {N:>10,}")
    
    # Warm-start from γ, with RNA split: 65% to mRNA, 35% to nRNA
    theta_rna_total = (1 - gamma) * N
    theta_mrna_each = 0.65 * theta_rna_total / n_mrna
    theta_nrna_each = 0.35 * theta_rna_total / n_nrna
    theta_init = np.zeros(n_c)
    theta_init[0] = gamma * N  # gDNA
    theta_init[1:1+n_mrna] = theta_mrna_each  # mRNA
    theta_init[1+n_mrna:] = theta_nrna_each  # nRNA
    
    for C in [1.0, 2.0, 5.0, 10.0, 100.0]:
        # Prior: gDNA gets γ*C, RNA gets (1-γ)*C split by coverage
        alphas = np.zeros(n_c)
        alphas[0] = gamma * C
        # RNA prior proportional to coverage
        rna_coverage = theta_init[1:].copy()
        rna_coverage /= rna_coverage.sum()
        alphas[1:] = (1 - gamma) * C * rna_coverage
        
        result = vbem(N, alphas, theta_init)
        mrna_total = result[1:1+n_mrna].sum()
        nrna_total = result[1+n_mrna:].sum()
        print(f"  C={C:>6.1f}: gDNA={result[0]:.4f} mRNA={mrna_total:.4f} "
              f"nRNA={nrna_total:.4f}")
    print()


print("=" * 70)
print("Experiment 3: Partially degenerate likelihoods")
print("=" * 70)
print()
print("gDNA has per-fragment LLR relative to RNA")
print("γ=0.10, N=2000, 2 components. Warm-start from γ.")
print()

def vbem_partial(N, alphas, theta_init, llr_gdna, n_iter=500):
    """VBEM with gDNA having different likelihood from RNA."""
    n_c = len(alphas)
    theta = theta_init.copy().astype(np.float64)
    lik = np.ones(n_c)
    lik[0] = np.exp(llr_gdna)  # gDNA relative likelihood
    
    for it in range(n_iter):
        psi = digamma(theta + alphas)
        psi -= psi.max()
        r = lik * np.exp(psi)
        r /= r.sum()
        theta_new = N * r
        if np.allclose(theta_new, theta, atol=1e-10):
            break
        theta = theta_new
    return theta / theta.sum()


gamma = 0.10
N = 2000
theta_init = np.array([gamma * N, (1 - gamma) * N])

for llr in [0.0, -0.005, -0.01, -0.05, -0.1, -0.5]:
    print(f"LLR={llr:+.3f}:")
    for C in [1.0, 2.0, 5.0, 10.0, 100.0]:
        alphas = np.array([gamma * C, (1 - gamma) * C])
        result = vbem_partial(N, alphas, theta_init, llr)
        print(f"  C={C:>6.1f}: θ_gDNA={result[0]:.6f}")
    print()


print("=" * 70)
print("Experiment 4: Average strand LLR at various SS values")
print("=" * 70)
print()

for SS in [0.50, 0.52, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 0.95, 1.00]:
    if SS == 1.0:
        avg_llr = float('-inf')
        print(f"SS={SS:.2f}: avg strand LLR = -inf (perfect separation)")
        continue
    if SS == 0.50:
        avg_llr = 0.0
        print(f"SS={SS:.2f}: avg strand LLR = 0.000 (degenerate)")
        continue
    # Average LLR for gDNA vs RNA (strand signal only)
    # Sense reads: log(0.5/SS), fraction SS of RNA
    # Antisense reads: log(0.5/(1-SS)), fraction (1-SS) of RNA
    llr_sense = np.log(0.5 / SS)
    llr_anti = np.log(0.5 / (1 - SS))
    avg_llr = SS * llr_sense + (1 - SS) * llr_anti
    # Total for N=2000
    total_llr = avg_llr * 2000
    print(f"SS={SS:.2f}: avg strand LLR = {avg_llr:.4f}/frag, "
          f"total for 2000 frags = {total_llr:.1f} nats")
