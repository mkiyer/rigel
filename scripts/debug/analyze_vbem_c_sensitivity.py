#!/usr/bin/env python3
"""Analyze C magnitude effect on VBEM (what rigel actually uses)."""
import numpy as np
from scipy.special import digamma


def vbem_degenerate(n_frags, alpha_1, alpha_2, n_iter=500):
    """Run VBEM with identical likelihoods for all components."""
    alpha = np.array([alpha_1, alpha_2])
    theta = np.array([0.5, 0.5]) * n_frags  # start equal
    for _ in range(n_iter):
        # VBEM E-step: use digamma of theta + alpha
        q = np.exp(digamma(theta + alpha))  # likelihoods are 1.0 for both
        q /= q.sum()
        # VBEM M-step: each component gets share of N
        theta = n_frags * q
    # Effective allocation
    return theta / theta.sum()


def vbem_partial(n_frags, alpha_1, alpha_2, llr_per_frag, n_iter=500):
    """VBEM with partially degenerate likelihoods."""
    alpha = np.array([alpha_1, alpha_2])
    lik_ratio = np.array([np.exp(llr_per_frag), 1.0])
    theta = np.array([0.5, 0.5]) * n_frags
    for _ in range(n_iter):
        q = lik_ratio * np.exp(digamma(theta + alpha))
        q /= q.sum()
        theta = n_frags * q
    return theta / theta.sum()


print("VBEM with degenerate likelihoods, γ=0.10, N=2000:")
print(f"  {'C':>8} {'α_gDNA':>8} {'α_RNA':>8} {'θ_gDNA':>10} {'gDNA frags':>12}")
for C in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0]:
    gamma = 0.10
    theta = vbem_degenerate(2000, gamma * C, (1 - gamma) * C)
    print(f"  {C:>8.1f} {gamma*C:>8.2f} {(1-gamma)*C:>8.2f} {theta[0]:>10.6f} {theta[0]*2000:>12.1f}")

print()
print("VBEM partially degenerate, γ=0.10, N=2000:")
for llr in [-0.01, -0.05, -0.1]:
    print(f"\n  LLR={llr} (gDNA slightly disfavored by strand signal):")
    print(f"  {'C':>8} {'θ_gDNA':>10} {'gDNA frags':>12} {'vs γ=10%':>10}")
    for C in [1.0, 2.0, 5.0, 10.0, 50.0, 100.0, 500.0]:
        gamma = 0.10
        theta = vbem_partial(2000, gamma * C, (1 - gamma) * C, llr)
        expected = 0.10 * 2000
        actual = theta[0] * 2000
        print(f"  {C:>8.1f} {theta[0]:>10.6f} {actual:>12.1f} {actual-expected:>+10.1f}")

print()
print("VBEM degenerate, γ=0.05, N=2000 (low gDNA):")
print(f"  {'C':>8} {'α_gDNA':>8} {'α_RNA':>8} {'θ_gDNA':>10} {'gDNA frags':>12}")
for C in [1.0, 2.0, 5.0, 100.0, 1000.0]:
    gamma = 0.05
    theta = vbem_degenerate(2000, gamma * C, (1 - gamma) * C)
    print(f"  {C:>8.1f} {gamma*C:>8.2f} {(1-gamma)*C:>8.2f} {theta[0]:>10.6f} {theta[0]*2000:>12.1f}")

print()
print("VBEM degenerate with ZERO γ, N=2000 (pure RNA):")
print(f"  {'C':>8} {'α_gDNA':>8} {'α_RNA':>8} {'θ_gDNA':>10} {'gDNA frags':>12}")
for C in [1.0, 2.0, 5.0, 100.0, 1000.0]:
    gamma = 0.0
    alpha_g = max(gamma * C, 1e-10)  # tiny epsilon to avoid log(0)
    alpha_r = max((1 - gamma) * C, 1e-10)
    theta = vbem_degenerate(2000, alpha_g, alpha_r)
    print(f"  {C:>8.1f} {alpha_g:>8.4f} {alpha_r:>8.2f} {theta[0]:>10.6f} {theta[0]*2000:>12.1f}")
