"""Theory checks for the rejection-criterion question (no data needed).

1. Minimum ATTAINABLE two-sided p-value under BetaBinom(n, a, a) -> the seed size a screen needs.
2. The Fisher-information weight for the intraclass correlation, vs the shipped pair-count weight.
3. Bias of the pooled MoM under contamination (od_s -> 1 for a pure-RNA seed).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import betabinom


def pmin_two_sided(n, a):
    """2*P(K=0) under BetaBinom(n,a,a) — the smallest attainable two-sided tail p-value."""
    return 2.0 * betabinom.pmf(0, n, a, a)


def pmin_closed(n, a):
    """Closed form for integer a: 2*B(a, n+a)/B(a,a)."""
    from scipy.special import betaln

    return 2.0 * np.exp(betaln(a, n + a) - betaln(a, a))


def n_needed(a, alpha):
    """Smallest n whose most extreme outcome can reach p <= alpha."""
    n = 2
    while pmin_closed(n, a) > alpha and n < 10_000_000:
        n = int(n * 1.05) + 1
    return n


def two_sided_p(k, n, a):
    """Two-sided equal-tail-mass p under symmetric BetaBinom (mean n/2)."""
    k = min(k, n - k)
    return float(2.0 * betabinom.cdf(k, n, a, a))


print("=" * 78)
print("1. MINIMUM ATTAINABLE two-sided p under BetaBinom(n,a,a)  (K=0, the most extreme seed)")
print("=" * 78)
print(f"{'n':>8} " + " ".join(f"{'a=%d' % a:>12}" for a in (2, 3, 4)))
for n in (2, 10, 100, 1000, 1523, 3162, 10000):
    print(f"{n:8d} " + " ".join(f"{pmin_closed(n, a):12.3e}" for a in (2, 3, 4)))
print()
print("closed form check vs scipy at n=1523, a=3:", pmin_closed(1523, 3), pmin_two_sided(1523, 3))
print("asymptote  a=2: 12/n^2 ; a=3: 120/n^3 ; a=4: 1680/n^4")
for a, c in ((2, 12), (3, 120), (4, 1680)):
    print(f"  a={a}: exact(1523)={pmin_closed(1523, a):.4e}  approx={c / 1523.0**a:.4e}")

print()
print("=" * 78)
print("2. SEED SIZE a screen would need, for several multiplicity levels")
print("=" * 78)
print(f"{'level alpha':>14} {'a=2':>10} {'a=3':>10} {'a=4':>10}   (min n with any power)")
for label, alpha in (
    ("1/160366", 1 / 160366),
    ("1/863246", 1 / 863246),
    ("0.05/863246", 0.05 / 863246),
    ("1/1000", 1e-3),
    ("1/100", 1e-2),
    ("0.05", 0.05),
):
    print(f"{label:>14} " + " ".join(f"{n_needed(a, alpha):10d}" for a in (2, 3, 4)))

print()
print("the vcap top seed  N=1523 sense=5:")
for a in (2, 3, 4):
    print(f"   a={a}:  p = {two_sided_p(5, 1523, a):.3e}   (p_min attainable = {pmin_closed(1523, a):.3e})")

print()
print("=" * 78)
print("3. WEIGHTS: shipped pair-count vs the model's own information")
print("=" * 78)
# per-seed MoM  od_s = [(K-n/2)^2 - n/4] / [n(n-1)/4]  ; Var(K) = n/4*(1+(n-1)r)
# Var(od_s) ~ 2*Var(K)^2 / (n(n-1)/4)^2 = 2*[(1+(n-1)r)/(n-1)]^2   (normal approx)
# optimal weight ∝ 1/Var ∝ [(n-1)/(1+(n-1)r)]^2 ; SATURATES at 1/r^2 as n->inf.
def w_opt(n, r):
    return ((n - 1.0) / (1.0 + (n - 1.0) * r)) ** 2


def w_ship(n):
    return n * (n - 1.0) / 4.0


print(f"{'n':>7} {'shipped w':>14} {'rel to n=2':>12} | " + " | ".join(
    f"r={r}: {'w_opt':>9} {'rel':>9}" for r in (0.05, 0.2)
))
for n in (2, 5, 10, 100, 1000, 1523):
    row = f"{n:7d} {w_ship(n):14.1f} {w_ship(n) / w_ship(2):12.3e} | "
    row += " | ".join(
        f"r={r}: {w_opt(n, r):9.3f} {w_opt(n, r) / w_opt(2, r):9.3e}" for r in (0.05, 0.2)
    )
    print(row)
print()
for r in (0.05, 0.2):
    over = (w_ship(1523) / w_ship(2)) / (w_opt(1523, r) / w_opt(2, r))
    print(f"  at r={r}: the shipped estimator OVER-weights an n=1523 seed vs an n=2 seed by {over:,.0f}x")

print()
print("normalised 'effective independent Beta draws' per seed  d(n,r) = r^2 * w_opt(n,r) in (0,1]:")
print(f"{'n':>7} " + " ".join(f"{'r=%.2f' % r:>10}" for r in (0.01, 0.03, 0.05, 0.1, 0.2)))
for n in (2, 3, 5, 10, 50, 100, 1000, 10**6):
    print(f"{n:7d} " + " ".join(f"{r * r * w_opt(n, r):10.4f}" for r in (0.01, 0.03, 0.05, 0.1, 0.2)))

print()
print("=" * 78)
print("4. CONTAMINATION BIAS of the pooled MoM")
print("=" * 78)
# a fully displaced seed (all fragments one strand): excess = (n/2)^2 - n/4 = n(n-1)/4 = scale
# => its own od_s is EXACTLY 1.  So od_hat = (1-pi)*od + pi*1 with pi the WEIGHT share.
for n in (10, 100, 1523):
    k = 0
    ex = (k - n * 0.5) ** 2 - n * 0.25
    sc = n * (n - 1.0) / 4.0
    print(f"  n={n:5d} all-antisense seed: excess/scale = {ex / sc:.6f}")
print("  => od_hat ≈ od_true + pi*(1 - od_true), pi = contaminated share of the pooling WEIGHT")
print("  => 1 percentage point of contaminated weight  ==  +0.01 of overdispersion.  Direct.")
