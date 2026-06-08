"""Study: normal-moment vs exact two-component BB-mixture strand likelihood.

Decides docs/em_strand/03 §2.1: at what fragment count n does the normal-moment
approximation's gdna_frac posterior diverge from the exact mixture marginal? Below that
threshold use exact; above, normal-moment.
"""

import numpy as np
from scipy.stats import betabinom, binom

KAPPA_RNA = 0.95
GRID = np.linspace(1e-4, 1 - 1e-4, 400)


def normal_moment_ll(kplus, n, f, od):
    p = 0.5 * f + KAPPA_RNA * (1 - f)
    var = max(n * p * (1 - p) + (n * f) ** 2 * 0.25 * od, 1e-9)
    return -0.5 * (kplus - n * p) ** 2 / var - 0.5 * np.log(var)


def exact_mixture_logpmf(kplus, n, f, od):
    a = b = (1e9 if od < 1e-9 else 0.5 * (1 - od) / od)
    total = 0.0
    for ng in range(n + 1):
        w = binom.pmf(ng, n, f)
        if w < 1e-12:
            continue
        conv = np.convolve(
            betabinom.pmf(np.arange(ng + 1), ng, a, b),
            binom.pmf(np.arange(n - ng + 1), n - ng, KAPPA_RNA),
        )
        if kplus < len(conv):
            total += w * conv[kplus]
    return np.log(max(total, 1e-300))


def median_frac(ll_fn, kplus, n, od):
    ll = np.array([ll_fn(kplus, n, f, od) for f in GRID])
    w = np.exp(ll - ll.max())
    w /= w.sum()
    return float(np.interp(0.5, np.cumsum(w), GRID))


print(f'{"n":>5} {"od":>5} {"tf":>4} {"k+/n":>9} {"normal":>8} {"exact":>8} {"d":>7}')
maxd = {}
for n in [8, 20, 50, 150, 500]:
    for od in [0.05, 0.1, 0.2]:
        for tf in [1.0, 0.5]:
            p = 0.5 * tf + KAPPA_RNA * (1 - tf)
            kplus = int(round(n * p))
            mn = median_frac(normal_moment_ll, kplus, n, od)
            ex = median_frac(exact_mixture_logpmf, kplus, n, od)
            maxd[n] = max(maxd.get(n, 0), abs(mn - ex))
            print(f"{n:>5} {od:>5.2f} {tf:>4.1f} {kplus:>4}/{n:<4} {mn:>8.3f} {ex:>8.3f} {mn-ex:>+7.3f}")
print("\nmax |d median gdna_frac| by n:", {k: round(v, 3) for k, v in maxd.items()})
print()
for od in [0.0, 0.05, 0.1, 0.2]:
    mn = median_frac(normal_moment_ll, 126, 152, od)
    ex = median_frac(exact_mixture_logpmf, 126, 152, od)
    print(f"126/152 od={od:.2f}: normal={mn:.3f}  exact={ex:.3f}")
