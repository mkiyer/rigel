"""Exact Beta-Binomial truncation machinery for the seed-screen bias derivation.

Symmetric BetaBinom(n, a, a) with mean n/2 and intraclass correlation od = 1/(2a+1).
"""
from __future__ import annotations

import numpy as np
from scipy.special import gammaln


def od_of_a(a: float) -> float:
    return 1.0 / (2.0 * a + 1.0)


def a_of_od(od: float) -> float:
    """a = (1-od)/(2 od); od=0 -> inf (Binomial limit)."""
    return np.inf if od <= 0 else 0.5 * (1.0 - od) / od


def bb_logpmf(n: int, a: float) -> np.ndarray:
    """log P(K=k) for k=0..n under symmetric BetaBinom(n,a,a). a=inf -> Binomial(n,1/2)."""
    k = np.arange(n + 1, dtype=np.float64)
    lchoose = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    if not np.isfinite(a):
        return lchoose - n * np.log(2.0)
    lb = gammaln(k + a) + gammaln(n - k + a) - gammaln(n + 2 * a)
    lb0 = gammaln(a) + gammaln(a) - gammaln(2 * a)
    return lchoose + lb - lb0


def bb_pmf(n: int, a: float) -> np.ndarray:
    lp = bb_logpmf(n, a)
    return np.exp(lp - lp.max()) / np.exp(lp - lp.max()).sum()


def two_sided_p(n: int, a_null: float) -> np.ndarray:
    """t[k] = P(|K - n/2| >= |k - n/2|) under BetaBinom(n, a_null, a_null). Exact, symmetric."""
    pmf = bb_pmf(n, a_null)
    k = np.arange(n + 1)
    d = np.abs(k - n / 2.0)
    # cumulative from the outside in, over distinct distance levels
    order = np.argsort(-d, kind="stable")
    csum = np.cumsum(pmf[order])
    t = np.empty(n + 1)
    ds = d[order]
    # for ties (k and n-k share d) both get the same inclusive tail
    i = 0
    while i <= n:
        j = i
        while j + 1 <= n and ds[j + 1] == ds[i]:
            j += 1
        t[order[i:j + 1]] = csum[j]
        i = j + 1
    return np.minimum(t, 1.0)


def keep_mask(n: int, a_null: float, alpha: float) -> np.ndarray:
    """True for k values KEPT (two-sided tail p >= alpha)."""
    return two_sided_p(n, a_null) >= alpha


def min_pvalue(n: int, a_null: float) -> float:
    """Smallest achievable two-sided p at size n = 2*P(K=0) = 2*B(a, a+n)/B(a,a).

    Closed form: a=2 -> 12/((n+2)(n+3)); a=3 -> 120/((n+3)(n+4)(n+5)).
    """
    lp0 = gammaln(a_null) + gammaln(n + a_null) - gammaln(n + 2 * a_null) - (
        gammaln(a_null) + gammaln(a_null) - gammaln(2 * a_null)
    )
    return float(2.0 * np.exp(lp0))


def seed_moments(n: int, od_true: float, a_null: float, alpha: float):
    """Exact per-seed truncated moments under true dispersion od_true.

    Returns (p_keep, E_excess_all, E_excess_kept_unnormalised, r) where
      E_excess_all = scale_n * od_true   (the untruncated expectation, exact)
      E_excess_kept_unnorm = E[excess * 1{keep}]
      r = E[excess|keep] / E[excess]     (the task's requested ratio)
    """
    scale = n * (n - 1) / 4.0
    a_t = a_of_od(od_true)
    pmf = bb_pmf(n, a_t)
    k = np.arange(n + 1, dtype=np.float64)
    dev2 = (k - n / 2.0) ** 2
    keep = keep_mask(n, a_null, alpha)
    p_keep = float(pmf[keep].sum())
    T = float((dev2[keep] * pmf[keep]).sum())
    e_kept = T - (n / 4.0) * p_keep          # E[excess * 1{keep}]
    e_all = scale * od_true                   # exact: Var - n/4 = scale*od
    r = np.nan if e_all == 0 else (e_kept / p_keep) / e_all if p_keep > 0 else np.nan
    return p_keep, e_all, e_kept, r
