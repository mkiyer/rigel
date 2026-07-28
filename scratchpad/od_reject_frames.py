"""WHICH REJECTION FRAME? Bonferroni / BH / Tarone / robust weighting, measured on the 4 real libs.

Uses the seed arrays dumped by od_reject_dump.py (exact production inputs).
Run: OMP_NUM_THREADS=1 python scratchpad/od_reject_frames.py
"""

from __future__ import annotations

import pickle
import sys

import numpy as np
from scipy.special import betaln
from scipy.stats import betabinom

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")

SEEDS = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
_EPS = 1e-12
CEIL = 0.2
PRIOR_OD = 1.0 / 29.0  # Beta(14,14)
PRIOR_W = 30.0


def ab_from(mu, rho):
    """BetaBinom (alpha, beta) for mean mu, intraclass correlation rho."""
    c = (1.0 - rho) / rho
    return mu * c, (1.0 - mu) * c


def pmin_two_sided(n, mu, rho):
    """Smallest attainable two-sided tail p: 2*min(P(K=0), P(K=n))  -> at mu=1/2, 2*P(K=0)."""
    a, b = ab_from(mu, rho)
    p0 = np.exp(betaln(a, n + b) - betaln(a, b))
    pn = np.exp(betaln(n + a, b) - betaln(a, b))
    return 2.0 * np.minimum(p0, pn)


def two_sided_p(k, n, mu, rho):
    """Two-sided p = 2*min(P(K<=k), P(K>=k)) clipped at 1 (vectorised over unique (k,n))."""
    a, b = ab_from(mu, rho)
    lo = betabinom.cdf(k, n, a, b)
    hi = betabinom.sf(k - 1, n, a, b)
    return np.minimum(1.0, 2.0 * np.minimum(lo, hi))


def load():
    with open(SEEDS, "rb") as fh:
        return pickle.load(fh)


def terms(sense, total, weight, kappa):
    node_mean = 0.5 * weight + kappa * (1.0 - weight)
    binv = total * node_mean * (1.0 - node_mean)
    excess = (sense - total * node_mean) ** 2 - binv
    nc = weight * total
    scale = np.maximum(nc * (nc - 1.0), 0.0) * 0.25
    return excess, scale, nc, node_mean


def shrink(mom, w_data, prior_od=PRIOR_OD, prior_w=PRIOR_W):
    tw = w_data + prior_w
    return float(np.clip((w_data * mom + prior_w * prior_od) / tw, 0.0, CEIL))


def w_opt(nc, rho):
    """Inverse-variance weight for the per-seed MoM od_s under intraclass correlation rho."""
    return ((nc - 1.0) / (1.0 + (nc - 1.0) * rho)) ** 2


def info_weighted_od(excess, scale, nc, rho0, iters=25):
    """Fixed point: od = sum(w*od_s)/sum(w), w = inverse variance evaluated at the current od."""
    ok = scale > 0
    ods = excess[ok] / scale[ok]
    n = nc[ok]
    rho = rho0
    for _ in range(iters):
        w = w_opt(n, max(rho, 1e-4))
        new = float(np.sum(w * ods) / np.sum(w))
        if abs(new - rho) < 1e-9:
            rho = new
            break
        rho = min(max(new, 0.0), CEIL)
    return rho, ods, n


def eff_draws(nc, rho):
    """Normalised information: rho^2 * w_opt in (0,1] — 'independent Beta draws' per seed."""
    return rho * rho * w_opt(nc, rho)


def tarone_threshold(pmin, alpha):
    """Tarone (1990): restrict the multiplicity to tests that COULD reach the threshold.

    K* = min{k : #(pmin <= alpha/k) <= k};  threshold = alpha/K*.
    """
    order = np.sort(pmin)
    m = order.size
    for k in range(1, m + 1):
        t = alpha / k
        if np.searchsorted(order, t, side="right") <= k:
            return t, k
    return alpha / m, m


def bh_threshold(p, q):
    ps = np.sort(p)
    m = ps.size
    crit = q * np.arange(1, m + 1) / m
    ok = np.nonzero(ps <= crit)[0]
    return (ps[ok[-1]], ok[-1] + 1) if ok.size else (0.0, 0)


def report(tag, sense, total, weight, kappa):
    excess, scale, nc, node_mean = terms(sense, total, weight, kappa)
    ok = (total > 0) & (nc > 0)
    excess, scale, nc, node_mean, sense, total = (
        excess[ok],
        scale[ok],
        nc[ok],
        node_mean[ok],
        sense[ok],
        total[ok],
    )
    m = excess.size
    pooled = float(np.sum(excess) / np.sum(scale))
    shipped = shrink(pooled, m)
    print(f"\n{'=' * 100}\n{tag}\n{'=' * 100}")
    print(
        f"  seeds={m:8d}  pairs=Sum n(n-1)/2 = {np.sum(scale) * 2:,.0f}   "
        f"median n_c={np.median(nc):.1f}  max n_c={nc.max():.0f}"
    )
    print(f"  pooled MoM = {pooled:.4f}   shipped (shrunk by seed COUNT, clipped) = {shipped:.4f}")

    # leverage
    o = np.argsort(-np.abs(excess))
    tot_ex = np.sum(excess)
    print(
        f"  LEVERAGE: top seed = {excess[o[0]] / tot_ex:7.2%} of net numerator "
        f"(n={total[o[0]]:.0f}, sense={sense[o[0]]:.0f}); top-10 = "
        f"{np.sum(excess[o[:10]]) / tot_ex:7.2%}; top-0.1% = "
        f"{np.sum(excess[o[: max(1, m // 1000)]]) / tot_ex:7.2%}"
    )

    # ---- p-values under the CEILING null Beta(2,2) and a=3 -----------------------------------
    ni = np.rint(nc).astype(np.int64)
    ki = np.rint(np.where(sense > total * node_mean, sense, sense)).astype(np.int64)
    ki = np.clip(np.rint(sense * np.where(total > 0, nc / np.maximum(total, _EPS), 1.0)), 0, ni)
    ki = ki.astype(np.int64)
    keep = ni >= 2
    for a in (2.0, 3.0):
        rho = 1.0 / (2.0 * a + 1.0)
        pm = np.ones(m)
        pm[keep] = pmin_two_sided(ni[keep], 0.5, rho)
        # unique (n,k) -> p
        key = ni.astype(np.int64) * 100003 + ki.astype(np.int64)
        uk, inv = np.unique(key, return_inverse=True)
        un = (uk // 100003).astype(np.int64)
        ukk = (uk % 100003).astype(np.int64)
        pu = np.where(un >= 2, two_sided_p(ukk, un, 0.5, rho), 1.0)
        p = pu[inv]

        print(f"\n  --- null = Beta({a:.0f},{a:.0f})   (od={rho:.4f}) ---")
        print(f"      seeds that CAN EVER reach p<1/m ({1 / m:.2e}): {int(np.sum(pm <= 1 / m)):d} of {m}")
        rows = []
        # (a) Bonferroni FWER 0.05 ; (c) expected-1-false ; BH q ; Tarone
        for lbl, t in (
            (f"Bonferroni a=0.05  t={0.05 / m:.2e}", 0.05 / m),
            (f"expect 1 false     t={1.0 / m:.2e}", 1.0 / m),
        ):
            rows.append((lbl, t))
        for q in (0.001, 0.05, 1.0):
            tb, kb = bh_threshold(p, q)
            rows.append((f"BH q={q:<5g}       t={tb:.2e}", tb))
        for alpha in (0.05, 1.0):
            tt, kk = tarone_threshold(pm, alpha)
            rows.append((f"TARONE a={alpha:<4g} (m_eff={kk})  t={tt:.2e}", tt))
        print(
            f"      {'criterion':46s} {'#rej':>7} {'%seeds':>8} {'%pairs':>8} "
            f"{'%numer':>8} {'od_pooled':>10} {'od_shrunk':>10}"
        )
        for lbl, t in rows:
            rej = p < t
            kept = ~rej
            if np.sum(scale[kept]) <= 0:
                continue
            newmom = float(np.sum(excess[kept]) / np.sum(scale[kept]))
            print(
                f"      {lbl:46s} {int(rej.sum()):7d} {rej.mean():8.3%} "
                f"{np.sum(scale[rej]) / np.sum(scale):8.2%} "
                f"{np.sum(excess[rej]) / tot_ex:8.2%} {newmom:10.4f} "
                f"{shrink(newmom, int(kept.sum())):10.4f}"
            )

    # ---- the ROBUST-WEIGHT alternative: no rejection at all ------------------------------------
    print("\n  --- (e) NO rejection: inverse-variance (information) weighting instead ---")
    for start in (PRIOR_OD, CEIL):
        rho, ods, n = info_weighted_od(excess, scale, nc, start)
        w = w_opt(n, max(rho, 1e-4))
        d = eff_draws(n, max(rho, 1e-4))
        lev = np.max(w) / np.sum(w)
        print(
            f"      start={start:.4f} -> od_infoweighted = {rho:.4f}   "
            f"Sum effective Beta-draws = {np.sum(d):10.1f}   max single-seed weight share = {lev:.2e}"
        )
        print(
            f"           prior weight 30 draws would then be "
            f"{30.0 / (np.sum(d) + 30.0):6.2%} of the posterior weight   "
            f"(shipped: {PRIOR_W / (m + PRIOR_W):.4%})"
        )
    return None


def main():
    data = load()
    for name in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        g = data[name]["gdna"]
        sense, total, weight, kappa, kw = g
        report(f"{name}   gDNA seeds   (kappa={kappa:.4f})", sense, total, weight, kappa)


if __name__ == "__main__":
    main()
