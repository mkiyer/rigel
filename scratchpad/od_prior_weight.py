"""DERIVING W — the strand-overdispersion prior's weight, in PAIRS.

Part A: MC-verify the exact sampling variance of the pooled MoM estimator.
        CLAIM: Var(od_mom | od=0) = 1 / I   EXACTLY, with I = SUM_s n_s(n_s-1)/2 (pairs).
        => one pair = one unit of Fisher-style precision on od; pairs are not an analogy.
Part B: information SATURATION at od > 0 — a seed's information caps at i_max(od) pairs.
Part C: the candidates for W, applied to the real / synthetic seed sets.

Run: OMP_NUM_THREADS=1 python scratchpad/od_prior_weight.py
"""
from __future__ import annotations

import pickle

import numpy as np
from scipy.stats import betabinom

SEEDS = "/Users/mkiyer/proj/rigel/scratchpad/od_seeds.pkl"
EPS = 1e-12


# ────────────────────────────── the estimator's own per-seed terms ─────────────────────────────
def terms(sense, total, weight, kappa, mu_c):
    """excess_s, scale_s exactly as `gdna_strand._fit_overdispersion` forms them."""
    if weight is None:
        node_mean = np.full(total.shape, kappa)
        frac = np.ones_like(total)
    else:
        node_mean = 0.5 * weight + kappa * (1.0 - weight)
        frac = weight
    binom = total * node_mean * (1.0 - node_mean)
    excess = (sense - total * node_mean) ** 2 - binom
    nc = frac * total
    scale = np.maximum(nc * (nc - 1.0), 0.0) * (mu_c * (1.0 - mu_c))
    ok = total > 0.0
    return excess[ok], scale[ok], nc[ok]


# ───────────────────────────────── PART A: the variance identity ───────────────────────────────
def part_a(rng):
    print("=" * 100)
    print("PART A — Var(od_mom | od = 0) = 1/I, with I = SUM n(n-1)/2.   [MC verification]")
    print("=" * 100)
    print("  algebra:  e_s = (K-n/2)^2 - n/4,  Var(e_s)|_{od=0} = mu4 - mu2^2 = n(n-1)/8 = pairs_s/4")
    print("            scale_s = n(n-1)/4 = pairs_s/2")
    print("            Var(od_mom) = SUM Var(e)/ (SUM scale)^2 = (I/4)/(I/2)^2 = 1/I\n")
    print(f"  {'seed sizes':>34s} {'I = pairs':>12s} {'1/sqrt(I)':>11s} {'MC sd':>11s} {'ratio':>8s}")
    trials = 40000
    for label, sizes in (
        ("100 x n=2", np.full(100, 2)),
        ("1000 x n=2", np.full(1000, 2)),
        ("100 x n=10", np.full(100, 10)),
        ("50 x n=50", np.full(50, 50)),
        ("1 x n=1000", np.full(1, 1000)),
        ("mixed 1..40 (500 seeds)", np.arange(1, 41).repeat(13)[:500]),
    ):
        n = sizes.astype(float)
        pairs = float(np.sum(n * (n - 1) / 2))
        scale_sum = float(np.sum(n * (n - 1) / 4))
        k = rng.binomial(sizes[None, :], 0.5, size=(trials, sizes.size))
        e = (k - n / 2) ** 2 - n / 4
        od = e.sum(axis=1) / scale_sum
        print(
            f"  {label:>34s} {pairs:12,.0f} {1 / np.sqrt(pairs):11.5f} {od.std():11.5f} "
            f"{od.std() * np.sqrt(pairs):8.3f}"
        )
    print("\n  ratio ~ 1.000 on every row  =>  the identity is exact, not asymptotic.")


# ───────────────────────────── PART B: information saturation at od>0 ──────────────────────────
def exact_var_excess(n, od):
    """Var[(K - n/2)^2 - n/4] for K ~ BetaBinom(n, a, a), od = 1/(2a+1).  Exact pmf."""
    if od <= 0:
        return n * (n - 1) / 8.0
    a = 0.5 * (1.0 - od) / od
    k = np.arange(n + 1)
    p = betabinom.pmf(k, n, a, a)
    e = (k - n / 2.0) ** 2 - n / 4.0
    m1 = float(np.sum(p * e))
    m2 = float(np.sum(p * e * e))
    return m2 - m1 * m1


def i_max(od):
    """Asymptotic per-seed information cap, in pairs:  (1+2od) / (2 od^2 (1-od))."""
    return (1.0 + 2.0 * od) / (2.0 * od * od * (1.0 - od))


def part_b():
    print("\n" + "=" * 100)
    print("PART B — at od > 0 a seed's information SATURATES: i_s(od) = scale_s^2 / Var(e_s)")
    print("=" * 100)
    print("  large-n limit:  Var(e) -> n^4 * Var(u^2),  u = p-1/2 ~ Beta(a,a)-1/2")
    print("                  Var(u^2) = od^2 (1-od) / (8 (1+2od))")
    print("                  => i_max(od) = (1+2od) / (2 od^2 (1-od))    [pairs, per seed]\n")
    ods = [0.0, 0.01, 0.0345, 0.1, 0.143, 0.2]
    print(f"  {'n':>7s} {'pairs':>12s} " + "".join(f"{f'od={o:g}':>12s}" for o in ods))
    for n in (2, 5, 10, 50, 200, 1523):
        row = []
        for o in ods:
            v = exact_var_excess(n, o)
            sc = n * (n - 1) / 4.0
            row.append(sc * sc / max(v, EPS))
        print(f"  {n:7d} {n * (n - 1) / 2:12,.0f} " + "".join(f"{r:12,.1f}" for r in row))
    print(f"  {'limit':>7s} {'inf':>12s} " + "".join(f"{i_max(o) if o > 0 else np.inf:12,.1f}" for o in ods))
    print("\n  At the ceiling od=0.2 one seed is worth <= 22 pairs no matter how large it is.")
    print("  At the current prior od=0.0345 the cap is ~465 pairs.")


# ────────────────────────────────── PART C: the candidates for W ───────────────────────────────
def shrink(od_mom, info, w, od0, ceil=0.2):
    return float(np.clip((info * od_mom + w * od0) / (info + w), 0.0, ceil))


def part_c():
    print("\n" + "=" * 100)
    print("PART C — the candidates for W, on the production seed sets")
    print("=" * 100)
    d = pickle.load(open(SEEDS, "rb"))
    ceil = 0.2
    print(
        "\n  W candidates (all in PAIRS):\n"
        "    W_pseudo(a)  = a(2a-1)          Beta(a,a) as 2a pseudo-fragments -> pairs   [a=14 -> 378]\n"
        "    W_unif       = 12/od_max^2      max-entropy (uniform) prior on [0, od_max]  [-> 300]\n"
        "    W_eb         = 1/var(od_mom across libraries)                               [empirical Bayes]\n"
    )
    hdr = (
        f"  {'sample':>40s} {'n_seed':>9s} {'I=pairs':>13s} {'I_sat':>13s} {'I_sand':>10s} "
        f"{'od_mom':>8s} {'shipped':>8s} {'W=378':>8s} {'W=300':>8s} {'sand300':>8s}"
    )
    rows = {}
    for comp, mu in (("gdna", 0.5), ("rna", None)):
        print(f"\n  ---- {comp.upper()} ----")
        print(hdr)
        for name, g in d.items():
            if comp not in g:
                continue
            s, t, w, kappa, _ = g[comp]
            mu_c = 0.5 if comp == "gdna" else kappa
            ex, sc, nc = terms(s, t, w if comp == "gdna" else None, kappa, mu_c)
            inf = sc > 0
            ex, sc, nc = ex[inf], sc[inf], nc[inf]
            if ex.size == 0 or sc.sum() <= 0:
                continue
            od_mom = float(ex.sum() / sc.sum())
            pairs = nc * (nc - 1.0) / 2.0
            I = float(pairs.sum())
            # saturation-capped information, evaluated self-consistently at the fitted od
            od_ref = min(max(od_mom, 1e-4), ceil)
            I_sat = float(np.minimum(pairs, i_max(od_ref)).sum())
            # sandwich (Huber-White / GMM) variance of the ratio estimator
            resid = ex - od_mom * sc
            I_sand = float(sc.sum() ** 2 / max((resid**2).sum(), EPS))
            n_nodes = int(np.sum(nc > 0))
            shipped = shrink(od_mom, n_nodes, 30.0, 1.0 / 29.0, ceil)
            rows[(comp, name)] = (od_mom, I, I_sat, I_sand, n_nodes)
            print(
                f"  {name[:40]:>40s} {n_nodes:9,d} {I:13,.0f} {I_sat:13,.0f} {I_sand:10,.1f} "
                f"{od_mom:8.4f} {shipped:8.4f} "
                f"{shrink(od_mom, I, 378.0, 1 / 29, ceil):8.4f} "
                f"{shrink(od_mom, I, 300.0, 0.1, ceil):8.4f} "
                f"{shrink(od_mom, I_sand, 300.0, 0.1, ceil):8.4f}"
            )

    # empirical-Bayes tau^2 from the between-library spread of od_mom (real libraries only)
    print("\n  ---- empirical Bayes: tau^2 from the between-library spread of od_mom ----")
    for comp in ("gdna", "rna"):
        v = [r[0] for (c, n), r in rows.items() if c == comp and n.startswith("REAL")]
        if len(v) < 2:
            continue
        v = np.array(v)
        tau2 = float(v.var(ddof=1))
        print(
            f"    {comp:>5s}  od_mom over {v.size} libraries = "
            f"[{', '.join(f'{x:.4f}' for x in v)}]  tau^2={tau2:.6f}  tau={np.sqrt(tau2):.4f}  "
            f"W_eb = 1/tau^2 = {1 / max(tau2, EPS):,.0f} pairs"
        )

    # what W implies operationally
    print("\n  ---- what a given W means: the data outvotes the prior at I = W pairs ----")
    print(f"    {'W':>8s} {'tau=1/sqrt(W)':>14s} {'seeds of n=2 needed':>21s} {'one seed of n=':>16s}")
    for w in (22.0, 30.0, 300.0, 378.0, 1e4, 1e6):
        n_single = 0.5 * (1 + np.sqrt(1 + 8 * w))
        print(f"    {w:8,.0f} {1 / np.sqrt(w):14.4f} {w:21,.0f} {n_single:16,.0f}")


def main():
    rng = np.random.default_rng(0)
    part_a(rng)
    part_b()
    part_c()


if __name__ == "__main__":
    main()
