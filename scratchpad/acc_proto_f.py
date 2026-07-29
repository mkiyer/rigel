"""EXPERIMENT F — gDNA FL 50 vs RNA FL 200: where does each object's information live?

The owner's concern: with a 4x FL contrast, short gDNA fragments cross far fewer
boundaries, so edges see almost no gDNA and the edge data is "sparse" for it.

Question 1: is that a BIAS or a SELECTION effect the effective length already prices?
Question 2: how much precision does each object class actually have for each component?
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260731)


def fl(mu, sd, n=1600):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


def edge_design(pmf):
    w = np.arange(pmf.size, dtype=np.float64)
    return np.array([np.sum(pmf * w), np.sum(pmf * w * w)])


def node_design(pmf, B):
    w = np.arange(pmf.size, dtype=np.float64)
    fit = np.maximum(0.0, B - w + 1.0)
    return np.array([np.sum(pmf * fit), np.sum(pmf * w * fit)])


def report(label, M, rho, scale, reps=600):
    """M columns = (gDNA, RNA); rows = (count eq, sum-length eq). rho in molecules/bp."""
    lam = M[0] * rho * scale                       # expected counts per component
    tot = lam.sum()
    fg_true = rho[0] / rho.sum()
    obs_share = lam[0] / max(tot, 1e-12)           # gDNA share of the OBSERVED events
    if tot < 1e-9 or abs(np.linalg.det(M)) < 1e-9:
        print(f"  {label:34s} events {tot:8.1f}  gDNA share {obs_share:6.3f}   SINGULAR/EMPTY")
        return
    est = []
    for _ in range(reps):
        ng, nr = RNG.poisson(lam[0]), RNG.poisson(lam[1])
        wg = RNG.choice(PG.size, size=ng, p=_ctx(PG, M, 0))
        wr = RNG.choice(PR.size, size=nr, p=_ctx(PR, M, 1))
        o = np.array([ng + nr, wg.sum() + wr.sum()]) / scale
        s = np.maximum(np.linalg.solve(M, o), 0.0)
        if s.sum() > 0:
            est.append(s[0] / s.sum())
    e = np.asarray(est)
    print(f"  {label:34s} events {tot:8.1f}  gDNA share {obs_share:6.3f}   "
          f"f_g {e.mean():.3f} +/- {e.std():.3f}  bias {e.mean()-fg_true:+.4f}")


CTX = {}


def _ctx(pmf, M, col):
    key = (id(pmf), CTX.get("mode"), CTX.get("B"))
    if key in CTX:
        return CTX[key]
    w = np.arange(pmf.size, dtype=np.float64)
    q = pmf * w if CTX["mode"] == "edge" else pmf * np.maximum(0.0, CTX["B"] - w + 1.0)
    q = q / q.sum() if q.sum() > 0 else pmf
    CTX[key] = q
    return q


PG, PR = fl(50, 12), fl(200, 50)

if __name__ == "__main__":
    print("gDNA FL 50 (sd 12)   RNA FL 200 (sd 50)\n")
    w = np.arange(PG.size, dtype=np.float64)
    print(f"  E_g[w] = {np.sum(PG*w):.1f}   E_r[w] = {np.sum(PR*w):.1f}"
          f"   -> an edge sees RNA {np.sum(PR*w)/np.sum(PG*w):.1f}x over-represented per molecule\n")

    # cfRNA-like mixture: gDNA is 15 % of MOLECULES
    rho = np.array([0.15, 0.85]) * 0.02
    fg_true = rho[0] / rho.sum()
    print(f"TRUE composition: f_g = {fg_true:.3f}  (molecules per bp)")
    print("SCALE: enough genome that each object class gets a realistic event count\n")

    CTX.clear()
    CTX["mode"] = "edge"
    CTX["B"] = None
    Me = np.column_stack([edge_design(PG), edge_design(PR)])
    print("EDGE  (one genomic point)")
    for scale, lbl in ((50, "thin  (~13 events)"), (500, "typical (~130)"), (5000, "deep (~1300)")):
        report(f"edge, {lbl}", Me, rho, scale)

    print("\nNODE, contained (by region length)")
    for B in (100, 147, 300, 1000, 5000):
        CTX.clear()
        CTX["mode"] = "node"
        CTX["B"] = B
        Mn = np.column_stack([node_design(PG, B), node_design(PR, B)])
        report(f"node L={B}", Mn, rho, 500)

    print("\nCOMPONENT SENSITIVITY  E_g / E_r  (>1 = gDNA-favouring object)")
    print(f"  edge                        {edge_design(PG)[0]/edge_design(PR)[0]:8.3f}")
    for B in (100, 147, 300, 1000, 5000):
        eg, er = node_design(PG, B)[0], node_design(PR, B)[0]
        r = eg / er if er > 1e-9 else float("inf")
        print(f"  node L={B:5d}                {r:8.3f}")
