"""EXPERIMENT D — does (count, total_length) identify gDNA vs RNA at a single object?

The owner's v3 proposal keeps the current node/edge architecture and stores, per object,
an integer COUNT and an integer TOTAL LENGTH.  Claim under test: those two numbers are
not just a density convenience -- they are a complete 2-component deconvolution, because
the two components have different fragment-length distributions.

EDGE (a genomic point p).  Fragments of length w crossing p have rate rho*w*f(w):
    E[flux]   = rho * E[w]
    E[sum w]  = rho * E[w^2]
NODE (contained in a region of length B).  Starts that fit: max(0, B-w+1):
    E[N]      = rho * sum_w f(w) * max(0, B-w+1)
    E[sum w]  = rho * sum_w f(w) * w * max(0, B-w+1)

Two equations, two unknowns, per object. No strand, no prior, no paths, no histogram.
"""

from __future__ import annotations

import numpy as np

RNG = np.random.default_rng(20260729)


def fl(mu, sd, n=1500):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[:10] = 0.0
    return p / p.sum()


def design_edge(pmf):
    w = np.arange(pmf.size, dtype=np.float64)
    return np.array([np.sum(pmf * w), np.sum(pmf * w * w)])


def design_node(pmf, B):
    w = np.arange(pmf.size, dtype=np.float64)
    fit = np.maximum(0.0, B - w + 1.0)
    return np.array([np.sum(pmf * fit), np.sum(pmf * w * fit)])


def solve2(obs, M):
    """min ||M rho - obs|| with rho >= 0 (2x2; direct then clip)."""
    try:
        r = np.linalg.solve(M, obs)
    except np.linalg.LinAlgError:
        return np.array([np.nan, np.nan])
    return np.maximum(r, 0.0)


def run(kind, gmu, rmu, B=None, rho_g=0.02, rho_r=0.05, reps=400, scale=3000.0):
    fg, fr = fl(gmu, gmu * 0.25), fl(rmu, rmu * 0.25)
    if kind == "edge":
        cg, cr = design_edge(fg), design_edge(fr)
    else:
        cg, cr = design_node(fg, B), design_node(fr, B)
    M = np.column_stack([cg, cr])            # rows = (count eq, sum-w eq)
    lam_g, lam_r = rho_g * cg[0] * scale, rho_r * cr[0] * scale
    fg_true = rho_g / (rho_g + rho_r)

    est = []
    for _ in range(reps):
        ng, nr = RNG.poisson(lam_g), RNG.poisson(lam_r)
        # draw the observed lengths from each component's length-biased-in-context law
        wg = RNG.choice(fg.size, size=ng, p=_ctx(fg, kind, B))
        wr = RNG.choice(fr.size, size=nr, p=_ctx(fr, kind, B))
        obs = np.array([ng + nr, wg.sum() + wr.sum()]) / scale
        r = solve2(obs, M)
        if r.sum() > 0:
            est.append(r[0] / r.sum())
    est = np.asarray(est)
    cond = np.linalg.cond(M)
    lbl = f"{kind:5s} " + (f"B={B:5d}" if B else "point ")
    print(f"  {lbl} gDNA {gmu:3.0f} / RNA {rmu:3.0f} | cond {cond:8.2f} | "
          f"f_g true {fg_true:.3f} est {est.mean():.3f} +/- {est.std():.3f} "
          f"| bias {est.mean()-fg_true:+.4f}")
    return cond


def _ctx(pmf, kind, B):
    """The length distribution of the fragments this object actually sees."""
    w = np.arange(pmf.size, dtype=np.float64)
    q = pmf * w if kind == "edge" else pmf * np.maximum(0.0, B - w + 1.0)
    s = q.sum()
    return q / s if s > 0 else pmf


if __name__ == "__main__":
    print("=" * 82)
    print("EDGE (flux, sum-of-lengths at a genomic point) -- region-size independent")
    print("=" * 82)
    for rmu in (200, 150, 120, 100, 80, 65, 60):
        run("edge", 60, rmu)

    print()
    print("=" * 82)
    print("NODE (contained count, sum-of-lengths) -- vs region length")
    print("=" * 82)
    for B in (100, 150, 300, 600, 1500, 5000):
        run("node", 60, 200, B=B)

    print()
    print("=" * 82)
    print("SENSITIVITY: what if the FL pmfs are MIS-SPECIFIED by 10 %?")
    print("=" * 82)
    fgood, fbad = fl(60, 15), fl(66, 16.5)
    fr = fl(200, 50)
    Mtrue = np.column_stack([design_edge(fgood), design_edge(fr)])
    Mbad = np.column_stack([design_edge(fbad), design_edge(fr)])
    rho = np.array([0.02, 0.05])
    obs = Mtrue @ rho
    r = solve2(obs, Mbad)
    print(f"  true f_g 0.286 -> with a 10 % gDNA-FL error: {r[0]/r.sum():.3f} "
          f"(bias {r[0]/r.sum()-0.286:+.4f})")
