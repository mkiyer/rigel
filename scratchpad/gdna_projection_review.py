"""PER-SCENARIO REVIEW of the PROJECTION — does it pull an observed density to the RIGHT gDNA level?

The projection is the endpoint: node's observed gDNA density `d = log10(g_hat/E)` -> anchor `mu*`, which is
what the solver gets pulled toward. This asks the only question that matters — **is mu*(d) the right
function of d?**

THE YARDSTICK. The projection is an estimator of the population conditional mean

        IDEAL(d) = E[ true log10 rho_g | observed log10 rho_g = d ]

so the honest test is to plot the projection's transfer curve `mu*(d)` against IDEAL(d) measured from the
oracle, on the same axis, per scenario. Three reference lines make the reading unambiguous:
  * IDENTITY      mu* = d           — do nothing (the projection is inert)
  * IDEAL         E[truth | d]      — the best any projection of d alone can do
  * mu*(d)        the candidate     — where it sits between them IS its score

`proj_bias` = mean(mu* - IDEAL) over live nodes, weighted by fragment mass (want ~0; negative = still
under-calling), and `proj_gain` = the fraction of the identity->IDEAL gap that the projection closes
(1.0 = perfect, 0.0 = inert, <0 = moves the wrong way).

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_projection_review.py [--suite ambig] [--out FIG.png]
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_projection import project_asym  # noqa: E402
from gdna_resurrect import load_recipe  # noqa: E402

GRID = L.GRID
PROJECTIONS = (
    ("symmetric h=0.15", lambda P, d: L.project(P, d, 0.15, "mean")),
    ("asymmetric", project_asym),
)


def node_view(s):
    """Live nodes with their observed / true gDNA log10-density and a mass weight.

    ⚠ BOTH densities are floored at the ONE-COUNT RESOLUTION WALL `1/E`, not at some epsilon. A node with
    zero gDNA has a true density of `-inf`; the honest statement the data supports is "gDNA <= 1 count
    here", i.e. `log10(1/E)`. This is the fit's own convention (`npmle._kde_density` centres kernels at
    `log(max(g,1))-log(E)`), so the metric and the object being measured agree.

    Flooring at 1e-9 instead puts every zero-gDNA node at -9..-13 decades and makes E[truth|obs] on the
    `none_*` family a meaningless -12, which reads out as a +8..+13 decade projection bias."""
    live = (s["eff"] > 1e-9) & (s["mass"] > 1e-12)
    eff = np.maximum(s["eff"], 1e-9)
    obs = np.log10(np.maximum(s["g_hat"], 1.0) / eff)
    tru = np.log10(np.maximum(s["G"], 1.0) / eff)
    return live, obs, tru, s["mass"]


def ideal_curve(obs, tru, w, edges):
    """IDEAL(d) = E[truth | obs = d], mass-weighted, binned. NaN where a bin has too little mass."""
    idx = np.digitize(obs, edges) - 1
    out = np.full(len(edges) - 1, np.nan)
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() >= 3 and w[m].sum() > 0:
            out[b] = float(np.average(tru[m], weights=w[m]))
    return out


def score(s, P, proj, edges):
    """proj_bias (mu* - IDEAL) and proj_gain (fraction of the identity->IDEAL gap closed), mass-weighted."""
    live, obs, tru, w = node_view(s)
    if P is None or live.sum() < 5:
        return np.nan, np.nan
    o, t, ww = obs[live], tru[live], w[live]
    mu = proj(P, o)
    # per-node IDEAL via the binned curve, so the score is against the achievable target, not the truth
    ic = ideal_curve(o, t, ww, edges)
    b = np.clip(np.digitize(o, edges) - 1, 0, len(ic) - 1)
    ideal = ic[b]
    ok = np.isfinite(ideal) & (ww > 0)
    if ok.sum() < 5:
        return np.nan, np.nan
    bias = float(np.average((mu - ideal)[ok], weights=ww[ok]))
    gap = ideal - o                                   # what the identity leaves on the table
    denom = float(np.average(np.abs(gap)[ok], weights=ww[ok]))
    gain = float(np.average((np.sign(gap) * (mu - o))[ok], weights=ww[ok])) / denom if denom > 1e-9 else np.nan
    return bias, gain


def review(suite="ambig", out=None):
    rec_fn = load_recipe()
    scen = sorted(L.load_scenarios(suite),
                  key=lambda s: (s["group"][0], s["group"][2], s["group"][1], s["group"][3]))
    edges = np.linspace(-5, 2.0, 29)
    mid = 0.5 * (edges[1:] + edges[:-1])

    print(f"=== {suite}: PROJECTION review — mu*(d) vs IDEAL(d) = E[truth | obs=d] ===")
    print(f"{'condition':50s} | " + " | ".join(f"{n:>22s}" for n, _ in PROJECTIONS))
    print(f"{'':50s} | " + " | ".join(f"{'bias':>10s} {'gain':>11s}" for _ in PROJECTIONS))
    agg = {n: ([], []) for n, _ in PROJECTIONS}
    rows = []
    for s in scen:
        P = rec_fn(s)
        line, cols = f"{s['cond']:50s} | ", []
        for name, pf in PROJECTIONS:
            b, g = score(s, P, pf, edges)
            cols.append(f"{b:+10.3f} {g:+11.2f}")
            if np.isfinite(b):
                agg[name][0].append(b)
                agg[name][1].append(g)
        print(line + " | ".join(cols))
        rows.append((s, P))
    print()
    for name, _ in PROJECTIONS:
        b, g = agg[name]
        print(f"  {name:22s} mean bias {np.mean(b):+.3f}   mean gain {np.nanmean(g):+.2f}   "
              f"(bias 0 = lands on IDEAL; gain 1.0 = closes the whole identity->IDEAL gap, 0 = inert)")

    if out is None:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(rows)
    ncol, nrow = 4, int(np.ceil(n / 4))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.4 * nrow), sharex=True, sharey=True)
    for ax, (s, P) in zip(axes.ravel(), rows):
        live, obs, tru, w = node_view(s)
        o, t, ww = obs[live], tru[live], w[live]
        ax.scatter(o, t, s=np.clip(ww / max(ww.max(), 1) * 60, 1, 60), c="0.72", lw=0, alpha=.5,
                   label="nodes (obs, truth)")
        ax.plot([-5, 2], [-5, 2], "k--", lw=0.9, label="IDENTITY (inert)")
        ic = ideal_curve(o, t, ww, edges)
        ax.plot(mid, ic, "C2", lw=2.2, label="IDEAL E[truth|obs]")
        if P is not None:
            for (name, pf), c in zip(PROJECTIONS, ("C0", "C3")):
                ax.plot(mid, pf(P, mid), c, lw=1.5, label=name)
        cap, dna, ss, nr = s["group"]
        ax.set_title(f"{cap} · {dna} · ss{ss} · nrna={nr}", fontsize=8)
        ax.set_xlim(-5, 2.0)
        ax.set_ylim(-5, 2.0)
        ax.tick_params(labelsize=7)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    axes.ravel()[0].legend(fontsize=6.5, loc="upper left")
    fig.suptitle(f"PROJECTION transfer function — every {suite} scenario. "
                 f"A good projection lies ON the green IDEAL curve, not on the dashed identity.", fontsize=13)
    fig.supxlabel("observed log10 rho_g  (pass-0)")
    fig.supylabel("projected mu*  /  true log10 rho_g")
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.985))
    fig.savefig(out, dpi=100)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default="ambig")
    ap.add_argument("--out",
                    default="/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_projection_review.png")
    a = ap.parse_args()
    review(a.suite, a.out)
