"""W3 — THE PROJECTION: the posterior form, and a measurement of variance honesty.

Scored against the only defensible target: IDEAL(d) = E[true log10 rho_g | observed = d], measured from the
oracle. That is the best any function of the observation alone can do, so `gain` = the fraction of the
identity->IDEAL gap a projection closes (1 = perfect, 0 = inert, <0 = moves the wrong way).

Established in `gdna_hyperprior_resurrection.md` §4b: the IDEAL curve is largely FLAT (the observation is
often uninformative, so the job is SHRINKAGE, not translation); the symmetric projection sits on the
identity (inert); the asymmetric one is identity + a constant and is negative on 19/32 conditions; the
POSTERIOR form with sigma = the node's own belief width is best on both averages and never harmful, with no
new constant.

This re-runs it on the W2 kNN landscape and adds the measurement that explains the residual:

  THE SIGMA-CALIBRATION.  gain is ~+0.2, not 1.0. Either the landscape is wrong or `sigma_obs` is
  mis-calibrated. Scale sigma by c and find the c* that maximises gain:
      c* ~ 1   => declared variances are honest IN THE UNITS THE PRIOR CONSUMES; the residual is the
                  landscape's, and W2 is the lever.
      c* >> 1  => pass-0 is over-confident by that factor — a direct measurement of variance honesty on a
                  new axis, feeding the existing z2 programme.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_w3_projection.py [--plot]
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_projection import project_asym  # noqa: E402
from gdna_projection_review import ideal_curve, node_view  # noqa: E402

GRID = L.GRID
LN10 = np.log(10.0)
EDGES = np.linspace(-5, 2.0, 29)
MID = 0.5 * (EDGES[1:] + EDGES[:-1])
KNN = 0.5                      # the W2 kernel
C_GRID = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)


def sigma_obs(s, c=1.0):
    """The node's OWN belief width in decades — sqrt(Var(log f_g))/ln10. Already computed by the solver;
    no new constant. `c` scales it for the calibration measurement."""
    v = np.nan_to_num(s["var"], nan=8.0, posinf=8.0)
    return c * np.sqrt(np.maximum(v, 0.0)) / LN10


def project_posterior(P, d, sig):
    """A proper Bayesian update: landscape prior P(mu) x node likelihood N(d; mu, sigma_obs). Floored at
    the grid's own resolution, which is the narrowest width this axis can represent (derived, not tuned)."""
    logP = np.log(np.maximum(P, 1e-30))
    sd = np.maximum(np.asarray(sig, float), L.GRID_H)[:, None]
    z = (np.asarray(d, float)[:, None] - GRID[None, :]) / sd
    lr = logP[None, :] - 0.5 * z * z
    lr -= lr.max(1, keepdims=True)
    r = np.exp(lr)
    r /= np.maximum(r.sum(1, keepdims=True), 1e-30)
    return (r * GRID[None, :]).sum(1)


def gain_bias(s, P, mu_fn):
    """(gain, bias) vs IDEAL, mass-weighted over live nodes."""
    live, obs, tru, wt = node_view(s)
    if P is None or live.sum() < 5:
        return np.nan, np.nan
    o, t, ww = obs[live], tru[live], wt[live]
    ic = ideal_curve(o, t, ww, EDGES)
    ideal = ic[np.clip(np.digitize(o, EDGES) - 1, 0, len(ic) - 1)]
    ok = np.isfinite(ideal) & (ww > 0)
    if ok.sum() < 5:
        return np.nan, np.nan
    mu = mu_fn(P, o, live)
    gap = ideal - o
    den = float(np.average(np.abs(gap)[ok], weights=ww[ok]))
    g = float(np.average((np.sign(gap) * (mu - o))[ok], weights=ww[ok])) / den if den > 1e-9 else np.nan
    return g, float(np.average((mu - ideal)[ok], weights=ww[ok]))


def main():
    projections = [
        ("symmetric h=0.15", lambda P, d, lv: L.project(P, d, 0.15, "mean")),
        ("asymmetric (3 consts)", lambda P, d, lv: project_asym(P, d)),
        ("POSTERIOR (sigma=belief)", None),   # needs the scenario; filled per-scenario below
    ]
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        print(f"=== {su}: projection vs IDEAL, on the W2 kNN landscape ===")
        print(f"{'projection':26s} {'mean gain':>10s} {'median':>8s} {'mean bias':>10s} {'negative on':>12s}")
        for name, fn in projections:
            gs, bs = [], []
            for s in scen:
                P = L.recipe(s, knn_scale=KNN)
                f = fn if fn is not None else (
                    lambda Pp, d, lv, _s=s: project_posterior(Pp, d, sigma_obs(_s)[lv]))
                g, b = gain_bias(s, P, f)
                if np.isfinite(g):
                    gs.append(g)
                    bs.append(b)
            a = np.array(gs)
            print(f"{name:26s} {a.mean():+10.3f} {np.median(a):+8.3f} {np.mean(bs):+10.3f} "
                  f"{int((a < 0).sum()):>7d}/{len(a)}")
        print()

    print("=== THE SIGMA-CALIBRATION: scale the node's declared width by c ===")
    print(f"{'suite':6s} {'stratum':16s} " + " ".join(f"{f'c={c}':>8s}" for c in C_GRID) + f" {'c*':>6s}")
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        for stratum, keep in (("all", lambda s: True),
                              ("unstranded", lambda s: s["group"][2] == "0.50"),
                              ("stranded", lambda s: s["group"][2] == "0.99"),
                              ("capture ON", lambda s: s["group"][0] == "ON")):
            means = []
            for c in C_GRID:
                gs = []
                for s in scen:
                    if not keep(s):
                        continue
                    P = L.recipe(s, knn_scale=KNN)
                    g, _ = gain_bias(s, P, lambda Pp, d, lv, _s=s, _c=c:
                                     project_posterior(Pp, d, sigma_obs(_s, _c)[lv]))
                    if np.isfinite(g):
                        gs.append(g)
                means.append(float(np.mean(gs)) if gs else np.nan)
            best = C_GRID[int(np.nanargmax(means))]
            print(f"{su:6s} {stratum:16s} " + " ".join(f"{m:8.3f}" for m in means) + f" {best:6.2f}")
    print("\n  c* ~ 1 => the declared widths are honest in the units the prior consumes.")
    print("  c* >> 1 => pass-0 is over-confident by that factor on this axis.")


def plot(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    scen = {s["cond"]: s for s in L.load_scenarios("ambig")}
    picks = ["gdna_gdna300_ss_0.50_nrna_none_capture_on", "gdna_gdna100_ss_0.50_nrna_none_capture_on",
             "gdna_none_ss_0.50_nrna_none_capture_on", "gdna_gdna300_ss_0.99_nrna_none_capture_on",
             "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
             "gdna_gdna300_ss_0.50_nrna_none_capture_off"]
    fig, axes = plt.subplots(2, 3, figsize=(20, 9))
    for ax, cond in zip(axes.ravel(), picks):
        s = scen[cond]
        P = L.recipe(s, knn_scale=KNN)
        live, obs, tru, wt = node_view(s)
        o, t, ww = obs[live], tru[live], wt[live]
        ax.scatter(o, t, s=np.clip(ww / max(ww.max(), 1) * 60, 1, 60), c="0.75", lw=0, alpha=.5)
        ax.plot([-5, 2], [-5, 2], "k--", lw=0.9, label="IDENTITY (inert)")
        ax.plot(MID, ideal_curve(o, t, ww, EDGES), "C2", lw=2.4, label="IDEAL E[truth|obs]")
        ax.plot(MID, L.project(P, MID, 0.15, "mean"), "C7", lw=1.3, label="symmetric h=0.15")
        ax.plot(MID, project_asym(P, MID), "C3", lw=1.3, label="asymmetric")
        # the posterior transfer at the median declared width for this scenario
        sm = float(np.median(sigma_obs(s)[live]))
        ax.plot(MID, project_posterior(P, MID, np.full_like(MID, sm)), "C0", lw=2.0,
                label=f"POSTERIOR (sigma={sm:.2f})")
        ax.set_title(f"{' · '.join(s['group'])}\n{cond}", fontsize=8)
        ax.set_xlim(-5, 2)
        ax.set_ylim(-5, 2)
        ax.tick_params(labelsize=7)
    axes.ravel()[0].legend(fontsize=7, loc="upper left")
    fig.suptitle("W3 — projection transfer functions on the W2 kNN landscape. "
                 "A good projection lies ON the green IDEAL curve.", fontsize=13)
    fig.supxlabel("observed log10 rho_g (pass-0)")
    fig.supylabel("projected mu* / true log10 rho_g")
    fig.tight_layout(rect=(0.01, 0.02, 1, 0.95))
    fig.savefig(out, dpi=110)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
    if "--plot" in sys.argv:
        plot("/Users/mkiyer/proj/rigel/docs/calibration/figures/gdna_w3_projection.png")
