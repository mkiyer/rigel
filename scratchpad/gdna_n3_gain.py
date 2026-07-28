"""N3 — does the FITTED nu deliver W3b's Student-t gain? Closes the loop.

W3b measured the projection gain for hand-picked nu (2 / 4 / 8) and said "fit it, do not pick it".
`gdna_n3_nu.py` fits it two ways: the oracle MLE says nu ~ 2.0-2.2; the ORACLE-FREE held-out predictive —
the only estimator that can run on cfRNA — says nu ~ 3.0, one step light-tailed, which is the same
direction of bias W2 already documented for held-out likelihood.

This scores the projection at both, so the cost of using the deployable estimator is a number rather than
an assumption. Same instrument as `gdna_w3_projection.py`: gain vs IDEAL(d) = E[true | observed = d].

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n3_gain.py
"""
from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_w3_projection import KNN, gain_bias, project_posterior, sigma_obs  # noqa: E402


def project_t(P, d, sig, nu):
    """The same posterior projection with a Student-t likelihood: a Gaussian marginalised over an UNCERTAIN
    width, which is the true state of knowledge about `sigma_declared` (W3b's quintile measurement: honest
    in the median, heavy in the tail). No new constant — `nu` is fitted, `sigma` is the node's own."""
    logP = np.log(np.maximum(P, 1e-30))
    sd = np.maximum(np.asarray(sig, float), L.GRID_H)[:, None]
    u = (np.asarray(d, float)[:, None] - L.GRID[None, :]) / sd
    lr = logP[None, :] - 0.5 * (nu + 1.0) * np.log1p(u * u / nu)
    lr -= lr.max(1, keepdims=True)
    r = np.exp(lr)
    r /= np.maximum(r.sum(1, keepdims=True), 1e-30)
    return (r * L.GRID[None, :]).sum(1)


def main():
    arms = [("Gaussian c=1 (W3 adopted)", None, 1.0), ("Gaussian c=5 (rejected)", None, 5.0)]
    arms += [(f"Student-t nu={v}", v, 1.0) for v in (2.0, 2.2, 3.0, 4.0, 8.0)]
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        print(f"=== {su}: projection gain on the W2 kNN landscape (n={len(scen)}) ===")
        print(f"{'arm':28s} {'mean gain':>10s} {'median':>8s} {'mean bias':>10s} {'made WORSE':>11s}")
        for name, nu, c in arms:
            gs, bs = [], []
            for s in scen:
                P = L.recipe(s, knn_scale=KNN)

                def fn(Pp, d, lv, _s=s, _nu=nu, _c=c):
                    sig = sigma_obs(_s, _c)[lv]
                    return (project_posterior(Pp, d, sig) if _nu is None
                            else project_t(Pp, d, sig, _nu))
                g, b = gain_bias(s, P, fn)
                if np.isfinite(g):
                    gs.append(g)
                    bs.append(b)
            a = np.array(gs)
            print(f"{name:28s} {a.mean():+10.3f} {np.median(a):+8.3f} {np.mean(bs):+10.3f} "
                  f"{int((a < 0).sum()):>6d}/{len(a)}")
        print()
    print("  nu = 2.2 is the oracle MLE; nu = 3.0 is the oracle-free held-out fit (the deployable one).")


if __name__ == "__main__":
    main()
