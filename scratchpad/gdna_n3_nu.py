"""N3 (W3b) — FIT the Student-t's nu, do not pick it.

W3b established that the right defence against pass-0's confidently-wrong TAIL is a heavier-tailed
likelihood, not a wider one (Gaussian c=5 harms 9/32 conditions; t(nu=2) gets 70 % of the gain with 4/32),
and set the task: `nu` follows from the residual's kurtosis by method of moments — the allowed pattern
(kappa, both strand overdispersions, omega_graft).

The residual is  z = (observed - true) / sigma_declared,  in decades, over the nodes psi applies the prior
to.  For a Student-t, excess kurtosis = 6/(nu-4)  =>  nu = 4 + 6/kurt_excess.

⚠ That estimator EXISTS ONLY FOR nu > 4.  This measures whether it is usable here at all, and if not, what
the defensible replacement is.  Reported alongside are two moment estimators that stay finite for every
nu > 0, so the answer does not depend on a moment the data may not have:

  * the LOG-VARIANCE estimator.  For z ~ t_nu, log(z^2) has variance psi'(1/2) + psi'(nu/2) exactly, for
    every nu > 0 (a t is a Gaussian over a chi-square scale, and log of that is a sum of two independent
    log-chi-squares).  Var(log z^2) is finite whatever the tail, so this inverts where kurtosis cannot.
  * the profile MLE on nu, as the arbiter both moment estimators are checked against.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_n3_nu.py
"""
from __future__ import annotations

import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import gammaln, polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import gdna_explore_lib as L  # noqa: E402
from gdna_projection_review import node_view  # noqa: E402

LN10 = np.log(10.0)


def residual_z(scen, keep=lambda s: True):
    """z = (observed - true)/sigma_declared over the live nodes of the kept conditions, with their mass."""
    zs, ws = [], []
    for s in scen:
        if not keep(s):
            continue
        live, obs, tru, wt = node_view(s)
        sig = np.sqrt(np.maximum(np.nan_to_num(s["var"], nan=np.inf, posinf=np.inf), 0.0)) / LN10
        ok = live & np.isfinite(sig) & (sig > L.GRID_H)
        zs.append(((obs - tru) / sig)[ok])
        ws.append(wt[ok])
    return np.concatenate(zs), np.concatenate(ws)


def nu_from_kurtosis(z):
    """MoM on the 4th moment: nu = 4 + 6/excess-kurtosis. Undefined for nu <= 4."""
    k = float(np.mean(z**4) / np.mean(z**2) ** 2 - 3.0)
    return k, (4.0 + 6.0 / k if k > 0 else np.nan)


def nu_from_logvar(z):
    """MoM on Var(log z^2) = psi'(1/2) + psi'(nu/2) — finite for EVERY nu > 0, so it inverts where the
    kurtosis estimator cannot exist."""
    lv = float(np.var(np.log(np.maximum(z**2, 1e-300))))
    target = lv - float(polygamma(1, 0.5))
    if target <= 0:
        return lv, np.inf
    try:
        return lv, 2.0 * brentq(lambda v: polygamma(1, v / 2.0) - target, 1e-3, 1e6)
    except ValueError:
        return lv, np.nan


def nu_mle(z):
    """Profile MLE of nu for z ~ t_nu with the declared scale held at 1 (the projection's own assumption)."""
    def nll(v):
        return -float(np.sum(gammaln((v + 1) / 2) - gammaln(v / 2) - 0.5 * np.log(np.pi * v)
                             - (v + 1) / 2 * np.log1p(z**2 / v)))
    grid = np.concatenate([np.linspace(0.2, 8, 40), np.linspace(8.5, 60, 30)])
    return float(grid[int(np.argmin([nll(v) for v in grid]))])


def nu_heldout(scen, keep=lambda s: True, grid=(0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 8, 12, 20, 40)):
    """The ORACLE-FREE estimator, and the only one that can run on cfRNA: fit the landscape on half the
    training nodes and score the OTHER half's observed density under  P (x) t_nu(sigma_i)  — the same
    held-out predictive protocol W2 used for the bandwidth. No truth anywhere.

    Halves are taken by alternating rank on the observed density, so both halves span the same support and
    the split introduces no randomness."""
    tot = np.zeros(len(grid))
    for s in scen:
        if not keep(s):
            continue
        mk = L.masks(s)
        sel = L.recipe_substrate(s, mk)
        if sel.sum() < 40:
            continue
        idx = np.flatnonzero(sel)
        d_all = np.log10(np.maximum(s["g_hat"][idx], 1.0) / np.maximum(s["eff"][idx], 1e-9))
        rank = np.argsort(np.argsort(d_all))
        for a, b in ((rank % 2 == 0, rank % 2 == 1), (rank % 2 == 1, rank % 2 == 0)):
            fit_sel = np.zeros_like(sel)
            fit_sel[idx[a]] = True
            P = L.recipe(s, sel=fit_sel, w=L.recipe_weights(s, fit_sel, mk), knn_scale=0.5)
            if P is None:
                continue
            j = idx[b]
            sig = np.sqrt(np.maximum(np.nan_to_num(s["var"][j], nan=np.inf, posinf=np.inf), 0.0)) / LN10
            ok = np.isfinite(sig) & (sig > L.GRID_H)
            if ok.sum() < 10:
                continue
            d = np.log10(np.maximum(s["g_hat"][j][ok], 1.0) / np.maximum(s["eff"][j][ok], 1e-9))
            u = (d[:, None] - L.GRID[None, :]) / sig[ok][:, None]
            for k, v in enumerate(grid):
                lt = (gammaln((v + 1) / 2) - gammaln(v / 2) - 0.5 * np.log(np.pi * v)
                      - (v + 1) / 2 * np.log1p(u**2 / v) - np.log(sig[ok])[:, None])
                m = lt.max(1, keepdims=True)
                tot[k] += float(np.sum(np.log(np.sum(np.maximum(P, 1e-30)[None, :] * np.exp(lt - m), 1))
                                       + m[:, 0]))
    return float(grid[int(np.argmax(tot))]), tot


def main():
    print("=== N3: fitting the Student-t nu from the residual z = (obs - true)/sigma_declared ===\n")
    print(f"{'suite':6s} {'stratum':16s} {'n':>8s} {'z==0':>6s} {'sd(z)':>7s} {'MAD-sd':>7s} "
          f"{'kurt':>8s} {'nu(kurt)':>9s} {'nu(logvar)':>11s} {'nu(MLE)':>8s} {'nu(HELD-OUT)':>13s}")
    for su in ("ambig", "quick"):
        scen = L.load_scenarios(su)
        for name, keep in (("all", lambda s: True),
                           ("unstranded", lambda s: s["group"][2] == "0.50"),
                           ("stranded", lambda s: s["group"][2] == "0.99"),
                           ("capture ON", lambda s: s["group"][0] in ("ON", "VSTRONG")),
                           ("capture OFF", lambda s: s["group"][0] == "OFF")):
            z, _ = residual_z(scen, keep)
            if z.size < 100:
                continue
            k, nk = nu_from_kurtosis(z)
            # the exact ties (obs == true, both at the one-count wall) send log z^2 to -inf and destroy the
            # log-variance moment; they are real successes, not outliers, so the estimator is reported on
            # z != 0 and flagged rather than silently repaired.
            _, nl = nu_from_logvar(z[z != 0])
            mad = float(np.median(np.abs(z - np.median(z))) * 1.4826)
            ho, _ = nu_heldout(scen, keep)
            print(f"{su:6s} {name:16s} {z.size:8d} {100 * (z == 0).mean():5.1f}% {np.std(z):7.2f} "
                  f"{mad:7.2f} {k:8.1f} {nk:9.2f} {nl:11.2f} {nu_mle(z):8.2f} {ho:13.1f}")
    print("\n  nu(kurt)     = 4 + 6/excess-kurtosis. ⛔ EXISTS ONLY FOR nu > 4, so it can never return the")
    print("                 empirical optimum nu ~ 2 whatever the data say — it saturates at its boundary.")
    print("  nu(logvar)   = from Var(log z^2) = psi'(1/2) + psi'(nu/2). Finite for every nu > 0, but ~24 % of")
    print("                 z are EXACTLY 0 (obs and true both at the one-count wall) and log z^2 diverges.")
    print("  nu(MLE)      = the arbiter — needs the ORACLE, so it cannot run on real data.")
    print("  nu(HELD-OUT) = fit the landscape on half the nodes, score the other half's observed density")
    print("                 under P (x) t_nu(sigma_i). NO TRUTH — the only estimator that runs on cfRNA.")


if __name__ == "__main__":
    main()
