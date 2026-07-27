"""Shared library for GOAL-1 landscape exploration (workflow agents import this). Loads the precomputed
substrate cache (gdna_cache_build.py) and provides fit methods, the DUAL metric (overall EMD + enriched-region
L1), mask helpers, and evaluate(). NO calibrate re-runs — pure numpy on the cached per-node arrays.

A "recipe" is a function  recipe(scen) -> (fitted density on GRID)  using the helpers below. evaluate(recipe)
runs it across the battery and returns per-group + aggregate metrics.

Fit methods (all additive / unit-mass, no EM competition):
  fit_kde(g,E,h)                 point-estimate Gaussian at log10(max(g,1)/E), fixed h. Floors at 1/E.
  fit_poisson(g,E,h_smooth)      additive Poisson P(g|ρE) (ZERO-NATIVE). Count sets the width.
  fit_precwidth(g,E,tau,h_floor) additive Gaussian at log10(g/E) with per-node width sqrt(tau^2+h_floor^2):
                                 PRECISION-AS-WIDTH — imprecise (high var) => broad => self-down-weighted,
                                 NO cutoff/percentile. tau = sqrt(Var(log f_g)) in decades.
  fit_poisln(g,E,tau)            additive Poisson-LOGNORMAL: Poisson count width ⊗ belief width tau. Both.
"""
from __future__ import annotations

import os
import pickle
from pathlib import Path

import numpy as np
from scipy.special import gammaln
from scipy.stats import wasserstein_distance

_EPS = 1e-12
LN10 = np.log(10.0)
GRID = np.linspace(-5.0, 2.5, 260)  # log10 rho_g
_lnat = GRID * LN10
# The caches live in the repo's untracked `scratchpad/` — NOT a session scratchpad. The original location was
# a per-session /private/tmp dir, which is volatile: the 2026-07-21 caches were lost when that session ended.
_CACHE = Path("/Users/mkiyer/proj/rigel/scratchpad/gdna_substrate_cache.pkl")


_SCR = _CACHE.parent
CACHES = {"ambig": _CACHE, "quick": _SCR / "gdna_substrate_cache_quick.pkl"}


def load_scenarios(which=None):
    """Load a substrate cache. `which` = a key in CACHES ('ambig'/'quick'), a path, or None (env
    RIGEL_GDNA_CACHE or the default ambig cache)."""
    if which in CACHES:
        p = CACHES[which]
    elif which:
        p = Path(which)
    else:
        p = Path(os.environ.get("RIGEL_GDNA_CACHE", str(_CACHE)))
    return pickle.loads(Path(p).read_bytes())


def _norm(d):
    return d / max(float(np.sum(d)), _EPS)


# ---------------- fit methods (g, E are per-node arrays of the SELECTED training nodes) ----------------
def fit_kde(g, E, h=0.15):
    if not len(g):
        return None
    a = np.log10(np.maximum(g, 1.0)) - np.log10(E)
    z = (GRID[:, None] - a[None, :]) / h
    return _norm((np.exp(-0.5 * z * z) / (h * np.sqrt(2 * np.pi))).sum(1))


def fit_poisson(g, E, h_smooth=0.0):
    if not len(g):
        return None
    lam = np.exp(_lnat)[None, :] * E[:, None]
    ll = g[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    pn /= np.maximum(pn.sum(1, keepdims=True), _EPS)
    d = pn.sum(0)
    if h_smooth > 0:
        k = np.exp(-0.5 * ((GRID[:, None] - GRID[None, :]) / h_smooth) ** 2)
        d = (k / k.sum(0, keepdims=True)) @ d
    return _norm(d)


def fit_precwidth(g, E, tau, h_floor=0.05):
    if not len(g):
        return None
    a = np.log10(np.maximum(g, 1.0)) - np.log10(E)
    hh = np.sqrt(np.maximum(tau, 0.0) ** 2 + h_floor**2)  # per-node width (decades)
    z = (GRID[:, None] - a[None, :]) / hh[None, :]
    return _norm((np.exp(-0.5 * z * z) / (hh[None, :] * np.sqrt(2 * np.pi))).sum(1))


def fit_poisln(g, E, tau, n_gh=7):
    """Additive Poisson-lognormal: marginalize log g ~ N(log ĝ, (tau·ln10)^2) against Poisson(g|ρE)."""
    if not len(g):
        return None
    tnat = np.maximum(tau, 0.0) * LN10  # tau is in decades; the Poisson count is in nat-log
    pure = float(np.max(tnat)) <= _EPS
    x, wq = (np.zeros(1), np.array([np.sqrt(np.pi)])) if pure else np.polynomial.hermite.hermgauss(n_gh)
    lwq = np.log(wq) - 0.5 * np.log(np.pi)
    gq = g[:, None] * np.exp(np.sqrt(2.0) * tnat[:, None] * x[None, :])  # (n,Q)
    lam = np.exp(_lnat)[None, None, :] * E[:, None, None]  # (n,1,G)
    lp = gq[:, :, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(gq[:, :, None] + 1.0)  # (n,Q,G)
    mx = lp.max((1, 2), keepdims=True)
    cell = np.log(np.sum(np.exp(lwq[None, :, None] + lp - mx), 1)) + mx[:, 0, :]  # (n,G)
    cell -= cell.max(1, keepdims=True)
    pn = np.exp(cell)
    pn /= np.maximum(pn.sum(1, keepdims=True), _EPS)
    return _norm(pn.sum(0))


# ---------------- mask helpers on a scenario dict ----------------
def masks(s):
    isr = s["is_region"]
    fp, fn = s["fp"], s["fn"]
    haveE = s["eff"] > 1e-9
    havem = s["mass"] > _EPS
    return dict(
        region=isr, boundary=~isr, ntype=s["ntype"], haveE=haveE, havemass=havem,
        single=fp ^ fn, gonly=~fp & ~fn, ambig=fp & fn,
        struct_zero=haveE & isr & ((s["ntype"] == 0) | (s["ntype"] == 1)) & ~havem,
        base=haveE & isr & havem & ~(fp & fn),  # region, has mass, not ambiguous
    )


def vpercentile(s, sel, p):
    """var threshold at the p-th percentile OF THE SELECTED nodes (a percentile subset)."""
    v = s["var"][sel]
    thr = np.percentile(v, p) if v.size else np.inf
    return sel & (s["var"] <= thr)


# ---------------- metrics ----------------
def oracle_landscape(s):
    m = (s["eff"] > 1e-9) & s["is_region"]
    return fit_poisson(s["G"][m], s["eff"][m])


def emd(fit, oracle):
    if fit is None or oracle is None:
        return np.nan
    return float(wasserstein_distance(GRID, GRID, fit, oracle))


def enriched_l1(fit, oracle):
    """L1 mass difference ABOVE the oracle's primary (depleted) mode — the enriched-region fidelity. Catches a
    fit that scores a good overall EMD while dropping a small high-density enriched mode."""
    if fit is None or oracle is None:
        return np.nan
    hi = GRID > GRID[int(np.argmax(oracle))]
    return float(np.sum(np.abs(fit[hi] - oracle[hi])))


def evaluate(recipe, scen=None, verbose=False):
    """Run recipe(scen)->density across the battery. Returns dict with mean/max EMD and enriched_l1, per-group
    breakdown, and (mean EMD, mean enr_l1). recipe returns a density on GRID (or None)."""
    scen = scen if scen is not None else load_scenarios()
    from collections import defaultdict
    g_emd, g_enr = defaultdict(list), defaultdict(list)
    for s in scen:
        orc = oracle_landscape(s)
        fit = recipe(s)
        g_emd[s["group"]].append(emd(fit, orc))
        g_enr[s["group"]].append(enriched_l1(fit, orc))
    groups = sorted(g_emd)
    gm_emd = {g: float(np.nanmean(g_emd[g])) for g in groups}
    gm_enr = {g: float(np.nanmean(g_enr[g])) for g in groups}
    out = dict(
        mean_emd=float(np.nanmean(list(gm_emd.values()))),
        max_emd=float(np.nanmax(list(gm_emd.values()))),
        mean_enr=float(np.nanmean(list(gm_enr.values()))),
        max_enr=float(np.nanmax(list(gm_enr.values()))),
        per_group_emd=gm_emd, per_group_enr=gm_enr,
    )
    if verbose:
        print(f"mean_emd={out['mean_emd']:.3f} max_emd={out['max_emd']:.3f}  "
              f"mean_enr={out['mean_enr']:.3f} max_enr={out['max_enr']:.3f}")
    return out


def project(P, d_obs, h_proj=0.15, readout="mean"):
    """THE projection (the endpoint): observed DNA log10-density d_obs -> anchor mu* onto landscape P via the
    sampling likelihood r_j ∝ P(mu_j)·N(d_obs; mu_j, h_proj). readout='mean' (responsibility mean, dilutes
    toward the tall depleted mode) | 'mode' (argmax — snaps to the nearest mode) | 'median'."""
    d_obs = np.asarray(d_obs, float)
    logP = np.log(np.maximum(np.asarray(P, float), 1e-30))
    z = (d_obs[:, None] - GRID[None, :]) / h_proj
    lr = logP[None, :] - 0.5 * z * z
    lr -= lr.max(1, keepdims=True)
    r = np.exp(lr)
    r /= np.maximum(r.sum(1, keepdims=True), _EPS)
    if readout == "mode":
        return GRID[np.argmax(r, axis=1)]
    if readout == "median":
        cw = np.cumsum(r, axis=1)
        return GRID[np.clip((cw < 0.5).sum(1), 0, len(GRID) - 1)]
    return (r * GRID[None, :]).sum(1)


def enriched_sensitivity(recipe, scen=None, h_proj=0.15, readout="mean", t_enr=0.0):
    """THE sensitivity metric. For each scenario: fit landscape=recipe(s), project every node's observed DNA
    density, and measure — over the TRULY-enriched nodes (oracle log10-density > t_enr) — how close mu* lands to
    the oracle. Returns per-scenario means aggregated: enr_abs_err = mean|mu*-oracle| (lower=better), enr_recovery
    = mean(mu*-observed) (want POSITIVE = pulled up toward the enriched mode), n_scen with an enriched population.
    Also a 'fabrication' canary: on zero-DNA/capOFF scenarios, how high mu* drifts (want ~0)."""
    scen = scen if scen is not None else load_scenarios()
    errs, recov, fab = [], [], []
    for s in scen:
        P = recipe(s)
        if P is None:
            continue
        g, E, G = s["g_hat"], s["eff"], s["G"]
        live = E > 1e-9
        obs = np.log10(np.maximum(g, 1e-9) / np.maximum(E, 1e-9))
        tru = np.log10(np.maximum(G, 1e-9) / np.maximum(E, 1e-9))
        enr = live & (G > 0) & (tru > t_enr)
        if enr.sum() >= 5:
            mu = project(P, obs[enr], h_proj, readout)
            errs.append(float(np.mean(np.abs(mu - tru[enr]))))
            recov.append(float(np.mean(mu - obs[enr])))
        if s["group"][1] == "none":  # zero-DNA: mu* over expressed nodes should stay low (fabrication canary)
            m = live & (s["mass"] > _EPS) & (s["ntype"] == 2)
            if m.sum() >= 5:
                fab.append(float(np.mean(project(P, obs[m], h_proj, readout))))
    return dict(
        enr_abs_err=float(np.mean(errs)) if errs else np.nan,
        enr_recovery=float(np.mean(recov)) if recov else np.nan,
        fabrication=float(np.mean(fab)) if fab else np.nan,
        n_enr_scen=len(errs),
    )


def enriched_sensitivity_suites(recipe, suites=("ambig", "quick"), h_proj=0.15, readout="mean", t_enr=0.0):
    """Cross-suite enriched-sensitivity. Returns per-suite + aggregate enr_abs_err (lower=better), enr_recovery
    (want POSITIVE — the mean upward pull of enriched nodes toward truth), fabrication (want ~low)."""
    per = {}
    for su in suites:
        try:
            per[su] = enriched_sensitivity(recipe, load_scenarios(su), h_proj, readout, t_enr)
        except FileNotFoundError:
            continue
    if not per:
        return None
    import numpy as _np
    return dict(
        per_suite=per,
        enr_abs_err=float(_np.nanmean([p["enr_abs_err"] for p in per.values()])),
        enr_recovery=float(_np.nanmean([p["enr_recovery"] for p in per.values()])),
        fabrication=float(_np.nanmean([p["fabrication"] for p in per.values()])),
    )


def evaluate_suites(recipe, suites=("ambig", "quick")):
    """Evaluate a recipe across MULTIPLE suites (generalization). Returns per-suite metrics + the cross-suite
    aggregate (mean over suites, and the WORST max_emd over suites — the robustness-across-datatypes number)."""
    per = {}
    for su in suites:
        try:
            per[su] = evaluate(recipe, load_scenarios(su))
        except FileNotFoundError:
            continue
    if not per:
        return None
    return dict(
        per_suite=per,
        mean_emd=float(np.mean([p["mean_emd"] for p in per.values()])),
        worst_max_emd=float(np.max([p["max_emd"] for p in per.values()])),
        mean_enr=float(np.mean([p["mean_enr"] for p in per.values()])),
        worst_max_enr=float(np.max([p["max_enr"] for p in per.values()])),
    )


if __name__ == "__main__":
    # verify: reproduce the v25 winner and show the precision-width no-cutoff candidate
    scen = load_scenarios()

    def r_v25(s):
        mk = masks(s)
        sel = vpercentile(s, mk["base"] & (mk["single"] | mk["gonly"]), 25) | mk["struct_zero"]
        return fit_poisson(s["g_hat"][sel], s["eff"][sel])

    def r_precwidth(s):  # NO cutoff: all region non-ambig + struct-zero, precision-as-width
        mk = masks(s)
        sel = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
        tau = np.sqrt(np.maximum(s["var"], 0.0)) / LN10  # decades
        return fit_precwidth(s["g_hat"][sel], s["eff"][sel], tau[sel])

    print("v25 winner        :", end=" ")
    evaluate(r_v25, scen, verbose=True)
    print("precwidth no-cutoff:", end=" ")
    evaluate(r_precwidth, scen, verbose=True)
