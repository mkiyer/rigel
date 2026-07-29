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
#: the log10-rho grid's own resolution (decades) — the smallest width this axis can represent, and therefore
#: the natural "minimal smoothing" for a point-mass reference. Derived from the grid, not asserted.
GRID_H = float(GRID[1] - GRID[0])


def live_region(s):
    """The default oracle substrate: live REGION nodes (what `oracle_landscape` has always used)."""
    return (s["eff"] > 1e-9) & s["is_region"]


def oracle_landscape(s, sel=None, knn_scale=0.0):
    """ORACLE density over `sel` (default: live region nodes), rendered with the SAME zero-native additive
    Poisson smoother the fits use — so a fit-vs-oracle comparison isolates the input COUNTS and weights
    rather than the smoother.

    ⚠⚠ **`knn_scale=0` MAKES THIS REFERENCE A COMB ABOVE log10 ρ ≈ 0, AND IT BIASED EVERY BANDWIDTH
    DECISION TAKEN AGAINST IT** (measured 2026-07-27, owner-spotted). The Poisson kernel's width on the log
    axis is `1/(√G·ln10)`, so it **shrinks as ρ^(−1/2)**: across this suite it collapses 41× (0.32 → 0.008
    decades) and drops below `GRID_H`, at which point a node cannot be rendered as anything but a delta in
    one cell. Above +1 decade **69 % of oracle nodes are such deltas**, and the reference's roughness (total
    variation of log P per decade) is **40.7** against the ~2–4 of a smooth unimodal bump.

    That is a *measurement* kernel applied to numbers that are the **TRUTH** — `G` is not an observation of
    `ρ`, it IS `ρ·E`, so there is no posterior width to render. What a reference needs instead is a
    **POPULATION** resolution: how finely a sample of this many nodes can resolve a density.

    **Pass `knn_scale` to get it** — the identical nearest-neighbour-spacing rule the fit uses
    (:func:`recipe`), so the two are rendered by one resolution law and the comparison still isolates the
    counts and the weights. The zero-count nodes keep their zero-native Poisson decay either way, which is
    the honest "ρ is anything below the resolution wall" (see :func:`oracle_pointmass`, and W1a).

    ⚠ Passing `sel` matters — the default is REGION-ONLY, so a fit trained partly on BOUNDARY nodes is
    otherwise scored against a mismatched population."""
    m = live_region(s) if sel is None else sel
    if knn_scale <= 0:
        return fit_poisson(s["G"][m], s["eff"][m])
    return _render(s["G"][m], s["eff"][m], np.ones(int(np.sum(m))), knn_scale)


def oracle_pointmass(s, sel=None):
    """The TRUE per-node density distribution at MINIMAL smoothing — the honest reference for WIDTH.

    Each node contributes unit mass at its true `log10(max(G,1)/E)` (the one-count resolution wall — a
    zero-gDNA node's honest statement is "gDNA <= 1 count here"), smoothed only by `GRID_H`, the axis's own
    resolution. Where `oracle_landscape` is the population CONVOLVED with the Poisson measurement width,
    this is the population itself: the difference between the two IS the smear."""
    m = live_region(s) if sel is None else sel
    E = np.maximum(s["eff"][m], _EPS)
    a = np.log10(np.maximum(s["G"][m], 1.0)) - np.log10(E)
    if not a.size:
        return None
    z = (GRID[:, None] - a[None, :]) / GRID_H
    return _norm(np.exp(-0.5 * z * z).sum(1))


def spread(P):
    """sd of a density on GRID, in decades — the effective width of a rendered landscape."""
    if P is None:
        return np.nan
    m = float((P * GRID).sum())
    return float(np.sqrt(max((P * (GRID - m) ** 2).sum(), 0.0)))


def modes(P, prominence=0.02):
    """Local maxima of a density on GRID as (location, height), prominence-filtered relative to its max.

    ⚠⚠ **NOT a structural metric — do not score landscapes by counting these.** The ground truth of this
    simulator is TWO components (gDNA is uniform; capture partitions nodes into captured / not-captured),
    and at `prominence=0.02` this counts noise wiggles: on one enriched-only landscape it reports 7 modes at
    0.02, 3 at 0.10 and **1 at 0.35**. Mode-recall / spurious-rate columns produced before 2026-07-27 are
    invalid for that reason. Use :func:`two_component` instead; keep this for eyeballing curves only."""
    if P is None:
        return []
    from scipy.signal import find_peaks
    pk, _ = find_peaks(P, prominence=prominence * float(np.max(P)))
    return [(float(GRID[i]), float(P[i])) for i in pk]


def two_component(P, split=None):
    """Score a landscape the shape the GROUND TRUTH actually has: a depleted component and an enriched one.

    gDNA is uniform, so the depleted level is a near point mass (measured IQR 0.02–0.04 decades on
    intergenic/intron); capture lifts the covered nodes into one broad enriched mode ~2.7 decades above it.
    Returns the location and width of each component and the enriched mass fraction — all robust to the
    wiggles that defeat mode counting.

    `split` defaults to the deepest valley of the density between its two extremes."""
    if P is None:
        return None
    c = np.cumsum(P)
    if split is None:
        lo_i, hi_i = int(np.argmax(P)), int(np.searchsorted(c, 0.995))
        seg = P[min(lo_i, hi_i):max(lo_i, hi_i) + 1]
        split = float(GRID[min(lo_i, hi_i) + int(np.argmin(seg))]) if seg.size > 2 else 0.0
    dep, enr = GRID <= split, GRID > split

    def _loc_wid(m):
        w = P[m]
        t = float(w.sum())
        if t <= 1e-12:
            return np.nan, np.nan
        g = GRID[m]
        mu = float((w * g).sum() / t)
        return mu, float(np.sqrt(max((w * (g - mu) ** 2).sum() / t, 0.0)))

    d_loc, d_wid = _loc_wid(dep)
    e_loc, e_wid = _loc_wid(enr)
    return dict(split=split, dep_loc=d_loc, dep_width=d_wid, enr_loc=e_loc, enr_width=e_wid,
                enr_mass=float(P[enr].sum()))


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


# ---------------- THE LANDSCAPE RECIPE (the Role-B candidate; was scripts/debug/gdna_landscape_recipe.py) --
# Promoted 2026-07-27 from a bare `def recipe(s)` that four call sites loaded with `exec(open(...).read())`.
# It is the seed of the production `GdnaLandscape` (docs/CARRY_FORWARD.md), so
# it needs to be importable, lintable and diffable. Behaviour is byte-identical to the exec'd version.

#: kernel/grid resolution floor as a log-rate variance — the library's own `fit_kde` bandwidth (0.15 decades).
#: A fixed reference SCALE on the rho axis, NOT a per-sample percentile.
_S0 = (0.15 * LN10) ** 2


def recipe_substrate(s, mk):
    """The training set, from the taxonomy (production plan §2.2) — no thresholds:

    * **circularity** → AMBIG is EXCLUDED (structural): it is the two-root ambiguity the prior exists to
      resolve. Measured 2026-07-27: admitting it wins only against an oracle that itself contains AMBIG;
      against the non-AMBIG population the prior is applied to it LOSES on both suites.
    * **identifiability** → the zero-count structural anchor is INCLUDED: "gDNA is absent here" is the
      strongest depletion evidence there is, and dropping it costs +0.32 / +0.60 EMD.
    * **precision** → handled ENTIRELY by the continuous weight (:func:`recipe_weights`), never by
      admission. A hard precision cutoff scores WORSE than ignoring precision altogether (+0.25 / +0.43 vs
      +0.14 / +0.14).
    * **geometry** → boundaries are INCLUDED (owner, 2026-07-27); include-vs-exclude is a standing ablation.

    ⛔ **The `tau_prior` boundary gate was REMOVED here on 2026-07-27** — three bare constants
    (0.15 / 0.70 / 0.3) and a cliff, clipping to its lower bound on 19/32 conditions, discarding boundaries
    that carried 12–16 % of boundary weight. Removing it is FREE: EMD −0.0034 (ambig) / −0.0255 (quick),
    mode recall and spurious rate unchanged. Do not reintroduce a single-constant admission cutoff."""
    region = (mk["base"] & (mk["single"] | mk["gonly"])) | mk["struct_zero"]
    boundary = mk["boundary"] & (s["eff"] > 1e-9) & (s["mass"] > _EPS) & (mk["single"] | mk["gonly"])
    return region | boundary


def recipe_weights(s, sel, mk):
    """Per-node RELIABILITY as MASS. `ref` is the sum of the two IRREDUCIBLE variance sources of the log
    gDNA rate — the node's own Poisson counting floor `1/max(g,1)` and the kernel resolution floor `S0` —
    so `w = ref/(v+ref)` is the irreducible share against the deconvolution AMBIGUITY `v`. Confident
    high-count nodes keep mass (v~S0 ⇒ w~0.5) so REAL enriched modes survive; give-up nodes (v >> ref)
    collapse toward 0. A zero-count structural node is the trusted zero-DNA anchor, w = 1.

    Exposed so the PRECISION treatment can be ablated independently of the SUBSTRATE."""
    g = s["g_hat"][sel]
    v = np.maximum(np.nan_to_num(s["var"][sel], nan=1e9, posinf=1e9), 0.0)
    ref = 1.0 / np.maximum(g, 1.0) + _S0
    return np.where(mk["struct_zero"][sel], 1.0, ref / (v + ref))


def _convolve(d, h):
    """Mass-preserving Gaussian convolution of an (unnormalised) density on GRID. `h <= 0` is a no-op."""
    if h <= 0:
        return d
    k = np.exp(-0.5 * ((GRID[:, None] - GRID[None, :]) / h) ** 2)
    return (k / np.maximum(k.sum(0, keepdims=True), _EPS)) @ d


def smooth(d, h):
    """Convolve a density on GRID with a Gaussian of width `h` decades — a GLOBAL population resolution.

    By linearity this equals widening every per-node kernel before summing, so it is the cheap way to
    express `h_i² = h_within,i² + h_pop²` when `h_pop` is constant.

    ⚠ **Measured 2026-07-27: a global `h_pop` is the WRONG instrument.** The oracle's depleted mode is
    itself sharp (peak 0.193; the fit at h=0 reaches 0.141, i.e. already UNDER-peaked), so global smoothing
    flattens the one region the fit gets right (0.026 at h=0.30, 7× too flat) while the sparse enriched
    tail — where the spurious modes actually live — needs MORE width. EMD and held-out likelihood both
    prefer h=0. Use :func:`recipe`'s ``adaptive=True`` instead."""
    if d is None or h <= 0:
        return d
    return _norm(_convolve(d, h))


def _poisson_kernels(g, E):
    """Per-node zero-native Poisson posterior on GRID, row-normalised to unit mass."""
    lam = np.exp(_lnat)[None, :] * E[:, None]
    ll = g[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    return pn / np.maximum(pn.sum(1, keepdims=True), _EPS)


def knn_widths(a, knn_scale, floor=None):
    """THE population-resolution law, in one place: ``h_i = scale · dist(a_i, k-th nearest neighbour)``,
    ``k = √n``. Spacing is the quantity that decides whether two kernels merge into one mode or stand apart
    as two spurious ones, so it — not the density value, and never a node's own measurement precision — is
    what sets the rendering width.

    `floor` defaults to :data:`GRID_H`, the axis's own resolution. That floor is **forced, not chosen**: no
    estimator on a grid of spacing `GRID_H` can represent a feature narrower than one cell, so a kernel
    below it is not a narrow density, it is a delta at the wrong height."""
    n = a.size
    kk = max(int(round(np.sqrt(n))), 2)
    srt = np.sort(a)
    idx = np.searchsorted(srt, a)
    lo, hi = np.maximum(idx - kk, 0), np.minimum(idx + kk, n - 1)
    h = float(knn_scale) * np.maximum(np.maximum(a - srt[lo], srt[hi] - a), GRID_H)
    return np.maximum(h, GRID_H if floor is None else float(floor))


def _render(g, E, w, knn_scale, n_bins=12, floor=None):
    """Sum weighted per-node kernels at the POPULATION resolution — the one rule the fit and the oracle
    reference both use, so a fit-vs-reference comparison isolates the counts and the weights rather than
    the smoother. Binned by width so the convolution stays a handful of matrix products."""
    if not len(g):
        return None
    pn = _poisson_kernels(g, E)
    a = np.log10(np.maximum(g, 1.0)) - np.log10(np.maximum(E, _EPS))
    h_i = knn_widths(a, knn_scale, floor)
    edges = np.quantile(h_i, np.linspace(0.0, 1.0, int(n_bins) + 1))
    out = np.zeros_like(GRID)
    for b in range(int(n_bins)):
        m = (h_i >= edges[b]) & (h_i <= edges[b + 1] if b == n_bins - 1 else h_i < edges[b + 1])
        if m.any():
            out += _convolve((w[m][:, None] * pn[m]).sum(0), float(np.mean(h_i[m])))
    tot = float(out.sum())
    return out / tot if tot > 0 and np.isfinite(tot) else None


def recipe(s, sel=None, w=None, h_pop=0.0, knn_scale=0.0, n_bins=12):
    """The 2026-07-21 'unified' landscape: one reliability-mass weight, zero-native additive Poisson.

    `sel` overrides the training substrate (default: :func:`recipe_substrate`), so substrate variants can be
    A/B'd **with the weighting held fixed** — dropping a node class and the weight at the same time is not a
    substrate experiment.

    Weight `w = ref/(v+ref)` with `ref = 1/max(g,1) + S0` — the two IRREDUCIBLE variance sources of the log
    gDNA rate (own Poisson counting floor + kernel resolution floor) over the deconvolution AMBIGUITY `v`.
    Flooring by S0 keeps confident high-count nodes (v~S0 ⇒ w~0.5) so REAL enriched modes survive, while
    give-up nodes (v >> ref) collapse to ~0 mass.

    ⚠ Its bandwidth is ENTIRELY per-node (the Poisson count), with no population resolution term — which is
    the measured overfitting defect (mode count rises as the training set shrinks). G1 of the production
    plan replaces this kernel; the substrate and weight are the parts worth keeping."""
    mk = masks(s)
    if sel is None:
        sel = recipe_substrate(s, mk)
    if not sel.any():
        return None
    g, E = s["g_hat"][sel], s["eff"][sel]
    w = recipe_weights(s, sel, mk) if w is None else np.asarray(w, float)
    lam = np.exp(_lnat)[None, :] * E[:, None]
    ll = g[:, None] * np.log(np.maximum(lam, _EPS)) - lam - gammaln(g[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll)
    pn /= np.maximum(pn.sum(1, keepdims=True), _EPS)
    d = (w[:, None] * pn).sum(0)
    ss = float(d.sum())
    if not (ss > 0 and np.isfinite(ss)):
        return None
    if knn_scale <= 0:
        return smooth(d / ss, float(h_pop))

    # ── ADAPTIVE population resolution by NEAREST-NEIGHBOUR SPACING ───────────────────────────────────
    # The quantity that decides whether two per-node kernels merge into one mode — or stand apart as two
    # spurious ones — is the SPACING between observations on the density axis, so that is what sets the
    # width: h_i = scale · (distance from node i to its k-th nearest neighbour), k = sqrt(n).
    #
    # It fixes the overfitting BY CONSTRUCTION and with no tuning: fewer training nodes ⇒ farther
    # neighbours ⇒ wider kernels. Measured (mode count vs training-set size, oracle = 3):
    #   h=0      4.0 → 9.6 modes as n shrinks 10×   (slope +5.6)
    #   kNN 0.5  3.0 → 2.0                          (slope −1.0, and exactly 3 at full data)
    #
    # ⚠ Two rejected alternatives, both measured: a GLOBAL h_pop (flattens the genuinely-sharp depleted
    # bulk to fix the sparse tail — EMD and held-out likelihood both prefer h=0), and ABRAMSON
    # sample-point adaptivity h_i ∝ f(a_i)^(−1/2), which is ~neutral because a spike IS high pilot density
    # at its own location and so keeps itself narrow.
    #
    # ⚠⚠ IT IS ALSO THE ONLY THING STANDING BETWEEN THIS ESTIMATOR AND A COMB (measured 2026-07-27).
    # The Poisson kernel's width on the log axis is `1/(√g·ln10)`, which SHRINKS as ρ^(−1/2): across this
    # suite it collapses 41× (0.32 → 0.008 decades) and goes below `GRID_H`, at which point 88 % of the
    # nodes above +1 decade — carrying 100 % of that band's mass — are deltas in a single cell. Roughness
    # (total variation of log P per decade) above +1: **h=0 scores 46.9, kNN 0.5 scores 5.1.**
    return _render(g, E, w, knn_scale, n_bins)


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
