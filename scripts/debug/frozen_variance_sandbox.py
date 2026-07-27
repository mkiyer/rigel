"""Frozen-variance sandbox — validate the count-zero-information fix (O2) before touching production.

Design: `docs/calibration/archive/node_prior_design.md` §2, §9 (the count-zero-info tradeoff).

THE CLAIM. `_mixture_strand_loglik` is a correct heteroscedastic Gaussian-BB, but its variance normalizer
(`-½·log var` + `(u-mean)²/var`) depends on the composition `f_g`, so when the mean channel degenerates
(κ→½) the *count magnitude* leaks a composition preference toward the variance-minimum `od_r/(od_g+od_r)`
— a count-zero-info violation. The fix FREEZES the variance's composition argument at a reference `f_ref`
(Fisher scoring): the LIVE mean stays in the numerator (the legitimate mean channel), the variance is
evaluated at `f_ref`. The count then sets precision only, never location.

We work on production's log-odds grid + `_posterior_median_fg` read-out (faithful to the shipped solver;
the log-odds grid resolves the f_g→0/1 vertices, so a wide posterior's tails are handled without the
hard-truncation artifact of a uniform-f_g grid). "current" = the imported production `_mixture_strand_loglik`
(byte-faithful); the sandbox reimplements only the frozen variant. Both paths share the SAME grid, Jacobian
measure, and prior, so any current-vs-frozen difference is the freeze alone.

Primary, grid-invariant evidence is the strand loglik's *f_g-dependence* at κ=½ (current > 0, frozen = 0).
Medians are the shipped read-out and are reported too, but the mechanism lives in the loglik shape.

Run inside the activated `rigel` env:  python scripts/debug/frozen_variance_sandbox.py
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.special import log_expit

from rigel.calibration.simplex import _mixture_strand_loglik
from rigel.calibration.simplex_logodds import _logodds_grid, _posterior_median_fg

_EPS = 1.0e-9
_LAM, _FG = _logodds_grid(4001, 10.0)  # production log-odds grid; _FG = σ(λ) ascending
_JAC = log_expit(_LAM) + log_expit(-_LAM)  # log σ'(λ) = log f_g + log(1−f_g); uniform-f measure


# --------------------------------------------------------------------------------------------------
# Strand loglik over the f_g grid for a single-strand (+) node: f_pos = 1−f_g, f_neg = 0.
# --------------------------------------------------------------------------------------------------
def strand_current(u_pos, n, kappa, od_g, od_r, fg=_FG):
    """The OLD (pre-B1) leaky likelihood — variance live in f_g — reproduced by passing the grid
    composition as the (now-required) freeze reference. Contrast `strand_frozen`."""
    fg = np.asarray(fg, float)
    grid_g, grid_p, grid_n = fg[None, :], (1.0 - fg)[None, :], np.zeros_like(fg)[None, :]
    return _mixture_strand_loglik(
        np.array([[float(u_pos)]]),
        np.array([[float(n)]]),
        grid_g,
        grid_p,
        grid_n,
        kappa,
        od_g,
        od_r,
        grid_g,
        grid_p,
        grid_n,
    )[0]


def strand_frozen(u_pos, n, kappa, od_g, od_r, f_ref, fg=_FG):
    """The fix: LIVE mean (`n·p(f_g)`), variance FROZEN at `f_ref` (constant over the grid)."""
    fg = np.asarray(fg, float)
    mean = n * (0.5 * fg + kappa * (1.0 - fg))  # single-strand +, f_neg = 0 — the live mean channel
    p_ref = 0.5 * f_ref + kappa * (1.0 - f_ref)
    rscale = kappa * (1.0 - kappa)
    var_ref = max(
        n * p_ref * (1.0 - p_ref)
        + (n * f_ref) ** 2 * 0.25 * od_g
        + (n * (1.0 - f_ref)) ** 2 * rscale * od_r,
        _EPS,
    )
    return -0.5 * (u_pos - mean) ** 2 / var_ref - 0.5 * np.log(var_ref)


# --------------------------------------------------------------------------------------------------
# Priors + posterior read-outs (production log-odds measure).
# --------------------------------------------------------------------------------------------------
def logprior_flat(fg=_FG):
    return np.zeros_like(np.asarray(fg, float))


def logprior_density_anchor(rho_obs, rho_bg, sigma2, fg=_FG):
    """The nascent≈0 prior (§4.2): prior on log ρ_g = log(f_g·ρ_obs) ~ 𝒩(log ρ_bg, σ²); lone-mode at
    f_g = clip(ρ_bg/ρ_obs, 0, 1) → 1 at background density, yielding above it."""
    fg = np.clip(np.asarray(fg, float), _EPS, 1.0)
    lr = np.log(fg * rho_obs) - np.log(rho_bg)
    return -0.5 * lr * lr / sigma2


def _post(strand, logprior):
    psi = strand + _JAC + logprior
    psi = psi - psi.max()
    w = np.exp(psi)
    return w / w.sum()


def median_fg(strand, logprior):
    return float(_posterior_median_fg(_post(strand, logprior)[None, :], _FG)[0])


def map_fg(strand, logprior):
    return float(_FG[int(np.argmax(strand + _JAC + logprior))])


def var_logfg(strand, logprior):
    w = _post(strand, logprior)
    lf = np.log(np.clip(_FG, _EPS, 1.0))
    m = float((w * lf).sum())
    return max(float((w * lf * lf).sum()) - m * m, 0.0)


# --------------------------------------------------------------------------------------------------
# Tests. Each returns (name, passed, detail).
# --------------------------------------------------------------------------------------------------
def test_A_count_zero_info_at_half():
    """κ=½, balanced counts, asymmetric od. THE reviewer's acceptance test. Grid-invariant mechanism:
    the frozen strand loglik is FLAT in f_g (no count→composition path); the current one is not. Read-out:
    the frozen median is invariant to n while the current median is pulled off the flat answer & drifts."""
    kappa, od_g, od_r = 0.5, 0.1, 0.3
    ns = [10, 100, 1000, 10000]
    range_c, range_f, med_c, med_f, var_f = [], [], [], [], []
    for n in ns:
        u = 0.5 * n  # p̂ = 0.5
        sc = strand_current(u, n, kappa, od_g, od_r)
        sf = strand_frozen(u, n, kappa, od_g, od_r, f_ref=1.0)
        range_c.append(float(np.ptp(sc)))
        range_f.append(float(np.ptp(sf)))
        med_c.append(median_fg(sc, logprior_flat()))
        med_f.append(median_fg(sf, logprior_flat()))
        var_f.append(var_logfg(sf, logprior_flat()))
    frozen_strand_flat = max(range_f) < 1.0e-9  # THE mechanism (grid-invariant)
    current_strand_not_flat = min(range_c) > 0.2  # current injects f_g info from the count alone
    frozen_median_invariant = (max(med_f) - min(med_f)) < 1.0e-2
    current_median_biased = abs(np.mean(med_c) - np.mean(med_f)) > 0.02  # leak vs the flat answer
    frozen_var_invariant = (max(var_f) - min(var_f)) < 1.0e-3  # κ=½ carries ZERO f_g info (correct)
    passed = (
        frozen_strand_flat
        and current_strand_not_flat
        and frozen_median_invariant
        and current_median_biased
        and frozen_var_invariant
    )
    detail = (
        f"n={ns}\n"
        f"      strand loglik f_g-range: current={[round(x, 2) for x in range_c]} (not-flat="
        f"{current_strand_not_flat}); frozen={[f'{x:.1e}' for x in range_f]} (flat={frozen_strand_flat})\n"
        f"      median: current={[round(x, 3) for x in med_c]} (biased-off-flat={current_median_biased}); "
        f"frozen={[round(x, 3) for x in med_f]} (invariant={frozen_median_invariant})\n"
        f"      frozen Var(log f_g)={[round(x, 3) for x in var_f]} (n-invariant={frozen_var_invariant} — "
        f"κ=½ carries zero f_g info; the count→precision role is Test D)"
    )
    return "A: count-zero-info @ κ=½ (reviewer's invariance)", passed, detail


def test_B_stranded_no_regression():
    """κ=0.8 stranded. With a REASONABLE reference (neutral ½ or near-truth) the frozen median agrees
    with current — the mean channel dominates. A far reference (f_ref=1 for a truly-RNA node) is reported
    to quantify the sensitivity (Test C studies it)."""
    kappa, od_g, od_r, n = 0.8, 0.2, 0.2, 2000
    fg_true = 0.3
    u = (kappa - (kappa - 0.5) * fg_true) * n  # p̂ = 0.71
    mc = median_fg(strand_current(u, n, kappa, od_g, od_r), logprior_flat())
    mfs = {
        fr: median_fg(strand_frozen(u, n, kappa, od_g, od_r, fr), logprior_flat())
        for fr in (0.3, 0.5, 1.0)
    }
    agree_reasonable = abs(mc - mfs[0.3]) < 0.02 and abs(mc - mfs[0.5]) < 0.02
    detail = (
        f"p̂={u / n:.3f}; current median={mc:.3f}; frozen "
        f"{{{', '.join(f'{k}:{v:.3f}' for k, v in mfs.items())}}}; "
        f"far-ref (1.0) shift={abs(mc - mfs[1.0]):.3f}"
    )
    return "B: stranded no-regression (κ=0.8, reasonable ref)", agree_reasonable, detail


def test_C_reference_insensitivity():
    """The reference is low-stakes for LOCATION at κ=½ (prior decides regardless of f_ref) and MODEST
    (not zero) at κ≫½ under overdispersion — the reviewer's #2, quantified honestly."""
    refs = [0.0, 0.25, 0.5, 0.75, 1.0]
    anchor = logprior_density_anchor(rho_obs=1.0, rho_bg=1.0, sigma2=0.5)  # background → f_g→1
    med_half = [median_fg(strand_frozen(0.5 * 500, 500, 0.5, 0.1, 0.3, fr), anchor) for fr in refs]
    n, kappa = 2000, 0.8
    u = (kappa - (kappa - 0.5) * 0.3) * n
    med_strand = [
        median_fg(strand_frozen(u, n, kappa, 0.2, 0.2, fr), logprior_flat()) for fr in refs
    ]
    half_spread = max(med_half) - min(med_half)
    strand_spread = max(med_strand) - min(med_strand)
    half_ok = half_spread < 1.0e-3  # κ=½: prior fully decides
    strand_ok = strand_spread < 0.05  # κ=0.8: modest, not negligible (validates prior-anchor rec)
    passed = half_ok and strand_ok
    detail = (
        f"f_ref={refs}\n"
        f"      κ=½ (anchor) median → {[round(x, 3) for x in med_half]} (spread={half_spread:.1e}, "
        f"prior-decides={half_ok})\n"
        f"      κ=0.8       median → {[round(x, 3) for x in med_strand]} (spread={strand_spread:.3f} "
        f"< 0.05={strand_ok}) — reference matters under overdispersion ⇒ use the prior anchor"
    )
    return "C: reference sensitivity (honest, prior-anchor claim)", passed, detail


def test_D_precision_monotonic_in_n():
    """Under the frozen variance the count sets PRECISION (the count-zero-info-compliant role). At LOW
    overdispersion the count drives it clearly; the posterior precision rises monotonically with n."""
    kappa, od_g, od_r, fg_true = 0.8, 0.01, 0.01, 0.4
    ns = [50, 200, 1000, 5000]
    prec = []
    for n in ns:
        u = (kappa - (kappa - 0.5) * fg_true) * n
        prec.append(
            1.0
            / max(var_logfg(strand_frozen(u, n, kappa, od_g, od_r, fg_true), logprior_flat()), _EPS)
        )
    monotone = all(prec[i + 1] > prec[i] for i in range(len(prec) - 1))
    detail = f"n={ns}; precision → {[round(x, 1) for x in prec]} (increasing={monotone})"
    return "D: precision monotonic in n (frozen, low od)", monotone, detail


def test_E_leak_exact_prediction():
    """The leak's falsifiable prediction, on the PURE strand term (no measure/prior): at κ=½ the current
    strand MODE = the variance-min od_r/(od_g+od_r); the frozen strand is flat (no mode)."""
    kappa, n = 0.5, 5000
    pairs = [(0.1, 0.3), (0.3, 0.1), (0.2, 0.2), (0.05, 0.4)]
    rows, ok = [], True
    for od_g, od_r in pairs:
        pred = od_r / (od_g + od_r)
        u = 0.5 * n
        sc = strand_current(u, n, kappa, od_g, od_r)
        mode_c = float(_FG[int(np.argmax(sc))])  # pure strand mode (no Jacobian/prior)
        flat_f = float(np.ptp(strand_frozen(u, n, kappa, od_g, od_r, 1.0))) < 1.0e-9
        hit = abs(mode_c - pred) < 0.02 and flat_f
        ok = ok and hit
        rows.append(
            f"od=({od_g},{od_r}) pred={pred:.3f} current-mode={mode_c:.3f} frozen-flat={flat_f} ok={hit}"
        )
    return "E: leak's exact var-min prediction (MODE)", ok, "\n      " + "\n      ".join(rows)


def test_F_median_vs_map_stability():
    """Why §6/§7 read the MEDIAN, not the MAP (880ac429's failure). At κ=0.499 with a FLAT prior and a
    tiny strand imbalance, the current-likelihood MAP tips while the median is stable. (With the nascent≈0
    anchor both stabilize — the prior is the other half of the fix.)"""
    kappa, n, od_g, od_r = 0.499, 500, 0.2, 0.2
    phats = (0.50, 0.51, 0.52, 0.53)
    med_flat, map_flat = [], []
    for phat in phats:
        sc = strand_current(phat * n, n, kappa, od_g, od_r)
        med_flat.append(median_fg(sc, logprior_flat()))
        map_flat.append(map_fg(sc, logprior_flat()))
    # with the anchor prior both should be pinned near f_g→1
    anchor = logprior_density_anchor(rho_obs=1.0, rho_bg=1.0, sigma2=0.5)
    med_anch = [median_fg(strand_current(p * n, n, kappa, od_g, od_r), anchor) for p in phats]
    med_flat_range = max(med_flat) - min(med_flat)
    map_flat_range = max(map_flat) - min(map_flat)
    anchor_stable = (max(med_anch) - min(med_anch)) < 0.05
    median_beats_map = med_flat_range < map_flat_range or map_flat_range < 1.0e-6
    passed = median_beats_map and anchor_stable
    detail = (
        f"κ=0.499, p̂∈{phats}\n"
        f"      flat prior:   median range={med_flat_range:.3f}  MAP range={map_flat_range:.3f} "
        f"(median≤MAP={median_beats_map})\n"
        f"      anchor prior: median={[round(x, 3) for x in med_anch]} (stable={anchor_stable})"
    )
    return "F: median vs MAP stability @ κ=0.499", passed, detail


def main():
    tests = [
        test_A_count_zero_info_at_half,
        test_B_stranded_no_regression,
        test_C_reference_insensitivity,
        test_D_precision_monotonic_in_n,
        test_E_leak_exact_prediction,
        test_F_median_vs_map_stability,
    ]
    print("=" * 92)
    print("FROZEN-VARIANCE SANDBOX — count-zero-information fix (O2)")
    print("=" * 92)
    n_pass = 0
    for t in tests:
        name, passed, detail = t()
        n_pass += int(passed)
        print(f"\n[{'PASS' if passed else 'FAIL'}] {name}\n      {detail}")
    print("\n" + "=" * 92)
    print(f"{n_pass}/{len(tests)} tests passed")
    print("=" * 92)
    return 0 if n_pass == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())
