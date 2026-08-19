"""Unit tests for the generic density-deconvolution primitive (intron factory = special case).

Phase 1 = the factor math: the NegBinom log-pmf, and the per-intron λ-factor's mode / precision / regimes.
The end-to-end fit-vs-oracle validation lives in scripts/scratch (cached ambig scenarios), not here.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import nbinom, poisson

from rigel.calibration.density_deconv import (
    GdnaBackground,
    _log_negbinom,
    density_lambda_factor,
    fit_gdna_background,
)

_GRID = np.linspace(1e-6, 1.0 - 1e-6, 400)  # dense f_g grid to locate the peak


def _bg(log_mu_bg, alpha=np.inf, size=1.0e6, informative=True):
    return GdnaBackground(
        log_mu_bg=float(np.log(log_mu_bg)) if log_mu_bg > 0 else -np.inf,
        alpha=float(alpha),
        size=float(size),
        n_regions=100,
        informative=informative,
    )


def _peak_fg(factor_row):
    return _GRID[int(np.argmax(factor_row))]


# ---- the NegBinom log-pmf primitive ----


def test_negbinom_matches_scipy_at_integer_g():
    """_log_negbinom (mean/size param) matches scipy.stats.nbinom at integer counts, finite α."""
    mu = np.array([5.0, 20.0, 100.0])
    r = 8.0  # size
    g = np.array([3.0, 18.0, 90.0])
    p = r / (r + mu)  # scipy nbinom uses (n=size, p); mean = n(1-p)/p
    expected = nbinom.logpmf(g, r, p)
    got = _log_negbinom(g, mu, r)
    assert np.allclose(got, expected, atol=1e-9)


def test_negbinom_poisson_limit():
    """α = ∞ is the exact Poisson log-pmf (the r→∞ limit, no Γln(∞) overflow)."""
    mu = np.array([2.0, 30.0, 300.0])
    g = np.array([2.0, 25.0, 310.0])
    got = _log_negbinom(g, mu, np.inf)
    assert np.allclose(got, poisson.logpmf(g, mu), atol=1e-9)


def test_negbinom_large_size_approaches_poisson():
    """A large-but-finite size converges to Poisson (continuity of the parameterization).

    Integer g only — scipy's poisson.logpmf is −inf off the integers, while _log_negbinom is continuous."""
    mu = np.full(5, 40.0)
    g = np.array([10.0, 25.0, 40.0, 60.0, 80.0])
    near = _log_negbinom(g, mu, 1.0e7)
    pois = poisson.logpmf(g, mu)
    assert np.max(np.abs(near - pois)) < 1e-2


# ---- the λ-factor: mode, regimes, precision ----


def test_factor_peaks_at_rho_bg_over_rho_obs():
    """The gDNA deconvolve lands at f_g = ρ_bg/ρ_obs = μ/C (the confident background deconvolve)."""
    rho_bg, E_g = 0.5, 1000.0
    C = np.array([2000.0])  # ρ_obs = 2.0 ⇒ f_g* = 0.25
    fac = density_lambda_factor(_bg(rho_bg), C, np.array([E_g]), _GRID)
    assert _peak_fg(fac[0]) == pytest.approx(rho_bg / (C[0] / E_g), abs=0.03)


def test_regime_pure_gdna_fg_near_one():
    """ρ_obs ≈ ρ_bg (no nascent) ⇒ f_g deconvolves near 1 (all gDNA)."""
    rho_bg, E_g = 0.5, 1000.0
    C = np.array([rho_bg * E_g])  # ρ_obs = ρ_bg
    fac = density_lambda_factor(_bg(rho_bg), C, np.array([E_g]), _GRID)
    assert _peak_fg(fac[0]) > 0.9


def test_regime_nascent_present_peels_to_background():
    """ρ_obs ≫ ρ_bg ⇒ f_g ≈ ρ_bg/ρ_obs ≪ 1: the excess is confidently nascent."""
    rho_bg, E_g = 0.5, 1000.0
    C = np.array([10.0 * rho_bg * E_g])  # 10× above background
    fac = density_lambda_factor(_bg(rho_bg), C, np.array([E_g]), _GRID)
    assert _peak_fg(fac[0]) == pytest.approx(0.1, abs=0.03)


def test_regime_dna_free_pins_low():
    """Σg = 0 (DNA-free): the posterior sits at ½/ΣE with honest width ⇒ the factor pulls f_g to ~0.

    ⭐ Through :func:`fit_gdna_background` rather than a hand-built background, because the DNA-free
    regime is exactly where the hand-built wall and the real fit used to diverge (the 2026-08-18
    detonation) — this regime test must exercise the real path."""
    bg = fit_gdna_background(np.zeros(200), np.full(200, 5_000.0))  # a genuinely empty pool
    C = np.array([2.0])  # a sparse intron, ρ_obs = 2/1000 = 2e-3
    fac = density_lambda_factor(bg, C, np.array([1000.0]), _GRID)
    assert _peak_fg(fac[0]) < 0.2


def test_precision_monotone_in_count():
    """Honest count-over-length: a higher-count intron deconvolves f_g more sharply (narrower factor)."""
    rho_bg, E_g = 0.5, 1000.0

    def width(C):
        fac = density_lambda_factor(_bg(rho_bg), np.array([C]), np.array([E_g]), _GRID)[0]
        w = np.exp(fac - fac.max())
        w /= w.sum()
        m = (w * _GRID).sum()
        return float(np.sqrt((w * (_GRID - m) ** 2).sum()))  # posterior std over the grid

    # scale C and E_g together so f_g* is fixed but the count grows
    lo = width(2000.0)
    hi = width(20000.0)  # 10× the count at the same ρ_obs
    assert hi < lo  # more count ⇒ sharper


def test_overdispersion_widens_the_peel():
    """A finite α (per-region CNV/mappability spread) widens the factor vs the Poisson (α=∞) case."""
    rho_bg, E_g, C = 0.5, 1000.0, np.array([2000.0])

    def width(alpha):
        fac = density_lambda_factor(_bg(rho_bg, alpha=alpha), C, np.array([E_g]), _GRID)[0]
        w = np.exp(fac - fac.max())
        w /= w.sum()
        m = (w * _GRID).sum()
        return float(np.sqrt((w * (_GRID - m) ** 2).sum()))

    assert width(3.0) > width(np.inf)


def test_non_informative_is_flat():
    """An empty background pool ⇒ a flat (all-zero) factor: the factory says nothing."""
    bg = _bg(0.0, size=0.5, informative=False)
    fac = density_lambda_factor(bg, np.array([100.0, 5.0]), np.array([1000.0, 1000.0]), _GRID)
    assert fac.shape == (2, 400)
    assert np.allclose(fac, 0.0)


def test_factor_shape_and_finiteness():
    """Multi-intron call returns (n, K) finite, each row offset to max 0."""
    bg = _bg(0.3)
    C = np.array([100.0, 3000.0, 12.0])
    Eg = np.array([800.0, 1500.0, 400.0])
    fac = density_lambda_factor(bg, C, Eg, _GRID)
    assert fac.shape == (3, _GRID.shape[0])
    assert np.all(np.isfinite(fac))
    assert np.allclose(fac.max(axis=1), 0.0)


# ---- ⭐⭐⭐ the SMOOTH background posterior — the 2026-08-18 zero-gDNA catastrophe, pinned ----
#
# The shipped fit branched: pooled MLE at Σg>0, a "resolution wall" mean(1/E) fallback at Σg=0. That
# wall is a mean of reciprocals and is OWNED by the smallest regions of the partition
# (TRAPS: a-mean-of-ratios-inherits-the-partition): measured on the ladder, ONE intergenic sliver
# (E=0.0074 bp) carried 35 % of it and the wall landed at 0.2985/bp on a library whose true background
# was EXACTLY 0 — so the intron factory confidently manufactured mu = 0.2985*E phantom gDNA and 80 % of
# all nascent intron mass was called gDNA. The repair is the conjugate posterior
# rho_bg ~ Gamma(Sumg + 1/2, SumE): ONE formula, no branch, no wall, honest width, and its Sumg >> 1
# limit is the shipped pooled rate exactly. These four tests are the falsification set for that repair.


def test_the_background_location_is_sliver_invariant():
    """⛔ THE BUG, PINNED: one fragment-length intergenic sliver must not move the background."""
    E = np.full(50, 10_000.0)
    g = np.zeros(50)
    clean = fit_gdna_background(g, E)
    slivered = fit_gdna_background(np.append(g, 0.0), np.append(E, 0.01))
    assert abs(clean.log_mu_bg - slivered.log_mu_bg) < 0.01, (
        clean.log_mu_bg,
        slivered.log_mu_bg,
    )


def test_the_background_is_smooth_through_zero_counts():
    """⛔ One observed fragment may move the location by ~(1+1/2)/(0+1/2) = 3x — never by orders."""
    E = np.full(50, 10_000.0)
    at0 = fit_gdna_background(np.zeros(50), E)
    g1 = np.zeros(50)
    g1[7] = 1.0
    at1 = fit_gdna_background(g1, E)
    assert abs(at1.log_mu_bg - at0.log_mu_bg) <= np.log(3.0) + 1e-9


def test_an_empty_pool_is_not_confident():
    """⛔ Σg=0 says "around 1/(2ΣE), and I genuinely do not know" — the factor must neither place the
    deconvolve away from ~0 on a dense intron NOR carry populated-pool precision. The shipped path read
    tau ~ 3e8 here and called 80 % of nascent gDNA."""
    from rigel.calibration.density_deconv import density_factor_precision

    E = np.full(500, 10_000.0)
    empty = fit_gdna_background(np.zeros(500), E)
    # ⛔ the information HALF of the contract, pinned directly: an empty region is not a unit of
    #   Fisher information, so 500 empty regions carry exactly the Jeffreys ½ — never 500. A break
    #   that restores ``Σg + n0`` here changes no factor peak (the location is already tiny), so the
    #   behavioural gates below cannot see it; this line is the one that fires.
    assert empty.size == pytest.approx(0.5), empty.size
    C, Eg = np.array([10_000.0]), np.array([5_000.0])
    fac = density_lambda_factor(empty, C, Eg, _GRID)
    assert _peak_fg(fac[0]) < 0.01
    lam = np.log(_GRID) - np.log1p(-_GRID)
    tau_empty = density_factor_precision(fac, lam)[0]
    full = fit_gdna_background(np.full(500, 500.0), E)  # same support, populated
    tau_full = density_factor_precision(density_lambda_factor(full, C, Eg, _GRID), lam)[0]
    assert tau_empty < 1e-2 * tau_full, (tau_empty, tau_full)


def test_the_populated_limit_is_the_pooled_rate():
    """⭐ Σg ≫ 1 reduces to the shipped pooled MLE exactly: ln((Σg+1/2)/ΣE) − ln(Σg/ΣE) ~ 1/(2Σg)."""
    g = np.full(100, 10_000.0)
    E = np.full(100, 20_000.0)
    bg = fit_gdna_background(g, E)
    assert bg.log_mu_bg == pytest.approx(np.log(g.sum() / E.sum()), abs=1e-6)
    assert bg.informative
