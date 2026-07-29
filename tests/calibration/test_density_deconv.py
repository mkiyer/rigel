"""Unit tests for the generic density-deconvolution primitive (intron factory = special case) (docs/CARRY_FORWARD.md).

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
)

_GRID = np.linspace(1e-6, 1.0 - 1e-6, 400)  # dense f_g grid to locate the peak


def _bg(log_mu_bg, alpha=np.inf, sg=1.0e6, n0=0.0, informative=True):
    return GdnaBackground(
        log_mu_bg=float(np.log(log_mu_bg)) if log_mu_bg > 0 else -np.inf,
        alpha=float(alpha),
        sg=float(sg),
        n0=float(n0),
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
    """The gDNA peel lands at f_g = ρ_bg/ρ_obs = μ/C (the confident background peel)."""
    rho_bg, E_g = 0.5, 1000.0
    C = np.array([2000.0])  # ρ_obs = 2.0 ⇒ f_g* = 0.25
    fac = density_lambda_factor(_bg(rho_bg), C, np.array([E_g]), _GRID)
    assert _peak_fg(fac[0]) == pytest.approx(rho_bg / (C[0] / E_g), abs=0.03)


def test_regime_pure_gdna_fg_near_one():
    """ρ_obs ≈ ρ_bg (no nascent) ⇒ f_g peels near 1 (all gDNA)."""
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
    """Σg = 0 (DNA-free): μ falls to the resolution wall ⇒ the factor pulls f_g toward ~0."""
    # wall ρ_res tiny; ρ_obs modest ⇒ f_g ≈ ρ_res/ρ_obs ≈ 0
    bg = _bg(1e-4, alpha=np.inf, sg=0.0, n0=500.0)  # Σg=0, wall at 1e-4
    C = np.array([2.0])  # a sparse intron, ρ_obs = 2/1000 = 2e-3
    fac = density_lambda_factor(bg, C, np.array([1000.0]), _GRID)
    assert _peak_fg(fac[0]) < 0.2


def test_precision_monotone_in_count():
    """Honest count-over-length: a higher-count intron peels f_g more sharply (narrower factor)."""
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
    bg = _bg(0.0, sg=0.0, n0=0.0, informative=False)
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
