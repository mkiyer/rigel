"""The generic density-deconvolution primitive, of which the intron factory is one special case.

What is gated here is the factor arithmetic — the kernel's, read through `native.transfer_rows` — the
NegBinom log-pmf against scipy's own, and the
per-intron λ-factor's mode, precision and limiting regimes — a confident background deconvolve, the
no-nascent and nascent-rich extremes, the sharpening with count, the widening with overdispersion,
and the flat factor an empty pool must produce. The background posterior itself is held at the end,
where the interesting failure lives: a fit that goes confident on a pool with no counts in it.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import nbinom, poisson

from rigel.calibration.density_deconv import GdnaBackground, fit_gdna_background
from rigel.native import transfer_rows as R

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


def _log_negbinom(g, mu, size):
    """The kernel's NegBinom cell (`native.transfer_rows.log_negbinom`), elementwise."""
    g, mu = np.broadcast_arrays(np.asarray(g, np.float64), np.asarray(mu, np.float64))
    return R.log_negbinom(np.ascontiguousarray(g), np.ascontiguousarray(mu), float(size))


def density_lambda_factor(bg, count, eff_g, fg_grid):
    """The kernel's factory rows (`native.transfer_rows.factory_rows`) on a fraction grid, every slot an
    intron: the row ``log NegBinom(f_g·C; ρ_bg·E, α_eff)``, max-normalised."""
    fg = np.asarray(fg_grid, np.float64)
    count = np.ascontiguousarray(count, np.float64)
    return R.factory_rows(
        np.ones(count.shape[0], bool),
        count,
        np.ascontiguousarray(eff_g, np.float64),
        float(bg.log_mu_bg),
        float(bg.alpha),
        float(bg.size),
        bool(bg.informative),
        np.log(fg) - np.log1p(-fg),
    )


def density_factor_precision(rows, lam):
    """The kernel's factor precision (`native.transfer_rows.factor_precision`): 1/Var_λ under each row."""
    return R.factor_precision(
        np.ascontiguousarray(rows, np.float64), np.ascontiguousarray(lam, np.float64)
    )


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

    Through :func:`fit_gdna_background` rather than a hand-built background, because the DNA-free
    regime is exactly where a hand-built wall and the real fit can diverge — this regime test must
    exercise the real path."""
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


# ---- the SMOOTH background posterior, and why it may not branch at zero counts ----
#
# A fit that branches — a pooled MLE at Sumg>0 and a "resolution wall" mean(1/E) fallback at Sumg=0 —
# puts a mean of reciprocals on the zero-count path, and that is OWNED by the smallest regions of the
# partition (TRAPS: a-mean-of-ratios-inherits-the-partition): a single sub-fragment-length intergenic
# sliver can carry a third of it, so the wall lands at a real per-base rate on a library whose true
# background is EXACTLY 0, the intron factory manufactures phantom gDNA from it, and most nascent
# intron mass is called gDNA. The conjugate posterior rho_bg ~ Gamma(Sumg + 1/2, SumE) is one formula
# with no branch and no wall, it has an honest width, and its Sumg >> 1 limit is the pooled rate
# exactly. These four tests are the falsification set for it.


def test_the_background_location_is_sliver_invariant():
    """One fragment-length intergenic sliver must not move the background."""
    E = np.full(50, 10_000.0)
    g = np.zeros(50)
    clean = fit_gdna_background(g, E)
    slivered = fit_gdna_background(np.append(g, 0.0), np.append(E, 0.01))
    assert abs(clean.log_mu_bg - slivered.log_mu_bg) < 0.01, (
        clean.log_mu_bg,
        slivered.log_mu_bg,
    )


def test_the_background_is_smooth_through_zero_counts():
    """One observed fragment may move the location by ~(1+1/2)/(0+1/2) = 3x — never by orders."""
    E = np.full(50, 10_000.0)
    at0 = fit_gdna_background(np.zeros(50), E)
    g1 = np.zeros(50)
    g1[7] = 1.0
    at1 = fit_gdna_background(g1, E)
    assert abs(at1.log_mu_bg - at0.log_mu_bg) <= np.log(3.0) + 1e-9


def test_an_empty_pool_is_not_confident():
    """Σg=0 says "around 1/(2ΣE), and I genuinely do not know" — the factor must neither place the
    deconvolve away from ~0 on a dense intron NOR carry populated-pool precision. A branching fit
    reads an enormous precision here and calls most nascent intron mass gDNA."""
    E = np.full(500, 10_000.0)
    empty = fit_gdna_background(np.zeros(500), E)
    # the information HALF of the contract, pinned directly: an empty region is not a unit of Fisher
    #   information, so 500 empty regions carry exactly the Jeffreys ½ — never 500. PERTURBATION:
    #   restoring ``Σg + n0`` here changes no factor peak (the location is already tiny), so the
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
    """Σg ≫ 1 reduces to the shipped pooled MLE exactly: ln((Σg+1/2)/ΣE) − ln(Σg/ΣE) ~ 1/(2Σg)."""
    g = np.full(100, 10_000.0)
    E = np.full(100, 20_000.0)
    bg = fit_gdna_background(g, E)
    assert bg.log_mu_bg == pytest.approx(np.log(g.sum() / E.sum()), abs=1e-6)
    assert bg.informative


def test_the_factor_precision_is_chunk_exact():
    """`density_factor_precision` read one row at a time equals the whole array's, to the bit: the
    precision is a per-row moment and must not depend on which rows share the call — the locus
    solve reads it per block. A BLAS matrix-vector product breaks this at a one-row call (the
    library dispatches a different kernel), so the moments are per-row sums instead."""
    rng = np.random.default_rng(11)
    lam = np.linspace(-10.0, 10.0, 60)
    rows = (
        -0.5
        * ((lam[None, :] - rng.normal(0.0, 3.0, 37)[:, None]) / rng.uniform(0.5, 4.0, 37)[:, None])
        ** 2
    )
    rows[5] = 0.0  # a flat row: no information, τ = 0 in every tiling
    whole = density_factor_precision(rows, lam)
    singles = np.concatenate(
        [density_factor_precision(rows[i : i + 1], lam) for i in range(rows.shape[0])]
    )
    halves = np.concatenate(
        [density_factor_precision(rows[:20], lam), density_factor_precision(rows[20:], lam)]
    )
    assert whole[5] == 0.0 and (whole > 0).sum() == 36
    assert np.array_equal(singles, whole), (
        f"one-row calls moved {int((singles != whole).sum())} rows"
    )
    assert np.array_equal(halves, whole)
