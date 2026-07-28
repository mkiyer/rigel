"""Unit tests for the population gDNA-density hyperprior (`calibration.gdna_landscape`).

One test per property the design actually rests on, so a regression names itself.
"""

import numpy as np
import pytest

from rigel.calibration.gdna_landscape import (
    GdnaLandscape,
    fit_gdna_landscape,
    knn_widths,
)

LN10 = np.log(10.0)


def _two_mode(n_dep=340, n_enr=60, eff=500.0, seed=0):
    """A library shaped like the truth: a uniform depleted level plus a capture-enriched minority."""
    rng = np.random.default_rng(seed)
    rho = np.concatenate([np.full(n_dep, 0.05), np.full(n_enr, 25.0)])
    e = np.full(n_dep + n_enr, eff)
    count = rng.poisson(rho * e).astype(float)
    return count, np.maximum(count * 3.0, 10.0), e, np.full(count.size, 0.3)


def _density(ls):
    d = np.exp(ls.logP - ls.logP.max())
    return d / d.sum()


def test_fit_is_deterministic():
    args = _two_mode()
    anchor = np.zeros(args[0].size, bool)
    a = fit_gdna_landscape(*args, anchor=anchor)
    b = fit_gdna_landscape(*args, anchor=anchor)
    assert np.array_equal(a.logP, b.logP) and np.array_equal(a.log_rho, b.log_rho)


def test_recovers_both_modes():
    """The estimator exists to keep a capture-enriched MINORITY — 15 % of nodes here — so a fit that
    smooths them into the depleted bulk is the failure mode, not a nuance."""
    count, mass, eff, var = _two_mode()
    ls = fit_gdna_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    x, d = ls.log_rho / LN10, _density(ls)
    assert d[x < 0].sum() > 0.5, "depleted bulk lost"
    assert d[x > 0.5].sum() > 0.05, "enriched minority competed away"
    assert abs(x[d[x > 0.5].argmax() + int((x <= 0.5).sum())] - np.log10(25.0)) < 0.5


def test_grid_is_the_domain_logprior_is_asked_about():
    """ψ evaluates at ρ_g = f_g·M/E with f_g ≤ 1, so the grid top is max(M/E) and the bottom is the deepest
    one-count resolution wall. Nothing outside is representable, so nothing outside is represented."""
    count, mass, eff, var = _two_mode()
    ls = fit_gdna_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    x = ls.log_rho / LN10
    assert x[0] == pytest.approx(np.min(-np.log10(eff)))
    assert x[-1] == pytest.approx(np.max(np.log10(mass) - np.log10(eff)))


def test_zero_count_anchor_is_native_and_low():
    """A node with no gDNA must say 'ρ is anything below the wall' — a downward decay, NOT an invented
    location at 1/E. It is also the depleted anchor, and dropping it costs +0.26/+0.61 EMD."""
    eff = np.full(200, 500.0)
    ls = fit_gdna_landscape(
        np.zeros(200), np.full(200, 1000.0), eff, np.full(200, np.inf), anchor=np.ones(200, bool)
    )
    d = _density(ls)
    assert d.argmax() < len(d) // 4, "zero-count mass must pile at the bottom of the grid"
    assert np.all(np.diff(d[d.argmax() :]) <= 1e-12), (
        "the zero-count kernel must decay monotonically"
    )


def test_reliability_weight_damps_but_never_deletes():
    """Precision belongs in a continuous weight; expressing it as an admission threshold is measurably
    worse than ignoring precision altogether. So an imprecise node must be DAMPED — strictly smaller
    contribution — and never DELETED, because the one node reporting a genuine enriched mode may be exactly
    the imprecise one."""
    count, mass, eff, var = _two_mode()
    # move a single node far from every other, so its contribution is separable on the axis
    count, var = count.copy(), var.copy()
    count[0], mass[0] = 5000.0, 6000.0
    at = np.log10(count[0] / eff[0]) * LN10
    anchor = np.zeros(count.size, bool)

    def mass_near(v):
        var[0] = v
        ls = fit_gdna_landscape(count, mass, eff, var, anchor=anchor)
        near = np.abs(ls.log_rho - at) < 0.4 * LN10
        return float(_density(ls)[near].sum())

    confident, vague = mass_near(0.01), mass_near(50.0)
    assert vague < confident, "an imprecise node must be damped"
    assert vague > 0.0, "and never deleted — a cutoff is what this design rejects"


def test_anchor_is_trusted_outright():
    """The zero-count structural anchor carries w = 1: its density is 0 for EVERY f_g, so its composition
    ambiguity is irrelevant and must not down-weight it."""
    count, mass, eff, var = _two_mode()
    var = var.copy()
    var[:20] = np.inf  # unidentified composition, but zero mass ⇒ still an exact density statement
    count, mass = count.copy(), mass.copy()
    count[:20], mass[:20] = 0.0, 0.0
    anchored = fit_gdna_landscape(count, mass, eff, var, anchor=np.arange(count.size) < 20)
    unanchored = fit_gdna_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    low = _density(anchored)[: len(anchored.logP) // 4].sum()
    low_un = _density(unanchored)[: len(unanchored.logP) // 4].sum()
    assert low > low_un, "the anchor must carry full weight into the depleted floor"


def test_knn_width_never_below_the_grid_step():
    """Forced by the axis: a kernel narrower than one cell is a delta at the wrong height, which is what
    turned the enriched half of the landscape into a comb (roughness 46.9 vs a smooth bump's 2-4)."""
    a = np.linspace(-3.0, 2.0, 500)
    step = 0.02
    assert (knn_widths(a, step) >= step).all()
    assert (knn_widths(np.zeros(50), step) == step).all()  # degenerate: all nodes coincident


def test_knn_width_is_the_exact_kth_nearest_neighbour_distance():
    """Not "the far edge of a 2k window" — that hands the WIDEST kernel in the fit to the most ISOLATED
    node, which is backwards and was measured as the false-positive channel on zero-gDNA libraries. Checked
    against brute force on samples deliberately built as a bulk plus outliers."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = int(rng.integers(8, 300))
        a = rng.normal(size=n) * np.where(rng.random(n) < 0.12, 6.0, 1.0)
        k = max(int(round(np.sqrt(n))), 2)
        if n <= k:
            continue
        d = np.abs(a[:, None] - a[None, :])
        d.sort(1)
        assert np.allclose(knn_widths(a, 0.0, 1.0), d[:, k])


def test_knn_width_widens_as_the_sample_thins():
    """The self-correcting property: fewer nodes ⇒ farther neighbours ⇒ wider kernels, with no tuning."""
    full = knn_widths(np.linspace(-3, 2, 1000), 1e-6)
    thin = knn_widths(np.linspace(-3, 2, 100), 1e-6)
    assert np.median(thin) > np.median(full)


def test_logprior_shape_and_strength():
    count, mass, eff, var = _two_mode()
    ls = fit_gdna_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    fg = np.linspace(0.01, 0.99, 7)
    lp = ls.logprior(fg, np.full(4, 1000.0), np.full(4, 500.0))
    assert lp.shape == (4, 7) and np.isfinite(lp).all()
    half = GdnaLandscape(ls.log_rho, ls.logP, ls.n_train, 0.5)
    assert np.allclose(half.logprior(fg, np.full(4, 1000.0), np.full(4, 500.0)), 0.5 * lp)
    zero = GdnaLandscape(ls.log_rho, ls.logP, ls.n_train, 0.0)
    assert np.all(zero.logprior(fg, np.full(4, 1000.0), np.full(4, 500.0)) == 0.0)


def test_logprior_tracks_the_node_mass():
    """ρ_g = f_g·M/E, so at fixed f_g a heavier node sits higher on the landscape."""
    count, mass, eff, var = _two_mode()
    ls = fit_gdna_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    fg = np.array([0.5])
    lo = ls.logprior(fg, np.array([100.0]), np.array([500.0]))
    hi = ls.logprior(fg, np.array([100000.0]), np.array([500.0]))
    assert not np.allclose(lo, hi)


def test_declines_gracefully_on_degenerate_input():
    assert (
        fit_gdna_landscape(
            np.array([1.0]),
            np.array([1.0]),
            np.array([0.0]),
            np.array([1.0]),
            anchor=np.array([False]),
        )
        is None
    )
    assert (
        fit_gdna_landscape(
            np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), anchor=np.zeros(0, bool)
        )
        is None
    )
