"""The population gDNA-density hyperprior (`calibration.landscape`), one property per test.

What the fit has to survive is a library whose truth is a uniform depleted level plus a
capture-enriched minority: the minority must stay a mode rather than being smoothed into the bulk,
the grid must span exactly what ψ can represent and nothing else, a zero-count region must decay
downward instead of inventing a location at 1/E, an imprecise region must be damped rather than
dropped, and the kernel widths must come from the local neighbour spacing so they widen on their own
as regions thin out. The census that splits the fit into the bulk and the enriched basins closes the file.
Each gate holds one property, so a regression names itself.
"""

import numpy as np
import pytest

import rigel.calibration.landscape as LS
from rigel.calibration.capture_efficiency import capture_efficiencies
from rigel.calibration.landscape import (
    DensityLandscape,
    LandscapeMode,
    census,
    fit_landscape,
    knn_widths,
    split_basins,
)
from rigel.config import CONSTANTS

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
    a = fit_landscape(*args, anchor=anchor)
    b = fit_landscape(*args, anchor=anchor)
    assert np.array_equal(a.logP, b.logP) and np.array_equal(a.log_rho, b.log_rho)


def test_the_row_tiles_move_no_bit(monkeypatch):
    """The kernels are built and summed a row tile at a time, and the tile's height is a resource setting
    (`simplex_logodds._block_rows`), so it must order no addition: a population spanning several tiles, fitted
    and refitted (the refit's E-step runs tile by tile), gives the same bits at a one-row tile, an odd tile, the
    shipped tile and one tile for the whole population. PERTURBATION: adding each tile's partial sum into the
    per-group totals, as the fit once did, fails here at every height."""
    rng = np.random.default_rng(5)
    n = 1600
    eff = rng.uniform(100.0, 3000.0, n)
    enriched = rng.random(n) < 0.15
    rho = np.where(enriched, 10.0 ** rng.uniform(-0.5, 1.0, n), 10.0 ** rng.uniform(-3.5, -2.5, n))
    count = rng.poisson(rho * eff).astype(float)
    mass = np.maximum(count * rng.uniform(1.5, 4.0, n), 1.0)
    var = rng.uniform(0.01, 0.9, n)
    anchor = (count == 0.0) & (rng.random(n) < 0.2)

    def fit(prev=None):
        return fit_landscape(count, mass, eff, var, anchor=anchor, prev=prev)

    shipped = LS._block_rows(CONSTANTS.landscape.grid_points, np.dtype(np.float64).itemsize)
    assert n > 2 * shipped, "the shipped tile must split the population, or the gate tests one tile"
    assert (count < 1.0).sum() > 100, "the refit's E-step must act on location-free kernels"
    first = fit()
    refit = fit(prev=first)
    for rows in (1, 7, 333, n):
        monkeypatch.setattr(LS, "_block_rows", lambda cells, itemsize, rows=rows: rows)
        for got, want in ((fit(), first), (fit(prev=first), refit)):
            for field in ("log_rho", "logP", "logP_free"):
                a, b = getattr(got, field), getattr(want, field)
                assert np.array_equal(a, b), (
                    f"rows {rows}, {field}: {int((a != b).sum())} cells moved"
                )


def test_recovers_both_modes():
    """The estimator exists to keep a capture-enriched MINORITY — 15 % of regions here — so a fit that
    smooths them into the depleted bulk is the failure mode, not a nuance."""
    count, mass, eff, var = _two_mode()
    ls = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    x, d = ls.log_rho / LN10, _density(ls)
    assert d[x < 0].sum() > 0.5, "depleted bulk lost"
    assert d[x > 0.5].sum() > 0.05, "enriched minority competed away"
    assert abs(x[d[x > 0.5].argmax() + int((x <= 0.5).sum())] - np.log10(25.0)) < 0.5


def test_grid_is_the_domain_logprior_is_asked_about():
    """ψ evaluates at ρ_g = f_g·M/E with f_g ≤ 1, so the grid top is max(M/E) and the bottom is the deepest
    one-count resolution wall. Nothing outside is representable, so nothing outside is represented."""
    count, mass, eff, var = _two_mode()
    ls = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    x = ls.log_rho / LN10
    assert x[0] == pytest.approx(np.min(-np.log10(eff)))
    assert x[-1] == pytest.approx(np.max(np.log10(mass) - np.log10(eff)))


def test_zero_count_anchor_is_native_and_low():
    """A region with no gDNA must say "ρ is anything below the wall" — a downward decay, NOT an
    invented location at 1/E. It is also the depleted anchor, so dropping it moves the whole fit."""
    eff = np.full(200, 500.0)
    ls = fit_landscape(
        np.zeros(200), np.full(200, 1000.0), eff, np.full(200, np.inf), anchor=np.ones(200, bool)
    )
    d = _density(ls)
    assert d.argmax() < len(d) // 4, "zero-count mass must pile at the bottom of the grid"
    assert np.all(np.diff(d[d.argmax() :]) <= 1e-12), (
        "the zero-count kernel must decay monotonically"
    )


def test_reliability_weight_damps_but_never_deletes():
    """Precision belongs in a continuous weight; expressing it as an admission threshold is measurably
    worse than ignoring precision altogether. So an imprecise region must be DAMPED — strictly smaller
    contribution — and never DELETED, because the one region reporting a genuine enriched mode may be exactly
    the imprecise one."""
    count, mass, eff, var = _two_mode()
    # move a single region far from every other, so its contribution is separable on the axis
    count, var = count.copy(), var.copy()
    count[0], mass[0] = 5000.0, 6000.0
    at = np.log10(count[0] / eff[0]) * LN10
    anchor = np.zeros(count.size, bool)

    def mass_near(v):
        var[0] = v
        ls = fit_landscape(count, mass, eff, var, anchor=anchor)
        near = np.abs(ls.log_rho - at) < 0.4 * LN10
        return float(_density(ls)[near].sum())

    confident, vague = mass_near(0.01), mass_near(50.0)
    assert vague < confident, "an imprecise region must be damped"
    assert vague > 0.0, "and never deleted — a cutoff is what this design rejects"


def test_anchor_is_trusted_outright():
    """The zero-count structural anchor carries w = 1: its density is 0 for EVERY f_g, so its composition
    ambiguity is irrelevant and must not down-weight it."""
    count, mass, eff, var = _two_mode()
    var = var.copy()
    var[:20] = np.inf  # unidentified composition, but zero mass ⇒ still an exact density statement
    count, mass = count.copy(), mass.copy()
    count[:20], mass[:20] = 0.0, 0.0
    anchored = fit_landscape(count, mass, eff, var, anchor=np.arange(count.size) < 20)
    unanchored = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    low = _density(anchored)[: len(anchored.logP) // 4].sum()
    low_un = _density(unanchored)[: len(unanchored.logP) // 4].sum()
    assert low > low_un, "the anchor must carry full weight into the depleted floor"


def test_knn_width_never_below_the_grid_step():
    """Forced by the axis: a kernel narrower than one cell is a delta at the wrong height, and the
    enriched half of the landscape becomes a comb of them rather than a bump."""
    a = np.linspace(-3.0, 2.0, 500)
    step = 0.02
    assert (knn_widths(a, step) >= step).all()
    assert (knn_widths(np.zeros(50), step) == step).all()  # degenerate: all regions coincident


def test_knn_width_is_the_exact_kth_nearest_neighbour_distance():
    """Not "the far boundary of a 2k window" — that hands the WIDEST kernel in the fit to the most
    ISOLATED region, which is backwards, and on a zero-gDNA library it is the channel that
    manufactures a false enriched mode. Checked against brute force on samples deliberately built as
    a bulk plus outliers."""
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
    """The self-correcting property: fewer regions ⇒ farther neighbours ⇒ wider kernels, with no tuning."""
    full = knn_widths(np.linspace(-3, 2, 1000), 1e-6)
    thin = knn_widths(np.linspace(-3, 2, 100), 1e-6)
    assert np.median(thin) > np.median(full)


def _arm(ls, fg, mass, eff):
    """ψ's fitted gDNA arm on a fraction grid, read through the kernel's own construction
    (`_psi_reference.gdna_arm`): the curve at ``log f_g + log M − log E`` per slot and cell."""
    from _psi_reference import gdna_arm

    fg = np.asarray(fg, np.float64)
    return gdna_arm(ls.log_rho, ls.logP, np.log(fg / (1.0 - fg)), mass, eff)


def test_the_arm_is_one_row_per_slot_and_finite():
    count, mass, eff, var = _two_mode()
    ls = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    lp = _arm(ls, np.linspace(0.01, 0.99, 7), np.full(4, 1000.0), np.full(4, 500.0))
    assert lp.shape == (4, 7) and np.isfinite(lp).all()


def test_the_arm_tracks_the_region_mass():
    """ρ_g = f_g·M/E, so at fixed f_g a heavier region sits higher on the landscape."""
    count, mass, eff, var = _two_mode()
    ls = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    lo = _arm(ls, [0.5], [100.0], [500.0])
    hi = _arm(ls, [0.5], [100000.0], [500.0])
    assert not np.allclose(lo, hi)


def test_the_kernels_arm_is_numpys_interpolation_of_the_curve_to_the_bit():
    """The arm ψ reads per cell (`native/psi_kernel.h`, ``Arm``) is the former ``(n, K)`` projection to the
    bit: ``np.interp`` of the curve at ``log σ(λ) + log M − log E`` with the ends held, the fraction, mass and
    opportunity clipped at the landscape's guard. Scored against numpy, a different implementation (TRAPS:
    a-test-that-redefines). PERTURBATION: an arm that extrapolated past the grid, or skipped a clip, fails here."""
    from scipy.special import expit

    from _psi_reference import gdna_arm

    from rigel.calibration.simplex_logodds import _logodds_grid

    count, mass, eff, var = _two_mode()
    ls = fit_landscape(count, mass, eff, var, anchor=np.zeros(count.size, bool))
    lam, _fg = _logodds_grid(101, 10.0)
    rng = np.random.default_rng(7)
    m = 40
    # two masses below the guard, one opportunity below it, and one heavy slot on a short opportunity that
    # lands above the curve's top: both held ends and every clip are exercised
    slot_mass = np.concatenate([rng.uniform(0.0, 3000.0, m - 3), [0.0, 1e-15, 1e6]])
    slot_eff = np.concatenate([rng.uniform(1.0, 5000.0, m - 2), [0.0, 1.0]])
    eps = 1e-12
    frac = np.clip(expit(lam), eps, 1.0 - eps)
    x = (
        np.log(frac)[None, :]
        + (np.log(np.maximum(slot_mass, eps)) - np.log(np.maximum(slot_eff, eps)))[:, None]
    )
    want = np.interp(x.ravel(), ls.log_rho, ls.logP, left=ls.logP[0], right=ls.logP[-1]).reshape(
        x.shape
    )
    got = gdna_arm(ls.log_rho, ls.logP, lam, slot_mass, slot_eff)
    assert np.array_equal(got, want), f"max |Δ| {np.abs(got - want).max():.3e}"
    assert (x < ls.log_rho[0]).any() and (x > ls.log_rho[-1]).any(), "the ends were not exercised"


def test_declines_gracefully_on_degenerate_input():
    assert (
        fit_landscape(
            np.array([1.0]),
            np.array([1.0]),
            np.array([0.0]),
            np.array([1.0]),
            anchor=np.array([False]),
        )
        is None
    )
    assert (
        fit_landscape(np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), anchor=np.zeros(0, bool))
        is None
    )


def test_the_landscape_keeps_its_own_density_beside_the_padded_one():
    """``logP_free`` is the rendered population before the one pseudo-region of ignorance, and ``logP`` is that
    density plus the uniform pseudo-mass and nothing else: one affine map takes the first to the second exactly."""
    count, mass, e, var = _two_mode()
    ls = fit_landscape(count, mass, e, var, anchor=np.zeros(count.size, bool))
    free, padded = np.exp(ls.logP_free), np.exp(ls.logP)
    assert free.sum() == pytest.approx(1.0, abs=1e-12)
    slope, intercept = np.polyfit(free, padded, 1)
    np.testing.assert_allclose(padded, slope * free + intercept, rtol=0, atol=1e-15)
    assert 0.0 < slope < 1.0 and intercept > 0.0


# ── the census: the bulk and the enriched basins ───────────────────────────────────────────────────────


def test_split_basins_names_the_largest_basin_depleted_and_returns_only_the_basins_above_it():
    """Depleted is the largest-mass basin; the basins returned beside it are those strictly above it, so a
    smaller basin BELOW the depleted one is never an enriched candidate, and a density whose largest basin
    is the top one has nothing above it."""
    below = LandscapeMode(log_rho=-7.0, basin_mass=0.10, lo=-9.0, hi=-6.0)
    base = LandscapeMode(log_rho=-4.0, basin_mass=0.60, lo=-6.0, hi=-2.0)
    mid = LandscapeMode(log_rho=0.0, basin_mass=0.30, lo=-2.0, hi=1.5)
    dep, above = split_basins((below, base, mid))
    assert dep is base and above == (mid,)
    small = LandscapeMode(log_rho=-4.0, basin_mass=0.20, lo=-6.0, hi=-2.0)
    top = LandscapeMode(log_rho=0.0, basin_mass=0.70, lo=-2.0, hi=1.5)
    dep2, above2 = split_basins((below, small, top))
    assert dep2 is top and above2 == ()


def _drawn(p) -> DensityLandscape:
    """A landscape whose density is ``p`` itself, on a unit grid."""
    p = np.asarray(p, dtype=np.float64)
    with np.errstate(divide="ignore"):
        lp = np.log(p / p.sum())
    return DensityLandscape(
        log_rho=np.arange(p.size, dtype=np.float64), logP=lp, n_train=p.size, logP_free=lp
    )


def test_the_basins_partition_the_grid():
    """Masses sum to one, adjacent basins share their bound and the basins span the grid, whatever the shape."""
    for p in (
        [0.70, 0.10, 0.05, 0.02, 0.03, 0.08, 0.02],
        [1, 3, 1, 2, 5, 2, 1, 4, 9],
        [5, 4, 3, 2, 1],
    ):
        modes = census(_drawn(p))
        assert sum(m.basin_mass for m in modes) == pytest.approx(1.0, rel=0, abs=1e-12)
        assert all(a.hi == b.lo for a, b in zip(modes, modes[1:]))
        assert modes[0].lo == 0.0 and modes[-1].hi == len(p) - 1


def test_a_density_highest_at_the_grid_edge_peaks_there():
    """The review's case: a floor of mass 0.70 at the grid's low edge and a small bump at x = 5. The edge is a
    maximum, so it is the bulk, and the bump is the enriched basin above it. PERTURBATION: reading interior
    maxima only makes the bump the one basin, holding everything, and the floor-versus-enriched contrast is
    gone."""
    modes = census(_drawn([0.70, 0.10, 0.05, 0.02, 0.03, 0.08, 0.02]))
    assert [m.log_rho for m in modes] == [0.0, 5.0]
    bulk, above = split_basins(modes)
    assert bulk.log_rho == 0.0 and [m.log_rho for m in above] == [5.0]
    assert bulk.basin_mass == pytest.approx(0.85)


def test_the_high_edge_peaks_too_and_a_flat_edge_does_not():
    """A density rising to the grid's top edge peaks there; an edge level with its neighbour is a plateau's end,
    no peak; a monotone density peaks once, at its arg-max."""
    assert [m.log_rho for m in census(_drawn([3, 1, 2, 5]))] == [0.0, 3.0]
    assert [m.log_rho for m in census(_drawn([1, 1, 2, 5, 3]))] == [3.0]
    assert [m.log_rho for m in census(_drawn([5, 4, 3, 2, 1]))] == [0.0]
    assert [m.log_rho for m in census(_drawn([1, 2, 3, 4, 5]))] == [4.0]


def test_a_plateau_peaks_at_its_left_end():
    assert [m.log_rho for m in census(_drawn([1, 3, 3, 3, 1]))] == [1.0]


def _walls_around_a_tail(seed=0, n_long=20_000, n_dep=200, n_short=600, n_tail=12):
    """A sparse human-like library, refitted through the loop's own E-step: the depleted bulk is 20,000
    long anchors (walls at 10^-5) with 200 located intron pieces at the depleted level; 600 SHORT anchors —
    20–60 bp intergenic slivers that sequenced nothing — whose resolution walls 1/E sit between 10^-1.8
    and 10^-1.3; and twelve located exon kernels near 10^-1.3, a tail and not a population."""
    rng = np.random.default_rng(seed)
    eff = np.concatenate(
        [
            np.full(n_long, 1e5),
            np.full(n_dep, 1e5),
            rng.uniform(20.0, 60.0, n_short),
            np.full(n_tail, 300.0),
        ]
    )
    count = np.concatenate(
        [
            np.zeros(n_long),
            rng.poisson(3e-5 * 1e5, n_dep) + 1.0,
            np.zeros(n_short),
            rng.poisson(0.05 * 300.0, n_tail) + 1.0,
        ]
    )
    mass = np.maximum(count * 2.0, 1.0)
    var = np.where(count > 0, 0.1, np.inf)
    anchor = count <= 0.0
    ls = None
    for _ in range(
        3
    ):  # calib_refit_iters: the E-step moves the walls' mass to the bulk, not their centres
        ls = fit_landscape(count, mass, eff, var, anchor=anchor, prev=ls)
    return ls


def test_walls_around_a_handful_of_exons_lift_nothing_and_the_exons_read_their_own_level():
    """cfRNA-A's shape, through the estimator itself: the walls' mass goes to the bulk in the E-step, and the
    twelve exons are a basin of a few parts in ten thousand. Read by the capture reader, nothing without evidence
    is lifted — a 40 bp sliver with no fragment reads the bulk to a part in a thousand, a 100 kb intron with none
    reads it exactly — and an exon at the tail's level, fifteen fragments on 300 bp, reads its own level: the
    small panel the located-mode gate would not read at all."""
    ls = _walls_around_a_tail()
    assert ls is not None
    _bulk, above = split_basins(census(ls))
    assert len(above) == 1 and above[0].basin_mass < 1e-3
    r = capture_efficiencies(
        ls, np.array([0.0, 0.0, 15.0]), np.array([40.0, 1e5, 300.0]), np.zeros(0), np.zeros(0)
    )
    assert abs(r.efficiency_region[0] - 1.0) < 1e-3
    assert abs(r.efficiency_region[1] - 1.0) < 1e-9
    assert abs(np.log(r.efficiency_region[2] * r.bulk_density / 0.05)) < 0.2
