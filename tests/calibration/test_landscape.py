"""The population gDNA-density hyperprior (`calibration.landscape`), one property per test.

What the fit has to survive is a library whose truth is a uniform depleted level plus a
capture-enriched minority: the minority must stay a mode rather than being smoothed into the bulk,
the grid must span exactly what ψ can represent and nothing else, a zero-count region must decay
downward instead of inventing a location at 1/E, an imprecise region must be damped rather than
dropped, and the kernel widths must come from the local neighbour spacing so they widen on their own
as regions thin out. The mode census that reads the ruler's reference off the fit closes the file.
Each gate holds one property, so a regression names itself.
"""

import numpy as np
import pytest

from rigel.calibration.landscape import (
    DensityLandscape,
    LandscapeMode,
    _census,
    fit_landscape,
    knn_widths,
    located_enriched_mode,
    split_basins,
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
    a = fit_landscape(*args, anchor=anchor)
    b = fit_landscape(*args, anchor=anchor)
    assert np.array_equal(a.logP, b.logP) and np.array_equal(a.log_rho, b.log_rho)


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


# ── the mode census: the ruler's reference ─────────────────────────────────────────────────────────


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


def _hand_landscape(peaks, sds, masses, members, walls=(), spread=0.1, lo=-16.0, hi=3.0, n=400):
    """A DensityLandscape rendered by hand — a mixture of Gaussians in natural-log density — with its
    kernels stated outright: per peak, ``members`` LOCATED kernels spread evenly across ±``spread`` nat
    of it, and ``walls`` as ``(centre, n)`` pairs of location-free kernels (anchors' resolution walls) —
    so a test controls where each basin sits, how much it holds, and how many located kernels stand
    behind it."""
    x = np.linspace(lo, hi, n)
    p = np.zeros_like(x)
    for mu, sd, m in zip(peaks, sds, masses, strict=True):
        p += m * np.exp(-0.5 * ((x - mu) / sd) ** 2) / sd
    p /= p.sum()
    centre = [
        np.linspace(mu - spread, mu + spread, k) for mu, k in zip(peaks, members, strict=True)
    ]
    located = [np.ones(k, dtype=bool) for k in members]
    for mu, k in walls:
        centre.append(np.linspace(mu - 0.1, mu + 0.1, k))
        located.append(np.zeros(k, dtype=bool))
    return DensityLandscape(
        log_rho=x,
        logP=np.log(p + 1e-300),
        n_train=int(sum(members) + sum(k for _, k in walls)),
        centre=np.concatenate(centre),
        located=np.concatenate(located),
    )


def test_a_unimodal_landscape_names_NO_enriched_mode():
    """Capture-OFF: one basin, nothing above the depleted mode, no reference, no contraction."""
    assert located_enriched_mode(_hand_landscape([-3.0], [0.3], [1000], [400])) is None


def test_a_bimodal_landscape_names_the_LOCATED_upper_mode_at_its_peak_with_its_members():
    """Capture-ON: the basin above the depleted one holding the most located kernels, resolved by them
    (60 members against k = √300 ≈ 17); the reference is its peak and the members are published."""
    ref = located_enriched_mode(_hand_landscape([-7.0, 0.5], [0.4, 0.25], [700, 300], [240, 60]))
    assert ref is not None and abs(ref.mode.log_rho - 0.5) < 0.05
    assert ref.n_members == 60


def test_a_basin_of_fewer_than_k_located_kernels_is_not_a_mode_however_narrow():
    """The blank contig's shadow transcription, a few false-positive fragments on slivers of support, or a
    sparse library's tail of a dozen measured exons: a basin whose located members number at most
    k = √n_located has no k-th neighbour inside itself and is no mode, whatever the density's cut made
    of it. PERTURBATION: the same landscape with the cluster grown past k IS a mode."""
    # 400 located kernels in all ⇒ k = 20; 20 members do not resolve the basin, 21 do
    assert (
        located_enriched_mode(_hand_landscape([-14.0, -5.0], [0.5, 0.2], [998, 2], [380, 20]))
        is None
    )
    grown = located_enriched_mode(_hand_landscape([-14.0, -5.0], [0.5, 0.2], [998, 2], [379, 21]))
    assert grown is not None and abs(grown.mode.log_rho - (-5.0)) < 0.05 and grown.n_members == 21
    # and the members are read at the POPULATION's k, not their own √21: the same 21 strewn across
    # ±1.5 nat put every member's 20th neighbour ~3 nat away, and that is no location either
    strewn = _hand_landscape([-14.0, -5.0], [0.5, 2.0], [998, 2], [379, 21], spread=1.5)
    assert located_enriched_mode(strewn) is None


def test_anchors_walls_are_not_members_and_cannot_locate_a_basin():
    """`ISSUES: the-ruler-reference-on-sparse-real-libraries`, the reader's half: a basin above the bulk
    packed with 20,000 anchors' resolution walls around twelve located kernels names no reference — the
    walls are not members, and twelve is below k = √n_located. PERTURBATION: admitting every centre as a
    member, walls included, reads the walls as a located mode of 20,012 members."""
    ls = _hand_landscape([-7.0, -1.5], [0.4, 0.2], [700, 300], [1000, 12], walls=[(-1.5, 20_000)])
    assert located_enriched_mode(ls) is None
    # and with the same twelve grown into a population the walls change nothing either way
    ls2 = _hand_landscape([-7.0, -1.5], [0.4, 0.2], [700, 300], [1000, 200], walls=[(-1.5, 20_000)])
    ref = located_enriched_mode(ls2)
    assert ref is not None and ref.n_members == 200


def test_the_enriched_mode_is_chosen_by_LOCATED_MEMBERS_above_the_depleted_one():
    """Two basins above the depleted one: the reference is the one holding the located population, not
    the one holding the most rendered mass — here the upper basin carries more density (walls render
    mass too) but thirty located kernels against the middle basin's 250, so the middle one is the
    reference. PERTURBATION: chosen by rendered mass, the upper basin is the candidate, and with thirty
    members against k = 31 it is no mode at all. A basin with no located member is nothing."""
    ls = _hand_landscape([-7.0, -1.0, 0.8], [0.4, 0.3, 0.1], [700, 100, 250], [700, 250, 30])
    ref = located_enriched_mode(ls)
    assert ref is not None and abs(ref.mode.log_rho - (-1.0)) < 0.05 and ref.n_members == 250
    empty = _hand_landscape([-7.0, 0.5], [0.4, 0.25], [700, 300], [700, 0])
    assert located_enriched_mode(empty) is None
    modes = _census(ls)  # the basins partition the grid: masses sum to one, adjacent bounds shared
    assert sum(m.basin_mass for m in modes) == pytest.approx(1.0, rel=0, abs=1e-9)
    assert all(a.hi == b.lo for a, b in zip(modes, modes[1:]))


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


def test_a_fitted_basin_of_walls_around_a_handful_of_located_kernels_is_NOT_a_mode():
    """LBX0190's defect through the estimator itself (the falsification test, verified failing on the
    shipped reader, which read a located mode of 612 members here): the walls' centres pack the basin
    above the bulk so densely that every rendered width there sits at the grid step. A wall is not a
    location: the basin's located members are twelve, fewer than k = √212, and it names no reference."""
    ls = _walls_around_a_tail()
    assert ls is not None
    assert int(ls.located.sum()) == 212
    assert located_enriched_mode(ls) is None
