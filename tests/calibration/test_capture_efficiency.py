"""``capture_efficiency``: an object's capture efficiency is its own gDNA density read against the bulk — the posterior
mean of its LOG density, saturated per enriched basin, under the fitted landscape as a two-group prior (the bulk basin
a point at its peak, the population's own density elsewhere), from its OWN count: a region's contained count on its
contained support, a boundary's crossing count on its crossing support. No mode, no gate, no constant.

The falsification test was written first and verified failing on the shipped ruler through its own fixture
(`ISSUES: ruler-multimapper-floor-caps-the-correction`): an unprobed exon holding no gDNA fragment read the
multimapper floor, 24× the depleted level. Here it reads the bulk. The rest pins the model's properties: an informative
object reads its own level, capture saturates per level, a one-basin landscape reads exactly 1, a lone level lifts an
object without evidence by at most ``e^{μΔ}``, a level below the bulk reads below 1, each object reads only its own
count, the density unit is free, the row tiles change nothing, and calibrate's wiring. A piece too short to contain a
fragment is read through its boundaries by the length, not here (``tests/calibration/test_capture_eff_length.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

import rigel.calibration.capture_efficiency as CE
from rigel.calibration.capture_efficiency import capture_efficiencies
from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_eff_length
from rigel.calibration.landscape import DensityLandscape, census, fit_landscape, split_basins

RHO_CAP, RHO_DEP = 1.0, 1e-3


def _pmf(mean=200.0, sd=50.0, lo=100, hi=300) -> np.ndarray:
    w = np.arange(hi + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((w - mean) / sd) ** 2)
    p[(w < lo) | (w > hi)] = 0.0
    return p / p.sum()


PMF = _pmf()
S_E = float(crossing_eff_length(PMF, np.array([UNBOUNDED_REACH]), np.array([UNBOUNDED_REACH]))[0])


def _support(lengths, pmf=PMF) -> np.ndarray:
    """The contained support ``E_f[(ℓ − w + 1)⁺]`` per piece — 0 for a piece shorter than every fragment."""
    w = np.arange(pmf.shape[0], dtype=np.float64)
    L = np.asarray(lengths, dtype=np.float64)[:, None]
    return (pmf[None, :] * np.maximum(L - w[None, :] + 1.0, 0.0)).sum(1)


def _fitted(rng_seed: int, n_dep: int, n_cap: int, n_anchor: int) -> DensityLandscape:
    rng = np.random.default_rng(rng_seed)
    n = n_dep + n_cap + n_anchor
    eff = np.full(n, 1000.0)
    count = np.concatenate(
        [
            rng.poisson(RHO_DEP * 1000.0, n_dep).astype(float),
            rng.poisson(RHO_CAP * 1000.0, n_cap).astype(float),
            np.zeros(n_anchor),
        ]
    )
    anchor = np.r_[np.zeros(n_dep + n_cap, bool), np.ones(n_anchor, bool)]
    ls = fit_landscape(count, np.maximum(count * 2.0, 1.0), eff, np.full(n, 0.1), anchor=anchor)
    assert ls is not None
    return ls


@pytest.fixture(scope="module")
def landscape():
    """A bimodal gDNA landscape fitted on a population: 900 depleted pieces at 10^-3 and 300 captured pieces at
    1, each a Poisson count on a 1 kb support, plus 200 zero-count anchors."""
    return _fitted(3, 900, 300, 200)


def _hand(levels, masses, widths=0.3, lo=1e-7, hi=100.0, n=300) -> DensityLandscape:
    """A landscape drawn by hand on a natural-log grid: Gaussian bumps in log density at ``levels`` carrying
    ``masses``; ``logP`` adds a uniform pseudo-mass as the fit does, ``logP_free`` is the bumps alone."""
    x = np.linspace(np.log(lo), np.log(hi), n)
    p = np.zeros(n)
    for level, m, w in zip(levels, masses, np.broadcast_to(widths, len(levels)), strict=True):
        b = np.exp(-0.5 * ((x - np.log(level)) / w) ** 2)
        p += m * b / b.sum()
    padded = p + 1.0 / (n * 1000.0)
    with np.errstate(divide="ignore"):
        free = np.log(p / p.sum())
    return DensityLandscape(
        log_rho=x, logP=np.log(padded / padded.sum()), n_train=1000, logP_free=free
    )


def _field(lengths, density, crossing_density):
    """Contained counts ``ρ_r S_r`` per region and crossing counts ``ρ_e S_E`` per boundary — each object's
    own count at its own density, which is all an efficiency reads."""
    S = _support(lengths)
    k_reg = np.asarray(density, dtype=np.float64) * S
    k_bnd = np.asarray(crossing_density, dtype=np.float64) * S_E
    return k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E)


def _dense_reference(ls: DensityLandscape, k, S) -> np.ndarray:
    """The reader written out once with no row tiles and its bulk term in closed form — the specification
    the tiled implementation must equal."""
    x = np.asarray(ls.log_rho, dtype=np.float64)
    bulk, above = split_basins(census(ls))
    xb = float(bulk.log_rho)
    in_bulk = (x >= bulk.lo) & (x <= bulk.hi)
    p = np.exp(ls.logP_free - np.max(ls.logP_free))
    p /= p.sum()
    pi0 = float(p[in_bulk].sum())
    alt = np.where(in_bulk, 0.0, p)
    tiny = np.finfo(float).tiny
    log_alt = np.where(alt > 0.0, np.log(np.maximum(alt, tiny)), -np.inf)
    cap = np.full(x.shape, np.inf)
    for m in above:
        cap[(x >= m.lo) & (x <= m.hi)] = m.log_rho
    floor = CE.CONSTANTS.calibration.rate_log_floor
    lam = np.exp(x)[None, :] * S[:, None]
    ll = np.where(S[:, None] > 0.0, k[:, None] * np.log(np.maximum(lam, floor)) - lam, 0.0)
    lam_b = np.exp(xb) * S
    ll_b = np.where(S > 0.0, k * np.log(np.maximum(lam_b, floor)) - lam_b, 0.0)
    a_null = np.log(pi0) + ll_b
    a_alt = log_alt[None, :] + ll
    top = np.maximum(a_null, a_alt.max(1))
    w_null, w_alt = np.exp(a_null - top), np.exp(a_alt - top[:, None])
    den = w_null + w_alt.sum(1)
    mean_level = (w_null * xb + w_alt @ np.minimum(x, cap)) / den
    return np.where(S > 0.0, np.exp(mean_level - xb), 1.0)


# --- the falsification test ---------------------------------------------------------------------------


def test_an_unprobed_exon_with_no_gdna_fragment_reads_the_bulk_not_a_floor(landscape):
    """A 1 kb exon holding no gDNA at all among depleted neighbours: its posterior sits at the bulk — within a
    factor of two of 1 — and nowhere near the multimapper floor's 24× the bulk."""
    k_reg, S, k_bnd, S_b = _field([5000, 1000, 5000], [RHO_DEP] * 3, [RHO_DEP] * 2)
    k_reg[1] = 0.0
    c = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b).efficiency_region
    assert abs(np.log(c[1])) < np.log(2.0), c[1]
    assert c[1] < 5.0


# --- the model's properties ---------------------------------------------------------------------------


def test_an_informative_object_reads_its_own_level(landscape):
    """At high depth the posterior is the object's own level against the bulk: 50 kb pieces at the captured level
    read the captured level over the bulk, at the depleted level 1, and at a third of the captured level — inside
    the enriched basin, below its peak — a third of it."""
    r = capture_efficiencies(
        landscape, *_field([50000] * 3, [RHO_CAP, RHO_DEP, RHO_CAP / 3.0], [RHO_DEP] * 2)
    )
    c, rb = r.efficiency_region, r.bulk_density
    assert abs(np.log(c[1])) < 0.2
    assert abs(np.log(c[0] * rb / RHO_CAP)) < 0.2
    assert abs(np.log(c[2] * rb * 3.0 / RHO_CAP)) < 0.2


def test_capture_saturates_at_its_own_basins_peak(landscape):
    """A density above its enriched basin's peak reads that peak: a 50 kb piece at three times the captured level
    reads exactly the saturated level, the most massive (here the only) enriched basin's peak over the bulk.
    PERTURBATION: with no saturation it reads its own level, three times higher."""
    r = capture_efficiencies(landscape, *_field([50000], [3.0 * RHO_CAP], []))
    assert r.saturation > 100.0
    assert r.efficiency_region[0] == pytest.approx(r.saturation, rel=1e-9)


def test_a_one_basin_landscape_reads_exactly_one_everywhere():
    """A library with no capture is one gDNA level: the prior is then a point, so every object reads exactly 1 —
    shallow or deep, at the bulk or far from it — and the saturated level is 1. No decision is taken."""
    ls = _fitted(4, 1400, 0, 200)
    _bulk, above = split_basins(census(ls))
    assert above == ()
    k_reg, S, k_bnd, S_b = _field(
        [50000, 5000, 300, 40], [RHO_DEP, RHO_CAP, 0.0, RHO_CAP], [RHO_CAP] * 3
    )
    r = capture_efficiencies(ls, k_reg, S, k_bnd, S_b)
    np.testing.assert_array_equal(r.efficiency_region, 1.0)
    np.testing.assert_array_equal(r.efficiency_boundary, 1.0)
    assert r.saturation == 1.0


def test_a_lone_level_lifts_an_object_without_evidence_by_at_most_e_to_the_mu_delta():
    """An object whose own count says nothing reads the prior. In the log readout a lone level of prior mass ``μ``
    at ``Δ`` nats above the bulk lifts it by at most ``e^{μΔ}``. PERTURBATION: the posterior mean of the density
    itself reads ``1 + μ(e^Δ − 1)`` — eleven times the bulk here, the zero-gDNA catastrophe of a linear reader."""
    mu, delta = 0.01, np.log(1000.0)
    ls = _hand([1e-3, 1.0], [1.0 - mu, mu], widths=0.15)
    r = capture_efficiencies(ls, np.zeros(1), np.full(1, 1e-12), np.zeros(0), np.zeros(0))
    e = float(r.efficiency_region[0])
    assert 1.0 <= e <= np.exp(mu * delta) * 1.01
    assert e < 0.2 * (1.0 + mu * (np.exp(delta) - 1.0))


def test_a_level_below_the_bulk_reads_below_one():
    """The bulk is the largest basin, not the lowest: a well-measured object at a level below it reads below 1,
    as its own count says."""
    ls = _hand([1e-5, 1e-3, 1.0], [0.05, 0.85, 0.10])
    S = np.full(1, 1e7)
    r = capture_efficiencies(ls, 1e-5 * S, S, np.zeros(0), np.zeros(0))
    assert abs(np.log(r.efficiency_region[0] / 1e-2)) < 0.2


@pytest.mark.parametrize("captured_mass", [0.4, 0.6])
def test_a_well_measured_capture_ratio_survives_the_captured_level_becoming_the_majority(
    captured_mass,
):
    """THE COUNTEREXAMPLE that refused the one-sided null (the capture-reference review, 2026-10-13). Two levels, 0.001
    and 0.1 fragments per bp, read by objects with a million positions each: the truth ratio is exactly 100. When the
    captured level holds 60% of the population it is the largest basin and becomes the unit, so the uncaptured object
    reads 0.01 and the captured one 1 — the ratio is still 100. PERTURBATION: absorbing every level at or below the
    bulk into the null (the one-sided prior) reads both at 1 and erases a ratio no amount of depth can restore."""
    x = np.log(np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0]))
    p = np.array([0.0, 1.0 - captured_mass, 0.0, captured_mass, 0.0])
    padded = p + 1.0 / (1000 * x.size)
    with np.errstate(divide="ignore"):
        ls = DensityLandscape(
            log_rho=x, logP=np.log(padded / padded.sum()), n_train=1000, logP_free=np.log(p)
        )
    S = np.array([1e6, 1e6])
    e = capture_efficiencies(
        ls, np.array([1e3, 1e5]), S, np.zeros(0), np.zeros(0)
    ).efficiency_region
    assert e[1] / e[0] == pytest.approx(100.0, rel=1e-9)
    assert (
        (e[0] == pytest.approx(1.0, rel=1e-9))
        if captured_mass < 0.5
        else (e[1] == pytest.approx(1.0, rel=1e-9))
    )


def test_a_boundarys_efficiency_is_its_own_crossing_counts_posterior(landscape):
    """A boundary reads its own crossing count on its crossing support, by the same posterior a region reads its
    contained count: crossings at the captured level read the captured level, at the depleted level 1, and every
    boundary reads exactly what a region holding its count on its support would."""
    k_reg, S, k_bnd, S_b = _field([5000] * 4, [RHO_CAP] * 4, [RHO_CAP, RHO_DEP, 0.5 * RHO_CAP])
    r = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b)
    cb, rb = r.efficiency_boundary, r.bulk_density
    assert abs(np.log(cb[0] * rb / RHO_CAP)) < 0.3
    assert abs(np.log(cb[1])) < 0.5
    as_regions = capture_efficiencies(
        landscape, k_bnd, S_b, np.zeros(0), np.zeros(0)
    ).efficiency_region
    np.testing.assert_array_equal(cb, as_regions)


def test_each_object_reads_only_its_own_count(landscape):
    """A region's efficiency is its contained count's and a boundary's its crossing count's: a depleted intron
    beside captured crossings stays at the bulk, and changing every boundary's count leaves every region's
    efficiency bit-identical (and the reverse) — the crossings are priced once, at the boundary that holds them.
    PERTURBATION: apportioning the crossings onto the pieces within reach lifts the intron and moves the regions
    with the boundaries."""
    lengths = [5000, 1000, 40, 1000, 5000]
    density = [RHO_DEP, RHO_DEP, RHO_CAP, RHO_DEP, RHO_DEP]
    k_reg, S, k_bnd, S_b = _field(
        lengths, density, [RHO_DEP, 0.5 * RHO_CAP, 0.5 * RHO_CAP, RHO_DEP]
    )
    r = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b)
    assert r.efficiency_region[1] < 2.0 and r.efficiency_region[3] < 2.0
    r2 = capture_efficiencies(landscape, k_reg, S, np.full(4, RHO_CAP) * S_E, S_b)
    np.testing.assert_array_equal(r2.efficiency_region, r.efficiency_region)
    assert not np.array_equal(r2.efficiency_boundary, r.efficiency_boundary)
    r3 = capture_efficiencies(landscape, np.full(5, RHO_CAP) * S, S, k_bnd, S_b)
    np.testing.assert_array_equal(r3.efficiency_boundary, r.efficiency_boundary)
    assert not np.array_equal(r3.efficiency_region, r.efficiency_region)


def test_an_object_with_no_support_reads_the_null(landscape):
    """A piece shorter than every fragment has no gDNA support and no likelihood: no gDNA fragment could sit there
    to witness its capture, so it reads exactly 1 — the null — whatever its neighbours hold."""
    lengths = [5000, 40, 5000]
    k1, S, kb, S_b = _field(lengths, [RHO_CAP] * 3, [RHO_CAP] * 2)
    k2, _, kb2, _ = _field(lengths, [RHO_DEP] * 3, [RHO_DEP] * 2)
    assert S[1] == 0.0
    c1 = capture_efficiencies(landscape, k1, S, kb, S_b).efficiency_region
    c2 = capture_efficiencies(landscape, k2, S, kb2, S_b).efficiency_region
    assert c1[1] == 1.0 and c2[1] == 1.0
    assert c1[0] > 100.0 and c2[0] < 2.0


def test_the_density_unit_is_free():
    """Rescaling the grid's densities by ``c`` and every support by ``1/c`` leaves every rate, and so every
    efficiency, unchanged; the bulk density scales by ``c``."""
    ls = _hand([1e-5, 1e-3, 0.2, 1.0], [0.03, 0.80, 0.05, 0.12])
    rng = np.random.default_rng(9)
    S = rng.uniform(10.0, 1e5, 400)
    k = np.floor(np.exp(rng.uniform(np.log(1e-5), np.log(3.0), 400)) * S)
    c = 37.0
    shifted = DensityLandscape(
        log_rho=ls.log_rho + np.log(c), logP=ls.logP, n_train=ls.n_train, logP_free=ls.logP_free
    )
    a = capture_efficiencies(ls, k, S, k[:50], S[:50])
    b = capture_efficiencies(shifted, k, S / c, k[:50], S[:50] / c)
    np.testing.assert_allclose(b.efficiency_region, a.efficiency_region, rtol=1e-9)
    np.testing.assert_allclose(b.efficiency_boundary, a.efficiency_boundary, rtol=1e-9)
    assert b.bulk_density == pytest.approx(c * a.bulk_density, rel=1e-12)
    assert b.saturation == pytest.approx(a.saturation, rel=1e-12)


def test_the_row_tiles_change_nothing_and_equal_the_dense_reader(monkeypatch):
    """The objects are read in row tiles so a million never exist against the grid at once: the tile size moves no
    efficiency beyond the last bits of the matrix product's summation order, and the tiled reader equals the reader
    written out densely with its bulk term in closed form."""
    ls = _hand([1e-5, 1e-3, 0.03, 0.2, 1.0], [0.02, 0.80, 0.03, 0.03, 0.12])
    rng = np.random.default_rng(12)
    S = np.where(rng.random(1000) < 0.05, 0.0, rng.uniform(1.0, 1e5, 1000))
    k = np.floor(np.exp(rng.uniform(np.log(1e-6), np.log(5.0), 1000)) * S)
    tiled = capture_efficiencies(ls, k, S, k[:10], S[:10]).efficiency_region
    np.testing.assert_allclose(tiled, _dense_reference(ls, k, S), rtol=1e-12)
    for rows in (1, 7, 333):
        monkeypatch.setattr(CE, "_block_rows", lambda n_grid, n_arrays, rows=rows: rows)
        np.testing.assert_allclose(
            capture_efficiencies(ls, k, S, k[:10], S[:10]).efficiency_region, tiled, rtol=1e-13
        )


def test_the_saturation_level_is_the_most_massive_enriched_basins_peak():
    """Of two enriched basins the most massive sets the saturated level, not the highest: the junction's bound is
    the level capture fills, and a lighter basin above it is a stronger probe's own level."""
    ls = _hand([1e-3, 0.3, 3.0], [0.80, 0.15, 0.05])
    bulk, above = split_basins(census(ls))
    heavy = max(above, key=lambda m: m.basin_mass)
    r = capture_efficiencies(ls, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0))
    assert len(above) == 2 and heavy.log_rho < max(m.log_rho for m in above)
    assert r.saturation == pytest.approx(np.exp(heavy.log_rho - bulk.log_rho), rel=1e-12)
    assert r.efficiency_region.shape == (0,) and r.efficiency_boundary.shape == (0,)


@pytest.mark.parametrize(
    "bad",
    [
        dict(count_region=np.array([np.nan, 1.0])),
        dict(eff_region=np.array([np.inf, 1.0])),
        dict(count_boundary=np.array([1.0, -np.inf])),
        dict(eff_region=np.array([1.0, 2.0, 3.0])),
    ],
)
def test_a_non_finite_or_misshapen_input_is_refused(landscape, bad):
    kw = dict(
        count_region=np.array([1.0, 2.0]),
        eff_region=np.array([100.0, 200.0]),
        count_boundary=np.array([0.0, 1.0]),
        eff_boundary=np.array([50.0, 50.0]),
    )
    kw.update(bad)
    with pytest.raises(ValueError, match="capture_efficiencies"):
        capture_efficiencies(landscape, **kw)


def test_calibrate_publishes_exactly_one_and_no_bulk_without_a_landscape():
    """The wiring: a solve with no refit fits no landscape, so the result carries no bulk density, a saturated
    level of 1 and efficiencies of exactly 1 on both axes; and it carries gDNA's conserved share at every boundary
    — its two flanking regions' shares on the gDNA pmf, at the unbounded reach of a template that does not end
    (`effective_length.conserved_cut_shares`)."""
    import sys

    sys.path.insert(0, "tests/calibration")
    from _synthetic import (
        make_gdna_fl_pmf,
        make_strand_models,
        make_synthetic_payload,
        make_synthetic_sj,
    )

    from rigel.calibration import calibrate
    from rigel.config import CalibrationConfig

    payload, ra = make_synthetic_payload()
    # gDNA fragments longer than the fixture's 100 bp regions, so one crosses both boundaries and its
    # conserved share differs from the crossing support
    pmf = np.zeros(251)
    pmf[150:251] = 1.0 / 101.0
    res = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=make_strand_models(0.95, 40),
        gdna_fl_pmf=pmf,
        rna_fl_pmf=make_gdna_fl_pmf(),
        config=CalibrationConfig(calib_refit_iters=0),
        sj=make_synthetic_sj(),
    )
    assert res.gdna_bulk_density is None
    assert res.gdna_capture_saturation == 1.0
    np.testing.assert_array_equal(res.gdna_capture_efficiency_region, 1.0)
    np.testing.assert_array_equal(res.gdna_capture_efficiency_boundary, 1.0)
    from rigel.calibration.effective_length import conserved_cut_shares
    from rigel.calibration.region_arrays import boundary_region_indices

    lo, hi = boundary_region_indices(np.asarray(ra.ref_id))
    length = np.asarray(ra.end, dtype=np.float64) - np.asarray(ra.start, dtype=np.float64)
    left, right = conserved_cut_shares(
        pmf, length[lo], length[hi], UNBOUNDED_REACH, UNBOUNDED_REACH
    )
    assert lo.size > 0
    np.testing.assert_array_equal(res.gdna_boundary_conserved_len, left + right)
    assert np.all(res.gdna_boundary_conserved_len < res.gdna_boundary_eff_len - 1.0)
