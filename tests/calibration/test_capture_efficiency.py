"""``capture_efficiency``: an object's capture efficiency is its own gDNA density read against the bulk.

It is the posterior mean of the object's log density under the fitted landscape as a two-group prior (the bulk basin a
point at its peak, the population's own density elsewhere), from the object's own count: a region's contained count on
its contained support, a boundary's crossing count on its crossing support. No mode, no gate, no constant. An unprobed
exon holding no gDNA fragment reads the bulk, not a multimapper floor. The tests pin the model's properties: an
informative object reads its own level, above its basin's peak as below it, a one-basin landscape reads exactly 1, a
lone level lifts an object without evidence by at most ``e^{μΔ}``, a level below the bulk reads below 1, each object
reads only its own count, the density unit is free, the row tiles change nothing, and calibrate's wiring. A piece too short to contain a fragment is read through its boundaries by the length
(``test_capture_eff_length.py``), not here.

The published reading is the reader's raised to one library-wide weight, the share of the one-strand unspliced mass
whose object would show a gDNA fragment of its own (``visible_share``, ``weighted``). The tests pin the share's form
(mass-weighted, the Poisson chance, both-strand objects silent), the power on every level alike, the hole it closes
(the reads that sit where no gDNA shows take the claim with them) and the locality it keeps (one probe in a library
whose gDNA shows keeps its reading).
"""

from __future__ import annotations

import numpy as np
import pytest
from _synthetic import contained_support, normal_pmf, spike_pmf

import rigel.calibration.capture_efficiency as CE
from rigel.calibration.capture_efficiency import capture_efficiencies, visible_share, weighted
from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_eff_length
from rigel.calibration.landscape import DensityLandscape, census, fit_landscape, split_basins

RHO_CAP, RHO_DEP = 1.0, 1e-3


PMF = normal_pmf()
S_E = float(crossing_eff_length(PMF, np.array([UNBOUNDED_REACH]), np.array([UNBOUNDED_REACH]))[0])


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
    """Contained counts ``ρ_r S_r`` per region and crossing counts ``ρ_e S_E`` per boundary: each object's own count at
    its own density, all an efficiency reads."""
    S = contained_support(lengths, PMF)
    k_reg = np.asarray(density, dtype=np.float64) * S
    k_bnd = np.asarray(crossing_density, dtype=np.float64) * S_E
    return k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E)


def _dense_reference(ls: DensityLandscape, k, S) -> np.ndarray:
    """The reader written out once with no row tiles and its bulk term in closed form: the specification the tiled
    implementation must equal."""
    x = np.asarray(ls.log_rho, dtype=np.float64)
    bulk, _above = split_basins(census(ls))
    xb = float(bulk.log_rho)
    in_bulk = (x >= bulk.lo) & (x <= bulk.hi)
    p = np.exp(ls.logP_free - np.max(ls.logP_free))
    p /= p.sum()
    pi0 = float(p[in_bulk].sum())
    alt = np.where(in_bulk, 0.0, p)
    tiny = np.finfo(float).tiny
    log_alt = np.where(alt > 0.0, np.log(np.maximum(alt, tiny)), -np.inf)
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
    mean_level = (w_null * xb + w_alt @ x) / den
    return np.where(S > 0.0, np.exp(mean_level - xb), 1.0)


# --- an unprobed exon with no gDNA ---------------------------------------------------------------------------


def test_an_unprobed_exon_with_no_gdna_fragment_reads_the_bulk_not_a_floor(landscape):
    """A 1 kb exon holding no gDNA among depleted neighbours: its posterior sits at the bulk, within a factor of two of
    1, far from a multimapper floor's 24× the bulk."""
    k_reg, S, k_bnd, S_b = _field([5000, 1000, 5000], [RHO_DEP] * 3, [RHO_DEP] * 2)
    k_reg[1] = 0.0
    c = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b).efficiency_region
    assert abs(np.log(c[1])) < np.log(2.0), c[1]
    assert c[1] < 5.0


# --- the model's properties ---------------------------------------------------------------------------


def test_an_informative_object_reads_its_own_level(landscape):
    """At high depth the posterior is the object's own level against the bulk: 50 kb pieces at the captured level read
    the captured level over the bulk, at the depleted level 1, and at a third of the captured level (inside the enriched
    basin, below its peak) a third of it."""
    r = capture_efficiencies(
        landscape, *_field([50000] * 3, [RHO_CAP, RHO_DEP, RHO_CAP / 3.0], [RHO_DEP] * 2)
    )
    c, rb = r.efficiency_region, r.bulk_density
    assert abs(np.log(c[1])) < 0.2
    assert abs(np.log(c[0] * rb / RHO_CAP)) < 0.2
    assert abs(np.log(c[2] * rb * 3.0 / RHO_CAP)) < 0.2


def test_a_density_above_its_basins_peak_reads_its_own_level(landscape):
    """The read-out is the object's own level on either side of its basin's peak: a 50 kb piece at one and a half times
    the captured level reads one and a half times what a piece at the captured level reads, to the grid's step, and
    above the saturated level. Read at its enriched basin's peak it would read what the piece at the peak reads, a ratio
    of 1, and every captured object above a peak would read low."""
    r = capture_efficiencies(landscape, *_field([50000] * 2, [RHO_CAP, 1.5 * RHO_CAP], [RHO_DEP]))
    c = r.efficiency_region
    step = float(np.diff(landscape.log_rho)[0])
    assert r.saturation > 100.0
    assert abs(np.log(c[1] / c[0] / 1.5)) < step
    assert c[1] > r.saturation


def test_a_one_basin_landscape_reads_exactly_one_everywhere():
    """A library with no capture is one gDNA level: the prior is a point, so every object reads exactly 1 (shallow or
    deep, at the bulk or far from it) and the saturated level is 1. No decision is taken."""
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


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_anchors_alone_are_one_basin_and_read_exactly_one(seed):
    """A library with no gDNA evidence trains its landscape on its zero-count anchors alone, on a grid that spans every
    slot the prior is read at, RNA-rich exons included, and refits three times (each refit's E-step places the
    location-free anchors where the previous landscape is). An anchor says only "below the resolution wall", so the
    density has one maximum: the census finds one basin, the saturated level is 1, and every object holding no gDNA
    fragment reads exactly 1, however long."""
    rng = np.random.default_rng(seed)
    n_anchor, n_exon = 3000, 2000
    eff_anchor = np.exp(rng.uniform(np.log(50.0), np.log(50_000.0), n_anchor))
    eff_exon = np.exp(rng.uniform(np.log(100.0), np.log(5_000.0), n_exon))
    mass_exon = eff_exon * np.exp(rng.uniform(np.log(1e-3), np.log(10.0), n_exon))
    domain = (np.r_[np.zeros(n_anchor), mass_exon], np.r_[eff_anchor, eff_exon])
    ls = None
    for _refit in range(3):
        ls = fit_landscape(
            np.zeros(n_anchor),
            np.zeros(n_anchor),
            eff_anchor,
            np.zeros(n_anchor),
            anchor=np.ones(n_anchor, dtype=bool),
            domain=domain,
            prev=ls,
        )
        assert ls is not None
        _bulk, above = split_basins(census(ls))
        assert above == ()
        r = capture_efficiencies(ls, np.zeros(n_exon), eff_exon, np.zeros(n_anchor), eff_anchor)
        assert r.saturation == 1.0
        np.testing.assert_array_equal(r.efficiency_region, 1.0)
        np.testing.assert_array_equal(r.efficiency_boundary, 1.0)


def test_a_lone_level_lifts_an_object_without_evidence_by_at_most_e_to_the_mu_delta():
    """An object whose own count says nothing reads the prior. In the log readout a lone level of prior mass ``μ`` at
    ``Δ`` nats above the bulk lifts it by at most ``e^{μΔ}``, where the posterior mean of the density itself would read
    ``1 + μ(e^Δ − 1)``, eleven times the bulk here, inventing capture at a zero-gDNA object."""
    mu, delta = 0.01, np.log(1000.0)
    ls = _hand([1e-3, 1.0], [1.0 - mu, mu], widths=0.15)
    r = capture_efficiencies(ls, np.zeros(1), np.full(1, 1e-12), np.zeros(0), np.zeros(0))
    e = float(r.efficiency_region[0])
    assert 1.0 <= e <= np.exp(mu * delta) * 1.01
    assert e < 0.2 * (1.0 + mu * (np.exp(delta) - 1.0))


def test_a_level_below_the_bulk_reads_below_one():
    """The bulk is the largest basin, not the lowest: a well-measured object at a level below it reads below 1, as its
    own count says."""
    ls = _hand([1e-5, 1e-3, 1.0], [0.05, 0.85, 0.10])
    S = np.full(1, 1e7)
    r = capture_efficiencies(ls, 1e-5 * S, S, np.zeros(0), np.zeros(0))
    assert abs(np.log(r.efficiency_region[0] / 1e-2)) < 0.2


@pytest.mark.parametrize("captured_mass", [0.4, 0.6])
def test_a_well_measured_capture_ratio_survives_a_captured_majority(
    captured_mass,
):
    """Two levels, 0.001 and 0.1 fragments per bp, read by objects with a million positions each: the true ratio is
    exactly 100. When the captured level holds 60% of the population it is the largest basin and becomes the unit, so
    the uncaptured object reads 0.01 and the captured one 1, still a ratio of 100. Absorbing every level at or below the
    bulk into the null (a one-sided prior) would read both at 1 and erase a ratio no depth can restore.
    """
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
    """A region's efficiency is its contained count's and a boundary's its crossing count's: a depleted intron beside
    captured crossings stays at the bulk, and changing every boundary's count leaves every region's efficiency
    bit-identical (and the reverse), so the crossings are priced once, at the boundary that holds them. Apportioning the
    crossings onto the pieces within reach would lift the intron and move the regions with the boundaries.
    """
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
    """A piece shorter than every fragment has no gDNA support and no likelihood (no gDNA fragment could sit there to
    witness its capture), so it reads exactly 1, the null, whatever its neighbours hold."""
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
    """The objects are read in row tiles so a million never exist against the grid at once, and the tile height is a
    resource setting, so it moves no bit: an object's moments are sums over its own row, in an order no tile changes.
    The tiled reader equals the dense reader with its bulk term in closed form; a level moment computed as a
    matrix–vector product would fail at a one-row tile.
    """
    ls = _hand([1e-5, 1e-3, 0.03, 0.2, 1.0], [0.02, 0.80, 0.03, 0.03, 0.12])
    rng = np.random.default_rng(12)
    S = np.where(rng.random(1000) < 0.05, 0.0, rng.uniform(1.0, 1e5, 1000))
    k = np.floor(np.exp(rng.uniform(np.log(1e-6), np.log(5.0), 1000)) * S)
    tiled = capture_efficiencies(ls, k, S, k[:10], S[:10]).efficiency_region
    np.testing.assert_allclose(tiled, _dense_reference(ls, k, S), rtol=1e-12)
    for rows in (1, 7, 333, k.size):
        monkeypatch.setattr(CE, "block_rows", lambda n_grid, n_arrays, rows=rows: rows)
        got = capture_efficiencies(ls, k, S, k[:10], S[:10]).efficiency_region
        assert np.array_equal(got, tiled), f"rows {rows}: {int((got != tiled).sum())} objects moved"


def test_the_reading_is_the_same_at_every_thread_count(monkeypatch):
    """The tiles are read on calibration's thread budget, each whole by one thread into its own rows, so the count moves
    no bit, and every tile is read: a tile no thread reached would leave its objects at 1, off the dense reader."""
    ls = _hand([1e-5, 1e-3, 0.03, 0.2, 1.0], [0.02, 0.80, 0.03, 0.03, 0.12])
    rng = np.random.default_rng(12)
    S = np.where(rng.random(1000) < 0.05, 0.0, rng.uniform(1.0, 1e5, 1000))
    k = np.floor(np.exp(rng.uniform(np.log(1e-6), np.log(5.0), 1000)) * S)
    monkeypatch.setattr(CE, "block_rows", lambda n_grid, n_arrays: 7)
    one = capture_efficiencies(ls, k, S, k[::-1], S[::-1], n_threads=1)
    np.testing.assert_allclose(one.efficiency_region, _dense_reference(ls, k, S), rtol=1e-12)
    for n_threads in (2, 5, 0):
        got = capture_efficiencies(ls, k, S, k[::-1], S[::-1], n_threads=n_threads)
        assert np.array_equal(got.efficiency_region, one.efficiency_region), n_threads
        assert np.array_equal(got.efficiency_boundary, one.efficiency_boundary), n_threads


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


# --- the weight: the reading fades as a whole where the library's gDNA cannot be seen --------------------------------


def test_the_visible_share_is_the_mass_weighted_poisson_chance_of_a_gdna_fragment():
    """The weight is ``Σ m (1 − e^{−k}) / Σ m`` over the objects that admit at most one RNA strand, ``k`` the gDNA booked
    on an object and ``m`` its unspliced count: the share of the reads whose object would show a gDNA fragment of its
    own. Weighted by the reads, not by the objects; the Poisson chance, not a cut at one fragment; and an object that
    admits both strands does not vote. Each of the three alternatives gives another number here (0.54, 0.90, 0.89)."""
    mass = np.array([100.0, 300.0, 600.0, 1000.0])
    gdna = np.array([0.0, 1.0, 40.0, 5.0])
    one_strand = np.array([True, True, True, False])
    expected = (300.0 * (1.0 - np.exp(-1.0)) + 600.0 * (1.0 - np.exp(-40.0))) / 1000.0
    assert visible_share(gdna, mass, one_strand) == pytest.approx(expected, rel=1e-14)


def test_an_object_admitting_both_rna_strands_does_not_vote():
    """A both-strand object's gDNA is where the landscape put it, not its own evidence (the landscape is not trained on
    it, for the same reason), so neither its booking nor its mass moves the weight."""
    one_strand = np.array([True, True, False])
    a = visible_share(np.array([0.2, 3.0, 0.0]), np.array([500.0, 500.0, 2000.0]), one_strand)
    b = visible_share(np.array([0.2, 3.0, 900.0]), np.array([500.0, 500.0, 1e5]), one_strand)
    assert a == b


def test_with_no_one_strand_mass_nothing_can_show_and_the_weight_is_zero():
    assert visible_share(np.array([5.0]), np.array([10.0]), np.array([False])) == 0.0
    assert visible_share(np.zeros(0), np.zeros(0), np.zeros(0, dtype=bool)) == 0.0


@pytest.mark.parametrize(
    "bad",
    [
        dict(gdna=np.array([np.nan, 1.0])),
        dict(mass=np.array([np.inf, 1.0])),
        dict(mass=np.array([1.0, 2.0, 3.0])),
        dict(one_strand=np.array([True])),
    ],
)
def test_a_non_finite_or_misshapen_share_input_is_refused(bad):
    kw = dict(
        gdna=np.array([0.5, 1.0]), mass=np.array([3.0, 4.0]), one_strand=np.array([True, True])
    )
    kw.update(bad)
    with pytest.raises(ValueError, match="visible_share"):
        visible_share(kw["gdna"], kw["mass"], kw["one_strand"])


def test_the_weight_raises_every_efficiency_and_the_saturation_alike(landscape):
    """``weighted`` raises the reader's answer to ``w``: every efficiency and the saturated level become ``x ** w``, the
    same fraction of each level's distance from the bulk in nats, and the bulk stays the unit. ``w = 1`` is the reader
    bit for bit; ``w = 0`` reads exactly 1 everywhere, the no-capture line."""
    r = capture_efficiencies(
        landscape,
        *_field([50000, 5000, 1000], [RHO_CAP, RHO_DEP, RHO_CAP / 3.0], [RHO_CAP, RHO_DEP]),
    )
    assert r.saturation > 100.0
    for w in (0.0, 0.37, 1.0):
        rw = weighted(r, w)
        np.testing.assert_array_equal(rw.efficiency_region, r.efficiency_region**w)
        np.testing.assert_array_equal(rw.efficiency_boundary, r.efficiency_boundary**w)
        assert rw.saturation == r.saturation**w
        assert rw.bulk_density == r.bulk_density
    one, none = weighted(r, 1.0), weighted(r, 0.0)
    np.testing.assert_array_equal(one.efficiency_region, r.efficiency_region)
    np.testing.assert_array_equal(one.efficiency_boundary, r.efficiency_boundary)
    assert one.saturation == r.saturation
    np.testing.assert_array_equal(none.efficiency_region, 1.0)
    np.testing.assert_array_equal(none.efficiency_boundary, 1.0)
    assert none.saturation == 1.0


@pytest.mark.parametrize("w", [-0.01, 1.01, np.nan, np.inf])
def test_a_weight_outside_the_unit_interval_is_refused(landscape, w):
    r = capture_efficiencies(landscape, *_field([5000], [RHO_CAP], [RHO_CAP]))
    with pytest.raises(ValueError, match="weighted"):
        weighted(r, w)


def test_the_reads_that_sit_where_no_gdna_shows_take_the_claim_with_them(landscape):
    """The hole between 0 and 1 % gDNA under capture: one probed exon holds enough gDNA to be read at the captured
    level, while the library's reads sit on exons holding none, which cannot be read and stay near the bulk; a
    transcript across them takes its lengths from a reading on some objects and none on others. Weighted by the share
    of the reads whose object shows gDNA (here under 0.1 %), the captured level fades from hundreds-fold to under 1 %."""
    k_reg, S, k_bnd, S_b = _field([5000] * 100, [0.0] * 99 + [RHO_CAP], [0.0] * 99)
    rna = np.r_[np.full(99, 2e5), 0.0]
    w = visible_share(np.r_[k_reg, k_bnd], np.r_[k_reg + rna, k_bnd], np.ones(199, dtype=bool))
    r = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b)
    out = weighted(r, w)
    assert w < 1e-3
    assert r.efficiency_region[-1] > 100.0
    assert out.efficiency_region[-1] < 1.01 and out.saturation < 1.01


def test_one_probe_keeps_its_reading_in_a_library_whose_gdna_shows(landscape):
    """Capture inference is local: one probe must not move the claim. The weight asks whether gDNA shows where the reads
    sit, not how many objects are enriched, so a single probed exon among exons that each show a few gDNA fragments
    (``k`` about 5, ``1 − e^{−k}`` about 0.99) keeps its reading to within a few percent. A weight scaled by the
    landscape's spread would read one probe in a hundred as weak capture and fade it."""
    k_reg, S, k_bnd, S_b = _field([5000] * 100, [RHO_DEP] * 99 + [RHO_CAP], [RHO_DEP] * 99)
    rna = np.r_[np.full(99, 2e5), 0.0]
    w = visible_share(np.r_[k_reg, k_bnd], np.r_[k_reg + rna, k_bnd], np.ones(199, dtype=bool))
    r = capture_efficiencies(landscape, k_reg, S, k_bnd, S_b)
    out = weighted(r, w)
    assert 0.98 < w < 1.0
    assert r.efficiency_region[-1] > 100.0
    assert out.efficiency_region[-1] > 0.9 * r.efficiency_region[-1]


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
    """A solve with no refit fits no landscape, so the result carries no bulk density, a saturated level of 1 and
    efficiencies of exactly 1 on both axes; and it carries gDNA's conserved share at every boundary, its two flanking
    regions' shares on the gDNA pmf at the unbounded reach of a template that does not end
    (`effective_length.conserved_cut_shares`)."""
    from _synthetic import make_strand_models, make_synthetic_payload, make_synthetic_sj

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
        rna_fl_pmf=spike_pmf(),
        config=CalibrationConfig(calib_refit_iters=0),
        sj=make_synthetic_sj(),
        boundary_flags=np.zeros(int(payload.n_boundaries), np.uint16),
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
