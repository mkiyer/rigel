"""``capture_efficiency``: an object's capture efficiency is the posterior mean of its clipped gDNA density
against the reference, under the fitted landscape, from its OWN count — a region's contained count on its
contained support, a boundary's crossing count on its crossing support — no floor, no constant, and nothing
apportioned between them.

The falsification test was written first and verified failing on the shipped ruler through its own fixture
(`ISSUES: ruler-multimapper-floor-caps-the-correction`): an unprobed exon holding no gDNA fragment read the
multimapper floor ``1/(C+1)``, 0.024 on 40 RNA fragments against a depleted level of 0.001 (+3.19 nat). Here
the posterior reads the depleted level from the population. The rest pins the mechanism's parts: the plug-in
limit at high depth, a boundary's own efficiency, each object reading only its own count, and calibrate's
wiring (efficiencies exactly 1 with no reference). A piece too short to contain a fragment is read through
its boundaries by the length, not here (``tests/calibration/test_capture_eff_length.py``).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.capture_efficiency import capture_efficiencies
from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_eff_length
from rigel.calibration.landscape import fit_landscape

RHO_REF, RHO_DEP = 1.0, 1e-3


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


@pytest.fixture(scope="module")
def landscape():
    """A bimodal gDNA landscape fitted on a population: 900 depleted pieces at 10^-3 and 300 captured
    pieces at the reference, each a Poisson count on a 1 kb support, plus 200 zero-count anchors."""
    rng = np.random.default_rng(3)
    eff = np.full(1400, 1000.0)
    count = np.concatenate(
        [
            rng.poisson(RHO_DEP * 1000.0, 900).astype(float),
            rng.poisson(RHO_REF * 1000.0, 300).astype(float),
            np.zeros(200),
        ]
    )
    anchor = np.r_[np.zeros(1200, bool), np.ones(200, bool)]
    ls = fit_landscape(count, np.maximum(count * 2.0, 1.0), eff, np.full(1400, 0.1), anchor=anchor)
    assert ls is not None
    return ls


def _field(lengths, density, crossing_density):
    """Contained counts ``ρ_r S_r`` per region and crossing counts ``ρ_e S_E`` per boundary — each object's
    own count at its own density, which is all an efficiency reads."""
    S = _support(lengths)
    k_reg = np.asarray(density, dtype=np.float64) * S
    k_bnd = np.asarray(crossing_density, dtype=np.float64) * S_E
    return k_reg, S, k_bnd, np.full(k_bnd.shape[0], S_E)


# --- the falsification test ---------------------------------------------------------------------------


def test_an_unprobed_exon_with_no_gdna_fragment_reads_the_depleted_level_not_a_floor(landscape):
    """A 1 kb exon holding no gDNA at all among depleted neighbours: its posterior sits at the
    depleted level, within a factor of two of 0.001 and nowhere near the floor's 0.024 or 1."""
    k_reg, S, k_bnd, S_b = _field([5000, 1000, 5000], [RHO_DEP] * 3, [RHO_DEP] * 2)
    k_reg[1] = 0.0
    c, _ = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, S_b)
    assert abs(np.log(c[1] / RHO_DEP)) < np.log(2.0), c[1]
    assert c[1] < 0.005


# --- the parts ----------------------------------------------------------------------------------------


def test_a_well_measured_piece_reads_its_plug_in(landscape):
    """At high depth the posterior is the plug-in: 50 kb pieces (fifty depleted fragments, thousands of
    captured ones) at the reference read 1, at the depleted level 0.001, and at a third of the reference
    a third."""
    density = [RHO_REF, RHO_DEP, RHO_REF / 3.0, RHO_REF]
    k_reg, S, k_bnd, S_b = _field([50000] * 4, density, [RHO_REF] * 3)
    c, _ = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, S_b)
    assert c[0] > 0.95 and c[3] > 0.95
    assert abs(np.log(c[1] / RHO_DEP)) < 0.2
    assert abs(np.log(c[2] * 3.0)) < 0.2


def test_every_efficiency_lies_in_the_unit_interval(landscape):
    lengths = [5000, 40, 1000, 7, 9, 300, 5000]
    k_reg, S, k_bnd, S_b = _field(lengths, np.full(7, RHO_REF), np.full(6, 2.0 * RHO_REF))
    c, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, S_b)
    assert np.all(c >= 0.0) and np.all(c <= 1.0)
    assert cb.shape == (6,) and np.all(cb >= 0.0) and np.all(cb <= 1.0)


def test_a_boundarys_efficiency_is_its_own_crossing_counts_posterior(landscape):
    """A boundary reads its own crossing count on its crossing support: crossings at the reference read
    1, at the depleted level the depleted level, at half the reference — the crossing fragments at a
    probed exon's edge beside a depleted intron — about half."""
    k_reg, S, k_bnd, S_b = _field([5000] * 4, [RHO_REF] * 4, [RHO_REF, RHO_DEP, 0.5 * RHO_REF])
    _, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, S_b)
    assert cb[0] > 0.9
    assert cb[1] < 0.02
    assert 0.3 < cb[2] < 0.7, cb


def test_each_object_reads_only_its_own_count(landscape):
    """A region's efficiency is its contained count's and a boundary's its crossing count's: a depleted
    intron beside captured crossings stays depleted, and changing every boundary's count leaves every
    region's efficiency bit-identical (and the reverse) — the crossings are priced once, at the boundary
    that holds them. PERTURBATION: apportioning the crossings onto the pieces within reach lifts the
    intron and moves the regions with the boundaries."""
    lengths = [5000, 1000, 40, 1000, 5000]
    density = [RHO_DEP, RHO_DEP, RHO_REF, RHO_DEP, RHO_DEP]
    k_reg, S, k_bnd, S_b = _field(
        lengths, density, [RHO_DEP, 0.5 * RHO_REF, 0.5 * RHO_REF, RHO_DEP]
    )
    c, cb = capture_efficiencies(landscape, RHO_REF, k_reg, S, k_bnd, S_b)
    assert c[1] < 0.01 and c[3] < 0.01, (c[1], c[3])
    c2, cb2 = capture_efficiencies(landscape, RHO_REF, k_reg, S, np.full(4, RHO_REF) * S_E, S_b)
    np.testing.assert_array_equal(c2, c)
    assert not np.array_equal(cb2, cb)
    c3, cb3 = capture_efficiencies(landscape, RHO_REF, np.full(5, RHO_REF) * S, S, k_bnd, S_b)
    np.testing.assert_array_equal(cb3, cb)
    assert not np.array_equal(c3, c)


def test_an_object_with_no_support_reads_the_population(landscape):
    """A piece shorter than every fragment has no contained support and no count: it reads the
    population's clipped mean, whatever its neighbours hold — no length multiplies it."""
    lengths = [5000, 40, 5000]
    k1, S, kb, S_b = _field(lengths, [RHO_REF, RHO_REF, RHO_REF], [RHO_REF, RHO_REF])
    k2, _, kb2, _ = _field(lengths, [RHO_DEP, RHO_DEP, RHO_DEP], [RHO_DEP, RHO_DEP])
    assert S[1] == 0.0
    c1, _ = capture_efficiencies(landscape, RHO_REF, k1, S, kb, S_b)
    c2, _ = capture_efficiencies(landscape, RHO_REF, k2, S, kb2, S_b)
    assert c1[1] == c2[1]
    assert c1[0] > 0.9 and c2[0] < 0.01


def test_calibrate_publishes_the_efficiencies_and_exactly_one_without_a_reference():
    """The wiring: a solve with no refit has no landscape and no reference, and the result's
    efficiencies are exactly 1 on both axes; and the result carries gDNA's conserved share at every
    boundary — its two flanking regions' shares on the gDNA pmf, at the unbounded reach of a template that
    does not end (`effective_length.conserved_cut_shares`)."""
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
    assert res.gdna_reference_density is None
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
