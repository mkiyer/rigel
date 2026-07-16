"""gDNA FL effective lengths: boundary = μ_FL, region = E_f[max(0, L−ℓ)]."""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import (
    boundary_eff_length,
    boundary_side_crossing_count_eff_length,
    boundary_side_eff_length,
    fl_mean,
    region_eff_length,
    spliced_side_eff_length,
)


def _brute_region(L: float, pmf: np.ndarray) -> float:
    """Σ_{ℓ≤L} (L−ℓ+1) f(ℓ) — the DISCRETE start-position count.

    NB this restates the closed form; it does not independently validate the geometry (the old version
    said "(L−ℓ)" and matched an implementation that was off by one, so the pair agreed while both were
    wrong). The geometry itself is enumerated in
    `test_region_eff_length_is_the_discrete_start_position_count`.
    """
    p = pmf / pmf.sum()
    return float(sum((L - ell + 1.0) * p[ell] for ell in range(p.shape[0]) if ell <= L))


def _spike(at: int, n: int = 801) -> np.ndarray:
    p = np.zeros(n, dtype=np.float64)
    p[at] = 1.0
    return p


def test_fl_mean_normalizes_unnormalized_counts():
    counts = np.zeros(801)
    counts[100] = 3.0  # raw counts, not a pmf
    counts[300] = 1.0
    # mean = (3*100 + 1*300) / 4 = 150
    assert fl_mean(counts) == 150.0


def test_boundary_eff_length_is_mu_fl_independent_of_regions():
    # THE regression against the old "tent" error: a boundary's gDNA exposure is
    # the FL mean, NOT min(L_left, L_right, …). With μ_FL = 300, the boundary
    # exposure is 300 regardless of the 400 bp / 200 bp neighbours.
    pmf = _spike(300)
    assert boundary_eff_length(pmf) == 300.0  # not capped at min(400, 200) = 200


def test_region_eff_length_spike():
    pmf = _spike(300)
    # Start-position COUNTS for a 300 bp fragment (the accumulator increments a contained fragment by +1,
    # so the divisor is the integer opportunity count `max(0, L-ell+1)`):
    #   L=400 -> 101 positions;  L=200 -> 0 (cannot be contained);  L=300 -> exactly 1 (it just fits).
    # The last cell is the regression: it previously asserted 0, encoding the off-by-one as expected
    # behaviour and floor-dividing every L==ell region to a ~1e9 density.
    np.testing.assert_allclose(
        region_eff_length(np.array([400.0, 200.0, 300.0]), pmf), [101.0, 0.0, 1.0]
    )


def test_region_eff_length_matches_brute_force():
    rng_pmf = np.zeros(801, dtype=np.float64)
    rng_pmf[50:801] = 1.0  # uniform over [50, 800]; mean = 425
    for L in (0.0, 100.0, 300.0, 425.0, 800.0, 1000.0):
        np.testing.assert_allclose(
            region_eff_length(np.array([L]), rng_pmf)[0], _brute_region(L, rng_pmf), rtol=1e-12
        )


def test_region_eff_length_large_L_is_L_plus_one_minus_mu():
    pmf = np.zeros(801, dtype=np.float64)
    pmf[200:401] = 1.0  # mean = 300
    mu = fl_mean(pmf)
    # For L well beyond the support, contained exposure → L + 1 − μ_FL (the +1 is the discrete count;
    # it is irrelevant at this scale — 4701 vs 4700 — and decisive only at L ≈ ell, where it is 1 vs 0).
    np.testing.assert_allclose(region_eff_length(np.array([5000.0]), pmf)[0], 5001.0 - mu)


def test_region_eff_length_zero_and_nonnegative():
    pmf = _spike(150)
    eff = region_eff_length(np.array([0.0, 50.0, 150.0]), pmf)
    assert np.all(eff >= 0.0)
    assert eff[0] == 0.0  # zero-length region: no contained exposure
    assert eff[1] == 0.0  # 50 bp region cannot contain a 150 bp fragment
    assert eff[2] == 1.0  # a 150 bp region fits a 150 bp fragment at exactly ONE start position


def test_boundary_face_density_length_matches_the_accumulator_deposit_rule():
    """The per-face DENSITY length must recover ρ from a boundary face's deposited mass.

    Regression guard for a systematic 2× frame error on every boundary↔region message: the face length
    was ``E[min(ℓ,R)]`` (the crossing COUNT) while the accumulator deposits only that face's *share* of a
    straddling fragment, so every message read ρ/2. Pinned against the reference accumulator's rule
    (``share = (slice_len/ℓ)/n_cross``, ``accumulator/00_design.md`` §4.3), which integrates to exactly
    ``min(ℓ,R)/2`` for EVERY R — the short-R cells below are what discriminate it from the alternatives
    (a one-sided half-triangle ``min(ℓ,R)²/(2ℓ)`` gives 25/6.25, not 50/25).
    """
    ell = 200
    fl = np.zeros(400)
    fl[ell] = 1.0
    R = np.array([600.0, 400.0, 200.0, 100.0, 50.0])
    # E_face = min(ℓ,R)/2, exact for all R
    expect = np.minimum(ell, R) / 2.0
    assert np.allclose(boundary_side_eff_length(fl, R), expect)
    # and it is exactly half the crossing-count length, which keeps the whole fragment
    assert np.allclose(
        boundary_side_eff_length(fl, R), 0.5 * boundary_side_crossing_count_eff_length(fl, R)
    )
    # the one-sided spliced divisor is a DIFFERENT quantity and must not be conflated with it
    assert not np.allclose(spliced_side_eff_length(fl, R), boundary_side_eff_length(fl, R))


def test_region_eff_length_is_the_discrete_start_position_count():
    """`region_eff_length` divides an INTEGER contained count (the accumulator adds exactly +1 for a
    fragment with all bp in one region, `accumulator/00_design.md` §4.1-4.2), so it must be the INTEGER
    opportunity count: a length-l fragment fits wholly in a length-L region at `max(0, L-l+1)` positions.

    Regression: this read `max(0, L-l)`, which is EXACTLY ZERO when L equals the shortest fragment — a
    division by zero, floored to _EPS by callers, producing ~1e9 densities. Measured on ambig_dense_10mb:
    211/1698 regions had eff <= 1e-9. Very short exons are deliberate in these suites (they should relay
    messages, not project a belief) — they need a small eff-length, never a zero one.
    """
    import numpy as np

    from rigel.calibration.effective_length import region_eff_length

    def enumerate_starts(L, ell):
        return sum(1 for s in range(0, L + 1) if s + ell <= L)

    for ell in (50, 100, 200):
        pmf = np.zeros(max(ell, 250) + 1)
        pmf[ell] = 1.0  # point-mass FL: eff-length must equal the raw enumeration
        for L in (1, 10, ell - 1, ell, ell + 1, 2 * ell, 5 * ell):
            got = float(region_eff_length(np.array([float(L)]), pmf)[0])
            want = float(enumerate_starts(L, ell))
            assert got == want, f"L={L} ell={ell}: eff={got} but {want} start positions exist"

    # THE REGRESSION ITSELF: L == the shortest fragment must give 1 opportunity, never 0.
    pmf = np.zeros(201)
    pmf[100] = 1.0
    assert float(region_eff_length(np.array([100.0]), pmf)[0]) == 1.0

    # A region SHORTER than every fragment genuinely contains none -> 0 is correct there.
    assert float(region_eff_length(np.array([99.0]), pmf)[0]) == 0.0

    # Short regions stay SMALL (so their density is imprecise and they relay) but never zero-and-dividing.
    short = region_eff_length(np.array([100.0, 120.0, 150.0]), pmf)
    long_ = region_eff_length(np.array([5000.0]), pmf)
    assert np.all(short > 0.0)
    assert np.all(short < 0.05 * long_[0]), "a short region must stay a small divisor"
