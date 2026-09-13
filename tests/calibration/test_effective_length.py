"""Effective lengths: ONE placement formula, per component, per frame.

An effective length is the expected number of admissible fragment START POSITIONS — the divisor that
turns an observed count into a start density — and there is one formula per frame and nothing else::

    contained   E_f[ (region_len − w + 1)+ ]
    crossing    E_f[ max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)) ]

The crossing formula covers BOTH boundary kinds and both components: mean fragment length is its
large-reach limit rather than a separate case, since gDNA's template is the chromosome and its
reaches are unbounded, while RNA's template ends where its transcript ends. Every test here
ENUMERATES integer start positions instead of restating the closed form, because a "brute force"
written from the same algebra as the implementation agrees with it while both are off by one.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import norm

from rigel.calibration import effective_length as el
from rigel.calibration.effective_length import (
    contained_eff_length,
    crossing_eff_length,
    fl_mean,
)

UNBOUNDED = el.UNBOUNDED_REACH


def _spike(at: int, n: int = 1301) -> np.ndarray:
    p = np.zeros(n, dtype=np.float64)
    p[at] = 1.0
    return p


def _normal_pmf(mean: float, sd: float, n: int = 1301) -> np.ndarray:
    p = np.diff(norm.cdf(np.arange(0, n + 1), mean, sd))
    p[0] = 0.0
    return p / p.sum()


# ---------------------------------------------------------------------------
# enumeration — the ground truth, independent of the closed form
# ---------------------------------------------------------------------------


def _enumerate_contained(region_len: int, w: int) -> int:
    """Integer starts placing a length-``w`` fragment wholly inside ``[0, region_len)``."""
    starts = np.arange(-w - 2, region_len + w + 2)
    return int(np.sum((starts >= 0) & (starts + w <= region_len)))


def _enumerate_crossing(w: int, reach_lo: int, reach_hi: int) -> int:
    """Integer placements of a length-``w`` molecule across a point with ``reach`` bases either side.

    The molecule occupies ``a`` bases to the left of the point and ``w − a`` to the right; it crosses iff
    both are ≥ 1, and it must FIT in what remains of its own template on each side.
    """
    return int(sum(1 for a in range(1, w) if a <= reach_lo and (w - a) <= reach_hi))


# ---------------------------------------------------------------------------
# contained
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("region_len", [1, 10, 60, 100, 150, 500, 2000])
@pytest.mark.parametrize("w", [1, 2, 59, 100, 101, 200, 501])
def test_contained_is_the_enumerated_start_count(region_len, w):
    got = contained_eff_length(np.array([float(region_len)]), _spike(w))[0]
    assert got == pytest.approx(float(_enumerate_contained(region_len, w)))


def test_contained_at_a_region_exactly_one_fragment_long_is_ONE_not_zero():
    """The ``+1`` is the discrete count of start positions, not a correction factor.

    Dropping it makes the divisor exactly 0 when a region is one fragment long — a division by zero
    that a floor turns into an absurd density on every short region of a fine partition, and short
    regions are a large share of one.
    """
    assert contained_eff_length(np.array([100.0]), _spike(100))[0] == pytest.approx(1.0)


def test_contained_beyond_the_pmf_support_is_region_plus_one_minus_mean():
    pmf = _normal_pmf(200.0, 50.0)
    got = contained_eff_length(np.array([5000.0]), pmf)[0]
    assert got == pytest.approx(5000.0 + 1.0 - fl_mean(pmf), rel=1e-9)


def test_contained_is_never_negative():
    assert float(contained_eff_length(np.array([0.0, 1.0, 5.0]), _spike(200)).min()) >= 0.0


# ---------------------------------------------------------------------------
# crossing — the one formula, both boundary kinds, both components
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("w", [2, 5, 40, 100, 200, 350])
@pytest.mark.parametrize("reach_lo,reach_hi", [(1, 1), (3, 500), (50, 50), (100, 150), (400, 400)])
def test_crossing_is_the_enumerated_placement_count(w, reach_lo, reach_hi):
    got = crossing_eff_length(_spike(w), np.array([float(reach_lo)]), np.array([float(reach_hi)]))[
        0
    ]
    assert got == pytest.approx(float(_enumerate_crossing(w, reach_lo, reach_hi)))


def _crossing_matrix(fl_pmf, reach_lo, reach_hi):
    """The brute force: every object against every fragment length, the four-way min materialised."""
    p = np.asarray(fl_pmf, np.float64) / np.sum(fl_pmf)
    lengths = np.arange(p.shape[0], dtype=np.float64)
    lo, hi = np.broadcast_arrays(np.asarray(reach_lo, float), np.asarray(reach_hi, float))
    lo_col, hi_col = lo.reshape(-1, 1), hi.reshape(-1, 1)
    placements = np.minimum(
        np.minimum(lengths - 1.0, np.minimum(lo_col, hi_col)), lo_col + hi_col - lengths + 1.0
    )
    return (np.maximum(placements, 0.0) @ p).reshape(lo.shape)


def test_crossing_closed_form_equals_the_matrix_brute_force_on_real_valued_reaches():
    """The closed form over the pmf's cumulative sums against the materialised ``(objects × lengths)``
    brute force, on real-valued reaches spanning below one base, inside the support, across its end and
    UNBOUNDED, on a pmf with mass at every length. PERTURBATION: any one of the three segment sums
    dropped, or a segment end off by one length, moves a reach in the middle of the support by more than
    the tolerance."""
    rng = np.random.default_rng(7)
    pmf = _normal_pmf(200.0, 50.0) + 1e-4  # mass at every length, so every segment end matters
    lo = np.concatenate(
        [rng.uniform(0.0, 600.0, 200), [0.0, 0.4, 1.0, 1.5, 349.0, 350.0, 351.0, UNBOUNDED]]
    )
    hi = np.concatenate(
        [rng.uniform(0.0, 600.0, 200), [500.0, 0.6, 1.0, 2.5, 349.0, 350.0, 2000.0, UNBOUNDED]]
    )
    got = crossing_eff_length(pmf, lo, hi)
    want = _crossing_matrix(pmf, lo, hi)
    np.testing.assert_allclose(got, want, rtol=1e-12, atol=1e-12)
    # non-integer reaches: the segment ends are floors, and the brute force agrees with them
    frac = crossing_eff_length(pmf, np.array([120.3, 120.9]), np.array([300.7, 300.1]))
    np.testing.assert_allclose(
        frac, _crossing_matrix(pmf, [120.3, 120.9], [300.7, 300.1]), rtol=1e-12
    )


def test_crossing_at_UNBOUNDED_reach_is_the_mean_length_minus_one():
    """gDNA's template is the chromosome, so it never tapers — and then the divisor is just ``mu − 1``.

    This is why mean fragment length is the large-reach LIMIT of the placement formula rather than a
    separate gDNA case: one formula, one code path, no branch on component.
    """
    pmf = _normal_pmf(200.0, 50.0)
    got = crossing_eff_length(pmf, np.array([UNBOUNDED]), np.array([UNBOUNDED]))[0]
    assert got == pytest.approx(fl_mean(pmf) - 1.0, rel=1e-9)


def test_crossing_is_SYMMETRIC_in_the_two_reaches():
    """A crossing point does not know which side is which; only the pair matters."""
    pmf = _normal_pmf(200.0, 50.0)
    a = crossing_eff_length(pmf, np.array([120.0]), np.array([300.0]))[0]
    b = crossing_eff_length(pmf, np.array([300.0]), np.array([120.0]))[0]
    assert a == pytest.approx(b)


def test_a_ZERO_reach_gives_ZERO_opportunity_not_a_floor():
    """An object with no opportunity for a component must emit NOTHING, never a floored division.

    A floored divisor turns "this component cannot be here" into a rate, and a slot that then defaults
    to all-gDNA seeds false gDNA into its neighbours through the messages. Zero is the correct answer
    here and must survive as zero.
    """
    pmf = _normal_pmf(200.0, 50.0)
    assert crossing_eff_length(pmf, np.array([0.0]), np.array([500.0]))[0] == 0.0
    assert crossing_eff_length(pmf, np.array([500.0]), np.array([0.0]))[0] == 0.0


def test_crossing_reproduces_the_MEASURED_taper_table():
    """An independent cross-check of the taper against values derived outside this file, RNA N(200,50).

    They mix two conventions and this test pins both: the first four are SYMMETRIC (both reaches = R),
    while the last, at a first exon, is ONE-SIDED — a first exon is short on one side and long on the
    other, and reading it as symmetric gives a different number entirely.
    """
    pmf = _normal_pmf(200.0, 50.0)

    def sym(R):
        return crossing_eff_length(pmf, np.array([float(R)]), np.array([float(R)]))[0]

    assert sym(200) == pytest.approx(160.1, abs=0.1)
    assert sym(147) == pytest.approx(87.8, abs=0.5)
    assert sym(100) == pytest.approx(19.6, abs=0.3)
    assert sym(550) == pytest.approx(199.0, abs=0.6)
    one_sided = crossing_eff_length(pmf, np.array([50.0]), np.array([UNBOUNDED]))[0]
    assert one_sided == pytest.approx(50.0, abs=0.1)


def test_the_taper_is_a_MULTI_FOLD_error_if_ignored():
    """Not a refinement: at a short reach the correct divisor is a fraction of the mean length, so using
    the mean blindly over-divides by many fold and under-reads the density by the same factor."""
    pmf = _normal_pmf(200.0, 50.0)
    unbounded = crossing_eff_length(pmf, np.array([UNBOUNDED]), np.array([UNBOUNDED]))[0]
    tapered = crossing_eff_length(pmf, np.array([100.0]), np.array([100.0]))[0]
    assert unbounded / tapered > 8.0


def test_crossing_is_vectorised_over_objects_and_agrees_elementwise():
    """The real call site passes a whole boundary axis at once; a broadcasting slip would mix objects."""
    pmf = _normal_pmf(200.0, 50.0)
    lo = np.array([50.0, 120.0, UNBOUNDED, 0.0])
    hi = np.array([500.0, 300.0, UNBOUNDED, 500.0])
    batch = crossing_eff_length(pmf, lo, hi)
    for i in range(lo.size):
        one = crossing_eff_length(pmf, lo[i : i + 1], hi[i : i + 1])[0]
        assert batch[i] == pytest.approx(one)


# ---------------------------------------------------------------------------
# one answer per question: the per-face divisors must stay gone
# ---------------------------------------------------------------------------


def test_the_THREE_OLD_DIVISORS_ARE_GONE():
    """A divisor for a deposit rule that no longer exists must not survive beside the one that does.

    ``boundary_side_eff_length`` (``E[min(l,R)]/2``), ``spliced_side_eff_length`` (``E[min^2/2l]``)
    and ``boundary_side_crossing_count_eff_length`` all divided a per-FACE quantity, and a contiguous
    boundary has no faces — it is a 0-bp boundary with one set of numbers. Two answers for one
    question is how an exact factor of 2 survives a whole file of assertions.
    """
    for dead in (
        "boundary_side_eff_length",
        "boundary_side_crossing_count_eff_length",
        "spliced_side_eff_length",
        "boundary_eff_length",
        "region_eff_length",
    ):
        assert not hasattr(el, dead), f"{dead} still exists"
