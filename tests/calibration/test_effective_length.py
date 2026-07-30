"""Effective lengths: ONE placements formula, per component, per frame.

    Derivation: ``docs/NODE_DENSITY_DERIVATION.md``   ·   Design: ``ACCUMULATOR_DESIGN.md`` §7

An effective length is the expected number of admissible fragment START POSITIONS — the divisor that
turns an observed count into a start density. There is one formula per frame and nothing else:

    contained   E_f[ (node_len − w + 1)+ ]
    crossing    E_f[ max(0, min(w−1, R_lo, R_hi, R_lo + R_hi − w + 1)) ]

⭐ **The crossing formula covers BOTH edge kinds and both components.** Mean fragment length is its
large-reach limit, not a separate case: gDNA's template is the chromosome, so its reaches are unbounded
and its divisor collapses to ``mu − 1``. RNA's template ends where its transcript ends.

⚠ Every test here enumerates integer start positions rather than restating the closed form. An earlier
version of this file computed the "brute force" from the same algebra as the implementation, so the pair
agreed while BOTH were off by one.
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


def _enumerate_contained(node_len: int, w: int) -> int:
    """Integer starts placing a length-``w`` fragment wholly inside ``[0, node_len)``."""
    starts = np.arange(-w - 2, node_len + w + 2)
    return int(np.sum((starts >= 0) & (starts + w <= node_len)))


def _enumerate_crossing(w: int, reach_lo: int, reach_hi: int) -> int:
    """Integer placements of a length-``w`` molecule across a point with ``reach`` bases either side.

    The molecule occupies ``a`` bases to the left of the point and ``w − a`` to the right; it crosses iff
    both are ≥ 1, and it must FIT in what remains of its own template on each side.
    """
    return int(sum(1 for a in range(1, w) if a <= reach_lo and (w - a) <= reach_hi))


# ---------------------------------------------------------------------------
# contained
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("node_len", [1, 10, 60, 100, 150, 500, 2000])
@pytest.mark.parametrize("w", [1, 2, 59, 100, 101, 200, 501])
def test_contained_is_the_enumerated_start_count(node_len, w):
    got = contained_eff_length(np.array([float(node_len)]), _spike(w))[0]
    assert got == pytest.approx(float(_enumerate_contained(node_len, w)))


def test_contained_at_a_node_exactly_one_fragment_long_is_ONE_not_zero():
    """⚠ The ``+1`` is the discrete count of start positions, not a correction factor.

    Dropping it makes the divisor exactly 0 when a node is one fragment long — a division by zero that
    was floored to an epsilon and produced densities of ~1e9 on 12.4 % of fine-partition nodes.
    """
    assert contained_eff_length(np.array([100.0]), _spike(100))[0] == pytest.approx(1.0)


def test_contained_beyond_the_pmf_support_is_node_plus_one_minus_mean():
    pmf = _normal_pmf(200.0, 50.0)
    got = contained_eff_length(np.array([5000.0]), pmf)[0]
    assert got == pytest.approx(5000.0 + 1.0 - fl_mean(pmf), rel=1e-9)


def test_contained_is_never_negative():
    assert float(contained_eff_length(np.array([0.0, 1.0, 5.0]), _spike(200)).min()) >= 0.0


# ---------------------------------------------------------------------------
# crossing — the one formula, both edge kinds, both components
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("w", [2, 5, 40, 100, 200, 350])
@pytest.mark.parametrize("reach_lo,reach_hi", [(1, 1), (3, 500), (50, 50), (100, 150), (400, 400)])
def test_crossing_is_the_enumerated_placement_count(w, reach_lo, reach_hi):
    got = crossing_eff_length(_spike(w), np.array([float(reach_lo)]), np.array([float(reach_hi)]))[
        0
    ]
    assert got == pytest.approx(float(_enumerate_crossing(w, reach_lo, reach_hi)))


def test_crossing_at_UNBOUNDED_reach_is_the_mean_length_minus_one():
    """gDNA's template is the chromosome, so it never tapers — and then the divisor is just ``mu − 1``.

    ⭐ This is why mean fragment length is the large-reach LIMIT of the placement formula rather than a
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

    ``CARRY_FORWARD.md`` §3 trap 23: a "no data" default of 100 % gDNA was actively seeding false gDNA
    into neighbouring exons. Zero is the correct answer here and must survive as zero.
    """
    pmf = _normal_pmf(200.0, 50.0)
    assert crossing_eff_length(pmf, np.array([0.0]), np.array([500.0]))[0] == 0.0
    assert crossing_eff_length(pmf, np.array([500.0]), np.array([0.0]))[0] == 0.0


def test_crossing_reproduces_the_MEASURED_taper_table():
    """⭐ An independent cross-check: ``CARRY_FORWARD.md`` §2's published table, RNA N(200,50).

    ⚠ That table mixes two conventions and this test pins both. The first four entries are SYMMETRIC
    (both reaches = R); the last, captioned "at a first exon", is ONE-SIDED — a first exon is short on
    one side and long on the other. Recomputing them here is what revealed that; the numbers agree to
    the precision the document quotes.
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
    """The real call site passes a whole edge axis at once; a broadcasting slip would mix objects."""
    pmf = _normal_pmf(200.0, 50.0)
    lo = np.array([50.0, 120.0, UNBOUNDED, 0.0])
    hi = np.array([500.0, 300.0, UNBOUNDED, 500.0])
    batch = crossing_eff_length(pmf, lo, hi)
    for i in range(lo.size):
        one = crossing_eff_length(pmf, lo[i : i + 1], hi[i : i + 1])[0]
        assert batch[i] == pytest.approx(one)


# ---------------------------------------------------------------------------
# what died
# ---------------------------------------------------------------------------


def test_the_THREE_OLD_DIVISORS_ARE_GONE():
    """The mass-era divisors described a deposit rule that no longer exists.

    ``boundary_side_eff_length`` (``E[min(l,R)]/2``), ``spliced_side_eff_length`` (``E[min^2/2l]``) and
    ``boundary_side_crossing_count_eff_length`` all divided a per-FACE quantity, and a contiguous edge no
    longer has faces — it is a 0-bp line with one set of numbers. Keeping them would leave two answers
    for one question, which is how an exact factor of 2 survived 29 tests (``CARRY_FORWARD.md`` §3
    trap 2).
    """
    for dead in (
        "boundary_side_eff_length",
        "boundary_side_crossing_count_eff_length",
        "spliced_side_eff_length",
        "boundary_eff_length",
        "region_eff_length",
    ):
        assert not hasattr(el, dead), f"{dead} still exists"
