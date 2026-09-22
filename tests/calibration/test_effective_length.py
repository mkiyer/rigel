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


# ---------------------------------------------------------------------------


def _enumerate_tau(L: int, pmf: np.ndarray) -> np.ndarray:
    """τ(x) by brute force: every start of every length, each start spread over its w bases."""
    p = np.asarray(pmf, np.float64) / np.sum(pmf)
    tau = np.zeros(L)
    for w in range(1, p.shape[0]):
        if p[w] == 0.0:
            continue
        for s in range(0, L - w + 1):
            tau[s : s + w] += p[w] / w
    return tau


@pytest.mark.parametrize("L", [1, 30, 99, 100, 250, 600, 998, 999, 1000, 1500])
def test_the_taper_interval_sums_are_the_enumerated_per_base_weights(L):
    """Long templates through the table, short ones base by base — both against the brute force. 600 and
    998 sit between one and two fragment lengths, where the table's collapse of the four-way min is
    invalid: a threshold at one fragment length passed every other length here."""
    from rigel.calibration.effective_length import base_taper

    pmf = _normal_pmf(200.0, 60.0, n=501)
    tap = base_taper(pmf)
    tau = _enumerate_tau(L, pmf)
    cum = np.concatenate([[0.0], np.cumsum(tau)])
    x0 = np.array([0, 0, L // 3, max(L - 7, 0)])
    x1 = np.array([L, min(5, L), min(L // 3 + 40, L), L])
    np.testing.assert_allclose(
        tap.interval_sums(x0, x1, L), cum[x1] - cum[x0], rtol=1e-10, atol=1e-12
    )


def test_the_taper_interval_sums_take_one_template_length_per_interval():
    """One call over intervals on many templates — long and short, in any order — is bit-identical to
    the per-template calls: the vectorised path is the same arithmetic grouped, and it is what makes a
    real library's 457,000 transcripts one pass rather than one loop."""
    from rigel.calibration.effective_length import base_taper

    tap = base_taper(_normal_pmf(200.0, 60.0, n=501))
    rng = np.random.default_rng(7)
    L = rng.choice([1, 30, 99, 250, 999, 1000, 1500, 20000], size=200)
    x0 = (rng.random(200) * L).astype(np.int64)
    x1 = x0 + (rng.random(200) * (L - x0)).astype(np.int64)
    expect = np.array([tap.interval_sums([a], [b], int(n))[0] for a, b, n in zip(x0, x1, L)])
    np.testing.assert_array_equal(tap.interval_sums(x0, x1, L), expect)


@pytest.mark.parametrize("L", [1, 60, 500, 3000])
def test_the_taper_partitions_the_fl_marginal_length_exactly(L):
    """Σ_x τ(x) over the whole template is Σ_w f(w)(L − w + 1)⁺ — the per-base frame is a partition of
    the start count, never a second definition of it."""
    from rigel.calibration.effective_length import base_taper

    pmf = _normal_pmf(200.0, 60.0, n=501)
    w = np.arange(pmf.shape[0], dtype=np.float64)
    fl = float((pmf * np.maximum(L - w + 1.0, 0.0)).sum())
    assert base_taper(pmf).interval_sums(np.array([0]), np.array([L]), L)[0] == pytest.approx(fl)


def _regions(lengths):
    import pandas as pd

    from rigel.calibration.region_arrays import RegionArrays

    lengths = np.asarray(lengths, dtype=np.int64)
    b = np.concatenate([[0], np.cumsum(lengths)])
    frame = pd.DataFrame(
        {
            "region_id": np.arange(lengths.size, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * lengths.size, dtype="string"),
            "start": b[:-1],
            "end": b[1:],
            "length": lengths,
            "signature": np.zeros(lengths.size, np.uint8),
        }
    )
    return RegionArrays.from_frame(frame, {"chr1": 0})


def _enumerate_base_shares(lengths, e, pmf):
    """Base-starts per piece by brute force: every crossing placement (w, a) at boundary e, each base of
    the fragment attributed to the piece it lies in."""
    p = np.asarray(pmf, np.float64) / np.sum(pmf)
    b = np.concatenate([[0], np.cumsum(lengths)])
    B = b[e + 1]
    shares = np.zeros(len(lengths))
    for w in range(2, p.shape[0]):
        if p[w] == 0.0:
            continue
        for a in range(1, w):
            for x in range(B - a, B - a + w):
                q = int(np.searchsorted(b, x, side="right") - 1)
                if 0 <= q < len(lengths):
                    shares[q] += p[w] / w
    return shares


def test_the_crossing_base_shares_are_the_enumerated_placements():
    """Tiny pieces beside long ones: the enumeration attributes every base of every placement to the
    piece it lies in, and the closed form must match piece by piece."""
    from rigel.calibration.effective_length import crossing_base_shares

    pmf = _normal_pmf(60.0, 15.0, n=121)
    lengths = [500, 25, 300, 7, 9, 400]
    E, Q, A = crossing_base_shares(_regions(lengths), pmf)
    for e in range(len(lengths) - 1):
        got = np.zeros(len(lengths))
        np.add.at(got, Q[E == e], A[E == e])
        np.testing.assert_allclose(
            got, _enumerate_base_shares(lengths, e, pmf), rtol=1e-9, atol=1e-12
        )


def test_the_crossing_base_shares_sum_to_the_crossing_opportunity():
    """Σ_q share_eq = E_f[w − 1] at every boundary whose pieces reach a fragment on both sides — the
    crossing divisor, partitioned over the bases it counts."""
    from rigel.calibration.effective_length import UNBOUNDED_REACH, crossing_base_shares

    pmf = _normal_pmf(200.0, 60.0, n=501)
    ra = _regions([5000, 40, 1000, 40, 300, 30, 8000])
    E, Q, A = crossing_base_shares(ra, pmf)
    S_e = float(
        crossing_eff_length(pmf, np.array([UNBOUNDED_REACH]), np.array([UNBOUNDED_REACH]))[0]
    )
    total = np.zeros(6)
    np.add.at(total, E, A)
    np.testing.assert_allclose(total, S_e, rtol=1e-10)
    # and the left side alone is half of it: the two sides are symmetric
    left = np.zeros(6)
    np.add.at(left, E[Q < E + 1], A[Q < E + 1])
    np.testing.assert_allclose(left, S_e / 2.0, rtol=1e-10)


def test_a_length_model_with_mass_at_zero_length_tapers_without_a_nan():
    """A smoothed real-library model put mass at ``w = 0`` and the human library's factors went NaN
    (found on LBX0588 the day the taper landed): a zero-length fragment covers no base and is dropped,
    so the taper equals the taper of the same pmf with that mass removed, on long and short templates."""
    from rigel.calibration.effective_length import base_taper

    pmf = _normal_pmf(200.0, 60.0, n=501)
    with_zero = pmf.copy()
    with_zero[0] = 0.05
    clean = base_taper(pmf)
    dirty = base_taper(with_zero)
    for L in (30, 1500):
        x0, x1 = np.array([0, L // 3]), np.array([L, L // 3 + 20])
        got = dirty.interval_sums(x0, x1, L)
        assert np.all(np.isfinite(got))
        # the pmf is normalised over all its mass and the zero-length part then dropped, so the positive
        # lengths carry their weights scaled by the total, 1/1.05 — a ratio of tapers is unmoved
        np.testing.assert_allclose(got, clean.interval_sums(x0, x1, L) / 1.05, rtol=1e-12)
