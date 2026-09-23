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
    conserved_cut_shares,
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


def _small_pmf() -> np.ndarray:
    """Uneven mass on widths 2..40: a fragment crosses several of the fixtures' 1–3 bp pieces at once, and the enumeration stays
    small enough to run every placement through the reference accumulator."""
    w = np.arange(41, dtype=np.float64)
    p = np.where(w >= 2.0, 1.0 + 0.5 * np.sin(w), 0.0)
    return p / p.sum()


def _deposit_every_placement(bounds, blocks, pmf, chromosome=False):
    """THE DEPOSIT RULE, per object: every placement of every width through the reference accumulator
    (`tests/native/_accumulator_reference.py`), weighted by the pmf. A template is its ``blocks`` on one
    reference cut at ``bounds``, the introns between blocks annotated; ``chromosome`` places every start on
    the reference instead (gDNA's template). Returns the partition and the per-region contained units,
    per-boundary mass and per-sj mass — nothing of the slice rule is restated here."""
    from native._accumulator_reference import Accumulator, DepositOutcome, Partition

    from rigel.types import Strand

    sj = [(0, blocks[k][1], blocks[k + 1][0], int(Strand.POS)) for k in range(len(blocks) - 1)]
    part = Partition.from_region_bounds([bounds], sj=sj)
    offs = np.cumsum([0] + [e - s for s, e in blocks])
    L = int(bounds[-1] - bounds[0]) if chromosome else int(offs[-1])
    region = np.zeros(part.n_regions)
    boundary = np.zeros(part.n_boundaries)
    junction = np.zeros(part.n_sj)
    for w in np.flatnonzero(pmf):
        acc = Accumulator(part, max_fragment_length=pmf.shape[0] + 10)
        for s in range(0, L - int(w) + 1):
            if chromosome:
                g0, g1, introns = bounds[0] + s, bounds[0] + s + int(w), ()
            else:
                k0 = int(np.searchsorted(offs, s, side="right")) - 1
                k1 = int(np.searchsorted(offs, s + w - 1, side="right")) - 1
                g0 = blocks[k0][0] + s - int(offs[k0])
                g1 = blocks[k1][0] + s + int(w) - int(offs[k1])
                introns = tuple((blocks[k][1], blocks[k + 1][0]) for k in range(k0, k1))
            out = acc.deposit(
                0,
                g0,
                g1,
                observed_introns=introns,
                align_strand=Strand.POS,
                sj_strand=Strand.POS if introns else Strand.NONE,
            )
            assert out is DepositOutcome.DEPOSITED, (w, s, out)
        t = acc.tally
        region += pmf[w] * t.region_contained_count.sum(1)
        boundary += pmf[w] * (t.boundary_unspliced_mass + t.boundary_spliced_mass)
        junction += pmf[w] * t.sj_mass.sum(1)
    return part, region, boundary, junction


def _template_objects(bounds, blocks):
    """The template's pieces (its blocks cut at ``bounds``) in its own coordinates: per piece its region
    and length, per cut between consecutive pieces the genomic position of each side and the bases of
    template left and right of it."""
    B = np.asarray(bounds)
    region, length = [], []
    for gs, ge in blocks:
        inner = [int(x) for x in B[(B > gs) & (B < ge)]]
        edges = [gs, *inner, ge]
        for a, b in zip(edges[:-1], edges[1:]):
            region.append(int(np.searchsorted(B, a, side="right")) - 1)
            length.append(b - a)
    length = np.asarray(length, dtype=np.float64)
    through = np.cumsum(length)
    return np.asarray(region), length, through[:-1], through[-1] - through[:-1]


# 1-, 2- and 3-bp pieces beside long ones, and two junctions whose exons are cut into pieces as well
_BOUNDS = [0, 10, 11, 13, 16, 40, 100, 101, 108, 150, 180, 181, 230, 260]
_BLOCKS = [(10, 40), (100, 150), (180, 230)]


def test_each_share_is_what_the_deposit_rule_gives_the_object():
    """THE GATE, per object. Every placement of a spliced template goes through the reference accumulator;
    each piece's contained units must be its contained share, and each cut's mass — at a boundary or at a
    junction — its conserved share (:func:`conserved_cut_shares` at the template's reach), to 1e-12. Per
    object, because a ``1/K`` split over the cuts a fragment crosses also conserves the total. Every other
    object holds nothing, and the whole is the fl-marginal length."""
    pmf = _small_pmf()
    part, region, boundary, junction = _deposit_every_placement(_BOUNDS, _BLOCKS, pmf)
    reg, length, lo, hi = _template_objects(_BOUNDS, _BLOCKS)
    np.testing.assert_allclose(region[reg], contained_eff_length(length, pmf), rtol=0, atol=1e-12)
    left, right = conserved_cut_shares(pmf, length[:-1], length[1:], lo, hi)
    share = left + right
    is_boundary = reg[1:] == reg[:-1] + 1  # a boundary's id is its left region's
    got = np.where(is_boundary, boundary[reg[:-1]], 0.0)
    got[~is_boundary] = junction  # the sj ids run in genomic order, as the junctions do here
    assert part.n_sj == int((~is_boundary).sum())
    np.testing.assert_allclose(got, share, rtol=0, atol=1e-12)
    others = np.setdiff1d(np.arange(region.size), reg)
    assert np.all(region[others] == 0.0)
    assert np.all(np.delete(boundary, reg[:-1][is_boundary]) == 0.0)
    L = length.sum()
    w = np.arange(pmf.shape[0], dtype=np.float64)
    fl = float((pmf * np.maximum(L - w + 1.0, 0.0)).sum())
    assert region.sum() + boundary.sum() + junction.sum() == pytest.approx(fl, rel=1e-12)


def test_the_gdna_shares_are_what_the_deposit_rule_gives_each_boundary():
    """gDNA's template is the chromosome: every start on it, through the reference accumulator. Away from
    the chromosome's ends each boundary's mass is :func:`conserved_cut_shares` of its two flanking regions
    at ``UNBOUNDED_REACH`` — ``½E_f[min(a, w − 1)] + ½E_f[min(b, w − 1)]`` — and each region's contained
    units its contained share, beside pieces of 1–3 bp where one fragment crosses up to seven boundaries.
    PERTURBATION: the crossing support ``E_f[w − 1]``, a start once per boundary its fragment crosses,
    overstates every boundary beside a piece shorter than a fragment."""
    pmf = _small_pmf()
    bounds = [0, 200, 201, 203, 206, 207, 230, 231, 260, 460]
    part, region, boundary, _ = _deposit_every_placement(bounds, [(200, 260)], pmf, chromosome=True)
    length = np.diff(np.asarray(bounds, dtype=np.float64))
    left, right = conserved_cut_shares(pmf, length[:-1], length[1:], UNBOUNDED, UNBOUNDED)
    np.testing.assert_allclose(boundary, left + right, rtol=0, atol=1e-12)
    np.testing.assert_allclose(region, contained_eff_length(length, pmf), rtol=0, atol=1e-12)
    incidence = float(crossing_eff_length(pmf, np.array([UNBOUNDED]), np.array([UNBOUNDED]))[0])
    short = (length[:-1] < 40) | (length[1:] < 40)
    assert short.sum() >= 6
    assert np.all(boundary[short] < incidence - 1e-9)


def test_the_conserved_share_at_unbounded_reach_closes():
    """Where neither reach binds the share is ``½E_f[min(a, w − 1)] + ½E_f[min(b, w − 1)]``, computed here
    by summing over the pmf, and beside two pieces longer than every fragment it is the crossing support
    ``E_f[w − 1]`` exactly — no fragment crosses a second cut, so the conserved frame and the incidence
    frame agree."""
    pmf = _normal_pmf(200.0, 60.0, n=501)
    w = np.arange(pmf.shape[0], dtype=np.float64)
    a = np.array([1.0, 7.0, 150.0, 499.0, 5000.0])
    b = np.array([3000.0, 2.0, 80.0, 500.0, 5000.0])
    left, right = conserved_cut_shares(pmf, a, b, UNBOUNDED, UNBOUNDED)
    half_min = lambda x: 0.5 * float((pmf * np.minimum(x, np.maximum(w - 1.0, 0.0))).sum())  # noqa: E731
    np.testing.assert_allclose(left, [half_min(x) for x in a], rtol=1e-12)
    np.testing.assert_allclose(right, [half_min(x) for x in b], rtol=1e-12)
    S_e = float(crossing_eff_length(pmf, np.array([UNBOUNDED]), np.array([UNBOUNDED]))[0])
    assert left[-1] + right[-1] == pytest.approx(S_e, rel=1e-12)


def test_a_reach_that_excludes_its_own_piece_is_refused():
    with pytest.raises(ValueError, match="reach"):
        conserved_cut_shares(
            _small_pmf(), np.array([5.0]), np.array([5.0]), np.array([4.0]), np.array([9.0])
        )


def test_a_length_model_with_mass_at_zero_length_shares_without_a_nan():
    """A smoothed real-library model put mass at ``w = 0`` (LBX0588): a zero-length fragment crosses no
    cut, so the shares are finite and equal those of the same pmf with that mass removed, rescaled by the
    normalisation."""
    pmf = _normal_pmf(200.0, 60.0, n=501)
    with_zero = pmf.copy()
    with_zero[0] = 0.05
    a, b = np.array([30.0, 1500.0]), np.array([7.0, 300.0])
    lo, hi = np.array([30.0, 4000.0]), np.array([900.0, 300.0])
    clean = np.add(*conserved_cut_shares(pmf, a, b, lo, hi))
    dirty = np.add(*conserved_cut_shares(with_zero, a, b, lo, hi))
    assert np.all(np.isfinite(dirty))
    np.testing.assert_allclose(dirty, clean / 1.05, rtol=1e-12)
