"""`assemble_priors`: the EM pseudocounts must be CONSERVED FRAGMENT COUNTS, not incidence sums.

``gdna_prior_count`` / ``rna_prior_count`` reach the EM as additive pseudocounts in fragment units,
but a fragment deposits on ``max(K, 1)`` objects, ``K`` being the contiguous boundaries it crosses —
so a sum of per-object masses is a fragment count only where every region is longer than every
fragment, and a large share of real regions is not. The region term needs no conversion, a contained
fragment depositing on exactly one region; only the crossing term is converted, by the accumulator's
own conserved ``mass / count`` share at that boundary. That share is ONE number for both components,
because the accumulator sees them mixed, and it is a compositional bias rather than an approximation:
rescaling both by the mixture's share conserves the locus TOTAL exactly while tilting the g:r SPLIT —
invisible to a total-mass gate (`TRAPS: conservation-misses-mis-attribution`), and equally invisible
on an equal-length substrate (`TRAPS: an-equal-length-panel-defeats-the-lift`). Until a per-component
share is built the biased value IS the specified value, so it is asserted exactly. Every mass below
is the deposition law enumerated through the SPECIFICATION, so a failure is a defect, never noise.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.priors import assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus

_UNB = np.full(1, UNBOUNDED_REACH)


def _point_pmf(mean_len: int, max_len: int = 512) -> np.ndarray:
    """A point-mass fragment-length pmf at ``mean_len``. Exact arithmetic, no distributional slack."""
    p = np.zeros(max_len, dtype=np.float64)
    p[mean_len] = 1.0
    return p


def _enumerate(region_len, w):
    """Deposit EVERY length-``w`` fragment at EVERY start position, through the SPECIFICATION itself.

    Returns ``(n_fragments, contained_count, crossing_count, conserved_mass)`` — the exact banks the
    accumulator would hold for a uniform field of unit density and a point-mass length.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
    from native._accumulator_reference import Accumulator, Partition

    region_bounds = np.concatenate(
        [[0.0], np.cumsum(np.asarray(region_len, dtype=np.float64))]
    ).astype(int)
    partition = Partition.from_region_bounds(
        [region_bounds.tolist()], region_types=[[0] * (len(region_bounds) - 1)]
    )
    acc = Accumulator(partition, max_fragment_length=10**6)
    n = 0
    for start in range(0, int(region_bounds[-1]) - int(w) + 1):
        acc.deposit(0, start, start + int(w))
        n += 1
    t = acc.tally
    return (
        n,
        np.asarray(t.region_contained_count, np.int64).sum(axis=1).astype(np.float64),
        np.asarray(t.boundary_unspliced_count, np.int64).sum(axis=1).astype(np.float64),
        np.asarray(t.boundary_unspliced_mass, np.float64),
    )


def _mass_per_crossing(region_len, rho_g, rho_r, pmf_g, pmf_r) -> np.ndarray:
    """The share the ACCUMULATOR itself would deposit on this tiling, by brute-force enumeration.

    Not derived from the answer these tests check. Every fragment of each component's length is
    deposited at every start position through the SPECIFICATION, and the share is read off its own
    conserved-mass bank as ``mass / count``. Enumerating every start at unit weight IS the analytic
    uniform field the masses above model — for a point-mass pmf the two agree exactly — so the fixture
    stays one self-consistent library rather than two models that have to be argued equal
    (``TRAPS: a-test-that-redefines``).

    ONE share for BOTH components, which is what ``assemble_priors`` uses. The accumulator sees the
    two populations mixed and cannot tell them apart, so the pooled share is
    ``(rho_g·m_g + rho_r·m_r) / (rho_g·c_g + rho_r·c_r)``. Where the two components share a mean length
    that is exact; where they do not it is an approximation, and these tests are where that shows.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
    from native._accumulator_reference import (
        Accumulator,
        Partition,
    )

    region_bounds = np.concatenate(
        [[0.0], np.cumsum(np.asarray(region_len, dtype=np.float64))]
    ).astype(int)
    n_boundaries = max(len(region_bounds) - 2, 0)
    mass = np.zeros(n_boundaries, dtype=np.float64)
    count = np.zeros(n_boundaries, dtype=np.float64)
    for rho, pmf in ((rho_g, pmf_g), (rho_r, pmf_r)):
        w = int(np.argmax(pmf))
        partition = Partition.from_region_bounds(
            [region_bounds.tolist()], region_types=[[0] * (len(region_bounds) - 1)]
        )
        acc = Accumulator(partition, max_fragment_length=max(1000, w + 1))
        for start in range(0, int(region_bounds[-1]) - w + 1):
            acc.deposit(0, start, start + w)
        t = acc.tally
        mass += rho * np.asarray(t.boundary_unspliced_mass, np.float64)
        count += rho * np.asarray(t.boundary_unspliced_count, np.int64).sum(axis=1).astype(
            np.float64
        )
    out = np.ones(n_boundaries, dtype=np.float64)
    np.divide(mass, count, out=out, where=count > 0)
    return out


def _uniform_library(region_len, rho_g, rho_r, pmf_g, pmf_r) -> CalibrationResult:
    """A CalibrationResult for ONE reference tiled by ``region_len``, under a UNIFORM field.

    EVERY BANK IS ENUMERATED THROUGH THE SPECIFICATION, not evaluated analytically. The masses and the
    conserved share must describe ONE population or the fixture is internally inconsistent, and that
    inconsistency is not hypothetical: taking the masses from the analytic infinite-chromosome forms
    ``rho·E[(len−w+1)+]`` and ``rho·E[w−1]`` while the share comes from a finite-reference enumeration
    puts the two well apart on the fine tiling, and it reads as a defect in ``assemble_priors``.

    For a point-mass pmf, enumerating every start at unit weight IS the analytic uniform field, so
    nothing is lost — the fixture is still exact, deterministic arithmetic with no distributional
    slack.

    The target is the CONSERVED FRAGMENT COUNT, which on a finite reference of span ``S`` is
    ``rho·(S−w+1)``, the number of fragments that fit. ``rho·span`` — what dividing out the
    opportunity and re-integrating over the genomic span produces — counts ``w−1`` start positions no
    fragment can occupy.
    """
    region_len = np.asarray(region_len, dtype=np.float64)
    n = region_len.shape[0]
    ne = max(n - 1, 0)
    # the OPPORTUNITY arrays stay analytic — they are the divisors `gdna_eff_len` contracts against,
    # and they are not a population statement
    a_g_region = contained_eff_length(region_len, pmf_g)
    a_r_region = contained_eff_length(region_len, pmf_r)
    a_g_boundary = np.full(ne, float(crossing_eff_length(pmf_g, _UNB, _UNB)[0]))
    a_r_boundary = np.full(ne, float(crossing_eff_length(pmf_r, _UNB, _UNB)[0]))
    # ...and every MASS is what the specification actually deposits on this tiling
    _n_g, cont_g, cross_g, _m_g = _enumerate(region_len, int(np.argmax(pmf_g)))
    _n_r, cont_r, cross_r, _m_r = _enumerate(region_len, int(np.argmax(pmf_r)))
    return CalibrationResult(
        count_gdna_region=rho_g * cont_g,
        count_rna_region=rho_r * cont_r,
        count_gdna_boundary=rho_g * cross_g,
        count_rna_boundary=rho_r * cross_r,
        count_rna_spliced_boundary=np.zeros(ne, dtype=np.float64),
        boundary_mass_per_crossing=_mass_per_crossing(region_len, rho_g, rho_r, pmf_g, pmf_r),
        count_rna_sj=np.zeros(0, dtype=np.float64),
        boundary_spliced_mass_per_crossing=np.ones_like(
            _mass_per_crossing(region_len, rho_g, rho_r, pmf_g, pmf_r)
        ),
        sj_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_region_eff_len=a_g_region,
        gdna_boundary_eff_len=a_g_boundary,
        rna_region_eff_len=a_r_region,
        rna_boundary_eff_len=a_r_boundary,
        gdna_frac_region=np.zeros_like(cont_g),
        rna_pos_frac_region=np.zeros_like(cont_g),
        rna_neg_frac_region=np.zeros_like(cont_g),
        gdna_frac_boundary=np.zeros_like(cross_g),
        rna_pos_frac_boundary=np.zeros_like(cross_g),
        rna_neg_frac_boundary=np.zeros_like(cross_g),
        gdna_density_global=rho_g,
        gdna_reference_density=None,
        gdna_reference_members=0,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_regions=n,
        n_boundaries=ne,
        n_sj=0,
        config=CalibrationConfig(),
    )


def _regions_tiling(region_len) -> RegionArrays:
    """One reference tiled contiguously from 0 by ``region_len``."""
    region_len = np.asarray(region_len, dtype=np.int64)
    ends = np.cumsum(region_len)
    starts = ends - region_len
    n = region_len.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        # exon on both strands: the region must survive the locus projection's intergenic drop
        signature=np.full(n, 0b0000_0011, dtype=np.uint8),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=region_len.astype(np.float64),
        ref_offsets=np.array([0, n], dtype=np.int32),
        n_refs=1,
    )


def _one_locus(span: int) -> list[MultiLocus]:
    return [
        MultiLocus(
            multi_locus_id=0,
            transcript_indices=np.array([], dtype=np.int32),
            unit_indices=np.array([], dtype=np.int32),
            gdna_span=span,
            loci=(Locus(ref="0", ref_id=0, start=0, end=span),),
        )
    ]


def _priors_for(tiling, rho_g, rho_r, pmf_g, pmf_r):
    span = int(np.sum(tiling))
    return assemble_priors(
        _uniform_library(tiling, rho_g, rho_r, pmf_g, pmf_r),
        _regions_tiling(tiling),
        _one_locus(span),
    )


def _truth_and_prediction(tiling, rho_g, rho_r, pmf_g, pmf_r):
    """``(truth_g, truth_r, pred_g, pred_r, q_g, q_r)`` for one tiling — from the SPECIFICATION only.

    Neither number is read back off ``assemble_priors`` (`TRAPS: a-test-that-redefines`). Both
    are re-derived from the reference accumulator's own banks:

    * **truth** — the conserved fragment count each component really deposited, ``SUM contained_c +
      SUM mass_c``, where ``mass_c`` is that component's OWN conserved-mass bank. On a finite reference
      of span ``S`` this is ``rho_c·(S − w_c + 1)``, the fragments that FIT, and the closed form is
      asserted against it separately.
    * **pred** — what a single POOLED ``q`` must produce: the same contained term, plus each
      component's own crossing COUNT rescaled by the MIXTURE's share. This is the pooling result
      stated as arithmetic, and its only content beyond ``truth`` is that ``q_c`` has been replaced by
      ``q_pooled``.

    ``q_c = mass_c / count_c`` per boundary is each component's own share, and ``pred == truth`` exactly
    where the two agree.
    """
    cont_g, cross_g, mass_g = _enumerate(tiling, int(np.argmax(pmf_g)))[1:]
    cont_r, cross_r, mass_r = _enumerate(tiling, int(np.argmax(pmf_r)))[1:]
    q_pooled = _mass_per_crossing(tiling, rho_g, rho_r, pmf_g, pmf_r)
    ones = np.ones_like(q_pooled)
    q_g = np.divide(mass_g, cross_g, out=ones.copy(), where=cross_g > 0)
    q_r = np.divide(mass_r, cross_r, out=ones.copy(), where=cross_r > 0)
    return (
        rho_g * (cont_g.sum() + mass_g.sum()),
        rho_r * (cont_r.sum() + mass_r.sum()),
        rho_g * (cont_g.sum() + (cross_g * q_pooled).sum()),
        rho_r * (cont_r.sum() + (cross_r * q_pooled).sum()),
        q_g,
        q_r,
    )


# --- T1: the conserved count ----------------------------------------------------------------------

# 1200 bp of reference, tiled four ways. The library is IDENTICAL in all four; only the bookkeeping
# grid moves. The 100 bp tiling is finer than the 200 bp RNA fragment, which is where most human
# regions live and where the raw incidence sum diverges hardest.
_SPAN = 1200
_MU_G, _MU_R = 50, 200
_TILINGS = {
    "coarse (1 x 1200)": [1200],
    "medium (3 x 400)": [400] * 3,
    "fine   (12 x 100)": [100] * 12,
    "ragged (mixed)": [37, 400, 63, 300, 1, 199, 200],
}
# The tilings where every region exceeds BOTH fragment lengths, so ``min(w−1, flank) == w−1`` on every
# boundary and ``q_g == q_r == 1``: the pooled share is then each component's own and the split is exact.
_SHARES_AGREE = ["coarse (1 x 1200)", "medium (3 x 400)"]
# Not 1e-9: the conserved-mass bank is fixed-point at 2^-32 per fragment, so a 1,001-fragment total
# carries ~2e-11 of relative rounding. Anything above 1e-10 here would be a real error.
_FIXED_POINT_RTOL = 1e-8


@pytest.mark.parametrize("name", list(_TILINGS))
def test_the_total_prior_is_the_true_fragment_count_on_every_tiling(name):
    """The conservation gate. The same physical library, re-tiled, deposits the same TOTAL — and it is
    the right total: the number of fragments that FIT, ``rho_g·(S−mu_g+1) + rho_r·(S−mu_r+1)``.

    A raw incidence sum grows as the tiling is refined, every new boundary adding a crossing to every
    fragment that spans it, so it fails here by construction; a ``rho_c·span_bp`` form fails too, by
    the ``w−1`` start positions no fragment can occupy.

    This gate is BLIND to the defect the next two tests exist for
    (`TRAPS: conservation-misses-mis-attribution`): rescaling both components by one pooled share
    conserves the total exactly while tilting the split, so on the fine tiling this passes while the
    gDNA side alone is a fifth low. Never read it as "the prior is right".
    """
    rho_g, rho_r = 0.03, 0.05
    p = _priors_for(_TILINGS[name], rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    total = float(p.gdna_prior_count[0] + p.rna_prior_count[0])
    fit = rho_g * (_SPAN - _MU_G + 1) + rho_r * (_SPAN - _MU_R + 1)
    assert total == pytest.approx(fit, rel=_FIXED_POINT_RTOL), (
        f"{name}: total {total:.9f} against {fit:.9f} fragments that fit"
    )
    # and it is NOT rho·span, which is the same on every tiling and a few per cent too big
    assert not np.isclose(total, (rho_g + rho_r) * _SPAN, rtol=1e-3)


@pytest.mark.parametrize("name", _SHARES_AGREE)
def test_each_component_is_its_true_fragment_count_where_the_two_shares_agree(name):
    """And where the pooled share IS each component's own (``q_g == q_r``), the SPLIT is exact too.

    Stronger than the total alone: a form uniformly wrong by a constant factor passes the
    conservation gate and fails this one. Restricted to two of the four tilings on purpose — that is
    the substrate condition an equal-length panel satisfies by construction, which is exactly why
    such a panel cannot see the bias (`TRAPS: an-equal-length-panel-defeats-the-lift`).
    """
    rho_g, rho_r = 0.03, 0.05
    tiling = _TILINGS[name]
    _, _, _, _, q_g, q_r = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R)
    )
    np.testing.assert_allclose(q_g, q_r, rtol=1e-12)  # the precondition, asserted not assumed
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    np.testing.assert_allclose(p.gdna_prior_count, [rho_g * (_SPAN - _MU_G + 1)], rtol=1e-9)
    np.testing.assert_allclose(p.rna_prior_count, [rho_r * (_SPAN - _MU_R + 1)], rtol=1e-9)


@pytest.mark.parametrize(
    ("name", "gdna_bias"), [("fine   (12 x 100)", -0.199), ("ragged (mixed)", -0.026)]
)
def test_the_split_carries_exactly_the_POOLED_SHARE_bias(name, gdna_bias):
    """Where the two shares disagree the split is wrong, and this pins the wrong value exactly.

    The accumulator cannot tell the two populations apart, so ``boundary_mass_per_crossing`` is the
    MIXTURE's share and both components are rescaled by it. With the library physically unchanged and
    only the bookkeeping grid moved, the split goes:

    ==================  ==========  ==========  ==========
    tiling              gDNA        RNA         total
    ==================  ==========  ==========  ==========
    coarse / medium     exact       exact       exact
    fine (12 x 100)     −19.9 %     +13.7 %     exact
    ragged (mixed)      −2.6 %      +1.8 %      exact
    ==================  ==========  ==========  ==========

    gDNA is the SHORTER component here (50 bp against 200), so its own ``q_g`` is the larger and the
    pooled share drags it DOWN. Reverse the lengths and the sign reverses — that is the sweep below.

    The biased value is the SPECIFIED value until a per-component ``q`` is built, so it is asserted
    to 1e-9 rather than tolerated with a loose bound, and the gDNA column above is asserted with it
    so the bias can neither drift nor silently vanish.
    """
    rho_g, rho_r = 0.03, 0.05
    tiling = _TILINGS[name]
    truth_g, truth_r, pred_g, pred_r, q_g, q_r = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R)
    )
    assert not np.allclose(q_g, q_r), "fixture no longer separates the two shares"
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    np.testing.assert_allclose(p.gdna_prior_count, [pred_g], rtol=1e-9)
    np.testing.assert_allclose(p.rna_prior_count, [pred_r], rtol=1e-9)
    # ...and the bias is real, in the recorded direction and of the recorded size
    assert float(p.gdna_prior_count[0] / truth_g - 1.0) == pytest.approx(gdna_bias, abs=5e-4)
    assert p.rna_prior_count[0] > truth_r  # the longer component absorbs what the shorter one lost


# --- T2: the length sweep -------------------------------------------------------------------------


@pytest.mark.parametrize("mu_g", [50, 100, 150, 200, 300, 400])
def test_the_prior_ratio_moves_with_the_length_ratio_by_exactly_the_pooled_share(mu_g):
    """The composition test: fixed true g:r, sweeping the two components' mean lengths against each
    other. The prior's ratio MOVES by a factor either way, and this pins where it lands.

    Swept in BOTH directions (``mu_g`` from 0.25x to 2x the RNA mean), because there is no rule that
    RNA is longer than gDNA and assuming one is how a tool overfits to cfRNA. The direction is the
    finding: the SHORTER component is under-called and the longer one over-called, because a longer
    fragment is censored harder by a 100 bp flank and so carries the smaller share.

    The distortion is NOT ``q_r/q_g`` at the locus level, and mixture-independence does not survive
    contact with contained mass. Only the CROSSING term passes through the share; the contained term
    is already a fragment count and is untouched, so the locus-level tilt is diluted by each
    component's contained fraction — and since the pooled share itself depends on the mixture, the
    dilution does too. That is why the expectation is recomputed per arm rather than stated as a
    constant, and why the gate is the direction plus the recomputed value.
    """
    rho_g, rho_r = 0.02, 0.06
    tiling = _TILINGS["fine   (12 x 100)"]
    _, _, pred_g, pred_r, _, _ = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(mu_g), _point_pmf(_MU_R)
    )
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(mu_g), _point_pmf(_MU_R))
    ratio = float(p.gdna_prior_count[0] / p.rna_prior_count[0])
    assert ratio == pytest.approx(pred_g / pred_r, rel=1e-9), (
        f"prior g:r is {ratio:.6f} at mu_g={mu_g}, not the pooled-share {pred_g / pred_r:.6f}"
    )
    true_ratio = rho_g / rho_r
    if mu_g < _MU_R:
        assert ratio < true_ratio, "the SHORTER component must be under-called"
    elif mu_g > _MU_R:
        assert ratio > true_ratio, "the LONGER component must be over-called"


def test_the_ratio_IS_exact_where_the_two_components_share_a_length():
    """At equal lengths the bias is exactly zero, which is why a panel built that way cannot measure it
    (`TRAPS: an-equal-length-panel-defeats-the-lift`). The ladder gives its two components equal
    fragment lengths by design, so this gate is the executable statement of the equal-length half and
    the gapped half needs a length-gap panel to be measured at all.

    Asserted at two very different mixtures, because "exact" here must not depend on the mixing
    ratio: when ``q_g == q_r`` the pooled share equals both regardless of ``phi``.
    """
    tiling = _TILINGS["fine   (12 x 100)"]
    for rho_g, rho_r in ((0.02, 0.06), (0.05, 0.01)):
        p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_R), _point_pmf(_MU_R))
        ratio = float(p.gdna_prior_count[0] / p.rna_prior_count[0])
        assert ratio == pytest.approx(rho_g / rho_r, rel=1e-12), (
            f"equal lengths must be unbiased; got {ratio:.9f} at rho {rho_g}/{rho_r}"
        )


# --- T4: zero opportunity emits nothing, never a floored division -----------------------------------


def test_zero_rna_opportunity_gives_zero_rna_prior():
    """An object with no opportunity for a component must emit NOTHING at zero precision, never a
    floored division. Every region here is shorter than one RNA fragment and the RNA crossing
    opportunity is zeroed, so the RNA support is identically 0.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4  # every region < 400 bp ⇒ contained_eff_length(RNA) == 0
    cal = _uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r)
    cal = _zero_rna_opportunity(cal)
    p = assemble_priors(cal, _regions_tiling(tiling), _one_locus(int(np.sum(tiling))))
    assert np.all(np.isfinite(p.rna_prior_count))
    np.testing.assert_allclose(p.rna_prior_count, [0.0])
    # and the gDNA side is EXACT here, not merely finite: with rho_r = 0 the library is one
    # component, so the pooled share IS the gDNA's own and there is nothing to bias the split.
    np.testing.assert_allclose(p.gdna_prior_count, [0.03 * (np.sum(tiling) - 20 + 1)], rtol=1e-9)


def test_mass_on_a_zero_opportunity_object_STILL_COUNTS_because_a_count_has_no_divisor():
    """Stray mass on a zero-opportunity object is KEPT by the prior, and this is half of a pair.

    ``mass > 0`` with ``support == 0`` is an ordinary configuration rather than a corner:
    ``contained_eff_length`` is exactly 0 wherever a region is shorter than that component's shortest
    fragment, which is a large minority of regions in a real index for both components, and the
    solver can still put mass there because ``f_g`` is an inference and not a fact.

    Dropping it would be right for a RATE — ``rho = SUM m / SUM S`` with mass in the numerator and no
    exposure in the denominator inflates without bound — but the prior no longer divides by anything.
    ``mass_c_region[r]`` is ``f_c(r)·contained_count[r]``, and a contained fragment deposits on
    exactly one region, so the mass IS the count and dropping it would silently lose fragments the
    accumulator really deposited. The inflation is structurally unreachable here rather than guarded,
    and the assertion says so: the prior is exactly the deposited 4 × 2.5.

    The guard is still LIVE where a divisor still lives, the eff-length, and that half is
    `test_priors.test_stray_mass_on_a_zero_opportunity_boundary_is_dropped_from_the_eff_len`, which
    perturbs it and measures the bp it holds back. Do not delete one without the other: alone, either
    one reads as a rule about the whole file.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4
    cal = _zero_rna_opportunity(_uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r))
    # the difference from the test above: put REAL mass on the zero-opportunity RNA objects
    stray = dataclasses.replace(cal, count_rna_region=np.full(4, 2.5))
    regions, loci = _regions_tiling(tiling), _one_locus(int(np.sum(tiling)))
    p = assemble_priors(stray, regions, loci)
    assert np.all(np.isfinite(p.rna_prior_count)), "a floored divisor produced a non-finite prior"
    np.testing.assert_allclose(p.rna_prior_count, [4 * 2.5], rtol=1e-12)
    # and the gDNA side, which DOES have opportunity everywhere, is untouched by the stray RNA mass
    np.testing.assert_allclose(p.gdna_prior_count, [0.03 * (np.sum(tiling) - 20 + 1)], rtol=1e-9)


def _zero_rna_opportunity(cal: CalibrationResult) -> CalibrationResult:
    """Remove every RNA crossing opportunity, so the RNA support is identically 0 on all objects."""
    return dataclasses.replace(
        cal,
        rna_boundary_eff_len=np.zeros_like(cal.rna_boundary_eff_len),
        count_rna_boundary=np.zeros_like(cal.count_rna_boundary),
    )
