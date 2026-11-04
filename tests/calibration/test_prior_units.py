"""`assemble_priors`: the gDNA count must be a conserved fragment count, not an incidence sum.

The EM reads ``gdna_count`` against its own count of the locus's fragments, so it must be in fragment units. A fragment
deposits on ``max(K, 1)`` objects (``K`` the contiguous boundaries it crosses), so a sum of per-object masses is a
fragment count only where every region is longer than every fragment, which many real regions are not. The region term
needs no conversion, since a contained fragment deposits on exactly one region; the crossing term is converted by the
accumulator's own conserved ``mass / count`` share at that boundary. That share is one number for both components,
because the accumulator sees them mixed, so on the gDNA count the mixture's share stands in for gDNA's own: a
compositional bias that a total over both components would conserve and hide, and that an equal-length substrate cannot
see. Until a per-component share exists, the biased value is the specified value and is asserted exactly. Every mass
here is the deposition law enumerated through the specification, so a failure is a defect, not noise.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    conserved_cut_shares,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.priors import assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.locus import Locus, MultiLocus
from native._accumulator_reference import Accumulator, Partition

_UNB = np.full(1, UNBOUNDED_REACH)


def _point_pmf(mean_len: int, max_len: int = 512) -> np.ndarray:
    """A point-mass fragment-length pmf at ``mean_len``. Exact arithmetic, no distributional slack."""
    p = np.zeros(max_len, dtype=np.float64)
    p[mean_len] = 1.0
    return p


def _enumerate(region_len, w):
    """Deposit every length-``w`` fragment at every start position, through the specification. Returns ``(n_fragments,
    contained_count, crossing_count, conserved_mass)``, the banks the accumulator would hold for a uniform field of unit
    density and a point-mass length.
    """

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
    """The share the accumulator itself would deposit on this tiling, by brute-force enumeration through the
    specification, read off its conserved-mass bank as ``mass / count``. For a point-mass pmf, every start at unit
    weight is the analytic uniform field the masses model, so the fixture is one self-consistent library.

    One share for both components, as ``assemble_priors`` uses: the accumulator cannot tell the populations apart, so
    the pooled share is ``(rho_g·m_g + rho_r·m_r) / (rho_g·c_g + rho_r·c_r)``, exact where the two components share a
    mean length and an approximation where they do not.
    """

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
    """A CalibrationResult for one reference tiled by ``region_len``, under a uniform field.

    Every bank is enumerated through the specification. The masses and the conserved share must describe one population:
    analytic infinite-chromosome masses (``rho·E[(len−w+1)+]``, ``rho·E[w−1]``) beside a finite-reference share would
    disagree on the fine tiling and read as a defect in ``assemble_priors``. For a point-mass pmf the enumeration is the
    analytic uniform field, so the fixture stays exact.

    The target is the conserved fragment count, ``rho·(S−w+1)`` on a finite reference of span ``S``, the fragments that
    fit; ``rho·span`` would count ``w−1`` start positions no fragment can occupy.
    """
    region_len = np.asarray(region_len, dtype=np.float64)
    n = region_len.shape[0]
    ne = max(n - 1, 0)
    # the opportunity arrays stay analytic: they are the divisors `gdna_eff_len` contracts against, not a population
    # statement
    a_g_region = contained_eff_length(region_len, pmf_g)
    a_g_boundary = np.full(ne, float(crossing_eff_length(pmf_g, _UNB, _UNB)[0]))
    # and every mass is what the specification deposits on this tiling
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
        gdna_boundary_conserved_len=np.add(
            *conserved_cut_shares(pmf_g, region_len[:-1], region_len[1:], _UNB, _UNB)
        ),
        gdna_frac_region=np.zeros_like(cont_g),
        rna_pos_frac_region=np.zeros_like(cont_g),
        rna_neg_frac_region=np.zeros_like(cont_g),
        gdna_frac_boundary=np.zeros_like(cross_g),
        rna_pos_frac_boundary=np.zeros_like(cross_g),
        rna_neg_frac_boundary=np.zeros_like(cross_g),
        gdna_density_global=rho_g,
        gdna_bulk_density=None,
        gdna_capture_saturation=1.0,
        gdna_capture_efficiency_region=np.ones(n),
        gdna_capture_efficiency_boundary=np.ones(ne),
        rna_sense_frac=0.9,
        n_regions=n,
        n_boundaries=ne,
        n_sj=0,
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
            span_bp=span,
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
    """``(truth_g, pred_g, q_g, q_r)`` for one tiling, re-derived from the reference accumulator's own banks rather than
    read back off ``assemble_priors``:

    * **truth**: the conserved fragment count gDNA deposited, ``SUM contained_g + SUM mass_g`` from gDNA's own
      conserved-mass bank; on a finite reference of span ``S`` this is ``rho_g·(S − w_g + 1)``, asserted separately.
    * **pred**: what a single pooled ``q`` produces: the same contained term plus gDNA's own crossing count rescaled by
      the mixture's share.

    ``q_c = mass_c / count_c`` per boundary is each component's own share, and ``pred == truth`` exactly where the two
    agree.
    """
    cont_g, cross_g, mass_g = _enumerate(tiling, int(np.argmax(pmf_g)))[1:]
    _cont_r, cross_r, mass_r = _enumerate(tiling, int(np.argmax(pmf_r)))[1:]
    q_pooled = _mass_per_crossing(tiling, rho_g, rho_r, pmf_g, pmf_r)
    ones = np.ones_like(q_pooled)
    q_g = np.divide(mass_g, cross_g, out=ones.copy(), where=cross_g > 0)
    q_r = np.divide(mass_r, cross_r, out=ones.copy(), where=cross_r > 0)
    return (
        rho_g * (cont_g.sum() + mass_g.sum()),
        rho_g * (cont_g.sum() + (cross_g * q_pooled).sum()),
        q_g,
        q_r,
    )


# --- T1: the conserved count ----------------------------------------------------------------------

# 1200 bp of reference, tiled four ways. The library is identical in all four; only the bookkeeping grid moves. The 100
# bp tiling is finer than the 200 bp RNA fragment, as most human regions are, and is where the raw incidence sum
# diverges most.
_SPAN = 1200
_MU_G, _MU_R = 50, 200
_TILINGS = {
    "coarse (1 x 1200)": [1200],
    "medium (3 x 400)": [400] * 3,
    "fine   (12 x 100)": [100] * 12,
    "ragged (mixed)": [37, 400, 63, 300, 1, 199, 200],
}
# The tilings where every region exceeds both fragment lengths, so ``min(w−1, flank) == w−1`` on every boundary and
# ``q_g == q_r == 1``: the pooled share is gDNA's own and the gDNA count is exact.
_SHARES_AGREE = ["coarse (1 x 1200)", "medium (3 x 400)"]


@pytest.mark.parametrize("name", _SHARES_AGREE)
def test_the_gdna_count_is_its_true_fragment_count_where_the_two_shares_agree(name):
    """Where the pooled share is gDNA's own (``q_g == q_r``) the gDNA count is exact: the gDNA fragments that fit,
    ``rho_g·(S−mu_g+1)``. This separates it from ``rho_g·span_bp``, which counts the ``w−1`` start positions no fragment
    can occupy. It cannot catch a raw incidence sum, since ``q`` is 1 on these two tilings; the finer tilings' tests do.
    """
    rho_g, rho_r = 0.03, 0.05
    tiling = _TILINGS[name]
    _, _, q_g, q_r = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R)
    )
    np.testing.assert_allclose(q_g, q_r, rtol=1e-12)  # the precondition, asserted not assumed
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    np.testing.assert_allclose(p.gdna_count, [rho_g * (_SPAN - _MU_G + 1)], rtol=1e-9)


@pytest.mark.parametrize(
    ("name", "gdna_bias"), [("fine   (12 x 100)", -0.199), ("ragged (mixed)", -0.026)]
)
def test_the_gdna_count_carries_exactly_the_pooled_share_bias(name, gdna_bias):
    """Where the two shares disagree the gDNA count is wrong, and this pins the wrong value exactly.
    ``boundary_mass_per_crossing`` is the mixture's share, and gDNA's crossings are rescaled by it. With the library
    unchanged and only the bookkeeping grid moved, the gDNA count goes:

    ==================  ==========
    tiling              gDNA
    ==================  ==========
    coarse / medium     exact fine (12 x 100)     −19.9 % ragged (mixed)      −2.6 %
    ==================  ==========

    gDNA is the shorter component here (50 bp against 200), so its own ``q_g`` is the larger and the pooled share drags
    it down; reversing the lengths reverses the sign (the sweep below). The biased value is the specified value until a
    per-component ``q`` exists, so it is asserted to 1e-9, along with the column above.
    """
    rho_g, rho_r = 0.03, 0.05
    tiling = _TILINGS[name]
    truth_g, pred_g, q_g, q_r = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R)
    )
    assert not np.allclose(q_g, q_r), "fixture no longer separates the two shares"
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    np.testing.assert_allclose(p.gdna_count, [pred_g], rtol=1e-9)
    # the bias is in the expected direction and of the expected size
    assert float(p.gdna_count[0] / truth_g - 1.0) == pytest.approx(gdna_bias, abs=5e-4)


# --- T2: the length sweep -------------------------------------------------------------------------


@pytest.mark.parametrize("mu_g", [50, 100, 150, 200, 300, 400])
def test_the_gdna_count_moves_with_the_length_ratio_by_exactly_the_pooled_share(mu_g):
    """Fixed true g:r, sweeping the two components' mean lengths against each other: the gDNA count moves either way,
    and this pins where it lands. Swept in both directions (``mu_g`` from 0.25x to 2x the RNA mean), since RNA need not
    be longer than gDNA. A shorter gDNA is under-counted and a longer one over-counted, because a longer fragment is
    censored harder by a 100 bp flank and carries the smaller share.

    The locus-level distortion is not ``q_pooled/q_g``: only the crossing term passes through the share, the contained
    term is already a fragment count, and the pooled share itself depends on the mixture. So the expectation is
    recomputed per arm, and the test checks the direction plus the recomputed value.
    """
    rho_g, rho_r = 0.02, 0.06
    tiling = _TILINGS["fine   (12 x 100)"]
    truth_g, pred_g, _, _ = _truth_and_prediction(
        tiling, rho_g, rho_r, _point_pmf(mu_g), _point_pmf(_MU_R)
    )
    p = _priors_for(tiling, rho_g, rho_r, _point_pmf(mu_g), _point_pmf(_MU_R))
    count = float(p.gdna_count[0])
    assert count == pytest.approx(pred_g, rel=1e-9), (
        f"gDNA count is {count:.6f} at mu_g={mu_g}, not the pooled-share {pred_g:.6f}"
    )
    if mu_g < _MU_R:
        assert count < truth_g, "a SHORTER gDNA must be under-counted"
    elif mu_g > _MU_R:
        assert count > truth_g, "a LONGER gDNA must be over-counted"


def test_the_gdna_count_is_exact_where_the_two_components_share_a_length():
    """At equal lengths the bias is exactly zero, so a panel built that way cannot measure it; the ladder gives its two
    components equal lengths by design, so measuring the gapped case needs a length-gap panel. Asserted at two very
    different mixtures: when ``q_g == q_r`` the pooled share equals both regardless of ``phi``.
    """
    tiling = _TILINGS["fine   (12 x 100)"]
    for rho_g, rho_r in ((0.02, 0.06), (0.05, 0.01)):
        p = _priors_for(tiling, rho_g, rho_r, _point_pmf(_MU_R), _point_pmf(_MU_R))
        count = float(p.gdna_count[0])
        fit = rho_g * (_SPAN - _MU_R + 1)
        assert count == pytest.approx(fit, rel=1e-12), (
            f"equal lengths must be unbiased; got {count:.9f} against {fit:.9f} at rho {rho_g}/{rho_r}"
        )


# --- T4: a component with no opportunity, and stray mass on it, never reach the gDNA count ---------


def test_zero_rna_opportunity_leaves_the_gdna_count_exact():
    """Every region is shorter than one RNA fragment and the RNA crossing count is zeroed; with ``rho_r = 0`` the
    library is one component, so the pooled share is gDNA's own and the count is exact.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4  # every region < 400 bp ⇒ contained_eff_length(RNA) == 0
    cal = _uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r)
    cal = _zero_rna_opportunity(cal)
    p = assemble_priors(cal, _regions_tiling(tiling), _one_locus(int(np.sum(tiling))))
    np.testing.assert_allclose(p.gdna_count, [0.03 * (np.sum(tiling) - 20 + 1)], rtol=1e-9)


def test_stray_rna_mass_on_a_zero_opportunity_object_never_reaches_the_gdna_count():
    """Calibration's RNA count does not enter the EM, so RNA mass on an object with no RNA opportunity must leave the
    gDNA count where it was. That configuration is ordinary: ``contained_eff_length`` is 0 wherever a region is shorter
    than the component's shortest fragment, and ``f_g`` is an inference. That a count keeps stray mass of its own
    component is the next test.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4
    cal = _zero_rna_opportunity(_uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r))
    # unlike the test above, real mass sits on the zero-opportunity RNA objects
    stray = dataclasses.replace(cal, count_rna_region=np.full(4, 2.5))
    regions, loci = _regions_tiling(tiling), _one_locus(int(np.sum(tiling)))
    p = assemble_priors(stray, regions, loci)
    # the gDNA side, which has opportunity everywhere, is untouched by the stray RNA mass
    np.testing.assert_allclose(p.gdna_count, [0.03 * (np.sum(tiling) - 20 + 1)], rtol=1e-9)


def test_gdna_mass_on_a_zero_opportunity_region_counts_as_a_count_has_no_divisor():
    """A count is a fragment count, not a density: gDNA mass on a region with zero contained support is still gDNA
    fragments in the locus and may not be dropped for want of a divisor. That the same stray mass never reaches
    ``gdna_eff_len`` is tested in ``test_priors.py``.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4
    cal = _uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r)
    regions, loci = _regions_tiling(tiling), _one_locus(int(np.sum(tiling)))
    base = assemble_priors(cal, regions, loci)
    stray = dataclasses.replace(
        cal,
        gdna_region_eff_len=np.zeros_like(cal.gdna_region_eff_len),
        count_gdna_region=np.asarray(cal.count_gdna_region, np.float64) + 2.5,
    )
    p = assemble_priors(stray, regions, loci)
    np.testing.assert_allclose(p.gdna_count - base.gdna_count, [4 * 2.5], rtol=1e-12)


def _zero_rna_opportunity(cal: CalibrationResult) -> CalibrationResult:
    """Remove every RNA crossing count."""
    return dataclasses.replace(cal, count_rna_boundary=np.zeros_like(cal.count_rna_boundary))
