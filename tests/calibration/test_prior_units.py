"""assemble_priors — the EM pseudocounts must be CONSERVED FRAGMENT COUNTS, not object-incidence sums.

⭐ **THE DEFECT THESE TESTS PIN.** ``gdna_prior_count`` / ``rna_prior_count`` are handed to the EM as
**additive pseudocounts in fragment units** — ``G = n_gdna + a_g`` in ``apply_grouped_prior_update`` —
where ``n_gdna`` is a count of fragments. But a fragment deposits on ``max(K, 1)`` objects, ``K`` being
the number of contiguous lines it crosses, so summing per-object masses does NOT give a fragment count::

    incidences(w) = max( 1 , (w-1)/s )        for a partition of spacing s

Counts are conserved exactly where every node is longer than every fragment, and become a
**length-weighted** count where they are not — and 56.7 % of human nodes are shorter than one 200 bp
fragment. Because the weighting is by length, it does not cancel between two components with different
mean lengths: measured on the chr22 pilot, gDNA deposits 1.031 incidences per fragment and RNA ≈1.17,
so the prior's g:r ratio under-calls gDNA by 13–19 %.

⭐⭐ **THE FIX, AND IT IS A READ-OUT RATHER THAN A DERIVATION.** The node term is already a fragment
count — a contained fragment deposits on exactly one node. Only the crossing term is converted, by the
accumulator's own conserved ``mass / count`` at that line::

    prior_c = SUM_locus share · [ mass_c_node[r] + SUM_{e owned by r} mass_c_edge[e] · q[e] ]

    q[e] = edge_mass_per_crossing = [ min(w−1,a) + min(w−1,b) ] / 2(w−1)    under a uniform field

⛔ **The predecessor rule was ``rho_c = SUM m / SUM A ; prior_c = rho_c · span_bp``**, and these tests
used to target it. It reached fragment units by dividing out the opportunity and re-integrating over the
genomic span, which on a finite reference of span ``S`` counts the ``w−1`` start positions no fragment
can occupy: the truth is ``rho·(S−w+1)``, not ``rho·S``. Both numbers appear below and they differ by
4.1 % on this fixture — small, systematic, and in fragment units it is simply wrong.

⛔⛔ **ONE ``q`` FOR TWO COMPONENTS IS A COMPOSITIONAL BIAS, AND THESE TESTS PIN IT AS A NUMBER RATHER
THAN TOLERATE IT.** The accumulator sees the two populations mixed, so ``q`` is the MIXTURE's share.
Rescaling both components by it conserves the locus TOTAL exactly to the fragment while tilting the g:r
SPLIT — which is why `test_the_total_prior_is_the_true_fragment_count_on_every_tiling` passes on every
tiling and `test_each_component_is_its_true_fragment_count_where_the_two_shares_agree` is restricted to
the tilings where ``q_g == q_r``. A total-mass gate cannot see this
(`TRAPS: conservation-misses-mis-attribution`), and neither can a substrate where the two components
share a length distribution (`TRAPS: an-equal-length-panel-defeats-the-lift`) — the sweep below reads
EXACTLY 1.000000 at ``mu_g == mu_r`` and 0.56–1.59× away from it. The repair is a per-component ``q``
and it is not built; until it is, the biased value is the specified value and is asserted to 1e-9.

⛔ **These tests are deterministic, not simulated.** Every mass below is the accumulator's deposition law
evaluated exactly through the SPECIFICATION, so a failure is a defect and never noise. That is a
deliberate departure from the plan's original end-to-end phrasing of T1/T2: rebuilding an index with
extra cuts also moves the transcript set, the reach and every effective length, so it would not isolate
the partition. The end-to-end conservation check against ``node_start_count`` is T3, and it lives in
`scripts/design/prior_units_check.py`.
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


def _enumerate(node_len, w):
    """Deposit EVERY length-``w`` fragment at EVERY start position, through the SPECIFICATION itself.

    Returns ``(n_fragments, contained_count, crossing_count, conserved_mass)`` — the exact banks the
    accumulator would hold for a uniform field of unit density and a point-mass length.
    """
    import sys
    from pathlib import Path as _P

    sys.path.insert(0, str(_P(__file__).resolve().parents[1]))
    from native._accumulator_reference import Accumulator, Partition

    cuts = np.concatenate([[0.0], np.cumsum(np.asarray(node_len, dtype=np.float64))]).astype(int)
    partition = Partition.from_cuts([cuts.tolist()], node_types=[[0] * (len(cuts) - 1)])
    acc = Accumulator(partition, max_fragment_length=10**6)
    n = 0
    for start in range(0, int(cuts[-1]) - int(w) + 1):
        acc.deposit(0, start, start + int(w))
        n += 1
    t = acc.tally
    return (
        n,
        np.asarray(t.node_contained_count, np.int64).sum(axis=1).astype(np.float64),
        np.asarray(t.edge_unspliced_count, np.int64).sum(axis=1).astype(np.float64),
        np.asarray(t.edge_unspliced_mass, np.float64),
    )


def _mass_per_crossing(node_len, rho_g, rho_r, pmf_g, pmf_r) -> np.ndarray:
    """The share the ACCUMULATOR itself would deposit on this tiling, by brute-force enumeration.

    ⛔ **Not derived from the answer these tests check.** Every fragment of each component's length is
    deposited at every start position through the SPECIFICATION, and the share is read off its own
    conserved-mass bank as ``mass / count``. Enumerating every start at unit weight IS the analytic
    uniform field the masses above model — for a point-mass pmf the two agree exactly — so the fixture
    stays one self-consistent library rather than two models that have to be argued equal
    (``TRAPS: a-test-that-redefines``).

    ⚠ **ONE share for BOTH components, which is what ``assemble_priors`` uses.** The accumulator sees the
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

    cuts = np.concatenate([[0.0], np.cumsum(np.asarray(node_len, dtype=np.float64))]).astype(int)
    n_edges = max(len(cuts) - 2, 0)
    mass = np.zeros(n_edges, dtype=np.float64)
    count = np.zeros(n_edges, dtype=np.float64)
    for rho, pmf in ((rho_g, pmf_g), (rho_r, pmf_r)):
        w = int(np.argmax(pmf))
        partition = Partition.from_cuts(
            [cuts.tolist()], node_types=[[0] * (len(cuts) - 1)]
        )
        acc = Accumulator(partition, max_fragment_length=max(1000, w + 1))
        for start in range(0, int(cuts[-1]) - w + 1):
            acc.deposit(0, start, start + w)
        t = acc.tally
        mass += rho * np.asarray(t.edge_unspliced_mass, np.float64)
        count += rho * np.asarray(t.edge_unspliced_count, np.int64).sum(axis=1).astype(np.float64)
    out = np.ones(n_edges, dtype=np.float64)
    np.divide(mass, count, out=out, where=count > 0)
    return out


def _uniform_library(node_len, rho_g, rho_r, pmf_g, pmf_r) -> CalibrationResult:
    """A CalibrationResult for ONE reference tiled by ``node_len``, under a UNIFORM field.

    ⭐⭐ **EVERY BANK IS ENUMERATED THROUGH THE SPECIFICATION, not evaluated analytically.** The masses
    and the conserved share must describe ONE population or the fixture is internally inconsistent —
    and that inconsistency is not hypothetical. The predecessor of this docstring used the analytic
    infinite-chromosome forms ``rho·E[(len−w+1)+]`` and ``rho·E[w−1]`` for the masses while the share
    came from a finite-reference enumeration. On the fine tiling the two disagreed by **10 %** and read
    as a defect in ``assemble_priors`` (2026-08-08).

    ⚠ For a point-mass pmf, enumerating every start at unit weight IS the analytic uniform field, so
    nothing is lost — the fixture is still exact, deterministic arithmetic with no distributional slack.

    ⛔ **And the target moved with the rule.** The old assembler produced ``rho·span``; the new one
    produces the CONSERVED FRAGMENT COUNT, which on a finite reference of span ``S`` is ``rho·(S−w+1)``
    — the number of fragments that fit. ``rho·span`` counted ``w−1`` start positions no fragment can
    occupy, and was the approximation the density conversion happened to produce.
    """
    node_len = np.asarray(node_len, dtype=np.float64)
    n = node_len.shape[0]
    ne = max(n - 1, 0)
    # the OPPORTUNITY arrays stay analytic — they are the divisors `gdna_eff_len` contracts against,
    # and they are not a population statement
    a_g_node = contained_eff_length(node_len, pmf_g)
    a_r_node = contained_eff_length(node_len, pmf_r)
    a_g_edge = np.full(ne, float(crossing_eff_length(pmf_g, _UNB, _UNB)[0]))
    a_r_edge = np.full(ne, float(crossing_eff_length(pmf_r, _UNB, _UNB)[0]))
    # ...and every MASS is what the specification actually deposits on this tiling
    _n_g, cont_g, cross_g, _m_g = _enumerate(node_len, int(np.argmax(pmf_g)))
    _n_r, cont_r, cross_r, _m_r = _enumerate(node_len, int(np.argmax(pmf_r)))
    return CalibrationResult(
        mass_gdna_node=rho_g * cont_g,
        mass_rna_node=rho_r * cont_r,
        mass_gdna_edge=rho_g * cross_g,
        mass_rna_edge=rho_r * cross_r,
        mass_rna_spliced_edge=np.zeros(ne, dtype=np.float64),
        edge_mass_per_crossing=_mass_per_crossing(node_len, rho_g, rho_r, pmf_g, pmf_r),
        mass_rna_junction=np.zeros(0, dtype=np.float64),
        edge_spliced_mass_per_crossing=np.ones_like(
            _mass_per_crossing(node_len, rho_g, rho_r, pmf_g, pmf_r)
        ),
        junction_mass_per_crossing=np.ones(0, dtype=np.float64),
        gdna_node_eff_len=a_g_node,
        gdna_edge_eff_len=a_g_edge,
        rna_node_eff_len=a_r_node,
        rna_edge_eff_len=a_r_edge,
        gdna_density_global=rho_g,
        rna_sense_frac=0.9,
        gdna_strand_overdispersion=0.05,
        rna_strand_overdispersion=0.05,
        n_nodes=n,
        n_edges=ne,
        n_junctions=0,
        config=CalibrationConfig(),
    )


def _regions_tiling(node_len) -> RegionArrays:
    """One reference tiled contiguously from 0 by ``node_len``."""
    node_len = np.asarray(node_len, dtype=np.int64)
    ends = np.cumsum(node_len)
    starts = ends - node_len
    n = node_len.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        # exon on both strands: the node must survive the locus projection's intergenic drop
        signature=np.full(n, 0b0000_0011, dtype=np.uint8),
        strand_class=np.zeros(n, dtype=np.int8),
        region_size_bp=node_len.astype(np.float64),
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

    ⛔ **Neither number is read back off ``assemble_priors``** (`TRAPS: a-test-that-redefines`). Both
    are re-derived from the reference accumulator's own banks:

    * **truth** — the conserved fragment count each component really deposited, ``SUM contained_c +
      SUM mass_c``, where ``mass_c`` is that component's OWN conserved-mass bank. On a finite reference
      of span ``S`` this is ``rho_c·(S − w_c + 1)``, the fragments that FIT, and the closed form is
      asserted against it separately.
    * **pred** — what a single POOLED ``q`` must produce: the same contained term, plus each
      component's own crossing COUNT rescaled by the MIXTURE's share. This is the pooling result
      stated as arithmetic, and its only content beyond ``truth`` is that ``q_c`` has been replaced by
      ``q_pooled``.

    ``q_c = mass_c / count_c`` per line is each component's own share, and ``pred == truth`` exactly
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
# grid moves. ⭐ The 100 bp tiling is finer than the 200 bp RNA fragment, which is the regime where
# 56.7 % of human nodes live and where the raw incidence sum diverges hardest.
_SPAN = 1200
_MU_G, _MU_R = 50, 200
_TILINGS = {
    "coarse (1 x 1200)": [1200],
    "medium (3 x 400)": [400] * 3,
    "fine   (12 x 100)": [100] * 12,
    "ragged (mixed)": [37, 400, 63, 300, 1, 199, 200],
}
# ⭐ The tilings where every node exceeds BOTH fragment lengths, so ``min(w−1, flank) == w−1`` on every
# line and ``q_g == q_r == 1``: the pooled share is then each component's own and the split is exact.
_SHARES_AGREE = ["coarse (1 x 1200)", "medium (3 x 400)"]
# ⚠ Not 1e-9: the conserved-mass bank is fixed-point at 2^-32 per fragment, so a 1,001-fragment total
# carries ~2e-11 of relative rounding. Anything above 1e-10 here would be a real error.
_FIXED_POINT_RTOL = 1e-8


@pytest.mark.parametrize("name", list(_TILINGS))
def test_the_total_prior_is_the_true_fragment_count_on_every_tiling(name):
    """⭐⭐ THE CONSERVATION GATE. The same physical library, re-tiled, deposits the same TOTAL — and it
    is the right total: the number of fragments that FIT, ``rho_g·(S−mu_g+1) + rho_r·(S−mu_r+1)``.

    The raw incidence sum grows as the tiling is refined (every new line adds a crossing to every
    fragment that spans it), so it fails here by construction. The retired ``rho_c·span_bp`` form fails
    too, by the ``w−1`` start positions no fragment can occupy — 4.1 % on this fixture.

    ⛔⛔ **AND THIS GATE IS BLIND TO THE DEFECT THE NEXT TWO TESTS EXIST FOR**
    (`TRAPS: conservation-misses-mis-attribution`). Rescaling both components by one pooled share
    conserves the total EXACTLY while tilting the split: on the fine tiling this passes at 1e-11 while
    the gDNA side alone is 19.9 % low. Never read this test as "the prior is right".
    """
    rho_g, rho_r = 0.03, 0.05
    p = _priors_for(_TILINGS[name], rho_g, rho_r, _point_pmf(_MU_G), _point_pmf(_MU_R))
    total = float(p.gdna_prior_count[0] + p.rna_prior_count[0])
    fit = rho_g * (_SPAN - _MU_G + 1) + rho_r * (_SPAN - _MU_R + 1)
    assert total == pytest.approx(fit, rel=_FIXED_POINT_RTOL), (
        f"{name}: total {total:.9f} against {fit:.9f} fragments that fit"
    )
    # ⛔ and it is NOT the retired rho·span, which is the same on every tiling and 4.1 % too big
    assert not np.isclose(total, (rho_g + rho_r) * _SPAN, rtol=1e-3)


@pytest.mark.parametrize("name", _SHARES_AGREE)
def test_each_component_is_its_true_fragment_count_where_the_two_shares_agree(name):
    """⭐ And where the pooled share IS each component's own (``q_g == q_r``), the SPLIT is exact too.

    ⚠ Stronger than the total alone — a form uniformly wrong by a constant factor passes the
    conservation gate and fails this one. ⛔ Restricted to two of the four tilings on purpose: this is
    the substrate condition an equal-length panel satisfies by construction, which is exactly why such
    a panel cannot see the bias (`TRAPS: an-equal-length-panel-defeats-the-lift`).
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
    """⛔⛔ WHERE THE TWO SHARES DISAGREE THE SPLIT IS WRONG, AND THIS PINS THE WRONG VALUE EXACTLY.

    The accumulator cannot tell the two populations apart, so ``edge_mass_per_crossing`` is the
    MIXTURE's share and both components are rescaled by it. The measured consequence, with the library
    physically unchanged and only the bookkeeping grid moved:

    ==================  ==========  ==========  ==========
    tiling              gDNA        RNA         total
    ==================  ==========  ==========  ==========
    coarse / medium     exact       exact       exact
    fine (12 x 100)     **−19.9 %** **+13.7 %** exact
    ragged (mixed)      **−2.6 %**  **+1.8 %**  exact
    ==================  ==========  ==========  ==========

    ⭐ gDNA is the SHORTER component here (50 bp against 200), so its own ``q_g`` is the larger and the
    pooled share drags it DOWN. Reverse the lengths and the sign reverses — that is the sweep below.

    ⚠ **The biased value is the SPECIFIED value until a per-component ``q`` is built**, so it is
    asserted to 1e-9 rather than tolerated with a loose bound. The recorded percentages are asserted
    too, so the bias cannot silently drift or silently vanish.
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
    """⭐ THE COMPOSITION TEST. Fixed true g:r; sweep the two components' mean lengths against each
    other. ⛔ **The prior's ratio MOVES, by 0.56× to 1.59×**, and this pins where it lands.

    ⛔ Swept in BOTH directions (``mu_g`` from 0.25x to 2x the RNA mean) — owner ruling: there is no rule
    that RNA is longer than gDNA, and assuming one is how a tool overfits to cfRNA. The direction is the
    finding: the SHORTER component is under-called and the longer one over-called, because a longer
    fragment is censored harder by a 100 bp flank and so carries the smaller share.

    ⚠ **The distortion is NOT ``q_r/q_g`` at the locus level, and the theorem's mixture-independence
    does not survive contact with contained mass.** Only the CROSSING term passes through the share; the
    contained term is already a fragment count and is untouched, so the locus-level tilt is diluted by
    each component's contained fraction — and since the pooled share itself depends on the mixture, the
    dilution does too. Measured at ``mu_g = 50``, where gDNA is 53 % contained: 0.837× at ``rho`` 0.02 /
    0.06 but 0.665× at 0.05 / 0.01, against a pure-crossing ``q_r/q_g`` of 0.5025. At ``mu_g = 100``,
    where gDNA is 1 % contained, both read ≈0.56 and the pure-crossing limit is nearly recovered.
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
    """⛔⛔ **AND AT EQUAL LENGTHS THE BIAS IS EXACTLY ZERO — which is why a panel built that way cannot
    measure it** (`TRAPS: an-equal-length-panel-defeats-the-lift`). The ladder's realised gDNA/RNA gap is
    +1.5–2.1 %; the flgap PAIR exists because of this line.

    ⭐ Asserted at two very different mixtures, because "exact" here must not depend on the mixing ratio:
    when ``q_g == q_r`` the pooled share equals both regardless of ``phi``.
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
    """⛔: an object with no opportunity for a component must emit NOTHING
    at zero precision — never a floored division. Every node here is shorter than one RNA fragment and
    the RNA crossing opportunity is zeroed, so the RNA support is identically 0.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4  # every node < 400 bp ⇒ contained_eff_length(RNA) == 0
    cal = _uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r)
    cal = _zero_rna_opportunity(cal)
    p = assemble_priors(cal, _regions_tiling(tiling), _one_locus(int(np.sum(tiling))))
    assert np.all(np.isfinite(p.rna_prior_count))
    np.testing.assert_allclose(p.rna_prior_count, [0.0])
    # ⭐ and the gDNA side is EXACT here, not merely finite: with rho_r = 0 the library is one
    # component, so the pooled share IS the gDNA's own and there is nothing to bias the split.
    np.testing.assert_allclose(p.gdna_prior_count, [0.03 * (np.sum(tiling) - 20 + 1)], rtol=1e-9)


def test_mass_on_a_zero_opportunity_object_STILL_COUNTS_because_a_count_has_no_divisor():
    """⛔⛔ **A DELIBERATE REVERSAL, AND HALF OF A PAIR.** This test used to be
    ``..._is_dropped_from_BOTH_sides`` and asserted **0.0**. It asserts the mass now, because the rule
    the drop existed for is gone from this path.

    ``mass > 0`` with ``support == 0`` is an ordinary configuration, not a corner: ``contained_eff_length``
    is exactly 0 wherever a node is shorter than that component's shortest fragment, which on the chr22
    pilot against its own measured pure pools is **21.7 % of nodes for RNA** and 18.7 % for gDNA. The
    solver can still put mass there — ``f_g`` is an inference, not a fact.

    ⭐ **What changed.** The drop existed because ``rho = SUM m / SUM S`` is a rate, and mass in the
    numerator with no exposure in the denominator inflates it — with ``mass / max(support, 1e-9)`` the
    inflation reaching ~1e9. **The prior no longer divides by anything.** ``mass_c_node[r]`` is
    ``f_c(r)·contained_count[r]`` and a contained fragment deposits on exactly one node, so the mass IS
    the count; dropping it would silently lose fragments the accumulator really deposited. The
    catastrophe is now structurally unreachable here rather than guarded, and the assertion says so: the
    prior is exactly the deposited 4 × 2.5, which is neither 0 nor any multiple of 1e9.

    ⛔ **The guard is still LIVE where a divisor still lives — the eff-length** — and that half is
    `test_priors.test_stray_mass_on_a_zero_opportunity_line_is_dropped_from_the_eff_len`, which
    perturbs it and measures the +19.97 / +44.93 bp it holds back. ⚠ Do not delete one without the
    other: alone, either one reads as a rule about the whole file.
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4
    cal = _zero_rna_opportunity(_uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r))
    # ⭐ the change from the test above: put REAL mass on the zero-opportunity RNA objects
    stray = dataclasses.replace(cal, mass_rna_node=np.full(4, 2.5))
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
        rna_edge_eff_len=np.zeros_like(cal.rna_edge_eff_len),
        mass_rna_edge=np.zeros_like(cal.mass_rna_edge),
    )
