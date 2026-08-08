"""assemble_priors — the EM pseudocounts must be FRAGMENT COUNTS, not object-incidence sums.

     (P1)

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

⭐ **THE FIX, AND WHY IT IS THE ONE-LINE STATEMENT OF THE MODEL.** Density is intensive; it is pooled as
a ratio of sums and then integrated over the span::

    rho_c = SUM_locus share·mass_c / SUM_locus share·A_c        <- ratio of SUMS, never mean of ratios
    prior_c = rho_c · span_bp                                    <- the SAME genomic span for both

``A_c`` is that component's own opportunity — ``contained_eff_length`` at a node, ``crossing_eff_length``
at a line — which is exactly what the accumulator's deposition law divides out. Under a uniform field
``mass_c = rho_c·A_c`` on every object, so ``prior_c`` is the true fragment count on ANY partition.

⚠ **The pooling is A-weighted, and that is an approximation the tests must not hide.** ``SUM m / SUM A``
is the support-weighted mean density, so ``rho·span`` is exact only where ``rho_c`` is uniform *within*
the locus. It is the same pooling `derive.gdna_density_global`'s
``rho_bg = Sum g / Sum E`` already use, and it is a strict improvement on the raw sum — but a locus with
a strong internal density gradient carries a second-order residual. Stated, not tested away.

⛔ **These tests are deterministic, not simulated.** Every mass below is the accumulator's deposition law
evaluated exactly, so a failure is a defect and never noise. That is a deliberate departure from the
plan's original end-to-end phrasing of T1/T2: rebuilding an index with extra cuts also moves the
transcript set, the reach and every effective length, so it would not isolate the partition. The
end-to-end conservation check against ``node_start_count`` is T3, and it lives in
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
    from native._accumulator_reference import INV_LENGTH_SCALE, Accumulator, Partition

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
        np.asarray(t.edge_unspliced_mass, np.uint64).astype(np.float64) / INV_LENGTH_SCALE,
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
        INV_LENGTH_SCALE,
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
        mass += rho * np.asarray(t.edge_unspliced_mass, np.float64) / INV_LENGTH_SCALE
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


# --- T1: partition invariance ---------------------------------------------------------------------

# 1200 bp of reference, tiled three ways. The library is IDENTICAL in all three; only the bookkeeping
# grid moves. ⭐ The 100 bp tiling is finer than the 200 bp RNA fragment, which is the regime where
# 56.7 % of human nodes live and where the raw sum diverges hardest.
_SPAN = 1200
_TILINGS = {
    "coarse (1 x 1200)": [1200],
    "medium (3 x 400)": [400] * 3,
    "fine   (12 x 100)": [100] * 12,
    "ragged (mixed)": [37, 400, 63, 300, 1, 199, 200],
}


@pytest.mark.parametrize("name", list(_TILINGS))
def test_prior_is_partition_invariant(name):
    """⭐ THE CORE INVARIANT. The same physical library, re-tiled, must give the same prior.

    The raw-sum form gives ``rho_c·(SUM A_c)``, which grows as the tiling is refined (every new line
    adds ``E_c[w−1]`` of opportunity), so it fails here by construction.
    """
    pmf_g, pmf_r = _point_pmf(50), _point_pmf(200)
    ref = _priors_for(_TILINGS["coarse (1 x 1200)"], 0.03, 0.05, pmf_g, pmf_r)
    got = _priors_for(_TILINGS[name], 0.03, 0.05, pmf_g, pmf_r)
    np.testing.assert_allclose(got.gdna_prior_count, ref.gdna_prior_count, rtol=1e-9)
    np.testing.assert_allclose(got.rna_prior_count, ref.rna_prior_count, rtol=1e-9)


@pytest.mark.parametrize("name", list(_TILINGS))
def test_prior_is_the_true_fragment_count(name):
    """And the invariant value is the RIGHT one: ``rho_c · span``, the fragments that started here.

    ⚠ Stronger than invariance alone — a form that was uniformly wrong by a constant factor would pass
    the invariance test and fail this one.
    """
    rho_g, rho_r = 0.03, 0.05
    p = _priors_for(_TILINGS[name], rho_g, rho_r, _point_pmf(50), _point_pmf(200))
    np.testing.assert_allclose(p.gdna_prior_count, [rho_g * _SPAN], rtol=1e-9)
    np.testing.assert_allclose(p.rna_prior_count, [rho_r * _SPAN], rtol=1e-9)


# --- T2: the length sweep -------------------------------------------------------------------------


@pytest.mark.parametrize("mu_g", [50, 100, 150, 200, 300, 400])
def test_prior_ratio_is_flat_in_the_length_ratio(mu_g):
    """⭐ THE COMPOSITION TEST. Fixed true g:r; sweep the two components' mean lengths against each
    other. The prior's ratio must not move.

    ⛔ Swept in BOTH directions (``mu_g`` from 0.25x to 2x the RNA mean) — owner ruling: there is no rule
    that RNA is longer than gDNA, and assuming one is how a tool overfits to cfRNA. The raw-sum form
    tilts by ``SUM A_g / SUM A_r``, which is monotone in ``mu_g − mu_r``, so it drifts across this sweep.
    """
    rho_g, rho_r = 0.02, 0.06
    p = _priors_for(_TILINGS["fine   (12 x 100)"], rho_g, rho_r, _point_pmf(mu_g), _point_pmf(200))
    ratio = float(p.gdna_prior_count[0] / p.rna_prior_count[0])
    assert ratio == pytest.approx(rho_g / rho_r, rel=1e-9), (
        f"prior g:r moved to {ratio:.6f} at mu_g={mu_g} against a true {rho_g / rho_r:.6f}"
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
    np.testing.assert_allclose(p.gdna_prior_count, [0.03 * np.sum(tiling)], rtol=1e-9)


def test_mass_on_a_zero_opportunity_object_is_dropped_from_BOTH_sides():
    """⭐ **THE PERTURBATION TEST (P1e).** The test above passes with a floored divisor, because its
    only zero-support object also has zero mass. This one does not.

    ``mass > 0`` with ``support == 0`` is an ordinary configuration, not a corner: ``contained_eff_length``
    is exactly 0 wherever a node is shorter than that component's shortest fragment, which on the chr22
    pilot against its own measured pure pools is **21.7 % of nodes for RNA** and 18.7 % for gDNA. The
    solver can still put mass there — ``f_g`` is an inference, not a fact.

    Two wrong answers this pins out:

    * ``mass / max(support, 1e-9)``  ⇒  a density of ~1e9. Trap 23, and how a "no data" default of
      100 % gDNA once seeded false gDNA into neighbouring exons.
    * mass kept in the numerator, support omitted from the denominator ⇒ ``rho`` inflated with no
      exposure to pay for it. **Both sides of a pooled rate, or neither.**
    """
    pmf_g, pmf_r = _point_pmf(20), _point_pmf(400)
    tiling = [50] * 4
    cal = _zero_rna_opportunity(_uniform_library(tiling, 0.03, 0.0, pmf_g, pmf_r))
    # ⭐ the change from the test above: put REAL mass on the zero-opportunity RNA objects
    stray = dataclasses.replace(cal, mass_rna_node=np.full(4, 2.5))
    regions, loci = _regions_tiling(tiling), _one_locus(int(np.sum(tiling)))
    p = assemble_priors(stray, regions, loci)
    assert np.all(np.isfinite(p.rna_prior_count)), "a floored divisor produced a non-finite prior"
    np.testing.assert_allclose(p.rna_prior_count, [0.0])
    # and the gDNA side, which DOES have opportunity everywhere, is untouched by the stray RNA mass
    np.testing.assert_allclose(p.gdna_prior_count, [0.03 * np.sum(tiling)], rtol=1e-9)


def _zero_rna_opportunity(cal: CalibrationResult) -> CalibrationResult:
    """Remove every RNA crossing opportunity, so the RNA support is identically 0 on all objects."""
    return dataclasses.replace(
        cal,
        rna_edge_eff_len=np.zeros_like(cal.rna_edge_eff_len),
        mass_rna_edge=np.zeros_like(cal.mass_rna_edge),
    )
