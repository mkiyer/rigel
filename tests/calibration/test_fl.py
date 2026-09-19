"""gDNA / RNA fragment-length laws: the pure pools, the smooth-EB build, and the two estimands.

Purity is the whole point of the pools: a length model is fitted only from populations known to be ONE
component, so nothing is ever estimated from the fragments it will later explain. The partition gates
hold the pool accessors exhaustive and disjoint, the build gates hold that a sparse pool shrinks
smoothly toward the global anchor with no threshold anywhere, and the accessor gates hold the raw
empirical views the QC report reads. Above them sit two estimands of the gDNA law: ``gdna_pmf``, the
UNIFORM-FRAME law that the chemistry makes and the opportunity and prior mathematics assume, and
``gdna_realized_pmf``, the LIBRARY-CENSUS law a sequenced gDNA fragment follows with capture selection
included, which is what the EM's per-fragment scorer conditions on. They coincide off capture and
split under it, and handing either consumer the other's estimand misassigns transcripts in bulk — so
the ROUTING is gated as hard as the estimator: geometry stays bit-identical whether or not the
realized law is computed, the realized law falls back to the uniform one exactly when it cannot be
estimated, and the on-target correction vanishes identically with no enrichment excess.

PERTURBATION (`TRAPS: perturb-every-gate`): leaking the realized law into ``gdna_pmf``, breaking the
fallback, and dropping the excess's ``-1`` each fire two gates. Removing the estimator's early
no-density guard fires NOTHING — at ``rho_off = 0`` every mass term multiplies to zero and the
mass-zero decline catches it downstream — so that guard is kept for its diagnostic string and its
near-unreachable unbracketed-fit branch, and this note stands in place of a gate that cannot fire.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.fl import (
    _fl_models_from_histograms,
    gdna_contained_fl_mass,
    gdna_fl_mass,
    rna_fl_mass,
    splash_fl_mass,
)
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    POOL_DNA_INTERGENIC,
    POOL_DNA_INTERGENIC_EXON,
    POOL_DNA_INTRONIC,
    POOL_DNA_INTRON_EXON,
    POOL_RNA_SPLICED,
)


def _spike(at: int, total: float, n: int = 1001) -> np.ndarray:
    c = np.zeros(n, dtype=np.float64)
    c[at] = total
    return c


def _pools(n_bins: int = 5) -> np.ndarray:
    """A ``pool_lengths`` block with a different value in every pool, so a wrong index cannot pass."""
    pools = np.zeros((N_FRAGMENT_POOLS, n_bins), dtype=np.int64)
    pools[POOL_DNA_INTERGENIC, 2] = 1
    pools[POOL_DNA_INTRONIC, 3] = 3
    pools[POOL_DNA_INTRON_EXON, 2] = 700
    pools[POOL_DNA_INTERGENIC_EXON, 3] = 900
    pools[POOL_RNA_SPLICED, 4] = 11
    return pools


def test_the_pool_indices_ARE_the_specifications_FragmentPool():
    """The three-way contract: the reference enum, the C++ enum and these constants are one axis.

    A silent disagreement here re-labels every pool — the gDNA model would be fitted from the RNA pool
    and nothing would look wrong. Checked against the executable specification itself, not a written-out
    list, for the same reason the payload schema test does.
    """
    from native._accumulator_reference import FragmentPool

    assert N_FRAGMENT_POOLS == len(FragmentPool)
    assert POOL_DNA_INTERGENIC == FragmentPool.DNA_INTERGENIC
    assert POOL_DNA_INTRONIC == FragmentPool.DNA_INTRONIC
    assert POOL_DNA_INTRON_EXON == FragmentPool.DNA_INTRON_EXON
    assert POOL_DNA_INTERGENIC_EXON == FragmentPool.DNA_INTERGENIC_EXON
    assert POOL_RNA_SPLICED == FragmentPool.RNA_SPLICED


def test_gdna_fl_mass_is_ALL_FOUR_gdna_pools():
    """All four, because the contained pair alone measures the SHORT HALF of one population.

    Under hybrid capture the surviving off-target gDNA sits beside a probe, and a fragment beside a
    probe *reaches* the exon boundary — so it stops being contained and becomes crossing. Fitting from
    the contained pair alone therefore reads short under capture, while all four, each divided by its
    own opportunity, stay close to truth in both regimes.

    This histogram is only meaningful paired with the matching divisor: the four pools tilt in
    opposite directions, so the raw sum is biased long (TRAPS: opposite-tilts-must-not-pool). Nothing
    in the tool consumes this sum on its own.
    """
    g = gdna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(g, [0.0, 0.0, 701.0, 903.0, 0.0])


def test_gdna_contained_fl_mass_is_the_pair_and_is_the_NO_DIVISOR_FALLBACK():
    """The honest fallback, not the convenient one.

    With no annotation offered there is no opportunity function, and pooling the four raw would be
    measurably WORSE than either the contained pair or the de-tilted four. So "no divisor" falls back to
    the pair — a model that is right off capture and knowably short under it — rather than to a sum that
    is wrong everywhere.
    """
    c = gdna_contained_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(c, [0.0, 0.0, 1.0, 3.0, 0.0])


def test_gdna_fl_mass_excludes_the_RNA_pool():
    """The circularity guard: the gDNA length model must never see a certified-RNA fragment."""
    pools = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
    pools[POOL_RNA_SPLICED, 1] = 12345
    assert float(gdna_fl_mass(SimpleNamespace(pool_lengths=pools)).sum()) == 0.0


def test_rna_fl_mass_is_the_ANNOTATED_SJ_pool_alone():
    """gDNA cannot be spliced, so an observed annotated sj certifies RNA — and only that pool
    does. ``sj_implicit`` fragments are already excluded by the accumulator, because a splice that was
    never sequenced is a product of the very model this pool is used to fit."""
    r = rna_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(r, [0.0, 0.0, 0.0, 0.0, 11.0])


def test_the_splash_pools_are_reachable_SEPARATELY_for_QC():
    """Named pools, not folded in — so 'is the off-target model mis-centred for the fragments that
    actually leak?' is answerable rather than assumed."""
    s = splash_fl_mass(SimpleNamespace(pool_lengths=_pools()))
    np.testing.assert_allclose(s, [0.0, 0.0, 700.0, 900.0, 0.0])


def test_the_two_COMPONENT_accessors_partition_every_pool():
    """Teeth on the partition itself: no pool double-counted, none unreachable.

    ``gdna_fl_mass`` and ``rna_fl_mass`` are the two COMPONENT accessors and they must be exhaustive and
    disjoint — a pool no accessor returns is silently discarded evidence, and one two accessors return is
    double-counted. Both are invisible in any single-accessor test.

    ``gdna_contained_fl_mass`` and ``splash_fl_mass`` are named SUBSETS of the gDNA side, not
    members of the partition; the test below pins that they tile it exactly.
    """
    payload = SimpleNamespace(pool_lengths=_pools())
    total = float(gdna_fl_mass(payload).sum()) + float(rna_fl_mass(payload).sum())
    assert total == float(_pools().sum()), "the two component accessors must be exhaustive"

    for pool in range(N_FRAGMENT_POOLS):
        one = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
        one[pool, 1] = 5
        p = SimpleNamespace(pool_lengths=one)
        reached = [float(gdna_fl_mass(p).sum()), float(rna_fl_mass(p).sum())]
        assert sorted(reached) == [0.0, 5.0], (
            f"pool {pool} reaches {reached}, expected exactly one component"
        )


def test_the_gdna_subsets_TILE_the_gdna_side_exactly():
    """contained + splash == all four, on every pool individually.

    This is what keeps "the crossing pools are reported separately" and "the crossing pools are fitted"
    from drifting apart: they are the same rows seen twice, not two different definitions.
    """
    payload = SimpleNamespace(pool_lengths=_pools())
    np.testing.assert_allclose(
        gdna_contained_fl_mass(payload) + splash_fl_mass(payload), gdna_fl_mass(payload)
    )
    for pool in range(N_FRAGMENT_POOLS):
        one = np.zeros((N_FRAGMENT_POOLS, 5), dtype=np.int64)
        one[pool, 1] = 5
        p = SimpleNamespace(pool_lengths=one)
        np.testing.assert_allclose(gdna_contained_fl_mass(p) + splash_fl_mass(p), gdna_fl_mass(p))


def test_build_fl_large_pool_is_empirical():
    glob = _spike(100, 1.0e6)
    gdna = _spike(300, 1.0e7)  # huge gDNA pool ⇒ dominated by its own evidence
    fl = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    assert int(np.argmax(fl.gdna_pmf)) == 300
    np.testing.assert_allclose(fl.gdna_pmf.sum(), 1.0)
    np.testing.assert_allclose(fl.global_pmf.sum(), 1.0)


def test_build_fl_empty_pool_collapses_to_global():
    glob = _spike(100, 1.0e6)
    fl = _fl_models_from_histograms(
        global_counts=glob,
        rna_counts=glob,
        gdna_counts=np.zeros(1001, dtype=np.float64),
        max_size=1000,
    )
    np.testing.assert_allclose(fl.gdna_pmf, fl.global_pmf)  # no gDNA evidence → global
    assert fl.n_gdna == 0.0


def test_build_fl_small_pool_shrinks_toward_global_no_cliff():
    glob = _spike(100, 1.0e6)  # global mode at 100
    gdna = _spike(300, 1.0)  # a single gDNA fragment at 300
    fl = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    # Smooth shrinkage (no threshold): with pool_total=1 ≪ ρ_ess=1000 the global
    # anchor (bin 100) dominates the lone gDNA fragment (bin 300).
    assert fl.gdna_pmf[100] > fl.gdna_pmf[300]
    assert fl.gdna_pmf[300] > 0.0  # but the empirical fragment still nudges its bin


def test_build_fl_just_below_5000_is_not_a_cliff():
    # 4999 vs 5001 differ only marginally — there is no GOOD/WEAK jump to land on.
    glob = _spike(100, 1.0e6)
    lo = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=_spike(300, 4999.0), max_size=1000
    )
    hi = _fl_models_from_histograms(
        global_counts=glob, rna_counts=glob, gdna_counts=_spike(300, 5001.0), max_size=1000
    )
    np.testing.assert_allclose(lo.gdna_pmf[300], hi.gdna_pmf[300], rtol=1e-3)


# ---------------------------------------------------------------------------
# Raw empirical counts + FragmentLengthModel accessors (QC report views)
# ---------------------------------------------------------------------------


def test_build_fl_stores_raw_aligned_counts():
    """Raw (unsmoothed) counts are stored aligned to max_size, totals = n_*."""
    glob = _spike(100, 500.0)
    rna = _spike(200, 300.0)
    gdna = _spike(300, 40.0)
    fl = _fl_models_from_histograms(
        global_counts=glob, rna_counts=rna, gdna_counts=gdna, max_size=1000, prior_ess=1000.0
    )
    assert fl.global_counts.shape == (1001,)
    assert int(np.argmax(fl.rna_counts)) == 200
    assert int(np.argmax(fl.gdna_counts)) == 300
    # Raw counts are the empirical evidence — NOT the EB-smoothed pmf.
    assert fl.rna_counts[200] == 300.0
    assert fl.gdna_counts[300] == 40.0
    assert fl.n_global == 500.0
    assert fl.n_rna == 300.0
    assert fl.n_gdna == 40.0


def test_build_fl_counts_overflow_folds_into_last_bin():
    """Counts beyond max_size fold into the overflow bin when aligned."""
    over = np.zeros(1201, dtype=np.float64)
    over[1150] = 7.0  # beyond max_size=1000
    fl = _fl_models_from_histograms(
        global_counts=_spike(100, 10.0),
        rna_counts=_spike(100, 10.0),
        gdna_counts=over,
        max_size=1000,
    )
    assert fl.gdna_counts.shape == (1001,)
    assert fl.gdna_counts[1000] == 7.0
    assert fl.n_gdna == 7.0


def test_accessors_return_empirical_models():
    """rna_model()/gdna_model()/global_model() expose the raw empirical FL."""
    fl = _fl_models_from_histograms(
        global_counts=_spike(100, 500.0),
        rna_counts=_spike(200, 300.0),
        gdna_counts=_spike(325, 40.0),
        max_size=1000,
        prior_ess=1000.0,
    )
    rna_m = fl.rna_model()
    gdna_m = fl.gdna_model()
    glob_m = fl.global_model()
    # Empirical modes/means track the raw counts (not smoothed toward global).
    assert rna_m.mode == 200
    assert gdna_m.mode == 325
    assert glob_m.mode == 100
    assert rna_m.mean == pytest.approx(200.0)
    # n_observations = pool total.
    assert rna_m.n_observations == 300
    assert gdna_m.n_observations == 40
    # to_dict() works (drives the summary QC report).
    assert gdna_m.to_dict()["summary"]["mode"] == 325


# ── the two estimands of the gDNA length law, and the routing guarantee between them ─────────


def _models(**kw):
    rng = np.random.default_rng(5)
    g = rng.random(101)
    r = rng.random(101)
    return _fl_models_from_histograms(
        global_counts=g + r, rna_counts=r, gdna_counts=g, max_size=100, **kw
    )


# ── the field exists and is never None: the scorer must be able to read it unconditionally ───────


def test_the_adjacent_pair_table_IS_the_reference_major_walk():
    """The census used to walk every adjacent region pair of every reference in Python. The table it
    walked is a property of the partition, so it is built once as arrays — and this asserts that the
    arrays are the SAME pairs in the SAME order, because the order is what makes the sums over them
    reproduce the loop's arithmetic bit for bit.

    Three references on purpose: one ordinary, one with a single region, which contributes no pair, and
    one whose boundary count disagrees with its region count, which the walk skipped rather than
    guessed at. PERTURBATION: dropping the usable mask emits pairs for the malformed reference, and
    ordering by anything but reference-major changes which sums land where.
    """
    from rigel.calibration.fl import _adjacent_pairs

    # ref 0: regions 0..3, boundaries 0..2 (3 = 4 - 1) — usable
    # ref 1: one region, no boundary — no pair
    # ref 2: regions 5..8, but only 2 boundaries where 3 are needed — skipped
    region_offsets = np.array([0, 4, 5, 9], dtype=np.int64)
    boundary_offsets = np.array([0, 3, 3, 5], dtype=np.int64)
    left, boundary = _adjacent_pairs(region_offsets, boundary_offsets)
    assert left.tolist() == [0, 1, 2]
    assert boundary.tolist() == [0, 1, 2]

    # and the plain case: two usable references, concatenated in reference order
    left, boundary = _adjacent_pairs(
        np.array([0, 3, 6], dtype=np.int64), np.array([0, 2, 4], dtype=np.int64)
    )
    assert left.tolist() == [0, 1, 3, 4]
    assert boundary.tolist() == [0, 1, 2, 3]
    # an empty partition is legal and gives empty arrays, not an error
    for empty in (np.array([0], np.int64), np.array([0, 1], np.int64)):
        got_left, got_boundary = _adjacent_pairs(empty, np.zeros(empty.size, np.int64))
        assert got_left.size == 0 and got_boundary.size == 0


def test_the_realized_law_is_always_present():
    m = _models()
    assert m.gdna_realized_pmf is not None
    assert m.gdna_realized_pmf.shape == m.gdna_pmf.shape


def test_without_a_realized_estimate_the_two_estimands_coincide_exactly():
    """The fallback is one law, byte-equal, so a consumer that reads the realized field on an
    off-capture or input-starved build gets the uniform law rather than an absent one."""
    m = _models()
    np.testing.assert_array_equal(m.gdna_realized_pmf, m.gdna_pmf)


def test_a_supplied_realized_histogram_is_shrunk_like_its_sibling():
    rng = np.random.default_rng(9)
    realized = rng.random(101) * 50
    m = _models(gdna_realized_counts=realized)
    assert not np.array_equal(m.gdna_realized_pmf, m.gdna_pmf)
    assert m.gdna_realized_pmf.min() > 0.0  # EB-shrunk toward the global anchor, like gdna_pmf
    assert m.gdna_realized_pmf.sum() == pytest.approx(1.0)


def test_the_uniform_law_is_bit_identical_with_and_without_the_realized_input():
    """The routing guarantee: computing the realized law must not perturb `gdna_pmf` by one ULP.
    Every geometry consumer reads that field, and geometry eating the wrong estimand is the large
    end-to-end regression this whole split exists to prevent."""
    rng = np.random.default_rng(9)
    a = _models()
    b = _models(gdna_realized_counts=rng.random(101) * 50)
    np.testing.assert_array_equal(a.gdna_pmf, b.gdna_pmf)
    np.testing.assert_array_equal(a.rna_pmf, b.rna_pmf)
    np.testing.assert_array_equal(a.global_pmf, b.global_pmf)


# ── the realized estimator's own invariants, on constructed payloads ─────────────────────────────


def _payload_fixture(boundary_excess: float):
    """A two-reference toy: intergenic/intron/exon regions with uniform-consistent contained counts,
    and intron|exon + intergenic|exon boundaries whose counts carry ``boundary_excess`` TIMES the
    uniform expectation. At 1.0 the boundaries say "no capture"."""
    from _fl_realized_fixture import build_fixture

    return build_fixture(boundary_excess)


def test_the_two_boundary_CLASSES_are_told_apart_by_what_flanks_the_exon():
    """A boundary's class is the NON-exon side: an exon against an INTRON is one estimand, an exon
    against anything else the other, and the two are inverted against different crossing pools. The
    fixture's symmetric form cannot see the difference — every off-target region carries exactly the
    uniform expectation, so both classes get the same weight and swapping their labels swaps two equal
    sums.

    So this breaks the symmetry: the intron carries RNA-side excess and the intergenic flanks do not,
    which pushes the intron-flanking class's weight below 1 while the intergenic one stays at it.
    PERTURBATION: swapping the two class labels fires this and nothing else in the file.
    """
    import dataclasses

    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=50.0)
    counts = np.array(payload.region_contained_count, dtype=np.float64, copy=True)
    counts[2] *= 40.0  # the INTRON alone reads far above the uniform field: RNA, not gDNA
    payload = dataclasses.replace(payload, region_contained_count=counts)

    _counts, _uniform_out, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.intron_exon_share < diag.intergenic_exon_share, (
        "the intron-flanking class carries the RNA-diluted weight; the intergenic one does not"
    )
    # not vacuous: the intergenic class is essentially undiluted, so the gap is the intron's doing
    assert diag.intergenic_exon_share > 0.99


def test_no_enrichment_excess_means_no_on_target_correction():
    """The closure property: boundaries consistent with the uniform field ⇒ the (eps−1)+ term is
    identically zero and the realized law is the sampled blend alone."""
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=1.0)
    counts, _uniform_out, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share == pytest.approx(0.0, abs=1e-9)


def test_enriched_boundaries_raise_the_on_target_share():
    """At 50x enrichment the fixture's own arithmetic puts the excess classes near
    ``rho·49·E_contained(200)·2`` of a ~2.6k total — about 0.27. Gate the ORDER, derived, not a guess:
    well clear of the no-excess case's exact 0, and the boundary share alive beside it."""
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=50.0)
    counts, _uniform_out, diag = _realized_gdna_counts(payload, opp, rl, rt, rna_pmf, uniform)
    assert diag.applied
    assert diag.ontarget_share > 0.2
    assert diag.boundary_share > 0.05


def test_the_on_target_correction_rises_smoothly_with_enrichment():
    """The capture SPECTRUM, not a switch: the correction is EXACTLY 0 with no excess and rises
    monotonically from there. Its weight is its OWN resolution — how far an exon's enrichment sits
    from 1 against its own sampling error — and NOT the strata-split weight `lam`, which asks whether
    the two strata's LAWS differ. That is a different question, and gating on it suppresses a
    correction that is genuinely resolved."""
    from rigel.calibration.fl import _realized_gdna_counts

    shares = [
        _realized_gdna_counts(*_payload_fixture(boundary_excess=x))[2].ontarget_share
        for x in (1.0, 2.0, 5.0, 10.0, 25.0, 50.0)
    ]
    assert shares[0] == 0.0
    assert all(b > a for a, b in zip(shares, shares[1:]))


def test_an_unresolved_enrichment_weight_collapses_with_its_precision():
    """The excess's own weight, gated at the MECHANISM rather than end to end, and the reason is worth
    recording. PERTURBATION: every way of starving the enrichment's precision in the fixture —
    thinning the boundaries, or the whole library at fixed enrichment — also trips an EARLIER guard,
    so the end-to-end share reads 0 whether or not this weight exists and cannot isolate it. What can
    be isolated is the weight itself, which is where the logic lives.

    The two arms below hold the same enrichment ratio on a well-populated crossing pool and on an
    almost empty one, and the weight must go to ~1 and to ~0 respectively. The ratio alone would fire
    identically at both.
    """
    from rigel.calibration.fl import _resolution_weight

    eps = 10.0
    strong = _resolution_weight((eps - 1.0) ** 2, eps * eps / 98.0)
    weak = _resolution_weight((eps - 1.0) ** 2, eps * eps / 1e-3)
    assert strong > 0.95
    assert weak < 0.01
    # and it is monotone in the evidence, with no step anywhere along it
    ws = [_resolution_weight((eps - 1.0) ** 2, eps * eps / n) for n in np.logspace(-4, 4, 200)]
    assert all(b >= a for a, b in zip(ws, ws[1:]))
    assert max(abs(b - a) for a, b in zip(ws, ws[1:])) < 0.05


def test_zero_gdna_declines_rather_than_fabricating_a_law():
    from rigel.calibration.fl import _realized_gdna_counts

    payload, opp, rl, rt, rna_pmf, uniform = _payload_fixture(boundary_excess=1.0)
    import dataclasses

    rc = np.zeros_like(np.asarray(payload.region_contained_count))
    pl = np.asarray(payload.pool_lengths, dtype=np.float64).copy()
    pl[:4] = 0.0
    starved = dataclasses.replace(payload, region_contained_count=rc, pool_lengths=pl)
    counts, _uniform_out, diag = _realized_gdna_counts(starved, opp, rl, rt, rna_pmf, uniform)
    assert counts is None and not diag.applied


# ── the CONVERGENCE law: no cliffs, and the two estimands merge when the split is unresolvable ───


def test_resolution_weight_is_a_smooth_signal_to_noise_ratio():
    """`lam = S/(S+N)`: 0 when the split is pure noise, 1 when noise vanishes, ½ at S = N, and
    MONOTONE in between. No threshold, so there is no value of the data at which behaviour jumps."""
    from rigel.calibration.fl import _resolution_weight

    assert _resolution_weight(0.0, 1.0) == 0.0
    assert _resolution_weight(1.0, 0.0) == 1.0
    assert _resolution_weight(1.0, 1.0) == pytest.approx(0.5)
    prev = -1.0
    for s in np.linspace(0.0, 10.0, 60):
        lam = _resolution_weight(float(s), 1.0)
        assert 0.0 <= lam <= 1.0 and lam >= prev
        prev = lam


def test_the_estimands_converge_when_the_boundary_stratum_is_starved():
    """The capture-OFF side: sparse boundary data must not let the two laws diverge. The difference is
    unmeasurable there, so the honest answer is that they agree."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(3)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)  # a real shift, but measured on almost nothing
    uni, real, lam = _couple_estimands(g_c, 1e6, g_b, 1e-3)
    assert lam < 1e-3
    # The guarantee is that the residual disagreement is BOUNDED BY lam times the split, not that it
    # is bitwise zero — asserting equality would be asserting more than the mathematics gives.
    assert float(np.abs(uni - real).sum()) <= lam * float(np.abs(g_b - g_c).sum()) + 1e-12


def test_the_estimands_converge_when_the_CONTAINED_stratum_is_starved():
    """The infinite-capture side: with no off-target data the chemistry law is not estimable, so it
    must borrow the only law that IS rather than hold a stale value."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(4)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    uni, real, lam = _couple_estimands(g_c, 1e-3, g_b, 1e6)
    assert lam < 1e-3
    assert float(np.abs(uni - real).sum()) <= lam * float(np.abs(g_b - g_c).sum()) + 1e-12
    # and the common law it converged on is the well-measured one, not the starved one
    assert float(np.abs(uni - g_b).sum()) < float(np.abs(uni - g_c).sum())


def test_a_well_measured_split_keeps_the_two_estimands_apart():
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(5)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    uni, real, lam = _couple_estimands(g_c, 1e7, g_b, 1e7)
    assert lam > 0.99
    # lam -> 1 recovers the uncoupled behaviour, to within (1 - lam) times the split
    assert float(np.abs(uni - g_c).sum()) <= (1.0 - lam) * float(np.abs(g_b - g_c).sum()) + 1e-12
    assert float(np.abs(uni - real).sum()) > 0.1


def test_there_is_no_cliff_anywhere_along_the_capture_spectrum():
    """The anti-cliff gate: sweep the boundary stratum's mass over eight orders of magnitude and
    assert the uniform law moves CONTINUOUSLY — a binary decline would show as a step."""
    from rigel.calibration.fl import _couple_estimands

    rng = np.random.default_rng(6)
    g_c = np.abs(rng.random(60)) + 0.01
    g_c /= g_c.sum()
    g_b = np.roll(g_c, 6)
    masses = np.logspace(-4, 4, 200)
    unis = np.array([_couple_estimands(g_c, 1e3, g_b, float(m))[0] for m in masses])
    steps = np.abs(np.diff(unis, axis=0)).sum(axis=1)
    assert steps.max() < 0.05, f"a step of {steps.max():.3f} is a cliff, not a fade"
