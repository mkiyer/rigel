"""Tests for the gDNA strand Beta-Binomial overdispersion fit — the AWAY-HALF estimator.

Core property: data generated with a known overdispersion is recovered by the estimator from pure seeds,
and contaminated seeds can only pull it DOWN or into the fallback, never up (the contamination gates
themselves live in ``test_gdna_strand_fit.py``).
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.gdna_strand import (
    _MAX_OVERDISPERSION,
    GdnaStrandModel,
    fit_gdna_strand_from_substrate,
    fit_gdna_strand_overdispersion,
    overdispersion_for_beta,
)
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.strand_likelihood import strand_loglik


def _beta_binom_regions(rng, n_regions, depth, overdispersion, mean=0.5):
    """Draw per-region (sense, total) from BetaBinom(depth, mean, overdispersion).

    Mirrors the simulator: each region has a shared latent sense rate p ~ Beta(a, b),
    then sense ~ Binomial(depth, p). With a symmetric mean=½, a = b = ½(1−od)/od.
    """
    conc = (1.0 - overdispersion) / overdispersion  # a + b
    a = mean * conc
    b = (1.0 - mean) * conc
    p = rng.beta(a, b, size=n_regions)
    total = np.full(n_regions, depth, dtype=np.float64)
    sense = rng.binomial(depth, p).astype(np.float64)
    return sense, total


@pytest.mark.parametrize("true_od", [0.01, 0.05, 0.10, 0.20])
def test_recovers_overdispersion_pure_gdna(true_od):
    """Pure-gDNA seeds: recovered overdispersion ≈ truth (from half the pairs — the away half)."""
    rng = np.random.default_rng(12345)
    sense, total = _beta_binom_regions(rng, n_regions=4000, depth=120, overdispersion=true_od)
    model = fit_gdna_strand_overdispersion(sense, total, rna_sense_frac=0.95)
    assert not model.fallback_used
    # MoM on 4000 regions: relative error should be small; absolute floor for tiny od.
    assert model.gdna_strand_overdispersion == pytest.approx(true_od, rel=0.20, abs=0.005)


def test_binomial_limit_recovers_near_zero():
    """Zero overdispersion (Binomial) → fit returns ≈ floor, not a spurious positive."""
    rng = np.random.default_rng(7)
    # Binomial: shared rate exactly ½ (no Beta spread).
    total = np.full(4000, 120, dtype=np.float64)
    sense = rng.binomial(120, 0.5, size=4000).astype(np.float64)
    model = fit_gdna_strand_overdispersion(sense, total, rna_sense_frac=0.95)
    assert model.gdna_strand_overdispersion < 0.01


@pytest.mark.parametrize("weight", [0.8, 0.5])
def test_uniform_contamination_never_inflates_the_fit(weight):
    """EVERY seed 20 % / 50 % stranded RNA at depth 200: each is pulled far onto the RNA side, so the
    away half is reached by noise alone — the fit falls back or sits at/below the truth. It is NEVER
    inflated, which is the lemma's one-sided guarantee; the predecessor read this RNA as gDNA spread."""
    rng = np.random.default_rng(99)
    true_od = 0.10
    kappa = 0.95
    depth = 200
    n_regions = 6000
    conc = (1.0 - true_od) / true_od
    a = 0.5 * conc
    sense = np.empty(n_regions)
    total = np.full(n_regions, depth, dtype=np.float64)
    for i in range(n_regions):
        n_g = rng.binomial(depth, weight)
        p_g = rng.beta(a, a)  # shared gDNA rate for this region
        sense[i] = rng.binomial(n_g, p_g) + rng.binomial(depth - n_g, kappa)
    model = fit_gdna_strand_overdispersion(sense, total, rna_sense_frac=kappa)
    # ⚠ Either the away half saw enough to answer (and it cannot be INFLATED), or it saw nothing at all and
    # says so — the ceiling fallback is "I measured nothing", never a claim about this population.
    assert model.fallback_used or model.gdna_strand_overdispersion <= true_od + 0.02


def test_thin_seed_fallback():
    """No gDNA strand signal (empty, or every seed on the RNA side) → the CEILING, no crash.

    ⭐ 2026-08-30 (owner ruling): no conjured constant survives anywhere in this fit. Having measured
    nothing, the estimator returns the widest dispersion the model admits, so the strand channel says
    NOTHING rather than something confident; ``calibrate`` then reconciles it against the RNA fit, which
    usually HAS measured something (``gdna_strand.reconcile_overdispersions``)."""
    empty = fit_gdna_strand_overdispersion(np.array([]), np.array([]), rna_sense_frac=0.95)
    assert empty.fallback_used
    assert empty.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)
    assert empty.information == 0.0

    rng = np.random.default_rng(1)
    total = np.full(10, 100.0)
    sense = rng.binomial(100, 0.95, size=10).astype(np.float64)
    # pure RNA at κ = 0.95: every seed sits on the RNA side of ½, so the away half is EMPTY
    no_gdna = fit_gdna_strand_overdispersion(sense, total, rna_sense_frac=0.95)
    assert no_gdna.fallback_used
    assert no_gdna.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)


def test_beta_concentration_roundtrip():
    """Model exposes the Beta(a, a) concentration consistent with a = ½(1−od)/od."""
    m = GdnaStrandModel(
        gdna_strand_overdispersion=0.1,
        information=1000.0,
        n_seed_regions=10,
        n_seed_fragments=1000,
        fallback_used=False,
    )
    assert m.beta_concentration() == pytest.approx(0.5 * (1.0 - 0.1) / 0.1)


# --- the substrate wrapper (Phase 2 extraction + fit) ---------------------------------------


def _view(pos, neg):
    """A :class:`PopulationView`-shaped stand-in: ONE genome-strand count array, nothing else.

    ⭐ ``mass_unspliced`` is gone with the fractional numerator — the accumulator deposits ``+1`` on
    every object the fragment touched, so ``count`` IS the mass and there is no second array to keep
    consistent with it.
    """
    return SimpleNamespace(
        count=np.stack(
            [np.asarray(pos, dtype=np.float64), np.asarray(neg, dtype=np.float64)], axis=1
        )
    )


def _mock_substrate(pos, neg, ts, count_evidence, observable):
    """Mock with only contained-REGION signal: the boundary axis carries ZERO counts ⇒ no boundary seeds.

    ⚠ One reference with ``n`` regions owns exactly ``n − 1`` boundaries, so the boundary arrays are sized from
    the region axis rather than left empty — an boundary axis inconsistent with its own ``ref_id`` is not a
    "no boundaries" fixture, it is a mis-shaped one.
    """
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n = pos.shape[0]
    n_boundaries = max(n - 1, 0)
    substrate = SimpleNamespace(
        region_contained=_view(pos, neg),
        boundary_unspliced=_view(np.zeros(n_boundaries), np.zeros(n_boundaries)),
    )
    region_arrays = SimpleNamespace(strand_class=np.asarray(ts), ref_id=np.zeros(n, dtype=np.int64))
    ce = np.asarray(count_evidence, dtype=np.float64)
    region_density = SimpleNamespace(
        count_evidence=ce,
        region_count_observable=np.asarray(observable, dtype=bool),
        boundary_count_observable=np.zeros(n_boundaries, dtype=bool),  # no boundary seeds
        density=np.zeros(n, dtype=np.float64),
    )
    return substrate, region_arrays, region_density


def test_wrapper_recovers_overdispersion_from_intron_seeds():
    """End-to-end of the wrapper: pure-gDNA INTRON seed regions (a + gene) → recovered od."""
    rng = np.random.default_rng(2024)
    true_od = 0.10
    sense, total = _beta_binom_regions(rng, n_regions=3000, depth=100, overdispersion=true_od)
    pos, neg = sense, total - sense
    n = len(total)
    substrate, region_arrays, region_density = _mock_substrate(
        pos, neg, np.full(n, TS_POS), count_evidence=total.copy(), observable=np.ones(n, bool)
    )
    model = fit_gdna_strand_from_substrate(
        substrate,
        region_arrays,
        region_count_observable=region_density.region_count_observable,
        boundary_count_observable=region_density.boundary_count_observable,
        rna_sense_frac=0.95,
    )
    assert not model.fallback_used
    assert model.n_seed_regions == n
    assert model.gdna_strand_overdispersion == pytest.approx(true_od, rel=0.25, abs=0.01)


def test_wrapper_excludes_ambig_and_non_observable():
    """AMBIG regions (no defined sense) and non-count-observable (exonic) regions are dropped."""
    n = 100
    pos = np.full(n, 50.0)
    neg = np.full(n, 50.0)
    # all AMBIG → excluded even though count-observable
    s_ambig, ra_ambig, nd_ambig = _mock_substrate(
        pos, neg, np.full(n, TS_AMBIG), np.full(n, 100.0), np.ones(n, bool)
    )
    assert (
        fit_gdna_strand_from_substrate(
            s_ambig,
            ra_ambig,
            region_count_observable=nd_ambig.region_count_observable,
            boundary_count_observable=nd_ambig.boundary_count_observable,
            rna_sense_frac=0.95,
        ).n_seed_regions
        == 0
    )
    # intergenic but non-observable (exonic) → excluded
    s_exon, ra_exon, nd_exon = _mock_substrate(
        pos, neg, np.full(n, TS_NONE), np.full(n, 100.0), np.zeros(n, bool)
    )
    assert fit_gdna_strand_from_substrate(
        s_exon,
        ra_exon,
        region_count_observable=nd_exon.region_count_observable,
        boundary_count_observable=nd_exon.boundary_count_observable,
        rna_sense_frac=0.95,
    ).fallback_used


# --- boundary seeds: ONE per boundary, not two per boundary -------------------------------------------


def _boundary_parts(signatures, boundary_pos, boundary_neg, ref_id=None):
    """A substrate + geometry-free ``RegionDensity`` over ``signatures``, on the region/boundary axes."""
    from rigel.calibration.density_model import count_observable_masks
    from rigel.calibration.signature import transcript_strand_class

    sig = np.asarray(signatures, dtype=np.uint8)
    n = sig.shape[0]
    ref = np.zeros(n, dtype=np.int64) if ref_id is None else np.asarray(ref_id, dtype=np.int64)
    region_obs, boundary_obs = count_observable_masks(sig, ref)
    substrate = SimpleNamespace(
        region_contained=_view(np.zeros(n), np.zeros(n)),
        boundary_unspliced=_view(boundary_pos, boundary_neg),
    )
    region_arrays = SimpleNamespace(
        signature=sig, strand_class=transcript_strand_class(sig.astype(np.int64)), ref_id=ref
    )
    region_density = SimpleNamespace(
        region_count_observable=region_obs,
        boundary_count_observable=boundary_obs,
        density=np.zeros(n),
    )
    return substrate, region_arrays, region_density


def test_boundary_seeds_emits_ONE_seed_per_boundary_not_two_per_boundary():
    """⭐ The S5.f collapse, and the whole point of the change.

    The predecessor emitted TWO seeds for one boundary — region ``r``'s right side and region ``r+1``'s
    left side — because the old accumulator split one crossing fragment's mass across the two flanks.
    They were the same physical crossing, counted twice into one pooled moment estimator, which
    inflates its apparent sample size by 2× and correlates every pair perfectly. A contiguous boundary is
    a 0-bp boundary with ONE count, so there is ONE seed.
    """
    from rigel.calibration.strand_deconv import boundary_seeds

    # exon+ | intron+ : one boundary, count-observable (no shared EXON bit), oriented POS.
    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_EXON_POS, BIT_INTRON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, total = boundary_seeds(
        substrate, region_arrays, region_density.boundary_count_observable
    )
    assert sense.shape == (1,)  # ⛔ was (2,) — the same crossing, twice
    np.testing.assert_allclose(sense, [70.0])
    np.testing.assert_allclose(total, [100.0])


def test_boundary_seed_sense_follows_the_flanking_transcript_strand():
    """A NEG-strand boundary orients to the NEG genome column; the sense count is not always ``pos``."""
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_INTRON_NEG, 0], boundary_pos=[70.0], boundary_neg=[30.0]
    )  # a NEG gene's edge, intergenic on the right
    sense, total = boundary_seeds(
        substrate, region_arrays, region_density.boundary_count_observable
    )
    np.testing.assert_allclose(sense, [30.0])
    np.testing.assert_allclose(total, [100.0])


def test_an_intergenic_flank_is_a_strand_WILDCARD():
    """Intergenic carries no transcript, so a gene boundary is oriented by its gene flank."""
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [0, BIT_INTRON_NEG], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    np.testing.assert_allclose(sense, [30.0])  # oriented NEG by the gene side


def test_an_opposite_strand_boundary_is_not_strand_observable():
    """``{POS, NEG}`` leaves 'sense' undefined, so the boundary cannot seed the fit at all."""
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_INTRON_POS, BIT_INTRON_NEG], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


def test_an_AMBIG_flank_cannot_seed():
    """⛔ Found by PERTURBATION, not by review: nothing pinned this rule.

    A flank carrying overlapping ± transcripts has no defined transcript sense, so neither genome
    column is "sense" and the boundary cannot seed a strand fit. ``boundary_strand_orientation`` used to carry
    an explicit ``~either_ambig`` guard — deleting it broke NOTHING, because ``TS_AMBIG`` is a fourth
    distinct value that the ``POS-or-NONE`` / ``NEG-or-NONE`` tests already exclude. The guard was
    dead code claiming to be the rule; this test is the rule.

    ⚠ The fixture uses INTRON bits on both strands, not exon bits: an AMBIG flank must still be
    count-observable, or the test would pass for the wrong reason.
    """
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_INTRON_POS | BIT_INTRON_NEG, BIT_INTRON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    assert region_density.boundary_count_observable[
        0
    ]  # count-observable — it fails on STRAND alone
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


def test_a_boundary_inside_one_exon_is_not_count_observable():
    """A shared exon bit means an exon-strand continues across the boundary, so unspliced MATURE RNA
    crosses it and its count is not gDNA."""
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_EXON_POS, BIT_EXON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


def test_boundary_seeds_never_straddle_a_reference():
    """Two single-region references own ZERO boundaries between them — nothing can leak across."""
    from rigel.calibration.strand_deconv import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_INTRON_POS, BIT_INTRON_POS], boundary_pos=[], boundary_neg=[], ref_id=[0, 1]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


# --- the deconv application (strand_likelihood) ---------------------------------------------


def test_strand_loglik_od_zero_is_binomial_variance():
    """od → 0 recovers the Binomial mixture variance exactly (regression-safe)."""
    grid = np.linspace(0.01, 0.99, 50)
    sense, antisense, kappa = 30.0, 70.0, 0.95
    ll0 = strand_loglik(grid, sense, antisense, kappa, gdna_strand_overdispersion=0.0)
    n = sense + antisense
    p = 0.5 * grid + kappa * (1.0 - grid)
    var = n * p * (1.0 - p)
    ref = -0.5 * (sense - n * p) ** 2 / var - 0.5 * np.log(var)
    np.testing.assert_allclose(ll0, ref)


def test_overdispersion_rescues_strand_skewed_gdna():
    """The 126/26 motivation: a noise-skewed pure-gDNA split reads as MORE gDNA at od>0.

    Under Binomial (od=0) a strong sense skew is implausible for symmetric gDNA and gets pulled
    toward RNA; the Beta-Binomial (od>0) tolerates the skew and keeps more of it as gDNA.
    """
    grid = np.linspace(1e-4, 1.0 - 1e-4, 400)

    def median_gdna_frac(od):
        ll = strand_loglik(grid, 126.0, 26.0, 0.95, gdna_strand_overdispersion=od)
        w = np.exp(ll - ll.max())
        w /= w.sum()
        return float(np.interp(0.5, np.cumsum(w), grid))

    assert median_gdna_frac(0.2) > median_gdna_frac(0.05) > median_gdna_frac(0.0)


# --- prior conversion + shrinkage (replaces the old min-region / significance gates) -----------


def test_overdispersion_for_beta_conversion():
    """od = 1/(2a+1): a=1→1/3 (uniform), a=2→0.2 (ceiling), a=3→1/7 (default prior)."""
    assert overdispersion_for_beta(1.0) == pytest.approx(1.0 / 3.0)
    assert overdispersion_for_beta(2.0) == pytest.approx(0.2)
    assert overdispersion_for_beta(3.0) == pytest.approx(1.0 / 7.0)
    # Beta(2,2) is the ceiling the estimator clamps to.
    assert _MAX_OVERDISPERSION == pytest.approx(0.2)


def test_no_shrinkage_the_fit_is_the_raw_moment():
    """⭐ 2026-08-29 (owner ruling): the gDNA fit carries NO location prior. Deep seeds follow the fit,
    and sparse seeds return the raw (noisy) moment inside the physical support — they are not pulled
    toward any constant, because every constant on offer was either conjured (Beta(14,14)) or measured to
    over-shrink where gDNA information is scarce (the RNA overdispersion, on the true-od arm)."""
    rng = np.random.default_rng(11)
    s, t = _beta_binom_regions(rng, n_regions=5000, depth=120, overdispersion=0.05)
    abundant = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.95)
    assert abundant.gdna_strand_overdispersion == pytest.approx(0.05, rel=0.20, abs=0.01)

    s, t = _beta_binom_regions(rng, n_regions=3, depth=400, overdispersion=0.05)
    deep_few = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.95)
    assert abs(deep_few.gdna_strand_overdispersion - 0.05) < 0.05

    s, t = _beta_binom_regions(rng, n_regions=40, depth=2, overdispersion=0.01)
    sparse = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.95)
    assert not sparse.fallback_used
    assert 0.0 <= sparse.gdna_strand_overdispersion <= _MAX_OVERDISPERSION


def test_overdispersion_clamped_to_ceiling():
    """A wildly overdispersed fit is capped at the Beta(2,2) ceiling (od=0.2)."""
    rng = np.random.default_rng(3)
    s, t = _beta_binom_regions(rng, n_regions=2000, depth=120, overdispersion=0.45)
    model = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.95)
    assert model.gdna_strand_overdispersion <= _MAX_OVERDISPERSION + 1e-12
    assert model.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION, abs=1e-9)


def test_a_component_that_measured_NOTHING_returns_the_ceiling():
    """No seed signal ⇒ the ceiling and ZERO information, which is what lets `reconcile_overdispersions`
    hand this component the OTHER one's measured value instead of a conjured constant."""
    m = fit_gdna_strand_overdispersion(np.array([]), np.array([]), rna_sense_frac=0.95)
    assert m.fallback_used
    assert m.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)
    assert m.information == 0.0
