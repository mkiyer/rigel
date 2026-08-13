"""Tests for the gDNA strand Beta-Binomial overdispersion fit (docs/em_strand/03 §4.1-4.3).

Core property: data generated with a known overdispersion is recovered by the estimator.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.gdna_strand import (
    _PRIOR_OVERDISPERSION,
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
    """§4.1 — pure-gDNA seeds (weight=1): recovered overdispersion ≈ truth."""
    rng = np.random.default_rng(12345)
    sense, total = _beta_binom_regions(rng, n_regions=4000, depth=120, overdispersion=true_od)
    weight = np.ones_like(total)  # pure gDNA → rna_sense_frac irrelevant
    model = fit_gdna_strand_overdispersion(sense, total, weight, rna_sense_frac=0.95)
    assert not model.fallback_used
    # MoM on 4000 regions: relative error should be small; absolute floor for tiny od.
    assert model.gdna_strand_overdispersion == pytest.approx(true_od, rel=0.20, abs=0.005)


def test_binomial_limit_recovers_near_zero():
    """Zero overdispersion (Binomial) → fit returns ≈ floor, not a spurious positive."""
    rng = np.random.default_rng(7)
    # Binomial: shared rate exactly ½ (no Beta spread).
    total = np.full(4000, 120, dtype=np.float64)
    sense = rng.binomial(120, 0.5, size=4000).astype(np.float64)
    model = fit_gdna_strand_overdispersion(sense, total, np.ones_like(total), rna_sense_frac=0.95)
    assert model.gdna_strand_overdispersion < 0.01


@pytest.mark.parametrize("weight", [0.8, 0.5])
def test_mixture_identifiability(weight):
    """§4.2 — seeds contaminated by stranded RNA (weight<1): overdispersion still recovered."""
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
    model = fit_gdna_strand_overdispersion(
        sense, total, np.full(n_regions, weight), rna_sense_frac=kappa
    )
    assert not model.fallback_used
    assert model.gdna_strand_overdispersion == pytest.approx(true_od, rel=0.30, abs=0.02)


def test_thin_seed_fallback():
    """§4.3 — no gDNA signal (all weight 0, or empty) → fallback to the PRIOR, no crash.

    ⭐ 2026-07-28: the no-argument fallback is now ``od₀`` (0.0345), not a hard ``0``. A hard zero is a
    claim of perfect Binomiality — the *most* confident strand likelihood the model can assert — from a
    library that supplied no evidence at all. Falling back to the near-binomial prior is the honest
    behaviour, and it is why the shrinkage is kept rather than replaced by a bare clip."""
    empty = fit_gdna_strand_overdispersion(
        np.array([]), np.array([]), np.array([]), rna_sense_frac=0.95
    )
    assert empty.fallback_used
    assert empty.gdna_strand_overdispersion == pytest.approx(_PRIOR_OVERDISPERSION)

    rng = np.random.default_rng(1)
    total = np.full(10, 100.0)
    sense = rng.binomial(100, 0.95, size=10).astype(np.float64)
    no_gdna = fit_gdna_strand_overdispersion(
        sense, total, np.zeros_like(total), rna_sense_frac=0.95
    )
    assert no_gdna.fallback_used
    assert no_gdna.gdna_strand_overdispersion == pytest.approx(_PRIOR_OVERDISPERSION)


def test_beta_concentration_roundtrip():
    """Model exposes the Beta(a, a) concentration consistent with a = ½(1−od)/od."""
    m = GdnaStrandModel(
        gdna_strand_overdispersion=0.1, n_seed_regions=10, n_seed_fragments=1000, fallback_used=False
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
    """Mock with only contained-REGION signal: the edge axis carries ZERO counts ⇒ no edge seeds.

    ⚠ One reference with ``n`` regions owns exactly ``n − 1`` lines, so the edge arrays are sized from
    the region axis rather than left empty — an edge axis inconsistent with its own ``ref_id`` is not a
    "no edges" fixture, it is a mis-shaped one.
    """
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n = pos.shape[0]
    n_edges = max(n - 1, 0)
    mass = pos + neg
    substrate = SimpleNamespace(
        region_contained=_view(pos, neg),
        edge_unspliced=_view(np.zeros(n_edges), np.zeros(n_edges)),
    )
    region_arrays = SimpleNamespace(strand_class=np.asarray(ts), ref_id=np.zeros(n, dtype=np.int64))
    # The seed weight is count_gdna_frac directly (was count_evidence/mass) — preserve the intent.
    ce = np.asarray(count_evidence, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(np.where(mass > 0.0, ce / mass, 0.0), 0.0, 1.0)
    region_density = SimpleNamespace(
        count_evidence=ce,
        count_gdna_frac=count_gdna_frac,
        region_count_observable=np.asarray(observable, dtype=bool),
        edge_count_observable=np.zeros(n_edges, dtype=bool),  # no edge seeds
        density=np.zeros(n, dtype=np.float64),
    )
    return substrate, region_arrays, region_density


def test_wrapper_recovers_overdispersion_from_intergenic_seeds():
    """End-to-end of the Phase-2 wrapper: pure-gDNA (intergenic) seed regions → recovered od."""
    rng = np.random.default_rng(2024)
    true_od = 0.10
    sense, total = _beta_binom_regions(rng, n_regions=3000, depth=100, overdispersion=true_od)
    pos, neg = sense, total - sense
    n = len(total)
    substrate, region_arrays, region_density = _mock_substrate(
        pos, neg, np.full(n, TS_NONE), count_evidence=total.copy(), observable=np.ones(n, bool)
    )  # count_evidence == mass ⇒ weight ≈ 1 (pure gDNA)
    model = fit_gdna_strand_from_substrate(
        substrate, region_arrays, region_density, rna_sense_frac=0.95
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
            s_ambig, ra_ambig, nd_ambig, rna_sense_frac=0.95
        ).n_seed_regions
        == 0
    )
    # intergenic but non-observable (exonic) → excluded
    s_exon, ra_exon, nd_exon = _mock_substrate(
        pos, neg, np.full(n, TS_NONE), np.full(n, 100.0), np.zeros(n, bool)
    )
    assert fit_gdna_strand_from_substrate(
        s_exon, ra_exon, nd_exon, rna_sense_frac=0.95
    ).fallback_used


# --- edge seeds: ONE per line, not two per boundary -------------------------------------------


def _edge_parts(signatures, edge_pos, edge_neg, ref_id=None):
    """A substrate + geometry-free ``RegionDensity`` over ``signatures``, on the region/edge axes."""
    from rigel.calibration.density_model import count_observable_masks
    from rigel.calibration.signature import transcript_strand_class

    sig = np.asarray(signatures, dtype=np.uint8)
    n = sig.shape[0]
    ref = np.zeros(n, dtype=np.int64) if ref_id is None else np.asarray(ref_id, dtype=np.int64)
    region_obs, edge_obs = count_observable_masks(sig, ref)
    substrate = SimpleNamespace(
        region_contained=_view(np.zeros(n), np.zeros(n)),
        edge_unspliced=_view(edge_pos, edge_neg),
    )
    region_arrays = SimpleNamespace(
        signature=sig, strand_class=transcript_strand_class(sig.astype(np.int64)), ref_id=ref
    )
    region_density = SimpleNamespace(
        count_gdna_frac=np.zeros(n),
        region_count_observable=region_obs,
        edge_count_observable=edge_obs,
        density=np.zeros(n),
    )
    return substrate, region_arrays, region_density


def test_edge_seeds_emits_ONE_seed_per_line_not_two_per_boundary():
    """⭐ The S5.f collapse, and the whole point of the change.

    The predecessor emitted TWO seeds for one seam — region ``r``'s right side and region ``r+1``'s
    left side — because the old accumulator split one crossing fragment's mass across the two flanks.
    They were the same physical crossing, counted twice into one pooled moment estimator, which
    inflates its apparent sample size by 2× and correlates every pair perfectly. A contiguous edge is
    a 0-bp line with ONE count, so there is ONE seed.
    """
    from rigel.calibration.strand_deconv import edge_seeds

    # intron+ | intron+ : one line, count-observable (no shared EXON bit), sense-POS on both flanks.
    substrate, region_arrays, region_density = _edge_parts(
        [BIT_INTRON_POS, BIT_INTRON_POS], edge_pos=[70.0], edge_neg=[30.0]
    )
    sense, total, weight = edge_seeds(substrate, region_arrays, region_density)
    assert sense.shape == (1,)  # ⛔ was (2,) — the same crossing, twice
    np.testing.assert_allclose(sense, [70.0])
    np.testing.assert_allclose(total, [100.0])
    np.testing.assert_allclose(weight, [1.0])


def test_edge_seed_sense_follows_the_flanking_transcript_strand():
    """A NEG-strand line orients to the NEG genome column; the sense count is not always ``pos``."""
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [BIT_INTRON_NEG, BIT_INTRON_NEG], edge_pos=[70.0], edge_neg=[30.0]
    )
    sense, total, _ = edge_seeds(substrate, region_arrays, region_density)
    np.testing.assert_allclose(sense, [30.0])
    np.testing.assert_allclose(total, [100.0])


def test_an_intergenic_flank_is_a_strand_WILDCARD():
    """Intergenic carries no transcript, so a gene-edge line is oriented by its gene flank."""
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [0, BIT_INTRON_NEG], edge_pos=[70.0], edge_neg=[30.0]
    )
    sense, _, _ = edge_seeds(substrate, region_arrays, region_density)
    np.testing.assert_allclose(sense, [30.0])  # oriented NEG by the gene side


def test_an_opposite_strand_line_is_not_strand_observable():
    """``{POS, NEG}`` leaves 'sense' undefined, so the line cannot seed the fit at all."""
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [BIT_INTRON_POS, BIT_INTRON_NEG], edge_pos=[70.0], edge_neg=[30.0]
    )
    sense, _, _ = edge_seeds(substrate, region_arrays, region_density)
    assert sense.shape == (0,)


def test_an_AMBIG_flank_cannot_seed():
    """⛔ Found by PERTURBATION, not by review: nothing pinned this rule.

    A flank carrying overlapping ± transcripts has no defined transcript sense, so neither genome
    column is "sense" and the line cannot seed a strand fit. ``edge_strand_orientation`` used to carry
    an explicit ``~either_ambig`` guard — deleting it broke NOTHING, because ``TS_AMBIG`` is a fourth
    distinct value that the ``POS-or-NONE`` / ``NEG-or-NONE`` tests already exclude. The guard was
    dead code claiming to be the rule; this test is the rule.

    ⚠ The fixture uses INTRON bits on both strands, not exon bits: an AMBIG flank must still be
    count-observable, or the test would pass for the wrong reason.
    """
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [BIT_INTRON_POS | BIT_INTRON_NEG, BIT_INTRON_POS], edge_pos=[70.0], edge_neg=[30.0]
    )
    assert region_density.edge_count_observable[0]  # count-observable — it fails on STRAND alone
    sense, _, _ = edge_seeds(substrate, region_arrays, region_density)
    assert sense.shape == (0,)


def test_a_line_inside_one_exon_is_not_count_observable():
    """A shared exon bit means an exon-strand continues across the line, so unspliced MATURE RNA
    crosses it and its count is not gDNA."""
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [BIT_EXON_POS, BIT_EXON_POS], edge_pos=[70.0], edge_neg=[30.0]
    )
    sense, _, _ = edge_seeds(substrate, region_arrays, region_density)
    assert sense.shape == (0,)


def test_edge_seeds_never_straddle_a_reference():
    """Two single-region references own ZERO lines between them — nothing can leak across."""
    from rigel.calibration.strand_deconv import edge_seeds

    substrate, region_arrays, region_density = _edge_parts(
        [BIT_INTRON_POS, BIT_INTRON_POS], edge_pos=[], edge_neg=[], ref_id=[0, 1]
    )
    sense, _, _ = edge_seeds(substrate, region_arrays, region_density)
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


def test_shrinkage_sparse_leans_on_prior_abundant_on_fit():
    """LOW-INFORMATION seeds → shrink toward the prior; high-information → follow the fitted MoM.

    ⭐ 2026-07-28: "sparse" is now measured in the right currency. Overdispersion is a correlation
    BETWEEN fragments, so the evidence unit is a PAIR — a seed of one fragment carries none of it. This
    test used to call 3 regions of depth 100 "sparse" and expect the prior to win; that is 14,850 pairs and
    the fit should win, which is exactly the confusion the information-weighted shrinkage fixes. Genuine
    sparsity is many SHALLOW seeds."""
    rng = np.random.default_rng(11)
    prior = overdispersion_for_beta(3.0)  # 1/7 ≈ 0.143
    weight = 909.0  # the derived prior weight, in information units

    # LOW information: 40 two-fragment seeds = 40 pairs, far below the prior's 909 ⇒ ≈ prior.
    s, t = _beta_binom_regions(rng, n_regions=40, depth=2, overdispersion=0.01)
    sparse = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert sparse.gdna_strand_overdispersion == pytest.approx(prior, abs=0.03)

    # HIGH information: 5000 regions at depth 120 = 35.7M pairs ⇒ the fit dominates.
    s, t = _beta_binom_regions(rng, n_regions=5000, depth=120, overdispersion=0.05)
    abundant = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert abundant.gdna_strand_overdispersion == pytest.approx(0.05, rel=0.20, abs=0.01)

    # ⭐ And the point of the units fix: 3 DEEP regions are NOT sparse — they must follow the fit.
    s, t = _beta_binom_regions(rng, n_regions=3, depth=400, overdispersion=0.05)
    deep_few = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert abs(deep_few.gdna_strand_overdispersion - 0.05) < abs(
        deep_few.gdna_strand_overdispersion - prior
    )


def test_overdispersion_clamped_to_ceiling():
    """A wildly overdispersed fit is capped at the Beta(2,2) ceiling (od=0.2)."""
    rng = np.random.default_rng(3)
    s, t = _beta_binom_regions(rng, n_regions=2000, depth=120, overdispersion=0.45)
    model = fit_gdna_strand_overdispersion(s, t, np.ones_like(t), rna_sense_frac=0.95)
    assert model.gdna_strand_overdispersion <= _MAX_OVERDISPERSION + 1e-12
    assert model.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION, abs=1e-9)


def test_fallback_returns_prior():
    """No seed signal (empty) ⇒ fall back to the prior overdispersion, not 0."""
    prior = overdispersion_for_beta(3.0)
    m = fit_gdna_strand_overdispersion(
        np.array([]),
        np.array([]),
        np.array([]),
        rna_sense_frac=0.95,
        prior_overdispersion=prior,
        prior_weight=30.0,
    )
    assert m.fallback_used
    assert m.gdna_strand_overdispersion == pytest.approx(prior)
