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
from rigel.calibration.signature import TS_AMBIG, TS_NONE
from rigel.calibration.strand_likelihood import strand_loglik


def _beta_binom_nodes(rng, n_nodes, depth, overdispersion, mean=0.5):
    """Draw per-node (sense, total) from BetaBinom(depth, mean, overdispersion).

    Mirrors the simulator: each node has a shared latent sense rate p ~ Beta(a, b),
    then sense ~ Binomial(depth, p). With a symmetric mean=½, a = b = ½(1−od)/od.
    """
    conc = (1.0 - overdispersion) / overdispersion  # a + b
    a = mean * conc
    b = (1.0 - mean) * conc
    p = rng.beta(a, b, size=n_nodes)
    total = np.full(n_nodes, depth, dtype=np.float64)
    sense = rng.binomial(depth, p).astype(np.float64)
    return sense, total


@pytest.mark.parametrize("true_od", [0.01, 0.05, 0.10, 0.20])
def test_recovers_overdispersion_pure_gdna(true_od):
    """§4.1 — pure-gDNA seeds (weight=1): recovered overdispersion ≈ truth."""
    rng = np.random.default_rng(12345)
    sense, total = _beta_binom_nodes(rng, n_nodes=4000, depth=120, overdispersion=true_od)
    weight = np.ones_like(total)  # pure gDNA → rna_sense_frac irrelevant
    model = fit_gdna_strand_overdispersion(sense, total, weight, rna_sense_frac=0.95)
    assert not model.fallback_used
    # MoM on 4000 nodes: relative error should be small; absolute floor for tiny od.
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
    n_nodes = 6000
    conc = (1.0 - true_od) / true_od
    a = 0.5 * conc
    sense = np.empty(n_nodes)
    total = np.full(n_nodes, depth, dtype=np.float64)
    for i in range(n_nodes):
        n_g = rng.binomial(depth, weight)
        p_g = rng.beta(a, a)  # shared gDNA rate for this node
        sense[i] = rng.binomial(n_g, p_g) + rng.binomial(depth - n_g, kappa)
    model = fit_gdna_strand_overdispersion(
        sense, total, np.full(n_nodes, weight), rna_sense_frac=kappa
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
        gdna_strand_overdispersion=0.1, n_seed_nodes=10, n_seed_fragments=1000, fallback_used=False
    )
    assert m.beta_concentration() == pytest.approx(0.5 * (1.0 - 0.1) / 0.1)


# --- the substrate wrapper (Phase 2 extraction + fit) ---------------------------------------


def _zero_view(n):
    z = np.zeros(n, dtype=np.float64)
    return SimpleNamespace(
        n_unspliced_pos=z.copy(),
        n_unspliced_neg=z.copy(),
        mass_unspliced=z.copy(),
        mass_spliced=z.copy(),
    )


def _mock_substrate(pos, neg, ts, count_evidence, observable):
    """Mock with only contained-region signal (boundary sides zeroed ⇒ no boundary seeds)."""
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n = pos.shape[0]
    mass = pos + neg
    contained = SimpleNamespace(n_unspliced_pos=pos, n_unspliced_neg=neg, mass_unspliced=mass)
    substrate = SimpleNamespace(contained=contained, left=_zero_view(n), right=_zero_view(n))
    region_arrays = SimpleNamespace(strand_class=np.asarray(ts), ref_id=np.zeros(n, dtype=np.int64))
    # New seed weight is count_gdna_frac directly (was count_evidence/mass) — preserve the intent.
    ce = np.asarray(count_evidence, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_gdna_frac = np.clip(np.where(mass > 0.0, ce / mass, 0.0), 0.0, 1.0)
    node_density = SimpleNamespace(
        count_evidence=ce,
        count_gdna_frac=count_gdna_frac,
        region_count_observable=np.asarray(observable, dtype=bool),
        boundary_count_observable=np.zeros(n, dtype=bool),  # no boundary seeds
        density=np.zeros(n, dtype=np.float64),
    )
    return substrate, region_arrays, node_density


def test_wrapper_recovers_overdispersion_from_intergenic_seeds():
    """End-to-end of the Phase-2 wrapper: pure-gDNA (intergenic) seed regions → recovered od."""
    rng = np.random.default_rng(2024)
    true_od = 0.10
    sense, total = _beta_binom_nodes(rng, n_nodes=3000, depth=100, overdispersion=true_od)
    pos, neg = sense, total - sense
    n = len(total)
    substrate, region_arrays, node_density = _mock_substrate(
        pos, neg, np.full(n, TS_NONE), count_evidence=total.copy(), observable=np.ones(n, bool)
    )  # count_evidence == mass ⇒ weight ≈ 1 (pure gDNA)
    model = fit_gdna_strand_from_substrate(
        substrate, region_arrays, node_density, np.zeros(n), rna_sense_frac=0.95
    )
    assert not model.fallback_used
    assert model.n_seed_nodes == n
    assert model.gdna_strand_overdispersion == pytest.approx(true_od, rel=0.25, abs=0.01)


def test_wrapper_excludes_ambig_and_non_observable():
    """AMBIG nodes (no defined sense) and non-count-observable (exonic) nodes are dropped."""
    n = 100
    pos = np.full(n, 50.0)
    neg = np.full(n, 50.0)
    eff = np.zeros(n)
    # all AMBIG → excluded even though count-observable
    s_ambig, ra_ambig, nd_ambig = _mock_substrate(
        pos, neg, np.full(n, TS_AMBIG), np.full(n, 100.0), np.ones(n, bool)
    )
    assert (
        fit_gdna_strand_from_substrate(
            s_ambig, ra_ambig, nd_ambig, eff, rna_sense_frac=0.95
        ).n_seed_nodes
        == 0
    )
    # intergenic but non-observable (exonic) → excluded
    s_exon, ra_exon, nd_exon = _mock_substrate(
        pos, neg, np.full(n, TS_NONE), np.full(n, 100.0), np.zeros(n, bool)
    )
    assert fit_gdna_strand_from_substrate(
        s_exon, ra_exon, nd_exon, eff, rna_sense_frac=0.95
    ).fallback_used


def test_boundary_side_seeds_extracts_observable_seam_sides():
    """boundary_side_seeds picks the count- & strand-observable seam sides with count weights."""
    from rigel.calibration.signature import TS_POS
    from rigel.calibration.strand_deconv import boundary_side_seeds

    def view(pos, neg):
        pos = np.asarray(pos, dtype=np.float64)
        neg = np.asarray(neg, dtype=np.float64)
        return SimpleNamespace(
            n_unspliced_pos=pos,
            n_unspliced_neg=neg,
            mass_unspliced=pos + neg,
            mass_spliced=np.zeros_like(pos),
        )

    # two adjacent POS regions, one ref; boundary (0,1) is a count-observable seam.
    left = view([0, 70], [0, 30])  # region 1's left side = right side of boundary (0,1)
    right = view([80, 0], [20, 0])  # region 0's right side = left side of boundary (0,1)
    substrate = SimpleNamespace(left=left, right=right, contained=None)
    region_arrays = SimpleNamespace(
        strand_class=np.array([TS_POS, TS_POS]), ref_id=np.array([0, 0])
    )
    node_density = SimpleNamespace(
        boundary_count_observable=np.array([True, False]), density=np.zeros(2)
    )
    sense, total, weight = boundary_side_seeds(
        substrate, region_arrays, node_density, np.array([100.0, 100.0])
    )
    assert sorted(np.round(sense).tolist()) == [70.0, 80.0]
    assert sorted(total.tolist()) == [100.0, 100.0]
    # weight = count_gdna_frac: a count-observable seam side reads its own crossing density
    # (mass/eff), so count_gdna_frac = clip((mass/eff)·eff/mass) = 1 (gDNA by signature).
    assert np.allclose(weight, 1.0)


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


# --- prior conversion + shrinkage (replaces the old min-node / significance gates) -----------


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
    test used to call 3 nodes of depth 100 "sparse" and expect the prior to win; that is 14,850 pairs and
    the fit should win, which is exactly the confusion the information-weighted shrinkage fixes. Genuine
    sparsity is many SHALLOW seeds."""
    rng = np.random.default_rng(11)
    prior = overdispersion_for_beta(3.0)  # 1/7 ≈ 0.143
    weight = 909.0  # the derived prior weight, in information units

    # LOW information: 40 two-fragment seeds = 40 pairs, far below the prior's 909 ⇒ ≈ prior.
    s, t = _beta_binom_nodes(rng, n_nodes=40, depth=2, overdispersion=0.01)
    sparse = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert sparse.gdna_strand_overdispersion == pytest.approx(prior, abs=0.03)

    # HIGH information: 5000 nodes at depth 120 = 35.7M pairs ⇒ the fit dominates.
    s, t = _beta_binom_nodes(rng, n_nodes=5000, depth=120, overdispersion=0.05)
    abundant = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert abundant.gdna_strand_overdispersion == pytest.approx(0.05, rel=0.20, abs=0.01)

    # ⭐ And the point of the units fix: 3 DEEP nodes are NOT sparse — they must follow the fit.
    s, t = _beta_binom_nodes(rng, n_nodes=3, depth=400, overdispersion=0.05)
    deep_few = fit_gdna_strand_overdispersion(
        s, t, np.ones_like(t), rna_sense_frac=0.95, prior_overdispersion=prior, prior_weight=weight
    )
    assert abs(deep_few.gdna_strand_overdispersion - 0.05) < abs(
        deep_few.gdna_strand_overdispersion - prior
    )


def test_overdispersion_clamped_to_ceiling():
    """A wildly overdispersed fit is capped at the Beta(2,2) ceiling (od=0.2)."""
    rng = np.random.default_rng(3)
    s, t = _beta_binom_nodes(rng, n_nodes=2000, depth=120, overdispersion=0.45)
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
