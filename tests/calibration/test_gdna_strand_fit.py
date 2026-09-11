"""The gDNA strand Beta-Binomial overdispersion fit, and everything that feeds it.

No structural class of objects can be assumed PURE gDNA — purity is a property of the annotation and of
the sample rather than of the genome — so the estimator trusts no class
(TRAPS: purity-is-a-property-of-the-annotation). It rests on a lemma instead: on a strand-specific
library, RNA of a seed's own gene can only pull that seed's sense fraction TOWARD κ, so orienting
``d = K − N/2`` such that RNA pulls ``d`` negative makes the moment excess ``(d² − N/4)`` an EVEN
function of ``d``, and over the AWAY half (``d > 0``, ties at weight ½) the pooled method-of-moments
ratio is unbiased for ρ_g whatever the seeds' RNA content — a contaminated seed reaches the away side
only by noise, with a small ``d``, biasing ρ̂ DOWN and never up. The seed set is every count- and
strand-observable GENIC object (intron regions; exon|intron and gene-edge boundaries) in whatever GTF
is supplied: an INTERGENIC seed has no gene strand to orient by and an AMBIG object has no defined
sense, so both are out, while unannotated ANTISENSE transcription pushes toward the away side and
inflates ρ̂ — a known limit rather than a hidden one. The blocks below gate that moment, the RNA
pooled moment, the reconciliation between the two components, the influence weighting, the seed
selectors, the conversion the strand likelihood consumes, and the end-to-end recovery of a planted
overdispersion through ``calibrate``. PERTURBATION: forcing the full two-sided moment back fires the
contamination gate.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from _synthetic import make_gdna_fl_pmf, make_strand_models

from rigel.calibration import calibrate
from rigel.calibration.gdna_strand import (
    _MAX_OVERDISPERSION,
    GdnaStrandModel,
    _null_information,
    away_half_moment,
    between_seed_variance,
    boundary_seeds,
    fit_gdna_strand_from_substrate,
    fit_gdna_strand_overdispersion,
    influence_weights,
    overdispersion_for_beta,
    reconcile_overdispersions,
    seed_participation,
)
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_AMBIG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.strand_likelihood import strand_loglik
from rigel.config import CalibrationConfig
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    AccumulatorPayload,
    DeferredFragments,
    GapCensus,
    ScanQC,
)

KAPPA = 0.01  # the RNA sense fraction on an R1-antisense library


def _bb(rng, n, mu, rho, size):
    if rho <= 0:
        return rng.binomial(n, mu, size)
    a, b = mu * (1 / rho - 1), (1 - mu) * (1 / rho - 1)
    return rng.binomial(n, rng.beta(a, b, size))


def _seeds(rng, *, n_pure=3000, od_true=0.05, contam=0.0, w_contam=None, antisense=False):
    """Pure BB(N, ½, od_true) seeds at N ~ U[3, 400], plus a fraction ``contam`` of seeds that ALSO hold RNA
    with gDNA share ``w`` (U(0,1) unless given), sense κ (or 1−κ for antisense)."""
    N = rng.integers(3, 400, n_pure)
    sense = _bb(rng, N, 0.5, od_true, n_pure).astype(float)
    total = N.astype(float)
    m = int(round(n_pure * contam))
    if m:
        Nc = rng.integers(3, 400, m)
        w = np.full(m, w_contam) if w_contam is not None else rng.uniform(0, 1, m)
        g = rng.binomial(Nc, w)
        k = (1 - KAPPA) if antisense else KAPPA
        sc = _bb(rng, g, 0.5, od_true, m) + rng.binomial(Nc - g, k)
        sense = np.concatenate([sense, sc.astype(float)])
        total = np.concatenate([total, Nc.astype(float)])
    return sense, total


def _away_sd(total, kappa=KAPPA):
    n = np.asarray(total, float)
    return 1.0 / np.sqrt(0.5 * _null_information(n, 0.25))


# ── the lemma ──────────────────────────────────────────────────────────────────────────────────


def test_pure_seeds_recover_the_truth():
    rng = np.random.default_rng(1)
    s, t = _seeds(rng, od_true=0.05)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert not m.fallback_used
    assert m.gdna_strand_overdispersion == pytest.approx(0.05, abs=3 * _away_sd(t) + 0.003)


@pytest.mark.parametrize("contam,w", [(0.3, None), (0.3, 0.1), (0.6, None)])
def test_same_strand_contamination_of_any_amount_does_not_move_the_fit(contam, w):
    """30 % of seeds; a nascent-dominated 30 % (w = 0.1); a MAJORITY (60 %) — the away half is unmoved."""
    rng = np.random.default_rng(2)
    s, t = _seeds(rng, od_true=0.05, contam=contam, w_contam=w)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert abs(m.gdna_strand_overdispersion - 0.05) < 3 * _away_sd(t) + 0.004


def test_perturbation_the_two_sided_moment_reads_contamination_as_overdispersion():
    """PERTURB THE FIX: the full (two-sided) moment on the same seeds — the previous premise — is pulled far
    above the truth, which is the defect this estimator exists to remove."""
    rng = np.random.default_rng(2)
    s, t = _seeds(rng, od_true=0.05, contam=0.3, w_contam=0.1)
    d = s - t / 2
    full = float(((d**2 - t / 4).sum()) / ((t * (t - 1) / 4).sum()))
    away = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA).gdna_strand_overdispersion
    assert full > 0.15 and abs(away - 0.05) < 0.02


def test_truth_zero_reads_zero_under_contamination():
    """The simulator's world: od 0 with nascent-laden seeds must fit ≈ 0 (clipped at the floor)."""
    rng = np.random.default_rng(3)
    s, t = _seeds(rng, od_true=0.0, contam=0.4)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert m.gdna_strand_overdispersion < 3 * _away_sd(t) + 0.002


def test_antisense_contamination_is_the_recorded_limit():
    """Unannotated ANTISENSE RNA pushes toward the away side and INFLATES ρ̂ — recorded, not hidden."""
    rng = np.random.default_rng(4)
    s, t = _seeds(rng, od_true=0.05, contam=0.1, antisense=True)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert m.gdna_strand_overdispersion > 0.05 + 3 * _away_sd(t)


def test_orientation_follows_kappa_on_an_r1_sense_library():
    """On a library where RNA reads land on the SENSE strand (κ ≈ 1), RNA pulls d POSITIVE and the away half
    is d < 0. The lemma is symmetric in κ ↔ 1−κ; the fit must be too."""
    rng = np.random.default_rng(5)
    N = rng.integers(3, 400, 2000)
    pure = _bb(rng, N, 0.5, 0.05, 2000).astype(float)
    Nc = rng.integers(3, 400, 600)
    g = rng.binomial(Nc, rng.uniform(0, 1, 600))
    contam = (_bb(rng, g, 0.5, 0.05, 600) + rng.binomial(Nc - g, 0.99)).astype(float)
    s, t = np.concatenate([pure, contam]), np.concatenate([N, Nc]).astype(float)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.99)
    assert abs(m.gdna_strand_overdispersion - 0.05) < 3 * _away_sd(t) + 0.004


def test_a_tie_enters_at_half_weight_exactly():
    """K = N/2 is a tie: it belongs to neither half, so it carries weight ½. The N = 2 case is exact —
    under BetaBinom(2, ½, ρ): P(K=1) = ½(1−ρ), P(K=0) = P(K=2) = (1+ρ)/4 — and at those frequencies the
    away half returns EXACTLY ρ, while full tie weight returns (3ρ−1)/(3−ρ) < 0 for ρ < ⅓ (clipped to 0)
    and zero tie weight returns 1. Watched: the tie perturbation fired NO other gate."""
    rho = 0.1
    n_pairs = 4000
    k0 = k2 = int(round(n_pairs * (1 + rho) / 4))
    k1 = n_pairs - k0 - k2
    sense = np.array([0.0] * k0 + [1.0] * k1 + [2.0] * k2)
    total = np.full(sense.shape, 2.0)
    od, _info = away_half_moment(sense, total, KAPPA)
    assert od == pytest.approx(rho, abs=1e-9)


def test_the_information_is_HALF_THE_TOTAL_PAIR_COUNT_not_half_the_away_halfs():
    """The double-halving. ``information`` must be ½ × the pair count of ALL seeds, because the away
    half's own pair count is ALREADY about half of that — halving it again understates the information 2×
    and overstates the standard error by √2.

    The property, not the implementation: ``Var(od_mom)|₀ = 2/P`` for the total pair count ``P``. Derivation
    — ``Var(e_s)|₀ = n(n−1)/8`` and ``E[a_s] = ½`` (the sign of a symmetric residual is independent of its
    size), so ``Var(num) = P/8``, ``E[den] = P/4`` and ``Var(num/den) = 2/P``. Asserted here BOTH ways: the
    closed form, and a Monte-Carlo null whose empirical sd must match ``1/sqrt(information)``.

    PERTURBATION: a gate asserting ``info == 0.5·_null_information(total[sense >= 1])`` restates the
    implementation's own expression, so it certifies the defect and fires on the repair instead
    (TRAPS: a-gate-that-restates-the-implementation)."""
    rng = np.random.default_rng(21)
    total = rng.integers(2, 40, 3000).astype(float)
    sense = rng.binomial(total.astype(int), 0.5).astype(float)
    _od, info = away_half_moment(sense, total, KAPPA)
    pairs = float((total * (total - 1.0) / 2.0).sum())
    assert info == pytest.approx(0.5 * pairs, rel=1e-12)  # NOT the away half's own pair count

    # the Monte-Carlo null: 1/sqrt(info) IS the standard error it claims to be
    n = np.full(2000, 5.0)
    draws = rng.binomial(5, 0.5, size=(500, 2000)).astype(float)
    ods = np.array([away_half_moment(draws[i], n, KAPPA)[0] for i in range(draws.shape[0])])
    _od0, info0 = away_half_moment(draws[0], n, KAPPA)
    assert ods.std() == pytest.approx(1.0 / np.sqrt(info0), rel=0.15)


def test_an_unstranded_library_uses_the_FULL_moment_and_is_not_forced_to_zero():
    """κ = ½ exactly is reachable and must not collapse the fit. ``rna_sense_frac`` is the posterior mean
    ``(n_same + 1)/(n_obs + 2)``, which is exactly ½ whenever ``2·n_same = n_obs`` — the modal outcome on an
    unstranded library. There ``sign(½ − κ) = 0``, every oriented residual is 0, every seed is a tie with
    excess ``−N/4``, and the away-half branch would return a hard ``od = 0`` — perfect Binomiality, the most
    confident strand likelihood assertable — with ``fallback_used`` still False.

    At κ = ½ the RNA is unstranded too, so its contamination is symmetric and the FULL two-sided moment is
    the right estimator: unbiased on pure seeds, and still one-sided under contamination (a same-mean RNA
    component can only shrink the excess). Both halves are asserted, plus the un-halved information."""
    rng = np.random.default_rng(22)
    s, t = _seeds(rng, od_true=0.05)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=0.5)
    assert not m.fallback_used
    assert m.gdna_strand_overdispersion == pytest.approx(0.05, abs=0.01)  # never a hard 0.0

    _od, info = away_half_moment(s, t, 0.5)
    assert info == pytest.approx(float((t * (t - 1.0) / 2.0).sum()), rel=1e-12)  # FULL, not halved

    # contamination at the same mean ½ can only pull it DOWN
    n_c = rng.integers(3, 400, 1500)
    g = rng.binomial(n_c, rng.uniform(0, 1, 1500))
    sc = _bb(rng, g, 0.5, 0.05, 1500) + rng.binomial(n_c - g, 0.5)
    mixed = fit_gdna_strand_overdispersion(
        np.concatenate([s, sc.astype(float)]),
        np.concatenate([t, n_c.astype(float)]),
        rna_sense_frac=0.5,
    )
    assert mixed.gdna_strand_overdispersion <= 0.05 + 0.005


def test_a_component_with_no_evidence_returns_the_ceiling_and_zero_information():
    """No conjured constant anywhere. With no pair the fit returns the CEILING — the widest strand
    likelihood the model admits, so the channel says nothing — and ZERO information, the signal
    `reconcile_overdispersions` uses to hand it the other component's measured value instead."""
    m = fit_gdna_strand_overdispersion(np.array([]), np.array([]), rna_sense_frac=KAPPA)
    assert m.fallback_used and m.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)
    assert (
        m.information == 0.0
    )  # zero information is what lets the OTHER component supply the value


# ── the seed SELECTORS: genic only; intergenic and AMBIG out ───────────────────────────────────


def _view(pos, neg):
    """A :class:`PopulationView`-shaped stand-in: ONE genome-strand count array, nothing else.

    There is no separate mass array: the accumulator deposits ``+1`` on every object the fragment
    touched, so ``count`` IS the mass and there is nothing to keep consistent with it.
    """
    return SimpleNamespace(
        count=np.stack(
            [np.asarray(pos, dtype=np.float64), np.asarray(neg, dtype=np.float64)], axis=1
        )
    )


def _parts(signatures, boundary_pos, boundary_neg, region_pos=None, region_neg=None):
    from rigel.calibration.density_model import count_observable_masks
    from rigel.calibration.signature import transcript_strand_class

    sig = np.asarray(signatures, dtype=np.uint8)
    n = sig.shape[0]
    ref = np.zeros(n, dtype=np.int64)
    region_obs, boundary_obs = count_observable_masks(sig, ref)
    rp = np.zeros(n) if region_pos is None else np.asarray(region_pos, float)
    rn = np.zeros(n) if region_neg is None else np.asarray(region_neg, float)
    substrate = SimpleNamespace(
        region_contained=_view(rp, rn), boundary_unspliced=_view(boundary_pos, boundary_neg)
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


def test_exon_intron_and_gene_edge_boundaries_both_seed():
    """intergenic | exon+ | intron+ : the gene edge AND the exon|intron boundary are seeds — no purity class
    decides membership; the away half handles what they hold."""
    substrate, ra, rd = _parts(
        [0, BIT_EXON_POS, BIT_INTRON_POS], boundary_pos=[60.0, 70.0], boundary_neg=[40.0, 30.0]
    )
    sense, total = boundary_seeds(substrate, ra, rd.boundary_count_observable)
    np.testing.assert_allclose(sense, [60.0, 70.0])
    np.testing.assert_allclose(total, [100.0, 100.0])


def test_an_ambig_flank_cannot_seed():
    """Annotated sense + antisense at one place has no defined sense: no orientation, no seed."""
    substrate, ra, rd = _parts(
        [BIT_INTRON_POS | BIT_INTRON_NEG, BIT_INTRON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, ra, rd.boundary_count_observable)
    assert sense.shape == (0,)


def test_intergenic_regions_are_not_seeds_but_intron_regions_are():
    """An intergenic region has no gene strand to orient by — the lemma cannot protect it — so it is OUT of
    the od fit; an intron region of a gene is IN, with no weight."""
    rng = np.random.default_rng(6)
    n = 40
    pos = rng.binomial(200, 0.5, n).astype(float)
    neg = 200 - pos
    sig = np.array([0] * 20 + [BIT_INTRON_POS] * 20, dtype=np.uint8)
    substrate, ra, rd = _parts(
        sig,
        boundary_pos=np.zeros(n - 1),
        boundary_neg=np.zeros(n - 1),
        region_pos=pos,
        region_neg=neg,
    )
    m = fit_gdna_strand_from_substrate(
        substrate,
        ra,
        region_count_observable=rd.region_count_observable,
        boundary_count_observable=rd.boundary_count_observable,
        rna_sense_frac=KAPPA,
    )
    assert m.n_seed_regions == 20  # the 20 intron regions; the 20 intergenic regions do not seed


def test_a_NEG_strand_genes_seed_orients_to_the_NEG_column():
    """PERTURBATION: the whole calibration suite passes with ``_region_seeds``' ``np.where(ts ==
    TS_NEG, neg, pos)`` flip deleted, so nothing but this pins it.

    Counts are stored by GENOME strand; the lemma needs TRANSCRIPT sense. Without the flip a minus-strand
    gene's residual is oriented backwards, so its RNA pulls ``d`` toward the AWAY side instead of away from
    it — the estimator's one-sided guarantee inverts for every NEG gene and contamination INFLATES od.
    Asserted twice: the seed's sense column, and the inversion itself on a contaminated NEG population."""
    n = 30
    pos = np.full(n, 70.0)
    neg = np.full(n, 30.0)
    substrate, ra, rd = _parts(
        np.full(n, BIT_INTRON_NEG, dtype=np.uint8),
        boundary_pos=np.zeros(n - 1),
        boundary_neg=np.zeros(n - 1),
        region_pos=pos,
        region_neg=neg,
    )
    from rigel.calibration.gdna_strand import _region_seeds

    sense, total = _region_seeds(substrate, ra, rd.region_count_observable)
    np.testing.assert_allclose(sense, neg)  # the NEG column is sense here, not `pos`
    np.testing.assert_allclose(total, pos + neg)

    # and the consequence: a NEG gene read with the POS column as "sense" inverts the guarantee
    rng = np.random.default_rng(23)
    s, t = _seeds(rng, od_true=0.0, contam=0.5, w_contam=0.3)
    right = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    wrong = fit_gdna_strand_overdispersion(
        t - s, t, rna_sense_frac=KAPPA
    )  # sense read off the wrong column
    assert right.gdna_strand_overdispersion < 0.01 < wrong.gdna_strand_overdispersion


# ── the QC facts: a clamp must announce itself, and concentration must be visible ───────────────


def test_the_participation_ratio_counts_EFFECTIVE_seeds():
    """``n_eff = (Σ|x|)²/Σx²`` — the seed count when every seed contributes alike, ~1 when one seed IS the
    estimate. Threshold-free: no ``k``, no cutoff, no constant. Both limits asserted exactly."""
    n = np.full(400, 9.0)
    even = np.full(400, 8.0)  # every seed the same distance from n/2 ⇒ identical contributions
    assert seed_participation(even, n, KAPPA) == pytest.approx(400.0, rel=1e-9)

    # one seed at depth 400 among 400 shallow ones: the deep seed's d² dwarfs the rest
    tot = np.concatenate([np.full(400, 9.0), [400.0]])
    sen = np.concatenate([np.full(400, 5.0), [400.0]])
    assert seed_participation(sen, tot, KAPPA) < 1.05

    assert np.isnan(seed_participation(np.array([]), np.array([]), KAPPA))


def test_a_CLAMPED_fit_announces_itself_and_keeps_the_raw_moment():
    """A value at the ceiling is a CLAMP, not a measurement, and real libraries do reach it. The model
    must carry both facts, or a reader takes the ceiling for a fit."""
    rng = np.random.default_rng(31)
    s, t = _seeds(rng, od_true=0.6, n_pure=2000)  # far above the ceiling
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert m.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)
    assert m.clamped_at_ceiling
    assert m.raw_overdispersion > _MAX_OVERDISPERSION

    s, t = _seeds(rng, od_true=0.05, n_pure=3000)
    ok = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert not ok.clamped_at_ceiling
    # ``raw`` is the ρ = 0 moment exactly — the pair-count estimator, before the influence weighting
    assert ok.raw_overdispersion == pytest.approx(away_half_moment(s, t, KAPPA)[0])
    assert ok.effective_seeds > 1.0

    fb = fit_gdna_strand_overdispersion(np.array([]), np.array([]), rna_sense_frac=KAPPA)
    assert fb.fallback_used and not fb.clamped_at_ceiling  # a fallback is not a clamp
    assert np.isnan(fb.raw_overdispersion)


def test_the_participation_ratio_is_WIRED_IN_and_uses_contribution_MAGNITUDE():
    """PERTURBATION, two holes: the preceding gates test the function and the model separately and
    leave both of these open.

    (a) the MODEL must publish the participation ratio, not the seed COUNT: with one dominant seed among
        401 the two differ 385-fold, and substituting the count fired nothing;
    (b) the ratio is over contribution MAGNITUDES. Dropping the ``|·|`` also fired nothing, because every
        earlier fixture had same-signed contributions. On a null population — where the signed sum nearly
        cancels by construction — the two answers are 713 and 0.55.
    """
    tot = np.concatenate([np.full(400, 9.0), [400.0]])
    sen = np.concatenate([np.full(400, 5.0), [400.0]])
    m = fit_gdna_strand_overdispersion(sen, tot, rna_sense_frac=KAPPA)
    assert m.n_seed_regions == 401
    assert m.effective_seeds < 0.05 * m.n_seed_regions  # the count would say 401

    rng = np.random.default_rng(32)
    s, t = _seeds(rng, od_true=0.0, n_pure=4000)
    assert seed_participation(s, t, KAPPA) > 400.0  # a signed sum collapses this to ~0.55

    # (c) — and it counts only the seeds the ESTIMATE uses. 400 fully-ANTISENSE deep seeds carry a huge
    # residual but sit on the RNA side at weight 0, so they are not evidence and must not dilute the
    # concentration; one sense-side seed IS the estimate. Ignoring the away-half weight reports 401.
    tot = np.full(401, 400.0)
    sen = np.concatenate([np.zeros(400), [400.0]])
    assert seed_participation(sen, tot, KAPPA) == pytest.approx(1.0, abs=1e-9)


# ── the influence weighting ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("od_true", [0.0, 0.01, 0.05, 0.2])
def test_the_influence_weighting_is_UNBIASED_at_a_known_overdispersion(od_true):
    """The weights depend only on ``n_s`` and ρ — never on a seed's own data — so the ratio keeps its
    expectation. Efficiency is what the weighting buys; correctness is not what it spends."""
    got = []
    for rep in range(12):
        n = np.random.default_rng(700 + rep).integers(2, 400, 2500)
        s = _bb(np.random.default_rng(900 + rep), n, 0.5, od_true, len(n)).astype(float)
        got.append(
            fit_gdna_strand_overdispersion(
                s, n.astype(float), rna_sense_frac=KAPPA
            ).gdna_strand_overdispersion
        )
    assert float(np.mean(got)) == pytest.approx(od_true, abs=0.006)


def test_the_between_seed_variance_IS_the_Beta_fourth_moment():
    """Pinned directly, because no behavioural gate can pin it. ``V∞(ρ) = 2ρ²(1−ρ)/(1+2ρ)`` is derived
    from the symmetric Beta's moments — ``Var(od̂_s | p) → Var(4u²) = 16(E[u⁴] − (E[u²])²)`` — and a
    wrong ``V∞`` costs only efficiency, never bias, so the estimator keeps working. PERTURBATION:
    dropping the ``(1−ρ)/(1+2ρ)`` factor fires NOTHING elsewhere, so the derivation is checked against
    the distribution itself, by Monte Carlo.
    """
    assert between_seed_variance(0.0) == 0.0
    for rho in (0.05, 0.2, 0.5):
        a = 0.5 * (1.0 - rho) / rho
        u = np.random.default_rng(int(1000 * rho)).beta(a, a, 2_000_000) - 0.5
        assert between_seed_variance(rho) == pytest.approx(float(np.var(4.0 * u * u)), rel=0.02)
    # and the shape the weights depend on: it must GROW with rho and saturate the deep-seed weight
    assert between_seed_variance(0.2) > between_seed_variance(0.05) > between_seed_variance(0.01)


def test_at_rho_zero_the_weighting_IS_the_pair_count_estimator():
    """``w_s = 1/(½ + c_s·V∞(0)) = 2`` — a constant that cancels. So the previous estimator is this one with
    ρ pinned at 0, and on a null population the fit and the ρ = 0 moment agree to the bisection's
    precision. A panel whose truth is ρ = 0 therefore cannot see this mechanism at all."""
    assert between_seed_variance(0.0) == 0.0
    c = np.array([0.5, 18.0, 39900.0])
    np.testing.assert_allclose(influence_weights(c, 0.0), 2.0)

    rng = np.random.default_rng(42)
    s, t = _seeds(rng, od_true=0.0, n_pure=3000)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert m.gdna_strand_overdispersion == pytest.approx(max(m.raw_overdispersion, 0.0), abs=1e-6)


def test_a_HANDFUL_OF_DEEP_SEEDS_CANNOT_DECIDE_THE_ANSWER():
    """The mechanism's own gate. 4,000 shallow seeds and 8 very deep ones at a true ρ = 0.2: the pair
    count hands each deep seed hundreds of times the weight of a shallow one, and the ρ = 0 moment is
    biased LOW because a deep seed's ``od̂_s`` is enormously variable. Weighting by ``1/Var(od̂_s|ρ)``
    caps a deep seed at ``1/V∞`` and recovers the truth, so the weighting is not only more robust but
    more ACCURATE, and that is the property to defend.

    The fixture is that lopsided because a real library can put most of its numerator on one seed
    (TRAPS: pair-count-weighting-lets-one-seed-decide).
    """
    raw, fitted = [], []
    for rep in range(8):
        rng = np.random.default_rng(950 + rep)
        n = np.concatenate([rng.integers(2, 20, 4000), np.full(8, 2000)])
        s = _bb(rng, n, 0.5, 0.2, len(n)).astype(float)
        m = fit_gdna_strand_overdispersion(s, n.astype(float), rna_sense_frac=KAPPA)
        raw.append(m.raw_overdispersion)
        fitted.append(m.gdna_strand_overdispersion)
    assert float(np.mean(raw)) < 0.17  # the pair-count moment is pulled off the truth
    assert float(np.mean(fitted)) == pytest.approx(0.2, abs=0.012)


def test_the_weighting_DE_CONCENTRATES_the_estimate():
    """One deep seed among 400 shallow ones: at ρ = 0 it IS the estimate (0.83, participation ~1); the fit
    demotes it, and both the value and the concentration move. Concentration is reported at the FIT's ρ, so
    ``effective_seeds`` describes the estimate that was actually made."""
    tot = np.concatenate([np.full(400, 9.0), [400.0]])
    sen = np.concatenate([np.full(400, 5.0), [400.0]])
    m = fit_gdna_strand_overdispersion(sen, tot, rna_sense_frac=KAPPA)
    assert m.raw_overdispersion > 0.8  # the deep seed alone
    assert m.gdna_strand_overdispersion < 0.05  # the population, once it cannot dominate
    assert m.effective_seeds > 3.0 * seed_participation(sen, tot, KAPPA)


def test_the_root_is_BRACKETED_so_bisection_always_terminates():
    """``g(ρ) = clip(moment(ρ)) − ρ`` has ``g(0) ≥ 0`` and ``g(ceiling) ≤ 0`` by construction, so a root
    always exists in ``[0, ceiling]`` and no iteration limit is asserted. The returned value is that root:
    re-evaluating the weighted moment at it reproduces it."""
    rng = np.random.default_rng(43)
    for od_true in (0.0, 0.03, 0.2, 0.6):
        s, t = _seeds(rng, od_true=od_true, n_pure=1500)
        m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
        od = m.gdna_strand_overdispersion
        assert 0.0 <= od <= _MAX_OVERDISPERSION
        g0 = np.clip(away_half_moment(s, t, KAPPA, overdispersion=0.0)[0], 0.0, _MAX_OVERDISPERSION)
        gc = np.clip(
            away_half_moment(s, t, KAPPA, overdispersion=_MAX_OVERDISPERSION)[0],
            0.0,
            _MAX_OVERDISPERSION,
        )
        assert g0 - 0.0 >= 0.0 and gc - _MAX_OVERDISPERSION <= 0.0  # the bracket
        if 0.0 < od < _MAX_OVERDISPERSION:  # an interior root reproduces itself
            back = np.clip(
                away_half_moment(s, t, KAPPA, overdispersion=od)[0], 0.0, _MAX_OVERDISPERSION
            )
            assert back == pytest.approx(od, abs=1e-9)


# ── the two components reconcile against EACH OTHER, not against a constant ─────────────────────


def test_the_WEAKER_component_shrinks_toward_the_BETTER_MEASURED_one():
    """The reference each component is pulled toward is a MEASUREMENT of the same library, and the
    weight is that measurement's own information, so there is no conjured target and no conjured
    weight anywhere in the reconciliation. The better-informed component does not move."""
    r, g = reconcile_overdispersions(0.010, 1e6, 0.130, 1e3)  # RNA far better measured
    assert r == pytest.approx(0.010)  # the strong one is untouched
    assert 0.010 < g < 0.011  # the weak one is pulled almost all the way to it
    r2, g2 = reconcile_overdispersions(0.010, 1e3, 0.130, 1e6)  # roles swapped
    assert g2 == pytest.approx(0.130)
    assert 0.129 < r2 < 0.130


def test_EQUALLY_MEASURED_components_do_not_move_at_all():
    """PERTURBATION: with the naive ``(I_w·od_w + I_s·od_s)/(I_w + I_s)`` blend, two EQUALLY
    well-measured components still have the weaker one dragged to their midpoint while the stronger
    does not move — asymmetric, and a claim neither measurement supports. Borrowing the DEFICIT
    ``(I_s − I_w)/I_s`` makes the borrow weight exactly 0 here.

    The difference is not academic: a real library's two components can differ by more than an order
    of magnitude with both well measured, and pooling them erases that."""
    r, g = reconcile_overdispersions(0.008, 1e6, 0.133, 1e6)
    assert r == pytest.approx(0.008) and g == pytest.approx(0.133)

    # and it is continuous in between: half the information ⇒ half the borrow
    _r, g_half = reconcile_overdispersions(0.008, 1e6, 0.133, 5e5)
    assert g_half == pytest.approx(0.5 * 0.133 + 0.5 * 0.008)


def test_a_component_with_no_evidence_takes_the_others_MEASURED_value():
    """Zero information ⇒ borrow outright. That is the case the deleted constant used to fill, and it is
    now filled by a number measured on the same library."""
    r, g = reconcile_overdispersions(float("nan"), 0.0, 0.075, 1e5)
    assert r == pytest.approx(0.075) and g == pytest.approx(0.075)
    r, g = reconcile_overdispersions(0.021, 1e5, float("nan"), 0.0)
    assert r == pytest.approx(0.021) and g == pytest.approx(0.021)


def test_with_NEITHER_measured_both_take_the_CEILING_and_coincide():
    """Any common value leaves the strand channel uninformative, because the composition term reads the
    DIFFERENCE between the two dispersions. The ceiling is the one constant already asserted, so this
    introduces none — and it errs toward 'the channel says nothing' rather than 'the channel is certain'."""
    r, g = reconcile_overdispersions(float("nan"), 0.0, float("nan"), 0.0)
    assert r == g == pytest.approx(_MAX_OVERDISPERSION)


def test_the_reconciled_pair_stays_inside_the_physical_support():
    """Both outputs are clipped, including when an input is out of range."""
    r, g = reconcile_overdispersions(0.9, 1e6, -0.3, 1e6)
    assert 0.0 <= r <= _MAX_OVERDISPERSION and 0.0 <= g <= _MAX_OVERDISPERSION
    assert r == pytest.approx(_MAX_OVERDISPERSION) and g == pytest.approx(0.0)


# ── the wrapper, the boundary seed selectors and the strand likelihood that reads the fit ─────


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
    # Either the away half saw enough to answer (and it cannot be INFLATED), or it saw nothing at all and
    # says so — the ceiling fallback is "I measured nothing", never a claim about this population.
    assert model.fallback_used or model.gdna_strand_overdispersion <= true_od + 0.02


def test_thin_seed_fallback():
    """No gDNA strand signal (empty, or every seed on the RNA side) → the CEILING, no crash.

    No conjured constant survives anywhere in this fit. Having measured nothing, the estimator
    returns the widest dispersion the model admits, so the strand channel says
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


def _mock_substrate(pos, neg, ts, count_evidence, observable):
    """Mock with only contained-REGION signal: the boundary axis carries ZERO counts ⇒ no boundary seeds.

    One reference with ``n`` regions owns exactly ``n − 1`` boundaries, so the boundary arrays are
    sized from the region axis rather than left empty — a boundary axis inconsistent with its own
    ``ref_id`` is not a "no boundaries" fixture, it is a mis-shaped one.
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
    """A substrate plus the two count-observability MASKS over ``signatures``, on the region/boundary axes."""
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
    """ONE seed per boundary, not two.

    Emitting region ``r``'s right side and region ``r+1``'s left side separately would put the same
    physical crossing into one pooled moment estimator twice, inflating its apparent sample size by 2×
    and correlating every pair perfectly. A contiguous boundary is a 0-bp boundary with ONE count, so
    there is ONE seed.
    """
    from rigel.calibration.gdna_strand import boundary_seeds

    # exon+ | intron+ : one boundary, count-observable (no shared EXON bit), oriented POS.
    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_EXON_POS, BIT_INTRON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, total = boundary_seeds(
        substrate, region_arrays, region_density.boundary_count_observable
    )
    assert sense.shape == (1,)  # (2,) would be the same crossing, twice
    np.testing.assert_allclose(sense, [70.0])
    np.testing.assert_allclose(total, [100.0])


def test_boundary_seed_sense_follows_the_flanking_transcript_strand():
    """A NEG-strand boundary orients to the NEG genome column; the sense count is not always ``pos``."""
    from rigel.calibration.gdna_strand import boundary_seeds

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
    from rigel.calibration.gdna_strand import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [0, BIT_INTRON_NEG], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    np.testing.assert_allclose(sense, [30.0])  # oriented NEG by the gene side


def test_an_opposite_strand_boundary_is_not_strand_observable():
    """``{POS, NEG}`` leaves 'sense' undefined, so the boundary cannot seed the fit at all."""
    from rigel.calibration.gdna_strand import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_INTRON_POS, BIT_INTRON_NEG], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


def test_an_AMBIG_flank_cannot_seed():
    """PERTURBATION: nothing else pins this rule.

    A flank carrying overlapping ± transcripts has no defined transcript sense, so neither genome
    column is "sense" and the boundary cannot seed a strand fit. ``boundary_strand_orientation`` used to carry
    an explicit ``~either_ambig`` guard — deleting it broke NOTHING, because ``TS_AMBIG`` is a fourth
    distinct value that the ``POS-or-NONE`` / ``NEG-or-NONE`` tests already exclude. The guard was
    dead code claiming to be the rule; this test is the rule.

    The fixture uses INTRON bits on both strands, not exon bits: an AMBIG flank must still be
    count-observable, or the test would pass for the wrong reason.
    """
    from rigel.calibration.gdna_strand import boundary_seeds

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
    from rigel.calibration.gdna_strand import boundary_seeds

    substrate, region_arrays, region_density = _boundary_parts(
        [BIT_EXON_POS, BIT_EXON_POS], boundary_pos=[70.0], boundary_neg=[30.0]
    )
    sense, _ = boundary_seeds(substrate, region_arrays, region_density.boundary_count_observable)
    assert sense.shape == (0,)


def test_boundary_seeds_never_straddle_a_reference():
    """Two single-region references own ZERO boundaries between them — nothing can leak across."""
    from rigel.calibration.gdna_strand import boundary_seeds

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
    """A noise-skewed pure-gDNA split reads as MORE gDNA at od>0.

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
    """The gDNA fit carries NO location prior. Deep seeds follow the fit and sparse seeds return the
    raw (noisy) moment inside the physical support; neither is pulled toward a constant, because every
    constant on offer is either conjured or over-shrinks exactly where gDNA information is scarce."""
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


# ── end to end: `calibrate` recovers a planted gDNA strand overdispersion ──────────────────────


_KAPPA = 0.95


def _consistent_strand_model(overdispersion, n_sj=800, depth=60, seed=7):
    """A strand model whose SPLICED seeds carry the SAME overdispersion as the planted gDNA one.

    A library has ONE technical process, and calibrate reconciles the two components against each
    other (`gdna_strand.reconcile_overdispersions`), so a fixture that plants `od_g = 0.2` beside a
    single junction sitting exactly at κ — which reads `od_r = 0` with real information — asserts a
    library that cannot exist, and the reconciliation correctly splits the difference. Planting the
    same dispersion in both components is what keeps this an end-to-end RECOVERY test.
    """
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels
    from rigel.types import Strand

    rng = np.random.default_rng(seed)
    if overdispersion > 0:
        conc = (1.0 - overdispersion) / overdispersion
        p = rng.beta(_KAPPA * conc, (1.0 - _KAPPA) * conc, size=n_sj)
    else:
        p = np.full(n_sj, _KAPPA)
    n_sense = rng.binomial(depth, p).astype(np.int64)
    table = SJStrandTable(
        ref_id=np.zeros(n_sj, dtype=np.int32),
        start=np.arange(n_sj, dtype=np.int64) * 1000,
        end=np.arange(n_sj, dtype=np.int64) * 1000 + 100,
        motif_strand=np.full(n_sj, int(Strand.POS), dtype=np.int8),
        n_sense=n_sense,
        n_antisense=np.full(n_sj, depth, dtype=np.int64) - n_sense,
    )
    return StrandModels(exonic_spliced=StrandModel.from_sj_table(table))


_STRAND_MODEL = make_strand_models(_KAPPA, 200)
_FRAG_LEN = 50  # the delta the gDNA pmf sits at; every fixture fragment is this long


def _intron_betabinom_payload(n_regions, depth, overdispersion, seed):
    """A 1-reference payload of ``n_regions`` intron(+) regions; contained gDNA ~ BetaBinom(½, od).

    The boundary axis is empty of COUNTS but not of ROWS. One reference with ``k`` regions owns
    ``k − 1`` boundaries, and the payload must carry them or the chain builder refuses it. Leaving the
    counts at zero is what makes this test isolate the CONTAINED-region seed arm: with no crossing
    fragments there are no boundary seeds, so a recovered overdispersion can only have come from the regions.
    """
    rng = np.random.default_rng(seed)
    a = 0.5 * (1.0 - overdispersion) / overdispersion if overdispersion > 0 else 1e9
    p = rng.beta(a, a, size=n_regions)
    pos = rng.binomial(depth, p)
    neg = depth - pos

    contained = np.stack([pos, neg], axis=1).astype(np.uint32)
    n_boundaries = n_regions - 1
    quantum = 1.0 / _FRAG_LEN

    def region_zeros(dtype):
        return np.zeros((n_regions, 2), dtype=dtype)

    def boundary_zeros(dtype):
        return np.zeros((n_boundaries, 2), dtype=dtype)

    def flat(rows, dtype):
        """A single-column bank — the length moments and the conserved mass carry no strand axis."""
        return np.zeros(rows, dtype=dtype)

    payload = AccumulatorPayload(
        region_bounds=np.arange(n_regions + 1, dtype=np.int64) * 100,
        ref_region_bound_offsets=np.array([0, n_regions + 1], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_boundaries], dtype=np.int64),
        ref_sj_offsets=np.array([0, 0], dtype=np.int64),
        region_contained_count=contained,
        # ONE column: the length moments carry no strand axis, so the two are summed.
        region_contained_inv_opportunity_sum=(
            contained.sum(axis=1).astype(np.uint64) * np.uint64(quantum)
        ),
        region_start_count=contained.astype(np.uint32),
        region_end_count=contained.astype(np.uint32),
        region_span_count=region_zeros(np.uint32),
        boundary_unspliced_count=boundary_zeros(np.uint32),
        boundary_unspliced_inv_length_sum=flat(n_boundaries, np.float64),
        # ONE value per boundary — the conserved mass has no strand axis. Zero here is a real state and
        # not a stub: this fixture deposits no crossings at all, so there is no mass to conserve.
        boundary_unspliced_mass=np.zeros(n_boundaries, dtype=np.float64),
        boundary_spliced_count=boundary_zeros(np.uint32),
        boundary_spliced_mass=np.zeros(n_boundaries, dtype=np.float64),
        sj_count=np.zeros((0, 2), dtype=np.uint32),
        sj_inv_length_sum=np.zeros(0, dtype=np.float64),
        sj_mass=np.zeros(0, dtype=np.float64),
        pool_lengths=np.zeros((N_FRAGMENT_POOLS, 201), dtype=np.int64),
        deposited_lengths=np.zeros(201, dtype=np.uint32),
        # Nothing was deferred here, and that is a real state, not a stub: this fixture has no
        # annotated intron in any mate gap. The two empty spellings live on the classes so a
        # hand-built payload cannot get the `[0]`-not-`[]` offset boundary wrong.
        deferred=DeferredFragments.empty(),
        gap_resolution=GapCensus.zeros(),
        qc=ScanQC(
            deposited=int(contained.sum()),
            dropped_too_long=0,
            dropped_empty=0,
            dropped_strand_undefined=0,
            deferred_undetermined_gap=0,
            unannotated_introns=0,
            contradictory_sj_strand=0,
            introns_absorbed=0,
        ),
        max_length=200,
        n_refs=1,
    )
    starts = np.arange(n_regions, dtype=np.int64) * 100
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n_regions, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * n_regions, dtype="string"),
            "start": starts,
            "end": starts + 100,
            "length": np.full(n_regions, 100, dtype=np.int64),
            "signature": np.full(
                n_regions, BIT_INTRON_POS, dtype=np.uint8
            ),  # intron ⇒ count-observable, genic
        }
    )
    return payload, RegionArrays.from_frame(region_df, {"chr1": 0})


def _calibrate(payload, ra, strand_model=None):
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model if strand_model is not None else _STRAND_MODEL,
        gdna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        rna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        config=CalibrationConfig(),
    )


@pytest.mark.parametrize("od_true", [0.05, 0.10, 0.20])
def test_calibrate_recovers_overdispersion(od_true):
    payload, ra = _intron_betabinom_payload(
        n_regions=400, depth=150, overdispersion=od_true, seed=int(od_true * 1000)
    )
    # ONE library, ONE technical process: the spliced seeds carry the same dispersion as the genomic
    # ones, so the two components agree and the reconciliation is a no-op — which is what makes this an
    # end-to-end RECOVERY test rather than a test of how far the components get pulled toward each other.
    result = _calibrate(payload, ra, _consistent_strand_model(od_true))
    assert result.gdna_strand_overdispersion == pytest.approx(od_true, rel=0.25, abs=0.02)
    assert result.rna_strand_overdispersion == pytest.approx(od_true, rel=0.30, abs=0.02)


def test_calibrate_binomial_gdna_floors_to_zero():
    """Non-overdispersed (50/50) gDNA → the identifiability gate floors od to 0 (Binomial)."""
    payload, ra = _intron_betabinom_payload(n_regions=400, depth=150, overdispersion=0.0, seed=7)
    assert _calibrate(payload, ra).gdna_strand_overdispersion < 0.02


def test_a_sj_free_library_calibrates(od_true=0.10):
    """``n_sj == 0`` is legal and must not be confused with "no sj flux": this payload's
    references are all single-region-signature intergenic, so the graph has no sj boundary at all.
    ``calibrate`` defaults to an empty sj axis and the result carries a length-0 array."""
    payload, ra = _intron_betabinom_payload(
        n_regions=400, depth=150, overdispersion=od_true, seed=11
    )
    result = _calibrate(payload, ra)
    assert result.n_sj == 0
    assert result.count_rna_sj.shape == (0,)
    assert result.n_boundaries == result.n_regions - 1
