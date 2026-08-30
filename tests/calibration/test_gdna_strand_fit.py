"""FALSIFICATION GATES for the gDNA strand-overdispersion fit — the AWAY-HALF estimator (landed 2026-08-29).

The premise every earlier fit rested on — that some structural class of objects is PURE gDNA — is a property
of the annotation and of the sample, not of the genome: pervasive transcription is real, the intergenic space
is whatever the user's GTF leaves over, and most genes are OFF in any one sample but nobody knows which. So
the estimator trusts no class.

THE LEMMA. On a strand-specific library RNA of a seed's own gene can only pull its sense fraction TOWARD κ.
Orient ``d = K − N/2`` so RNA pulls ``d`` negative. Under pure gDNA ``d`` is symmetric about 0 and the moment
excess ``(d² − N/4)`` is an EVEN function of ``d``, so over the AWAY half (``d > 0``, ties at weight ½) the
pooled method-of-moments ratio is unbiased for ρ_g for EVERY distribution of the seeds' RNA content — a
contaminated seed reaches the away side only by noise, with a small ``d``, biasing ρ̂ DOWN, never up.
The seed set is every count- and strand-observable GENIC object (intron regions; exon|intron AND gene-edge
boundaries) in whatever GTF is supplied. INTERGENIC seeds have no gene strand to orient by and are OUT; AMBIG
(annotated sense + antisense) objects have no defined sense and are OUT. Known limit, recorded not hidden:
unannotated ANTISENSE transcription pushes toward the away side and inflates ρ̂.

Each property below was watched FAILING against the previous fit before this one landed, and each carries its
perturbation: force the full (two-sided) moment back and the contamination gate fires.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from rigel.calibration.gdna_strand import (
    reconcile_overdispersions,
    between_seed_variance,
    influence_weights,
    seed_participation,
    _MAX_OVERDISPERSION,
    _null_information,
    away_half_moment,
    fit_gdna_strand_from_substrate,
    fit_gdna_strand_overdispersion,
)
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)
from rigel.calibration.strand_deconv import boundary_seeds

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
    """⛔ THE DOUBLE-HALVING. ``information`` must be ½ × the pair count of ALL seeds, because the away
    half's own pair count is ALREADY about half of that — halving it again understates the information 2×
    and overstates the standard error by √2.

    The property, not the implementation: ``Var(od_mom)|₀ = 2/P`` for the total pair count ``P``. Derivation
    — ``Var(e_s)|₀ = n(n−1)/8`` and ``E[a_s] = ½`` (the sign of a symmetric residual is independent of its
    size), so ``Var(num) = P/8``, ``E[den] = P/4`` and ``Var(num/den) = 2/P``. Asserted here BOTH ways: the
    closed form, and a Monte-Carlo null whose empirical sd must match ``1/sqrt(information)``.

    ⚠ The gate this replaces asserted ``info == 0.5·_null_information(total[sense >= 1])`` — the
    implementation's own expression restated — so it certified the defect and fired on the repair."""
    rng = np.random.default_rng(21)
    total = rng.integers(2, 40, 3000).astype(float)
    sense = rng.binomial(total.astype(int), 0.5).astype(float)
    _od, info = away_half_moment(sense, total, KAPPA)
    pairs = float((total * (total - 1.0) / 2.0).sum())
    assert info == pytest.approx(0.5 * pairs, rel=1e-12)  # ⛔ NOT the away half's own pair count

    # the Monte-Carlo null: 1/sqrt(info) IS the standard error it claims to be
    n = np.full(2000, 5.0)
    draws = rng.binomial(5, 0.5, size=(500, 2000)).astype(float)
    ods = np.array([away_half_moment(draws[i], n, KAPPA)[0] for i in range(draws.shape[0])])
    _od0, info0 = away_half_moment(draws[0], n, KAPPA)
    assert ods.std() == pytest.approx(1.0 / np.sqrt(info0), rel=0.15)


def test_an_unstranded_library_uses_the_FULL_moment_and_is_not_forced_to_zero():
    """⛔ κ = ½ EXACTLY IS REACHABLE and must not collapse the fit. ``rna_sense_frac`` is the posterior mean
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
    assert m.gdna_strand_overdispersion == pytest.approx(0.05, abs=0.01)  # ⛔ was a hard 0.0

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
    """⛔ No conjured constant anywhere. With no pair the fit returns the CEILING — the widest strand
    likelihood the model admits, so the channel says nothing — and ZERO information, the signal
    `reconcile_overdispersions` uses to hand it the other component's measured value instead."""
    m = fit_gdna_strand_overdispersion(np.array([]), np.array([]), rna_sense_frac=KAPPA)
    assert m.fallback_used and m.gdna_strand_overdispersion == pytest.approx(_MAX_OVERDISPERSION)
    assert (
        m.information == 0.0
    )  # zero information is what lets the OTHER component supply the value


# ── the seed SELECTORS: genic only; intergenic and AMBIG out ───────────────────────────────────


def _view(pos, neg):
    return SimpleNamespace(count=np.stack([np.asarray(pos, float), np.asarray(neg, float)], axis=1))


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
    """⛔ LOAD-BEARING AND PREVIOUSLY UNPINNED (found by adversarial review, 2026-08-29): the whole
    calibration suite passed with ``_region_seeds``' ``np.where(ts == TS_NEG, neg, pos)`` flip deleted.

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
    np.testing.assert_allclose(sense, neg)  # ⛔ the NEG column is sense here, not `pos`
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
    """⛔ A value at the ceiling is a CLAMP, not a measurement — measured on 2 of 4 real cfRNA libraries
    (raw 0.675 and 0.974 against a ceiling of 0.2). The model must carry both facts, or a reader takes the
    ceiling for a fit."""
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
    """⛔ TWO HOLES FOUND BY PERTURBATION, NOT BY REVIEW — the preceding gates tested the function and the
    model separately and left both of these open:

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
    assert m.effective_seeds < 0.05 * m.n_seed_regions  # ⛔ the count would say 401

    rng = np.random.default_rng(32)
    s, t = _seeds(rng, od_true=0.0, n_pure=4000)
    assert seed_participation(s, t, KAPPA) > 400.0  # ⛔ a signed sum collapses this to ~0.55

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
    """⛔ PINNED DIRECTLY, because no behavioural gate can pin it. ``V∞(ρ) = 2ρ²(1−ρ)/(1+2ρ)`` is derived
    from the symmetric Beta's moments — ``Var(od̂_s | p) → Var(4u²) = 16(E[u⁴] − (E[u²])²)`` — and a WRONG
    ``V∞`` costs only efficiency, never bias, so the estimator keeps working and every other gate stays
    green (measured: dropping the ``(1−ρ)/(1+2ρ)`` factor fired NOTHING). The derivation is therefore
    checked against the distribution itself, by Monte Carlo.
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
    ρ pinned at 0, and on a null population the fit and the ρ = 0 moment agree to the bisection's precision.
    ⭐ This is why neither simulated panel moved: their truth is ρ = 0."""
    assert between_seed_variance(0.0) == 0.0
    c = np.array([0.5, 18.0, 39900.0])
    np.testing.assert_allclose(influence_weights(c, 0.0), 2.0)

    rng = np.random.default_rng(42)
    s, t = _seeds(rng, od_true=0.0, n_pure=3000)
    m = fit_gdna_strand_overdispersion(s, t, rna_sense_frac=KAPPA)
    assert m.gdna_strand_overdispersion == pytest.approx(max(m.raw_overdispersion, 0.0), abs=1e-6)


def test_a_HANDFUL_OF_DEEP_SEEDS_CANNOT_DECIDE_THE_ANSWER():
    """⭐ THE MECHANISM'S OWN GATE. 4,000 shallow seeds and 8 very deep ones at a true ρ = 0.2: the pair
    count gives the 8 deep seeds ~500× the weight each, and the ρ = 0 moment is biased LOW (~0.147) because
    a deep seed's ``od̂_s`` is enormously variable. Weighting by ``1/Var(od̂_s|ρ)`` caps a deep seed at
    ``1/V∞`` and recovers the truth. ⛔ The weighting is therefore not only more robust but more ACCURATE,
    and that is the property to defend.

    The real-data analogue this stands in for: one seed carrying 77.8 % of a library's numerator.
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
    """⭐ THE REPLACEMENT FOR `Beta(14,14)` (owner, 2026-08-30). The reference is a MEASUREMENT of the same
    library and the weight is that measurement's own information, so both conjured constants — the target
    0.0345 and its derived weight ~909 — are gone. The better-informed component does not move."""
    r, g = reconcile_overdispersions(0.010, 1e6, 0.130, 1e3)  # RNA far better measured
    assert r == pytest.approx(0.010)  # the strong one is untouched
    assert 0.010 < g < 0.011  # the weak one is pulled almost all the way to it
    r2, g2 = reconcile_overdispersions(0.010, 1e3, 0.130, 1e6)  # roles swapped
    assert g2 == pytest.approx(0.130)
    assert 0.129 < r2 < 0.130


def test_EQUALLY_MEASURED_components_do_not_move_at_all():
    """⛔ Found by perturbing my own gate, not by review: with the naive
    ``(I_w·od_w + I_s·od_s)/(I_w + I_s)`` blend, two EQUALLY well-measured components still had the weaker
    one dragged to their midpoint while the stronger did not move — asymmetric, and a claim neither
    measurement supports. Borrowing the DEFICIT ``(I_s − I_w)/I_s`` makes the borrow weight exactly 0 here.

    The difference is not academic: on one cfRNA library RNA measures 0.008 against gDNA 0.133, and pooling
    them would erase a real, well-measured difference."""
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
    """⭐ Any common value leaves the strand channel uninformative, because the composition term reads the
    DIFFERENCE between the two dispersions. The ceiling is the one constant already asserted, so this
    introduces none — and it errs toward 'the channel says nothing' rather than 'the channel is certain'."""
    r, g = reconcile_overdispersions(float("nan"), 0.0, float("nan"), 0.0)
    assert r == g == pytest.approx(_MAX_OVERDISPERSION)


def test_the_reconciled_pair_stays_inside_the_physical_support():
    """Both outputs are clipped, including when an input is out of range."""
    r, g = reconcile_overdispersions(0.9, 1e6, -0.3, 1e6)
    assert 0.0 <= r <= _MAX_OVERDISPERSION and 0.0 <= g <= _MAX_OVERDISPERSION
    assert r == pytest.approx(_MAX_OVERDISPERSION) and g == pytest.approx(0.0)
