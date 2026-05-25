"""Tests for the v4 gDNA density model."""

from __future__ import annotations

import numpy as np
import pytest
from scipy.stats import nbinom

from rigel.calibration.density_model import (
    DENSITY_FALLBACK_PRIOR_ALPHA,
    DENSITY_MIN_EFF_LENGTH,
    DENSITY_PHI_FLOOR,
    DENSITY_PRIOR_MIN_CV,
    FLAG_FALLBACK_USED,
    FLAG_HIGH_TAIL_TENSION,
    FLAG_LOW_BOUNDARY_OPPORTUNITY,
    FLAG_LOW_CONTAINED_OPPORTUNITY,
    FLAG_NON_ANCHOR,
    FLAG_PRIOR_DOMINATED,
    PRIOR_FAMILY_ALL,
    PRIOR_FAMILY_DETERMINISTIC_ZERO,
    PRIOR_FAMILY_FALLBACK_BROAD,
    compute_beta_cap,
    fit_density_evidence,
    fit_gamma_prior,
    select_rho_ref,
)
from rigel.calibration.density_model import (  # private math helpers are tested as sentinels
    _broad_fallback_prior,
    _compute_predictive_diagnostics,
    _assign_region_priors,
)
from rigel.calibration.density_observation import DensityObservation


def _obs(
    contained: np.ndarray,
    contained_leff: np.ndarray,
    *,
    boundary: np.ndarray | None = None,
    boundary_leff: np.ndarray | None = None,
    intergenic: np.ndarray | None = None,
    intron: np.ndarray | None = None,
) -> DensityObservation:
    C = np.asarray(contained, dtype=np.float32)
    Lc = np.asarray(contained_leff, dtype=np.float64)
    B = np.zeros_like(C, dtype=np.float32) if boundary is None else np.asarray(boundary, dtype=np.float32)
    Lb = (
        np.zeros_like(Lc, dtype=np.float64)
        if boundary_leff is None
        else np.asarray(boundary_leff, dtype=np.float64)
    )
    if intergenic is None:
        intergenic = np.ones(C.size, dtype=bool)
    if intron is None:
        intron = np.zeros(C.size, dtype=bool)
    inter = np.asarray(intergenic, dtype=bool)
    intr = np.asarray(intron, dtype=bool)
    anchor = inter | intr
    return DensityObservation(
        contained_count=C,
        boundary_left_count=(B / 2.0).astype(np.float32),
        boundary_right_count=(B / 2.0).astype(np.float32),
        boundary_count=B,
        observed_compatible_count=(C + B).astype(np.float32),
        contained_leff=Lc,
        boundary_left_leff=Lb / 2.0,
        boundary_right_leff=Lb / 2.0,
        boundary_leff=Lb,
        anchor_intergenic=inter,
        anchor_intron=intr,
        is_anchor=anchor,
        spliced_count=np.zeros_like(C, dtype=np.float32),
        region_length=Lc.copy(),
    )


def _prior(name: str, mean: float, *, alpha: float = 25.0, status: str = "ok"):
    beta = alpha / mean if mean > 0.0 else 1.0
    return fit_gamma_prior(
        np.full(25, mean * 1000.0, dtype=np.float64),
        np.full(25, 1000.0, dtype=np.float64),
        family=name,
        beta_cap=float("inf"),
        min_regions=20,
    ).__class__(
        family=name,
        alpha=alpha,
        beta=beta,
        mean_density=mean,
        phi=1.0 / alpha,
        beta_raw=beta,
        beta_cap=float("inf"),
        cap_applied=False,
        n_regions=25,
        n_fragments=mean * 25_000.0,
        eff_length=25_000.0,
        residual_sum=0.0,
        poisson_variance_sum=0.0,
        extra_variance_basis_sum=0.0,
        trim_upper=0.0,
        n_trimmed=0,
        trimmed_mu_fraction=0.0,
        fallback_depth=0,
        fit_status=status,
        pearson_chi2_trimmed=0.0,
    )


def test_fit_gamma_prior_recovers_mean_density_on_poisson_anchors() -> None:
    rng = np.random.default_rng(1)
    rho = 0.02
    L = np.full(500, 1000.0, dtype=np.float64)
    C = rng.poisson(rho * L).astype(np.float64)

    prior = fit_gamma_prior(C, L, family="INTERGENIC", min_regions=20)

    assert prior.fit_status == "ok"
    assert prior.mean_density == pytest.approx(rho, rel=0.08)
    assert prior.n_regions == 500


def test_fit_gamma_prior_reports_finite_phi_on_overdispersed_anchors() -> None:
    rng = np.random.default_rng(2)
    rho = 0.02
    L = np.full(800, 1000.0, dtype=np.float64)
    rates = rng.gamma(shape=5.0, scale=rho / 5.0, size=L.size)
    C = rng.poisson(rates * L).astype(np.float64)

    prior = fit_gamma_prior(C, L, family="INTRON", beta_cap=float("inf"), min_regions=20)

    assert np.isfinite(prior.phi)
    assert prior.phi > DENSITY_PHI_FLOOR


def test_zero_phi_uses_floor_instead_of_infinite_alpha() -> None:
    L = np.full(50, 1000.0, dtype=np.float64)
    C = np.full(50, 10.0, dtype=np.float64)

    prior = fit_gamma_prior(
        C,
        L,
        family="INTERGENIC",
        beta_cap=float("inf"),
        min_regions=20,
    )

    assert prior.phi == pytest.approx(DENSITY_PHI_FLOOR)
    assert np.isfinite(prior.alpha)
    assert np.isfinite(prior.beta_raw)


def test_pearson_trim_does_not_preferentially_discard_high_mu_rows() -> None:
    rng = np.random.default_rng(3)
    L = np.exp(rng.uniform(np.log(10.0), np.log(10_000.0), size=6000))
    rho = 0.05
    C = rng.poisson(rho * L).astype(np.float64)

    prior = fit_gamma_prior(C, L, family="ALL", beta_cap=float("inf"), min_regions=20)

    mu = (C.sum() / L.sum()) * L
    raw_score = (C - mu) ** 2
    k_keep = int(np.ceil((1.0 - 0.05) * raw_score.size))
    raw_keep = np.argpartition(raw_score, k_keep - 1)[:k_keep]
    kept = np.zeros(raw_score.size, dtype=bool)
    kept[raw_keep] = True
    raw_trimmed_mu_fraction = float(mu[~kept].sum() / mu.sum())

    assert abs(prior.trimmed_mu_fraction - 0.05) < abs(raw_trimmed_mu_fraction - 0.05)
    assert raw_trimmed_mu_fraction > 0.10


def test_pearson_trim_resists_small_nrna_contaminated_anchor_fraction() -> None:
    rng = np.random.default_rng(4)
    L = np.full(1000, 1000.0, dtype=np.float64)
    rho = 0.02
    C = rng.poisson(rho * L).astype(np.float64)
    contaminated = C.copy()
    outlier_idx = rng.choice(C.size, size=20, replace=False)
    contaminated[outlier_idx] = 50.0 * rho * L[outlier_idx]

    robust = fit_gamma_prior(
        contaminated,
        L,
        family="INTRON",
        beta_cap=float("inf"),
        trim_upper=0.05,
        min_regions=20,
    )
    untrimmed = fit_gamma_prior(
        contaminated,
        L,
        family="INTRON",
        beta_cap=float("inf"),
        trim_upper=0.0,
        min_regions=20,
    )

    assert robust.beta_raw > 10.0 * untrimmed.beta_raw
    assert robust.pearson_chi2_trimmed < untrimmed.pearson_chi2_trimmed


def test_compute_beta_cap_and_depth_invariance() -> None:
    assert compute_beta_cap(0.01) == pytest.approx(40_000.0)
    assert np.isinf(compute_beta_cap(0.0))

    C = np.full(50, 10.0, dtype=np.float64)
    L = np.full(50, 1000.0, dtype=np.float64)
    p1 = fit_gamma_prior(C, L, family="ALL", min_regions=20)
    p2 = fit_gamma_prior(C * 10.0, L * 10.0, family="ALL", min_regions=20)

    assert p1.mean_density == pytest.approx(p2.mean_density)
    assert p1.beta_cap == pytest.approx(p2.beta_cap)


def test_beta_cap_preserves_mean_and_enforces_prior_cv_floor() -> None:
    C = np.full(50, 10.0, dtype=np.float64)
    L = np.full(50, 100.0, dtype=np.float64)

    prior = fit_gamma_prior(C, L, family="ALL", beta_cap=10.0, min_regions=20)

    assert prior.cap_applied
    assert prior.alpha / prior.beta == pytest.approx(prior.mean_density)
    assert 1.0 / np.sqrt(prior.alpha) >= DENSITY_PRIOR_MIN_CV


def test_select_rho_ref_priority_order() -> None:
    inter = _prior("INTERGENIC", 0.01, alpha=10.0)
    intron = _prior("INTRON", 0.03, alpha=30.0)
    all_prior = _prior("ALL", 0.02, alpha=40.0)

    assert select_rho_ref({"ALL": all_prior, "INTERGENIC": inter}) == pytest.approx(
        (0.02, "ALL")
    )

    rho, source = select_rho_ref(
        {"ALL": _prior("ALL", 0.02, status="failed"), "INTERGENIC": inter, "INTRON": intron}
    )
    assert source == "WEIGHTED_FAMILIES"
    assert rho == pytest.approx((10.0 * 0.01 + 30.0 * 0.03) / 40.0)

    assert select_rho_ref({"INTRON": intron}) == pytest.approx((0.03, "WEIGHTED_FAMILIES"))
    assert select_rho_ref({}) == (0.0, "ZERO")


def test_depth2_fallback_prior_shape_is_weak_gamma() -> None:
    rho_ref = 0.0125
    prior = _broad_fallback_prior(rho_ref)

    assert prior.alpha == pytest.approx(DENSITY_FALLBACK_PRIOR_ALPHA)
    assert prior.beta == pytest.approx(prior.alpha / rho_ref)
    assert prior.mean_density == pytest.approx(rho_ref)
    assert 1.0 / np.sqrt(prior.alpha) == pytest.approx(1.0)


def test_deterministic_zero_density_branch() -> None:
    observation = _obs(
        np.zeros(25, dtype=np.float32),
        np.full(25, 1000.0, dtype=np.float64),
    )

    evidence = fit_density_evidence(observation)

    assert np.all(evidence.relative_exposure == 1.0)
    assert np.all(evidence.rho_post == 0.0)
    assert np.all(evidence.prior_family == PRIOR_FAMILY_DETERMINISTIC_ZERO)
    assert np.all((evidence.flags & (FLAG_FALLBACK_USED | FLAG_PRIOR_DOMINATED)) != 0)
    assert evidence.to_summary_dict()["rho_ref_source"] == "ZERO"


def test_local_update_zero_boundary_leaves_posterior_equal_prior() -> None:
    C = np.concatenate([np.full(25, 10.0), np.array([0.0])])
    Lc = np.concatenate([np.full(25, 1000.0), np.array([500.0])])
    boundary = np.zeros_like(C)
    Lb = np.zeros_like(C)
    inter = np.concatenate([np.ones(25, dtype=bool), np.array([False])])
    observation = _obs(C, Lc, boundary=boundary, boundary_leff=Lb, intergenic=inter)

    evidence = fit_density_evidence(observation)

    assert evidence.rho_post[-1] == pytest.approx(evidence.rho_ref)
    assert evidence.relative_exposure[-1] == pytest.approx(1.0)


def test_high_boundary_count_increases_local_density() -> None:
    C = np.concatenate([np.full(25, 10.0), np.array([0.0])])
    Lc = np.concatenate([np.full(25, 1000.0), np.array([500.0])])
    boundary = np.concatenate([np.zeros(25), np.array([100.0])])
    Lb = np.concatenate([np.zeros(25), np.array([1000.0])])
    inter = np.concatenate([np.ones(25, dtype=bool), np.array([False])])
    observation = _obs(C, Lc, boundary=boundary, boundary_leff=Lb, intergenic=inter)

    evidence = fit_density_evidence(observation)

    assert evidence.rho_post[-1] > evidence.rho_ref
    assert evidence.relative_exposure[-1] > 1.0


def test_nb_predictive_mean_and_ppf_match_scipy() -> None:
    diag = _compute_predictive_diagnostics(
        contained_count=np.array([5.0]),
        boundary_count=np.array([3.0]),
        contained_leff=np.array([20.0]),
        boundary_leff=np.array([7.0]),
        alpha_prior=np.array([2.0]),
        beta_prior=np.array([10.0]),
        confidence=0.95,
    )

    alpha_post = 5.0
    beta_post = 17.0
    p_nb = beta_post / (beta_post + 20.0)
    mean_c = alpha_post / beta_post * 20.0

    assert diag.mean_c[0] == pytest.approx(mean_c)
    assert diag.p_nb[0] == pytest.approx(p_nb)
    assert diag.upper_unbounded[0] == pytest.approx(3.0 + nbinom.ppf(0.95, alpha_post, p_nb))


def test_fractional_count_convention_floors_contained_for_sf_only() -> None:
    kwargs = dict(
        boundary_count=np.array([3.7]),
        contained_leff=np.array([20.0]),
        boundary_leff=np.array([7.0]),
        alpha_prior=np.array([2.0]),
        beta_prior=np.array([10.0]),
        confidence=0.95,
    )
    d1 = _compute_predictive_diagnostics(contained_count=np.array([5.1]), **kwargs)
    d2 = _compute_predictive_diagnostics(contained_count=np.array([5.9]), **kwargs)
    d3 = _compute_predictive_diagnostics(contained_count=np.array([6.0]), **kwargs)

    assert d1.alpha_post[0] == pytest.approx(5.7)
    assert d1.tail_probability[0] == pytest.approx(d2.tail_probability[0])
    assert d1.expected_tail_count[0] == pytest.approx(d2.expected_tail_count[0])
    assert d1.tail_probability[0] != pytest.approx(d3.tail_probability[0])


def test_closed_form_stop_loss_matches_truncated_sum() -> None:
    diag = _compute_predictive_diagnostics(
        contained_count=np.array([4.0]),
        boundary_count=np.array([1.0]),
        contained_leff=np.array([3.0]),
        boundary_leff=np.array([4.0]),
        alpha_prior=np.array([2.0]),
        beta_prior=np.array([8.0]),
        confidence=0.95,
    )
    c = 4
    k_max = int(nbinom.ppf(1.0 - 1.0e-14, diag.alpha_post[0], diag.p_nb[0])) + 10
    k = np.arange(c + 1, k_max + 1)
    pmf = nbinom.pmf(k, diag.alpha_post[0], diag.p_nb[0])
    truncated = float(np.sum((k - c) * pmf))

    assert diag.expected_tail_count[0] == pytest.approx(truncated, abs=1.0e-9)


def test_tension_flag_fires_only_for_surprisingly_low_observed_contained_count() -> None:
    C = np.concatenate([np.full(25, 10.0), np.array([10.0, 0.0])])
    Lc = np.concatenate([np.full(25, 1000.0), np.array([1000.0, 5000.0])])
    inter = np.concatenate([np.ones(25, dtype=bool), np.array([False, False])])
    observation = _obs(C, Lc, intergenic=inter)

    evidence = fit_density_evidence(observation)

    assert (evidence.flags[-2] & FLAG_HIGH_TAIL_TENSION) == 0
    assert (evidence.flags[-1] & FLAG_HIGH_TAIL_TENSION) != 0


def test_flag_bits_set_on_handcrafted_region() -> None:
    C = np.concatenate([np.full(25, 10.0), np.array([0.0])])
    Lc = np.concatenate([np.full(25, 1000.0), np.array([DENSITY_MIN_EFF_LENGTH / 2.0])])
    boundary = np.zeros_like(C)
    Lb = np.zeros_like(C)
    inter = np.concatenate([np.ones(25, dtype=bool), np.array([False])])
    observation = _obs(C, Lc, boundary=boundary, boundary_leff=Lb, intergenic=inter)

    evidence = fit_density_evidence(observation)
    flag = int(evidence.flags[-1])

    assert flag & FLAG_LOW_CONTAINED_OPPORTUNITY
    assert flag & FLAG_LOW_BOUNDARY_OPPORTUNITY
    assert flag & FLAG_PRIOR_DOMINATED
    assert flag & FLAG_NON_ANCHOR
    assert flag & FLAG_FALLBACK_USED


def test_sparse_family_falls_back_to_all_then_broad_prior() -> None:
    C = np.concatenate([np.full(5, 10.0), np.full(25, 12.0), np.array([0.0])])
    Lc = np.full(C.size, 1000.0)
    inter = np.concatenate([np.ones(5, dtype=bool), np.zeros(26, dtype=bool)])
    intron = np.concatenate([np.zeros(5, dtype=bool), np.ones(25, dtype=bool), [False]])
    observation = _obs(C, Lc, intergenic=inter, intron=intron)

    evidence = fit_density_evidence(observation)
    assert np.all(evidence.fallback_depth[:5] == 1)
    assert np.all(evidence.prior_family[:5] == PRIOR_FAMILY_ALL)

    broad_alpha, broad_beta, broad_family, broad_depth = _assign_region_priors(
        observation,
        {"INTRON": evidence.priors["INTRON"]},
        rho_ref=evidence.rho_ref,
    )
    assert np.all(broad_depth[:5] == 2)
    assert np.all(broad_family[:5] == PRIOR_FAMILY_FALLBACK_BROAD)
    assert np.all(broad_alpha[:5] == DENSITY_FALLBACK_PRIOR_ALPHA)
    assert np.all(broad_beta[:5] == pytest.approx(DENSITY_FALLBACK_PRIOR_ALPHA / evidence.rho_ref))