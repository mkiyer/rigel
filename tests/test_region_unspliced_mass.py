"""Tests for PR 03 RegionUnsplicedMass + BackgroundDensity dataclasses.

See docs/newcalib/pr03_impl_plan_v3.md Section 9 for the full 18-case spec.
This file lands incrementally; cases marked with ``pytest.mark.skip(reason="Step N")``
are unblocked as the implementation order in Section 11 progresses.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.background_model import BackgroundModel
from rigel.calibration.calibration_iteration import (
    FLAG_M_CLIPPED_TO_TOTAL,
    FLAG_M_IMPUTED_FROM_BACKGROUND,
    METHOD_BACKGROUND_FALLBACK,
    METHOD_BOUNDARY,
    METHOD_STRAND,
    BackgroundDensity,
    RegionUnsplicedMass,
)


def _make_mass(
    *,
    total: list[float],
    gdna: list[float],
    region_size: list[float] | None = None,
    counts: list[int] | None = None,
    method: list[int] | None = None,
    precision: list[float] | None = None,
    flags: list[int] | None = None,
) -> RegionUnsplicedMass:
    n = len(total)
    total_arr = np.asarray(total, dtype=np.float64)
    gdna_arr = np.asarray(gdna, dtype=np.float64)
    rna_arr = total_arr - gdna_arr  # exact in float64
    return RegionUnsplicedMass(
        total_mass=total_arr,
        gdna_mass=gdna_arr,
        rna_mass=rna_arr,
        region_size_bp=np.asarray(region_size or [100.0] * n, dtype=np.float64),
        unspliced_counts=np.asarray(counts or [1] * n, dtype=np.uint64),
        method=np.asarray(method or [METHOD_STRAND] * n, dtype=np.uint8),
        precision=np.asarray(precision or [1.0] * n, dtype=np.float64),
        flags=np.asarray(flags or [0] * n, dtype=np.uint16),
    )


# ---------------------------------------------------------------------------
# Tier-hierarchy correctness
# ---------------------------------------------------------------------------


def test_01_mass_conservation_exact_float64() -> None:
    """Case 1: gdna_mass + rna_mass == total_mass exactly in float64 (no tolerance)."""
    # Mix Tier 1/2/3 sized regions with non-trivial values.
    mass = _make_mass(
        total=[0.0, 10.0, 123.456789012345, 1e6, 1e-9],
        gdna=[0.0, 3.5, 50.0, 999.9999, 5e-10],
    )
    # Exact equality, not allclose.
    assert np.array_equal(mass.gdna_mass + mass.rna_mass, mass.total_mass)


def test_05_sentinel_empty_region() -> None:
    """Case 5: T_r == 0 and unspliced_counts == 0 -> M_r == R_r == 0."""
    mass = _make_mass(
        total=[0.0, 0.0],
        gdna=[0.0, 0.0],
        counts=[0, 0],
        method=[METHOD_BACKGROUND_FALLBACK, METHOD_STRAND],
        flags=[FLAG_M_IMPUTED_FROM_BACKGROUND, 0],
        precision=[0.0, 1.0],
    )
    assert mass.gdna_mass.tolist() == [0.0, 0.0]
    assert mass.rna_mass.tolist() == [0.0, 0.0]


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------


def test_16_dtypes_match_plan_contract() -> None:
    """Case 16: primary tensors float64; unspliced_counts uint64; method uint8; flags uint16."""
    mass = _make_mass(total=[1.0, 2.0], gdna=[0.5, 1.0])
    assert mass.total_mass.dtype == np.float64
    assert mass.gdna_mass.dtype == np.float64
    assert mass.rna_mass.dtype == np.float64
    assert mass.region_size_bp.dtype == np.float64
    assert mass.precision.dtype == np.float64
    assert mass.unspliced_counts.dtype == np.uint64
    assert mass.method.dtype == np.uint8
    assert mass.flags.dtype == np.uint16

    bg = BackgroundDensity(
        rho0_mean=0.01,
        alpha0=1.0,
        beta0=100.0,
        log_dispersion=np.log(2.0),
        n_effective_regions=10.0,
        n_regions_in_pool=5,
        method_histogram=(3, 2, 0),
        fit_status="converged",
    )
    assert isinstance(bg.rho0_mean, float)
    assert isinstance(bg.alpha0, float)
    assert isinstance(bg.beta0, float)
    assert isinstance(bg.log_dispersion, float)
    assert isinstance(bg.n_effective_regions, float)
    assert isinstance(bg.n_regions_in_pool, int)


def test_18_no_capture_latent_fields() -> None:
    """Case 18: no PR 03 dataclass exposes p_captured / gamma_r / capture_enrichment_target."""
    mass = _make_mass(total=[1.0], gdna=[0.5])
    forbidden = ("p_captured", "gamma_r", "capture_enrichment_target")
    for name in forbidden:
        assert not hasattr(mass, name)
    bg = BackgroundDensity(
        rho0_mean=0.01,
        alpha0=1.0,
        beta0=100.0,
        log_dispersion=0.5,
        n_effective_regions=0.0,
        n_regions_in_pool=0,
        method_histogram=(0, 0, 0),
        fit_status="fallback_bootstrap",
    )
    for name in forbidden:
        assert not hasattr(bg, name)


# ---------------------------------------------------------------------------
# Bootstrap handoff (Test 13 -- lands together with dataclass since it has no
# dependencies on the M_r builder or estimator).
# ---------------------------------------------------------------------------


def test_13_bootstrap_handoff() -> None:
    """Case 13: BackgroundDensity.from_bootstrap copies model fields and sets fallback status."""
    model = BackgroundModel(
        rho_off_alpha=2.5,
        rho_off_beta=250.0,
        rho_off_mean=2.5 / 250.0,
        seed_mask=np.zeros(3, dtype=bool),
        top_t_exclusion_mask=np.zeros(3, dtype=bool),
        n_seed_regions=0,
        n_fragments=1.0,
        eff_length=100.0,
        fit_status="ok",
        flags=np.zeros(3, dtype=np.uint16),
    )
    bg = BackgroundDensity.from_bootstrap(model)
    assert bg.rho0_mean == pytest.approx(model.rho_off_mean)
    assert bg.alpha0 == pytest.approx(model.rho_off_alpha)
    assert bg.beta0 == pytest.approx(model.rho_off_beta)
    assert bg.n_effective_regions == 0.0
    assert bg.n_regions_in_pool == 0
    assert bg.method_histogram == (0, 0, 0)
    assert bg.fit_status == "fallback_bootstrap"
    assert bg.log_dispersion > 0.0


# ---------------------------------------------------------------------------
# Validation guards
# ---------------------------------------------------------------------------


def test_post_init_rejects_mass_overflow() -> None:
    with pytest.raises(ValueError, match="gdna_mass must not exceed total_mass"):
        _make_mass(total=[1.0], gdna=[2.0])


def test_post_init_rejects_invalid_method() -> None:
    with pytest.raises(ValueError, match="method must be one of"):
        _make_mass(
            total=[1.0],
            gdna=[0.5],
            method=[99],
        )


def test_post_init_rejects_nonpositive_region_size_bp() -> None:
    with pytest.raises(ValueError, match="region_size_bp must be strictly positive"):
        _make_mass(total=[1.0], gdna=[0.5], region_size=[0.0])


def test_post_init_background_density_rejects_bad_status() -> None:
    with pytest.raises(ValueError, match="fit_status must be one of"):
        BackgroundDensity(
            rho0_mean=0.01,
            alpha0=1.0,
            beta0=100.0,
            log_dispersion=0.1,
            n_effective_regions=0.0,
            n_regions_in_pool=0,
            method_histogram=(0, 0, 0),
            fit_status="weird",
        )


def test_flag_bits_disjoint_from_legacy_flags() -> None:
    """Sanity: the new region-level flags do not collide numerically with each other."""
    assert FLAG_M_IMPUTED_FROM_BACKGROUND != FLAG_M_CLIPPED_TO_TOTAL
    assert FLAG_M_IMPUTED_FROM_BACKGROUND.bit_count() == 1
    assert FLAG_M_CLIPPED_TO_TOTAL.bit_count() == 1


def test_methods_are_distinct_positive_uint8() -> None:
    """Sanity: METHOD_STRAND / METHOD_BOUNDARY / METHOD_BACKGROUND_FALLBACK are distinct."""
    values = {METHOD_STRAND, METHOD_BOUNDARY, METHOD_BACKGROUND_FALLBACK}
    assert len(values) == 3
    for v in values:
        assert 0 < v <= 255


# ---------------------------------------------------------------------------
# Three-tier M_r hierarchy (build_region_unspliced_mass)
# ---------------------------------------------------------------------------


from types import SimpleNamespace  # noqa: E402

from rigel.calibration.calibration_iteration import build_region_unspliced_mass  # noqa: E402


def _bd(rho0: float = 0.01) -> BackgroundDensity:
    return BackgroundDensity(
        rho0_mean=rho0,
        alpha0=1.0,
        beta0=1.0 / rho0,
        log_dispersion=0.5,
        n_effective_regions=0.0,
        n_regions_in_pool=0,
        method_histogram=(0, 0, 0),
        fit_status="fallback_bootstrap",
    )


def _obs(total: list[float]) -> SimpleNamespace:
    arr = np.asarray(total, dtype=np.float64)
    return SimpleNamespace(observed_compatible_count=arr)


def _empty_local(n: int) -> SimpleNamespace:
    return SimpleNamespace(
        alpha_excess=np.zeros(n, dtype=np.float64),
        beta_excess=np.zeros(n, dtype=np.float64),
    )


def _local(alpha: list[float], beta: list[float]) -> SimpleNamespace:
    return SimpleNamespace(
        alpha_excess=np.asarray(alpha, dtype=np.float64),
        beta_excess=np.asarray(beta, dtype=np.float64),
    )


def _sweep(mu: list[float]) -> SimpleNamespace:
    return SimpleNamespace(mu_sweep=np.asarray(mu, dtype=np.float64))


def _strand(
    *,
    n: int,
    contained_mean: list[float] | None = None,
    contained_reliability: list[float] | None = None,
    contained_precision: list[float] | None = None,
    left_mean: list[float] | None = None,
    left_reliability: list[float] | None = None,
    left_precision: list[float] | None = None,
    right_mean: list[float] | None = None,
    right_reliability: list[float] | None = None,
    right_precision: list[float] | None = None,
    flags: list[int] | None = None,
) -> SimpleNamespace:
    def _arr(values: list[float] | None) -> np.ndarray:
        return np.asarray(values if values is not None else [0.0] * n, dtype=np.float64)

    return SimpleNamespace(
        contained_mean=_arr(contained_mean),
        contained_reliability=_arr(contained_reliability),
        contained_precision=_arr(contained_precision),
        boundary_left_mean=_arr(left_mean),
        boundary_left_reliability=_arr(left_reliability),
        boundary_left_precision=_arr(left_precision),
        boundary_right_mean=_arr(right_mean),
        boundary_right_reliability=_arr(right_reliability),
        boundary_right_precision=_arr(right_precision),
        flags=np.asarray(flags or [0] * n, dtype=np.uint16),
    )


def test_02_tier1_strand_uses_reliability_weighted_sum() -> None:
    """Case 2: strong strand evidence -> METHOD_STRAND, weighted-sum M_r, precision > 0."""
    n = 1
    strand = _strand(
        n=n,
        contained_mean=[10.0],
        contained_reliability=[0.8],
        contained_precision=[0.9],
        left_mean=[2.0],
        left_reliability=[0.5],
        left_precision=[0.7],
        right_mean=[1.0],
        right_reliability=[0.1],
        right_precision=[0.6],
    )
    mass = build_region_unspliced_mass(
        _obs([100.0]),
        region_size_bp=np.asarray([1000.0], dtype=np.float64),
        unspliced_counts=np.asarray([50], dtype=np.uint64),
        strand_channels=strand,
        local_posterior=_empty_local(n),
        sweep=_sweep([0.0]),
        background_density=_bd(),
    )
    assert mass.method.tolist() == [METHOD_STRAND]
    expected = 10.0 * 0.8 + 2.0 * 0.5 + 1.0 * 0.1
    assert mass.gdna_mass[0] == pytest.approx(expected)
    # precision = max channel-product = max(0.9*0.8, 0.7*0.5, 0.6*0.1) = 0.72
    assert mass.precision[0] == pytest.approx(0.72)
    assert mass.rna_mass[0] == pytest.approx(100.0 - expected)
    assert mass.flags[0] & FLAG_M_IMPUTED_FROM_BACKGROUND == 0


def test_03_tier2_boundary_imputation() -> None:
    """Case 3: no strand contrast but boundary excess > 0 -> METHOD_BOUNDARY."""
    mass = build_region_unspliced_mass(
        _obs([100.0]),
        region_size_bp=np.asarray([1000.0], dtype=np.float64),
        unspliced_counts=np.asarray([10], dtype=np.uint64),
        strand_channels=None,
        local_posterior=_local(alpha=[3.0], beta=[2.0]),
        sweep=_sweep([42.0]),
        background_density=_bd(),
    )
    assert mass.method.tolist() == [METHOD_BOUNDARY]
    # mu_sweep clipped to [0, T_r=100] -> 42.0
    assert mass.gdna_mass[0] == pytest.approx(42.0)
    # precision proxy is alpha_excess + beta_excess (PR03 deviation: see plan note).
    assert mass.precision[0] == pytest.approx(5.0)
    assert mass.rna_mass[0] == pytest.approx(58.0)
    assert mass.flags[0] & FLAG_M_IMPUTED_FROM_BACKGROUND == 0


def test_04_tier3_background_fallback() -> None:
    """Case 4: neither strand nor boundary -> METHOD_BACKGROUND_FALLBACK + flag bit."""
    n = 1
    rho0 = 0.05
    region_bp = 1000.0
    total = 100.0
    mass = build_region_unspliced_mass(
        _obs([total]),
        region_size_bp=np.asarray([region_bp], dtype=np.float64),
        unspliced_counts=np.asarray([7], dtype=np.uint64),
        strand_channels=None,
        local_posterior=_empty_local(n),
        sweep=_sweep([0.0]),
        background_density=_bd(rho0=rho0),
    )
    assert mass.method.tolist() == [METHOD_BACKGROUND_FALLBACK]
    expected = min(rho0 * region_bp, total)
    assert mass.gdna_mass[0] == pytest.approx(expected)
    assert mass.precision[0] == 0.0
    assert mass.flags[0] & FLAG_M_IMPUTED_FROM_BACKGROUND


def test_04b_tier3_background_clips_to_total_when_dense() -> None:
    """Sub-case: rho0 * region_bp > T_r forces clipping and the clip flag fires."""
    n = 1
    mass = build_region_unspliced_mass(
        _obs([5.0]),  # T_r small
        region_size_bp=np.asarray([1000.0], dtype=np.float64),
        unspliced_counts=np.asarray([1], dtype=np.uint64),
        strand_channels=None,
        local_posterior=_empty_local(n),
        sweep=_sweep([0.0]),
        background_density=_bd(rho0=0.1),  # rho0 * bp = 100 > 5
    )
    assert mass.gdna_mass[0] == 5.0
    assert mass.rna_mass[0] == 0.0
    assert mass.flags[0] & FLAG_M_CLIPPED_TO_TOTAL
    assert mass.flags[0] & FLAG_M_IMPUTED_FROM_BACKGROUND


def test_04c_force_zero_gdna_overrides_strand_and_boundary_evidence() -> None:
    """Deterministic-zero density evidence keeps prior mass at D=0."""
    n = 1
    mass = build_region_unspliced_mass(
        _obs([100.0]),
        region_size_bp=np.asarray([1000.0], dtype=np.float64),
        unspliced_counts=np.asarray([50], dtype=np.uint64),
        strand_channels=_strand(
            n=n,
            contained_mean=[80.0],
            contained_reliability=[1.0],
            contained_precision=[1.0],
        ),
        local_posterior=_local(alpha=[10.0], beta=[10.0]),
        sweep=_sweep([75.0]),
        background_density=_bd(rho0=0.05),
        force_zero_gdna=True,
    )
    assert mass.method.tolist() == [METHOD_BACKGROUND_FALLBACK]
    assert mass.gdna_mass[0] == 0.0
    assert mass.rna_mass[0] == 100.0
    assert mass.precision[0] == 0.0
    assert mass.flags[0] & FLAG_M_IMPUTED_FROM_BACKGROUND


def test_06_tier_promotion_with_added_evidence() -> None:
    """Case 6: same region; tier moves BACKGROUND_FALLBACK -> BOUNDARY -> STRAND."""
    n = 1
    total = 100.0
    region_bp = np.asarray([1000.0], dtype=np.float64)
    counts = np.asarray([10], dtype=np.uint64)
    obs = _obs([total])
    bd = _bd(rho0=0.05)

    # 1) Background fallback baseline.
    no_evidence = build_region_unspliced_mass(
        obs,
        region_size_bp=region_bp,
        unspliced_counts=counts,
        strand_channels=None,
        local_posterior=_empty_local(n),
        sweep=_sweep([0.0]),
        background_density=bd,
    )
    assert no_evidence.method.tolist() == [METHOD_BACKGROUND_FALLBACK]

    # 2) Add boundary excess -> Tier 2 wins over Tier 3.
    boundary_only = build_region_unspliced_mass(
        obs,
        region_size_bp=region_bp,
        unspliced_counts=counts,
        strand_channels=None,
        local_posterior=_local(alpha=[2.0], beta=[1.0]),
        sweep=_sweep([35.0]),
        background_density=bd,
    )
    assert boundary_only.method.tolist() == [METHOD_BOUNDARY]
    assert boundary_only.gdna_mass[0] == pytest.approx(35.0)

    # 3) Add strand reliability -> Tier 1 wins over Tier 2.
    strand_promoted = build_region_unspliced_mass(
        obs,
        region_size_bp=region_bp,
        unspliced_counts=counts,
        strand_channels=_strand(
            n=n,
            contained_mean=[20.0],
            contained_reliability=[0.5],
            contained_precision=[1.0],
        ),
        local_posterior=_local(alpha=[2.0], beta=[1.0]),
        sweep=_sweep([35.0]),
        background_density=bd,
    )
    assert strand_promoted.method.tolist() == [METHOD_STRAND]
    assert strand_promoted.gdna_mass[0] == pytest.approx(10.0)  # 20 * 0.5


def test_mass_conservation_all_tiers_in_one_call() -> None:
    """Case 1 extension: a single call exercising all three tiers preserves M+R==T exactly."""
    # Region 0: strand-only; Region 1: boundary-only; Region 2: background fallback.
    n = 3
    totals = np.asarray([50.0, 80.0, 30.0], dtype=np.float64)
    strand = _strand(
        n=n,
        contained_mean=[25.0, 0.0, 0.0],
        contained_reliability=[0.6, 0.0, 0.0],
        contained_precision=[0.8, 0.0, 0.0],
    )
    mass = build_region_unspliced_mass(
        _obs(totals.tolist()),
        region_size_bp=np.asarray([500.0, 1000.0, 1500.0], dtype=np.float64),
        unspliced_counts=np.asarray([10, 5, 1], dtype=np.uint64),
        strand_channels=strand,
        local_posterior=_local(alpha=[0.0, 1.0, 0.0], beta=[0.0, 1.0, 0.0]),
        sweep=_sweep([0.0, 40.0, 0.0]),
        background_density=_bd(rho0=0.02),
    )
    assert mass.method.tolist() == [
        METHOD_STRAND,
        METHOD_BOUNDARY,
        METHOD_BACKGROUND_FALLBACK,
    ]
    assert np.array_equal(mass.gdna_mass + mass.rna_mass, totals)


# ---------------------------------------------------------------------------
# Step 3: estimate_background_density (Section 6.3-6.7)
# ---------------------------------------------------------------------------


from rigel.calibration.calibration_iteration import (  # noqa: E402
    estimate_background_density,
)


def _bd_prev(
    rho0: float = 0.01, alpha0: float = 1.0, beta0: float | None = None
) -> BackgroundDensity:
    """Previous BackgroundDensity that triggers the 'real fit' branch (not prior_only)."""
    return BackgroundDensity(
        rho0_mean=rho0,
        alpha0=alpha0,
        beta0=beta0 if beta0 is not None else alpha0 / rho0,
        log_dispersion=0.5,
        n_effective_regions=10.0,
        n_regions_in_pool=10,
        method_histogram=(5, 3, 2),
        fit_status="iterating",
    )


def _mass_from(
    *,
    gdna: list[float],
    region_bp: list[float],
    counts: list[int],
    method: list[int],
    precision: list[float] | None = None,
    total: list[float] | None = None,
    flags: list[int] | None = None,
) -> RegionUnsplicedMass:
    g = np.asarray(gdna, dtype=np.float64)
    t = np.asarray(total if total is not None else [v * 2 + 1.0 for v in gdna], dtype=np.float64)
    g = np.minimum(g, t)
    r = t - g
    return RegionUnsplicedMass(
        total_mass=t,
        gdna_mass=g,
        rna_mass=r,
        region_size_bp=np.asarray(region_bp, dtype=np.float64),
        unspliced_counts=np.asarray(counts, dtype=np.uint64),
        method=np.asarray(method, dtype=np.uint8),
        precision=np.asarray(precision if precision is not None else [1.0] * len(gdna), dtype=np.float64),
        flags=np.asarray(flags or [0] * len(gdna), dtype=np.uint16),
    )


def test_07_empty_pool_carries_previous_unchanged() -> None:
    """Case 7: zero usable regions -> fallback_bootstrap, previous fields preserved."""
    prev = _bd_prev(rho0=0.02, alpha0=2.0, beta0=100.0)
    # All regions are METHOD_BACKGROUND_FALLBACK -> excluded from pool.
    mass = _mass_from(
        gdna=[5.0, 3.0],
        region_bp=[1000.0, 2000.0],
        counts=[10, 5],
        method=[METHOD_BACKGROUND_FALLBACK, METHOD_BACKGROUND_FALLBACK],
    )
    bg = estimate_background_density(
        mass,
        p_unexpressed=np.asarray([1.0, 1.0]),
        previous=prev,
    )
    assert bg.fit_status == "fallback_bootstrap"
    assert bg.n_regions_in_pool == 0
    assert bg.n_effective_regions == 0.0
    assert bg.rho0_mean == pytest.approx(prev.rho0_mean)
    assert bg.alpha0 == pytest.approx(prev.alpha0)
    assert bg.beta0 == pytest.approx(prev.beta0)
    assert bg.log_dispersion == pytest.approx(prev.log_dispersion)
    assert bg.method_histogram == (0, 0, 2)


def test_08_zero_count_regions_excluded_from_pool() -> None:
    """Case 8: METHOD_STRAND region with unspliced_counts==0 is dropped from pool."""
    prev = _bd_prev(rho0=0.02)
    mass = _mass_from(
        gdna=[5.0, 3.0],
        region_bp=[1000.0, 1000.0],
        counts=[0, 4],  # first region excluded
        method=[METHOD_STRAND, METHOD_BOUNDARY],
    )
    bg = estimate_background_density(mass, p_unexpressed=np.full(2, 1.0), previous=prev)
    assert bg.n_regions_in_pool == 1
    # Only the second region contributes -> rho_hat ~ 3/1000 (after shrinkage + damping).
    # With damping=0.5 the log-mean is between log(0.02) and log(~0.003).
    assert 0.002 < bg.rho0_mean < 0.02


def test_09_bayesian_shrinkage_handles_zero_mass() -> None:
    """Case 9: a strand region with gdna==0 must NOT pin log_center at -inf."""
    prev = _bd_prev(rho0=0.01)
    mass = _mass_from(
        gdna=[0.0, 5.0],          # first region is genuinely clean
        region_bp=[1000.0, 1000.0],
        counts=[8, 8],
        method=[METHOD_STRAND, METHOD_STRAND],
        total=[200.0, 200.0],
    )
    bg = estimate_background_density(mass, p_unexpressed=np.full(2, 1.0), previous=prev)
    assert np.isfinite(bg.rho0_mean) and bg.rho0_mean > 0.0
    assert np.isfinite(bg.log_dispersion)
    # Should not collapse to rho_floor.
    assert bg.rho0_mean > 1e-10


def test_10_huber_bounds_single_region_pull() -> None:
    """Case 10: one mega-enriched outlier moves rho0 by a BOUNDED multiplicative factor."""
    prev = _bd_prev(rho0=0.01)
    # 9 typical regions at ~0.01 density + 1 outlier at ~10x.
    typical_gdna = [10.0] * 9 + [10000.0]
    region_bp = [1000.0] * 10
    counts = [100] * 10
    mass = _mass_from(
        gdna=typical_gdna,
        region_bp=region_bp,
        counts=counts,
        method=[METHOD_STRAND] * 10,
        total=[20000.0] * 10,
    )
    bg = estimate_background_density(
        mass,
        p_unexpressed=np.full(10, 1.0),
        previous=prev,
        damping=1.0,  # disable damping to expose the Huber-only effect
    )
    # Outlier alone would push log_center by log(1000) ~ 6.9; Huber cap with
    # huber_k=1.5 and weighted MAD limits the shift dramatically.
    # We assert the resulting rho0 is below the naive (unclipped) center.
    naive_unclipped = float(np.exp(np.mean(np.log(np.asarray(typical_gdna) / np.asarray(region_bp)))))
    assert bg.rho0_mean < naive_unclipped
    # And still positive / finite.
    assert bg.rho0_mean > 0.0 and np.isfinite(bg.rho0_mean)


def test_11_damping_halves_log_jump() -> None:
    """Case 11: damping=0.5 -> log_center sits halfway between log(prev) and log(target)."""
    prev = _bd_prev(rho0=0.01)
    # Pool consistently points at ~0.1 (10x current).
    mass = _mass_from(
        gdna=[100.0, 100.0, 100.0],
        region_bp=[1000.0, 1000.0, 1000.0],
        counts=[50, 50, 50],
        method=[METHOD_STRAND, METHOD_STRAND, METHOD_STRAND],
        total=[200.0, 200.0, 200.0],
    )
    bg_full = estimate_background_density(
        mass, p_unexpressed=np.full(3, 1.0), previous=prev, damping=1.0
    )
    bg_half = estimate_background_density(
        mass, p_unexpressed=np.full(3, 1.0), previous=prev, damping=0.5
    )
    log_full = np.log(bg_full.rho0_mean)
    log_half = np.log(bg_half.rho0_mean)
    log_prev = np.log(prev.rho0_mean)
    # log_half should be ~ (log_full + log_prev) / 2.
    assert log_half == pytest.approx(0.5 * (log_full + log_prev), abs=1e-9)


def test_12_gamma_view_matches_rho0_mean() -> None:
    """Case 12 / Section 6.4: returned (alpha0, beta0) satisfies alpha0/beta0 == rho0_mean."""
    prev = _bd_prev(rho0=0.05)
    mass = _mass_from(
        gdna=[10.0, 20.0, 5.0],
        region_bp=[500.0, 1000.0, 2000.0],
        counts=[3, 7, 12],
        method=[METHOD_STRAND, METHOD_BOUNDARY, METHOD_STRAND],
        total=[100.0, 100.0, 100.0],
    )
    bg = estimate_background_density(mass, p_unexpressed=np.full(3, 1.0), previous=prev)
    assert bg.alpha0 > 0.0 and bg.beta0 > 0.0
    assert bg.alpha0 / bg.beta0 == pytest.approx(bg.rho0_mean, rel=1e-9)


def test_14_fit_status_converged_below_tol() -> None:
    """Case 14: small log_jump -> fit_status=='converged'."""
    prev = _bd_prev(rho0=0.01)
    # Pool already at ~0.01 -> next iteration barely moves.
    mass = _mass_from(
        gdna=[10.0, 10.0, 10.0, 10.0, 10.0],
        region_bp=[1000.0] * 5,
        counts=[100] * 5,
        method=[METHOD_STRAND] * 5,
        total=[200.0] * 5,
    )
    bg = estimate_background_density(
        mass,
        p_unexpressed=np.full(5, 1.0),
        previous=prev,
        damping=1.0,
        convergence_log_tol=0.1,  # loose so the test is robust to shrinkage
    )
    assert bg.fit_status == "converged"


def test_15_fit_status_prior_only_after_bootstrap() -> None:
    """Case 15: first real fit when previous is the bootstrap stub."""
    prev = BackgroundDensity.from_bootstrap.__func__(BackgroundDensity, None) if False else None  # noqa: E501
    # Simulate a bootstrap-stub previous explicitly:
    prev = BackgroundDensity(
        rho0_mean=0.01,
        alpha0=1.0,
        beta0=100.0,
        log_dispersion=float(np.log(10.0)),
        n_effective_regions=0.0,
        n_regions_in_pool=0,
        method_histogram=(0, 0, 0),
        fit_status="fallback_bootstrap",
    )
    mass = _mass_from(
        gdna=[10.0, 20.0],
        region_bp=[1000.0, 2000.0],
        counts=[5, 10],
        method=[METHOD_STRAND, METHOD_BOUNDARY],
        total=[100.0, 100.0],
    )
    bg = estimate_background_density(mass, p_unexpressed=np.full(2, 1.0), previous=prev)
    assert bg.fit_status == "prior_only"
    assert bg.n_regions_in_pool == 2


def test_method_histogram_counts_all_regions() -> None:
    """The histogram reports total method counts, not just pool members."""
    prev = _bd_prev(rho0=0.01)
    mass = _mass_from(
        gdna=[5.0, 6.0, 7.0, 8.0],
        region_bp=[1000.0] * 4,
        counts=[0, 5, 3, 2],  # first dropped from pool, others stay
        method=[
            METHOD_STRAND,
            METHOD_STRAND,
            METHOD_BOUNDARY,
            METHOD_BACKGROUND_FALLBACK,
        ],
    )
    bg = estimate_background_density(mass, p_unexpressed=np.full(4, 1.0), previous=prev)
    assert bg.method_histogram == (2, 1, 1)
    # Pool excludes counts==0 strand region and the background-fallback region.
    assert bg.n_regions_in_pool == 2


def test_n_effective_regions_equals_weight_sum() -> None:
    """n_effective_regions == sum of pool weights (precision * counts * p_unx)."""
    prev = _bd_prev(rho0=0.01)
    mass = _mass_from(
        gdna=[5.0, 10.0],
        region_bp=[1000.0, 1000.0],
        counts=[4, 6],
        method=[METHOD_STRAND, METHOD_BOUNDARY],
        precision=[0.5, 1.0],
    )
    p_unx = np.asarray([0.8, 0.2])
    bg = estimate_background_density(mass, p_unexpressed=p_unx, previous=prev)
    expected = 0.5 * 4.0 * 0.8 + 1.0 * 6.0 * 0.2
    assert bg.n_effective_regions == pytest.approx(expected)
