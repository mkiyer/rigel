"""Tests for ``fuse_density_and_strand`` (PR 07 Phase 3 — FMA fusion engine)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.density_model import DensityEvidence
from rigel.calibration.fusion import (
    FUSED_CLIPPED_TO_TOTAL,
    FUSED_DEGENERATE,
    FUSED_STRAND_STRUCTURAL_ABSENT,
    FusedRegionGdnaEvidence,
    fuse_density_and_strand,
)
from rigel.calibration.strand_evidence import StrandGdnaEvidence


# ---------------------------------------------------------------------------
# Synthetic constructors
# ---------------------------------------------------------------------------


def _density_evidence(
    *,
    mu: np.ndarray,
    var: np.ndarray,
    info: np.ndarray | None = None,
    applicable: np.ndarray | None = None,
) -> DensityEvidence:
    n = int(mu.size)
    if info is None:
        info = np.where(var > 1e-12, 1.0 / np.maximum(var, 1e-12), 0.0)
    if applicable is None:
        applicable = np.ones(n, dtype=bool)
    return DensityEvidence(
        rho_post=np.zeros(n, dtype=np.float64),
        relative_exposure=np.ones(n, dtype=np.float64),
        mean_unbounded=np.asarray(mu, dtype=np.float64),
        upper_unbounded=np.asarray(mu, dtype=np.float64),
        prior_family=np.zeros(n, dtype=np.uint8),
        fallback_depth=np.zeros(n, dtype=np.uint8),
        flags=np.zeros(n, dtype=np.uint8),
        confidence=0.95,
        priors={},
        rho_ref=1.0,
        rho_ref_source="TEST",
        variance_unbounded=np.asarray(var, dtype=np.float64),
        information=np.asarray(info, dtype=np.float64),
        applicable=np.asarray(applicable, dtype=bool),
    )


def _strand_evidence(
    *,
    Ts: np.ndarray,
    Ta: np.ndarray,
    p_r1_sense: float = 0.99,
    global_info_scale: float = 1.0,
    region_info_gain: np.ndarray | None = None,
    applicable: np.ndarray | None = None,
    support_total: np.ndarray | None = None,
) -> StrandGdnaEvidence:
    Ts = np.asarray(Ts, dtype=np.float64)
    Ta = np.asarray(Ta, dtype=np.float64)
    n = int(Ts.size)
    T = Ts + Ta
    if region_info_gain is None:
        # Saturate so g_lib · g_region · applicable == g_lib (≈ kappa_eff = kappa_lib).
        region_info_gain = np.ones(n, dtype=np.float64)
    if applicable is None:
        applicable = np.ones(n, dtype=bool)
    if support_total is None:
        support_total = np.maximum(np.ceil(T).astype(np.uint64), np.uint64(1))
    kappa_lib = min(p_r1_sense, 1.0 - p_r1_sense)
    return StrandGdnaEvidence(
        n_sense=Ts,
        n_anti=Ta,
        n_total=T,
        support_total=np.asarray(support_total, dtype=np.uint64),
        eff_support=Ts + Ta,
        kappa=np.full(n, kappa_lib, dtype=np.float64),
        applicable=np.asarray(applicable, dtype=bool),
        structural_absent=~np.asarray(applicable, dtype=bool),
        p_r1_sense=float(p_r1_sense),
        global_info_scale=float(global_info_scale),
        region_info_gain=np.asarray(region_info_gain, dtype=np.float64),
    )


def _zeros_supports(n: int) -> dict:
    z = np.zeros(n, dtype=np.uint32)
    return {"n_contained": z, "n_left": z, "n_right": z}


# ---------------------------------------------------------------------------
# Edge-case limit tests (the heart of the FMA correctness story)
# ---------------------------------------------------------------------------


class TestLimits:
    def test_pi_zero_kappa_zero_yields_T_anti(self):
        """Strict (π=0, κ=0) limit: D̂_r = T_anti (NOT 2·T_anti)."""
        T_sense = np.array([5.0])
        T_anti = np.array([3.0])
        dens = _density_evidence(mu=np.array([0.0]), var=np.array([1.0]))
        # κ_lib = 0 requires p_r1_sense = 1.0 exactly.
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=1.0)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.gdna_mass[0] == pytest.approx(T_anti[0], abs=1e-12)

    def test_pi_positive_kappa_zero_approaches_two_T_anti(self):
        """Monotone recovery toward 2·T_anti as π_r grows under κ → 0."""
        T_sense = np.array([5.0])
        T_anti = np.array([3.0])
        # Build a density evidence with a sweeping prior π = μ/T.
        T_total = T_sense + T_anti  # 8.0
        prev = T_anti[0] - 1e-9
        D_hats = []
        for pi in [0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]:
            mu = pi * T_total
            dens = _density_evidence(mu=mu, var=np.array([1.0]))
            strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=1.0)
            out = fuse_density_and_strand(
                density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
            )
            assert out.gdna_mass[0] >= prev - 1e-12
            prev = out.gdna_mass[0]
            D_hats.append(float(out.gdna_mass[0]))
        # At π=0.5: q_sense = 0.25/(0.25+0.5) = 1/3; q_anti = 1.
        # D̂ = 5·(1/3) + 3·1 = 4.667, vs 2·T_anti = 6. So we expect to
        # have moved well above T_anti=3 but not yet at 2·T_anti.
        assert D_hats[-1] > T_anti[0]
        assert D_hats[-1] < 2.0 * T_anti[0]
        # Closed-form check at π=0.5.
        assert D_hats[-1] == pytest.approx(5.0 / 3.0 + 3.0, rel=1e-12)

    def test_kappa_half_unstranded_equals_mu_density(self):
        """κ=0.5 (p_r1_sense=0.5): D̂_r = T_r · π_r = μ_density."""
        T_sense = np.array([5.0])
        T_anti = np.array([3.0])
        mu = np.array([2.0])
        dens = _density_evidence(mu=mu, var=np.array([1.0]))
        # With p_r1_sense=0.5, global_info_scale=0 → g=0 → κ_eff=0.5.
        strand = _strand_evidence(
            Ts=T_sense, Ta=T_anti, p_r1_sense=0.5, global_info_scale=0.0,
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.gdna_mass[0] == pytest.approx(mu[0], abs=1e-12)

    def test_pi_one_yields_T(self):
        """π_r → 1 (μ_density ≥ T): D̂_r = T_r (gDNA-dominated)."""
        T_sense = np.array([5.0])
        T_anti = np.array([3.0])
        T_total = T_sense + T_anti
        dens = _density_evidence(mu=T_total, var=np.array([1.0]))
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=1.0)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.gdna_mass[0] == pytest.approx(T_total[0], abs=1e-12)

    def test_structurally_absent_uses_density_only(self):
        """applicable=False: q_sense = q_anti = π_r, so D̂_r = T_r · π_r."""
        T_sense = np.array([5.0])
        T_anti = np.array([3.0])
        mu = np.array([2.0])
        dens = _density_evidence(mu=mu, var=np.array([1.0]))
        strand = _strand_evidence(
            Ts=T_sense, Ta=T_anti, p_r1_sense=0.99, applicable=np.array([False]),
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.gdna_mass[0] == pytest.approx(mu[0], abs=1e-12)
        # And the structural-absent flag is set.
        assert out.flags[0] & FUSED_STRAND_STRUCTURAL_ABSENT

    def test_T_zero_yields_D_zero_and_no_nan(self):
        T_sense = np.array([0.0])
        T_anti = np.array([0.0])
        dens = _density_evidence(mu=np.array([0.0]), var=np.array([1.0]))
        strand = _strand_evidence(
            Ts=T_sense, Ta=T_anti, p_r1_sense=0.99,
            support_total=np.array([0], dtype=np.uint64),
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.gdna_mass[0] == 0.0
        for arr in (out.gdna_mass, out.rna_mass, out.q_sense, out.q_anti,
                    out.density_weight, out.strand_weight):
            assert np.all(np.isfinite(arr))
        assert out.flags[0] & FUSED_DEGENERATE

    def test_pi_zero_kappa_eff_one_guarded(self):
        """π=0 ∧ κ_eff=1 would naively yield 0/0 in q_anti denom. Guarded."""
        # κ_lib = 0 ⇒ denom_a = 0.5·0 + 1·(1-0) = 1, fine. So craft
        # the truly degenerate case: p_r1_sense=1, π=0 → denom_a=1.
        # The actual 0/0 case is π=0 with both denoms = 0, which never
        # happens when 0 ≤ κ_eff ≤ 1. Test the boundary still finite.
        T_sense = np.array([2.0])
        T_anti = np.array([1.0])
        dens = _density_evidence(mu=np.array([0.0]), var=np.array([1.0]))
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=1.0)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert np.all(np.isfinite(out.gdna_mass))
        assert np.all(np.isfinite(out.variance))

    def test_N_total_zero_var_strand_zero(self):
        """No fragments touching: strand contributes zero precision."""
        T_sense = np.array([0.0])
        T_anti = np.array([0.0])
        # But density has nonzero info (e.g. via prior).
        dens = _density_evidence(mu=np.array([1.0]), var=np.array([2.0]))
        strand = _strand_evidence(
            Ts=T_sense, Ta=T_anti, p_r1_sense=0.99,
            support_total=np.array([0], dtype=np.uint64),
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert out.strand_information[0] == 0.0
        assert out.density_weight[0] == pytest.approx(1.0, abs=1e-12)
        assert out.strand_weight[0] == pytest.approx(0.0, abs=1e-12)
        assert out.variance[0] == pytest.approx(2.0, rel=1e-12)


# ---------------------------------------------------------------------------
# Output schema + invariants
# ---------------------------------------------------------------------------


class TestSchema:
    def test_dataclass_shape_and_dtype(self):
        n = 5
        T_sense = np.full(n, 4.0)
        T_anti = np.full(n, 2.0)
        dens = _density_evidence(
            mu=np.full(n, 1.0), var=np.full(n, 1.5),
        )
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(n),
        )
        assert isinstance(out, FusedRegionGdnaEvidence)
        for arr in (out.total_mass, out.gdna_mass, out.rna_mass, out.variance,
                    out.pi_prior, out.q_sense, out.q_anti,
                    out.density_information, out.strand_information,
                    out.density_weight, out.strand_weight):
            assert arr.dtype == np.float64
            assert arr.shape == (n,)
        for arr in (out.n_contained, out.n_left, out.n_right):
            assert arr.dtype == np.uint32
        assert out.flags.dtype == np.uint16

    def test_total_equals_sum(self):
        T_sense = np.array([5.0, 1.0, 0.0])
        T_anti = np.array([3.0, 1.0, 0.0])
        dens = _density_evidence(mu=np.array([2.0, 1.0, 0.0]), var=np.full(3, 1.0))
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(3),
        )
        np.testing.assert_allclose(out.gdna_mass + out.rna_mass, out.total_mass, atol=1e-12)

    def test_weights_sum_to_one_or_zero(self):
        n = 4
        dens = _density_evidence(mu=np.full(n, 1.0), var=np.full(n, 1.5))
        strand = _strand_evidence(
            Ts=np.full(n, 3.0), Ta=np.full(n, 1.0),
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(n),
        )
        s = out.density_weight + out.strand_weight
        # Where any information exists, weights sum to 1.
        assert np.all(np.isclose(s, 1.0, atol=1e-12) | np.isclose(s, 0.0, atol=1e-12))

    def test_clipped_flag_when_naive_overshoots(self):
        """If naive D̂ > T (e.g. pathological strand fit), CLIPPED flag set."""
        # Construct π=1 with T_anti only -> D̂ = T_anti·1 + T_sense·(0.5/(0.5)) = T_sense + T_anti = T.
        # Tweak to overshoot: set μ > T to ensure π=1 and the formula sits at T.
        # FMA naturally clamps π; to force overshoot we need q_sense > 1 which can't happen.
        # So this scenario primarily verifies the flag is *not* set under normal inputs.
        T_sense = np.array([2.0])
        T_anti = np.array([1.0])
        dens = _density_evidence(mu=np.array([0.5]), var=np.array([1.0]))
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=0.95)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        assert not (out.flags[0] & FUSED_CLIPPED_TO_TOTAL)
        assert 0.0 <= out.gdna_mass[0] <= out.total_mass[0]


# ---------------------------------------------------------------------------
# Continuity (no thresholds → smooth in every input)
# ---------------------------------------------------------------------------


class TestContinuity:
    def _sweep_gdna(self, p_grid: np.ndarray) -> np.ndarray:
        out = np.empty(p_grid.size)
        T_sense = np.array([5.0])
        T_anti = np.array([2.0])
        for i, p in enumerate(p_grid):
            # Build strand with p_r1_sense and a fixed library-confidence
            # to avoid coupling g_lib to the sweep variable.
            strand = _strand_evidence(
                Ts=T_sense, Ta=T_anti, p_r1_sense=float(p),
                global_info_scale=1.0,
            )
            dens = _density_evidence(mu=np.array([1.0]), var=np.array([1.0]))
            o = fuse_density_and_strand(
                density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
            )
            out[i] = o.gdna_mass[0]
        return out

    def test_smooth_in_p_r1_sense(self):
        """Doubling grid resolution must shrink max first-difference."""
        coarse = self._sweep_gdna(np.linspace(0.5, 1.0, 51))
        fine = self._sweep_gdna(np.linspace(0.5, 1.0, 101))
        jc = float(np.max(np.abs(np.diff(coarse))))
        jf = float(np.max(np.abs(np.diff(fine))))
        assert jf < 0.75 * jc

    def test_smooth_in_mu_density(self):
        T_sense = np.array([5.0])
        T_anti = np.array([2.0])
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=0.95)

        def sweep(npts):
            mu_grid = np.linspace(0.0, T_sense[0] + T_anti[0], npts)
            vals = np.empty(npts)
            for i, mu in enumerate(mu_grid):
                dens = _density_evidence(mu=np.array([mu]), var=np.array([1.0]))
                o = fuse_density_and_strand(
                    density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
                )
                vals[i] = o.gdna_mass[0]
            return vals

        coarse = sweep(51)
        fine = sweep(101)
        # Monotone increasing in μ_density.
        assert np.all(np.diff(coarse) >= -1e-12)
        # Smooth.
        jc = float(np.max(np.abs(np.diff(coarse))))
        jf = float(np.max(np.abs(np.diff(fine))))
        assert jf < 0.75 * jc

    def test_smooth_in_T_anti(self):
        """Sweep T_anti from 0 to large; D̂ continuous, no jump at T_anti=0."""
        T_sense = np.array([3.0])
        dens = _density_evidence(mu=np.array([1.0]), var=np.array([1.0]))

        def sweep(npts):
            ta_grid = np.linspace(0.0, 5.0, npts)
            vals = np.empty(npts)
            for i, ta in enumerate(ta_grid):
                strand = _strand_evidence(
                    Ts=T_sense, Ta=np.array([ta]), p_r1_sense=0.99,
                )
                o = fuse_density_and_strand(
                    density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
                )
                vals[i] = o.gdna_mass[0]
            return vals

        coarse = sweep(51)
        fine = sweep(101)
        jc = float(np.max(np.abs(np.diff(coarse))))
        jf = float(np.max(np.abs(np.diff(fine))))
        assert jf < 0.75 * jc


# ---------------------------------------------------------------------------
# Sentinels for downstream phases
# ---------------------------------------------------------------------------


class TestPR06Sentinel:
    """Strong strand-specific + zero density prior → D̂ ≈ T_anti (no `2·T_anti`)."""

    def test_anti_intron_like_scenario(self):
        # Region with mostly antisense unspliced mass (intronic gDNA contam suspected).
        T_sense = np.array([1.0])
        T_anti = np.array([8.0])
        # Density prior strictly zero — no boundary evidence.
        dens = _density_evidence(mu=np.array([0.0]), var=np.array([0.0]),
                                 info=np.array([0.0]))
        strand = _strand_evidence(Ts=T_sense, Ta=T_anti, p_r1_sense=1.0)
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        # FMA produces T_anti (it does not impute matched sense gDNA without π_r > 0).
        assert out.gdna_mass[0] == pytest.approx(T_anti[0], abs=1e-12)
        # And density info is zero → LOW_INFORMATION flag if strand also zero.
        # (Here strand info is non-zero from antisense fragments, so no flag.)


class TestMultiMapperDilution:
    """``g_region`` shrinkage propagates into ``kappa_eff`` toward 0.5.

    Note: the absolute strand_weight in diluted regions is not bounded
    below the well-supported regions because Var(D\u0302) scales as (T/N)^2,
    making absolute variance shrink faster than g_region penalizes
    precision. The proper dilution sentinel lives in Phase 7 with
    multi-mapper simulation; here we only verify the schema invariant
    that g_region actively shifts ``q_sense, q_anti`` toward ``\u03c0``.
    """

    def test_low_g_region_pushes_q_toward_pi(self):
        T_sense = np.array([0.05])
        T_anti = np.array([0.10])
        T_total = T_sense + T_anti  # 0.15
        mu = np.array([0.05])       # \u03c0 = 1/3
        dens = _density_evidence(mu=mu, var=np.array([1.0]))
        # Tiny g_region (diluted region).
        strand = _strand_evidence(
            Ts=T_sense, Ta=T_anti, p_r1_sense=0.99,
            region_info_gain=np.array([T_total[0] / (T_total[0] + 4.0)]),
        )
        out = fuse_density_and_strand(
            density_evidence=dens, strand_evidence=strand, **_zeros_supports(1),
        )
        pi = mu[0] / T_total[0]
        # Both q's should be close to \u03c0 because kappa_eff was pushed near 0.5.
        # Generous tolerance: g_region ~ 0.036 so kappa_eff ~ 0.482, q's near \u03c0
        # but not exact.
        assert abs(out.q_sense[0] - pi) < 0.15
        assert abs(out.q_anti[0] - pi) < 0.15
