"""RNA strand-balance fit: κ_rna from StrandModel, ρ_r_bb by method-of-moments."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.strand_balance import (
    _BB_FLOOR,
    _RHO_R_BB_FALLBACK,
    StrandBalance,
    fit_strand_balance,
)
from rigel.calibration.substrate import CalibrationSubstrate, SubstrateView


def _strand_model(p_r1_sense: float):
    return SimpleNamespace(p_r1_sense=p_r1_sense)


def _empty_view(r: int) -> SubstrateView:
    z = np.zeros(r, dtype=np.int64)
    zf = np.zeros(r, dtype=np.float64)
    return SubstrateView(z.copy(), z.copy(), z.copy(), z.copy(), zf.copy(), zf.copy())


def _substrate_with_contained_spliced(sense, antisense) -> CalibrationSubstrate:
    """Build a substrate carrying spliced (sense, antisense) only in the contained view."""
    sense = np.asarray(sense, dtype=np.int64)
    antisense = np.asarray(antisense, dtype=np.int64)
    r = len(sense)
    z = np.zeros(r, dtype=np.int64)
    zf = np.zeros(r, dtype=np.float64)
    contained = SubstrateView(
        n_unspliced_pos=z.copy(),
        n_unspliced_neg=z.copy(),
        n_spliced_sense=sense,
        n_spliced_antisense=antisense,
        mass_unspliced=zf.copy(),
        mass_spliced=zf.copy(),
    )
    return CalibrationSubstrate(
        n_regions=r,
        L_eff=np.full(r, 100.0),
        ts_class=z.astype(np.int8),
        contained=contained,
        left=_empty_view(r),
        right=_empty_view(r),
    )


def test_kappa_rna_from_strand_model_clamped():
    sub = _substrate_with_contained_spliced([8, 16], [2, 4])
    # Perfectly stranded library (p_r1_sense == 1.0) is clamped into (0, 1).
    sb = fit_strand_balance(sub, _strand_model(1.0))
    assert sb.kappa_rna == 1.0 - _BB_FLOOR
    sb2 = fit_strand_balance(sub, _strand_model(0.8))
    assert sb2.kappa_rna == 0.8


def test_at_mean_data_is_underdispersed_floored():
    # Every observation sits exactly at κ = 0.8 → residual variance below the
    # binomial expectation → ρ floored.
    sub = _substrate_with_contained_spliced([8, 16, 80], [2, 4, 20])
    sb = fit_strand_balance(sub, _strand_model(0.8))
    assert not sb.fallback_used
    assert sb.n_observations == 3
    assert sb.rho_r_bb == _BB_FLOOR


def test_overdispersed_data_yields_large_rho():
    # All-sense and all-antisense observations at κ=0.5 → maximal overdispersion.
    sub = _substrate_with_contained_spliced([10, 0, 10, 0], [0, 10, 0, 10])
    sb = fit_strand_balance(sub, _strand_model(0.5))
    assert not sb.fallback_used
    assert sb.rho_r_bb > 0.5
    assert sb.n_observations == 4
    assert sb.n_total_reads == 40


def test_fallback_no_spliced_observations():
    sub = _substrate_with_contained_spliced([0, 0], [0, 0])  # no spliced reads
    sb = fit_strand_balance(sub, _strand_model(0.9))
    assert sb.fallback_used
    assert sb.rho_r_bb == _RHO_R_BB_FALLBACK
    assert sb.kappa_rna == 0.9
    assert sb.n_observations == 0


def test_fallback_single_observation():
    # One observation: overdispersion is not estimable (needs >= 2).
    sub = _substrate_with_contained_spliced([7], [3])
    sb = fit_strand_balance(sub, _strand_model(0.7))
    assert sb.fallback_used
    assert sb.rho_r_bb == _RHO_R_BB_FALLBACK
    assert sb.n_observations == 1


def test_pools_all_three_views():
    # Spliced observations in contained + left + right views are all pooled.
    r = 2
    z = np.zeros(r, dtype=np.int64)
    zf = np.zeros(r, dtype=np.float64)

    def view(sense, anti):
        return SubstrateView(
            z.copy(),
            z.copy(),
            np.asarray(sense, dtype=np.int64),
            np.asarray(anti, dtype=np.int64),
            zf.copy(),
            zf.copy(),
        )

    sub = CalibrationSubstrate(
        n_regions=r,
        L_eff=np.full(r, 100.0),
        ts_class=z.astype(np.int8),
        contained=view([10, 0], [0, 0]),  # 1 obs (region 0)
        left=view([0, 8], [0, 2]),  # 1 obs (region 1)
        right=view([5, 0], [5, 0]),  # 1 obs (region 0)
    )
    sb = fit_strand_balance(sub, _strand_model(0.5))
    assert sb.n_observations == 3
    assert sb.n_total_reads == 10 + 10 + 10


def test_returns_strand_balance_type():
    sub = _substrate_with_contained_spliced([8, 16], [2, 4])
    sb = fit_strand_balance(sub, _strand_model(0.8))
    assert isinstance(sb, StrandBalance)
    assert 0.0 < sb.kappa_rna < 1.0
    assert 0.0 < sb.rho_r_bb < 1.0
