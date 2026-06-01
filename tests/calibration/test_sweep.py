"""AMBIG boundary↔region sweep: alternating-chain exposure imputation (D7 leg 2)."""

from __future__ import annotations

import numpy as np

from rigel.calibration.estep import Allocation
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, TS_POS
from rigel.calibration.substrate import CalibrationSubstrate, SubstrateView
from rigel.calibration.sweep import _internal_boundary_nodes, sweep_ambig_exposure


def _view(unspliced, spliced=None) -> SubstrateView:
    u = np.asarray(unspliced, dtype=np.int64)
    r = u.shape[0]
    z = np.zeros(r, dtype=np.int64)
    zf = np.zeros(r, dtype=np.float64)
    s = z.copy() if spliced is None else np.asarray(spliced, dtype=np.int64)
    return SubstrateView(u, z.copy(), s, z.copy(), zf.copy(), zf.copy())


def _alloc(m_g_unspl, pi_g) -> Allocation:
    m_g_unspl = np.asarray(m_g_unspl, dtype=np.float64)
    z = np.zeros_like(m_g_unspl)
    return Allocation(
        m_g=z.copy(),
        m_d=z.copy(),
        m_g_unspl=m_g_unspl,
        m_d_unspl=z.copy(),
        pi_g=np.asarray(pi_g, np.float64),
        k_sense=np.zeros_like(m_g_unspl, dtype=np.int64),
    )


def _region_arrays(ts_class) -> RegionArrays:
    r = len(ts_class)
    return RegionArrays(
        ref_id=np.zeros(r, dtype=np.int32),
        start=np.arange(r, dtype=np.int64) * 100,
        end=(np.arange(r, dtype=np.int64) + 1) * 100,
        signature=np.zeros(r, dtype=np.uint8),
        ts_class=np.asarray(ts_class, dtype=np.int8),
        region_size_bp=np.full(r, 100.0),
        ref_offsets=np.array([0, r], dtype=np.int32),
        order=np.arange(r, dtype=np.int64),
        n_refs=1,
    )


def _substrate(ts_class, left, right) -> CalibrationSubstrate:
    r = len(ts_class)
    z = np.zeros(r, dtype=np.int64)
    zf = np.zeros(r, dtype=np.float64)
    contained = SubstrateView(z.copy(), z.copy(), z.copy(), z.copy(), zf.copy(), zf.copy())
    return CalibrationSubstrate(
        n_regions=r,
        L_eff=np.full(r, 100.0),
        ts_class=np.asarray(ts_class, dtype=np.int8),
        contained=contained,
        left=left,
        right=right,
    )


def test_internal_boundary_nodes_clean():
    ts = [TS_POS, TS_AMBIG, TS_POS]
    ra = _region_arrays(ts)
    left = _view([0, 20, 20])  # R0 terminal, R1 left, R2 left
    right = _view([20, 20, 0])  # R0 right, R1 right, R2 terminal
    alloc_left = _alloc([0, 0, 0], [0.5, 0.5, 1.0])
    alloc_right = _alloc([0, 0, 0], [1.0, 0.5, 0.5])
    gdna_flux, weight = _internal_boundary_nodes(ra.ref_id, alloc_left, alloc_right, left, right)
    # boundary 0: 0.5(1.0*20 + 0.5*20)=15; weight=(20/21)*1*0.75; boundary at last region = 0.
    np.testing.assert_allclose(gdna_flux, [15.0, 15.0, 0.0])
    np.testing.assert_allclose(weight, [(20 / 21) * 0.75, (20 / 21) * 0.75, 0.0])


def test_internal_boundary_weight_drops_for_spliced_and_low_purity():
    ts = [TS_POS, TS_POS]
    ra = _region_arrays(ts)
    left = _view([0, 20], spliced=[0, 100])  # boundary carries heavy spliced (RNA)
    right = _view([20, 0], spliced=[100, 0])
    alloc_left = _alloc([0, 0], [0.0, 0.02])  # low gDNA purity at the boundary
    alloc_right = _alloc([0, 0], [0.02, 0.0])
    _, weight = _internal_boundary_nodes(ra.ref_id, alloc_left, alloc_right, left, right)
    assert weight[0] < 0.01  # spliced + low purity ⇒ choked conduit


def _sweep(
    ts,
    a_reg,
    region_eff_len,
    *,
    left,
    right,
    alloc_left,
    alloc_right,
    base_omega,
    mu_fl=50.0,
    rho_0=0.1,
    exposure_dispersion=1.0,
):
    ra = _region_arrays(ts)
    sub = _substrate(ts, left, right)
    alloc_contained = _alloc(a_reg, np.full(len(ts), 0.5))
    base_lov = np.full(len(ts), 1.0)
    return sweep_ambig_exposure(
        sub,
        ra,
        alloc_contained=alloc_contained,
        alloc_left=alloc_left,
        alloc_right=alloc_right,
        region_eff_len=np.asarray(region_eff_len, np.float64),
        mu_fl=mu_fl,
        rho_0=rho_0,
        exposure_dispersion=exposure_dispersion,
        base_omega=np.asarray(base_omega, np.float64),
        base_log_omega_var=base_lov,
    )


def test_ambig_inherits_neighbour_density():
    # R0, R2 decodable with gDNA density 0.2 (= ρ_0·ω with ω=2); R1 AMBIG, no
    # contained gDNA, flanked by clean boundaries → ω imputed toward the neighbours.
    ts = [TS_POS, TS_AMBIG, TS_POS]
    res = _sweep(
        ts,
        a_reg=[10, 0, 10],
        region_eff_len=[50, 50, 50],
        left=_view([0, 20, 20]),
        right=_view([20, 20, 0]),
        alloc_left=_alloc([0, 0, 0], [0.5, 0.5, 1.0]),
        alloc_right=_alloc([0, 0, 0], [1.0, 0.5, 0.5]),
        base_omega=[2.0, 1.0, 2.0],
    )
    assert res.n_ambig == 1
    np.testing.assert_allclose(res.omega[[0, 2]], [2.0, 2.0])  # decodable rows untouched
    assert res.omega[1] > 1.0  # imputed up from the ω=1 fallback...
    assert 1.7 < res.omega[1] < 1.95  # ...toward the neighbours' ω=2


def test_ambig_falls_back_when_isolated_and_choked():
    # Tiny AMBIG region (no contained exposure), boundaries RNA-choked → no
    # evidence reaches it → ω → 1 (the global-fallback prior mean).
    ts = [TS_POS, TS_AMBIG, TS_POS]
    res = _sweep(
        ts,
        a_reg=[10, 0, 10],
        region_eff_len=[50, 0, 50],
        left=_view([0, 20, 20], spliced=[0, 100, 100]),
        right=_view([20, 20, 0], spliced=[100, 100, 0]),
        alloc_left=_alloc([0, 0, 0], [0.01, 0.01, 0.01]),
        alloc_right=_alloc([0, 0, 0], [0.01, 0.01, 0.01]),
        base_omega=[2.0, 1.0, 2.0],
    )
    np.testing.assert_allclose(res.omega[1], 1.0, atol=1e-3)


def test_no_ambig_is_a_noop():
    ts = [TS_POS, TS_POS]
    base = np.array([1.3, 0.7])
    res = _sweep(
        ts,
        a_reg=[10, 10],
        region_eff_len=[50, 50],
        left=_view([0, 20]),
        right=_view([20, 0]),
        alloc_left=_alloc([0, 0], [1.0, 1.0]),
        alloc_right=_alloc([0, 0], [1.0, 1.0]),
        base_omega=base,
    )
    assert res.n_ambig == 0
    np.testing.assert_array_equal(res.omega, base)
