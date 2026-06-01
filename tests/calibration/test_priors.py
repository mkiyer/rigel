"""assemble_priors — CalibrationResult → per-locus EM prior (PR 6 §I.2)."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _calibration(*, mass_d_cont, mass_d_left, mass_d_right, omega, rho_0=0.01) -> CalibrationResult:
    md_c = np.asarray(mass_d_cont, dtype=np.float64)
    n = md_c.shape[0]
    z = np.zeros(n, dtype=np.float64)
    return CalibrationResult(
        mass_g_contained=z.copy(),
        mass_d_contained=md_c,
        mass_g_left=z.copy(),
        mass_d_left=np.asarray(mass_d_left, dtype=np.float64),
        mass_g_right=z.copy(),
        mass_d_right=np.asarray(mass_d_right, dtype=np.float64),
        omega=np.asarray(omega, dtype=np.float64),
        log_omega_var=np.ones(n, dtype=np.float64),
        rho_0=rho_0,
        exposure_dispersion=0.2,
        rho_d_bb=0.01,
        kappa_rna=0.9,
        rho_r_bb=0.01,
        n_iterations=0,
        converged=True,
        mass_change_history=np.empty(0, dtype=np.float64),
        n_regions=n,
        config=CalibrationConfig(),
    )


def _regions(starts, ends) -> RegionArrays:
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    n = starts.shape[0]
    return RegionArrays(
        ref_id=np.zeros(n, dtype=np.int32),
        start=starts,
        end=ends,
        signature=np.zeros(n, dtype=np.uint8),
        ts_class=np.zeros(n, dtype=np.int8),
        region_size_bp=(ends - starts).astype(np.float64),
        ref_offsets=np.array([0, n], dtype=np.int32),
        order=np.arange(n, dtype=np.int64),
        n_refs=1,
    )


def _ml(locus_id, blocks) -> MultiLocus:
    loci = tuple(Locus(ref=str(rid), ref_id=rid, start=s, end=e) for rid, s, e in blocks)
    return MultiLocus(
        multi_locus_id=locus_id,
        transcript_indices=np.array([], dtype=np.int32),
        unit_indices=np.array([], dtype=np.int32),
        gdna_span=sum(e - s for _, s, e in blocks),
        loci=loci,
    )


def test_conservation_single_locus():
    # 3 regions fully inside one locus → φ = 1 each. D_r = cont+left+right,
    # e_r = ω·L_phys. The conservation identity (doc 04 §6.3, PR 6 §I.2):
    #   Σ(α_rna + α_gdna) = w·Σ_r [D_r + ρ_0·e_r].
    cal = _calibration(
        mass_d_cont=[2.0, 3.0, 4.0],
        mass_d_left=[0.5, 0.5, 0.5],
        mass_d_right=[0.5, 0.5, 0.5],
        omega=[1.0, 2.0, 1.5],
        rho_0=0.01,
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])], prior_weight=1.0)
    # D_r = [3, 4, 5] → 12; e_r = ω·100 = [100, 200, 150] → 450.
    np.testing.assert_allclose(priors.alpha_rna_add, [12.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [0.01 * 450.0])  # 4.5
    np.testing.assert_allclose(priors.gdna_eff_len, [450.0])
    np.testing.assert_allclose(
        priors.alpha_rna_add[0] + priors.alpha_gdna_add[0], 12.0 + 0.01 * 450.0
    )


def test_prior_weight_scales_alphas_not_eff_len():
    cal = _calibration(
        mass_d_cont=[10.0], mass_d_left=[0.0], mass_d_right=[0.0], omega=[1.0], rho_0=0.1
    )
    ra = _regions([0], [100])
    ml = [_ml(0, [(0, 0, 100)])]
    base = assemble_priors(cal, ra, ml, prior_weight=1.0)
    doubled = assemble_priors(cal, ra, ml, prior_weight=2.0)
    np.testing.assert_allclose(doubled.alpha_rna_add, 2.0 * base.alpha_rna_add)
    np.testing.assert_allclose(doubled.alpha_gdna_add, 2.0 * base.alpha_gdna_add)
    # gdna_eff_len is a geometric quantity — NOT scaled by prior_weight.
    np.testing.assert_allclose(doubled.gdna_eff_len, base.gdna_eff_len)
    np.testing.assert_allclose(base.gdna_eff_len, [100.0])


def test_region_split_between_two_loci():
    # One region [0,100) straddling two adjacent loci ([0,50), [50,100)) →
    # overlap shares 0.5/0.5; each locus gets half the region's mass.
    cal = _calibration(
        mass_d_cont=[10.0], mass_d_left=[0.0], mass_d_right=[0.0], omega=[1.0], rho_0=0.1
    )
    ra = _regions([0], [100])
    priors = assemble_priors(
        cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])], prior_weight=1.0
    )
    np.testing.assert_allclose(priors.alpha_rna_add, [5.0, 5.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [0.1 * 50.0, 0.1 * 50.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_zero_fragment_region_still_gives_gdna_prior():
    # No RNA mass, but the modeled gDNA prior ω·ρ_0·L is nonzero (doc 04 §6.4).
    cal = _calibration(
        mass_d_cont=[0.0], mass_d_left=[0.0], mass_d_right=[0.0], omega=[1.0], rho_0=0.02
    )
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
    np.testing.assert_allclose(priors.alpha_rna_add, [0.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [2.0])  # 0.02·100
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _calibration(
        mass_d_cont=[5.0, 99.0],
        mass_d_left=[0.0, 0.0],
        mass_d_right=[0.0, 0.0],
        omega=[1.0, 1.0],
        rho_0=0.01,
    )
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
    np.testing.assert_allclose(priors.alpha_rna_add, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.alpha_gdna_add, [0.01 * 100.0])


def test_empty_multiloci_returns_empty():
    cal = _calibration(mass_d_cont=[1.0], mass_d_left=[0.0], mass_d_right=[0.0], omega=[1.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [], prior_weight=1.0)
    assert priors.alpha_rna_add.shape == (0,)
    assert priors.alpha_gdna_add.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _calibration(
        mass_d_cont=[1.0, 2.0], mass_d_left=[0.0, 0.0], mass_d_right=[0.0, 0.0], omega=[1.0, 1.0]
    )
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
