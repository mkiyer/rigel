"""assemble_priors — acyclic CalibrationResult → per-locus EM prior."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.priors import assemble_priors
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig
from rigel.locus import Locus, MultiLocus


def _calibration(*, mass_g, mass_d, gdna_exposure_len, rho_0=0.01) -> CalibrationResult:
    """Build a result with all mass in the contained node (left/right = 0).

    ``mass_g`` / ``mass_d`` are the per-region deconvolved gDNA / RNA mass;
    ``gdna_exposure_len`` is the gDNA component's exposure-weighted length. ``rho_0``
    is a required field but is **not** read by ``assemble_priors`` (it is already
    baked into the deconvolved mass).
    """
    mg = np.asarray(mass_g, dtype=np.float64)
    md = np.asarray(mass_d, dtype=np.float64)
    el = np.asarray(gdna_exposure_len, dtype=np.float64)
    n = mg.shape[0]
    z = np.zeros(n, dtype=np.float64)
    o = np.ones(n, dtype=np.float64)
    return CalibrationResult(
        mass_g_contained=mg,
        mass_d_contained=md,
        mass_g_left=z.copy(),
        mass_d_left=z.copy(),
        mass_g_right=z.copy(),
        mass_d_right=z.copy(),
        omega_contained=o.copy(),  # not read by assemble_priors
        omega_left=o.copy(),
        omega_right=o.copy(),
        gdna_exposure_len=el,
        rho_0=rho_0,
        kappa_rna=0.9,
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


def test_single_locus_sums_region_nodes():
    # 3 regions fully inside one locus → φ = 1 each. alpha_gdna = Σ G_r,
    # alpha_rna = Σ D_r, gdna_eff_len = Σ E_r.
    cal = _calibration(
        mass_g=[1.0, 2.0, 1.5],
        mass_d=[3.0, 4.0, 5.0],
        gdna_exposure_len=[100.0, 200.0, 150.0],
    )
    ra = _regions([0, 100, 200], [100, 200, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 300)])], prior_weight=1.0)
    np.testing.assert_allclose(priors.alpha_rna_add, [12.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [4.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [450.0])


def test_prior_weight_scales_alphas_not_eff_len():
    cal = _calibration(mass_g=[10.0], mass_d=[10.0], gdna_exposure_len=[100.0])
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
    # overlap shares 0.5/0.5; each locus gets half of every projected quantity.
    cal = _calibration(mass_g=[5.0], mass_d=[10.0], gdna_exposure_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(
        cal, ra, [_ml(0, [(0, 0, 50)]), _ml(1, [(0, 50, 100)])], prior_weight=1.0
    )
    np.testing.assert_allclose(priors.alpha_rna_add, [5.0, 5.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [2.5, 2.5])
    np.testing.assert_allclose(priors.gdna_eff_len, [50.0, 50.0])


def test_evidence_free_region_gives_zero_gdna_prior():
    # Acyclic: no observed gDNA ⇒ zero gDNA pseudocount (no ρ_0·ω·L hallucination).
    # The gDNA eff-len floors at 1.0 (the EM's own default).
    cal = _calibration(mass_g=[0.0], mass_d=[0.0], gdna_exposure_len=[0.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
    np.testing.assert_allclose(priors.alpha_rna_add, [0.0])
    np.testing.assert_allclose(priors.alpha_gdna_add, [0.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [1.0])  # floored


def test_intergenic_region_dropped():
    # Region 1 ([200,300)) overlaps no locus → its mass is dropped, not allocated.
    cal = _calibration(
        mass_g=[1.0, 50.0],
        mass_d=[5.0, 99.0],
        gdna_exposure_len=[100.0, 100.0],
    )
    ra = _regions([0, 200], [100, 300])
    priors = assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
    np.testing.assert_allclose(priors.alpha_rna_add, [5.0])  # the intergenic 99 is gone
    np.testing.assert_allclose(priors.alpha_gdna_add, [1.0])
    np.testing.assert_allclose(priors.gdna_eff_len, [100.0])


def test_empty_multiloci_returns_empty():
    cal = _calibration(mass_g=[1.0], mass_d=[1.0], gdna_exposure_len=[100.0])
    ra = _regions([0], [100])
    priors = assemble_priors(cal, ra, [], prior_weight=1.0)
    assert priors.alpha_rna_add.shape == (0,)
    assert priors.alpha_gdna_add.shape == (0,)
    assert priors.gdna_eff_len.shape == (0,)


def test_region_count_mismatch_raises():
    cal = _calibration(mass_g=[1.0, 2.0], mass_d=[1.0, 2.0], gdna_exposure_len=[10.0, 20.0])
    ra = _regions([0], [100])  # 1 region vs calibration's 2
    with pytest.raises(ValueError):
        assemble_priors(cal, ra, [_ml(0, [(0, 0, 100)])], prior_weight=1.0)
