"""calibrate(): the acyclic single-pass calibrator — schema + structural invariants.

These pin the *mechanics* (a valid result, mass conservation per node, bounded masses,
sane exposure, κ_rna provenance, the confidence knob). Converged *biology* (paralog
rescue, exon→RNA) needs realistic data and is covered by the scenario suite.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from _synthetic import make_gdna_fl_pmf, make_synthetic_payload

from rigel.calibration import calibrate
from rigel.calibration.result import CalibrationResult
from rigel.config import CalibrationConfig


def _run(config=None):
    payload, ra = make_synthetic_payload()
    strand_model = SimpleNamespace(p_r1_sense=0.95, n_observations=40)
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=make_gdna_fl_pmf(),
        rna_fl_pmf=make_gdna_fl_pmf(),  # a valid pmf; these tests pin mechanics, not the splice fraction
        config=config or CalibrationConfig(),
    )


def test_returns_valid_result():
    result = _run()
    assert isinstance(result, CalibrationResult)
    assert result.n_regions == 3


def test_mass_conserved_per_node():
    # Each node's gDNA + RNA mass equals the node's total accumulator mass.
    result = _run()
    np.testing.assert_allclose(
        result.mass_gdna_contained + result.mass_rna_contained, [15.0, 26.0, 15.0]
    )
    np.testing.assert_allclose(result.mass_gdna_left + result.mass_rna_left, [0.0, 3.0, 1.5])
    np.testing.assert_allclose(result.mass_gdna_right + result.mass_rna_right, [2.0, 4.5, 0.0])


def test_masses_bounded():
    result = _run()
    for g, tot in (
        (result.mass_gdna_contained, np.array([15.0, 26.0, 15.0])),
        (result.mass_gdna_left, np.array([0.0, 3.0, 1.5])),
        (result.mass_gdna_right, np.array([2.0, 4.5, 0.0])),
    ):
        assert np.all(g >= -1e-9)
        assert np.all(g <= tot + 1e-9)


def test_gdna_strand_overdispersion_populated():
    # The fitted gDNA strand overdispersion is surfaced on the result and in range [0, 1).
    result = _run()
    assert 0.0 <= result.gdna_strand_overdispersion < 1.0


def test_rna_strand_overdispersion_populated():
    # The fitted RNA strand overdispersion is surfaced and clamped to the Beta(2,2) ceiling (0.2).
    result = _run()
    assert 0.0 <= result.rna_strand_overdispersion <= 0.2


def test_density_and_geom_len_sane():
    result = _run()
    assert np.isfinite(result.gdna_density_global) and result.gdna_density_global >= 0.0
    assert np.all(np.isfinite(result.gdna_geom_len))
    assert np.all(result.gdna_geom_len >= 0.0)
    assert 0.0 <= result.rna_sense_frac <= 1.0


def test_kappa_matches_strand_balance():
    # κ_rna is the posterior-predictive strand fit; the calibrator passes it through.
    from rigel.calibration.strand_balance import fit_strand_balance

    sb = fit_strand_balance(SimpleNamespace(p_r1_sense=0.95, n_observations=40))
    assert _run().rna_sense_frac == sb.rna_sense_frac


def _total_g(r):
    return float(r.mass_gdna_contained.sum() + r.mass_gdna_left.sum() + r.mass_gdna_right.sum())


def test_deconv_quantile_default_is_noop():
    # q=0.5 ⇒ Φ⁻¹=0 ⇒ no shift: bit-identical to the implicit default (the variance is consumed
    # only when q≠0.5). Pins the "default-preserving" property of the Phase-2 quantile knob.
    base = _run()
    half = _run(CalibrationConfig(gdna_deconv_quantile=0.5))
    np.testing.assert_array_equal(half.mass_gdna_contained, base.mass_gdna_contained)
    np.testing.assert_array_equal(half.mass_gdna_left, base.mass_gdna_left)
    np.testing.assert_array_equal(half.mass_gdna_right, base.mass_gdna_right)


def test_deconv_quantile_shifts_split_toward_gdna():
    # The FP-rate quantile g(q)=clip(center+Φ⁻¹(q)·σ) reports a higher posterior quantile of each
    # node's gDNA fraction as q rises: q>0.5 (FP-averse) cannot decrease the total deconvolved gDNA
    # mass vs neutral (0.5), and q<0.5 cannot increase it — monotone in q (widening, never sharpening).
    g_lo = _total_g(_run(CalibrationConfig(gdna_deconv_quantile=0.05)))
    g_mid = _total_g(_run(CalibrationConfig(gdna_deconv_quantile=0.5)))
    g_hi = _total_g(_run(CalibrationConfig(gdna_deconv_quantile=0.95)))
    assert g_lo <= g_mid + 1e-9
    assert g_mid <= g_hi + 1e-9
