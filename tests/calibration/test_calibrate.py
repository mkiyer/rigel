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


def test_exposure_and_density_sane():
    result = _run()
    assert np.isfinite(result.rho_0) and result.rho_0 >= 0.0
    for w in (result.omega_contained, result.omega_left, result.omega_right):
        assert np.all(np.isfinite(w))
        assert np.all(w >= 0.0)
    assert np.all(result.gdna_geom_len >= 0.0)
    assert 0.0 <= result.kappa_rna <= 1.0


def test_kappa_matches_strand_balance():
    # κ_rna is the posterior-predictive strand fit; the calibrator passes it through.
    from rigel.calibration.strand_balance import fit_strand_balance

    sb = fit_strand_balance(SimpleNamespace(p_r1_sense=0.95, n_observations=40))
    assert _run().kappa_rna == sb.kappa_rna


def test_confidence_knob_shifts_split_toward_gdna():
    # confidence is a posterior quantile: a higher quantile (0.9) reports a higher gDNA
    # fraction than the median (0.5), so the total deconvolved gDNA mass cannot decrease.
    def total_g(r):
        return float(r.mass_gdna_contained.sum() + r.mass_gdna_left.sum() + r.mass_gdna_right.sum())

    assert total_g(_run(CalibrationConfig(confidence=0.9))) >= total_g(
        _run(CalibrationConfig(confidence=0.5))
    ) - 1e-9
