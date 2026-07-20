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
    assert np.all(np.isfinite(result.gdna_region_eff_len))
    assert np.all(result.gdna_region_eff_len >= 0.0)
    assert np.all(np.isfinite(result.gdna_boundary_len))
    assert np.all(result.gdna_boundary_len >= 0.0)
    assert 0.0 <= result.rna_sense_frac <= 1.0


def test_kappa_matches_strand_balance():
    # κ_rna is the posterior-predictive strand fit; the calibrator passes it through.
    from rigel.calibration.strand_balance import fit_strand_balance

    sb = fit_strand_balance(SimpleNamespace(p_r1_sense=0.95, n_observations=40))
    assert _run().rna_sense_frac == sb.rna_sense_frac


def test_intron_factory_flag_runs_and_conserves_mass():
    # The gDNA intron factory flag is safe end-to-end: a valid result, mass still conserved per node.
    import dataclasses

    result = _run(dataclasses.replace(CalibrationConfig(), intron_factory=True))
    assert isinstance(result, CalibrationResult)
    np.testing.assert_allclose(
        result.mass_gdna_contained + result.mass_rna_contained, [15.0, 26.0, 15.0]
    )


def test_intron_factory_noop_without_introns():
    # Correct scoping: with no INTRON regions (the synthetic is +exon/−exon/intergenic) the factory is a
    # graceful no-op — byte-identical to the flag-off calibration.
    import dataclasses

    off = _run(CalibrationConfig())
    on = _run(dataclasses.replace(CalibrationConfig(), intron_factory=True))
    np.testing.assert_array_equal(off.mass_gdna_contained, on.mass_gdna_contained)
    np.testing.assert_array_equal(off.mass_rna_contained, on.mass_rna_contained)
