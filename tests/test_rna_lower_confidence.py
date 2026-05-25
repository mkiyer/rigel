"""Phase 2 (strand-deconv v5): rna_lower_confidence config knob.

These tests assert that the single statistical user knob introduced by
``docs/fineregions/strand_model_impl_plan_v5.md`` Phase 2 is plumbed
end-to-end (config -> CLI -> orchestrator -> calibration result summary)
without affecting any downstream behaviour. Phase 2 is intentionally inert
beyond the validation layer; the value is still unused by the math.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.cli import (
    _PARAM_SPECS,
    _build_pipeline_config,
    _build_quant_defaults,
    _resolve_quant_args,
    build_parser,
)
from rigel.config import CalibrationConfig


# ---------------------------------------------------------------------------
# Defaults & validation on the dataclass
# ---------------------------------------------------------------------------


def test_rna_lower_confidence_default_is_095():
    cfg = CalibrationConfig()
    assert cfg.rna_lower_confidence == 0.95


@pytest.mark.parametrize("value", [0.5, 0.95, 0.999])
def test_rna_lower_confidence_accepts_in_range(value):
    cfg = CalibrationConfig(rna_lower_confidence=value)
    assert cfg.rna_lower_confidence == value


@pytest.mark.parametrize("value", [0.49, 0.0, -1.0, 1.0, 1.5, float("inf")])
def test_rna_lower_confidence_validation_rejects_out_of_range(value):
    with pytest.raises(ValueError, match="rna_lower_confidence"):
        CalibrationConfig(rna_lower_confidence=value)


# ---------------------------------------------------------------------------
# Orchestrator-level validation mirrors the config (defence in depth)
# ---------------------------------------------------------------------------


def test_orchestrator_rejects_out_of_range_rna_lower_confidence():
    from rigel.calibration._orchestrator import calibrate

    class _StubIndex:
        region_df = None

    with pytest.raises(ValueError, match="rna_lower_confidence"):
        calibrate(
            index=_StubIndex(),
            payload=None,  # type: ignore[arg-type]
            scan_trained=None,  # type: ignore[arg-type]
            rna_lower_confidence=1.0,
        )


# ---------------------------------------------------------------------------
# CLI round-trip
# ---------------------------------------------------------------------------


_QUANT_REQ = ["quant", "--bam", "x.bam", "--index", "idx", "-o", "out"]


def test_cli_rna_lower_confidence_default_is_none_before_resolution():
    parser = build_parser()
    args = parser.parse_args(_QUANT_REQ)
    assert args.rna_lower_confidence is None


def test_cli_rna_lower_confidence_round_trips_to_config():
    parser = build_parser()
    args = parser.parse_args([*_QUANT_REQ, "--rna-lower-confidence", "0.99"])
    assert args.rna_lower_confidence == pytest.approx(0.99)

    _resolve_quant_args(args, _build_quant_defaults())
    config = _build_pipeline_config(args, seed=0, sj_strand_tag="XS")
    assert config.calibration.rna_lower_confidence == pytest.approx(0.99)


def test_cli_rna_lower_confidence_default_resolves_to_095():
    parser = build_parser()
    args = parser.parse_args(_QUANT_REQ)
    _resolve_quant_args(args, _build_quant_defaults())
    config = _build_pipeline_config(args, seed=0, sj_strand_tag="XS")
    assert config.calibration.rna_lower_confidence == 0.95


def test_param_spec_registers_rna_lower_confidence():
    dests = {spec.cli_dest for spec in _PARAM_SPECS}
    assert "rna_lower_confidence" in dests
    targets = {spec.config_path for spec in _PARAM_SPECS}
    assert "calibration.rna_lower_confidence" in targets


# ---------------------------------------------------------------------------
# Summary echo
# ---------------------------------------------------------------------------


def test_calibration_result_summary_echoes_rna_lower_confidence(monkeypatch):
    """``to_summary_dict`` advertises the knob under ``calibration_config``."""
    from rigel.calibration._diagnostics import Diagnostics
    from rigel.calibration._result import CalibrationResult
    from rigel.calibration.density_model import DensityEvidence
    from rigel.calibration.exposure import RegionExposure
    from rigel.calibration.fl import FLModels
    from rigel.calibration.strand_deconv import RegionGdnaEstimate

    monkeypatch.setattr(FLModels, "to_summary_dict", lambda self: {}, raising=True)
    monkeypatch.setattr(Diagnostics, "to_summary_dict", lambda self: {}, raising=True)

    fl_models = FLModels.__new__(FLModels)
    diagnostics = Diagnostics.__new__(Diagnostics)
    density_evidence = DensityEvidence(
        rho_post=np.zeros(0, dtype=np.float64),
        relative_exposure=np.zeros(0, dtype=np.float64),
        mean_unbounded=np.zeros(0, dtype=np.float64),
        upper_unbounded=np.zeros(0, dtype=np.float64),
        prior_family=np.zeros(0, dtype=np.uint8),
        fallback_depth=np.zeros(0, dtype=np.uint8),
        flags=np.zeros(0, dtype=np.uint8),
        confidence=0.95,
        priors={},
        rho_ref=0.0,
        rho_ref_source="ZERO",
    )

    result = CalibrationResult(
        density_evidence=density_evidence,
        fl_models=fl_models,
        diagnostics=diagnostics,
        region_gdna=RegionGdnaEstimate(
            n_total=np.zeros(0, dtype=np.float32),
            mean_count=np.zeros(0, dtype=np.float32),
            upper_count=np.zeros(0, dtype=np.float32),
            rna_lower_count=np.zeros(0, dtype=np.float32),
            precision=np.zeros(0, dtype=np.float32),
            flags=np.zeros(0, dtype=np.uint8),
            kappa_d=2.0,
            kappa_d_n_seed_regions=0,
            kappa_d_n_exon_self_training=0,
            p_r1_sense=0.5,
            rna_lower_confidence=0.97,
        ),
        region_exposure=RegionExposure.uniform(0),
        rna_lower_confidence=0.97,
    )
    summary = result.to_summary_dict()
    assert "calibration_config" in summary
    assert summary["calibration_config"] == {"rna_lower_confidence": 0.97}
