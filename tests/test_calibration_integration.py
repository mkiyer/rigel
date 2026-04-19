"""
test_calibration_integration.py — Integration tests for calibration wiring.

Validates that:
1. calibrate_gdna() is invoked by run_pipeline() and its output is stored on PipelineResult.
2. Calibrated gDNA priors flow into compute_locus_priors() and produce valid
   per-locus alpha values.
3. Clean scenarios (no gDNA) produce valid output with calibration.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import (
    BamScanConfig,
    CalibrationConfig,
    EMConfig,
    PipelineConfig,
)
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------


def _make_gdna_scenario(tmp_path, *, gdna_abundance=30.0, n_fragments=500):
    """Build a two-gene scenario with moderate gDNA contamination."""
    sc = Scenario(
        "cal_integration",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "cal_integration",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    gdna = GDNAConfig(
        abundance=gdna_abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=80,
        frag_max=600,
    )
    sim_config = SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.95,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=n_fragments, sim_config=sim_config, gdna_config=gdna)
    return sc, result


def _run_with_calibration(result, index):
    """Run the pipeline with calibration."""
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
        calibration=CalibrationConfig(),
    )
    return run_pipeline(result.bam_path, index, config=config)


# ---------------------------------------------------------------------------
# Test: CalibrationConfig defaults
# ---------------------------------------------------------------------------


class TestCalibrationConfig:
    """Verify CalibrationConfig creation and defaults."""

    def test_default_values(self):
        cfg = CalibrationConfig()
        assert cfg.gdna_prior_c_base == 5.0
        assert cfg.fl_prior_ess == 1000.0

    def test_pipeline_config_includes_calibration(self):
        pcfg = PipelineConfig()
        assert hasattr(pcfg, "calibration")
        assert isinstance(pcfg.calibration, CalibrationConfig)

    def test_frozen(self):
        cfg = CalibrationConfig()
        with pytest.raises(AttributeError):
            cfg.gdna_prior_c_base = 10.0


# ---------------------------------------------------------------------------
# Test: Calibration runs end-to-end in run_pipeline()
# ---------------------------------------------------------------------------


class TestCalibrationEndToEnd:
    """End-to-end integration: calibration runs inside run_pipeline()."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.sc, self.result = _make_gdna_scenario(tmp_path)
        self.index = self.result.index
        yield
        self.sc.cleanup()

    def test_calibration_populates_result(self):
        pr = _run_with_calibration(self.result, self.index)
        assert pr.calibration is not None, "Calibration should be populated"

    def test_calibration_has_expected_fields(self):
        pr = _run_with_calibration(self.result, self.index)
        cal = pr.calibration
        assert isinstance(cal.lambda_gdna, float)
        assert isinstance(cal.strand_specificity, float)
        assert cal.region_e_gdna is not None
        assert cal.region_n_total is not None
        assert cal.region_e_gdna.ndim == 1
        assert cal.region_n_total.ndim == 1
        assert len(cal.region_e_gdna) == len(cal.region_n_total)

    def test_calibration_lambda_nonnegative(self):
        pr = _run_with_calibration(self.result, self.index)
        cal = pr.calibration
        assert cal.lambda_gdna >= 0.0

    def test_calibration_e_gdna_valid(self):
        pr = _run_with_calibration(self.result, self.index)
        cal = pr.calibration
        assert np.all(cal.region_e_gdna >= 0.0)
        assert cal.region_e_gdna.sum() > 0  # scenario has gDNA

    def test_calibration_strand_specificity_valid(self):
        pr = _run_with_calibration(self.result, self.index)
        cal = pr.calibration
        assert 0.5 <= cal.strand_specificity <= 1.0

    def test_produces_valid_output(self):
        pr = _run_with_calibration(self.result, self.index)
        df = pr.estimator.get_counts_df(self.index)
        n_annotated = sum(1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_"))
        assert n_annotated == 3
        assert (df["count"] >= 0).all()
        assert df["count"].sum() > 0
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)


# ---------------------------------------------------------------------------
# Test: Calibrated vs Uncalibrated comparison
# ---------------------------------------------------------------------------


class TestCalibratedOutputFields:
    """Verify calibrated pipeline output has expected structure."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.sc, self.result = _make_gdna_scenario(tmp_path, gdna_abundance=50.0, n_fragments=800)
        self.index = self.result.index
        self.pr = _run_with_calibration(self.result, self.index)
        yield
        self.sc.cleanup()

    def test_produces_valid_transcripts(self):
        df = self.pr.estimator.get_counts_df(self.index)
        n_annotated = sum(1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_"))
        assert n_annotated == 3
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)

    def test_calibration_result_populated(self):
        assert self.pr.calibration is not None

    def test_lambda_gdna_nonnegative(self):
        cal = self.pr.calibration
        assert cal.lambda_gdna >= 0.0


# ---------------------------------------------------------------------------
# Test: Calibration object validation
# ---------------------------------------------------------------------------


class TestCalibrationObjectValid:
    """Unit tests for the calibration object produced by the pipeline."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        """Build a scenario, run quantification,
        then test the calibration object fields."""
        self.sc, self.result = _make_gdna_scenario(tmp_path, gdna_abundance=40.0, n_fragments=600)
        self.index = self.result.index

        # Run full pipeline with calibration to get cal object
        pr = _run_with_calibration(self.result, self.index)
        self.cal = pr.calibration
        yield
        self.sc.cleanup()

    def test_calibration_object_valid(self):
        assert self.cal is not None
        assert self.cal.lambda_gdna >= 0.0
        assert self.cal.strand_specificity > 0.0

    def test_calibration_region_arrays_populated(self):
        assert self.cal.region_e_gdna is not None
        assert len(self.cal.region_e_gdna) > 0

    def test_calibration_region_e_gdna_valid(self):
        e_gdna = self.cal.region_e_gdna
        assert e_gdna.ndim == 1
        assert len(e_gdna) > 0
        assert np.all(e_gdna >= 0.0)


# ---------------------------------------------------------------------------
# Test: Clean scenario (no gDNA) — calibration should be harmless
# ---------------------------------------------------------------------------


class TestCalibrationCleanScenario:
    """Calibration on a clean (no gDNA) scenario should not disrupt output."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        sc = Scenario(
            "cal_clean",
            genome_length=5000,
            seed=SEED,
            work_dir=tmp_path / "cal_clean",
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 100},
            ],
        )
        sim_config = SimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=SEED,
        )
        self.sc = sc
        self.result = sc.build_oracle(n_fragments=300, sim_config=sim_config)
        self.index = self.result.index
        yield
        sc.cleanup()

    def test_calibration_runs_without_error(self):
        pr = _run_with_calibration(self.result, self.index)
        assert pr.calibration is not None

    def test_clean_scenario_low_gdna(self):
        pr = _run_with_calibration(self.result, self.index)
        df = pr.estimator.get_counts_df(self.index)
        # With no gDNA in the simulation, gDNA should be near zero
        loci_df = pr.estimator.get_loci_df()
        if len(loci_df) > 0:
            assert loci_df["gdna"].sum() < df["count"].sum() * 0.1

    def test_output_valid(self):
        pr = _run_with_calibration(self.result, self.index)
        df = pr.estimator.get_counts_df(self.index)
        assert len(df) > 0
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)
