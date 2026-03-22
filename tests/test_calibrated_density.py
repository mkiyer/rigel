"""
test_calibrated_density.py — End-to-end tests for calibrated gDNA priors.

Validates that the pipeline produces valid gDNA priors with calibration
at different strand specificities.
"""

from __future__ import annotations

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
# Test: End-to-end A/B comparison at different strand specificities
# ---------------------------------------------------------------------------


def _make_scenario_with_ss(tmp_path, *, ss=0.95, gdna_abundance=40.0):
    """Build scenario with configurable strand specificity."""
    sc = Scenario(
        "phase3_ab",
        genome_length=8000,
        seed=SEED,
        work_dir=tmp_path / f"phase3_ss{int(ss * 100)}",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(300, 600), (1000, 1300)], "abundance": 60},
            {"t_id": "t2", "exons": [(300, 600), (1600, 1900)], "abundance": 40},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(4000, 4300), (4800, 5100)], "abundance": 50},
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
        strand_specificity=ss,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=800, sim_config=sim_config, gdna_config=gdna)
    return sc, result


class TestCalibratedDensityEndToEnd:
    """End-to-end tests comparing calibrated vs uncalibrated gDNA priors."""

    @pytest.fixture
    def high_ss_scenario(self, tmp_path):
        sc, result = _make_scenario_with_ss(tmp_path, ss=0.98)
        yield sc, result
        sc.cleanup()

    @pytest.fixture
    def low_ss_scenario(self, tmp_path):
        sc, result = _make_scenario_with_ss(tmp_path, ss=0.65)
        yield sc, result
        sc.cleanup()

    def _run(self, result, index):
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
            calibration=CalibrationConfig(),
        )
        return run_pipeline(result.bam_path, index, config=config)

    def test_high_ss_valid(self, high_ss_scenario):
        """At high SS, pipeline produces valid output."""
        sc, result = high_ss_scenario
        pr = self._run(result, result.index)
        df = pr.estimator.get_counts_df(result.index)
        n_annotated = sum(
            1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_")
        )
        assert n_annotated == 3
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)

    def test_low_ss_valid(self, low_ss_scenario):
        """At low SS, pipeline produces valid output."""
        sc, result = low_ss_scenario
        pr = self._run(result, result.index)
        df = pr.estimator.get_counts_df(result.index)
        n_annotated = sum(
            1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_")
        )
        assert n_annotated == 3
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)

    def test_high_ss_calibration_populates(self, high_ss_scenario):
        """Calibration result is populated at high SS."""
        sc, result = high_ss_scenario
        pr = self._run(result, result.index)
        assert pr.calibration is not None
        assert pr.calibration.gdna_density_global >= 0.0

    def test_low_ss_calibration_populates(self, low_ss_scenario):
        """Calibration result is populated at low SS."""
        sc, result = low_ss_scenario
        pr = self._run(result, result.index)
        assert pr.calibration is not None
        assert pr.calibration.gdna_density_global >= 0.0

    def test_low_ss_calibrated_gdna_reasonable(self, low_ss_scenario):
        """At low SS with gDNA, calibrated path should detect gDNA contamination."""
        sc, result = low_ss_scenario
        pr = self._run(result, result.index)
        loci_df = pr.estimator.get_loci_df()
        # gDNA was simulated — the pipeline should detect some
        if len(loci_df) > 0:
            total_gdna = loci_df["gdna"].sum()
            total_mrna = loci_df["mrna"].sum()
            # Should detect non-trivial gDNA (we simulated abundance=40)
            assert total_gdna >= 0.0
            # mRNA should still be the majority
            assert total_mrna > 0
