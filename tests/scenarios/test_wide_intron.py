"""Scenario: Wide-intron fragment-length discrimination.

Exons (1000,2000) and (3000,4000) create a 1000 bp intron gap.
Validates that per-candidate fragment lengths correctly penalise
nRNA candidates for junction-spanning fragments, preventing false
nRNA allocation from spliced mRNA reads.
"""

import pytest
from hulkrna.sim import Scenario

from .conftest import (
    GDNA_LEVELS, SIM_SEED, build_and_run,
    assert_alignment, assert_accountability, assert_transcript_accuracy,
    assert_negative_control, assert_gdna_accuracy, assert_nrna_detected,
)


class TestWideIntronInsertPenalty:

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("wide_intron", genome_length=6000, seed=SIM_SEED,
                       work_dir=tmp_path / "wide_intron")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(1000, 2000), (3000, 4000)],
             "abundance": 100},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(5000, 5300)], "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_baseline_no_false_nrna(self, scenario):
        """Pure mRNA: wide intron must not produce false nRNA counts."""
        bench = build_and_run(scenario,
                              scenario_name="wide_intron_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=5)
        assert_negative_control(bench)
        assert bench.n_nrna_pipeline <= 3, (
            f"False nRNA detected: {bench.n_nrna_pipeline:.0f} "
            f"(expected <= 3 for pure mRNA scenario)"
        )

    @pytest.mark.parametrize("nrna", [30, 70],
                             ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_sweep(self, scenario, nrna):
        """With real nRNA present, pipeline detects it without siphoning mRNA."""
        bench = build_and_run(scenario, nrna_abundance=nrna,
                              scenario_name=f"wide_intron_nrna_{nrna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(scenario, gdna_abundance=gdna,
                              scenario_name=f"wide_intron_gdna_{gdna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=5)
        else:
            assert_gdna_accuracy(bench, gdna)
