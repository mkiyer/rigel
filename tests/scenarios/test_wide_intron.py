"""Scenario: Wide-intron fragment-length discrimination.

Exons (1000,2000) and (3000,4000) create a 1000 bp intron gap.
Validates that per-candidate fragment lengths correctly penalise
nascent RNA candidates for sj-spanning fragments, preventing false
nascent RNA allocation from spliced annotated RNA reads.
"""

import pytest
from rigel.sim import Scenario

from _scenarios import (
    GDNA_LEVELS,
    SIM_SEED,
    build_and_run,
    assert_alignment,
    assert_accountability,
    assert_transcript_accuracy,
    assert_negative_control,
    assert_gdna_accuracy,
    assert_nascent_detected,
)


class TestWideIntronInsertPenalty:
    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "wide_intron", genome_length=6000, seed=SIM_SEED, work_dir=tmp_path / "wide_intron"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(1000, 2000), (3000, 4000)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(5000, 5300)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_baseline_no_false_nascent(self, scenario):
        """Pure annotated RNA: wide intron must not produce false nascent RNA counts."""
        bench = build_and_run(scenario, scenario_name="wide_intron_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=5)
        assert_negative_control(bench)
        assert bench.n_nascent_pipeline <= 3, (
            f"False nascent RNA detected: {bench.n_nascent_pipeline:.0f} "
            f"(expected <= 3 for pure annotated RNA scenario)"
        )

    @pytest.mark.parametrize("nascent", [30, 70], ids=[f"nascent_{n}" for n in [30, 70]])
    def test_nascent_sweep(self, scenario, nascent):
        """With real nascent RNA present, pipeline detects it without siphoning annotated RNA."""
        bench = build_and_run(
            scenario, nascent_abundance=nascent, scenario_name=f"wide_intron_nascent_{nascent}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        assert_nascent_detected(bench, nascent)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(
            scenario, gdna_abundance=gdna, scenario_name=f"wide_intron_gdna_{gdna}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=5)
        else:
            assert_gdna_accuracy(bench, gdna)
