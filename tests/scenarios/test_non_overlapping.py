"""Scenario: Two non-overlapping genes at distant loci.

g1 (+, spliced) and g2 (−, single-exon) — no spatial overlap,
no ambiguity between expressed genes.
"""

import pytest
from hulkrna.sim import Scenario

from .conftest import (
    GDNA_LEVELS, STRAND_LEVELS, STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability, assert_transcript_accuracy,
    assert_negative_control, assert_gdna_accuracy, assert_nrna_detected,
)


class TestNonOverlappingGenes:

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("non_overlapping", genome_length=10000, seed=SIM_SEED,
                       work_dir=tmp_path / "non_overlapping")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4400)], "abundance": 100},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(8000, 8300)], "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = build_and_run(scenario, scenario_name="nonoverlap_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=3)
        assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(scenario, gdna_abundance=gdna,
                              scenario_name=f"nonoverlap_gdna_{gdna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [30, 70],
                             ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_sweep(self, scenario, nrna):
        bench = build_and_run(scenario, nrna_abundance=nrna,
                              scenario_name=f"nonoverlap_nrna_{nrna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = build_and_run(scenario, strand_specificity=ss,
                              scenario_name=f"nonoverlap_ss_{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            assert_transcript_accuracy(bench, max_abs_diff=40)

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario, gdna_abundance=gdna,
            nrna_abundance=nrna, strand_specificity=ss,
            scenario_name=f"nonoverlap_g{gdna}_n{nrna}_s{int(ss*100)}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna,
                                 strand_specificity=ss)
