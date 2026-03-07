"""Scenario: Single-exon unspliced gene.

The hardest case for gDNA separation — unspliced RNA is physically
identical to gDNA overlapping the gene region.  A separate multi-exon
helper gene provides splice junctions for strand-model training.
"""

import pytest
from rigel.sim import Scenario

from .conftest import (
    GDNA_LEVELS, STRAND_LEVELS, STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability, assert_transcript_accuracy,
    assert_negative_control, assert_gdna_accuracy,
)


class TestSingleExon:

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("single_exon", genome_length=8000, seed=SIM_SEED,
                       work_dir=tmp_path / "single_exon")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(2500, 3000), (3500, 4000)], "abundance": 50},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = build_and_run(scenario, scenario_name="single_exon_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=3)
        assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(scenario, gdna_abundance=gdna,
                              scenario_name=f"single_exon_gdna_{gdna}")
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
        """Nascent RNA on a single-exon gene is identical to mRNA."""
        bench = build_and_run(scenario, nrna_abundance=nrna,
                              scenario_name=f"single_exon_nrna_{nrna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = build_and_run(scenario, strand_specificity=ss,
                              scenario_name=f"single_exon_ss_{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            assert_transcript_accuracy(bench, max_abs_diff=55)

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario, gdna_abundance=gdna,
            nrna_abundance=nrna, strand_specificity=ss,
            scenario_name=f"single_exon_g{gdna}_n{nrna}_s{int(ss*100)}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna,
                                 strand_specificity=ss)
