"""Scenario: Overlapping genes on opposite strands (antisense).

g1 (+) and g2 (−) overlap spatially — the most challenging
scenario for strand-based disambiguation.
"""

import pytest
from rigel.sim import Scenario

from .conftest import (
    GDNA_LEVELS, STRAND_LEVELS, STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability,
    assert_negative_control, assert_gdna_accuracy,
)


class TestOverlappingAntisense:

    def _make_scenario(self, tmp_path, g1_abundance, g2_abundance,
                       name_suffix=""):
        sc = Scenario("overlap_anti" + name_suffix, genome_length=8000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / ("overlap_anti" + name_suffix))
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": g1_abundance},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
             "abundance": g2_abundance},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(5500, 5800)], "abundance": 0},
        ])
        return sc

    @pytest.mark.parametrize("fold_change", [1, 4, 16],
                             ids=["fc_1", "fc_4", "fc_16"])
    def test_abundance_sweep(self, tmp_path, fold_change):
        g1_abund, g2_abund = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, g1_abund, g2_abund)
        try:
            bench = build_and_run(sc, n_fragments=1000,
                                  scenario_name=f"anti_fc_{fold_change}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            max_overlap_loss = 0.35
            min_expected = bench.total_expected * (1.0 - max_overlap_loss)
            assert bench.total_rna_observed >= min_expected
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100, f"_gdna_{gdna}")
        try:
            bench = build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                  scenario_name=f"anti_gdna_{gdna}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                max_loss = 0.35
                min_exp = bench.total_expected * (1.0 - max_loss)
                assert bench.total_rna_observed >= min_exp
            else:
                assert_gdna_accuracy(bench, gdna, max_rel_err=1.0)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100, f"_ss_{ss}")
        try:
            bench = build_and_run(sc, strand_specificity=ss, n_fragments=1000,
                                  scenario_name=f"anti_ss_{ss}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, strand_specificity=ss)
            max_loss = (1.0 - ss) + 0.35
            min_exp = bench.total_expected * (1.0 - max_loss)
            assert bench.total_rna_observed >= min_exp
            if ss >= 0.95:
                for ta in bench.transcripts:
                    if ta.t_id == "t_ctrl":
                        continue
                    assert ta.rel_error < 0.40
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_g{gdna}_n{nrna}_s{int(ss*100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna,
                nrna_abundance=nrna, strand_specificity=ss,
                n_fragments=1000,
                scenario_name=f"anti_stress_{gdna}_{nrna}_{int(ss*100)}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna,
                                     strand_specificity=ss)
        finally:
            sc.cleanup()
