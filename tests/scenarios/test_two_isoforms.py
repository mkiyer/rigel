"""Scenario: Two isoforms of a single gene.

t1 (major): 3 exons — shared_5p, middle, shared_3p.
t2 (minor): 2 exons — shared_5p, shared_3p (skips middle exon).
"""

import pytest
from hulkrna.sim import Scenario

from .conftest import (
    GDNA_LEVELS, STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability,
    assert_negative_control, assert_gdna_accuracy,
)


class TestTwoIsoforms:

    ABUNDANCE_RATIOS = [1, 4, 16]

    def _make_scenario(self, tmp_path, major_abundance, minor_abundance,
                       name_suffix=""):
        sc = Scenario("two_isoforms" + name_suffix, genome_length=6000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / ("two_isoforms" + name_suffix))
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": major_abundance},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300)],
             "abundance": minor_abundance},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(4500, 4800)], "abundance": 0},
        ])
        return sc

    @pytest.mark.parametrize("fold_change", ABUNDANCE_RATIOS,
                             ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS])
    def test_abundance_sweep(self, tmp_path, fold_change):
        major, minor = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, major, minor)
        try:
            bench = build_and_run(sc, n_fragments=1000,
                                  scenario_name=f"iso_fc_{fold_change}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(
                bench.total_expected, abs=5)
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 10, f"_gdna_{gdna}")
        try:
            bench = build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                  scenario_name=f"iso_gdna_{gdna}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(
                    bench.total_expected, abs=5)
            else:
                assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 10,
                                  f"_g{gdna}_n{nrna}_s{int(ss*100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna,
                nrna_abundance=nrna, strand_specificity=ss,
                n_fragments=1000,
                scenario_name=f"iso_stress_{gdna}_{nrna}_{int(ss*100)}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna,
                                     strand_specificity=ss)
        finally:
            sc.cleanup()
