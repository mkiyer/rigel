"""Paralog multimapping scenarios (aligner-dependent).

These scenarios depend on a real aligner (minimap2) to produce
multimapped reads (NH>1) from duplicated genome sequences.
The oracle BAM simulator currently writes NH=1 for all reads and
cannot model this behavior.

Scenario 6: Identical paralogs — all RNA reads are multimappers.
Scenario 7: Distinguishable paralogs — shared exon is multimapped,
            unique exons anchor the EM.
"""

import shutil

import pytest

from rigel.sim import Scenario

from .conftest import (
    GDNA_LEVELS, STRAND_LEVELS, STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability,
    assert_negative_control, assert_gdna_accuracy,
)

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Scenario 6: Paralog multimapping (identical sequences)
# =====================================================================


class TestParalogMultimapping:
    """Two genes with identical exonic sequences + negative control.

    All RNA reads are multimappers (NH=2).  The EM must distribute
    counts evenly when abundances are equal.
    """

    def _make_scenario(self, tmp_path, g1_abund, g2_abund, *,
                       spliced=False, name_suffix=""):
        sc = Scenario("paralogs" + name_suffix, genome_length=12000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / ("paralogs" + name_suffix))
        if spliced:
            sc.add_gene("g1", "+", [
                {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
                 "abundance": g1_abund},
            ])
            sc.add_gene("g2", "+", [
                {"t_id": "t2", "exons": [(5000, 5300), (5700, 6000)],
                 "abundance": g2_abund},
            ])
            sc.genome.edit(5000, sc.genome[500:1500])
        else:
            sc.add_gene("g1", "+", [
                {"t_id": "t1", "exons": [(500, 1000)],
                 "abundance": g1_abund},
            ])
            sc.add_gene("g2", "+", [
                {"t_id": "t2", "exons": [(5000, 5500)],
                 "abundance": g2_abund},
            ])
            sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(8000, 8300), (8700, 9000)],
             "abundance": 50},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
        ])
        return sc

    def test_equal_unspliced(self, tmp_path):
        """Equal-abundance unspliced paralogs -> ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100)
        try:
            bench = build_and_run(sc, include_multimap=True,
                                  scenario_name="paralogs_eq_unspliced")
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(
                bench.total_expected, abs=25)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            assert abs(t1.observed - t2.observed) < total * 0.15 + 5
        finally:
            sc.cleanup()

    def test_equal_spliced(self, tmp_path):
        """Equal-abundance spliced paralogs -> ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100, spliced=True)
        try:
            bench = build_and_run(sc, include_multimap=True,
                                  scenario_name="paralogs_eq_spliced")
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(
                bench.total_expected, abs=10)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            assert abs(t1.observed - t2.observed) < total * 0.15 + 5
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_gdna_{gdna}")
        try:
            # Use n_fragments=3000 (6× default) so the EM has enough
            # symmetric shared-multimapper signal to overcome asymmetric
            # gDNA noise at gene boundaries.  With identical paralogs,
            # gDNA fragments that extend into unique flanking sequence
            # map to only one paralog, creating stochastic warm-start
            # asymmetry that SQUAREM amplifies at small N.
            bench = build_and_run(sc, n_fragments=3000,
                                  gdna_abundance=gdna,
                                  include_multimap=True,
                                  scenario_name=f"paralogs_gdna_{gdna}")
            assert_alignment(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            if total > 10:
                tol = 0.30 if gdna > 0 else 0.20
                assert abs(t1.observed - t2.observed) < total * tol + 5
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_ss_{ss}")
        try:
            bench = build_and_run(sc, strand_specificity=ss,
                                  include_multimap=True,
                                  scenario_name=f"paralogs_ss_{ss}")
            assert_alignment(bench)
            assert_negative_control(bench, strand_specificity=ss)
            if ss >= 0.95:
                assert bench.total_rna_observed == pytest.approx(
                    bench.total_expected, abs=55)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            if total > 10:
                assert abs(t1.observed - t2.observed) < total * 0.20 + 5
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_g{gdna}_n{nrna}_s{int(ss*100)}")
        try:
            bench = build_and_run(sc, gdna_abundance=gdna,
                                  nrna_abundance=nrna, strand_specificity=ss,
                                  include_multimap=True,
                                  scenario_name=f"paralogs_stress_{gdna}_{nrna}_{int(ss*100)}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna,
                                     strand_specificity=ss)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 7: Distinguishable paralogs (shared + unique exons)
# =====================================================================


class TestDistinguishableParalogs:
    """Two genes sharing one exon but with distinct second exons + control.

    Reads from the shared exon are multimappers; reads from unique
    exons anchor the EM to correctly resolve shared-exon reads.
    """

    def _make_scenario(self, tmp_path, g1_abund, g2_abund,
                       name_suffix=""):
        sc = Scenario("dist_paralogs" + name_suffix, genome_length=12000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / ("dist_paralogs" + name_suffix))
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": g1_abund},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5300), (5700, 5900)],
             "abundance": g2_abund},
        ])
        sc.genome.edit(5000, sc.genome[500:800])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
        ])
        return sc

    @pytest.mark.parametrize("fold_change", [1, 4, 16],
                             ids=["fc_1", "fc_4", "fc_16"])
    def test_abundance_sweep(self, tmp_path, fold_change):
        g1_abund, g2_abund = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, g1_abund, g2_abund,
                                  f"_fc_{fold_change}")
        try:
            bench = build_and_run(sc, n_fragments=1000,
                                  include_multimap=True,
                                  scenario_name=f"dist_fc_{fold_change}")
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(
                bench.total_expected, abs=15)
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_gdna_{gdna}")
        try:
            bench = build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                  include_multimap=True,
                                  scenario_name=f"dist_gdna_{gdna}")
            assert_alignment(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(
                    bench.total_expected, abs=10)
            else:
                assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_g{gdna}_n{nrna}_s{int(ss*100)}")
        try:
            bench = build_and_run(sc, gdna_abundance=gdna,
                                  nrna_abundance=nrna, strand_specificity=ss,
                                  n_fragments=1000, include_multimap=True,
                                  scenario_name=f"dist_stress_{gdna}_{nrna}_{int(ss*100)}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna,
                                     strand_specificity=ss)
        finally:
            sc.cleanup()
