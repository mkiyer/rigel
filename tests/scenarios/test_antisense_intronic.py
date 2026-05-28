"""Scenario: Antisense single-exon gene embedded in multi-exon gene intron.

Topology (PVT1 / LINC02912 model):

    g1 (POS): exons [(1000,1200), (2000,2200), (5000,5600)]
       → intronic spans: 1200–2000 and 2200–5000
    g2 (NEG): single exon [(3000,3800)] inside g1's second intron

Key properties exercised:
 - g2 is single-exon → nRNA prior must be zeroed
 - Intronic g1 nRNA fragments overlap g2 exon → gene-ambiguous
 - Strand model must prevent wrong-strand nRNA assignment
 - Multiple isoform variants stress EM disambiguation
"""

import pytest
from rigel.sim import Scenario

from .conftest import (
    STRAND_LEVELS, NRNA_LEVELS,
    STRESS_COMBOS, STRESS_IDS,
    SIM_SEED, build_and_run,
    assert_alignment, assert_accountability, assert_negative_control, assert_nrna_detected,
)


# =====================================================================
# Single-isoform: one POS multi-exon + one NEG single-exon
# =====================================================================


class TestAntisenseIntronicOverlap:
    """g1(+) multi-exon with g2(−) single-exon in its intron."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("anti_intron", genome_length=8000, seed=SIM_SEED,
                       work_dir=tmp_path / "anti_intron")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
             "abundance": 100},
        ])
        # g2: single-exon antisense inside g1 intron; abundance=0
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_mRNA_only_no_leak(self, scenario):
        """Pure mRNA from g1: zero nRNA, zero T2 counts."""
        bench = build_and_run(scenario, n_fragments=2000,
                              scenario_name="anti_intron_mRNA")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 1, (
            f"T2 received {t2.observed:.0f} mRNA counts (expected 0)"
        )
        assert bench.n_nrna_pipeline <= 3, (
            f"nRNA leak: {bench.n_nrna_pipeline:.0f} (expected ~0)"
        )

    @pytest.mark.parametrize("nrna", NRNA_LEVELS,
                             ids=[f"nrna_{n}" for n in NRNA_LEVELS])
    def test_nrna_sweep(self, scenario, nrna):
        """nRNA from g1 creates intronic reads overlapping g2 exon.

        T2 must still get zero assignments; all nRNA goes to T1.
        """
        bench = build_and_run(scenario, nrna_abundance=nrna,
                              n_fragments=2000,
                              scenario_name=f"anti_intron_nrna_{nrna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, (
            f"T2 mRNA leak: {t2.observed:.0f} (expected 0)"
        )
        if nrna > 0:
            assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep_with_nrna(self, scenario, ss):
        """nRNA + reduced strand specificity."""
        bench = build_and_run(scenario, nrna_abundance=50,
                              strand_specificity=ss, n_fragments=2000,
                              scenario_name=f"anti_intron_nrna50_ss{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        # Allow some T2 leak at low SS
        max_t2 = 5 if ss >= 0.9 else 20
        assert t2.observed <= max_t2, (
            f"T2 mRNA leak: {t2.observed:.0f} at SS={ss} (limit={max_t2})"
        )


# =====================================================================
# Multi-isoform: 5 POS isoforms + NEG single-exon
# =====================================================================


class TestAntisenseIntronicMultiIsoform:
    """5 g1 isoforms creates isoform ambiguity atop gene ambiguity."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("anti_intron_iso", genome_length=12000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / "anti_intron_iso")
        sc.add_gene("g1", "+", [
            {"t_id": "t1a",
             "exons": [(1000, 1200), (2000, 2200), (8000, 8600)],
             "abundance": 30},
            {"t_id": "t1b",
             "exons": [(1000, 1200), (3000, 3400), (8000, 8600)],
             "abundance": 25},
            {"t_id": "t1c",
             "exons": [(1000, 1200), (2000, 2200), (5000, 5400), (8000, 8600)],
             "abundance": 20},
            {"t_id": "t1d",
             "exons": [(1000, 1200), (3000, 3400), (5000, 5400), (8000, 8600)],
             "abundance": 15},
            {"t_id": "t1e",
             "exons": [(1000, 1200), (8000, 8600)],
             "abundance": 10},
        ])
        # g2: single-exon antisense within g1 intron
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(4000, 4400)],
             "abundance": 0},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(10500, 10800)],
             "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_multi_isoform_nrna(self, scenario):
        """nRNA from multi-isoform g1: T2 must get zero."""
        bench = build_and_run(scenario, nrna_abundance=40,
                              n_fragments=2000,
                              scenario_name="anti_iso_nrna40")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, (
            f"T2 mRNA leak: {t2.observed:.0f}"
        )

    def test_multi_isoform_nrna_low_ss(self, scenario):
        """nRNA + SS=0.9 + multiple isoforms."""
        bench = build_and_run(scenario, nrna_abundance=40,
                              strand_specificity=0.9, n_fragments=2000,
                              scenario_name="anti_iso_nrna40_s90")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=0.9)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 10, (
            f"T2 mRNA leak: {t2.observed:.0f} at SS=0.9 (limit=10)"
        )

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS,
                             ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario, gdna_abundance=gdna,
            nrna_abundance=nrna, strand_specificity=ss,
            n_fragments=2000,
            scenario_name=f"anti_iso_stress_{gdna}_{nrna}_{int(ss*100)}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna,
                                 strand_specificity=ss)


# =====================================================================
# Multi-exon antisense: nRNA prior NOT zeroed, strand-only defense
# =====================================================================


class TestAntisenseIntronicMultiExonT2:
    """g2 is multi-exon → its nRNA prior is NOT zeroed.

    The strand model ALONE must prevent wrong-strand nRNA assignment.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("anti_intron_t2me", genome_length=10000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / "anti_intron_t2me")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(1000, 1200), (2000, 2200), (7000, 7600)],
             "abundance": 100},
        ])
        # g2: multi-exon antisense, abundance=0
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(3000, 3400), (4500, 4900)],
             "abundance": 0},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(9000, 9300)],
             "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_nrna_multiexon_t2(self, scenario):
        """nRNA from T1 with multi-exon T2: strand model must hold."""
        bench = build_and_run(scenario, nrna_abundance=50,
                              n_fragments=2000,
                              scenario_name="anti_me_nrna50")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, (
            f"T2 mRNA leak with multi-exon T2: {t2.observed:.0f}"
        )

    def test_nrna_multiexon_t2_low_ss(self, scenario):
        """nRNA + multi-exon T2 + SS=0.65 (worst case)."""
        bench = build_and_run(scenario, nrna_abundance=50,
                              strand_specificity=0.65, n_fragments=2000,
                              scenario_name="anti_me_nrna50_s65")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=0.65)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        # At SS=0.65 we expect some leakage
        max_t2 = 50
        assert t2.observed <= max_t2, (
            f"T2 mRNA leak at SS=0.65: {t2.observed:.0f} (limit={max_t2})"
        )
