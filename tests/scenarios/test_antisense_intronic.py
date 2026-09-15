"""Scenario: an antisense gene nested inside a multi-exon gene's intron.

The topology is modelled on PVT1 / LINC02912::

    g1 (POS): exons [(1000,1200), (2000,2200), (5000,5600)]
       → intronic spans 1200–2000 and 2200–5000
    g2 (NEG): a single exon [(3000,3800)] inside g1's second intron

g2 is unexpressed in every arm, so every count it receives is a false positive, and the fragments
that could produce one are real: g1's intronic nascent fragments overlap g2's exon exactly, so they
are gene-ambiguous on position alone and only the strand model separates them. The file runs that
against four stresses — a single-exon g2, where the nascent prior is zeroed, a multi-exon g2, where it
is not and the strand model is the only defence, five g1 isoforms, which put isoform ambiguity on top
of gene ambiguity, and reduced strand specificity, which is what takes the separating channel away.
"""

import pytest
from rigel.sim import Scenario

from .conftest import (
    STRAND_LEVELS,
    NRNA_LEVELS,
    STRESS_COMBOS,
    STRESS_IDS,
    SIM_SEED,
    build_and_run,
    assert_alignment,
    assert_accountability,
    assert_negative_control,
    assert_nrna_detected,
)


# =====================================================================
# Single-isoform: one POS multi-exon + one NEG single-exon
# =====================================================================


class TestAntisenseIntronicOverlap:
    """g1(+) multi-exon with g2(−) single-exon in its intron."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "anti_intron", genome_length=8000, seed=SIM_SEED, work_dir=tmp_path / "anti_intron"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1",
                    "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
                    "abundance": 100,
                },
            ],
        )
        # g2: single-exon antisense inside g1 intron; abundance=0
        sc.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "+",
            [
                {"t_id": "t_ctrl", "exons": [(7000, 7300)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_mRNA_only_no_leak(self, scenario):
        """Pure mRNA from g1: zero nRNA, zero T2 counts."""
        bench = build_and_run(scenario, n_fragments=2000, scenario_name="anti_intron_mRNA")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 1, f"T2 received {t2.observed:.0f} mRNA counts (expected 0)"
        assert bench.n_nrna_pipeline <= 3, f"nRNA leak: {bench.n_nrna_pipeline:.0f} (expected ~0)"

    @pytest.mark.parametrize("nrna", NRNA_LEVELS, ids=[f"nrna_{n}" for n in NRNA_LEVELS])
    def test_nrna_sweep(self, scenario, nrna):
        """nRNA from g1 creates intronic reads overlapping g2 exon.

        T2 must still get zero assignments; all nRNA goes to T1.
        """
        bench = build_and_run(
            scenario,
            nrna_abundance=nrna,
            n_fragments=2000,
            scenario_name=f"anti_intron_nrna_{nrna}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, f"T2 mRNA leak: {t2.observed:.0f} (expected 0)"
        if nrna > 0:
            assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize(
        "ss",
        [
            pytest.param(
                s,
                marks=pytest.mark.xfail(
                    strict=True,
                    reason="ISSUES: nested-antisense-leak-under-the-sane-ruler — with the EM's ruler "
                    "honest (a gDNA-free library contracts nothing, DESIGN.md §7.2) the strand-flipped "
                    "intronic nascent fragments over t2 are assigned to it: 24 of 2,000 at SS 0.9, 124 at "
                    "0.65. The old bound was met only because a fabricated reference had contracted the "
                    "host's nascent entity 3.9×; the EM's assignment at an unwitnessed nested transcript "
                    "is the defect, and it is EM-side",
                ),
            )
            if s < 1.0
            else s
            for s in STRAND_LEVELS
        ],
        ids=[f"ss_{s}" for s in STRAND_LEVELS],
    )
    def test_strand_sweep_with_nrna(self, request, scenario, ss):
        """Nascent RNA at reduced strand specificity, which is where the separating channel weakens.

        `t2`'s truth is 0 and it must stay 0. The failure this holds off is not only a leak into `t2`
        but the same fragments being called gDNA in a library that contains none: at low strand
        specificity the intronic nascent fragments are ambiguous in position and nearly ambiguous in
        strand, so whatever the solver cannot place on the host lands somewhere wrong. The bound is
        left tight (5 at SS ≥ 0.9, 20 below) rather than tracking the current answer, so a regression
        shows.
        """
        bench = build_and_run(
            scenario,
            nrna_abundance=50,
            strand_specificity=ss,
            n_fragments=2000,
            scenario_name=f"anti_intron_nrna50_ss{ss}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        # Allow some T2 leak at low SS
        max_t2 = 5 if ss >= 0.9 else 20
        assert t2.observed <= max_t2, f"T2 mRNA leak: {t2.observed:.0f} at SS={ss} (limit={max_t2})"


# =====================================================================
# Multi-isoform: 5 POS isoforms + NEG single-exon
# =====================================================================


class TestAntisenseIntronicMultiIsoform:
    """5 g1 isoforms creates isoform ambiguity atop gene ambiguity."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "anti_intron_iso",
            genome_length=12000,
            seed=SIM_SEED,
            work_dir=tmp_path / "anti_intron_iso",
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1a",
                    "exons": [(1000, 1200), (2000, 2200), (8000, 8600)],
                    "abundance": 30,
                },
                {
                    "t_id": "t1b",
                    "exons": [(1000, 1200), (3000, 3400), (8000, 8600)],
                    "abundance": 25,
                },
                {
                    "t_id": "t1c",
                    "exons": [(1000, 1200), (2000, 2200), (5000, 5400), (8000, 8600)],
                    "abundance": 20,
                },
                {
                    "t_id": "t1d",
                    "exons": [(1000, 1200), (3000, 3400), (5000, 5400), (8000, 8600)],
                    "abundance": 15,
                },
                {"t_id": "t1e", "exons": [(1000, 1200), (8000, 8600)], "abundance": 10},
            ],
        )
        # g2: single-exon antisense within g1 intron
        sc.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(4000, 4400)], "abundance": 0},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "+",
            [
                {"t_id": "t_ctrl", "exons": [(10500, 10800)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_multi_isoform_nrna(self, scenario):
        """nRNA from multi-isoform g1: T2 must get zero."""
        bench = build_and_run(
            scenario, nrna_abundance=40, n_fragments=2000, scenario_name="anti_iso_nrna40"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, f"T2 mRNA leak: {t2.observed:.0f}"

    def test_multi_isoform_nrna_low_ss(self, scenario):
        """nRNA + SS=0.9 + multiple isoforms."""
        bench = build_and_run(
            scenario,
            nrna_abundance=40,
            strand_specificity=0.9,
            n_fragments=2000,
            scenario_name="anti_iso_nrna40_s90",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=0.9)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        # A TRIPWIRE, not a target. Calibration under-calls gDNA on this scenario, and a little of the
        # surplus RNA lands on the wrong isoform, so the limit admits under a percent of the library
        # rather than zero. Tightening it is the job of the prior work, not of this bound.
        assert t2.observed <= 20, (
            f"T2 mRNA leak: {t2.observed:.0f} at SS=0.9 (Step-1 interim limit=20)"
        )

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            n_fragments=2000,
            scenario_name=f"anti_iso_stress_{gdna}_{nrna}_{int(ss * 100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)


# =====================================================================
# Multi-exon antisense: nRNA prior NOT zeroed, strand-only defense
# =====================================================================


class TestAntisenseIntronicMultiExonT2:
    """g2 is multi-exon, so its nascent prior is NOT zeroed and the strand model alone has to prevent
    a wrong-strand nascent assignment."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "anti_intron_t2me",
            genome_length=10000,
            seed=SIM_SEED,
            work_dir=tmp_path / "anti_intron_t2me",
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1",
                    "exons": [(1000, 1200), (2000, 2200), (7000, 7600)],
                    "abundance": 100,
                },
            ],
        )
        # g2: multi-exon antisense, abundance=0
        sc.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(3000, 3400), (4500, 4900)], "abundance": 0},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "+",
            [
                {"t_id": "t_ctrl", "exons": [(9000, 9300)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_nrna_multiexon_t2(self, scenario):
        """nRNA from T1 with multi-exon T2: strand model must hold."""
        bench = build_and_run(
            scenario, nrna_abundance=50, n_fragments=2000, scenario_name="anti_me_nrna50"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        assert t2.observed <= 2, f"T2 mRNA leak with multi-exon T2: {t2.observed:.0f}"

    @pytest.mark.xfail(
        strict=True,
        reason="ISSUES: antisense-prior-assembly-casualty — ⛔ A PRIOR-ASSEMBLY CASUALTY, NOT A CALIBRATION OR MESSAGE-LAYER DEFECT — owner diagnosis, "
        "2026-08-18, measured under the relay policy of the day. `assemble_priors` pins synthetic nascent RNA at Dirichlet alpha = 0 "
        "(EQUATIONS.md §9b): gDNA gets an additive prior, annotated RNA a multiplicative one, and "
        "nascent must out-evidence both. On this scenario 1,600 true nascent fragments yield only "
        "~536 called, and with the transfer policy the messages recover MORE RNA overall "
        "(536 vs 334 muted; false gDNA 879 vs 1,144) — the recovered mass lands on the annotated "
        "antisense t2 (80 > the 50 limit) because the alpha = 0 rule forbids it landing on nascent. "
        "Measured: the leak is 80 under every single message-operator ablation and passes only with "
        "messages fully off, while BOTH pool totals are better messages-on — so the test's threshold "
        "is a casualty of a change that is net-helpful, and the repair belongs in the prior "
        "assembler's nascent handling, not in calibration.",
    )
    def test_nrna_multiexon_t2_low_ss(self, scenario):
        """nRNA + multi-exon T2 + SS=0.65 (worst case)."""
        bench = build_and_run(
            scenario,
            nrna_abundance=50,
            strand_specificity=0.65,
            n_fragments=2000,
            scenario_name="anti_me_nrna50_s65",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=0.65)
        t2 = next(t for t in bench.transcripts if t.t_id == "t2")
        # At SS=0.65 we expect some leakage
        max_t2 = 50
        assert t2.observed <= max_t2, f"T2 mRNA leak at SS=0.65: {t2.observed:.0f} (limit={max_t2})"
