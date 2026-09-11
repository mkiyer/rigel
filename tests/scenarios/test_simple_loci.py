"""End-to-end oracle scenarios for the simple loci: two isoforms of one gene, a
single-exon unspliced gene, a spliced two-exon gene, and two genes at distant
non-overlapping loci. Each block runs the whole pipeline on simulated data and
sweeps abundance ratio, gDNA contamination, nascent RNA and strand specificity,
checking alignment rate, fragment accountability, per-transcript accuracy and a
silent negative-control transcript.
"""

import pytest
from rigel.sim import Scenario

from .conftest import (
    GDNA_LEVELS,
    STRAND_LEVELS,
    STRESS_COMBOS,
    STRESS_IDS,
    SIM_SEED,
    build_and_run,
    assert_alignment,
    assert_accountability,
    assert_transcript_accuracy,
    assert_negative_control,
    assert_gdna_accuracy,
    assert_nrna_detected,
)


# ── two isoforms of a single gene ──────────────────────────────────────


class TestTwoIsoforms:
    """t1 (major) is three exons — shared_5p, middle, shared_3p; t2 (minor) is
    two exons — shared_5p, shared_3p — and so skips the middle exon."""

    ABUNDANCE_RATIOS = [1, 4, 16]

    def _make_scenario(self, tmp_path, major_abundance, minor_abundance, name_suffix=""):
        sc = Scenario(
            "two_isoforms" + name_suffix,
            genome_length=6000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("two_isoforms" + name_suffix),
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1",
                    "exons": [(200, 500), (1000, 1300), (2000, 2300)],
                    "abundance": major_abundance,
                },
                {"t_id": "t2", "exons": [(200, 500), (2000, 2300)], "abundance": minor_abundance},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(4500, 4800)], "abundance": 0},
            ],
        )
        return sc

    @pytest.mark.parametrize(
        "fold_change", ABUNDANCE_RATIOS, ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS]
    )
    def test_abundance_sweep(self, tmp_path, fold_change):
        major, minor = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, major, minor)
        try:
            bench = build_and_run(sc, n_fragments=1000, scenario_name=f"iso_fc_{fold_change}")
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=5)
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 10, f"_gdna_{gdna}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, n_fragments=1000, scenario_name=f"iso_gdna_{gdna}"
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=5)
            else:
                assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 10, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                scenario_name=f"iso_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)
        finally:
            sc.cleanup()


# ── a single-exon unspliced gene ───────────────────────────────────────


class TestSingleExon:
    """The hardest case for gDNA separation: unspliced RNA is physically
    identical to gDNA overlapping the gene region. A separate multi-exon helper
    gene provides the splice junctions the strand model trains on."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "single_exon", genome_length=8000, seed=SIM_SEED, work_dir=tmp_path / "single_exon"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(2500, 3000), (3500, 4000)], "abundance": 50},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(6500, 6800)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = build_and_run(scenario, scenario_name="single_exon_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=10)
        assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(
            scenario, gdna_abundance=gdna, scenario_name=f"single_exon_gdna_{gdna}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=10)
        else:
            assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [30, 70], ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_sweep(self, scenario, nrna):
        """Nascent RNA on a single-exon gene is identical to mRNA."""
        bench = build_and_run(
            scenario, nrna_abundance=nrna, scenario_name=f"single_exon_nrna_{nrna}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)

    @pytest.mark.parametrize("ss", STRAND_LEVELS, ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = build_and_run(scenario, strand_specificity=ss, scenario_name=f"single_exon_ss_{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            assert_transcript_accuracy(bench, max_abs_diff=55)

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"single_exon_g{gdna}_n{nrna}_s{int(ss * 100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)


# ── a spliced two-exon gene ────────────────────────────────────────────


class TestSplicedGene:
    """Splice junctions provide strong mRNA anchors, making gDNA/RNA separation
    easier than the unspliced case."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "spliced_gene", genome_length=5000, seed=SIM_SEED, work_dir=tmp_path / "spliced_gene"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(3500, 3800)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = build_and_run(scenario, scenario_name="spliced_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=3)
        assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(scenario, gdna_abundance=gdna, scenario_name=f"spliced_gdna_{gdna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [30, 70], ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_sweep(self, scenario, nrna):
        """Spliced gene: nRNA produces intronic reads detectable by pipeline."""
        bench = build_and_run(scenario, nrna_abundance=nrna, scenario_name=f"spliced_nrna_{nrna}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS, ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = build_and_run(scenario, strand_specificity=ss, scenario_name=f"spliced_ss_{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        max_loss = (1.0 - ss) + 0.20
        t1 = next(t for t in bench.transcripts if t.t_id == "t1")
        assert t1.observed >= t1.expected * (1.0 - max_loss), (
            f"t1: {t1.observed:.0f} too low vs {t1.expected} at ss={ss}"
        )

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"spliced_g{gdna}_n{nrna}_s{int(ss * 100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)


# ── two non-overlapping genes at distant loci ──────────────────────────


class TestNonOverlappingGenes:
    """g1 (+, spliced) and g2 (−, single-exon) have no spatial overlap, so there
    is no ambiguity between the expressed genes."""

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "non_overlapping",
            genome_length=10000,
            seed=SIM_SEED,
            work_dir=tmp_path / "non_overlapping",
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g2",
            "-",
            [
                {"t_id": "t2", "exons": [(4000, 4400)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "+",
            [
                {"t_id": "t_ctrl", "exons": [(8000, 8300)], "abundance": 0},
            ],
        )
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = build_and_run(scenario, scenario_name="nonoverlap_base")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_transcript_accuracy(bench, max_abs_diff=3)
        assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = build_and_run(
            scenario, gdna_abundance=gdna, scenario_name=f"nonoverlap_gdna_{gdna}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [30, 70], ids=[f"nrna_{n}" for n in [30, 70]])
    def test_nrna_sweep(self, scenario, nrna):
        bench = build_and_run(
            scenario, nrna_abundance=nrna, scenario_name=f"nonoverlap_nrna_{nrna}"
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench)
        assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS, ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = build_and_run(scenario, strand_specificity=ss, scenario_name=f"nonoverlap_ss_{ss}")
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            assert_transcript_accuracy(bench, max_abs_diff=40)

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = build_and_run(
            scenario,
            gdna_abundance=gdna,
            nrna_abundance=nrna,
            strand_specificity=ss,
            scenario_name=f"nonoverlap_g{gdna}_n{nrna}_s{int(ss * 100)}",
        )
        assert_alignment(bench)
        assert_accountability(bench)
        assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)
