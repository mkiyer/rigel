"""Aligner-dependent scenarios: a multimapping fragment is counted once per physical molecule throughout the
pipeline, never once per alignment hit, and the EM resolves the paralog cases sensibly.

Each genome duplicates some region so reads map to two or more places: intergenic duplicates, an intronic region
with an intergenic copy, exonic paralogs identical in sequence, a processed pseudogene beside its spliced parent,
mixed exonic-plus-intergenic hits, a three-way hit, all of them in one genome, and paralogs made distinguishable by a
unique exon. Every class pairs its counting claim with a fragment-budget conservation assertion and a stress case
combining gDNA, nascent RNA and low strand specificity.
"""

import shutil
from functools import partial

import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario

import _scenarios
from _scenarios import (
    GDNA_LEVELS,
    PIPELINE_SEED,
    SIM_SEED,
    STRAND_LEVELS,
    STRESS_COMBOS,
    STRESS_IDS,
    assert_accountability,
)

#: Built through FASTQ and a real alignment. An aligner loses reads and leaves orientation ambiguity an oracle BAM does
#: not, so the alignment rate, the negative control's strand reference and the gDNA error bound are looser here.
build_and_run = partial(_scenarios.build_and_run, aligned=True)
assert_alignment = partial(_scenarios.assert_alignment, min_rate=0.70)
assert_negative_control = _scenarios.assert_negative_control
assert_gdna_accuracy = partial(_scenarios.assert_gdna_accuracy, max_rel_err=0.55)

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Scenario 1: All-intergenic multimappers
# =====================================================================


class TestIntergenicMultimappers:
    """Duplicated intergenic regions — gDNA reads multimap (NH=2).

    Genome layout (20 kb):
        [2000-2500]  intergenic duplicate A
        [6000-6200]  gene g1 (+ strand, single exon, abundance=100)
        [10000-10500] intergenic duplicate B (identical to A)
        [14000-14300][14700-15000] gene g_helper (+ strand, spliced)
        [17000-17300] gene g_ctrl (- strand, abundance=0)

    gDNA fragments landing in either duplicate region produce NH=2.
    The pipeline must count each such molecule exactly once in
    n_intergenic, not N times.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "ig_mm" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("ig_mm" + name_suffix),
        )
        # Real gene far from the duplicated regions
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(6000, 6200)], "abundance": 100},
            ],
        )
        # Helper gene for strand training
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(14000, 14300), (14700, 15000)], "abundance": 80},
            ],
        )
        # Negative control
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
            ],
        )
        # Duplicate a 500bp intergenic region → multimapping gDNA
        sc.genome.edit(10000, sc.genome[2000:2500])
        return sc

    def test_baseline_no_gdna(self, tmp_path):
        """Without gDNA, no intergenic multimappers appear.
        Fragment budget must be exact."""
        sc = self._make_scenario(tmp_path, "_baseline")
        try:
            bench = build_and_run(
                sc,
                include_multimap=True,
                scenario_name="ig_mm_baseline",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
        finally:
            sc.cleanup()

    def test_gdna_budget_conservation(self, tmp_path):
        """With gDNA, multimapping intergenic fragments must be
        counted once.  The total predicted must equal n_fragments."""
        sc = self._make_scenario(tmp_path, "_gdna_budget")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="ig_mm_gdna_budget",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=50)
        finally:
            sc.cleanup()

    def test_high_gdna_budget(self, tmp_path):
        """Heavy gDNA contamination — many intergenic multimappers."""
        sc = self._make_scenario(tmp_path, "_high_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=200,
                n_fragments=2000,
                include_multimap=True,
                scenario_name="ig_mm_high_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                include_multimap=True,
                scenario_name=f"ig_mm_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 2: Intronic + intergenic multimappers
# =====================================================================


class TestIntronicIntergenicMultimappers:
    """Gene with intronic region duplicated in intergenic space.

    Genome layout (20 kb):
        [1000-1300][3000-3300] gene g1 (+ strand, spliced, large intron)
                               intron spans [1300-3000] = 1700bp
        [8000-9700]  intergenic region identical to [1300-3000]
                     (duplicated intronic sequence)
        [14000-14300][14700-15000] gene g_helper (+ strand, spliced)
        [17000-17300] gene g_ctrl (- strand, abundance=0)

    Nascent fragments from g1's intron can align to both:
      (a) the intronic region [1300-3000] → resolved to g1's synthetic span
      (b) the intergenic copy [8000-9700] → unresolved (intergenic)

    If any hit resolves, the fragment is not counted intergenic.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "intr_ig" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("intr_ig" + name_suffix),
        )
        # Gene with large intron → nascent reads come from the intron
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(1000, 1300), (3000, 3300)], "abundance": 100},
            ],
        )
        # Helper for strand training
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(14000, 14300), (14700, 15000)], "abundance": 80},
            ],
        )
        # Negative control
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
            ],
        )
        # Duplicate g1's intron to an intergenic region
        intronic_seq = sc.genome[1300:3000]
        sc.genome.edit(8000, intronic_seq)
        return sc

    def test_nrna_budget(self, tmp_path):
        """Nascent fragments from g1 multimap to the intergenic copy.
        Fragment budget must be conserved."""
        sc = self._make_scenario(tmp_path, "_nrna")
        try:
            bench = build_and_run(
                sc,
                nrna_abundance=50,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="intr_ig_nrna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_nrna_plus_gdna(self, tmp_path):
        """Nascent RNA and gDNA combined.  Both contribute multimappers
        to the intronic/intergenic duplicate."""
        sc = self._make_scenario(tmp_path, "_nrna_gdna")
        try:
            bench = build_and_run(
                sc,
                nrna_abundance=50,
                gdna_abundance=30,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="intr_ig_nrna_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"intr_ig_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 3: Exonic paralogs with strict budget assertion
# =====================================================================


class TestExonicParalogBudget:
    """Identical single-exon paralogs: every RNA read is NH=2.

    Same structure as TestParalogMultimapping but focused on
    fragment budget conservation rather than EM accuracy.

    Genome layout (16 kb):
        [500-1000]  gene g1 exon (+ strand)
        [5000-5500] gene g2 exon (identical to g1)
        [8000-8300][8700-9000] gene g_helper (+ strand, spliced)
        [12000-12300] gene g_ctrl (- strand, abundance=0)
    """

    def _make_scenario(self, tmp_path, g1_abund=100, g2_abund=100, name_suffix=""):
        sc = Scenario(
            "para_budget" + name_suffix,
            genome_length=16000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("para_budget" + name_suffix),
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(500, 1000)], "abundance": g1_abund},
            ],
        )
        sc.add_gene(
            "g2",
            "+",
            [
                {"t_id": "t2", "exons": [(5000, 5500)], "abundance": g2_abund},
            ],
        )
        # Copy g1 sequence to g2 → identical paralogs
        sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 80},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(12000, 12300)], "abundance": 0},
            ],
        )
        return sc

    def test_budget_no_gdna(self, tmp_path):
        """RNA multimappers alone: budget must be exact."""
        sc = self._make_scenario(tmp_path, name_suffix="_budget_no_gdna")
        try:
            bench = build_and_run(
                sc,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="para_budget_no_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            # EM should split roughly 50/50
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            if total > 10:
                assert abs(t1.observed - t2.observed) < total * 0.20 + 5
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """RNA multimappers + gDNA: budget conserved."""
        sc = self._make_scenario(tmp_path, name_suffix="_budget_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="para_budget_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_budget_asymmetric(self, tmp_path):
        """Asymmetric paralog abundances (4:1) — budget conserved."""
        sc = self._make_scenario(tmp_path, g1_abund=100, g2_abund=25, name_suffix="_budget_asym")
        try:
            bench = build_and_run(
                sc,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="para_budget_asym",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, name_suffix=f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"para_budget_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 4: Pseudogene — spliced parent + intronless retrocopy
# =====================================================================


class TestPseudogene:
    """Spliced parent gene with processed pseudogene retrocopy.

    Genome layout (20 kb):
        [1000-1300][2000-2300] gene g_parent (+ strand, two exons,
                               spliced RNA of 600bp)
        [8000-8600]  pseudogene g_pseudo (+ strand, single exon,
                     sequence = concatenation of parent exons)
        [14000-14300][14700-15000] gene g_helper
        [17000-17300] gene g_ctrl (- strand, abundance=0)

    RNA reads from g_parent span exon1+exon2 (spliced).  When
    aligned, they match:
      (a) g_parent at the spliced locus (with intron gap in CIGAR)
      (b) g_pseudo contiguously (no splice)

    This is the classic processed pseudogene scenario.  The aligner
    produces NH=2, and the EM must use splice evidence to prefer
    the parent gene.
    """

    def _make_scenario(self, tmp_path, parent_abund=100, pseudo_abund=0, name_suffix=""):
        sc = Scenario(
            "pseudogene" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("pseudogene" + name_suffix),
        )
        # Parent gene: two exons with a 700bp intron
        sc.add_gene(
            "g_parent",
            "+",
            [
                {
                    "t_id": "t_parent",
                    "exons": [(1000, 1300), (2000, 2300)],
                    "abundance": parent_abund,
                },
            ],
        )
        # Pseudogene: single exon = concatenated parent exons
        sc.add_gene(
            "g_pseudo",
            "+",
            [
                {"t_id": "t_pseudo", "exons": [(8000, 8600)], "abundance": pseudo_abund},
            ],
        )
        # Copy concatenated exon sequences to pseudogene location
        exon1 = sc.genome[1000:1300]
        exon2 = sc.genome[2000:2300]
        sc.genome.edit(8000, exon1 + exon2)

        # Helper gene
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(14000, 14300), (14700, 15000)], "abundance": 80},
            ],
        )
        # Negative control
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
            ],
        )
        return sc

    def test_parent_only(self, tmp_path):
        """Only parent gene expressed — pseudogene silent.
        Budget must be conserved.  The EM may split counts between
        parent and pseudo since most reads don't span the splice
        sj, but the total should match expected."""
        sc = self._make_scenario(
            tmp_path, parent_abund=100, pseudo_abund=0, name_suffix="_parent_only"
        )
        try:
            bench = build_and_run(
                sc,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="pseudo_parent_only",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            # Total across parent+pseudo should match expected
            t_parent = next(t for t in bench.transcripts if t.t_id == "t_parent")
            t_pseudo = next(t for t in bench.transcripts if t.t_id == "t_pseudo")
            combined = t_parent.observed + t_pseudo.observed
            assert combined == pytest.approx(t_parent.expected, abs=20), (
                f"Combined parent+pseudo={combined:.0f} vs expected={t_parent.expected}"
            )
        finally:
            sc.cleanup()

    def test_both_expressed(self, tmp_path):
        """Both parent and pseudogene expressed — budget conserved."""
        sc = self._make_scenario(tmp_path, parent_abund=100, pseudo_abund=50, name_suffix="_both")
        try:
            bench = build_and_run(
                sc,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="pseudo_both",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_parent_with_gdna(self, tmp_path):
        """Parent gene + gDNA contamination — budget conserved."""
        sc = self._make_scenario(tmp_path, parent_abund=100, pseudo_abund=0, name_suffix="_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="pseudo_parent_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path,
            parent_abund=100,
            pseudo_abund=0,
            name_suffix=f"_g{gdna}_n{nrna}_s{int(ss * 100)}",
        )
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"pseudo_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 5: Exonic + intergenic multimappers
# =====================================================================


class TestExonicIntergenicMultimappers:
    """Gene exon duplicated in intergenic space.

    Genome layout (20 kb):
        [1000-1300][2000-2300] gene g1 (+ strand, spliced)
        [7000-7300]  intergenic duplicate of exon 1 [1000-1300]
        [12000-12300][12700-13000] gene g_helper
        [17000-17300] gene g_ctrl (- strand, abundance=0)

    RNA reads from g1 exon 1 align to:
      (a) the gene itself (resolved → the EM)
      (b) the intergenic copy at 7000-7300 (unresolved → intergenic)

    If hit (a) resolves, the intergenic hit (b) does not also
    increment n_intergenic.  The molecule enters the EM once.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "exon_ig" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("exon_ig" + name_suffix),
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(1000, 1300), (2000, 2300)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(12000, 12300), (12700, 13000)], "abundance": 80},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
            ],
        )
        # Duplicate exon 1 to intergenic space
        sc.genome.edit(7000, sc.genome[1000:1300])
        return sc

    def test_budget_no_gdna(self, tmp_path):
        """RNA reads multimap to exonic + intergenic.  Budget exact."""
        sc = self._make_scenario(tmp_path, "_budget")
        try:
            bench = build_and_run(
                sc,
                n_fragments=1000,
                include_multimap=True,
                scenario_name="exon_ig_budget",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            # g1 should get most of its expected counts despite multimap
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            assert t1.observed >= t1.expected * 0.50, (
                f"g1 under-quantified: {t1.observed:.0f} vs expected {t1.expected}"
            )
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """Exonic+intergenic multimap + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, "_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="exon_ig_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"exon_ig_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 6: Triple-hit multimapper (3 loci)
# =====================================================================


class TestTripleMultimapper:
    """Single-exon gene duplicated to two other locations → NH=3.

    Genome layout (24 kb):
        [1000-1500]  gene g1 exon (+ strand, abundance=100)
        [6000-6500]  gene g2 exon (identical to g1)
        [11000-11500] gene g3 exon (identical to g1)
        [16000-16300][16700-17000] gene g_helper
        [20000-20300] gene g_ctrl (- strand, abundance=0)

    Every RNA read from g1/g2/g3 aligns to all three loci (NH=3).
    The EM must handle 3-way multimappers and the budget must not
    inflate to 3× the true count.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "triple_mm" + name_suffix,
            genome_length=24000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("triple_mm" + name_suffix),
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(1000, 1500)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g2",
            "+",
            [
                {"t_id": "t2", "exons": [(6000, 6500)], "abundance": 100},
            ],
        )
        sc.add_gene(
            "g3",
            "+",
            [
                {"t_id": "t3", "exons": [(11000, 11500)], "abundance": 100},
            ],
        )
        # Make all three loci identical
        seq = sc.genome[1000:1500]
        sc.genome.edit(6000, seq)
        sc.genome.edit(11000, seq)

        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(16000, 16300), (16700, 17000)], "abundance": 80},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(20000, 20300)], "abundance": 0},
            ],
        )
        return sc

    def test_budget_three_way(self, tmp_path):
        """3-way multimapper — budget must be exact, not 3×."""
        sc = self._make_scenario(tmp_path, "_3way")
        try:
            bench = build_and_run(
                sc,
                n_fragments=1500,
                include_multimap=True,
                scenario_name="triple_mm_3way",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            # All three should get approximately equal counts
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            t3 = next(t for t in bench.transcripts if t.t_id == "t3")
            total = t1.observed + t2.observed + t3.observed
            if total > 20:
                avg = total / 3
                for t in [t1, t2, t3]:
                    assert abs(t.observed - avg) < avg * 0.30 + 5, (
                        f"{t.t_id}={t.observed:.0f} too far from mean={avg:.0f}"
                    )
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """3-way + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, "_3way_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=2000,
                include_multimap=True,
                scenario_name="triple_mm_3way_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1500,
                include_multimap=True,
                scenario_name=f"triple_mm_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 7: Mixed multimapper kitchen sink
# =====================================================================


class TestMixedMultimapperKitchenSink:
    """Multiple duplication events in a single genome.

    Genome layout (30 kb):
        [1000-1300][2000-2300]  gene g_spliced (+ strand, two exons)
        [5000-5300]   intergenic dup of exon 1 [1000-1300]
        [8000-8500]   gene g_single_a (+ strand, single exon)
        [12000-12500] gene g_single_b (+ strand, identical to g_single_a)
        [16000-16600] pseudogene g_retro (+ strand, single exon =
                      concat of g_spliced exons)
        [20000-20300][20700-21000] gene g_helper
        [25000-25300] gene g_ctrl (- strand, abundance=0)

    This creates multiple overlapping multimapper populations:
      - g_spliced exon1 reads → NH=2 (gene + intergenic dup)
      - g_single_a/b reads → NH=2 (identical paralogs)
      - g_spliced spliced reads → NH=2 (gene + pseudogene)
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "kitchen_sink" + name_suffix,
            genome_length=30000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("kitchen_sink" + name_suffix),
        )
        # Spliced gene
        sc.add_gene(
            "g_spliced",
            "+",
            [
                {"t_id": "t_spliced", "exons": [(1000, 1300), (2000, 2300)], "abundance": 100},
            ],
        )
        # Identical single-exon paralogs
        sc.add_gene(
            "g_single_a",
            "+",
            [
                {"t_id": "t_single_a", "exons": [(8000, 8500)], "abundance": 80},
            ],
        )
        sc.add_gene(
            "g_single_b",
            "+",
            [
                {"t_id": "t_single_b", "exons": [(12000, 12500)], "abundance": 80},
            ],
        )
        # Pseudogene (retrocopy of g_spliced)
        sc.add_gene(
            "g_retro",
            "+",
            [
                {"t_id": "t_retro", "exons": [(16000, 16600)], "abundance": 0},
            ],
        )
        # Helper
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(20000, 20300), (20700, 21000)], "abundance": 80},
            ],
        )
        # Negative control
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(25000, 25300)], "abundance": 0},
            ],
        )

        # Dup 1: copy exon 1 of g_spliced to intergenic
        sc.genome.edit(5000, sc.genome[1000:1300])
        # Dup 2: make g_single_b identical to g_single_a
        sc.genome.edit(12000, sc.genome[8000:8500])
        # Dup 3: pseudogene = concatenated exons of g_spliced
        exon1 = sc.genome[1000:1300]
        exon2 = sc.genome[2000:2300]
        sc.genome.edit(16000, exon1 + exon2)

        return sc

    def test_budget_baseline(self, tmp_path):
        """Multiple multimap classes in one genome — budget exact."""
        sc = self._make_scenario(tmp_path, "_baseline")
        try:
            bench = build_and_run(
                sc,
                n_fragments=2000,
                include_multimap=True,
                scenario_name="kitchen_sink_baseline",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """Kitchen sink + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, "_gdna")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=50,
                n_fragments=2500,
                include_multimap=True,
                scenario_name="kitchen_sink_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=2000,
                include_multimap=True,
                scenario_name=f"kitchen_sink_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_a_locus_holding_expressed_rna_is_not_handed_a_gdna_share_of_one(self, tmp_path):
        """The prior's gDNA share is calibration's gDNA count of a locus over calibration's count of every fragment
        it deposited there, the multimappers on both sides. Here nearly every fragment of the spliced gene and of the
        two paralogs multimaps, and gDNA multimaps with them: read over the locus's unique fragments alone, the count
        exceeds what it is divided by, the share clips at one, and loci whose genes are expressed get an RNA
        pseudocount of exactly zero under a prior as strong as their data."""
        sc = self._make_scenario(tmp_path, "_share")
        try:
            result = sc.build(
                n_fragments=3000,
                sim_config=_scenarios.sim_config(strand_specificity=1.0),
                gdna_config=_scenarios.gdna_config(20),
                nrna_abundance=0,
            )
            config = PipelineConfig(
                em=EMConfig(seed=PIPELINE_SEED),
                scan=BamScanConfig(sj_strand_tag="ts", include_multimap=True),
            )
            pr = run_pipeline(result.bam_path, result.index, config=config)
            t_ids = result.index.t_df["t_id"].tolist()
            locus_of = pr.estimator.locus_id_per_transcript
            by_locus = {int(r["locus_id"]): r for r in pr.estimator.locus_results}
            for t_id in ("t_spliced", "t_single_a", "t_single_b", "t_helper"):
                row = by_locus[int(locus_of[t_ids.index(t_id)])]
                strength = row["count_gdna_prior"] + row["count_rna_prior"]
                assert strength > 0.0, f"the locus of {t_id} has no prior"
                assert row["count_rna_prior"] > 0.0, (
                    f"the locus of {t_id}, an expressed gene, has an RNA pseudocount of 0 "
                    f"beside a gDNA pseudocount of {row['count_gdna_prior']:.1f}"
                )
        finally:
            sc.cleanup()


# ── paralogs with identical exonic sequence, every RNA read a multimapper ──


#: The identical-paralog gates run the EM to convergence. On two sequence-identical templates the split is the EM's
#: own degeneracy: it drifts from the even start toward a vertex, slowly, and at the shipped iteration budget it can
#: stop on the way, which is where the solver stops, not a property of the model. The shape is asserted where the
#: model decides it.
CONVERGED_EM = EMConfig(seed=PIPELINE_SEED, iterations=100_000, convergence_delta=1e-12)


def assert_identical_paralogs(t1, t2, *, tol: float) -> None:
    """What the EM can be held to on two sequence-identical templates: the total, and the shape.

    The two are not separable (every fragment inside them multimaps to both), so the split is an open
    EM degeneracy: bimodal rather than noisy, landing even or at a vertex, decided by the symmetry of the
    two loci's priors, not by depth. Calibration draws each multimapper to one placement, so the two
    priors differ by that draw's noise and the degenerate EM leaves the knife-edge for a vertex. The
    shape is asserted: an even split within ``tol`` or a clean vertex, and nothing in between. The
    tolerance must stay as it is, and no seed or depth may be moved until a mode comes out even, because
    that tunes the gate to green. `TestDistinguishableParalogs` is where the split is identifiable.
    """
    total = t1.observed + t2.observed
    if total <= 10:
        return
    even = abs(t1.observed - t2.observed) < total * tol + 5
    vertex = min(t1.observed, t2.observed) == 0.0 and max(t1.observed, t2.observed) == total
    assert even or vertex, (
        f"identical paralogs split {t1.observed} / {t2.observed}: neither even within {tol:.0%} nor a vertex"
    )


class TestParalogMultimapping:
    """Two genes with identical exonic sequences + negative control.

    All RNA reads are multimappers (NH=2).  The EM must distribute
    counts evenly when abundances are equal.
    """

    def _make_scenario(self, tmp_path, g1_abund, g2_abund, *, spliced=False, name_suffix=""):
        sc = Scenario(
            "paralogs" + name_suffix,
            genome_length=12000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("paralogs" + name_suffix),
        )
        if spliced:
            sc.add_gene(
                "g1",
                "+",
                [
                    {"t_id": "t1", "exons": [(500, 800), (1200, 1500)], "abundance": g1_abund},
                ],
            )
            sc.add_gene(
                "g2",
                "+",
                [
                    {"t_id": "t2", "exons": [(5000, 5300), (5700, 6000)], "abundance": g2_abund},
                ],
            )
            sc.genome.edit(5000, sc.genome[500:1500])
        else:
            sc.add_gene(
                "g1",
                "+",
                [
                    {"t_id": "t1", "exons": [(500, 1000)], "abundance": g1_abund},
                ],
            )
            sc.add_gene(
                "g2",
                "+",
                [
                    {"t_id": "t2", "exons": [(5000, 5500)], "abundance": g2_abund},
                ],
            )
            sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
            ],
        )
        return sc

    def test_equal_unspliced(self, tmp_path):
        """Equal-abundance unspliced paralogs -> ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100)
        try:
            bench = build_and_run(
                sc,
                include_multimap=True,
                scenario_name="paralogs_eq_unspliced",
                em_config=CONVERGED_EM,
            )
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=25)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            assert_identical_paralogs(t1, t2, tol=0.15)
        finally:
            sc.cleanup()

    def test_equal_spliced(self, tmp_path):
        """Equal-abundance spliced paralogs -> ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100, spliced=True)
        try:
            bench = build_and_run(
                sc,
                include_multimap=True,
                scenario_name="paralogs_eq_spliced",
                em_config=CONVERGED_EM,
            )
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=10)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            assert_identical_paralogs(t1, t2, tol=0.15)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100, name_suffix=f"_gdna_{gdna}")
        try:
            # n_fragments is 6× the default, so the EM has plenty of symmetric shared-multimapper
            # signal. That does not make the split between two identical paralogs reliable: the
            # collapse is bimodal, not noisy, and depth does not decide which way it goes; the total
            # is right either way. The tie is broken only by gDNA fragments extending into unique
            # flanking sequence, which are short and rare, so this is an EM degeneracy, not a
            # simulator defect (`assert_identical_paralogs`).
            bench = build_and_run(
                sc,
                n_fragments=3000,
                gdna_abundance=gdna,
                include_multimap=True,
                scenario_name=f"paralogs_gdna_{gdna}",
                em_config=CONVERGED_EM,
            )
            assert_alignment(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            # The total: with both paralog exons counted at their density by calibration's multimapper
            # second pass, the locus priors hold the gDNA and the EM's total runs near the truth instead
            # of absorbing the gDNA as RNA. It is a draw under `assignment_mode="sample"`, so it is held
            # to a band, not pinned.
            assert t1.observed + t2.observed == pytest.approx(
                t1.expected + t2.expected, abs=0.20 * (t1.expected + t2.expected) + 10
            )
            assert_identical_paralogs(t1, t2, tol=0.30 if gdna > 0 else 0.20)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS, ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100, name_suffix=f"_ss_{ss}")
        try:
            bench = build_and_run(
                sc,
                strand_specificity=ss,
                include_multimap=True,
                scenario_name=f"paralogs_ss_{ss}",
                em_config=CONVERGED_EM,
            )
            assert_alignment(bench)
            assert_negative_control(bench, strand_specificity=ss)
            if ss >= 0.95:
                assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=55)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            assert_identical_paralogs(t1, t2, tol=0.20)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, 100, 100, name_suffix=f"_g{gdna}_n{nrna}_s{int(ss * 100)}"
        )
        try:
            # 6× the default depth, matching the gDNA sweep above: at heavy gDNA on identical
            # paralogs the RNA pool is tiny and mostly multimapped, so the default depth leaves too
            # few unique spliced reads for the strand model to be identifiable at all.
            bench = build_and_run(
                sc,
                n_fragments=3000,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                include_multimap=True,
                scenario_name=f"paralogs_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)
        finally:
            sc.cleanup()


# ── paralogs made distinguishable by a unique exon ──────────────────


class TestDistinguishableParalogs:
    """Two genes sharing one exon but with distinct second exons + control.

    Reads from the shared exon are multimappers; reads from unique
    exons anchor the EM to correctly resolve shared-exon reads.
    """

    def _make_scenario(self, tmp_path, g1_abund, g2_abund, name_suffix=""):
        sc = Scenario(
            "dist_paralogs" + name_suffix,
            genome_length=12000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("dist_paralogs" + name_suffix),
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(500, 800), (1200, 1500)], "abundance": g1_abund},
            ],
        )
        sc.add_gene(
            "g2",
            "+",
            [
                {"t_id": "t2", "exons": [(5000, 5300), (5700, 5900)], "abundance": g2_abund},
            ],
        )
        sc.genome.edit(5000, sc.genome[500:800])
        # Spliced helper gene anchors calibration (provides splice signal
        # so the algebraic fallback doesn't misclassify RNA as gDNA).
        sc.add_gene(
            "g_helper",
            "+",
            [
                {"t_id": "t_helper", "exons": [(8000, 8300), (8700, 9000)], "abundance": 50},
            ],
        )
        sc.add_gene(
            "g_ctrl",
            "-",
            [
                {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
            ],
        )
        return sc

    @pytest.mark.parametrize("fold_change", [1, 4, 16], ids=["fc_1", "fc_4", "fc_16"])
    def test_abundance_sweep(self, tmp_path, fold_change):
        g1_abund, g2_abund = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, g1_abund, g2_abund, f"_fc_{fold_change}")
        try:
            bench = build_and_run(
                sc, n_fragments=1000, include_multimap=True, scenario_name=f"dist_fc_{fold_change}"
            )
            assert_alignment(bench)
            assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=15)
            if fold_change > 1:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS, ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100, f"_gdna_{gdna}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"dist_gdna_{gdna}",
            )
            assert_alignment(bench)
            assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(bench.total_expected, abs=10)
            else:
                assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc,
                gdna_abundance=gdna,
                nrna_abundance=nrna,
                strand_specificity=ss,
                n_fragments=1000,
                include_multimap=True,
                scenario_name=f"dist_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench, gdna_abundance=gdna, strand_specificity=ss)
        finally:
            sc.cleanup()
