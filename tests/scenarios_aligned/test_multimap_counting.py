"""Multimapper counting scenarios (aligner-dependent).

These tests verify that multimapping fragments are counted ONCE per
physical molecule throughout the pipeline — never inflated by
the number of alignment hits.

Scenario overview:

1. **All-intergenic multimappers** — Duplicated intergenic regions
   cause gDNA/random reads to map to 2+ intergenic locations.
   Tests n_intergenic counting (once per molecule, not per hit).

2. **Intronic + intergenic multimappers** — A gene with large introns
   sits in a region that is duplicated elsewhere (intergenic).
   nRNA reads can map to both the gene intron AND the intergenic
   duplicate.  Tests mixed resolved/unresolved hit handling.

3. **Exonic paralogs** — Two genes with identical exonic sequence.
   All mRNA reads are NH=2.  Tests EM distribution.
   (Supplements test_paralogs.py with budget-focused assertions.)

4. **Pseudogene scenario** — A spliced gene with an intronless
   retrocopy elsewhere in the genome.  mRNA reads from the parent
   gene align to both the spliced parent (exonic, resolved) AND
   the processed pseudogene (single-exon, resolved).  Tests that
   the fragment budget is conserved and the EM resolves correctly.

5. **Mixed multimapper: exonic + intergenic** — One gene overlaps a
   region that is duplicated in an intergenic zone.  Some reads
   map to the gene (resolved) and to the intergenic copy
   (unresolved).  Tests that intergenic is NOT counted when
   a resolved hit exists.

Each test class includes:
- A baseline test with no gDNA, no nRNA
- A fragment budget conservation assertion (assert_accountability)
- Stress tests with combined gDNA + nRNA + low strand specificity
"""

import shutil

import pytest

from rigel.sim import Scenario

from .conftest import (
    STRESS_COMBOS,
    STRESS_IDS,
    SIM_SEED,
    build_and_run,
    assert_alignment,
    assert_accountability,
    assert_negative_control,
)

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
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(6000, 6200)], "abundance": 100},
        ])
        # Helper gene for strand training
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(14000, 14300), (14700, 15000)],
             "abundance": 80},
        ])
        # Negative control
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
        ])
        # Duplicate a 500bp intergenic region → multimapping gDNA
        sc.genome.edit(10000, sc.genome[2000:2500])
        return sc

    def test_baseline_no_gdna(self, tmp_path):
        """Without gDNA, no intergenic multimappers appear.
        Fragment budget must be exact."""
        sc = self._make_scenario(tmp_path, "_baseline")
        try:
            bench = build_and_run(
                sc, include_multimap=True,
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
                sc, gdna_abundance=50, n_fragments=1000,
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
                sc, gdna_abundance=200, n_fragments=2000,
                include_multimap=True,
                scenario_name="ig_mm_high_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, include_multimap=True,
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

    nRNA fragments from g1's intron can align to both:
      (a) the intronic region [1300-3000] → resolved to g1 nRNA
      (b) the intergenic copy [8000-9700] → unresolved (intergenic)

    With the bugfix, if ANY hit resolves, intergenic is NOT counted.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "intr_ig" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("intr_ig" + name_suffix),
        )
        # Gene with large intron → nRNA reads come from intron
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(1000, 1300), (3000, 3300)],
             "abundance": 100},
        ])
        # Helper for strand training
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(14000, 14300), (14700, 15000)],
             "abundance": 80},
        ])
        # Negative control
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
        ])
        # Duplicate g1's intron to an intergenic region
        intronic_seq = sc.genome[1300:3000]
        sc.genome.edit(8000, intronic_seq)
        return sc

    def test_nrna_budget(self, tmp_path):
        """nRNA fragments from g1 multimap to the intergenic copy.
        Fragment budget must be conserved."""
        sc = self._make_scenario(tmp_path, "_nrna")
        try:
            bench = build_and_run(
                sc, nrna_abundance=50, n_fragments=1000,
                include_multimap=True,
                scenario_name="intr_ig_nrna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_nrna_plus_gdna(self, tmp_path):
        """Combined nRNA + gDNA.  Both contribute multimappers
        to the intronic/intergenic duplicate."""
        sc = self._make_scenario(tmp_path, "_nrna_gdna")
        try:
            bench = build_and_run(
                sc, nrna_abundance=50, gdna_abundance=30,
                n_fragments=1500, include_multimap=True,
                scenario_name="intr_ig_nrna_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=1000,
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
    """Identical single-exon paralogs — all mRNA reads are NH=2.

    Same structure as TestParalogMultimapping but focused on
    fragment budget conservation rather than EM accuracy.

    Genome layout (16 kb):
        [500-1000]  gene g1 exon (+ strand)
        [5000-5500] gene g2 exon (identical to g1)
        [8000-8300][8700-9000] gene g_helper (+ strand, spliced)
        [12000-12300] gene g_ctrl (- strand, abundance=0)
    """

    def _make_scenario(self, tmp_path, g1_abund=100, g2_abund=100,
                       name_suffix=""):
        sc = Scenario(
            "para_budget" + name_suffix,
            genome_length=16000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("para_budget" + name_suffix),
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": g1_abund},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5500)], "abundance": g2_abund},
        ])
        # Copy g1 sequence to g2 → identical paralogs
        sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(8000, 8300), (8700, 9000)],
             "abundance": 80},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(12000, 12300)], "abundance": 0},
        ])
        return sc

    def test_budget_no_gdna(self, tmp_path):
        """Pure mRNA multimappers — budget must be exact."""
        sc = self._make_scenario(tmp_path, name_suffix="_budget_no_gdna")
        try:
            bench = build_and_run(
                sc, n_fragments=1000, include_multimap=True,
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
        """mRNA multimappers + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, name_suffix="_budget_gdna")
        try:
            bench = build_and_run(
                sc, gdna_abundance=50, n_fragments=1500,
                include_multimap=True,
                scenario_name="para_budget_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_budget_asymmetric(self, tmp_path):
        """Asymmetric paralog abundances (4:1) — budget conserved."""
        sc = self._make_scenario(tmp_path, g1_abund=100, g2_abund=25,
                                 name_suffix="_budget_asym")
        try:
            bench = build_and_run(
                sc, n_fragments=1000, include_multimap=True,
                scenario_name="para_budget_asym",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, name_suffix=f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=1000,
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
                               spliced mRNA of 600bp)
        [8000-8600]  pseudogene g_pseudo (+ strand, single exon,
                     sequence = concatenation of parent exons)
        [14000-14300][14700-15000] gene g_helper
        [17000-17300] gene g_ctrl (- strand, abundance=0)

    mRNA reads from g_parent span exon1+exon2 (spliced).  When
    aligned, they match:
      (a) g_parent at the spliced locus (with intron gap in CIGAR)
      (b) g_pseudo contiguously (no splice)

    This is the classic processed pseudogene scenario.  The aligner
    produces NH=2, and the EM must use splice evidence to prefer
    the parent gene.
    """

    def _make_scenario(self, tmp_path, parent_abund=100, pseudo_abund=0,
                       name_suffix=""):
        sc = Scenario(
            "pseudogene" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("pseudogene" + name_suffix),
        )
        # Parent gene: two exons with a 700bp intron
        sc.add_gene("g_parent", "+", [
            {"t_id": "t_parent",
             "exons": [(1000, 1300), (2000, 2300)],
             "abundance": parent_abund},
        ])
        # Pseudogene: single exon = concatenated parent exons
        sc.add_gene("g_pseudo", "+", [
            {"t_id": "t_pseudo",
             "exons": [(8000, 8600)],
             "abundance": pseudo_abund},
        ])
        # Copy concatenated exon sequences to pseudogene location
        exon1 = sc.genome[1000:1300]
        exon2 = sc.genome[2000:2300]
        sc.genome.edit(8000, exon1 + exon2)

        # Helper gene
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(14000, 14300), (14700, 15000)],
             "abundance": 80},
        ])
        # Negative control
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
        ])
        return sc

    def test_parent_only(self, tmp_path):
        """Only parent gene expressed — pseudogene silent.
        Budget must be conserved.  The EM may split counts between
        parent and pseudo since most reads don't span the splice
        junction, but the total should match expected."""
        sc = self._make_scenario(tmp_path, parent_abund=100, pseudo_abund=0,
                                 name_suffix="_parent_only")
        try:
            bench = build_and_run(
                sc, n_fragments=1000, include_multimap=True,
                scenario_name="pseudo_parent_only",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            # Total across parent+pseudo should match expected
            t_parent = next(
                t for t in bench.transcripts if t.t_id == "t_parent")
            t_pseudo = next(
                t for t in bench.transcripts if t.t_id == "t_pseudo")
            combined = t_parent.observed + t_pseudo.observed
            assert combined == pytest.approx(
                t_parent.expected, abs=20), (
                f"Combined parent+pseudo={combined:.0f} vs "
                f"expected={t_parent.expected}"
            )
        finally:
            sc.cleanup()

    def test_both_expressed(self, tmp_path):
        """Both parent and pseudogene expressed — budget conserved."""
        sc = self._make_scenario(tmp_path, parent_abund=100, pseudo_abund=50,
                                 name_suffix="_both")
        try:
            bench = build_and_run(
                sc, n_fragments=1500, include_multimap=True,
                scenario_name="pseudo_both",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    def test_parent_with_gdna(self, tmp_path):
        """Parent gene + gDNA contamination — budget conserved."""
        sc = self._make_scenario(tmp_path, parent_abund=100, pseudo_abund=0,
                                 name_suffix="_gdna")
        try:
            bench = build_and_run(
                sc, gdna_abundance=50, n_fragments=1500,
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
            tmp_path, parent_abund=100, pseudo_abund=0,
            name_suffix=f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=1000,
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

    mRNA reads from g1 exon 1 align to:
      (a) the gene itself (resolved → mRNA/EM)
      (b) the intergenic copy at 7000-7300 (unresolved → intergenic)

    With the bugfix, if hit (a) resolves, the intergenic hit (b) must
    NOT also increment n_intergenic.  The molecule enters the EM once.
    """

    def _make_scenario(self, tmp_path, name_suffix=""):
        sc = Scenario(
            "exon_ig" + name_suffix,
            genome_length=20000,
            seed=SIM_SEED,
            work_dir=tmp_path / ("exon_ig" + name_suffix),
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(1000, 1300), (2000, 2300)],
             "abundance": 100},
        ])
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(12000, 12300), (12700, 13000)],
             "abundance": 80},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(17000, 17300)], "abundance": 0},
        ])
        # Duplicate exon 1 to intergenic space
        sc.genome.edit(7000, sc.genome[1000:1300])
        return sc

    def test_budget_no_gdna(self, tmp_path):
        """mRNA reads multimap to exonic + intergenic.  Budget exact."""
        sc = self._make_scenario(tmp_path, "_budget")
        try:
            bench = build_and_run(
                sc, n_fragments=1000, include_multimap=True,
                scenario_name="exon_ig_budget",
            )
            assert_alignment(bench)
            assert_accountability(bench)
            assert_negative_control(bench)
            # g1 should get most of its expected counts despite multimap
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            assert t1.observed >= t1.expected * 0.50, (
                f"g1 under-quantified: {t1.observed:.0f} vs "
                f"expected {t1.expected}"
            )
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """Exonic+intergenic multimap + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, "_gdna")
        try:
            bench = build_and_run(
                sc, gdna_abundance=50, n_fragments=1500,
                include_multimap=True,
                scenario_name="exon_ig_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=1000,
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
    """Single-exon gene duplicated to TWO other locations → NH=3.

    Genome layout (24 kb):
        [1000-1500]  gene g1 exon (+ strand, abundance=100)
        [6000-6500]  gene g2 exon (identical to g1)
        [11000-11500] gene g3 exon (identical to g1)
        [16000-16300][16700-17000] gene g_helper
        [20000-20300] gene g_ctrl (- strand, abundance=0)

    Every mRNA read from g1/g2/g3 aligns to all three loci (NH=3).
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
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(1000, 1500)], "abundance": 100},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(6000, 6500)], "abundance": 100},
        ])
        sc.add_gene("g3", "+", [
            {"t_id": "t3", "exons": [(11000, 11500)], "abundance": 100},
        ])
        # Make all three loci identical
        seq = sc.genome[1000:1500]
        sc.genome.edit(6000, seq)
        sc.genome.edit(11000, seq)

        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(16000, 16300), (16700, 17000)],
             "abundance": 80},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(20000, 20300)], "abundance": 0},
        ])
        return sc

    def test_budget_three_way(self, tmp_path):
        """3-way multimapper — budget must be exact, not 3×."""
        sc = self._make_scenario(tmp_path, "_3way")
        try:
            bench = build_and_run(
                sc, n_fragments=1500, include_multimap=True,
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
                        f"{t.t_id}={t.observed:.0f} too far from "
                        f"mean={avg:.0f}"
                    )
        finally:
            sc.cleanup()

    def test_budget_with_gdna(self, tmp_path):
        """3-way + gDNA — budget conserved."""
        sc = self._make_scenario(tmp_path, "_3way_gdna")
        try:
            bench = build_and_run(
                sc, gdna_abundance=50, n_fragments=2000,
                include_multimap=True,
                scenario_name="triple_mm_3way_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=1500,
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
        sc.add_gene("g_spliced", "+", [
            {"t_id": "t_spliced",
             "exons": [(1000, 1300), (2000, 2300)],
             "abundance": 100},
        ])
        # Identical single-exon paralogs
        sc.add_gene("g_single_a", "+", [
            {"t_id": "t_single_a", "exons": [(8000, 8500)],
             "abundance": 80},
        ])
        sc.add_gene("g_single_b", "+", [
            {"t_id": "t_single_b", "exons": [(12000, 12500)],
             "abundance": 80},
        ])
        # Pseudogene (retrocopy of g_spliced)
        sc.add_gene("g_retro", "+", [
            {"t_id": "t_retro", "exons": [(16000, 16600)],
             "abundance": 0},
        ])
        # Helper
        sc.add_gene("g_helper", "+", [
            {"t_id": "t_helper",
             "exons": [(20000, 20300), (20700, 21000)],
             "abundance": 80},
        ])
        # Negative control
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(25000, 25300)], "abundance": 0},
        ])

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
                sc, n_fragments=2000, include_multimap=True,
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
                sc, gdna_abundance=50, n_fragments=2500,
                include_multimap=True,
                scenario_name="kitchen_sink_gdna",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", STRESS_COMBOS, ids=STRESS_IDS)
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(
            tmp_path, f"_g{gdna}_n{nrna}_s{int(ss * 100)}")
        try:
            bench = build_and_run(
                sc, gdna_abundance=gdna, nrna_abundance=nrna,
                strand_specificity=ss, n_fragments=2000,
                include_multimap=True,
                scenario_name=f"kitchen_sink_stress_{gdna}_{nrna}_{int(ss * 100)}",
            )
            assert_alignment(bench)
            assert_accountability(bench)
        finally:
            sc.cleanup()
