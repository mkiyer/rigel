"""
End-to-end integration tests for hulkrna's 3-pool EM pipeline.

Every scenario includes a negative control transcript (t_ctrl, abundance=0)
to track false positives.  Parameter sweeps cover gDNA contamination,
nascent RNA fraction, and strand specificity (SS > 0.6 required).

Tests sweep three independent axes:
  1. gDNA contamination: [0, 5, 20, 50, 100]
  2. Nascent RNA fraction: [0.0, 0.1, 0.3, 0.5, 0.7]
  3. Strand specificity: [0.65, 0.8, 0.9, 0.95, 1.0]

Requirements: minimap2 and samtools must be available in PATH.
"""

import logging
import shutil

import numpy as np
import pytest

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark


logger = logging.getLogger(__name__)

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Parameter grids
# =====================================================================

GDNA_LEVELS = [0, 5, 20, 50, 100]
STRAND_LEVELS = [0.65, 0.8, 0.9, 0.95, 1.0]
NRNA_FRACTIONS = [0.0, 0.1, 0.3, 0.5, 0.7]
ABUNDANCE_RATIOS = [1, 2, 4, 8, 16]

N_FRAGMENTS = 500
SIM_SEED = 42
PIPELINE_SEED = 42


# =====================================================================
# Helpers
# =====================================================================


def _sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity, seed=seed,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


def _build_and_run(scenario, *, n_fragments=N_FRAGMENTS,
                   gdna_abundance=0, strand_specificity=1.0,
                   nrna_fraction=0.0, include_multimap=False,
                   scenario_name=""):
    sim_config = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)
    result = scenario.build(
        n_fragments=n_fragments, sim_config=sim_config,
        gdna_config=gdna, nrna_fraction=nrna_fraction,
    )
    pr = run_pipeline(
        result.bam_path, result.index,
        sj_strand_tag="ts", seed=PIPELINE_SEED,
        include_multimap=include_multimap,
    )
    bench = run_benchmark(result, pr, scenario_name=scenario_name)
    logger.info("\n%s", bench.summary())
    return bench


def _assert_alignment(bench, min_rate=0.70):
    assert bench.n_fragments > 0, "No fragments entered the pipeline"
    assert bench.alignment_rate > min_rate, (
        f"Low alignment rate: {bench.alignment_rate:.1%}"
    )


def _assert_accountability(bench, tolerance=5):
    """Every fragment must be accounted for (mRNA + nRNA + gDNA + chimeric).

    Multimapper alignments are counted individually in n_fragments but
    the EM distributes 1.0 per read *group*, so we must subtract the
    extra duplicate alignments.
    """
    sd = bench.stats_dict
    mm_extra = sd.get("n_multimapper_alignments", 0) - sd.get(
        "n_multimapper_groups", 0
    )
    n_gated_out = sd.get("n_gated_out", 0)
    effective_fragments = bench.n_fragments - mm_extra
    total = (
        bench.total_observed + bench.n_nrna_pipeline
        + bench.n_gdna_pipeline + bench.n_chimeric
        + n_gated_out
    )
    assert abs(total - effective_fragments) <= tolerance, (
        f"Accountability gap: {abs(total - effective_fragments):.0f} "
        f"(total={total:.0f}, eff_frags={effective_fragments}, "
        f"mm_extra={mm_extra})"
    )


def _assert_transcript_accuracy(bench, max_abs_diff=3,
                                 exclude_ids=("t_ctrl", "t_helper")):
    """Per-transcript mRNA accuracy (skips controls/helpers)."""
    for ta in bench.transcripts:
        if ta.t_id in exclude_ids:
            continue
        assert ta.abs_diff <= max_abs_diff, (
            f"{ta.t_id}: expected={ta.expected}, "
            f"observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
        )


def _assert_negative_control(bench, ctrl_id="t_ctrl", *,
                              gdna_abundance=0,
                              strand_specificity=1.0):
    """Assert negative-control transcript has near-zero mRNA counts.

    Tolerance scales with gDNA contamination (gDNA landing on the
    control gene may leak past the shadow) and low strand specificity
    (weakens strand-based discrimination).
    """
    ctrl = next(t for t in bench.transcripts if t.t_id == ctrl_id)
    max_fp = 5
    if gdna_abundance > 0:
        max_fp += min(gdna_abundance, 60)
    if strand_specificity < 0.9:
        ss_gap = 0.9 - strand_specificity
        # SS penalty + gDNA×SS interaction (low SS amplifies leakage)
        max_fp += round(ss_gap * 200 + gdna_abundance * ss_gap * 8)
    assert ctrl.observed <= max_fp, (
        f"Negative control {ctrl_id}: {ctrl.observed:.0f} counts "
        f"(limit={max_fp}, gdna={gdna_abundance}, "
        f"ss={strand_specificity})"
    )


def _assert_gdna_accuracy(bench, gdna_abundance, max_rel_err=0.55):
    if gdna_abundance == 0:
        return
    for ta in bench.transcripts:
        if ta.t_id in ("t_ctrl", "t_helper"):
            continue
        assert ta.observed <= ta.expected + bench.n_gdna_expected + 5, (
            f"{ta.t_id}: observed={ta.observed:.0f} exceeds "
            f"RNA({ta.expected}) + gDNA({bench.n_gdna_expected})"
        )
    if bench.n_gdna_expected > 10:
        rel_err = bench.gdna_abs_diff / bench.n_gdna_expected
        assert rel_err < max_rel_err, (
            f"gDNA rel error: pipeline={bench.n_gdna_pipeline:.0f}, "
            f"expected={bench.n_gdna_expected}, rel={rel_err:.2f}"
        )


def _assert_nrna_detected(bench, nrna_fraction, max_rel_err=0.80):
    """When nRNA is present, pipeline should detect some.

    Skips assertion when nRNA signal is negligible or when gDNA
    contamination overwhelms the nRNA signal (making intronic nRNA
    indistinguishable from intronic gDNA).
    """
    if nrna_fraction == 0 or bench.n_nrna_expected <= 5:
        return
    # When gDNA >> nRNA, the intronic nRNA signal is buried
    if bench.n_gdna_expected > 3 * bench.n_nrna_expected:
        return
    assert bench.n_nrna_pipeline > 0, (
        f"No nRNA detected but expected {bench.n_nrna_expected:.0f}"
    )


# =====================================================================
# Scenario 1: Single-exon gene (unspliced)
# =====================================================================


class TestSingleExon:
    """Single unspliced gene + zero-expression negative control.

    The hardest case for gDNA separation: unspliced RNA is physically
    identical to gDNA overlapping the gene region.

    A separate multi-exon helper gene provides splice junctions for
    strand-model training — mirroring real data where spliced genes
    are always present alongside unspliced ones.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("single_exon", genome_length=8000, seed=SIM_SEED,
                       work_dir=tmp_path / "single_exon")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        # Multi-exon helper: provides spliced reads for strand training
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
        bench = _build_and_run(scenario, scenario_name="single_exon_base")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_transcript_accuracy(bench, max_abs_diff=3)
        _assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               scenario_name=f"single_exon_gdna_{gdna}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, scenario, nrna):
        """Nascent RNA on a single-exon gene is identical to mRNA."""
        bench = _build_and_run(scenario, nrna_fraction=nrna,
                               scenario_name=f"single_exon_nrna_{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = _build_and_run(scenario, strand_specificity=ss,
                               scenario_name=f"single_exon_ss_{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            _assert_transcript_accuracy(bench, max_abs_diff=55)

    @pytest.mark.parametrize("gdna,nrna", [
        (20, 0.3), (50, 0.5), (100, 0.3),
    ], ids=["gdna20_nrna30", "gdna50_nrna50", "gdna100_nrna30"])
    def test_gdna_nrna(self, scenario, gdna, nrna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna,
                               scenario_name=f"single_exon_g{gdna}_n{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        _assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna, strand_specificity=ss,
                               scenario_name=f"single_exon_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna,
                                  strand_specificity=ss)


# =====================================================================
# Scenario 2: Single spliced gene (multi-exon)
# =====================================================================


class TestSplicedGene:
    """Spliced two-exon gene + zero-expression negative control.

    Splice junctions provide a strong mRNA anchor, making gDNA/RNA
    separation easier than the unspliced case.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("spliced_gene", genome_length=5000, seed=SIM_SEED,
                       work_dir=tmp_path / "spliced_gene")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(3500, 3800)], "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        bench = _build_and_run(scenario, scenario_name="spliced_base")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_transcript_accuracy(bench, max_abs_diff=3)
        _assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               scenario_name=f"spliced_gdna_{gdna}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, scenario, nrna):
        """Spliced gene: nRNA produces intronic reads detectable by pipeline."""
        bench = _build_and_run(scenario, nrna_fraction=nrna,
                               scenario_name=f"spliced_nrna_{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench)
        _assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = _build_and_run(scenario, strand_specificity=ss,
                               scenario_name=f"spliced_ss_{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, strand_specificity=ss)
        max_loss = (1.0 - ss) + 0.20
        t1 = next(t for t in bench.transcripts if t.t_id == "t1")
        assert t1.observed >= t1.expected * (1.0 - max_loss), (
            f"t1: {t1.observed:.0f} too low vs {t1.expected} at ss={ss}"
        )

    @pytest.mark.parametrize("gdna,nrna", [
        (20, 0.3), (50, 0.5), (100, 0.3),
    ], ids=["gdna20_nrna30", "gdna50_nrna50", "gdna100_nrna30"])
    def test_gdna_nrna(self, scenario, gdna, nrna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna,
                               scenario_name=f"spliced_g{gdna}_n{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        _assert_gdna_accuracy(bench, gdna)
        _assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna, strand_specificity=ss,
                               scenario_name=f"spliced_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna,
                                  strand_specificity=ss)


# =====================================================================
# Scenario 3: Multiple non-overlapping genes
# =====================================================================


class TestNonOverlappingGenes:
    """Two non-overlapping expressed genes + negative control.

    g1 (+, spliced) and g2 (−, single-exon) at distant loci.
    No spatial overlap → no ambiguity between expressed genes.
    """

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
        bench = _build_and_run(scenario,
                               scenario_name="nonoverlap_base")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_transcript_accuracy(bench, max_abs_diff=3)
        _assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               scenario_name=f"nonoverlap_gdna_{gdna}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        if gdna == 0:
            _assert_transcript_accuracy(bench, max_abs_diff=3)
        else:
            _assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, scenario, nrna):
        bench = _build_and_run(scenario, nrna_fraction=nrna,
                               scenario_name=f"nonoverlap_nrna_{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench)
        _assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = _build_and_run(scenario, strand_specificity=ss,
                               scenario_name=f"nonoverlap_ss_{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, strand_specificity=ss)
        if ss >= 0.95:
            _assert_transcript_accuracy(bench, max_abs_diff=40)

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, scenario, gdna, nrna, ss):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna, strand_specificity=ss,
                               scenario_name=f"nonoverlap_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna,
                                  strand_specificity=ss)


# =====================================================================
# Scenario 4: Two isoforms of a single gene
# =====================================================================


class TestTwoIsoforms:
    """Two isoforms sharing 5' and 3' exons + negative control gene.

    t1 (major): 3 exons — shared_5p, middle, shared_3p.
    t2 (minor): 2 exons — shared_5p, shared_3p (skips middle exon).
    t_ctrl: single exon on opposite strand, abundance=0.
    """

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
            bench = _build_and_run(sc, n_fragments=1000,
                                   scenario_name=f"iso_fc_{fold_change}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
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
            bench = _build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                   scenario_name=f"iso_gdna_{gdna}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(
                    bench.total_expected, abs=5)
            else:
                _assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, tmp_path, nrna):
        sc = self._make_scenario(tmp_path, 100, 10, f"_nrna_{int(nrna*100)}")
        try:
            bench = _build_and_run(sc, nrna_fraction=nrna, n_fragments=1000,
                                   scenario_name=f"iso_nrna_{int(nrna*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
            _assert_nrna_detected(bench, nrna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 10, f"_ss_{ss}")
        try:
            bench = _build_and_run(sc, strand_specificity=ss, n_fragments=1000,
                                   scenario_name=f"iso_ss_{ss}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, strand_specificity=ss)
            max_loss = (1.0 - ss) + 0.20
            min_expected = bench.total_expected * (1.0 - max_loss)
            assert bench.total_rna_observed >= min_expected
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 10,
                                  f"_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   nrna_fraction=nrna,
                                   strand_specificity=ss, n_fragments=1000,
                                   scenario_name=f"iso_stress_{gdna}_{int(nrna*100)}_{int(ss*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna,
                                      strand_specificity=ss)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 5: Overlapping genes on opposite strands (antisense)
# =====================================================================


class TestOverlappingAntisense:
    """Two overlapping antisense genes + distant negative control.

    g1 (+) and g2 (−) overlap spatially — the most challenging
    scenario for strand-based disambiguation.  t_ctrl is a distant
    zero-expression gene.
    """

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

    @pytest.mark.parametrize("fold_change", ABUNDANCE_RATIOS,
                             ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS])
    def test_abundance_sweep(self, tmp_path, fold_change):
        g1_abund, g2_abund = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, g1_abund, g2_abund)
        try:
            bench = _build_and_run(sc, n_fragments=1000,
                                   scenario_name=f"anti_fc_{fold_change}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
            # Cross-gene antisense may inflate gDNA — allow up to 35% loss
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
            bench = _build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                   scenario_name=f"anti_gdna_{gdna}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                max_loss = 0.35
                min_exp = bench.total_expected * (1.0 - max_loss)
                assert bench.total_rna_observed >= min_exp
            else:
                _assert_gdna_accuracy(bench, gdna, max_rel_err=1.0)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, tmp_path, nrna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_nrna_{int(nrna*100)}")
        try:
            bench = _build_and_run(sc, nrna_fraction=nrna, n_fragments=1000,
                                   scenario_name=f"anti_nrna_{int(nrna*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100, f"_ss_{ss}")
        try:
            bench = _build_and_run(sc, strand_specificity=ss, n_fragments=1000,
                                   scenario_name=f"anti_ss_{ss}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, strand_specificity=ss)
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

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   nrna_fraction=nrna, strand_specificity=ss,
                                   n_fragments=1000,
                                   scenario_name=f"anti_stress_{gdna}_{int(nrna*100)}_{int(ss*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna,
                                      strand_specificity=ss)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 5b: Contained antisense transcript
# =====================================================================


class TestContainedAntisense:
    """Antisense gene fully contained within another gene's intron.

    g1 (+) has exons (100,600) and (3000,3500) with a large intron
    spanning 600-3000.  g2 (-) has exons (1000,1500) and (2000,2500),
    sitting entirely within g1's intron.  This geometry means ALL of
    g2's exonic sequence is intronic for g1, making it the hardest
    case for exon-only overlap filtering.

    t_ctrl is a distant zero-expression negative control.
    """

    def _make_scenario(self, tmp_path, g1_abundance, g2_abundance,
                       name_suffix=""):
        sc = Scenario("contained_anti" + name_suffix, genome_length=10000,
                       seed=SIM_SEED,
                       work_dir=tmp_path / ("contained_anti" + name_suffix))
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(100, 600), (3000, 3500)],
             "abundance": g1_abundance},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(1000, 1500), (2000, 2500)],
             "abundance": g2_abundance},
        ])
        sc.add_gene("g_ctrl", "+", [
            {"t_id": "t_ctrl", "exons": [(7500, 7800)], "abundance": 0},
        ])
        return sc

    @pytest.mark.parametrize("a_abund,b_abund", [
        (100, 10), (100, 50), (100, 100), (50, 100), (10, 100),
    ], ids=["10to1", "2to1", "1to1", "1to2", "1to10"])
    def test_abundance_ratio(self, tmp_path, a_abund, b_abund):
        sc = self._make_scenario(tmp_path, a_abund, b_abund,
                                  f"_a{a_abund}_b{b_abund}")
        try:
            bench = _build_and_run(sc, n_fragments=1000,
                                   scenario_name=f"contained_anti_{a_abund}v{b_abund}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
            # Both genes should retain most of their RNA
            max_loss = 0.40
            min_expected = bench.total_expected * (1.0 - max_loss)
            assert bench.total_rna_observed >= min_expected, (
                f"RNA loss too high: observed={bench.total_rna_observed:.0f}, "
                f"expected>={min_expected:.0f} (total_exp={bench.total_expected})"
            )
            # When one gene dominates, it should have more counts
            if a_abund > b_abund:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t1.observed > t2.observed, (
                    f"t1 should dominate: t1={t1.observed:.0f}, t2={t2.observed:.0f}"
                )
            elif b_abund > a_abund:
                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")
                assert t2.observed > t1.observed, (
                    f"t2 should dominate: t1={t1.observed:.0f}, t2={t2.observed:.0f}"
                )
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, tmp_path, gdna):
        sc = self._make_scenario(tmp_path, 100, 100, f"_gdna_{gdna}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                   scenario_name=f"contained_anti_gdna_{gdna}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                max_loss = 0.40
                min_exp = bench.total_expected * (1.0 - max_loss)
                assert bench.total_rna_observed >= min_exp
            else:
                _assert_gdna_accuracy(bench, gdna, max_rel_err=1.0)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("nrna", NRNA_FRACTIONS,
                             ids=[f"nrna_{int(f*100)}" for f in NRNA_FRACTIONS])
    def test_nrna_sweep(self, tmp_path, nrna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_nrna_{int(nrna*100)}")
        try:
            bench = _build_and_run(sc, nrna_fraction=nrna, n_fragments=1000,
                                   scenario_name=f"contained_anti_nrna_{int(nrna*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100, f"_ss_{ss}")
        try:
            bench = _build_and_run(sc, strand_specificity=ss, n_fragments=1000,
                                   scenario_name=f"contained_anti_ss_{ss}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, strand_specificity=ss)
            max_loss = (1.0 - ss) + 0.40
            min_exp = bench.total_expected * (1.0 - max_loss)
            assert bench.total_rna_observed >= min_exp
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   nrna_fraction=nrna, strand_specificity=ss,
                                   n_fragments=1000,
                                   scenario_name=f"contained_anti_stress_{gdna}_{int(nrna*100)}_{int(ss*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna,
                                      strand_specificity=ss)
        finally:
            sc.cleanup()


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
        # Spliced helper gene to train the strand model (paralogs
        # alone are all multi-gene so strand model gets 0 obs).
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
        """Equal-abundance unspliced paralogs → ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100)
        try:
            bench = _build_and_run(sc, include_multimap=True,
                                   scenario_name="paralogs_eq_unspliced")
            _assert_alignment(bench)
            _assert_negative_control(bench)
            assert bench.total_rna_observed == pytest.approx(
                bench.total_expected, abs=10)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            assert abs(t1.observed - t2.observed) < total * 0.15 + 5
        finally:
            sc.cleanup()

    def test_equal_spliced(self, tmp_path):
        """Equal-abundance spliced paralogs → ~50/50 split."""
        sc = self._make_scenario(tmp_path, 100, 100, spliced=True)
        try:
            bench = _build_and_run(sc, include_multimap=True,
                                   scenario_name="paralogs_eq_spliced")
            _assert_alignment(bench)
            _assert_negative_control(bench)
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
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   include_multimap=True,
                                   scenario_name=f"paralogs_gdna_{gdna}")
            _assert_alignment(bench)
            _assert_negative_control(bench, gdna_abundance=gdna)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total = t1.observed + t2.observed
            if total > 10:
                tol = 0.45 if gdna > 0 else 0.20
                assert abs(t1.observed - t2.observed) < total * tol + 5
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, tmp_path, nrna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_nrna_{int(nrna*100)}")
        try:
            bench = _build_and_run(sc, nrna_fraction=nrna,
                                   include_multimap=True,
                                   scenario_name=f"paralogs_nrna_{int(nrna*100)}")
            _assert_alignment(bench)
            _assert_negative_control(bench)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_ss_{ss}")
        try:
            bench = _build_and_run(sc, strand_specificity=ss,
                                   include_multimap=True,
                                   scenario_name=f"paralogs_ss_{ss}")
            _assert_alignment(bench)
            _assert_negative_control(bench, strand_specificity=ss)
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

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (50, 0.5, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g50n50s65"])
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  name_suffix=f"_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   nrna_fraction=nrna, strand_specificity=ss,
                                   include_multimap=True,
                                   scenario_name=f"paralogs_stress_{gdna}_{int(nrna*100)}_{int(ss*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna,
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
        sc.genome.edit(5000, sc.genome[500:800])  # copy exon 1
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(9500, 9800)], "abundance": 0},
        ])
        return sc

    @pytest.mark.parametrize("fold_change", ABUNDANCE_RATIOS,
                             ids=[f"fc_{r}" for r in ABUNDANCE_RATIOS])
    def test_abundance_sweep(self, tmp_path, fold_change):
        g1_abund, g2_abund = 100, 100 / fold_change
        sc = self._make_scenario(tmp_path, g1_abund, g2_abund,
                                  f"_fc_{fold_change}")
        try:
            bench = _build_and_run(sc, n_fragments=1000,
                                   include_multimap=True,
                                   scenario_name=f"dist_fc_{fold_change}")
            _assert_alignment(bench)
            _assert_negative_control(bench)
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
            bench = _build_and_run(sc, gdna_abundance=gdna, n_fragments=1000,
                                   include_multimap=True,
                                   scenario_name=f"dist_gdna_{gdna}")
            _assert_alignment(bench)
            _assert_negative_control(bench, gdna_abundance=gdna)
            if gdna == 0:
                assert bench.total_rna_observed == pytest.approx(
                    bench.total_expected, abs=10)
            else:
                _assert_gdna_accuracy(bench, gdna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, tmp_path, nrna):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_nrna_{int(nrna*100)}")
        try:
            bench = _build_and_run(sc, nrna_fraction=nrna, n_fragments=1000,
                                   include_multimap=True,
                                   scenario_name=f"dist_nrna_{int(nrna*100)}")
            _assert_alignment(bench)
            _assert_negative_control(bench)
            _assert_nrna_detected(bench, nrna)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, tmp_path, ss):
        sc = self._make_scenario(tmp_path, 100, 100, f"_ss_{ss}")
        try:
            bench = _build_and_run(sc, strand_specificity=ss, n_fragments=1000,
                                   include_multimap=True,
                                   scenario_name=f"dist_ss_{ss}")
            _assert_alignment(bench)
            _assert_negative_control(bench, strand_specificity=ss)
        finally:
            sc.cleanup()

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (50, 0.5, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g50n50s65"])
    def test_stress(self, tmp_path, gdna, nrna, ss):
        sc = self._make_scenario(tmp_path, 100, 100,
                                  f"_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        try:
            bench = _build_and_run(sc, gdna_abundance=gdna,
                                   nrna_fraction=nrna, strand_specificity=ss,
                                   n_fragments=1000, include_multimap=True,
                                   scenario_name=f"dist_stress_{gdna}_{int(nrna*100)}_{int(ss*100)}")
            _assert_alignment(bench)
            _assert_accountability(bench)
            _assert_negative_control(bench, gdna_abundance=gdna,
                                      strand_specificity=ss)
        finally:
            sc.cleanup()


# =====================================================================
# Scenario 8: Two-exon gene with single-exon negative control
# =====================================================================
#
# Comprehensive evaluation: gene 1 is a spliced two-exon gene
# (abundance=100), gene 2 is a single-exon zero-expression control.
# Sweeps all three axes (gDNA, nRNA, SS) plus interactions.


class TestTwoExonWithControl:
    """Two-exon spliced gene + single-exon zero-expression control.

    Comprehensive parameter sweep for the 3-pool EM: mature RNA,
    nascent RNA, and gDNA.

    g1/t1: two exons, + strand, abundance=100
    g_ctrl/t_ctrl: single exon, − strand, abundance=0
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario("two_exon_ctrl", genome_length=5000, seed=SIM_SEED,
                       work_dir=tmp_path / "two_exon_ctrl")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 400), (900, 1100)],
             "abundance": 100},
        ])
        sc.add_gene("g_ctrl", "-", [
            {"t_id": "t_ctrl", "exons": [(3500, 3800)],
             "abundance": 0},
        ])
        yield sc
        sc.cleanup()

    def test_baseline(self, scenario):
        """Pure mature RNA: all counts to t1, none to t_ctrl."""
        bench = _build_and_run(scenario, scenario_name="ctrl_baseline")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_transcript_accuracy(bench, max_abs_diff=5)
        _assert_negative_control(bench)

    @pytest.mark.parametrize("gdna", GDNA_LEVELS,
                             ids=[f"gdna_{g}" for g in GDNA_LEVELS])
    def test_gdna_sweep(self, scenario, gdna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               scenario_name=f"ctrl_gdna_{gdna}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)
        if gdna > 0:
            _assert_gdna_accuracy(bench, gdna)

    @pytest.mark.parametrize("nrna", [0.1, 0.3, 0.5, 0.7],
                             ids=[f"nrna_{int(f*100)}" for f in
                                  [0.1, 0.3, 0.5, 0.7]])
    def test_nrna_sweep(self, scenario, nrna):
        bench = _build_and_run(scenario, nrna_fraction=nrna,
                               scenario_name=f"ctrl_nrna_{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench)
        _assert_nrna_detected(bench, nrna)

    @pytest.mark.parametrize("ss", STRAND_LEVELS,
                             ids=[f"ss_{s}" for s in STRAND_LEVELS])
    def test_strand_sweep(self, scenario, ss):
        bench = _build_and_run(scenario, strand_specificity=ss,
                               scenario_name=f"ctrl_ss_{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, strand_specificity=ss)
        t1 = next(t for t in bench.transcripts if t.t_id == "t1")
        max_loss = (1.0 - ss) + 0.15
        assert t1.observed >= t1.expected * (1.0 - max_loss), (
            f"t1: {t1.observed:.0f} too low vs {t1.expected} at ss={ss}"
        )

    @pytest.mark.parametrize("gdna,nrna", [
        (20, 0.3), (50, 0.3), (100, 0.3), (20, 0.5), (50, 0.5),
    ], ids=["gdna20_nrna30", "gdna50_nrna30", "gdna100_nrna30",
            "gdna20_nrna50", "gdna50_nrna50"])
    def test_gdna_nrna_interaction(self, scenario, gdna, nrna):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna,
                               scenario_name=f"ctrl_gdna{gdna}_nrna{int(nrna*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna)

    @pytest.mark.parametrize("gdna,ss", [
        (20, 0.95), (50, 0.9), (20, 0.65), (100, 0.8),
    ], ids=["gdna20_ss95", "gdna50_ss90", "gdna20_ss65", "gdna100_ss80"])
    def test_gdna_strand_interaction(self, scenario, gdna, ss):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               strand_specificity=ss,
                               scenario_name=f"ctrl_gdna{gdna}_ss{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna,
                                  strand_specificity=ss)

    @pytest.mark.parametrize("nrna,ss", [
        (0.3, 0.9), (0.5, 0.8), (0.3, 0.65), (0.7, 0.65),
    ], ids=["nrna30_ss90", "nrna50_ss80", "nrna30_ss65", "nrna70_ss65"])
    def test_nrna_strand_interaction(self, scenario, nrna, ss):
        bench = _build_and_run(scenario, nrna_fraction=nrna,
                               strand_specificity=ss,
                               scenario_name=f"ctrl_nrna{int(nrna*100)}_ss{ss}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, strand_specificity=ss)

    @pytest.mark.parametrize("gdna,nrna,ss", [
        (20, 0.3, 0.95), (50, 0.3, 0.9), (20, 0.5, 0.8),
        (50, 0.5, 0.65), (100, 0.3, 0.65),
    ], ids=["g20n30s95", "g50n30s90", "g20n50s80",
            "g50n50s65", "g100n30s65"])
    def test_combined_stress(self, scenario, gdna, nrna, ss):
        bench = _build_and_run(scenario, gdna_abundance=gdna,
                               nrna_fraction=nrna, strand_specificity=ss,
                               scenario_name=f"ctrl_g{gdna}_n{int(nrna*100)}_s{int(ss*100)}")
        _assert_alignment(bench)
        _assert_accountability(bench)
        _assert_negative_control(bench, gdna_abundance=gdna,
                                  strand_specificity=ss)
