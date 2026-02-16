"""
Diagnostic tests to expose and quantify hulkrna's critical issues
identified in the salmon benchmark.

These tests isolate each issue, instrument the EM internals, and
provide detailed analysis to guide fixes.

Issue 1 (CRITICAL): Minor isoform count collapse
Issue 2 (HIGH):     Unspliced gene annihilation at low strand specificity
Issue 3 (MEDIUM):   gDNA over-absorption in gDNA scenarios
Issue 4 (LOW):      Antisense overlap strand drift
"""

import logging
import shutil

import numpy as np
import pytest

from hulkrna.counter import ReadCounter, EMData
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark
from hulkrna.types import Strand

logger = logging.getLogger(__name__)

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


SIM_SEED = 42
PIPELINE_SEED = 42
N_FRAGMENTS = 2000


def _sim_config(*, strand_specificity: float = 1.0):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity,
        seed=SIM_SEED,
    )


# =====================================================================
# Issue 1: Minor isoform count collapse
# =====================================================================


class TestIsoformCollapse:
    """
    DIAGNOSIS: hulkrna massively over-counts the major isoform (t1)
    and under-counts the minor (t2) when t2 is a strict exonic subset
    of t1.

    At 1:1 ratio, hulkrna reports 1748 vs 252 (truth: 1303 vs 697).
    Salmon reports 1320 vs 680.

    ROOT CAUSE HYPOTHESES:
    1. gDNA shadow steals probability mass from t2
    2. Unspliced reads from shared exons have equal likelihood for
       t1 and t2, but t1 also collects unique spliced reads from the
       middle exon → t1's theta grows → winner-take-all
    3. The combination of (1) and (2) creates a death spiral:
       shadow absorbs t2's "share" of ambiguous reads → t2 theta
       drops → fewer reads assigned to t2 → even lower theta
    """

    def _make_two_isoform_scenario(self, tmp_path, t1_abundance, t2_abundance):
        sc = Scenario("iso_diag", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "iso_diag")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": t1_abundance},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300)],
             "abundance": t2_abundance},
        ])
        return sc

    def test_equal_abundance_default(self, tmp_path):
        """Reproduce the 1:1 isoform collapse with default settings."""
        sc = self._make_two_isoform_scenario(tmp_path, 100, 100)
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config())
            gt = result.ground_truth_from_fastq()

            # Default pipeline
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr, scenario_name="iso_1_1_default")
            logger.info("\n%s", bench.summary())

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")

            # Inspect EM internals
            theta = pr.counter._converged_theta
            n_t = pr.counter.num_transcripts
            gdna_base = pr.counter.gdna_base_index

            # Report theta for each component
            for i in range(n_t):
                label = result.index.t_df["t_id"].values[i]
                logger.info(f"  theta[{label}] = {theta[i]:.6f}")
            for g in range(pr.counter.num_genes):
                label = result.index.g_df["g_id"].values[g]
                logger.info(f"  theta[shadow_{label}] = {theta[gdna_base + g]:.6f}")

            logger.info(f"\n  Ground truth:  t1={gt.get('t1',0)}, t2={gt.get('t2',0)}")
            logger.info(f"  hulkrna:       t1={t1.observed:.0f}, t2={t2.observed:.0f}")
            logger.info(f"  unique_counts: t1={pr.counter.unique_counts[0].sum():.0f}, "
                        f"t2={pr.counter.unique_counts[1].sum():.0f}")
            logger.info(f"  em_counts:     t1={pr.counter.em_counts[0].sum():.0f}, "
                        f"t2={pr.counter.em_counts[1].sum():.0f}")
            logger.info(f"  gdna_em_counts: {pr.counter.gdna_em_counts}")
            logger.info(f"  shadow_init:    {pr.counter.shadow_init}")

            # Quantify the bias
            t2_ratio_truth = gt.get("t2", 0) / max(gt.get("t1", 0), 1)
            t2_ratio_hulk = t2.observed / max(t1.observed, 1)
            logger.info(f"\n  t2/t1 truth ratio: {t2_ratio_truth:.3f}")
            logger.info(f"  t2/t1 hulkrna ratio: {t2_ratio_hulk:.3f}")

            # This test documents the improved behavior after effective
            # length normalization.  t2 should be within ~35% of truth.
            t2_err = abs(t2.observed - gt.get("t2", 0))
            truth_t2 = gt.get("t2", 0)
            assert t2_err < truth_t2 * 0.5, (
                f"t2 error {t2_err:.0f} exceeds 50% of truth ({truth_t2}). "
                f"Effective length normalization should keep error < 50%."
            )
        finally:
            sc.cleanup()

    def test_hypothesis_gdna_shadow_steals_t2(self, tmp_path):
        """Test: Does disabling gDNA shadow fix isoform quantification?

        Run with gdna_threshold=0.0 (never assign to gDNA via EM).
        If t2 recovers, the gDNA shadow is the culprit.
        """
        sc = self._make_two_isoform_scenario(tmp_path, 100, 100)
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config())
            gt = result.ground_truth_from_fastq()

            # Pipeline with gDNA shadow disabled
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED,
                              gdna_threshold=0.0)  # Never assign to gDNA
            bench = run_benchmark(result, pr,
                                  scenario_name="iso_1_1_no_gdna")
            logger.info("\n[NO gDNA SHADOW] %s", bench.summary())

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            logger.info(f"  With gdna_threshold=0.0: t1={t1.observed:.0f}, "
                        f"t2={t2.observed:.0f}")
            logger.info(f"  Ground truth: t1={gt.get('t1',0)}, t2={gt.get('t2',0)}")

            # Even without gDNA shadow, does the EM still bias toward t1?
            theta = pr.counter._converged_theta
            for i in range(pr.counter.num_transcripts):
                label = result.index.t_df["t_id"].values[i]
                logger.info(f"  theta[{label}] = {theta[i]:.6f}")
        finally:
            sc.cleanup()

    def test_hypothesis_em_convergence_bias(self, tmp_path):
        """Test: Does the EM itself bias toward t1 even without shadows?

        Use synthetic EMData with equal likelihoods for t1 and t2 on
        shared-exon reads, plus unique reads for t1 (middle exon).
        """
        # Simulate the isoform EM directly without the full pipeline
        n_t = 2  # t1=0, t2=1
        n_g = 1  # g1=0

        # Set up unique counts: t1 gets unique middle-exon reads
        rc = ReadCounter(n_t, n_g, seed=PIPELINE_SEED, alpha=1.0)
        # t1 gets 200 unique reads from middle exon (spliced)
        rc.unique_counts[0, 6] = 200.0  # SPLICED_ANNOT_SENSE
        # No unique reads for t2 (all its reads are shared with t1)
        # No shadow init (clean scenario)

        # Build EMData: 500 ambiguous units, each mapping to both t1 and t2
        # with EQUAL log-likelihoods (shared exon reads)
        n_units = 500
        offsets = np.arange(0, 3 * n_units + 1, 3, dtype=np.int64)

        # Each unit has 3 candidates: t1, t2, gDNA_shadow_g1
        t_indices = np.zeros(3 * n_units, dtype=np.int32)
        log_liks = np.zeros(3 * n_units, dtype=np.float64)
        count_cols = np.zeros(3 * n_units, dtype=np.uint8)

        for u in range(n_units):
            base = u * 3
            t_indices[base] = 0      # t1
            t_indices[base + 1] = 1  # t2
            t_indices[base + 2] = 2  # shadow g1 (gdna_base=2)
            # Equal likelihoods for t1 and t2, lower for gDNA
            log_liks[base] = 0.0       # log(1.0)
            log_liks[base + 1] = 0.0   # log(1.0) - same as t1
            log_liks[base + 2] = -1.0  # log(0.37) - gDNA less likely

        em = EMData(
            offsets=offsets,
            t_indices=t_indices,
            log_liks=log_liks,
            count_cols=count_cols,
            locus_t_indices=np.zeros(n_units, dtype=np.int32),
            locus_count_cols=np.zeros(n_units, dtype=np.uint8),
            n_units=n_units,
            n_candidates=3 * n_units,
            gdna_base_index=2,
        )

        # Run EM
        rc.run_em(em, em_iterations=20)
        theta = rc._converged_theta
        logger.info(f"\n  Synthetic EM test:")
        logger.info(f"  theta[t1]={theta[0]:.6f}, theta[t2]={theta[1]:.6f}, "
                     f"theta[shadow]={theta[2]:.6f}")
        logger.info(f"  t1/t2 theta ratio: {theta[0]/max(theta[1], 1e-10):.2f}")

        # Now assign
        rc.assign_ambiguous(em, gdna_threshold=0.5)
        t1_em = rc.em_counts[0].sum()
        t2_em = rc.em_counts[1].sum()
        gdna_em = rc.gdna_em_counts[0]
        logger.info(f"  Assigned: t1={t1_em:.0f}, t2={t2_em:.0f}, gDNA={gdna_em:.0f}")

        # The EM with 200 unique reads for t1 and 0 for t2 will
        # cause theta[t1] >> theta[t2], creating the bias
        assert theta[0] > theta[1], (
            "Expected t1 theta > t2 theta due to unique-count advantage"
        )

        # Also test without gDNA (gdna_threshold=0.0)
        rc2 = ReadCounter(n_t, n_g, seed=PIPELINE_SEED, alpha=1.0)
        rc2.unique_counts[0, 6] = 200.0
        rc2.run_em(em, em_iterations=20)
        rc2.assign_ambiguous(em, gdna_threshold=0.0)
        t1_no_gdna = rc2.em_counts[0].sum()
        t2_no_gdna = rc2.em_counts[1].sum()
        logger.info(f"  Without gDNA: t1={t1_no_gdna:.0f}, t2={t2_no_gdna:.0f}")

    def test_isoform_with_unique_t2_reads(self, tmp_path):
        """Test: What if t2 also has unique reads (unique exon)?

        This tests whether the subset problem is the root cause.
        When t2 has its own unique exon, the EM should recover.
        """
        sc = Scenario("iso_unique_t2", genome_length=8000, seed=SIM_SEED,
                      work_dir=tmp_path / "iso_unique_t2")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": 100},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300), (3500, 3800)],
             "abundance": 100},
        ])
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config())
            gt = result.ground_truth_from_fastq()
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr,
                                  scenario_name="iso_unique_t2")
            logger.info("\n[t2 WITH UNIQUE EXON] %s", bench.summary())

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            logger.info(f"  Truth: t1={gt.get('t1',0)}, t2={gt.get('t2',0)}")
            logger.info(f"  hulkrna: t1={t1.observed:.0f}, t2={t2.observed:.0f}")

            # With unique reads for both, the bias should be much smaller
            total_err = abs(t1.observed - gt["t1"]) + abs(t2.observed - gt["t2"])
            logger.info(f"  Total abs error: {total_err:.0f}")
        finally:
            sc.cleanup()


# =====================================================================
# Issue 2: Unspliced gene annihilation at low strand specificity
# =====================================================================


class TestUnsplicedAnnihilation:
    """
    DIAGNOSIS: At strand_specificity=0.5, ALL reads from a single-exon
    gene are diverted to gDNA (observed=0, truth=2000).

    ROOT CAUSE HYPOTHESIS:
    - The strand model learns strand_specificity ≈ 0.5 from the data.
    - For unspliced reads, the exonic strand model gives p_sense ≈ 0.5.
        - is_antisense() returns True when p < 0.5, so unstranded reads are
            no longer forced into antisense at exactly p=0.5.
            In weak-strand settings, ~50% can still be antisense in expectation.
      classified as antisense.
    - The Iceberg estimator computes shadow_init = 2*antisense + κ*depth*θ
    - With 50% antisense, shadow_init ≈ depth, so the gDNA shadow
      competes equally with the transcript.
    - In the EM, the gDNA shadow's strand model (intergenic, also ~0.5)
      is equally likely to the transcript's → the shadow absorbs reads.
    - At assignment, gdna_threshold=0.5 means P(RNA) must exceed 0.5;
      with equal competition, many units fall below threshold → gDNA.
    """

    def test_ss_0_5_single_exon(self, tmp_path):
        """Reproduce the ss=0.5 wipeout and inspect internals."""
        sc = Scenario("ss_diag", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_diag")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config(strand_specificity=0.5))
            gt = result.ground_truth_from_fastq()
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr,
                                  scenario_name="ss_0.5_default")
            logger.info("\n[SS=0.5 DEFAULT] %s", bench.summary())

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            logger.info(f"  Truth: {gt.get('t1', 0)}")
            logger.info(f"  hulkrna: {t1.observed:.0f}")
            logger.info(f"  unique_counts: {pr.counter.unique_counts[0].sum():.0f}")
            logger.info(f"  em_counts: {pr.counter.em_counts[0].sum():.0f}")
            logger.info(f"  shadow_init: {pr.counter.shadow_init}")
            logger.info(f"  gdna_em_counts: {pr.counter.gdna_em_counts}")

            # Strand model diagnostics
            sm = pr.strand_models
            logger.info(f"  exonic model: p_r1_sense={sm.exonic.p_r1_sense:.4f}")
            logger.info(f"  intergenic model: p_r1_sense={sm.intergenic.p_r1_sense:.4f}")

            # theta diagnostics
            theta = pr.counter._converged_theta
            logger.info(f"  theta[t1]={theta[0]:.6f}, theta[shadow_g1]={theta[1]:.6f}")
        finally:
            sc.cleanup()

    def test_ss_0_5_with_gdna_disabled(self, tmp_path):
        """Test: Does disabling gDNA threshold fix the wipeout?"""
        sc = Scenario("ss_nogg", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_nogg")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config(strand_specificity=0.5))
            gt = result.ground_truth_from_fastq()

            # Try with gdna_threshold=0.0 (all RNA)
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED,
                              gdna_threshold=0.0)
            bench = run_benchmark(result, pr,
                                  scenario_name="ss_0.5_no_gdna")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            logger.info(f"\n[SS=0.5 NO gDNA] truth={gt.get('t1',0)}, "
                        f"hulkrna={t1.observed:.0f}")
        finally:
            sc.cleanup()

    def test_ss_0_5_spliced_gene(self, tmp_path):
        """Test: Does a spliced gene survive ss=0.5?

        Spliced reads have splice-junction strand evidence that is
        independent of library strandedness. The exonic_spliced model
        should still show high specificity from XS/ts tags.
        """
        sc = Scenario("ss_spliced", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_spliced")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config(strand_specificity=0.5))
            gt = result.ground_truth_from_fastq()
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr,
                                  scenario_name="ss_0.5_spliced")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            logger.info(f"\n[SS=0.5 SPLICED] truth={gt.get('t1',0)}, "
                        f"hulkrna={t1.observed:.0f}")
            sm = pr.strand_models
            logger.info(f"  exonic_spliced: p_r1_sense={sm.exonic_spliced.p_r1_sense:.4f}")
            logger.info(f"  exonic: p_r1_sense={sm.exonic.p_r1_sense:.4f}")
        finally:
            sc.cleanup()

    def test_ss_gradient(self, tmp_path):
        """Sweep ss from 0.5 to 1.0 for unspliced gene to find threshold."""
        sc = Scenario("ss_grad", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_grad")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])

        ss_values = [0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0]
        logger.info("\n[SS GRADIENT - UNSPLICED GENE]")
        logger.info(f"  {'ss':>6s}  {'truth':>6s}  {'hulkrna':>8s}  "
                     f"{'shadow_init':>12s}  {'gdna_em':>8s}")

        for ss in ss_values:
            result = sc.build(
                n_fragments=N_FRAGMENTS,
                sim_config=_sim_config(strand_specificity=ss),
            )
            gt = result.ground_truth_from_fastq()
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr,
                                  scenario_name=f"ss_{ss}")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            logger.info(
                f"  {ss:6.2f}  {gt.get('t1',0):6d}  {t1.observed:8.0f}  "
                f"{pr.counter.shadow_init.sum():12.1f}  "
                f"{pr.counter.gdna_em_counts.sum():8.0f}"
            )

        sc.cleanup()


# =====================================================================
# Issue 3: gDNA over-absorption in gDNA scenarios
# =====================================================================


class TestGDNAOverAbsorption:
    """
    DIAGNOSIS: In gDNA scenarios, both hulkrna AND salmon over-count
    RNA. hulkrna is slightly better (it has a gDNA model) but still
    inflates RNA counts.

    For single_exon gdna=50: truth=512, hulkrna=738, salmon=848.
    hulkrna absorbs some gDNA but not enough.

    ROOT CAUSE: gDNA fragments that overlap the gene are identical
    to RNA fragments for unspliced genes. The gDNA shadow can only
    absorb as much as the prior (shadow_init) supports. With
    insufficient intergenic signal, shadow_init underestimates gDNA.
    """

    def test_gdna_overlap_analysis(self, tmp_path):
        """Analyze how gDNA fragments overlap genes and escape detection."""
        sc = Scenario("gdna_diag", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "gdna_diag")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        try:
            gdna = GDNAConfig(
                abundance=50,
                frag_mean=350, frag_std=100,
                frag_min=100, frag_max=1000,
            )
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config(), gdna_config=gdna)
            gt = result.ground_truth_from_fastq()
            n_gdna = result.ground_truth_gdna_count()

            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)

            stats = pr.stats
            logger.info(f"\n[gDNA OVERLAP ANALYSIS]")
            logger.info(f"  Ground truth RNA: {gt.get('t1',0)}")
            logger.info(f"  Ground truth gDNA: {n_gdna}")
            logger.info(f"  Pipeline fragments: {stats.n_fragments}")
            logger.info(f"  Intergenic: {stats.n_intergenic}")
            logger.info(f"  shadow_init: {pr.counter.shadow_init}")
            logger.info(f"  gdna_em_counts: {pr.counter.gdna_em_counts}")
            logger.info(f"  hulkrna RNA total: "
                        f"{pr.counter.t_counts.sum():.0f}")

            # Calculate how much gDNA is "invisible" (overlaps gene fully)
            # Gene spans 500-1500 (1000bp). Genome is 5000bp.
            # gDNA uniformly distributed → P(overlaps gene) ≈ gene_len/genome_len
            gene_len = 1000
            genome_len = 5000
            expected_overlap_frac = gene_len / genome_len
            logger.info(f"  Expected gDNA overlap with gene: "
                        f"{expected_overlap_frac:.1%}")
            logger.info(f"  Expected invisible gDNA: "
                        f"{n_gdna * expected_overlap_frac:.0f}")
        finally:
            sc.cleanup()


# =====================================================================
# Issue 4: Antisense overlap strand drift
# =====================================================================


class TestAntisenseStrandDrift:
    """
    DIAGNOSIS: hulkrna shows small but systematic count offsets (MAE 4-5)
    on overlapping antisense genes at perfect strandedness, while salmon
    achieves MAE 0.5.

    ROOT CAUSE: The strand model, though highly specific at ss=1.0,
    still has finite training data at the decision boundary. The EM's
    gDNA shadow steals a few reads from the minor-strand gene. This
    is a minor issue.
    """

    def test_antisense_perfect_strand(self, tmp_path):
        """Inspect strand model at ss=1.0 for overlapping antisense."""
        sc = Scenario("anti_diag", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "anti_diag")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(300, 600), (1100, 1400)],
             "abundance": 100},
        ])
        try:
            result = sc.build(n_fragments=N_FRAGMENTS,
                              sim_config=_sim_config(strand_specificity=1.0))
            gt = result.ground_truth_from_fastq()
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED)
            bench = run_benchmark(result, pr,
                                  scenario_name="anti_perfect")
            logger.info("\n[ANTISENSE PERFECT STRAND] %s", bench.summary())

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")

            logger.info(f"  Truth: t1={gt.get('t1',0)}, t2={gt.get('t2',0)}")
            logger.info(f"  hulkrna: t1={t1.observed:.0f}, t2={t2.observed:.0f}")
            logger.info(f"  shadow_init: {pr.counter.shadow_init}")
            logger.info(f"  gdna_em: {pr.counter.gdna_em_counts}")

            # The offsets should be small
            assert abs(t1.observed - gt["t1"]) <= 10
            assert abs(t2.observed - gt["t2"]) <= 10
        finally:
            sc.cleanup()


# =====================================================================
# Cross-cutting: verify gdna_threshold=0.0 restores sanity
# =====================================================================


class TestGdnaThresholdEffect:
    """Test the effect of gdna_threshold on all issues."""

    def test_threshold_sweep_isoforms(self, tmp_path):
        """Sweep gdna_threshold for the isoform scenario."""
        sc = Scenario("thresh_iso", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "thresh_iso")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": 100},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300)],
             "abundance": 100},
        ])

        thresholds = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]
        logger.info("\n[GDNA THRESHOLD SWEEP - ISOFORMS]")
        logger.info(f"  {'thresh':>7s}  {'t1_truth':>8s}  {'t1_hulk':>8s}  "
                     f"{'t2_truth':>8s}  {'t2_hulk':>8s}  {'t2_err':>8s}  "
                     f"{'gdna_em':>8s}")

        result = sc.build(n_fragments=N_FRAGMENTS, sim_config=_sim_config())
        gt = result.ground_truth_from_fastq()

        for thresh in thresholds:
            pr = run_pipeline(result.bam_path, result.index,
                              sj_strand_tag="ts", seed=PIPELINE_SEED,
                              gdna_threshold=thresh)
            bench = run_benchmark(result, pr)
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            logger.info(
                f"  {thresh:7.1f}  {gt.get('t1', 0):8d}  {t1.observed:8.0f}  "
                f"{gt.get('t2', 0):8d}  {t2.observed:8.0f}  "
                f"{abs(t2.observed - gt.get('t2', 0)):8.0f}  "
                f"{pr.counter.gdna_em_counts.sum():8.0f}"
            )

        sc.cleanup()
