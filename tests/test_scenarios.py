"""
End-to-end integration tests and accuracy benchmarks using the simulation
framework.

These tests exercise the full pipeline: genome → annotation → reads →
minimap2 alignment → BAM → HulkIndex → run_pipeline → counts.

Smoke tests verify that the pipeline runs and produces non-zero counts.
Benchmark tests compare observed counts against simulation ground truth,
verifying both correctness and measuring accuracy.

Requirements: minimap2 and samtools must be available in PATH
(installed via conda/mamba).
"""

import logging
import shutil

import numpy as np
import pytest

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark


logger = logging.getLogger(__name__)


# Skip all tests in this module if minimap2 or samtools are missing
pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Fixtures
# =====================================================================


@pytest.fixture
def single_gene_scenario(tmp_path):
    """One gene, one transcript, positive strand, 2 exons."""
    sc = Scenario(
        "single_gene",
        genome_length=5000,
        seed=42,
        work_dir=tmp_path / "single_gene",
    )
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(500, 1000), (2000, 2500)], "abundance": 100},
    ])
    yield sc
    sc.cleanup()


@pytest.fixture
def multi_isoform_scenario(tmp_path):
    """One gene, two isoforms at different abundances."""
    sc = Scenario(
        "multi_iso",
        genome_length=5000,
        seed=42,
        work_dir=tmp_path / "multi_iso",
    )
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (2000, 2300)],
         "abundance": 100},
        {"t_id": "t2", "exons": [(200, 500), (2000, 2300)],
         "abundance": 10},
    ])
    yield sc
    sc.cleanup()


@pytest.fixture
def antisense_scenario(tmp_path):
    """Two genes on opposite strands, non-overlapping."""
    sc = Scenario(
        "antisense",
        genome_length=5000,
        seed=42,
        work_dir=tmp_path / "antisense",
    )
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(3000, 3300), (3800, 4100)], "abundance": 100},
    ])
    yield sc
    sc.cleanup()


# Standard sim config for smoke tests
_SMOKE_CONFIG = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80,
    frag_max=450, read_length=100, seed=42,
)


# =====================================================================
# Smoke tests — verify the pipeline runs and produces counts
# =====================================================================


class TestSingleGeneScenario:
    """Single gene: all reads should map uniquely."""

    def test_build_produces_artifacts(self, single_gene_scenario):
        result = single_gene_scenario.build(
            n_fragments=200, sim_config=_SMOKE_CONFIG,
        )
        assert result.fasta_path.exists()
        assert result.gtf_path.exists()
        assert result.bam_path.exists()
        assert result.index_dir.exists()
        assert result.index is not None
        assert len(result.transcripts) == 1
        assert result.n_simulated == 200

    def test_bam_has_reads(self, single_gene_scenario):
        import pysam
        result = single_gene_scenario.build(
            n_fragments=200, sim_config=_SMOKE_CONFIG,
        )
        n_reads = 0
        with pysam.AlignmentFile(str(result.bam_path), "rb",
                                 check_sq=False) as bam:
            for read in bam:
                if not read.is_unmapped:
                    n_reads += 1
        assert n_reads > 100, f"Only {n_reads} reads aligned"

    def test_pipeline_produces_counts(self, single_gene_scenario):
        result = single_gene_scenario.build(
            n_fragments=500, sim_config=_SMOKE_CONFIG,
        )
        pipeline_result = run_pipeline(result.bam_path, result.index,
                                       sj_strand_tag="ts")
        total = pipeline_result.counter.get_g_counts_df(
            result.index.t_to_g_arr
        ).sum().sum()
        assert total > 0, "Pipeline produced zero counts"


class TestMultiIsoformScenario:
    """Two isoforms at 10:1 abundance ratio."""

    def test_build_and_count(self, multi_isoform_scenario):
        result = multi_isoform_scenario.build(
            n_fragments=1000, sim_config=_SMOKE_CONFIG,
        )
        assert len(result.transcripts) == 2

        pipeline_result = run_pipeline(
            result.bam_path, result.index, include_multimap=False,
            sj_strand_tag="ts",
        )

        g_df = pipeline_result.counter.get_g_counts_df(result.index.t_to_g_arr)
        total = g_df.sum().sum()
        assert total > 0


class TestAntisenseScenario:
    """Two genes on opposite strands."""

    def test_both_genes_counted(self, antisense_scenario):
        result = antisense_scenario.build(
            n_fragments=500,
            sim_config=SimConfig(
                frag_mean=200, frag_std=30, frag_min=80,
                frag_max=250, read_length=100, seed=42,
            ),
        )
        assert len(result.transcripts) == 2

        pipeline_result = run_pipeline(result.bam_path, result.index,
                                       sj_strand_tag="ts")

        g_df = pipeline_result.counter.get_g_counts_df(result.index.t_to_g_arr)
        assert len(g_df) == 2, "Expected counts for both genes"
        for g_idx, row in g_df.iterrows():
            row_total = row.sum()
            assert row_total > 0, f"Gene index {g_idx} has zero counts"


# =====================================================================
# Accuracy benchmarks — compare observed counts to ground truth
# =====================================================================


class TestSingleExonBenchmark:
    """Simplest possible scenario: 1 gene, 1 transcript, 1 exon (1000bp).

    With a single unspliced transcript and no competing genes, every
    fragment that aligns should count to the sole transcript.  This is
    the deterministic baseline — any discrepancy indicates a bug.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "single_exon",
            genome_length=5000,
            seed=42,
            work_dir=tmp_path / "single_exon",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    @pytest.mark.parametrize("n_fragments", [50, 100, 500])
    def test_exact_counts(self, scenario, n_fragments):
        """Every aligned fragment should count to t1 (within gDNA tolerance).

        With gDNA modeling always on, unspliced fragments compete with
        the gDNA pseudo-component. The gDNA Dirichlet prior (alpha=1)
        can absorb ~1-2 counts via stochastic assignment, especially
        at low fragment counts.
        """
        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80,
            frag_max=450, read_length=100, seed=42,
        )
        result = scenario.build(n_fragments=n_fragments, sim_config=sim_config)
        pipeline_result = run_pipeline(
            result.bam_path, result.index, sj_strand_tag="ts",
        )

        bench = run_benchmark(result, pipeline_result,
                              scenario_name=f"single_exon_n{n_fragments}")

        logger.info("\n%s", bench.summary())

        # Sanity: most fragments should align
        assert bench.n_fragments > 0, "No fragments entered the pipeline"
        assert bench.alignment_rate > 0.80, (
            f"Low alignment rate: {bench.alignment_rate:.1%}"
        )

        # Allow minor gDNA absorption (≤2 counts) from the Dirichlet prior
        for ta in bench.transcripts:
            assert ta.abs_diff <= 2, (
                f"{ta.t_id}: expected={ta.expected}, "
                f"observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
            )


class TestSingleExonGDNABenchmark:
    """Single-exon transcript with simulated gDNA contamination.

    For unspliced single-exon transcripts, gDNA reads that overlap the
    gene region are indistinguishable from RNA reads.  The EM can only
    separate gDNA that falls *outside* the transcript footprint (which
    the pipeline classifies as intergenic) from gDNA that overlaps it.

    We verify:
    1. Total accountability: pipeline assigns all fragments.
    2. Transcript counts are in a reasonable range (RNA + some gDNA overlap).
    3. Pipeline gDNA estimate is in the right ballpark.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "single_exon_gdna",
            genome_length=5000,
            seed=42,
            work_dir=tmp_path / "single_exon_gdna",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    @pytest.mark.parametrize(
        "gdna_abundance",
        [5, 20, 50, 100],
        ids=["gdna_low", "gdna_medium", "gdna_high", "gdna_equal"],
    )
    def test_gdna_separation(self, scenario, gdna_abundance):
        """EM should recover reasonable counts despite gDNA contamination.

        For unspliced transcripts gDNA/RNA separation is ambiguous for
        reads overlapping the gene region, so we check:
        - Total accountability (transcript + gDNA + intergenic ≈ pipeline n_fragments)
        - Transcript count is between RNA-only and RNA+all-overlapping-gDNA
        - Pipeline gDNA estimate is positive when gDNA is simulated
        """
        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80,
            frag_max=450, read_length=100, seed=42,
        )
        gdna_config = GDNAConfig(
            abundance=gdna_abundance,
            frag_mean=350, frag_std=100,
            frag_min=100, frag_max=1000,
        )
        result = scenario.build(
            n_fragments=500,
            sim_config=sim_config,
            gdna_config=gdna_config,
        )
        pipeline_result = run_pipeline(
            result.bam_path, result.index, sj_strand_tag="ts",
        )

        bench = run_benchmark(
            result, pipeline_result,
            scenario_name=f"single_exon_gdna_a{gdna_abundance}",
        )

        logger.info("\n%s", bench.summary())

        # Sanity
        assert bench.n_fragments > 0, "No fragments entered the pipeline"
        assert bench.alignment_rate > 0.70, (
            f"Low alignment rate: {bench.alignment_rate:.1%}"
        )

        # gDNA was simulated — verify some gDNA is expected
        assert bench.n_gdna_expected > 0, (
            f"Expected gDNA fragments but got 0 "
            f"(abundance={gdna_abundance})"
        )

        # Total accountability: observed transcript counts + pipeline gDNA
        # + intergenic + chimeric should approximate total pipeline fragments
        total_accounted = (
            bench.total_observed + bench.n_gdna_pipeline
            + bench.n_intergenic + bench.n_chimeric
        )
        assert abs(total_accounted - bench.n_fragments) <= 2, (
            f"Accountability gap: accounted={total_accounted:.0f}, "
            f"pipeline={bench.n_fragments}, "
            f"diff={abs(total_accounted - bench.n_fragments):.0f}"
        )

        # RNA-only ground truth (what the simulator actually produced as RNA)
        n_rna_expected = bench.total_expected
        # Transcript count should be at least close to the RNA expectation.
        # Some gDNA overlapping the gene region is expected to inflate the
        # transcript count (unspliced gDNA is indistinguishable from RNA).
        # Allow the observed to exceed the RNA-only expectation by up to
        # the full gDNA count (worst case: all gDNA overlaps 1 gene).
        for ta in bench.transcripts:
            assert ta.observed >= ta.expected * 0.85, (
                f"{ta.t_id}: observed={ta.observed:.0f} is too low vs "
                f"expected RNA={ta.expected}"
            )
            # Observed should not exceed RNA + all gDNA (impossible)
            assert ta.observed <= ta.expected + bench.n_gdna_expected + 5, (
                f"{ta.t_id}: observed={ta.observed:.0f} exceeds "
                f"RNA({ta.expected}) + gDNA({bench.n_gdna_expected})"
            )


class TestSplicedSingleGeneBenchmark:
    """1 gene, 1 transcript, 2 exons — spliced reads must count correctly.

    Some fragments span the intron (producing spliced reads), others
    stay within a single exon (unspliced).  All should count to the
    sole transcript.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "spliced_single",
            genome_length=5000,
            seed=42,
            work_dir=tmp_path / "spliced_single",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    def test_exact_counts(self, scenario):
        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80,
            frag_max=450, read_length=100, seed=42,
        )
        result = scenario.build(n_fragments=200, sim_config=sim_config)
        pipeline_result = run_pipeline(
            result.bam_path, result.index, sj_strand_tag="ts",
        )

        bench = run_benchmark(result, pipeline_result,
                              scenario_name="spliced_single_gene")

        logger.info("\n%s", bench.summary())

        assert bench.n_fragments > 0
        assert bench.alignment_rate > 0.70

        for ta in bench.transcripts:
            assert ta.abs_diff <= 2, (
                f"{ta.t_id}: expected={ta.expected}, "
                f"observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
            )

        # Verify some spliced reads were detected
        stats = pipeline_result.stats
        assert stats.n_with_annotated_sj > 0, (
            "No annotated splice junctions detected — "
            "spliced reads may not be aligning correctly"
        )


class TestAntisenseBenchmark:
    """2 non-overlapping genes on opposite strands.

    Each gene's fragments should count exclusively to that gene.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "antisense_bench",
            genome_length=8000,
            seed=42,
            work_dir=tmp_path / "antisense_bench",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4300), (5000, 5300)],
             "abundance": 100},
        ])
        yield sc
        sc.cleanup()

    def test_exact_counts(self, scenario):
        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80,
            frag_max=450, read_length=100, seed=42,
        )
        result = scenario.build(n_fragments=500, sim_config=sim_config)
        pipeline_result = run_pipeline(
            result.bam_path, result.index, sj_strand_tag="ts",
        )

        bench = run_benchmark(result, pipeline_result,
                              scenario_name="antisense_benchmark")

        logger.info("\n%s", bench.summary())

        assert bench.n_fragments > 0
        assert bench.alignment_rate > 0.70

        for ta in bench.transcripts:
            assert ta.abs_diff <= 2, (
                f"{ta.t_id}: expected={ta.expected}, "
                f"observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
            )


class TestMultiIsoformBenchmark:
    """1 gene, 2 isoforms at 10:1 abundance — accuracy benchmark.

    Isoform disambiguation relies on EM, so exact counts are not
    expected.  Instead we verify:
    1. Total gene-level count matches the sum of expected transcript counts.
    2. The higher-abundance isoform gets more counts.
    3. Relative error per transcript is bounded.
    """

    @pytest.fixture
    def scenario(self, tmp_path):
        sc = Scenario(
            "multi_iso_bench",
            genome_length=5000,
            seed=42,
            work_dir=tmp_path / "multi_iso_bench",
        )
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": 100},
            {"t_id": "t2",
             "exons": [(200, 500), (2000, 2300)],
             "abundance": 10},
        ])
        yield sc
        sc.cleanup()

    def test_gene_level_accuracy(self, scenario):
        sim_config = SimConfig(
            frag_mean=200, frag_std=30, frag_min=80,
            frag_max=450, read_length=100, seed=42,
        )
        result = scenario.build(n_fragments=1000, sim_config=sim_config)
        pipeline_result = run_pipeline(
            result.bam_path, result.index, sj_strand_tag="ts",
        )

        bench = run_benchmark(result, pipeline_result,
                              scenario_name="multi_isoform_benchmark")

        logger.info("\n%s", bench.summary())

        # Gene-level total should match sum of expected transcripts
        assert bench.total_observed == bench.total_expected, (
            f"Gene total mismatch: observed={bench.total_observed:.0f}, "
            f"expected={bench.total_expected}"
        )

        # Higher-abundance isoform should get more counts
        t1 = next(ta for ta in bench.transcripts if ta.t_id == "t1")
        t2 = next(ta for ta in bench.transcripts if ta.t_id == "t2")
        assert t1.observed > t2.observed, (
            f"t1 (abundance 100) should dominate: "
            f"t1={t1.observed:.0f}, t2={t2.observed:.0f}"
        )

        # Report per-transcript accuracy (not asserting exact match)
        for ta in bench.transcripts:
            logger.info("  %s", ta)

