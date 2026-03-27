"""Diagnostic scenarios: isoform collapse and unspliced low-strand.

These are focused regression tests for specific algorithm behaviors
rather than full sweep-based scenarios.
"""


from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, run_benchmark

from .conftest import sim_config, SIM_SEED, PIPELINE_SEED


class TestIsoformCollapse:
    """Minor isoform count collapse when t2 is a strict exonic subset of t1.

    Validates that effective-length normalization keeps t2 error < 50%.
    Also verifies recovery when t2 has its own unique exon.
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

    def test_equal_abundance(self, tmp_path):
        """1:1 isoform ratio: t2 should be within 50% of truth."""
        sc = self._make_two_isoform_scenario(tmp_path, 100, 100)
        try:
            result = sc.build_oracle(n_fragments=2000,
                                     sim_config=sim_config())
            gt = result.ground_truth_auto()
            pr = run_pipeline(
                result.bam_path, result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="ts"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name="iso_1_1")

            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            t2_err = abs(t2.observed - gt.get("t2", 0))
            truth_t2 = gt.get("t2", 0)
            assert t2_err < truth_t2 * 0.5, (
                f"t2 error {t2_err:.0f} exceeds 50% of truth ({truth_t2})"
            )
        finally:
            sc.cleanup()

    def test_unique_t2_exon_recovers(self, tmp_path):
        """When t2 has its own unique exon, the EM should recover."""
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
            result = sc.build_oracle(n_fragments=2000,
                                     sim_config=sim_config())
            gt = result.ground_truth_auto()
            pr = run_pipeline(
                result.bam_path, result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="ts"),
                ),
            )
            bench = run_benchmark(result, pr,
                                  scenario_name="iso_unique_t2")

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total_err = abs(t1.observed - gt["t1"]) + abs(t2.observed - gt["t2"])
            assert total_err < 0.5 * (gt["t1"] + gt["t2"]), (
                f"Total error {total_err:.0f} too high with unique exons"
            )
        finally:
            sc.cleanup()


class TestUnsplicedLowStrand:
    """At low SS, unspliced genes may lose reads to gDNA shadow.

    Spliced genes survive because splice-junction strand evidence
    is independent of library strandedness.
    """

    def test_ss_gradient_unspliced(self, tmp_path):
        """Sweep SS for unspliced gene: pipeline should retain reads at SS>=0.9."""
        sc = Scenario("ss_grad", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_grad")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
        ])

        for ss in [0.65, 0.9, 1.0]:
            result = sc.build_oracle(
                n_fragments=2000,
                sim_config=sim_config(strand_specificity=ss),
            )
            gt = result.ground_truth_auto()
            pr = run_pipeline(
                result.bam_path, result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="ts"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name=f"ss_{ss}")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            if ss >= 0.9:
                assert t1.observed >= gt.get("t1", 0) * 0.5, (
                    f"t1 too low at ss={ss}: {t1.observed:.0f} "
                    f"vs truth {gt.get('t1', 0)}"
                )

        sc.cleanup()

    def test_spliced_survives_low_ss(self, tmp_path):
        """Spliced gene retains reads even at SS=0.65."""
        sc = Scenario("ss_spliced", genome_length=5000, seed=SIM_SEED,
                      work_dir=tmp_path / "ss_spliced")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(200, 500), (1000, 1300)],
             "abundance": 100},
        ])
        try:
            result = sc.build_oracle(
                n_fragments=2000,
                sim_config=sim_config(strand_specificity=0.65),
            )
            gt = result.ground_truth_auto()
            pr = run_pipeline(
                result.bam_path, result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="ts"),
                ),
            )
            bench = run_benchmark(result, pr,
                                  scenario_name="ss_0.65_spliced")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            assert t1.observed >= gt.get("t1", 0) * 0.3, (
                f"Spliced gene lost too many reads at ss=0.65: "
                f"{t1.observed:.0f} vs truth {gt.get('t1', 0)}"
            )
        finally:
            sc.cleanup()
