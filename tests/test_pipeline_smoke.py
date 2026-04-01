"""Smoke test for the full quantification pipeline.

Runs ``run_pipeline`` end-to-end on a minimal oracle scenario and
validates that output DataFrames have the expected schema, non-zero
counts, and internally consistent totals.  This catches import/interface
breakage that component-level unit tests might miss.
"""

import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, SimConfig

SEED = 42


def _make_scenario(tmp_path, *, n_fragments=500):
    """Build a two-gene oracle scenario with spliced and unspliced signal."""
    sc = Scenario(
        "smoke",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "smoke",
    )
    # Gene on + strand: two isoforms sharing an exon (triggers EM)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
    # Gene on - strand: single transcript
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(2500, 2700), (3000, 3200)], "abundance": 50},
        ],
    )
    sim_config = SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SEED,
    )
    return sc, sc.build_oracle(n_fragments=n_fragments, sim_config=sim_config)


def _run(result, index):
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    return run_pipeline(result.bam_path, index, config=config)


class TestPipelineSmoke:
    """End-to-end smoke tests for ``run_pipeline``."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.sc, self.result = _make_scenario(tmp_path)
        self.index = self.result.index
        self.pr = _run(self.result, self.index)
        yield
        self.sc.cleanup()

    # -- transcript-level output --

    def test_transcript_df_schema(self):
        df = self.pr.estimator.get_counts_df(self.index)
        required = {
            "transcript_id",
            "gene_id",
            "count",
            "tpm",
        }
        assert required.issubset(df.columns), (
            f"Missing columns: {required - set(df.columns)}"
        )

    def test_transcript_df_rows(self):
        df = self.pr.estimator.get_counts_df(self.index)
        n_annotated = sum(
            1 for tid in df["transcript_id"] if not tid.startswith("RIGEL_NRNA_")
        )
        assert n_annotated == 3, f"Expected 3 annotated transcripts, got {n_annotated}"

    def test_transcript_positive_mrna(self):
        df = self.pr.estimator.get_counts_df(self.index)
        assert (df["count"] >= 0).all()
        assert df["count"].sum() > 0

    # -- gene-level output --

    def test_gene_df_schema(self):
        df = self.pr.estimator.get_gene_counts_df(self.index)
        required = {"gene_id", "count", "tpm"}
        assert required.issubset(df.columns), (
            f"Missing columns: {required - set(df.columns)}"
        )

    def test_gene_df_rows(self):
        df = self.pr.estimator.get_gene_counts_df(self.index)
        assert len(df) == 2, f"Expected 2 genes, got {len(df)}"

    # -- loci output --

    def test_loci_df_schema(self):
        df = self.pr.estimator.get_loci_df()
        required = {"locus_id", "mrna", "gdna"}
        assert required.issubset(df.columns), (
            f"Missing columns: {required - set(df.columns)}"
        )

    def test_loci_df_nonempty(self):
        df = self.pr.estimator.get_loci_df()
        assert len(df) > 0

    # -- pipeline stats --

    def test_stats_fragment_count(self):
        assert self.pr.stats.n_fragments > 0

    # -- internal consistency --

    def test_tpm_sums_to_million(self):
        df = self.pr.estimator.get_counts_df(self.index)
        assert df["tpm"].sum() == pytest.approx(1e6, rel=1e-3)
