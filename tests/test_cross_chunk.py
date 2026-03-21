"""Cross-chunk regression test for fused_score_buffer.

Verifies that splitting fragments across multiple buffer chunks
produces identical quantification output to processing in a single
chunk.  This catches chunk-boundary bugs in the C++ scoring and
EM data preparation path.
"""

import numpy as np
import pandas as pd
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, SimConfig

SEED = 42
N_FRAGS = 600  # enough to fill several tiny chunks and exercise EM


def _make_scenario(tmp_path):
    """Two overlapping genes → ambiguous fragments → EM required."""
    sc = Scenario(
        "crosschunk",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "crosschunk",
    )
    # Two isoforms on + strand sharing first exon → some ambiguous frags
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 60},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 40},
        ],
    )
    # Second gene on - strand
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
    return sc, sc.build_oracle(n_fragments=N_FRAGS, sim_config=sim_config)


def _run_with_chunk_size(result, index, chunk_size):
    config = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional"),
        scan=BamScanConfig(sj_strand_tag="auto", chunk_size=chunk_size),
    )
    return run_pipeline(result.bam_path, index, config=config)


class TestCrossChunkRegression:
    """Ensure chunk-size has no effect on quantification output."""

    @pytest.fixture(autouse=True)
    def _setup(self, tmp_path):
        self.sc, self.result = _make_scenario(tmp_path)
        self.index = self.result.index
        # Baseline: all fragments in one chunk
        self.pr_big = _run_with_chunk_size(
            self.result, self.index, chunk_size=N_FRAGS + 100
        )
        # Test: many tiny chunks (forces chunk-boundary splits)
        self.pr_tiny = _run_with_chunk_size(
            self.result, self.index, chunk_size=10
        )
        yield
        self.sc.cleanup()

    def _compare_df(self, df_big, df_tiny, sort_col):
        df_big = df_big.sort_values(sort_col).reset_index(drop=True)
        df_tiny = df_tiny.sort_values(sort_col).reset_index(drop=True)
        pd.testing.assert_frame_equal(df_big, df_tiny, atol=1e-10, rtol=1e-12)

    def test_transcript_counts_match(self):
        df_big = self.pr_big.estimator.get_counts_df(self.index)
        df_tiny = self.pr_tiny.estimator.get_counts_df(self.index)
        self._compare_df(df_big, df_tiny, "transcript_id")

    def test_gene_counts_match(self):
        df_big = self.pr_big.estimator.get_gene_counts_df(self.index)
        df_tiny = self.pr_tiny.estimator.get_gene_counts_df(self.index)
        self._compare_df(df_big, df_tiny, "gene_id")

    def test_loci_match(self):
        df_big = self.pr_big.estimator.get_loci_df()
        df_tiny = self.pr_tiny.estimator.get_loci_df()
        self._compare_df(df_big, df_tiny, "locus_id")

    def test_stats_match(self):
        assert self.pr_big.stats.n_fragments == self.pr_tiny.stats.n_fragments
        assert (
            self.pr_big.stats.deterministic_unambig_units
            == self.pr_tiny.stats.deterministic_unambig_units
        )
        assert (
            self.pr_big.stats.em_routed_ambig_same_strand_units
            == self.pr_tiny.stats.em_routed_ambig_same_strand_units
        )

    def test_gdna_em_count_match(self):
        np.testing.assert_allclose(
            self.pr_big.estimator.gdna_em_count,
            self.pr_tiny.estimator.gdna_em_count,
            atol=1e-10,
            rtol=1e-12,
        )

    def test_nrna_em_count_match(self):
        np.testing.assert_allclose(
            self.pr_big.estimator.nrna_em_count,
            self.pr_tiny.estimator.nrna_em_count,
            atol=1e-10,
            rtol=1e-12,
        )
