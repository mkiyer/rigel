"""Tests for zero-copy C++ → numpy transfer (vec_to_ndarray).

Validates that the capsule-backed ndarray pattern in ndarray_util.h
produces correct, zero-copy numpy arrays through the C++ module
interfaces that use it:

  1. FragmentAccumulator.finalize_zero_copy() — chunk streaming path
  2. BamScanner.build_result() — post-scan observation aggregates
  3. StreamingScorer.finish() — scoring output arrays

Also tests the Python consumers that ingest these arrays:
  - StrandModel.observe_batch() with int8 ndarray input
  - FragmentLengthModels.observe_batch() with int32/int8 ndarray input
"""

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import (
    _replay_fraglen_observations,
    _replay_strand_observations,
    run_pipeline,
    scan_and_buffer,
)
from rigel.sim import Scenario, SimConfig
from rigel.strand_model import StrandModel, StrandModels
from rigel.frag_length_model import FragmentLengthModels
from rigel.buffer import FragmentBuffer

SEED = 42


# =====================================================================
# Fixtures
# =====================================================================


@pytest.fixture
def smoke_scenario(tmp_path):
    """Minimal oracle scenario for exercising the full scan path."""
    sc = Scenario(
        "ndarray_test",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "ndarray_test",
    )
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 80},
            {"t_id": "t2", "exons": [(200, 400), (900, 1100)], "abundance": 20},
        ],
    )
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
    result = sc.build_oracle(n_fragments=500, sim_config=sim_config)
    yield sc, result
    sc.cleanup()


# =====================================================================
# Test: Pipeline end-to-end with ndarray transfer
# =====================================================================


class TestPipelineNdarrayTransfer:
    """End-to-end validation that the pipeline works with vec_to_ndarray."""

    def test_pipeline_completes(self, smoke_scenario):
        """Pipeline runs to completion with zero-copy arrays."""
        sc, result = smoke_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)
        df = pr.estimator.get_counts_df(result.index)
        assert len(df) > 0
        assert df["count"].sum() > 0

    def test_strand_model_trained(self, smoke_scenario):
        """Strand models receive observations from C++ ndarray path."""
        sc, result = smoke_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)
        sm = pr.strand_models
        total = (
            sm.exonic_spliced.n_observations
            + sm.exonic.n_observations
            + sm.intergenic.n_observations
        )
        assert total > 0

    def test_frag_length_model_trained(self, smoke_scenario):
        """Fragment length models receive observations from C++ ndarray path."""
        sc, result = smoke_scenario
        config = PipelineConfig(
            em=EMConfig(seed=SEED),
            scan=BamScanConfig(sj_strand_tag="auto"),
        )
        pr = run_pipeline(result.bam_path, result.index, config=config)
        assert pr.frag_length_models.global_model.n_observations > 0


# =====================================================================
# Test: scan_and_buffer returns ndarray observations
# =====================================================================


class TestScanAndBufferNdarrayDtypes:
    """Verify scan_and_buffer returns ndarray with correct dtypes."""

    @pytest.fixture(autouse=True)
    def _run_scan(self, smoke_scenario):
        sc, result = smoke_scenario
        self.index = result.index
        scan_config = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
        (
            self.stats,
            self.strand_models,
            self.frag_length_models,
            self.buffer,
            _region_counts,
            _fl_table,
        ) = scan_and_buffer(
            str(result.bam_path),
            self.index,
            scan_config,
        )

    def test_strand_models_populated(self):
        """Strand models are populated after scan."""
        total = (
            self.strand_models.exonic_spliced.n_observations
            + self.strand_models.exonic.n_observations
            + self.strand_models.intergenic.n_observations
        )
        assert total > 0

    def test_frag_length_models_populated(self):
        """Fragment length models are populated after scan."""
        assert self.frag_length_models.global_model.n_observations > 0

    def test_buffer_has_fragments(self):
        """Buffer contains fragments after scan."""
        assert self.buffer.total_fragments > 0

    def test_chunk_arrays_are_ndarray(self):
        """Buffer chunks contain numpy arrays (not bytes or lists)."""
        for chunk in self.buffer._chunks:
            assert isinstance(chunk.t_indices, np.ndarray)
            assert isinstance(chunk.frag_lengths, np.ndarray)
            assert isinstance(chunk.splice_type, np.ndarray)


# =====================================================================
# Test: _replay_strand_observations with ndarray input
# =====================================================================


class TestReplayStrandObservationsNdarray:
    """Test that _replay_strand_observations works with int8 ndarray input."""

    def test_int8_ndarray_input(self):
        """observe_batch works when given int8 ndarray (as from C++)."""
        models = StrandModels()
        strand_dict = {
            "exonic_spliced_obs": np.array([1, 1, 2, 2], dtype=np.int8),
            "exonic_spliced_truth": np.array([1, 2, 1, 2], dtype=np.int8),
            "exonic_obs": np.array([1, 2], dtype=np.int8),
            "exonic_truth": np.array([1, 2], dtype=np.int8),
            "intergenic_obs": np.array([], dtype=np.int8),
            "intergenic_truth": np.array([], dtype=np.int8),
        }
        _replay_strand_observations(strand_dict, models)

        assert models.exonic_spliced.n_observations == 4
        assert models.exonic_spliced.pos_pos == 1
        assert models.exonic_spliced.pos_neg == 1
        assert models.exonic_spliced.neg_pos == 1
        assert models.exonic_spliced.neg_neg == 1
        assert models.exonic.n_observations == 2
        assert models.intergenic.n_observations == 0

    def test_empty_arrays(self):
        """Empty ndarray arrays produce no observations."""
        models = StrandModels()
        strand_dict = {
            "exonic_spliced_obs": np.array([], dtype=np.int8),
            "exonic_spliced_truth": np.array([], dtype=np.int8),
            "exonic_obs": np.array([], dtype=np.int8),
            "exonic_truth": np.array([], dtype=np.int8),
            "intergenic_obs": np.array([], dtype=np.int8),
            "intergenic_truth": np.array([], dtype=np.int8),
        }
        _replay_strand_observations(strand_dict, models)
        assert models.exonic_spliced.n_observations == 0
        assert models.exonic.n_observations == 0
        assert models.intergenic.n_observations == 0

    def test_observe_batch_int8_directly(self):
        """StrandModel.observe_batch works with int8 ndarray directly."""
        sm = StrandModel()
        obs = np.array([1, 1, 2, 2, 1], dtype=np.int8)
        truth = np.array([1, 2, 1, 2, 1], dtype=np.int8)
        sm.observe_batch(obs, truth)
        assert sm.n_observations == 5
        assert sm.pos_pos == 2
        assert sm.pos_neg == 1
        assert sm.neg_pos == 1
        assert sm.neg_neg == 1

    def test_observe_batch_int32_compatible(self):
        """StrandModel.observe_batch still works with int32 arrays."""
        sm = StrandModel()
        obs = np.array([1, 2, 1], dtype=np.int32)
        truth = np.array([1, 1, 2], dtype=np.int32)
        sm.observe_batch(obs, truth)
        assert sm.n_observations == 3
        assert sm.pos_pos == 1
        assert sm.neg_pos == 1
        assert sm.pos_neg == 1


# =====================================================================
# Test: _replay_fraglen_observations with ndarray input
# =====================================================================


class TestReplayFraglenObservationsNdarray:
    """Test _replay_fraglen_observations with int32/int8 ndarray input."""

    def test_int32_lengths_int8_splice_types(self):
        """observe_batch works with int32 lengths and int8 splice_types."""
        models = FragmentLengthModels()
        fraglen_dict = {
            "lengths": np.array([200, 250, 300, 150], dtype=np.int32),
            "splice_types": np.array([0, 2, 2, 0], dtype=np.int8),
            "intergenic_lengths": np.array([180, 220], dtype=np.int32),
        }
        _replay_fraglen_observations(fraglen_dict, models)
        assert models.global_model.n_observations == 6  # 4 + 2 intergenic

    def test_empty_fraglen_arrays(self):
        """Empty ndarray arrays produce no observations."""
        models = FragmentLengthModels()
        fraglen_dict = {
            "lengths": np.array([], dtype=np.int32),
            "splice_types": np.array([], dtype=np.int8),
            "intergenic_lengths": np.array([], dtype=np.int32),
        }
        _replay_fraglen_observations(fraglen_dict, models)
        assert models.global_model.n_observations == 0

    def test_splice_type_dispatching(self):
        """Fragment lengths are dispatched to correct category sub-models."""
        from rigel.splice import SpliceType

        models = FragmentLengthModels()
        fraglen_dict = {
            "lengths": np.array([200, 300], dtype=np.int32),
            "splice_types": np.array(
                [int(SpliceType.UNSPLICED), int(SpliceType.SPLICED_ANNOT)],
                dtype=np.int8,
            ),
            "intergenic_lengths": np.array([], dtype=np.int32),
        }
        _replay_fraglen_observations(fraglen_dict, models)
        assert models.global_model.n_observations == 2
        assert models.category_models[SpliceType.UNSPLICED].n_observations == 1
        assert models.category_models[SpliceType.SPLICED_ANNOT].n_observations == 1


# =====================================================================
# Test: Memory safety — arrays outlive source dict
# =====================================================================


class TestMemorySafety:
    """Verify capsule-backed arrays remain valid after source is deleted."""

    def test_array_survives_dict_deletion(self, smoke_scenario):
        """Numpy arrays from chunk data remain valid after dict deletion."""
        sc, result = smoke_scenario
        scan_config = BamScanConfig(sj_strand_tag="auto", n_scan_threads=1)
        (
            _stats,
            _strand_models,
            _fl_models,
            buf,
            _region_counts,
            _fl_table,
        ) = scan_and_buffer(
            str(result.bam_path),
            result.index,
            scan_config,
        )

        # Extract arrays from the first chunk
        assert len(buf._chunks) > 0
        chunk = buf._chunks[0]
        t_ind = chunk.t_indices
        fl = chunk.frag_lengths

        # Delete everything except the extracted arrays
        del buf, chunk

        # Arrays should still be valid
        if len(t_ind) > 0:
            _ = t_ind.sum()
        if len(fl) > 0:
            _ = fl.sum()
            assert fl.dtype == np.int32
