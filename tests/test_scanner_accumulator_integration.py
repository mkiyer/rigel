"""End-to-end check that BamScanner populates an AccumulatorPayload.

Builds a tiny oracle scenario, calls ``scan_and_buffer`` directly, and
verifies that the returned ``AccumulatorPayload`` has the expected shape
derived from the index's region partition and that fragments deposited
non-trivial fractional mass into at least one region or boundary.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.splice_graph import build_node_partition_arrays
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import AccumulatorPayload
from rigel.sim import ReadSimConfig, Scenario

SEED = 1234


@pytest.fixture
def oracle(tmp_path):
    sc = Scenario(
        "acc_int",
        genome_length=5000,
        seed=SEED,
        work_dir=tmp_path / "acc_int",
    )
    sc.add_gene(
        "g1",
        "+",
        [{"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 50}],
    )
    sc.add_gene(
        "g2",
        "-",
        [{"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 30}],
    )
    sim_config = ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SEED,
    )
    result = sc.build_oracle(n_fragments=200, sim_config=sim_config)
    yield result
    sc.cleanup()


def _scan(result) -> AccumulatorPayload:
    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    _, _, _, _, payload = scan_and_buffer(str(result.bam_path), result.index, config.scan)
    return payload


class TestScannerAccumulatorIntegration:
    def test_payload_is_populated(self, oracle):
        payload = _scan(oracle)
        assert payload is not None
        assert isinstance(payload, AccumulatorPayload)
        assert payload.n_strand_columns == 2

    def test_payload_shape_matches_index_partition(self, oracle):
        index = oracle.index
        boundaries, ref_pos_offsets, region_types = build_node_partition_arrays(index)
        payload = _scan(oracle)

        np.testing.assert_array_equal(payload.boundaries, boundaries)
        np.testing.assert_array_equal(payload.ref_pos_offsets, ref_pos_offsets)
        assert payload.n_refs == len(index.ref_names)
        assert region_types.shape == (payload.r_total,)  # one type per region

        # Per-ref: a ref with k boundary positions contributes (k-1) regions
        # and k boundary objects; an empty ref contributes (0, 0).
        diffs = np.diff(ref_pos_offsets)
        expected_regions = int(np.sum(np.maximum(diffs - 1, 0)))
        expected_boundaries = int(np.sum(np.where(diffs > 0, diffs, 0)))
        assert payload.r_total == expected_regions
        assert payload.b_obj_total == expected_boundaries

    def test_fl_pools_emitted(self, oracle):
        # PR 4c: the scan emits the gDNA FL pools (set_regions passes region_types
        # + fl_max_size), so the payload carries a (N_FL_POOLS, fl_max_size+1) grid.
        from rigel.calibration.fl import gdna_fl_mass
        from rigel.calibration.fl import N_FL_POOLS

        payload = _scan(oracle)
        assert payload.fl_pool_mass is not None
        assert payload.fl_max_size > 0
        assert payload.fl_pool_mass.shape == (N_FL_POOLS, payload.fl_max_size + 1)
        # gDNA pool aggregation is well-formed (non-negative, right length).
        gdna = gdna_fl_mass(payload)
        assert gdna.shape == (payload.fl_max_size + 1,)
        assert float(gdna.sum()) >= 0.0

    def test_at_least_some_mass_deposited(self, oracle):
        payload = _scan(oracle)
        # Either contained counts, boundary mass, or flux must be non-zero.
        total = (
            int(payload.region_contained.sum())
            + int(payload.boundary_flux_left.sum())
            + int(payload.boundary_flux_right.sum())
            + float(payload.boundary_mass_left.sum())
            + float(payload.boundary_mass_right.sum())
        )
        assert total > 0, "scanner did not deposit any fractional mass"
