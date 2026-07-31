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
        """The payload's three axes must be exactly what the index's partition implies.

        ⭐ ``cuts`` are the CUT POSITIONS; a reference with ``k`` cuts owns ``k − 1`` nodes and
        ``k − 2`` interior lines. The predecessor counted ``k`` boundary objects per reference — the
        ``k − 1`` interiors plus two data-free terminals — which is the axis S5.f retired.
        """
        index = oracle.index
        cuts, ref_cut_offsets, node_types = build_node_partition_arrays(index)
        payload = _scan(oracle)

        np.testing.assert_array_equal(payload.cut_positions, cuts)
        np.testing.assert_array_equal(payload.ref_cut_offsets, ref_cut_offsets)
        assert payload.n_refs == len(index.ref_names)
        assert node_types.shape == (payload.n_nodes,)  # one type per node

        diffs = np.diff(ref_cut_offsets)
        expected_nodes = int(np.sum(np.maximum(diffs - 1, 0)))
        expected_edges = int(np.sum(np.maximum(diffs - 2, 0)))
        assert payload.n_nodes == expected_nodes
        assert payload.n_edges == expected_edges
        # ⚠ E = N − (non-empty refs), stated a second way: the two derivations must agree.
        n_live_refs = int(np.sum(diffs > 1))
        assert payload.n_edges == payload.n_nodes - n_live_refs

    def test_fl_pools_emitted(self, oracle):
        """The scan emits the FIVE PURE fragment-length pools, binned at the same L as every other
        bank. ⚠ The pools are integer counts on a ``(N_FRAGMENT_POOLS, max_length + 1)`` grid — there
        is no ``fl_pool_mass`` and no separate ``fl_max_size``, because a pool is a histogram of the
        same molecule length the accumulator deposits by."""
        from rigel.calibration.fl import gdna_fl_mass
        from rigel.scan_payload import N_FRAGMENT_POOLS

        payload = _scan(oracle)
        assert payload.pool_lengths.shape == (N_FRAGMENT_POOLS, payload.max_length + 1)
        # gDNA pool aggregation is well-formed (non-negative, right length).
        gdna = gdna_fl_mass(payload)
        assert gdna.shape == (payload.max_length + 1,)
        assert float(gdna.sum()) >= 0.0

    def test_at_least_some_mass_deposited(self, oracle):
        payload = _scan(oracle)
        # ⭐ ONE tally answers this now: node_start_count is incremented once per ACCEPTED fragment, so
        # its total IS the deposit count. The predecessor had to add five arrays across two dtypes
        # because mass was fractional and carried separately from the integer flux.
        assert int(np.asarray(payload.node_start_count).sum()) > 0, "scanner deposited nothing"
        assert int(payload.qc.deposited) == int(np.asarray(payload.node_start_count).sum())
