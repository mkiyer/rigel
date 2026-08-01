"""End-to-end: calibrate() recovers a planted gDNA strand overdispersion (docs/em_strand/03 §4).

Builds a payload of many intergenic (pure-gDNA) seed regions whose contained sense/antisense
counts are drawn from BetaBinom(½, od_true), runs the real calibrator, and checks the fitted
gdna_strand_overdispersion. This exercises the full Phase-2 path: substrate → count clue →
seed extraction → fit → result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from _synthetic import INV_LENGTH_SCALE, make_gdna_fl_pmf, make_strand_models

from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import CalibrationConfig
from rigel.scan_payload import N_FRAGMENT_POOLS, AccumulatorPayload, ScanQC

_STRAND_MODEL = make_strand_models(0.95, 200)
_FRAG_LEN = 50  # the delta the gDNA pmf sits at; every fixture fragment is this long


def _intergenic_betabinom_payload(n_nodes, depth, overdispersion, seed):
    """A 1-reference payload of ``n_nodes`` intergenic nodes; contained gDNA ~ BetaBinom(½, od).

    ⚠ **The edge axis is empty of COUNTS but not of ROWS.** One reference with ``k`` nodes owns
    ``k − 1`` lines, and the payload must carry them or the chain builder refuses it. Leaving the
    counts at zero is what makes this test isolate the CONTAINED-node seed arm: with no crossing
    fragments there are no edge seeds, so a recovered overdispersion can only have come from the nodes.
    """
    rng = np.random.default_rng(seed)
    a = 0.5 * (1.0 - overdispersion) / overdispersion if overdispersion > 0 else 1e9
    p = rng.beta(a, a, size=n_nodes)
    pos = rng.binomial(depth, p)
    neg = depth - pos

    contained = np.stack([pos, neg], axis=1).astype(np.uint32)
    n_edges = n_nodes - 1
    quantum = (2 * INV_LENGTH_SCALE + _FRAG_LEN) // (2 * _FRAG_LEN)

    def node_zeros(dtype):
        return np.zeros((n_nodes, 2), dtype=dtype)

    def edge_zeros(dtype):
        return np.zeros((n_edges, 2), dtype=dtype)

    payload = AccumulatorPayload(
        cut_positions=np.arange(n_nodes + 1, dtype=np.int64) * 100,
        ref_cut_offsets=np.array([0, n_nodes + 1], dtype=np.int64),
        ref_node_offsets=np.array([0, n_nodes], dtype=np.int64),
        ref_edge_offsets=np.array([0, n_edges], dtype=np.int64),
        ref_sj_offsets=np.array([0, 0], dtype=np.int64),
        node_contained_count=contained,
        node_contained_inv_length_sum=(contained.astype(np.uint64) * np.uint64(quantum)),
        node_contained_length_sum=(contained.astype(np.uint64) * np.uint64(_FRAG_LEN)),
        node_spanning_count=node_zeros(np.uint32),
        node_spanning_inv_length_sum=node_zeros(np.uint64),
        node_spanning_length_sum=node_zeros(np.uint64),
        node_start_count=contained.sum(axis=1).astype(np.uint32),
        edge_unspliced_count=edge_zeros(np.uint32),
        edge_unspliced_inv_length_sum=edge_zeros(np.uint64),
        edge_unspliced_length_sum=edge_zeros(np.uint64),
        edge_spliced_count=edge_zeros(np.uint32),
        edge_spliced_inv_length_sum=edge_zeros(np.uint64),
        edge_spliced_length_sum=edge_zeros(np.uint64),
        sj_count=np.zeros((0, 2), dtype=np.uint32),
        sj_inv_length_sum=np.zeros((0, 2), dtype=np.uint64),
        sj_length_sum=np.zeros((0, 2), dtype=np.uint64),
        pool_lengths=np.zeros((N_FRAGMENT_POOLS, 201), dtype=np.int64),
        deposited_lengths=np.zeros(201, dtype=np.uint32),
        qc=ScanQC(
            deposited=int(contained.sum()),
            dropped_too_long=0,
            dropped_empty=0,
            dropped_strand_undefined=0,
            deferred_undetermined_gap=0,
            unannotated_introns=0,
            contradictory_sj_strand=0,
            introns_absorbed=0,
        ),
        max_length=200,
        n_refs=1,
    )
    starts = np.arange(n_nodes, dtype=np.int64) * 100
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n_nodes, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * n_nodes, dtype="string"),
            "start": starts,
            "end": starts + 100,
            "length": np.full(n_nodes, 100, dtype=np.int64),
            "signature": np.zeros(n_nodes, dtype=np.uint8),  # intergenic ⇒ count-observable
        }
    )
    return payload, RegionArrays.from_frame(region_df, {"chr1": 0})


def _calibrate(payload, ra):
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=_STRAND_MODEL,
        gdna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        rna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        config=CalibrationConfig(),
    )


@pytest.mark.parametrize("od_true", [0.05, 0.10, 0.20])
def test_calibrate_recovers_overdispersion(od_true):
    payload, ra = _intergenic_betabinom_payload(
        n_nodes=400, depth=150, overdispersion=od_true, seed=int(od_true * 1000)
    )
    result = _calibrate(payload, ra)
    assert result.gdna_strand_overdispersion == pytest.approx(od_true, rel=0.25, abs=0.02)


def test_calibrate_binomial_gdna_floors_to_zero():
    """Non-overdispersed (50/50) gDNA → the identifiability gate floors od to 0 (Binomial)."""
    payload, ra = _intergenic_betabinom_payload(n_nodes=400, depth=150, overdispersion=0.0, seed=7)
    assert _calibrate(payload, ra).gdna_strand_overdispersion < 0.02


def test_a_junction_free_library_calibrates(od_true=0.10):
    """⚠ ``n_sj == 0`` is legal and must not be confused with "no junction flux": this payload's
    references are all single-node-signature intergenic, so the graph has no junction edge at all.
    ``calibrate`` defaults to an empty junction axis and the result carries a length-0 array."""
    payload, ra = _intergenic_betabinom_payload(
        n_nodes=400, depth=150, overdispersion=od_true, seed=11
    )
    result = _calibrate(payload, ra)
    assert result.n_junctions == 0
    assert result.mass_rna_junction.shape == (0,)
    assert result.n_edges == result.n_nodes - 1
