"""Hand-built synthetic payloads for calibration unit tests.

A plain helper module (not ``conftest.py``) to avoid colliding with the top-level ``conftest`` that
several tests import by name.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.scan_payload import N_FRAGMENT_POOLS, AccumulatorPayload, ScanQC

#: The fixed-point scale the accumulator deposits ``round(SCALE / placements)`` at. Decoded exactly once,
#: in ``CalibrationSubstrate.from_payload``; a test that needs the raw integer builds it with this.
INV_LENGTH_SCALE = 1 << 32


def make_synthetic_payload() -> tuple[AccumulatorPayload, RegionArrays]:
    """A 1-reference, 3-node payload + aligned :class:`RegionArrays`, with every bank distinct.

    chr1 is cut at 0/100/200/300, so it owns **3 nodes and 2 contiguous edges** — the axes are off by one
    per reference, and a fixture that used the same length for both would hide an axis mix-up. One
    junction edge exists so the third axis is non-trivial.

    Nodes: n0 +exon, n1 −exon, n2 intergenic. Every population gets its own values, and no two banks
    share a value, so a consumer reading the wrong one cannot pass by coincidence::

        node_contained_count   n0 [10, 2]   n1 [1, 20]   n2 [7, 8]
        node_spanning_count    n0 [3, 0]    n1 [0, 5]    n2 [1, 6]
        edge_unspliced_count   e0 [4, 1]    e1 [2, 3]
        edge_spliced_count     e0 [0, 0]    e1 [6, 0]
        sj_count               j0 [9, 4]
    """
    n_nodes, n_edges, n_sj = 3, 2, 1

    def bank(rows, values, dtype):
        return np.asarray(values, dtype=dtype).reshape(rows, 2)

    contained = bank(n_nodes, [[10, 2], [1, 20], [7, 8]], np.uint32)
    spanning = bank(n_nodes, [[3, 0], [0, 5], [1, 6]], np.uint32)
    unspliced = bank(n_edges, [[4, 1], [2, 3]], np.uint32)
    spliced = bank(n_edges, [[0, 0], [6, 0]], np.uint32)
    sj = bank(n_sj, [[9, 4]], np.uint32)

    def inv(counts, placements):
        """The fixed-point sum a bank of ``counts`` fragments at one placement count would deposit.

        ⚠ Rounds HALF AWAY FROM ZERO, exactly as ``_accumulator_reference.inv_length_quantum`` does —
        not floor. A fixture that quantised differently from the accumulator would make the decode test
        assert the fixture's own arithmetic rather than the substrate's.
        """
        quantum = (2 * INV_LENGTH_SCALE + placements) // (2 * placements)
        return (np.asarray(counts, np.uint64) * np.uint64(quantum)).astype(np.uint64)

    def lengths(counts, length):
        return (np.asarray(counts, np.uint64) * np.uint64(length)).astype(np.uint64)

    payload = AccumulatorPayload(
        cut_positions=np.array([0, 100, 200, 300], dtype=np.int64),
        ref_cut_offsets=np.array([0, 4], dtype=np.int64),
        ref_node_offsets=np.array([0, n_nodes], dtype=np.int64),
        ref_edge_offsets=np.array([0, n_edges], dtype=np.int64),
        ref_sj_offsets=np.array([0, n_sj], dtype=np.int64),
        node_contained_count=contained,
        node_contained_inv_length_sum=inv(contained, 50),
        node_contained_length_sum=lengths(contained, 50),
        node_spanning_count=spanning,
        node_spanning_inv_length_sum=inv(spanning, 40),
        node_spanning_length_sum=lengths(spanning, 40),
        node_start_count=np.array([11, 12, 13], dtype=np.uint32),
        edge_unspliced_count=unspliced,
        edge_unspliced_inv_length_sum=inv(unspliced, 25),
        edge_unspliced_length_sum=lengths(unspliced, 26),
        edge_spliced_count=spliced,
        edge_spliced_inv_length_sum=inv(spliced, 20),
        edge_spliced_length_sum=lengths(spliced, 21),
        sj_count=sj,
        sj_inv_length_sum=inv(sj, 10),
        sj_length_sum=lengths(sj, 11),
        pool_lengths=np.zeros((N_FRAGMENT_POOLS, 201), dtype=np.int64),
        qc=ScanQC(
            deposited=36,
            dropped_too_long=0,
            dropped_empty=0,
            dropped_strand_undefined=0,
            dropped_ambiguous_path=0,
            unannotated_introns=0,
            contradictory_sj_strand=0,
            sj_implicit_fragments=0,
            introns_absorbed=0,
        ),
        max_length=200,
        n_refs=1,
    )
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(3, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * 3, dtype="string"),
            "start": np.array([0, 100, 200], dtype=np.int64),
            "end": np.array([100, 200, 300], dtype=np.int64),
            "length": np.array([100, 100, 100], dtype=np.int64),
            "signature": np.array([BIT_EXON_POS, BIT_EXON_NEG, 0], dtype=np.uint8),
        }
    )
    return payload, RegionArrays.from_frame(region_df, {"chr1": 0})


def make_gdna_fl_pmf(mean: int = 50, max_size: int = 200) -> np.ndarray:
    """A spiked gDNA fragment-length pmf for calibration tests (μ_FL = ``mean``)."""
    pmf = np.zeros(max_size + 1, dtype=np.float64)
    pmf[mean] = 1.0
    return pmf


def make_strand_models(p_r1_sense: float, n_observations: int, n_junctions: int = 1):
    """A real :class:`StrandModels` with a chosen κ and observation count.

    The calibrator now reads BOTH halves of the RNA strand Beta-Binomial from the per-junction SJ
    strand table — κ as its marginal, the overdispersion as its spread — so a unit fixture must
    supply a real table rather than duck-type two scalars. Observations are spread evenly over
    ``n_junctions`` motif-POS junctions, sense/antisense split to give exactly ``p_r1_sense``.
    """
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels
    from rigel.types import Strand

    per = n_observations // n_junctions
    rem = n_observations - per * n_junctions
    depth = np.full(n_junctions, per, dtype=np.int64)
    if n_junctions:
        depth[0] += rem
    n_sense = np.rint(depth * float(p_r1_sense)).astype(np.int64)
    # Repair rounding so the marginal is EXACTLY the requested rate where it can be.
    want = int(round(n_observations * float(p_r1_sense)))
    if n_junctions and n_sense.sum() != want:
        n_sense[0] = np.clip(n_sense[0] + (want - n_sense.sum()), 0, depth[0])
    table = SJStrandTable(
        ref_id=np.zeros(n_junctions, dtype=np.int32),
        start=np.arange(n_junctions, dtype=np.int64) * 1000,
        end=np.arange(n_junctions, dtype=np.int64) * 1000 + 100,
        motif_strand=np.full(n_junctions, int(Strand.POS), dtype=np.int8),
        n_sense=n_sense,
        n_antisense=depth - n_sense,
    )
    return StrandModels(exonic_spliced=StrandModel.from_sj_table(table))
