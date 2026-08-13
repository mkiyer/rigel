"""Hand-built synthetic payloads for calibration unit tests.

A plain helper module (not ``conftest.py``) to avoid colliding with the top-level ``conftest`` that
several tests import by name.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    AccumulatorPayload,
    DeferredFragments,
    GapCensus,
    ScanQC,
)

#: ⭐ ONE NUMERIC CONVENTION: a COUNT is an integer, a FRACTION is float64. The accumulator deposits
#: ``1/placements`` directly — there is no scale and nothing to decode, so this fixture builds the same
#: number the accumulator would.


def make_synthetic_payload() -> tuple[AccumulatorPayload, RegionArrays]:
    """A 1-reference, 3-node payload + aligned :class:`RegionArrays`, with every bank distinct.

    chr1 is cut at 0/100/200/300, so it owns **3 nodes and 2 contiguous edges** — the axes are off by one
    per reference, and a fixture that used the same length for both would hide an axis mix-up. One
    junction edge exists so the third axis is non-trivial.

    Nodes: n0 +exon, n1 −exon, n2 intergenic. Every population gets its own values, and no two banks
    share a value, so a consumer reading the wrong one cannot pass by coincidence::

        node_contained_count   n0 [10, 2]   n1 [1, 20]   n2 [7, 8]
        edge_unspliced_count   e0 [4, 1]    e1 [2, 3]
        edge_spliced_count     e0 [0, 0]    e1 [6, 0]
        sj_count               j0 [9, 4]
    """
    n_nodes, n_edges, n_sj = 3, 2, 1

    def bank(rows, values, dtype):
        return np.asarray(values, dtype=dtype).reshape(rows, 2)

    contained = bank(n_nodes, [[10, 2], [1, 20], [7, 8]], np.uint32)
    unspliced = bank(n_edges, [[4, 1], [2, 3]], np.uint32)
    spliced = bank(n_edges, [[0, 0], [6, 0]], np.uint32)
    sj = bank(n_sj, [[9, 4]], np.uint32)

    def inv(counts, placements):
        """The fixed-point sum a bank of ``counts`` fragments at one placement count would deposit.

        ⚠ ONE column — the two strands are SUMMED, because the length moments carry no strand axis.
        """
        return np.asarray(counts, np.float64).sum(axis=1) / placements

    def lengths(counts, length):
        """⚠ ONE column, for the same reason :func:`inv` is."""
        return (np.asarray(counts, np.uint64).sum(axis=1) * np.uint64(length)).astype(np.uint64)

    def mass(counts, per_crossing=2):
        """The conserved mass a bank of ``counts`` crossings would deposit, at ``1/per_crossing`` each.

        ⚠ ONE value per edge — the two strand columns are SUMMED, because the mass has no strand axis.
        ⭐ And it is a plausible state rather than an arbitrary array: a crossing's share is at most one
        fragment, so ``mass <= count`` must hold. ``per_crossing = 2`` is the value for a fragment that
        crosses two lines, which is what the multi-line geometry this fixture describes produces.
        """
        return np.asarray(counts, np.float64).sum(axis=1) / per_crossing

    payload = AccumulatorPayload(
        cut_positions=np.array([0, 100, 200, 300], dtype=np.int64),
        ref_cut_offsets=np.array([0, 4], dtype=np.int64),
        ref_node_offsets=np.array([0, n_nodes], dtype=np.int64),
        ref_edge_offsets=np.array([0, n_edges], dtype=np.int64),
        ref_sj_offsets=np.array([0, n_sj], dtype=np.int64),
        node_contained_count=contained,
        node_contained_inv_opportunity_sum=inv(contained, 50),
        node_start_count=np.array([11, 12, 13], dtype=np.uint32),
        edge_unspliced_count=unspliced,
        edge_unspliced_inv_length_sum=inv(unspliced, 25),
        edge_unspliced_mass=mass(unspliced),
        edge_spliced_count=spliced,
        edge_spliced_mass=mass(spliced),
        sj_count=sj,
        sj_inv_length_sum=inv(sj, 10),
        # ⭐ TWO COLUMNS, and they SUM to what the one-column bank held (1.3), so every
        # expectation downstream of the substrate's fold is unchanged by the strand split.
        sj_mass=np.asarray(sj, np.float64) / 10.0,
        pool_lengths=np.zeros((N_FRAGMENT_POOLS, 201), dtype=np.int64),
        deposited_lengths=np.zeros(201, dtype=np.uint32),
        # ⚠ Nothing was deferred here, and that is a real state, not a stub: this fixture has no
        # annotated intron in any mate gap. The two empty spellings live on the classes so a
        # hand-built payload cannot get the `[0]`-not-`[]` offset boundary wrong.
        deferred=DeferredFragments.empty(),
        gap_resolution=GapCensus.zeros(),
        qc=ScanQC(
            deposited=36,
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


def make_synthetic_junctions():
    """The :class:`JunctionGeometry` matching :func:`make_synthetic_payload`'s one junction edge.

    ⚠ **It must exist, and match.** The payload declares ``n_sj = 1``; ``calibrate`` refuses a junction
    axis of a different length, because an axis addressing a different graph would place every splice
    on the wrong line and nothing downstream would fault on it.

    The junction runs ``node 0 → node 2``, i.e. it splices OVER node 1 — the only shape a 3-node
    reference admits, and the one that makes the donor and acceptor two DIFFERENT lines (edge 0 and
    edge 1). A ``0 → 1`` junction would put both endpoints on the same line and hide an
    endpoint mix-up.
    """
    from rigel.calibration.splice_graph import JunctionGeometry
    from rigel.types import Strand

    return JunctionGeometry(
        src_node=np.array([0], dtype=np.int64),
        dst_node=np.array([2], dtype=np.int64),
        strand=np.array([int(Strand.POS)], dtype=np.int8),
        reach_lo=np.array([100.0]),
        reach_hi=np.array([100.0]),
    )


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


# ---------------------------------------------------------------------------
# The S5.e chain fixture — hand-built numbers on the node / edge / junction axes.
# ---------------------------------------------------------------------------


def delta_pmf(length: int, size: int | None = None) -> np.ndarray:
    """A fragment-length pmf concentrated at one length — the simplest oracle for a divisor."""
    p = np.zeros((size or length) + 1, dtype=np.float64)
    p[length] = 1.0
    return p


def make_chain_parts(
    signatures,
    *,
    node_size_bp=1000.0,
    node_pos=0.0,
    node_neg=0.0,
    edge_pos=0.0,
    edge_neg=0.0,
    edge_spliced=0.0,
    junctions=None,
    gdna_fl=None,
    rna_fl=None,
    ref_names=None,
):
    """A chain + substrate + geometry + statics over ``signatures``, on the S5.e axes.

    ⭐ **The axes are off by one per reference and that is the point of the helper**: a reference with
    ``k`` nodes owns ``k`` node rows and ``k − 1`` contiguous-edge rows, with **no terminal slots**. Every
    per-object argument is broadcast, so a test states only the numbers it cares about.

    ``junctions`` is a list of ``(src_node, dst_node, strand, reach_lo, reach_hi, count)``; each becomes a
    row on the junction axis and is placed on the lines it leaves and enters.

    Returns ``SimpleNamespace(chain, substrate, region_arrays, geometry, statics)``.
    """
    from types import SimpleNamespace

    from rigel.calibration.node_chain import build_node_chain
    from rigel.calibration.node_geometry import build_node_geometry, build_node_statics
    from rigel.calibration.signature import transcript_strand_class
    from rigel.calibration.splice_graph import JunctionGeometry

    sig = np.asarray(signatures, dtype=np.uint8)
    n_nodes = sig.shape[0]
    ref_names = list(ref_names or ["chr1"] * n_nodes)
    ref_id = np.array(
        [{n: i for i, n in enumerate(dict.fromkeys(ref_names))}[r] for r in ref_names]
    )
    rno = np.zeros(int(ref_id.max()) + 2, dtype=np.int64)
    np.cumsum(np.bincount(ref_id, minlength=rno.shape[0] - 1), out=rno[1:])
    reo = np.zeros_like(rno)
    np.cumsum(np.maximum(np.diff(rno) - 1, 0), out=reo[1:])
    n_edges = int(reo[-1])

    def pair(pos, neg, n):
        return np.stack(
            [
                np.broadcast_to(np.asarray(pos, float), (n,)),
                np.broadcast_to(np.asarray(neg, float), (n,)),
            ],
            axis=1,
        ).copy()

    j = list(junctions or [])
    junction_geometry = JunctionGeometry(
        src_node=np.array([x[0] for x in j], dtype=np.int64),
        dst_node=np.array([x[1] for x in j], dtype=np.int64),
        strand=np.array([x[2] for x in j], dtype=np.int8),
        reach_lo=np.array([float(x[3]) for x in j]),
        reach_hi=np.array([float(x[4]) for x in j]),
    )
    gdna_pmf = delta_pmf(200, 400) if gdna_fl is None else gdna_fl
    rna_pmf = delta_pmf(100, 400) if rna_fl is None else rna_fl
    node_sz = np.broadcast_to(np.asarray(node_size_bp, float), (n_nodes,)).copy()

    node_counts = pair(node_pos, node_neg, n_nodes)
    edge_counts = pair(edge_pos, edge_neg, n_edges)
    substrate = SimpleNamespace(
        node_contained=SimpleNamespace(count=node_counts),
        edge_unspliced=SimpleNamespace(count=edge_counts),
        edge_spliced=SimpleNamespace(count=pair(edge_spliced, 0.0, n_edges)),
        junction=SimpleNamespace(
            count=np.array([[float(x[5]), 0.0] for x in j]).reshape(len(j), 2)
        ),
    )
    region_arrays = SimpleNamespace(
        signature=sig,
        strand_class=transcript_strand_class(sig.astype(np.int64)),
        region_size_bp=node_sz,
        ref_id=ref_id,
        ref_offsets=rno,
        n_regions=n_nodes,
    )
    chain = build_node_chain(rno, reo)
    geometry = build_node_geometry(
        chain,
        substrate,
        region_arrays,
        junction_geometry,
        gdna_pmf,
        rna_pmf,
    )
    return SimpleNamespace(
        chain=chain,
        substrate=substrate,
        region_arrays=region_arrays,
        geometry=geometry,
        statics=build_node_statics(chain, region_arrays),
    )
