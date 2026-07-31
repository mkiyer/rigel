"""Substrate conservation on a REAL scanned payload — the three axes, end to end.

⛔ **This file was rewritten, not ported (S5.f).** Its predecessor asserted that the region→boundary
map "attributes exactly the non-terminal boundary mass, no double count, no loss" — a property of the
``k + 1`` boundary axis, whose two data-free terminal slots were the only thing that could lose mass.
That axis does not exist: a reference with ``k`` nodes owns ``k − 1`` lines and **every one of them has
a node on both sides**, so there is no terminal to drop and no face to double. A test of a model that
no longer exists can only be rewritten from scratch (the call S5.c made on ``test_message_frames.py``
and S5.e made on three ``test_bp_solver`` tests).

⭐ **What is kept is the thing that made it worth having: it runs on a REAL SCAN**, not a hand-built
fixture, so the payload, the index and the geometry must agree with each other rather than with a
fixture author. The invariants are re-derived from the FRAGMENT BUFFER and the INDEX — sources
independent of the accumulator (`CARRY_FORWARD.md` §3 trap 1: a validator that calls the builder's own
helper validates nothing).
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays, edge_node_indices
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario

SEED = 4321


@pytest.fixture(scope="module")
def scanned():
    import tempfile
    from pathlib import Path

    work = Path(tempfile.mkdtemp())
    sc = Scenario("subcons", genome_length=6000, seed=SEED, work_dir=work / "subcons")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(300, 600), (900, 1200)], "abundance": 60}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3300), (3700, 4000)], "abundance": 40}])
    # ⚠ gDNA is not decoration here. Without it every fragment is mature RNA, which either sits inside
    # one exon node or JUMPS — so ``edge_unspliced``, the one population calibration actually
    # deconvolves, is identically zero and a "conservation" test would pass over an empty bank.
    from rigel.sim.reads import GDNAConfig

    result = sc.build_oracle(
        n_fragments=300,
        sim_config=ReadSimConfig(frag_mean=200, frag_std=30, seed=SEED),
        gdna_config=GDNAConfig(abundance=400.0, frag_mean=250, frag_std=60),
    )
    config = PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="auto"))
    _, _, _, buffer, payload = scan_and_buffer(str(result.bam_path), result.index, config.scan)
    ra = RegionArrays.from_frame(result.index.nodes_df, result.index.ref_name_to_id)
    yield payload, ra, buffer, result.index
    sc.cleanup()


def test_one_node_start_per_ACCEPTED_fragment(scanned):
    """⭐ THE invariant. ``node_start_count`` is incremented once, at the node holding the fragment's
    first base, for every fragment the accumulator accepts — so its total IS the accepted count.

    ⚠ Checked against the payload's own QC tally, which is written on the SAME line of the deposit but
    is a separate counter: they can only agree if every accepted fragment reached both.
    """
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    assert int(np.asarray(sub.node_start_count).sum()) == int(payload.qc.deposited)
    assert int(payload.qc.deposited) > 0, "a scan that deposited nothing proves nothing"


def test_every_buffered_fragment_is_ACCEPTED_or_DROPPED_FOR_A_NAMED_REASON(scanned):
    """The accepted count plus every drop reason must account for the whole buffer — an independent
    source, filled by the scanner's own fragment grouping rather than by the accumulator."""
    payload, _ra, buffer, _index = scanned
    qc = payload.qc
    dropped = (
        qc.dropped_too_long
        + qc.dropped_empty
        + qc.dropped_strand_undefined
        + qc.dropped_ambiguous_path
    )
    assert qc.deposited + dropped == buffer.total_fragments


def test_the_edge_axis_is_N_minus_the_NONEMPTY_references(scanned):
    """``E = N − (references that own at least one node)``, re-derived from the INDEX rather than from
    the payload's own offsets. A reference with one node owns no line; a reference with none owns
    nothing at all — neither is a special case in the formula, and both were terminal slots before."""
    payload, ra, _buffer, index = scanned
    nodes_per_ref = np.diff(np.asarray(ra.ref_offsets, dtype=np.int64))
    expected_edges = int(np.maximum(nodes_per_ref - 1, 0).sum())
    assert int(payload.n_edges) == expected_edges
    assert int(payload.n_nodes) == len(index.nodes_df)


def test_contained_deposits_never_exceed_the_accepted_fragments(scanned):
    """A fragment lies wholly inside at most one node, so the contained bank cannot out-count the
    fragments. ⚠ It is an inequality, not an identity: a fragment that crosses any line is contained
    in nothing and contributes only to the edge banks."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    contained = int(np.asarray(sub.node_contained.count).sum())
    assert 0 < contained <= int(np.asarray(sub.node_start_count).sum())


def test_no_line_STRADDLES_a_reference(scanned):
    """Every contiguous edge joins two nodes of the same reference — the invariant the ``k + 1`` axis
    could not state, because its terminal slots had a node on one side only."""
    _payload, ra, _buffer, _index = scanned
    lo, hi = edge_node_indices(np.asarray(ra.ref_id))
    ref = np.asarray(ra.ref_id)
    np.testing.assert_array_equal(ref[lo], ref[hi])
    np.testing.assert_array_equal(hi, lo + 1)


def test_the_two_edge_banks_are_DISJOINT_populations(scanned):
    """``edge_unspliced`` and ``edge_spliced`` are different molecules at the same line — crossed
    contiguously having spliced NOWHERE, versus having spliced ELSEWHERE — so a fragment lands in
    exactly one of them and neither is a subset of the other. Pinned here on real data as a
    non-degeneracy check: if the scan produced both, the two banks cannot be the same array."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    unspliced = np.asarray(sub.edge_unspliced.count)
    spliced = np.asarray(sub.edge_spliced.count)
    assert unspliced.sum() > 0, "the fixture must exercise the unspliced bank"
    assert not np.array_equal(unspliced, spliced)


def test_the_junction_axis_matches_the_payloads_own(scanned):
    """The substrate's junction population is exactly ``n_sj`` rows — the third axis, independent of
    the other two, and the one a consumer must not size from ``n_nodes`` or ``n_edges``."""
    payload, ra, _buffer, _index = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    assert sub.n_junctions == int(payload.n_sj)
    assert np.asarray(sub.junction.count).shape == (int(payload.n_sj), 2)
    assert np.asarray(sub.junction.count).sum() > 0, "the fixture must exercise the junction axis"
