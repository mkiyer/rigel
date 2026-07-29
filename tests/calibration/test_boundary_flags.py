"""The splice graph's structural flags, re-indexed onto the ACCUMULATOR's boundary axis (plan W1c).

    Plan: docs/accumulator/06_implementation_plan.md §3 (W1c)   ·   Graph: docs/index/00_splice_graph_design.md §2.3

W2 needs to ask, at a seam, *"is this a transcript terminus?"* and *"is this a splice site, on which
flank?"* — questions the 4-bit signature is structurally blind to. The graph answers them, but in its
own index space, and **the two spaces are off by one per reference in opposite directions**: a
reference with ``k`` nodes has ``k − 1`` contiguous edges and ``k + 1`` accumulator boundary slots.

Getting that mapping wrong shifts every flag by one seam — which is invisible in aggregate, would
sail through a bit-identity gate (nothing reads the flags yet), and would silently corrupt W2. So the
mapping is what this module tests, and it tests it against a **really scanned payload**, not against
a re-derivation of the same arithmetic.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.splice_graph import (
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_POS,
    FLAG_TES_POS,
    FLAG_TSS_POS,
    build_boundary_flags_array,
    build_node_partition_arrays,
    is_splice_site,
    is_terminus,
)
from rigel.types import Strand

from conftest import build_test_index


#: t0 splices [400,700); t1 ENDS at 700, where t0's intron ends — so position 700 is an ACCEPTOR
#: *and* a TES, the case §2.3 exists for. t2 on chr2 keeps the per-reference offsets honest.
GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t651\t700\t.\t+\t.\tgene_id "g2"; transcript_id "t1";
chr2\ttest\texon\t301\t500\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
chr2\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g3"; transcript_id "t2";
"""

REFS = {"chr1": 1500, "chr2": 1500}


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF, name="w1c_flags", refs=REFS)


def _slot_positions(index):
    """Genomic position of every accumulator boundary slot, in slot order."""
    positions, _offsets, _types = build_node_partition_arrays(index)
    return positions


def test_flags_are_one_per_accumulator_boundary_slot(index):
    """The array must address the payload's boundary axis with a single index — same length, same
    per-reference offsets as the cut array the scanner is handed."""
    flags = build_boundary_flags_array(index)
    positions, offsets, _types = build_node_partition_arrays(index)
    assert flags.dtype == np.uint16
    assert flags.shape == positions.shape
    assert int(offsets[-1]) == flags.shape[0]


def test_reference_terminals_are_zero(index):
    """⚠ Not padding — a terminal is genuinely flag-less. It is not a transition, no edge exists
    there, and it never carried a deposit either."""
    flags = build_boundary_flags_array(index)
    _positions, offsets, _types = build_node_partition_arrays(index)
    for f in range(len(index.ref_names)):
        lo, hi = int(offsets[f]), int(offsets[f + 1])
        if hi == lo:
            continue
        assert flags[lo] == 0, "the reference-start terminal carries a flag"
        assert flags[hi - 1] == 0, "the reference-end terminal carries a flag"
        assert flags[lo + 1 : hi - 1].any(), "the fixture must set SOME interior flag"


def test_each_slot_carries_the_flags_of_its_own_genomic_position(index):
    """⭐ THE mapping assertion, checked position by position rather than index by index.

    Every contiguous edge sits at a known genomic coordinate; every boundary slot sits at a known
    genomic coordinate. Matching them up by COORDINATE cannot be fooled by an off-by-one in the
    index arithmetic, which is the failure this whole module exists to catch.
    """
    flags = build_boundary_flags_array(index)
    positions = _slot_positions(index)
    nodes, edges = index.nodes_df, index.edges_df

    contiguous = edges[edges["kind"] == 0]
    end_of_src = nodes["end"].to_numpy(np.int64)[contiguous["src"].to_numpy(np.int64)]
    want = dict(zip(end_of_src.tolist(), contiguous["flags"].to_numpy(np.uint16).tolist()))

    got: dict[int, int] = {}
    _positions, offsets, _types = build_node_partition_arrays(index)
    for f in range(len(index.ref_names)):
        lo, hi = int(offsets[f]), int(offsets[f + 1])
        for b in range(lo + 1, hi - 1):  # interior slots only; terminals are covered above
            got[int(positions[b])] = int(flags[b])
    assert got == want


def test_the_terminus_and_splice_site_predicates(index):
    """Position 700 is BOTH a TES (t1 ends) and an ACCEPTOR (t0's intron ends) — the case the 4-bit
    signature cannot represent, and the reason the two predicates are independent bits."""
    flags = build_boundary_flags_array(index)
    positions = _slot_positions(index)
    at = {int(p): int(f) for p, f in zip(positions, flags)}

    assert at[700] & FLAG_TES_POS and at[700] & FLAG_ACCEPTOR_POS
    assert at[400] & FLAG_DONOR_POS  # t0's intron starts
    assert at[200] & FLAG_TSS_POS  # t0's 5' end

    arr = np.array([at[700], at[400], at[200]], dtype=np.uint16)
    np.testing.assert_array_equal(is_terminus(arr, Strand.POS), [True, False, True])
    np.testing.assert_array_equal(is_splice_site(arr, Strand.POS), [True, True, False])
    # the fixture is + only, so nothing may leak onto the − strand
    assert not is_terminus(arr, Strand.NEG).any()
    assert not is_splice_site(arr, Strand.NEG).any()


def test_flags_align_with_a_REALLY_SCANNED_payload(tmp_path):
    """⭐ The end-to-end statement: the array indexes the payload the C++ accumulator actually
    produced, not a re-derivation of the arithmetic that produced it.

    ``NodeStatics`` then places each slot's bits on the right chain node — asserted through the
    chain rather than assumed, because that is the second place an off-by-one could hide.
    """
    from rigel.calibration.node_chain import BOUNDARY, build_node_chain
    from rigel.calibration.node_geometry import build_node_statics
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer
    from rigel.sim import ReadSimConfig, Scenario

    sc = Scenario("w1c_flags", genome_length=5000, seed=11, work_dir=tmp_path / "w1c")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(200, 400), (600, 800)], "abundance": 50}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(2500, 2700), (3000, 3200)], "abundance": 30}])
    result = sc.build_oracle(
        n_fragments=200,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=11,
        ),
    )
    index = result.index
    _stats, _sm, _flm, _buf, payload = scan_and_buffer(
        str(result.bam_path), index, BamScanConfig(sj_strand_tag="auto")
    )

    flags = build_boundary_flags_array(index)
    assert flags.shape[0] == int(payload.ref_boundary_offsets[-1]), (
        "the flags array does not address the payload's boundary axis"
    )

    ra = RegionArrays.from_index(index)
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    statics = build_node_statics(
        chain,
        CalibrationSubstrate.from_payload(payload, ra),
        BoundarySubstrate.from_payload(payload),
        ra,
        flags,
    )
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    assert not statics.boundary_flags[kind != BOUNDARY].any(), (
        "a region node carries boundary flags"
    )
    np.testing.assert_array_equal(
        statics.boundary_flags[kind == BOUNDARY], flags[idx[kind == BOUNDARY]]
    )
    assert statics.boundary_flags.any(), "the scenario must set SOME flag, or this asserts nothing"
    sc.cleanup()


def test_a_wrong_length_is_refused(index):
    """A silently mis-sized array would shift every flag by one seam — invisible in aggregate, and
    exactly what W2 must not inherit."""
    from rigel.calibration.node_chain import build_node_chain
    from rigel.calibration.node_geometry import build_node_statics
    from rigel.calibration.region_arrays import RegionArrays

    ra = RegionArrays.from_index(index)
    rro = np.asarray(ra.ref_offsets, np.int64)
    rbo = rro + np.arange(rro.shape[0], dtype=np.int64)
    chain = build_node_chain(rro, rbo)

    class _Sub:
        pass

    with pytest.raises(ValueError, match="boundary slot"):
        build_node_statics.__wrapped__ if False else None
        build_node_statics(chain, _Sub(), _Sub(), ra, np.zeros(int(rbo[-1]) + 1, dtype=np.uint16))
