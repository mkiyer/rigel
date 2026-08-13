"""The splice graph's structural flags on the accumulator's CONTIGUOUS-EDGE axis.

    Graph: · (S5.e)

The solver needs to ask, at a line, *"is this a transcript terminus?"* and *"is this a splice site, on
which flank?"* — questions the 4-bit signature is structurally blind to. The graph answers them, and
until S5.e it answered them in a **different index space**: a reference with ``k`` regions has ``k − 1``
contiguous edges but the old accumulator had ``k + 1`` boundary slots, so the two ran off by one per
reference *in opposite directions*.

⭐ **That mismatch is gone.** There are no terminal slots — a contiguous edge is the line BETWEEN two
adjacent regions, and there is no such line before the first or after the last — so the flags array and
the payload's edge axis are the same axis, with no padding to align. What is left to test is that the
flags land on the right LINE, and it is tested by **genomic coordinate** against a really scanned
payload, never by re-deriving the same index arithmetic.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.splice_graph import (
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_POS,
    FLAG_TES_POS,
    FLAG_TSS_POS,
    build_edge_flags_array,
    build_region_partition_arrays,
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
    return build_test_index(tmp_path_factory, GTF, name="edge_flags", refs=REFS)


def _edge_positions(index) -> np.ndarray:
    """Genomic position of every contiguous edge, in edge order.

    A reference contributing ``c`` cuts owns ``c − 1`` regions and ``c − 2`` interior lines, and line
    ``e`` sits at cut ``e + 1`` — so the interior cuts, per reference, ARE the edge coordinates.
    """
    positions, cut_offsets, _types = build_region_partition_arrays(index)
    out = []
    for f in range(len(cut_offsets) - 1):
        lo, hi = int(cut_offsets[f]), int(cut_offsets[f + 1])
        if hi - lo >= 2:
            out.append(positions[lo + 1 : hi - 1])
    return np.concatenate(out) if out else np.zeros(0, np.int64)


def test_there_is_EXACTLY_ONE_ENTRY_PER_LINE_and_no_padding(index):
    """⭐ The shape change S5.d/S5.e make. A reference with ``k`` regions contributes ``k − 1`` entries,
    never ``k + 1``: the two data-free terminals existed only so every region had an object on each
    side, and an edge is defined by having a region on BOTH sides."""
    flags = build_edge_flags_array(index)
    assert flags.dtype == np.uint16
    n_regions = len(index.regions_df)
    n_refs_with_regions = index.regions_df["ref_name"].nunique()
    assert flags.shape[0] == n_regions - n_refs_with_regions
    assert flags.shape == _edge_positions(index).shape


def test_each_entry_carries_the_flags_of_its_own_GENOMIC_POSITION(index):
    """⭐ THE mapping assertion, checked position by position rather than index by index — matching by
    COORDINATE cannot be fooled by an off-by-one in the index arithmetic, which is the failure this
    module exists to catch."""
    flags = build_edge_flags_array(index)
    positions = _edge_positions(index)
    regions, edges = index.regions_df, index.edges_df

    contiguous = edges[edges["kind"] == 0]
    end_of_src = regions["end"].to_numpy(np.int64)[contiguous["src"].to_numpy(np.int64)]
    want = dict(zip(end_of_src.tolist(), contiguous["flags"].to_numpy(np.uint16).tolist()))
    got = {int(p): int(f) for p, f in zip(positions, flags)}
    assert got == want
    assert any(want.values()), "the fixture must set SOME flag, or this asserts nothing"


def test_the_terminus_and_splice_site_predicates(index):
    """Position 700 is BOTH a TES (t1 ends) and an ACCEPTOR (t0's intron ends) — the case the 4-bit
    signature cannot represent, and the reason the two predicates are independent bits."""
    flags = build_edge_flags_array(index)
    at = {int(p): int(f) for p, f in zip(_edge_positions(index), flags)}

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

    ``RegionStatics`` then places each edge's bits on the right chain SLOT — asserted through the chain
    rather than assumed, because that is the second place an off-by-one could hide.
    """
    from rigel.calibration.region_chain import EDGE, build_region_chain
    from rigel.calibration.region_geometry import build_region_statics
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import BamScanConfig
    from rigel.pipeline import scan_and_buffer
    from rigel.sim import ReadSimConfig, Scenario

    sc = Scenario("edge_flags", genome_length=5000, seed=11, work_dir=tmp_path / "ef")
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
    _stats, _sm, _buf, payload = scan_and_buffer(
        str(result.bam_path), index, BamScanConfig(sj_strand_tag="auto")
    )

    flags = build_edge_flags_array(index)
    assert flags.shape[0] == int(payload.ref_edge_offsets[-1]), (
        "the flags array does not address the payload's contiguous-edge axis"
    )

    ra = RegionArrays.from_index(index)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_edge_offsets)
    statics = build_region_statics(chain, ra, flags)
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, np.int64)
    assert not statics.edge_flags[kind != EDGE].any(), "a REGION slot carries edge flags"
    np.testing.assert_array_equal(statics.edge_flags[kind == EDGE], flags[idx[kind == EDGE]])
    assert statics.edge_flags.any(), "the scenario must set SOME flag, or this asserts nothing"
    sc.cleanup()


def test_a_wrong_length_is_refused(index):
    """A silently mis-sized array would shift every flag by one line — invisible in aggregate, and
    exactly what the solver must not inherit. ⚠ The old ``k+1`` shape is the specific wrong length
    most likely to be handed in during the transition, so that is the one used here."""
    from rigel.calibration.region_chain import build_region_chain
    from rigel.calibration.region_geometry import build_region_statics
    from rigel.calibration.region_arrays import RegionArrays

    ra = RegionArrays.from_index(index)
    rno = np.asarray(ra.ref_offsets, np.int64)
    reo = np.zeros_like(rno)
    np.cumsum(np.maximum(np.diff(rno) - 1, 0), out=reo[1:])
    chain = build_region_chain(rno, reo)

    class _Sub:
        pass

    old_shape = rno + np.arange(rno.shape[0], dtype=np.int64)  # the retired k+1 boundary axis
    with pytest.raises(ValueError, match="one per contiguous edge"):
        build_region_statics(chain, ra, np.zeros(int(old_shape[-1]), dtype=np.uint16))
