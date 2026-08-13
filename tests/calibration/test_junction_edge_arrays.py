"""The junction boundaries, re-indexed onto the accumulator's flat CUT axis as a CSR.



The deposit rule must answer one question per observed intron, inside the BAM-scan hot loop: *is this
intron an annotated junction, and if so which boundary?* Getting it wrong is not a rounding error — an
intron that fails to match deposits nothing at all (an unannotated intron credits no junction boundary), so
a mis-keyed table silently deletes the entire spliced-RNA signal.

The table is keyed by the **donor cut index**, which is the index the deposit has already computed while
locating the lines its path crosses. ``_lookup`` below is written to mirror that inner loop exactly, so
these tests double as the API's worked example.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.splice_graph import (
    EDGE_KIND_JUNCTION,
    build_junction_edge_arrays,
    build_region_partition_arrays,
)
from rigel.types import Strand

from conftest import build_test_index
from native._accumulator_reference import Partition


#: chr1 t0 splices 400->700 and 900->1200; t1 shares t0's DONOR at 400 but lands at 1000, so cut 400
#: has fan-out 2 — the alternative-3'-splice-site case a naive one-junction-per-cut table would drop.
#: chr2 carries a NEG-strand junction, which keeps the per-reference offsets and the strand honest.
#: ⭐ chr3 is a NESTED pair with the OUTER intron on the minus strand: intron [400,1400) NEG encloses
#: [600,1000) POS. It is the fixture that makes the slot ordering falsifiable. Everywhere else in this
#: module — and on both real indexes — donor order, acceptor order and strand order happen to agree, so
#: every permutation of the sort key produces the identical answer and the contract test proves nothing.
#: Nesting breaks donor-vs-acceptor (400 < 600 but 1400 > 1000) and the strand assignment breaks
#: donor-vs-strand (the smaller donor is the NEG one). Verified: three separate key permutations each turn
#: ``test_the_csr_slot_order_matches_the_reference_accumulator`` red only because chr3 is here.
GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t1201\t1400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr2\ttest\texon\t301\t500\t.\t-\t.\tgene_id "g2"; transcript_id "t2";
chr2\ttest\texon\t801\t1000\t.\t-\t.\tgene_id "g2"; transcript_id "t2";
chr3\ttest\texon\t201\t400\t.\t-\t.\tgene_id "g5"; transcript_id "t5";
chr3\ttest\texon\t1401\t1600\t.\t-\t.\tgene_id "g5"; transcript_id "t5";
chr3\ttest\texon\t201\t600\t.\t+\t.\tgene_id "g6"; transcript_id "t6";
chr3\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "g6"; transcript_id "t6";
"""

REFS = {"chr1": 2000, "chr2": 2000, "chr3": 2000}

#: ⭐ A STRAND-COINCIDENT PAIR: two genes on opposite strands whose intron coordinates are byte-identical.
#: It is the ONLY configuration in which the two builders' slot orderings can differ, so it is the
#: discriminating case for the ordering contract — and the index warns about it, correctly, which is why it
#: gets its own fixture instead of polluting every test above with the warning.
COINCIDENT_GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g3"; transcript_id "t3";
chr1\ttest\texon\t801\t1000\t.\t+\t.\tgene_id "g3"; transcript_id "t3";
chr1\ttest\texon\t201\t400\t.\t-\t.\tgene_id "g4"; transcript_id "t4";
chr1\ttest\texon\t801\t1000\t.\t-\t.\tgene_id "g4"; transcript_id "t4";
"""


@pytest.fixture(scope="module")
def index(tmp_path_factory):
    return build_test_index(tmp_path_factory, GTF, name="s1_junctions", refs=REFS)


@pytest.fixture(scope="module")
def coincident_index(tmp_path_factory):
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        return build_test_index(
            tmp_path_factory, COINCIDENT_GTF, name="s1_coincident", refs={"chr1": 2000}
        )


def _cut_index(cuts, offsets, ref_id, position):
    """The flat cut index of ``position`` on reference ``ref_id``, or -1 if it is not a cut."""
    lo, hi = int(offsets[ref_id]), int(offsets[ref_id + 1])
    k = lo + int(np.searchsorted(cuts[lo:hi], position))
    return k if k < hi and int(cuts[k]) == position else -1


def _lookup(arrays, boundary_left, boundary_right):
    """What the deposit's inner loop does: scan the donor's CSR slice for the acceptor.

    Returns the **junction-boundary id, which is the CSR slot itself**, paired with its strand — or ``None``.
    One to three iterations at human scale (measured mean fan-out over donor cuts: 1.31, max 25).

    ⚠ It returns ``k``, not ``edge_row[k]``. ``edge_row`` is the key for joining back to ``edges_df`` and
    is not an index into any junction bank; see :class:`JunctionEdgeArrays`.
    """
    if boundary_left < 0 or boundary_right < 0:
        return None
    for k in range(int(arrays.offsets[boundary_left]), int(arrays.offsets[boundary_left + 1])):
        if int(arrays.boundary_right[k]) == boundary_right:
            return k, int(arrays.strand[k])
    return None


def test_the_csr_addresses_the_flat_cut_axis(index):
    """One slot per cut, and the totals close against the boundary table."""
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    n_junction = int((index.edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_JUNCTION).sum())
    assert arrays.offsets.shape == (cuts.shape[0] + 1,)
    assert int(arrays.offsets[0]) == 0
    assert int(arrays.offsets[-1]) == n_junction
    assert arrays.boundary_right.shape == arrays.edge_row.shape == arrays.strand.shape
    assert arrays.boundary_right.shape == (n_junction,)
    assert int(offsets[-1]) == cuts.shape[0]


def test_every_annotated_intron_is_found_at_its_LEFT_BOUNDARY(index):
    """The four annotated introns, looked up the way the deposit will look them up."""
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    expected = [
        (0, 400, 700, Strand.POS),  # t0 intron 1
        (0, 900, 1200, Strand.POS),  # t0 intron 2
        (0, 400, 1000, Strand.POS),  # t1, sharing t0's donor
        (1, 500, 800, Strand.NEG),  # t2 on chr2
    ]
    for ref_id, start, end, strand in expected:
        hit = _lookup(
            arrays,
            _cut_index(cuts, offsets, ref_id, start),
            _cut_index(cuts, offsets, ref_id, end),
        )
        assert hit is not None, f"intron [{start},{end}) on ref {ref_id} not found"
        slot, got_strand = hit
        assert 0 <= slot < arrays.boundary_right.shape[0]
        row = index.edges_df.iloc[int(arrays.edge_row[slot])]  # the JOIN, not the id
        assert int(row["kind"]) == EDGE_KIND_JUNCTION
        assert got_strand == int(strand)


def test_a_shared_donor_keeps_BOTH_junctions(index):
    """⭐ Alternative 3' splice site: cut 400 on chr1 is the donor of two distinct junctions.

    A table storing one junction per cut would silently drop one of them, and the loss would be
    invisible — the dropped intron would simply be treated as unannotated.
    """
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    donor = _cut_index(cuts, offsets, 0, 400)
    lo, hi = int(arrays.offsets[donor]), int(arrays.offsets[donor + 1])
    assert hi - lo == 2
    landed = sorted(int(cuts[arrays.boundary_right[k]]) for k in range(lo, hi))
    assert landed == [700, 1000]


def test_a_cut_that_is_not_a_LEFT_BOUNDARY_has_an_empty_slice(index):
    """Most cuts are not left_boundaries — measured 70.4 % on both the toy and the human annotation. The slice
    must be empty rather than absent, so the deposit needs no special case."""
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    for ref_id, position in ((0, 200), (0, 1400), (1, 300)):  # TSS / TES, never a donor
        cut = _cut_index(cuts, offsets, ref_id, position)
        assert cut >= 0, f"{position} should be a cut on ref {ref_id}"
        assert int(arrays.offsets[cut + 1]) - int(arrays.offsets[cut]) == 0


def test_an_unannotated_intron_does_not_match(index):
    """A coordinate pair that is not an annotated junction must miss, even when both of its endpoints
    happen to be cuts — that is what routes it to the unspliced channel (design §4.1)."""
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    # 700 and 900 are both cuts on chr1, but [700,900) is an EXON, not an intron
    assert (
        _lookup(arrays, _cut_index(cuts, offsets, 0, 700), _cut_index(cuts, offsets, 0, 900))
        is None
    )
    # and a position that is not a cut at all
    assert _cut_index(cuts, offsets, 0, 401) == -1


def test_the_csr_round_trips_to_the_boundary_table(index):
    """Re-derive the junction set from the CSR and compare with ``edges_df`` — the two agree only if the
    region-id → cut-index shift is right on every reference, which is the one thing that can silently
    break when a reference has no regions."""
    arrays = build_junction_edge_arrays(index)
    cuts, offsets, _ = build_region_partition_arrays(index)
    donor = np.repeat(np.arange(cuts.shape[0]), np.diff(arrays.offsets))
    from_csr = np.stack(
        [cuts[donor], cuts[arrays.boundary_right], arrays.strand.astype(np.int64)], axis=1
    )

    boundaries = index.edges_df
    junction = boundaries["kind"].to_numpy(np.uint8) == EDGE_KIND_JUNCTION
    regions = index.regions_df
    src, dst = boundaries["src"].to_numpy(np.int64)[junction], boundaries["dst"].to_numpy(np.int64)[junction]
    from_boundaries = np.stack(
        [
            regions["end"].to_numpy(np.int64)[src],  # the intron starts where src ENDS
            regions["start"].to_numpy(np.int64)[dst],  # and ends where dst BEGINS
            boundaries["strand"].to_numpy(np.int8)[junction].astype(np.int64),
        ],
        axis=1,
    )
    order = lambda a: a[np.lexsort(a.T[::-1])]  # noqa: E731
    assert np.array_equal(order(from_csr), order(from_boundaries))


def _reference_partition(index):
    """The same graph, built through the reference accumulator's OWN constructor.

    Deliberately independent of :func:`build_junction_edge_arrays`: this route names each junction by its
    genomic ``(ref, intron_start, intron_end, strand)`` and lets ``Partition.from_cuts`` resolve both
    endpoints with its own ``searchsorted``, where the builder walks region ids and applies a per-reference
    ``cut_base − region_base`` shift. Only the *definition* of a junction is shared.
    """
    cuts, cut_offsets, region_types = build_region_partition_arrays(index)
    boundaries, regions = index.edges_df, index.regions_df
    junction = boundaries["kind"].to_numpy(np.uint8) == EDGE_KIND_JUNCTION
    src = boundaries["src"].to_numpy(np.int64)[junction]
    dst = boundaries["dst"].to_numpy(np.int64)[junction]
    strand = boundaries["strand"].to_numpy(np.int8)[junction]
    region_end, region_start = regions["end"].to_numpy(np.int64), regions["start"].to_numpy(np.int64)
    ref_of_region = regions["ref_name"].to_numpy()
    ref_id = {name: i for i, name in enumerate(index.ref_names)}

    junctions = [
        # the intron starts where src ENDS and ends where dst BEGINS
        (ref_id[ref_of_region[s]], int(region_end[s]), int(region_start[d]), int(st))
        for s, d, st in zip(src, dst, strand)
    ]
    n_refs = len(index.ref_names)
    return Partition.from_cuts(
        [cuts[cut_offsets[r] : cut_offsets[r + 1]] for r in range(n_refs)],
        # a reference contributing c cuts owns c-1 regions, so r earlier references own cut_offsets[r]-r
        region_types=[
            region_types[cut_offsets[r] - r : cut_offsets[r + 1] - r - 1] for r in range(n_refs)
        ],
        junctions=junctions,
    )


@pytest.mark.parametrize("fixture", ["index", "coincident_index"])
def test_the_csr_slot_order_matches_the_reference_accumulator(fixture, request):
    """⛔ THE CONTRACT: the junction-boundary id IS the CSR slot, so both builders must emit the slots in the
    SAME order — otherwise every junction row permutes and the native build's byte-identity gate compares
    two different labellings of one graph.

    ⚠ This had never been tested. The two orderings disagreed once during S2 — ``(acceptor, donor)``
    against ``(strand, acceptor, donor)`` — and nothing would have caught it: the spec matrix exercises
    only ``Partition.from_cuts``, and the real-data shim builds its ``Partition`` straight from
    ``build_junction_edge_arrays``, so the two sorts were never compared to each other.

    ⚠ **What this test does and does not cover**, measured by perturbing the builder's key:

    * promoting ``strand`` to the primary key, or swapping donor/acceptor priority → **caught**, but only
      because of chr3's nested pair (see ``GTF``).
    * *dropping* ``strand`` from the builder's key → **NOT caught, and cannot be.** ⭐ ``edges_df`` emits a
      strand-coincident pair POS-before-NEG whichever order the GTF lists the genes in, and ``np.lexsort``
      is stable, so both routes start from an already-correct tie order. The builder's ``strand`` key is
      therefore *defensive* — keep it, because it makes the contract explicit instead of resting on
      ``edges_df``'s internal sort, but do not expect this test to defend it.
    * ``from_cuts``'s ``strand`` key **is** load-bearing, since a caller may pass junctions in any order,
      and it is pinned by ``test_a_junction_id_is_a_function_of_the_PARTITION_not_of_argument_order`` in
      the spec matrix.
    """
    index = request.getfixturevalue(fixture)
    arrays = build_junction_edge_arrays(index)
    reference = _reference_partition(index)

    assert np.array_equal(reference.sj_offsets, arrays.offsets)
    assert np.array_equal(reference.sj_boundary_right, arrays.boundary_right)
    assert np.array_equal(reference.sj_strand, arrays.strand)


def test_a_strand_coincident_pair_is_two_distinct_slots(coincident_index):
    """Two genes on opposite strands sharing their intron coordinates exactly.

    Both must survive as separate junction boundaries, in the same CSR slice and ordered by strand — that
    adjacency is the only thing that makes the ordering contract above falsifiable, since every other
    junction is already separated by its donor or its acceptor.
    """
    arrays = build_junction_edge_arrays(coincident_index)
    cuts, offsets, _ = build_region_partition_arrays(coincident_index)
    donor = _cut_index(cuts, offsets, 0, 400)
    lo, hi = int(arrays.offsets[donor]), int(arrays.offsets[donor + 1])
    assert hi - lo == 2, "the strand-coincident pair collapsed to one junction boundary"
    assert [int(cuts[arrays.boundary_right[k]]) for k in range(lo, hi)] == [800, 800]
    assert [int(arrays.strand[k]) for k in range(lo, hi)] == [int(Strand.POS), int(Strand.NEG)]
