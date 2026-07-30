"""The v8 SPLICE GRAPH test matrix — G1-G18, I1-I12, P1-P5.

    Design: docs/CARRY_FORWARD.md §7   ·   Plan: docs/CARRY_FORWARD.md

The graph replaces the v7 region/boundary partition. Its ONE behavioural change is that adjacent
equal-signature segments are no longer merged, which is what preserves transcript termini — 59.5 % of
human TSS/TES fall strictly inside a merged v7 region today.

⚠ **The migration gates P2/P2′/P3 are GONE with the partition they compared against.** P3 — "merging
adjacent equal-signature nodes reproduces regions.feather" — was the only independent check on the
signature computation, and it is replaced by **I3b**: ``validate_graph`` recomputes every node's
signature from its MIDPOINT by direct interval containment, a different algorithm from the builder's
cumulative-difference sweep over the same interval sets. **I13** does the same for the flags. Both
are proven to fire.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from rigel.index import TranscriptIndex

from rigel.calibration.splice_graph import (
    EDGE_KIND_CONTIGUOUS,
    EDGE_KIND_JUNCTION,
    FLAG_ACCEPTOR_NEG,
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_NEG,
    FLAG_DONOR_POS,
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
    build_splice_graph,
    validate_graph,
)
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_POS
from rigel.transcript import Interval, Transcript
from rigel.types import Strand

REF = {"chr1": 2000}


def _tx(exons, strand=Strand.POS, ref="chr1", t_id=None, **kw):
    return Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(a, b) for a, b in exons],
        t_id=t_id or f"t{len(exons)}_{exons[0][0]}",
        **kw,
    )


def _graph(transcripts, reflen=None):
    reflen = reflen or REF
    nodes, edges = build_splice_graph(transcripts, reflen)
    with warnings.catch_warnings():
        # Several cases here deliberately build strand-coincident junctions (G18, the per-strand reach
        # case). validate_graph warns on those because they are biologically impossible; the warning
        # itself is asserted by test_strand_coincident_junctions_warn_but_still_work, so it is muted
        # here to keep the suite warning-clean rather than swallowed anywhere it is not already tested.
        warnings.simplefilter("ignore", RuntimeWarning)
        validate_graph(nodes, edges, reflen, transcripts=transcripts)  # I1-I12 on every case
    return nodes, edges


def _nodes(nodes):
    return list(zip(nodes["start"].tolist(), nodes["end"].tolist()))


def _junctions(edges):
    j = edges[edges["kind"] == EDGE_KIND_JUNCTION]
    return list(zip(j["src"].tolist(), j["dst"].tolist(), j["strand"].tolist()))


def _flags_at(nodes, edges, pos):
    """The flags on the contiguous edge at genomic position ``pos``."""
    c = edges[edges["kind"] == EDGE_KIND_CONTIGUOUS]
    ends = nodes["end"].to_numpy()
    hit = [int(f) for s, f in zip(c["src"].tolist(), c["flags"].tolist()) if ends[s] == pos]
    assert len(hit) == 1, f"expected one contiguous edge at {pos}, found {len(hit)}"
    return hit[0]


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# G1-G18 — the toy-GTF matrix. One graph assertion per structural case.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_G1_single_exon_transcript():
    """3 nodes (upstream / exon / downstream), 2 contiguous edges, 0 junctions, TSS+TES on the two."""
    n, e = _graph([_tx([(500, 800)])])
    assert _nodes(n) == [(0, 500), (500, 800), (800, 2000)]
    assert int((e["kind"] == EDGE_KIND_CONTIGUOUS).sum()) == 2
    assert _junctions(e) == []
    assert _flags_at(n, e, 500) & FLAG_TSS_POS  # + strand: 5' end is the low edge
    assert _flags_at(n, e, 800) & FLAG_TES_POS


def test_G2_two_exon_transcript():
    """One junction edge donor→acceptor; the intron node exists and is not on the transcript's path."""
    n, e = _graph([_tx([(500, 700), (1200, 1500)])])
    assert _nodes(n) == [(0, 500), (500, 700), (700, 1200), (1200, 1500), (1500, 2000)]
    assert _junctions(e) == [(1, 3, Strand.POS)]
    assert n["signature"][2] & BIT_INTRON_POS  # the skipped intron node is annotated as intron


def test_G3_three_exon_transcript():
    n, e = _graph([_tx([(200, 400), (700, 900), (1200, 1400)])])
    assert len(_junctions(e)) == 2


def test_G4_alternative_TSS_inside_another_exon():
    """⭐ THE case v7's merge deletes: an interior cut whose two sides share a signature."""
    n, e = _graph([_tx([(400, 1000)], t_id="long"), _tx([(600, 1000)], t_id="short")])
    assert (600, 1000) in _nodes(n), "the alternative-TSS cut at 600 was merged away"
    i = _nodes(n).index((400, 600))
    assert n["signature"][i] == n["signature"][i + 1], "both sides carry the same signature"
    assert _flags_at(n, e, 600) & FLAG_TSS_POS


def test_G5_alternative_TES_inside_another_exon():
    n, e = _graph([_tx([(400, 1000)], t_id="long"), _tx([(400, 700)], t_id="short")])
    assert (400, 700) in _nodes(n)
    assert _flags_at(n, e, 700) & FLAG_TES_POS


def test_G6_position_is_both_terminus_and_junction():
    """⭐ TES_s AND DONOR_s on one edge — the case the 4-bit signature is structurally blind to."""
    n, e = _graph([_tx([(400, 800)], t_id="ends"), _tx([(200, 800), (1200, 1500)], t_id="splices")])
    f = _flags_at(n, e, 800)
    assert f & FLAG_TES_POS and f & FLAG_DONOR_POS


def test_G7_exon_skipping():
    """Two junction edges out of A's last node; the undirected graph now has a cycle."""
    n, e = _graph(
        [
            _tx([(200, 400), (700, 900), (1200, 1400)], t_id="abc"),
            _tx([(200, 400), (1200, 1400)], t_id="ac"),
        ]
    )
    a_last = _nodes(n).index((200, 400))
    assert sum(1 for s, _d, _st in _junctions(e) if s == a_last) == 2


def test_G8_mutually_exclusive_exons():
    n, e = _graph(
        [
            _tx([(100, 300), (500, 700), (1500, 1700)], t_id="m1"),
            _tx([(100, 300), (900, 1100), (1500, 1700)], t_id="m2"),
        ]
    )
    assert len(_junctions(e)) == 4
    ns = _nodes(n)
    # no edge joins the two alternative exons directly
    assert (ns.index((500, 700)), ns.index((900, 1100)), Strand.POS) not in _junctions(e)


def test_G9_retained_intron():
    """The intron node carries exon_s AND intron_s; the junction edge spans it."""
    n, e = _graph(
        [
            _tx([(300, 500), (800, 1000)], t_id="spliced"),
            _tx([(300, 1000)], t_id="retained"),
        ]
    )
    i = _nodes(n).index((500, 800))
    assert n["signature"][i] & BIT_EXON_POS and n["signature"][i] & BIT_INTRON_POS
    assert len(_junctions(e)) == 1


def test_G10_overlapping_transcripts_opposite_strands():
    """AMBIG signature nodes; junction edges keep distinct strand; no strand leakage in the flags."""
    n, e = _graph(
        [
            _tx([(200, 500), (900, 1200)], strand=Strand.POS, t_id="p"),
            _tx([(300, 600), (1000, 1300)], strand=Strand.NEG, t_id="m"),
        ]
    )
    i = _nodes(n).index((300, 500))
    assert n["signature"][i] & BIT_EXON_POS and n["signature"][i] & BIT_EXON_NEG
    strands = {st for _s, _d, st in _junctions(e)}
    assert strands == {Strand.POS, Strand.NEG}
    # the + transcript's TSS at 200 must not set the − bit
    assert _flags_at(n, e, 200) & FLAG_TSS_POS
    assert not _flags_at(n, e, 200) & FLAG_TSS_NEG


def test_G11_nested_transcript_inside_another_intron():
    n, e = _graph(
        [
            _tx([(100, 300), (1500, 1700)], t_id="outer"),
            _tx([(700, 900)], t_id="inner"),
        ]
    )
    assert (700, 900) in _nodes(n)
    assert len(_junctions(e)) == 1


def test_G12_shared_exon_endpoint_across_transcripts():
    """Exactly one cut, one edge — the cut set dedups."""
    n, e = _graph([_tx([(400, 800)], t_id="a"), _tx([(400, 800)], t_id="b")])
    assert _nodes(n) == [(0, 400), (400, 800), (800, 2000)]


def test_G13_one_bp_node():
    """A node of length 1 is emitted and walkable — human has 15,687 of them."""
    n, _e = _graph([_tx([(500, 700)], t_id="a"), _tx([(501, 700)], t_id="b")])
    assert (500, 501) in _nodes(n)
    assert int(n["length"].min()) == 1


def test_G14_bookended_exons_no_intron():
    """Adjacent exons with a zero-length intron: a contiguous edge, NO junction edge."""
    n, e = _graph([_tx([(400, 700), (700, 1000)])])
    assert _junctions(e) == []
    assert 700 in n["end"].tolist()


def test_G15_transcript_at_reference_edges():
    """No zero-length node and no duplicate cut when an exon touches 0 or ref_length."""
    n, _e = _graph([_tx([(0, 300), (1700, 2000)])])
    assert int(n["length"].min()) > 0
    assert _nodes(n)[0][0] == 0 and _nodes(n)[-1][1] == 2000
    assert len(set(_nodes(n))) == len(_nodes(n))


def test_G16_reference_with_no_transcripts():
    n, e = _graph([], reflen={"chr1": 2000})
    assert _nodes(n) == [(0, 2000)]
    assert len(e) == 0


def test_G17_two_references():
    """Node/edge ids contiguous per reference; NO cross-reference edge."""
    reflen = {"chr1": 1500, "chr2": 1500}
    n, e = _graph(
        [_tx([(200, 400), (900, 1100)], ref="chr1"), _tx([(300, 500)], ref="chr2")], reflen=reflen
    )
    ref = n["ref_name"].astype(str).to_numpy()
    for name in reflen:
        rows = np.flatnonzero(ref == name)
        assert np.array_equal(rows, np.arange(rows[0], rows[0] + rows.size))
    assert np.all(ref[e["src"].to_numpy()] == ref[e["dst"].to_numpy()])


def test_G18_coincident_opposite_strand_junctions():
    """Two junction edges at the same donor/acceptor, distinct strand. Zero occurrences in GENCODE
    (verified on the real index), so this exists to prove it WORKS, not because it fires."""
    n, e = _graph(
        [
            _tx([(300, 500), (900, 1100)], strand=Strand.POS, t_id="p"),
            _tx([(300, 500), (900, 1100)], strand=Strand.NEG, t_id="m"),
        ]
    )
    js = _junctions(e)
    assert len(js) == 2
    assert {st for _s, _d, st in js} == {Strand.POS, Strand.NEG}
    assert len({(s, d) for s, d, _st in js}) == 1  # same endpoints


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# The REACH columns (plan C1/C2/C3) — not in the original matrix, added because the schema is new.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def _reach(edges, src, dst, kind, strand=None):
    row = edges[(edges["src"] == src) & (edges["dst"] == dst) & (edges["kind"] == kind)]
    if strand is not None:
        row = row[row["strand"] == strand]
    assert len(row) == 1
    r = row.iloc[0]
    return (
        int(r["reach_lo_pos"]),
        int(r["reach_hi_pos"]),
        int(r["reach_lo_neg"]),
        int(r["reach_hi_neg"]),
    )


def test_reach_on_a_junction_is_the_exonic_length_either_side():
    """The owner's worked example: TSS 500, first exon [500,550), junction at 550 ⇒ reach_lo = 50."""
    n, e = _graph([_tx([(500, 550), (1000, 1300)])])
    ns = _nodes(n)
    lo, hi, nlo, nhi = _reach(e, ns.index((500, 550)), ns.index((1000, 1300)), EDGE_KIND_JUNCTION)
    assert (lo, hi) == (50, 300)  # 50 exonic bases before the intron, 300 after
    assert (nlo, nhi) == (0, 0)  # nothing on the − strand


def test_reach_is_maximal_over_isoforms_independently_per_side():
    """D2: a position open on ANY isoform is open. The two isoforms disagree on BOTH sides here."""
    n, e = _graph(
        [
            _tx([(400, 600), (1000, 1100)], t_id="short"),  # lo 200, hi 100
            _tx([(300, 600), (1000, 1400)], t_id="long"),  # lo 300, hi 400
        ]
    )
    ns = _nodes(n)
    lo, hi, _, _ = _reach(e, ns.index((400, 600)), ns.index((1000, 1100)), EDGE_KIND_JUNCTION)
    assert (lo, hi) == (300, 400)


def test_reach_is_per_strand_and_does_not_mix():
    """⚠ plan C2: the two strands have different reaches at the SAME endpoints and must not be conflated.

    Both transcripts splice 200→600, so this is also a G18 coincident-junction pair: two edges sharing
    ``(src, dst)`` and differing only in strand. Each must carry its OWN reach in its OWN columns and
    zero in the other strand's — a strand-agnostic maximum would give the + junction the − transcript's
    1000 bp downstream reach and over-state its mature opportunity 10-fold.
    """
    n, e = _graph(
        [
            _tx([(0, 200), (600, 700)], strand=Strand.POS, t_id="p"),  # 200 before, 100 after
            _tx([(0, 200), (600, 1600)], strand=Strand.NEG, t_id="m"),  # 200 before, 1000 after
        ]
    )
    ns = _nodes(n)
    src, dst = ns.index((0, 200)), ns.index((600, 700))
    assert _reach(e, src, dst, EDGE_KIND_JUNCTION, Strand.POS) == (200, 100, 0, 0)
    assert _reach(e, src, dst, EDGE_KIND_JUNCTION, Strand.NEG) == (0, 0, 200, 1000)


def test_reach_on_a_contiguous_edge_inside_an_exon():
    """⭐ plan C1/Q1: reaches live on CONTIGUOUS edges too, which is where the taper near a TES bites."""
    n, e = _graph([_tx([(400, 1000)], t_id="a"), _tx([(700, 1000)], t_id="b")])
    ns = _nodes(n)
    i = ns.index((400, 700))
    lo, hi, _, _ = _reach(e, i, i + 1, EDGE_KIND_CONTIGUOUS)
    assert (lo, hi) == (300, 300)  # transcript "a": 300 exonic bases either side of position 700


def test_contiguous_reach_is_NONZERO_INSIDE_AN_INTRON():
    """⭐ The reach on a CONTIGUOUS edge is the genomic distance to the transcript's span ends.

    Nascent RNA is an ordinary transcript spanning its gene, so RNA opportunity inside an intron is
    real — the intron is where nascent RNA lives. An exonic reach would report 0 here and so declare
    zero RNA opportunity across every intron in the genome.

    Fixture: one + transcript, exons [500,700) and [1200,1500), so span [500,1500) and intron
    [700,1200). Both interior interfaces sit on that span.
    """
    n, e = _graph([_tx([(500, 700), (1200, 1500)])])
    ns = _nodes(n)
    at_donor = ns.index((500, 700))  # the edge at 700, the intron's low end
    at_acceptor = ns.index((700, 1200))  # the edge at 1200, the intron's high end
    assert _reach(e, at_donor, at_donor + 1, EDGE_KIND_CONTIGUOUS) == (200, 800, 0, 0)
    assert _reach(e, at_acceptor, at_acceptor + 1, EDGE_KIND_CONTIGUOUS) == (700, 300, 0, 0)


def test_reach_is_zero_outside_a_span_and_on_a_strand_with_no_transcript():
    """A reach of 0 is meaningful, not a sentinel — but it now means *no RNA of any form* here.

    Two ways to earn a zero, both asserted: the outward side of a span's own edge (nothing continues
    past a transcript end), and every position on a strand carrying no transcript at all.
    """
    n, e = _graph([_tx([(500, 700), (1200, 1500)], strand=Strand.POS)])
    ns = _nodes(n)
    at_tss = ns.index((0, 500))  # the edge at 500 — the span's low edge
    at_tes = ns.index((1200, 1500))  # the edge at 1500 — the span's high edge
    assert _reach(e, at_tss, at_tss + 1, EDGE_KIND_CONTIGUOUS) == (0, 1000, 0, 0)
    assert _reach(e, at_tes, at_tes + 1, EDGE_KIND_CONTIGUOUS) == (1000, 0, 0, 0)


def test_a_SYNTHETIC_nrna_span_is_excluded_from_everything():
    """A manufactured nRNA span row contributes no cut, no flag and no reach.

    ⚠ This test used to build its span with ``is_nrna=True`` alone and assert the *opposite* — that
    the span cuts the partition but sets no flags. That encoded a real bug (see the test below):
    on a NON-synthetic row ``is_nrna`` means "this real transcript is single-exon", not "this is a
    manufactured span". The rows that must be excluded are the SYNTHETIC ones, and ``~is_synthetic``
    excludes exactly them.
    """
    real = _tx([(400, 600), (1000, 1200)], t_id="real")
    span = _tx([(300, 1500)], t_id="span", is_nrna=True, is_synthetic=True)
    n, e = _graph([real, span])
    cuts = set(n["start"].tolist()) | set(n["end"].tolist())
    assert not ({300, 1500} & cuts), "a synthetic span must contribute NO cut"
    assert set(n["end"].tolist()[:-1]) == {400, 600, 1000, 1200}


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# I1-I12 — the validators must FIRE. A validator that cannot fail is not a validator.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def sample_graph():
    txs = [
        _tx([(200, 400), (700, 900), (1200, 1400)], t_id="a"),
        _tx([(250, 400), (1200, 1600)], strand=Strand.NEG, t_id="b"),
        _tx([(1700, 1900)], t_id="c"),
    ]
    return (*build_splice_graph(txs, REF), txs)


def test_I_all_hold_on_a_valid_graph(sample_graph):
    n, e, txs = sample_graph
    validate_graph(n, e, REF, transcripts=txs)


@pytest.mark.parametrize(
    "inv,mutate",
    [
        ("I1 tiling", lambda n, e: n.__setitem__("end", n["end"].mask(n.index == 0, 999))),
        ("I2 node_id", lambda n, e: n.__setitem__("node_id", n["node_id"] + 1)),
        (
            "I3 signature",
            lambda n, e: n.__setitem__("signature", n["signature"].mask(n.index == 0, 99)),
        ),
        ("I8 src<dst", lambda n, e: e.__setitem__("dst", e["dst"].mask(e.index == 0, 0))),
        ("I12 edge_id", lambda n, e: e.__setitem__("edge_id", e["edge_id"] + 5)),
    ],
)
def test_I_validators_fire_when_violated(sample_graph, inv, mutate):
    n, e, _txs = sample_graph
    n, e = n.copy(), e.copy()
    mutate(n, e)
    with pytest.raises(ValueError):
        validate_graph(n, e, REF)


def test_I11_fires_when_a_junction_edge_is_missing(sample_graph):
    n, e, txs = sample_graph
    e = e[e["kind"] != EDGE_KIND_JUNCTION].reset_index(drop=True)
    e["edge_id"] = np.arange(len(e), dtype=np.int64)
    with pytest.raises(ValueError, match="I11"):
        validate_graph(n, e, REF, transcripts=txs)


def test_I4_fires_when_an_interface_is_not_an_annotation_event(sample_graph):
    """The strongest structural statement: interior interfaces are EXACTLY the annotation events."""
    n, e, txs = sample_graph
    with pytest.raises(ValueError, match="I4"):
        validate_graph(n, e, REF, transcripts=txs + [_tx([(1000, 1050)], t_id="unbuilt")])


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# P1-P5 — properties. P2′ and P3 are THE migration gates.
# ═══════════════════════════════════════════════════════════════════════════════════════════════

_RANDOM_CASES = [
    [_tx([(100, 300)])],
    [_tx([(100, 300), (600, 900)]), _tx([(150, 300), (600, 800)])],
    [_tx([(0, 400)], strand=Strand.NEG), _tx([(200, 600), (1000, 1200)])],
    [_tx([(100, 200), (200, 300), (900, 1000)]), _tx([(100, 1000)], t_id="ri")],
    [_tx([(500, 501)]), _tx([(500, 502)])],
    [
        _tx([(100, 400), (800, 1200)], strand=s, t_id=f"x{i}")
        for i, s in enumerate((Strand.POS, Strand.NEG))
    ],
]


@pytest.mark.parametrize("txs", _RANDOM_CASES, ids=[f"case{i}" for i in range(len(_RANDOM_CASES))])
def test_P1_invariants_hold_on_every_case(txs):
    """I1-I13 on every random case — including I3b (the signature, recomputed from each node's
    midpoint) and I13 (the flags ARE the events), both of which need the transcripts."""
    n, e = build_splice_graph(txs, REF)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # case5 is strand-coincident by design
        validate_graph(n, e, REF, transcripts=txs)


def test_I3b_FIRES_on_a_corrupted_signature():
    """A validator that cannot fail is worthless — this is the proof it can."""
    txs = [_tx([(300, 500), (900, 1100)])]
    n, e = build_splice_graph(txs, REF)
    n.loc[0, "signature"] = np.uint8(BIT_EXON_POS)  # node 0 is intergenic
    with pytest.raises(ValueError, match="I3.*recomputed"):
        validate_graph(n, e, REF, transcripts=txs)


def test_a_single_exon_transcript_is_a_REAL_transcript_with_REAL_termini():
    """⛔ THE BUG plan F1's second filter caused, pinned.

    F1 specified the flags and reaches as ``~is_synthetic & ~is_nrna``, reasoning that "an nRNA span's
    ends are not real transcript termini". True — but on a NON-synthetic row ``is_nrna`` does not mean
    "manufactured span": it means **the transcript is single-exon, so mature ≡ nascent**. Measured on
    the human annotation, all **26,475** ``is_nrna & ~is_synthetic`` rows have ``n_exons == 1`` and
    none is a ``RIGEL_NRNA_*`` row, so the extra clause deleted **52,104 distinct real terminus
    positions** — the exact visibility the whole v8 partition exists to buy.
    """
    single = _tx([(300, 1500)], t_id="single_exon", is_nrna=True)  # is_synthetic is False
    txs = [_tx([(400, 600), (1000, 1200)], t_id="real"), single]
    n, e = build_splice_graph(txs, REF)
    validate_graph(n, e, REF, transcripts=txs)
    assert {300, 1500} <= set(n["start"].tolist()) | set(n["end"].tolist())
    assert _flags_at(n, e, 300) & FLAG_TSS_POS, "a real single-exon transcript's 5' end IS a TSS"
    assert _flags_at(n, e, 1500) & FLAG_TES_POS, "...and its 3' end IS a TES"


def test_reach_covers_a_single_exon_transcript():
    """The same row is also a real MATURE molecule, so it carries reach — 1,200 exonic bases split
    at whichever interior interface falls inside it."""
    single = _tx([(300, 1500)], t_id="single_exon", is_nrna=True)
    n, e = _graph([single, _tx([(700, 800)], t_id="cutter")])
    ns = _nodes(n)
    i = ns.index((300, 700))
    lo, hi, _ln, _hn = _reach(e, i, i + 1, EDGE_KIND_CONTIGUOUS)
    assert (lo, hi) == (400, 800), "400 exonic bases before position 700, 800 after"


def test_P4_determinism_build_twice_is_identical():
    txs = _RANDOM_CASES[3]
    a_n, a_e = build_splice_graph(txs, REF)
    b_n, b_e = build_splice_graph(list(reversed(txs)), REF)
    assert a_n.equals(b_n), "node table depends on transcript ORDER"
    assert a_e.equals(b_e), "edge table depends on transcript ORDER"


def test_P5_every_transcript_walks_on_a_realistic_multilocus_case():
    """I11 over a denser case: overlapping loci, both strands, a retained intron and a 1 bp node."""
    txs = [
        _tx([(100, 300), (600, 800), (1100, 1300)], t_id="a"),
        _tx([(100, 800), (1100, 1300)], t_id="a_ri"),
        _tx([(150, 300), (1100, 1200)], t_id="a_alt"),
        _tx([(400, 700), (900, 1000)], strand=Strand.NEG, t_id="b"),
        _tx([(1500, 1501)], t_id="tiny"),
        _tx([(1500, 1900)], strand=Strand.NEG, t_id="c"),
    ]
    n, e = build_splice_graph(txs, REF)
    validate_graph(n, e, REF, transcripts=txs)
    assert int(n["length"].min()) == 1


def test_edge_rows_are_sorted_by_src_kind_dst():
    """I12 as a contract, not just a validator: out-edges of a node are contiguous ⇒ CSR is one
    searchsorted, which every downstream consumer depends on."""
    n, e = build_splice_graph(_RANDOM_CASES[3] + _RANDOM_CASES[1], REF)
    key = list(zip(e["src"].tolist(), e["kind"].tolist(), e["dst"].tolist()))
    assert key == sorted(key)
    assert e["edge_id"].tolist() == list(range(len(e)))


def test_flags_are_not_mutually_exclusive():
    """All four bit classes can co-occur at one position on one strand — the case the signature
    cannot represent, and the reason the flags exist at all."""
    txs = [
        _tx([(200, 600)], t_id="ends_at_600"),
        _tx([(600, 900)], t_id="starts_at_600"),
        _tx([(300, 600), (900, 1100)], t_id="donates_at_600"),
        _tx([(100, 200), (600, 800)], t_id="accepts_at_600"),
    ]
    n, e = _graph(txs)
    f = _flags_at(n, e, 600)
    for bit, name in (
        (FLAG_TES_POS, "TES"),
        (FLAG_TSS_POS, "TSS"),
        (FLAG_DONOR_POS, "DONOR"),
        (FLAG_ACCEPTOR_POS, "ACCEPTOR"),
    ):
        assert f & bit, f"{name}_POS not set at position 600"
    assert not f & (FLAG_TSS_NEG | FLAG_TES_NEG | FLAG_DONOR_NEG | FLAG_ACCEPTOR_NEG)


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# T-D1 / integration — byte-identical rebuilds, and the artifact actually lands on disk.
# ═══════════════════════════════════════════════════════════════════════════════════════════════

_INTEGRATION_GTF = """\
chr1\ttest\texon\t201\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t701\t900\t.\t+\t.\tgene_id "g1"; transcript_id "t0";
chr1\ttest\texon\t251\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t701\t1000\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\ttest\texon\t1301\t1500\t.\t-\t.\tgene_id "g2"; transcript_id "t2";
chr1\ttest\texon\t1701\t1900\t.\t-\t.\tgene_id "g2"; transcript_id "t2";
"""


def test_TD1_rebuilds_are_byte_identical(tmp_path_factory):
    """Determinism is an INVARIANT, not an observation: the output is a pure function of
    (transcripts, ref_lengths) — np.unique sorts, ids come from position, edges from an explicit
    total order. No dict iteration, no hashing, no parallel reduction."""
    from conftest import build_test_index

    a = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_det_a")
    b = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_det_b")
    for fname in ("nodes.feather", "edges.feather"):
        pa = Path(a.index_dir) / fname
        pb = Path(b.index_dir) / fname
        assert pa.exists(), f"{fname} was not written by the index build"
        assert pa.read_bytes() == pb.read_bytes(), f"{fname} is not byte-identical across rebuilds"


def test_index_build_writes_and_loads_the_graph(tmp_path_factory):
    """R1: the graph lands on disk, reloads, and re-validates."""
    from conftest import build_test_index

    idx = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_int")
    assert idx.nodes_df is not None and idx.edges_df is not None
    validate_graph(idx.nodes_df, idx.edges_df, idx.ref_lengths)
    # I6 — ONE edge per DISTINCT (donor, acceptor, strand). t0 and t1 share intron [400,700) despite
    # different exon extents, so three intron INSTANCES dedup to two junction EDGES.
    assert int((idx.edges_df["kind"] == EDGE_KIND_JUNCTION).sum()) == 2


def test_graph_is_REQUIRED_at_load(tmp_path_factory):
    """⚠ INVERTED at W1b, deliberately. Through W1a the graph was OPTIONAL — nothing read it, so an
    index without one stayed fully usable and the arm was reversible. W1b made it THE partition the
    scanner deposits into, so an index without it cannot serve a scan, and loading one anyway is the
    worst available failure: calibration running on the merged v7 geometry while the caller believes
    it is on v8. It must raise, and say what to do.
    """
    from conftest import build_test_index

    idx = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_opt")
    for fname in ("nodes.feather", "edges.feather"):
        (Path(idx.index_dir) / fname).unlink()
        with pytest.raises(RuntimeError, match=r"(?s)splice graph.*[Rr]ebuild"):
            TranscriptIndex.load(idx.index_dir)


def test_strand_coincident_junctions_warn_but_still_work():
    """⚠ Biologically impossible: splice motifs are non-palindromic, so the same (donor, acceptor)
    cannot be a valid intron on both strands (a GT..AG intron reverse-complements to CT..AC). Measured:
    ZERO in GENCODE. One in a GTF means the ANNOTATION is wrong.

    The graph must (a) still be CORRECT — two distinct edges, kept apart by the strand in the sort key —
    and (b) SAY SO, because silence would let a bad annotation propagate into every downstream density.
    Warning rather than raising, so that this very case stays testable.
    """
    txs = [
        _tx([(300, 500), (900, 1100)], strand=Strand.POS, t_id="p"),
        _tx([(300, 500), (900, 1100)], strand=Strand.NEG, t_id="m"),
    ]
    nodes, edges = build_splice_graph(txs, REF)
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        validate_graph(nodes, edges, REF, transcripts=txs)
    js = _junctions(edges)
    assert len(js) == 2 and {st for _s, _d, st in js} == {Strand.POS, Strand.NEG}


def test_no_warning_on_a_biologically_normal_annotation():
    txs = [_tx([(300, 500), (900, 1100)], strand=Strand.POS), _tx([(1300, 1400), (1600, 1700)])]
    nodes, edges = build_splice_graph(txs, REF)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning fails the test
        validate_graph(nodes, edges, REF, transcripts=txs)
