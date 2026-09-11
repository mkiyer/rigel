"""The splice graph: the build matrix, the validators, the reach, and the sj CSR the deposit reads.

Adjacent equal-signature segments are not merged, which is what preserves transcript termini: under
a merging rule most human TSS/TES fall strictly inside a region and become invisible. The build
matrix walks the annotation shapes one at a time; the invariant block asserts that each validator
fires when its property is violated, the signature and flag validators recomputing from a region's
midpoint by interval containment — a different algorithm from the builder's sweep, so neither can
confirm the other's mistake. The last block gates the sj boundaries re-indexed onto the
accumulator's flat region-bound axis as a CSR, whose one question inside the scan's hot loop is
whether an observed intron is an annotated sj and which boundary it is: an intron that fails to
match deposits nothing, so a mis-keyed table silently deletes the whole spliced-RNA signal.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from rigel.index import TranscriptIndex

from rigel.calibration.splice_graph import (
    EDGE_KIND_CONTIGUOUS,
    EDGE_KIND_SJ,
    FLAG_ACCEPTOR_NEG,
    FLAG_ACCEPTOR_POS,
    FLAG_DONOR_NEG,
    FLAG_DONOR_POS,
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
    build_region_partition_arrays,
    build_sj_arrays,
    build_splice_graph,
    validate_graph,
)
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_POS
from rigel.transcript import Interval, Transcript
from rigel.types import Strand

from _index_builder import build_test_index
from native._accumulator_reference import Partition


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
    regions, boundaries = build_splice_graph(transcripts, reflen)
    with warnings.catch_warnings():
        # Several cases here deliberately build strand-coincident sj (G18, the per-strand reach
        # case). validate_graph warns on those because they are biologically impossible; the warning
        # itself is asserted by test_strand_coincident_sj_warn_but_still_work, so it is muted
        # here to keep the suite warning-clean rather than swallowed anywhere it is not already tested.
        warnings.simplefilter("ignore", RuntimeWarning)
        validate_graph(regions, boundaries, reflen, transcripts=transcripts)  # I1-I12 on every case
    return regions, boundaries


def _regions(regions):
    return list(zip(regions["start"].tolist(), regions["end"].tolist()))


def _sj(boundaries):
    j = boundaries[boundaries["kind"] == EDGE_KIND_SJ]
    return list(zip(j["src"].tolist(), j["dst"].tolist(), j["strand"].tolist()))


def _flags_at(regions, boundaries, pos):
    """The flags on the contiguous boundary at genomic position ``pos``."""
    c = boundaries[boundaries["kind"] == EDGE_KIND_CONTIGUOUS]
    ends = regions["end"].to_numpy()
    hit = [int(f) for s, f in zip(c["src"].tolist(), c["flags"].tolist()) if ends[s] == pos]
    assert len(hit) == 1, f"expected one contiguous boundary at {pos}, found {len(hit)}"
    return hit[0]


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# G1-G18 — the toy-GTF matrix. One graph assertion per structural case.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def test_G1_single_exon_transcript():
    """3 regions (upstream / exon / downstream), 2 contiguous boundaries, 0 sj, TSS+TES on the two."""
    n, e = _graph([_tx([(500, 800)])])
    assert _regions(n) == [(0, 500), (500, 800), (800, 2000)]
    assert int((e["kind"] == EDGE_KIND_CONTIGUOUS).sum()) == 2
    assert _sj(e) == []
    assert _flags_at(n, e, 500) & FLAG_TSS_POS  # + strand: 5' end is the low boundary
    assert _flags_at(n, e, 800) & FLAG_TES_POS


def test_G2_two_exon_transcript():
    """One sj boundary donor→acceptor; the intron region exists and is not on the transcript's path."""
    n, e = _graph([_tx([(500, 700), (1200, 1500)])])
    assert _regions(n) == [(0, 500), (500, 700), (700, 1200), (1200, 1500), (1500, 2000)]
    assert _sj(e) == [(1, 3, Strand.POS)]
    assert n["signature"][2] & BIT_INTRON_POS  # the skipped intron region is annotated as intron


def test_G3_three_exon_transcript():
    n, e = _graph([_tx([(200, 400), (700, 900), (1200, 1400)])])
    assert len(_sj(e)) == 2


def test_G4_alternative_TSS_inside_another_exon():
    """The case a merging partition deletes: an interior bound whose two sides share a signature."""
    n, e = _graph([_tx([(400, 1000)], t_id="long"), _tx([(600, 1000)], t_id="short")])
    assert (600, 1000) in _regions(n), "the alternative-TSS region_bound at 600 was merged away"
    i = _regions(n).index((400, 600))
    assert n["signature"][i] == n["signature"][i + 1], "both sides carry the same signature"
    assert _flags_at(n, e, 600) & FLAG_TSS_POS


def test_G5_alternative_TES_inside_another_exon():
    n, e = _graph([_tx([(400, 1000)], t_id="long"), _tx([(400, 700)], t_id="short")])
    assert (400, 700) in _regions(n)
    assert _flags_at(n, e, 700) & FLAG_TES_POS


def test_G6_position_is_both_terminus_and_sj():
    """TES_s AND DONOR_s on one boundary — the case the 4-bit signature is structurally blind to."""
    n, e = _graph([_tx([(400, 800)], t_id="ends"), _tx([(200, 800), (1200, 1500)], t_id="splices")])
    f = _flags_at(n, e, 800)
    assert f & FLAG_TES_POS and f & FLAG_DONOR_POS


def test_G7_exon_skipping():
    """Two sj boundaries out of A's last region; the undirected graph now has a cycle."""
    n, e = _graph(
        [
            _tx([(200, 400), (700, 900), (1200, 1400)], t_id="abc"),
            _tx([(200, 400), (1200, 1400)], t_id="ac"),
        ]
    )
    a_last = _regions(n).index((200, 400))
    assert sum(1 for s, _d, _st in _sj(e) if s == a_last) == 2


def test_G8_mutually_exclusive_exons():
    n, e = _graph(
        [
            _tx([(100, 300), (500, 700), (1500, 1700)], t_id="m1"),
            _tx([(100, 300), (900, 1100), (1500, 1700)], t_id="m2"),
        ]
    )
    assert len(_sj(e)) == 4
    ns = _regions(n)
    # no boundary joins the two alternative exons directly
    assert (ns.index((500, 700)), ns.index((900, 1100)), Strand.POS) not in _sj(e)


def test_G9_retained_intron():
    """The intron region carries exon_s AND intron_s; the sj boundary spans it."""
    n, e = _graph(
        [
            _tx([(300, 500), (800, 1000)], t_id="spliced"),
            _tx([(300, 1000)], t_id="retained"),
        ]
    )
    i = _regions(n).index((500, 800))
    assert n["signature"][i] & BIT_EXON_POS and n["signature"][i] & BIT_INTRON_POS
    assert len(_sj(e)) == 1


def test_G10_overlapping_transcripts_opposite_strands():
    """AMBIG signature regions; sj boundaries keep distinct strand; no strand leakage in the flags."""
    n, e = _graph(
        [
            _tx([(200, 500), (900, 1200)], strand=Strand.POS, t_id="p"),
            _tx([(300, 600), (1000, 1300)], strand=Strand.NEG, t_id="m"),
        ]
    )
    i = _regions(n).index((300, 500))
    assert n["signature"][i] & BIT_EXON_POS and n["signature"][i] & BIT_EXON_NEG
    strands = {st for _s, _d, st in _sj(e)}
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
    assert (700, 900) in _regions(n)
    assert len(_sj(e)) == 1


def test_G12_shared_exon_endpoint_across_transcripts():
    """Exactly one region_bound, one boundary — the region_bound set dedups."""
    n, e = _graph([_tx([(400, 800)], t_id="a"), _tx([(400, 800)], t_id="b")])
    assert _regions(n) == [(0, 400), (400, 800), (800, 2000)]


def test_G13_one_bp_region():
    """A region of length 1 is emitted and walkable — human has many thousands of them."""
    n, _e = _graph([_tx([(500, 700)], t_id="a"), _tx([(501, 700)], t_id="b")])
    assert (500, 501) in _regions(n)
    assert int(n["length"].min()) == 1


def test_G14_bookended_exons_no_intron():
    """Adjacent exons with a zero-length intron: a contiguous boundary, NO sj boundary."""
    n, e = _graph([_tx([(400, 700), (700, 1000)])])
    assert _sj(e) == []
    assert 700 in n["end"].tolist()


def test_G15_transcript_at_reference_boundaries():
    """No zero-length region and no duplicate region_bound when an exon touches 0 or ref_length."""
    n, _e = _graph([_tx([(0, 300), (1700, 2000)])])
    assert int(n["length"].min()) > 0
    assert _regions(n)[0][0] == 0 and _regions(n)[-1][1] == 2000
    assert len(set(_regions(n))) == len(_regions(n))


def test_G16_reference_with_no_transcripts():
    n, e = _graph([], reflen={"chr1": 2000})
    assert _regions(n) == [(0, 2000)]
    assert len(e) == 0


def test_G17_two_references():
    """Region/boundary ids contiguous per reference; NO cross-reference boundary."""
    reflen = {"chr1": 1500, "chr2": 1500}
    n, e = _graph(
        [_tx([(200, 400), (900, 1100)], ref="chr1"), _tx([(300, 500)], ref="chr2")], reflen=reflen
    )
    ref = n["ref_name"].astype(str).to_numpy()
    for name in reflen:
        rows = np.flatnonzero(ref == name)
        assert np.array_equal(rows, np.arange(rows[0], rows[0] + rows.size))
    assert np.all(ref[e["src"].to_numpy()] == ref[e["dst"].to_numpy()])


def test_G18_coincident_opposite_strand_sj():
    """Two sj boundaries at the same donor/acceptor, distinct strand. There are no occurrences in
    GENCODE, so this exists to prove the case works and not because it fires."""
    n, e = _graph(
        [
            _tx([(300, 500), (900, 1100)], strand=Strand.POS, t_id="p"),
            _tx([(300, 500), (900, 1100)], strand=Strand.NEG, t_id="m"),
        ]
    )
    js = _sj(e)
    assert len(js) == 2
    assert {st for _s, _d, st in js} == {Strand.POS, Strand.NEG}
    assert len({(s, d) for s, d, _st in js}) == 1  # same endpoints


# ═══════════════════════════════════════════════════════════════════════════════════════════════
# The REACH columns — the per-strand exonic and genomic distances each boundary carries.
# ═══════════════════════════════════════════════════════════════════════════════════════════════


def _reach(boundaries, src, dst, kind, strand=None):
    row = boundaries[
        (boundaries["src"] == src) & (boundaries["dst"] == dst) & (boundaries["kind"] == kind)
    ]
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


def test_reach_on_a_sj_is_the_exonic_length_either_side():
    """The owner's worked example: TSS 500, first exon [500,550), sj at 550 ⇒ reach_lo = 50."""
    n, e = _graph([_tx([(500, 550), (1000, 1300)])])
    ns = _regions(n)
    lo, hi, nlo, nhi = _reach(e, ns.index((500, 550)), ns.index((1000, 1300)), EDGE_KIND_SJ)
    assert (lo, hi) == (50, 300)  # 50 exonic bases before the intron, 300 after
    assert (nlo, nhi) == (0, 0)  # nothing on the − strand


def test_reach_is_maximal_over_isoforms_independently_per_side():
    """A position open on ANY isoform is open. The two isoforms disagree on BOTH sides here."""
    n, e = _graph(
        [
            _tx([(400, 600), (1000, 1100)], t_id="short"),  # lo 200, hi 100
            _tx([(300, 600), (1000, 1400)], t_id="long"),  # lo 300, hi 400
        ]
    )
    ns = _regions(n)
    lo, hi, _, _ = _reach(e, ns.index((400, 600)), ns.index((1000, 1100)), EDGE_KIND_SJ)
    assert (lo, hi) == (300, 400)


def test_reach_is_per_strand_and_does_not_mix():
    """The two strands have different reaches at the same endpoints and must not be conflated.

    Both transcripts splice 200→600, so this is also a G18 coincident-sj pair: two boundaries sharing
    ``(src, dst)`` and differing only in strand. Each must carry its OWN reach in its OWN columns and
    zero in the other strand's — a strand-agnostic maximum would give the + sj the − transcript's
    1000 bp downstream reach and over-state its mature opportunity 10-fold.
    """
    n, e = _graph(
        [
            _tx([(0, 200), (600, 700)], strand=Strand.POS, t_id="p"),  # 200 before, 100 after
            _tx([(0, 200), (600, 1600)], strand=Strand.NEG, t_id="m"),  # 200 before, 1000 after
        ]
    )
    ns = _regions(n)
    src, dst = ns.index((0, 200)), ns.index((600, 700))
    assert _reach(e, src, dst, EDGE_KIND_SJ, Strand.POS) == (200, 100, 0, 0)
    assert _reach(e, src, dst, EDGE_KIND_SJ, Strand.NEG) == (0, 0, 200, 1000)


def test_reach_on_a_contiguous_boundary_inside_an_exon():
    """Reaches live on contiguous boundaries too, which is where the taper near a TES bites."""
    n, e = _graph([_tx([(400, 1000)], t_id="a"), _tx([(700, 1000)], t_id="b")])
    ns = _regions(n)
    i = ns.index((400, 700))
    lo, hi, _, _ = _reach(e, i, i + 1, EDGE_KIND_CONTIGUOUS)
    assert (lo, hi) == (300, 300)  # transcript "a": 300 exonic bases either side of position 700


def test_contiguous_reach_is_NONZERO_INSIDE_AN_INTRON():
    """The reach on a CONTIGUOUS boundary is the genomic distance to the transcript's span ends.

    Nascent RNA is an ordinary transcript spanning its gene, so RNA opportunity inside an intron is
    real — the intron is where nascent RNA lives. An exonic reach would report 0 here and so declare
    zero RNA opportunity across every intron in the genome.

    Fixture: one + transcript, exons [500,700) and [1200,1500), so span [500,1500) and intron
    [700,1200). Both interior interfaces sit on that span.
    """
    n, e = _graph([_tx([(500, 700), (1200, 1500)])])
    ns = _regions(n)
    at_left_boundary = ns.index((500, 700))  # the boundary at 700, the intron's low end
    at_right_boundary = ns.index((700, 1200))  # the boundary at 1200, the intron's high end
    assert _reach(e, at_left_boundary, at_left_boundary + 1, EDGE_KIND_CONTIGUOUS) == (
        200,
        800,
        0,
        0,
    )
    assert _reach(e, at_right_boundary, at_right_boundary + 1, EDGE_KIND_CONTIGUOUS) == (
        700,
        300,
        0,
        0,
    )


def test_reach_is_zero_outside_a_span_and_on_a_strand_with_no_transcript():
    """A reach of 0 is meaningful, not a sentinel — but it now means *no RNA of any form* here.

    Two ways to earn a zero, both asserted: the outward side of a span's own boundary (nothing continues
    past a transcript end), and every position on a strand carrying no transcript at all.
    """
    n, e = _graph([_tx([(500, 700), (1200, 1500)], strand=Strand.POS)])
    ns = _regions(n)
    at_tss = ns.index((0, 500))  # the boundary at 500 — the span's low boundary
    at_tes = ns.index((1200, 1500))  # the boundary at 1500 — the span's high boundary
    assert _reach(e, at_tss, at_tss + 1, EDGE_KIND_CONTIGUOUS) == (0, 1000, 0, 0)
    assert _reach(e, at_tes, at_tes + 1, EDGE_KIND_CONTIGUOUS) == (1000, 0, 0, 0)


def test_a_SYNTHETIC_nrna_span_is_excluded_from_everything():
    """A manufactured nRNA span row contributes no region_bound, no flag and no reach.

    The filter is ``~is_synthetic`` alone. On a non-synthetic row ``is_nrna`` means "this real
    transcript is single-exon", not "this is a manufactured span", so using it as part of the
    exclusion deletes real transcripts (TRAPS: nrna-does-not-mean-synthetic); the test below is what
    that costs.
    """
    real = _tx([(400, 600), (1000, 1200)], t_id="real")
    span = _tx([(300, 1500)], t_id="span", is_nrna=True, is_synthetic=True)
    n, e = _graph([real, span])
    region_bounds = set(n["start"].tolist()) | set(n["end"].tolist())
    assert not ({300, 1500} & region_bounds), "a synthetic span must contribute NO region_bound"
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
        ("I2 region_id", lambda n, e: n.__setitem__("region_id", n["region_id"] + 1)),
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


def test_I11_fires_when_a_sj_edge_is_missing(sample_graph):
    n, e, txs = sample_graph
    e = e[e["kind"] != EDGE_KIND_SJ].reset_index(drop=True)
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
    """I1-I13 on every random case — including I3b (the signature, recomputed from each region's
    midpoint) and I13 (the flags ARE the events), both of which need the transcripts."""
    n, e = build_splice_graph(txs, REF)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # case5 is strand-coincident by design
        validate_graph(n, e, REF, transcripts=txs)


def test_I3b_FIRES_on_a_corrupted_signature():
    """A validator that cannot fail is worthless — this is the proof it can."""
    txs = [_tx([(300, 500), (900, 1100)])]
    n, e = build_splice_graph(txs, REF)
    n.loc[0, "signature"] = np.uint8(BIT_EXON_POS)  # region 0 is intergenic
    with pytest.raises(ValueError, match="I3.*recomputed"):
        validate_graph(n, e, REF, transcripts=txs)


def test_a_single_exon_transcript_is_a_REAL_transcript_with_REAL_termini():
    """What the second filter costs (TRAPS: nrna-does-not-mean-synthetic).

    Specifying the flags and reaches as ``~is_synthetic & ~is_nrna`` reasons that an nRNA span's ends
    are not real transcript termini. True — but on a non-synthetic row ``is_nrna`` does not mean
    "manufactured span": it means the transcript is single-exon, so mature is nascent. On the human
    annotation every such row is single-exon and none is a manufactured one, so the extra clause
    deletes tens of thousands of real terminus positions — the exact visibility the unmerged
    partition exists to buy.
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
    ns = _regions(n)
    i = ns.index((300, 700))
    lo, hi, _ln, _hn = _reach(e, i, i + 1, EDGE_KIND_CONTIGUOUS)
    assert (lo, hi) == (400, 800), "400 exonic bases before position 700, 800 after"


def test_P4_determinism_build_twice_is_identical():
    txs = _RANDOM_CASES[3]
    a_n, a_e = build_splice_graph(txs, REF)
    b_n, b_e = build_splice_graph(list(reversed(txs)), REF)
    assert a_n.equals(b_n), "region table depends on transcript ORDER"
    assert a_e.equals(b_e), "boundary table depends on transcript ORDER"


def test_P5_every_transcript_walks_on_a_realistic_multilocus_case():
    """I11 over a denser case: overlapping loci, both strands, a retained intron and a 1 bp region."""
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
    """I12 as a contract, not just a validator: out-boundaries of a region are contiguous ⇒ CSR is one
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
    (transcripts, ref_lengths) — np.unique sorts, ids come from position, boundaries from an explicit
    total order. No dict iteration, no hashing, no parallel reduction."""
    from _index_builder import build_test_index

    a = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_det_a")
    b = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_det_b")
    for fname in ("regions.feather", "edges.feather"):
        pa = Path(a.index_dir) / fname
        pb = Path(b.index_dir) / fname
        assert pa.exists(), f"{fname} was not written by the index build"
        assert pa.read_bytes() == pb.read_bytes(), f"{fname} is not byte-identical across rebuilds"


def test_index_build_writes_and_loads_the_graph(tmp_path_factory):
    """R1: the graph lands on disk, reloads, and re-validates."""
    from _index_builder import build_test_index

    idx = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_int")
    assert idx.regions_df is not None and idx.edges_df is not None
    validate_graph(idx.regions_df, idx.edges_df, idx.ref_lengths)
    # I6 — ONE boundary per DISTINCT (donor, acceptor, strand). t0 and t1 share intron [400,700) despite
    # different exon extents, so three intron INSTANCES dedup to two sj BOUNDARIES.
    assert int((idx.edges_df["kind"] == EDGE_KIND_SJ).sum()) == 2


def test_graph_is_REQUIRED_at_load(tmp_path_factory):
    """The graph is the partition the scanner deposits into, so an index without one cannot serve a
    scan. Loading one anyway is the worst available failure: calibration running on a merged
    geometry while the caller believes it is on the unmerged one. It must raise, and say what to do.
    """
    from _index_builder import build_test_index

    idx = build_test_index(tmp_path_factory, _INTEGRATION_GTF, genome_size=2000, name="sg_opt")
    for fname in ("regions.feather", "edges.feather"):
        (Path(idx.index_dir) / fname).unlink()
        with pytest.raises(RuntimeError, match=r"(?s)splice graph.*[Rr]ebuild"):
            TranscriptIndex.load(idx.index_dir)


def test_strand_coincident_sj_warn_but_still_work():
    """Biologically impossible: splice motifs are non-palindromic, so the same (donor, acceptor)
    cannot be a valid intron on both strands (a GT..AG intron reverse-complements to CT..AC), and
    there are none in GENCODE. One in a GTF means the annotation is wrong.

    The graph must (a) still be CORRECT — two distinct boundaries, kept apart by the strand in the sort key —
    and (b) SAY SO, because silence would let a bad annotation propagate into every downstream density.
    Warning rather than raising, so that this very case stays testable.
    """
    txs = [
        _tx([(300, 500), (900, 1100)], strand=Strand.POS, t_id="p"),
        _tx([(300, 500), (900, 1100)], strand=Strand.NEG, t_id="m"),
    ]
    regions, boundaries = build_splice_graph(txs, REF)
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        validate_graph(regions, boundaries, REF, transcripts=txs)
    js = _sj(boundaries)
    assert len(js) == 2 and {st for _s, _d, st in js} == {Strand.POS, Strand.NEG}


def test_no_warning_on_a_biologically_normal_annotation():
    txs = [_tx([(300, 500), (900, 1100)], strand=Strand.POS), _tx([(1300, 1400), (1600, 1700)])]
    regions, boundaries = build_splice_graph(txs, REF)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any warning fails the test
        validate_graph(regions, boundaries, REF, transcripts=txs)


# ── the sj boundaries as a CSR on the accumulator's flat region-bound axis ────────────────────


#: chr1 t0 splices 400->700 and 900->1200; t1 shares t0's DONOR at 400 but lands at 1000, so region_bound 400
#: has fan-out 2 — the alternative-3'-splice-site case a naive one-sj-per-region_bound table would drop.
#: chr2 carries a NEG-strand sj, which keeps the per-reference offsets and the strand honest.
#: chr3 is a NESTED pair with the OUTER intron on the minus strand: intron [400,1400) NEG encloses
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

#: A STRAND-COINCIDENT PAIR: two genes on opposite strands whose intron coordinates are byte-identical.
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
    return build_test_index(tmp_path_factory, GTF, name="s1_sj", refs=REFS)


@pytest.fixture(scope="module")
def coincident_index(tmp_path_factory):
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        return build_test_index(
            tmp_path_factory, COINCIDENT_GTF, name="s1_coincident", refs={"chr1": 2000}
        )


def _region_bound_index(region_bounds, offsets, ref_id, position):
    """The flat region_bound index of ``position`` on reference ``ref_id``, or -1 if it is not a region_bound."""
    lo, hi = int(offsets[ref_id]), int(offsets[ref_id + 1])
    k = lo + int(np.searchsorted(region_bounds[lo:hi], position))
    return k if k < hi and int(region_bounds[k]) == position else -1


def _lookup(arrays, boundary_left, boundary_right):
    """What the deposit's inner loop does: scan the donor's CSR slice for the acceptor.

    Returns the sj-boundary id, which is the CSR slot itself, paired with its strand — or ``None``.
    One to three iterations at human scale: a donor bound's fan-out is close to one.

    It returns ``k``, not ``edge_row[k]``. ``edge_row`` is the key for joining back to ``edges_df`` and
    is not an index into any sj bank; see :class:`SpliceJunctionArrays`.
    """
    if boundary_left < 0 or boundary_right < 0:
        return None
    for k in range(int(arrays.offsets[boundary_left]), int(arrays.offsets[boundary_left + 1])):
        if int(arrays.boundary_right[k]) == boundary_right:
            return k, int(arrays.strand[k])
    return None


def test_the_csr_addresses_the_flat_region_bound_axis(index):
    """One slot per region_bound, and the totals close against the boundary table."""
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    n_sj = int((index.edges_df["kind"].to_numpy(np.uint8) == EDGE_KIND_SJ).sum())
    assert arrays.offsets.shape == (region_bounds.shape[0] + 1,)
    assert int(arrays.offsets[0]) == 0
    assert int(arrays.offsets[-1]) == n_sj
    assert arrays.boundary_right.shape == arrays.edge_row.shape == arrays.strand.shape
    assert arrays.boundary_right.shape == (n_sj,)
    assert int(offsets[-1]) == region_bounds.shape[0]


def test_every_annotated_intron_is_found_at_its_LEFT_BOUNDARY(index):
    """The four annotated introns, looked up the way the deposit will look them up."""
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    expected = [
        (0, 400, 700, Strand.POS),  # t0 intron 1
        (0, 900, 1200, Strand.POS),  # t0 intron 2
        (0, 400, 1000, Strand.POS),  # t1, sharing t0's donor
        (1, 500, 800, Strand.NEG),  # t2 on chr2
    ]
    for ref_id, start, end, strand in expected:
        hit = _lookup(
            arrays,
            _region_bound_index(region_bounds, offsets, ref_id, start),
            _region_bound_index(region_bounds, offsets, ref_id, end),
        )
        assert hit is not None, f"intron [{start},{end}) on ref {ref_id} not found"
        slot, got_strand = hit
        assert 0 <= slot < arrays.boundary_right.shape[0]
        row = index.edges_df.iloc[int(arrays.edge_row[slot])]  # the JOIN, not the id
        assert int(row["kind"]) == EDGE_KIND_SJ
        assert got_strand == int(strand)


def test_a_shared_donor_keeps_BOTH_sj(index):
    """Alternative 3' splice site: region_bound 400 on chr1 is the donor of two distinct sj.

    A table storing one sj per region_bound would silently drop one of them, and the loss would be
    invisible — the dropped intron would simply be treated as unannotated.
    """
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    donor = _region_bound_index(region_bounds, offsets, 0, 400)
    lo, hi = int(arrays.offsets[donor]), int(arrays.offsets[donor + 1])
    assert hi - lo == 2
    landed = sorted(int(region_bounds[arrays.boundary_right[k]]) for k in range(lo, hi))
    assert landed == [700, 1000]


def test_a_region_bound_that_is_not_a_LEFT_BOUNDARY_has_an_empty_slice(index):
    """Most region_bounds are not left_boundaries — measured 70.4 % on both the toy and the human annotation. The slice
    must be empty rather than absent, so the deposit needs no special case."""
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    for ref_id, position in ((0, 200), (0, 1400), (1, 300)):  # TSS / TES, never a donor
        region_bound = _region_bound_index(region_bounds, offsets, ref_id, position)
        assert region_bound >= 0, f"{position} should be a region_bound on ref {ref_id}"
        assert int(arrays.offsets[region_bound + 1]) - int(arrays.offsets[region_bound]) == 0


def test_an_unannotated_intron_does_not_match(index):
    """A coordinate pair that is not an annotated sj must miss, even when both of its endpoints
    happen to be region bounds — the miss is what routes the fragment to the unspliced channel."""
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    # 700 and 900 are both region_bounds on chr1, but [700,900) is an EXON, not an intron
    assert (
        _lookup(
            arrays,
            _region_bound_index(region_bounds, offsets, 0, 700),
            _region_bound_index(region_bounds, offsets, 0, 900),
        )
        is None
    )
    # and a position that is not a region_bound at all
    assert _region_bound_index(region_bounds, offsets, 0, 401) == -1


def test_the_csr_round_trips_to_the_boundary_table(index):
    """Re-derive the sj set from the CSR and compare with ``edges_df`` — the two agree only if the
    region-id → region_bound-index shift is right on every reference, which is the one thing that can silently
    break when a reference has no regions."""
    arrays = build_sj_arrays(index)
    region_bounds, offsets, _ = build_region_partition_arrays(index)
    donor = np.repeat(np.arange(region_bounds.shape[0]), np.diff(arrays.offsets))
    from_csr = np.stack(
        [
            region_bounds[donor],
            region_bounds[arrays.boundary_right],
            arrays.strand.astype(np.int64),
        ],
        axis=1,
    )

    boundaries = index.edges_df
    sj = boundaries["kind"].to_numpy(np.uint8) == EDGE_KIND_SJ
    regions = index.regions_df
    src, dst = boundaries["src"].to_numpy(np.int64)[sj], boundaries["dst"].to_numpy(np.int64)[sj]
    from_boundaries = np.stack(
        [
            regions["end"].to_numpy(np.int64)[src],  # the intron starts where src ENDS
            regions["start"].to_numpy(np.int64)[dst],  # and ends where dst BEGINS
            boundaries["strand"].to_numpy(np.int8)[sj].astype(np.int64),
        ],
        axis=1,
    )
    order = lambda a: a[np.lexsort(a.T[::-1])]  # noqa: E731
    assert np.array_equal(order(from_csr), order(from_boundaries))


def _reference_partition(index):
    """The same graph, built through the reference accumulator's OWN constructor.

    Deliberately independent of :func:`build_sj_arrays`: this route names each sj by its
    genomic ``(ref, intron_start, intron_end, strand)`` and lets ``Partition.from_region_bounds`` resolve both
    endpoints with its own ``searchsorted``, where the builder walks region ids and applies a per-reference
    ``region_bound_base − region_base`` shift. Only the *definition* of a sj is shared.
    """
    region_bounds, region_bound_offsets, region_types = build_region_partition_arrays(index)
    boundaries, regions = index.edges_df, index.regions_df
    sj = boundaries["kind"].to_numpy(np.uint8) == EDGE_KIND_SJ
    src = boundaries["src"].to_numpy(np.int64)[sj]
    dst = boundaries["dst"].to_numpy(np.int64)[sj]
    strand = boundaries["strand"].to_numpy(np.int8)[sj]
    region_end, region_start = (
        regions["end"].to_numpy(np.int64),
        regions["start"].to_numpy(np.int64),
    )
    ref_of_region = regions["ref_name"].to_numpy()
    ref_id = {name: i for i, name in enumerate(index.ref_names)}

    sj = [
        # the intron starts where src ENDS and ends where dst BEGINS
        (ref_id[ref_of_region[s]], int(region_end[s]), int(region_start[d]), int(st))
        for s, d, st in zip(src, dst, strand)
    ]
    n_refs = len(index.ref_names)
    return Partition.from_region_bounds(
        [
            region_bounds[region_bound_offsets[r] : region_bound_offsets[r + 1]]
            for r in range(n_refs)
        ],
        # a reference contributing c region_bounds owns c-1 regions, so r earlier references own region_bound_offsets[r]-r
        region_types=[
            region_types[region_bound_offsets[r] - r : region_bound_offsets[r + 1] - r - 1]
            for r in range(n_refs)
        ],
        sj=sj,
    )


@pytest.mark.parametrize("fixture", ["index", "coincident_index"])
def test_the_csr_slot_order_matches_the_reference_accumulator(fixture, request):
    """THE CONTRACT: the sj-boundary id IS the CSR slot, so both builders must emit the slots in the
    SAME order — otherwise every sj row permutes and the native build's byte-identity gate compares
    two different labellings of one graph.

    This had never been tested. The two orderings disagreed once during S2 — ``(acceptor, donor)``
    against ``(strand, acceptor, donor)`` — and nothing would have caught it: the spec matrix exercises
    only ``Partition.from_region_bounds``, and the real-data shim builds its ``Partition`` straight from
    ``build_sj_arrays``, so the two sorts were never compared to each other.

    What this test does and does not cover, measured by perturbing the builder's key:

    * promoting ``strand`` to the primary key, or swapping donor/acceptor priority → caught, but only
      because of chr3's nested pair (see ``GTF``).
    * *dropping* ``strand`` from the builder's key → NOT caught, and cannot be. ``edges_df`` emits a
      strand-coincident pair POS-before-NEG whichever order the GTF lists the genes in, and ``np.lexsort``
      is stable, so both routes start from an already-correct tie order. The builder's ``strand`` key is
      therefore *defensive* — keep it, because it makes the contract explicit instead of resting on
      ``edges_df``'s internal sort, but do not expect this test to defend it.
    * ``from_region_bounds``'s ``strand`` key is load-bearing, since a caller may pass sj in any order,
      and it is pinned by ``test_a_sj_id_is_a_function_of_the_PARTITION_not_of_argument_order`` in
      the spec matrix.
    """
    index = request.getfixturevalue(fixture)
    arrays = build_sj_arrays(index)
    reference = _reference_partition(index)

    assert np.array_equal(reference.sj_offsets, arrays.offsets)
    assert np.array_equal(reference.sj_boundary_right, arrays.boundary_right)
    assert np.array_equal(reference.sj_strand, arrays.strand)


def test_a_strand_coincident_pair_is_two_distinct_slots(coincident_index):
    """Two genes on opposite strands sharing their intron coordinates exactly.

    Both must survive as separate sj boundaries, in the same CSR slice and ordered by strand — that
    adjacency is the only thing that makes the ordering contract above falsifiable, since every other
    sj is already separated by its donor or its acceptor.
    """
    arrays = build_sj_arrays(coincident_index)
    region_bounds, offsets, _ = build_region_partition_arrays(coincident_index)
    donor = _region_bound_index(region_bounds, offsets, 0, 400)
    lo, hi = int(arrays.offsets[donor]), int(arrays.offsets[donor + 1])
    assert hi - lo == 2, "the strand-coincident pair collapsed to one sj boundary"
    assert [int(region_bounds[arrays.boundary_right[k]]) for k in range(lo, hi)] == [800, 800]
    assert [int(arrays.strand[k]) for k in range(lo, hi)] == [int(Strand.POS), int(Strand.NEG)]
