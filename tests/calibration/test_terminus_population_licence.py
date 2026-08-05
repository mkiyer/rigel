"""⭐⭐⭐ THE POPULATION HALF OF THE COMPOSITION LICENCE — "is the source measuring the same thing I am?"

The reframe between two adjacent objects is a COMPOSITION IMPUTATION (`EQUATIONS.md` §3.5): it delivers
the source's density SHARE applied to the destination's observed total, so it is exact only where the two
share a composition. Attributing the density discrepancy between them to hybrid-capture enrichment is
therefore licensed only if they are measuring the same RNA population — otherwise enrichment and a
difference in the measured population are indistinguishable, and the imputation reads one as the other.

⭐ **The structural fact.** An EDGE is a single genomic position and it counts the fragments spanning it
CONTIGUOUSLY, so the transcripts it can see are exactly those continuous across it::

    T(EDGE)  =  T(NODE_left)  ∩  T(NODE_right)

A transcript whose body begins at the EDGE is in the right flank and not in the EDGE; one that ends there
is in the left flank and not in the EDGE. So a transcript TERMINUS at the EDGE is precisely what makes one
flank's population strictly larger, and the licence is an **equality** of population per (EDGE, side)
pair — not a containment, because ``phi_R`` too high and too low both corrupt ``phi_g``.

⛔⛔ **AND THE PREDICATE MUST BE WRITTEN IN GENOMIC TERMS, WHICH IS THE POINT OF HALF THIS FILE.** TSS/TES
is transcript-relative and the strand flips it: a ``+`` transcript's body extends UP from its TSS, a ``−``
transcript's UP from its TES. So::

    the RIGHT flank gains RNA at an EDGE  <=>  a transcript's genomic LOW  end is there
                                          <=>  FLAG_TSS_POS  or  FLAG_TES_NEG
    the LEFT  flank gains RNA at an EDGE  <=>  a transcript's genomic HIGH end is there
                                          <=>  FLAG_TES_POS  or  FLAG_TSS_NEG

⭐ The mirror gate below is what pins that: two indexes with **identical geometry on opposite strands**
must give the **same** flank answer, while the TSS bit alone points at the opposite EDGEs. A predicate
spelled "is there a TSS here" passes on one index and fails on its mirror.

⛔ **TERMINI ONLY.** A DONOR/ACCEPTOR EDGE also changes the population — RNA splices out or in — but
there the flux is MEASURED (``junction_count``) and the graft and the peel exist to route it. A terminus
has no flux to measure: a transcript simply begins. That is the derived boundary between the two
treatments, and extending this licence to splice sites is a separate experiment (`ROADMAP.md` §1 step 4).

⭐ Measured on `toy_harness.py --spec nested_exons`, where every EDGE carries a terminus: the gene's
mass-weighted ``|Δf_g|`` goes **0.2264 → 0.0541**, and its strand mirror ``nested_exons_neg``
**0.2145 → 0.0535**. ⚠ Both are BELOW the 0.0760 that deleting the relay's mass pin outright reaches, so
this is not a slower route to the same place: withholding the licence stops the reframe as well as the pin.

⚠ **Perturbation-verified** (`scratchpad/perturb_b.sh`), and one of these found a real hole:

    swap the two genomic-end masks                        5 gates
    TSS/TES bits instead of genomic ends                  3 — including the MIRROR gate, on the − index only
    swap the sides in the pair algebra                    1
    drop the test from the RELAY only (twin drift)        ⛔ **0 — until
                                                          `test_the_RELAY_honours_the_population_test_TOO`
                                                          was written for it.** Every other gate here reads
                                                          the COMBINE's capture and cannot see `_relay`
    drop the test from the COMBINE only (twin drift)      2
    the backward direction not shifted through `right`    2
    extend the masks to DONOR/ACCEPTOR (out of scope)     1, and only incidentally (a `worst_objects`
                                                          instrument row) — so
                                                          `test_a_SPLICE_SITE_alone_breaks_no_population_here`
                                                          exists to pin the scope directly
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.node_geometry import (
    build_node_statics,
    init_beliefs,
    terminus_flank_gain,
)
from rigel.calibration.splice_graph import (
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
    build_edge_flags_array,
    build_node_partition_arrays,
)

from conftest import build_test_index

#: ⭐⭐ THE OWNER'S `nested_exons` SPEC, as a real annotation: three NESTED single-exon transcripts, so
#: the gene is 7 NODEs and 6 EDGEs with no intron anywhere and the RNA population is a strict staircase.
#: Every EDGE therefore carries a terminus, and which FLANK it implicates is known by hand:
#:
#:      @1,000  @2,000  @3,000   a transcript's genomic LOW  end  -> the RIGHT flank gains
#:      @8,000  @9,000  @10,000  a transcript's genomic HIGH end  -> the LEFT  flank gains
#:
#: ⚠ GTF is 1-based inclusive, so ``[1000, 10000)`` is written ``1001 10000``.
_EXONS = ((1000, 10000), (2000, 9000), (3000, 8000))
_LOW_ENDS = {1000, 2000, 3000}
_HIGH_ENDS = {8000, 9000, 10000}
REFS = {"chr1": 12000}


def _gtf(strand: str) -> str:
    return "".join(
        f'chr1\ttest\texon\t{lo + 1}\t{hi}\t.\t{strand}\t.\tgene_id "g1"; transcript_id "t{i}";\n'
        for i, (lo, hi) in enumerate(_EXONS)
    )


@pytest.fixture(scope="module")
def pos_index(tmp_path_factory):
    return build_test_index(tmp_path_factory, _gtf("+"), name="nested_pos", refs=REFS)


@pytest.fixture(scope="module")
def neg_index(tmp_path_factory):
    """⭐ The STRAND MIRROR — the same three nested transcripts on the minus strand."""
    return build_test_index(tmp_path_factory, _gtf("-"), name="nested_neg", refs=REFS)


def _edge_positions(index) -> np.ndarray:
    """Genomic position of every contiguous edge, in edge order (the interior cuts, per reference)."""
    positions, cut_offsets, _types = build_node_partition_arrays(index)
    out = []
    for f in range(len(cut_offsets) - 1):
        lo, hi = int(cut_offsets[f]), int(cut_offsets[f + 1])
        if hi - lo >= 2:
            out.append(positions[lo + 1 : hi - 1])
    return np.concatenate(out) if out else np.zeros(0, np.int64)


# ── 1. THE PREDICATE, against a SECOND ALGORITHM ────────────────────────────────────────────────────


@pytest.mark.parametrize("which", ("pos", "neg"))
def test_the_flank_that_gains_is_the_one_the_transcript_SPANS_say_it_is(
    which, pos_index, neg_index
):
    """⭐⭐ **THE SET ALGEBRA, re-derived by a DIFFERENT ALGORITHM and compared.** The predicate reads
    structural BITS off the graph; this recomputes ``T(EDGE) = T(left) ∩ T(right)`` from the transcript
    spans this file declared by hand, and the two can then only agree by both being right — the same
    discipline as the splice graph's own I3b/I13.

    ⚠ Run on BOTH strands from one body of code, because the span arithmetic is strand-free: which flank
    gains is a fact about coordinates, and if the flag-derived answer needed a strand branch the two
    would disagree on one of the two indexes."""
    index = pos_index if which == "pos" else neg_index
    flags = build_edge_flags_array(index)
    right_gains, left_gains = terminus_flank_gain(flags)
    pos = _edge_positions(index)
    assert pos.shape == flags.shape

    for e, p in enumerate(pos):
        # the SECOND algorithm: which transcripts occupy the bases either side of this position, and
        # which span it contiguously. Half-open spans, from the hand-written table above.
        t_left = {i for i, (lo, hi) in enumerate(_EXONS) if lo <= p - 1 < hi}
        t_right = {i for i, (lo, hi) in enumerate(_EXONS) if lo <= p < hi}
        t_edge = t_left & t_right
        assert bool(right_gains[e]) == (t_edge != t_right), (
            f"@{p}: right_gains={bool(right_gains[e])} but T(EDGE)={sorted(t_edge)} vs "
            f"T(right)={sorted(t_right)}"
        )
        assert bool(left_gains[e]) == (t_edge != t_left), (
            f"@{p}: left_gains={bool(left_gains[e])} but T(EDGE)={sorted(t_edge)} vs "
            f"T(left)={sorted(t_left)}"
        )


@pytest.mark.parametrize("which", ("pos", "neg"))
def test_the_flank_answer_is_the_hand_written_one(which, pos_index, neg_index):
    """⭐ The same statement as a literal table, so a failure names the coordinate rather than a set."""
    index = pos_index if which == "pos" else neg_index
    right_gains, left_gains = terminus_flank_gain(build_edge_flags_array(index))
    pos = _edge_positions(index)
    assert set(pos[right_gains].tolist()) == _LOW_ENDS
    assert set(pos[left_gains].tolist()) == _HIGH_ENDS


def test_THE_STRAND_MIRROR_gives_the_SAME_flanks_while_the_TSS_BIT_POINTS_THE_OTHER_WAY(
    pos_index, neg_index
):
    """⛔⛔ **THE TRAP, made executable — the single biggest hazard in this predicate.** The two indexes
    are the same three nested transcripts at the same coordinates on opposite strands, so the flank that
    carries RNA the EDGE cannot see is the SAME one in both: it is a fact about which coordinate direction
    a transcript's body extends into.

    ⭐ And the second half is what makes it a trap rather than a coincidence: the TSS bit alone points at
    the **opposite** EDGEs between the two indexes. So a predicate written in TSS/TES terms would be
    exactly inverted on one of them — silently, since both are valid annotations."""
    p_flags = build_edge_flags_array(pos_index)
    n_flags = build_edge_flags_array(neg_index)
    p_pos, n_pos = _edge_positions(pos_index), _edge_positions(neg_index)
    assert np.array_equal(p_pos, n_pos), "the mirror indexes are not the same partition"

    p_right, p_left = terminus_flank_gain(p_flags)
    n_right, n_left = terminus_flank_gain(n_flags)
    assert np.array_equal(p_right, n_right), (
        f"the RIGHT-gaining EDGEs differ between the mirrors: {p_pos[p_right]} vs {n_pos[n_right]}"
    )
    assert np.array_equal(p_left, n_left), (
        f"the LEFT-gaining EDGEs differ between the mirrors: {p_pos[p_left]} vs {n_pos[n_left]}"
    )

    # …and the naive TSS/TES reading is INVERTED between them, which is why the above is not free.
    tss_p = set(p_pos[(p_flags & FLAG_TSS_POS) != 0].tolist())
    tss_n = set(n_pos[(n_flags & FLAG_TSS_NEG) != 0].tolist())
    assert tss_p == _LOW_ENDS and tss_n == _HIGH_ENDS, (
        f"the TSS bits are not where this gate assumes (+ {sorted(tss_p)}, − {sorted(tss_n)})"
    )
    assert tss_p != tss_n, "the mirrors' TSS bits agree, so this gate proves nothing about the trap"
    tes_p = set(p_pos[(p_flags & FLAG_TES_POS) != 0].tolist())
    tes_n = set(n_pos[(n_flags & FLAG_TES_NEG) != 0].tolist())
    assert tes_p == _HIGH_ENDS and tes_n == _LOW_ENDS


def test_a_flagless_edge_breaks_no_population():
    """⭐ The empty case, stated so a future all-zeros plumbing bug cannot read as "everything matches".
    ``edge_flags`` is 0 when no graph was supplied, and 0 must mean *no terminus here*, so both flanks
    match — which is the permissive answer and therefore the one that has to be deliberate."""
    right_gains, left_gains = terminus_flank_gain(np.zeros(4, np.uint16))
    assert not right_gains.any() and not left_gains.any()


def test_a_SPLICE_SITE_alone_breaks_no_population_here():
    """⛔ **THE SCOPE, pinned so it cannot be widened by accident.** A DONOR/ACCEPTOR EDGE does change the
    RNA population — a transcript splices out or in — but there the flux is MEASURED (``junction_count``)
    and the graft and the peel exist to route it, so routing it a second time as a withheld licence would
    double-count. This predicate must therefore read the four TERMINUS bits and nothing else.

    ⚠ Extending it to splice sites is a separate experiment with its own measurement: deleting the
    junction channel was measured *better* under capture and *worse* off it (`ROADMAP.md` §1 step 4)."""
    from rigel.calibration.splice_graph import (
        FLAG_ACCEPTOR_NEG,
        FLAG_ACCEPTOR_POS,
        FLAG_DONOR_NEG,
        FLAG_DONOR_POS,
    )

    splice_only = np.array(
        [FLAG_DONOR_POS, FLAG_DONOR_NEG, FLAG_ACCEPTOR_POS, FLAG_ACCEPTOR_NEG], np.uint16
    )
    right_gains, left_gains = terminus_flank_gain(splice_only)
    assert not right_gains.any() and not left_gains.any(), (
        "a splice site is breaking a population in this predicate; DONOR/ACCEPTOR are out of scope and "
        "belong to the graft and the peel"
    )
    # …and a terminus sharing the position with a splice site is still seen (the measured 0.99 % case).
    both = np.array([FLAG_TSS_POS | FLAG_DONOR_POS], np.uint16)
    assert terminus_flank_gain(both)[0].all()


def test_the_two_masks_are_independent_and_an_edge_can_carry_BOTH():
    """⚠ They are not complements. A transcript can end where another begins, so one EDGE can break the
    population on BOTH sides — and then neither flank is licensed. A predicate built as "which side?"
    rather than as two independent bits would have to pick one."""
    both = np.array([FLAG_TSS_POS | FLAG_TES_POS, FLAG_TES_NEG | FLAG_TSS_NEG], np.uint16)
    right_gains, left_gains = terminus_flank_gain(both)
    assert right_gains.all() and left_gains.all()


# ── 2. THE WIRING — the licence per STEP, on a chain the solver actually walks ──────────────────────


def _flagged_chain(flags_by_edge):
    """The `nested_exons` chain as a synthetic payload, with ``edge_flags`` injected per EDGE.

    ⭐ STRANDED (κ = 0.95 at the solve), because the population test only *decides* anything where the
    source could otherwise lend a composition: on an unstranded chain with no junction nothing is
    licensed to begin with and the new conjunct is inert by construction.
    """
    from rigel.calibration.signature import BIT_EXON_POS

    from _synthetic import make_chain_parts

    parts = make_chain_parts(
        [0] + [BIT_EXON_POS] * 5 + [0],
        node_size_bp=1000.0,
        node_pos=[100.0, 900.0, 900.0, 900.0, 900.0, 900.0, 100.0],
        node_neg=[100.0, 50.0, 50.0, 50.0, 50.0, 50.0, 100.0],
        edge_pos=[20.0, 60.0, 60.0, 60.0, 60.0, 20.0],
        edge_neg=[20.0, 10.0, 10.0, 10.0, 10.0, 20.0],
    )
    parts.statics = build_node_statics(
        parts.chain, parts.region_arrays, np.asarray(flags_by_edge, np.uint16)
    )
    return parts


def _hop_licence(parts):
    """Every (destination slot, source slot) hop the COMBINE transports, with its licence.

    Reads ``_capture['_pin']``, which publishes the gDNA scale rule per hop — ``lend``, ``r`` and the
    scale ``r_g`` actually applied — so the gate reads the rule off a real solve rather than re-deriving
    it. The two appended entries are the left-source and right-source messages, in that order.
    """
    from rigel.calibration.bp_solver import node_sweep

    cap = {}
    node_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=200),
        parts.region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10_000.0,
        n_gdna_obs=10_000.0,
        n_grid=200,
        _capture=cap,
    )
    out = []
    for d in cap["_pin"]:
        src, valid = np.asarray(d["src"], np.int64), np.asarray(d["valid"], bool)
        lend, r, r_g = (np.asarray(d[k]) for k in ("lend", "r", "r_g"))
        for dst in range(src.shape[0]):
            if valid[dst]:
                out.append(
                    {
                        "dst": dst,
                        "src": int(src[dst]),
                        "lend": bool(lend[dst]),
                        "r": float(r[dst]),
                        "r_g": float(r_g[dst]),
                    }
                )
    return out


#: chain N E N E N E N E N E N E N — slot 2k is node k, slot 2k+1 is edge k.
#: The six EDGEs are slots 1,3,5,7,9,11; the `nested_exons` flag pattern puts a genomic LOW end on the
#: first three and a genomic HIGH end on the last three.
_NESTED_FLAGS = (
    FLAG_TSS_POS,
    FLAG_TSS_POS,
    FLAG_TSS_POS,
    FLAG_TES_POS,
    FLAG_TES_POS,
    FLAG_TES_POS,
)
_NO_FLAGS = (0,) * 6


def test_a_terminus_UNLICENSES_the_step_into_the_flank_that_gains_rna():
    """⭐⭐⭐ **THE CHANGE, as a per-step statement.** At an EDGE where a transcript's genomic LOW end
    sits, the step between that EDGE and its RIGHT flank crosses a population difference, so the
    composition may not be imputed across it and the gDNA level must cross UNSCALED (``r_g == 1``). At an
    EDGE carrying a genomic HIGH end, the same for its LEFT flank. Both directions of each pair, because
    the populations differ symmetrically — the question is about the pair, not about who is talking.

    ⚠ **And the steps a terminus does NOT touch must be unaffected**, which is the half that stops the
    conjunct being a blanket. On this chain each EDGE breaks exactly one of its two sides, so every EDGE
    contributes one broken pair and one intact pair."""
    licence = {(h["dst"], h["src"]): h for h in _hop_licence(_flagged_chain(_NESTED_FLAGS))}
    # slot 2k = node k, slot 2k+1 = edge k
    broken, intact = [], []
    for e, bits in enumerate(_NESTED_FLAGS):
        edge = 2 * e + 1
        gains_right = bits == FLAG_TSS_POS
        far = edge + 1 if gains_right else edge - 1
        near = edge - 1 if gains_right else edge + 1
        broken += [(edge, far), (far, edge)]
        intact += [(edge, near), (near, edge)]
    assert all(p in licence for p in broken + intact), (
        "the chain's hops are not what this gate assumes"
    )

    for dst, src in broken:
        h = licence[(dst, src)]
        assert not h["lend"], (
            f"step {src} → {dst} crosses a transcript terminus into the flank that gains RNA, yet it is "
            f"still licensed to lend a composition (r_g {h['r_g']:.4g}, r {h['r']:.4g})"
        )
        assert h["r_g"] == 1.0, f"step {src} → {dst}: the gDNA level was scaled by {h['r_g']:.6g}"
    assert any(licence[p]["lend"] for p in intact), (
        "no step survives the population test — the conjunct has become a blanket refusal, and this "
        "chain has an intact pair at every EDGE"
    )


def _ambig_chain(flags_by_edge):
    """``intergenic | exon+ | exon AMBIG | exon+ | intergenic``, STRANDED, flags injected per EDGE.

    ⭐ **Built for READABILITY of the relay's own output.** An exact statement about a relayed level needs
    a destination with no own gDNA precision — otherwise the fuse mixes its belief into the running value
    — and on a stranded chain every single-strand genic slot earns some. An **AMBIG** exon is the
    exception: both strands are admissible, so it is G3, holds ``{0,0,1}`` at MAX variance, and has
    precision 0 while NOT being structurally pure gDNA. Its flanking EDGE is single-strand (the ``+``
    neighbour has no ``−`` transcript), so the EDGE does earn evidence and can lend a composition.
    """
    from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS

    from _synthetic import make_chain_parts

    parts = make_chain_parts(
        [0, BIT_EXON_POS, BIT_EXON_POS | BIT_EXON_NEG, BIT_EXON_POS, 0],
        node_size_bp=1000.0,
        node_pos=[100.0, 900.0, 500.0, 900.0, 100.0],
        node_neg=[100.0, 50.0, 500.0, 50.0, 100.0],
        edge_pos=[20.0, 60.0, 60.0, 20.0],
        edge_neg=[20.0, 10.0, 10.0, 20.0],
    )
    parts.statics = build_node_statics(
        parts.chain, parts.region_arrays, np.asarray(flags_by_edge, np.uint16)
    )
    return parts


def _relay_levels(parts):
    """``(pg_own, forward running gDNA level)`` per slot — the RELAY's own output, before the combine."""
    from rigel.calibration.bp_solver import node_sweep

    cap = {}
    node_sweep(
        parts.chain,
        parts.statics,
        parts.geometry,
        init_beliefs(parts.chain, parts.geometry, parts.statics, rna_sense_frac=0.95, n_grid=200),
        parts.region_arrays,
        rna_sense_frac=0.95,
        n_rna_obs=10_000.0,
        n_gdna_obs=10_000.0,
        n_grid=200,
        _capture=cap,
    )
    st = cap["_uni_static"]
    return np.asarray(st["pg_own"], float), np.asarray(st["fwd_g"], float)


def test_the_RELAY_honours_the_population_test_TOO():
    """⭐⭐⭐ **THE TWIN-DRIFT GATE, and it is not decoration: without it, deleting the population test
    from ``_relay`` alone fires NOTHING in the whole calibration suite (measured).** ``_relay`` and
    ``_transport`` are two hand-maintained copies of one transform — the DO-NOT-MERGE note in `bp_solver`
    exists for exactly this — and every other gate in this file reads ``_capture['_pin']``, which is the
    COMBINE's. This one reads ``fwd_g``, the relay's running level, which the combine has not yet touched.

    The chain is ``intergenic | exon+ | exon AMBIG | exon+ | intergenic`` with a genomic LOW end on the
    EDGE at slot 3, so the step from that EDGE into the AMBIG node at slot 4 crosses a population
    difference. The AMBIG destination has no own gDNA precision, so the relayed claim survives the fuse
    unchanged and the statement is exact: the level must cross UNSCALED.

    ⚠ Two-sided against the same chain with the flags cleared, where the step IS licensed and the level
    must move — otherwise the gate would pass on a fixture where the two branches happen to agree."""
    edge, dst = 3, 4  # chain N E N E N E N E N → EDGE 1 is slot 3, its RIGHT flank is slot 4
    flags = [
        0,
        FLAG_TSS_POS,
        0,
        0,
    ]  # a transcript's genomic LOW end at EDGE 1 ⇒ the right flank gains
    own_f, lvl_f = _relay_levels(_ambig_chain(flags))
    own_p, lvl_p = _relay_levels(_ambig_chain(_NO_FLAGS[:4]))
    assert own_f[dst] == 0.0 and own_p[dst] == 0.0, (
        f"slot {dst} has own gDNA precision ({own_f[dst]:.4g}), so the relayed claim is not readable "
        "there and this fixture has drifted"
    )
    assert lvl_f[dst] == lvl_f[edge], (
        f"the relay scaled the level {lvl_f[edge]:.9g} → {lvl_f[dst]:.9g} across a step that gains RNA "
        "on the destination side, so `_relay` is not applying the population test"
    )
    assert lvl_p[dst] != lvl_p[edge], (
        f"with the flags cleared the same step left the level at {lvl_p[dst]:.9g}, so this chain cannot "
        "tell the population test from a step that was never licensed anyway"
    )


def test_the_conjunct_is_INERT_when_no_edge_carries_a_terminus():
    """⭐ **THE OTHER SIDE: with the flags all zero the licence must be exactly what it was.** Same
    chain, same numbers, ``edge_flags = 0`` — so every hop's ``lend`` is decided by the supply test
    alone. This is what says the new conjunct adds a restriction rather than changing the old one, and it
    is also the state a caller that supplies no graph gets."""
    flagged = {(h["dst"], h["src"]): h for h in _hop_licence(_flagged_chain(_NESTED_FLAGS))}
    plain = {(h["dst"], h["src"]): h for h in _hop_licence(_flagged_chain(_NO_FLAGS))}
    assert set(flagged) == set(plain)
    assert any(plain[p]["lend"] and not flagged[p]["lend"] for p in plain), (
        "the flags changed no hop's licence, so this fixture cannot see the conjunct at all"
    )
    assert all(plain[p]["lend"] or not flagged[p]["lend"] for p in plain), (
        "a hop became licensed BECAUSE of a terminus flag — the conjunct is meant to remove licences"
    )
