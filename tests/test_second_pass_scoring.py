"""⭐ THE P2 GATE — the second pass's score DISCRIMINATES, on loci where the truth is built in.

     (the score) and its P2 · Measurement:
    Subject: ``rigel.second_pass.score_held_fragments``, which shipped in `src/` with **no test at all**.

The gate P2 asks for is *"on a hand-built locus where the truth is known, the correct hypothesis takes the
larger share"*. ⛔ **The hard part is not the assertion, it is the fixture.** records why
this was deliberately not written against the four-fragment smoke fixture: it had no depth, every score
came out uniform, and the gate would have been green over a scorer that decided nothing.

⭐ **HOW rho IS ISOLATED WITHOUT A THRESHOLD — the mirror.** Loci 1 and 2 are the *same geometry* with the
deep sj swapped. Every hypothesis has the same implied ``L`` at both, so ``f`` and the strand term
are **identical between them** — measured, not assumed (:func:`test_the_two_arms_differ_ONLY_in_rho`). So
if the winner flips between the two loci, only ``rho`` can have flipped it. That is a derivation; a
"score(A) > 3 * score(B)" assertion would have been a tuned constant.

The five arms
-------------

======  ==========================================  ==========================================
arm 1   locus 1, the **wide** sj is deep      the wide hypothesis must take the larger share
arm 2   locus 2, the **narrow** sj is deep    ⭐ the mirror — the answer must FLIP
arm 3   locus 3, the gap is deeply crossed          ∅ must win. ⛔ **This arm is D-6**, and it
        **contiguously** and nothing splices it     failed until ∅'s evidence set was corrected
arm 4   locus 4, the gap's **donor** end is deeply  ∅ must LOSE — a molecule crossing the gap is
6       crossed and its acceptor never; locus 6 is  present at both ends, and nothing was seen at
        the mirror                                  one. ⛔ Added because half-fixes to D-6 that
                                                    kept ONE boundary passed arms 1–3
arm 5   locus 5, two sj at **equal** depth   the length term decides, so the hypothesis the
        and different implied lengths               anchor supports must win. ⛔ Added because
                                                    dropping ``f`` entirely passed arms 1–4
arm 7   locus 7, two **opposite-strand** hypotheses  scored twice on one payload, at an R1-sense and
        of equal width and equal depth              an R1-antisense library — the winner must flip.
                                                    ⛔ Added because dropping ``s`` passed arms 1–6
======  ==========================================  ==========================================

⚠ **The strand model here is deliberately prior-only** (``p_r1_sense = 0.5``): no read carries an ``XS``
tag, so there is no strand information and 0.5 is the correct answer rather than a degenerate one. That
is on purpose — it holds the strand term constant across each locus's two spliced hypotheses so this
module tests ``rho`` and ``f``, and `tests/native/test_gap_hypothesis_strand.py` already gates the strand
behaviour on its own. The scan warns about it, and the warning is correct.

⛔ **The fixture goes through the SCAN.** 's trap: driving ``FragmentResolver`` directly
leaves ``t_strand_arr_`` empty and every hypothesis's implied strand silently reads ``NONE``.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer

GENOME = 30_000
M, N = 0, 3

#: The six loci. Three MIRROR PAIRS: 1/2 swap which sj is deep, and 4/6 swap which END of the
#: gap is crossed. Every non-firing perturbation this module found was closed by adding the mirror.
L1, L2, L3, L4, L5, L6, L7 = 1000, 4000, 7000, 10_000, 13_000, 16_000, 19_000

#: Depths. ⚠ Fixture quantities, not algorithm constants: the only thing that matters is that DEEP and
#: SHALLOW are far enough apart to separate and that both are **non-zero**, so the gate tests
# discrimination rather than the elimination-by-zero that measured separately.
DEEP, SHALLOW, CONTIGUOUS = 40, 4, 8

#: ⭐ Ballast sets the LENGTH pmf and touches no locus. `build_fl_models` does **not** smooth the global
#: anchor (`_normalized`, not `_smooth_eb`), so ``global_pmf[L]`` is exactly 0 unless some deposited
#: fragment had that very length — and every ∅ hypothesis here is 600 bp. Separating the pmf's source
#: from the loci is what stops a length coincidence from doing rho's work.
#:
#: ⭐ ``(half, count)``, so ``L = 2 * half``. Adding mass at one length rescales every ``global_pmf``
#: entry by a common factor, which cancels in a score normalised within a record — so the ballast sets
#: length RATIOS and disturbs nothing else.
#:
#: ⛔ **Arm 5's two lengths (300 and 500) are reachable by NO locus**, because every depth fragment's
#: length comes from ``_SHAPES`` and those are 200/400/600. So arm 5's 3:1 ratio is a property of this
#: table alone and cannot drift when a locus is added — which is exactly how it broke once, when arm 7's
#: depth fragments quietly added mass at 200.
BALLAST = ((100, 60), (150, 60), (200, 60), (250, 180), (300, 60))

#: The implied lengths. Arms 1–4 put a 400 bp gap in a 600 bp span; arm 5's introns are 300 and 100 bp.
L_WIDE, L_NARROW, L_GENOMIC = 200, 400, 600
L_ARM5_WIDE, L_ARM5_NARROW = 300, 500


def _gtf() -> str:
    rows = []
    for base, wide, narrow in ((L1, "A1", "B1"), (L2, "A2", "B2")):
        rows += [
            f'chr1\tt\texon\t{base + 1}\t{base + 600}\t.\t+\t.\tgene_id "g{wide}"; transcript_id "t{wide}";\n',
            f'chr1\tt\texon\t{base + 1001}\t{base + 1200}\t.\t+\t.\tgene_id "g{wide}"; transcript_id "t{wide}";\n',
            f'chr1\tt\texon\t{base + 1}\t{base + 700}\t.\t+\t.\tgene_id "g{narrow}"; transcript_id "t{narrow}";\n',
            f'chr1\tt\texon\t{base + 901}\t{base + 1200}\t.\t+\t.\tgene_id "g{narrow}"; transcript_id "t{narrow}";\n',
        ]
    # ⭐ Loci 3 and 4 have ONE isoform on purpose, so the intron's endpoints are ADJACENT region_bounds and the
    # intron spans exactly one region. That is D-6's case: the boundaries strictly between them are an empty set.
    for base, gene in ((L3, "C"), (L4, "D"), (L6, "G")):
        rows += [
            f'chr1\tt\texon\t{base + 1}\t{base + 600}\t.\t+\t.\tgene_id "g{gene}"; transcript_id "t{gene}";\n',
            f'chr1\tt\texon\t{base + 1001}\t{base + 1200}\t.\t+\t.\tgene_id "g{gene}"; transcript_id "t{gene}";\n',
        ]
    # Locus 5: two isoforms whose introns differ in WIDTH (300 vs 100), so their implied lengths are 300
    # and 500 — bins only the ballast fills. Both sj get the same depth, so only ``f`` can decide.
    rows += [
        f'chr1\tt\texon\t{L5 + 1}\t{L5 + 650}\t.\t+\t.\tgene_id "gE"; transcript_id "tE";\n',
        f'chr1\tt\texon\t{L5 + 951}\t{L5 + 1200}\t.\t+\t.\tgene_id "gE"; transcript_id "tE";\n',
        f'chr1\tt\texon\t{L5 + 1}\t{L5 + 700}\t.\t+\t.\tgene_id "gF"; transcript_id "tF";\n',
        f'chr1\tt\texon\t{L5 + 801}\t{L5 + 1200}\t.\t+\t.\tgene_id "gF"; transcript_id "tF";\n',
    ]
    # Locus 7: opposite-strand isoforms whose introns have the SAME width (200 bp), so the two
    # hypotheses tie on implied length AND — given equal depth — on rho. Only the strand can separate
    # them. ⚠ Their coordinates differ, so this is not D-5's strand-coincident case and the index is quiet.
    rows += [
        f'chr1\tt\texon\t{L7 + 1}\t{L7 + 600}\t.\t+\t.\tgene_id "gH"; transcript_id "tH";\n',
        f'chr1\tt\texon\t{L7 + 801}\t{L7 + 1200}\t.\t+\t.\tgene_id "gH"; transcript_id "tH";\n',
        f'chr1\tt\texon\t{L7 + 1}\t{L7 + 700}\t.\t-\t.\tgene_id "gI"; transcript_id "tI";\n',
        f'chr1\tt\texon\t{L7 + 901}\t{L7 + 1200}\t.\t-\t.\tgene_id "gI"; transcript_id "tI";\n',
    ]
    return "".join(rows)


def _read(qname, pos, cigar, mate_pos, is_r1):
    import pysam

    a = pysam.AlignedSegment()
    a.query_name = qname
    a.reference_id = 0
    a.reference_start = pos
    a.mapping_quality = 60
    a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10)
    a.cigar = cigar
    n = sum(length for op, length in cigar if op == M)
    a.query_sequence = "A" * n
    a.query_qualities = pysam.qualitystring_to_array("I" * n)
    a.next_reference_id = 0
    a.next_reference_start = mate_pos
    a.set_tags([("NH", 1, "i")])
    return a


#: Three block layouts, so the depth fragments span the same lengths the ballast does and no hypothesis
#: length is reachable only through one kind of evidence.
_SHAPES = ((100, 50, 50), (200, 100, 100), (300, 150, 150))


def _spliced(qname, donor, acceptor, k):
    """A fully-sequenced fragment OBSERVING ``(donor, acceptor)``. No mate gap, so it deposits."""
    a, b, c = _SHAPES[k % len(_SHAPES)]
    return [
        _read(qname, donor - a, [(M, a), (N, acceptor - donor), (M, b)], acceptor + b, True),
        _read(qname, acceptor + b, [(M, c)], donor - a, False),
    ]


def _contiguous(qname, start, half):
    """A fully-sequenced UNSPLICED fragment: it crosses every boundary strictly inside its span."""
    return [
        _read(qname, start, [(M, half)], start + half, True),
        _read(qname, start + half, [(M, half)], start, False),
    ]


def _ambiguous(qname, base):
    """The held fragment: ``[base+500, base+600)`` and ``[base+1000, base+1100)``, gap unsequenced."""
    return [
        _read(qname, base + 500, [(M, 100)], base + 1000, True),
        _read(qname, base + 1000, [(M, 100)], base + 500, False),
    ]


@pytest.fixture(scope="module")
def scored(tmp_path_factory):
    """Scan the fixture once and score its held fragments. Returns everything the gates read."""
    import pysam

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.splice_graph import (
        build_sj_arrays,
        build_region_partition_arrays,
    )
    from rigel.index import TranscriptIndex
    from rigel.second_pass import score_held_fragments

    base = tmp_path_factory.mktemp("second_pass_scoring")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(">chr1\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n")
    pysam.faidx(str(fasta))
    gtf.write_text(_gtf())
    index_dir = base / "idx"
    TranscriptIndex.build(str(fasta), str(gtf), str(index_dir), write_tsv=False)
    index = TranscriptIndex.load(str(index_dir))

    reads = []
    for i, (half, count) in enumerate(BALLAST):
        for k in range(count):
            reads += _contiguous(f"ballast_{i}_{k}", 23_000 + (k % 20) * 100, half)

    # arms 1 and 2: the same locus twice, with the deep sj swapped
    for locus, deep_offset, shallow_offset in (
        (L1, (600, 1000), (700, 900)),
        (L2, (700, 900), (600, 1000)),
    ):
        for k in range(DEEP):
            reads += _spliced(
                f"deep_{locus}_{k}", locus + deep_offset[0], locus + deep_offset[1], k
            )
        for k in range(SHALLOW):
            reads += _spliced(
                f"shallow_{locus}_{k}", locus + shallow_offset[0], locus + shallow_offset[1], k
            )
        for k in range(CONTIGUOUS):
            reads += _contiguous(f"contig_{locus}_{k}", locus + 500, 300)
        reads += _ambiguous(f"ambig_{locus}", locus)

    # arm 3: ∅ is the truth — the gap is crossed contiguously DEEP and spliced only SHALLOW
    for k in range(DEEP):
        reads += _contiguous(f"contig_{L3}_{k}", L3 + 500, 300)
    for k in range(SHALLOW):
        reads += _spliced(f"spliced_{L3}_{k}", L3 + 600, L3 + 1000, k)
    reads += _ambiguous(f"ambig_{L3}", L3)

    # arm 4: the DONOR end of the gap is deeply crossed and the ACCEPTOR end never is. These fragments
    # span [L4+400, L4+800) — over the boundary at L4+600, and stopping well short of the one at L4+1000.
    for k in range(DEEP):
        reads += _contiguous(f"contig_{L4}_{k}", L4 + 400, 200)
    for k in range(SHALLOW):
        reads += _spliced(f"spliced_{L4}_{k}", L4 + 600, L4 + 1000, k)
    reads += _ambiguous(f"ambig_{L4}", L4)

    # arm 6: the MIRROR of arm 4 — the ACCEPTOR end is deeply crossed and the DONOR end never is.
    # These span [L6+800, L6+1200), over the boundary at L6+1000 and starting well past the one at L6+600.
    for k in range(DEEP):
        reads += _contiguous(f"contig_{L6}_{k}", L6 + 800, 200)
    for k in range(SHALLOW):
        reads += _spliced(f"spliced_{L6}_{k}", L6 + 600, L6 + 1000, k)
    reads += _ambiguous(f"ambig_{L6}", L6)

    # arm 5: both sj at EQUAL depth — the same k range, so the deposited length multiset and
    # therefore the deposited density are identical. Only the implied length differs.
    for k in range(DEEP):
        reads += _spliced(f"wide_{L5}_{k}", L5 + 650, L5 + 950, k)
        reads += _spliced(f"narrow_{L5}_{k}", L5 + 700, L5 + 800, k)
    reads += _ambiguous(f"ambig_{L5}", L5)

    # arm 7: two opposite-strand hypotheses, equal width and equal depth — only strand can separate them
    for k in range(DEEP):
        reads += _spliced(f"plus_{L7}_{k}", L7 + 600, L7 + 800, k)
        reads += _spliced(f"minus_{L7}_{k}", L7 + 700, L7 + 900, k)
    reads += _ambiguous(f"ambig_{L7}", L7)

    bam = str(base / "fixture.bam")
    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": GENOME}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)

    _stats, strand_model, _buffer, payload = scan_and_buffer(
        bam, index, BamScanConfig(sj_strand_tag="XS", total_threads=1)
    )
    _, _, region_types = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)
    fl_models = build_fl_models(payload)

    def rescore(rna_sense_frac: float):
        """⭐ Re-score the SAME payload at a different library sense fraction. Isolates the strand term
        exactly: nothing else in the score can move, because nothing else is a function of it."""
        return score_held_fragments(
            payload,
            fl_models=fl_models,
            rna_sense_frac=rna_sense_frac,
            region_types=region_types,
            sj=sj,
        )

    return payload, rescore(strand_model.p_r1_sense), rescore


def _shares(payload, result, base: int) -> dict[tuple, float]:
    """The held record at ``base`` as ``{implied introns: normalised score}``; ``()`` is ∅."""
    d = payload.deferred
    matches = [i for i in range(d.n_fragments) if int(d.start[i]) == base + 500]
    assert len(matches) == 1, f"expected one held record at {base + 500}, got {matches}"
    i = matches[0]
    return {
        tuple(tuple(p) for p in d.hypothesis_introns_of(h).tolist()): float(result.score[h])
        for h in range(int(d.hypothesis_offsets[i]), int(d.hypothesis_offsets[i + 1]))
    }


def _terms(payload, result, base: int) -> dict[tuple, tuple]:
    """The same record's ``(L, f, strand)`` per hypothesis — everything except ``rho``."""
    d = payload.deferred
    i = next(i for i in range(d.n_fragments) if int(d.start[i]) == base + 500)
    t = result.terms
    return {
        tuple(tuple(p) for p in d.hypothesis_introns_of(h).tolist()): (
            int(t.length[h]),
            float(t.length_likelihood[h]),
            float(t.strand[h]),
        )
        for h in range(int(d.hypothesis_offsets[i]), int(d.hypothesis_offsets[i + 1]))
    }


# ── arms 1 and 2: the deep sj wins, and the answer FLIPS when the depth moves ────────────────


def test_the_STRAND_term_decides_when_rho_and_LENGTH_both_tie(scored):
    """⛔ **Arm 7.** Dropping the strand term entirely passed arms 1–6, because every locus there offers
    hypotheses of one strand and the term cancels. This locus offers two.

    ``tH`` (+) and ``tI`` (−) imply introns of the **same width** at different coordinates, observed by
    the same number of fragments — so implied length and rho both tie and only ``s`` is left.

    ⭐ **Scored twice on ONE payload**, at an R1-sense and an R1-antisense library. Nothing but the strand
    term is a function of ``rna_sense_frac``, so a flipped winner isolates it exactly — the same mirror
    argument arms 1/2 use for rho, with no second fixture.

    ⚠ **The direction is the point**. ``rna_sense_frac`` is ``P(align_strand agrees |
    RNA)``, and on a real dUTP cfRNA library it is ≈ 0.01 — so the hypothesis that **disagrees** with
    ``align_strand`` is the likely one. ⛔ A scorer written as "agreement ⇒ multiply by
    ``rna_sense_frac``" is backwards on every real library, and this gate is what says so.
    """
    payload, _default, rescore = scored
    plus, minus = ((L7 + 600, L7 + 800),), ((L7 + 700, L7 + 900),)
    sense = _shares(payload, rescore(0.99), L7)
    antisense = _shares(payload, rescore(0.01), L7)
    assert sense[plus] > sense[minus], (
        f"the held fragment aligned to +, and at rna_sense_frac=0.99 agreement is likely, so the + "
        f"hypothesis must win; got {sense}"
    )
    assert antisense[minus] > antisense[plus], (
        f"at rna_sense_frac=0.01 — a dUTP library, which is what real cfRNA is — DISAGREEMENT is the "
        f"likely case, so the − hypothesis must win; got {antisense}. An unmoved winner means the "
        f"strand term is not being used, and a winner that moved the other way means it is inverted."
    )

    # ⛔ AND ∅'s OWN TERM MUST BE SYMMETRIC HERE. This locus offers candidates on BOTH strands, so there is
    # no locus orientation for an unspliced molecule to be sense or antisense TO — and a rule that picked
    # one anyway would make the answer depend on which annotation boundary was read first, which is D-5's defect
    # wearing a different hat.
    d = payload.deferred
    i = next(j for j in range(d.n_fragments) if int(d.start[j]) == L7 + 500)
    result = rescore(0.99)
    for h in range(int(d.hypothesis_offsets[i]), int(d.hypothesis_offsets[i + 1])):
        if len(d.hypothesis_introns_of(h)) == 0:
            assert float(result.terms.strand[h]) == 0.5, (
                f"at a both-strands locus the genomic candidate has no sense direction, so its strand term "
                f"must be 0.5 by symmetry; got {float(result.terms.strand[h])}"
            )


def test_the_DEEP_sj_takes_the_larger_share(scored):
    """⭐ Arm 1. Locus 1's wide intron carries ``DEEP`` observed sj fragments and its narrow rival
    ``SHALLOW``; the held fragment is compatible with both. The wide hypothesis must win.

    ⚠ Both rivals have **non-zero** flux by construction, so this is a test of discrimination and not of
    the elimination-by-zero measured. That distinction is the whole reason for
    ``SHALLOW`` being 4 rather than 0.
    """
    shares = _shares(scored[0], scored[1], L1)
    wide, narrow = ((1600, 2000),), ((1700, 1900),)
    assert shares[wide] > shares[narrow], (
        f"the DEEP sj must take the larger share; got wide={shares[wide]:.4f} "
        f"narrow={shares[narrow]:.4f}. The two differ only in how many fragments were observed "
        f"crossing them."
    )
    assert shares[wide] == max(shares.values()), (
        f"the deep sj must be the winner; got {shares}"
    )


def test_MOVING_the_depth_FLIPS_the_answer(scored):
    """⭐ Arm 2, the mirror — and the reason no threshold is needed anywhere in this module.

    Locus 2 is locus 1's geometry with the deep sj moved to the *narrow* intron. Every implied
    length is unchanged, so ``f`` and the strand term are unchanged
    (:func:`test_the_two_arms_differ_ONLY_in_rho` proves it). A flipped winner therefore isolates ``rho``
    **exactly**, where a "wins by more than X" assertion would have been a constant chosen after the fact.
    """
    first = _shares(scored[0], scored[1], L1)
    second = _shares(scored[0], scored[1], L2)
    wide1, narrow1 = ((1600, 2000),), ((1700, 1900),)
    wide2, narrow2 = ((4600, 5000),), ((4700, 4900),)
    assert first[wide1] > first[narrow1] and second[narrow2] > second[wide2], (
        f"moving the depth from the wide intron to the narrow one must flip the winner; got "
        f"locus1={first} locus2={second}"
    )


def test_the_two_arms_differ_ONLY_in_rho(scored):
    """⛔ The mirror argument is a claim about the fixture, so it is CHECKED rather than asserted in prose.

    If ``f`` or the strand term differed between the two loci, a flipped winner would no longer isolate
    ``rho`` and :func:`test_MOVING_the_depth_FLIPS_the_answer` would be green for a reason it does not
    state.
    """
    first, second = _terms(scored[0], scored[1], L1), _terms(scored[0], scored[1], L2)
    pairs = ((((1600, 2000),), ((4600, 5000),)), (((1700, 1900),), ((4700, 4900),)), ((), ()))
    for a, b in pairs:
        assert first[a] == second[b], (
            f"hypothesis {a} at locus 1 has (L, f, strand) = {first[a]} but its mirror {b} at locus 2 "
            f"has {second[b]}. The two loci must differ in rho and in nothing else."
        )


# ── arm 3: the genomic path, when the gap is genuinely crossed ─────────────────────────────────────


def test_a_DEEPLY_CROSSED_gap_is_won_by_the_GENOMIC_hypothesis(scored):
    """⛔ **Arm 3 — this is D-6, and it FAILED when it was written**.

    Locus 3's gap is crossed contiguously by ``DEEP`` fully-sequenced fragments and spliced by only
    ``SHALLOW``, so the evidence says the molecule is genomic. ⭐ Its intron spans **exactly one region** —
    the locus has one isoform, so the intron's endpoints are adjacent region_bounds — which is precisely the
    configuration where the shipped ``_boundaries_inside`` returned an EMPTY evidence set and handed ∅ a
    structural ``rho = 0``.

    ⚠ The deposit rule is what settles the right set, not taste: a boundary is crossed iff it is strictly
    inside a *segment*, so the boundaries distinguishing ∅ from a path splicing ``[a, b)`` are those at region_bounds
    ``a <= c <= b`` — **endpoints included**, and both endpoints are always region_bounds.
    """
    shares = _shares(scored[0], scored[1], L3)
    genomic, spliced = (), ((7600, 8000),)
    assert shares[genomic] > shares[spliced], (
        f"the gap is crossed contiguously {DEEP} times and spliced only {SHALLOW} times, so ∅ must take "
        f"the larger share; got genomic={shares[genomic]:.4f} spliced={shares[spliced]:.4f}. A "
        f"genomic share of exactly 0 means ∅'s contiguous-boundary evidence set was empty — D-6."
    )


@pytest.mark.parametrize("base,crossed_end", [(L4, "donor"), (L6, "acceptor")])
def test_the_GENOMIC_hypothesis_needs_evidence_at_BOTH_ENDS_of_the_gap(scored, base, crossed_end):
    """⛔ **Arm 4 — the half of D-6 that arms 1–3 could not see.**

    A perturbation keeping only the **donor** boundary of the distinguishing set passed every other gate in
    this module, because at locus 3 both ends carry the same deep coverage and either one alone answers.
    Loci 4 and 6 separate them, and they are a **mirror pair**: locus 4 crosses only the donor boundary and
    locus 6 only the acceptor, so a rule that consults either end alone fails one of the two.

    ⭐ A molecule that crosses the gap contiguously is present at *both* of its ends, so the scarcest
    object on the path bounds it — which is what ``min`` aggregation means (D-1's bottleneck reading).
    With no fragment ever seen crossing the acceptor boundary, there is no evidence for a contiguous crossing
    however deep the donor side is, and ∅ must lose to the sj that does have flux.

    ⚠ ∅'s zero here is a **correct** zero, and the owner's D-3 ruling is what makes it stand: the score
    keeps no fallback, so a hypothesis with no evidence is unselectable rather than floored.
    """
    shares = _shares(scored[0], scored[1], base)
    genomic, spliced = (), ((base + 600, base + 1000),)
    empty_end = base + 1000 if crossed_end == "donor" else base + 600
    assert shares[spliced] > shares[genomic], (
        f"nothing was ever observed crossing the {'acceptor' if crossed_end == 'donor' else 'donor'} "
        f"boundary at {empty_end}, so ∅ has no evidence for a contiguous crossing; got "
        f"genomic={shares[genomic]:.4f} spliced={shares[spliced]:.4f}. A winning ∅ means only ONE end "
        f"of the gap was consulted."
    )


def test_when_rho_TIES_the_LENGTH_term_decides(scored):
    """⛔ **Arm 5.** Dropping ``f`` from the score entirely passed arms 1–4, because rho is deliberately
    decisive there. This locus removes rho from the contest.

    Both sj are observed by the same number of fragments over the same block layouts, so their
    deposited densities are **identical** — checked below, not assumed. What differs is the implied
    length: 300 bp for the wide intron and 500 bp for the narrow one, and the ballast puts three times
    the mass at 500. ⭐ Neither bin is reachable by any locus's depth fragments, so that 3:1 is structural
    rather than a count that happens to come out right. The narrow hypothesis must win, through ``f``.
    """
    payload, result = scored[0], scored[1]
    shares = _shares(payload, result, L5)
    wide, narrow = ((L5 + 650, L5 + 950),), ((L5 + 700, L5 + 800),)
    d = payload.deferred
    i = next(j for j in range(d.n_fragments) if int(d.start[j]) == L5 + 500)
    rho = {
        tuple(tuple(p) for p in d.hypothesis_introns_of(h).tolist()): float(result.terms.density[h])
        for h in range(int(d.hypothesis_offsets[i]), int(d.hypothesis_offsets[i + 1]))
    }
    assert rho[wide] == rho[narrow], (
        f"the two sj were built with identical depth and identical block layouts, so their "
        f"densities must be equal; got {rho[wide]!r} and {rho[narrow]!r}. Without that this gate is not "
        f"isolating the length term."
    )
    assert shares[narrow] > shares[wide], (
        f"rho ties, so the anchor's greater mass at L={L_ARM5_NARROW} than at L={L_ARM5_WIDE} must "
        f"decide it; got narrow={shares[narrow]:.4f} wide={shares[wide]:.4f}"
    )


# ── non-vacuity ────────────────────────────────────────────────────────────────────────────────────


def test_the_fixture_actually_DECIDES_something(scored):
    """⛔ The trap this module exists to avoid. records that the P2 gate was **not** written
    against the smoke fixture because every score there came out uniform — a green gate over a scorer that
    decided nothing. Uniformity is therefore a failure condition here, checked directly."""
    payload, result = scored[0], scored[1]
    assert payload.deferred.n_fragments == 7, (
        f"the seven ambiguous fragments, and only those, must be held; got "
        f"{payload.deferred.n_fragments}"
    )
    # ⭐ Arm 7's locus is EXPECTED to tie here, and that is its whole premise: its two candidates differ
    # only in strand, and this fixture carries no strand information, so `p_r1_sense` is the neutral 0.5 and
    # `s` is 0.5 either way. ⚠ `test_the_STRAND_term_decides_when_rho_and_LENGTH_both_tie` is what breaks
    # that tie, by rescoring the same payload at 0.99 and at 0.01.
    assert result.n_undecided == 1, (
        f"{result.n_undecided} records had two or more candidates tied for the lead. Exactly one is "
        f"expected — arm 7's, which ties by construction at a neutral strand fraction. More than that "
        f"means a gate above is comparing equal numbers."
    )
    for base in (L1, L2, L3, L4, L5, L6):
        shares = _shares(payload, result, base)
        uniform = 1.0 / len(shares)
        assert max(shares.values()) > uniform, (
            f"the record at {base + 500} is uniform ({shares}); the evidence separated nothing"
        )
    assert payload.qc.dropped_too_long == 0 and payload.qc.dropped_empty == 0


def test_the_length_term_has_support_at_EVERY_hypothesis_length(scored):
    """⚠ ``build_fl_models`` does not smooth the global anchor, so ``global_pmf[L]`` is exactly zero
    unless a deposited fragment had that length — and a zero there would kill ∅ for a reason that has
    nothing to do with the density. That is what the ballast is for, and this checks it worked."""
    payload, result = scored[0], scored[1]
    for length in (L_WIDE, L_ARM5_WIDE, L_NARROW, L_ARM5_NARROW, L_GENOMIC):
        assert payload.deposited_lengths[length] > 0, (
            f"no deposited fragment has L = {length}, so the anchor has no mass there and every "
            f"hypothesis at that length scores 0 through f alone. Anchor support: "
            f"{np.nonzero(payload.deposited_lengths)[0].tolist()}"
        )
    assert np.all(result.terms.length_likelihood > 0.0), (
        "every hypothesis in this fixture must have positive length likelihood, or the gates above are "
        "measuring the length term's support rather than the density"
    )
