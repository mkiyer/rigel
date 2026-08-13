"""⭐ STRAND, in the gap-hypothesis enumeration — what pins it, what leaves it open, what contradicts it.

     and D-5
    Owner ruling, 2026-08-02, recorded verbatim below.

> *"Splice junctions are stranded and asymmetric. And so if we detect one splice junction in a fragment,
> we know this fragment strand. So any fragment with a splice junction, immediately, we can constrain the
> gap [hypotheses] to that strand. Fragments without a splice junction are unspliced and could be either
> strand, but the library strand specificity could constrain this."*

Two behaviours follow, and **neither was tested anywhere** before this module. They are not new — the
audit found the first already working — but an untested behaviour is one a refactor is free to delete.

⛔ **THESE GATES MUST GO THROUGH THE SCAN.** Driving ``FragmentResolver`` directly leaves
``t_strand_arr_`` empty, so every hypothesis's implied strand silently reads ``NONE`` and a strand gate
written that way passes against nothing. Measured during the audit; it cost an hour. Everything here is
read off ``AccumulatorPayload`` after a real scan of a real BAM.

⚠ **The index warns about the fixture, and that is expected.** Two of the three gene pairs annotate the
same ``(donor, acceptor)`` on both strands, which the index correctly calls biologically impossible — a
``GT..AG`` intron reverse-complements to ``CT..AC``. That is exactly the configuration these gates need,
and the warning is the index doing its job.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.types import Strand

GENOME = 8_000

#: ⭐ Three pairs, each isolating one behaviour.
#:
#: **tP/tM** share the observed junction ``(1200,1400)`` on OPPOSITE strands and cover the same sequenced
#: blocks, so neither annotation nor overlap can separate them — only the sequenced **motif** can. Their
#: gap introns differ in width (400 vs 200 bp), so the deposited ``L`` says which one was believed.
#:
#: **tQ/tR** imply the SAME gap intron on opposite strands: the D-5 case, where the two paths are
#: genuinely one path but their strands disagree.
GTF = (
    # tP (+): introns (1200,1400) observed, (1600,2000) in the gap  -> cuts 400
    'chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    'chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    'chr1\ttest\texon\t2001\t2200\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    # tM (-): same observed junction, but its gap intron is (1700,1900) -> cuts 200
    'chr1\ttest\texon\t1001\t1200\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    'chr1\ttest\texon\t1401\t1700\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    'chr1\ttest\texon\t1901\t2200\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    # tQ (+) and tR (-): the SAME gap intron (3200,3600) on both strands — D-5
    'chr1\ttest\texon\t3001\t3200\t.\t+\t.\tgene_id "gQ"; transcript_id "tQ";\n'
    'chr1\ttest\texon\t3601\t3800\t.\t+\t.\tgene_id "gQ"; transcript_id "tQ";\n'
    'chr1\ttest\texon\t3001\t3200\t.\t-\t.\tgene_id "gR"; transcript_id "tR";\n'
    'chr1\ttest\texon\t3601\t3800\t.\t-\t.\tgene_id "gR"; transcript_id "tR";\n'
)

#: The four molecules, and what each one is FOR. Extents are the merged block span.
#:
#:   name      reads                                   extent        the question
#:   pinned_p  100M 200N 150M @1100  +  100M @2050     [1100,2150)   does a + motif pin it to tP?
#:   pinned_m  the same reads with a - motif           [1100,2150)   ... and a - motif to tM?
#:   open      100M @1450            +  100M @2050     [1450,2150)   no motif -> BOTH strands offered
#:   same_path 100M @3050            +  100M @3650     [3050,3750)   D-5: one path, two strands
L_PINNED_P = (2150 - 1100) - 200 - 400  # 450 — tP's wide gap intron was cut
L_PINNED_M = (2150 - 1100) - 200 - 200  # 650 — tM's narrow one was cut


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    """Scan the fixture BAM once and hand back the payload."""
    import pysam

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("gap_strand")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(">chr1\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n")
    pysam.faidx(str(fasta))
    gtf.write_text(GTF)
    index_dir = base / "idx"
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        TranscriptIndex.build(str(fasta), str(gtf), str(index_dir), write_tsv=False)
    with pytest.warns(RuntimeWarning, match="strand-coincident"):
        index = TranscriptIndex.load(str(index_dir))

    header = {"HD": {"VN": "1.6", "SO": "queryname"}, "SQ": [{"SN": "chr1", "LN": GENOME}]}
    M, N = 0, 3

    def read(qname, pos, cigar, mate_pos, is_r1, motif=None):
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
        tags = [("NH", 1, "i")]
        # ⭐ XS is the aligner's GENOMIC MOTIF strand. It is what pins the fragment: a GT..AG intron is
        # on +, its reverse complement CT..AC on -, and no annotation is consulted to decide it.
        if motif is not None:
            tags.append(("XS", motif, "A"))
        a.set_tags(tags)
        return a

    spliced = [(M, 100), (N, 200), (M, 150)]
    reads = [
        read("pinned_p", 1100, spliced, 2050, True, motif="+"),
        read("pinned_p", 2050, [(M, 100)], 1100, False),
        read("pinned_m", 1100, spliced, 2050, True, motif="-"),
        read("pinned_m", 2050, [(M, 100)], 1100, False),
        read("open", 1450, [(M, 100)], 2050, True),
        read("open", 2050, [(M, 100)], 1450, False),
        read("same_path", 3050, [(M, 100)], 3650, True),
        read("same_path", 3650, [(M, 100)], 3050, False),
    ]
    bam_path = str(base / "gap_strand.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam_path, bam_path)

    _stats, _strand, _buffer, payload = scan_and_buffer(
        bam_path, index, BamScanConfig(sj_strand_tag="XS", total_threads=1)
    )
    return payload


def _record(payload, start: int):
    """The one held record whose clipped extent begins at ``start``."""
    matches = [
        i for i in range(payload.deferred.n_fragments) if int(payload.deferred.start[i]) == start
    ]
    assert len(matches) == 1, f"expected exactly one held record starting at {start}, got {matches}"
    return matches[0]


def _hypotheses(payload, i: int):
    """Record ``i``'s hypotheses as ``(implied introns, implied strand)``."""
    d = payload.deferred
    return [
        ([tuple(p) for p in d.hypothesis_introns_of(h).tolist()], int(d.hypothesis_sj_strand[h]))
        for h in range(int(d.hypothesis_offsets[i]), int(d.hypothesis_offsets[i + 1]))
    ]


def _lengths(payload) -> dict[int, int]:
    nz = np.nonzero(payload.deposited_lengths)[0]
    return {int(i): int(payload.deposited_lengths[i]) for i in nz}


# ── an OBSERVED junction pins the strand ───────────────────────────────────────────────────────────


def test_an_OBSERVED_junction_PINS_the_gap_hypotheses_to_its_own_strand(scanned):
    """⭐ The owner's ruling, and the audit found it already true — now it is gated.

    ``tP`` and ``tM`` cover the same sequenced blocks and share the observed junction's coordinates, so
    nothing separates them but the **motif**. A ``+`` motif must leave only tP's 400 bp gap intron on the
    table and a ``-`` motif only tM's 200 bp one — and the deposited ``L`` says which was believed,
    without needing to read the hypothesis set at all.

    ⛔ Without the pin both transcripts would contribute, the fragment would have two surviving
    hypotheses, and it would be **deferred instead of deposited** — so this gate fires on the count as
    well as on the lengths.
    """
    lengths = _lengths(scanned)
    assert lengths.get(L_PINNED_P) == 1, (
        f"expected a fragment at L={L_PINNED_P} (tP's 400 bp gap intron cut, pinned by the + motif); "
        f"deposited lengths are {lengths}"
    )
    assert lengths.get(L_PINNED_M) == 1, (
        f"expected a fragment at L={L_PINNED_M} (tM's 200 bp gap intron cut, pinned by the - motif); "
        f"deposited lengths are {lengths}"
    )
    # Both pinned fragments RESOLVED; only `open` and `same_path` are held.
    assert scanned.gap_resolution.gap_resolved_spliced == 2
    assert scanned.deferred.n_fragments == 2


# ── an UNSPLICED fragment leaves both strands open ─────────────────────────────────────────────────


def test_an_UNSPLICED_fragment_offers_BOTH_STRANDS(scanned):
    """⭐ The other half of the ruling: nothing sequenced pins the strand, so nothing may narrow it.

    ⚠ **This is the case that makes the second pass's strand term necessary** rather than decorative.
    The two spliced hypotheses here differ in strand *and* in implied length; drop the strand term and
    they are separated by length alone.
    """
    hypotheses = _hypotheses(scanned, _record(scanned, 1450))
    introns = {tuple(path[0]) for path, _strand in hypotheses if path}
    assert introns == {(1600, 2000), (1700, 1900)}, (
        f"both strands' gap introns must be offered; got {introns}"
    )
    strands = {strand for path, strand in hypotheses if path}
    assert strands == {int(Strand.POS), int(Strand.NEG)}, (
        f"the two spliced hypotheses must carry OPPOSITE implied strands; got {strands}. A single "
        f"strand here means the enumeration narrowed on something that was never sequenced."
    )
    assert any(not path for path, _strand in hypotheses), (
        "the genomic hypothesis must also be present — no annotated junction was observed, so the gap "
        "may be real template"
    )


# ── D-5: one path, two strands ─────────────────────────────────────────────────────────────────────


def test_ONE_PATH_claimed_by_BOTH_STRANDS_is_marked_AMBIGUOUS(scanned):
    """⛔ **D-5.** ``tQ`` (+) and ``tR`` (−) imply the *same* intron ``(3200,3600)``.

    Grouping by path is right — it **is** one path, and one hypothesis is the correct count. ⛔ But the
    hypothesis carries ONE ``sj_strand``, and taking the first supporter's silently asserts a strand the
    evidence does not support: swap the two GTF boundaries and the answer flips.

    ⭐ ``AMBIGUOUS`` is what that state is called everywhere else in this codebase — the fragment-level
    ``sj_strand`` uses it for exactly this, "contradictory evidence rather than missing evidence", and
    ``deposit`` already refuses to credit a junction on it. Reusing the value keeps one vocabulary.

    ⚠ Unreachable on human data — **0 of 404,168** junction coordinates are annotated on both strands,
    and the index warns that it is biologically impossible. Fixed anyway because it is three boundaries and
    the alternative is an answer that depends on GTF boundary order.
    """
    hypotheses = _hypotheses(scanned, _record(scanned, 3050))
    spliced = [(path, strand) for path, strand in hypotheses if path]
    assert len(spliced) == 1, (
        f"tQ and tR imply the same intron, so they are ONE path and must group into one hypothesis; "
        f"got {spliced}"
    )
    (path, strand) = spliced[0]
    assert path == [(3200, 3600)]
    assert strand == int(Strand.AMBIGUOUS), (
        f"the merged hypothesis reports strand {strand}, but its two supporters disagree. Taking the "
        f"first supporter's strand makes the answer depend on GTF boundary order."
    )
    assert (
        len(
            scanned.deferred.supporting_t_of(
                int(scanned.deferred.hypothesis_offsets[_record(scanned, 3050)])
            )
        )
        == 2
    ), "both transcripts must remain recorded as supporters of the merged path"


def test_the_fixture_reaches_every_branch_it_claims_to(scanned):
    """⚠ Non-vacuity. Each gate above reads one record or one length bin; if the scan quietly produced
    something else, several of them could pass for the wrong reason."""
    assert scanned.qc.deposited == 2, "the two pinned fragments, and only those, must deposit"
    assert scanned.deferred.n_fragments == 2, (
        "the unspliced and the same-path fragments must be held"
    )
    assert scanned.gap_resolution.deferred == 2
    assert scanned.qc.dropped_too_long == 0 and scanned.qc.dropped_empty == 0


# ⛔ A PERTURBATION THAT DOES NOT FIRE, AND WHY — recorded rather than left as a hole.
#
# Z3 removed `!certified_rna` from `enumerate_gap_hypotheses`'s unspliced-hypothesis condition and
# **nothing here failed**. The cause is not a weak fixture: the nRNA SHADOWS are in the candidate set.
# Measured on the `open` fragment — `t_inds = [0, 1, 4, 5]`, where 4 and 5 are
# `RIGEL_NRNA_chr1_{1,2}_1000_2200`, both `is_synthetic`. A shadow is single-exon, so it implies nothing
# in the gap, so `any_candidate_implies_nothing` is already true and ∅ is emitted without the clause.
#
# ⭐ That is 's own ruling arriving by a different route — *"the nascent shadow IS
# the ∅ hypothesis"* — so the two mechanisms agree rather than conflict:
#
#   * SPLICED fragment: an observed CIGAR-N intron falls inside the shadow's single exon, so the shadow
#     cannot explain the read and drops out of `t_inds`. ∅ is then correctly absent — measured in
#     `test_gap_introns_are_cut.py`, where the certified-RNA `mixed` fragment has exactly ONE hypothesis.
#   * UNSPLICED fragment: the shadow survives, implies nothing, and supplies ∅ — which is what
#     §2's table requires ("no annotated junction ⇒ ∅ always").
#
# ⚠ So `!certified_rna` is REDUNDANT here, not wrong, and it is kept: it states the rule directly instead
# of depending on the shadow mechanism continuing to exist. A single-exon gene has no separate shadow row
# at all (`is_nrna` on a non-synthetic row means mature ≡ nascent), and that is the case the clause covers
# on its own. ⛔ Do not delete it on the strength of this perturbation.
