"""⭐ `L` excludes an annotated intron in the mate gap, whatever the fragment's splice type.

    Spec: ``docs/SPEC_GAP_INTRONS.md``   ·   Cause and evidence: ``docs/JUNCTION_OPPORTUNITY.md`` §4

`docs/FRAGMENT_LENGTH_AUDIT.md`'s C2 left the tool with ONE definition of fragment length — the
accumulator's ``L``, the total length of the fragment's own path. C0 proved that definition correct
*given its inputs*. This module is about an input that was **incomplete**: implicit-splice detection ran
only on fragments the resolver had already called ``SPLICE_UNSPLICED``, so a fragment carrying an
observed CIGAR-N splice never had its unsequenced mate gap examined and kept that intron inside ``L``.

⚠ **Why a hand-written BAM and not the simulator.** Every number below is a fragment length asserted to
the base pair, and the point of the gate is that ``L`` is *exactly* the molecule. A simulated library
gives a distribution; three reads written by hand give an answer. ``deposited_lengths`` is C1's
unconditional histogram — one bin per deposited fragment, indexed by ``L`` — so a four-fragment BAM makes
``L`` directly readable with no model, no fit and no tolerance.

⛔ **The trap this module also pins** (spec §1): the gap finder walks consecutive aligned blocks and
emits every hole, so a CIGAR-N intron is a "hole" too. Dropping only gaps that EXACTLY match an observed
intron is what stops the detector substituting a *different* annotated intron that happens to lie within
the ±K anchor tolerance — which would merge into one wider interval and make ``L`` too SHORT.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import POOL_RNA_SPLICED
from rigel.splice import SpliceType, census_field

GENOME = 70_000

#: g6/t_mixed — three exons, so a fragment can carry an observed splice AND a gap intron.
#: g7/t_two_{a,b} — two isoforms whose gap introns DISAGREE, which is the deferral case.
#: (GTF is 1-based inclusive; the 0-based half-open geometry each line produces is in the comment.)
GTF = (
    # t_mixed exons (62000,62200) (62400,62600) (62800,63000); introns (62200,62400) (62600,62800)
    'chr1\ttest\texon\t62001\t62200\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62401\t62600\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62801\t63000\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    # t_two_a exons (66000,66200) (66400,66600) (66800,67000); gap intron (66600,66800)
    'chr1\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr1\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr1\ttest\texon\t66801\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    # t_two_b — same first two exons, but its gap intron is (66600,66700): a DIFFERENT L
    'chr1\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr1\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr1\ttest\texon\t66701\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
)

# ── the four fragments, and the L each one must produce ────────────────────────────────────────────
#
#   name      read1 CIGAR at its start          read2       extent      cut                     L
#   pure      62100  100M 200N 100M             62500 100M  [62100,62600)  (62200,62400)        300
#   mixed     62100  100M 200N 100M             62900 100M  [62100,63000)  + (62600,62800)      500
#   near      62100  102M 198N 100M             62900 100M  [62100,63000)  (62202,62400)
#                                                                          + (62600,62800)      502
#   ambig     66100  100M 200N 100M             66900 100M  — the two isoforms disagree —    REJECTED
#
#: ``pure`` is the control that must NOT move: an observed splice and no mate gap at all.
L_PURE = 300
#: ``mixed`` is the population the spec exists to find: observed splice PLUS a gap intron.
L_MIXED = 500
#: ``near`` is the §1 trap: its observed intron is 2 bp inside the annotated donor, so the ANNOTATED
#: intron lies within K=3 of the observed gap. Substituting it would cut 200 bp instead of 198 → 500.
L_NEAR = 502


@pytest.fixture(scope="module")
def gap_intron_index(tmp_path_factory):
    import pysam

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("gap_introns")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(">chr1\n" + "\n".join(["N" * 80] * (GENOME // 80)) + "\n")
    pysam.faidx(str(fasta))
    gtf.write_text(GTF)
    index_dir = base / "idx"
    TranscriptIndex.build(str(fasta), str(gtf), str(index_dir), write_tsv=False)
    return base, TranscriptIndex.load(str(index_dir))


@pytest.fixture(scope="module")
def gap_intron_bam(gap_intron_index):
    """One BAM holding the four fragments above, name-sorted and NH-tagged."""
    import pysam

    base, _ = gap_intron_index
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": GENOME}],
    }

    def _read(qname, pos, cigar, mate_pos, is_r1):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = 60
        # Proper FR pair: R1 forward with a reverse mate, R2 reverse. `build_fragment` flips R2's
        # strand, so both mates land in ONE (ref, strand) group and merge into a single extent.
        a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10)
        a.cigar = cigar
        n_query = sum(length for op, length in cigar if op == 0)
        a.query_sequence = "A" * n_query
        a.query_qualities = pysam.qualitystring_to_array("I" * n_query)
        a.next_reference_id = 0
        a.next_reference_start = mate_pos
        a.set_tags([("NH", 1, "i")])
        return a

    M, N = 0, 3
    reads = [
        # pure — observed splice, mates ABUT so there is no mate gap to search
        _read("pure", 62100, [(M, 100), (N, 200), (M, 100)], 62500, True),
        _read("pure", 62500, [(M, 100)], 62100, False),
        # mixed — observed splice AND an annotated intron inside the [62500,62900) mate gap
        _read("mixed", 62100, [(M, 100), (N, 200), (M, 100)], 62900, True),
        _read("mixed", 62900, [(M, 100)], 62100, False),
        # near — the observed splice lands 2 bp inside the annotated donor
        _read("near", 62100, [(M, 102), (N, 198), (M, 100)], 62900, True),
        _read("near", 62900, [(M, 100)], 62100, False),
        # ambig — two isoforms imply different introns in the [66500,66900) gap
        _read("ambig", 66100, [(M, 100), (N, 200), (M, 100)], 66900, True),
        _read("ambig", 66900, [(M, 100)], 66100, False),
    ]
    bam_path = str(base / "gap_introns.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam_path, bam_path)
    return bam_path


@pytest.fixture(scope="module")
def payload(gap_intron_index, gap_intron_bam):
    _, index = gap_intron_index
    _, _, _, payload = scan_and_buffer(gap_intron_bam, index, BamScanConfig(sj_strand_tag="auto"))
    return payload


def _lengths(payload) -> dict[int, int]:
    """``deposited_lengths`` as ``{L: count}`` over its non-zero bins."""
    nz = np.nonzero(payload.deposited_lengths)[0]
    return {int(i): int(payload.deposited_lengths[i]) for i in nz}


# ── U1 · the headline ──────────────────────────────────────────────────────────────────────────────


def test_U1_L_excludes_BOTH_the_observed_intron_and_the_one_in_the_mate_gap(payload):
    """⭐ The spec's U1, read off the tally rather than argued about.

    ``mixed`` spans [62100,63000) — 900 bp of genome — and is a 500 bp molecule: 200 bp of it was
    spliced out and sequenced as CIGAR-N, another 200 bp was spliced out and never sequenced at all.
    Before this work only the first cut was made and the fragment measured 700 bp.
    """
    assert _lengths(payload) == {L_PURE: 1, L_MIXED: 1, L_NEAR: 1}, (
        "the deposited fragment lengths are not the four molecules this BAM contains. 700 bp means the "
        "mate-gap intron was never cut (the SPLICE_UNSPLICED gate); 500 in place of 502 means a NEAR "
        "match was substituted for the observed gap — SPEC_GAP_INTRONS.md §0 and §1"
    )
    assert payload.qc.deposited == 3


def test_U1b_the_unspliced_control_is_untouched(payload):
    """``pure`` has an observed splice and NO mate gap: nothing about it may move."""
    assert _lengths(payload)[L_PURE] == 1


def test_the_two_intron_lists_are_DISJOINT_so_the_union_absorbs_nothing(payload):
    """⭐ The union of observed and implied introns is a union, not a merge.

    ``introns_absorbed`` counts introns the accumulator had to normalise away because they overlapped or
    abutted. The §1 filter's claim is that the implied list can only hold holes the CIGAR did not
    explain, so unioning the two lists can never produce an overlap — and this counter is where that
    claim shows up. ⚠ It is also what the adapter's sort protects: the de-duplication before the deposit
    is an adjacent-pair comparison, so an unsorted union would pass a duplicate through to be absorbed
    here instead of skipped there.
    """
    assert payload.qc.introns_absorbed == 0


# ── U3 · the near-match trap ───────────────────────────────────────────────────────────────────────


def test_U3_a_near_match_does_not_shorten_L(payload):
    """⛔ ``near``'s observed intron is [62202,62400); the ANNOTATED one is [62200,62400).

    They are 2 bp apart and K is 3, so a filter written on *overlap* rather than exact equality lets the
    annotated intron be emitted for a gap the CIGAR already explained. It then normalises together with
    the observed one into [62200,62400) and cuts 200 bp where the molecule lost 198 — ``L`` too SHORT.
    """
    assert L_NEAR in _lengths(payload), (
        f"expected the near-match fragment at L={L_NEAR}; L={L_MIXED} means the annotated intron was "
        f"substituted for the observed one"
    )


# ── U5 · the deposit is refused when the candidates disagree ───────────────────────────────────────


def test_U5_a_mixed_fragment_with_disagreeing_candidates_is_REJECTED_and_COUNTED(payload):
    """⛔ ``ambig`` carries an observed splice, so it never reached this test before.

    ``t_two_a`` puts [66600,66800) in the mate gap and ``t_two_b`` puts [66600,66700) there — a 100 bp
    difference in ``L`` for one molecule. There is no answer to deposit, and picking either is picking a
    fragment length at random. ⚠ The rejection must be COUNTED: `path_ambiguous` sizes the population a
    second pass can recover, and a silent return would erase it.
    """
    assert payload.qc.dropped_ambiguous_path == 1, (
        "the two isoforms imply different introns in the same unsequenced gap, so L is undetermined "
        "even though the fragment's OTHER splice was observed — SPEC_GAP_INTRONS.md §2.3"
    )


# ── D1 · a mixed fragment leaves the pure RNA length pool ──────────────────────────────────────────


def test_D1_only_the_fragment_with_a_FULLY_OBSERVED_path_stays_in_the_RNA_pool(payload):
    """⭐ Owner decision D1. The pool is a **length** histogram used to FIT the fragment-length model,
    and a length partly inferred from the annotation is a product of the model it would help fit.

    ⚠ Two-sided on purpose. ``pure``'s splice was sequenced end to end and it must stay; ``mixed`` and
    ``near`` each owe 200 bp of their ``L`` to an intron nobody read, and they must leave. A one-sided
    version of this test passes on an accumulator that simply stopped filling the pool.
    """
    rna = payload.pool_lengths[POOL_RNA_SPLICED]
    assert {int(i): int(rna[i]) for i in np.nonzero(rna)[0]} == {L_PURE: 1}, (
        "the pure-RNA length pool must hold exactly the fragment whose whole path was observed"
    )
    assert payload.qc.sj_implicit_fragments == 2, (
        "`sj_implicit` means 'this L depends on an intron that was never sequenced', which is now true "
        "of a fragment that ALSO has an observed splice"
    )


# ── U4 · the classification does not move ──────────────────────────────────────────────────────────


def test_U4_the_splice_census_is_unchanged_by_gap_intron_detection(
    gap_intron_index, gap_intron_bam
):
    """Detection is unconditional; the ``SPLICE_IMPLICIT`` **promotion** stays unspliced-only.

    ``splice_type`` feeds scoring, the buffer, the strand training and — since C2 — ``rigel report``'s
    census. All four fragments here carry an observed CIGAR-N splice, so all four must be counted as
    observed splices and **none** as implicit.
    """
    _, index = gap_intron_index
    stats, _, _, _ = scan_and_buffer(gap_intron_bam, index, BamScanConfig(sj_strand_tag="auto"))

    def census(stype):
        return getattr(stats, census_field(stype))

    assert census(SpliceType.SPLICED_ANNOT) == 3, "pure, mixed and ambig are annotated splices"
    assert census(SpliceType.SPLICED_UNANNOT) == 1, (
        "near's observed intron is 2 bp off the annotation"
    )
    assert census(SpliceType.SPLICED_IMPLICIT) == 0, (
        "no fragment here may be RE-LABELLED implicit: this work is about L, not classification"
    )
    assert census(SpliceType.UNSPLICED) == 0
