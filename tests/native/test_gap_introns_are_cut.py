"""⭐ An annotated intron in the mate gap is found on EVERY fragment, and never a near-miss for it.

     (which supersedes)
    Cause and evidence:

's TRAPS: pure-and-length-censored left the tool with ONE definition of fragment length — the
accumulator's ``L``, the total length of the fragment's own path. TRAPS: two-divisors-opposite-sign proved that definition correct
*given its inputs*. This module is about an input that was **incomplete**: implicit-splice detection ran
only on fragments the resolver had already called ``SPLICE_UNSPLICED``, so a fragment carrying an
observed CIGAR-N splice never had its unsequenced mate gap examined and kept that intron inside ``L``.

⚠ **Why a hand-written BAM and not the simulator.** Every number below is a fragment length asserted to
the base pair, and the point of the gate is that ``L`` is *exactly* the molecule. A simulated library
gives a distribution; four reads written by hand give an answer. ``deposited_lengths`` is TRAPS: a-purity-filter-is-a-length-filter's
unconditional histogram — one bin per deposited fragment, indexed by ``L`` — so a four-fragment BAM makes
``L`` directly readable with no model, no fit and no tolerance.

⛔ **The trap this module also pins** (spec §1): the gap finder walks consecutive aligned blocks and
emits every hole, so a CIGAR-N intron is a "hole" too. Dropping only gaps that EXACTLY match an observed
intron is what stops the detector substituting a *different* annotated intron that happens to lie within
the ±K anchor tolerance — which would merge into one wider interval and make ``L`` too SHORT.

⭐ **WHAT MOVED AT S1, AND WHY IT IS THE SAME GATE.** made the ACCUMULATOR the
arbiter: a fragment arrives with its hypothesis SET, and if more than one hypothesis survives its ``L`` is
undetermined and it is **held whole in the side buffer** rather than deposited. Two of the four fragments
here now take that route, so their ``L`` is no longer readable from ``deposited_lengths`` — it is readable
from the hypotheses in the bank, which is a **stronger** statement: the old gate said the deposited number
was 502, this one says the enumeration produced exactly the intron set that yields 502 and never the
near-miss that yields 500. ⛔ Nothing about the near-match trap is relaxed; it is asserted on coordinates
instead of on a length.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import POOL_RNA_SPLICED
from rigel.splice import SpliceType, census_field

GENOME = 70_000

#: ⭐ **TWO CONTIGS, and the second one is load-bearing.** ``g7`` — the disagreeing-isoform gene, the only
#: one that produces a deferral by structure — lives on **chr2**, so the side buffer must hold a record
#: stamped ``ref = 1``.
#:
#: ⚠ Measured, not supposed: with everything on one contig, a perturbation giving **every** reference the id
#: ``0`` passed the entire suite, 1860 tests. A single-contig fixture cannot tell a correct reference stamp
#: from a constant, and the second pass replays each held record onto that reference's cut axis — so a
#: constant stamp would drain chr2's coordinates onto chr1's partition. That is the defect the predecessor
#: adapter actually had (`multi_reference_bam`'s docstring in ``test_scanner_accumulator_integration.py``).
#:
#: g6/t_mixed — three exons, so a fragment can carry an observed splice AND a gap intron.
#: g7/t_two_{a,b} — two isoforms whose gap introns DISAGREE, which is the deferral case.
#: (GTF is 1-based inclusive; the 0-based half-open geometry each line produces is in the comment.)
GTF = (
    # t_mixed exons (62000,62200) (62400,62600) (62800,63000); introns (62200,62400) (62600,62800)
    'chr1\ttest\texon\t62001\t62200\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62401\t62600\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62801\t63000\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    # t_two_a exons (66000,66200) (66400,66600) (66800,67000); gap intron (66600,66800)  — ON CHR2
    'chr2\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr2\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr2\ttest\texon\t66801\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    # t_two_b — same first two exons, but its gap intron is (66600,66700): a DIFFERENT L
    'chr2\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr2\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr2\ttest\texon\t66701\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
)

#: Which reference each gene sits on, so the gates can say which cut axis a record must replay onto.
REF_CHR1, REF_CHR2 = 0, 1

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

#: ⭐ ``near`` and ``ambig`` are DEFERRED, not deposited — for two different reasons, and both are the
#: arbitration rule working rather than a loss:
#:
#: * ``near``'s observed CIGAR-N junction is 2 bp off the annotation, so the fragment is
#:   ``SPLICED_UNANNOT`` and **not certified RNA**. The unspliced (genomic) hypothesis therefore survives
#:   alongside the spliced one — the molecule may be gDNA or nascent with the gap as real template — and
# two survivors mean ``L`` is undetermined: ∅ is available whenever no
#:   *annotated* junction was sequenced.
#: * ``ambig``'s two isoforms imply DIFFERENT introns in the same gap. gDNA cannot be spliced and its
#:   junction IS annotated, so the molecule is certified RNA and only the structure is open.
_DEFERRED = 2


def _covered_length(start: int, end: int, introns) -> int:
    """``L`` by integer SET ARITHMETIC — a different algorithm from the accumulator's segment walk.

    ⚠ Deliberately not the reference's own ``_hypothesis_length``: a
    validator that calls the implementation's helper validates nothing, and the whole point here is that a
    hypothesis's intron set yields 502 rather than 500.
    """
    covered = set(range(start, end))
    for a, b in introns:
        covered -= set(range(int(a), int(b)))
    return len(covered)


def _hypotheses_of(payload, i: int) -> list[tuple[list[tuple[int, int]], int]]:
    """Deferred record ``i``'s hypotheses as ``(implied introns, implied L)``, observed introns included."""
    deferred = payload.deferred
    observed = [tuple(pair) for pair in deferred.observed_introns_of(i).tolist()]
    out = []
    for h in range(int(deferred.hypothesis_offsets[i]), int(deferred.hypothesis_offsets[i + 1])):
        implied = [tuple(pair) for pair in deferred.hypothesis_introns_of(h).tolist()]
        out.append(
            (
                implied,
                _covered_length(int(deferred.start[i]), int(deferred.end[i]), observed + implied),
            )
        )
    return out


@pytest.fixture(scope="module")
def gap_intron_index(tmp_path_factory):
    import pysam

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("gap_introns")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    contig = "\n".join(["N" * 80] * (GENOME // 80))
    fasta.write_text(f">chr1\n{contig}\n>chr2\n{contig}\n")
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
        "SQ": [{"SN": "chr1", "LN": GENOME}, {"SN": "chr2", "LN": GENOME}],
    }

    def _read(qname, pos, cigar, mate_pos, is_r1, ref_id=REF_CHR1):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = ref_id
        a.reference_start = pos
        a.mapping_quality = 60
        # Proper FR pair: R1 forward with a reverse mate, R2 reverse. `build_fragment` flips R2's
        # strand, so both mates land in ONE (ref, strand) group and merge into a single extent.
        a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10)
        a.cigar = cigar
        n_query = sum(length for op, length in cigar if op == 0)
        a.query_sequence = "A" * n_query
        a.query_qualities = pysam.qualitystring_to_array("I" * n_query)
        a.next_reference_id = ref_id
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
        # ambig — two isoforms imply different introns in the [66500,66900) gap. ⭐ ON CHR2, so the held
        # record must be stamped ref = 1: a constant stamp would replay these coordinates onto chr1.
        _read("ambig", 66100, [(M, 100), (N, 200), (M, 100)], 66900, True, ref_id=REF_CHR2),
        _read("ambig", 66900, [(M, 100)], 66100, False, ref_id=REF_CHR2),
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

    ⚠ ``mixed``'s junction IS annotated, so it is certified RNA, the unspliced hypothesis is ruled out, and
    its single surviving hypothesis DEPOSITS. ``near`` and ``ambig`` are held instead — U3 and U5 below.
    """
    assert _lengths(payload) == {L_PURE: 1, L_MIXED: 1}, (
        "the deposited fragment lengths are not the molecules this BAM determines. 700 bp for `mixed` "
        "means the mate-gap intron was never cut (the old SPLICE_UNSPLICED gate) — "
    )
    assert payload.qc.deposited == 2
    # ⭐ Conservation: nothing is discarded. Four fragments in, two deposited and two held.
    assert payload.qc.deferred_undetermined_gap == payload.deferred.n_fragments == _DEFERRED
    assert payload.qc.dropped_too_long == payload.qc.dropped_empty == 0


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

    ⭐ Asserted on the ENUMERATED COORDINATES, which is stronger than asserting a deposited length. The
    fragment is now held (its junction is unannotated, so the genomic hypothesis survives too), and what the
    bank must show is that the gap the CIGAR explained was left alone: the observed list is exactly the
    fragment's own 198 bp intron, and the one implied intron is the one in the MATE gap and nothing else.
    """
    (record,) = [
        i for i in range(payload.deferred.n_fragments) if int(payload.deferred.start[i]) == 62_100
    ]
    observed = [tuple(pair) for pair in payload.deferred.observed_introns_of(record).tolist()]
    assert observed == [(62_202, 62_400)], (
        f"the observed intron list is {observed}; the annotated [62200,62400) appearing here means a NEAR "
        f"match was substituted for the intron the CIGAR actually stated"
    )
    hypotheses = _hypotheses_of(payload, record)
    assert [introns for introns, _ in hypotheses] == [[(62_600, 62_800)], []], (
        "the spliced hypothesis must imply exactly the MATE-GAP intron; a second implied intron near "
        "62200 is the near-match substitution, and an empty spliced path is the gap intron going unfound"
    )
    assert [length for _, length in hypotheses] == [L_NEAR, 702], (
        f"the hypotheses must imply L={L_NEAR} (spliced) and 702 (genomic). L={L_MIXED} for the spliced "
        f"one means the annotated intron was substituted for the observed one"
    )


# ── U5 · the deposit is refused when the candidates disagree ───────────────────────────────────────


def test_U5_a_mixed_fragment_with_disagreeing_candidates_is_DEFERRED_and_HELD(payload):
    """⛔ ``ambig`` carries an observed splice, so it never reached this test before.

    ``t_two_a`` puts [66600,66800) in the mate gap and ``t_two_b`` puts [66600,66700) there — a 100 bp
    difference in ``L`` for one molecule. There is no answer to deposit, and picking either is picking a
    fragment length at random.

    ⚠ **Not "dropped".** The fragment is retained in full, with both paths and the transcript that supports
    each, so the second pass can score them against a fragment-length distribution the first pass did not
    have. A counter alone would size the population; the bank is what makes it recoverable.
    """
    (record,) = [
        i for i in range(payload.deferred.n_fragments) if int(payload.deferred.start[i]) == 66_100
    ]
    hypotheses = _hypotheses_of(payload, record)
    assert [introns for introns, _ in hypotheses] == [[(66_600, 66_800)], [(66_600, 66_700)]], (
        "both isoforms' implied intron sets must be held — the second pass cannot choose between answers "
        "it was not given"
    )
    assert [length for _, length in hypotheses] == [500, 600], (
        "a 100 bp difference in L for one molecule"
    )
    # ⛔ Its own junction IS annotated, so the molecule is certified RNA: the genomic hypothesis is dead and
    # the open question is purely WHICH STRUCTURE. `near`, whose junction is 2 bp off the annotation, is the
    # other subclass. The two arms are counted apart, and that is what the census is for.
    assert payload.gap_resolution.gap_deferred_which_introns == 1
    assert payload.gap_resolution.gap_deferred_rna_or_gdna == 1
    assert payload.gap_resolution.gap_deferred_both == 0
    assert payload.gap_resolution.deferred == payload.qc.deferred_undetermined_gap == _DEFERRED
    # Every supporting transcript is recorded: the second pass weights a path by their abundance.
    supporting = [
        payload.deferred.supporting_t_of(h).tolist()
        for h in range(
            int(payload.deferred.hypothesis_offsets[record]),
            int(payload.deferred.hypothesis_offsets[record + 1]),
        )
    ]
    assert all(len(t) == 1 for t in supporting), (
        f"each path has exactly one supporter here: {supporting}"
    )
    assert supporting[0] != supporting[1], "and they are DIFFERENT transcripts"


def test_EVERY_HELD_RECORD_IS_STAMPED_WITH_ITS_OWN_REFERENCE(payload):
    """⛔ **The stamp the drain replays onto, and nothing else in the suite was checking it.**

    Each held record carries the reference it came from, because the second pass re-enters ``deposit`` on
    *that* reference's accumulator — and a coordinate is only meaningful against its own cut axis.

    ⚠ **Measured.** A perturbation giving every reference in the ``AccumulatorSet`` the id ``0`` passed the
    entire suite, 1860 tests, because every fixture was single-contig or deferred only on reference 0. A
    constant is indistinguishable from a correct value until two references both hold something. ``g7``
    lives on chr2 for exactly this reason.

    ⚠ The extent check is the half that matters most: ``[66100,67000)`` is a perfectly plausible interval on
    chr1 too, so "the record is stamped 1" and "the record's coordinates lie inside reference 1" are
    different statements and a wrong stamp fails the second one against the real partition.
    """
    deferred = payload.deferred
    stamped = sorted(int(deferred.ref[i]) for i in range(deferred.n_fragments))
    assert stamped == [REF_CHR1, REF_CHR2], (
        f"the held records are stamped {stamped}; `near` is on chr1 and `ambig` is on chr2, so a constant "
        f"stamp — or a swapped one — shows up here and nowhere else"
    )
    cuts, offsets = payload.cut_positions, payload.ref_cut_offsets
    for i in range(deferred.n_fragments):
        ref = int(deferred.ref[i])
        lo, hi = int(offsets[ref]), int(offsets[ref + 1])
        assert hi > lo, f"record {i} names reference {ref}, which has no cuts at all"
        assert (
            int(cuts[lo]) <= int(deferred.start[i]) < int(deferred.end[i]) <= int(cuts[hi - 1])
        ), (
            f"record {i} spans [{int(deferred.start[i])},{int(deferred.end[i])}) but reference {ref} "
            f"covers [{int(cuts[lo])},{int(cuts[hi - 1])}) — the drain would replay it off the axis"
        )


# ── the RNA pool is keyed on DETERMINACY, not on provenance ────────────────────────────────────────


def test_the_RNA_pool_holds_every_fragment_whose_PATH_IS_DETERMINED(payload):
    """⭐ **TRAPS: a-variance-cannot-fix-a-bias IS DELETED, and this is the two-sided gate that replaces it**.

    TRAPS: a-variance-cannot-fix-a-bias barred a fragment whose splice was *inferred* rather than sequenced, on the grounds that a length
    partly inferred from the annotation is a product of the model the pool is used to fit. The purity
    argument was real; its price was measured and it is larger than what it buys. On the chr22 pilot the
    pool reads **+0.67 % mean / +2.40 % sd** against truth under determinacy and **−9.58 % / −22.46 %**
    under provenance, because barring inferred lengths preferentially bars the fragments whose mates sit
    far apart. ⛔ **A purity filter on a length pool is a length filter.**

    ⚠ What replaces it is stronger, not weaker: a fragment reaches the pool only when **exactly one
    hypothesis survived**, so its ``L`` is not in doubt at all — however it was arrived at.

    ⚠ Two-sided on purpose. ``pure`` (splice sequenced end to end) and ``mixed`` (200 bp of its ``L``
    inferred, but only one possible path) must BOTH be in; ``near`` and ``ambig`` must be out, and they are
    out because they were **deferred**, not because the pool stopped filling. A one-sided version of this
    test passes on an accumulator that emits no pool at all.
    """
    rna = payload.pool_lengths[POOL_RNA_SPLICED]
    assert {int(i): int(rna[i]) for i in np.nonzero(rna)[0]} == {L_PURE: 1, L_MIXED: 1}, (
        "the pure-RNA length pool must hold every fragment whose path was DETERMINED — including the one "
        "whose gap intron was inferred, because with one surviving hypothesis its L is not in doubt"
    )
    assert int(rna.sum()) == payload.qc.deposited, (
        "every deposited fragment here used an annotated junction, so the pool and the deposit total "
        "coincide; a smaller pool means a fragment was silently barred"
    )


# ── U4 · the classification does not move ──────────────────────────────────────────────────────────


def test_U4_the_splice_census_is_unchanged_by_gap_intron_detection(
    gap_intron_index, gap_intron_bam
):
    """Detection is unconditional; the ``SPLICE_IMPLICIT`` **promotion** stays unspliced-only.

    ``splice_type`` feeds scoring, the buffer, the strand training and — since TRAPS: pure-and-length-censored — ``rigel report``'s
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
