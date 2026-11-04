"""An annotated intron in the mate gap is found on every fragment, and never a near-miss for it.

The tool has one definition of fragment length, the accumulator's ``L`` (the length of the fragment's own path), and
this module checks one of its inputs. Implicit-splice detection runs on every fragment, not only those the resolver
called unspliced: a fragment with an observed CIGAR-N splice still has an unsequenced mate gap, and an annotated intron
inside it belongs outside ``L``. The gap finder walks consecutive aligned blocks and emits every hole (a CIGAR-N intron
is a hole too), dropping gaps that exactly match an observed intron. A near-miss junction never reaches the gap search,
because the scan withholds a hit whose junctions no annotated transcript holds, so every candidate of a usable hit holds
its observed junctions exactly.

The four fragments are written by hand because every number asserted is a length to the base pair. Two deposit and their
``L`` is read off ``deposited_lengths``; the other two are held, so their test reads the enumerated hypothesis
coordinates, which also shows the enumeration produced exactly the intended intron set.
"""

from __future__ import annotations

import numpy as np
import pytest
from _index_builder import write_fasta

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import POOL_RNA_SPLICED
from rigel.splice import SpliceType, census_field

GENOME = 70_000

#: Two contigs, and the second matters: ``g7``, the disagreeing-isoform gene and the only one that defers by structure,
#: lives on chr2, so the deferred bank must hold a record stamped ``ref = 1``. A single-contig fixture cannot tell a
#: correct reference stamp from a constant, and the second pass replays each held record onto its reference's coordinate
#: axis.
#:
#: g6/t_mixed: three exons, so a fragment can carry an observed splice and a gap intron.
#: g7/t_two_{a,b}: two isoforms whose gap introns disagree, the deferral case.
#: (GTF is 1-based inclusive; the comment on each row gives its 0-based half-open geometry.)
GTF = (
    # t_mixed exons (62000,62200) (62400,62600) (62800,63000); introns (62200,62400) (62600,62800)
    'chr1\ttest\texon\t62001\t62200\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62401\t62600\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    'chr1\ttest\texon\t62801\t63000\t.\t+\t.\tgene_id "g6"; transcript_id "t_mixed";\n'
    # t_two_a exons (66000,66200) (66400,66600) (66800,67000); gap intron (66600,66800), on chr2
    'chr2\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr2\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    'chr2\ttest\texon\t66801\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_a";\n'
    # t_two_b: same first two exons, but its gap intron is (66600,66700): a different L
    'chr2\ttest\texon\t66001\t66200\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr2\ttest\texon\t66401\t66600\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
    'chr2\ttest\texon\t66701\t67000\t.\t+\t.\tgene_id "g7"; transcript_id "t_two_b";\n'
)

#: Which reference each gene sits on, so the gates can say which coordinate axis a record replays onto.
REF_CHR1, REF_CHR2 = 0, 1

# ── the four fragments, and the L each one must produce ────────────────────────────────────────────
#
#   name      read1 CIGAR at its start          read2       extent      introns cut               L
#   pure      62100  100M 200N 100M             62500 100M  [62100,62600)  (62200,62400)        300
#   mixed     62100  100M 200N 100M             62900 100M  [62100,63000)  + (62600,62800)      500
#   implicit  62100  100M                       62400 100M  [62100,62500)  (62200,62400) or none  200 or 400
#   ambig     66100  100M 200N 100M             66900 100M  (the two isoforms disagree)    held
#
#: ``pure`` is the control that must not move: an observed splice and no mate gap at all.
L_PURE = 300
#: ``mixed`` is the population this module is about: an observed splice plus a gap intron.
L_MIXED = 500

#: ``implicit`` and ``ambig`` are deferred, not deposited, for two different reasons:
#:
#: * ``implicit`` sequenced no splice, so it is not certified RNA. Its mate gap is exactly ``t_mixed``'s first intron,
#:   and the unspliced (genomic) hypothesis survives beside the spliced one: the molecule may be gDNA, or RNA unspliced
#:   there, so ``L`` is undetermined.
#: * ``ambig``'s two isoforms imply different introns in the same gap. Its sequenced sj is annotated and gDNA cannot be
#:   spliced, so the molecule is certified RNA and only the structure is open.
_DEFERRED = 2


def _covered_length(start: int, end: int, introns) -> int:
    """``L`` by integer set arithmetic, a different algorithm from the accumulator's segment walk. Deliberately not the
    specification's ``_hypothesis_length``: a validator that calls the implementation's helper would validate nothing.
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

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("gap_introns")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    write_fasta(fasta, {"chr1": GENOME, "chr2": GENOME})
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
        # Proper FR pair: R1 forward with a reverse mate, R2 reverse. `build_fragment` flips R2's strand, so both mates
        # land in one (ref, strand) group and merge into a single extent.
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
        # pure: observed splice, mates abut so there is no mate gap to search
        _read("pure", 62100, [(M, 100), (N, 200), (M, 100)], 62500, True),
        _read("pure", 62500, [(M, 100)], 62100, False),
        # mixed: observed splice and an annotated intron inside the [62500,62900) mate gap
        _read("mixed", 62100, [(M, 100), (N, 200), (M, 100)], 62900, True),
        _read("mixed", 62900, [(M, 100)], 62100, False),
        # implicit: no splice sequenced; the [62200,62400) mate gap is exactly t_mixed's first intron
        _read("implicit", 62100, [(M, 100)], 62400, True),
        _read("implicit", 62400, [(M, 100)], 62100, False),
        # ambig: two isoforms imply different introns in the [66500,66900) gap. On chr2, so the held record must be
        # stamped ref = 1.
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
    """``mixed`` spans [62100,63000), 900 bp of genome, and is a 500 bp molecule: 200 bp was spliced out and sequenced
    as CIGAR-N, another 200 bp spliced out and never sequenced. Cutting only the sequenced intron would give 700 bp. Its
    sj is annotated, so it is certified RNA, the unspliced hypothesis is ruled out, and its single surviving hypothesis
    deposits. ``implicit`` and ``ambig`` are held (U3 and U5).
    """
    assert _lengths(payload) == {L_PURE: 1, L_MIXED: 1}, (
        "the deposited fragment lengths are not the molecules this BAM determines. 700 bp for `mixed` "
        "means the mate-gap intron was never cut (the old SPLICE_UNSPLICED gate) — "
    )
    assert payload.qc.deposited == 2
    # Conservation: nothing is discarded. Four fragments in, two deposited and two held.
    assert payload.qc.deferred_undetermined_gap == payload.deferred.n_fragments == _DEFERRED
    assert payload.qc.dropped_too_long == payload.qc.dropped_empty == 0


def test_U1b_the_unspliced_control_is_untouched(payload):
    """``pure`` has an observed splice and no mate gap: nothing about it may move."""
    assert _lengths(payload)[L_PURE] == 1


# ── U3 · an unspliced pair keeps both readings of its gap ─────────────────────────────────────────


def test_U3_an_unspliced_pair_keeps_both_readings_of_its_gap(payload):
    """``implicit`` sequenced no splice; its mate gap [62200,62400) is exactly ``t_mixed``'s first intron. The gap
    search finds the annotated intron (the spliced reading, L=200); nothing certifies the molecule as RNA, so the
    genomic reading (L=400, no intron) stays beside it and the record is held for the second pass.
    """
    (record,) = [
        i for i in range(payload.deferred.n_fragments) if int(payload.deferred.start[i]) == 62_100
    ]
    assert payload.deferred.observed_introns_of(record).tolist() == []
    hypotheses = _hypotheses_of(payload, record)
    assert hypotheses == [([(62_200, 62_400)], 200), ([], 400)], (
        "the spliced reading must cut exactly the gap's annotated intron, and the genomic one must stay"
    )


# ── U5 · the deposit is refused when the candidates disagree ───────────────────────────────────────


def test_U5_a_mixed_fragment_with_disagreeing_candidates_is_DEFERRED_and_HELD(payload):
    """``ambig`` carries an observed splice and a disputed mate gap: ``t_two_a`` puts [66600,66800) in it and
    ``t_two_b`` puts [66600,66700), a 100 bp difference in ``L`` for one molecule. Picking either would pick a fragment
    length at random, so the fragment is held in full with both paths for the second pass to score against a
    fragment-length distribution the first pass did not have.
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
    # Its own sj is annotated, so the molecule is certified RNA: the genomic hypothesis is dead and only the structure
    # is open. `implicit`, which sequenced no splice, is the other subclass; the census counts the two apart.
    assert payload.gap_resolution.gap_deferred_which_introns == 1
    assert payload.gap_resolution.gap_deferred_rna_or_gdna == 1
    assert payload.gap_resolution.gap_deferred_both == 0
    assert payload.gap_resolution.deferred == payload.qc.deferred_undetermined_gap == _DEFERRED


def test_EVERY_HELD_RECORD_IS_STAMPED_WITH_ITS_OWN_REFERENCE(payload):
    """Each held record carries its source reference, because the second pass re-enters ``deposit`` on that reference's
    accumulator and a coordinate means something only on its own axis. A constant stamp is indistinguishable from a
    correct one until two references both hold records, which is why ``g7`` lives on chr2. The extent check matters
    most: ``[66100,67000)`` is a plausible interval on chr1 too, so a wrong stamp fails only when the coordinates are
    checked against the real partition.
    """
    deferred = payload.deferred
    stamped = sorted(int(deferred.ref[i]) for i in range(deferred.n_fragments))
    assert stamped == [REF_CHR1, REF_CHR2], (
        f"the held records are stamped {stamped}; `implicit` is on chr1 and `ambig` is on chr2, so a constant "
        f"stamp — or a swapped one — shows up here and nowhere else"
    )
    region_bounds, offsets = payload.region_bounds, payload.ref_region_bound_offsets
    for i in range(deferred.n_fragments):
        ref = int(deferred.ref[i])
        lo, hi = int(offsets[ref]), int(offsets[ref + 1])
        assert hi > lo, f"record {i} names reference {ref}, which has no region_bounds at all"
        assert (
            int(region_bounds[lo])
            <= int(deferred.start[i])
            < int(deferred.end[i])
            <= int(region_bounds[hi - 1])
        ), (
            f"record {i} spans [{int(deferred.start[i])},{int(deferred.end[i])}) but reference {ref} "
            f"covers [{int(region_bounds[lo])},{int(region_bounds[hi - 1])}) — the drain would replay it off the axis"
        )


# ── the RNA pool is keyed on determinacy, not on provenance ────────────────────────────────────────


def test_the_RNA_pool_holds_every_fragment_whose_PATH_IS_DETERMINED(payload):
    """Pool membership is determinacy, not provenance, checked from both sides.

    Barring fragments whose splice was inferred rather than sequenced would preferentially bar those whose mates sit far
    apart, shifting the pool's own length distribution. Instead a fragment reaches the pool only when exactly one
    hypothesis survived, so its ``L`` is not in doubt however it was derived. ``pure`` (splice sequenced end to end) and
    ``mixed`` (200 bp of its ``L`` inferred, one possible path) must both be in; ``implicit`` and ``ambig`` must be out
    because they were deferred. A one-sided test would pass on an accumulator that emits no pool at all.
    """
    rna = payload.pool_lengths[POOL_RNA_SPLICED]
    assert {int(i): int(rna[i]) for i in np.nonzero(rna)[0]} == {L_PURE: 1, L_MIXED: 1}, (
        "the pure-RNA length pool must hold every fragment whose path was DETERMINED — including the one "
        "whose gap intron was inferred, because with one surviving hypothesis its L is not in doubt"
    )
    assert int(rna.sum()) == payload.qc.deposited, (
        "every deposited fragment here used an annotated sj, so the pool and the deposit total "
        "coincide; a smaller pool means a fragment was silently barred"
    )


# ── U4 · the classification does not move ──────────────────────────────────────────────────────────


def test_U4_the_splice_census_is_unchanged_by_gap_intron_detection(
    gap_intron_index, gap_intron_bam
):
    """Detection is unconditional; the ``SPLICE_IMPLICIT`` label stays unspliced-only.

    ``splice_type`` feeds scoring, the deferred bank, the strand training and ``rigel report``'s census. The three
    fragments with an observed CIGAR-N splice count as observed splices, though ``mixed`` and ``ambig`` also hold an
    intron in their mate gap; ``implicit``, which sequenced no splice, is the one implicit splice.
    """
    _, index = gap_intron_index
    stats, _, _, _ = scan_and_buffer(gap_intron_bam, index, BamScanConfig(sj_strand_tag="auto"))

    def census(stype):
        return getattr(stats, census_field(stype))

    assert census(SpliceType.SPLICED_ANNOT) == 3, "pure, mixed and ambig are annotated splices"
    assert census(SpliceType.SPLICED_IMPLICIT) == 1, (
        "only `implicit` may be implicit: an observed splice is never RE-LABELLED, this work is about L"
    )
    assert census(SpliceType.UNSPLICED) == 0
