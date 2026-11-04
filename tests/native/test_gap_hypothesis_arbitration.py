"""The gap hypothesis: the paths a fragment's unsequenced mate gap could take through the annotation, how the
accumulator arbitrates between them, and what strand each carries.

A gap may hold no intron, one or several, and which cannot be observed. Deciding is a likelihood question that needs a
fragment-length distribution the first pass does not have, so the first pass enumerates: one surviving hypothesis
deposits, more are held for the second pass. The first block tests that arbitration against the reference accumulator
(the compatible paths and their lengths, the deferral, the single-transcript deposit, the genomic hypothesis as the
longest and as the empty path, the span limit, conservation, and the deferred queue's flattening). The second block
tests strand.
"""

from __future__ import annotations

import numpy as np
import pytest
from _index_builder import write_fasta

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.types import Strand

from ._accumulator_reference import (
    Accumulator,
    DepositOutcome,
    GapHypothesis,
    Partition,
)


#: One reference, with a bound at every coordinate the worked example names, so both paths land on
#: real regions.
REGION_BOUNDS = [0, 1000, 2000, 3000, 3050, 4000, 5000, 6000]
TYPES = [0, 2, 1, 2, 1, 2, 0]


def _acc(**kw):
    return Accumulator(Partition.from_region_bounds([REGION_BOUNDS], region_types=[TYPES]), **kw)


# ── the worked example, to the base pair ──────────────────────────────────────────────────────────
#
#     TA  exons (1000,2000)  (3000,3050)  (4000,5000)   introns (2000,3000) (3050,4000)
#     TB  exons (1000,2000)               (4000,5000)   intron  (2000,4000)
#
#     fragment  blocks [1800,1950)                [4050,4200)      span 2400
#                      ==========|~~ unsequenced ~~|=========
#
#: TA's path crosses an exon, (3000,3050), that no read touched. A gap hypothesis is a path through the
#: annotation, not "an intron", which is why the first pass cannot resolve it.
TA = GapHypothesis(((2000, 3000), (3050, 4000)), sj_strand=Strand.POS)
TB = GapHypothesis(((2000, 4000),), sj_strand=Strand.POS)
L_TA = 2400 - (1000 + 950)  # 450
L_TB = 2400 - 2000  # 400


def test_the_two_compatible_paths_have_the_lengths_the_owner_computed():
    """Checked apart from the arbitration: a bank holding the right fragment with wrong hypothesis lengths still defers
    it, and the second pass would then choose between two wrong answers."""
    assert (L_TA, L_TB) == (450, 400)
    acc = _acc()
    for path, expected in ((TA, L_TA), (TB, L_TB)):
        length, _introns = acc._hypothesis_length(1800, 4200, (), path)
        assert length == expected


def test_TWO_COMPATIBLE_PATHS_are_BUFFERED_and_neither_is_chosen():
    """450 bp and 400 bp are both real molecules and nothing sequenced says which."""
    acc = _acc()
    outcome = acc.deposit(0, 1800, 4200, hypotheses=(TA, TB))
    assert outcome is DepositOutcome.DEFERRED
    t = acc.tally
    assert int(t.region_start_count.sum()) == 0, "an undetermined fragment locates nowhere"
    assert int(t.deposited_lengths.sum()) == 0, "and has no length to bin"
    assert int(t.pool_lengths.sum()) == 0
    assert len(t.deferred) == 1
    held = t.deferred[0].placements[0]
    # Stored whole, with both paths: the second pass cannot choose between answers it was not given.
    assert held.hypotheses == (TA, TB)
    # Both are spliced and neither is the unspliced hypothesis: the molecule is certified RNA and only the structure is
    # open.
    assert t.gap_resolution["gap_deferred_which_introns"] == 1


def test_the_SAME_fragment_deposits_when_only_ONE_transcript_is_compatible():
    """The ambiguity is a property of the annotation, not of the fragment: drop TB from the candidate set and the same
    molecule has one answer, so it deposits at TA's 450 bp.
    """
    acc = _acc()
    assert acc.deposit(0, 1800, 4200, hypotheses=(TA,)) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(np.nonzero(t.deposited_lengths)[0][0]) == L_TA
    assert not t.deferred
    assert t.gap_resolution["gap_resolved_spliced"] == 1


def test_the_GENOMIC_hypothesis_is_ALWAYS_the_LONGEST():
    """A spliced hypothesis cuts bases the genomic one keeps, so ``L_spliced <= L_genomic`` for any extent, observed
    introns and implied path. The single arbitration filter is ``L <= max_fragment_length``, so if the genomic path
    survives every spliced path survives too, and the survivors are never exactly ``{genomic}`` while a spliced path was
    offered.

    That is why :class:`GapResolution` has no ``RESOLVED_UNSPLICED``: no fragment could enter it. A future filter that
    broke the ordering fails here. The test is randomised over the coordinate space because the claim is universal;
    introns are drawn unsorted and may be reversed, zero-length, overlapping, nested or outside the fragment, every
    configuration ``_normalise_introns`` handles.
    """
    import random

    acc = _acc()
    rng = random.Random(20260801)
    checked = 0
    for _ in range(20_000):
        start, end = sorted(rng.sample(range(0, 6001), 2))
        if end == start:
            continue
        observed = tuple(tuple(rng.sample(range(0, 6001), 2)) for _ in range(rng.randint(0, 2)))
        genomic, _introns = acc._hypothesis_length(start, end, observed, GapHypothesis())
        for _ in range(rng.randint(1, 3)):
            spliced_path = GapHypothesis(
                tuple(tuple(rng.sample(range(0, 6001), 2)) for _ in range(rng.randint(1, 3)))
            )
            spliced, _introns = acc._hypothesis_length(start, end, observed, spliced_path)
            checked += 1
            assert spliced <= genomic, (
                f"[{start},{end}) observed={observed} path={spliced_path.introns}: the spliced path "
                f"L={spliced} exceeds the genomic L={genomic}. Cutting bases cannot lengthen a molecule, "
                f"so this breaks the ordering GapResolution's shape depends on."
            )
    assert checked > 15_000, f"only {checked} pairs compared; the sweep is not doing its job"


def test_the_genomic_hypothesis_is_the_EMPTY_path_and_needs_no_flag():
    """Cutting nothing means the gap is real template: the molecule is gDNA, or RNA unspliced there, the same span. So
    "could this be genomic?" needs no separate flag, and a synthetic span is not a separate candidate: it is this
    hypothesis.

    The limit is raised above the 2400 bp span on purpose. At the default 1000 the unspliced hypothesis is ruled out by
    length and the fragment deposits, which the next test covers.
    """
    acc = _acc(max_fragment_length=3000)
    outcome = acc.deposit(0, 1800, 4200, hypotheses=(TA, GapHypothesis()))
    assert outcome is DepositOutcome.DEFERRED
    # The genomic path against exactly one spliced path: the open question is RNA or gDNA, the composition question
    # calibration answers.
    assert acc.tally.gap_resolution["gap_deferred_rna_or_gdna"] == 1
    assert acc.tally.gap_resolution["gap_deferred_which_introns"] == 0


def test_a_fragment_with_NO_gap_hypothesis_is_not_in_the_umbrella_at_all():
    """A fragment that never had a question to answer must not be counted as having answered one, or the umbrella's
    denominator becomes the whole library."""
    acc = _acc()
    assert acc.deposit(0, 1800, 2400) is DepositOutcome.DEPOSITED
    assert sum(acc.tally.gap_resolution.values()) == 0


# ── the span rule is the max_fragment_length filter applied to the unspliced hypothesis ─────────────────────


def test_a_span_over_the_limit_RULES_OUT_the_genomic_hypothesis():
    """ "If the genomic span exceeds max_fragment_length, it is RNA" is not a separate rule: the unspliced hypothesis's
    ``L`` is that span, so the ordinary hypothesis filter deletes it. Span 2400 against a limit of 1000 leaves TA alone,
    and it deposits at 450.
    """
    acc = _acc(max_fragment_length=1000)
    assert acc.deposit(0, 1800, 4200, hypotheses=(TA, GapHypothesis())) is DepositOutcome.DEPOSITED
    t = acc.tally
    assert int(np.nonzero(t.deposited_lengths)[0][0]) == L_TA
    assert t.gap_resolution["gap_resolved_spliced"] == 1
    assert not t.deferred


def test_under_the_limit_the_unspliced_hypothesis_SURVIVES_and_the_fragment_is_DEFERRED():
    """The other side of the same rule: a short unspliced fragment with an intron in its gap is compatible with DNA and
    with RNA, so it is always deferred."""
    acc = _acc(max_fragment_length=1000)
    # span 900, and a 400 bp intron inside the gap -> unspliced L = 900, spliced L = 500. Both possible.
    outcome = acc.deposit(
        0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())
    )
    assert outcome is DepositOutcome.DEFERRED
    assert acc.tally.gap_resolution["gap_deferred_rna_or_gdna"] == 1


# ── conservation: nothing offered may be lost ──────────────────────────────────────────────────────


def test_EVERY_OFFERED_FRAGMENT_IS_ACCOUNTED_FOR_EXACTLY_ONCE():
    """``deposited + deferred + dropped_* == offered``, exactly.

    Until the second pass drains the bank the tally is thinner by design (the ambiguous mass is retained, not yet
    deposited), so conservation, not calibration accuracy, is what shows nothing was lost. Every branch is exercised,
    since a conservation identity over a single route is satisfied by any bookkeeping.
    """
    acc = _acc(max_fragment_length=1000)
    offered = 0

    def offer(*args, **kw):
        nonlocal offered
        offered += 1
        return acc.deposit(*args, **kw)

    assert offer(0, 1800, 2400) is DepositOutcome.DEPOSITED
    assert offer(0, 1800, 4200, hypotheses=(TA,)) is DepositOutcome.DEPOSITED
    assert offer(0, 1800, 4200, hypotheses=(TA, TB)) is DepositOutcome.DEFERRED
    assert offer(0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())) is (
        DepositOutcome.DEFERRED
    )
    assert offer(0, 1800, 4200) is DepositOutcome.TOO_LONG
    assert offer(0, 1800, 2400, align_strand=Strand.AMBIGUOUS) is DepositOutcome.STRAND_UNDEFINED
    assert offer(0, 9000, 9500) is DepositOutcome.EMPTY

    qc = acc.tally.qc
    assert sum(qc[outcome.value] for outcome in DepositOutcome) == offered
    # The deferred queue holds the fragments the counter claims.
    assert len(acc.tally.deferred) == qc["deferred_undetermined_gap"] == 2
    assert int(acc.tally.region_start_count.sum()) == qc["deposited"] == 2


def test_the_gap_resolution_SUBCLASSES_CLOSE_against_the_umbrella_and_the_deferred_queue():
    """Two identities: the subclasses partition the umbrella, and the three ``DEFERRED_*`` subclasses are exactly what
    was held.
    """
    acc = _acc(max_fragment_length=1000)
    acc.deposit(0, 1800, 4200, hypotheses=(TA,))  # determined, spliced
    acc.deposit(0, 1800, 4200, hypotheses=(TA, TB))  # deferred, path only
    acc.deposit(
        0, 1800, 2700, hypotheses=(GapHypothesis(((2100, 2500),)), GapHypothesis())
    )  # deferred, component
    acc.deposit(
        0,
        1800,
        2700,
        hypotheses=(
            GapHypothesis(((2100, 2500),)),
            GapHypothesis(((2200, 2600),)),
            GapHypothesis(),
        ),
    )
    acc.deposit(0, 1800, 2400)  # no hypothesis — not in the umbrella

    census = acc.tally.gap_resolution
    assert census == {
        "gap_resolved_spliced": 1,
        "gap_deferred_rna_or_gdna": 1,
        "gap_deferred_which_introns": 1,
        "gap_deferred_both": 1,
    }
    deferred = sum(v for k, v in census.items() if k.startswith("gap_deferred"))
    assert deferred == acc.tally.qc["deferred_undetermined_gap"] == len(acc.tally.deferred)


# ── the flattened deferred queue is what the payload carries, and it is specified in the reference ──────────


#: The two deferred fragments the flattening tests use. ``SHORT`` sorts before ``LONG`` (same start, smaller end), so
#: the expected arrays below are in that order, not in deposit order.
LONG = (0, 1800, 4200, (TA, TB))
SHORT = (0, 1800, 2700, (GapHypothesis(((2100, 2500),)), GapHypothesis()))


def test_the_deferred_queue_FLATTENS_to_a_CSR_that_round_trips():
    """The reference stores records for readability; the payload carries arrays. The flattening is specified in the
    reference and checked here, so the two are one representation."""
    acc = _acc(max_fragment_length=1000)
    for ref, start, end, hypotheses in (LONG, SHORT):
        acc.deposit(ref, start, end, hypotheses=hypotheses)

    arrays = acc.tally.deferred_arrays()
    assert arrays["start"].tolist() == [1800, 1800]
    assert arrays["end"].tolist() == [2700, 4200]
    # Offsets are cumulative and start at 0, so an empty deferred queue is [0] and never [].
    assert arrays["hypothesis_offsets"].tolist() == [0, 2, 4]
    # SHORT's spliced path has one intron and its genomic path none; then TA has two and TB one.
    assert arrays["hypothesis_intron_offsets"].tolist() == [0, 1, 1, 3, 4]
    assert arrays["hypothesis_introns"].tolist() == [
        2100, 2500,               # SHORT, spliced
                                  # SHORT, genomic: no introns, which is what makes it genomic
        2000, 3000, 3050, 4000,   # TA
        2000, 4000,               # TB
    ]  # fmt: skip
    assert Accumulator(
        Partition.from_region_bounds([REGION_BOUNDS], region_types=[TYPES])
    ).tally.deferred_arrays()["hypothesis_offsets"].tolist() == [0]


def test_the_FLATTENED_queue_DOES_NOT_DEPEND_ON_DEPOSIT_ORDER():
    """Every other bank is a per-object sum, so uint32 counts merge exactly and float64 fractions to a derived tolerance
    in any chunk order. The deferred queue is a list: concatenating per-worker queues would give different bytes at 1,
    2, 4 and 8 workers with identical contents. Sorting on each record's own content is the canonical form, which is why
    the sort lives in the reference and not in the exporter.
    """
    forward, backward = _acc(max_fragment_length=1000), _acc(max_fragment_length=1000)
    for ref, start, end, hypotheses in (LONG, SHORT):
        forward.deposit(ref, start, end, hypotheses=hypotheses)
    for ref, start, end, hypotheses in (SHORT, LONG):
        backward.deposit(ref, start, end, hypotheses=hypotheses)

    assert [f.placements[0].end for f in forward.tally.deferred] != [
        f.placements[0].end for f in backward.tally.deferred
    ], "the two accumulators must actually differ in deposit order, or this passes vacuously"
    a, b = forward.tally.deferred_arrays(), backward.tally.deferred_arrays()
    assert a.keys() == b.keys()
    for name in a:
        assert np.array_equal(a[name], b[name]), f"{name} depends on the order fragments arrived in"


# ── strand in the enumeration: what pins it, what leaves it open, what contradicts it ───────────
#
# Splice junctions are stranded, so one detected in a fragment gives its strand and the gap hypotheses are constrained
# to it. An unspliced fragment could be either strand, though the library's strand specificity may constrain it. These
# tests go through the scan: driving ``FragmentResolver`` directly leaves ``t_strand_arr_`` empty, so every implied
# strand would read ``NONE`` and a strand test would pass against nothing. The index warns about the fixture, as
# expected: two of the three gene pairs annotate the same (donor, acceptor) on both strands, which the index calls
# biologically impossible (GT..AG reverse-complements to CT..AC), and that configuration is what these tests need.


GENOME = 8_000

#: Three pairs, each isolating one behaviour.
#:
#: tP/tM share the observed sj ``(1200,1400)`` on opposite strands and cover the same sequenced blocks, so only the
#: sequenced motif can separate them. Their gap introns differ in width (400 vs 200 bp), so the deposited ``L`` says
#: which one was believed.
#:
#: tQ/tR imply the same gap intron on opposite strands: one path whose two strands disagree.
GTF = (
    # tP (+): introns (1200,1400) observed, (1600,2000) in the gap  -> cuts 400
    'chr1\ttest\texon\t1001\t1200\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    'chr1\ttest\texon\t1401\t1600\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    'chr1\ttest\texon\t2001\t2200\t.\t+\t.\tgene_id "gP"; transcript_id "tP";\n'
    # tM (-): same observed sj, but its gap intron is (1700,1900) -> cuts 200
    'chr1\ttest\texon\t1001\t1200\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    'chr1\ttest\texon\t1401\t1700\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    'chr1\ttest\texon\t1901\t2200\t.\t-\t.\tgene_id "gM"; transcript_id "tM";\n'
    # tQ (+) and tR (-): the same gap intron (3200,3600) on both strands, one path, two strands
    'chr1\ttest\texon\t3001\t3200\t.\t+\t.\tgene_id "gQ"; transcript_id "tQ";\n'
    'chr1\ttest\texon\t3601\t3800\t.\t+\t.\tgene_id "gQ"; transcript_id "tQ";\n'
    'chr1\ttest\texon\t3001\t3200\t.\t-\t.\tgene_id "gR"; transcript_id "tR";\n'
    'chr1\ttest\texon\t3601\t3800\t.\t-\t.\tgene_id "gR"; transcript_id "tR";\n'
)

#: The four molecules, and what each one tests. Extents are the merged block span.
#:
#:   name      reads                                   extent        the question
#:   pinned_p  100M 200N 150M @1100  +  100M @2050     [1100,2150)   does a + motif pin it to tP?
#:   pinned_m  the same reads with a - motif           [1100,2150)   ... and a - motif to tM?
#:   open      100M @1450            +  100M @2050     [1450,2150)   no motif -> both strands offered
#:   same_path 100M @3050            +  100M @3650     [3050,3750)   one path, two strands
L_PINNED_P = (2150 - 1100) - 200 - 400  # 450: tP's wide gap intron was cut
L_PINNED_M = (2150 - 1100) - 200 - 200  # 650: tM's narrow one was cut


@pytest.fixture(scope="module")
def scanned(tmp_path_factory):
    """Scan the fixture BAM once and hand back the payload."""
    import pysam

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("gap_strand")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    write_fasta(fasta, {"chr1": GENOME})
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
        # XS is the aligner's genomic motif strand, and it pins the fragment: a GT..AG intron is on +, its reverse
        # complement CT..AC on -, and no annotation is consulted.
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


# ── an observed sj pins the strand ───────────────────────────────────────────────────────────


def test_an_OBSERVED_sj_PINS_the_gap_hypotheses_to_its_own_strand(scanned):
    """``tP`` and ``tM`` cover the same sequenced blocks and share the observed sj's coordinates, so only the motif
    separates them. A ``+`` motif leaves only tP's 400 bp gap intron and a ``-`` motif only tM's 200 bp one, and the
    deposited ``L`` says which was believed. Without the pin both transcripts survive and the fragment is deferred, so
    the count fails as well as the lengths.
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
    # Both pinned fragments resolved; only `open` and `same_path` are held.
    assert scanned.gap_resolution.gap_resolved_spliced == 2
    assert scanned.deferred.n_fragments == 2


# ── an unspliced fragment leaves both strands open ─────────────────────────────────────────────────


def test_an_UNSPLICED_fragment_offers_BOTH_STRANDS(scanned):
    """Nothing sequenced pins the strand, so nothing may narrow it. The two spliced hypotheses here differ in strand and
    in implied length; without the second pass's strand term they would be separated by length alone.
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
        "the genomic hypothesis must also be present — no annotated sj was observed, so the gap "
        "may be real template"
    )


# ── one path claimed by two strands ────────────────────────────────────────────────────────────────


def test_ONE_PATH_claimed_by_BOTH_STRANDS_is_marked_AMBIGUOUS(scanned):
    """``tQ`` (+) and ``tR`` (−) imply the same intron ``(3200,3600)``.

    It is one path, so one hypothesis is the correct count, but the hypothesis carries one ``sj_strand``, and taking the
    first supporter's would assert a strand the evidence does not support (swapping the two GTF lines would flip it).
    ``AMBIGUOUS`` is the codebase's value for contradictory evidence: the fragment-level ``sj_strand`` uses it the same
    way, and ``deposit`` refuses to credit a sj on it. Human annotation never puts one sj coordinate pair on both
    strands, but the case is handled so the answer cannot depend on GTF row order.
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
        f"first supporter's strand makes the answer depend on GTF line order."
    )


def test_the_fixture_reaches_every_branch_it_claims_to(scanned):
    """Each test above reads one record or one length bin; if the scan produced something else, several could pass for
    the wrong reason."""
    assert scanned.qc.deposited == 2, "the two pinned fragments, and only those, must deposit"
    assert scanned.deferred.n_fragments == 2, (
        "the unspliced and the same-path fragments must be held"
    )
    assert scanned.gap_resolution.deferred == 2
    assert scanned.qc.dropped_too_long == 0 and scanned.qc.dropped_empty == 0


# `!certified_rna` in `enumerate_gap_hypotheses`'s unspliced-hypothesis condition is redundant whenever a synthetic span
# is a candidate: the span is single-exon, implies nothing in the gap, and so supplies the genomic path itself (and a
# spliced fragment's observed intron excludes the span, leaving the genomic path absent). The clause stays because it
# states the rule directly and alone covers a single-exon gene, which has no synthetic span; no test here can
# distinguish it.
