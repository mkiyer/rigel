"""⭐ P3 — THE DRAIN: one tally path, conservation, and byte-identity with the specification.

     (the draw), §6 (the drain) · Phase: its P3
    Specification: ``tests/native/_accumulator_reference.py`` — ``Accumulator.drain``

The second pass replays each held fragment with **one** chosen hypothesis. §6.1 makes a strong claim about
how, and the claim is the thing worth gating:

> The drain re-enters ``Accumulator::deposit`` with the chosen hypothesis **alone** — a set of size one, so
> arbitration is degenerate and it deposits (or is rejected by the ordinary rules). There is no second
> deposit implementation, no duplicated crossing logic, and byte-identity with the specification is
> preserved for free.

⭐ **"For free" is checkable, and this module checks it directly**: a fragment drained with hypothesis ``h``
must produce a byte-identically equal tally to the same fragment offered ``hypotheses=(h,)`` in the first
place (:func:`test_draining_a_choice_equals_depositing_it_directly`). If that holds, the drain has no
tally semantics of its own to test — which is the point of building it this way.

⛔ **The drain does NOT extend ``GapCensus``, and that is structural.** The census has no
``gap_resolved_unspliced`` class because pass-one arbitration cannot produce one — a spliced path always
region_bounds bases the genomic path keeps, so the genomic path can never be the sole survivor (S1 checked this over
200,000 random hypothesis sets before deleting the class). ⭐ But the drain *chooses*, and it can choose ∅.
So a naive drain would grow the census by however many draws happened to pick a spliced path while chosen-∅
fragments vanished from it entirely — a census that depends on the RNG.
:func:`test_the_drain_leaves_the_ARBITRATION_census_alone` is that gate.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.scan_payload import ADDITIVE_AXES
from rigel.types import Strand

from ._accumulator_reference import (
    Accumulator as ReferenceAccumulator,
)
from ._accumulator_reference import (
    GapHypothesis,
    Partition,
    Tally,
)


GENOME_REGION_BOUNDS = [0, 1000, 1600, 1700, 1900, 2000, 2200, 4000]
#: Two annotated introns inside the gap, so an unspliced fragment spanning it has three hypotheses.
SJ = ((0, 1600, 2000, int(Strand.POS)), (0, 1700, 1900, int(Strand.NEG)))
MAX_LENGTH = 1000


def _partition() -> Partition:
    return Partition.from_region_bounds([GENOME_REGION_BOUNDS], sj=SJ)


def _hypotheses() -> tuple[GapHypothesis, ...]:
    """The three competing explanations of the gap ``[1600, 2000)``: two introns and the genomic path."""
    return (
        GapHypothesis(introns=((1600, 2000),), sj_strand=int(Strand.POS)),
        GapHypothesis(introns=((1700, 1900),), sj_strand=int(Strand.NEG)),
        GapHypothesis(),
    )


def _offer_held(accumulator: ReferenceAccumulator, n: int = 1) -> None:
    """Offer ``n`` fragments that span the gap unspliced, so each is HELD with three hypotheses."""
    for _ in range(n):
        outcome = accumulator.deposit(
            0,
            1500,
            2100,
            align_strand=int(Strand.POS),
            sj_strand=int(Strand.NONE),
            hypotheses=_hypotheses(),
        )
        assert outcome.value == "deferred_undetermined_gap", (
            f"the fixture must HOLD this fragment or there is nothing to drain; got {outcome}"
        )


def _fresh() -> ReferenceAccumulator:
    return ReferenceAccumulator(_partition(), max_fragment_length=MAX_LENGTH)


def _tally_fields() -> list[str]:
    return [f.name for f in dataclasses.fields(Tally)]


def _compare(left: Tally, right: Tally) -> list[str]:
    """Field names on which two tallies differ. ⚠ Read off ``dataclasses.fields`` so a new channel is
    compared automatically rather than when someone remembers to add it here."""
    differing = []
    for name in _tally_fields():
        a, b = getattr(left, name), getattr(right, name)
        if isinstance(a, np.ndarray):
            same = a.shape == b.shape and a.dtype == b.dtype and np.array_equal(a, b)
        elif name == "deferred":
            same = len(a) == len(b)
        else:
            same = a == b
        if not same:
            differing.append(name)
    return differing


# ── §6.1 · ONE TALLY PATH ──────────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("choice", [0, 1, 2])
def test_draining_a_choice_equals_depositing_it_directly(choice):
    """⭐ **§6.1's claim, checked rather than argued.** Draining a held fragment with hypothesis ``choice``
    must give byte-identically the tally that offering only that hypothesis would have given.

    ⚠ Not a tautology, because the two routes differ in everything except the deposit: one goes through
    arbitration with three hypotheses, a canonical sort and a replay; the other deposits once. If they
    agree on every ``Tally`` field then the drain genuinely has no tally semantics of its own.

    ⛔ The one field that must NOT agree is ``qc``: the drained route counted a `deferred` on the way in.
    That is checked separately in :func:`test_the_drain_consumes_the_bank_and_conserves`, and excluded
    here by name so this gate cannot pass by comparing nothing.
    """
    drained = _fresh()
    _offer_held(drained)
    drained.drain([choice])

    direct = _fresh()
    outcome = direct.deposit(
        0,
        1500,
        2100,
        align_strand=int(Strand.POS),
        sj_strand=int(Strand.NONE),
        hypotheses=(_hypotheses()[choice],),
    )
    assert outcome.value == "deposited"

    differing = set(_compare(drained.tally, direct.tally))
    # `gap_resolution` differs by construction: the direct route arbitrated one hypothesis and the drained
    # route arbitrated three, and the drain deliberately restores pass one's census.
    assert differing <= {"gap_resolution"}, (
        f"a drained deposit must be byte-identical to a direct one on every tally channel; differs on "
        f"{sorted(differing)}"
    )
    assert "region_start_count" not in differing and "deposited_lengths" not in differing


def test_the_drain_consumes_the_bank_and_conserves():
    """⭐ §6.2's conservation, exactly, and the bank empty afterwards."""
    accumulator = _fresh()
    _offer_held(accumulator, n=7)
    held_before = accumulator.tally.qc["deferred_undetermined_gap"]
    assert held_before == 7

    counters = accumulator.drain([0, 1, 2, 0, 1, 2, 0])

    assert counters["offered"] == 7
    dropped = sum(v for k, v in counters.items() if k.startswith("dropped_"))
    assert counters["deposited"] + dropped == counters["offered"], (
        f"every offered fragment must deposit or be rejected exactly once; got {counters}"
    )
    assert counters["chose_genomic"] + counters["chose_spliced"] == 7
    assert len(accumulator.tally.deferred) == 0, "the drain consumes the bank"
    # ⭐ The held counter goes to 0 WITH the bank. The two must describe one population — the payload
    # refuses a bank that disagrees with its counter at the door, and a drained payload gets no exception.
    assert accumulator.tally.qc["deferred_undetermined_gap"] == 0
    assert int(accumulator.tally.region_start_count.sum()) == accumulator.tally.qc["deposited"]
    assert int(accumulator.tally.deposited_lengths.sum()) == accumulator.tally.qc["deposited"]


def test_the_drain_leaves_the_ARBITRATION_census_alone():
    """⛔ The RNG-dependent census. ``_record_gap_resolution`` sends a size-one SPLICED set to
    ``RESOLVED_SPLICED`` and returns early on an all-unspliced one, so a drain that let it run would move
    the census by however many draws picked a spliced path — and chosen-∅ fragments would vanish from it.

    ⭐ The census must read the same whatever was chosen, and the drain's own axis carries the difference.
    """
    censuses = []
    for choices in ([0, 0, 0], [2, 2, 2], [0, 1, 2]):
        accumulator = _fresh()
        _offer_held(accumulator, n=3)
        before = dict(accumulator.tally.gap_resolution)
        counters = accumulator.drain(choices)
        censuses.append(dict(accumulator.tally.gap_resolution))
        assert counters["census_before"] == before
        assert counters["chose_genomic"] == sum(1 for c in choices if c == 2)

    assert censuses[0] == censuses[1] == censuses[2], (
        f"the arbitration census must not depend on what the drain chose; got {censuses}"
    )
    # The three deferred_* go to 0 with the bank; `gap_resolved_spliced` is pass one's and stays.
    assert censuses[0]["gap_deferred_rna_or_gdna"] == 0
    assert censuses[0]["gap_deferred_which_introns"] == 0
    assert censuses[0]["gap_deferred_both"] == 0


def test_a_choice_vector_of_the_wrong_length_is_refused():
    """⚠ The choices index the canonical queue. A length mismatch means the producer of the scores and the
    consumer of the draw walked different queues, which is a wrong answer that looks entirely plausible."""
    accumulator = _fresh()
    _offer_held(accumulator, n=3)
    with pytest.raises(ValueError, match="choices"):
        accumulator.drain([0, 1])


# ── §5.1 · THE DRAW, and §5.2 · reproducibility ────────────────────────────────────────────────────


def test_the_draw_respects_the_scores_and_is_reproducible():
    """⭐ §5.1: one multinomial draw per record. §5.2: one seed, one stream, in canonical order.

    A degenerate score vector is the sharp test — a hypothesis with all the mass must always be chosen,
    and one with none never. ⚠ That is checkable without a threshold, unlike "roughly 60/40".
    """
    from rigel.second_pass import choose_hypotheses

    n, k = 200, 3
    offsets = np.arange(n + 1, dtype=np.int64) * k
    score = np.zeros(n * k)
    score[1::k] = 1.0  # all the mass on local index 1, in every record

    payload = _FakeDeferred(offsets, n)
    picked = choose_hypotheses(_FakeScores(score), payload, seed=7)
    assert np.array_equal(picked, np.ones(n, np.int64)), (
        "a hypothesis holding all the posterior mass must always be the one drawn"
    )

    # Uniform scores: the draw must use all three, and repeat exactly at one seed.
    uniform = _FakeScores(np.full(n * k, 1.0 / k))
    first = choose_hypotheses(uniform, payload, seed=7)
    again = choose_hypotheses(uniform, payload, seed=7)
    other = choose_hypotheses(uniform, payload, seed=8)
    assert np.array_equal(first, again), "one seed must give one answer"
    assert set(first.tolist()) == {0, 1, 2}, (
        f"the draw must reach every hypothesis; got {set(first)}"
    )
    assert not np.array_equal(first, other), "a different seed must give a different draw"
    assert first.min() >= 0 and first.max() < k, "every choice must be inside its own record's run"


def test_the_draw_never_leaks_into_the_NEXT_records_run():
    """⛔ A draw that overshoots its run's cumulative total lands on the NEXT run's first slot — another
    fragment's hypothesis, silently, and every downstream number stays plausible.

    ⚠ **The adversarial case is CONSTRUCTED, not sampled.** With properly normalised scores the overshoot
    needs a uniform within one ULP of the run total, which random draws will never produce — so a gate that
    just drew a lot of uniforms would pass whether or not the guard existed (measured: it does). Feeding
    scores that sum to **0.5** per run makes every draw above 0.5 overshoot, which is half of them.

    ⭐ So the property gated here is the function's contract rather than one arithmetic accident: a choice
    is always inside its own record's run, whatever the caller's normalisation.
    """
    from rigel.second_pass import choose_hypotheses

    n, k = 400, 3
    offsets = np.arange(n + 1, dtype=np.int64) * k
    undersummed = _FakeScores(np.full(n * k, 0.5 / k))
    picked = choose_hypotheses(undersummed, _FakeDeferred(offsets, n), seed=3)
    assert picked.min() >= 0 and picked.max() < k, (
        f"a choice escaped its run: got range [{picked.min()}, {picked.max()}] for runs of {k}. Scores "
        f"summing to less than 1 must still never index another fragment's hypotheses."
    )

    # And the degenerate run: one hypothesis can only ever be index 0, however the draw falls.
    single = np.arange(n + 1, dtype=np.int64)
    assert np.array_equal(
        choose_hypotheses(_FakeScores(np.ones(n)), _FakeDeferred(single, n), seed=3),
        np.zeros(n, np.int64),
    )


@dataclasses.dataclass
class _FakeScores:
    """Just the one field :func:`choose_hypotheses` reads. ⚠ Deliberately not a real `HeldScores`: the
    draw must depend on the scores and the offsets and on nothing else, and this makes that checkable."""

    score: np.ndarray


class _FakeDeferred:
    """The two payload attributes the draw reaches: ``deferred.hypothesis_offsets`` and ``n_fragments``."""

    def __init__(self, offsets: np.ndarray, n: int) -> None:
        self.deferred = dataclasses.make_dataclass("_D", ["hypothesis_offsets", "n_fragments"])(
            offsets, n
        )


# ── non-vacuity ────────────────────────────────────────────────────────────────────────────────────


def test_the_fixture_reaches_all_three_hypotheses_and_they_differ():
    """⚠ If the three hypotheses gave the same ``L`` and touched the same objects, every gate above would
    pass over a distinction that does not exist."""
    lengths = set()
    for choice in range(3):
        accumulator = _fresh()
        _offer_held(accumulator)
        accumulator.drain([choice])
        nz = np.nonzero(accumulator.tally.deposited_lengths)[0]
        assert nz.size == 1
        lengths.add(int(nz[0]))
    assert lengths == {200, 400, 600}, (
        f"the three hypotheses must imply three different L (600 span, minus a 400 bp intron, a 200 bp "
        f"one, or nothing); got {sorted(lengths)}"
    )


def test_every_additive_channel_is_named_in_ADDITIVE_AXES():
    """⛔ The drain adds only what ``ADDITIVE_AXES`` names. A channel the accumulator fills but the table
    omits would go quietly short by exactly the drained fragments — and `pool_lengths` has no
    externally-checkable sum, so nothing else would notice."""
    named = {name for name, _axis in ADDITIVE_AXES}
    array_fields = {
        f.name
        for f in dataclasses.fields(Tally)
        if isinstance(getattr(_fresh().tally, f.name), np.ndarray)
    }
    assert array_fields <= named, (
        f"these accumulator array channels are not in ADDITIVE_AXES and so would not be drained: "
        f"{sorted(array_fields - named)}"
    )


# ══ THE PAYLOAD PATH — score, draw, drain, on a real scan across TWO contigs ════════════════════════
#
# ⛔ TWO CONTIGS ON PURPOSE. S1 found that a constant `ref` stamp and a constant `AccumulatorSet` id each
# passed all 1860 tests, because every fixture was single-contig or deferred only on reference 0. The drain
# rebuilds one accumulator per reference from the payload and slices the sj CSR per reference — so
# reference 0's slot base is 0 and every arithmetic error there is invisible. chr2's is not.

DRAIN_GENOME = 6_000
_M, _N = 0, 3


def _drain_gtf() -> str:
    rows = []
    for ref in ("chr1", "chr2"):
        rows += [
            f'{ref}\tt\texon\t1001\t1600\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
            f'{ref}\tt\texon\t2001\t2200\t.\t+\t.\tgene_id "gW_{ref}"; transcript_id "tW_{ref}";\n',
            f'{ref}\tt\texon\t1001\t1700\t.\t+\t.\tgene_id "gN_{ref}"; transcript_id "tN_{ref}";\n',
            f'{ref}\tt\texon\t1901\t2200\t.\t+\t.\tgene_id "gN_{ref}"; transcript_id "tN_{ref}";\n',
        ]
    return "".join(rows)


@pytest.fixture(scope="module")
def scanned_two_contigs(tmp_path_factory):
    """Scan a two-contig fixture; return everything the payload-path gates read."""
    import pysam

    from rigel.calibration.splice_graph import (
        build_sj_arrays,
        build_region_partition_arrays,
    )
    from rigel.config import BamScanConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import scan_and_buffer

    base = tmp_path_factory.mktemp("drain_two_contigs")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(
        "".join(
            f">{ref}\n" + "\n".join(["N" * 80] * (DRAIN_GENOME // 80)) + "\n"
            for ref in ("chr1", "chr2")
        )
    )
    pysam.faidx(str(fasta))
    gtf.write_text(_drain_gtf())
    index = TranscriptIndex.load(
        str(
            TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
            or base / "idx"
        )
    )

    def read(qname, ref_id, pos, cigar, mate_pos, is_r1):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = ref_id
        a.reference_start = pos
        a.mapping_quality = 60
        a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10)
        a.cigar = cigar
        n = sum(length for op, length in cigar if op == _M)
        a.query_sequence = "A" * n
        a.query_qualities = pysam.qualitystring_to_array("I" * n)
        a.next_reference_id = ref_id
        a.next_reference_start = mate_pos
        a.set_tags([("NH", 1, "i")])
        return a

    reads = []
    for ref_id in (0, 1):
        # depth on both sj, so neither reads zero and the draw has something to weigh
        for k in range(12):
            for donor, acceptor in ((1600, 2000), (1700, 1900)):
                reads += [
                    read(
                        f"sp_{ref_id}_{donor}_{k}",
                        ref_id,
                        donor - 200,
                        [(_M, 200), (_N, acceptor - donor), (_M, 100)],
                        acceptor + 100,
                        True,
                    ),
                    read(
                        f"sp_{ref_id}_{donor}_{k}",
                        ref_id,
                        acceptor + 100,
                        [(_M, 100)],
                        donor - 200,
                        False,
                    ),
                ]
        # the ambiguous fragments — three per contig, so both banks are non-empty
        for k in range(3):
            reads += [
                read(f"amb_{ref_id}_{k}", ref_id, 1500, [(_M, 100)], 2000, True),
                read(f"amb_{ref_id}_{k}", ref_id, 2000, [(_M, 100)], 1500, False),
            ]

    bam = str(base / "drain.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": DRAIN_GENOME}, {"SN": "chr2", "LN": DRAIN_GENOME}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)

    _s, strand_model, _b, payload = scan_and_buffer(
        bam, index, BamScanConfig(sj_strand_tag="XS", total_threads=1)
    )
    _, _, region_types = build_region_partition_arrays(index)
    return payload, region_types, build_sj_arrays(index), strand_model


def _drained(scanned_two_contigs, seed=11):
    from rigel.calibration.fl import build_fl_models
    from rigel.second_pass import choose_hypotheses, drain, score_held_fragments

    payload, region_types, sj, strand_model = scanned_two_contigs
    scores = score_held_fragments(
        payload,
        fl_models=build_fl_models(payload),
        rna_sense_frac=strand_model.p_r1_sense,
        region_types=region_types,
        sj=sj,
    )
    choices = choose_hypotheses(scores, payload, seed=seed)
    return payload, drain(payload, choices, region_types=region_types, sj=sj), choices


def test_the_fixture_holds_fragments_on_BOTH_contigs(scanned_two_contigs):
    """⚠ Non-vacuity, and it is the whole reason this fixture exists. If everything deferred on reference 0,
    every per-reference slice below would be exercised only where its base offset is 0."""
    payload = scanned_two_contigs[0]
    refs = set(payload.deferred.ref.tolist())
    assert refs == {0, 1}, f"both contigs must hold fragments; got {refs}"
    assert payload.deferred.n_fragments == 6


def test_the_drained_payload_CONSERVES_and_the_bank_is_empty(scanned_two_contigs):
    """⭐ §6.2, at the payload. The identity, the emptied bank, and both externally-checkable sums."""
    before, after, _choices = _drained(scanned_two_contigs)

    assert after.drain is not None and after.drain.conserved
    assert after.drain.offered == before.qc.deferred_undetermined_gap
    assert after.deferred.n_fragments == 0
    assert after.qc.deferred_undetermined_gap == 0
    assert after.qc.deposited == before.qc.deposited + after.drain.deposited
    # ⭐ The two invariants that survive the per-reference placement arithmetic in `_gather_delta`.
    assert int(after.region_start_count.sum()) == after.qc.deposited
    assert int(after.deposited_lengths.sum()) == after.qc.deposited
    # ⛔ Pass one's payload must be untouched — the delta is a separate object, which is what makes the
    # drain's contribution to every channel a subtraction rather than a rerun (§6.3).
    assert before.deferred.n_fragments == 6
    assert before.drain is None


def test_the_drain_credits_the_sj_on_the_RIGHT_CONTIG(scanned_two_contigs):
    """⛔ The gate for the per-reference sj slice. A drained spliced choice must credit a sj slot
    **on its own contig**, and chr2's slots are not zero-based — which is what makes the arithmetic visible.

    ⚠ Two failures this catches, both silent: sj never installed (every slot stays at pass one's
    value, because an observed intron with no table reads as unannotated), and a slot base taken as 0 for
    every reference (chr2's traffic lands on chr1's sj).
    """
    before, after, choices = _drained(scanned_two_contigs)
    delta = after.sj_count.astype(np.int64) - before.sj_count.astype(np.int64)
    n_spliced = after.drain.chose_spliced
    assert n_spliced > 0, "the fixture must have drawn at least one spliced hypothesis to gate this"
    assert int(delta.sum()) == n_spliced, (
        f"every spliced choice credits exactly one sj slot; {n_spliced} were chosen but the sj "
        f"banks moved by {int(delta.sum())}. Zero means the sj table was never installed."
    )
    # Which contig did the movement land on? Slots are reference-major, so compare against ref_sj_offsets.
    moved = np.flatnonzero(delta.sum(axis=1))
    chr2_slot_base = int(after.ref_sj_offsets[1])
    assert chr2_slot_base > 0, (
        "the fixture must give chr2 a non-zero slot base or this gate is vacuous"
    )
    held_refs = before.deferred.ref
    spliced_on_chr2 = sum(
        1
        for i in range(before.deferred.n_fragments)
        if int(held_refs[i]) == 1
        and len(
            before.deferred.hypothesis_introns_of(
                int(before.deferred.hypothesis_offsets[i]) + int(choices[i])
            )
        )
    )
    if spliced_on_chr2:
        assert np.any(moved >= chr2_slot_base), (
            f"{spliced_on_chr2} chr2 fragments chose a spliced path, but every credited slot is below "
            f"chr2's base {chr2_slot_base} — chr2's sj traffic landed on chr1's slots."
        )


def test_the_drain_is_REPRODUCIBLE_and_the_seed_is_what_moves_it(scanned_two_contigs):
    """⭐ §5.2. Same payload, same seed → byte-identical on every channel. ⚠ The drain never sees a thread —
    it runs after the worker merge, over a canonically sorted bank — so this is structural; the gate exists
    to keep it that way."""
    _b, first, choices_a = _drained(scanned_two_contigs, seed=11)
    _b, again, choices_b = _drained(scanned_two_contigs, seed=11)
    assert np.array_equal(choices_a, choices_b)
    for name, _axis in ADDITIVE_AXES:
        assert np.array_equal(getattr(first, name), getattr(again, name)), (
            f"{name} differs between two drains at one seed"
        )
    assert first.drain == again.drain


def test_a_payload_cannot_be_drained_TWICE(scanned_two_contigs):
    """⛔ The drain consumes the bank, so a second one would deposit nothing and double the bookkeeping."""
    from rigel.second_pass import drain

    _before, after, choices = _drained(scanned_two_contigs)
    payload, region_types, sj, _strand = scanned_two_contigs
    with pytest.raises(ValueError, match="already been drained"):
        drain(after, np.zeros(0, np.int64), region_types=region_types, sj=sj)
    assert len(choices) == 6


def test_the_DRAINED_payload_is_byte_identical_at_1_2_4_8_WORKERS(tmp_path_factory):
    """⭐ **§5.2's gate, composed through the drain** — the one P3 names.

    S1 and S2.1 already establish that the *bank* is worker-independent, and the drain runs after the
    worker merge over that bank's canonical order, so this holds structurally rather than by luck. ⚠ The
    gate exists because "structurally" is a property of the current shape: the moment anything drains
    per-worker, or consumes the queue in append order, this is the test that says so.

    ⛔ Compared on every additive channel AND on the choice vector, because identical tallies from
    different choices would mean the draw had stopped depending on the queue order.

    ⚠ **Where the teeth actually are, stated honestly**: the canonical sort itself is gated upstream —
    S1's perturbation X1 (the sort does nothing) fails 5 tests in the parity and worker-determinism
    modules. This gate adds the composition: that scoring, drawing and draining on top of that bank
    introduces no new order dependence of their own.
    """
    import pysam

    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.splice_graph import (
        build_sj_arrays,
        build_region_partition_arrays,
    )
    from rigel.config import BamScanConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import scan_and_buffer
    from rigel.second_pass import choose_hypotheses, drain, score_held_fragments

    base = tmp_path_factory.mktemp("drain_workers")
    fasta, gtf = base / "g.fa", base / "a.gtf"
    fasta.write_text(
        "".join(
            f">{ref}\n" + "\n".join(["N" * 80] * (DRAIN_GENOME // 80)) + "\n"
            for ref in ("chr1", "chr2")
        )
    )
    pysam.faidx(str(fasta))
    gtf.write_text(_drain_gtf())
    TranscriptIndex.build(str(fasta), str(gtf), str(base / "idx"), write_tsv=False)
    index = TranscriptIndex.load(str(base / "idx"))

    def read(qname, ref_id, pos, cigar, mate_pos, is_r1):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.reference_id = ref_id
        a.reference_start = pos
        a.mapping_quality = 60
        a.flag = 0x1 | 0x2 | (0x40 | 0x20 if is_r1 else 0x80 | 0x10)
        a.cigar = cigar
        n = sum(length for op, length in cigar if op == _M)
        a.query_sequence = "A" * n
        a.query_qualities = pysam.qualitystring_to_array("I" * n)
        a.next_reference_id = ref_id
        a.next_reference_start = mate_pos
        a.set_tags([("NH", 1, "i")])
        return a

    # ⚠ Enough fragments that several chunks and several workers really are in play.
    reads = []
    for ref_id in (0, 1):
        for k in range(60):
            for donor, acceptor in ((1600, 2000), (1700, 1900)):
                reads += [
                    read(
                        f"sp_{ref_id}_{donor}_{k}",
                        ref_id,
                        donor - 200,
                        [(_M, 200), (_N, acceptor - donor), (_M, 100)],
                        acceptor + 100,
                        True,
                    ),
                    read(
                        f"sp_{ref_id}_{donor}_{k}",
                        ref_id,
                        acceptor + 100,
                        [(_M, 100)],
                        donor - 200,
                        False,
                    ),
                ]
        for k in range(40):
            reads += [
                read(f"amb_{ref_id}_{k}", ref_id, 1500, [(_M, 100)], 2000, True),
                read(f"amb_{ref_id}_{k}", ref_id, 2000, [(_M, 100)], 1500, False),
            ]
    bam = str(base / "workers.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "queryname"},
        "SQ": [{"SN": "chr1", "LN": DRAIN_GENOME}, {"SN": "chr2", "LN": DRAIN_GENOME}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.sort("-n", "-o", bam, bam)

    _, _, region_types = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)

    results = {}
    for threads in (1, 2, 4, 8):
        _s, strand_model, _b, payload = scan_and_buffer(
            bam,
            index,
            # ⚠ A SMALL chunk so several chunks and several workers really are in flight. At the default
            # 1 M the whole fixture is one chunk and the gate would pass whatever the worker count said.
            BamScanConfig(sj_strand_tag="XS", total_threads=threads, fragments_per_chunk=32),
        )
        scores = score_held_fragments(
            payload,
            fl_models=build_fl_models(payload),
            rna_sense_frac=strand_model.p_r1_sense,
            region_types=region_types,
            sj=sj,
        )
        choices = choose_hypotheses(scores, payload, seed=5)
        results[threads] = (
            choices,
            drain(payload, choices, region_types=region_types, sj=sj),
        )

    assert results[1][1].drain.offered == 80, (
        f"the fixture must hold 80 fragments across two contigs; got {results[1][1].drain.offered}"
    )
    reference_choices, reference_payload = results[1]
    for threads in (2, 4, 8):
        choices, drained = results[threads]
        assert np.array_equal(choices, reference_choices), (
            f"the draw differs at {threads} workers; the deferred bank's canonical order is what makes "
            f"'the i-th record' mean one thing, so this is that order having become observable."
        )
        for name, _axis in ADDITIVE_AXES:
            assert np.array_equal(getattr(drained, name), getattr(reference_payload, name)), (
                f"{name} differs between 1 and {threads} workers after the drain"
            )
        assert drained.drain == reference_payload.drain
        assert drained.qc == reference_payload.qc


# ══ THE COMBINATION RULE — a factor that is zero for EVERY candidate must not annihilate the others ══
#
# ⛔ THE BUG THIS GATES, found 2026-08-03 by scoring the pilot against the simulator's per-fragment truth.
# The score is a PRODUCT of three factors, normalised within each fragment's candidate set. If one factor
# is zero for *every* candidate the product is zero everywhere, the normalisation cannot run, and the whole
# record falls back to a uniform coin toss — discarding the other two factors, which usually DO decide.
#
# Measured on `gdna_none_ss_0.99_capture_off`: of the 10 records that fell back to uniform, the length term
# alone would have picked the CORRECT candidate 8 times out of the 8 it could decide — 100 %. The coin got
# them right by chance instead. Both of P4's over-ceiling fragments came from exactly these records, where
# the length term had already scored the impossible answer at zero.
#
# ⭐ THE RULE, AND WHY IT NEEDS NO CONSTANT. The score is normalised within the candidate set, so a factor
# that takes the SAME value for every candidate cancels and cannot affect the answer. Zero is the one value
# where the arithmetic loses that property: instead of cancelling it destroys the product. So an all-zero
# factor is treated as what it is — **uninformative** — and dropped from the product for that fragment.
#
# ⛔ IT DOES NOT TOUCH THE PARTIAL-ZERO CASE, and that is the point. A factor that is zero for SOME
# candidates and positive for others is highly informative: the zero says "no evidence for this path". That
# stays decisive, which is the owner's D-3 ruling ("no fallback") left exactly as it was.


def _score_one(rho, f, s):
    """Score ONE synthetic record with the given per-candidate factors. Returns the normalised scores.

    ⚠ Drives the real ``score_held_fragments`` combination step through a payload whose factors are forced,
    rather than reimplementing the arithmetic — a second copy of the rule would be a second rule.
    """
    from rigel.second_pass import combine_factors

    return combine_factors(
        np.asarray(rho, np.float64), np.asarray(f, np.float64), np.asarray(s, np.float64)
    )


def test_an_ALL_ZERO_factor_is_UNINFORMATIVE_and_does_not_annihilate_the_others():
    """⭐ The bug. Every candidate has zero local traffic, but the length term separates them cleanly."""
    scores, undecided = _score_one(rho=[0.0, 0.0], f=[0.0001, 0.0], s=[0.99, 1.0])
    assert not undecided, "the length term decides this record; it is not undecided"
    assert scores[0] > scores[1], (
        f"the candidate with 100x the length plausibility must win; got {scores}. A factor that is zero "
        f"for EVERY candidate carries no information and must drop out, not zero the product."
    )
    assert scores[1] == 0.0, "and the impossible length keeps zero posterior"


def test_a_PARTIAL_zero_stays_DECISIVE():
    """⛔ The other half, and it must not change. A zero for SOME candidates is evidence, not an absence of
    evidence — it says this path has no support. The owner's D-3 ruling keeps it decisive."""
    scores, undecided = _score_one(rho=[1.0, 0.0], f=[0.5, 0.5], s=[1.0, 1.0])
    assert not undecided
    assert scores[0] == 1.0 and scores[1] == 0.0, (
        f"a partial zero must remain decisive; got {scores}. Softening it would reopen D-3, which the "
        f"owner ruled needs no fallback."
    )


def test_a_factor_that_is_CONSTANT_and_positive_already_cancels():
    """⚠ The premise the rule rests on, checked rather than assumed: normalisation makes any constant
    factor irrelevant, so dropping an all-zero one is consistent rather than special-cased."""
    with_it, _ = _score_one(rho=[2.0, 6.0], f=[0.3, 0.3], s=[1.0, 1.0])
    without, _ = _score_one(rho=[2.0, 6.0], f=[1.0, 1.0], s=[1.0, 1.0])
    assert np.allclose(with_it, without), (
        f"a constant positive factor must already cancel under normalisation; got {with_it} vs {without}"
    )


def test_a_WEAKER_factor_cannot_VETO_what_a_STRONGER_one_left_standing():
    """⭐ Factors are applied strongest-evidence first, and a weaker one only refines within the survivors.

    Traffic favours candidate 0 and the strand term candidate 1, with the length term allowing both. Under
    a blind product these zero each other out and the record becomes a coin toss. Under the ordered rule the
    strand term — which rests on the library's whole spliced population — narrows to candidate 1, and
    traffic is then flat-zero **among the survivors**, so it is uninformative and says nothing.

    ⛔ This is what removes the "irreducible contradiction" case entirely: each factor either narrows a
    non-empty set or is skipped, so the product can never collapse to zero everywhere.
    """
    scores, undecided = _score_one(rho=[1.0, 0.0], f=[0.5, 0.5], s=[0.0, 1.0])
    assert not undecided
    assert np.allclose(scores, [0.0, 1.0]), (
        f"the strand term excluded candidate 0; traffic is zero for the one survivor and so is "
        f"uninformative, not a veto. got {scores}"
    )


def test_ALL_THREE_factors_absent_is_a_genuine_coin_toss():
    """⚠ And the floor: no evidence of any kind means uniform, which is what it meant before."""
    scores, undecided = _score_one(rho=[0.0, 0.0, 0.0], f=[0.0, 0.0, 0.0], s=[0.0, 0.0, 0.0])
    assert undecided
    assert np.allclose(scores, [1 / 3, 1 / 3, 1 / 3])


def test_an_IMPOSSIBLE_LENGTH_cannot_be_CHOSEN_however_much_traffic_favours_it():
    """⭐ **The owner's question, gated.** *"A fragment length of 739 should be exceedingly unlikely under
    the first-pass RNA length distribution — essentially nonexistent. How could traffic possibly overcome
    that?"* It cannot, and this is where that is enforced.

    ⛔ **No cutoff and no constant.** ``f = 0`` already means *"no fragment of this length was observed
    anywhere in the library"* — the same statement ``max_fragment_length`` makes, read off the measured
    distribution instead of a round number. Measured on the pilot, the RNA pmf's support ends at **713 bp**,
    which is exactly the library's true longest molecule, so a candidate at 739 is excluded by evidence
    rather than by a threshold anybody chose.

    ⚠ The case is real: pilot record 155262 had traffic of 11.25 behind candidates at L = 739 and 1024 and
    the length term behind the true one at L = 352. Traffic three orders of magnitude larger does not buy a
    molecule the library does not contain.
    """
    from rigel.second_pass import combine_factors

    # Traffic likes candidates 0 and 1; the length term rules both out and likes candidate 2.
    scores, undecided = combine_factors(
        np.array([11.25, 11.25, 0.0]),
        np.array([0.0, 0.0, 0.0014]),
        np.array([1.0, 1.0, 1.0]),
    )
    assert not undecided, (
        "the length term excludes both impossible candidates and traffic is then flat-zero among the one "
        "survivor, so it is uninformative — this is decidable, not a coin toss"
    )
    assert scores[2] == 1.0 and scores[0] == 0.0 and scores[1] == 0.0, (
        f"a length the library does not contain cannot be chosen however much traffic favours it; "
        f"got {scores}"
    )

    # ⚠ And when the length term rules out EVERYTHING it cannot narrow anything: uniform over all.
    scores, undecided = combine_factors(
        np.array([0.0, 0.0]), np.array([0.0, 0.0]), np.array([1.0, 1.0])
    )
    assert undecided and np.allclose(scores, [0.5, 0.5])


def test_SURVIVORS_are_weighted_by_LIKELIHOOD_not_tossed_for():
    """⛔ **The pilot's record 155262, and the reason "restrict then toss a fair coin" is not good enough.**

    The length term leaves two candidates possible — one at ``9.9e-06`` and the true one at ``1.4e-03``,
    a **143-fold** difference — and traffic is zero for both, so it is uninformative among them. The answer
    must be weighted by the surviving factor, not drawn uniformly from its support.

    ⚠ An earlier fix restricted the fallback to possible lengths and then flipped a fair coin, which picked
    the 143× less likely candidate half the time. Narrowing the draw and *weighting* it are different
    things, and only the second is using the likelihood.
    """
    scores, undecided = _score_one(
        rho=[11.247, 11.247, 0.0, 0.0, 0.010],
        f=[0.0, 0.0, 9.9236e-06, 1.4225e-03, 0.0],
        s=[1.0, 1.0, 1.0, 1.0, 1.0],
    )
    assert not undecided
    assert scores.argmax() == 3, f"the 143x more likely candidate must win; got {scores}"
    assert scores[3] / scores[2] == pytest.approx(1.4225e-03 / 9.9236e-06, rel=1e-9), (
        f"and it must win BY THAT RATIO, not by a coin toss over the two survivors; got {scores}"
    )
    assert scores[0] == 0.0 and scores[1] == 0.0 and scores[4] == 0.0


# ══ THE ∅ STRAND TERM — gDNA is biologically 50/50 ════════════════════════════════════════════════
#
# ⛔ THE GAP THIS CLOSES, owner-derived 2026-08-03. The genomic candidate carried `s = 1.0` in both
# orientations, which is not a probability at all — under a hypothesis the two orientations must sum to 1,
# and 1.0 twice sums to 2.
#
# ⭐ THE DERIVATION. Let `t` be the strand a candidate implies, `a` the strand the fragment aligned to, and
# `p = P(a == t | RNA)` the library's directional sense fraction (≈ 0.01 on dUTP).
#
#     H_spliced   used a sj; gDNA cannot splice, so RNA on strand t:   p  or  1 - p
#     H_genomic   crossed contiguously; the discriminating component is gDNA, which is DOUBLE-STRANDED
#                 and therefore has no sense direction at all:              0.5, either orientation
#
# ⭐ ∅ also covers unspliced RNA — but unspliced RNA would give the same p / (1-p) as the spliced candidate
# and cancel, contributing nothing. The only part of ∅ that can separate it is its gDNA component, whose
# orientation likelihood is a biological constant rather than anything to be fitted.
#
# ⛔ AND A GLOBAL MIXTURE MARGINAL IS WORSE THAN USELESS, which the algebra shows and measurement confirmed.
# The orientation discrimination with ANY constant c for ∅ is [c/p] / [c/(1-p)] = (1-p)/p = 98.0 — c cancels
# — so 0.5 and 1.0 discriminate identically and 0.5 is simply the one that is a probability. But an
# orientation-DEPENDENT ∅ term q/(1-q) from the library-wide genic marginal (q = 0.1825 measured on a
# gdna100 pilot) gives q(1-p)/p(1-q) = 21.9: it moves ∅ in the SAME direction as the spliced term and
# destroys 78 % of the signal. A global value says nothing about one fragment, whose gene may be silent
# (pure gDNA) or highly expressed.


def test_gDNA_is_STRAND_SYMMETRIC_in_both_orientations():
    """⭐ The biological fact the ∅ term rests on: double-stranded DNA has no sense direction."""
    from rigel.second_pass import strand_terms

    for align in (int(Strand.POS), int(Strand.NEG)):
        _spliced, genomic = strand_terms(
            align=align, implied_strand=int(Strand.POS), rna_sense_frac=0.0101
        )
        assert genomic == 0.5, (
            f"gDNA must be 50/50 in either orientation; got {genomic}. A flat 1.0 is not a probability — "
            f"the two orientations would sum to 2."
        )


def test_the_RNA_UNLIKELY_orientation_favours_the_GENOMIC_candidate():
    """⭐ The owner's model, and note the direction is protocol-dependent so nothing here hard-codes it.

    On a dUTP library ``p = 0.01``, so a fragment aligned CO-oriented with the candidate's strand is the
    RNA-unlikely case: RNA would produce that orientation only 1 % of the time while DNA produces it half
    the time, so ∅ is favoured 50:1. Counter-oriented, RNA is favoured — but only about 2:1.
    """
    from rigel.second_pass import strand_terms

    p = 0.0101
    co_spliced, co_genomic = strand_terms(
        align=int(Strand.POS), implied_strand=int(Strand.POS), rna_sense_frac=p
    )
    assert co_spliced == pytest.approx(p) and co_genomic == 0.5
    assert co_genomic / co_spliced == pytest.approx(0.5 / p), (
        "the RNA-unlikely orientation must favour ∅ by 0.5/p"
    )

    ct_spliced, ct_genomic = strand_terms(
        align=int(Strand.NEG), implied_strand=int(Strand.POS), rna_sense_frac=p
    )
    assert ct_spliced == pytest.approx(1.0 - p) and ct_genomic == 0.5
    assert ct_spliced > ct_genomic, "the RNA-likely orientation must favour RNA, and only mildly"
    assert ct_spliced / ct_genomic == pytest.approx((1.0 - p) / 0.5)


def test_the_ORIENTATION_DISCRIMINATION_is_exactly_the_LIBRARY_SPECIFICITY_ODDS():
    """⛔ **The regression gate against reintroducing a fitted marginal.**

    The whole information content of this term is the odds ratio between the two orientations, and it must
    equal ``(1 - p) / p`` — the library's own strand-specificity odds — exactly. ⭐ That is what makes ∅'s
    constant irrelevant to the discrimination and what a global mixture marginal destroys: measured,
    ``q = 0.1825`` collapses 98.0 to 21.9.
    """
    from rigel.second_pass import strand_terms

    p = 0.0101
    co_s, co_g = strand_terms(
        align=int(Strand.POS), implied_strand=int(Strand.POS), rna_sense_frac=p
    )
    ct_s, ct_g = strand_terms(
        align=int(Strand.NEG), implied_strand=int(Strand.POS), rna_sense_frac=p
    )
    discrimination = (co_g / co_s) / (ct_g / ct_s)
    assert discrimination == pytest.approx((1.0 - p) / p), (
        f"the orientation discrimination must be the library specificity odds {(1 - p) / p:.2f}; got "
        f"{discrimination:.2f}. A value below it means ∅'s term has become orientation-dependent — which "
        f"is what a fitted mixture marginal does, and it cost 78 % of the signal when measured."
    )


def test_an_UNSTRANDED_library_gives_a_NEUTRAL_strand_term():
    """⚠ At ``p = 0.5`` there is no strand information, so every candidate scores 0.5 and the factor
    cancels — which `combine_factors` then skips as uninformative."""
    from rigel.second_pass import strand_terms

    for align in (int(Strand.POS), int(Strand.NEG)):
        spliced, genomic = strand_terms(
            align=align, implied_strand=int(Strand.POS), rna_sense_frac=0.5
        )
        assert spliced == pytest.approx(0.5) and genomic == pytest.approx(0.5)


def test_an_AMBIGUOUS_implied_strand_says_NOTHING():
    """⚠ D-5's case: one path claimed by both strands names no orientation to compare against, so the
    spliced candidate cannot be scored on strand either and must fall back to symmetric."""
    from rigel.second_pass import strand_terms

    spliced, genomic = strand_terms(
        align=int(Strand.POS), implied_strand=int(Strand.AMBIGUOUS), rna_sense_frac=0.0101
    )
    assert spliced == 0.5 and genomic == 0.5
