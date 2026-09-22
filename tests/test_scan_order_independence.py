"""One BAM gives one answer whatever the scan's thread count — for the fragment buffer the EM reads
and for the accumulator's tally alike. The first block gates the pipeline. The scanner streams
finalized chunks from its worker threads in completion order, so the row a fragment occupies depends
on which worker got there first; the EM's ``assignment_mode="sample"`` draws one categorical sample
per unit from a per-locus RNG, so a permutation of the units permutes which fragment gets which draw,
and units in different equivalence classes have different posteriors. The remedy is an ordering
rather than a seed: ``frag_id`` is assigned by the single reader thread in BAM order and is a stable
identity, so a locus's units are sorted by it once, where ``MultiLocus.unit_indices`` is built, and
every per-locus array scattered by that index list inherits the canonical order. The second block
gates the tally. Its uint32 count banks merge exactly whatever order the chunks arrived in, because
integer addition is associative; its float64 fraction and mass banks agree to within a tolerance
derived from the deposit count, because the per-worker merge re-associates their sums.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline, scan_and_buffer
from rigel.scan_payload import AccumulatorPayload
from rigel.sim import ReadSimConfig, Scenario


# ── The FRAGMENT BUFFER — the pipeline's answer at every scan thread count ──────────────────────
#
# The remedy is an ordering by identity rather than a per-fragment seed, and the reason matters.
# Deriving the draw from a hash of the fragment's own content would also be
# order-independent, and would be wrong: identical fragments would hash identically and so all draw
# identically, turning a 60/40 posterior into 100/0 for every group of duplicates. Ordering by
# identity and keeping one stream per locus preserves the multinomial spread, which is the whole
# point of the sampling mode.


SEED = 7

#: Thread counts spanning single-chunk and many-chunk scans. A low pair alone passes on a small BAM
#: even without the ordering: the gate needs a count high enough that the scanner really does emit
#: several chunks concurrently, which is what makes the row order vary.
THREAD_COUNTS = (1, 4, 16)


@pytest.fixture(scope="module")
def oracle(tmp_path_factory):
    """A locus with two isoforms sharing an exon, so units land in DIFFERENT equivalence classes.

    That is what gives the gate teeth. Fragments that are all interchangeable would be permuted with no
    effect on the answer — the assignment only moves when units with *different* posteriors swap draws.
    ``t1``/``t2`` share their first exon and differ in the second, so a fragment in the shared region is
    ambiguous while one in the unique region is not.
    """
    work = tmp_path_factory.mktemp("scan_order")
    scenario = Scenario("scan_order", genome_length=8000, seed=SEED, work_dir=work / "sim")
    scenario.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(1000, 1400), (2000, 2400)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1400), (2200, 2400)], "abundance": 60},
        ],
    )
    scenario.add_gene("g2", "-", [{"t_id": "t3", "exons": [(4000, 4600)], "abundance": 40}])
    result = scenario.build_oracle(
        n_fragments=3000,
        gdna_fraction=0.2,
        sim_config=ReadSimConfig(
            frag_mean=220,
            frag_std=50,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=0.99,
            seed=SEED,
        ),
    )
    yield result
    scenario.cleanup()


def _counts(oracle, threads: int, mode: str = "sample") -> np.ndarray:
    config = PipelineConfig(
        em=EMConfig(seed=1234, assignment_mode=mode),
        scan=BamScanConfig(sj_strand_tag="auto", total_threads=threads),
    )
    estimator = run_pipeline(oracle.bam_path, oracle.index, config=config).estimator
    return np.asarray(estimator.t_counts, dtype=np.float64)


#: Two repeats per thread count, not one. The permutation is random, so a single pair can agree by luck
#: — two runs at the SAME thread count can disagree while three different counts happen to match. Six
#: runs against one baseline is what makes the gate reliable rather than occasionally lucky.
REPEATS = 2


@pytest.mark.parametrize("mode", ["sample", "map", "fractional"])
def test_the_answer_IS_THE_SAME_AT_EVERY_SCAN_THREAD_COUNT(oracle, mode):
    """The contract, for all three assignment modes: one BAM and one seed give one answer.

    Byte-identical across thread counts AND across repeats, which is strictly stronger than
    "reproducible at the default" and the only version that survives someone changing the default.

    All three modes, because each fails a different way. ``sample`` moves by whole counts, since a
    permutation changes which fragment draws which sample. ``map`` and ``fractional`` are immune to
    that — they read each unit's posterior alone — but ``fractional`` scatters float posteriors into
    shared accumulators, so a permutation reorders the summation and the answer drifts by ULPs. One
    ordering fixes both; asserting only the first leaves the second to be rediscovered as "flaky".
    """
    runs = [_counts(oracle, threads, mode) for threads in THREAD_COUNTS for _ in range(REPEATS)]
    assert runs[0].sum() > 0, "nothing was assigned; this fixture cannot detect a permutation"
    for i, other in enumerate(runs[1:], start=1):
        assert np.array_equal(runs[0], other), (
            f"{mode}: run {i} of {len(runs)} disagrees with the first by up to "
            f"{np.abs(runs[0] - other).max()}. The answer must not depend on the order the scanner "
            f"happened to fill the fragment buffer in."
        )


def test_THE_FIXTURE_REALLY_DOES_REORDER_THE_BUFFER(oracle):
    """Non-vacuity, and it is not optional here.

    The contract above passes trivially if the scan hands the buffer back in the same row order every
    time — there is then no permutation to be immune to, and the gate reads as coverage it does not
    have.

    Retried, because the reordering is a RACE: a gate that asserts a race fired on the first try goes
    red for the wrong reason on a loaded machine.
    """
    from rigel.pipeline import scan_and_buffer

    def frag_id_order(threads: int) -> list[int]:
        _stats, _strand, buffer, _payload = scan_and_buffer(
            str(oracle.bam_path),
            oracle.index,
            BamScanConfig(sj_strand_tag="auto", total_threads=threads),
        )
        return [int(i) for chunk in buffer.iter_chunks() for i in chunk.frag_id]

    serial = frag_id_order(THREAD_COUNTS[0])
    parallel_orders = [frag_id_order(THREAD_COUNTS[-1]) for _ in range(3)]
    assert all(sorted(p) == sorted(serial) for p in parallel_orders), (
        "the scans saw different fragments, so this is not a pure reordering"
    )
    assert any(p != serial for p in parallel_orders), (
        f"three {THREAD_COUNTS[-1]}-thread scans all returned the serial row order, so this module is "
        f"not exercising the permutation it exists to be immune to — raise the fragment count or lower "
        f"BamScanConfig.fragments_per_chunk until it does"
    )


def test_frag_id_IS_AN_IDENTITY_which_is_what_makes_it_a_legal_sort_key(oracle):
    """The assumption the ordering rests on, pinned where the ordering is gated.

    Ordering by ``frag_id`` canonicalises the units ONLY IF no two units share one. If ids could repeat,
    ``argsort(kind="stable")`` would break the tie by the unit's buffer row, silently reintroducing the
    exact dependence the ordering removes, on exactly the duplicates it matters for.

    Uniqueness, not contiguity. The EM's units are a SUBSET of the scanned fragments, so the ids have
    gaps — which is fine, because the sort needs a total order rather than a dense one.
    """
    import numpy as np_

    from rigel import locus as locus_module

    captured: dict[str, np_.ndarray] = {}
    original = locus_module.build_multi_loci

    def capture(em_data, index):
        captured["frag_ids"] = np_.asarray(em_data.frag_ids).copy()
        return original(em_data, index)

    locus_module.build_multi_loci = capture
    try:
        _counts(oracle, threads=THREAD_COUNTS[-1])
    finally:
        locus_module.build_multi_loci = original

    frag_ids = captured["frag_ids"]
    assert frag_ids.size > 0, "no EM units were built, so this proves nothing"
    assert np_.unique(frag_ids).size == frag_ids.size, (
        f"{frag_ids.size - np_.unique(frag_ids).size} EM units share a frag_id. The canonical unit "
        f"order is a sort on that id, so duplicates fall back to buffer-row order and the pipeline "
        f"stops being reproducible for precisely the fragments where it matters most."
    )


# ── The SPLIT — how the budget is divided between decompression and workers ─────────────────────


def test_the_thread_split_is_derived_from_the_budget_and_an_explicit_one_still_wins():
    """The split is a RATIO, not a count: one decompression thread keeps about eight scan workers fed,
    so the budget is divided by that ratio rather than by a fixed reservation.

    Measured on the 18.6M-fragment library, scan seconds by budget and decompression threads — 4: (0)
    32.0, (1) 41.2; 8: (1) 25.8, (2) 26.2, (0) 28.3, (4) 34.5; 16: (2) 17.6, (1) 19.6, (4) 20.1 — so the
    derived cell is the best one at every budget measured, and the old fixed 4 was the worst at all
    three. What this gate holds is the ARITHMETIC: every thread is spent, none is spent twice, at least
    one worker always runs, and an explicit request still overrides. The answer itself cannot move with
    the split — that is the tally's own gate below, and the reason is that every bank is a sum of
    integers.
    """
    for total, expect in (
        (1, (1, 0)),
        (2, (2, 0)),
        (4, (4, 0)),
        (8, (7, 1)),
        (16, (14, 2)),
        (64, (56, 8)),
    ):
        got = BamScanConfig(total_threads=total).resolved_scan_threads()
        assert got == expect, f"budget {total} split {got}, expected {expect}"
        workers, bgzf = got
        assert workers + bgzf == total and workers >= 1
    # an explicit request is honoured, and still capped to leave a worker
    assert BamScanConfig(total_threads=8, bgzf_threads=4).resolved_scan_threads() == (4, 4)
    assert BamScanConfig(total_threads=8, bgzf_threads=0).resolved_scan_threads() == (8, 0)
    assert BamScanConfig(total_threads=2, bgzf_threads=9).resolved_scan_threads() == (1, 1)


# ── The TALLY — the same BAM gives the same accumulator payload at any worker count ─────────────
#
# The count banks are uint32 integers and integer addition is associative, so they merge bit for
# bit. The fraction and mass banks are float64 and the per-worker merge re-associates their sums,
# so they agree only to within the float64 representation. Testing ``Accumulator.merge_from``
# directly — which the parity module does — is not enough, because this exercises the SCANNER's
# worker path: each worker builds its own ``AccumulatorSet`` from the scanner's members, and the
# chunk-to-worker split is data-dependent. A worker whose set was constructed differently, with the
# sj installed on the template but not on the copies say, is invisible to a merge test on one
# accumulator and shows up here. The fixture is named for the property it has to have: every bank of
# the tally must receive something, or a bit-identity gate over an all-zero array passes for the
# wrong reason.


TALLY_SEED = 20260730


#: Derived from the machine, never fitted.
EPS = float(np.finfo(np.float64).eps)


@pytest.fixture
def every_bank_oracle(tmp_path):
    """A scenario built so that EVERY bank of the tally receives something.

    The obvious scenario does not. Two single-isoform genes deposit into the regions and the sj
    boundaries and leave both contiguous-boundary banks identically zero: every region boundary is an
    exon boundary, so a mature fragment either fits inside an exon (contained) or splices across the gap
    (sj boundary) and never has bases on both sides of a boundary. A bit-identity gate over an all-zero
    array passes for the wrong reason, and this project has already had one arm with zero rows report
    itself fully identical. So:

    * ``t2`` starts at 500, INSIDE ``t1``'s first exon. That makes 500 a region boundary, and a ``t1``
      fragment spanning it has bases on both sides — a contiguous crossing. If that fragment also uses
      the sj, it is a crossing in the SPLICED bank, which is a channel of its own.
    * ``t4`` ends at 650, which cuts a 50 bp region out of ``[500, 700)``. Nothing else here can be
      SPANNED: spanning needs one segment covering a region whole, so at 200 bp regions and 220 bp
      fragments it essentially never happens, and mature RNA can never span the region before a sj at
      all — it has no base past the exon end, it splices there. A 50 bp region is spanned by gDNA.
    * ``gdna_fraction`` puts genomic fragments in the intronic and intergenic regions, which is the
      unspliced bank and the two contained gDNA length pools.

    Every one of those three was added because the assertion below found the array empty. That is the
    fixture doing its job: each is a population the tally has, and a determinism gate that never sees
    one is not testing it.
    """
    scenario = Scenario(
        "worker_determinism",
        genome_length=6000,
        seed=TALLY_SEED,
        work_dir=tmp_path / "worker_determinism",
    )
    scenario.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(200, 700), (1200, 1700)], "abundance": 60},
            {"t_id": "t2", "exons": [(500, 700), (1200, 1700)], "abundance": 25},
            {"t_id": "t4", "exons": [(200, 650), (1200, 1700)], "abundance": 15},
        ],
    )
    scenario.add_gene(
        "g2",
        "-",
        [{"t_id": "t3", "exons": [(3000, 3400), (4000, 4400)], "abundance": 40}],
    )
    # 300 bp mean at 12,000 fragments, and BOTH numbers are load-bearing. The RESOLVED-gap arm needs a
    # fragment whose two 100 bp mates straddle the sj with the intron strictly inside the gap, so it
    # needs length ≥ ~300 on a 1,000 nt transcript — and the simulator models the effective length, so
    # a 600 bp fragment has 401 placements on that transcript against an 80 bp fragment's 921 and the
    # long tail is genuinely suppressed. A shorter mean or a smaller budget empties that arm entirely,
    # and the assertion below then holds over nothing.
    result = scenario.build_oracle(
        n_fragments=12000,
        gdna_fraction=0.3,
        sim_config=ReadSimConfig(
            frag_mean=300,
            frag_std=60,
            frag_min=80,
            frag_max=600,
            read_length=100,
            strand_specificity=1.0,
            seed=TALLY_SEED,
        ),
    )
    yield result
    scenario.cleanup()


def _tally(result, n_workers: int) -> AccumulatorPayload:
    """Scan at ``n_workers`` and return the accumulator payload."""
    _, _, _, payload = scan_and_buffer(
        str(result.bam_path),
        result.index,
        BamScanConfig(sj_strand_tag="auto", total_threads=n_workers),
    )
    assert payload is not None, "the scan produced no accumulator payload"
    return payload


def test_the_tally_is_bit_identical_at_1_2_4_and_8_workers(every_bank_oracle):
    baseline = _tally(every_bank_oracle, 1)

    # A tally of zeros would satisfy bit-identity trivially, so the baseline has to be shown to contain
    # something first. A bit-identity gate in this project has lied in exactly this way: an arm with zero
    # rows scored fully identical, because the comparison looped over the empty arm's rows.
    for key, why in [
        ("region_start_count", "nothing was deposited at all"),
        ("region_contained_count", "no fragment fitted inside a region"),
        (
            "boundary_unspliced_count",
            "no unspliced contiguous crossing — the mixture being deconvolved",
        ),
        ("boundary_spliced_count", "no SPLICED contiguous crossing — the certified-RNA channel"),
        ("sj_count", "no annotated sj was used"),
        ("pool_lengths", "no fragment entered a length pool"),
    ]:
        assert int(getattr(baseline, key).sum()) > 0, f"{key}: {why}"

    # Read off the payload's own fields, so a bank added later joins this gate automatically.
    array_keys = [
        f.name
        for f in dataclasses.fields(AccumulatorPayload)
        if isinstance(getattr(baseline, f.name), np.ndarray)
    ]
    assert len(array_keys) >= 16, (
        f"only {len(array_keys)} arrays on the payload; the gate is too narrow"
    )
    # The side buffer is the one bank this is not free for, and it is read the same way — off its own
    # fields, so an array added to it joins this gate too. Every other bank is a sum, compared below
    # exactly (counts) or within the derived tolerance (float64). The deferred queue is a LIST:
    # concatenating per-worker queues gives a different byte sequence at 1, 2, 4 and 8 workers with
    # identical CONTENTS, so the C++ export sorts on the record's own content — and this is what
    # says it does.
    deferred_keys = [f.name for f in dataclasses.fields(baseline.deferred)]
    assert baseline.deferred.n_fragments > 0, (
        "no fragment was deferred, so the canonical sort is compared only against an empty bank — which "
        "is bit-identical for free. The fixture must produce an undetermined gap."
    )

    n_deposits = int(np.asarray(baseline.region_start_count).sum())
    assert n_deposits > 0, "nothing was deposited, so this comparison could not have differed"

    for n_workers in (2, 4, 8):
        other = _tally(every_bank_oracle, n_workers)
        assert other.graph_hash == baseline.graph_hash, f"{n_workers} workers: different index"
        for key in array_keys:
            expected = np.asarray(getattr(baseline, key))
            actual = np.asarray(getattr(other, key))
            assert actual.dtype == expected.dtype, f"{n_workers} workers: {key} dtype"
            # One convention, two standards, and that is the point. A COUNT is an integer and integer
            # addition IS associative, so a count bank must be bit-identical at any worker count and is
            # asserted with no tolerance at all. A FRACTION is float64 and float addition is NOT
            # associative: the per-worker merge re-associates the sum, so those banks agree only to the
            # representation. The budget is derived from the deposit count, never fitted, and a count
            # bank that needed it would fail here.
            if expected.dtype == np.float64:
                # The budget is `n_deposits * EPS`: each cell accumulates at most one addition per
                # deposited fragment, and re-associating `n` round-to-nearest additions moves the sum
                # by at most `n` ulp. DERIVED from the fragment count — not the array size, which is
                # unrelated, and not a number anyone chose.
                assert np.allclose(actual, expected, rtol=n_deposits * EPS, atol=0.0), (
                    f"{n_workers} workers: {key} differs by MORE than the float64 representation — "
                    f"that is a merge defect, not re-association"
                )
            else:
                assert np.array_equal(actual, expected), (
                    f"{n_workers} workers: {key} is not bit-identical to the single-worker run — "
                    f"{int(np.count_nonzero(actual != expected))} of {expected.size} cells differ. "
                    f"This bank is an INTEGER count and integer addition is associative, so a "
                    f"difference here is a real merge defect."
                )
        for key in deferred_keys:
            expected = getattr(baseline.deferred, key)
            actual = getattr(other.deferred, key)
            assert actual.dtype == expected.dtype, f"{n_workers} workers: deferred.{key} dtype"
            assert np.array_equal(actual, expected), (
                f"{n_workers} workers: deferred.{key} is not bit-identical to the single-worker run. "
                f"The queue's ORDER is observable — this is the canonical sort failing, not a tally bug."
            )
        assert other.qc == baseline.qc, f"{n_workers} workers: qc"
        assert other.gap_resolution == baseline.gap_resolution, (
            f"{n_workers} workers: gap_resolution"
        )


def test_the_deferred_bank_holds_the_fragments_ITS_COUNTER_CLAIMS(every_bank_oracle):
    """The fragment census on a real scan. ``deposited + deferred + dropped_* == offered`` is worth
    nothing if the deferred term is a number with no fragments behind it.

    Two identities, and they are different statements: the bank holds as many records as the counter
    says, and the umbrella census's three ``gap_deferred_*`` partition exactly that population. The
    payload refuses both at its door, so this test's real content is that the fixture REACHES the
    population at all — a conservation identity over an empty term is satisfied by any bookkeeping.
    """
    payload = _tally(every_bank_oracle, 1)
    assert payload.deferred.n_fragments == payload.qc.deferred_undetermined_gap > 0
    assert payload.gap_resolution.deferred == payload.qc.deferred_undetermined_gap
    assert payload.gap_resolution.gap_resolved_spliced > 0, (
        "no fragment's gap was RESOLVED, so this fixture only exercises the deferral arm"
    )

    # Every record replays: two or more hypotheses, an extent inside its own reference, and a strand.
    runs = np.diff(payload.deferred.hypothesis_offsets)
    assert int(runs.min()) >= 2
    assert np.all(payload.deferred.end > payload.deferred.start)
    assert np.all(payload.deferred.ref < payload.n_refs)


def test_the_start_count_invariant_holds_on_a_real_scan(every_bank_oracle):
    """``sum(region_start_count) == deposited`` — the accumulator's one non-tautological invariant.

    An identity that can only be evaluated by re-running the deposit is satisfied by a deliberately
    broken replay, whatever the crossings contain. This one is checkable against a number the deposit
    counts independently, which is what makes it worth asserting.
    """
    tally = _tally(every_bank_oracle, 1)
    assert int(tally.region_start_count.sum()) == tally.qc.deposited
