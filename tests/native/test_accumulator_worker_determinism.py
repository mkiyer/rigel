"""The same BAM must give the same tally at any worker count — bit for bit.

     step 8

⭐ THIS IS NEWLY ACHIEVABLE, and it is the reason every channel is an integer. The float accumulator this
replaced differed by only ~3.7e-7 per cell between worker counts — and that propagated to a **~2.6 %
difference in the calibration output**, so the same BAM gave different answers on different machines.
Integer addition is associative, so the per-worker merge is now exact.

⚠ WHY IT IS NOT ENOUGH TO TEST ``Accumulator.merge_from`` DIRECTLY (the parity module does that). This
exercises the *scanner's* worker path: each worker builds its own ``AccumulatorSet`` from the scanner's
members, and the chunk→worker split is data-dependent. A worker whose set was constructed differently —
junctions installed on the template but not on the copies, say — is invisible to a merge test on one
accumulator and shows up here.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from rigel.config import BamScanConfig
from rigel.pipeline import scan_and_buffer
from rigel.scan_payload import AccumulatorPayload
from rigel.sim import ReadSimConfig, Scenario


SEED = 20260730


@pytest.fixture
def oracle(tmp_path):
    """A scenario built so that **every** bank of the tally receives something.

    ⚠ The obvious scenario does not. Two single-isoform genes deposit into the nodes and the junction
    edges and leave **both contiguous-edge banks identically zero**: every cut is an exon boundary, so a
    mature fragment either fits inside an exon (contained) or splices across the gap (junction edge), and
    it never has bases on both sides of a line. A bit-identity gate over an all-zero array passes for the
    wrong reason — this project has already had one report "32/32 IDENTICAL" on an arm with zero rows. So:

    * ``t2`` starts at 500, **inside** ``t1``'s first exon. That makes 500 a cut, and a ``t1`` fragment
      spanning it has bases on both sides — a contiguous crossing. If that fragment also uses the
      junction, it is a crossing in the SPLICED bank, which is the channel the old design merged away.
    * ``t4`` ends at 650, which cuts a **50 bp** node out of ``[500, 700)``. Nothing else here can be
      SPANNED: spanning needs one segment covering a node whole, so at 200 bp nodes and 220 bp fragments
      it essentially never happens, and mature RNA can never span the node before a junction at all —
      it has no base past the exon end, it splices there. A 50 bp node is spanned by ordinary gDNA.
    * ``gdna_fraction`` puts genomic fragments in the intronic and intergenic nodes, which is the
      unspliced bank and the two pure gDNA length pools.

    ⚠ Every one of those three was added because the assertion below found the array empty. That is the
    fixture doing its job: each is a population the tally has, and a determinism gate that never sees one
    is not testing it.
    """
    scenario = Scenario(
        "worker_determinism",
        genome_length=6000,
        seed=SEED,
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
    # ⚠ 300 bp mean at 12,000 fragments, and BOTH numbers are load-bearing. The RESOLVED-gap arm needs
    # a fragment whose two 100 bp mates straddle the junction with the intron strictly inside the gap,
    # so it needs length ≥ ~300 on a 1,000 nt transcript — and the simulator now models the effective
    # length, so a 600 bp fragment has 401 placements on that transcript against an 80 bp fragment's
    # 921 and the long tail is genuinely suppressed (a 2.8x tilt). At the old 220/4,000 the arm read
    # exactly 0 and the assertion below said so. At 300/12,000 it reads 51, which is a margin rather
    # than a coincidence.
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
            seed=SEED,
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


def test_the_tally_is_bit_identical_at_1_2_4_and_8_workers(oracle):
    baseline = _tally(oracle, 1)

    # A tally of zeros would satisfy bit-identity trivially, so the baseline has to be shown to contain
    # something first. ⚠ A bit-identity gate in this project has already lied in exactly this way: an arm
    # with ZERO rows scored "32/32 IDENTICAL" because the comparison looped over the empty arm's rows.
    for key, why in [
        ("node_start_count", "nothing was deposited at all"),
        ("node_contained_count", "no fragment fitted inside a node"),
        (
            "edge_unspliced_count",
            "no unspliced contiguous crossing — the mixture being deconvolved",
        ),
        ("edge_spliced_count", "no SPLICED contiguous crossing — the certified-RNA channel"),
        ("sj_count", "no annotated junction was used"),
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
    # ⭐ THE SIDE BUFFER IS THE ONE BANK THIS IS NOT FREE FOR, and it is read the same way — off its own
    # fields, so an array added to it joins this gate too. Every other bank is a sum of integers and
    # integer addition is associative, so a per-worker merge is exact whatever order the chunks arrived in.
    # ⛔ The deferred queue is a LIST. Concatenating per-worker queues gives a different byte sequence at 1,
    # 2, 4 and 8 workers with identical CONTENTS, so the C++ export sorts on the record's own content — and
    # this is what says it does.
    deferred_keys = [f.name for f in dataclasses.fields(baseline.deferred)]
    assert baseline.deferred.n_fragments > 0, (
        "no fragment was deferred, so the canonical sort is compared only against an empty bank — which "
        "is bit-identical for free. The fixture must produce an undetermined gap."
    )

    for n_workers in (2, 4, 8):
        other = _tally(oracle, n_workers)
        assert other.graph_hash == baseline.graph_hash, f"{n_workers} workers: different index"
        for key in array_keys:
            expected = np.asarray(getattr(baseline, key))
            actual = np.asarray(getattr(other, key))
            assert actual.dtype == expected.dtype, f"{n_workers} workers: {key} dtype"
            assert np.array_equal(actual, expected), (
                f"{n_workers} workers: {key} is not bit-identical to the single-worker run — "
                f"{int(np.count_nonzero(actual != expected))} of {expected.size} cells differ"
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


def test_the_deferred_bank_holds_the_fragments_ITS_COUNTER_CLAIMS(oracle):
    """⭐ **S1's gate, on a real scan.** ``deposited + deferred + dropped_* == offered`` is worth nothing if
    the deferred term is a number with no fragments behind it.

    Two identities, and they are different statements: the bank holds as many records as the counter says,
    and the umbrella census's three ``gap_deferred_*`` partition exactly that population. ⚠ The payload
    refuses both at its door, so this test's real content is that the fixture REACHES the population at all
    — a conservation identity over an empty term is satisfied by any bookkeeping.
    """
    payload = _tally(oracle, 1)
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


def test_the_start_count_invariant_holds_on_a_real_scan(oracle):
    """``sum(node_start_count) == deposited`` — the accumulator's one non-tautological invariant.

    The three "conservation identities" it replaced could only be evaluated by re-running the deposit, so a
    deliberately broken replay satisfied all three while 91 % of the crossings were junk. This one is
    checkable against a number the deposit counts independently.
    """
    tally = _tally(oracle, 1)
    assert int(tally.node_start_count.sum()) == tally.qc.deposited
