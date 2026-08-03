"""⭐ THE SAME BAM MUST GIVE THE SAME ANSWER, whatever the scan's thread count.

    Cause and measurement:, the "scan-order independence" entry
    Related: `tests/native/test_accumulator_worker_determinism.py` — the TALLY's version of this

The accumulator's tally is order-independent for free: every bank is a sum of integers and integer
addition is associative, so a per-worker merge is exact whatever order the chunks arrived in. ⛔ **The
FRAGMENT BUFFER is not.** The scanner streams finalized chunks from worker threads in *completion* order,
so the row a fragment occupies depends on which worker got there first.

That was invisible for as long as nothing downstream cared about row order. The EM's
``assignment_mode="sample"`` cares: it draws one categorical sample per unit from a per-locus RNG, so a
permutation of the units permutes which fragment gets which draw — and units in different equivalence
classes have different posteriors, so the assignment genuinely changes.

**Measured before the fix**, one 3,000-fragment BAM, seed fixed, three runs each:

===================  ==================  ===========================================
scan ``total_threads``  chunks emitted   seeded ``sample`` reproduces?
===================  ==================  ===========================================
1                    1                   ✅
4                    1                   ✅ — one chunk, so nothing to reorder
8                    ~6                  ⛔ counts differ by up to **22**
16 (the default)     ~6                  ⛔ counts differ by up to **3**
===================  ==================  ===========================================

⚠ **The 1- and 4-thread passes were an artifact of scale, not a property.** They emitted a single chunk
because the BAM is small. A real library emits many chunks at any thread count above one, so the shipped
default was not reproducible on real data at all.

⭐ **THE FIX IS AN ORDERING, NOT A SEED.** ``frag_id`` is assigned by the single reader thread in BAM
order and is a stable identity — measured: the ``frag_id → fragment`` mapping is byte-identical across
scans while the buffer's row order is not. So a locus's units are ordered by ``frag_id``, once, where
``MultiLocus.unit_indices`` is built; every per-locus array is scattered by that index list, so all of
them inherit the canonical order and no consumer has to know.

⛔ **WHY NOT SEED EACH DRAW FROM THE FRAGMENT INSTEAD.** Deriving the draw from a hash of the fragment's
own content would also be order-independent — and would be **wrong**. Identical fragments would hash
identically and so all draw identically, turning a 60/40 posterior into 100/0 for every group of
duplicates. Ordering by identity and keeping one stream per locus preserves the multinomial spread, which
is the whole point of the sampling mode. (Owner ruling, D-D.)
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario


SEED = 7

#: Thread counts spanning single-chunk and many-chunk scans. ⚠ 1 and 4 alone would pass on a small BAM
#: before the fix — the gate needs a count high enough that the scanner really does emit several chunks
#: concurrently, which is what makes the row order vary.
THREAD_COUNTS = (1, 4, 16)


@pytest.fixture(scope="module")
def oracle(tmp_path_factory):
    """A locus with two isoforms sharing an exon, so units land in DIFFERENT equivalence classes.

    ⚠ That is what gives the gate teeth. Fragments that are all interchangeable would be permuted with no
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


#: ⚠ Two repeats per thread count, not one. The permutation is *random*, so a single pair can agree by
#: luck — measured: two 16-thread runs disagreed while 1 vs 4 vs 16 happened to match on the same tree.
#: Six runs compared against one baseline is what makes the gate reliable rather than occasionally lucky.
REPEATS = 2


@pytest.mark.parametrize("mode", ["sample", "map", "fractional"])
def test_the_answer_IS_THE_SAME_AT_EVERY_SCAN_THREAD_COUNT(oracle, mode):
    """⭐ **THE CONTRACT, for all three assignment modes:** one BAM and one seed give one answer.

    Byte-identical, across thread counts *and* across repeats — which is a strictly stronger statement
    than "reproducible at the default", and the only version that survives someone changing the default.

    ⛔ **All three modes, and each fails a different way.** ``sample`` moves by whole counts, because a
    permutation changes which fragment draws which sample. ``map`` and ``fractional`` are immune to that —
    they read each unit's posterior alone — but ``fractional`` scatters float posteriors into shared
    accumulators, so a permutation reorders the summation and the answer drifts by ULPs. One ordering
    fixes both; asserting only the first would leave the second to be rediscovered later as "flaky".
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
    """⚠ Non-vacuity, and it is not optional here.

    The contract above passes trivially if the scan hands the buffer back in the same row order every
    time — there is then no permutation to be immune to, and the gate reads as coverage it does not have.

    ⚠ Retried, because the reordering is a RACE. Measured at 6/6 on this fixture, but a gate that asserts
    a race fired on the first try is a gate that goes red for the wrong reason on a loaded machine.
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
    """⛔ The assumption the whole fix rests on, pinned where the fix is gated.

    Ordering by ``frag_id`` canonicalises the units **only if no two units share one**. If ids could
    repeat, ``argsort(kind="stable")`` would break the tie by the unit's buffer row — silently
    reintroducing the exact dependence the ordering removes, on exactly the duplicates it matters for.

    ⚠ Uniqueness, not contiguity. The EM's units are a SUBSET of the scanned fragments (measured on this
    fixture: 2,118 units from 3,000 fragments), so the ids have gaps — which is fine, because the sort
    only needs a total order, not a dense one.
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
