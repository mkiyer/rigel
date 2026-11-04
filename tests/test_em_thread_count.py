"""The EM gives one answer whatever its thread count.

A locus whose E-step work reaches a thread's fair share runs its E-step on every thread. Each equivalence class's
column sums go to a slot of their own and are added in class order on every path, so the thread count changes only
which thread runs a class. Sums reduced in thread order would group differently at each thread count, and SQUAREM
turns those ulps into whole isoforms.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario

SEED = 19

#: Thread counts at which a split of each class into row ranges per thread lands at a different row.
THREAD_COUNTS = (2, 3, 4, 16)


@pytest.fixture(scope="module")
def one_locus(tmp_path_factory):
    """One gene of three isoforms that share exons, under gDNA: one locus, so at any thread count above one it is
    the whole of the work and its E-step runs threaded, with classes deep enough that a split into row ranges per
    thread would land at a different row for every thread count."""
    work = tmp_path_factory.mktemp("em_threads")
    scenario = Scenario("em_threads", genome_length=9000, seed=SEED, work_dir=work / "sim")
    scenario.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(1000, 1600), (2400, 3000), (4000, 4800)], "abundance": 100},
            {"t_id": "t2", "exons": [(1000, 1600), (4000, 4800)], "abundance": 70},
            {"t_id": "t3", "exons": [(1300, 1600), (2400, 3000), (4000, 4500)], "abundance": 40},
        ],
    )
    result = scenario.build_oracle(
        n_fragments=40000,
        gdna_fraction=0.3,
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


def _solve(one_locus, threads: int):
    config = PipelineConfig(
        em=EMConfig(n_threads=threads, assignment_mode="fractional"),
        scan=BamScanConfig(sj_strand_tag="auto", total_threads=1),
        emit_locus_stats=True,
    )
    return run_pipeline(one_locus.bam_path, one_locus.index, config=config).estimator


@pytest.fixture(scope="module")
def single_thread(one_locus):
    return _solve(one_locus, 1)


@pytest.mark.parametrize("threads", THREAD_COUNTS)
def test_the_em_answer_is_the_same_at_every_thread_count(one_locus, single_thread, threads):
    threaded = _solve(one_locus, threads)
    assert [s["estep_threads_used"] for s in threaded.locus_stats] == [threads], (
        "the locus must run its E-step on every thread, or this gate holds nothing"
    )
    reference = np.asarray(single_thread.t_counts)
    answer = np.asarray(threaded.t_counts)
    moved = np.abs(answer - reference).max()
    assert np.array_equal(answer, reference), (
        f"{threads} EM threads moved a transcript count by {moved!r} against one thread"
    )
    assert threaded.gdna_em_count == single_thread.gdna_em_count
