"""The ZN tag holds a fragment's EM candidate count whole, from the EM through the annotation table to the BAM.

A BAM integer tag holds an int32, and a unit's candidates are counted in an int, so a count above int16's 32,767 is
written as it is: never clipped at 32,767 by the EM, never wrapped negative by the table.
"""

import numpy as np
import pysam
import pytest

from rigel.annotate import AnnotationTable
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.estimator import AbundanceEstimator
from rigel.pipeline import run_pipeline
from rigel.sim import ReadSimConfig, Scenario

#: A candidate count int16 cannot hold: above its largest value, 32,767.
N_WIDE = 40_000
SEED = 7


def test_the_em_counts_every_candidate_of_a_unit():
    """One unit over N_WIDE transcripts, no gDNA and no intergenic place: its count is N_WIDE."""
    est = AbundanceEstimator(
        num_transcripts=N_WIDE,
        em_config=EMConfig(assignment_mode="fractional", seed=0, n_threads=1),
    )
    partition = (
        np.array([0, N_WIDE], dtype=np.int64),
        np.arange(N_WIDE, dtype=np.int32),
        np.zeros(N_WIDE, dtype=np.float32),
        np.ones(N_WIDE, dtype=np.float32),
        np.zeros(N_WIDE, dtype=np.uint8),
        np.ones(1, dtype=np.uint8),
        np.full(1, -np.inf, dtype=np.float32),
        np.full(1, -np.inf, dtype=np.float32),
    )
    *_, n_candidates = est.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.arange(N_WIDE, dtype=np.int32)],
        gdna_prior_count=np.zeros(1, dtype=np.float64),
        emit_assignments=True,
    )
    assert n_candidates.tolist() == [N_WIDE]


def test_the_annotation_table_holds_the_count_whole():
    table = AnnotationTable.create(1)
    for _ in range(2):  # the second row grows the table
        table.add_batch(
            frag_ids=np.zeros(1, dtype=np.int64),
            best_tids=np.zeros(1, dtype=np.int32),
            best_gids=np.zeros(1, dtype=np.int32),
            tx_flags=np.zeros(1, dtype=np.uint8),
            posteriors=np.ones(1, dtype=np.float32),
            frag_classes=np.zeros(1, dtype=np.int8),
            n_candidates=np.full(1, N_WIDE, dtype=np.int32),
        )
    assert table.n_candidates[: table.size].tolist() == [N_WIDE, N_WIDE]


@pytest.fixture
def scenario(tmp_path):
    sc = Scenario("zn_tag", genome_length=4000, seed=SEED, work_dir=tmp_path / "sc")
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(300, 700), (1100, 1500)], "abundance": 60},
            {"t_id": "t2", "exons": [(300, 700), (1200, 1500)], "abundance": 40},
        ],
    )
    result = sc.build_oracle(
        n_fragments=200,
        sim_config=ReadSimConfig(
            frag_mean=200,
            frag_std=30,
            frag_min=80,
            frag_max=450,
            read_length=100,
            strand_specificity=1.0,
            seed=SEED,
        ),
    )
    yield result
    sc.cleanup()


def test_the_bam_tag_holds_the_count_whole(scenario, tmp_path, monkeypatch):
    """The EM reports N_WIDE candidates for every unit it solves; each such fragment's ZN tag must read N_WIDE back,
    through the run's own path from the EM's output to the table and the writer."""
    solve = AbundanceEstimator.run_batch_locus_em_partitioned

    def wide(self, *args, **kwargs):
        *out, n_candidates = solve(self, *args, **kwargs)
        return (*out, np.full(len(n_candidates), N_WIDE, dtype=np.int32))

    monkeypatch.setattr(AbundanceEstimator, "run_batch_locus_em_partitioned", wide)
    annotated = tmp_path / "annotated.bam"
    run_pipeline(
        scenario.bam_path,
        scenario.index,
        config=PipelineConfig(
            em=EMConfig(seed=SEED, n_threads=1),
            scan=BamScanConfig(sj_strand_tag="auto"),
            annotated_bam_path=annotated,
        ),
    )
    with pysam.AlignmentFile(str(annotated), "rb", check_sq=False) as bam:
        zn = {r.get_tag("ZN") for r in bam.fetch(until_eof=True)}
    assert N_WIDE in zn
    # beside the EM's units: a deterministic fragment's one candidate, and none for a fragment no pool holds
    assert zn <= {0, 1, N_WIDE}
