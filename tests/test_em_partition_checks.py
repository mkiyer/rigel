"""The EM's input is checked once, at its bindings: a locus partition or a candidate CSR the EM could not read is
refused there (``batch_locus_em_partitioned``, ``connected_components``), so the EM's loops read every candidate
without a guard. Each refusal here was a silent skip or an out-of-range read before."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.native import connected_components


def _locus(candidates: list[list[int]]) -> list[np.ndarray]:
    """A partition whose units hold these candidates, each an even log-likelihood and the first count column."""
    offsets = np.cumsum([0] + [len(c) for c in candidates]).astype(np.int64)
    t_idx = np.concatenate([np.asarray(c, np.int32) for c in candidates] or [np.zeros(0, np.int32)])
    n, units = len(t_idx), len(candidates)
    return [
        offsets,
        t_idx,
        np.zeros(n, np.float32),
        np.ones(n, np.float32),
        np.zeros(n, np.uint8),
        np.zeros(units, np.uint8),
        np.full(units, -1.0, np.float32),
        np.full(units, -np.inf, np.float32),
    ]


def _run(part: list[np.ndarray], transcripts: list[int]):
    est = AbundanceEstimator(4, em_config=EMConfig(n_threads=1, seed=7))
    return est.run_batch_locus_em_partitioned(
        [tuple(part)],
        [np.asarray(transcripts, np.int32)],
        np.array([5.0]),
        rna_prior_count=np.array([3.0]),
        intergenic_rate=20.0,
    )


def _with(index: int, value: np.ndarray) -> list[np.ndarray]:
    part = _locus([[0, 1]] * 3)
    part[index] = value
    return part


def test_a_partition_the_em_can_read_runs():
    _run(_locus([[0, 1]] * 3), [0, 1])


@pytest.mark.parametrize(
    "part, transcripts, said",
    [
        (
            _with(1, np.array([0, 1, 0, 3, 0, 1], np.int32)),
            [0, 1],
            "transcript 3, which is not one of its",
        ),
        (_with(4, np.full(6, 200, np.uint8)), [0, 1], "count column 200"),
        (_locus([]), [0, 1], "holds no unit"),
        (_locus([[0, 1]] * 3), [0, 9], "transcript 9, outside the 4"),
        (_with(2, np.zeros(5, np.float32)), [0, 1], "candidate array whose length"),
        (_with(6, np.zeros(2, np.float32)), [0, 1], "unit array whose length"),
        (_with(0, np.array([0, 2, 1, 6], np.int64)), [0, 1], "offsets that fall"),
    ],
)
def test_a_partition_the_em_cannot_read_is_refused(part, transcripts, said):
    with pytest.raises(ValueError, match=said):
        _run(part, transcripts)


def test_a_candidate_csr_naming_a_transcript_past_the_index_is_refused():
    with pytest.raises(ValueError, match="candidate transcript 7 is outside the 4"):
        connected_components(np.array([0, 2], np.int64), np.array([0, 7], np.int32), np.int32(4))
