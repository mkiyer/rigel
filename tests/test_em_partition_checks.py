"""The EM's input is checked once, at its bindings: a locus partition or a candidate CSR the EM could not read is
refused there (``batch_locus_em_partitioned``, ``connected_components``), so the EM's loops read every candidate
without a guard. Each refusal here was a silent skip or an out-of-range read before."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from _em_harness import em_partition

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
        [em_partition(*part)],
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


def _run_loci(parts: list, transcripts: list, gdna_prior: list[float], rna_prior: list[float]):
    est = AbundanceEstimator(4, em_config=EMConfig(n_threads=1, seed=7))
    return est.run_batch_locus_em_partitioned(
        [em_partition(*p) for p in parts],
        [np.asarray(t, np.int32) for t in transcripts],
        np.asarray(gdna_prior),
        rna_prior_count=np.asarray(rna_prior),
        intergenic_rate=20.0,
    )


def test_a_transcript_two_loci_hold_is_refused():
    """Loci are disjoint: two naming one transcript would add into its output row from two workers at once."""
    part = _locus([[0, 1]] * 3)
    with pytest.raises(ValueError, match="names transcript 0, which locus 0 holds"):
        _run_loci([part, part], [[0, 1], [0, 1]], [5.0, 5.0], [3.0, 3.0])


def test_a_locus_naming_a_transcript_twice_is_refused():
    with pytest.raises(ValueError, match="names transcript 0 twice"):
        _run(_locus([[0, 1]] * 3), [0, 0, 1])


@pytest.mark.parametrize("gdna, rna", [(-1.0, 3.0), (5.0, float("nan")), (float("inf"), 3.0)])
def test_a_locus_prior_count_that_is_not_a_count_is_refused(gdna, rna):
    """A prior count is a count of pseudo-fragments; one that is negative or not finite was read as 0 before."""
    with pytest.raises(ValueError, match="prior count is negative or not finite"):
        _run_loci([_locus([[0, 1]] * 3)], [[0, 1]], [gdna], [rna])


def _named(part: list[np.ndarray], **drop_or_rename) -> SimpleNamespace:
    """``part`` as a named partition, with each ``old=new`` renaming a field (``new=None`` drops it)."""
    fields = dict(vars(em_partition(*part)))
    for old, new in drop_or_rename.items():
        value = fields.pop(old)
        if new is not None:
            fields[new] = value
    return SimpleNamespace(**fields)


def test_the_em_reads_each_array_by_its_name_not_its_position():
    """Four float arrays cross the boundary, two per candidate and two per unit, and the same-typed ones are
    interchangeable by shape and dtype: read by position, a reordering passes every check and swaps their meanings.
    Read by name, the order a partition lists its fields in is not read at all, and the field names are the contract."""
    part = _locus([[0, 1]] * 3)
    part[3] = np.array(
        [1.0, 30.0, 1.0, 30.0, 1.0, 30.0], np.float32
    )  # coverage weights unlike the log-likelihoods
    part[7] = np.full(3, -2.0, np.float32)  # an intergenic reading unlike the gDNA one
    reference = _run_named(_named(part))
    shuffled = SimpleNamespace(**dict(reversed(list(vars(_named(part)).items()))))
    for got, want in zip(_run_named(shuffled), reference):
        np.testing.assert_array_equal(got, want)


@pytest.mark.parametrize(
    "rename, said",
    [
        ({"gdna_excluded": None}, "has no array gdna_excluded"),
        ({"gdna_excluded": "is_spliced"}, "has no array gdna_excluded"),
        ({"intergenic_log_liks": None}, "has no array intergenic_log_liks"),
    ],
)
def test_a_partition_missing_a_named_array_is_refused(rename, said):
    with pytest.raises(ValueError, match=said):
        _run_named(_named(_locus([[0, 1]] * 3), **rename))


def test_a_positional_partition_is_refused():
    """A tuple carries no names, so nothing in it says which array is which."""
    with pytest.raises(ValueError, match="has no array offsets"):
        _run_named(tuple(_locus([[0, 1]] * 3)))


def _run_named(partition):
    est = AbundanceEstimator(4, em_config=EMConfig(n_threads=1, seed=7))
    total, rna, gdna = est.run_batch_locus_em_partitioned(
        [partition],
        [np.asarray([0, 1], np.int32)],
        np.array([5.0]),
        rna_prior_count=np.array([3.0]),
        intergenic_rate=20.0,
    )
    return (
        np.asarray(total),
        rna,
        gdna,
        est.em_counts.copy(),
        np.asarray(est.gdna_intergenic_em_count),
    )
