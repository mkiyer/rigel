"""The EM's per-unit class masses: each unit's posterior on the annotated transcripts, the synthetic spans, the
locus's gDNA and the intergenic term (`run_batch_locus_em_partitioned(emit_class_mass=True)`).

A paired proper score reads a source forecast per fragment, and the annotated BAM holds only the winner's class and
its float32 posterior. The masses are the posterior the assignment itself keeps or draws from, so they must be a
distribution per unit, add up to the pools the EM books, sit in the per-unit annotation's row order, and change no
count. The fixture holds every class: two loci, one with an annotated transcript, a synthetic span, gDNA candidates
and units with an intergenic place, and one of two annotated transcripts.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator
from rigel.native import EM_CLASSES, EM_WINNER_INTERGENIC

ANN, NAS, GDNA, IG = (
    EM_CLASSES.index(name) for name in ("annotated", "nascent", "gdna", "intergenic")
)
IS_SYNTHETIC = np.array([False, True, False, False])
#: The bound a sum of n posteriors in [0, 1] keeps from its exact value: gamma(n) = n·u / (1 − n·u).
U = np.finfo(np.float64).eps / 2


def _gamma(n: int) -> float:
    return n * U / (1 - n * U)


def _locus(rows):
    """A partition from ``(candidates, log_liks, gdna_log_lik, intergenic_log_lik)`` per unit."""
    offsets = np.cumsum([0] + [len(c) for c, *_ in rows]).astype(np.int64)
    t_idx = np.concatenate([np.asarray(c, np.int32) for c, *_ in rows])
    ll = np.concatenate([np.asarray(lk, np.float32) for _, lk, *_ in rows])
    n_cand = len(t_idx)
    return (
        offsets,
        t_idx,
        ll,
        np.ones(n_cand, np.float32),
        np.zeros(n_cand, np.uint8),
        np.zeros(len(rows), np.uint8),
        np.array([g for *_, g, _ in rows], np.float32),
        np.array([i for *_, i in rows], np.float32),
    )


def _fixture():
    a = (
        [([0, 1], (0.0, -0.4), -1.0, -np.inf)] * 30
        + [([0, 1], (-1.5, 0.0), -1.0, -np.inf)] * 12
        + [([0], (0.0,), -0.5, -np.inf)] * 10
        + [([0, 1], (0.0, 0.0), -0.2, 0.0)] * 6
    )
    b = [([2, 3], (0.0, -1.2), -np.inf, -np.inf)] * 20 + [
        ([2, 3], (-0.6, 0.0), -np.inf, -np.inf)
    ] * 9
    return (
        [_locus(a), _locus(b)],
        [np.array([0, 1], np.int32), np.array([2, 3], np.int32)],
        [58, 29],
    )


def _run(mode: str, *, emit_class_mass: bool):
    est = AbundanceEstimator(
        4,
        em_config=EMConfig(assignment_mode=mode, n_threads=1, seed=7, iterations=300),
        is_synthetic=IS_SYNTHETIC,
    )
    parts, loci, _ = _fixture()
    out = est.run_batch_locus_em_partitioned(
        parts,
        loci,
        np.array([5.0, 0.0]),
        rna_prior_count=np.array([3.0, 2.0]),
        intergenic_rate=20.0,
        emit_assignments=True,
        emit_class_mass=emit_class_mass,
    )
    return est, out


@pytest.mark.parametrize("mode", ["fractional", "sample"])
def test_each_units_masses_are_a_distribution_in_the_annotations_row_order(mode):
    _, out = _run(mode, emit_class_mass=True)
    _, _, _, winner_tid, winner_post, n_candidates, mass = out
    n_units = sum(_fixture()[2])
    assert mass.shape == (n_units, len(EM_CLASSES)) and mass.dtype == np.float64
    assert (mass >= 0.0).all()
    np.testing.assert_allclose(mass.sum(axis=1), 1.0, rtol=0, atol=_gamma(3 * len(EM_CLASSES) + 2))
    # every class is reached, so a class booked into another's column would show
    assert (mass.sum(axis=0) > 1.0).all(), mass.sum(axis=0)
    # the winner's class holds at least the winner's posterior (float32): the rows are the annotation's rows
    cls = np.where(
        winner_tid >= 0,
        np.where(IS_SYNTHETIC[np.maximum(winner_tid, 0)], NAS, ANN),
        np.where(winner_tid == EM_WINNER_INTERGENIC, IG, GDNA),
    )
    held = mass[np.arange(n_units), cls]
    assert (held >= winner_post.astype(np.float64) * (1 - np.finfo(np.float32).eps)).all()
    # and a unit's candidate count is the annotation's: the intergenic units have a third candidate
    assert (n_candidates[52:58] == 4).all() and (n_candidates[:42] == 3).all()


def test_under_fractional_assignment_the_columns_are_the_pools():
    """The class columns add up to what the EM booked: the annotated transcripts' and the spans' counts, each locus's
    gDNA exactly (one gDNA candidate per unit, added in the same order), and the intergenic term's mass."""
    est, out = _run("fractional", emit_class_mass=True)
    total_gdna, _, locus_gdna, *_, mass = out
    counts = est.em_counts.sum(axis=1)
    n = mass.shape[0]
    np.testing.assert_allclose(
        mass[:, ANN].sum(), counts[~IS_SYNTHETIC].sum(), rtol=_gamma(n), atol=0
    )
    np.testing.assert_allclose(
        mass[:, NAS].sum(), counts[IS_SYNTHETIC].sum(), rtol=_gamma(n), atol=0
    )
    first = np.cumsum([0] + _fixture()[2])
    for li in range(len(locus_gdna)):
        booked = 0.0
        for p in mass[first[li] : first[li + 1], GDNA]:  # in unit order, as the solver adds them
            booked += float(p)
        assert booked == locus_gdna[li]
    assert mass[:, IG].sum() == pytest.approx(est.gdna_intergenic_em_count, rel=_gamma(n))
    assert mass[:, GDNA].sum() == pytest.approx(total_gdna, rel=_gamma(n))


@pytest.mark.parametrize("mode", ["fractional", "sample"])
def test_asking_for_the_masses_changes_no_count(mode):
    off_est, off = _run(mode, emit_class_mass=False)
    on_est, on = _run(mode, emit_class_mass=True)
    assert len(on) == len(off) + 1
    for a, b in zip(off, on[: len(off)]):
        np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(off_est.em_counts, on_est.em_counts)
    np.testing.assert_array_equal(off_est.posterior_mean(), on_est.posterior_mean())
    assert off_est.gdna_intergenic_em_count == on_est.gdna_intergenic_em_count


def test_no_masses_are_made_unless_asked():
    _, out = _run("fractional", emit_class_mass=False)
    assert len(out) == 6


def _stats_and_counts(parts):
    est = AbundanceEstimator(
        4,
        em_config=EMConfig(assignment_mode="fractional", n_threads=1, seed=7, iterations=300),
        is_synthetic=IS_SYNTHETIC,
    )
    _, loci, _ = _fixture()
    est.run_batch_locus_em_partitioned(
        parts,
        loci,
        np.array([5.0, 0.0]),
        rna_prior_count=np.array([3.0, 2.0]),
        intergenic_rate=20.0,
        emit_locus_stats=True,
    )
    return sum(s["ec_total_elements"] for s in est.locus_stats), est.em_counts.sum(axis=1)


def test_units_with_identical_rows_share_one_row_of_the_e_step():
    """The fixture's 87 units hold 6 distinct log-likelihood rows, gDNA and intergenic terms included (4 in the first
    locus's 3 candidate sets, 2 in the second's one), so the E-step reads 16 entries, not 228. Its answer is the one
    the units give as distinct rows: nudging each unit's first log-likelihood by its own few parts in 1e7 merges
    nothing and moves no count beyond the nudge, while a merged row that lost its units' weight would move whole
    fragments."""
    parts, _, _ = _fixture()
    entries, merged = _stats_and_counts(parts)
    assert entries == 16
    nudged = []
    for offsets, t_idx, ll, *rest in parts:
        ll = ll.copy()
        ll[offsets[:-1]] += np.arange(len(offsets) - 1, dtype=np.float32) * np.float32(2e-7)
        nudged.append((offsets, t_idx, ll, *rest))
    entries_distinct, distinct = _stats_and_counts(nudged)
    assert entries_distinct == 228
    np.testing.assert_allclose(merged, distinct, rtol=1e-3, atol=1e-3)
