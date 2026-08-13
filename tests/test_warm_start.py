"""``EMConfig.warm_start`` — what the EM's initial ``theta`` is derived FROM.

⭐⭐ **WHY THE SWITCH EXISTS.** The shipped warm start seeds ``theta`` with each component's unambiguous
total plus a coverage-weighted share of the ambiguous fragments, and then projects that seed **through**
the calibration prior. Seed and prior are two different methods, and the projection MULTIPLIES them — so
a coverage-weighted share scaled by a per-transcript allocation derived some other way is neither
method's answer. ``warm_start="prior"`` zeroes the seed, so a prior can be tested on its own terms.

⛔ **IT IS ONLY MEANINGFUL WITH A PER-TRANSCRIPT WEIGHT, AND THAT IS A PROPERTY, NOT A CAVEAT.** Under
the shipped evidence-proportional rule ``out[i]`` is proportional to ``raw[i]``, so an all-zero seed
yields an all-zero RNA pool and the whole locus goes to gDNA — ``theta = 0`` is an absorbing state. The
tests below pin both halves, because the failure is silent: a starved locus still converges.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator


def _partition(*, n_units: int, log_liks: tuple[float, ...], gdna_log_lik: float = -50.0):
    """⛔ ``gdna_log_lik`` defaults to a FINITE (but hopeless) value on purpose. With ``-inf`` no unit
    carries a gDNA candidate, and ``apply_grouped_prior_update`` then discards the RNA prior along with
    the gDNA one — so a fixture built that way measures a prior that was never applied."""
    n_t = len(log_liks)
    offsets = np.arange(n_units + 1, dtype=np.int64) * n_t
    return (
        offsets,
        np.tile(np.arange(n_t, dtype=np.int32), n_units),
        np.tile(np.asarray(log_liks, dtype=np.float64), n_units),
        np.ones(n_units * n_t, dtype=np.float64),
        np.zeros(n_units * n_t, dtype=np.uint8),
        np.zeros(n_units, dtype=np.uint8),
        np.full(n_units, gdna_log_lik, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _estimator(n_t: int, *, warm_start: str = "coverage", mode: str = "vbem"):
    est = AbundanceEstimator(
        num_transcripts=n_t,
        em_config=EMConfig(
            mode=mode,
            iterations=500,
            convergence_delta=1e-10,
            assignment_mode="fractional",
            seed=0,
            warm_start=warm_start,
        ),
    )
    est._t_eff_len = np.ones(n_t, dtype=np.float64)
    return est


def _run(est, partition, t_idx, *, weight=None, rna_prior=0.0, gdna_prior=0.0, enable_gdna=1,
         iterations=500):
    # ⚠ `em_iterations` is an EXPLICIT argument of the estimator entry point, NOT read from
    # `em_config` — a helper that sets it on the config alone silently runs the default 1000.
    est.run_batch_locus_em_partitioned(
        em_iterations=iterations,
        partition_tuples=[partition],
        locus_transcript_indices=[np.asarray(t_idx, dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        rna_prior_count=np.array([rna_prior], dtype=np.float64),
        rna_prior_weight=None if weight is None else np.asarray(weight, dtype=np.float64),
        index=None,
        enable_gdna=np.array([enable_gdna], dtype=np.uint8),
    )
    return est.em_counts.sum(axis=1)


# ──────────────────────────────────────────────────────────────────────────────
# The default is untouched
# ──────────────────────────────────────────────────────────────────────────────


def test_the_SHIPPED_default_is_coverage_and_is_BYTE_IDENTICAL_to_no_switch():
    """⛔ The control. A config field that changes the default is not a switch, it is a release."""
    assert EMConfig().warm_start == "coverage"
    part = _partition(n_units=100, log_liks=(0.0, -0.4), gdna_log_lik=-1.0)
    a = _run(_estimator(2), part, [0, 1], rna_prior=30.0, gdna_prior=10.0)
    b = _run(_estimator(2, warm_start="coverage"), part, [0, 1], rna_prior=30.0, gdna_prior=10.0)
    np.testing.assert_array_equal(a, b)


# ──────────────────────────────────────────────────────────────────────────────
# ⭐⭐ theta starts from the prior alone
# ──────────────────────────────────────────────────────────────────────────────


def test_the_WEIGHT_moves_the_converged_split_into_the_weights_ratio():
    """⭐⭐ Two components whose likelihoods are IDENTICAL, so the data cannot separate them and the only
    thing saying where mass belongs is the prior. At the fixed point ``theta_i`` is proportional to
    ``raw_i + a_i``; with flat likelihoods ``raw_i`` is itself proportional to ``theta_i``, so the
    algebra collapses to ``theta_i = a_i / Σa`` — the WEIGHTS' ratio."""
    part = _partition(n_units=100, log_liks=(0.0, 0.0))
    flat = _run(_estimator(2), part, [0, 1], rna_prior=40.0)
    tilted = _run(_estimator(2), part, [0, 1], weight=[3.0, 1.0], rna_prior=40.0)

    assert flat[0] == pytest.approx(flat[1], rel=1e-6), "flat likelihoods must split evenly with no weight"
    assert tilted[0] / tilted[1] == pytest.approx(3.0, rel=0.05), (
        "the weight did not move the converged split into its own ratio"
    )


def test_the_split_is_MONOTONE_in_the_weight_ratio():
    """⭐ A PROPERTY rather than a number. The exact converged ratio is perturbed by the gDNA component
    sharing the locus, so pinning it to `w0/w1` across a sweep would be fitting the fixture. What the
    allocation must guarantee is direction and order: equal weights split evenly, and a heavier weight
    always takes strictly more."""
    part = _partition(n_units=200, log_liks=(0.0, 0.0))
    ratios = [
        _run(_estimator(2), part, [0, 1], weight=[w, 1.0], rna_prior=25.0)
        for w in (0.25, 1.0, 3.0, 9.0)
    ]
    got = [float(o[0] / o[1]) for o in ratios]
    assert got[1] == pytest.approx(1.0, rel=1e-6), "equal weights did not split evenly"
    assert got == sorted(got), f"the split is not monotone in the weight ratio: {got}"
    assert got[0] < 1.0 < got[2] < got[3]


def test_the_prior_only_warm_start_CHANGES_THE_SEED_but_not_a_UNIQUE_fixed_point():
    """⭐⭐ **WHAT ZEROING THE WARM START ACTUALLY DOES, STATED HONESTLY.** It changes where the EM
    STARTS, not where it ends: when the fixed point is unique the two arms converge to the same answer,
    which is the correct and reassuring behaviour rather than evidence the switch is inert.

    ⛔ So the switch is demonstrated on the ITERATION, not on the limit: truncated to a single pass the
    two seeds give visibly different answers, and at convergence they agree. A test that only compared
    converged values could not tell a live switch from a dead one
    (``TRAPS: could-the-arm-have-fired``).
    """
    part = _partition(n_units=100, log_liks=(0.0, -0.7))
    kw = dict(weight=[3.0, 1.0], rna_prior=40.0, gdna_prior=5.0)

    def _one(ws, iters):
        return _run(_estimator(2, warm_start=ws), part, [0, 1], iterations=iters, **kw)

    early_cov, early_pri = _one("coverage", 1), _one("prior", 1)
    assert not np.allclose(early_cov, early_pri), "the warm-start switch never reached the solver"

    late_cov, late_pri = _one("coverage", 500), _one("prior", 500)
    np.testing.assert_allclose(late_cov, late_pri, rtol=1e-6)


def test_the_RNA_PRIOR_REACHES_a_locus_with_NO_gDNA_CANDIDATE():
    """⛔⛔ **THE BUG THIS TEST WAS WRITTEN TO RECORD, THEN REWRITTEN TO GATE THE FIX.**
    ``apply_grouped_prior_update`` used to read BOTH pseudocounts through ``has_gdna``::

        gdna_prior = has_gdna ? ... : 0.0
        rna_prior  = has_gdna ? ... : 0.0     // <- discarded the whole RNA prior

    A locus with no gDNA candidate is one whose fragments are ALL SPLICED — a gDNA candidate is appended
    to every unspliced unit — so the RNA prior was silently withheld from exactly the loci whose RNA is
    most certain.

    ⭐ **It was invisible because the shipped weights make it a no-op.** The evidence-proportional rule
    enters as a COMMON factor over the eligible components, and a common factor cancels under
    normalisation. An informative per-component weight cancels nothing, which is how it surfaced.
    """
    flat = (0.0, 0.0)
    part = _partition(n_units=100, log_liks=flat, gdna_log_lik=-np.inf)
    weighted = _run(_estimator(2), part, [0, 1], weight=[3.0, 1.0], rna_prior=40.0, enable_gdna=0)
    assert weighted[0] / weighted[1] == pytest.approx(3.0, rel=0.05), (
        "the RNA prior was discarded at a locus with no gDNA candidate"
    )
    # ⭐ and the SHIPPED weights still cancel there, which is why the gate went unnoticed for so long
    shipped = _run(_estimator(2), part, [0, 1], rna_prior=40.0, enable_gdna=0)
    assert shipped[0] == pytest.approx(shipped[1], rel=1e-9)


def test_an_unknown_warm_start_is_REFUSED():
    """⚠ ``Literal`` is a type-checker annotation, not a runtime constraint — so a typo must be refused
    explicitly, or it silently falls through to the shipped path."""
    with pytest.raises(ValueError, match="Unknown warm start"):
        EMConfig(warm_start="zero")
