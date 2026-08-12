"""Production-path tests for the native batch locus EM wrapper."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.config import EMConfig
from rigel.estimator import AbundanceEstimator


def _partition(
    *,
    n_units: int,
    log_liks: tuple[float, ...],
    gdna_log_lik: float = -np.inf,
    is_spliced: bool = False,
):
    n_t = len(log_liks)
    offsets = np.arange(n_units + 1, dtype=np.int64) * n_t
    return (
        offsets,
        np.tile(np.arange(n_t, dtype=np.int32), n_units),
        np.tile(np.asarray(log_liks, dtype=np.float64), n_units),
        np.ones(n_units * n_t, dtype=np.float64),
        np.zeros(n_units * n_t, dtype=np.uint8),
        np.full(n_units, 1 if is_spliced else 0, dtype=np.uint8),
        np.full(n_units, gdna_log_lik, dtype=np.float64),
        np.zeros(n_units, dtype=np.int32),
        np.zeros(n_units, dtype=np.uint8),
    )


def _estimator(
    n_t: int, *, mode: str = "vbem", gdna_em_llr_bias: float = 0.0
) -> AbundanceEstimator:
    est = AbundanceEstimator(
        num_transcripts=n_t,
        em_config=EMConfig(
            mode=mode,
            iterations=200,
            convergence_delta=1e-9,
            assignment_mode="fractional",
            seed=0,
            gdna_em_llr_bias=gdna_em_llr_bias,
        ),
    )
    est._t_eff_len = np.ones(n_t, dtype=np.float64)
    return est


def test_equal_rna_likelihoods_split_evenly():
    n_units = 100
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(0.0, 0.0))],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
    )

    assert total_gdna == pytest.approx(0.0)
    assert locus_gdna[0] == pytest.approx(0.0)
    assert locus_rna[0] == pytest.approx(n_units)
    np.testing.assert_allclose(est.em_counts.sum(axis=1), [50.0, 50.0], atol=1e-6)


def test_enabled_gdna_component_absorbs_likelihood_mass_without_prior_count():
    n_units = 80
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -1.0), gdna_log_lik=0.0)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert total_gdna > 0.8 * n_units
    assert locus_gdna[0] == pytest.approx(total_gdna)
    assert locus_rna[0] + locus_gdna[0] == pytest.approx(n_units)


def test_gdna_effective_length_downweights_gdna_component():
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)

    est_short = _estimator(1, mode="map")
    _total_short, _rna_short, gdna_short = est_short.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    est_long = _estimator(1, mode="map")
    _total_long, _rna_long, gdna_long = est_long.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([100.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert gdna_long[0] < gdna_short[0]


def test_gdna_em_llr_bias_favors_gdna_assignment():
    # An ambiguous fragment (transcript and gDNA equally likely, equal eff-lens):
    # neutral (bias=0) splits ~evenly; a positive gDNA LLR bias shifts mass toward
    # gDNA (the FP-aversion knob — fewer gDNA->RNA leaks, more RNA->gDNA siphons);
    # a negative bias shifts toward RNA. Monotone in the bias.
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)
    kw = dict(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )
    _t, _r, gdna_neg = _estimator(1, gdna_em_llr_bias=-2.0).run_batch_locus_em_partitioned(**kw)
    _t, _r, gdna_neutral = _estimator(1, gdna_em_llr_bias=0.0).run_batch_locus_em_partitioned(**kw)
    _t, _r, gdna_pos = _estimator(1, gdna_em_llr_bias=2.0).run_batch_locus_em_partitioned(**kw)

    assert gdna_neg[0] < gdna_neutral[0] < gdna_pos[0]
    # neutral splits an exactly-tied 1-transcript-vs-gDNA fragment ~50/50
    assert gdna_neutral[0] == pytest.approx(n_units / 2, rel=0.05)
    # a +2-nat (~7.4:1 odds) bias pushes the gDNA share well above half
    assert gdna_pos[0] > 0.8 * n_units


def test_aggregate_rna_prior_reduces_gdna_share_without_isoform_floor():
    n_units = 100
    partition = _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0)

    est_unprior = _estimator(1, mode="map")
    _total_unprior, _rna_unprior, gdna_unprior = est_unprior.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    est_rna = _estimator(1, mode="map")
    _total_rna, _rna_rna, gdna_rna = est_rna.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        rna_prior_count=np.array([10.0], dtype=np.float64),
        index=None,
        gdna_eff_len=np.array([1.0], dtype=np.float64),
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert gdna_rna[0] < gdna_unprior[0]


def test_grouped_priors_inactive_without_structural_gdna_candidate():
    n_units = 40
    est = _estimator(1, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[
            _partition(n_units=n_units, log_liks=(0.0,), gdna_log_lik=0.0, is_spliced=True)
        ],
        locus_transcript_indices=[np.array([0], dtype=np.int32)],
        gdna_prior_count=np.array([100.0], dtype=np.float64),
        rna_prior_count=np.array([100.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([1], dtype=np.uint8),
    )

    assert total_gdna == pytest.approx(0.0)
    assert locus_gdna[0] == pytest.approx(0.0)
    assert locus_rna[0] == pytest.approx(n_units)


def test_assignment_outputs_follow_partition_units():
    n_units = 5
    est = _estimator(2, mode="map")
    result = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(0.0, -4.0))],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([0.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
        emit_assignments=True,
    )

    total_gdna, _locus_rna, _locus_gdna, winner_tid, winner_post, n_candidates = result
    assert total_gdna == pytest.approx(0.0)
    np.testing.assert_array_equal(winner_tid, np.zeros(n_units, dtype=np.int32))
    assert np.all(winner_post > 0.98)
    np.testing.assert_array_equal(n_candidates, np.full(n_units, 2, dtype=np.int16))


# ---------------------------------------------------------------------------
# Structural gDNA eligibility: ``gdna_prior_count`` no longer gates the gDNA
# component, and the compatibility ``enable_gdna`` array is not a modeling gate
# — native derives candidate availability from the partition itself.
# ---------------------------------------------------------------------------


def test_compat_enable_false_does_not_disable_structural_candidate():
    """A compatibility ``enable_gdna=False`` input is ignored by native v3."""
    n_units = 50
    est = _estimator(2, mode="map")
    total_gdna, _rna, _g = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -1.0), gdna_log_lik=-0.5)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([100.0], dtype=np.float64),
        index=None,
        enable_gdna=np.array([0], dtype=np.uint8),
    )

    assert total_gdna > 0.5 * n_units, (
        "native should derive gDNA availability from finite unspliced candidates, "
        f"not the compatibility enable_gdna array; got total_gdna={total_gdna}"
    )


def test_default_enable_gdna_inferred_from_partition():
    """When ``enable_gdna`` is None, the wrapper computes it from the partition
    (any unspliced unit with a finite gDNA log-lik ⇒ enabled).
    """
    n_units = 20
    locus_t_lists = [np.array([0, 1], dtype=np.int32)]

    # All-spliced partition ⇒ no per-unit gDNA candidate ⇒ enable=0.
    spliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=True
    )
    # Unspliced partition with finite gDNA log-liks ⇒ enable=1.
    unspliced_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-0.5, is_spliced=False
    )
    # Unspliced but non-finite gDNA log-lik ⇒ enable=0.
    nogdna_part = _partition(
        n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-np.inf, is_spliced=False
    )

    # All three with gdna_prior_count=0 and enable_gdna omitted (inferred).
    # Only the unspliced+finite case produces gDNA assignments.
    g_spl, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [spliced_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )
    g_uns, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [unspliced_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )
    g_no, _, _ = _estimator(2, mode="map").run_batch_locus_em_partitioned(
        [nogdna_part],
        locus_t_lists,
        np.zeros(1),
        index=None,
    )

    assert g_spl == 0.0, "spliced partition has no gDNA candidates"
    assert g_no == 0.0, "non-finite gDNA log-liks ⇒ no gDNA candidates"
    assert g_uns > 0.0, (
        "unspliced+finite gDNA log-lik must enable component even when gdna_prior_count == 0"
    )


def test_positive_gdna_prior_produces_finite_outputs():
    """A positive gDNA prior count produces finite, conserved outputs."""
    n_units = 30
    est = _estimator(2, mode="map")
    total_gdna, locus_rna, locus_gdna = est.run_batch_locus_em_partitioned(
        partition_tuples=[_partition(n_units=n_units, log_liks=(-1.0, -2.0), gdna_log_lik=-1.0)],
        locus_transcript_indices=[np.array([0, 1], dtype=np.int32)],
        gdna_prior_count=np.array([5.0], dtype=np.float64),
        index=None,
    )

    # Outputs must be finite and sum to total fragments.
    assert np.isfinite(total_gdna)
    assert np.isfinite(locus_rna).all() and np.isfinite(locus_gdna).all()
    assert total_gdna >= 0.0
    # Conservation up to assignment fractional remainder.
    assigned = float(locus_rna[0]) + float(locus_gdna[0])
    assert assigned == pytest.approx(30.0, rel=1e-9)


# ── SYNTHETIC NASCENT ENTITIES GET NO RNA PRIOR ──────────────────────────────────────────────────
#
# A synthetic nascent entity is a shadow span the INDEX manufactured; no annotation asserts it exists.
# The null hypothesis is therefore that it is ABSENT, and it earns mass only from fragments the data
# cannot explain any other way. These gates pin that behaviour and, just as importantly, pin the two
# things it must NOT do: kill an entity the data supports, and disturb a locus that has none.


class _StubIndex:
    """The minimum `run_batch_locus_em_partitioned` reads: `t_df["is_synthetic"]`."""

    def __init__(self, flags):
        import pandas as pd

        self.t_df = pd.DataFrame({"is_synthetic": np.asarray(flags, dtype=bool)})


def _run(est, partition, t_idx, *, index=None, rna_prior=0.0, gdna_prior=0.0, enable_gdna=1):
    return est.run_batch_locus_em_partitioned(
        partition_tuples=[partition],
        locus_transcript_indices=[np.asarray(t_idx, dtype=np.int32)],
        gdna_prior_count=np.array([gdna_prior], dtype=np.float64),
        rna_prior_count=np.array([rna_prior], dtype=np.float64),
        index=index,
        enable_gdna=np.array([enable_gdna], dtype=np.uint8),
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_locus_with_no_synthetic_component_is_BIT_IDENTICAL(mode):
    """⛔ The control. The rule must be invisible where it does not apply — and BIT-identical, not
    close: the C++ skips the annotated/synthetic split entirely when the mask is empty, and that is
    what makes every other number here attributable (TRAPS: byte-identity-gate)."""
    outs = []
    for index in (None, _StubIndex([False, False])):
        est = _estimator(2, mode=mode)
        _run(est, _partition(n_units=100, log_liks=(0.0, -0.5), gdna_log_lik=-1.0), [0, 1],
             index=index, rna_prior=50.0, gdna_prior=20.0)
        outs.append(est.em_counts.sum(axis=1).copy())
    np.testing.assert_array_equal(outs[0], outs[1])


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_SHARED_ONLY_synthetic_entity_LOSES_mass_to_the_annotated_one(mode):
    """⭐ The zombie. Every fragment it holds is equally well explained by the annotated transcript,
    so it has no evidence of its own and the prior no longer props it up."""
    part = _partition(n_units=200, log_liks=(0.0, 0.0), gdna_log_lik=-2.0)
    counts = {}
    for tag, index in (("shipped", None), ("gated", _StubIndex([False, True]))):
        est = _estimator(2, mode=mode)
        _run(est, part, [0, 1], index=index, rna_prior=100.0, gdna_prior=10.0)
        counts[tag] = est.em_counts.sum(axis=1).copy()
    assert counts["gated"][1] < counts["shipped"][1], (
        f"the synthetic component did not lose mass: {counts}"
    )
    assert counts["gated"][0] > counts["shipped"][0], (
        f"the annotated component did not gain what the synthetic lost: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_a_synthetic_entity_the_DATA_SUPPORTS_still_survives(mode):
    """⛔⛔ THE GATE THAT STOPS THIS BECOMING 'KILL ALL NASCENT RNA'. A rule that zeroes everything
    passes the zombie test above; only this one separates the two. The synthetic component is the
    strictly better explanation here, and it must keep the mass."""
    est = _estimator(2, mode=mode)
    _run(est, _partition(n_units=200, log_liks=(-8.0, 0.0), gdna_log_lik=-8.0), [0, 1],
         index=_StubIndex([False, True]), rna_prior=100.0, gdna_prior=10.0)
    counts = est.em_counts.sum(axis=1)
    assert counts[1] > 0.9 * counts.sum(), (
        f"a synthetic entity with decisive likelihood support was suppressed anyway: {counts}"
    )


@pytest.mark.parametrize("mode", ["map", "vbem"])
def test_an_ALL_SYNTHETIC_locus_takes_NO_rna_prior(mode):
    """⚠ The degenerate locus of the derivation: no annotated component means no eligible recipient,
    so the RNA prior is 0 and the pool must outcompete gDNA unaided. It must not silently fall back
    to handing the prior to the synthetic components after all."""
    part = _partition(n_units=100, log_liks=(0.0, 0.0), gdna_log_lik=0.0)
    res = {}
    for tag, rna_prior in (("no_prior", 0.0), ("big_prior", 500.0)):
        est = _estimator(2, mode=mode)
        _total, _rna, gdna = _run(est, part, [0, 1], index=_StubIndex([True, True]),
                                  rna_prior=rna_prior, gdna_prior=10.0)
        res[tag] = (est.em_counts.sum(axis=1).copy(), float(gdna[0]))
    np.testing.assert_allclose(res["no_prior"][0], res["big_prior"][0], rtol=1e-9, atol=1e-9)
    assert res["no_prior"][1] == pytest.approx(res["big_prior"][1], rel=1e-9), (
        "the RNA prior reached an all-synthetic locus and moved the gDNA split"
    )


def test_is_synthetic_is_read_NOT_is_nrna():
    """⛔ A single-exon annotated transcript carries `is_nrna = True` and is simultaneously the
    nascent and the mature form of a REAL gene. It must keep its prior. The estimator must key on
    `is_synthetic` alone, so an index carrying only `is_nrna` changes nothing."""
    import pandas as pd

    class _NrnaOnlyIndex:
        def __init__(self):
            self.t_df = pd.DataFrame({"is_nrna": [False, True]})

    est = AbundanceEstimator(num_transcripts=2, em_config=EMConfig(mode="map"))
    np.testing.assert_array_equal(est._t_is_synthetic(_NrnaOnlyIndex(), 2), np.zeros(0, np.uint8))
    np.testing.assert_array_equal(est._t_is_synthetic(_StubIndex([False, False]), 2),
                                  np.zeros(0, np.uint8))
    np.testing.assert_array_equal(est._t_is_synthetic(_StubIndex([False, True]), 2),
                                  np.array([0, 1], np.uint8))


# ⛔⛔ THE ONE INVARIANT THIS DESIGN RESTS ON HAS NO TEST, AND CANNOT HAVE ONE FROM PYTHON TODAY.
#
# `apply_grouped_prior_update` guarantees, FOR A GIVEN raw_counts VECTOR, that the RNA components sum
# to `rna_count + rna_prior` — so withholding the prior from synthetic components redistributes mass
# strictly WITHIN the RNA pool. It is asserted only in a C++ comment.
#
# ⚠ It is a PER-M-STEP algebraic identity, NOT an end-to-end one, and conflating the two has now
# produced a wrong test twice. End to end the gDNA total legitimately DEPENDS on `rna_prior`: setting
# the gDNA:RNA split is the prior's whole purpose (`gdna_total = gdna_count + gdna_prior`,
# `rna_total = rna_count + rna_prior`), and a larger RNA prior shifts theta, hence the E-step, hence
# the next iteration's `gdna_count`. A test asserting "the gDNA total is invariant to rna_prior" is
# therefore asserting something FALSE BY DESIGN — verified: it fails on the shipped no-synthetic
# configuration too.
#
# ⭐ The function is `static` in em_solver.cpp, so nothing can call it directly. Exposing it behind a
# test-only binding — the executable-specification pattern this repo already uses for the accumulator
# — is the way to gate the identity, and it is a prerequisite for the per-transcript prior work
# (a per-transcript vector makes the identity harder, not easier, to hold).
