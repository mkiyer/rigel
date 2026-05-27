"""Tests for :mod:`rigel.calibration.strand_deconv` (Phase 3 of v5 plan)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.special import betaln, gammaln, logsumexp

from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.scan_payload import FL_HIST_N_BINS, CalibrationScanPayload
from rigel.calibration.signature import (
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_CONTAINED,
    N_CHANNELS,
    N_FL_POOLS,
    N_SIGNATURES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
    pack_signature,
)
from rigel.calibration.strand_deconv import (
    DNA_STRAND_UNIFORM_PRIOR,
    FLAG_INELIGIBLE,
    FLAG_KAPPA_FALLBACK,
    FLAG_NEAR_UNSTRANDED,
    FLAG_RELIABILITY_APPROX,
    MAX_EXACT_POSTERIOR_N,
    _log_normal_interval_prob,
    build_strand_region_counts,
    deconvolve_regions_by_strand,
    estimate_kappa_d,
    screen_no_rna_exons,
    strand_log_likelihood_d_grid_minor_beta_binom,
)
from rigel.calibration.strand_summary import StrandSummary


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


def _summary(p_r1_sense: float, n_observations: int) -> StrandSummary:
    if n_observations == 0:
        return StrandSummary.uninformative()
    n_same = int(round(float(p_r1_sense) * int(n_observations)))
    assert n_same / int(n_observations) == pytest.approx(float(p_r1_sense), abs=1e-12)
    return StrandSummary(
        p_r1_sense=float(p_r1_sense),
        n_observations=int(n_observations),
        n_same=n_same,
        n_opposite=int(n_observations) - n_same,
    )


def _build_arrays(
    rows: list[tuple[str, int, int, int]],
    channel_assignments: dict[int, dict[tuple[int, int, int], float]],
) -> tuple[RegionArrays, PayloadArrays]:
    """rows = list of (ref, start, end, signature)
    channel_assignments = {row_idx: {(compartment, splice, strand): mass}}
    """
    df = pd.DataFrame(rows, columns=["ref_name", "start", "end", "signature"])
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    df.index = df["region_id"].to_numpy()

    n_regions = len(rows)
    region_counts = np.zeros((n_regions, N_CHANNELS), dtype=np.float32)
    for row_idx, channels in channel_assignments.items():
        for (comp, splice, strand), mass in channels.items():
            region_counts[row_idx, channel_index(comp, splice, strand)] = mass

    signature_mass = np.zeros(N_SIGNATURES, dtype=np.float64)
    for r, sig in enumerate(df["signature"].to_numpy()):
        signature_mass[int(sig)] += float(region_counts[r].sum(dtype=np.float64))

    fl_pool_mass = np.zeros((N_FL_POOLS, FL_HIST_N_BINS), dtype=np.float64)
    payload = CalibrationScanPayload(
        region_counts=region_counts,
        channel_mass=region_counts.sum(axis=0, dtype=np.float64),
        signature_mass=signature_mass,
        fl_pool_mass=fl_pool_mass,
        fl_pool_total=fl_pool_mass.sum(axis=1),
        n_observed=int(region_counts.sum(dtype=np.float64)),
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_excluded_strand_ambig=0,
        n_unobserved=0,
        n_unannotated_ref=0,
        n_fl_unavailable=0,
        resolver_splicing_anchor_tolerance=0,
        n_regions=n_regions,
    )
    region_arrays = RegionArrays.from_region_df(df, {"chr1": 0})
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    return region_arrays, payload_arrays


def _single_region(
    signature: int,
    pos: float,
    neg: float,
    *,
    spliced_pos: float = 0.0,
    spliced_neg: float = 0.0,
) -> tuple[RegionArrays, PayloadArrays]:
    rows = [("chr1", 0, 1000, signature)]
    channels: dict[tuple[int, int, int], float] = {}
    if pos:
        channels[(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS)] = pos
    if neg:
        channels[(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG)] = neg
    if spliced_pos:
        channels[(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_POS)] = spliced_pos
    if spliced_neg:
        channels[(COMPARTMENT_CONTAINED, SPLICE_SPLICED, CHANNEL_STRAND_NEG)] = spliced_neg
    return _build_arrays(rows, {0: channels})


# ---------------------------------------------------------------------------
# build_strand_region_counts
# ---------------------------------------------------------------------------


def test_observations_ts_pos_and_ts_neg_basic():
    rows = [
        ("chr1", 0, 1000, pack_signature(exon_pos=True)),  # TS_POS
        ("chr1", 1000, 2000, pack_signature(intron_neg=True)),  # TS_NEG
        ("chr1", 2000, 3000, pack_signature(exon_pos=True, exon_neg=True)),  # TS_AMBIG
        ("chr1", 3000, 4000, pack_signature()),  # TS_NONE
    ]
    channels = {
        0: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 80.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 20.0,
        },
        1: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 30.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 70.0,
        },
        2: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 40.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 40.0,
        },
        3: {(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 25.0},
    }
    region_arrays, payload_arrays = _build_arrays(rows, channels)

    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.9)

    # Region table is sorted by (ref, start); these inputs already are.
    # TS_POS: sense=POS=80, antisense=NEG=20.
    assert counts.k_sense[0] == pytest.approx(80.0)
    assert counts.k_antisense[0] == pytest.approx(20.0)
    assert counts.n_total[0] == pytest.approx(100.0)
    assert counts.eligible[0]
    # TS_NEG: sense=NEG=70, antisense=POS=30.
    assert counts.k_sense[1] == pytest.approx(70.0)
    assert counts.k_antisense[1] == pytest.approx(30.0)
    assert counts.eligible[1]
    # TS_AMBIG and TS_NONE: ineligible, zero folded counts.
    assert not counts.eligible[2]
    assert counts.k_sense[2] == 0.0 and counts.k_antisense[2] == 0.0
    assert not counts.eligible[3]
    assert counts.k_sense[3] == 0.0 and counts.k_antisense[3] == 0.0
    assert counts.p_r1_sense == pytest.approx(0.9)


def test_protocol_symmetry():
    # Same physical distribution observed under R1-sense (p=0.9) vs R1-antisense (p=0.1).
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=60.0, neg=40.0
    )
    counts_sense = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.9)
    est_sense = deconvolve_regions_by_strand(
        counts_sense,
        kappa_d=10.0,
        strand_summary=_summary(0.9, 10000),
    )

    # For the antisense protocol the SAME molecules produce reflected reads:
    # POS and NEG channel counts swap.
    region_arrays2, payload_arrays2 = _single_region(
        pack_signature(exon_pos=True), pos=40.0, neg=60.0
    )
    counts_anti = build_strand_region_counts(region_arrays2, payload_arrays2, p_r1_sense=0.1)
    est_anti = deconvolve_regions_by_strand(
        counts_anti,
        kappa_d=10.0,
        strand_summary=_summary(0.1, 10000),
    )

    assert est_sense.mean_count[0] == pytest.approx(est_anti.mean_count[0], abs=1e-3)
    assert est_sense.rna_lower_count[0] == pytest.approx(est_anti.rna_lower_count[0], abs=1.0)
    assert est_sense.upper_count[0] == pytest.approx(est_anti.upper_count[0], abs=1.0)


# ---------------------------------------------------------------------------
# Deconvolution
# ---------------------------------------------------------------------------


def test_near_unstranded_is_conservative():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=60.0, neg=40.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.5)
    est = deconvolve_regions_by_strand(
        counts,
        kappa_d=10.0,
        strand_summary=StrandSummary.uninformative(),
    )

    assert (est.flags[0] & FLAG_NEAR_UNSTRANDED) != 0
    assert est.rna_lower_count[0] == 0.0
    assert est.mean_count[0] == pytest.approx(est.n_total[0])
    assert est.upper_count[0] == pytest.approx(est.n_total[0])
    assert est.precision[0] == 0.0


def test_exact_posterior_handles_perfect_strand_specificity():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=80.0, neg=20.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=1.0)

    est = deconvolve_regions_by_strand(counts, kappa_d=1.0e4, strand_summary=_summary(1.0, 1000))

    assert np.isfinite(est.mean_count[0])
    assert np.isfinite(est.rna_lower_count[0])
    assert est.mean_count[0] < est.n_total[0]
    assert est.rna_lower_count[0] > 0.0


def test_fractional_exact_outputs_are_clamped_to_n_total():
    # The exact posterior rounds fractional accumulator mass to an integer grid.
    # Output summaries still need to stay on the original fractional-count scale.
    region_arrays, payload_arrays = _single_region(pack_signature(exon_pos=True), pos=1.6, neg=0.0)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)

    est = deconvolve_regions_by_strand(
        counts,
        kappa_d=1.0e6,
        strand_summary=_summary(0.99, 10000),
    )

    assert 0.0 <= est.mean_count[0] <= est.n_total[0]
    assert 0.0 <= est.upper_count[0] <= est.n_total[0]
    assert 0.0 <= est.rna_lower_count[0] <= est.n_total[0]


def _manual_log_beta_binom(k: np.ndarray, n: int, alpha: float, beta: float) -> np.ndarray:
    k_arr = np.asarray(k, dtype=np.float64)
    if n <= 0:
        return np.where(k_arr == 0, 0.0, -np.inf)
    log_c = gammaln(n + 1) - gammaln(k_arr + 1) - gammaln(n - k_arr + 1)
    return log_c + betaln(k_arr + alpha, n - k_arr + beta) - betaln(alpha, beta)


def test_minor_beta_binom_mixed_likelihood_matches_direct_convolution():
    n = 20
    k_minor = 4
    d_grid = np.array([0, 3, 8, 20], dtype=np.int64)
    kappa = 9.0
    alpha_q = 2.0
    beta_q = 80.0

    observed = strand_log_likelihood_d_grid_minor_beta_binom(
        k_minor,
        n,
        d_grid,
        kappa_d=kappa,
        minor_rate_alpha=alpha_q,
        minor_rate_beta=beta_q,
    )

    expected = []
    dna_alpha = kappa / 2.0
    for d in d_grid:
        r = n - int(d)
        j_lo = max(0, k_minor - int(d))
        j_hi = min(r, k_minor)
        js = np.arange(j_lo, j_hi + 1)
        terms = _manual_log_beta_binom(js, r, alpha_q, beta_q) + _manual_log_beta_binom(
            k_minor - js,
            int(d),
            dna_alpha,
            dna_alpha,
        )
        expected.append(float(logsumexp(terms)))

    np.testing.assert_allclose(observed, expected, atol=1e-12)


def test_beta_binom_slab_mean_uses_same_likelihood_as_log_p_mixed():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=24.0, neg=6.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.95, 1000)
    kappa = 12.0

    estimate = deconvolve_regions_by_strand(counts, kappa_d=kappa, strand_summary=summary)

    n = 30
    k_minor = 6
    d_grid = np.arange(n + 1, dtype=np.int64)
    log_l = strand_log_likelihood_d_grid_minor_beta_binom(
        k_minor,
        n,
        d_grid,
        kappa_d=kappa,
        minor_rate_alpha=summary.minor_rate_alpha,
        minor_rate_beta=summary.minor_rate_beta,
    )
    log_norm = float(logsumexp(log_l))
    posterior = np.exp(log_l - log_norm)
    expected_mean = float(np.sum(d_grid * posterior))
    expected_log_p_mixed = log_norm - np.log(n + 1)

    assert estimate.mean_count[0] == pytest.approx(expected_mean, abs=1e-5)
    assert estimate.log_p_mixed[0] == pytest.approx(expected_log_p_mixed, abs=1e-5)


def test_reliability_is_vectorized_continuous_and_not_thresholded():
    summary = _summary(0.95, 1000)
    weights = []
    for k_minor in range(1, 8):
        region_arrays, payload_arrays = _single_region(
            pack_signature(exon_pos=True), pos=50.0 - k_minor, neg=float(k_minor)
        )
        counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
        estimate = deconvolve_regions_by_strand(counts, kappa_d=20.0, strand_summary=summary)
        weights.append(float(estimate.reliability[0]))

    assert all(0.0 < weight < 1.0 for weight in weights)
    assert all(next_weight >= weight for weight, next_weight in zip(weights, weights[1:]))
    assert max(np.diff(weights)) < 0.75


def test_near_unstranded_beta_binom_reliability_is_inactive():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=25.0, neg=25.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.5001)
    summary = _summary(0.5001, 10000)

    estimate = deconvolve_regions_by_strand(counts, kappa_d=10.0, strand_summary=summary)

    assert estimate.reliability[0] == pytest.approx(0.0)
    assert (estimate.flags[0] & FLAG_NEAR_UNSTRANDED) != 0


def test_pure_rna_beta_binom_reliability_is_low():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=99.0, neg=1.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)
    summary = _summary(0.99, 10_000)

    estimate = deconvolve_regions_by_strand(counts, kappa_d=20.0, strand_summary=summary)

    assert 0.0 < estimate.reliability[0] < 0.1


def test_strong_mixed_beta_binom_reliability_is_high():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=50.0, neg=50.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)
    summary = _summary(0.99, 10_000)

    estimate = deconvolve_regions_by_strand(counts, kappa_d=100.0, strand_summary=summary)

    assert estimate.reliability[0] > 0.99


def test_small_strand_training_set_lowers_reliability():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=70.0, neg=30.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)
    loose_summary = _summary(0.99, 100)
    tight_summary = _summary(0.99, 10_000)

    loose = deconvolve_regions_by_strand(counts, kappa_d=50.0, strand_summary=loose_summary)
    tight = deconvolve_regions_by_strand(counts, kappa_d=50.0, strand_summary=tight_summary)

    assert loose.reliability[0] < tight.reliability[0]


def test_large_count_beta_binom_reliability_uses_approximation():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=180.0, neg=70.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)
    summary = _summary(0.99, 10_000)

    estimate = deconvolve_regions_by_strand(counts, kappa_d=50.0, strand_summary=summary)

    assert counts.n_total[0] > MAX_EXACT_POSTERIOR_N
    assert np.isfinite(estimate.reliability[0])
    assert (estimate.flags[0] & FLAG_RELIABILITY_APPROX) != 0


def test_log_normal_interval_prob_remains_finite_in_deep_tails():
    log_probs = _log_normal_interval_prob(
        np.array([-9.0, 8.0], dtype=np.float64),
        np.array([-8.0, 9.0], dtype=np.float64),
    )

    assert np.all(np.isfinite(log_probs))
    assert np.all(log_probs < 0.0)


# ---------------------------------------------------------------------------
# Exon screen
# ---------------------------------------------------------------------------


def test_exon_screen_accepts_balanced_no_spliced_exon():
    # With a tight (high-kappa) gDNA strand prior, a balanced 100/100 unspliced
    # exon should be accepted at confidence 0.95.
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=100.0, neg=100.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(counts, region_arrays, payload_arrays, kappa_d_seed=1.0e4)
    assert bool(accepted[0])


def test_exon_screen_rejects_clear_rna_exon():
    region_arrays, payload_arrays = _single_region(pack_signature(exon_pos=True), pos=95.0, neg=5.0)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(counts, region_arrays, payload_arrays, kappa_d_seed=10.0)
    assert not bool(accepted[0])


def test_exon_screen_rejects_spliced_exon():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True),
        pos=50.0,
        neg=50.0,
        spliced_pos=5.0,
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(counts, region_arrays, payload_arrays, kappa_d_seed=10.0)
    assert not bool(accepted[0])


def test_exon_screen_requires_both_strand_channels():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=100.0, neg=0.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(counts, region_arrays, payload_arrays, kappa_d_seed=10.0)
    assert not bool(accepted[0])


# ---------------------------------------------------------------------------
# kappa_d refit
# ---------------------------------------------------------------------------


def _seed_plus_exons_arrays() -> tuple[RegionArrays, PayloadArrays]:
    """Four tightly balanced seed regions (→ large kappa) plus three balanced exons."""
    rows = [
        ("chr1", 0, 500, pack_signature()),
        ("chr1", 500, 1000, pack_signature()),
        ("chr1", 1000, 1500, pack_signature(intron_pos=True)),
        ("chr1", 1500, 2000, pack_signature(intron_neg=True)),
        ("chr1", 2000, 2500, pack_signature(exon_pos=True)),
        ("chr1", 2500, 3000, pack_signature(exon_pos=True)),
        ("chr1", 3000, 3500, pack_signature(exon_neg=True)),
    ]
    # Seeds chosen with residual variance only modestly above binomial, so
    # MoM kappa is large (tight strand balance). Exons are similarly balanced.
    channels = {
        0: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 52.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 48.0,
        },
        1: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 49.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 51.0,
        },
        2: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 51.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 49.0,
        },
        3: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 48.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 52.0,
        },
        4: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 100.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 100.0,
        },
        5: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 95.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 105.0,
        },
        6: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 103.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 97.0,
        },
    }
    return _build_arrays(rows, channels)


def test_kappa_refit_uses_accepted_exons():
    region_arrays, payload_arrays = _seed_plus_exons_arrays()
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.95, 10_000)

    est = estimate_kappa_d(region_arrays, payload_arrays, counts, summary)
    assert est.n_seed_regions == 4
    assert est.n_exon_self_training >= 1
    # Refit balance.n_regions should reflect seeds + accepted exons.
    assert est.balance.n_regions == est.n_seed_regions + est.n_exon_self_training


def test_kappa_soft_fallback_is_not_reported_as_uniform_fallback():
    region_arrays, payload_arrays = _seed_plus_exons_arrays()
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.95, 10_000)

    est = estimate_kappa_d(region_arrays, payload_arrays, counts, summary)

    assert est.balance.fallback_used
    assert est.balance.fallback_reason == "residual variance <= binomial expectation"
    assert not est.fallback_used
    assert est.kappa > DNA_STRAND_UNIFORM_PRIOR


def test_kappa_hard_fallback_skips_exon_self_training_even_with_exons():
    rows = [
        ("chr1", 0, 500, pack_signature()),
        ("chr1", 500, 1000, pack_signature(exon_pos=True)),
    ]
    channels = {
        0: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 10.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 10.0,
        },
        1: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 100.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 100.0,
        },
    }
    region_arrays, payload_arrays = _build_arrays(rows, channels)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.95, 10_000)

    est = estimate_kappa_d(region_arrays, payload_arrays, counts, summary)

    assert est.fallback_used
    assert est.kappa == pytest.approx(DNA_STRAND_UNIFORM_PRIOR)
    assert est.n_seed_regions == 1
    assert est.n_exon_self_training == 0
    assert est.balance.n_regions == 1


def test_kappa_fallback_when_unidentifiable():
    # Single intergenic region only — below MIN_REGIONS_FOR_STRAND_MOM,
    # forcing the symmetric-beta-binomial seed estimate to fall back.
    rows = [("chr1", 0, 500, pack_signature())]
    channels = {
        0: {
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_POS): 10.0,
            (COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, CHANNEL_STRAND_NEG): 10.0,
        },
    }
    region_arrays, payload_arrays = _build_arrays(rows, channels)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.95, 10_000)

    est = estimate_kappa_d(region_arrays, payload_arrays, counts, summary)
    assert est.fallback_used
    assert est.kappa == pytest.approx(DNA_STRAND_UNIFORM_PRIOR)
    assert est.n_exon_self_training == 0

    # Plumbing through deconvolution should propagate the fallback flag.
    region_arrays2, payload_arrays2 = _single_region(
        pack_signature(exon_pos=True), pos=80.0, neg=20.0
    )
    counts2 = build_strand_region_counts(region_arrays2, payload_arrays2, p_r1_sense=0.95)
    deconv = deconvolve_regions_by_strand(
        counts2,
        kappa_d=est.kappa,
        kappa_d_fallback_used=True,
        strand_summary=_summary(0.95, 10_000),
    )
    assert (deconv.flags[0] & FLAG_KAPPA_FALLBACK) != 0


def test_estimate_kappa_requires_matching_strand_summary():
    region_arrays, payload_arrays = _seed_plus_exons_arrays()
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = _summary(0.90, 10_000)

    with pytest.raises(ValueError, match="strand_summary.p_r1_sense"):
        estimate_kappa_d(
            region_arrays,
            payload_arrays,
            counts,
            summary,
        )


# ---------------------------------------------------------------------------
# Ineligible / sanity
# ---------------------------------------------------------------------------


def test_ineligible_regions_carry_zero_estimates():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True, exon_neg=True), pos=50.0, neg=50.0
    )  # TS_AMBIG -> ineligible
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.9)
    est = deconvolve_regions_by_strand(counts, kappa_d=5.0, strand_summary=_summary(0.9, 10000))
    assert (est.flags[0] & FLAG_INELIGIBLE) != 0
    assert est.rna_lower_count[0] == 0.0
    assert est.precision[0] == 0.0
