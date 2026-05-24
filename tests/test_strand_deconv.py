"""Tests for :mod:`rigel.calibration.strand_deconv` (Phase 3 of v5 plan)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

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
    MAX_EXACT_POSTERIOR_N,
    _exact_posterior_R,
    _summarize_exact,
    _summarize_normal,
    build_strand_region_counts,
    deconvolve_regions_by_strand,
    estimate_kappa_d,
    screen_no_rna_exons,
)
from rigel.calibration.strand_summary import StrandSummary


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------


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
    est_sense = deconvolve_regions_by_strand(counts_sense, kappa_d=10.0, rna_lower_confidence=0.95)

    # For the antisense protocol the SAME molecules produce reflected reads:
    # POS and NEG channel counts swap.
    region_arrays2, payload_arrays2 = _single_region(
        pack_signature(exon_pos=True), pos=40.0, neg=60.0
    )
    counts_anti = build_strand_region_counts(region_arrays2, payload_arrays2, p_r1_sense=0.1)
    est_anti = deconvolve_regions_by_strand(counts_anti, kappa_d=10.0, rna_lower_confidence=0.95)

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
    est = deconvolve_regions_by_strand(counts, kappa_d=10.0, rna_lower_confidence=0.95)

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

    est = deconvolve_regions_by_strand(counts, kappa_d=1.0e4, rna_lower_confidence=0.95)

    assert np.isfinite(est.mean_count[0])
    assert np.isfinite(est.rna_lower_count[0])
    assert est.mean_count[0] < est.n_total[0]
    assert est.rna_lower_count[0] > 0.0


def test_fractional_exact_outputs_are_clamped_to_n_total():
    # The exact posterior rounds fractional accumulator mass to an integer grid.
    # Output summaries still need to stay on the original fractional-count scale.
    region_arrays, payload_arrays = _single_region(pack_signature(exon_pos=True), pos=1.6, neg=0.0)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.99)

    est = deconvolve_regions_by_strand(counts, kappa_d=1.0e6, rna_lower_confidence=0.50)

    assert 0.0 <= est.mean_count[0] <= est.n_total[0]
    assert 0.0 <= est.upper_count[0] <= est.n_total[0]
    assert 0.0 <= est.rna_lower_count[0] <= est.n_total[0]


def test_rna_lower_confidence_monotone_in_bound():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=80.0, neg=20.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)

    bounds = [
        deconvolve_regions_by_strand(counts, kappa_d=5.0, rna_lower_confidence=c).rna_lower_count[0]
        for c in (0.50, 0.80, 0.95, 0.99)
    ]
    # As we demand higher confidence the lower bound on R must not increase.
    for a, b in zip(bounds, bounds[1:]):
        assert b <= a + 1e-6, f"non-monotone bound: {bounds}"


def test_exact_and_normal_agree_at_boundary():
    # At n == MAX_EXACT_POSTERIOR_N, with a tight gDNA strand-balance (kappa large
    # → BetaBinomial ≈ Binomial(0.5)) and a high-contrast protocol, both the
    # exact and the normal-approximation summaries should agree closely. Moderate
    # kappa values produce strongly asymmetric posteriors that the normal approx
    # cannot capture — a known limitation of evaluating Var at R̂ (see Phase 3 log).
    n = MAX_EXACT_POSTERIOR_N
    p = 0.99
    kappa = 1.0e4
    # Mix so that ~half the reads are RNA and ~half are gDNA. Under high contrast
    # k_sense ≈ 0.5*D + p*R so for a 50/50 mix k_sense ≈ 0.5*n*(0.5+p).
    r_true = 0.5 * n
    k_sense = int(round(0.5 * (n - r_true) + p * r_true))

    posterior = _exact_posterior_R(k_sense, n, kappa_d=kappa, p_r1_sense=p)
    r_hat_e, r_lower_e, _r_upper_e, sd_e = _summarize_exact(
        posterior, rna_lower_confidence=0.95, n_int=n
    )
    r_hat_n, r_lower_n, _r_upper_n, sd_n = _summarize_normal(
        float(k_sense),
        float(n),
        kappa_d=kappa,
        p_r1_sense=p,
        rna_lower_confidence=0.95,
    )

    assert abs(r_hat_e - r_hat_n) < 0.05 * n
    assert abs(sd_e - sd_n) < 0.20 * max(sd_n, 1.0)
    assert abs(r_lower_e - r_lower_n) < 0.10 * n


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
    accepted = screen_no_rna_exons(
        counts, region_arrays, payload_arrays, kappa_d_seed=1.0e4, rna_lower_confidence=0.95
    )
    assert bool(accepted[0])


def test_exon_screen_rejects_clear_rna_exon():
    region_arrays, payload_arrays = _single_region(pack_signature(exon_pos=True), pos=95.0, neg=5.0)
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(
        counts, region_arrays, payload_arrays, kappa_d_seed=10.0, rna_lower_confidence=0.95
    )
    assert not bool(accepted[0])


def test_exon_screen_rejects_spliced_exon():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True),
        pos=50.0,
        neg=50.0,
        spliced_pos=5.0,
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(
        counts, region_arrays, payload_arrays, kappa_d_seed=10.0, rna_lower_confidence=0.95
    )
    assert not bool(accepted[0])


def test_exon_screen_requires_both_strand_channels():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True), pos=100.0, neg=0.0
    )
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    accepted = screen_no_rna_exons(
        counts, region_arrays, payload_arrays, kappa_d_seed=10.0, rna_lower_confidence=0.95
    )
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
    summary = StrandSummary(p_r1_sense=0.95, n_observations=10_000)

    est = estimate_kappa_d(
        region_arrays,
        payload_arrays,
        counts,
        summary,
        rna_lower_confidence=0.95,
    )
    assert est.n_seed_regions == 4
    assert est.n_exon_self_training >= 1
    # Refit balance.n_regions should reflect seeds + accepted exons.
    assert est.balance.n_regions == est.n_seed_regions + est.n_exon_self_training


def test_kappa_soft_fallback_is_not_reported_as_uniform_fallback():
    region_arrays, payload_arrays = _seed_plus_exons_arrays()
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = StrandSummary(p_r1_sense=0.95, n_observations=10_000)

    est = estimate_kappa_d(
        region_arrays,
        payload_arrays,
        counts,
        summary,
        rna_lower_confidence=0.95,
    )

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
    summary = StrandSummary(p_r1_sense=0.95, n_observations=10_000)

    est = estimate_kappa_d(
        region_arrays,
        payload_arrays,
        counts,
        summary,
        rna_lower_confidence=0.95,
    )

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
    summary = StrandSummary(p_r1_sense=0.95, n_observations=10_000)

    est = estimate_kappa_d(
        region_arrays,
        payload_arrays,
        counts,
        summary,
        rna_lower_confidence=0.95,
    )
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
        rna_lower_confidence=0.95,
        kappa_d_fallback_used=True,
    )
    assert (deconv.flags[0] & FLAG_KAPPA_FALLBACK) != 0


def test_estimate_kappa_requires_matching_strand_summary():
    region_arrays, payload_arrays = _seed_plus_exons_arrays()
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.95)
    summary = StrandSummary(p_r1_sense=0.90, n_observations=10_000)

    with pytest.raises(ValueError, match="strand_summary.p_r1_sense"):
        estimate_kappa_d(
            region_arrays,
            payload_arrays,
            counts,
            summary,
            rna_lower_confidence=0.95,
        )


# ---------------------------------------------------------------------------
# Ineligible / sanity
# ---------------------------------------------------------------------------


def test_ineligible_regions_carry_zero_estimates():
    region_arrays, payload_arrays = _single_region(
        pack_signature(exon_pos=True, exon_neg=True), pos=50.0, neg=50.0
    )  # TS_AMBIG -> ineligible
    counts = build_strand_region_counts(region_arrays, payload_arrays, p_r1_sense=0.9)
    est = deconvolve_regions_by_strand(counts, kappa_d=5.0, rna_lower_confidence=0.95)
    assert (est.flags[0] & FLAG_INELIGIBLE) != 0
    assert est.rna_lower_count[0] == 0.0
    assert est.precision[0] == 0.0
