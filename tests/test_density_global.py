"""Tests for ``rigel.calibration.density_global.compute_global_densities``.

Pins the three locked formulas (INTERGENIC, INTRON, EXON-INTRON), the
FL-PMF-weighted containment effective length
:math:`\\sum_\\ell h(\\ell)\\,\\max(0, L - \\ell + 1)`, and the
schema/validation contract.

See ``docs/calibration/density_eff_length_fix_2026_05.md``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._exposure import boundary_crossing_exposure
from rigel.calibration._kappa import KAPPA_DEFAULT
from rigel.calibration._orient import ORIENT_OPP, ORIENT_SAME, ORIENT_UNINF, StrandSummary
from rigel.calibration.density_global import (
    GlobalDensityTable,
    GlobalGdnaDensity,
    compute_global_densities,
    l_eff_contained,
)
from rigel.calibration.regions import RegionStrand, RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.frag_length_model import FragmentLengthModel


_MASK_EXON = 0b001
_MASK_INTRON = 0b010
_MASK_INTERGENIC = 0b100


def _delta_fl(length: int, *, max_size: int = 1024) -> FragmentLengthModel:
    """FL model with a sharp peak at ``length``.

    The 10_000 counts dominate Laplace smoothing but do not eliminate
    it; tests use :func:`_oracle_leff` to compute the exact expected
    L_eff under the same smoothed PMF that the production code sees.
    """
    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _oracle_leff(spans, gdna_fl: FragmentLengthModel) -> np.ndarray:
    """Ground-truth L_eff using the FL model's own eCDF — the
    contract test compares ``compute_global_densities`` against this.
    """
    return gdna_fl.compute_all_transcript_eff_lens(np.asarray(spans, dtype=np.int64), min_value=0.0)


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------


def _make_region_df(rows: list[dict]) -> pd.DataFrame:
    """Build a region DataFrame from row dicts (with sensible defaults)."""
    df = pd.DataFrame(rows)
    df["region_id"] = np.arange(len(df), dtype=np.int64)
    for col in ("start", "end"):
        df[col] = df[col].astype(np.int64)
    df["type"] = df["type"].astype(np.uint8)
    df.index = df["region_id"].to_numpy()
    return df


def _empty_payload(n_regions: int) -> dict:
    return {
        "global_counts": np.zeros(MASK_N_STATES, dtype=np.int64),
        "per_region_counts": np.zeros((n_regions, MASK_N_STATES), dtype=np.int64),
        "fl_hist": np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        "u_left": np.zeros(n_regions, dtype=np.int64),
        "u_right": np.zeros(n_regions, dtype=np.int64),
        "intron_counts_by_orient": np.zeros((n_regions, 3), dtype=np.int64),
        "u_left_by_orient": np.zeros((n_regions, 3), dtype=np.int64),
        "u_right_by_orient": np.zeros((n_regions, 3), dtype=np.int64),
        "n_observed": 0,
        "n_excluded_multimap": 0,
        "n_excluded_chimera": 0,
        "n_excluded_artifact": 0,
        "n_unobserved": 0,
        "n_unannotated_ref": 0,
    }


def _wrap_payload(d: dict) -> CalibrationScanPayload:
    n_obs = int(d["global_counts"].sum())
    d["n_observed"] = n_obs
    if not np.array_equal(
        d["intron_counts_by_orient"].sum(axis=1),
        d["per_region_counts"][:, _MASK_INTRON],
    ):
        d["intron_counts_by_orient"][:, :] = 0
        d["intron_counts_by_orient"][:, ORIENT_UNINF] = d["per_region_counts"][:, _MASK_INTRON]
    if not np.array_equal(d["u_left_by_orient"].sum(axis=1), d["u_left"]):
        d["u_left_by_orient"][:, :] = 0
        d["u_left_by_orient"][:, ORIENT_UNINF] = d["u_left"]
    if not np.array_equal(d["u_right_by_orient"].sum(axis=1), d["u_right"]):
        d["u_right_by_orient"][:, :] = 0
        d["u_right_by_orient"][:, ORIENT_UNINF] = d["u_right"]
    # Mirror n_obs into fl_hist[mask=1] bin 0 so the validator's
    # `sum(fl_hist) == n_observed` check passes — unused by M4.
    d["fl_hist"][_MASK_EXON, 0] = n_obs
    return CalibrationScanPayload.from_scan_dict(d, n_total=n_obs)


# ---------------------------------------------------------------------------
# Hand-counted scenarios
# ---------------------------------------------------------------------------


class TestHandCounted:
    def test_three_region_scenario(self):
        """One INTERGENIC, one INTRON, one EXON region with hand counts."""
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 1500,
                    "type": int(RegionType.EXON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": True,
                },
                {
                    "ref_name": "chr1",
                    "start": 1500,
                    "end": 2500,
                    "type": int(RegionType.INTRON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(3)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 50  # INTERGENIC
        d["per_region_counts"][2, _MASK_INTRON] = 20  # INTRON
        d["u_right"][1] = 7  # right edge of EXON
        d["global_counts"][_MASK_INTERGENIC] = 50
        d["global_counts"][_MASK_INTRON] = 20
        d["global_counts"][_MASK_EXON | _MASK_INTRON] = 7  # boundary fragments

        payload = _wrap_payload(d)
        fl_mean = 200
        gdna_fl = _delta_fl(fl_mean)
        out = compute_global_densities(df, payload, gdna_fl=gdna_fl)

        leff_ig = float(_oracle_leff([1000], gdna_fl)[0])
        assert out.intergenic.n_fragments == 50
        assert out.intergenic.eff_length_bp == pytest.approx(leff_ig)
        assert out.intergenic.rho == pytest.approx(50.0 / leff_ig)
        leff_in = leff_ig
        assert out.intron.n_fragments == 20
        assert out.intron.eff_length_bp == pytest.approx(leff_in)
        assert out.intron.rho == pytest.approx(20.0 / leff_in)
        # EXON-INTRON: per-side denominator is B_cross =
        # Σ h(ℓ) max(ℓ - 1, 0). Numerator = 7 * 1; denominator = 1 * B_cross.
        b_cross = boundary_crossing_exposure(gdna_fl)
        assert out.exon_intron.n_fragments == 7
        assert out.exon_intron.eff_length_bp == pytest.approx(b_cross)
        assert out.exon_intron.rho == pytest.approx(7.0 / b_cross)
        assert out.exon_intron.n_regions_used == 1

    def test_pure_mrna_library(self):
        """All counts in EXON_ONLY mask → all three gDNA densities zero."""
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 500,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 500,
                    "end": 1000,
                    "type": int(RegionType.EXON),
                    "boundary_flux_left": True,
                    "boundary_flux_right": True,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "type": int(RegionType.INTRON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(3)
        d["per_region_counts"][1, _MASK_EXON] = 100
        d["global_counts"][_MASK_EXON] = 100

        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl=_delta_fl(350))
        assert out.intergenic.rho == 0.0
        assert out.intron.rho == 0.0
        assert out.exon_intron.rho == 0.0
        # Empty numerators → kappa fallback.
        for d_ in (out.intergenic, out.intron, out.exon_intron):
            assert d_.kappa.value == KAPPA_DEFAULT
            assert d_.kappa.fallback_used

    def test_intron_only_library(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "type": int(RegionType.INTRON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(2)
        d["per_region_counts"][0, _MASK_INTRON] = 30
        d["per_region_counts"][1, _MASK_INTRON] = 70
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)
        out = compute_global_densities(df, payload, gdna_fl=gdna_fl)
        assert out.intergenic.rho == 0.0
        assert out.intron.n_fragments == 100
        leff = float(_oracle_leff([1000], gdna_fl)[0]) * 2
        assert out.intron.rho == pytest.approx(100.0 / leff)
        assert out.exon_intron.rho == 0.0


class TestExonIntronEligibility:
    def test_all_flags_false_zero_density(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 500,
                    "type": int(RegionType.EXON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["u_left"][0] = 99
        d["u_right"][0] = 99
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl=_delta_fl(200))
        assert out.exon_intron.rho == 0.0
        assert out.exon_intron.n_regions_used == 0
        assert out.exon_intron.eff_length_bp == 0.0

    def test_single_side_eligible(self):
        """bf_left=True, bf_right=False → u_R must NOT contribute."""
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 500,
                    "type": int(RegionType.EXON),
                    "boundary_flux_left": True,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["u_left"][0] = 5
        d["u_right"][0] = 999  # ineligible side; must be ignored
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(200)
        out = compute_global_densities(df, payload, gdna_fl=gdna_fl)
        b_cross = boundary_crossing_exposure(gdna_fl)
        assert out.exon_intron.n_fragments == 5
        assert out.exon_intron.eff_length_bp == pytest.approx(b_cross)
        assert out.exon_intron.rho == pytest.approx(5.0 / b_cross)


class TestStrandAwareDensity:
    def test_uninformative_strand_uses_uncorrected_count_exposure(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTRON] = 100
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 90
        d["intron_counts_by_orient"][0, ORIENT_OPP] = 10
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)

        uncorrected = compute_global_densities(df, payload, gdna_fl=gdna_fl)
        out = compute_global_densities(
            df,
            payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary.uninformative(),
        )

        assert out.intron.rho == pytest.approx(uncorrected.intron.rho)
        assert out.intron.rho_uncorrected == pytest.approx(uncorrected.intron.rho)
        assert not out.intron.strand_active
        assert out.intron.n_fragments_estimated == pytest.approx(100.0)

    def test_near_unstranded_sampling_noise_uses_uncorrected_count_exposure(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTRON] = 100
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 90
        d["intron_counts_by_orient"][0, ORIENT_OPP] = 10
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)

        uncorrected = compute_global_densities(df, payload, gdna_fl=gdna_fl)
        out = compute_global_densities(
            df,
            payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=0.5009, n_observations=324_000),
        )

        assert out.intron.rho == pytest.approx(uncorrected.intron.rho)
        assert out.intron.n_fragments_estimated == pytest.approx(100.0)
        assert not out.intron.strand_active

    def test_perfect_stranded_pure_rna_intron_zeroes_corrected_density(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTRON] = 100
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 100
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)

        out = compute_global_densities(
            df,
            payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
        )

        leff = float(_oracle_leff([1000], gdna_fl)[0])
        assert out.intron.strand_active
        assert out.intron.n_fragments == 100
        assert out.intron.rho_uncorrected == pytest.approx(100.0 / leff)
        assert out.intron.rho == 0.0
        assert out.intron.n_fragments_estimated == 0.0
        assert out.intron.strand_correction_fragments == pytest.approx(-100.0)

    def test_perfect_stranded_pure_gdna_intron_matches_uncorrected_density(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTRON] = 100
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 50
        d["intron_counts_by_orient"][0, ORIENT_OPP] = 50
        d["global_counts"][_MASK_INTRON] = 100
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)

        out = compute_global_densities(
            df,
            payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
        )

        leff = float(_oracle_leff([1000], gdna_fl)[0])
        assert out.intron.strand_active
        assert out.intron.rho == pytest.approx(100.0 / leff)
        assert out.intron.n_fragments_estimated == pytest.approx(100.0)

    def test_r1_antisense_mirror_is_invariant(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        gdna_fl = _delta_fl(300)

        sense_d = _empty_payload(1)
        sense_d["per_region_counts"][0, _MASK_INTRON] = 100
        sense_d["intron_counts_by_orient"][0, ORIENT_SAME] = 80
        sense_d["intron_counts_by_orient"][0, ORIENT_OPP] = 20
        sense_d["global_counts"][_MASK_INTRON] = 100
        sense_payload = _wrap_payload(sense_d)

        antisense_d = _empty_payload(1)
        antisense_d["per_region_counts"][0, _MASK_INTRON] = 100
        antisense_d["intron_counts_by_orient"][0, ORIENT_SAME] = 20
        antisense_d["intron_counts_by_orient"][0, ORIENT_OPP] = 80
        antisense_d["global_counts"][_MASK_INTRON] = 100
        antisense_payload = _wrap_payload(antisense_d)

        sense = compute_global_densities(
            df,
            sense_payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
        )
        antisense = compute_global_densities(
            df,
            antisense_payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=0.0, n_observations=1000),
        )

        assert sense.intron.rho == pytest.approx(antisense.intron.rho)
        assert sense.intron.n_fragments_estimated == pytest.approx(
            antisense.intron.n_fragments_estimated
        )

    def test_negative_moments_aggregate_before_final_clip(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(2)
        d["per_region_counts"][:, _MASK_INTRON] = [100, 40]
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 100
        d["intron_counts_by_orient"][1, ORIENT_OPP] = 40
        d["global_counts"][_MASK_INTRON] = 140
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(300)

        out = compute_global_densities(
            df,
            payload,
            gdna_fl=gdna_fl,
            strand_summary=StrandSummary(p_r1_sense=0.9, n_observations=1000),
        )

        leff = float(_oracle_leff([1000], gdna_fl)[0])
        # signed_strand_contrast = 0.8, so raw moments are -25 and 90.
        assert out.intron.rho == pytest.approx(65.0 / (2.0 * leff))
        assert out.intron.n_fragments_estimated == pytest.approx(65.0)

    def test_ambiguous_rows_are_diagnostics_not_corrected_numerator(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 2000,
                    "type": int(RegionType.INTRON),
                    "strand": int(RegionStrand.AMBIG),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(2)
        d["per_region_counts"][:, _MASK_INTRON] = [100, 200]
        d["intron_counts_by_orient"][0, ORIENT_SAME] = 100
        d["intron_counts_by_orient"][1, ORIENT_UNINF] = 200
        d["global_counts"][_MASK_INTRON] = 300
        payload = _wrap_payload(d)

        out = compute_global_densities(
            df,
            payload,
            gdna_fl=_delta_fl(300),
            strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
        )

        assert out.intron.strand_active
        assert out.intron.rho == 0.0
        assert out.intron.n_fragments_estimated == 0.0
        assert out.intron.n_rows_eligible == 2
        assert out.intron.n_strand_informative_regions == 1
        assert out.intron.strand_informative_exposure_fraction == pytest.approx(0.5)
        assert out.intron.n_uninf_observed == 200

    def test_boundary_flux_uses_orientation_payload(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 500,
                    "type": int(RegionType.EXON),
                    "strand": int(RegionStrand.POS),
                    "boundary_flux_left": True,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["u_left"][0] = 100
        d["u_left_by_orient"][0, ORIENT_SAME] = 100
        d["global_counts"][_MASK_EXON | _MASK_INTRON] = 100
        payload = _wrap_payload(d)

        out = compute_global_densities(
            df,
            payload,
            gdna_fl=_delta_fl(300),
            strand_summary=StrandSummary(p_r1_sense=1.0, n_observations=1000),
        )

        assert out.exon_intron.strand_active
        assert out.exon_intron.rho == 0.0
        assert out.exon_intron.n_fragments == 100
        assert out.exon_intron.n_fragments_estimated == 0.0


class TestEmptyPerType:
    def test_no_intron_regions(self):
        """Single-exon-genome: no INTRON rows at all → intron.rho == 0."""
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
                {
                    "ref_name": "chr1",
                    "start": 1000,
                    "end": 1500,
                    "type": int(RegionType.EXON),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(2)
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl=_delta_fl(200))
        assert out.intron.rho == 0.0
        assert out.intron.n_regions_used == 0
        assert out.intron.eff_length_bp == 0.0


class TestLeffGeometryLock:
    """The most important test in this file: pins the FL-PMF
    containment formula — ``l_eff_contained`` is a wrapper around
    :meth:`FragmentLengthModel.compute_all_transcript_eff_lens` with
    ``min_value=0.0`` (no salmon transcript-eff floor).
    """

    def test_l_eff_contained_matches_fl_model(self):
        gdna_fl = _delta_fl(200)
        spans = np.array([1000, 200, 50, 0])
        out = l_eff_contained(spans, gdna_fl)
        oracle = _oracle_leff(spans, gdna_fl)
        np.testing.assert_allclose(out, oracle)

    def test_l_eff_contained_floors_at_zero(self):
        # Span of 0 must give L_eff == 0 (not the salmon 1.0 floor).
        gdna_fl = _delta_fl(200)
        out = l_eff_contained(np.array([0]), gdna_fl)
        assert out[0] == 0.0

    def test_l_eff_contained_two_point_distribution(self):
        # Half-and-half FL between 100 and 300; for span L the
        # expected L_eff is approximately 0.5*(L-99) + 0.5*(L-299).
        # The exact value (with Laplace smoothing) is the same on both
        # sides via _oracle_leff; this test pins that the production
        # path agrees with the FL model's own computation.
        counts = np.zeros(1025, dtype=np.float64)
        counts[100] = 5_000.0
        counts[300] = 5_000.0
        gdna_fl = FragmentLengthModel.from_counts(counts, max_size=1024)
        spans = np.array([1000])
        out = l_eff_contained(spans, gdna_fl)
        oracle = _oracle_leff(spans, gdna_fl)
        np.testing.assert_allclose(out, oracle)

    def test_pinned_through_orchestrator(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 1000,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 10
        d["global_counts"][_MASK_INTERGENIC] = 10
        payload = _wrap_payload(d)
        gdna_fl = _delta_fl(200)
        out = compute_global_densities(df, payload, gdna_fl=gdna_fl)
        expected = float(_oracle_leff([1000], gdna_fl)[0])
        assert out.intergenic.eff_length_bp == pytest.approx(expected)
        assert out.intergenic.rho == pytest.approx(10.0 / expected)


# ---------------------------------------------------------------------------
# Schema / validation
# ---------------------------------------------------------------------------


class TestValidation:
    def test_region_count_mismatch(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 100,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(2)  # mismatched
        payload = _wrap_payload(d)
        with pytest.raises(ValueError, match="Rebuild the index"):
            compute_global_densities(df, payload, gdna_fl=_delta_fl(200))

    def test_gdna_fl_mean_must_be_positive(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 100,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        payload = _wrap_payload(d)
        # Empty FL counts → mean is 0 → error.
        empty_fl = FragmentLengthModel.from_counts(np.zeros(1025, dtype=np.float64), max_size=1024)
        with pytest.raises(ValueError, match="gdna_fl.mean"):
            compute_global_densities(df, payload, gdna_fl=empty_fl)

    def test_to_summary_dict_round_trip(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 100,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        d["per_region_counts"][0, _MASK_INTERGENIC] = 5
        d["global_counts"][_MASK_INTERGENIC] = 5
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl=_delta_fl(200))
        assert isinstance(out, GlobalDensityTable)
        sd = out.to_summary_dict()
        assert sd["gdna_fl_mean"] == pytest.approx(200.0, rel=2e-2)
        for key in ("INTERGENIC", "INTRON", "EXON-INTRON"):
            block = sd[key]
            for col in (
                "rho",
                "n_fragments",
                "eff_length_bp",
                "n_regions_used",
                "n_fragments_estimated",
                "n_rows_eligible",
                "strand_active",
                "rho_uncorrected",
                "strand_correction_fragments",
                "n_strand_informative_regions",
                "strand_informative_exposure_fraction",
                "n_uninf_observed",
                "kappa",
                "kappa_n_regions",
                "kappa_fallback",
                "kappa_reason",
            ):
                assert col in block, f"{key} missing {col}"

    def test_returns_globalgdnadensity_per_type(self):
        df = _make_region_df(
            [
                {
                    "ref_name": "chr1",
                    "start": 0,
                    "end": 100,
                    "type": int(RegionType.INTERGENIC),
                    "boundary_flux_left": False,
                    "boundary_flux_right": False,
                },
            ]
        )
        d = _empty_payload(1)
        payload = _wrap_payload(d)
        out = compute_global_densities(df, payload, gdna_fl=_delta_fl(200))
        for d_ in (out.intergenic, out.intron, out.exon_intron):
            assert isinstance(d_, GlobalGdnaDensity)
        assert out.intergenic.type == "INTERGENIC"
        assert out.intron.type == "INTRON"
        assert out.exon_intron.type == "EXON-INTRON"
