"""Tests for ``assemble_priors`` (PriorTable orchestrator)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._kappa import KappaEstimate
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.locus_prior import (
    LocusGdnaEstimate,
    PriorTable,
    assemble_multilocus_prior,
    assemble_priors,
)
from rigel.calibration.regions import RegionType
from rigel.calibration.scan_payload import (
    FL_HIST_N_BINS,
    MASK_N_STATES,
    CalibrationScanPayload,
)
from rigel.locus import Locus, MultiLocus


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------


def _kappa_zero() -> KappaEstimate:
    return KappaEstimate(value=0.0, n_regions=1, fallback_used=False, fallback_reason="")


def _delta_fl(length: int, *, max_size: int = 1024):
    from rigel.frag_length_model import FragmentLengthModel

    counts = np.zeros(max_size + 1, dtype=np.float64)
    counts[length] = 10_000.0
    return FragmentLengthModel.from_counts(counts, max_size=max_size)


def _gdt_zero(fl_mean: int = 200) -> GlobalDensityTable:
    return GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa_zero(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa_zero(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa_zero(),
        ),
        gdna_fl=_delta_fl(fl_mean),
    )


def _fake_index(
    region_rows: list[tuple[str, int, int, int, bool, bool]],
    transcripts: list[tuple[str, int, int]],  # (ref, start, end)
):
    """Build a SimpleNamespace mimicking the bits of TranscriptIndex
    that ``assemble_priors`` reads."""
    refs = sorted({r[0] for r in region_rows} | {t[0] for t in transcripts})
    ref_name_to_id = {name: i for i, name in enumerate(refs)}
    n = len(region_rows)
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n, dtype=np.int64),
            "ref_name": [r[0] for r in region_rows],
            "start": np.array([r[1] for r in region_rows], dtype=np.int64),
            "end": np.array([r[2] for r in region_rows], dtype=np.int64),
            "type": np.array([r[3] for r in region_rows], dtype=np.uint8),
            "strand": np.zeros(n, dtype=np.uint8),
            "tx_pos_bp": np.zeros(n, dtype=np.int64),
            "tx_neg_bp": np.zeros(n, dtype=np.int64),
            "exon_pos_bp": np.zeros(n, dtype=np.int64),
            "exon_neg_bp": np.zeros(n, dtype=np.int64),
            "boundary_flux_left": np.array([r[4] for r in region_rows], dtype=bool),
            "boundary_flux_right": np.array([r[5] for r in region_rows], dtype=bool),
        }
    )
    t_df = pd.DataFrame(
        {
            "ref": pd.Categorical([t[0] for t in transcripts], categories=refs),
            "start": np.array([t[1] for t in transcripts], dtype=np.int64),
            "end": np.array([t[2] for t in transcripts], dtype=np.int64),
        }
    )
    return SimpleNamespace(
        region_df=region_df,
        t_df=t_df,
        ref_name_to_id=ref_name_to_id,
        # Mirror TranscriptIndex.ref_lengths (insertion-ordered dict aligned
        # with ref_name_to_id). Use a length larger than any region/transcript
        # extent so the Phase-4 intergenic flank clamp is a no-op for tests.
        ref_lengths={name: 1_000_000 for name in refs},
    )


def _make_payload(
    n_regions: int,
    counts_intergenic: list[int] | None = None,
    counts_intron: list[int] | None = None,
    u_left: list[int] | None = None,
    u_right: list[int] | None = None,
) -> CalibrationScanPayload:
    per_region = np.zeros((n_regions, MASK_N_STATES), dtype=np.int64)
    if counts_intergenic is not None:
        per_region[:, 0b100] = np.array(counts_intergenic, dtype=np.int64)
    if counts_intron is not None:
        per_region[:, 0b010] = np.array(counts_intron, dtype=np.int64)
    intron_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    intron_by_orient[:, 2] = per_region[:, 0b010]
    u_left_arr = np.array(u_left if u_left is not None else [0] * n_regions, dtype=np.int64)
    u_right_arr = np.array(u_right if u_right is not None else [0] * n_regions, dtype=np.int64)
    u_left_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    u_right_by_orient = np.zeros((n_regions, 3), dtype=np.int64)
    u_left_by_orient[:, 2] = u_left_arr
    u_right_by_orient[:, 2] = u_right_arr
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=per_region,
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=u_left_arr,
        u_right=u_right_arr,
        intron_counts_by_orient=intron_by_orient,
        u_left_by_orient=u_left_by_orient,
        u_right_by_orient=u_right_by_orient,
        n_observed=0,
        n_excluded_multimap=0,
        n_excluded_chimera=0,
        n_excluded_artifact=0,
        n_unobserved=0,
        n_unannotated_ref=0,
    )


def _make_em(
    locus_t_indices: np.ndarray,
    *,
    n_units: int | None = None,
) -> SimpleNamespace:
    if n_units is None:
        n_units = int(locus_t_indices.size + 32)
    return SimpleNamespace(
        locus_t_indices=locus_t_indices,
        # Phase\u00a02: enable_gdna_for_multilocus consults these directly.
        # Default: every unit is unspliced with finite gDNA log-lik
        # \u21d2 enable_gdna == 1 unless the test overrides.
        is_spliced=np.zeros(n_units, dtype=np.uint8),
        gdna_log_liks=np.full(n_units, -1.0, dtype=np.float64),
    )


def _ml_single(
    multi_locus_id: int,
    transcript_indices: list[int],
    unit_indices: list[int],
    locus: Locus,
) -> MultiLocus:
    return MultiLocus(
        multi_locus_id=multi_locus_id,
        transcript_indices=np.array(transcript_indices, dtype=np.int32),
        unit_indices=np.array(unit_indices, dtype=np.int32),
        gdna_span=int(locus.span),
        loci=(locus,),
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def _make_estimate(
    locus: Locus,
    n_obs: int,
    *,
    n_gdna_intergenic: float = 0.0,
    n_gdna_intron: float = 0.0,
    n_gdna_boundary_observed: float = 0.0,
    n_gdna_exon_only: float = 0.0,
) -> LocusGdnaEstimate:
    """Test helper: construct a ``LocusGdnaEstimate`` with consistent totals."""
    n_gdna = n_gdna_intergenic + n_gdna_intron + n_gdna_boundary_observed + n_gdna_exon_only
    pi = (n_gdna / n_obs) if n_obs > 0 else 0.0
    return LocusGdnaEstimate(
        locus=locus,
        n_obs=n_obs,
        n_gdna_intergenic=n_gdna_intergenic,
        n_gdna_intron=n_gdna_intron,
        n_gdna_boundary_observed=n_gdna_boundary_observed,
        n_gdna_exon_only=n_gdna_exon_only,
        n_gdna=n_gdna,
        pi_gdna=min(1.0, max(0.0, pi)),
        rho_loco=(0.0, 0.0, 0.0),
        leff_loco=(0.0, 0.0, 0.0),
        n_eligible_boundaries=0,
        n_boundary_events=n_gdna_boundary_observed,
        nrna_active=False,
        fallback_flags=0,
    )


def test_assemble_multilocus_prior_single_locus_matches_estimate():
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = _make_estimate(locus, n_obs=10, n_gdna_intergenic=3.0)
    ml = _ml_single(0, [0], [10], locus)
    mlp = assemble_multilocus_prior(
        ml,
        (est,),
        gdna_prior_count=3.0,
    )
    assert mlp.n_obs == 10
    assert mlp.n_gdna == pytest.approx(3.0)
    assert mlp.n_rna == pytest.approx(7.0)
    assert mlp.pi_gdna == pytest.approx(0.3)
    assert mlp.gdna_prior_count == pytest.approx(3.0)
    assert mlp.per_locus == (est,)


def test_assemble_multilocus_prior_two_locus_aggregates():
    loc0 = Locus(ref="chr1", ref_id=0, start=0, end=500)
    loc1 = Locus(ref="chr1", ref_id=0, start=1000, end=1500)
    est0 = _make_estimate(loc0, n_obs=10, n_gdna_intergenic=2.0)
    est1 = _make_estimate(loc1, n_obs=10, n_gdna_intergenic=4.0)
    ml = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([0, 1], dtype=np.int32),
        unit_indices=np.array([10, 11], dtype=np.int32),
        gdna_span=1000,
        loci=(loc0, loc1),
    )
    mlp = assemble_multilocus_prior(
        ml,
        (est0, est1),
        gdna_prior_count=6.0,
    )
    assert mlp.n_obs == 20
    assert mlp.n_gdna == pytest.approx(6.0)
    assert mlp.pi_gdna == pytest.approx(0.3)


def test_assemble_priors_global_only_zero_density():
    """Phase 2: global-only prior is zero when global rho == 0 in every branch.

    The diagnostic ``pi_gdna`` is still populated from the legacy
    locoregional path for the per-locus dataframe, but it is no longer
    consumed by the EM prior.
    """
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[5])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], locus)
    em = _make_em(np.zeros(10, dtype=np.int32))
    gdt = _gdt_zero(fl_mean=200)

    pt = assemble_priors(
        multi_loci=[ml],
        em_data=em,
        index=index,
        payload=payload,
        global_densities=gdt,
        gdna_fl=_delta_fl(200),
    )
    assert isinstance(pt, PriorTable)
    # Diagnostic pi_gdna still computed from the legacy locoregional path.
    assert pt.multi_locus_priors[0].pi_gdna == pytest.approx(0.5)
    # Canonical prior: rho == 0 everywhere => eta_g == 0.
    assert pt.gdna_prior_count[0] == pytest.approx(0.0)
    assert pt.gdna_eff_len.shape == (1,)
    assert pt.gdna_eff_len[0] >= 1.0
    # Eligibility decoupled from prior strength: every unit is unspliced
    # with finite gdna_log_lik in the fixture, so enable_gdna == 1
    # even though eta_g == 0.
    assert pt.enable_gdna[0] == 1


def test_assemble_priors_global_only_positive_density():
    """gDNA prior count is positive when global gDNA density is positive."""
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[10])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], [0] * 10, locus)
    em = _make_em(np.zeros(10, dtype=np.int32))
    gdt = GlobalDensityTable(
        intergenic=GlobalGdnaDensity(
            type="INTERGENIC",
            rho=1e-3,
            n_fragments=10,
            eff_length_bp=10_000.0,
            n_regions_used=1,
            kappa=_kappa_zero(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa_zero(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON",
            rho=0.0,
            n_fragments=0,
            eff_length_bp=0.0,
            n_regions_used=0,
            kappa=_kappa_zero(),
        ),
        gdna_fl=_delta_fl(200),
    )
    pt = assemble_priors(
        multi_loci=[ml],
        em_data=em,
        index=index,
        payload=payload,
        global_densities=gdt,
    )
    # eta_g must be strictly positive: rho_ig > 0 and L_ig > 0.
    assert pt.gdna_prior_count[0] > 0.0
    assert pt.gdna_eff_len[0] >= 1.0


def test_assemble_priors_uniform_regional_exposure_matches_baseline():
    """With uniform regional exposure (default), weighted ≡ unweighted L_g."""
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[5])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], list(range(10)), locus)
    em = _make_em(np.zeros(10, dtype=np.int32))
    gdt = _gdt_zero(fl_mean=200)

    pt = assemble_priors(
        multi_loci=[ml],
        em_data=em,
        index=index,
        payload=payload,
        global_densities=gdt,
        gdna_fl=_delta_fl(200),
    )
    assert pt.gdna_eff_len_unweighted.shape == (1,)
    assert pt.gdna_prior_count_em.shape == (1,)
    assert pt.gdna_em_exposure_weight.shape == (1,)
    assert pt.gdna_eff_len[0] == pytest.approx(pt.gdna_eff_len_unweighted[0])
    assert pt.gdna_prior_count_em[0] == pytest.approx(pt.gdna_prior_count[0])
    assert pt.gdna_em_exposure_weight[0] == pytest.approx(1.0)


def test_assemble_priors_regional_exposure_attenuates_weighted_l_g():
    """An explicit regional exposure with sub-unit weights must shrink
    ``gdna_eff_len`` strictly below the unweighted baseline."""
    from rigel.calibration._arrays import RegionArrays
    from rigel.calibration._regional_exposure import RegionalGdnaExposure

    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[5])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], list(range(10)), locus)
    em = _make_em(np.zeros(10, dtype=np.int32))
    gdt = _gdt_zero(fl_mean=200)
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    weights = np.array([0.1], dtype=np.float64)
    regional = RegionalGdnaExposure(
        rho_hat=np.array([1e-4]),
        log_weight=np.log(weights),
        weight=weights,
        mode="regional",
        rho_ref=1e-3,
        n_at_floor=0,
        per_class={},
        ref_offsets=region_arrays.ref_offsets.copy(),
        ref_id=region_arrays.ref_id.copy(),
        start=region_arrays.start.copy(),
        end=region_arrays.end.copy(),
    )
    pt = assemble_priors(
        multi_loci=[ml],
        em_data=em,
        index=index,
        payload=payload,
        global_densities=gdt,
        gdna_fl=_delta_fl(200),
        regional_exposure=regional,
    )
    # Weight 0.1 attenuates weighted L_g below unweighted.
    assert pt.gdna_eff_len[0] < pt.gdna_eff_len_unweighted[0]
    assert pt.gdna_em_exposure_weight[0] == pytest.approx(0.1)
    # gdna_prior_count (canonical eta_g) unchanged by regional weighting.
    assert pt.gdna_prior_count[0] == pytest.approx(0.0)  # gdt has rho=0
    assert pt.gdna_prior_count_em[0] == pytest.approx(0.0)

