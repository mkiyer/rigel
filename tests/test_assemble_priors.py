"""Tests for ``assemble_priors`` (PriorTable orchestrator)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from rigel.calibration._kappa import KappaEstimate
from rigel.calibration.density_global import GlobalDensityTable, GlobalGdnaDensity
from rigel.calibration.locus_prior import (
    C_BASE_DEFAULT,
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
            type="INTERGENIC", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        intron=GlobalGdnaDensity(
            type="INTRON", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        exon_intron=GlobalGdnaDensity(
            type="EXON-INTRON", rho=0.0, n_fragments=0, eff_length_bp=0.0,
            n_regions_used=0, kappa=_kappa_zero(),
        ),
        gdna_fl=_delta_fl(fl_mean),
    )


def _fake_index(
    region_rows: list[tuple[str, int, int, int, bool, bool]],
    transcripts: list[tuple[str, int, int]],   # (ref, start, end)
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
    return CalibrationScanPayload(
        global_counts=np.zeros(MASK_N_STATES, dtype=np.int64),
        per_region_counts=per_region,
        fl_hist=np.zeros((MASK_N_STATES, FL_HIST_N_BINS), dtype=np.int64),
        u_left=np.array(u_left if u_left is not None else [0] * n_regions, dtype=np.int64),
        u_right=np.array(u_right if u_right is not None else [0] * n_regions, dtype=np.int64),
        n_observed=0, n_excluded_multimap=0, n_excluded_chimera=0,
        n_excluded_artifact=0, n_unobserved=0, n_unannotated_ref=0,
    )


def _make_em(locus_t_indices: np.ndarray) -> SimpleNamespace:
    return SimpleNamespace(locus_t_indices=locus_t_indices)


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

def test_assemble_multilocus_prior_single_locus_matches_estimate():
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    est = LocusGdnaEstimate(
        locus=locus, n_obs=10,
        n_gdna_intergenic=3.0, n_gdna_intron=0.0, n_gdna_exon_intron=0.0,
        n_gdna=3.0, pi_gdna=0.3,
        rho_loco=(0.0, 0.0, 0.0), leff_loco=(0.0, 0.0, 0.0),
        n_eligible_boundaries=0, fallback_flags=0,
    )
    ml = _ml_single(0, [0], [10], locus)
    mlp = assemble_multilocus_prior(ml, (est,))
    assert mlp.n_obs == 10
    assert mlp.n_gdna == pytest.approx(3.0)
    assert mlp.n_rna == pytest.approx(7.0)
    assert mlp.pi_gdna == pytest.approx(0.3)
    assert mlp.per_locus == (est,)


def test_assemble_multilocus_prior_two_locus_aggregates():
    loc0 = Locus(ref="chr1", ref_id=0, start=0, end=500)
    loc1 = Locus(ref="chr1", ref_id=0, start=1000, end=1500)
    est0 = LocusGdnaEstimate(
        locus=loc0, n_obs=10, n_gdna_intergenic=2.0, n_gdna_intron=0.0,
        n_gdna_exon_intron=0.0, n_gdna=2.0, pi_gdna=0.2,
        rho_loco=(0.0,)*3, leff_loco=(0.0,)*3,
        n_eligible_boundaries=0, fallback_flags=0,
    )
    est1 = LocusGdnaEstimate(
        locus=loc1, n_obs=10, n_gdna_intergenic=4.0, n_gdna_intron=0.0,
        n_gdna_exon_intron=0.0, n_gdna=4.0, pi_gdna=0.4,
        rho_loco=(0.0,)*3, leff_loco=(0.0,)*3,
        n_eligible_boundaries=0, fallback_flags=0,
    )
    ml = MultiLocus(
        multi_locus_id=0,
        transcript_indices=np.array([0, 1], dtype=np.int32),
        unit_indices=np.array([10, 11], dtype=np.int32),
        gdna_span=1000, loci=(loc0, loc1),
    )
    mlp = assemble_multilocus_prior(ml, (est0, est1))
    assert mlp.n_obs == 20
    assert mlp.n_gdna == pytest.approx(6.0)
    assert mlp.pi_gdna == pytest.approx(0.3)


def test_assemble_priors_alpha_scaling():
    """alpha_gdna == c_base · pi_gdna; alpha_rna == c_base · (1 - pi_gdna)."""
    # Single intergenic locus, 10 obs, 5 gdna fragments ⇒ pi = 0.5.
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[5])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], locus)
    em = _make_em(np.zeros(10, dtype=np.int32))   # all anchor on transcript 0
    gdt = _gdt_zero(fl_mean=200)

    pt = assemble_priors(
        multi_loci=[ml], em_data=em, index=index,
        payload=payload, global_densities=gdt,
        gdna_fl=_delta_fl(200), c_base=10.0,
    )
    assert isinstance(pt, PriorTable)
    # n_gdna = ρ_loco · L_eff_contained = 5 (κ=0 ⇒ cancels);
    # n_obs = 10 ⇒ pi = 0.5.
    assert pt.multi_locus_priors[0].pi_gdna == pytest.approx(0.5)
    assert pt.alpha_gdna[0] == pytest.approx(10.0 * 0.5)
    assert pt.alpha_rna[0] == pytest.approx(10.0 * 0.5)
    assert pt.c_base_value == 10.0


def test_assemble_priors_prior_weight_rna_shape():
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800), ("chr1", 200, 900)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[0])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0, 1], [0], locus)
    em = _make_em(np.array([0], dtype=np.int32))
    pt = assemble_priors(
        multi_loci=[ml], em_data=em, index=index,
        payload=payload, global_densities=_gdt_zero(),
    )
    pwr = pt.prior_weight_rna
    assert isinstance(pwr, list)
    assert len(pwr) == 1
    # n_t = 2 ⇒ length = 3 (2 transcripts + 1 gDNA).
    assert pwr[0].shape == (3,)
    assert pwr[0].dtype == np.float32
    assert np.all(pwr[0] == 1.0)


def test_assemble_priors_default_c_base():
    index = _fake_index(
        region_rows=[("chr1", 0, 1000, int(RegionType.INTERGENIC), False, False)],
        transcripts=[("chr1", 100, 800)],
    )
    payload = _make_payload(n_regions=1, counts_intergenic=[10])
    locus = Locus(ref="chr1", ref_id=0, start=0, end=1000)
    ml = _ml_single(0, [0], [0] * 10, locus)
    em = _make_em(np.zeros(10, dtype=np.int32))
    pt = assemble_priors(
        multi_loci=[ml], em_data=em, index=index,
        payload=payload, global_densities=_gdt_zero(),
    )
    # Default c_base.
    assert pt.c_base_value == C_BASE_DEFAULT
