"""Top-level calibration orchestrator for the v6 clean cutover."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._result import CalibrationResult, build_calibration_result
from .fl import (
    POOL_EB_PRIOR_ESS,
    POOL_SCORING_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    build_fl_models,
)
from .scan_payload import CalibrationScanPayload
from .strand_deconv import (
    build_compartment_strand_counts,
    build_strand_region_counts,
    deconvolve_compartments_by_strand,
    estimate_kappa_d,
)
from .strand_summary import STRAND_CONTRAST_NUMERICAL_FLOOR, StrandSummary

if TYPE_CHECKING:
    from ..frag_length_model import FragmentLengthModels
    from ..index import TranscriptIndex


__all__ = ["calibrate"]

_INTERNAL_GDNA_DENSITY_CI: float = 0.95


def _strand_summary_identifiable(
    strand_summary: StrandSummary,
    *,
    confidence: float = 0.99,
) -> bool:
    effective_min = max(
        STRAND_CONTRAST_NUMERICAL_FLOOR,
        strand_summary.signed_strand_contrast_margin(confidence=confidence),
    )
    return abs(strand_summary.signed_strand_contrast) >= effective_min


def calibrate(
    *,
    index: "TranscriptIndex",
    payload: CalibrationScanPayload,
    scan_trained: "FragmentLengthModels",
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
    fl_scoring_prior_ess: float = POOL_SCORING_PRIOR_ESS,
    pool_quality_good: int = POOL_QUALITY_GOOD_THRESHOLD,
    pool_quality_weak: int = POOL_QUALITY_WEAK_THRESHOLD,
    strand_summary: StrandSummary,
    density_min_eff_length: float = 1.0,
    background_trim_fraction: float = 0.01,
    max_calibration_passes: int = 5,
) -> CalibrationResult:
    """Run the v6 calibration stages that feed locus EM."""
    if density_min_eff_length < 0.0:
        raise ValueError(
            f"calibrate: density_min_eff_length must be >= 0; got {density_min_eff_length}."
        )
    if not (0.0 <= float(background_trim_fraction) < 1.0):
        raise ValueError(
            "calibrate: background_trim_fraction must be in [0, 1); "
            f"got {background_trim_fraction}."
        )
    if int(max_calibration_passes) < 1:
        raise ValueError(
            "calibrate: max_calibration_passes must be >= 1; "
            f"got {max_calibration_passes}."
        )
    if index.region_df is None:
        raise RuntimeError(
            "Index has no region table. Rebuild the index "
            "(rigel index --fasta ... --gtf ...). Older indexes "
            "are not supported."
        )

    fl_models = build_fl_models(
        global_counts=extract_global_counts(scan_trained),
        rna_counts=extract_rna_counts(scan_trained),
        gdna_counts=extract_gdna_counts(payload),
        max_size=scan_trained.max_size,
        prior_ess=fl_prior_ess,
        scoring_prior_ess=fl_scoring_prior_ess,
        good_threshold=pool_quality_good,
        weak_threshold=pool_quality_weak,
    )

    from ._arrays import PayloadArrays, RegionArrays
    from .background_model import fit_background_model
    from .boundaries import build_boundary_table
    from .calibration_iteration import run_calibration_iteration
    from .density_model import fit_density_evidence
    from .density_observation import build_density_observation
    from .region_count_ledger import build_region_count_ledger

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    observation = build_density_observation(region_arrays, ledger, fl_models.gdna)
    density_evidence = fit_density_evidence(
        observation,
        confidence=_INTERNAL_GDNA_DENSITY_CI,
        min_eff_length=float(density_min_eff_length),
    )
    force_zero_gdna_mass = density_evidence.rho_ref_source == "ZERO"

    strand_usable = _strand_summary_identifiable(strand_summary)

    strand_counts = build_strand_region_counts(
        region_arrays,
        payload_arrays,
        p_r1_sense=strand_summary.p_r1_sense,
    )
    kappa_d = estimate_kappa_d(
        region_arrays,
        payload_arrays,
        strand_counts,
        strand_summary,
    )
    compartment_counts = build_compartment_strand_counts(
        region_arrays,
        payload_arrays,
        p_r1_sense=strand_summary.p_r1_sense,
    )
    strand_channels = deconvolve_compartments_by_strand(
        compartment_counts,
        kappa_d=kappa_d.kappa,
        strand_summary=strand_summary,
    )
    calibration_strand_channels = strand_channels if strand_usable else None

    background = fit_background_model(
        observation,
        calibration_strand_channels,
        top_t_fraction=float(background_trim_fraction),
        min_eff_length=float(density_min_eff_length),
    )
    boundaries = build_boundary_table(
        region_arrays,
        ledger,
        observation.boundary_left_leff,
    )
    region_calibration = run_calibration_iteration(
        region_arrays,
        observation,
        boundaries,
        background,
        strand_channels=calibration_strand_channels,
        max_calibration_passes=int(max_calibration_passes),
        confidence=_INTERNAL_GDNA_DENSITY_CI,
        # PR 03: feed integer unspliced counts to enable RegionUnsplicedMass /
        # BackgroundDensity construction inside the calibration loop.
        unspliced_counts=ledger.unspliced_support,
        force_zero_gdna_mass=force_zero_gdna_mass,
    )

    return build_calibration_result(
        payload=payload,
        scan_trained=scan_trained,
        fl_models=fl_models,
        fl_prior_ess=fl_prior_ess,
        fl_scoring_prior_ess=fl_scoring_prior_ess,
        region_signature=region_arrays.signature,
        region_calibration=region_calibration,
        strand_channels=strand_channels,
    )
