"""Top-level calibration orchestrator."""

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
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    build_fl_models,
)
from .scan_payload import CalibrationScanPayload
from .exposure import RegionExposure
from .strand_deconv import (
    build_strand_region_counts,
    deconvolve_regions_by_strand,
    estimate_kappa_d,
)
from .strand_summary import StrandSummary

if TYPE_CHECKING:
    from ..frag_length_model import FragmentLengthModels
    from ..index import TranscriptIndex


__all__ = ["calibrate"]


def calibrate(
    *,
    index: "TranscriptIndex",
    payload: CalibrationScanPayload,
    scan_trained: "FragmentLengthModels",
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
    pool_quality_good: int = POOL_QUALITY_GOOD_THRESHOLD,
    pool_quality_weak: int = POOL_QUALITY_WEAK_THRESHOLD,
    strand_summary: StrandSummary | None = None,
    rna_lower_confidence: float = 0.95,
) -> CalibrationResult:
    """Run the calibration stages that are live before locus EM."""
    if not (0.5 <= rna_lower_confidence < 1.0):
        raise ValueError(
            f"calibrate: rna_lower_confidence must be in [0.5, 1.0); got {rna_lower_confidence}."
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
        good_threshold=pool_quality_good,
        weak_threshold=pool_quality_weak,
    )

    from ._arrays import PayloadArrays, RegionArrays
    from .density_global import compute_global_densities

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    global_densities = compute_global_densities(
        region_arrays,
        payload_arrays,
        gdna_fl=fl_models.gdna,
    )
    if strand_summary is None:
        strand_summary = StrandSummary.uninformative()

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
        rna_lower_confidence=rna_lower_confidence,
    )
    region_gdna = deconvolve_regions_by_strand(
        strand_counts,
        kappa_d=kappa_d.kappa,
        rna_lower_confidence=rna_lower_confidence,
        kappa_d_n_seed_regions=kappa_d.n_seed_regions,
        kappa_d_n_exon_self_training=kappa_d.n_exon_self_training,
        kappa_d_fallback_used=kappa_d.fallback_used,
    )
    region_exposure = RegionExposure.uniform(R=int(region_gdna.n_total.size))

    return build_calibration_result(
        payload=payload,
        scan_trained=scan_trained,
        global_densities=global_densities,
        fl_models=fl_models,
        fl_prior_ess=fl_prior_ess,
        region_signature=region_arrays.signature,
        region_gdna=region_gdna,
        region_exposure=region_exposure,
        rna_lower_confidence=rna_lower_confidence,
    )
