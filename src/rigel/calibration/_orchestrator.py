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
    strand_summary: StrandSummary | None = None,  # noqa: ARG001 - reserved for Phase 3
    rna_lower_confidence: float = 0.95,  # noqa: ARG001 - accepted for Phase 2/3 plumbing
) -> CalibrationResult:
    """Run the calibration stages that are live before locus EM."""
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

    return build_calibration_result(
        payload=payload,
        scan_trained=scan_trained,
        global_densities=global_densities,
        fl_models=fl_models,
        fl_prior_ess=fl_prior_ess,
        region_signature=region_arrays.signature,
    )
