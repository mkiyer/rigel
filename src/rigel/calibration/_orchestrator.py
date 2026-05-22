"""rigel.calibration._orchestrator \u2014 top-level ``calibrate(...)``.

Under the fractional cutover this function builds the FL models, the
initial contained-region global gDNA density estimates, and a uniform
regional exposure scaffold. Downstream prior assembly is still a stub
that raises :class:`FractionalCutoverPending` (see ``locus_prior.py``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._regional_exposure import RegionalGdnaExposure
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
    strand_summary: StrandSummary | None = None,  # noqa: ARG001 \u2014 reserved
    regional_exposure_enabled: bool = True,  # noqa: ARG001 \u2014 reserved
    regional_exposure_reference_quantile: float = 0.95,  # noqa: ARG001
    resolver_splicing_anchor_tolerance: int = 0,
) -> CalibrationResult:
    """Run the fractional-era calibration pipeline.

    During the cutover this is intentionally limited to:

    1. Validating the index carries a ``region_df`` (still needed for
       region geometry in downstream stages).
    2. Building the FL models from the fractional payload + scanner-trained
       global model.
    3. Estimating contained intergenic/intronic global gDNA density and
       strand variance diagnostics.
    4. Returning a :class:`CalibrationResult` with an empty prior table.

    The CLI is expected to consume the FL models and the diagnostics
    block from the result, then fail fast on any attempt to assemble
    priors (``with_priors``) or run the EM stage.
    """
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

    # Build a degenerate uniform regional exposure off the sorted region
    # arrays so result-consumers that introspect ``regional_exposure``
    # see a well-formed (R,)-shape object instead of None. Cheap and
    # documented as the cutover-default behaviour.
    from ._arrays import PayloadArrays, RegionArrays  # local import: avoid load-time cycle
    from .density_global import compute_global_densities

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    global_densities = compute_global_densities(
        region_arrays,
        payload_arrays,
        gdna_fl=fl_models.gdna,
    )
    regional_exposure = RegionalGdnaExposure.uniform(region_arrays)

    return build_calibration_result(
        payload=payload,
        scan_trained=scan_trained,
        global_densities=global_densities,
        fl_models=fl_models,
        fl_prior_ess=fl_prior_ess,
        regional_exposure=regional_exposure,
        resolver_splicing_anchor_tolerance=int(resolver_splicing_anchor_tolerance),
        region_signature=region_arrays.signature,
    )
