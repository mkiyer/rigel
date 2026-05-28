"""Rigel gDNA calibration public surface."""

from ._diagnostics import Diagnostics
from ._exposure import l_eff_contained
from .background_model import BackgroundModel, fit_background_model
from .boundaries import BoundaryTable, build_boundary_table, validate_boundary_table
from .boundary_model import (
    BoundaryLocalPosterior,
    build_boundary_local_posterior,
    predict_contained_gdna_from_excess,
)
from .boundary_sweep import (
    BoundarySweepResult,
    compute_boundary_transfer_weight,
    run_boundary_sweep,
)
from .calibration_iteration import (
    CalibrationStepResult,
    RegionCalibration,
    calibration_e_step,
    calibration_m_step,
    run_calibration_iteration,
)
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._orchestrator import calibrate
from ._result import CalibrationResult, build_calibration_result
from .latent_states import (
    N_STATES,
    STATE_EXPRESSED,
    STATE_IS_EXPRESSED,
    STATE_UNEXPRESSED,
    STATE_NAMES,
    build_logbf_expression,
    build_logbf_strand,
    build_state_log_prior,
    build_state_log_tensor,
    build_state_tensor,
    normalize_state_log_tensor,
)
from .prior import (
    ComponentExposureTable,
    EMInputTable,
    PriorTable,
    assemble_em_inputs,
    assemble_priors,
)
from .strand_balance import StrandBalanceEstimate, estimate_strand_balance
from .strand_deconv import (
    CompartmentStrandCounts,
    RegionGdnaChannelEstimate,
    build_compartment_strand_counts,
    deconvolve_compartments_by_strand,
)
from .fl import (
    POOL_EB_PRIOR_ESS,
    POOL_SCORING_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    FLModels,
    Quality,
    build_fl_models,
)
from .strand_summary import StrandSummary

__all__ = [
    "calibrate",
    "CalibrationResult",
    "build_calibration_result",
    "BoundaryTable",
    "build_boundary_table",
    "validate_boundary_table",
    "BackgroundModel",
    "fit_background_model",
    "BoundaryLocalPosterior",
    "build_boundary_local_posterior",
    "predict_contained_gdna_from_excess",
    "BoundarySweepResult",
    "compute_boundary_transfer_weight",
    "run_boundary_sweep",
    "RegionCalibration",
    "CalibrationStepResult",
    "calibration_e_step",
    "calibration_m_step",
    "run_calibration_iteration",
    "STATE_UNEXPRESSED",
    "STATE_EXPRESSED",
    "N_STATES",
    "STATE_IS_EXPRESSED",
    "STATE_NAMES",
    "build_state_log_prior",
    "build_logbf_expression",
    "build_logbf_strand",
    "build_state_log_tensor",
    "normalize_state_log_tensor",
    "build_state_tensor",
    "CompartmentStrandCounts",
    "RegionGdnaChannelEstimate",
    "build_compartment_strand_counts",
    "deconvolve_compartments_by_strand",
    "PriorTable",
    "ComponentExposureTable",
    "EMInputTable",
    "assemble_em_inputs",
    "assemble_priors",
    "StrandBalanceEstimate",
    "estimate_strand_balance",
    "l_eff_contained",
    "StrandSummary",
    "FLModels",
    "Quality",
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "POOL_EB_PRIOR_ESS",
    "POOL_SCORING_PRIOR_ESS",
    "build_fl_models",
    "Diagnostics",
    "extract_global_counts",
    "extract_rna_counts",
    "extract_gdna_counts",
]
