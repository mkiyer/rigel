"""Rigel gDNA calibration (v6).

Public entry point is :func:`calibrate`, which assembles a
:class:`CalibrationResult` carrying:

* :class:`FLModels` — empirical-Bayes shrunk RNA / gDNA / global
  fragment-length distributions, ready for scoring.
* :class:`GlobalDensityTable` — per-region eligible counts and the
  global gDNA density used by per-locus prior assembly.
* :class:`PriorTable` — back-filled by
  :meth:`CalibrationResult.with_priors` after multi-locus
  construction in :func:`rigel.pipeline.quant_from_buffer`.
"""

from ._diagnostics import Diagnostics
from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._kappa import (
    KAPPA_DEFAULT,
    KAPPA_MAX,
    KAPPA_MIN,
    KappaEstimate,
    estimate_kappa,
)
from ._orient import StrandSummary
from ._orchestrator import calibrate
from ._result import (
    CalibrationResult,
    build_calibration_result,
    build_multi_locus_prior_df,
    build_per_locus_gdna_df,
)
from .density_global import (
    GlobalDensityTable,
    GlobalGdnaDensity,
    compute_global_densities,
    l_eff_contained,
)
from .density_loco import shrink_to_loco
from .fl import (
    POOL_EB_PRIOR_ESS,
    POOL_QUALITY_GOOD_THRESHOLD,
    POOL_QUALITY_WEAK_THRESHOLD,
    FLModels,
    Quality,
    build_fl_models,
)
from .locus_prior import (
    FLAG_EXON_INTRON_NO_ELIGIBLE,
    FLAG_INTERGENIC_ZERO_LEFF,
    FLAG_INTRON_ZERO_LEFF,
    FLAG_PI_CLIPPED,
    ExpectedGdnaPriorParts,
    LocusGdnaEstimate,
    MultiLocusPrior,
    PriorTable,
    assemble_multilocus_prior,
    assemble_priors,
    build_prior_weight_rna,
    enable_gdna_for_multilocus,
    estimate_locus_gdna,
    expected_gdna_count_global,
)

__all__ = [
    # Top-level orchestrator
    "calibrate",
    # CalibrationResult + helpers
    "CalibrationResult",
    "build_calibration_result",
    "build_multi_locus_prior_df",
    "build_per_locus_gdna_df",
    # Global gDNA densities + κ
    "GlobalDensityTable",
    "GlobalGdnaDensity",
    "compute_global_densities",
    "l_eff_contained",
    "KappaEstimate",
    "estimate_kappa",
    "KAPPA_DEFAULT",
    "KAPPA_MIN",
    "KAPPA_MAX",
    "StrandSummary",
    # EM prior_weight_rna helper
    "build_prior_weight_rna",
    # Locoregional shrinkage + per-locus / MultiLocus priors
    "shrink_to_loco",
    "LocusGdnaEstimate",
    "MultiLocusPrior",
    "PriorTable",
    "estimate_locus_gdna",
    "expected_gdna_count_global",
    "enable_gdna_for_multilocus",
    "ExpectedGdnaPriorParts",
    "assemble_multilocus_prior",
    "assemble_priors",
    "FLAG_INTERGENIC_ZERO_LEFF",
    "FLAG_INTRON_ZERO_LEFF",
    "FLAG_EXON_INTRON_NO_ELIGIBLE",
    "FLAG_PI_CLIPPED",
    # Fragment-length models + diagnostics
    "FLModels",
    "Quality",
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "POOL_EB_PRIOR_ESS",
    "build_fl_models",
    "Diagnostics",
    "extract_global_counts",
    "extract_rna_counts",
    "extract_gdna_counts",
]
