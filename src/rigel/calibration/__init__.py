"""Rigel gDNA calibration.

The legacy SRD v1 surface (``calibrate_gdna``, ``CalibrationResult``,
``GdnaFlQuality``) is kept here while the v6 surface is built up
incrementally per ``docs/calibration/calibration_v6_plan.md``.  The v6
exports below are additive; they will replace the v1 names in M8c.
"""

from ._kappa import (
    KAPPA_DEFAULT,
    KAPPA_MAX,
    KAPPA_MIN,
    KappaEstimate,
    estimate_kappa,
)
from ._result import CalibrationResult, GdnaFlQuality
from ._simple import calibrate_gdna
from .density_global import (
    GlobalDensityTable,
    GlobalGdnaDensity,
    compute_global_densities,
    l_eff_overlap,
)
from .density_loco import shrink_to_loco
from .locus_prior import (
    C_BASE_DEFAULT,
    FLAG_EXON_INTRON_NO_ELIGIBLE,
    FLAG_INTERGENIC_ZERO_LEFF,
    FLAG_INTRON_ZERO_LEFF,
    FLAG_PI_CLIPPED,
    LocusGdnaEstimate,
    MultiLocusPrior,
    PriorTable,
    assemble_multilocus_prior,
    assemble_priors,
    build_prior_weight_rna,
    estimate_locus_gdna,
)


__all__ = [
    # SRD v1 (legacy; deleted in M8c)
    "CalibrationResult",
    "GdnaFlQuality",
    "calibrate_gdna",
    # M4 — global gDNA densities + κ
    "GlobalDensityTable",
    "GlobalGdnaDensity",
    "compute_global_densities",
    "l_eff_overlap",
    "KappaEstimate",
    "estimate_kappa",
    "KAPPA_DEFAULT",
    "KAPPA_MIN",
    "KAPPA_MAX",
    # M5 — EM prior_weight_rna helper
    "build_prior_weight_rna",
    # M6 — locoregional shrinkage + per-Locus / MultiLocus priors
    "shrink_to_loco",
    "LocusGdnaEstimate",
    "MultiLocusPrior",
    "PriorTable",
    "estimate_locus_gdna",
    "assemble_multilocus_prior",
    "assemble_priors",
    "C_BASE_DEFAULT",
    "FLAG_INTERGENIC_ZERO_LEFF",
    "FLAG_INTRON_ZERO_LEFF",
    "FLAG_EXON_INTRON_NO_ELIGIBLE",
    "FLAG_PI_CLIPPED",
]
