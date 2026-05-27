"""Projection of adaptive calibration evidence onto MultiLocus EM inputs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..config import EMConfig
from ._arrays import RegionArrays
from ._exposure import bp_weighted_mean_exposure_over_blocks, gdna_eff_len_for_loci
from .adaptive_prior import (
    MAX_ESS,
    PRIOR_BIAS_APPLIED,
    PRIOR_ESS_CAPPED,
    PRIOR_NO_UNSPLICED_MASS,
    PRIOR_STRUCTURAL_GATED,
    compute_adaptive_prior,
)

if TYPE_CHECKING:  # pragma: no cover - annotations only.
    from ..index import TranscriptIndex
    from ..locus import MultiLocus
    from ..scored_fragments import ScoredFragments
    from ._result import CalibrationResult


__all__ = [
    "PriorTable",
    "assemble_priors",
    "enable_gdna_for_multilocus",
]

_INT64_MIN: int = -(2**63)


@dataclass(frozen=True, slots=True)
class PriorTable:
    """Per-MultiLocus adaptive prior inputs and diagnostics."""

    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray
    gdna_eff_len: np.ndarray
    gdna_eff_len_unweighted: np.ndarray
    gdna_em_exposure_weight: np.ndarray
    enable_gdna: np.ndarray
    prior_locus_weight: np.ndarray
    prior_shrink_weight: np.ndarray
    prior_n_local_gdna: np.ndarray
    prior_n_local_rna: np.ndarray
    prior_n_other_gdna: np.ndarray
    prior_n_other_rna: np.ndarray
    prior_ess_final: np.ndarray
    prior_rna_share_v5: np.ndarray
    prior_rna_share_final: np.ndarray
    prior_unspliced_total: np.ndarray
    prior_flags: np.ndarray
    prior_region_weight: np.ndarray
    n_regions_touched: np.ndarray
    n_units_used_for_diagnostics: np.ndarray
    multi_locus_region_mass: np.ndarray
    partial_coverage_region_mass: np.ndarray
    global_n_gdna: float = 0.0
    global_n_rna: float = 0.0
    unallocated_unspliced_count: float = 0.0
    unallocated_weighted_unspliced_count: float = 0.0
    prior_policy_name: str = "entropy_dirichlet_v5_v6"
    prior_max_ess: float = MAX_ESS
    rna_call_bias: float = 0.5

    def to_summary_dict(self) -> dict[str, object]:
        """Return compact JSON-safe adaptive-prior diagnostics."""
        share_shift = np.asarray(self.prior_rna_share_final, dtype=np.float64) - np.asarray(
            self.prior_rna_share_v5,
            dtype=np.float64,
        )
        share_shift = share_shift[np.isfinite(share_shift)]
        if share_shift.size == 0:
            share_shift = np.array([0.0], dtype=np.float64)
        return {
            "name": str(self.prior_policy_name),
            "max_ess": float(self.prior_max_ess),
            "rna_call_bias": float(self.rna_call_bias),
            "global_counts": {
                "gdna": float(self.global_n_gdna),
                "rna": float(self.global_n_rna),
            },
            "region_weight": _summary_stats(self.prior_region_weight),
            "locus_weight": _summary_stats(self.prior_locus_weight),
            "shrink_weight": _summary_stats(self.prior_shrink_weight),
            "prior_ess_final": _summary_stats(self.prior_ess_final),
            "rna_share_shift_p50": float(np.quantile(share_shift, 0.50)),
            "rna_share_shift_p90": float(np.quantile(share_shift, 0.90)),
            "n_loci_total": int(self.alpha_gdna_add.size),
            "n_loci_with_prior_mass": int(np.count_nonzero(self.prior_ess_final > 0.0)),
            "n_loci_structural_gated": int(
                np.count_nonzero((self.prior_flags & PRIOR_STRUCTURAL_GATED) != 0)
            ),
            "n_loci_no_unspliced_mass": int(
                np.count_nonzero((self.prior_flags & PRIOR_NO_UNSPLICED_MASS) != 0)
            ),
            "n_loci_ess_capped": int(
                np.count_nonzero((self.prior_flags & PRIOR_ESS_CAPPED) != 0)
            ),
            "n_loci_bias_applied": int(
                np.count_nonzero((self.prior_flags & PRIOR_BIAS_APPLIED) != 0)
            ),
            "flag_histogram": _prior_flag_histogram(self.prior_flags),
            "n_regions_touched": _summary_stats(self.n_regions_touched),
            "gdna_eff_len": _summary_stats(self.gdna_eff_len),
            "gdna_em_exposure_weight": _summary_stats(self.gdna_em_exposure_weight),
            "multi_locus_region_mass": _summary_stats(self.multi_locus_region_mass),
            "partial_coverage_region_mass": _summary_stats(self.partial_coverage_region_mass),
            "unallocated_unspliced_count": float(self.unallocated_unspliced_count),
            "unallocated_weighted_unspliced_count": float(
                self.unallocated_weighted_unspliced_count
            ),
        }


def assemble_priors(
    *,
    multi_loci: list["MultiLocus"],
    em_data: "ScoredFragments",
    index: "TranscriptIndex",
    calibration: "CalibrationResult",
    em_config: EMConfig | None = None,
) -> PriorTable:
    """Assemble adaptive grouped priors, denominators, and native EM eligibility."""
    if em_config is None:
        em_config = EMConfig()
    n_loci = len(multi_loci)
    if n_loci == 0:
        return _empty_prior_table(float(em_config.rna_call_bias))
    if index.region_df is None:
        raise RuntimeError(
            "assemble_priors: index has no region table. Rebuild the index with current rigel."
        )
    region_calibration = getattr(calibration, "region_calibration", None)
    if region_calibration is None:
        raise ValueError("assemble_priors: calibration.region_calibration is required.")

    prior_mass = getattr(region_calibration, "prior_mass", None)
    if prior_mass is None:
        raise ValueError("assemble_priors: calibration.region_calibration.prior_mass is required.")

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    has_gdna_candidate = np.array(
        [enable_gdna_for_multilocus(locus, em_data) for locus in multi_loci],
        dtype=bool,
    )
    adaptive = compute_adaptive_prior(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        p_states=np.asarray(region_calibration.p_states, dtype=np.float64),
        unspliced_total=np.asarray(prior_mass.unspliced_total, dtype=np.float64),
        gdna_unspliced_mean=np.asarray(prior_mass.gdna_unspliced_mean, dtype=np.float64),
        rna_unspliced_mean=np.asarray(prior_mass.rna_unspliced_mean, dtype=np.float64),
        has_gdna_candidate=has_gdna_candidate,
        rna_call_bias=float(em_config.rna_call_bias),
        max_ess=MAX_ESS,
    )

    gdna_eff_len_unweighted = np.zeros(n_loci, dtype=np.float64)
    gdna_eff_len = np.zeros(n_loci, dtype=np.float64)
    gdna_em_exposure_weight = np.ones(n_loci, dtype=np.float64)
    for locus in multi_loci:
        locus_id = int(locus.multi_locus_id)
        unweighted = gdna_eff_len_for_loci(
            locus.loci,
            index.ref_lengths,
            calibration.fl_models.gdna,
        )
        exposure_weight = bp_weighted_mean_exposure_over_blocks(
            blocks=[(loc.ref_id, loc.start, loc.end) for loc in locus.loci],
            region_arrays=region_arrays,
            exposure=region_calibration,
        )
        gdna_eff_len_unweighted[locus_id] = float(unweighted)
        gdna_em_exposure_weight[locus_id] = float(exposure_weight)
        gdna_eff_len[locus_id] = float(max(float(unweighted) * float(exposure_weight), 1.0))

    n_units_used_for_diagnostics = np.array(
        [_count_diagnostic_units(locus, em_data) for locus in multi_loci],
        dtype=np.int64,
    )

    return PriorTable(
        alpha_gdna_add=adaptive.alpha_gdna_add,
        alpha_rna_add=adaptive.alpha_rna_add,
        gdna_eff_len=gdna_eff_len,
        gdna_eff_len_unweighted=gdna_eff_len_unweighted,
        gdna_em_exposure_weight=gdna_em_exposure_weight,
        enable_gdna=has_gdna_candidate.astype(np.uint8),
        prior_locus_weight=adaptive.locus_weight,
        prior_shrink_weight=adaptive.shrink_weight,
        prior_n_local_gdna=adaptive.n_local[:, 0],
        prior_n_local_rna=adaptive.n_local[:, 1],
        prior_n_other_gdna=adaptive.n_other[:, 0],
        prior_n_other_rna=adaptive.n_other[:, 1],
        prior_ess_final=adaptive.ess_final,
        prior_rna_share_v5=adaptive.rna_share_v5,
        prior_rna_share_final=adaptive.rna_share_final,
        prior_unspliced_total=adaptive.locus_unspliced,
        prior_flags=adaptive.flags,
        prior_region_weight=adaptive.region_weight,
        n_regions_touched=adaptive.n_regions_touched,
        n_units_used_for_diagnostics=n_units_used_for_diagnostics,
        multi_locus_region_mass=adaptive.multi_locus_region_mass,
        partial_coverage_region_mass=adaptive.partial_coverage_region_mass,
        global_n_gdna=float(adaptive.global_counts[0]),
        global_n_rna=float(adaptive.global_counts[1]),
        unallocated_unspliced_count=float(adaptive.unallocated_unspliced),
        unallocated_weighted_unspliced_count=float(adaptive.unallocated_weighted_unspliced),
        rna_call_bias=float(em_config.rna_call_bias),
    )


def enable_gdna_for_multilocus(locus: "MultiLocus", em_data: "ScoredFragments") -> bool:
    """Return True iff a locus has an unspliced unit with finite gDNA likelihood."""
    units = np.asarray(locus.unit_indices, dtype=np.int64)
    if units.size == 0:
        return False
    is_spliced = np.asarray(em_data.is_spliced, dtype=bool)
    gdna_log_liks = np.asarray(em_data.gdna_log_liks, dtype=np.float64)
    valid_units = units[(units >= 0) & (units < is_spliced.size) & (units < gdna_log_liks.size)]
    if valid_units.size == 0:
        return False
    return bool(np.any((~is_spliced[valid_units]) & np.isfinite(gdna_log_liks[valid_units])))


def _count_diagnostic_units(locus: "MultiLocus", em_data: "ScoredFragments") -> int:
    midpoint = getattr(em_data, "genomic_midpoint", None)
    if midpoint is None:
        return 0
    units = np.asarray(locus.unit_indices, dtype=np.int64)
    mid = np.asarray(midpoint, dtype=np.int64)
    valid = units[(units >= 0) & (units < mid.size)]
    if valid.size == 0:
        return 0
    return int(np.sum(mid[valid] != _INT64_MIN))


def _empty_prior_table(rna_call_bias: float = 0.5) -> PriorTable:
    return PriorTable(
        alpha_gdna_add=np.zeros(0, dtype=np.float64),
        alpha_rna_add=np.zeros(0, dtype=np.float64),
        gdna_eff_len=np.zeros(0, dtype=np.float64),
        gdna_eff_len_unweighted=np.zeros(0, dtype=np.float64),
        gdna_em_exposure_weight=np.zeros(0, dtype=np.float64),
        enable_gdna=np.zeros(0, dtype=np.uint8),
        prior_locus_weight=np.zeros(0, dtype=np.float64),
        prior_shrink_weight=np.zeros(0, dtype=np.float64),
        prior_n_local_gdna=np.zeros(0, dtype=np.float64),
        prior_n_local_rna=np.zeros(0, dtype=np.float64),
        prior_n_other_gdna=np.zeros(0, dtype=np.float64),
        prior_n_other_rna=np.zeros(0, dtype=np.float64),
        prior_ess_final=np.zeros(0, dtype=np.float64),
        prior_rna_share_v5=np.zeros(0, dtype=np.float64),
        prior_rna_share_final=np.zeros(0, dtype=np.float64),
        prior_unspliced_total=np.zeros(0, dtype=np.float64),
        prior_flags=np.zeros(0, dtype=np.uint16),
        prior_region_weight=np.zeros(0, dtype=np.float64),
        n_regions_touched=np.zeros(0, dtype=np.int32),
        n_units_used_for_diagnostics=np.zeros(0, dtype=np.int64),
        multi_locus_region_mass=np.zeros(0, dtype=np.float64),
        partial_coverage_region_mass=np.zeros(0, dtype=np.float64),
        rna_call_bias=float(rna_call_bias),
    )


def _prior_flag_histogram(values: np.ndarray) -> dict[str, int]:
    flags = np.asarray(values, dtype=np.uint16)
    return {
        "PRIOR_NO_UNSPLICED_MASS": int(np.count_nonzero(flags & PRIOR_NO_UNSPLICED_MASS)),
        "PRIOR_STRUCTURAL_GATED": int(np.count_nonzero(flags & PRIOR_STRUCTURAL_GATED)),
        "PRIOR_ESS_CAPPED": int(np.count_nonzero(flags & PRIOR_ESS_CAPPED)),
        "PRIOR_BIAS_APPLIED": int(np.count_nonzero(flags & PRIOR_BIAS_APPLIED)),
    }


def _summary_stats(values: np.ndarray) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    if arr.size == 0:
        return {"min": 0.0, "p50": 0.0, "p95": 0.0, "max": 0.0, "mean": 0.0}
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        finite = np.array([0.0], dtype=np.float64)
    return {
        "min": float(np.min(finite)),
        "p50": float(np.quantile(finite, 0.50)),
        "p95": float(np.quantile(finite, 0.95)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
    }