"""Projection of v6 regional calibration evidence onto MultiLocus EM inputs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..config import EMConfig
from ._arrays import RegionArrays
from ._exposure import bp_weighted_mean_exposure_over_blocks, gdna_eff_len_for_loci

if TYPE_CHECKING:  # pragma: no cover - annotations only.
    from ..index import TranscriptIndex
    from ..locus import MultiLocus
    from ..scored_fragments import ScoredFragments
    from ._result import CalibrationResult


__all__ = [
    "COUNT_ALLOC_GEOMETRY",
    "COUNT_ALLOC_MIDPOINT",
    "COUNT_ALLOC_SPAN",
    "GroupedPriorCounts",
    "PriorTable",
    "_compute_grouped_prior_counts",
    "assemble_priors",
    "enable_gdna_for_multilocus",
]


COUNT_ALLOC_GEOMETRY: int = 1
COUNT_ALLOC_MIDPOINT: int = 2
COUNT_ALLOC_SPAN: int = 3

_INT64_MIN: int = -(2**63)
_EPS: float = 1.0e-12


@dataclass(frozen=True, slots=True)
class GroupedPriorCounts:
    """Paired additive grouped-prior counts and budget diagnostics."""

    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray
    budget_raw: np.ndarray
    budget: np.ndarray
    gdna_share_raw: np.ndarray
    gdna_share_biased: np.ndarray


@dataclass(frozen=True, slots=True)
class PriorTable:
    """Per-MultiLocus gDNA prior inputs and diagnostics."""

    gdna_prior_count_em: np.ndarray
    gdna_expected_count: np.ndarray
    rna_expected_count: np.ndarray
    prior_unspliced_total: np.ndarray
    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray
    prior_budget_raw: np.ndarray
    prior_budget: np.ndarray
    prior_gdna_share_raw: np.ndarray
    prior_gdna_share_biased: np.ndarray
    prior_mass_conservation_error: np.ndarray
    prior_allocated_fraction: np.ndarray
    gdna_prior_density: np.ndarray
    gdna_upper_count: np.ndarray
    gdna_eff_len: np.ndarray
    gdna_eff_len_unweighted: np.ndarray
    gdna_em_exposure_weight: np.ndarray
    enable_gdna: np.ndarray
    n_regions_touched: np.ndarray
    n_units_used_for_diagnostics: np.ndarray
    count_allocation_mode: np.ndarray
    count_allocation_fallback: np.ndarray
    multi_locus_region_mass: np.ndarray
    partial_coverage_region_mass: np.ndarray
    unallocated_expected_count: float = 0.0
    unallocated_rna_expected_count: float = 0.0
    unallocated_unspliced_count: float = 0.0
    unallocated_upper_count: float = 0.0

    def to_summary_dict(self) -> dict[str, object]:
        """Return compact JSON-safe prior diagnostics."""
        return {
            "n_loci": int(self.gdna_prior_count_em.size),
            "sum_gdna_prior_count_em": float(np.sum(self.gdna_prior_count_em)),
            "sum_gdna_expected_count": float(np.sum(self.gdna_expected_count)),
            "sum_rna_expected_count": float(np.sum(self.rna_expected_count)),
            "sum_prior_unspliced_total": float(np.sum(self.prior_unspliced_total)),
            "sum_alpha_gdna_add": float(np.sum(self.alpha_gdna_add)),
            "sum_alpha_rna_add": float(np.sum(self.alpha_rna_add)),
            "sum_gdna_upper_count": float(np.sum(self.gdna_upper_count)),
            "unallocated_expected_count": float(self.unallocated_expected_count),
            "unallocated_rna_expected_count": float(self.unallocated_rna_expected_count),
            "unallocated_unspliced_count": float(self.unallocated_unspliced_count),
            "unallocated_upper_count": float(self.unallocated_upper_count),
            "n_loci_enable_gdna_true": int(np.sum(self.enable_gdna != 0)),
            "n_regions_touched": _summary_stats(self.n_regions_touched),
            "prior_budget": _summary_stats(self.prior_budget),
            "prior_gdna_share_biased": _summary_stats(self.prior_gdna_share_biased),
            "prior_mass_conservation_error": _summary_stats(self.prior_mass_conservation_error),
            "prior_allocated_fraction": _summary_stats(self.prior_allocated_fraction),
            "gdna_prior_density": _summary_stats(self.gdna_prior_density),
            "gdna_eff_len": _summary_stats(self.gdna_eff_len),
            "gdna_em_exposure_weight": _summary_stats(self.gdna_em_exposure_weight),
            "multi_locus_region_mass": _summary_stats(self.multi_locus_region_mass),
            "partial_coverage_region_mass": _summary_stats(self.partial_coverage_region_mass),
            "count_allocation_mode_histogram": _uint_histogram(self.count_allocation_mode),
            "count_allocation_fallback_histogram": _uint_histogram(self.count_allocation_fallback),
        }


def _compute_grouped_prior_counts(
    *,
    gdna_expected: np.ndarray,
    rna_expected: np.ndarray,
    em_config: EMConfig,
) -> GroupedPriorCounts:
    """Compute bounded paired additive pseudocounts from projected mass."""
    gdna = np.maximum(np.asarray(gdna_expected, dtype=np.float64), 0.0)
    rna = np.maximum(np.asarray(rna_expected, dtype=np.float64), 0.0)
    if gdna.shape != rna.shape:
        raise ValueError(
            "_compute_grouped_prior_counts: gdna_expected and rna_expected shapes must match; "
            f"got {gdna.shape} and {rna.shape}."
        )

    total = gdna + rna
    share_raw = np.divide(gdna, total, out=np.zeros_like(gdna), where=total > _EPS)

    balanced_budget = 2.0 * np.minimum(gdna, rna)
    edge_budget = np.minimum(total, float(em_config.aggregate_prior_edge_count))
    budget_raw = np.maximum(balanced_budget, edge_budget)
    budget_raw = np.where(total > _EPS, budget_raw, 0.0)
    budget_raw = np.minimum(budget_raw, float(em_config.aggregate_prior_max_count))
    budget = float(em_config.aggregate_prior_strength) * budget_raw
    budget = np.minimum(budget, float(em_config.aggregate_prior_max_count))

    eps = 1.0e-15
    clipped = np.clip(share_raw, eps, 1.0 - eps)
    logits = np.log(clipped) - np.log1p(-clipped) + float(em_config.gdna_prior_logit_bias)
    share_biased = np.empty_like(logits)
    pos = logits >= 0.0
    share_biased[pos] = 1.0 / (1.0 + np.exp(-logits[pos]))
    exp_logits = np.exp(logits[~pos])
    share_biased[~pos] = exp_logits / (1.0 + exp_logits)
    share_biased = np.where(total > _EPS, share_biased, 0.0)

    alpha_gdna = budget * share_biased
    alpha_rna = budget * (1.0 - share_biased)
    alpha_gdna = np.where(budget > 0.0, alpha_gdna, 0.0)
    alpha_rna = np.where(budget > 0.0, alpha_rna, 0.0)

    return GroupedPriorCounts(
        alpha_gdna_add=alpha_gdna.astype(np.float64, copy=False),
        alpha_rna_add=alpha_rna.astype(np.float64, copy=False),
        budget_raw=budget_raw.astype(np.float64, copy=False),
        budget=budget.astype(np.float64, copy=False),
        gdna_share_raw=share_raw.astype(np.float64, copy=False),
        gdna_share_biased=share_biased.astype(np.float64, copy=False),
    )


def assemble_priors(
    *,
    multi_loci: list["MultiLocus"],
    em_data: "ScoredFragments",
    index: "TranscriptIndex",
    calibration: "CalibrationResult",
    em_config: EMConfig | None = None,
) -> PriorTable:
    """Assemble gDNA prior counts, denominators, and native EM eligibility."""
    if em_config is None:
        em_config = EMConfig()
    n_loci = len(multi_loci)
    if n_loci == 0:
        return _empty_prior_table()
    if index.region_df is None:
        raise RuntimeError(
            "assemble_priors: index has no region table. Rebuild the index with current rigel."
        )
    region_calibration = getattr(calibration, "region_calibration", None)
    if region_calibration is None:
        raise ValueError("assemble_priors: calibration.region_calibration is required.")

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    prior_mass = getattr(region_calibration, "prior_mass", None)
    if prior_mass is None:
        raise ValueError("assemble_priors: calibration.region_calibration.prior_mass is required.")
    region_mean = np.asarray(prior_mass.gdna_unspliced_mean, dtype=np.float64)
    region_rna = np.asarray(prior_mass.rna_unspliced_mean, dtype=np.float64)
    region_unspliced = np.asarray(prior_mass.unspliced_total, dtype=np.float64)
    region_upper = np.asarray(region_calibration.upper_gdna, dtype=np.float64)
    if (
        region_mean.shape != region_arrays.start.shape
        or region_rna.shape != region_arrays.start.shape
        or region_unspliced.shape != region_arrays.start.shape
        or region_upper.shape != region_arrays.start.shape
    ):
        raise ValueError("assemble_priors: RegionCalibration arrays are not aligned to index regions.")

    gdna_expected = np.zeros(n_loci, dtype=np.float64)
    rna_expected = np.zeros(n_loci, dtype=np.float64)
    prior_unspliced_total = np.zeros(n_loci, dtype=np.float64)
    gdna_upper = np.zeros(n_loci, dtype=np.float64)
    n_regions_touched = np.zeros(n_loci, dtype=np.int32)
    multi_locus_region_mass = np.zeros(n_loci, dtype=np.float64)
    partial_coverage_region_mass = np.zeros(n_loci, dtype=np.float64)

    (
        unallocated_expected,
        unallocated_upper,
        unallocated_rna,
        unallocated_unspliced,
    ) = _allocate_counts_by_geometry(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        region_mean=region_mean,
        region_rna=region_rna,
        region_unspliced=region_unspliced,
        region_upper=region_upper,
        gdna_expected=gdna_expected,
        rna_expected=rna_expected,
        prior_unspliced_total=prior_unspliced_total,
        gdna_upper=gdna_upper,
        n_regions_touched=n_regions_touched,
        multi_locus_region_mass=multi_locus_region_mass,
        partial_coverage_region_mass=partial_coverage_region_mass,
    )

    grouped = _compute_grouped_prior_counts(
        gdna_expected=gdna_expected,
        rna_expected=rna_expected,
        em_config=em_config,
    )

    gdna_eff_len_unweighted = np.zeros(n_loci, dtype=np.float64)
    gdna_eff_len = np.zeros(n_loci, dtype=np.float64)
    gdna_em_exposure_weight = np.ones(n_loci, dtype=np.float64)
    for locus in multi_loci:
        lid = int(locus.multi_locus_id)
        unweighted = gdna_eff_len_for_loci(
            locus.loci,
            index.ref_lengths,
            calibration.fl_models.gdna,
        )
        weight = bp_weighted_mean_exposure_over_blocks(
            blocks=[(loc.ref_id, loc.start, loc.end) for loc in locus.loci],
            region_arrays=region_arrays,
            exposure=region_calibration,
        )
        gdna_eff_len_unweighted[lid] = float(unweighted)
        gdna_em_exposure_weight[lid] = float(weight)
        gdna_eff_len[lid] = float(max(float(unweighted) * float(weight), 1.0))

    candidate_enable = np.array(
        [enable_gdna_for_multilocus(locus, em_data) for locus in multi_loci],
        dtype=bool,
    )
    enable_gdna = candidate_enable.astype(np.uint8)
    alpha_gdna_add = np.where(candidate_enable, grouped.alpha_gdna_add, 0.0).astype(
        np.float64,
        copy=False,
    )
    alpha_rna_add = np.where(candidate_enable, grouped.alpha_rna_add, 0.0).astype(
        np.float64,
        copy=False,
    )
    prior_budget_raw = np.where(candidate_enable, grouped.budget_raw, 0.0).astype(
        np.float64,
        copy=False,
    )
    prior_budget = np.where(candidate_enable, grouped.budget, 0.0).astype(
        np.float64,
        copy=False,
    )
    prior_gdna_share_biased = np.where(candidate_enable, grouped.gdna_share_biased, 0.0).astype(
        np.float64,
        copy=False,
    )
    n_units_used_for_diagnostics = np.array(
        [_count_diagnostic_units(locus, em_data) for locus in multi_loci],
        dtype=np.int64,
    )
    prior_mass_conservation_error = np.abs(
        gdna_expected + rna_expected - prior_unspliced_total
    ).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        prior_allocated_fraction = np.where(
            prior_unspliced_total > 0.0,
            (gdna_expected + rna_expected) / prior_unspliced_total,
            1.0,
        )
        gdna_prior_density = np.where(
            gdna_eff_len > 0.0,
            alpha_gdna_add / gdna_eff_len,
            0.0,
        )

    return PriorTable(
        gdna_prior_count_em=alpha_gdna_add.copy(),
        gdna_expected_count=gdna_expected,
        rna_expected_count=rna_expected,
        prior_unspliced_total=prior_unspliced_total,
        alpha_gdna_add=alpha_gdna_add,
        alpha_rna_add=alpha_rna_add,
        prior_budget_raw=prior_budget_raw,
        prior_budget=prior_budget,
        prior_gdna_share_raw=grouped.gdna_share_raw,
        prior_gdna_share_biased=prior_gdna_share_biased,
        prior_mass_conservation_error=prior_mass_conservation_error,
        prior_allocated_fraction=prior_allocated_fraction.astype(np.float64, copy=False),
        gdna_prior_density=gdna_prior_density.astype(np.float64, copy=False),
        gdna_upper_count=gdna_upper,
        gdna_eff_len=gdna_eff_len,
        gdna_eff_len_unweighted=gdna_eff_len_unweighted,
        gdna_em_exposure_weight=gdna_em_exposure_weight,
        enable_gdna=enable_gdna,
        n_regions_touched=n_regions_touched,
        n_units_used_for_diagnostics=n_units_used_for_diagnostics,
        count_allocation_mode=np.full(n_loci, COUNT_ALLOC_GEOMETRY, dtype=np.uint8),
        count_allocation_fallback=np.zeros(n_loci, dtype=np.uint8),
        multi_locus_region_mass=multi_locus_region_mass,
        partial_coverage_region_mass=partial_coverage_region_mass,
        unallocated_expected_count=float(unallocated_expected),
        unallocated_rna_expected_count=float(unallocated_rna),
        unallocated_unspliced_count=float(unallocated_unspliced),
        unallocated_upper_count=float(unallocated_upper),
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


def _allocate_counts_by_geometry(
    *,
    region_arrays: RegionArrays,
    multi_loci: list["MultiLocus"],
    region_mean: np.ndarray,
    region_rna: np.ndarray,
    region_unspliced: np.ndarray,
    region_upper: np.ndarray,
    gdna_expected: np.ndarray,
    rna_expected: np.ndarray,
    prior_unspliced_total: np.ndarray,
    gdna_upper: np.ndarray,
    n_regions_touched: np.ndarray,
    multi_locus_region_mass: np.ndarray,
    partial_coverage_region_mass: np.ndarray,
) -> tuple[float, float, float, float]:
    blocks_by_ref: dict[int, list[tuple[int, int, int]]] = {}
    for locus in multi_loci:
        lid = int(locus.multi_locus_id)
        for block in locus.loci:
            if int(block.end) <= int(block.start):
                continue
            blocks_by_ref.setdefault(int(block.ref_id), []).append(
                (int(block.start), int(block.end), lid)
            )
    for blocks in blocks_by_ref.values():
        blocks.sort(key=lambda item: (item[0], item[1], item[2]))

    unallocated_expected = 0.0
    unallocated_upper = 0.0
    unallocated_rna = 0.0
    unallocated_unspliced = 0.0
    for ref_id in range(region_arrays.n_refs):
        blocks = blocks_by_ref.get(ref_id)
        lo = int(region_arrays.ref_offsets[ref_id])
        hi = int(region_arrays.ref_offsets[ref_id + 1])
        if not blocks:
            unallocated_expected += float(np.sum(region_mean[lo:hi]))
            unallocated_upper += float(np.sum(region_upper[lo:hi]))
            unallocated_rna += float(np.sum(region_rna[lo:hi]))
            unallocated_unspliced += float(np.sum(region_unspliced[lo:hi]))
            continue

        block_starts = np.array([b[0] for b in blocks], dtype=np.int64)
        for region_idx in range(lo, hi):
            r_start = int(region_arrays.start[region_idx])
            r_end = int(region_arrays.end[region_idx])
            r_len = max(r_end - r_start, 0)
            if r_len <= 0:
                continue

            candidate_hi = int(np.searchsorted(block_starts, r_end, side="left"))
            raw_by_locus: dict[int, float] = {}
            for b_start, b_end, lid in blocks[:candidate_hi]:
                if b_end <= r_start:
                    continue
                overlap = min(b_end, r_end) - max(b_start, r_start)
                if overlap <= 0:
                    continue
                raw_by_locus[lid] = raw_by_locus.get(lid, 0.0) + float(overlap) / float(r_len)

            raw_total = float(sum(raw_by_locus.values()))
            if raw_total <= 0.0:
                unallocated_expected += float(region_mean[region_idx])
                unallocated_upper += float(region_upper[region_idx])
                unallocated_rna += float(region_rna[region_idx])
                unallocated_unspliced += float(region_unspliced[region_idx])
                continue

            touches_multiple = len(raw_by_locus) > 1
            partial_coverage = raw_total < 1.0 - 1.0e-9
            for lid, raw_share in raw_by_locus.items():
                share = float(raw_share / raw_total)
                expected_share = share * float(region_mean[region_idx])
                rna_share = share * float(region_rna[region_idx])
                unspliced_share = share * float(region_unspliced[region_idx])
                upper_share = share * float(region_upper[region_idx])
                gdna_expected[lid] += expected_share
                rna_expected[lid] += rna_share
                prior_unspliced_total[lid] += unspliced_share
                gdna_upper[lid] += upper_share
                n_regions_touched[lid] += 1
                if touches_multiple:
                    multi_locus_region_mass[lid] += expected_share
                if partial_coverage:
                    partial_coverage_region_mass[lid] += expected_share

    return unallocated_expected, unallocated_upper, unallocated_rna, unallocated_unspliced


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


def _empty_prior_table() -> PriorTable:
    return PriorTable(
        gdna_prior_count_em=np.zeros(0, dtype=np.float64),
        gdna_expected_count=np.zeros(0, dtype=np.float64),
        rna_expected_count=np.zeros(0, dtype=np.float64),
        prior_unspliced_total=np.zeros(0, dtype=np.float64),
        alpha_gdna_add=np.zeros(0, dtype=np.float64),
        alpha_rna_add=np.zeros(0, dtype=np.float64),
        prior_budget_raw=np.zeros(0, dtype=np.float64),
        prior_budget=np.zeros(0, dtype=np.float64),
        prior_gdna_share_raw=np.zeros(0, dtype=np.float64),
        prior_gdna_share_biased=np.zeros(0, dtype=np.float64),
        prior_mass_conservation_error=np.zeros(0, dtype=np.float64),
        prior_allocated_fraction=np.zeros(0, dtype=np.float64),
        gdna_prior_density=np.zeros(0, dtype=np.float64),
        gdna_upper_count=np.zeros(0, dtype=np.float64),
        gdna_eff_len=np.zeros(0, dtype=np.float64),
        gdna_eff_len_unweighted=np.zeros(0, dtype=np.float64),
        gdna_em_exposure_weight=np.zeros(0, dtype=np.float64),
        enable_gdna=np.zeros(0, dtype=np.uint8),
        n_regions_touched=np.zeros(0, dtype=np.int32),
        n_units_used_for_diagnostics=np.zeros(0, dtype=np.int64),
        count_allocation_mode=np.zeros(0, dtype=np.uint8),
        count_allocation_fallback=np.zeros(0, dtype=np.uint8),
        multi_locus_region_mass=np.zeros(0, dtype=np.float64),
        partial_coverage_region_mass=np.zeros(0, dtype=np.float64),
    )


def _uint_histogram(values: np.ndarray) -> dict[int, int]:
    arr = np.asarray(values, dtype=np.uint8)
    return {int(v): int(np.sum(arr == v)) for v in np.unique(arr)}


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
