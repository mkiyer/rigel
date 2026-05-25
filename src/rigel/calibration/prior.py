"""Projection of fused regional gDNA evidence onto MultiLocus EM inputs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

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
    "PriorTable",
    "assemble_priors",
    "enable_gdna_for_multilocus",
]


COUNT_ALLOC_GEOMETRY: int = 1
COUNT_ALLOC_MIDPOINT: int = 2
COUNT_ALLOC_SPAN: int = 3

_INT64_MIN: int = -(2**63)
_GDNA_ENABLE_UPPER_EPS: float = 1.0e-6


@dataclass(frozen=True, slots=True)
class PriorTable:
    """Per-MultiLocus gDNA prior inputs and diagnostics."""

    gdna_prior_count_em: np.ndarray
    gdna_expected_count: np.ndarray
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
    unallocated_upper_count: float = 0.0

    def to_summary_dict(self) -> dict[str, object]:
        """Return compact JSON-safe prior diagnostics."""
        return {
            "n_loci": int(self.gdna_prior_count_em.size),
            "sum_gdna_prior_count_em": float(np.sum(self.gdna_prior_count_em)),
            "sum_gdna_expected_count": float(np.sum(self.gdna_expected_count)),
            "sum_gdna_upper_count": float(np.sum(self.gdna_upper_count)),
            "unallocated_expected_count": float(self.unallocated_expected_count),
            "unallocated_upper_count": float(self.unallocated_upper_count),
            "n_loci_enable_gdna_true": int(np.sum(self.enable_gdna != 0)),
            "n_regions_touched": _summary_stats(self.n_regions_touched),
            "gdna_eff_len": _summary_stats(self.gdna_eff_len),
            "gdna_em_exposure_weight": _summary_stats(self.gdna_em_exposure_weight),
            "multi_locus_region_mass": _summary_stats(self.multi_locus_region_mass),
            "partial_coverage_region_mass": _summary_stats(self.partial_coverage_region_mass),
            "count_allocation_mode_histogram": _uint_histogram(self.count_allocation_mode),
            "count_allocation_fallback_histogram": _uint_histogram(self.count_allocation_fallback),
        }


def assemble_priors(
    *,
    multi_loci: list["MultiLocus"],
    em_data: "ScoredFragments",
    index: "TranscriptIndex",
    calibration: "CalibrationResult",
) -> PriorTable:
    """Assemble gDNA prior counts, denominators, and eligibility for EM."""
    n_loci = len(multi_loci)
    if n_loci == 0:
        return _empty_prior_table()
    if index.region_df is None:
        raise RuntimeError(
            "assemble_priors: index has no region table. Rebuild the index with current rigel."
        )
    if calibration.fused_region_gdna is None:
        raise ValueError("assemble_priors: calibration.fused_region_gdna is required.")

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    fused = calibration.fused_region_gdna
    if np.asarray(fused.mean_count).shape != region_arrays.start.shape:
        raise ValueError("assemble_priors: fused region evidence is not aligned to index regions.")

    gdna_expected = np.zeros(n_loci, dtype=np.float64)
    gdna_upper = np.zeros(n_loci, dtype=np.float64)
    n_regions_touched = np.zeros(n_loci, dtype=np.int32)
    multi_locus_region_mass = np.zeros(n_loci, dtype=np.float64)
    partial_coverage_region_mass = np.zeros(n_loci, dtype=np.float64)

    unallocated_expected, unallocated_upper = _allocate_counts_by_geometry(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        fused_mean=np.asarray(fused.mean_count, dtype=np.float64),
        fused_upper=np.asarray(fused.upper_count, dtype=np.float64),
        gdna_expected=gdna_expected,
        gdna_upper=gdna_upper,
        n_regions_touched=n_regions_touched,
        multi_locus_region_mass=multi_locus_region_mass,
        partial_coverage_region_mass=partial_coverage_region_mass,
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
            exposure=calibration.region_exposure,
        )
        gdna_eff_len_unweighted[lid] = float(unweighted)
        gdna_em_exposure_weight[lid] = float(weight)
        gdna_eff_len[lid] = float(max(float(unweighted) * float(weight), 1.0))

    candidate_enable = np.array(
        [enable_gdna_for_multilocus(locus, em_data) for locus in multi_loci],
        dtype=bool,
    )
    evidence_enable = gdna_upper > _GDNA_ENABLE_UPPER_EPS
    enable_gdna = (candidate_enable & evidence_enable).astype(np.uint8)
    n_units_used_for_diagnostics = np.array(
        [_count_diagnostic_units(locus, em_data) for locus in multi_loci],
        dtype=np.int64,
    )

    return PriorTable(
        gdna_prior_count_em=gdna_expected.copy(),
        gdna_expected_count=gdna_expected,
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
    fused_mean: np.ndarray,
    fused_upper: np.ndarray,
    gdna_expected: np.ndarray,
    gdna_upper: np.ndarray,
    n_regions_touched: np.ndarray,
    multi_locus_region_mass: np.ndarray,
    partial_coverage_region_mass: np.ndarray,
) -> tuple[float, float]:
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
    for ref_id in range(region_arrays.n_refs):
        blocks = blocks_by_ref.get(ref_id)
        lo = int(region_arrays.ref_offsets[ref_id])
        hi = int(region_arrays.ref_offsets[ref_id + 1])
        if not blocks:
            unallocated_expected += float(np.sum(fused_mean[lo:hi]))
            unallocated_upper += float(np.sum(fused_upper[lo:hi]))
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
                unallocated_expected += float(fused_mean[region_idx])
                unallocated_upper += float(fused_upper[region_idx])
                continue

            touches_multiple = len(raw_by_locus) > 1
            partial_coverage = raw_total < 1.0 - 1.0e-9
            for lid, raw_share in raw_by_locus.items():
                share = float(raw_share / raw_total)
                expected_share = share * float(fused_mean[region_idx])
                upper_share = share * float(fused_upper[region_idx])
                gdna_expected[lid] += expected_share
                gdna_upper[lid] += upper_share
                n_regions_touched[lid] += 1
                if touches_multiple:
                    multi_locus_region_mass[lid] += expected_share
                if partial_coverage:
                    partial_coverage_region_mass[lid] += expected_share

    return unallocated_expected, unallocated_upper


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
