"""Entropy-weighted grouped EM priors for the adaptive v5/v6 cutover."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ._arrays import RegionArrays
from .latent_states import N_STATES, STATE_IS_EXPRESSED

if TYPE_CHECKING:  # pragma: no cover - annotations only.
    from ..locus import MultiLocus


__all__ = [
    "EPS_MASS",
    "MAX_ESS",
    "PRIOR_BIAS_APPLIED",
    "PRIOR_ESS_CAPPED",
    "PRIOR_NO_UNSPLICED_MASS",
    "PRIOR_STRUCTURAL_GATED",
    "AdaptivePriorResult",
    "compute_adaptive_prior",
]


MAX_ESS: float = 3000.0
EPS_MASS: float = 1.0e-12

PRIOR_NO_UNSPLICED_MASS: np.uint16 = np.uint16(0x1)
PRIOR_STRUCTURAL_GATED: np.uint16 = np.uint16(0x2)
PRIOR_ESS_CAPPED: np.uint16 = np.uint16(0x4)
PRIOR_BIAS_APPLIED: np.uint16 = np.uint16(0x8)

_PROB_TOL: float = 1.0e-8


@dataclass(frozen=True, slots=True)
class AdaptivePriorResult:
    """Outputs of the v5 entropy-Dirichlet prior plus the v6 split dial."""

    alpha_gdna_add: np.ndarray
    alpha_rna_add: np.ndarray
    n_local: np.ndarray
    n_other: np.ndarray
    locus_weight: np.ndarray
    shrink_weight: np.ndarray
    global_counts: np.ndarray
    region_weight: np.ndarray
    locus_unspliced: np.ndarray
    n_regions_touched: np.ndarray
    multi_locus_region_mass: np.ndarray
    partial_coverage_region_mass: np.ndarray
    unallocated_unspliced: float
    unallocated_weighted_unspliced: float
    rna_share_v5: np.ndarray
    rna_share_final: np.ndarray
    ess_final: np.ndarray
    flags: np.ndarray


@dataclass(frozen=True, slots=True)
class _ProjectionResult:
    values: dict[str, np.ndarray]
    unallocated: dict[str, float]
    n_regions_touched: np.ndarray
    multi_locus_region_mass: np.ndarray
    partial_coverage_region_mass: np.ndarray


def compute_adaptive_prior(
    *,
    region_arrays: RegionArrays,
    multi_loci: list["MultiLocus"],
    p_states: np.ndarray,
    unspliced_total: np.ndarray,
    has_gdna_candidate: np.ndarray,
    rna_call_bias: float = 0.5,
    max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    """Compute v5/v6 grouped additive priors for every MultiLocus.

    The output arrays are indexed by ``multi_locus_id``. ``multi_loci`` may be
    provided in any order, but ids must be contiguous in ``0..L-1``.
    """
    n_loci = len(multi_loci)
    _validate_locus_ids(multi_loci, n_loci)
    p, unspliced, has_candidate, bias, cap_max = _validate_inputs(
        region_arrays=region_arrays,
        n_loci=n_loci,
        p_states=p_states,
        unspliced_total=unspliced_total,
        has_gdna_candidate=has_gdna_candidate,
        rna_call_bias=rna_call_bias,
        max_ess=max_ess,
    )

    region_weight = _entropy_weight(p)
    rna_mask = np.asarray(STATE_IS_EXPRESSED, dtype=bool)
    gdna_mask = ~rna_mask

    q_gdna = np.sum(p[:, gdna_mask], axis=1, dtype=np.float64)
    q_rna = np.sum(p[:, rna_mask], axis=1, dtype=np.float64)
    q_total = q_gdna + q_rna
    q_gdna = np.divide(q_gdna, q_total, out=np.zeros_like(q_gdna), where=q_total > 0.0)
    q_rna = np.divide(q_rna, q_total, out=np.zeros_like(q_rna), where=q_total > 0.0)

    weighted_unspliced = unspliced * region_weight
    n_gdna_region = weighted_unspliced * q_gdna
    n_rna_region = weighted_unspliced * q_rna

    projection = _project_to_loci(
        region_arrays=region_arrays,
        multi_loci=multi_loci,
        n_loci=n_loci,
        region_values={
            "unspliced": unspliced,
            "n_gdna": n_gdna_region,
            "n_rna": n_rna_region,
            "weighted_unspliced": weighted_unspliced,
        },
        diagnostic_mass_key="unspliced",
    )
    locus_unspliced = projection.values["unspliced"]
    n_local_gdna = projection.values["n_gdna"]
    n_local_rna = projection.values["n_rna"]
    weighted_locus_unspliced = projection.values["weighted_unspliced"]

    locus_weight = np.divide(
        weighted_locus_unspliced,
        locus_unspliced,
        out=np.zeros_like(weighted_locus_unspliced),
        where=locus_unspliced > EPS_MASS,
    )
    locus_weight = np.clip(locus_weight, 0.0, 1.0)
    shrink_weight = 1.0 - locus_weight

    global_counts = np.array(
        [float(np.sum(n_local_gdna, dtype=np.float64)), float(np.sum(n_local_rna, dtype=np.float64))],
        dtype=np.float64,
    )
    n_other_gdna = np.maximum(global_counts[0] - n_local_gdna, 0.0)
    n_other_rna = np.maximum(global_counts[1] - n_local_rna, 0.0)

    alpha_gdna = n_local_gdna + shrink_weight * n_other_gdna
    alpha_rna = n_local_rna + shrink_weight * n_other_rna

    total_before_cap = alpha_gdna + alpha_rna
    cap = np.minimum(locus_unspliced, cap_max)
    capped = total_before_cap > cap
    scale = np.divide(
        cap,
        total_before_cap,
        out=np.ones_like(total_before_cap),
        where=capped & (total_before_cap > 0.0),
    )
    alpha_gdna = alpha_gdna * scale
    alpha_rna = alpha_rna * scale

    has_mass = locus_unspliced > EPS_MASS
    enabled = has_candidate & has_mass
    alpha_gdna = np.where(enabled, alpha_gdna, 0.0).astype(np.float64, copy=False)
    alpha_rna = np.where(enabled, alpha_rna, 0.0).astype(np.float64, copy=False)

    total_v5 = alpha_gdna + alpha_rna
    rna_share_v5 = np.divide(
        alpha_rna,
        total_v5,
        out=np.zeros_like(total_v5),
        where=total_v5 > 0.0,
    )

    bias_applied = np.zeros(n_loci, dtype=bool)
    if bias != 0.5:
        shifted_share = _apply_rna_call_bias(rna_share_v5, bias)
        has_prior = total_v5 > 0.0
        alpha_gdna = np.where(has_prior, total_v5 * (1.0 - shifted_share), alpha_gdna)
        alpha_rna = np.where(has_prior, total_v5 * shifted_share, alpha_rna)
        bias_applied = has_prior

    ess_final = alpha_gdna + alpha_rna
    rna_share_final = np.divide(
        alpha_rna,
        ess_final,
        out=np.zeros_like(ess_final),
        where=ess_final > 0.0,
    )

    flags = np.zeros(n_loci, dtype=np.uint16)
    flags |= np.where(~has_mass, PRIOR_NO_UNSPLICED_MASS, np.uint16(0)).astype(np.uint16)
    flags |= np.where(~has_candidate, PRIOR_STRUCTURAL_GATED, np.uint16(0)).astype(np.uint16)
    flags |= np.where(capped, PRIOR_ESS_CAPPED, np.uint16(0)).astype(np.uint16)
    flags |= np.where(bias_applied, PRIOR_BIAS_APPLIED, np.uint16(0)).astype(np.uint16)

    return AdaptivePriorResult(
        alpha_gdna_add=alpha_gdna,
        alpha_rna_add=alpha_rna,
        n_local=np.column_stack((n_local_gdna, n_local_rna)).astype(np.float64, copy=False),
        n_other=np.column_stack((n_other_gdna, n_other_rna)).astype(np.float64, copy=False),
        locus_weight=locus_weight.astype(np.float64, copy=False),
        shrink_weight=shrink_weight.astype(np.float64, copy=False),
        global_counts=global_counts,
        region_weight=region_weight.astype(np.float64, copy=False),
        locus_unspliced=locus_unspliced.astype(np.float64, copy=False),
        n_regions_touched=projection.n_regions_touched.astype(np.int32, copy=False),
        multi_locus_region_mass=projection.multi_locus_region_mass.astype(np.float64, copy=False),
        partial_coverage_region_mass=projection.partial_coverage_region_mass.astype(
            np.float64,
            copy=False,
        ),
        unallocated_unspliced=float(projection.unallocated["unspliced"]),
        unallocated_weighted_unspliced=float(projection.unallocated["weighted_unspliced"]),
        rna_share_v5=rna_share_v5.astype(np.float64, copy=False),
        rna_share_final=rna_share_final.astype(np.float64, copy=False),
        ess_final=ess_final.astype(np.float64, copy=False),
        flags=flags,
    )


def _validate_locus_ids(multi_loci: list["MultiLocus"], n_loci: int) -> None:
    ids = sorted(int(locus.multi_locus_id) for locus in multi_loci)
    if ids != list(range(n_loci)):
        raise ValueError("multi_loci must have contiguous multi_locus_id values 0..L-1.")


def _validate_inputs(
    *,
    region_arrays: RegionArrays,
    n_loci: int,
    p_states: np.ndarray,
    unspliced_total: np.ndarray,
    has_gdna_candidate: np.ndarray,
    rna_call_bias: float,
    max_ess: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    n_regions = int(region_arrays.start.shape[0])
    p = np.asarray(p_states, dtype=np.float64)
    if p.shape != (n_regions, N_STATES):
        raise ValueError(f"p_states must have shape ({n_regions}, {N_STATES}); got {p.shape}.")
    if not np.all(np.isfinite(p)):
        raise ValueError("p_states must contain only finite values.")
    if np.any(p < -_PROB_TOL):
        raise ValueError("p_states must be nonnegative within floating-point tolerance.")
    p = np.maximum(p, 0.0)
    row_sum = np.sum(p, axis=1, dtype=np.float64)
    if np.any(row_sum <= 0.0):
        raise ValueError("every p_states row must have positive probability mass.")
    p = np.divide(p, row_sum[:, None], out=np.zeros_like(p), where=row_sum[:, None] > 0.0)

    unspliced = np.asarray(unspliced_total, dtype=np.float64)
    if unspliced.shape != (n_regions,):
        raise ValueError(
            f"unspliced_total must have shape ({n_regions},); got {unspliced.shape}."
        )
    if not np.all(np.isfinite(unspliced)):
        raise ValueError("unspliced_total must contain only finite values.")
    if np.any(unspliced < 0.0):
        raise ValueError("unspliced_total must be nonnegative.")

    has_candidate = np.asarray(has_gdna_candidate, dtype=bool)
    if has_candidate.shape != (n_loci,):
        raise ValueError(
            f"has_gdna_candidate must have shape ({n_loci},); got {has_candidate.shape}."
        )

    bias = float(rna_call_bias)
    if not math.isfinite(bias) or not (0.0 < bias < 1.0):
        raise ValueError(f"rna_call_bias must be finite and in (0, 1); got {rna_call_bias!r}.")

    cap_max = float(max_ess)
    if not math.isfinite(cap_max) or cap_max <= 0.0:
        raise ValueError(f"max_ess must be finite and > 0; got {max_ess!r}.")

    return p, unspliced, has_candidate, bias, cap_max


def _entropy_weight(p_states: np.ndarray) -> np.ndarray:
    p = np.asarray(p_states, dtype=np.float64)
    xlogx = np.zeros_like(p)
    positive = p > 0.0
    xlogx[positive] = p[positive] * np.log(p[positive])
    entropy = -np.sum(xlogx, axis=1, dtype=np.float64)
    weight = 1.0 - entropy / math.log(float(N_STATES))
    return np.clip(weight, 0.0, 1.0)


def _apply_rna_call_bias(rna_share: np.ndarray, rna_call_bias: float) -> np.ndarray:
    two_w = 2.0 * float(rna_call_bias)
    two_omw = 2.0 * (1.0 - float(rna_call_bias))
    numerator = two_w * rna_share
    denominator = numerator + two_omw * (1.0 - rna_share)
    return np.divide(numerator, denominator, out=np.zeros_like(rna_share), where=denominator > 0.0)


def _project_to_loci(
    *,
    region_arrays: RegionArrays,
    multi_loci: list["MultiLocus"],
    n_loci: int,
    region_values: dict[str, np.ndarray],
    diagnostic_mass_key: str,
) -> _ProjectionResult:
    n_regions = int(region_arrays.start.shape[0])
    values = {
        name: _validate_region_value(name, value, n_regions) for name, value in region_values.items()
    }
    if diagnostic_mass_key not in values:
        raise ValueError(f"diagnostic_mass_key {diagnostic_mass_key!r} is not in region_values.")

    out = {name: np.zeros(n_loci, dtype=np.float64) for name in values}
    unallocated = {name: 0.0 for name in values}
    n_regions_touched = np.zeros(n_loci, dtype=np.int32)
    multi_locus_region_mass = np.zeros(n_loci, dtype=np.float64)
    partial_coverage_region_mass = np.zeros(n_loci, dtype=np.float64)

    blocks_by_ref: dict[int, list[tuple[int, int, int]]] = {}
    for multi_locus in multi_loci:
        locus_id = int(multi_locus.multi_locus_id)
        for block in multi_locus.loci:
            start = int(block.start)
            end = int(block.end)
            if end <= start:
                continue
            blocks_by_ref.setdefault(int(block.ref_id), []).append((start, end, locus_id))
    for blocks in blocks_by_ref.values():
        blocks.sort(key=lambda item: (item[0], item[1], item[2]))

    for ref_id in range(int(region_arrays.n_refs)):
        lo = int(region_arrays.ref_offsets[ref_id])
        hi = int(region_arrays.ref_offsets[ref_id + 1])
        blocks = blocks_by_ref.get(ref_id)
        if not blocks:
            for name, array in values.items():
                unallocated[name] += float(np.sum(array[lo:hi], dtype=np.float64))
            continue

        block_starts = np.array([block[0] for block in blocks], dtype=np.int64)
        for region_idx in range(lo, hi):
            region_start = int(region_arrays.start[region_idx])
            region_end = int(region_arrays.end[region_idx])
            region_len = region_end - region_start
            if region_len <= 0:
                continue

            candidate_hi = int(np.searchsorted(block_starts, region_end, side="left"))
            raw_by_locus: dict[int, float] = {}
            for block_start, block_end, locus_id in blocks[:candidate_hi]:
                if block_end <= region_start:
                    continue
                overlap = min(block_end, region_end) - max(block_start, region_start)
                if overlap <= 0:
                    continue
                raw_by_locus[locus_id] = raw_by_locus.get(locus_id, 0.0) + (
                    float(overlap) / float(region_len)
                )

            raw_total = float(sum(raw_by_locus.values()))
            if raw_total <= 0.0:
                for name, array in values.items():
                    unallocated[name] += float(array[region_idx])
                continue

            touches_multiple_loci = len(raw_by_locus) > 1
            partial_coverage = raw_total < 1.0 - 1.0e-9
            diagnostic_mass = float(values[diagnostic_mass_key][region_idx])
            for locus_id, raw_share in raw_by_locus.items():
                share = float(raw_share / raw_total)
                for name, array in values.items():
                    out[name][locus_id] += share * float(array[region_idx])
                n_regions_touched[locus_id] += 1
                if touches_multiple_loci:
                    multi_locus_region_mass[locus_id] += share * diagnostic_mass
                if partial_coverage:
                    partial_coverage_region_mass[locus_id] += share * diagnostic_mass

    return _ProjectionResult(
        values=out,
        unallocated=unallocated,
        n_regions_touched=n_regions_touched,
        multi_locus_region_mass=multi_locus_region_mass,
        partial_coverage_region_mass=partial_coverage_region_mass,
    )


def _validate_region_value(name: str, value: np.ndarray, n_regions: int) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (n_regions,):
        raise ValueError(f"region value {name!r} must have shape ({n_regions},); got {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"region value {name!r} must contain only finite values.")
    return array