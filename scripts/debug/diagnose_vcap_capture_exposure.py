"""Audit VCaP hybrid-capture regional gDNA exposure calibration.

This script is intentionally read-only with respect to production code. It
rescans the VCaP mixed BAM to recover the calibration payload, rebuilds the
current v4.3 regional exposure field, and then compares alternative shrinkage
and reference-density choices for the same evidence.

Run:
    conda activate rigel && python scripts/debug/diagnose_vcap_capture_exposure.py
"""

from __future__ import annotations

import argparse
import gc
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from rigel.calibration import calibrate
from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration._exposure import (
    boundary_crossing_exposure,
    footprint_exposure_weight,
)
from rigel.calibration._kappa import estimate_kappa
from rigel.calibration._orient import ORIENT_OPP, ORIENT_SAME, StrandSummary
from rigel.calibration._regional_exposure import LOG_A_FLOOR, RegionalGdnaExposure
from rigel.calibration.density_global import (
    _gdna_count_moment,
    _strand_identifiable_rows,
    l_eff_contained,
    strand_correction_usable,
)
from rigel.calibration.regions import RegionType
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer


DEFAULT_BAM = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam")
DEFAULT_INDEX = Path("/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index")
DEFAULT_QUANT_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/v4_3_with_mm")
DEFAULT_OUT_DIR = Path("results/vcap_capture_exposure_diagnostics_2026-05-19")

CLASS_INTERGENIC = "INTERGENIC"
CLASS_INTRON = "INTRON"
CLASS_EXON = "EXON-INTRON"
CLASS_ORDER = (CLASS_INTERGENIC, CLASS_INTRON, CLASS_EXON)


@dataclass(frozen=True)
class KappaMoment:
    label: str
    n_regions: int
    sum_y: float
    sum_e: float
    rho: float
    sum_mu: float
    sum_mu2: float
    excess: float
    raw_alpha: float | None
    clipped_alpha: float
    fallback_used: bool
    fallback_reason: str
    beta_from_clipped_alpha: float | None
    beta_from_raw_alpha: float | None
    median_e: float
    q90_e: float
    q99_e: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", default=str(DEFAULT_BAM))
    parser.add_argument("--index", default=str(DEFAULT_INDEX))
    parser.add_argument("--quant-dir", default=str(DEFAULT_QUANT_DIR))
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--bgzf-threads", type=int, default=4)
    parser.add_argument("--splicing-anchor-tolerance", type=int, default=3)
    parser.add_argument("--top-transcripts", type=int, default=50_000)
    return parser.parse_args()


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float, fallback: float) -> float:
    values = np.asarray(values, dtype=np.float64).ravel()
    weights = np.asarray(weights, dtype=np.float64).ravel()
    keep = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not keep.any():
        return float(fallback)
    values = values[keep]
    weights = weights[keep]
    order = np.argsort(values, kind="stable")
    values = values[order]
    weights = weights[order]
    cumulative = np.cumsum(weights)
    target = float(q) * float(cumulative[-1])
    idx = int(np.searchsorted(cumulative, target, side="left"))
    idx = min(idx, values.size - 1)
    return float(values[idx])


def unweighted_quantiles(values: np.ndarray, quantiles: Iterable[float]) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {f"q{int(q * 1000):03d}": float("nan") for q in quantiles}
    return {f"q{int(q * 1000):03d}": float(np.quantile(values, q)) for q in quantiles}


def weighted_quantiles(values: np.ndarray, weights: np.ndarray, fallback: float = 0.0) -> dict[str, float]:
    qs = (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99, 0.999)
    return {f"q{int(q * 1000):03d}": weighted_quantile(values, weights, q, fallback) for q in qs}


def compute_e_y(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    gdna_fl,
    strand_summary: StrandSummary,
    splicing_anchor_tolerance: int,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], dict[str, object]]:
    spans = (region_arrays.end - region_arrays.start).astype(np.int64, copy=False)
    b_cross = boundary_crossing_exposure(
        gdna_fl,
        splicing_anchor_tolerance=splicing_anchor_tolerance,
    )

    type_arr = region_arrays.type
    is_intergenic = type_arr == int(RegionType.INTERGENIC)
    is_intron = type_arr == int(RegionType.INTRON)
    is_exon = type_arr == int(RegionType.EXON)

    E = np.zeros(region_arrays.start.size, dtype=np.float64)
    contained_mask = is_intergenic | is_intron
    E[contained_mask] = l_eff_contained(spans[contained_mask], gdna_fl)

    bf_left = region_arrays.bf_left.astype(bool, copy=False)
    bf_right = region_arrays.bf_right.astype(bool, copy=False)
    eligible_left = (is_exon & bf_left).astype(np.int64)
    eligible_right = (is_exon & bf_right).astype(np.int64)
    sides = eligible_left + eligible_right
    E[is_exon] = sides[is_exon].astype(np.float64) * float(b_cross)

    Y = np.zeros(region_arrays.start.size, dtype=np.float64)
    Y[is_intergenic] = payload_arrays.intergenic_per_region[is_intergenic].astype(
        np.float64, copy=False
    )

    strand_active = strand_correction_usable(strand_summary)
    signed_strand_contrast = float(strand_summary.signed_strand_contrast)
    identifiable = _strand_identifiable_rows(region_arrays.strand)

    intron_raw = payload_arrays.intron_by_orient.sum(axis=1).astype(np.float64, copy=False)
    if strand_active:
        same = payload_arrays.intron_by_orient[:, ORIENT_SAME]
        opp = payload_arrays.intron_by_orient[:, ORIENT_OPP]
        corrected = _gdna_count_moment(same, opp, signed_strand_contrast=signed_strand_contrast)
        intron_y = np.where(identifiable, np.maximum(corrected, 0.0), intron_raw)
    else:
        intron_y = intron_raw
    Y[is_intron] = intron_y[is_intron]

    ul = payload_arrays.u_left_by_orient
    ur = payload_arrays.u_right_by_orient
    exon_raw = (eligible_left * ul.sum(axis=1) + eligible_right * ur.sum(axis=1)).astype(
        np.float64, copy=False
    )
    if strand_active:
        same_ex = eligible_left * ul[:, ORIENT_SAME] + eligible_right * ur[:, ORIENT_SAME]
        opp_ex = eligible_left * ul[:, ORIENT_OPP] + eligible_right * ur[:, ORIENT_OPP]
        corrected_ex = _gdna_count_moment(
            same_ex, opp_ex, signed_strand_contrast=signed_strand_contrast
        )
        exon_y = np.where(identifiable, np.maximum(corrected_ex, 0.0), exon_raw)
    else:
        exon_y = exon_raw
    Y[is_exon] = exon_y[is_exon]

    class_masks = {
        CLASS_INTERGENIC: is_intergenic,
        CLASS_INTRON: is_intron,
        CLASS_EXON: is_exon,
    }
    diagnostics = {
        "b_cross": float(b_cross),
        "strand_active": bool(strand_active),
        "signed_strand_contrast": float(signed_strand_contrast),
        "eligible_exon_sides": {
            "zero": int(np.sum(is_exon & (sides == 0))),
            "one": int(np.sum(is_exon & (sides == 1))),
            "two": int(np.sum(is_exon & (sides == 2))),
        },
        "raw_counts": {
            "intron_raw": float(intron_raw[is_intron].sum()),
            "intron_corrected": float(intron_y[is_intron].sum()),
            "exon_boundary_raw": float(exon_raw[is_exon].sum()),
            "exon_boundary_corrected": float(exon_y[is_exon].sum()),
        },
    }
    return E, Y, class_masks, diagnostics


def moment_for(label: str, Y: np.ndarray, E: np.ndarray, mask: np.ndarray) -> KappaMoment:
    valid = mask & np.isfinite(Y) & np.isfinite(E) & (E > 0.0)
    if not valid.any():
        return KappaMoment(
            label=label,
            n_regions=0,
            sum_y=0.0,
            sum_e=0.0,
            rho=0.0,
            sum_mu=0.0,
            sum_mu2=0.0,
            excess=0.0,
            raw_alpha=None,
            clipped_alpha=float("nan"),
            fallback_used=True,
            fallback_reason="no valid regions",
            beta_from_clipped_alpha=None,
            beta_from_raw_alpha=None,
            median_e=float("nan"),
            q90_e=float("nan"),
            q99_e=float("nan"),
        )

    y = Y[valid].astype(np.float64, copy=False)
    e = E[valid].astype(np.float64, copy=False)
    sum_y = float(y.sum())
    sum_e = float(e.sum())
    rho = sum_y / sum_e if sum_e > 0.0 else 0.0
    mu = rho * e
    excess = float(np.sum((y - mu) ** 2) - np.sum(mu))
    sum_mu2 = float(np.sum(mu * mu))
    raw_alpha = sum_mu2 / excess if excess > 0.0 else None
    est = estimate_kappa(y, e, rho)
    clipped_alpha = float(est.value)
    beta_from_clipped = clipped_alpha / rho if rho > 0.0 else None
    beta_from_raw = raw_alpha / rho if raw_alpha is not None and rho > 0.0 else None
    return KappaMoment(
        label=label,
        n_regions=int(valid.sum()),
        sum_y=sum_y,
        sum_e=sum_e,
        rho=float(rho),
        sum_mu=float(mu.sum()),
        sum_mu2=sum_mu2,
        excess=excess,
        raw_alpha=None if raw_alpha is None else float(raw_alpha),
        clipped_alpha=clipped_alpha,
        fallback_used=bool(est.fallback_used),
        fallback_reason=str(est.fallback_reason),
        beta_from_clipped_alpha=beta_from_clipped,
        beta_from_raw_alpha=beta_from_raw,
        median_e=float(np.quantile(e, 0.50)),
        q90_e=float(np.quantile(e, 0.90)),
        q99_e=float(np.quantile(e, 0.99)),
    )


def rho_hat_for(Y: np.ndarray, E: np.ndarray, rho_global: float, kappa_len: float) -> np.ndarray:
    valid = E > 0.0
    rho_hat = np.full(Y.size, rho_global, dtype=np.float64)
    denom = E[valid] + float(kappa_len)
    rho_hat[valid] = np.where(denom > 0.0, (Y[valid] + kappa_len * rho_global) / denom, rho_global)
    return rho_hat


def make_exposure(
    region_arrays: RegionArrays,
    rho_hat: np.ndarray,
    rho_global: float,
    rho_ref: float,
    kappa_len: float,
) -> RegionalGdnaExposure:
    valid = np.isfinite(rho_hat) & (rho_hat > 0.0) & (rho_ref > 0.0)
    log_weight = np.zeros_like(rho_hat, dtype=np.float64)
    raw = np.full_like(rho_hat, LOG_A_FLOOR, dtype=np.float64)
    raw[valid] = np.log(rho_hat[valid]) - float(np.log(rho_ref))
    raw = np.minimum(raw, 0.0)
    log_weight[:] = np.maximum(raw, LOG_A_FLOOR)
    weight = np.exp(log_weight)
    return RegionalGdnaExposure(
        rho_hat=rho_hat.copy(),
        log_weight=log_weight,
        weight=weight,
        mode="regional",
        rho_ref=float(rho_ref),
        n_at_floor=int(np.sum(log_weight <= LOG_A_FLOOR + 1e-15)),
        per_class={},
        rho_global=float(rho_global),
        kappa_alpha_global=float(kappa_len * rho_global),
        kappa_opportunity_bp=float(kappa_len),
        ref_offsets=np.asarray(region_arrays.ref_offsets, dtype=np.int32).copy(),
        ref_id=np.asarray(region_arrays.ref_id, dtype=np.int32).copy(),
        start=np.asarray(region_arrays.start, dtype=np.int64).copy(),
        end=np.asarray(region_arrays.end, dtype=np.int64).copy(),
    )


def summarize_field(
    label: str,
    rho_hat: np.ndarray,
    weight: np.ndarray,
    E: np.ndarray,
    Y: np.ndarray,
    class_masks: dict[str, np.ndarray],
    rho_global: float,
    rho_ref: float,
    kappa_len: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    valid = E > 0.0
    for class_label in ("ALL", *CLASS_ORDER):
        mask = valid if class_label == "ALL" else valid & class_masks[class_label]
        if not mask.any():
            continue
        row = {
            "field": label,
            "class": class_label,
            "kappa_len": float(kappa_len),
            "rho_global": float(rho_global),
            "rho_ref": float(rho_ref),
            "n_regions": int(mask.sum()),
            "sum_y": float(Y[mask].sum()),
            "sum_e": float(E[mask].sum()),
            "mean_rho_hat_e_weighted": float(np.sum(rho_hat[mask] * E[mask]) / np.sum(E[mask])),
            "mean_weight_e_weighted": float(np.sum(weight[mask] * E[mask]) / np.sum(E[mask])),
            "n_at_floor": int(np.sum(weight[mask] <= np.exp(LOG_A_FLOOR) * (1.0 + 1e-12))),
            "n_at_one": int(np.sum(weight[mask] >= 1.0 - 1e-12)),
        }
        for prefix, values in (("rho", rho_hat[mask]), ("weight", weight[mask])):
            qs = weighted_quantiles(values, E[mask], fallback=0.0)
            for key, value in qs.items():
                row[f"{prefix}_{key}"] = value
        rows.append(row)
    return rows


def reference_table(
    rho_hat: np.ndarray,
    E: np.ndarray,
    class_masks: dict[str, np.ndarray],
    rho_global: float,
) -> dict[str, float]:
    valid = (E > 0.0) & np.isfinite(rho_hat)
    refs: dict[str, float] = {}
    for q in (0.95, 0.99, 0.999):
        refs[f"global_q{int(q * 1000):03d}"] = weighted_quantile(
            rho_hat[valid], E[valid], q, rho_global
        )
    refs["global_max"] = float(np.nanmax(rho_hat[valid]))
    for class_label in CLASS_ORDER:
        mask = valid & class_masks[class_label]
        if not mask.any():
            continue
        for q in (0.50, 0.95, 0.99, 0.999):
            refs[f"{class_label.lower()}_q{int(q * 1000):03d}"] = weighted_quantile(
                rho_hat[mask], E[mask], q, rho_global
            )
        refs[f"{class_label.lower()}_max"] = float(np.nanmax(rho_hat[mask]))
    return refs


def select_transcript_indices(index: TranscriptIndex, quant_dir: Path, top_n: int) -> np.ndarray:
    quant_path = quant_dir / "quant.feather"
    if not quant_path.exists():
        return np.arange(min(index.num_transcripts, top_n), dtype=np.int64)
    quant_df = pd.read_feather(quant_path, columns=["transcript_id", "count", "is_nrna"])
    quant_df = quant_df.sort_values("count", ascending=False, kind="stable").head(top_n)
    tid_to_index = pd.Series(np.arange(index.num_transcripts), index=index.t_df["t_id"])
    mapped = tid_to_index.reindex(quant_df["transcript_id"]).dropna().to_numpy(dtype=np.int64)
    return np.unique(mapped)


def subset_transcript_weights(
    index: TranscriptIndex,
    exposure: RegionalGdnaExposure,
    t_indices: np.ndarray,
) -> pd.DataFrame:
    is_nrna = index.t_df["is_nrna"].to_numpy(dtype=bool)
    starts = index.t_df["start"].to_numpy(dtype=np.int64)
    ends = index.t_df["end"].to_numpy(dtype=np.int64)
    ref_ids = np.asarray(index.t_to_ref_arr, dtype=np.int32)
    transcript_ids = index.t_df["t_id"].to_numpy()
    exon_offsets, exon_starts, exon_ends, _ = index.build_exon_csr()

    rows = []
    for t_idx in t_indices:
        t = int(t_idx)
        ref_id = int(ref_ids[t])
        if bool(is_nrna[t]):
            blocks = [(ref_id, int(starts[t]), int(ends[t]))]
        else:
            lo = int(exon_offsets[t])
            hi = int(exon_offsets[t + 1])
            if hi > lo:
                blocks = [
                    (ref_id, int(exon_starts[pos]), int(exon_ends[pos]))
                    for pos in range(lo, hi)
                ]
            else:
                blocks = [(ref_id, int(starts[t]), int(ends[t]))]
        rows.append(
            {
                "t_idx": t,
                "transcript_id": str(transcript_ids[t]),
                "is_nrna": bool(is_nrna[t]),
                "weight": footprint_exposure_weight(blocks, exposure),
            }
        )
    return pd.DataFrame(rows)


def summarize_transcript_weights(label: str, weights_df: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for class_label, mask in (
        ("mRNA", ~weights_df["is_nrna"].to_numpy(dtype=bool)),
        ("nRNA", weights_df["is_nrna"].to_numpy(dtype=bool)),
        ("ALL", np.ones(len(weights_df), dtype=bool)),
    ):
        sub = weights_df.loc[mask, "weight"].to_numpy(dtype=np.float64)
        if sub.size == 0:
            continue
        row = {
            "field": label,
            "component_class": class_label,
            "n_transcripts": int(sub.size),
            "mean": float(np.mean(sub)),
            "min": float(np.min(sub)),
            "max": float(np.max(sub)),
        }
        row.update(unweighted_quantiles(sub, (0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)))
        rows.append(row)
    return rows


def boundary_side_summary(
    region_arrays: RegionArrays,
    payload_arrays: PayloadArrays,
    E: np.ndarray,
    Y: np.ndarray,
    class_masks: dict[str, np.ndarray],
) -> list[dict[str, object]]:
    is_exon = class_masks[CLASS_EXON]
    sides = region_arrays.bf_left.astype(np.int64) + region_arrays.bf_right.astype(np.int64)
    groups = {
        "all_eligible": is_exon & (sides > 0),
        "internal_two_sides": is_exon & (sides == 2),
        "terminal_one_side": is_exon & (sides == 1),
        "ineligible_zero_sides": is_exon & (sides == 0),
        "left_only": is_exon & (sides == 1) & region_arrays.bf_left,
        "right_only": is_exon & (sides == 1) & region_arrays.bf_right,
    }
    rows = []
    for label, mask in groups.items():
        if not mask.any():
            continue
        sum_e = float(E[mask].sum())
        sum_y = float(Y[mask].sum())
        rows.append(
            {
                "boundary_group": label,
                "n_regions": int(mask.sum()),
                "sum_y": sum_y,
                "sum_e": sum_e,
                "rho": sum_y / sum_e if sum_e > 0.0 else float("nan"),
                "raw_u_left": float(payload_arrays.u_left[mask].sum()),
                "raw_u_right": float(payload_arrays.u_right[mask].sum()),
            }
        )
    return rows


def safe_float(obj):
    if isinstance(obj, dict):
        return {str(k): safe_float(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [safe_float(v) for v in obj]
    if isinstance(obj, tuple):
        return [safe_float(v) for v in obj]
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[vcap-capture] bam       : {args.bam}")
    print(f"[vcap-capture] index     : {args.index}")
    print(f"[vcap-capture] quant-dir : {args.quant_dir}")
    print(f"[vcap-capture] out-dir   : {out_dir}")

    index = TranscriptIndex.load(args.index)
    scan_cfg = BamScanConfig(
        total_threads=int(args.threads),
        bgzf_threads=int(args.bgzf_threads),
        splicing_anchor_tolerance=int(args.splicing_anchor_tolerance),
    )
    stats, strand_models, fl_models_scan, buffer, payload = scan_and_buffer(args.bam, index, scan_cfg)
    del buffer
    gc.collect()
    if payload is None:
        raise RuntimeError("Calibration payload is missing; rebuild the index with regions.feather.")

    strand_summary = StrandSummary.from_model(strand_models.exonic_spliced)
    calibration = calibrate(
        index=index,
        payload=payload,
        scan_trained=fl_models_scan,
        strand_summary=strand_summary,
        regional_exposure_enabled=True,
    )

    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(payload, region_arrays)
    E, Y, class_masks, evidence_diagnostics = compute_e_y(
        region_arrays,
        payload_arrays,
        calibration.fl_models.gdna,
        strand_summary,
        int(args.splicing_anchor_tolerance),
    )
    valid = E > 0.0
    rho_global = float(Y[valid].sum() / E[valid].sum())

    kappa_rows = [moment_for("ALL", Y, E, valid)]
    kappa_rows.extend(moment_for(label, Y, E, class_masks[label]) for label in CLASS_ORDER)
    kappa_df = pd.DataFrame([asdict(row) for row in kappa_rows])
    kappa_df.to_csv(out_dir / "kappa_moments.tsv", sep="\t", index=False)

    all_moment = kappa_rows[0]
    current_kappa_len = float(calibration.regional_exposure.kappa_opportunity_bp)
    beta_global_len = (
        float(all_moment.beta_from_clipped_alpha)
        if all_moment.beta_from_clipped_alpha is not None
        else current_kappa_len
    )
    beta_raw_global_len = (
        float(all_moment.beta_from_raw_alpha)
        if all_moment.beta_from_raw_alpha is not None
        else beta_global_len
    )
    median_e_len = float(all_moment.median_e) if np.isfinite(all_moment.median_e) else beta_global_len
    one_frag_len = 1.0 / rho_global if rho_global > 0.0 else beta_global_len
    five_frag_len = 5.0 / rho_global if rho_global > 0.0 else beta_global_len

    kappa_candidates = {
        "production_opportunity_length": current_kappa_len,
        "raw_alpha_over_rho_global": beta_raw_global_len,
        "clipped_alpha_over_rho_global": beta_global_len,
        "one_prior_fragment_over_rho": one_frag_len,
        "five_prior_fragments_over_rho": five_frag_len,
        "median_valid_E": median_e_len,
    }

    reference_rows = []
    field_rows = []
    transcript_rows = []
    t_indices = select_transcript_indices(index, Path(args.quant_dir), int(args.top_transcripts))

    field_specs: list[tuple[str, str]] = [
        ("production_opportunity_length", "global_q950"),
        ("production_opportunity_length", "global_q990"),
        ("production_opportunity_length", "exon-intron_q990"),
        ("clipped_alpha_over_rho_global", "global_q950"),
        ("clipped_alpha_over_rho_global", "exon-intron_q950"),
        ("clipped_alpha_over_rho_global", "exon-intron_q990"),
        ("clipped_alpha_over_rho_global", "exon-intron_max"),
        ("five_prior_fragments_over_rho", "exon-intron_q990"),
    ]

    refs_by_kappa: dict[str, dict[str, float]] = {}
    for kappa_label, kappa_len in kappa_candidates.items():
        rho_hat = rho_hat_for(Y, E, rho_global, float(kappa_len))
        refs = reference_table(rho_hat, E, class_masks, rho_global)
        refs_by_kappa[kappa_label] = refs
        for ref_label, ref_value in refs.items():
            reference_rows.append(
                {
                    "kappa_label": kappa_label,
                    "kappa_len": float(kappa_len),
                    "reference_label": ref_label,
                    "rho_ref": float(ref_value),
                    "rho_ref_over_current": float(ref_value / calibration.regional_exposure.rho_ref),
                    "rho_ref_over_rho_global": float(ref_value / rho_global),
                }
            )

    for kappa_label, ref_label in field_specs:
        kappa_len = float(kappa_candidates[kappa_label])
        rho_hat = rho_hat_for(Y, E, rho_global, kappa_len)
        refs = refs_by_kappa[kappa_label]
        if ref_label not in refs:
            continue
        rho_ref = refs[ref_label]
        exposure = make_exposure(region_arrays, rho_hat, rho_global, rho_ref, kappa_len)
        label = f"{kappa_label}__{ref_label}"
        field_rows.extend(
            summarize_field(
                label,
                rho_hat,
                exposure.weight,
                E,
                Y,
                class_masks,
                rho_global,
                rho_ref,
                kappa_len,
            )
        )
        weights_df = subset_transcript_weights(index, exposure, t_indices)
        weights_df["field"] = label
        transcript_rows.extend(summarize_transcript_weights(label, weights_df))

    references_df = pd.DataFrame(reference_rows)
    references_df.to_csv(out_dir / "reference_candidates.tsv", sep="\t", index=False)
    pd.DataFrame(field_rows).to_csv(out_dir / "field_candidate_summaries.tsv", sep="\t", index=False)
    pd.DataFrame(transcript_rows).to_csv(
        out_dir / "top_transcript_weight_summaries.tsv", sep="\t", index=False
    )
    pd.DataFrame(boundary_side_summary(region_arrays, payload_arrays, E, Y, class_masks)).to_csv(
        out_dir / "exon_boundary_side_summary.tsv", sep="\t", index=False
    )

    zero_shrink_rows = []
    e_probe_values = sorted(
        {
            1.0,
            float(np.quantile(E[valid], 0.10)),
            float(np.quantile(E[valid], 0.50)),
            float(np.quantile(E[valid], 0.90)),
            float(np.quantile(E[valid], 0.99)),
        }
    )
    for kappa_label, kappa_len in kappa_candidates.items():
        for e_value in e_probe_values:
            zero_shrink_rows.append(
                {
                    "kappa_label": kappa_label,
                    "kappa_len": float(kappa_len),
                    "E": float(e_value),
                    "rho_hat_zero_over_rho_global": float(kappa_len / (e_value + kappa_len)),
                }
            )
    pd.DataFrame(zero_shrink_rows).to_csv(out_dir / "zero_count_shrinkage.tsv", sep="\t", index=False)

    summary = {
        "inputs": {
            "bam": str(args.bam),
            "index": str(args.index),
            "quant_dir": str(args.quant_dir),
            "splicing_anchor_tolerance": int(args.splicing_anchor_tolerance),
        },
        "scan": {
            "n_fragments": int(stats.n_fragments),
            "n_read_names": int(stats.n_read_names),
            "payload_n_observed": int(payload.n_observed),
        },
        "strand_summary": {
            "n_observations": int(strand_summary.n_observations),
            "p_r1_sense": float(strand_summary.p_r1_sense),
            "signed_strand_contrast": float(strand_summary.signed_strand_contrast),
            "margin_99": float(strand_summary.signed_strand_contrast_margin(confidence=0.99)),
            "usable": bool(strand_correction_usable(strand_summary)),
        },
        "current_regional_exposure": calibration.regional_exposure.to_summary_dict(),
        "global_densities": calibration.global_densities.to_summary_dict(),
        "evidence_diagnostics": evidence_diagnostics,
        "rho_global_recomputed": float(rho_global),
        "kappa_candidates": kappa_candidates,
        "outputs": {
            "kappa_moments": str(out_dir / "kappa_moments.tsv"),
            "reference_candidates": str(out_dir / "reference_candidates.tsv"),
            "field_candidate_summaries": str(out_dir / "field_candidate_summaries.tsv"),
            "top_transcript_weight_summaries": str(out_dir / "top_transcript_weight_summaries.tsv"),
            "zero_count_shrinkage": str(out_dir / "zero_count_shrinkage.tsv"),
            "exon_boundary_side_summary": str(out_dir / "exon_boundary_side_summary.tsv"),
        },
    }
    with (out_dir / "summary.json").open("w") as handle:
        json.dump(safe_float(summary), handle, indent=2, sort_keys=True)

    print("[vcap-capture] wrote diagnostics:")
    for path in summary["outputs"].values():
        print(f"  {path}")
    print(f"  {out_dir / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())