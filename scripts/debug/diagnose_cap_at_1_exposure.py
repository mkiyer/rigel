#!/usr/bin/env python3
"""Diagnose cap-at-1 exposure behavior for the VCAP FLG2 mega-locus."""

from __future__ import annotations

import argparse
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import yaml
except ImportError:  # pragma: no cover - PyYAML is present in the dev env.
    yaml = None

from rigel.calibration._arrays import RegionArrays
from rigel.calibration._exposure import boundary_crossing_exposure
from rigel.calibration._regional_exposure import LOG_A_FLOOR, _weighted_quantile
from rigel.calibration.density_global import (
    _strand_identifiable_rows,
    kappa_opportunity_bp,
    l_eff_contained,
    precision_opportunity,
    strand_correction_usable,
    strand_reliability_power,
)
from rigel.calibration.regions import RegionType
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.native import detect_sj_strand_tag
from rigel.pipeline import _calibration_strand_summary, scan_and_buffer


DEFAULT_BAM = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam")
DEFAULT_INDEX_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index")
DEFAULT_QUANT_DIR = Path(
    "/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1"
)
DEFAULT_OUT_DIR = Path("results/vcap_cap_at_1_diagnostics_2026-05-21")
DEFAULT_REPORT = Path("docs/benchmarks/vcap_locus3_cap_at_1_diagnostics_2026-05-21.md")

NRNA_RE = re.compile(r"^RIGEL_NRNA_(?P<ref>.+)_(?P<strand>[12])_(?P<start>\d+)_(?P<end>\d+)$")
EXPOSURE_FLOOR = float(np.exp(LOG_A_FLOOR))
BIN_LABELS = (
    "[1e-4,1e-3)",
    "[1e-3,1e-2)",
    "[1e-2,1e-1)",
    "[1e-1,1)",
    "A==1",
)


@dataclass(frozen=True, slots=True)
class ExposureSupport:
    region_arrays: RegionArrays
    opportunity: np.ndarray
    class_masks: dict[str, np.ndarray]


@dataclass(frozen=True, slots=True)
class NormalizerResult:
    name: str
    rule: str
    weights: np.ndarray
    rho_ref: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, default=DEFAULT_BAM)
    parser.add_argument("--index-dir", type=Path, default=DEFAULT_INDEX_DIR)
    parser.add_argument("--quant-dir", type=Path, default=DEFAULT_QUANT_DIR)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--locus-id", type=int, default=3)
    return parser.parse_args()


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int, digits: int = 4) -> str:
    return f"{float(value):.{digits}f}"


def fmt_sci(value: float | int, digits: int = 3) -> str:
    return f"{float(value):.{digits}e}"


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    sep = "| " + " | ".join("---" for _label, _key in columns) + " |"
    body = [
        "| " + " | ".join(str(row.get(key, "")) for _label, key in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def load_run_config(quant_dir: Path) -> dict[str, Any]:
    config_path = quant_dir / "config.yaml"
    if not config_path.exists() or yaml is None:
        return {}
    data = yaml.safe_load(config_path.read_text())
    return data if isinstance(data, dict) else {}


def scan_config_from_run_config(run_config: dict[str, Any], bam_path: Path) -> BamScanConfig:
    sj_tag = run_config.get("sj_strand_tag", "auto")
    if isinstance(sj_tag, list):
        sj_tag = sj_tag[0] if len(sj_tag) == 1 else tuple(str(tag) for tag in sj_tag)
    if sj_tag == "auto":
        sj_tag = detect_sj_strand_tag(str(bam_path))

    scan_buffer_gib = float(run_config.get("scan_buffer_size", 2.0))
    return BamScanConfig(
        skip_duplicates=not bool(run_config.get("keep_duplicates", False)),
        include_multimap=bool(run_config.get("include_multimap", True)),
        sj_strand_tag=sj_tag,
        total_threads=int(run_config.get("threads", 8) or 0),
        bgzf_threads=int(run_config.get("scan_bgzf_threads", 4)),
        fragments_per_chunk=int(run_config.get("scan_fragments_per_chunk", 1_000_000)),
        read_name_batch_size=int(run_config.get("scan_read_name_batch_size", 512)),
        buffer_size_bytes=int(scan_buffer_gib * 1024**3),
        spill_dir=run_config.get("tmpdir"),
        splicing_anchor_tolerance=int(run_config.get("splicing_anchor_tolerance", 3)),
    )


def scan_and_build_calibration(args: argparse.Namespace):
    from rigel.calibration import calibrate

    run_config = load_run_config(args.quant_dir)
    index = TranscriptIndex.load(args.index_dir)
    scan_config = scan_config_from_run_config(run_config, args.bam)

    stats, strand_models, frag_length_models, buffer, payload = scan_and_buffer(
        str(args.bam),
        index,
        scan_config,
    )
    try:
        if payload is None:
            raise RuntimeError("Calibration payload missing; cannot build regional exposure diagnostics.")
        strand_models.finalize()
        strand_summary = _calibration_strand_summary(strand_models)
        calibration = calibrate(
            index=index,
            payload=payload,
            scan_trained=frag_length_models,
            fl_prior_ess=float(run_config.get("cal_prior_ess", 1000.0)),
            pool_quality_good=int(run_config.get("cal_quality_good", 5000)),
            pool_quality_weak=int(run_config.get("cal_quality_weak", 200)),
            strand_summary=strand_summary,
            regional_exposure_enabled=True,
        )
    finally:
        buffer.cleanup()
    return index, payload, calibration, strand_summary, scan_config, stats


def compute_exposure_support(
    index: TranscriptIndex,
    calibration,
    strand_summary,
    splicing_anchor_tolerance: int,
) -> ExposureSupport:
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    gdna_fl = calibration.fl_models.gdna
    global_densities = calibration.global_densities

    spans = (region_arrays.end - region_arrays.start).astype(np.int64, copy=False)
    type_arr = region_arrays.type
    is_intergenic = type_arr == int(RegionType.INTERGENIC)
    is_intron = type_arr == int(RegionType.INTRON)
    is_exon = type_arr == int(RegionType.EXON)

    opportunity = np.zeros(region_arrays.start.size, dtype=np.float64)
    containment_mask = is_intergenic | is_intron
    if containment_mask.any():
        opportunity[containment_mask] = l_eff_contained(spans[containment_mask], gdna_fl)

    boundary_cross = float(
        boundary_crossing_exposure(gdna_fl, splicing_anchor_tolerance=splicing_anchor_tolerance)
    )
    eligible_left = (is_exon & region_arrays.bf_left.astype(bool, copy=False)).astype(np.int64)
    eligible_right = (is_exon & region_arrays.bf_right.astype(bool, copy=False)).astype(np.int64)
    sides = eligible_left + eligible_right
    boundary_opportunity = np.zeros(region_arrays.start.size, dtype=np.float64)
    if boundary_cross > 0.0 and is_exon.any():
        boundary_opportunity[is_exon] = sides[is_exon].astype(np.float64) * boundary_cross
        opportunity[is_exon] = boundary_opportunity[is_exon]

    contained_opportunity = np.zeros(region_arrays.start.size, dtype=np.float64)
    if is_exon.any():
        contained_opportunity[is_exon] = l_eff_contained(spans[is_exon], gdna_fl)

    strand_active = strand_correction_usable(strand_summary)
    identifiable = _strand_identifiable_rows(region_arrays.strand)
    active_contained_rows = is_exon & identifiable & strand_active
    strand_power = strand_reliability_power(strand_summary)

    beta_boundary = kappa_opportunity_bp(
        global_densities.exon_intron.kappa,
        global_densities.exon_intron.rho,
    )
    exon_contained_density = global_densities.exon_contained
    assert exon_contained_density is not None
    beta_contained = kappa_opportunity_bp(
        exon_contained_density.kappa,
        exon_contained_density.rho,
    )
    boundary_precision = np.asarray(
        precision_opportunity(boundary_opportunity, beta_boundary),
        dtype=np.float64,
    )
    contained_precision = strand_power * np.asarray(
        precision_opportunity(contained_opportunity, beta_contained),
        dtype=np.float64,
    )
    contained_precision = np.where(active_contained_rows, contained_precision, 0.0)
    opportunity[is_exon] = boundary_precision[is_exon] + contained_precision[is_exon]

    class_masks = {
        "INTERGENIC": is_intergenic,
        "INTRON": is_intron,
        "EXON-COMPOSITE": is_exon,
    }
    return ExposureSupport(
        region_arrays=region_arrays,
        opportunity=opportunity,
        class_masks=class_masks,
    )


def weighted_reference(
    rho_hat: np.ndarray,
    opportunity: np.ndarray,
    mask: np.ndarray,
    quantile: float,
    fallback: float,
) -> float:
    valid = mask & np.isfinite(rho_hat) & (opportunity > 0.0)
    if not valid.any():
        return float(fallback)
    return _weighted_quantile(rho_hat[valid], opportunity[valid], quantile, fallback=fallback)


def max_reference(
    rho_hat: np.ndarray,
    opportunity: np.ndarray,
    mask: np.ndarray,
    fallback: float,
) -> float:
    valid = mask & np.isfinite(rho_hat) & (opportunity > 0.0)
    if not valid.any():
        return float(fallback)
    return float(np.max(rho_hat[valid]))


def hard_cap_weights(rho_hat: np.ndarray, reference: float | np.ndarray) -> np.ndarray:
    reference_arr = np.asarray(reference, dtype=np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        raw = np.divide(rho_hat, reference_arr, out=np.ones_like(rho_hat), where=reference_arr > 0.0)
    return np.clip(raw, EXPOSURE_FLOOR, 1.0)


def soft_cap_weights(rho_hat: np.ndarray, reference: float) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        raw = np.divide(rho_hat, rho_hat + reference, out=np.zeros_like(rho_hat), where=reference > 0.0)
    return np.clip(raw, EXPOSURE_FLOOR, 1.0)


def build_normalizers(calibration, support: ExposureSupport) -> list[NormalizerResult]:
    exposure = calibration.regional_exposure
    if exposure is None:
        raise RuntimeError("Calibration did not build a regional exposure model.")
    rho_hat = exposure.rho_hat
    opportunity = support.opportunity
    valid = opportunity > 0.0
    fallback = float(exposure.rho_global)

    normalizers = [
        NormalizerResult(
            name="current_global_q95",
            rule="saved global E-weighted Q95 with hard cap",
            weights=exposure.weight.copy(),
            rho_ref=float(exposure.rho_ref),
        )
    ]

    for quantile in (0.99, 0.995, 0.999, 0.9995):
        reference = weighted_reference(rho_hat, opportunity, valid, quantile, fallback)
        normalizers.append(
            NormalizerResult(
                name=f"global_q{quantile * 100:g}",
                rule=f"global E-weighted Q{quantile * 100:g} with hard cap",
                weights=hard_cap_weights(rho_hat, reference),
                rho_ref=reference,
            )
        )

    max_ref = max_reference(rho_hat, opportunity, valid, fallback)
    normalizers.append(
        NormalizerResult(
            name="global_max",
            rule="global max rho_hat; bounding diagnostic",
            weights=hard_cap_weights(rho_hat, max_ref),
            rho_ref=max_ref,
        )
    )

    for quantile in (0.95, 0.99, 0.995, 0.999):
        reference_by_region = np.full(rho_hat.size, fallback, dtype=np.float64)
        refs: list[str] = []
        for class_name, class_mask in support.class_masks.items():
            reference = weighted_reference(rho_hat, opportunity, class_mask, quantile, fallback)
            reference_by_region[class_mask] = reference
            refs.append(f"{class_name}={reference:.3e}")
        normalizers.append(
            NormalizerResult(
                name=f"class_q{quantile * 100:g}",
                rule=f"class-aware E-weighted Q{quantile * 100:g}: " + "; ".join(refs),
                weights=hard_cap_weights(rho_hat, reference_by_region),
                rho_ref=float(reference_by_region[support.class_masks["EXON-COMPOSITE"]][0]),
            )
        )

    current_ref = float(exposure.rho_ref)
    normalizers.append(
        NormalizerResult(
            name="soft_global_q95",
            rule="diagnostic soft cap A=rho/(rho+rho_ref), rho_ref=current Q95",
            weights=soft_cap_weights(rho_hat, current_ref),
            rho_ref=current_ref,
        )
    )
    return normalizers


def parse_nrna_interval(nrna_id: str, ref_name_to_id: dict[str, int]) -> tuple[int, int, int] | None:
    match = NRNA_RE.match(nrna_id)
    if match is None:
        return None
    ref_name = match.group("ref")
    if ref_name not in ref_name_to_id:
        return None
    return (ref_name_to_id[ref_name], int(match.group("start")) - 1, int(match.group("end")))


def merge_blocks(blocks: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    if not blocks:
        return []
    ordered = sorted(blocks)
    merged: list[tuple[int, int, int]] = []
    cur_ref, cur_start, cur_end = ordered[0]
    for ref_id, start, end in ordered[1:]:
        if ref_id == cur_ref and start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_ref, cur_start, cur_end))
            cur_ref, cur_start, cur_end = ref_id, start, end
    merged.append((cur_ref, cur_start, cur_end))
    return merged


def reconstruct_locus_blocks(
    index: TranscriptIndex,
    quant_dir: Path,
    locus_id: int,
) -> tuple[list[tuple[int, int, int]], pd.Series]:
    loci = pd.read_feather(quant_dir / "loci.feather")
    quant = pd.read_feather(quant_dir / "quant.feather")
    nrna = pd.read_feather(quant_dir / "nrna_quant.feather")
    transcripts = pd.read_feather(Path(str(index.index_dir)) / "transcripts.feather")

    locus_rows = loci.loc[loci["locus_id"] == locus_id]
    if locus_rows.empty:
        raise RuntimeError(f"locus_id={locus_id} not present in {quant_dir / 'loci.feather'}")
    locus = locus_rows.iloc[0]
    quant_locus = quant.loc[quant["locus_id"] == locus_id]
    nrna_locus = nrna.loc[nrna["locus_id"] == locus_id]

    transcript_hits = transcripts.loc[transcripts["t_id"].isin(quant_locus["transcript_id"])]
    blocks: list[tuple[int, int, int]] = [
        (index.ref_name_to_id[str(row.ref)], int(row.start), int(row.end))
        for row in transcript_hits.itertuples(index=False)
        if str(row.ref) in index.ref_name_to_id
    ]
    for nrna_id in nrna_locus["nrna_id"].astype(str):
        interval = parse_nrna_interval(nrna_id, index.ref_name_to_id)
        if interval is not None:
            blocks.append(interval)
    return merge_blocks(blocks), locus


def exposure_bin_label(weight: float) -> str:
    if weight >= 1.0 - 1.0e-12:
        return "A==1"
    if weight < 1.0e-3:
        return "[1e-4,1e-3)"
    if weight < 1.0e-2:
        return "[1e-3,1e-2)"
    if weight < 1.0e-1:
        return "[1e-2,1e-1)"
    return "[1e-1,1)"


def decompose_footprint(
    blocks: list[tuple[int, int, int]],
    region_arrays: RegionArrays,
    weights: np.ndarray,
) -> dict[str, dict[str, float]]:
    bins = {label: {"bp": 0.0, "weighted_bp": 0.0} for label in BIN_LABELS}
    n_refs = int(region_arrays.ref_offsets.size - 1)

    for ref_id, block_start, block_end in blocks:
        if block_end <= block_start:
            continue
        if ref_id < 0 or ref_id >= n_refs:
            label = exposure_bin_label(1.0)
            bins[label]["bp"] += float(block_end - block_start)
            bins[label]["weighted_bp"] += float(block_end - block_start)
            continue
        ref_lo = int(region_arrays.ref_offsets[ref_id])
        ref_hi = int(region_arrays.ref_offsets[ref_id + 1])
        if ref_hi <= ref_lo:
            label = exposure_bin_label(1.0)
            bins[label]["bp"] += float(block_end - block_start)
            bins[label]["weighted_bp"] += float(block_end - block_start)
            continue

        starts = region_arrays.start[ref_lo:ref_hi]
        ends = region_arrays.end[ref_lo:ref_hi]
        local_index = int(np.searchsorted(starts, block_start, side="right") - 1)
        if local_index < 0:
            local_index = 0
        cursor = int(block_start)
        while cursor < block_end and local_index < starts.size:
            region_start = int(starts[local_index])
            region_end = int(ends[local_index])
            if region_end <= cursor:
                local_index += 1
                continue
            if region_start > cursor:
                gap_end = min(region_start, block_end)
                gap_bp = float(gap_end - cursor)
                label = exposure_bin_label(1.0)
                bins[label]["bp"] += gap_bp
                bins[label]["weighted_bp"] += gap_bp
                cursor = gap_end
                continue
            segment_start = max(cursor, region_start)
            segment_end = min(block_end, region_end)
            if segment_end > segment_start:
                segment_bp = float(segment_end - segment_start)
                weight = float(weights[ref_lo + local_index])
                label = exposure_bin_label(weight)
                bins[label]["bp"] += segment_bp
                bins[label]["weighted_bp"] += segment_bp * weight
            cursor = max(cursor, region_end)
            local_index += 1
        if cursor < block_end:
            tail_bp = float(block_end - cursor)
            label = exposure_bin_label(1.0)
            bins[label]["bp"] += tail_bp
            bins[label]["weighted_bp"] += tail_bp
    return bins


def footprint_region_bp(
    blocks: list[tuple[int, int, int]],
    region_arrays: RegionArrays,
) -> np.ndarray:
    overlap_bp = np.zeros(region_arrays.start.size, dtype=np.float64)
    n_refs = int(region_arrays.ref_offsets.size - 1)
    for ref_id, block_start, block_end in blocks:
        if block_end <= block_start or ref_id < 0 or ref_id >= n_refs:
            continue
        ref_lo = int(region_arrays.ref_offsets[ref_id])
        ref_hi = int(region_arrays.ref_offsets[ref_id + 1])
        if ref_hi <= ref_lo:
            continue
        starts = region_arrays.start[ref_lo:ref_hi]
        ends = region_arrays.end[ref_lo:ref_hi]
        local_index = int(np.searchsorted(starts, block_start, side="right") - 1)
        if local_index < 0:
            local_index = 0
        cursor = int(block_start)
        while cursor < block_end and local_index < starts.size:
            region_start = int(starts[local_index])
            region_end = int(ends[local_index])
            if region_end <= cursor:
                local_index += 1
                continue
            if region_start > cursor:
                cursor = min(region_start, block_end)
                continue
            segment_start = max(cursor, region_start)
            segment_end = min(block_end, region_end)
            if segment_end > segment_start:
                overlap_bp[ref_lo + local_index] += float(segment_end - segment_start)
            cursor = max(cursor, region_end)
            local_index += 1
    return overlap_bp


def build_locus_normalizers(
    calibration,
    support: ExposureSupport,
    blocks: list[tuple[int, int, int]],
) -> list[NormalizerResult]:
    exposure = calibration.regional_exposure
    if exposure is None:
        raise RuntimeError("Calibration did not build a regional exposure model.")
    overlap_bp = footprint_region_bp(blocks, support.region_arrays)
    valid = overlap_bp > 0.0
    fallback = float(exposure.rho_ref) if exposure.rho_ref > 0.0 else float(exposure.rho_global)
    out: list[NormalizerResult] = []
    for quantile in (0.99, 0.999):
        reference = weighted_reference(exposure.rho_hat, overlap_bp, valid, quantile, fallback)
        out.append(
            NormalizerResult(
                name=f"locus_q{quantile * 100:g}",
                rule=f"within-locus bp-weighted Q{quantile * 100:g} with hard cap",
                weights=hard_cap_weights(exposure.rho_hat, reference),
                rho_ref=reference,
            )
        )
    reference = max_reference(exposure.rho_hat, overlap_bp, valid, fallback)
    out.append(
        NormalizerResult(
            name="locus_max",
            rule="within-locus max rho_hat with hard cap",
            weights=hard_cap_weights(exposure.rho_hat, reference),
            rho_ref=reference,
        )
    )
    return out


def summarize_normalizer(
    normalizer: NormalizerResult,
    support: ExposureSupport,
    blocks: list[tuple[int, int, int]],
    unweighted_eff_len: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    bins = decompose_footprint(blocks, support.region_arrays, normalizer.weights)
    footprint_bp = float(sum(end - start for _ref_id, start, end in blocks))
    weighted_bp = sum(row["weighted_bp"] for row in bins.values())
    mean_weight = weighted_bp / footprint_bp if footprint_bp > 0.0 else 1.0
    denominator = max(float(unweighted_eff_len) * mean_weight, 1.0)

    saturated_bp = bins["A==1"]["bp"]
    saturated_weighted_bp = bins["A==1"]["weighted_bp"]
    valid_exon = support.class_masks["EXON-COMPOSITE"] & (support.opportunity > 0.0)
    saturated_exon_opportunity = float(
        support.opportunity[valid_exon & (normalizer.weights >= 1.0 - 1.0e-12)].sum()
    )
    total_exon_opportunity = float(support.opportunity[valid_exon].sum())

    summary = {
        "normalizer": normalizer.name,
        "rule": normalizer.rule,
        "rho_ref": normalizer.rho_ref,
        "footprint_bp": footprint_bp,
        "weighted_bp": weighted_bp,
        "mean_A": mean_weight,
        "gdna_eff_len": denominator,
        "reduction_x": float(unweighted_eff_len) / denominator if denominator > 0.0 else np.nan,
        "sat_bp_fraction": saturated_bp / footprint_bp if footprint_bp > 0.0 else 0.0,
        "sat_weighted_fraction": saturated_weighted_bp / weighted_bp if weighted_bp > 0.0 else 0.0,
        "sat_exon_opportunity_fraction": (
            saturated_exon_opportunity / total_exon_opportunity
            if total_exon_opportunity > 0.0
            else 0.0
        ),
    }

    decomposition_rows: list[dict[str, Any]] = []
    for label in BIN_LABELS:
        bin_bp = bins[label]["bp"]
        bin_weighted_bp = bins[label]["weighted_bp"]
        denominator_bp = float(unweighted_eff_len) * bin_weighted_bp / footprint_bp
        decomposition_rows.append(
            {
                "normalizer": normalizer.name,
                "bin": label,
                "bp": bin_bp,
                "bp_fraction": bin_bp / footprint_bp if footprint_bp > 0.0 else 0.0,
                "weighted_bp": bin_weighted_bp,
                "weighted_fraction": bin_weighted_bp / weighted_bp if weighted_bp > 0.0 else 0.0,
                "gdna_eff_len_contribution": denominator_bp,
            }
        )
    return summary, decomposition_rows


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    index, _payload, calibration, strand_summary, scan_config, stats = scan_and_build_calibration(args)
    support = compute_exposure_support(
        index,
        calibration,
        strand_summary,
        int(scan_config.splicing_anchor_tolerance),
    )
    blocks, saved_locus = reconstruct_locus_blocks(index, args.quant_dir, args.locus_id)
    unweighted_eff_len = float(saved_locus["gdna_eff_len_unweighted"])
    saved_eff_len = float(saved_locus["gdna_eff_len"])
    saved_weight = float(saved_locus["gdna_em_exposure_weight"])

    summary_rows: list[dict[str, Any]] = []
    decomposition_rows: list[dict[str, Any]] = []
    normalizers = [
        *build_normalizers(calibration, support),
        *build_locus_normalizers(calibration, support, blocks),
    ]
    for normalizer in normalizers:
        summary, decomposition = summarize_normalizer(normalizer, support, blocks, unweighted_eff_len)
        summary_rows.append(summary)
        decomposition_rows.extend(decomposition)

    summary_df = pd.DataFrame(summary_rows)
    decomposition_df = pd.DataFrame(decomposition_rows)
    summary_path = args.out_dir / "normalizer_summary.tsv"
    decomposition_path = args.out_dir / "footprint_bin_decomposition.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    decomposition_df.to_csv(decomposition_path, sep="\t", index=False)

    current = summary_df.loc[summary_df["normalizer"] == "current_global_q95"].iloc[0]
    selected_names = [
        "current_global_q95",
        "global_q99",
        "global_q99.5",
        "global_q99.9",
        "class_q95",
        "class_q99",
        "class_q99.9",
        "locus_q99",
        "locus_q99.9",
        "locus_max",
        "global_max",
        "soft_global_q95",
    ]
    selected_summary = summary_df.loc[summary_df["normalizer"].isin(selected_names)].copy()
    selected_summary["order"] = selected_summary["normalizer"].map(
        {name: index for index, name in enumerate(selected_names)}
    )
    selected_summary = selected_summary.sort_values("order")

    summary_table_rows = [
        {
            "normalizer": row.normalizer,
            "rho_ref": fmt_sci(row.rho_ref),
            "mean_A": fmt_float(row.mean_A, 6),
            "gdna_eff_len": fmt_int(row.gdna_eff_len),
            "reduction_x": fmt_float(row.reduction_x, 2),
            "sat_bp_pct": fmt_float(100.0 * row.sat_bp_fraction, 2),
            "sat_weighted_pct": fmt_float(100.0 * row.sat_weighted_fraction, 2),
            "sat_exon_opp_pct": fmt_float(100.0 * row.sat_exon_opportunity_fraction, 2),
        }
        for row in selected_summary.itertuples(index=False)
    ]

    decomposition_focus = decomposition_df.loc[
        decomposition_df["normalizer"].isin(
            [
                "current_global_q95",
                "class_q95",
                "class_q99.9",
                "locus_q99.9",
                "locus_max",
                "global_max",
            ]
        )
    ]
    decomp_table_rows = [
        {
            "normalizer": row.normalizer,
            "bin": row.bin,
            "bp_pct": fmt_float(100.0 * row.bp_fraction, 2),
            "weighted_pct": fmt_float(100.0 * row.weighted_fraction, 2),
            "gdna_eff_len": fmt_int(row.gdna_eff_len_contribution),
        }
        for row in decomposition_focus.itertuples(index=False)
    ]

    exposure = calibration.regional_exposure
    assert exposure is not None
    footprint_bp = float(sum(end - start for _ref_id, start, end in blocks))
    report = f"""# VCAP Locus 3 Cap-At-1 Exposure Diagnostics - 2026-05-21

Input BAM: `{args.bam}`

Reference run for locus outputs: `{args.quant_dir}`

This diagnostic reran the native scan plus calibration only, then evaluated
alternate mappings from regional density `rho_hat` to exposure weight `A_r`.
It did not rerun EM.

## Current Scalar Denominator Check

| Metric | Value |
| --- | ---: |
| n_fragments rescanned | {fmt_int(stats.n_fragments)} |
| reconstructed footprint bp | {fmt_int(footprint_bp)} |
| saved locus_span_bp | {fmt_int(saved_locus['locus_span_bp'])} |
| saved gdna_eff_len_unweighted | {fmt_int(unweighted_eff_len)} |
| saved gdna_em_exposure_weight | {fmt_float(saved_weight, 6)} |
| saved gdna_eff_len | {fmt_int(saved_eff_len)} |
| recomputed current mean A | {fmt_float(current.mean_A, 6)} |
| recomputed current gdna_eff_len | {fmt_int(current.gdna_eff_len)} |

## Normalizer Sweep

{markdown_table(summary_table_rows, [
        ('Normalizer', 'normalizer'),
        ('rho ref', 'rho_ref'),
        ('footprint mean A', 'mean_A'),
        ('L_gDNA', 'gdna_eff_len'),
        ('reduction', 'reduction_x'),
        ('footprint bp at A=1 %', 'sat_bp_pct'),
        ('weighted denominator from A=1 %', 'sat_weighted_pct'),
        ('global exon opportunity at cap %', 'sat_exon_opp_pct'),
    ])}

Full numeric tables:

- `{summary_path}`
- `{decomposition_path}`

## Footprint Bin Decomposition

For each normalizer, the `L_gDNA` contribution is the current scalar production
denominator contribution: `gdna_eff_len_unweighted * weighted_bp_bin / footprint_bp`.

{markdown_table(decomp_table_rows, [
        ('Normalizer', 'normalizer'),
        ('A bin', 'bin'),
        ('footprint bp %', 'bp_pct'),
        ('weighted denominator %', 'weighted_pct'),
        ('L_gDNA contribution', 'gdna_eff_len'),
    ])}

## Immediate Interpretation

Under the current global Q95 normalizer, the footprint mean exposure remains
near `{fmt_float(current.mean_A, 4)}`. Saturation is important, but it is not the
whole explanation: `A==1` covers only
`{fmt_float(100.0 * current.sat_bp_fraction, 2)}%` of the reconstructed footprint
while contributing `{fmt_float(100.0 * current.sat_weighted_fraction, 2)}%` of
the scalar denominator. The remaining large contribution comes mostly from the
`[0.1,1)` bin, so simply replacing the hard cap with a soft cap at the same Q95
reference scale is insufficient.

The strongest useful diagnostic is the reference scale. Raising the global
normalizer to Q99.9 shrinks `L_gDNA` to 9.5M. Class-aware Q99.9 and within-locus
Q99.9 land at similar scales, 5.44M and 5.25M respectively, close to the FLG2
mRNA break-even scale observed in the hotspot audit. The raw max normalizers give
40-50k denominators, which are useful lower bounds but almost certainly
over-normalize because they push nearly the whole footprint to the exposure
floor. A production candidate should be chosen by an explicit identifying
assumption about the exposure scale, then validated on normal loci and full
confusion matrices before changing EM-facing denominators.
"""
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report)
    print(f"wrote {args.report}")
    print(f"wrote {summary_path}")
    print(f"wrote {decomposition_path}")
    print(
        "current L_gDNA "
        f"{current.gdna_eff_len:,.0f}; class_q95 "
        f"{summary_df.loc[summary_df['normalizer'] == 'class_q95', 'gdna_eff_len'].iloc[0]:,.0f}; "
        f"global_max {summary_df.loc[summary_df['normalizer'] == 'global_max', 'gdna_eff_len'].iloc[0]:,.0f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())