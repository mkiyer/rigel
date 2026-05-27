"""Deep diagnostics for hybrid-capture calibration failures.

This script is intentionally file-based so it can be run safely from Copilot/VS Code
without heredoc stdin. It reruns a single condition in-process to expose the full
calibration arrays, then projects simulator oracle origins onto Rigel calibration
regions for per-region truth/error analysis.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

from rigel.annotate import AF_GDNA_BIT, AF_MRNA_BIT, AF_NRNA_BIT
from rigel.calibration._arrays import PayloadArrays, RegionArrays
from rigel.calibration.boundaries import build_boundary_table
from rigel.calibration.density_observation import build_density_observation
from rigel.calibration.latent_states import STATE_NAMES
from rigel.calibration.region_count_ledger import build_region_count_ledger
from rigel.config import BamScanConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline
from rigel.sim.truth import parse_origin
from rigel.types import STRAND_NEG


REGION_TYPE_LABELS = {
    0: "intergenic",
    1: "intron",
    2: "exon",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sim-base",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"),
        help="Simulation suite base directory.",
    )
    parser.add_argument(
        "--condition",
        default="gdna_high_ss_0.50_nrna_none_capture_on",
        help="Condition name to diagnose.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to <sim-base>/diagnostics/<condition>.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def resolve_suite_path(sim_base: Path, value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return sim_base / path


def find_condition(manifest: dict, name: str) -> dict:
    for condition in manifest.get("conditions", []):
        if condition.get("name") == name:
            return condition
    raise ValueError(f"Condition {name!r} not found in manifest.")


def region_type_label(value: object) -> str:
    try:
        return REGION_TYPE_LABELS.get(int(value), f"type_{int(value)}")
    except Exception:
        return str(value)


def transcript_blocks(
    origin,
    t_id_to_index: dict[str, int],
    t_index_to_exons: dict[int, np.ndarray],
    t_index_to_strand: dict[int, int],
) -> tuple[str | None, list[tuple[int, int]]]:
    if origin.transcript_id is None or origin.start is None or origin.end is None:
        return None, []
    t_index = t_id_to_index.get(origin.transcript_id)
    if t_index is None:
        return None, []
    exons = t_index_to_exons.get(t_index)
    if exons is None or len(exons) == 0:
        return None, []

    frag_start = int(origin.start)
    frag_end = int(origin.end)
    tx_len = int(np.sum(exons[:, 1] - exons[:, 0]))
    if int(t_index_to_strand[t_index]) == STRAND_NEG:
        frag_start, frag_end = tx_len - frag_end, tx_len - frag_start

    blocks: list[tuple[int, int]] = []
    consumed = 0
    for exon_start, exon_end in exons:
        exon_start = int(exon_start)
        exon_end = int(exon_end)
        exon_len = exon_end - exon_start
        exon_tx_start = consumed
        exon_tx_end = consumed + exon_len
        overlap_start = max(frag_start, exon_tx_start)
        overlap_end = min(frag_end, exon_tx_end)
        if overlap_start < overlap_end:
            offset_start = overlap_start - exon_tx_start
            offset_end = overlap_end - exon_tx_start
            blocks.append((exon_start + offset_start, exon_start + offset_end))
        consumed += exon_len
        if consumed >= frag_end:
            break

    t_row = t_index_to_row[t_index]
    return str(t_row["ref"]), blocks


def overlap_regions(
    ref_name: str | None,
    blocks: list[tuple[int, int]],
    ref_name_to_id: dict[str, int],
    region_arrays: RegionArrays,
) -> dict[int, int]:
    if ref_name is None or ref_name not in ref_name_to_id:
        return {}
    ref_id = int(ref_name_to_id[ref_name])
    ref_start = int(region_arrays.ref_offsets[ref_id])
    ref_end = int(region_arrays.ref_offsets[ref_id + 1])
    starts = region_arrays.start
    ends = region_arrays.end

    overlaps: dict[int, int] = defaultdict(int)
    local_ends = ends[ref_start:ref_end]
    for block_start, block_end in blocks:
        block_start = int(block_start)
        block_end = int(block_end)
        if block_end <= block_start:
            continue
        local_idx = int(np.searchsorted(local_ends, block_start, side="right"))
        region_idx = ref_start + local_idx
        while region_idx < ref_end and int(starts[region_idx]) < block_end:
            ov_start = max(block_start, int(starts[region_idx]))
            ov_end = min(block_end, int(ends[region_idx]))
            if ov_start < ov_end:
                overlaps[region_idx] += ov_end - ov_start
            region_idx += 1
    return dict(overlaps)


def origin_blocks(
    qname: str,
    t_id_to_index: dict[str, int],
    t_index_to_exons: dict[int, np.ndarray],
    t_index_to_strand: dict[int, int],
) -> tuple[str | None, str, list[tuple[int, int]]]:
    origin = parse_origin(qname)
    if origin.kind == "gdna":
        if origin.ref is None or origin.start is None or origin.end is None:
            return None, origin.kind, []
        return str(origin.ref), origin.kind, [(int(origin.start), int(origin.end))]
    ref_name, blocks = transcript_blocks(
        origin,
        t_id_to_index,
        t_index_to_exons,
        t_index_to_strand,
    )
    return ref_name, origin.kind, blocks


def assignment_from_zf(zf: int | None) -> str:
    if zf is None:
        return "missing"
    flags = int(zf)
    if flags & AF_GDNA_BIT:
        return "gdna"
    if flags & (AF_MRNA_BIT | AF_NRNA_BIT):
        return "rna"
    return "unresolved"


def add_truth_projection(
    annotated_bam: Path,
    region_df: pd.DataFrame,
    ref_name_to_id: dict[str, int],
    region_arrays: RegionArrays,
    t_id_to_index: dict[str, int],
    t_index_to_exons: dict[int, np.ndarray],
    t_index_to_strand: dict[int, int],
) -> tuple[pd.DataFrame, dict[str, object]]:
    n_regions = len(region_df)
    true_mass = {
        "gdna": np.zeros(n_regions, dtype=np.float64),
        "rna": np.zeros(n_regions, dtype=np.float64),
    }
    major_counts = {
        "gdna": np.zeros(n_regions, dtype=np.float64),
        "rna": np.zeros(n_regions, dtype=np.float64),
    }
    assignment_region_counts: dict[str, np.ndarray] = {
        "true_gdna_assigned_gdna_major": np.zeros(n_regions, dtype=np.float64),
        "true_gdna_assigned_rna_major": np.zeros(n_regions, dtype=np.float64),
        "true_gdna_assigned_unresolved_major": np.zeros(n_regions, dtype=np.float64),
        "true_rna_assigned_gdna_major": np.zeros(n_regions, dtype=np.float64),
        "true_rna_assigned_rna_major": np.zeros(n_regions, dtype=np.float64),
        "true_rna_assigned_unresolved_major": np.zeros(n_regions, dtype=np.float64),
    }
    confusion = Counter()
    seen: set[str] = set()
    skipped = Counter()

    with pysam.AlignmentFile(str(annotated_bam), "rb", check_sq=False) as bam:
        for read in bam:
            qname = read.query_name
            if qname in seen:
                continue
            seen.add(qname)
            try:
                ref_name, origin_kind, blocks = origin_blocks(
                    qname,
                    t_id_to_index,
                    t_index_to_exons,
                    t_index_to_strand,
                )
            except Exception:
                skipped["parse_error"] += 1
                continue
            truth_kind = "gdna" if origin_kind == "gdna" else "rna"
            overlaps = overlap_regions(ref_name, blocks, ref_name_to_id, region_arrays)
            if not overlaps:
                skipped[f"no_region_{truth_kind}"] += 1
                continue
            total_bp = float(sum(overlaps.values()))
            if total_bp <= 0.0:
                skipped[f"zero_bp_{truth_kind}"] += 1
                continue

            for region_idx, bp in overlaps.items():
                true_mass[truth_kind][region_idx] += float(bp) / total_bp
            major_region = max(overlaps.items(), key=lambda item: item[1])[0]
            major_counts[truth_kind][major_region] += 1.0

            try:
                zf = read.get_tag("ZF")
            except KeyError:
                zf = None
            assigned_kind = assignment_from_zf(zf)
            confusion[(truth_kind, assigned_kind)] += 1
            key = f"true_{truth_kind}_assigned_{assigned_kind}_major"
            if key in assignment_region_counts:
                assignment_region_counts[key][major_region] += 1.0

    out = pd.DataFrame(
        {
            "true_gdna_mass": true_mass["gdna"],
            "true_rna_mass": true_mass["rna"],
            "true_gdna_major_count": major_counts["gdna"],
            "true_rna_major_count": major_counts["rna"],
        }
    )
    for name, values in assignment_region_counts.items():
        out[name] = values
    out["true_total_mass"] = out["true_gdna_mass"] + out["true_rna_mass"]
    out["true_gdna_fraction"] = safe_div(out["true_gdna_mass"], out["true_total_mass"])

    summary = {
        "n_fragments_projected": int(len(seen) - sum(skipped.values())),
        "n_qnames_seen": int(len(seen)),
        "skipped": {str(key): int(value) for key, value in skipped.items()},
        "confusion": {f"{k[0]}->{k[1]}": int(v) for k, v in confusion.items()},
    }
    return out, summary


def safe_div(numer, denom):
    n = np.asarray(numer, dtype=np.float64)
    d = np.asarray(denom, dtype=np.float64)
    return np.divide(n, d, out=np.zeros_like(n, dtype=np.float64), where=d != 0.0)


def load_bed12_probe_blocks(path: Path | None) -> list[tuple[str, int, int]]:
    if path is None or not path.exists():
        return []
    blocks: list[tuple[str, int, int]] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            chrom_start = int(fields[1])
            chrom_end = int(fields[2])
            if len(fields) >= 12:
                block_count = int(fields[9])
                sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
                for size, rel_start in zip(sizes[:block_count], starts[:block_count], strict=False):
                    start = chrom_start + rel_start
                    end = start + size
                    if end > start:
                        blocks.append((chrom, start, end))
            elif chrom_end > chrom_start:
                blocks.append((chrom, chrom_start, chrom_end))
    return blocks


def probe_overlap_table(
    probe_blocks: list[tuple[str, int, int]],
    n_regions: int,
    ref_name_to_id: dict[str, int],
    region_arrays: RegionArrays,
) -> pd.DataFrame:
    overlap_bp = np.zeros(n_regions, dtype=np.float64)
    probe_block_count = np.zeros(n_regions, dtype=np.float64)
    for ref_name, start, end in probe_blocks:
        overlaps = overlap_regions(ref_name, [(start, end)], ref_name_to_id, region_arrays)
        for region_idx, bp in overlaps.items():
            overlap_bp[region_idx] += float(bp)
            probe_block_count[region_idx] += 1.0
    return pd.DataFrame(
        {
            "probe_overlap_bp": overlap_bp,
            "probe_block_count": probe_block_count,
            "has_probe_overlap": overlap_bp > 0.0,
        }
    )


def prepare_region_table(
    index: TranscriptIndex,
    result,
    condition_info: dict,
    sim_base: Path,
    annotated_bam: Path,
) -> tuple[pd.DataFrame, dict[str, object]]:
    region_arrays = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload_arrays = PayloadArrays.from_payload(result.calibration_payload, region_arrays)
    ledger = build_region_count_ledger(payload_arrays)
    observation = build_density_observation(region_arrays, ledger, result.calibration.fl_models.gdna)
    boundaries = build_boundary_table(region_arrays, ledger, observation.boundary_left_leff)

    calibration = result.calibration.region_calibration
    background = result.calibration.background_model
    local = result.calibration.boundary_local
    sweep = result.calibration.boundary_sweep
    p_states = np.asarray(calibration.p_states, dtype=np.float64)

    sorted_region_df = index.region_df.iloc[region_arrays.order].reset_index(drop=True).copy()
    table = sorted_region_df[
        [
            "region_id",
            "ref_name",
            "start",
            "end",
            "length",
            "signature",
            "type",
            "strand",
            "intron_pos",
            "intron_neg",
            "exon_pos",
            "exon_neg",
            "boundary_kind_left",
            "boundary_kind_right",
        ]
    ].copy()
    table["region_type"] = table["type"].map(region_type_label)
    table["midpoint"] = (table["start"].astype(float) + table["end"].astype(float)) / 2.0

    table["contained_unspliced_count"] = observation.contained_count
    table["boundary_unspliced_count"] = observation.boundary_count
    table["observed_compatible_count"] = observation.observed_compatible_count
    table["spliced_count"] = observation.spliced_count
    table["ts_class"] = region_arrays.ts_class
    table["contained_unspliced_pos"] = ledger.contained_unspliced_pos
    table["contained_unspliced_neg"] = ledger.contained_unspliced_neg
    table["boundary_left_unspliced_pos"] = ledger.boundary_left_unspliced_pos
    table["boundary_left_unspliced_neg"] = ledger.boundary_left_unspliced_neg
    table["boundary_right_unspliced_pos"] = ledger.boundary_right_unspliced_pos
    table["boundary_right_unspliced_neg"] = ledger.boundary_right_unspliced_neg
    table["contained_spliced_pos"] = ledger.contained_spliced_pos
    table["contained_spliced_neg"] = ledger.contained_spliced_neg
    table["boundary_left_spliced_pos"] = ledger.boundary_left_spliced_pos
    table["boundary_left_spliced_neg"] = ledger.boundary_left_spliced_neg
    table["boundary_right_spliced_pos"] = ledger.boundary_right_spliced_pos
    table["boundary_right_spliced_neg"] = ledger.boundary_right_spliced_neg
    table["contained_leff"] = observation.contained_leff
    table["boundary_leff"] = observation.boundary_leff
    table["observed_density"] = safe_div(
        observation.observed_compatible_count,
        np.maximum(observation.contained_leff, 0.0),
    )
    table["contained_density"] = safe_div(
        observation.contained_count,
        np.maximum(observation.contained_leff, 0.0),
    )
    table["is_anchor"] = observation.is_anchor
    table["anchor_intergenic"] = observation.anchor_intergenic
    table["anchor_intron"] = observation.anchor_intron

    for idx, name in enumerate(STATE_NAMES):
        table[f"p_state_{name}"] = p_states[:, idx]
    table["state_name"] = [STATE_NAMES[int(i)] for i in np.argmax(p_states, axis=1)]
    table["p_captured"] = calibration.p_captured
    table["p_expressed"] = calibration.p_expressed
    table["mu_gdna"] = calibration.mu_gdna
    table["upper_gdna"] = calibration.upper_gdna
    table["rna_lower"] = calibration.rna_lower
    table["prior_total"] = calibration.prior_mass.unspliced_total
    table["prior_gdna"] = calibration.prior_mass.gdna_unspliced_mean
    table["prior_rna"] = calibration.prior_mass.rna_unspliced_mean
    table["estimated_gdna_fraction"] = safe_div(table["prior_gdna"], table["prior_total"])
    table["A_r"] = calibration.A_r
    table["gamma_r"] = calibration.gamma_r
    table["calibration_flags"] = calibration.flags
    table["background_seed"] = background.seed_mask
    table["background_top_t_excluded"] = background.top_t_exclusion_mask
    table["background_flags"] = background.flags
    table["boundary_local_alpha"] = local.alpha_excess
    table["boundary_local_beta"] = local.beta_excess
    table["boundary_local_mu"] = local.mu_local
    table["boundary_sweep_alpha"] = sweep.alpha_excess
    table["boundary_sweep_beta"] = sweep.beta_excess
    table["boundary_sweep_mu"] = sweep.mu_sweep
    table["boundary_sweep_upper"] = sweep.upper_sweep

    strand_channels = result.calibration.strand_channels
    if strand_channels is None:
        for name in [
            "contained_mean",
            "contained_upper",
            "contained_rna_lower",
            "contained_precision",
            "boundary_left_mean",
            "boundary_left_upper",
            "boundary_left_rna_lower",
            "boundary_left_precision",
            "boundary_right_mean",
            "boundary_right_upper",
            "boundary_right_rna_lower",
            "boundary_right_precision",
            "flags",
        ]:
            table[f"strand_{name}"] = np.nan
    else:
        table["strand_contained_mean"] = strand_channels.contained_mean
        table["strand_contained_upper"] = strand_channels.contained_upper
        table["strand_contained_rna_lower"] = strand_channels.contained_rna_lower
        table["strand_contained_precision"] = strand_channels.contained_precision
        table["strand_boundary_left_mean"] = strand_channels.boundary_left_mean
        table["strand_boundary_left_upper"] = strand_channels.boundary_left_upper
        table["strand_boundary_left_rna_lower"] = strand_channels.boundary_left_rna_lower
        table["strand_boundary_left_precision"] = strand_channels.boundary_left_precision
        table["strand_boundary_right_mean"] = strand_channels.boundary_right_mean
        table["strand_boundary_right_upper"] = strand_channels.boundary_right_upper
        table["strand_boundary_right_rna_lower"] = strand_channels.boundary_right_rna_lower
        table["strand_boundary_right_precision"] = strand_channels.boundary_right_precision
        table["strand_flags"] = strand_channels.flags
        table["strand_total_mean"] = (
            table["strand_contained_mean"]
            + table["strand_boundary_left_mean"]
            + table["strand_boundary_right_mean"]
        )
        table["strand_total_rna_lower"] = (
            table["strand_contained_rna_lower"]
            + table["strand_boundary_left_rna_lower"]
            + table["strand_boundary_right_rna_lower"]
        )

    probe_path = resolve_suite_path(
        sim_base,
        condition_info.get("capture_probe_bed") or condition_info.get("capture_probe_panel"),
    )
    probe_table = probe_overlap_table(
        load_bed12_probe_blocks(probe_path),
        len(table),
        index.ref_name_to_id,
        region_arrays,
    )
    table = pd.concat([table, probe_table], axis=1)

    truth_table, truth_summary = add_truth_projection(
        annotated_bam,
        table,
        index.ref_name_to_id,
        region_arrays,
        t_id_to_index,
        t_index_to_exons,
        t_index_to_strand,
    )
    table = pd.concat([table, truth_table], axis=1)
    table["deconv_gdna_error"] = table["prior_gdna"] - table["true_gdna_mass"]
    table["deconv_rna_error"] = table["prior_rna"] - table["true_rna_mass"]
    table["absolute_deconv_error"] = np.abs(table["deconv_gdna_error"])
    table["true_gdna_density"] = safe_div(table["true_gdna_mass"], table["contained_leff"])
    table["true_rna_density"] = safe_div(table["true_rna_mass"], table["contained_leff"])
    table["gdna_to_rna_rate_major"] = safe_div(
        table["true_gdna_assigned_rna_major"],
        table["true_gdna_major_count"],
    )

    diagnostics = {
        "truth_projection": truth_summary,
        "n_probe_blocks": int(len(load_bed12_probe_blocks(probe_path))),
        "probe_path": None if probe_path is None else str(probe_path),
        "boundary_count": int(boundaries.n_boundaries),
    }
    return table, diagnostics


def weighted_corr(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 2:
        return 0.0
    if float(np.std(x[mask])) == 0.0 or float(np.std(y[mask])) == 0.0:
        return 0.0
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def aggregate_groups(table: pd.DataFrame) -> pd.DataFrame:
    grouped = table.assign(
        probe_bin=np.where(table["has_probe_overlap"], "probe_overlap", "no_probe"),
        captured_bin=np.where(table["p_captured"] >= 0.5, "p_captured>=0.5", "p_captured<0.5"),
    )
    keys = ["probe_bin", "captured_bin", "state_name", "region_type"]
    cols = [
        "true_gdna_mass",
        "true_rna_mass",
        "prior_gdna",
        "prior_rna",
        "observed_compatible_count",
        "contained_leff",
        "true_gdna_assigned_rna_major",
        "true_gdna_assigned_gdna_major",
        "true_rna_assigned_gdna_major",
        "true_rna_assigned_rna_major",
    ]
    out = grouped.groupby(keys, dropna=False)[cols].sum().reset_index()
    out.insert(0, "n_regions", grouped.groupby(keys, dropna=False).size().to_numpy())
    out["true_gdna_fraction"] = safe_div(out["true_gdna_mass"], out["true_gdna_mass"] + out["true_rna_mass"])
    out["estimated_gdna_fraction"] = safe_div(out["prior_gdna"], out["prior_gdna"] + out["prior_rna"])
    out["true_gdna_assigned_rna_rate"] = safe_div(
        out["true_gdna_assigned_rna_major"],
        out["true_gdna_assigned_rna_major"] + out["true_gdna_assigned_gdna_major"],
    )
    return out.sort_values(["true_gdna_assigned_rna_major", "true_gdna_mass"], ascending=False)


def summarize(table: pd.DataFrame, result, metrics_row: dict[str, object] | None) -> dict[str, object]:
    background = result.calibration.background_model
    calibration = result.calibration.region_calibration
    seed = np.asarray(table["background_seed"], dtype=bool)
    has_probe = np.asarray(table["has_probe_overlap"], dtype=bool)
    captured = np.asarray(table["p_captured"] >= 0.5, dtype=bool)
    leff = np.asarray(table["contained_leff"], dtype=np.float64)

    def density(mask: np.ndarray, col: str) -> float:
        return float(safe_div(table.loc[mask, col].sum(), leff[mask].sum())) if mask.any() else 0.0

    summary = {
        "rho_off": float(background.rho_off_mean),
        "background_seed_regions": int(background.n_seed_regions),
        "background_seed_observed_density": density(seed, "contained_unspliced_count"),
        "background_seed_true_gdna_density": density(seed, "true_gdna_mass"),
        "background_seed_true_rna_density": density(seed, "true_rna_mass"),
        "probe_true_gdna_density": density(has_probe, "true_gdna_mass"),
        "nonprobe_true_gdna_density": density(~has_probe, "true_gdna_mass"),
        "captured_true_gdna_density": density(captured, "true_gdna_mass"),
        "noncaptured_true_gdna_density": density(~captured, "true_gdna_mass"),
        "probe_estimated_gdna_density": density(has_probe, "prior_gdna"),
        "nonprobe_estimated_gdna_density": density(~has_probe, "prior_gdna"),
        "captured_estimated_gdna_density": density(captured, "prior_gdna"),
        "noncaptured_estimated_gdna_density": density(~captured, "prior_gdna"),
        "true_gdna_mass_total": float(table["true_gdna_mass"].sum()),
        "true_rna_mass_total": float(table["true_rna_mass"].sum()),
        "prior_gdna_total": float(table["prior_gdna"].sum()),
        "prior_rna_total": float(table["prior_rna"].sum()),
        "true_gdna_assigned_rna_major": float(table["true_gdna_assigned_rna_major"].sum()),
        "true_gdna_assigned_gdna_major": float(table["true_gdna_assigned_gdna_major"].sum()),
        "true_rna_assigned_gdna_major": float(table["true_rna_assigned_gdna_major"].sum()),
        "true_rna_assigned_rna_major": float(table["true_rna_assigned_rna_major"].sum()),
        "prior_vs_truth_gdna_corr": weighted_corr(table["true_gdna_mass"], table["prior_gdna"]),
        "mu_vs_truth_gdna_corr": weighted_corr(table["true_gdna_mass"], table["mu_gdna"]),
        "n_regions": int(len(table)),
        "n_probe_overlap_regions": int(has_probe.sum()),
        "n_p_captured_ge_0_5_regions": int(captured.sum()),
        "capture_enrichment_target": float(calibration.capture_enrichment_target),
        "gamma_r_p95": float(np.quantile(np.asarray(table["gamma_r"], dtype=np.float64), 0.95)),
        "A_r_p95": float(np.quantile(np.asarray(table["A_r"], dtype=np.float64), 0.95)),
        "strand_contained_mean_total": float(np.nansum(table.get("strand_contained_mean", 0.0))),
        "strand_total_mean_total": float(np.nansum(table.get("strand_total_mean", 0.0))),
        "strand_total_rna_lower_total": float(
            np.nansum(table.get("strand_total_rna_lower", 0.0))
        ),
        "strand_contained_precision_p95": float(
            np.nanquantile(np.asarray(table.get("strand_contained_precision", np.nan)), 0.95)
        )
        if "strand_contained_precision" in table.columns
        else float("nan"),
    }
    if metrics_row:
        for key, value in metrics_row.items():
            if isinstance(value, (np.integer, np.floating)):
                summary[f"metrics_{key}"] = float(value)
            elif isinstance(value, (int, float, str)):
                summary[f"metrics_{key}"] = value
    return summary


def write_json(path: Path, data: dict[str, object]) -> None:
    with path.open("w") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)
        handle.write("\n")


def save_scatter_truth_vs_est(table: pd.DataFrame, out_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, y_col, title in [
        (axes[0], "prior_gdna", "Prior gDNA mass vs oracle"),
        (axes[1], "mu_gdna", "Calibration mu_gDNA vs oracle"),
    ]:
        colors = np.where(table["has_probe_overlap"], "tab:red", "tab:blue")
        ax.scatter(
            table["true_gdna_mass"],
            table[y_col],
            c=colors,
            s=20 + 80 * safe_div(table["observed_compatible_count"], table["observed_compatible_count"].max()),
            alpha=0.75,
            linewidths=0.3,
            edgecolors="black",
        )
        max_val = max(float(table["true_gdna_mass"].max()), float(table[y_col].max()), 1.0)
        ax.plot([0, max_val], [0, max_val], color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("Oracle true gDNA mass by region")
        ax.set_ylabel(y_col)
        ax.set_title(title)
        ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "truth_vs_estimated_gdna.png", dpi=180)
    plt.close(fig)


def save_fraction_scatter(table: pd.DataFrame, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    scatter = ax.scatter(
        table["true_gdna_fraction"],
        table["estimated_gdna_fraction"],
        c=table["p_captured"],
        cmap="viridis",
        s=30 + 90 * safe_div(table["observed_compatible_count"], table["observed_compatible_count"].max()),
        alpha=0.85,
        linewidths=0.3,
        edgecolors="black",
    )
    ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Oracle gDNA fraction by region")
    ax.set_ylabel("Calibration prior gDNA fraction")
    ax.set_title("Per-region deconvolution fraction")
    ax.grid(alpha=0.2)
    fig.colorbar(scatter, ax=ax, label="p_captured")
    fig.tight_layout()
    fig.savefig(out_dir / "region_gdna_fraction_deconvolution.png", dpi=180)
    plt.close(fig)


def save_genome_map(table: pd.DataFrame, out_dir: Path) -> None:
    x = np.arange(len(table), dtype=np.int64)
    labels = table["region_id"].astype(str).to_numpy()
    fig, axes = plt.subplots(4, 1, figsize=(15, 11), sharex=True)

    axes[0].bar(x, table["true_rna_mass"], color="tab:green", alpha=0.75, label="truth RNA")
    axes[0].bar(
        x,
        table["true_gdna_mass"],
        bottom=table["true_rna_mass"],
        color="tab:purple",
        alpha=0.75,
        label="truth gDNA",
    )
    axes[0].set_ylabel("Oracle mass")
    axes[0].legend(loc="upper right")

    axes[1].bar(x, table["prior_rna"], color="tab:green", alpha=0.55, label="prior RNA")
    axes[1].bar(
        x,
        table["prior_gdna"],
        bottom=table["prior_rna"],
        color="tab:purple",
        alpha=0.55,
        label="prior gDNA",
    )
    axes[1].plot(x, table["mu_gdna"], color="black", linewidth=1, label="mu_gDNA")
    axes[1].set_ylabel("Calibration mass")
    axes[1].legend(loc="upper right")

    axes[2].plot(x, table["p_captured"], color="tab:orange", label="p_captured")
    axes[2].plot(x, table["p_expressed"], color="tab:blue", label="p_expressed")
    axes[2].fill_between(
        x,
        0,
        np.where(table["has_probe_overlap"], 1.0, 0.0),
        color="tab:red",
        alpha=0.15,
        label="probe overlap",
    )
    axes[2].set_ylabel("State probability")
    axes[2].set_ylim(-0.05, 1.05)
    axes[2].legend(loc="upper right")

    axes[3].bar(x, table["true_gdna_assigned_rna_major"], color="tab:red", alpha=0.8, label="true gDNA -> RNA")
    axes[3].bar(x, -table["true_rna_assigned_gdna_major"], color="tab:blue", alpha=0.8, label="true RNA -> gDNA")
    axes[3].axhline(0, color="black", linewidth=0.8)
    axes[3].set_ylabel("Major-region errors")
    axes[3].set_xlabel("Calibration region order")
    axes[3].legend(loc="upper right")

    step = max(1, len(x) // 20)
    axes[3].set_xticks(x[::step])
    axes[3].set_xticklabels(labels[::step], rotation=90, fontsize=7)
    for ax in axes:
        ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(out_dir / "genome_region_calibration_map.png", dpi=180)
    plt.close(fig)


def save_background_density(table: pd.DataFrame, result, out_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5.5))
    seed = np.asarray(table["background_seed"], dtype=bool)
    probe = np.asarray(table["has_probe_overlap"], dtype=bool)
    colors = np.where(seed, "tab:green", np.where(probe, "tab:red", "tab:gray"))
    ax.scatter(
        table["contained_leff"],
        table["contained_density"],
        c=colors,
        s=35,
        alpha=0.8,
        linewidths=0.3,
        edgecolors="black",
    )
    ax.axhline(result.calibration.background_model.rho_off_mean, color="black", linestyle="--", label="rho_off")
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=0.01)
    ax.set_xlabel("Contained gDNA effective length")
    ax.set_ylabel("Observed contained unspliced density")
    ax.set_title("Background fit seed regions")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_dir / "background_density_fit.png", dpi=180)
    plt.close(fig)


def save_fl_plot(truth_fl_path: Path, result, out_dir: Path) -> pd.DataFrame:
    truth = pd.read_csv(truth_fl_path, sep="\t")
    max_len = 450
    lengths = np.arange(max_len + 1)
    models = {
        "model_rna": result.calibration.fl_models.rna.pmf[: max_len + 1],
        "model_gdna": result.calibration.fl_models.gdna.pmf[: max_len + 1],
        "model_rna_scoring": result.calibration.fl_models.rna_scoring.pmf[: max_len + 1],
        "model_gdna_scoring": result.calibration.fl_models.gdna_scoring.pmf[: max_len + 1],
    }
    model_df = pd.DataFrame({"fragment_length": lengths, **models})
    truth_pivot = truth.pivot_table(
        index="fragment_length",
        columns="kind",
        values="fraction",
        aggfunc="sum",
        fill_value=0.0,
    ).reset_index()
    merged = model_df.merge(truth_pivot, on="fragment_length", how="left").fillna(0.0)

    fig, ax = plt.subplots(figsize=(10, 5.5))
    if "mrna" in merged.columns:
        ax.plot(merged["fragment_length"], merged["mrna"], color="tab:green", label="truth mRNA")
    if "gdna" in merged.columns:
        ax.plot(merged["fragment_length"], merged["gdna"], color="tab:purple", label="truth gDNA")
    ax.plot(merged["fragment_length"], merged["model_rna"], color="darkgreen", linestyle="--", label="cal RNA FL")
    ax.plot(merged["fragment_length"], merged["model_gdna"], color="indigo", linestyle="--", label="cal gDNA FL")
    ax.set_xlim(0, max_len)
    ax.set_xlabel("Fragment length")
    ax.set_ylabel("Probability")
    ax.set_title("Post-capture truth FL vs learned FL models")
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(out_dir / "fragment_length_truth_vs_model.png", dpi=180)
    plt.close(fig)
    return merged


def markdown_table(df: pd.DataFrame, max_rows: int = 10) -> str:
    view = df.head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: f"{x:.4g}")
    columns = [str(col) for col in view.columns]
    rows = [[str(value) for value in row] for row in view.to_numpy()]
    widths = [len(col) for col in columns]
    for row in rows:
        for idx, value in enumerate(row):
            widths[idx] = max(widths[idx], len(value))

    def fmt(row: list[str]) -> str:
        cells = [value.ljust(widths[idx]) for idx, value in enumerate(row)]
        return "| " + " | ".join(cells) + " |"

    header = fmt(columns)
    separator = "| " + " | ".join("-" * width for width in widths) + " |"
    body = [fmt(row) for row in rows]
    return "\n".join([header, separator, *body])


def write_summary_md(
    out_dir: Path,
    condition: str,
    summary: dict[str, object],
    diagnostics: dict[str, object],
    top_misfit: pd.DataFrame,
    group_summary: pd.DataFrame,
) -> None:
    lines: list[str] = []
    lines.append(f"# Calibration diagnostic: {condition}")
    lines.append("")
    lines.append("## High-level calibration")
    lines.append("")
    lines.append(f"- Learned rho_off: {summary['rho_off']:.6g}")
    lines.append(
        "- Background seed true gDNA density: "
        f"{summary['background_seed_true_gdna_density']:.6g}; "
        f"observed seed density: {summary['background_seed_observed_density']:.6g}"
    )
    lines.append(
        "- Probe-overlap true gDNA density: "
        f"{summary['probe_true_gdna_density']:.6g}; non-probe true gDNA density: "
        f"{summary['nonprobe_true_gdna_density']:.6g}"
    )
    lines.append(
        "- Probe-overlap estimated gDNA density: "
        f"{summary['probe_estimated_gdna_density']:.6g}; non-probe estimated gDNA density: "
        f"{summary['nonprobe_estimated_gdna_density']:.6g}"
    )
    lines.append(
        "- Prior-vs-oracle gDNA correlation: "
        f"{summary['prior_vs_truth_gdna_corr']:.4f}; mu-vs-oracle correlation: "
        f"{summary['mu_vs_truth_gdna_corr']:.4f}"
    )
    lines.append(
        "- Major-region true gDNA assigned RNA: "
        f"{summary['true_gdna_assigned_rna_major']:.0f}; true gDNA assigned gDNA: "
        f"{summary['true_gdna_assigned_gdna_major']:.0f}"
    )
    lines.append(
        "- Capture enrichment target: "
        f"{summary['capture_enrichment_target']:.4g}; gamma_r p95: {summary['gamma_r_p95']:.4g}"
    )
    lines.append(
        "- Strand deconv contained gDNA mean total: "
        f"{summary['strand_contained_mean_total']:.4g}; total compartment gDNA mean: "
        f"{summary['strand_total_mean_total']:.4g}; RNA lower total: "
        f"{summary['strand_total_rna_lower_total']:.4g}"
    )
    lines.append(
        "- Strand contained precision p95: "
        f"{summary['strand_contained_precision_p95']:.4g}"
    )
    lines.append("")
    lines.append("## Projection diagnostics")
    lines.append("")
    lines.append(f"- Qnames seen: {diagnostics['truth_projection']['n_qnames_seen']}")
    lines.append(f"- Fragments projected: {diagnostics['truth_projection']['n_fragments_projected']}")
    lines.append(f"- Skipped: `{diagnostics['truth_projection']['skipped']}`")
    lines.append(f"- Probe panel: `{diagnostics.get('probe_path')}`")
    lines.append("")
    lines.append("## Top misfit regions")
    lines.append("")
    lines.append(markdown_table(top_misfit, max_rows=12))
    lines.append("")
    lines.append("## Grouped error summary")
    lines.append("")
    lines.append(markdown_table(group_summary, max_rows=16))
    lines.append("")
    lines.append("## Plot files")
    lines.append("")
    for name in [
        "truth_vs_estimated_gdna.png",
        "region_gdna_fraction_deconvolution.png",
        "genome_region_calibration_map.png",
        "background_density_fit.png",
        "fragment_length_truth_vs_model.png",
    ]:
        lines.append(f"- {name}")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


def read_metrics_row(sim_base: Path, condition: str) -> dict[str, object] | None:
    path = sim_base / "condition_metrics.tsv"
    if not path.exists():
        return None
    df = pd.read_csv(path, sep="\t")
    row = df.loc[df["condition"] == condition]
    if row.empty:
        return None
    return row.iloc[0].to_dict()


t_index_to_row: dict[int, pd.Series] = {}
t_id_to_index: dict[str, int] = {}
t_index_to_exons: dict[int, np.ndarray] = {}
t_index_to_strand: dict[int, int] = {}


def initialize_transcript_maps(index: TranscriptIndex) -> None:
    global t_index_to_row, t_id_to_index, t_index_to_exons, t_index_to_strand
    t_df = index.t_df.reset_index(drop=True)
    t_index_to_row = {int(row["t_index"]): row for _, row in t_df.iterrows()}
    t_id_to_index = {str(row["t_id"]): int(row["t_index"]) for _, row in t_df.iterrows()}
    t_index_to_strand = {int(row["t_index"]): int(row["strand"]) for _, row in t_df.iterrows()}
    t_index_to_exons = {
        int(row["t_index"]): np.asarray(index.get_exon_intervals(int(row["t_index"])), dtype=np.int64)
        for _, row in t_df.iterrows()
    }


def main() -> None:
    args = parse_args()
    sim_base = args.sim_base.resolve()
    condition_dir = sim_base / args.condition
    out_dir = args.out_dir or (sim_base / "diagnostics" / args.condition)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = load_json(sim_base / "manifest.json")
    condition_info = find_condition(manifest, args.condition)
    index = TranscriptIndex.load(sim_base / "rigel_index")
    initialize_transcript_maps(index)

    bam_path = condition_dir / "sim_oracle.bam"
    annotated_bam = condition_dir / "annotated.bam"
    truth_fl_path = condition_dir / "truth_fragment_lengths.tsv"
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)
    if not annotated_bam.exists():
        raise FileNotFoundError(annotated_bam)

    config = PipelineConfig(
        scan=BamScanConfig(sj_strand_tag="auto"),
        annotated_bam_path=None,
        emit_locus_stats=False,
    )
    print(f"[diagnose] Running pipeline for {args.condition} to expose calibration arrays")
    result = run_pipeline(str(bam_path), index, config)
    if result.calibration is None or result.calibration_payload is None:
        raise RuntimeError("Pipeline did not produce calibration arrays.")

    print("[diagnose] Projecting oracle truth and assignments onto calibration regions")
    region_table, diagnostics = prepare_region_table(
        index,
        result,
        condition_info,
        sim_base,
        annotated_bam,
    )
    group_summary = aggregate_groups(region_table)
    top_misfit_cols = [
        "region_id",
        "ref_name",
        "start",
        "end",
        "region_type",
        "has_probe_overlap",
        "state_name",
        "p_captured",
        "p_expressed",
        "true_gdna_mass",
        "true_rna_mass",
        "prior_gdna",
        "prior_rna",
        "mu_gdna",
        "deconv_gdna_error",
        "true_gdna_assigned_rna_major",
        "true_gdna_assigned_gdna_major",
        "contained_unspliced_count",
        "spliced_count",
        "contained_leff",
        "boundary_sweep_mu",
        "gamma_r",
        "contained_unspliced_pos",
        "contained_unspliced_neg",
        "contained_spliced_pos",
        "contained_spliced_neg",
        "strand_contained_mean",
        "strand_contained_rna_lower",
        "strand_contained_precision",
        "strand_total_mean",
        "strand_total_rna_lower",
    ]
    top_misfit = region_table.sort_values("absolute_deconv_error", ascending=False)[top_misfit_cols]

    print("[diagnose] Writing TSV and JSON outputs")
    region_table.to_csv(out_dir / "per_region_calibration_truth.tsv", sep="\t", index=False)
    group_summary.to_csv(out_dir / "grouped_error_summary.tsv", sep="\t", index=False)
    top_misfit.to_csv(out_dir / "top_misfit_regions.tsv", sep="\t", index=False)
    seed_cols = [
        "region_id",
        "ref_name",
        "start",
        "end",
        "region_type",
        "background_seed",
        "background_top_t_excluded",
        "has_probe_overlap",
        "contained_unspliced_count",
        "true_gdna_mass",
        "true_rna_mass",
        "contained_leff",
        "contained_density",
        "true_gdna_density",
        "p_captured",
        "state_name",
    ]
    region_table.loc[region_table["is_anchor"], seed_cols].to_csv(
        out_dir / "background_candidate_regions.tsv",
        sep="\t",
        index=False,
    )
    metrics_row = read_metrics_row(sim_base, args.condition)
    summary = summarize(region_table, result, metrics_row)
    write_json(out_dir / "diagnostic_summary.json", {"summary": summary, "diagnostics": diagnostics})

    print("[diagnose] Generating plots")
    save_scatter_truth_vs_est(region_table, out_dir)
    save_fraction_scatter(region_table, out_dir)
    save_genome_map(region_table, out_dir)
    save_background_density(region_table, result, out_dir)
    fl_table = save_fl_plot(truth_fl_path, result, out_dir)
    fl_table.to_csv(out_dir / "fragment_length_truth_vs_model.tsv", sep="\t", index=False)

    write_summary_md(out_dir, args.condition, summary, diagnostics, top_misfit, group_summary)
    print(f"[diagnose] Done. Diagnostics written to {out_dir}")


if __name__ == "__main__":
    main()