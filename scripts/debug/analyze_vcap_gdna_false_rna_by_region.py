#!/usr/bin/env python3
"""Analyze VCaP gDNA->RNA errors by calibration region class."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08

TYPE_INTERGENIC = 0
TYPE_INTRON = 1
TYPE_EXON = 2
INTERVAL_EXON = 0
INTERVAL_TRANSCRIPT = 1

MASK_EXON = 1
MASK_INTRON = 2
MASK_INTERGENIC = 4

MASK_LABELS = {
    0: "NONE",
    MASK_EXON: "EXON_ONLY",
    MASK_INTRON: "INTRON_ONLY",
    MASK_EXON | MASK_INTRON: "EXON_INTRON",
    MASK_INTERGENIC: "INTERGENIC_ONLY",
    MASK_EXON | MASK_INTERGENIC: "EXON_INTERGENIC",
    MASK_INTRON | MASK_INTERGENIC: "INTRON_INTERGENIC",
    MASK_EXON | MASK_INTRON | MASK_INTERGENIC: "EXON_INTRON_INTERGENIC",
}

TYPE_LABELS = {
    -1: "NONE",
    TYPE_INTERGENIC: "INTERGENIC",
    TYPE_INTRON: "INTRON",
    TYPE_EXON: "EXON",
}

STRAND_LABELS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True, help="Rigel annotated BAM.")
    parser.add_argument("--index-dir", type=Path, required=True, help="Rigel index directory.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--report", type=Path, required=True, help="Markdown report path.")
    parser.add_argument("--rna-flowcell", action="append", default=[], help="RNA flowcell ID.")
    parser.add_argument("--gdna-flowcell", action="append", default=[], help="gDNA flowcell ID.")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--splicing-anchor-tolerance", type=int, default=3)
    parser.add_argument("--top-regions", type=int, default=100)
    return parser.parse_args()


@dataclass(frozen=True, slots=True)
class RefRegions:
    region_ids: np.ndarray
    starts: np.ndarray
    ends: np.ndarray
    types: np.ndarray


class RegionLookup:
    def __init__(self, region_df: pd.DataFrame) -> None:
        self.regions_by_ref: dict[str, RefRegions] = {}
        for ref, sub in region_df.groupby("ref_name", sort=False):
            ordered = sub.sort_values("start", kind="stable")
            self.regions_by_ref[str(ref)] = RefRegions(
                region_ids=ordered["region_id"].to_numpy(np.int64),
                starts=ordered["start"].to_numpy(np.int64),
                ends=ordered["end"].to_numpy(np.int64),
                types=ordered["type"].to_numpy(np.uint8),
            )

    def overlaps(self, ref_name: str, start: int, end: int) -> Iterable[tuple[int, int]]:
        ref_regions = self.regions_by_ref.get(ref_name)
        if ref_regions is None or end <= start:
            return ()
        starts = ref_regions.starts
        ends = ref_regions.ends
        hi = int(np.searchsorted(starts, end, side="left"))
        lo = int(np.searchsorted(ends, start, side="right"))
        if lo >= hi:
            return ()
        out: list[tuple[int, int]] = []
        for pos in range(lo, hi):
            overlap = min(end, int(ends[pos])) - max(start, int(starts[pos]))
            if overlap > 0:
                out.append((int(ref_regions.region_ids[pos]), int(overlap)))
        return out


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def flowcell_from_qname(query_name: str) -> str:
    parts = query_name.split(":")
    return parts[2] if len(parts) >= 3 else query_name


def source_map_from_args(args: argparse.Namespace) -> dict[str, str]:
    source_map: dict[str, str] = {}
    for flowcell in args.rna_flowcell:
        source_map[flowcell] = "rna"
    for flowcell in args.gdna_flowcell:
        if source_map.get(flowcell) == "rna":
            raise ValueError(f"Flowcell {flowcell} assigned to both RNA and gDNA")
        source_map[flowcell] = "gdna"
    return source_map


def predicted_detail_from_zf(zf_value: int | None) -> str:
    if zf_value is None:
        return "unresolved"
    if zf_value & AF_GDNA_BIT:
        return "gdna"
    if zf_value & AF_NRNA_BIT:
        return "nrna"
    if zf_value & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def binary_pool(detail: str) -> str:
    if detail in {"mrna", "nrna"}:
        return "rna"
    return detail


def type_mask(region_type: int) -> int:
    return 1 << (2 - int(region_type))


def broad_mask_label(mask: int) -> str:
    if mask & MASK_EXON:
        return "EXON_CONTAINING"
    if mask & MASK_INTRON:
        return "INTRON_ONLY"
    if mask & MASK_INTERGENIC:
        return "INTERGENIC_ONLY"
    return "NONE"


def length_bin(length: int) -> str:
    if length <= 100:
        return "1-100"
    if length <= 250:
        return "101-250"
    if length <= 500:
        return "251-500"
    if length <= 1000:
        return "501-1k"
    if length <= 5000:
        return "1k-5k"
    if length <= 10000:
        return "5k-10k"
    if length <= 50000:
        return "10k-50k"
    return ">50k"


def choose_representative(records: list[pysam.AlignedSegment]) -> pysam.AlignedSegment | None:
    for read in records:
        if read.is_paired and read.is_read1:
            return read
    return records[0] if records else None


def primary_records(records: list[pysam.AlignedSegment]) -> list[pysam.AlignedSegment]:
    return [read for read in records if not read.is_secondary and not read.is_supplementary]


def fragment_blocks(records: list[pysam.AlignedSegment]) -> list[tuple[str, int, int]]:
    blocks: list[tuple[str, int, int]] = []
    for read in records:
        if read.is_unmapped or read.reference_name is None:
            continue
        for start, end in read.get_blocks():
            if end > start:
                blocks.append((str(read.reference_name), int(start), int(end)))
    return blocks


def region_class_for_fragment(
    lookup: RegionLookup,
    region_types: np.ndarray,
    blocks: list[tuple[str, int, int]],
    min_overlap: int,
) -> tuple[int, int, int, dict[int, int]]:
    overlap_by_region: dict[int, int] = defaultdict(int)
    for ref_name, start, end in blocks:
        for region_id, overlap in lookup.overlaps(ref_name, start, end):
            overlap_by_region[region_id] += overlap

    mask = 0
    dominant_region_id = -1
    dominant_overlap = 0
    for region_id, overlap in overlap_by_region.items():
        if overlap < min_overlap:
            continue
        region_type = int(region_types[region_id])
        mask |= type_mask(region_type)
        if overlap > dominant_overlap:
            dominant_region_id = region_id
            dominant_overlap = overlap

    dominant_type = int(region_types[dominant_region_id]) if dominant_region_id >= 0 else -1
    return mask, dominant_region_id, dominant_type, dict(overlap_by_region)


def iter_fragment_groups(bam: pysam.AlignmentFile) -> Iterable[list[pysam.AlignedSegment]]:
    current_qname: str | None = None
    group: list[pysam.AlignedSegment] = []
    for read in bam.fetch(until_eof=True):
        qname = read.query_name
        if current_qname is None:
            current_qname = qname
        if qname != current_qname:
            yield group
            group = []
            current_qname = qname
        group.append(read)
    if group:
        yield group


def load_inputs(index_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    regions = pd.read_feather(index_dir / "regions.feather")
    transcripts = pd.read_feather(index_dir / "transcripts.feather")
    intervals = pd.read_feather(index_dir / "intervals.feather")
    regions["ref_name"] = regions["ref_name"].astype(str)
    intervals["ref"] = intervals["ref"].astype(str)
    transcripts["ref"] = transcripts["ref"].astype(str)
    return regions, transcripts, intervals


def scan_bam(args: argparse.Namespace, region_df: pd.DataFrame) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    lookup = RegionLookup(region_df)
    region_types = region_df.sort_values("region_id", kind="stable")["type"].to_numpy(np.uint8)
    n_regions = len(region_df)
    min_overlap = max(int(args.splicing_anchor_tolerance), 1)

    record_counts = Counter()
    mask_counts = Counter()
    dominant_counts = Counter()
    dominant_region_total = np.zeros(n_regions, dtype=np.int64)
    dominant_region_gdna_total = np.zeros(n_regions, dtype=np.int64)
    dominant_region_gdna_false_rna = np.zeros(n_regions, dtype=np.int64)
    dominant_region_gdna_false_mrna = np.zeros(n_regions, dtype=np.int64)
    dominant_region_gdna_false_nrna = np.zeros(n_regions, dtype=np.int64)
    dominant_region_rna_total = np.zeros(n_regions, dtype=np.int64)
    dominant_region_rna_false_gdna = np.zeros(n_regions, dtype=np.int64)
    exact_mask_gdna_total = Counter()
    exact_mask_gdna_false_rna = Counter()

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for group in iter_fragment_groups(bam):
            record_counts["qname_groups_seen"] += 1
            primaries = primary_records(group)
            if not primaries:
                record_counts["groups_no_primary"] += 1
                continue
            representative = choose_representative(primaries)
            if representative is None:
                record_counts["groups_no_representative"] += 1
                continue

            record_counts["fragments_counted"] += 1
            true_source = source_map.get(flowcell_from_qname(representative.query_name), "unknown")
            detail = predicted_detail_from_zf(get_tag(representative, "ZF", None))
            predicted_binary = binary_pool(detail)
            blocks = fragment_blocks(primaries)
            mask, dominant_region_id, dominant_type, _overlaps = region_class_for_fragment(
                lookup, region_types, blocks, min_overlap
            )
            mask_label = MASK_LABELS.get(mask, f"MASK_{mask}")
            broad_label = broad_mask_label(mask)
            dominant_label = TYPE_LABELS.get(dominant_type, str(dominant_type))

            mask_counts[(true_source, predicted_binary, detail, mask_label, broad_label)] += 1
            dominant_counts[(true_source, predicted_binary, detail, dominant_label)] += 1
            if true_source == "gdna":
                exact_mask_gdna_total[mask_label] += 1
            if dominant_region_id >= 0:
                dominant_region_total[dominant_region_id] += 1
                if true_source == "gdna":
                    dominant_region_gdna_total[dominant_region_id] += 1
                elif true_source == "rna":
                    dominant_region_rna_total[dominant_region_id] += 1

            if true_source == "gdna" and predicted_binary == "rna":
                exact_mask_gdna_false_rna[mask_label] += 1
                if dominant_region_id >= 0:
                    dominant_region_gdna_false_rna[dominant_region_id] += 1
                    if detail == "mrna":
                        dominant_region_gdna_false_mrna[dominant_region_id] += 1
                    elif detail == "nrna":
                        dominant_region_gdna_false_nrna[dominant_region_id] += 1
            if true_source == "rna" and predicted_binary == "gdna" and dominant_region_id >= 0:
                dominant_region_rna_false_gdna[dominant_region_id] += 1

    return {
        "source_map": source_map,
        "record_counts": record_counts,
        "mask_counts": mask_counts,
        "dominant_counts": dominant_counts,
        "dominant_region_total": dominant_region_total,
        "dominant_region_gdna_total": dominant_region_gdna_total,
        "dominant_region_gdna_false_rna": dominant_region_gdna_false_rna,
        "dominant_region_gdna_false_mrna": dominant_region_gdna_false_mrna,
        "dominant_region_gdna_false_nrna": dominant_region_gdna_false_nrna,
        "dominant_region_rna_total": dominant_region_rna_total,
        "dominant_region_rna_false_gdna": dominant_region_rna_false_gdna,
        "exact_mask_gdna_total": exact_mask_gdna_total,
        "exact_mask_gdna_false_rna": exact_mask_gdna_false_rna,
        "min_overlap": min_overlap,
    }


def build_interval_lookup(intervals: pd.DataFrame) -> dict[str, pd.DataFrame]:
    by_ref: dict[str, pd.DataFrame] = {}
    for ref, sub in intervals.groupby("ref", sort=False):
        by_ref[str(ref)] = sub.sort_values("start", kind="stable").reset_index(drop=True)
    return by_ref


def interval_complexity(
    region: pd.Series,
    intervals_by_ref: dict[str, pd.DataFrame],
    transcripts: pd.DataFrame,
) -> dict[str, Any]:
    ref = str(region["ref_name"])
    start = int(region["start"])
    end = int(region["end"])
    sub = intervals_by_ref.get(ref)
    if sub is None or sub.empty:
        return {
            "exon_interval_count": 0,
            "exon_transcript_count": 0,
            "exon_gene_count": 0,
            "transcript_span_count": 0,
            "transcript_span_gene_count": 0,
        }
    starts = sub["start"].to_numpy(np.int64)
    hi = int(np.searchsorted(starts, end, side="left"))
    if hi <= 0:
        overlapping = sub.iloc[0:0]
    else:
        candidates = sub.iloc[:hi]
        overlapping = candidates[candidates["end"] > start]
    exon_hits = overlapping[overlapping["interval_type"] == INTERVAL_EXON]
    span_hits = overlapping[overlapping["interval_type"] == INTERVAL_TRANSCRIPT]
    exon_t = pd.Index(exon_hits["t_index"].astype(int).unique())
    span_t = pd.Index(span_hits["t_index"].astype(int).unique())
    valid_exon_t = exon_t[(exon_t >= 0) & (exon_t < len(transcripts))]
    valid_span_t = span_t[(span_t >= 0) & (span_t < len(transcripts))]
    exon_gene_count = int(transcripts.iloc[valid_exon_t]["g_index"].nunique()) if len(valid_exon_t) else 0
    span_gene_count = int(transcripts.iloc[valid_span_t]["g_index"].nunique()) if len(valid_span_t) else 0
    return {
        "exon_interval_count": int(len(exon_hits)),
        "exon_transcript_count": int(len(valid_exon_t)),
        "exon_gene_count": exon_gene_count,
        "transcript_span_count": int(len(valid_span_t)),
        "transcript_span_gene_count": span_gene_count,
    }


def build_rows(
    args: argparse.Namespace,
    region_df: pd.DataFrame,
    transcripts: pd.DataFrame,
    intervals: pd.DataFrame,
    scan: dict[str, Any],
) -> dict[str, list[dict[str, Any]]]:
    mask_rows = []
    for (true_source, predicted_binary, detail, mask_label, broad_label), count in sorted(
        scan["mask_counts"].items()
    ):
        mask_rows.append(
            {
                "true_source": true_source,
                "predicted_binary": predicted_binary,
                "predicted_detail": detail,
                "region_mask": mask_label,
                "broad_region": broad_label,
                "count": count,
            }
        )

    dominant_rows = []
    for (true_source, predicted_binary, detail, dominant_type), count in sorted(
        scan["dominant_counts"].items()
    ):
        dominant_rows.append(
            {
                "true_source": true_source,
                "predicted_binary": predicted_binary,
                "predicted_detail": detail,
                "dominant_region_type": dominant_type,
                "count": count,
            }
        )

    total_gdna = sum(scan["exact_mask_gdna_total"].values())
    total_false = sum(scan["exact_mask_gdna_false_rna"].values())
    false_mask_rows = []
    for mask_label in sorted(set(scan["exact_mask_gdna_total"]) | set(scan["exact_mask_gdna_false_rna"])):
        gdna_total = scan["exact_mask_gdna_total"].get(mask_label, 0)
        false_count = scan["exact_mask_gdna_false_rna"].get(mask_label, 0)
        false_mask_rows.append(
            {
                "region_mask": mask_label,
                "broad_region": broad_mask_label(next((k for k, v in MASK_LABELS.items() if v == mask_label), 0)),
                "gdna_source_total": gdna_total,
                "gdna_false_rna": false_count,
                "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
                "fraction_of_gdna_source": gdna_total / total_gdna if total_gdna else 0.0,
                "fraction_of_false_rna": false_count / total_false if total_false else 0.0,
            }
        )
    false_mask_rows.sort(key=lambda row: -row["gdna_false_rna"])

    region_by_id = region_df.sort_values("region_id", kind="stable").reset_index(drop=True)
    top_region_ids = np.argsort(scan["dominant_region_gdna_false_rna"])[::-1]
    top_region_ids = [int(rid) for rid in top_region_ids if scan["dominant_region_gdna_false_rna"][rid] > 0]
    top_region_ids = top_region_ids[: args.top_regions]
    intervals_by_ref = build_interval_lookup(intervals)

    top_region_rows = []
    for rid in top_region_ids:
        region = region_by_id.iloc[rid]
        length = int(region["end"] - region["start"])
        complexity = interval_complexity(region, intervals_by_ref, transcripts)
        gdna_total = int(scan["dominant_region_gdna_total"][rid])
        false_count = int(scan["dominant_region_gdna_false_rna"][rid])
        top_region_rows.append(
            {
                "region_id": rid,
                "ref_name": str(region["ref_name"]),
                "start_1based": int(region["start"]) + 1,
                "end": int(region["end"]),
                "length_bp": length,
                "region_type": TYPE_LABELS.get(int(region["type"]), str(int(region["type"]))),
                "region_strand": STRAND_LABELS.get(int(region["strand"]), str(int(region["strand"]))),
                "boundary_flux_left": bool(region["boundary_flux_left"]),
                "boundary_flux_right": bool(region["boundary_flux_right"]),
                "tx_pos_bp": int(region["tx_pos_bp"]),
                "tx_neg_bp": int(region["tx_neg_bp"]),
                "exon_pos_bp": int(region["exon_pos_bp"]),
                "exon_neg_bp": int(region["exon_neg_bp"]),
                "gdna_source_total": gdna_total,
                "gdna_false_rna": false_count,
                "gdna_false_mrna": int(scan["dominant_region_gdna_false_mrna"][rid]),
                "gdna_false_nrna": int(scan["dominant_region_gdna_false_nrna"][rid]),
                "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
                "rna_source_total": int(scan["dominant_region_rna_total"][rid]),
                "rna_false_gdna": int(scan["dominant_region_rna_false_gdna"][rid]),
                **complexity,
            }
        )

    bin_counter = Counter()
    for rid, false_count in enumerate(scan["dominant_region_gdna_false_rna"]):
        gdna_total = int(scan["dominant_region_gdna_total"][rid])
        if gdna_total <= 0:
            continue
        region = region_by_id.iloc[rid]
        region_type = TYPE_LABELS.get(int(region["type"]), str(int(region["type"])))
        length = int(region["end"] - region["start"])
        key = (
            region_type,
            length_bin(length),
            STRAND_LABELS.get(int(region["strand"]), str(int(region["strand"]))),
            bool(region["boundary_flux_left"]),
            bool(region["boundary_flux_right"]),
        )
        bin_counter[(*key, "gdna_source_total")] += gdna_total
        bin_counter[(*key, "gdna_false_rna")] += int(false_count)
        bin_counter[(*key, "regions")] += 1

    bin_rows = []
    for key_and_metric, value in list(bin_counter.items()):
        *_key, metric = key_and_metric
        if metric != "regions":
            continue
        region_type, len_bin, strand, bfl, bfr = _key
        gdna_total = bin_counter[(region_type, len_bin, strand, bfl, bfr, "gdna_source_total")]
        false_count = bin_counter[(region_type, len_bin, strand, bfl, bfr, "gdna_false_rna")]
        bin_rows.append(
            {
                "region_type": region_type,
                "length_bin": len_bin,
                "region_strand": strand,
                "boundary_flux_left": bfl,
                "boundary_flux_right": bfr,
                "regions_with_gdna": value,
                "gdna_source_total": gdna_total,
                "gdna_false_rna": false_count,
                "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
            }
        )
    bin_rows.sort(key=lambda row: -row["gdna_false_rna"])

    return {
        "mask_rows": mask_rows,
        "dominant_rows": dominant_rows,
        "false_mask_rows": false_mask_rows,
        "top_region_rows": top_region_rows,
        "bin_rows": bin_rows,
    }


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def fmt_int(value: int) -> str:
    return f"{value:,}"


def pct(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * value:.2f}%"


def markdown_table(records: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not records:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _ in columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for record in records:
        body.append("| " + " | ".join(str(record.get(key, "")) for _, key in columns) + " |")
    return "\n".join([header, sep, *body])


def report_false_mask_rows(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    out = []
    for row in rows:
        out.append(
            {
                "mask": row["region_mask"],
                "broad": row["broad_region"],
                "gdna_total": fmt_int(int(row["gdna_source_total"])),
                "false": fmt_int(int(row["gdna_false_rna"])),
                "rate": pct(float(row["gdna_false_rna_rate"])),
                "false_frac": pct(float(row["fraction_of_false_rna"])),
            }
        )
    return out


def report_top_region_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "region": f"{row['ref_name']}:{fmt_int(int(row['start_1based']))}-{fmt_int(int(row['end']))}",
                "type": row["region_type"],
                "len": fmt_int(int(row["length_bp"])),
                "strand": row["region_strand"],
                "false": fmt_int(int(row["gdna_false_rna"])),
                "rate": pct(float(row["gdna_false_rna_rate"])),
                "mrna": fmt_int(int(row["gdna_false_mrna"])),
                "nrna": fmt_int(int(row["gdna_false_nrna"])),
                "exon_tx": fmt_int(int(row["exon_transcript_count"])),
                "exon_genes": fmt_int(int(row["exon_gene_count"])),
                "span_genes": fmt_int(int(row["transcript_span_gene_count"])),
                "boundary": f"{int(bool(row['boundary_flux_left']))}/{int(bool(row['boundary_flux_right']))}",
            }
        )
    return out


def report_bin_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "type": row["region_type"],
                "len_bin": row["length_bin"],
                "strand": row["region_strand"],
                "boundary": f"{int(bool(row['boundary_flux_left']))}/{int(bool(row['boundary_flux_right']))}",
                "regions": fmt_int(int(row["regions_with_gdna"])),
                "gdna_total": fmt_int(int(row["gdna_source_total"])),
                "false": fmt_int(int(row["gdna_false_rna"])),
                "rate": pct(float(row["gdna_false_rna_rate"])),
            }
        )
    return out


def write_report(args: argparse.Namespace, scan: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    total_false = sum(scan["exact_mask_gdna_false_rna"].values())
    total_gdna = sum(scan["exact_mask_gdna_total"].values())
    exon_false = sum(
        row["gdna_false_rna"]
        for row in rows["false_mask_rows"]
        if "EXON" in row["region_mask"]
    )
    intron_only_false = scan["exact_mask_gdna_false_rna"].get("INTRON_ONLY", 0)
    intergenic_only_false = scan["exact_mask_gdna_false_rna"].get("INTERGENIC_ONLY", 0)
    top20_false = sum(int(row["gdna_false_rna"]) for row in rows["top_region_rows"][:20])
    top100_false = sum(int(row["gdna_false_rna"]) for row in rows["top_region_rows"][:100])

    text = f"""# VCaP gDNA-to-RNA Error By Region - 2026-05-19

BAM: `{args.bam}`

Index: `{args.index_dir}`

This analysis classifies each counted fragment by overlap with Rigel's calibration
`regions.feather` partition. It uses both mates' primary aligned blocks, applies
minimum region overlap `q={scan['min_overlap']}` bp, and then counts the hard Rigel
winner from the fragment-level `ZF` tag.

Fragments counted: {fmt_int(int(scan['record_counts']['fragments_counted']))}

True gDNA fragments: {fmt_int(total_gdna)}

True gDNA -> RNA calls: {fmt_int(total_false)} ({pct(total_false / total_gdna if total_gdna else None)})

## Region Headline

- EXON-containing masks explain {fmt_int(int(exon_false))} gDNA->RNA calls ({pct(exon_false / total_false if total_false else None)} of false RNA).
- INTRON-only explains {fmt_int(int(intron_only_false))} calls ({pct(intron_only_false / total_false if total_false else None)}).
- INTERGENIC-only explains {fmt_int(int(intergenic_only_false))} calls ({pct(intergenic_only_false / total_false if total_false else None)}).

## Exact Region Mask

{markdown_table(report_false_mask_rows(rows['false_mask_rows']), [
        ('Region mask', 'mask'),
        ('Broad class', 'broad'),
        ('gDNA total', 'gdna_total'),
        ('gDNA->RNA', 'false'),
        ('False rate', 'rate'),
        ('False RNA frac', 'false_frac'),
    ])}

## Top Dominant Regions

Dominant region means the qualified region with the largest aligned-block overlap
for the fragment. Boundary is `left/right` from `regions.feather`.

{markdown_table(report_top_region_rows(rows['top_region_rows'], 40), [
        ('Region', 'region'),
        ('Type', 'type'),
        ('Len bp', 'len'),
        ('Strand', 'strand'),
        ('gDNA->RNA', 'false'),
        ('Rate', 'rate'),
        ('mRNA', 'mrna'),
        ('nRNA', 'nrna'),
        ('Exon tx', 'exon_tx'),
        ('Exon genes', 'exon_genes'),
        ('Span genes', 'span_genes'),
        ('Boundary', 'boundary'),
    ])}

Top dominant-region concentration:

- Top 20 regions: {fmt_int(top20_false)} ({pct(top20_false / total_false if total_false else None)} of false RNA)
- Top 100 regions: {fmt_int(top100_false)} ({pct(top100_false / total_false if total_false else None)} of false RNA)

## Dominant Region Bins

{markdown_table(report_bin_rows(rows['bin_rows'], 50), [
        ('Type', 'type'),
        ('Len bin', 'len_bin'),
        ('Strand', 'strand'),
        ('Boundary', 'boundary'),
        ('Regions', 'regions'),
        ('gDNA total', 'gdna_total'),
        ('gDNA->RNA', 'false'),
        ('False rate', 'rate'),
    ])}

## Artifacts

- Exact mask confusion: `{args.out_dir / 'regional_confusion_by_mask.tsv'}`
- Dominant type confusion: `{args.out_dir / 'dominant_region_type_confusion.tsv'}`
- gDNA false RNA by mask: `{args.out_dir / 'gdna_false_rna_by_mask.tsv'}`
- Top dominant regions: `{args.out_dir / 'dominant_region_gdna_false_rna.tsv'}`
- Dominant region bins: `{args.out_dir / 'dominant_region_bins.tsv'}`
- Summary JSON: `{args.out_dir / 'summary.json'}`
"""
    args.report.write_text(text)


def write_outputs(
    args: argparse.Namespace,
    scan: dict[str, Any],
    rows: dict[str, list[dict[str, Any]]],
) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "regional_confusion_by_mask.tsv",
        rows["mask_rows"],
        ["true_source", "predicted_binary", "predicted_detail", "region_mask", "broad_region", "count"],
    )
    write_tsv(
        args.out_dir / "dominant_region_type_confusion.tsv",
        rows["dominant_rows"],
        ["true_source", "predicted_binary", "predicted_detail", "dominant_region_type", "count"],
    )
    write_tsv(
        args.out_dir / "gdna_false_rna_by_mask.tsv",
        rows["false_mask_rows"],
        [
            "region_mask",
            "broad_region",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
            "fraction_of_gdna_source",
            "fraction_of_false_rna",
        ],
    )
    write_tsv(
        args.out_dir / "dominant_region_gdna_false_rna.tsv",
        rows["top_region_rows"],
        [
            "region_id",
            "ref_name",
            "start_1based",
            "end",
            "length_bp",
            "region_type",
            "region_strand",
            "boundary_flux_left",
            "boundary_flux_right",
            "tx_pos_bp",
            "tx_neg_bp",
            "exon_pos_bp",
            "exon_neg_bp",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_mrna",
            "gdna_false_nrna",
            "gdna_false_rna_rate",
            "rna_source_total",
            "rna_false_gdna",
            "exon_interval_count",
            "exon_transcript_count",
            "exon_gene_count",
            "transcript_span_count",
            "transcript_span_gene_count",
        ],
    )
    write_tsv(
        args.out_dir / "dominant_region_bins.tsv",
        rows["bin_rows"],
        [
            "region_type",
            "length_bin",
            "region_strand",
            "boundary_flux_left",
            "boundary_flux_right",
            "regions_with_gdna",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
        ],
    )
    summary = {
        "bam": str(args.bam),
        "index_dir": str(args.index_dir),
        "source_map": scan["source_map"],
        "record_counts": dict(scan["record_counts"]),
        "min_overlap": scan["min_overlap"],
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    write_report(args, scan, rows)


def main() -> None:
    args = parse_args()
    region_df, transcripts, intervals = load_inputs(args.index_dir)
    scan = scan_bam(args, region_df)
    rows = build_rows(args, region_df, transcripts, intervals, scan)
    write_outputs(args, scan, rows)
    print(f"counted {scan['record_counts']['fragments_counted']:,} fragments")
    print(f"gDNA->RNA by exact region mask written to {args.out_dir}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()
