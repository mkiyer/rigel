#!/usr/bin/env python3
"""Focused diagnosis for the chr1 FLG2/antisense gDNA false-RNA hotspot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08

STRAND_LABEL = {0: ".", 1: "+", 2: "-", 3: "?"}
INTERVAL_LABEL = {0: "EXON", 1: "TRANSCRIPT", 2: "INTERGENIC"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bam",
        type=Path,
        default=Path(
            "/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/"
            "exon_strand_deconv_v1/annotated.bam"
        ),
    )
    parser.add_argument(
        "--quant-dir",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1"),
    )
    parser.add_argument(
        "--index-dir",
        type=Path,
        default=Path("/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index"),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/flg2_hotspot_diagnostics_2026-05-21"),
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("docs/benchmarks/flg2_hotspot_diagnostics_2026-05-21.md"),
    )
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--start", type=int, default=152_350_001, help="1-based inclusive")
    parser.add_argument("--end", type=int, default=152_360_000, help="1-based inclusive")
    parser.add_argument("--rna-flowcell", action="append", default=["C6EL5ANXX"])
    parser.add_argument("--gdna-flowcell", action="append", default=["H7MFFDSXY"])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sample-limit", type=int, default=80)
    return parser.parse_args()


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any = None) -> Any:
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
    return "rna" if detail in {"mrna", "nrna"} else detail


def is_representative(read: pysam.AlignedSegment) -> bool:
    if read.is_secondary or read.is_supplementary:
        return False
    return (not read.is_paired) or read.is_read1


def fragment_span(read: pysam.AlignedSegment) -> tuple[int, int]:
    start = read.reference_start if read.reference_start >= 0 else 0
    end = read.reference_end if read.reference_end is not None else start + max(read.query_length or 1, 1)
    if read.is_paired and read.next_reference_id == read.reference_id and read.template_length:
        if read.template_length > 0:
            start = read.reference_start
            end = read.reference_start + read.template_length
        else:
            start = read.next_reference_start
            end = read.next_reference_start - read.template_length
    if end < start:
        start, end = end, start
    return max(start, 0), max(end, start + 1)


def overlaps(start: int, end: int, region_start: int, region_end: int) -> bool:
    return start < region_end and end > region_start


def overlap_bp(start: int, end: int, region_start: int, region_end: int) -> int:
    return max(0, min(end, region_end) - max(start, region_start))


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged: list[tuple[int, int]] = []
    current_start, current_end = ordered[0]
    for interval_start, interval_end in ordered[1:]:
        if interval_start <= current_end:
            current_end = max(current_end, interval_end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = interval_start, interval_end
    merged.append((current_start, current_end))
    return merged


def subtract_intervals(
    intervals: list[tuple[int, int]],
    masks: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    pieces: list[tuple[int, int]] = []
    merged_masks = merge_intervals(masks)
    for interval_start, interval_end in merge_intervals(intervals):
        cursor = interval_start
        for mask_start, mask_end in merged_masks:
            if mask_end <= cursor:
                continue
            if mask_start >= interval_end:
                break
            if mask_start > cursor:
                pieces.append((cursor, min(mask_start, interval_end)))
            cursor = max(cursor, mask_end)
            if cursor >= interval_end:
                break
        if cursor < interval_end:
            pieces.append((cursor, interval_end))
    return [(start, end) for start, end in pieces if end > start]


def overlaps_any(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    return any(overlaps(start, end, interval_start, interval_end) for interval_start, interval_end in intervals)


def parse_nrna_region(target: str) -> str:
    if not target.startswith("RIGEL_NRNA_"):
        return ""
    try:
        prefix, strand_value, start_value, end_value = target.rsplit("_", 3)
    except ValueError:
        return ""
    chrom = prefix.removeprefix("RIGEL_NRNA_")
    strand_label = {"1": "+", "2": "-"}.get(strand_value, strand_value)
    return f"{chrom}:{int(start_value):,}-{int(end_value):,}({strand_label})"


def target_label(read: pysam.AlignedSegment, predicted_detail: str) -> str:
    transcript = str(get_tag(read, "ZT", "."))
    gene = str(get_tag(read, "ZR", "."))
    gene_id = str(get_tag(read, "ZG", "."))
    if predicted_detail == "nrna":
        region = parse_nrna_region(transcript)
        return f"nRNA {region or transcript}"
    if predicted_detail == "mrna":
        return f"mRNA {gene if gene not in {'', '.'} else gene_id} {transcript}"
    return predicted_detail


def compact(value: str, max_len: int = 90) -> str:
    return value if len(value) <= max_len else value[: max_len - 3] + "..."


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int | None, digits: int = 3) -> str:
    if value is None:
        return "n/a"
    return f"{float(value):.{digits}f}"


def fmt_pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * float(value):.{digits}f}%"


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    sep = "| " + " | ".join("---" for _label, _key in columns) + " |"
    body = ["| " + " | ".join(str(row.get(key, "")) for _label, key in columns) + " |" for row in rows]
    return "\n".join([header, sep, *body])


def load_context(args: argparse.Namespace) -> dict[str, Any]:
    transcript_df = pd.read_feather(args.index_dir / "transcripts.feather")
    interval_df = pd.read_feather(args.index_dir / "intervals.feather")
    region_df = pd.read_feather(args.index_dir / "regions.feather")
    quant_df = pd.read_feather(args.quant_dir / "quant.feather")
    nrna_df = pd.read_feather(args.quant_dir / "nrna_quant.feather")
    loci_df = pd.read_feather(args.quant_dir / "loci.feather")
    summary = json.loads((args.quant_dir / "summary.json").read_text())
    return {
        "transcripts": transcript_df,
        "intervals": interval_df,
        "regions": region_df,
        "quant": quant_df,
        "nrna": nrna_df,
        "loci": loci_df,
        "summary": summary,
    }


def build_feature_intervals(
    args: argparse.Namespace,
    context: dict[str, Any],
) -> dict[str, Any]:
    start0 = args.start - 1
    end0 = args.end
    transcript_df: pd.DataFrame = context["transcripts"]
    interval_df: pd.DataFrame = context["intervals"]

    window_transcripts = transcript_df[
        (transcript_df["ref"] == args.chrom)
        & (transcript_df["end"] > start0)
        & (transcript_df["start"] < end0)
    ].copy()
    annotated_positive = window_transcripts[
        (~window_transcripts["is_synthetic"].astype(bool))
        & (window_transcripts["strand"] == 1)
    ].copy()
    flg2 = transcript_df[transcript_df["t_id"] == "ENST00000388718.5"].copy()
    flg2_t_indices = set(int(value) for value in flg2["t_index"].tolist())

    flg2_exon_rows = interval_df[
        (interval_df["t_index"].isin(flg2_t_indices))
        & (interval_df["interval_type"] == 0)
    ]
    flg2_exons = merge_intervals(
        [(int(row.start), int(row.end)) for row in flg2_exon_rows.itertuples(index=False)]
    )

    positive_t_indices = set(int(value) for value in annotated_positive["t_index"].tolist())
    positive_rows = interval_df[interval_df["t_index"].isin(positive_t_indices)]
    positive_span_rows = positive_rows[positive_rows["interval_type"] == 1]
    positive_exon_rows = positive_rows[positive_rows["interval_type"] == 0]
    positive_spans = merge_intervals(
        [(int(row.start), int(row.end)) for row in positive_span_rows.itertuples(index=False)]
    )
    positive_exons = merge_intervals(
        [(int(row.start), int(row.end)) for row in positive_exon_rows.itertuples(index=False)]
    )
    positive_introns = subtract_intervals(positive_spans, positive_exons)
    positive_exons_outside_window = subtract_intervals(positive_exons, [(start0, end0)])
    positive_introns_outside_window = subtract_intervals(positive_introns, [(start0, end0)])
    broad_intervals = merge_intervals(positive_spans + [(start0, end0)] + flg2_exons)
    broad_start = min(start for start, _end in broad_intervals)
    broad_end = max(end for _start, end in broad_intervals)

    return {
        "window_start0": start0,
        "window_end0": end0,
        "broad_start0": broad_start,
        "broad_end0": broad_end,
        "window_transcripts": window_transcripts,
        "annotated_positive": annotated_positive,
        "flg2_exons": flg2_exons,
        "positive_spans": positive_spans,
        "positive_exons": positive_exons,
        "positive_introns": positive_introns,
        "positive_exons_outside_window": positive_exons_outside_window,
        "positive_introns_outside_window": positive_introns_outside_window,
    }


def read1_strand(read: pysam.AlignedSegment) -> str:
    return "-" if read.is_reverse else "+"


def inferred_rna_strand_from_r1(read: pysam.AlignedSegment) -> str:
    # The VCaP run inferred R1-antisense. For an RNA molecule, R1 aligns to
    # the opposite genomic strand of the transcript.
    return "+" if read.is_reverse else "-"


def scan_bam(args: argparse.Namespace, features: dict[str, Any]) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    window_start = int(features["window_start0"])
    window_end = int(features["window_end0"])
    broad_start = int(features["broad_start0"])
    broad_end = int(features["broad_end0"])
    feature_intervals = {
        "hotspot_window": [(window_start, window_end)],
        "flg2_exons": features["flg2_exons"],
        "antisense_exons_outside_window": features["positive_exons_outside_window"],
        "antisense_introns_outside_window": features["positive_introns_outside_window"],
    }

    records = Counter()
    window_confusion = Counter()
    window_orientation = Counter()
    window_targets = Counter()
    false_targets = Counter()
    false_category = Counter()
    feature_counts = defaultdict(Counter)
    sample_rows: list[dict[str, Any]] = []
    fragment_lengths: list[int] = []

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen"] += 1
            if not is_representative(read):
                records["nonrepresentative_skipped"] += 1
                continue
            if read.reference_name != args.chrom:
                continue
            span_start, span_end = fragment_span(read)
            if not overlaps(span_start, span_end, broad_start, broad_end):
                continue
            records["broad_fragments"] += 1

            flowcell = flowcell_from_qname(read.query_name)
            true_source = source_map.get(flowcell, "unknown")
            predicted_detail = predicted_detail_from_zf(get_tag(read, "ZF"))
            predicted_pool = binary_pool(predicted_detail)
            read_strand = read1_strand(read)
            rna_strand = inferred_rna_strand_from_r1(read)
            zc_value = str(get_tag(read, "ZC", "missing"))
            zs_value = str(get_tag(read, "ZS", "missing"))
            target = target_label(read, predicted_detail)
            template_len = abs(int(read.template_length)) if read.template_length else span_end - span_start

            for feature_name, intervals in feature_intervals.items():
                if overlaps_any(span_start, span_end, intervals):
                    feature_counts[feature_name][(true_source, predicted_pool, predicted_detail)] += 1
                    feature_counts[feature_name][(true_source, "read1", read_strand)] += 1
                    feature_counts[feature_name][(true_source, "rna_strand", rna_strand)] += 1

            if not overlaps(span_start, span_end, window_start, window_end):
                continue

            records["window_fragments"] += 1
            window_confusion[(true_source, predicted_pool, predicted_detail)] += 1
            window_orientation[(true_source, predicted_detail, read_strand, rna_strand)] += 1
            window_targets[(true_source, predicted_detail, target)] += 1
            fragment_lengths.append(template_len)

            if true_source == "gdna" and predicted_pool == "rna":
                false_targets[(predicted_detail, target)] += 1
                false_category[(predicted_detail, zc_value, zs_value, read_strand, rna_strand)] += 1
                if len(sample_rows) < args.sample_limit:
                    sample_rows.append(
                        {
                            "qname": read.query_name,
                            "span": f"{args.chrom}:{span_start + 1}-{span_end}",
                            "predicted_detail": predicted_detail,
                            "target": target,
                            "read1_strand": read_strand,
                            "rna_strand_if_rna": rna_strand,
                            "template_length": template_len,
                            "cigar": read.cigarstring or "",
                            "nh": get_tag(read, "NH", ""),
                            "nm": get_tag(read, "NM", ""),
                            "zw": get_tag(read, "ZW", ""),
                            "zc": zc_value,
                            "zs": zs_value,
                            "zt": get_tag(read, "ZT", ""),
                            "zg": get_tag(read, "ZG", ""),
                            "zr": get_tag(read, "ZR", ""),
                        }
                    )

    return {
        "records": records,
        "window_confusion": window_confusion,
        "window_orientation": window_orientation,
        "window_targets": window_targets,
        "false_targets": false_targets,
        "false_category": false_category,
        "feature_counts": feature_counts,
        "sample_rows": sample_rows,
        "fragment_lengths": fragment_lengths,
    }


def build_transcript_rows(context: dict[str, Any], features: dict[str, Any]) -> list[dict[str, Any]]:
    window_start = int(features["window_start0"])
    window_end = int(features["window_end0"])
    interval_df: pd.DataFrame = context["intervals"]
    quant_df: pd.DataFrame = context["quant"]
    nrna_df: pd.DataFrame = context["nrna"]
    rows: list[dict[str, Any]] = []

    quant_by_tid = quant_df.set_index("transcript_id", drop=False)
    nrna_by_id = nrna_df.set_index("nrna_id", drop=False)

    for transcript in features["window_transcripts"].sort_values(
        ["is_synthetic", "strand", "start", "end", "t_id"]
    ).itertuples(index=False):
        transcript_id = str(transcript.t_id)
        t_index = int(transcript.t_index)
        exon_rows = interval_df[
            (interval_df["t_index"] == t_index)
            & (interval_df["interval_type"] == 0)
        ]
        exon_overlap = int(
            sum(
                overlap_bp(int(row.start), int(row.end), window_start, window_end)
                for row in exon_rows.itertuples(index=False)
            )
        )
        quant_count = None
        quant_em = None
        tpm = None
        if transcript_id in quant_by_tid.index:
            quant_row = quant_by_tid.loc[transcript_id]
            quant_count = float(quant_row["count"])
            quant_em = float(quant_row["count_em"])
            tpm = float(quant_row["tpm"])
        elif transcript_id in nrna_by_id.index:
            nrna_row = nrna_by_id.loc[transcript_id]
            quant_count = float(nrna_row["count"])
            quant_em = float(nrna_row["count"])
            tpm = float(nrna_row["tpm"])
        rows.append(
            {
                "t_index": t_index,
                "transcript_id": transcript_id,
                "gene_name": str(transcript.g_name),
                "strand": STRAND_LABEL.get(int(transcript.strand), str(transcript.strand)),
                "start_1based": int(transcript.start) + 1,
                "end": int(transcript.end),
                "length": int(transcript.length),
                "n_exons": int(transcript.n_exons),
                "is_synthetic": bool(transcript.is_synthetic),
                "span_overlap_bp": overlap_bp(int(transcript.start), int(transcript.end), window_start, window_end),
                "exon_overlap_bp": exon_overlap,
                "count": quant_count,
                "count_em": quant_em,
                "tpm": tpm,
            }
        )
    rows.sort(key=lambda row: (-int(row["exon_overlap_bp"]), bool(row["is_synthetic"]), row["strand"], row["start_1based"]))
    return rows


def build_region_rows(args: argparse.Namespace, context: dict[str, Any]) -> list[dict[str, Any]]:
    start0 = args.start - 1
    end0 = args.end
    region_df: pd.DataFrame = context["regions"]
    rows = []
    region_subset = region_df[
        (region_df["ref_name"] == args.chrom)
        & (region_df["end"] > start0)
        & (region_df["start"] < end0)
    ].copy()
    for region in region_subset.itertuples(index=False):
        rows.append(
            {
                "region_id": int(region.region_id),
                "start_1based": int(region.start) + 1,
                "end": int(region.end),
                "type": int(region.type),
                "strand": STRAND_LABEL.get(int(region.strand), str(region.strand)),
                "tx_pos_bp": int(region.tx_pos_bp),
                "tx_neg_bp": int(region.tx_neg_bp),
                "exon_pos_bp": int(region.exon_pos_bp),
                "exon_neg_bp": int(region.exon_neg_bp),
                "boundary_left": bool(region.boundary_flux_left),
                "boundary_right": bool(region.boundary_flux_right),
            }
        )
    return rows


def counter_rows(counter: Counter, columns: list[str], total: int | None = None) -> list[dict[str, Any]]:
    rows = []
    for key, count in counter.most_common():
        row = dict(zip(columns, key, strict=True))
        row["count"] = int(count)
        if total:
            row["fraction"] = count / total
        rows.append(row)
    return rows


def feature_rows(scan: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for feature_name, counter in scan["feature_counts"].items():
        totals = Counter()
        for key, count in counter.items():
            if len(key) == 3 and key[1] not in {"read1", "rna_strand"}:
                true_source, _pool, _detail = key
                totals[true_source] += count
        for key, count in counter.items():
            if len(key) != 3 or key[1] in {"read1", "rna_strand"}:
                continue
            true_source, predicted_pool, predicted_detail = key
            rows.append(
                {
                    "feature": feature_name,
                    "true_source": true_source,
                    "predicted_pool": predicted_pool,
                    "predicted_detail": predicted_detail,
                    "count": int(count),
                    "feature_true_source_fraction": count / totals[true_source]
                    if totals[true_source]
                    else 0.0,
                }
            )
    rows.sort(key=lambda row: (row["feature"], row["true_source"], -row["count"]))
    return rows


def orientation_rows(scan: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    totals = Counter()
    for (true_source, predicted_detail, _read1, _rna_strand), count in scan[
        "window_orientation"
    ].items():
        totals[(true_source, predicted_detail)] += count
    for (true_source, predicted_detail, read1, rna_strand), count in scan[
        "window_orientation"
    ].most_common():
        denominator = totals[(true_source, predicted_detail)]
        rows.append(
            {
                "true_source": true_source,
                "predicted_detail": predicted_detail,
                "read1_strand": read1,
                "rna_strand_if_rna": rna_strand,
                "count": int(count),
                "fraction": count / denominator if denominator else 0.0,
            }
        )
    return rows


def local_math_rows(context: dict[str, Any], scan: dict[str, Any]) -> list[dict[str, Any]]:
    quant_df: pd.DataFrame = context["quant"]
    nrna_df: pd.DataFrame = context["nrna"]
    loci_df: pd.DataFrame = context["loci"]
    summary = context["summary"]
    locus = loci_df[loci_df["locus_id"] == 3].iloc[0]
    gdna_alpha = float(locus["gdna"] + locus["gdna_prior_count_em"] + 0.5)
    gdna_eff_len = float(locus["gdna_eff_len"])

    fragment_lengths = scan["fragment_lengths"]
    median_fragment_length = int(np.median(fragment_lengths)) if fragment_lengths else 300
    rna_hist = summary["fragment_length"]["rna"]["histogram"]["values"]
    gdna_hist = summary["fragment_length"]["gdna"]["histogram"]["values"]
    has_fl_histograms = bool(rna_hist) and bool(gdna_hist)

    def empirical_log_prob(values: list[float], fragment_length: int) -> float | None:
        if not values or fragment_length <= 0 or fragment_length > len(values):
            return None
        total = float(sum(values))
        return math.log((float(values[fragment_length - 1]) + 1.0) / (total + len(values)))

    log_fl_gdna = empirical_log_prob(gdna_hist, median_fragment_length)
    log_fl_rna = empirical_log_prob(rna_hist, median_fragment_length)
    fl_log_odds = (
        float(log_fl_gdna - log_fl_rna)
        if has_fl_histograms and log_fl_gdna is not None and log_fl_rna is not None
        else None
    )
    log_strand_gdna = math.log(0.5)
    log_strand_rna_good = math.log(0.999249)

    rows: list[dict[str, Any]] = []

    def add_component(label: str, count: float, effective_length: float) -> None:
        alpha = float(count + 0.5)
        density_log_odds = math.log(gdna_alpha) - math.log(gdna_eff_len) - math.log(alpha) + math.log(effective_length)
        ll_log_odds = (fl_log_odds or 0.0) + (log_strand_gdna - log_strand_rna_good)
        rows.append(
            {
                "component": label,
                "component_count": count,
                "component_eff_len": effective_length,
                "gdna_alpha_approx": gdna_alpha,
                "gdna_eff_len": gdna_eff_len,
                "median_fragment_length": median_fragment_length,
                "density_log_odds_gdna_vs_component": density_log_odds,
                "fl_log_odds_gdna_vs_rna": fl_log_odds,
                "fl_log_odds_available": has_fl_histograms,
                "strand_log_odds_gdna_vs_good_rna": log_strand_gdna - log_strand_rna_good,
                "approx_total_log_odds_gdna_vs_component": density_log_odds + ll_log_odds,
            }
        )

    flg2 = quant_df[quant_df["transcript_id"] == "ENST00000388718.5"]
    if not flg2.empty:
        row = flg2.iloc[0]
        add_component("mRNA FLG2 ENST00000388718.5", float(row["count"]), float(row["em_effective_length"]))

    for nrna_id in [
        "RIGEL_NRNA_chr1_1_152335378_152364080",
        "RIGEL_NRNA_chr1_1_152168124_152332686",
    ]:
        hit = nrna_df[nrna_df["nrna_id"] == nrna_id]
        if not hit.empty:
            row = hit.iloc[0]
            add_component(f"nRNA {nrna_id}", float(row["count"]), float(row["em_effective_length"]))
    return rows


def report_rows(records: list[dict[str, Any]], columns: list[tuple[str, str]], max_rows: int = 20) -> str:
    return markdown_table(records[:max_rows], columns)


def write_outputs(
    args: argparse.Namespace,
    context: dict[str, Any],
    features: dict[str, Any],
    scan: dict[str, Any],
) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    transcript_rows = build_transcript_rows(context, features)
    region_rows = build_region_rows(args, context)
    window_confusion_rows = counter_rows(
        scan["window_confusion"], ["true_source", "predicted_pool", "predicted_detail"]
    )
    target_rows = counter_rows(
        scan["false_targets"], ["predicted_detail", "target"],
        total=sum(scan["false_targets"].values()),
    )
    category_rows = counter_rows(
        scan["false_category"],
        ["predicted_detail", "zc", "zs", "read1_strand", "rna_strand_if_rna"],
        total=sum(scan["false_category"].values()),
    )
    feature_summary_rows = feature_rows(scan)
    orientation_summary_rows = orientation_rows(scan)
    math_rows = local_math_rows(context, scan)

    write_tsv(args.out_dir / "overlapping_transcripts.tsv", transcript_rows, list(transcript_rows[0].keys()))
    write_tsv(args.out_dir / "calibration_regions.tsv", region_rows, list(region_rows[0].keys()))
    write_tsv(args.out_dir / "window_confusion.tsv", window_confusion_rows, list(window_confusion_rows[0].keys()))
    write_tsv(args.out_dir / "false_rna_targets.tsv", target_rows, list(target_rows[0].keys()))
    write_tsv(args.out_dir / "false_rna_categories.tsv", category_rows, list(category_rows[0].keys()))
    write_tsv(args.out_dir / "feature_summary.tsv", feature_summary_rows, list(feature_summary_rows[0].keys()))
    write_tsv(args.out_dir / "orientation_summary.tsv", orientation_summary_rows, list(orientation_summary_rows[0].keys()))
    write_tsv(args.out_dir / "local_math.tsv", math_rows, list(math_rows[0].keys()))
    if scan["sample_rows"]:
        write_tsv(args.out_dir / "sample_false_rna_fragments.tsv", scan["sample_rows"], list(scan["sample_rows"][0].keys()))

    summary = {
        "bam": str(args.bam),
        "window": {"chrom": args.chrom, "start_1based": args.start, "end": args.end},
        "records": dict(scan["records"]),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True))

    write_report(
        args,
        context,
        features,
        scan,
        transcript_rows,
        region_rows,
        window_confusion_rows,
        target_rows,
        category_rows,
        feature_summary_rows,
        orientation_summary_rows,
        math_rows,
    )


def write_report(
    args: argparse.Namespace,
    context: dict[str, Any],
    features: dict[str, Any],
    scan: dict[str, Any],
    transcript_rows: list[dict[str, Any]],
    region_rows: list[dict[str, Any]],
    window_confusion_rows: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
    category_rows: list[dict[str, Any]],
    feature_summary_rows: list[dict[str, Any]],
    orientation_summary_rows: list[dict[str, Any]],
    math_rows: list[dict[str, Any]],
) -> None:
    loci_df: pd.DataFrame = context["loci"]
    locus = loci_df[loci_df["locus_id"] == 3].iloc[0]
    gdna_total = sum(
        int(row["count"])
        for row in window_confusion_rows
        if row["true_source"] == "gdna"
    )
    gdna_false_rna = sum(
        int(row["count"])
        for row in window_confusion_rows
        if row["true_source"] == "gdna" and row["predicted_pool"] == "rna"
    )
    rna_total = sum(
        int(row["count"])
        for row in window_confusion_rows
        if row["true_source"] == "rna"
    )
    median_fragment_length = int(np.median(scan["fragment_lengths"])) if scan["fragment_lengths"] else 0

    rendered_math_rows = []
    for row in math_rows:
        rendered_math_rows.append(
            {
                "component": compact(str(row["component"]), 44),
                "count": fmt_int(row["component_count"]),
                "Le": fmt_int(row["component_eff_len"]),
                "density_log_odds": fmt_float(row["density_log_odds_gdna_vs_component"]),
                "fl_log_odds": fmt_float(row["fl_log_odds_gdna_vs_rna"]),
                "strand_log_odds": fmt_float(row["strand_log_odds_gdna_vs_good_rna"]),
                "total_log_odds": fmt_float(row["approx_total_log_odds_gdna_vs_component"]),
            }
        )

    rendered_feature_rows = []
    for row in feature_summary_rows:
        rendered_feature_rows.append(
            {
                "feature": row["feature"],
                "truth": row["true_source"],
                "pred": row["predicted_detail"],
                "count": fmt_int(row["count"]),
                "frac": fmt_pct(row["feature_true_source_fraction"]),
            }
        )

    rendered_orientation_rows = []
    for row in orientation_summary_rows:
        rendered_orientation_rows.append(
            {
                "truth": row["true_source"],
                "pred": row["predicted_detail"],
                "r1": row["read1_strand"],
                "rna_strand": row["rna_strand_if_rna"],
                "count": fmt_int(row["count"]),
                "frac": fmt_pct(row["fraction"]),
            }
        )

    rendered_targets = [
        {
            "pred": row["predicted_detail"],
            "target": compact(str(row["target"]), 55),
            "count": fmt_int(row["count"]),
            "frac": fmt_pct(row["fraction"]),
        }
        for row in target_rows
    ]
    rendered_categories = [
        {
            "pred": row["predicted_detail"],
            "zc": row["zc"],
            "zs": row["zs"],
            "r1": row["read1_strand"],
            "rna_strand": row["rna_strand_if_rna"],
            "count": fmt_int(row["count"]),
            "frac": fmt_pct(row["fraction"]),
        }
        for row in category_rows
    ]
    rendered_transcripts = [
        {
            "transcript": compact(str(row["transcript_id"]), 38),
            "gene": row["gene_name"],
            "strand": row["strand"],
            "kind": "synthetic nRNA" if row["is_synthetic"] else "annotated",
            "exon_bp": fmt_int(row["exon_overlap_bp"]),
            "span_bp": fmt_int(row["span_overlap_bp"]),
            "count": fmt_int(row["count"] or 0),
            "tpm": fmt_float(row["tpm"] or 0.0),
        }
        for row in transcript_rows
    ]
    rendered_regions = [
        {
            "region": f"{fmt_int(row['start_1based'])}-{fmt_int(row['end'])}",
            "type": row["type"],
            "strand": row["strand"],
            "tx+": fmt_int(row["tx_pos_bp"]),
            "tx-": fmt_int(row["tx_neg_bp"]),
            "exon+": fmt_int(row["exon_pos_bp"]),
            "exon-": fmt_int(row["exon_neg_bp"]),
        }
        for row in region_rows
    ]

    text = f"""# FLG2 Hotspot EM Failure Diagnosis - 2026-05-21

Window: `{args.chrom}:{args.start:,}-{args.end:,}`

BAM: `{args.bam}`

The hotspot is solved inside MultiLocus/locus `3`, a chromosome-scale component with {fmt_int(locus['n_transcripts'])} transcript components, {fmt_int(locus['n_nrna_entities'])} nRNA entities, {fmt_int(locus['n_em_fragments'])} EM fragments, and a gDNA effective length of {fmt_int(locus['gdna_eff_len'])} bp. The same component contains {fmt_int(locus['gdna'])} assigned gDNA fragments, but its gDNA mass is spread over a very large effective opportunity.

## Local Confusion

Fragments overlapping the 10 kb window: {fmt_int(scan['records']['window_fragments'])}. True RNA-source fragments in this window are only {fmt_int(rna_total)}, while true gDNA-source fragments are {fmt_int(gdna_total)}. Rigel assigns {fmt_int(gdna_false_rna)} true-gDNA fragments to RNA ({fmt_pct(gdna_false_rna / gdna_total if gdna_total else None)}).

{markdown_table(window_confusion_rows, [('Truth', 'true_source'), ('Pred pool', 'predicted_pool'), ('Detail', 'predicted_detail'), ('Count', 'count')])}

## Local Strand Split

`read1_strand` is the genomic strand of R1. Since this run inferred an R1-antisense protocol, `rna_strand_if_rna` is the transcript strand implied if the fragment were RNA. In the true-gDNA population, the two orientations are both present, which lets the EM explain the window as two RNA components: negative-strand FLG2 mRNA plus positive-strand antisense/nRNA.

{report_rows(rendered_orientation_rows, [('Truth', 'truth'), ('Pred', 'pred'), ('R1', 'r1'), ('RNA strand if RNA', 'rna_strand'), ('Count', 'count'), ('Fraction', 'frac')], 16)}

## False-RNA Targets

{report_rows(rendered_targets, [('Pred', 'pred'), ('Target', 'target'), ('Count', 'count'), ('False-RNA frac', 'frac')], 12)}

## False-RNA Categories

{report_rows(rendered_categories, [('Pred', 'pred'), ('ZC', 'zc'), ('ZS', 'zs'), ('R1', 'r1'), ('RNA strand if RNA', 'rna_strand'), ('Count', 'count'), ('False-RNA frac', 'frac')], 16)}

## Overlapping Transcript Models

The reference index presents FLG2, positive-strand annotated antisense transcripts, and synthetic positive/negative nRNA spans as simultaneous candidates in this window.

{report_rows(rendered_transcripts, [('Transcript', 'transcript'), ('Gene', 'gene'), ('Strand', 'strand'), ('Kind', 'kind'), ('Exon bp in window', 'exon_bp'), ('Span bp in window', 'span_bp'), ('Run count', 'count'), ('TPM', 'tpm')], 28)}

## Calibration Regions In The Window

{report_rows(rendered_regions, [('Region', 'region'), ('Type mask', 'type'), ('Strand', 'strand'), ('tx+', 'tx+'), ('tx-', 'tx-'), ('exon+', 'exon+'), ('exon-', 'exon-')], 20)}

## Nearby Antisense Evidence Check

This table compares the hotspot to the positive-strand antisense transcript span outside the hotspot. It is a read-level way to test the browser observation that the antisense transcript model has little support outside the FLG2-overlapping exon/intron region.

{report_rows(rendered_feature_rows, [('Feature', 'feature'), ('Truth', 'truth'), ('Pred', 'pred'), ('Count', 'count'), ('Within-truth frac', 'frac')], 28)}

## Approximate EM Log-Odds

For one unspliced fragment, the native VBEM E-step compares components roughly as:

```text
log posterior score_k = log_lik_k(fragment) + digamma(alpha_k) - log(effective_length_k) + constant
```

So for gDNA versus a transcript component:

```text
log odds(gDNA / RNA_k)
  ~= [log_lik_gDNA - log_lik_RNA_k]
   + [log(alpha_gDNA) - log(L_gDNA)]
   - [log(alpha_k) - log(L_k)]
```

Using the converged locus outputs as an approximation, the median fragment length in the window is {fmt_int(median_fragment_length)} bp. The saved summary for this run does not include the RNA/gDNA FL histograms, so the table below excludes the FL likelihood ratio and shows the density plus strand terms. Negative total log-odds means the RNA component wins for a strand-compatible fragment before any additional FL-model effect.

{markdown_table(rendered_math_rows, [('Component', 'component'), ('Count', 'count'), ('L_eff', 'Le'), ('Density log-odds', 'density_log_odds'), ('FL log-odds', 'fl_log_odds'), ('Strand log-odds', 'strand_log_odds'), ('Total log-odds', 'total_log_odds')])}

## Interpretation

The browser intuition is right: the local data look like gDNA, but the model can explain the same 50/50 strand mixture as the sum of two strand-specific RNA components. The key mathematical issue is that Rigel's gDNA competitor is one component for the entire MultiLocus, with `L_gDNA={fmt_int(locus['gdna_eff_len'])}` bp in this run. FLG2 has an EM effective length on the order of kilobases, and the synthetic antisense nRNA spans are still far shorter than the MultiLocus gDNA opportunity.

That creates a strong per-base density penalty against gDNA. Even though the locus has many gDNA-assigned fragments overall, the gDNA abundance must be divided by the huge MultiLocus effective length before it competes with FLG2 or the antisense nRNA spans. For a strand-compatible FLG2 read, RNA also gets a strand-likelihood advantage over gDNA: the RNA component gets probability near 1 on its matching orientation, while gDNA pays a 0.5 strand-symmetry factor. The opposite read orientation is not evidence against RNA globally because a positive-strand nRNA/antisense component can absorb it.

This means the failure is both a likelihood-model problem and a prior/geometry problem:

1. Likelihood problem: there is no coupled local test that says `positive and negative unspliced reads at the same exonic interval should be explained by one symmetric gDNA source unless the positive-strand transcript has independent support elsewhere`.
2. Prior/effective-length problem: the gDNA component is normalized over the whole chromosome-scale MultiLocus, not over the local captured opportunity that generated these reads.
3. Model-structure problem: nRNA components are allowed to claim intronic overlap locally without a penalty for absent support in their other exons/introns.

## What We Need

The direct fix is not just another scalar kappa tweak. We need local evidence coupling in the model. The most actionable design is a local gDNA-versus-ambiguous-unspliced guard:

1. For an exon-contained unspliced pileup, aggregate both orientations over a small genomic segment before hard assignment or before EM prior construction.
2. Require independent RNA support for the antisense/nRNA explanation: splice junctions, transcript-unique exons, or reads in other introns/exons of the same nascent span.
3. Give gDNA a local captured-opportunity effective length for that segment, not the entire chromosome-scale MultiLocus denominator.
4. Penalize nRNA components whose only support is a local overlap with a high-depth opposite-strand mature exon and whose remaining span has near-zero support.

Artifacts are in `{args.out_dir}`.
"""
    args.report.write_text(text)


def main() -> int:
    args = parse_args()
    context = load_context(args)
    features = build_feature_intervals(args, context)
    scan = scan_bam(args, features)
    write_outputs(args, context, features, scan)
    print(f"counted {scan['records']['window_fragments']:,} fragments in window")
    print(f"wrote {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())