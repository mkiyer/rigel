#!/usr/bin/env python3
"""Deep per-fragment analysis of selected no-MM VCaP gDNA false-RNA regions."""

from __future__ import annotations

import argparse
import csv
import json
import math
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
LOG_HALF = math.log(0.5)
OVERHANG_LOG_PENALTY = math.log(0.1)
EPS = 1e-12


@dataclass(frozen=True)
class RegionSpec:
    region_id: str
    chrom: str
    start: int
    end: int
    label: str
    reason: str


DEFAULT_REGIONS = [
    RegionSpec(
        "chr11_malat1_65500k_65510k",
        "chr11",
        65_500_000,
        65_510_000,
        "MALAT1 / TALAM1 neighborhood",
        "Top no-MM 10kb window: 3,822 gDNA->RNA calls, 90.9% local false-RNA rate.",
    ),
    RegionSpec(
        "chrx_ar_67540k_67550k",
        "chrX",
        67_540_000,
        67_550_000,
        "AR exon-rich window",
        "Second no-MM 10kb window: 3,270 gDNA->RNA calls, mostly mRNA AR.",
    ),
    RegionSpec(
        "chr2_long_nrna_178530k_178620k",
        "chr2",
        178_530_000,
        178_620_000,
        "Long nRNA span around chr2:178.53-178.62 Mb",
        "Contiguous severe nRNA false-positive band with little local RNA support.",
    ),
    RegionSpec(
        "chr1_flg2_152350k_152360k",
        "chr1",
        152_350_000,
        152_360_000,
        "FLG2 / FLG repeat-like exon window",
        "High false-RNA rate with nearly no RNA-source fragments in the window.",
    ),
    RegionSpec(
        "chr9_znf462_106920k_106930k",
        "chr9",
        106_920_000,
        106_930_000,
        "ZNF462 expressed exon window",
        "High mRNA false-positive count in a locally expressed gene.",
    ),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--quant-dir", type=Path, required=True)
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sample-per-region", type=int, default=12)
    return parser.parse_args()


def flowcell_from_qname(query_name: str) -> str:
    parts = query_name.split(":")
    if len(parts) >= 3:
        return parts[2]
    return query_name


def source_map_from_args(args: argparse.Namespace) -> dict[str, str]:
    source_map: dict[str, str] = {}
    for flowcell in args.rna_flowcell:
        source_map[flowcell] = "rna"
    for flowcell in args.gdna_flowcell:
        source_map[flowcell] = "gdna"
    return source_map


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def predicted_pool_from_zf(zf_value: int | None) -> str:
    if zf_value is None:
        return "unresolved"
    if zf_value & AF_GDNA_BIT:
        return "gdna"
    if zf_value & AF_NRNA_BIT:
        return "nrna"
    if zf_value & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def binary_pool(pool: str) -> str:
    if pool in {"mrna", "nrna"}:
        return "rna"
    return pool


def primary_records(records: list[pysam.AlignedSegment]) -> list[pysam.AlignedSegment]:
    return [record for record in records if not record.is_secondary and not record.is_supplementary]


def representative_record(records: list[pysam.AlignedSegment]) -> pysam.AlignedSegment | None:
    primaries = primary_records(records)
    if not primaries:
        return None
    for record in primaries:
        if record.is_paired and record.is_read1:
            return record
    return primaries[0]


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


def aligned_bases(record: pysam.AlignedSegment) -> int:
    total = 0
    for op, length in record.cigartuples or []:
        if op in {0, 2, 7, 8}:
            total += int(length)
    return total


def fragment_aligned_bases(records: Iterable[pysam.AlignedSegment]) -> int:
    return sum(aligned_bases(record) for record in primary_records(list(records)))


def overlapping_region(chrom: str, span_start: int, regions: list[RegionSpec]) -> RegionSpec | None:
    for region in regions:
        if region.chrom == chrom and region.start <= span_start < region.end:
            return region
    return None


def parse_nrna_target(target: str) -> tuple[str, int, int, str] | None:
    if not target.startswith("RIGEL_NRNA_"):
        return None
    try:
        prefix, strand_code, start_s, end_s = target.rsplit("_", 3)
    except ValueError:
        return None
    chrom = prefix.removeprefix("RIGEL_NRNA_")
    strand = {"1": "+", "2": "-"}.get(strand_code, strand_code)
    return chrom, int(start_s), int(end_s), strand


def format_region(chrom: str, start0: int, end: int) -> str:
    return f"{chrom}:{start0 + 1:,}-{end:,}"


def fmt_int(value: int) -> str:
    return f"{value:,}"


def pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * value:.{digits}f}%"


def fmt_float(value: float | None, digits: int = 3) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(records: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not records:
        return "_No rows._"
    lines = ["| " + " | ".join(label for label, _ in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for record in records:
        lines.append("| " + " | ".join(str(record.get(key, "")) for _, key in columns) + " |")
    return "\n".join(lines)


def iter_query_groups(bam: pysam.AlignmentFile) -> Iterable[tuple[str, list[pysam.AlignedSegment]]]:
    current_name = None
    current_records: list[pysam.AlignedSegment] = []
    for record in bam.fetch(until_eof=True):
        if current_name is None:
            current_name = record.query_name
        if record.query_name != current_name:
            yield current_name, current_records
            current_name = record.query_name
            current_records = [record]
        else:
            current_records.append(record)
    if current_name is not None:
        yield current_name, current_records


def init_region_stats() -> dict[str, Any]:
    return {
        "source_pred": Counter(),
        "false_categories": Counter(),
        "false_targets": Counter(),
        "false_zl": Counter(),
        "false_zw_sum": 0.0,
        "false_zw_count": 0,
        "false_zw_ge_90": 0,
        "false_frag_len_sum": 0,
        "false_frag_len_count": 0,
    }


def sample_row_from_group(
    region: RegionSpec,
    qname: str,
    records: list[pysam.AlignedSegment],
    read: pysam.AlignedSegment,
    span_start: int,
    span_end: int,
    true_source: str,
    predicted_detail: str,
) -> dict[str, Any]:
    primaries = primary_records(records)
    read1 = next((record for record in primaries if record.is_read1), read)
    read2 = next((record for record in primaries if record.is_read2), None)
    zt = str(get_tag(read, "ZT", "."))
    zg = str(get_tag(read, "ZG", "."))
    zr = str(get_tag(read, "ZR", "."))
    return {
        "region_id": region.region_id,
        "qname": qname,
        "true_source": true_source,
        "predicted_detail": predicted_detail,
        "chrom": read.reference_name,
        "fragment_start_1based": span_start + 1,
        "fragment_end": span_end,
        "fragment_length": span_end - span_start,
        "read1_start_1based": read1.reference_start + 1 if read1.reference_start >= 0 else -1,
        "read1_end": read1.reference_end or -1,
        "read1_cigar": read1.cigarstring or "",
        "read1_reverse": read1.is_reverse,
        "read2_start_1based": read2.reference_start + 1 if read2 and read2.reference_start >= 0 else -1,
        "read2_end": (read2.reference_end or -1) if read2 else -1,
        "read2_cigar": (read2.cigarstring or "") if read2 else "",
        "read2_reverse": read2.is_reverse if read2 else "",
        "template_length": read.template_length,
        "mapq": read.mapping_quality,
        "nh": get_tag(read, "NH", ""),
        "nm": get_tag(read, "NM", ""),
        "as": get_tag(read, "AS", ""),
        "md": get_tag(read, "MD", ""),
        "zf": get_tag(read, "ZF", ""),
        "zw": get_tag(read, "ZW", ""),
        "zc": get_tag(read, "ZC", ""),
        "zs": get_tag(read, "ZS", ""),
        "zn": get_tag(read, "ZN", ""),
        "zl": get_tag(read, "ZL", ""),
        "zi": get_tag(read, "ZI", ""),
        "zt": zt,
        "zg": zg,
        "zr": zr,
        "aligned_bases_total": fragment_aligned_bases(records),
    }


def scan_bam(args: argparse.Namespace, regions: list[RegionSpec]) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    record_counts = Counter()
    region_stats = {region.region_id: init_region_stats() for region in regions}
    samples_by_region: dict[str, list[dict[str, Any]]] = defaultdict(list)
    fl_counts = {"rna": np.ones(1001, dtype=np.float64), "gdna": np.ones(1001, dtype=np.float64)}

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for qname, records in iter_query_groups(bam):
            record_counts["query_groups_seen"] += 1
            read = representative_record(records)
            if read is None:
                record_counts["groups_without_primary"] += 1
                continue
            if read.is_paired and not read.is_read1:
                record_counts["groups_without_read1"] += 1
            record_counts["fragments_counted"] += 1

            true_source = source_map.get(flowcell_from_qname(qname), "unknown")
            zf_value = get_tag(read, "ZF", None)
            predicted_detail = predicted_pool_from_zf(zf_value)
            predicted_binary = binary_pool(predicted_detail)
            chrom = read.reference_name or "*"
            span_start, span_end = fragment_span(read)
            frag_len = span_end - span_start

            if true_source in fl_counts and 1 <= frag_len <= 1000:
                nh = get_tag(read, "NH", 1)
                if nh == 1:
                    fl_counts[true_source][frag_len] += 1.0

            region = overlapping_region(chrom, span_start, regions)
            if region is None:
                continue

            stats = region_stats[region.region_id]
            stats["source_pred"][(true_source, predicted_binary, predicted_detail)] += 1
            if true_source == "gdna" and predicted_binary == "rna":
                record_counts["selected_gdna_false_rna"] += 1
                zc = str(get_tag(read, "ZC", ""))
                zs = str(get_tag(read, "ZS", ""))
                zt = str(get_tag(read, "ZT", "."))
                zr = str(get_tag(read, "ZR", "."))
                target = zt if predicted_detail == "nrna" else f"{zr}:{zt}"
                zl = get_tag(read, "ZL", -1)
                zw = get_tag(read, "ZW", None)
                stats["false_categories"][(predicted_detail, zc, zs)] += 1
                stats["false_targets"][(predicted_detail, target)] += 1
                stats["false_zl"][zl] += 1
                if isinstance(zw, float):
                    stats["false_zw_sum"] += zw
                    stats["false_zw_count"] += 1
                    if zw >= 0.9:
                        stats["false_zw_ge_90"] += 1
                stats["false_frag_len_sum"] += frag_len
                stats["false_frag_len_count"] += 1
                sample = sample_row_from_group(
                    region, qname, records, read, span_start, span_end, true_source, predicted_detail
                )
                samples_by_region[region.region_id].append(sample)

    for source in fl_counts:
        fl_counts[source] /= fl_counts[source].sum()

    return {
        "record_counts": record_counts,
        "region_stats": region_stats,
        "samples_by_region": samples_by_region,
        "fl_probs": fl_counts,
        "source_map": source_map,
    }


def load_context(args: argparse.Namespace) -> dict[str, Any]:
    quant = pd.read_feather(args.quant_dir / "quant.feather")
    nrna_quant = pd.read_feather(args.quant_dir / "nrna_quant.feather")
    loci = pd.read_feather(args.quant_dir / "loci.feather")
    tx = pd.read_feather(args.index_dir / "transcripts.feather")
    intervals = pd.read_feather(args.index_dir / "intervals.feather")
    exon_intervals = intervals[intervals["interval_type"] == 0]
    exon_by_t: dict[int, list[tuple[int, int]]] = {}
    for t_index, group in exon_intervals.groupby("t_index", sort=False):
        exon_by_t[int(t_index)] = list(zip(group["start"].astype(int), group["end"].astype(int)))
    return {
        "quant_by_tid": {str(row.transcript_id): row for row in quant.itertuples(index=False)},
        "nrna_by_id": {str(row.nrna_id): row for row in nrna_quant.itertuples(index=False)},
        "loci_by_id": {int(row.locus_id): row for row in loci.itertuples(index=False)},
        "tx_by_index": {int(row.t_index): row for row in tx.itertuples(index=False)},
        "exon_by_t": exon_by_t,
    }


def cigar_blocks(start0: int, cigar: str) -> list[tuple[int, int]]:
    if start0 < 0 or not cigar:
        return []
    blocks: list[tuple[int, int]] = []
    ref_pos = start0
    length_chars: list[str] = []
    for char in cigar:
        if char.isdigit():
            length_chars.append(char)
            continue
        if not length_chars:
            continue
        length = int("".join(length_chars))
        length_chars = []
        if char in {"M", "=", "X"}:
            blocks.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif char in {"D", "N"}:
            ref_pos += length
        elif char in {"S", "H", "I", "P"}:
            continue
    return blocks


def sample_blocks(sample: dict[str, Any]) -> list[tuple[int, int]]:
    blocks = cigar_blocks(int(sample["read1_start_1based"]) - 1, str(sample["read1_cigar"]))
    if int(sample["read2_start_1based"]) > 0:
        blocks.extend(cigar_blocks(int(sample["read2_start_1based"]) - 1, str(sample["read2_cigar"])))
    return blocks


def block_overlap(blocks: list[tuple[int, int]], exons: list[tuple[int, int]]) -> int:
    if not exons:
        return 0
    total = 0
    for block_start, block_end in blocks:
        for exon_start, exon_end in exons:
            if exon_end <= block_start:
                continue
            if exon_start >= block_end:
                break
            total += max(0, min(block_end, exon_end) - max(block_start, exon_start))
    return total


def target_strand(sample: dict[str, Any], context: dict[str, Any]) -> str | None:
    zt = str(sample["zt"])
    parsed = parse_nrna_target(zt)
    if parsed is not None:
        return parsed[3]
    try:
        zi = int(sample["zi"])
    except (TypeError, ValueError):
        return None
    tx_row = context["tx_by_index"].get(zi)
    if tx_row is None:
        return None
    return "+" if int(tx_row.strand) == 1 else "-"


def protocol_consistent(read1_reverse: bool, strand: str | None) -> bool | None:
    if strand is None:
        return None
    if strand == "+":
        return bool(read1_reverse)
    if strand == "-":
        return not bool(read1_reverse)
    return None


def sigmoid(value: float) -> float:
    if value >= 0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


def enrich_samples(
    scan: dict[str, Any],
    context: dict[str, Any],
    regions: list[RegionSpec],
) -> list[dict[str, Any]]:
    samples = []
    fl_probs = scan["fl_probs"]
    log_p_sense = math.log(0.000751)
    log_p_antisense = math.log(0.999249)

    for region in regions:
        for sample in sorted(
            scan["samples_by_region"].get(region.region_id, []),
            key=lambda row: float(row.get("zw") or 0.0),
            reverse=True,
        ):
            frag_len = int(sample["fragment_length"])
            fl_idx = min(max(frag_len, 0), 1000)
            fl_log_ratio = math.log(fl_probs["rna"][fl_idx]) - math.log(fl_probs["gdna"][fl_idx])
            strand = target_strand(sample, context)
            consistent = protocol_consistent(bool(sample["read1_reverse"]), strand)
            if consistent is None:
                strand_log_ratio = None
            else:
                strand_log_ratio = (log_p_antisense if consistent else log_p_sense) - LOG_HALF

            predicted_detail = str(sample["predicted_detail"])
            zt = str(sample["zt"])
            try:
                zi = int(sample["zi"])
            except (TypeError, ValueError):
                zi = -1

            aligned_total = int(sample["aligned_bases_total"])
            exon_overlap = aligned_total
            overhang_bp = 0
            target_count = None
            target_eff_len = None
            target_locus = None
            if predicted_detail == "mrna":
                exons = context["exon_by_t"].get(zi, [])
                exon_overlap = block_overlap(sample_blocks(sample), exons)
                qrow = context["quant_by_tid"].get(zt)
                if qrow is not None:
                    target_count = float(qrow.count)
                    target_eff_len = float(qrow.effective_length)
                    target_locus = int(qrow.locus_id)
            elif predicted_detail == "nrna":
                nrow = context["nrna_by_id"].get(zt)
                if nrow is not None:
                    target_count = float(nrow.count)
                    target_eff_len = float(nrow.effective_length)
                    target_locus = int(nrow.locus_id)
                else:
                    qrow = context["quant_by_tid"].get(zt)
                    if qrow is not None:
                        target_count = float(qrow.count)
                        target_eff_len = float(qrow.effective_length)
                        target_locus = int(qrow.locus_id)

            if predicted_detail == "mrna":
                overhang_bp = max(0, aligned_total - exon_overlap)
            else:
                overhang_bp = 0
            overhang_log_ratio = overhang_bp * OVERHANG_LOG_PENALTY

            gdna_count = None
            gdna_eff_len = None
            abundance_log_ratio = None
            two_pool_log_odds = None
            two_pool_p_rna = None
            locus_row = context["loci_by_id"].get(target_locus) if target_locus is not None else None
            if locus_row is not None:
                gdna_count = float(locus_row.gdna)
                gdna_eff_len = float(locus_row.gdna_eff_len)
                if target_count is not None and target_eff_len is not None:
                    abundance_log_ratio = math.log((target_count + EPS) / target_eff_len) - math.log(
                        (gdna_count + EPS) / gdna_eff_len
                    )
            if abundance_log_ratio is not None and strand_log_ratio is not None:
                two_pool_log_odds = (
                    fl_log_ratio + strand_log_ratio + overhang_log_ratio + abundance_log_ratio
                )
                two_pool_p_rna = sigmoid(two_pool_log_odds)

            sample.update(
                {
                    "target_strand": strand or "",
                    "protocol_consistent": consistent if consistent is not None else "",
                    "fl_log_rna_over_gdna_empirical": fl_log_ratio,
                    "strand_log_rna_over_gdna": strand_log_ratio,
                    "exon_overlap_bp": exon_overlap,
                    "overhang_bp": overhang_bp,
                    "overhang_log_rna_over_gdna": overhang_log_ratio,
                    "target_count": target_count,
                    "target_eff_len": target_eff_len,
                    "target_locus_id": target_locus,
                    "gdna_count_in_target_locus": gdna_count,
                    "gdna_eff_len_in_target_locus": gdna_eff_len,
                    "abundance_log_rna_over_gdna": abundance_log_ratio,
                    "two_pool_log_odds_rna_over_gdna_approx": two_pool_log_odds,
                    "two_pool_p_rna_approx": two_pool_p_rna,
                }
            )
            samples.append(sample)
    return samples


def build_region_rows(scan: dict[str, Any], regions: list[RegionSpec]) -> list[dict[str, Any]]:
    rows = []
    for region in regions:
        stats = scan["region_stats"][region.region_id]
        source_pred = stats["source_pred"]
        gdna_total = sum(count for (source, _binary, _detail), count in source_pred.items() if source == "gdna")
        rna_total = sum(count for (source, _binary, _detail), count in source_pred.items() if source == "rna")
        gdna_false_rna = sum(
            count for (source, binary, _detail), count in source_pred.items() if source == "gdna" and binary == "rna"
        )
        gdna_pred_gdna = sum(
            count for (source, binary, _detail), count in source_pred.items() if source == "gdna" and binary == "gdna"
        )
        gdna_pred_unresolved = sum(
            count
            for (source, binary, _detail), count in source_pred.items()
            if source == "gdna" and binary == "unresolved"
        )
        pred_mrna = sum(
            count
            for (source, _binary, detail), count in source_pred.items()
            if source == "gdna" and detail == "mrna"
        )
        pred_nrna = sum(
            count
            for (source, _binary, detail), count in source_pred.items()
            if source == "gdna" and detail == "nrna"
        )
        zw_mean = stats["false_zw_sum"] / stats["false_zw_count"] if stats["false_zw_count"] else None
        frag_len_mean = (
            stats["false_frag_len_sum"] / stats["false_frag_len_count"]
            if stats["false_frag_len_count"]
            else None
        )
        rows.append(
            {
                "region_id": region.region_id,
                "label": region.label,
                "browser_region": format_region(region.chrom, region.start, region.end),
                "chrom": region.chrom,
                "start_0based": region.start,
                "end": region.end,
                "reason": region.reason,
                "gdna_source_total": gdna_total,
                "rna_source_total": rna_total,
                "gdna_false_rna": gdna_false_rna,
                "gdna_false_rna_rate": gdna_false_rna / gdna_total if gdna_total else 0.0,
                "gdna_pred_mrna": pred_mrna,
                "gdna_pred_nrna": pred_nrna,
                "gdna_pred_gdna": gdna_pred_gdna,
                "gdna_pred_unresolved": gdna_pred_unresolved,
                "false_zw_mean": zw_mean,
                "false_zw_ge_90_rate": stats["false_zw_ge_90"] / stats["false_zw_count"] if stats["false_zw_count"] else 0.0,
                "false_fragment_length_mean": frag_len_mean,
                "top_category": "; ".join(
                    f"{detail}/{zc}/{zs}:{count}"
                    for (detail, zc, zs), count in stats["false_categories"].most_common(3)
                ),
                "top_target": "; ".join(
                    f"{detail}:{target}:{count}"
                    for (detail, target), count in stats["false_targets"].most_common(3)
                ),
            }
        )
    return rows


def build_counter_rows(scan: dict[str, Any], regions: list[RegionSpec], key: str) -> list[dict[str, Any]]:
    rows = []
    for region in regions:
        counter = scan["region_stats"][region.region_id][key]
        for parts, count in counter.most_common():
            row = {"region_id": region.region_id, "count": count}
            if key == "false_categories":
                detail, zc, zs = parts
                row.update({"predicted_detail": detail, "zc": zc, "zs": zs})
            elif key == "false_targets":
                detail, target = parts
                row.update({"predicted_detail": detail, "target": target})
            elif key == "false_zl":
                row.update({"locus_id": parts})
            rows.append(row)
    return rows


def median_of(values: list[float]) -> float | None:
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return None
    return float(np.median(np.asarray(finite_values, dtype=np.float64)))


def build_likelihood_summary_rows(samples: list[dict[str, Any]], regions: list[RegionSpec]) -> list[dict[str, Any]]:
    rows = []
    for region in regions:
        region_samples = [sample for sample in samples if sample["region_id"] == region.region_id]
        for label, subset in [
            ("all_false_rna", region_samples),
            (
                "ambiguous_unspliced",
                [
                    sample
                    for sample in region_samples
                    if str(sample.get("zc", "")).startswith("ambig") and sample.get("zs") == "unspliced"
                ],
            ),
        ]:
            if not subset:
                continue
            p2_values = [float(sample["two_pool_p_rna_approx"]) for sample in subset if sample["two_pool_p_rna_approx"] is not None]
            rows.append(
                {
                    "region_id": region.region_id,
                    "subset": label,
                    "n": len(subset),
                    "median_zw": median_of([float(sample["zw"]) for sample in subset]),
                    "median_fl_log_rna_over_gdna": median_of(
                        [float(sample["fl_log_rna_over_gdna_empirical"]) for sample in subset]
                    ),
                    "median_strand_log_rna_over_gdna": median_of(
                        [
                            float(sample["strand_log_rna_over_gdna"])
                            for sample in subset
                            if sample["strand_log_rna_over_gdna"] is not None
                        ]
                    ),
                    "median_abundance_log_rna_over_gdna": median_of(
                        [
                            float(sample["abundance_log_rna_over_gdna"])
                            for sample in subset
                            if sample["abundance_log_rna_over_gdna"] is not None
                        ]
                    ),
                    "median_two_pool_p_rna_approx": median_of(p2_values),
                    "p2_ge_90_rate": sum(value >= 0.90 for value in p2_values) / len(p2_values)
                    if p2_values
                    else None,
                }
            )
    return rows


def report_region_rows(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "region": row["browser_region"],
            "label": row["label"],
            "gdna_total": fmt_int(int(row["gdna_source_total"])),
            "false": fmt_int(int(row["gdna_false_rna"])),
            "rate": pct(float(row["gdna_false_rna_rate"])),
            "mrna": fmt_int(int(row["gdna_pred_mrna"])),
            "nrna": fmt_int(int(row["gdna_pred_nrna"])),
            "rna_total": fmt_int(int(row["rna_source_total"])),
            "zw": fmt_float(row["false_zw_mean"]),
            "zw90": pct(float(row["false_zw_ge_90_rate"])),
        }
        for row in rows
    ]


def report_sample_rows(samples: list[dict[str, Any]], region_id: str, limit: int) -> list[dict[str, str]]:
    out = []
    region_samples = [row for row in samples if row["region_id"] == region_id]
    dominant = [
        row
        for row in region_samples
        if str(row.get("zc", "")).startswith("ambig") and row.get("zs") == "unspliced"
    ]
    selected = sorted(dominant or region_samples, key=lambda row: float(row.get("zw") or 0.0), reverse=True)[
        :limit
    ]
    for sample in selected:
        out.append(
            {
                "qname": sample["qname"],
                "span": format_region(sample["chrom"], int(sample["fragment_start_1based"]) - 1, int(sample["fragment_end"])),
                "pred": sample["predicted_detail"],
                "target": sample["zr"] if sample["predicted_detail"] == "mrna" else sample["zt"],
                "tlen": sample["template_length"],
                "cigars": f"{sample['read1_cigar']}/{sample['read2_cigar']}",
                "nh": sample["nh"],
                "nm": sample["nm"],
                "zc": sample["zc"],
                "zs": sample["zs"],
                "zw": fmt_float(float(sample["zw"])),
                "fl_lr": fmt_float(sample["fl_log_rna_over_gdna_empirical"]),
                "strand_lr": fmt_float(sample["strand_log_rna_over_gdna"]),
                "abund_lr": fmt_float(sample["abundance_log_rna_over_gdna"]),
                "p2": pct(sample["two_pool_p_rna_approx"]),
            }
        )
    return out


def write_outputs(
    args: argparse.Namespace,
    regions: list[RegionSpec],
    scan: dict[str, Any],
    region_rows: list[dict[str, Any]],
    category_rows: list[dict[str, Any]],
    target_rows: list[dict[str, Any]],
    locus_rows: list[dict[str, Any]],
    likelihood_rows: list[dict[str, Any]],
    sample_rows: list[dict[str, Any]],
) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "selected_regions_summary.tsv",
        region_rows,
        [
            "region_id",
            "label",
            "browser_region",
            "chrom",
            "start_0based",
            "end",
            "reason",
            "gdna_source_total",
            "rna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
            "gdna_pred_mrna",
            "gdna_pred_nrna",
            "gdna_pred_gdna",
            "gdna_pred_unresolved",
            "false_zw_mean",
            "false_zw_ge_90_rate",
            "false_fragment_length_mean",
            "top_category",
            "top_target",
        ],
    )
    write_tsv(
        args.out_dir / "region_false_categories.tsv",
        category_rows,
        ["region_id", "predicted_detail", "zc", "zs", "count"],
    )
    write_tsv(
        args.out_dir / "region_false_targets.tsv",
        target_rows,
        ["region_id", "predicted_detail", "target", "count"],
    )
    write_tsv(args.out_dir / "region_false_loci.tsv", locus_rows, ["region_id", "locus_id", "count"])
    write_tsv(
        args.out_dir / "likelihood_term_summary.tsv",
        likelihood_rows,
        [
            "region_id",
            "subset",
            "n",
            "median_zw",
            "median_fl_log_rna_over_gdna",
            "median_strand_log_rna_over_gdna",
            "median_abundance_log_rna_over_gdna",
            "median_two_pool_p_rna_approx",
            "p2_ge_90_rate",
        ],
    )
    sample_fields = [
        "region_id",
        "qname",
        "true_source",
        "predicted_detail",
        "chrom",
        "fragment_start_1based",
        "fragment_end",
        "fragment_length",
        "read1_start_1based",
        "read1_end",
        "read1_cigar",
        "read1_reverse",
        "read2_start_1based",
        "read2_end",
        "read2_cigar",
        "read2_reverse",
        "template_length",
        "mapq",
        "nh",
        "nm",
        "as",
        "md",
        "zf",
        "zw",
        "zc",
        "zs",
        "zn",
        "zl",
        "zi",
        "zt",
        "zg",
        "zr",
        "aligned_bases_total",
        "target_strand",
        "protocol_consistent",
        "fl_log_rna_over_gdna_empirical",
        "strand_log_rna_over_gdna",
        "exon_overlap_bp",
        "overhang_bp",
        "overhang_log_rna_over_gdna",
        "target_count",
        "target_eff_len",
        "target_locus_id",
        "gdna_count_in_target_locus",
        "gdna_eff_len_in_target_locus",
        "abundance_log_rna_over_gdna",
        "two_pool_log_odds_rna_over_gdna_approx",
        "two_pool_p_rna_approx",
    ]
    write_tsv(args.out_dir / "per_fragment_examples.tsv", sample_rows, sample_fields)
    with (args.out_dir / "selected_regions.bed").open("w") as handle:
        for region in regions:
            handle.write(f"{region.chrom}\t{region.start}\t{region.end}\t{region.region_id}\n")
    summary = {
        "bam": str(args.bam),
        "quant_dir": str(args.quant_dir),
        "index_dir": str(args.index_dir),
        "record_counts": dict(scan["record_counts"]),
        "regions": [region.__dict__ for region in regions],
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    write_report(args, regions, region_rows, likelihood_rows, sample_rows)


def write_report(
    args: argparse.Namespace,
    regions: list[RegionSpec],
    region_rows: list[dict[str, Any]],
    likelihood_rows: list[dict[str, Any]],
    sample_rows: list[dict[str, Any]],
) -> None:
    text = f"""# VCaP No-MM gDNA -> RNA Deep Region Analysis - 2026-05-18

BAM: `{args.bam}`

This report focuses on five no-multimap regions with severe gDNA-source fragments miscalled as RNA. Coordinates are 1-based closed for browser use in tables; BED artifacts are 0-based half-open.

Important limitation: the annotated BAM stores the winner (`ZT`/`ZF`) and winner posterior (`ZW`), but not the full candidate log-likelihood vector or runner-up gDNA posterior. The per-fragment likelihood table therefore reconstructs the decision drivers from available run outputs: empirical RNA/gDNA fragment-length likelihood ratio from this BAM, strand likelihood ratio from the run's R1-antisense model, target abundance/effective-length versus locus gDNA abundance/effective-length, and the observed `ZW`. This is enough to identify the mechanism, but adding runner-up gDNA posterior tags would make this exact.

## Selected Regions

{markdown_table(report_region_rows(region_rows), [
        ('Region', 'region'),
        ('Label', 'label'),
        ('gDNA fragments', 'gdna_total'),
        ('gDNA->RNA', 'false'),
        ('False rate', 'rate'),
        ('mRNA', 'mrna'),
        ('nRNA', 'nrna'),
        ('RNA fragments', 'rna_total'),
        ('Mean ZW', 'zw'),
        ('ZW >= .90', 'zw90'),
    ])}

## Per-Fragment Examples

The example rows below are chosen from the dominant ambiguous-unspliced false-positive category for each region when available. Columns `fl lr`, `strand lr`, and `abund lr` are log RNA-over-gDNA terms. `p2` is a two-pool RNA probability reconstructed from those terms and should be read as diagnostic, not as the exact Rigel posterior.

## Likelihood-Term Summary

{markdown_table([
        {
            'region': next(region.label for region in regions if region.region_id == row['region_id']),
            'subset': row['subset'],
            'n': fmt_int(int(row['n'])),
            'zw': fmt_float(row['median_zw']),
            'fl': fmt_float(row['median_fl_log_rna_over_gdna']),
            'strand': fmt_float(row['median_strand_log_rna_over_gdna']),
            'abund': fmt_float(row['median_abundance_log_rna_over_gdna']),
            'p2': pct(row['median_two_pool_p_rna_approx']),
            'p2_90': pct(row['p2_ge_90_rate']),
        }
        for row in likelihood_rows
        if row['subset'] == 'ambiguous_unspliced'
    ], [
        ('Region', 'region'),
        ('n', 'n'),
        ('Median ZW', 'zw'),
        ('Median fl lr', 'fl'),
        ('Median strand lr', 'strand'),
        ('Median abund lr', 'abund'),
        ('Median p2', 'p2'),
        ('p2 >= .90', 'p2_90'),
    ])}
"""
    for region in regions:
        text += f"""
### {region.label}

Browser region: `{format_region(region.chrom, region.start, region.end)}`

{region.reason}

{markdown_table(report_sample_rows(sample_rows, region.region_id, args.sample_per_region), [
        ('QNAME', 'qname'),
        ('Span', 'span'),
        ('Pred', 'pred'),
        ('Target', 'target'),
        ('TLEN', 'tlen'),
        ('CIGARs', 'cigars'),
        ('NH', 'nh'),
        ('NM', 'nm'),
        ('ZC', 'zc'),
        ('ZS', 'zs'),
        ('ZW', 'zw'),
        ('fl lr', 'fl_lr'),
        ('strand lr', 'strand_lr'),
        ('abund lr', 'abund_lr'),
        ('p2', 'p2'),
    ])}
"""

    text += f"""
## Mechanism

1. These fragments are mostly ordinary gDNA-compatible paired-end alignments: `NH=1`, contiguous `150M/150M`, low `NM`, and `ZS=unspliced`.
2. The RNA candidate often gets a strand likelihood advantage of about +0.69 log units over gDNA because the library is highly stranded while the gDNA likelihood is strand-symmetric (`log(0.5)`). This is appropriate for true RNA, but dangerous when the geometry is otherwise gDNA-compatible.
3. The fragment-length term is usually small compared with the abundance/effective-length term. Many false-positive fragments have lengths that are plausible under both RNA and gDNA models.
4. The decisive term is usually `abund lr`: the selected transcript or nRNA target has much larger abundance per effective length than the locus gDNA component. Once an ambiguous unspliced fragment has an RNA candidate, EM can make RNA the high-posterior winner even though the fragment has no splice-junction or RNA-unique evidence.
5. The low-RNA regions show self-seeding: gDNA fragments that overlap long nRNA spans can build enough nRNA mass locally to pull additional gDNA fragments into nRNA, even when RNA-source reads are scarce.

## Fix Direction

1. For diagnostic calling, treat `ZS=unspliced` plus `ZC=ambig_same_strand/ambig_opp_strand` as weak RNA evidence. Require a positive splice/unique-exon signal or a large RNA-over-gDNA margin before a hard RNA call.
2. Add a conservative hard-call policy that defaults ambiguous unspliced RNA/gDNA competition to gDNA unless runner-up gDNA posterior is safely low.
3. Emit runner-up pool posterior or `log_posterior_rna - log_posterior_gdna` in annotated BAMs. `ZW` tells us the winner can be confident, but not how close gDNA was.
4. Add local guards for long nRNA spans: if a region has little local RNA-source evidence and many gDNA-compatible unspliced reads, prevent nRNA self-seeding from becoming strong positive RNA evidence.

## Artifacts

- Region summary: `{args.out_dir / 'selected_regions_summary.tsv'}`
- Region BED: `{args.out_dir / 'selected_regions.bed'}`
- Per-fragment examples: `{args.out_dir / 'per_fragment_examples.tsv'}`
- Likelihood-term summary: `{args.out_dir / 'likelihood_term_summary.tsv'}`
- Category counts: `{args.out_dir / 'region_false_categories.tsv'}`
- Target counts: `{args.out_dir / 'region_false_targets.tsv'}`
- Locus counts: `{args.out_dir / 'region_false_loci.tsv'}`
"""
    args.report.write_text(text)


def main() -> None:
    args = parse_args()
    regions = DEFAULT_REGIONS
    scan = scan_bam(args, regions)
    context = load_context(args)
    samples = enrich_samples(scan, context, regions)
    region_rows = build_region_rows(scan, regions)
    category_rows = build_counter_rows(scan, regions, "false_categories")
    target_rows = build_counter_rows(scan, regions, "false_targets")
    locus_rows = build_counter_rows(scan, regions, "false_zl")
    likelihood_rows = build_likelihood_summary_rows(samples, regions)
    write_outputs(
        args, regions, scan, region_rows, category_rows, target_rows, locus_rows, likelihood_rows, samples
    )
    print(f"counted {scan['record_counts']['fragments_counted']:,} fragments")
    print(f"selected gDNA->RNA fragments: {scan['record_counts']['selected_gdna_false_rna']:,}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()