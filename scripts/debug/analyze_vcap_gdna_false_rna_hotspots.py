#!/usr/bin/env python3
"""Hotspot analysis for gDNA-source fragments miscalled as RNA by Rigel."""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True, help="Rigel annotated BAM to scan.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--report", type=Path, required=True, help="Markdown report path.")
    parser.add_argument(
        "--rna-flowcell",
        action="append",
        default=[],
        help="Flowcell ID from RNA source. May be repeated.",
    )
    parser.add_argument(
        "--gdna-flowcell",
        action="append",
        default=[],
        help="Flowcell ID from gDNA source. May be repeated.",
    )
    parser.add_argument("--window-size", type=int, default=10_000, help="Window size in bp.")
    parser.add_argument("--threads", type=int, default=4, help="BGZF decompression threads.")
    parser.add_argument("--sample-per-locus", type=int, default=4, help="Example fragments to keep per locus.")
    parser.add_argument("--top-n", type=int, default=30, help="Rows to show in report tables.")
    return parser.parse_args()


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


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
        previous = source_map.get(flowcell)
        if previous is not None and previous != "gdna":
            raise ValueError(f"Flowcell {flowcell} was assigned to both RNA and gDNA")
        source_map[flowcell] = "gdna"
    return source_map


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


def binary_pool(predicted_pool: str) -> str:
    if predicted_pool in {"mrna", "nrna"}:
        return "rna"
    return predicted_pool


def is_fragment_representative(read: pysam.AlignedSegment) -> bool:
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_paired:
        return read.is_read1
    return True


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


def parse_nrna_region(target: str) -> str:
    if not target.startswith("RIGEL_NRNA_"):
        return ""
    try:
        prefix, strand, start, end = target.rsplit("_", 3)
    except ValueError:
        return ""
    chrom = prefix.removeprefix("RIGEL_NRNA_")
    strand_label = {"1": "+", "2": "-"}.get(strand, strand)
    return f"{chrom}:{int(start):,}-{int(end):,}({strand_label})"


def target_label(read: pysam.AlignedSegment, predicted_detail: str) -> str:
    zt = str(get_tag(read, "ZT", "."))
    zg = str(get_tag(read, "ZG", "."))
    zr = str(get_tag(read, "ZR", "."))
    if predicted_detail == "nrna":
        region = parse_nrna_region(zt)
        if region:
            return sys.intern(f"nRNA {region}")
        return sys.intern(f"nRNA {zt}")
    if predicted_detail == "mrna":
        gene = zr if zr not in {"", "."} else zg
        return sys.intern(f"mRNA {gene} {zt}")
    return sys.intern(predicted_detail)


def compact_target(target: str, max_len: int = 80) -> str:
    if len(target) <= max_len:
        return target
    return target[: max_len - 3] + "..."


class OnlineStats:
    __slots__ = ("count", "sum", "min", "max", "ge_50", "ge_90", "ge_99")

    def __init__(self) -> None:
        self.count = 0
        self.sum = 0.0
        self.min = float("inf")
        self.max = float("-inf")
        self.ge_50 = 0
        self.ge_90 = 0
        self.ge_99 = 0

    def add(self, value: float | None) -> None:
        if value is None:
            return
        self.count += 1
        self.sum += value
        self.min = min(self.min, value)
        self.max = max(self.max, value)
        if value >= 0.50:
            self.ge_50 += 1
        if value >= 0.90:
            self.ge_90 += 1
        if value >= 0.99:
            self.ge_99 += 1

    def mean(self) -> float | None:
        if self.count == 0:
            return None
        return self.sum / self.count


def pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * value:.{digits}f}%"


def fmt_int(value: int) -> str:
    return f"{value:,}"


def fmt_float(value: float | None, digits: int = 3) -> str:
    if value is None:
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
    header = "| " + " | ".join(label for label, _ in columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = []
    for record in records:
        body.append("| " + " | ".join(str(record.get(key, "")) for _, key in columns) + " |")
    return "\n".join([header, sep, *body])


def sample_from_read(
    read: pysam.AlignedSegment,
    true_source: str,
    predicted_detail: str,
    span_start: int,
    span_end: int,
    target: str,
) -> dict[str, Any]:
    return {
        "qname": read.query_name,
        "true_source": true_source,
        "predicted_detail": predicted_detail,
        "chrom": read.reference_name,
        "fragment_start_1based": span_start + 1,
        "fragment_end": span_end,
        "read_start_1based": (read.reference_start + 1) if read.reference_start >= 0 else -1,
        "read_end": read.reference_end or -1,
        "cigar": read.cigarstring or "",
        "template_length": read.template_length,
        "mapq": read.mapping_quality,
        "is_reverse": read.is_reverse,
        "mate_is_reverse": read.mate_is_reverse,
        "nh": get_tag(read, "NH", ""),
        "nm": get_tag(read, "NM", ""),
        "as": get_tag(read, "AS", ""),
        "zf": get_tag(read, "ZF", ""),
        "zw": get_tag(read, "ZW", ""),
        "zc": get_tag(read, "ZC", ""),
        "zs": get_tag(read, "ZS", ""),
        "zl": get_tag(read, "ZL", ""),
        "zt": get_tag(read, "ZT", ""),
        "zg": get_tag(read, "ZG", ""),
        "zr": get_tag(read, "ZR", ""),
        "target_label": target,
    }


def update_counter_span(
    spans: dict[Any, list[Any]],
    key: Any,
    chrom: str,
    span_start: int,
    span_end: int,
) -> None:
    if key not in spans:
        spans[key] = [chrom, span_start, span_end]
        return
    entry = spans[key]
    if entry[0] != chrom:
        entry[0] = "mixed"
    entry[1] = min(entry[1], span_start)
    entry[2] = max(entry[2], span_end)


def scan_bam(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    record_counts = Counter()
    global_counts = Counter()
    category_counts = Counter()
    window_counts = defaultdict(Counter)
    window_targets = defaultdict(Counter)
    window_categories = defaultdict(Counter)
    locus_counts = defaultdict(Counter)
    locus_spans: dict[Any, list[Any]] = {}
    locus_targets = defaultdict(Counter)
    locus_zw = defaultdict(OnlineStats)
    target_counts = defaultdict(Counter)
    target_spans: dict[Any, list[Any]] = {}
    target_zw = defaultdict(OnlineStats)
    samples_by_locus: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    high_conf_heap: list[tuple[float, int, dict[str, Any]]] = []
    sequence_id = 0

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            record_counts["records_seen"] += 1
            if read.is_secondary:
                record_counts["secondary_skipped"] += 1
                continue
            if read.is_supplementary:
                record_counts["supplementary_skipped"] += 1
                continue
            if read.is_paired and not read.is_read1:
                record_counts["read2_skipped"] += 1
                continue
            if not is_fragment_representative(read):
                record_counts["other_skipped"] += 1
                continue

            record_counts["fragments_counted"] += 1
            flowcell = flowcell_from_qname(read.query_name)
            true_source = source_map.get(flowcell, "unknown")
            zf_value = get_tag(read, "ZF", None)
            predicted_detail = predicted_pool_from_zf(zf_value)
            predicted_binary = binary_pool(predicted_detail)
            span_start, span_end = fragment_span(read)
            chrom = read.reference_name or "*"
            window_start = (span_start // args.window_size) * args.window_size
            window_key = (chrom, window_start)

            global_counts[(true_source, predicted_binary, predicted_detail)] += 1
            if true_source in {"gdna", "rna"}:
                window_counts[window_key][f"{true_source}_total"] += 1
                window_counts[window_key][f"{true_source}_pred_{predicted_binary}"] += 1

            if true_source != "gdna" or predicted_binary != "rna":
                continue

            record_counts["gdna_false_rna"] += 1
            zc_value = sys.intern(str(get_tag(read, "ZC", "missing")))
            zs_value = sys.intern(str(get_tag(read, "ZS", "missing")))
            zw_value = get_tag(read, "ZW", None)
            try:
                zw_float = float(zw_value)
            except (TypeError, ValueError):
                zw_float = None
            locus_id = get_tag(read, "ZL", -999999)
            locus_key = (chrom, locus_id)
            category_key = (predicted_detail, zc_value, zs_value)
            target = target_label(read, predicted_detail)

            category_counts[category_key] += 1
            window_counts[window_key]["gdna_false_rna"] += 1
            window_counts[window_key][f"false_{predicted_detail}"] += 1
            window_targets[window_key][target] += 1
            window_categories[window_key][category_key] += 1

            locus_counts[locus_key]["false_total"] += 1
            locus_counts[locus_key][f"false_{predicted_detail}"] += 1
            locus_counts[locus_key][f"zc_{zc_value}"] += 1
            locus_counts[locus_key][f"zs_{zs_value}"] += 1
            locus_targets[locus_key][target] += 1
            locus_zw[locus_key].add(zw_float)
            update_counter_span(locus_spans, locus_key, chrom, span_start, span_end)

            target_counts[target]["false_total"] += 1
            target_counts[target][f"false_{predicted_detail}"] += 1
            target_counts[target][f"zc_{zc_value}"] += 1
            target_counts[target][f"zs_{zs_value}"] += 1
            target_zw[target].add(zw_float)
            update_counter_span(target_spans, target, chrom, span_start, span_end)

            sample = sample_from_read(
                read, true_source, predicted_detail, span_start, span_end, target
            )
            sample["locus_key"] = f"{chrom}:{locus_id}"
            sample["locus_id"] = locus_id
            if len(samples_by_locus[locus_key]) < args.sample_per_locus:
                samples_by_locus[locus_key].append(sample)
            if zw_float is not None:
                sequence_id += 1
                heap_item = (zw_float, sequence_id, sample)
                if len(high_conf_heap) < 100:
                    heapq.heappush(high_conf_heap, heap_item)
                elif zw_float > high_conf_heap[0][0]:
                    heapq.heapreplace(high_conf_heap, heap_item)

    return {
        "source_map": source_map,
        "record_counts": record_counts,
        "global_counts": global_counts,
        "category_counts": category_counts,
        "window_counts": window_counts,
        "window_targets": window_targets,
        "window_categories": window_categories,
        "locus_counts": locus_counts,
        "locus_spans": locus_spans,
        "locus_targets": locus_targets,
        "locus_zw": locus_zw,
        "target_counts": target_counts,
        "target_spans": target_spans,
        "target_zw": target_zw,
        "samples_by_locus": samples_by_locus,
        "high_conf_samples": [item[2] for item in sorted(high_conf_heap, reverse=True)],
    }


def top_count(counter: Counter) -> tuple[Any, int]:
    if not counter:
        return "", 0
    return counter.most_common(1)[0]


def build_rows(scan: dict[str, Any], window_size: int) -> dict[str, list[dict[str, Any]]]:
    gdna_total = sum(
        count
        for (source, _binary, _detail), count in scan["global_counts"].items()
        if source == "gdna"
    )
    false_total = scan["record_counts"]["gdna_false_rna"]

    category_rows = []
    for (detail, zc_value, zs_value), count in scan["category_counts"].most_common():
        category_rows.append(
            {
                "predicted_detail": detail,
                "zc": zc_value,
                "zs": zs_value,
                "count": count,
                "fraction_of_false_rna": count / false_total if false_total else 0.0,
                "fraction_of_gdna_source": count / gdna_total if gdna_total else 0.0,
            }
        )

    window_rows = []
    for (chrom, window_start), counts in scan["window_counts"].items():
        false_count = counts.get("gdna_false_rna", 0)
        if false_count <= 0:
            continue
        gdna_source_total = counts.get("gdna_total", 0)
        rna_source_total = counts.get("rna_total", 0)
        top_target, top_target_count = top_count(scan["window_targets"].get((chrom, window_start), Counter()))
        top_category, top_category_count = top_count(
            scan["window_categories"].get((chrom, window_start), Counter())
        )
        if isinstance(top_category, tuple):
            top_category_label = "/".join(str(part) for part in top_category)
        else:
            top_category_label = ""
        window_rows.append(
            {
                "chrom": chrom,
                "window_start_1based": window_start + 1,
                "window_end": window_start + window_size,
                "gdna_source_total": gdna_source_total,
                "gdna_false_rna": false_count,
                "gdna_false_rna_rate": false_count / gdna_source_total if gdna_source_total else 0.0,
                "gdna_pred_mrna": counts.get("gdna_pred_rna", 0)
                - counts.get("false_nrna", 0),
                "gdna_pred_nrna": counts.get("false_nrna", 0),
                "gdna_pred_gdna": counts.get("gdna_pred_gdna", 0),
                "gdna_pred_unresolved": counts.get("gdna_pred_unresolved", 0),
                "rna_source_total": rna_source_total,
                "rna_pred_rna": counts.get("rna_pred_rna", 0),
                "top_target": top_target,
                "top_target_count": top_target_count,
                "top_category": top_category_label,
                "top_category_count": top_category_count,
            }
        )
    window_rows.sort(key=lambda row: (-row["gdna_false_rna"], -row["gdna_false_rna_rate"]))

    locus_entries = []
    for locus_key, counts in scan["locus_counts"].items():
        chrom, locus_id = locus_key
        span = scan["locus_spans"].get(locus_key, ["", 0, 0])
        zw = scan["locus_zw"][locus_key]
        top_target, top_target_count = top_count(scan["locus_targets"].get(locus_key, Counter()))
        false_count = counts["false_total"]
        row = (
            {
                "locus_key": f"{chrom}:{locus_id}",
                "locus_id": locus_id,
                "chrom": span[0],
                "observed_start_1based": span[1] + 1,
                "observed_end": span[2],
                "gdna_false_rna": false_count,
                "fraction_of_false_rna": false_count / false_total if false_total else 0.0,
                "pred_mrna": counts.get("false_mrna", 0),
                "pred_nrna": counts.get("false_nrna", 0),
                "ambig_same": counts.get("zc_ambig_same_strand", 0),
                "ambig_opp": counts.get("zc_ambig_opp_strand", 0),
                "unambig": counts.get("zc_unambig", 0),
                "unspliced": counts.get("zs_unspliced", 0),
                "mean_zw": zw.mean(),
                "max_zw": zw.max if zw.count else None,
                "zw_ge_90_rate": zw.ge_90 / zw.count if zw.count else 0.0,
                "top_target": top_target,
                "top_target_count": top_target_count,
            }
        )
        locus_entries.append((false_count, locus_key, row))
    locus_entries.sort(key=lambda entry: -entry[0])
    locus_rows = [entry[2] for entry in locus_entries]

    target_rows = []
    for target, counts in scan["target_counts"].items():
        span = scan["target_spans"].get(target, ["", 0, 0])
        zw = scan["target_zw"][target]
        false_count = counts["false_total"]
        target_rows.append(
            {
                "target_label": target,
                "chrom": span[0],
                "observed_start_1based": span[1] + 1,
                "observed_end": span[2],
                "gdna_false_rna": false_count,
                "fraction_of_false_rna": false_count / false_total if false_total else 0.0,
                "pred_mrna": counts.get("false_mrna", 0),
                "pred_nrna": counts.get("false_nrna", 0),
                "ambig_same": counts.get("zc_ambig_same_strand", 0),
                "ambig_opp": counts.get("zc_ambig_opp_strand", 0),
                "unspliced": counts.get("zs_unspliced", 0),
                "mean_zw": zw.mean(),
                "max_zw": zw.max if zw.count else None,
                "zw_ge_90_rate": zw.ge_90 / zw.count if zw.count else 0.0,
            }
        )
    target_rows.sort(key=lambda row: -row["gdna_false_rna"])

    sample_rows = []
    top_locus_keys = [entry[1] for entry in locus_entries[:20]]
    for locus_key in top_locus_keys:
        for sample in scan["samples_by_locus"].get(locus_key, []):
            row = {**sample}
            sample_rows.append(row)

    high_conf_rows = []
    for sample in scan["high_conf_samples"][:50]:
        high_conf_rows.append(sample)

    return {
        "category_rows": category_rows,
        "window_rows": window_rows,
        "locus_rows": locus_rows,
        "target_rows": target_rows,
        "sample_rows": sample_rows,
        "high_conf_rows": high_conf_rows,
    }


def write_outputs(args: argparse.Namespace, scan: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    write_tsv(
        args.out_dir / "gdna_false_rna_categories.tsv",
        rows["category_rows"],
        [
            "predicted_detail",
            "zc",
            "zs",
            "count",
            "fraction_of_false_rna",
            "fraction_of_gdna_source",
        ],
    )
    write_tsv(
        args.out_dir / "gdna_false_rna_windows.tsv",
        rows["window_rows"],
        [
            "chrom",
            "window_start_1based",
            "window_end",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
            "gdna_pred_mrna",
            "gdna_pred_nrna",
            "gdna_pred_gdna",
            "gdna_pred_unresolved",
            "rna_source_total",
            "rna_pred_rna",
            "top_target",
            "top_target_count",
            "top_category",
            "top_category_count",
        ],
    )
    write_tsv(
        args.out_dir / "gdna_false_rna_loci.tsv",
        rows["locus_rows"],
        [
            "locus_id",
            "locus_key",
            "chrom",
            "observed_start_1based",
            "observed_end",
            "gdna_false_rna",
            "fraction_of_false_rna",
            "pred_mrna",
            "pred_nrna",
            "ambig_same",
            "ambig_opp",
            "unambig",
            "unspliced",
            "mean_zw",
            "max_zw",
            "zw_ge_90_rate",
            "top_target",
            "top_target_count",
        ],
    )
    write_tsv(
        args.out_dir / "gdna_false_rna_targets.tsv",
        rows["target_rows"],
        [
            "target_label",
            "chrom",
            "observed_start_1based",
            "observed_end",
            "gdna_false_rna",
            "fraction_of_false_rna",
            "pred_mrna",
            "pred_nrna",
            "ambig_same",
            "ambig_opp",
            "unspliced",
            "mean_zw",
            "max_zw",
            "zw_ge_90_rate",
        ],
    )
    sample_fields = [
        "locus_id",
        "locus_key",
        "qname",
        "true_source",
        "predicted_detail",
        "chrom",
        "fragment_start_1based",
        "fragment_end",
        "read_start_1based",
        "read_end",
        "cigar",
        "template_length",
        "mapq",
        "is_reverse",
        "mate_is_reverse",
        "nh",
        "nm",
        "as",
        "zf",
        "zw",
        "zc",
        "zs",
        "zl",
        "zt",
        "zg",
        "zr",
        "target_label",
    ]
    write_tsv(args.out_dir / "sample_false_rna_fragments.tsv", rows["sample_rows"], sample_fields)
    high_conf_fields = [field for field in sample_fields if field != "locus_id"]
    write_tsv(
        args.out_dir / "high_confidence_false_rna_fragments.tsv",
        rows["high_conf_rows"],
        high_conf_fields,
    )

    summary = {
        "bam": str(args.bam),
        "source_map": scan["source_map"],
        "record_counts": dict(scan["record_counts"]),
        "window_size": args.window_size,
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    write_report(args, scan, rows)


def report_category_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "detail": row["predicted_detail"],
                "zc": row["zc"],
                "zs": row["zs"],
                "count": fmt_int(int(row["count"])),
                "false_frac": pct(float(row["fraction_of_false_rna"])),
                "gdna_frac": pct(float(row["fraction_of_gdna_source"])),
            }
        )
    return out


def report_window_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "region": f"{row['chrom']}:{fmt_int(int(row['window_start_1based']))}-{fmt_int(int(row['window_end']))}",
                "false": fmt_int(int(row["gdna_false_rna"])),
                "rate": pct(float(row["gdna_false_rna_rate"])),
                "gdna_total": fmt_int(int(row["gdna_source_total"])),
                "mrna": fmt_int(int(row["gdna_pred_mrna"])),
                "nrna": fmt_int(int(row["gdna_pred_nrna"])),
                "top_target": compact_target(str(row["top_target"])),
                "top_cat": compact_target(str(row["top_category"]), 45),
            }
        )
    return out


def report_locus_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "locus": row["locus_key"],
                "region": f"{row['chrom']}:{fmt_int(int(row['observed_start_1based']))}-{fmt_int(int(row['observed_end']))}",
                "false": fmt_int(int(row["gdna_false_rna"])),
                "false_frac": pct(float(row["fraction_of_false_rna"])),
                "mrna": fmt_int(int(row["pred_mrna"])),
                "nrna": fmt_int(int(row["pred_nrna"])),
                "same": fmt_int(int(row["ambig_same"])),
                "opp": fmt_int(int(row["ambig_opp"])),
                "mean_zw": fmt_float(row["mean_zw"]),
                "zw90": pct(float(row["zw_ge_90_rate"])),
                "target": compact_target(str(row["top_target"])),
            }
        )
    return out


def report_local_rna_rows(rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    bins = [
        ("0", lambda value: value == 0),
        ("1-100", lambda value: 0 < value <= 100),
        ("101-1,000", lambda value: 100 < value <= 1000),
        (">1,000", lambda value: value > 1000),
    ]
    total_false = sum(int(row["gdna_false_rna"]) for row in rows)
    out = []
    for label, predicate in bins:
        selected = [row for row in rows if predicate(int(row["rna_source_total"]))]
        false_count = sum(int(row["gdna_false_rna"]) for row in selected)
        gdna_total = sum(int(row["gdna_source_total"]) for row in selected)
        out.append(
            {
                "rna_source_bin": label,
                "windows": fmt_int(len(selected)),
                "false_rna": fmt_int(false_count),
                "false_frac": pct(false_count / total_false if total_false else None),
                "weighted_false_rate": pct(false_count / gdna_total if gdna_total else None),
            }
        )
    return out


def report_target_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        out.append(
            {
                "target": compact_target(str(row["target_label"])),
                "region": f"{row['chrom']}:{fmt_int(int(row['observed_start_1based']))}-{fmt_int(int(row['observed_end']))}",
                "false": fmt_int(int(row["gdna_false_rna"])),
                "false_frac": pct(float(row["fraction_of_false_rna"])),
                "mrna": fmt_int(int(row["pred_mrna"])),
                "nrna": fmt_int(int(row["pred_nrna"])),
                "same": fmt_int(int(row["ambig_same"])),
                "opp": fmt_int(int(row["ambig_opp"])),
                "mean_zw": fmt_float(row["mean_zw"]),
            }
        )
    return out


def report_sample_rows(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, str]]:
    out = []
    for row in rows[:top_n]:
        locus = row.get("locus_key", "")
        out.append(
            {
                "locus": locus,
                "qname": compact_target(str(row["qname"]), 42),
                "region": f"{row['chrom']}:{fmt_int(int(row['fragment_start_1based']))}-{fmt_int(int(row['fragment_end']))}",
                "detail": row["predicted_detail"],
                "cigar": row["cigar"],
                "tlen": row["template_length"],
                "nh": row["nh"],
                "nm": row["nm"],
                "zw": row["zw"],
                "zc": row["zc"],
                "zs": row["zs"],
                "target": compact_target(str(row["target_label"]), 55),
            }
        )
    return out


def write_report(args: argparse.Namespace, scan: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    gdna_total = sum(
        count
        for (source, _binary, _detail), count in scan["global_counts"].items()
        if source == "gdna"
    )
    false_total = scan["record_counts"]["gdna_false_rna"]
    pred_mrna = sum(
        count
        for (source, _binary, detail), count in scan["global_counts"].items()
        if source == "gdna" and detail == "mrna"
    )
    pred_nrna = sum(
        count
        for (source, _binary, detail), count in scan["global_counts"].items()
        if source == "gdna" and detail == "nrna"
    )
    unspliced_false = sum(row["unspliced"] for row in rows["locus_rows"])
    ambig_false = sum(row["ambig_same"] + row["ambig_opp"] for row in rows["locus_rows"])
    top10_locus_false = sum(row["gdna_false_rna"] for row in rows["locus_rows"][:10])
    top10_window_false = sum(row["gdna_false_rna"] for row in rows["window_rows"][:10])

    text = f"""# VCaP gDNA False-RNA Hotspot Analysis - 2026-05-16

BAM: `{args.bam}`

Truth source is derived from query-name flowcell ID. This analysis focuses on gDNA-source fragments (`{', '.join(args.gdna_flowcell)}`) that Rigel assigned to an RNA pool (`mRNA` or `nRNA`). Counting rule is one primary read1 record per fragment.

Fragments counted: {fmt_int(scan['record_counts']['fragments_counted'])}

gDNA-source fragments: {fmt_int(gdna_total)}

gDNA -> RNA false positives: {fmt_int(false_total)} ({pct(false_total / gdna_total if gdna_total else None)} of gDNA-source fragments)

False-positive split: {fmt_int(pred_mrna)} mRNA ({pct(pred_mrna / false_total if false_total else None)} of false RNA) and {fmt_int(pred_nrna)} nRNA ({pct(pred_nrna / false_total if false_total else None)} of false RNA).

## What Dominates

The false-positive RNA problem is concentrated in unspliced, genic ambiguous fragments. Across false RNA calls, {pct(unspliced_false / false_total if false_total else None)} are `ZS=unspliced` and {pct(ambig_false / false_total if false_total else None)} are either `ambig_same_strand` or `ambig_opp_strand`.

The top 10 EM loci account for {pct(top10_locus_false / false_total if false_total else None)} of all gDNA -> RNA false positives. The top 10 genomic {args.window_size:,} bp windows account for {pct(top10_window_false / false_total if false_total else None)}. So this is not only a few pathological coordinates; it is a recurrent ambiguous-region failure mode, with some strong hotspots.

## Error Categories

{markdown_table(report_category_rows(rows['category_rows'], args.top_n), [
        ('Pred detail', 'detail'),
        ('ZC', 'zc'),
        ('ZS', 'zs'),
        ('Count', 'count'),
        ('False RNA frac', 'false_frac'),
        ('gDNA source frac', 'gdna_frac'),
    ])}

## Top 10kb Genomic Windows

{markdown_table(report_window_rows(rows['window_rows'], args.top_n), [
        ('Region', 'region'),
        ('False RNA', 'false'),
        ('False rate', 'rate'),
        ('gDNA total', 'gdna_total'),
        ('mRNA', 'mrna'),
        ('nRNA', 'nrna'),
        ('Top target', 'top_target'),
        ('Top category', 'top_cat'),
    ])}

## Local RNA Evidence Stratification

False RNA calls are not confined to windows with abundant RNA-source reads. This matters
because it separates two mechanisms: real local RNA expression can pull ambiguous gDNA
into RNA, but ambiguous gDNA can also self-seed RNA components in windows with little or
no RNA evidence.

{markdown_table(report_local_rna_rows(rows['window_rows']), [
        ('RNA-source fragments in 10kb window', 'rna_source_bin'),
        ('Windows', 'windows'),
        ('False RNA', 'false_rna'),
        ('False RNA frac', 'false_frac'),
        ('Weighted false rate', 'weighted_false_rate'),
    ])}

## Top EM Loci

The region span is the observed fragment span among false-positive fragments in the locus, not the full locus extent. Components such as `chr2:3`, `chr1:3`, and `chr7:3` are chromosome-scale mega-components rather than precise biological loci, so the 10kb windows and target tables above are the sharper coordinates for assay triage.

{markdown_table(report_locus_rows(rows['locus_rows'], args.top_n), [
        ('Locus', 'locus'),
        ('Observed region', 'region'),
        ('False RNA', 'false'),
        ('False frac', 'false_frac'),
        ('mRNA', 'mrna'),
        ('nRNA', 'nrna'),
        ('Same', 'same'),
        ('Opp', 'opp'),
        ('Mean ZW', 'mean_zw'),
        ('ZW >= .90', 'zw90'),
        ('Top target', 'target'),
    ])}

## Top Assigned RNA Targets

{markdown_table(report_target_rows(rows['target_rows'], args.top_n), [
        ('Target', 'target'),
        ('Observed region', 'region'),
        ('False RNA', 'false'),
        ('False frac', 'false_frac'),
        ('mRNA', 'mrna'),
        ('nRNA', 'nrna'),
        ('Same', 'same'),
        ('Opp', 'opp'),
        ('Mean ZW', 'mean_zw'),
    ])}

## Representative False-Positive Fragments From Top Loci

{markdown_table(report_sample_rows(rows['sample_rows'], min(args.top_n, 40)), [
        ('Locus', 'locus'),
        ('QNAME', 'qname'),
        ('Region', 'region'),
        ('Detail', 'detail'),
        ('CIGAR', 'cigar'),
        ('TLEN', 'tlen'),
        ('NH', 'nh'),
        ('NM', 'nm'),
        ('ZW', 'zw'),
        ('ZC', 'zc'),
        ('ZS', 'zs'),
        ('Target', 'target'),
    ])}

## High-Confidence False-Positive Examples

These examples have the largest `ZW` winner weights among gDNA-source fragments called RNA.

{markdown_table(report_sample_rows(rows['high_conf_rows'], min(args.top_n, 40)), [
        ('QNAME', 'qname'),
        ('Region', 'region'),
        ('Detail', 'detail'),
        ('CIGAR', 'cigar'),
        ('TLEN', 'tlen'),
        ('NH', 'nh'),
        ('NM', 'nm'),
        ('ZW', 'zw'),
        ('ZC', 'zc'),
        ('ZS', 'zs'),
        ('Target', 'target'),
    ])}

## Interpretation

1. The dominant false positives are not spliced RNA-like evidence. They are ordinary contiguous gDNA alignments (`ZS=unspliced`, usually `NH=1`) that overlap expressed annotated transcript or nRNA spans.
2. Ambiguous same-strand and opposite-strand labels mean the fragment has RNA candidates but no decisive splice-junction or transcript-unique evidence. In these cases, the EM can let RNA abundance in a large component overcome the generic gDNA candidate.
3. Local expression is only part of the failure. Many false-positive windows have little or no RNA-source support, meaning gDNA fragments themselves can self-seed an RNA component when the model has no conservative guard for ambiguous unspliced evidence.
4. Many false positives have high `ZW`, so this is not just low-confidence noise. The model is often confident after EM because the RNA component has enough mass to explain genic unspliced fragments.
5. A conservative diagnostic mode should treat ambiguous unspliced genic fragments from gDNA-compatible regions as weak RNA evidence unless they are supported by splice junctions, transcript-unique exonic geometry, or local RNA-only context.

## Recommended Changes

1. Add an RNA-call guard for `ZC in {{ambig_same_strand, ambig_opp_strand}}` and `ZS=unspliced`: require a stronger RNA-vs-gDNA posterior margin before emitting mRNA/nRNA as the hard winner.
2. Make gDNA a first-class competitor for ambiguous genic unspliced fragments even when locus RNA mass is high; practically, this means a conservative prior or likelihood floor for gDNA in these channels.
3. Split RNA evidence tiers in diagnostic output: splice-supported RNA, transcript-unique exonic RNA, and ambiguous-unspliced RNA. The last tier should not be used as strong positive RNA evidence in assays without additional support.
4. Preserve enough annotation in BAM output to inspect the runner-up gDNA posterior or RNA/gDNA posterior ratio. `ZW` alone reports winner confidence, but not the rejected gDNA probability.
5. Consider a diagnostic operating mode that defaults ambiguous unspliced RNA/gDNA ties to gDNA, trading RNA sensitivity for much lower false-positive RNA.
6. Break or regularize chromosome-scale mega-components for diagnostic calling. Ambiguous unspliced fragments should be judged against local evidence, not allowed to borrow RNA mass from distant transcripts in the same connected component.

## Artifacts

- Error categories: `{args.out_dir / 'gdna_false_rna_categories.tsv'}`
- Hotspot windows: `{args.out_dir / 'gdna_false_rna_windows.tsv'}`
- Hotspot EM loci: `{args.out_dir / 'gdna_false_rna_loci.tsv'}`
- Assigned RNA targets: `{args.out_dir / 'gdna_false_rna_targets.tsv'}`
- Representative fragments: `{args.out_dir / 'sample_false_rna_fragments.tsv'}`
- High-confidence examples: `{args.out_dir / 'high_confidence_false_rna_fragments.tsv'}`
- Summary JSON: `{args.out_dir / 'summary.json'}`
"""
    args.report.write_text(text)


def main() -> None:
    args = parse_args()
    scan = scan_bam(args)
    rows = build_rows(scan, args.window_size)
    write_outputs(args, scan, rows)
    print(f"counted {scan['record_counts']['fragments_counted']:,} fragments")
    print(f"gDNA -> RNA false positives: {scan['record_counts']['gdna_false_rna']:,}")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()