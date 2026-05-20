#!/usr/bin/env python3
"""Evaluate local splice-support guards for ambiguous unspliced RNA calls."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08
SPLICE_SUPPORTED = {"spliced_annot", "spliced_implicit", "spliced_unannot"}
BINARY_ORDER = ["rna", "gdna", "unresolved"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--window-size", type=int, default=10_000)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--support-threshold", type=int, action="append", default=[0, 10, 100])
    return parser.parse_args()


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


def is_representative(read: pysam.AlignedSegment) -> bool:
    if read.is_secondary or read.is_supplementary:
        return False
    if read.is_paired:
        return read.is_read1
    return True


def window_key(read: pysam.AlignedSegment, window_size: int) -> tuple[str, int]:
    chrom = read.reference_name or "*"
    start = read.reference_start if read.reference_start >= 0 else 0
    return chrom, (start // window_size) * window_size


def is_ambig_unspliced_rna(detail: str, zc: str, zs: str) -> bool:
    return detail in {"mrna", "nrna"} and zs == "unspliced" and zc in {
        "ambig_same_strand",
        "ambig_opp_strand",
    }


def collect_support(args: argparse.Namespace) -> tuple[dict[tuple[str, int], Counter[str]], Counter[str]]:
    support = defaultdict(Counter)
    records = Counter()
    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen_pass1"] += 1
            if not is_representative(read):
                records["records_skipped_pass1"] += 1
                continue
            detail = predicted_detail_from_zf(get_tag(read, "ZF", None))
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))
            if detail in {"mrna", "nrna"} and zs in SPLICE_SUPPORTED:
                key = window_key(read, args.window_size)
                support[key]["splice_supported"] += 1
                if zc == "unambig":
                    support[key]["unambig_splice_supported"] += 1
    return support, records


def should_guard(
    detail: str,
    zc: str,
    zs: str,
    support_counts: Counter[str],
    guard_kind: str,
    threshold: int,
) -> bool:
    if not is_ambig_unspliced_rna(detail, zc, zs):
        return False
    if guard_kind == "any_splice_le":
        return support_counts.get("splice_supported", 0) <= threshold
    if guard_kind == "unambig_splice_le":
        return support_counts.get("unambig_splice_supported", 0) <= threshold
    if guard_kind == "nrna_unambig_splice_le":
        return detail == "nrna" and support_counts.get("unambig_splice_supported", 0) <= threshold
    if guard_kind == "mrna_unambig_splice_le":
        return detail == "mrna" and support_counts.get("unambig_splice_supported", 0) <= threshold
    return False


def scan_with_guards(args: argparse.Namespace, support: dict[tuple[str, int], Counter[str]]) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    thresholds = sorted(set(args.support_threshold))
    guard_kinds = [
        "baseline",
        "any_splice_le",
        "unambig_splice_le",
        "nrna_unambig_splice_le",
        "mrna_unambig_splice_le",
    ]
    confusion = {(kind, threshold): Counter() for kind in guard_kinds for threshold in thresholds}
    guarded = Counter()
    records = Counter()

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen_pass2"] += 1
            if not is_representative(read):
                records["records_skipped_pass2"] += 1
                continue
            records["fragments_counted"] += 1
            true_source = source_map.get(flowcell_from_qname(read.query_name), "unknown")
            if true_source not in {"rna", "gdna"}:
                continue
            detail = predicted_detail_from_zf(get_tag(read, "ZF", None))
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))
            key = window_key(read, args.window_size)
            support_counts = support.get(key, Counter())
            for kind in guard_kinds:
                for threshold in thresholds:
                    adjusted_detail = detail
                    if kind != "baseline" and should_guard(
                        detail, zc, zs, support_counts, kind, threshold
                    ):
                        adjusted_detail = "gdna"
                        guarded[(kind, threshold, true_source, detail, zc)] += 1
                    confusion[(kind, threshold)][(true_source, binary_pool(adjusted_detail))] += 1
    return {
        "source_map": source_map,
        "thresholds": thresholds,
        "guard_kinds": guard_kinds,
        "confusion": confusion,
        "guarded": guarded,
        "records": records,
    }


def rate(count: int, total: int) -> float:
    return count / total if total else 0.0


def build_rows(result: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    metric_rows = []
    for kind in result["guard_kinds"]:
        for threshold in result["thresholds"]:
            counts = result["confusion"][(kind, threshold)]
            rna_total = sum(counts[("rna", pred)] for pred in BINARY_ORDER)
            gdna_total = sum(counts[("gdna", pred)] for pred in BINARY_ORDER)
            metric_rows.append(
                {
                    "guard_kind": kind,
                    "support_threshold": threshold,
                    "rna_total": rna_total,
                    "gdna_total": gdna_total,
                    "rna_to_rna": counts[("rna", "rna")],
                    "rna_to_gdna": counts[("rna", "gdna")],
                    "rna_to_unresolved": counts[("rna", "unresolved")],
                    "gdna_to_rna": counts[("gdna", "rna")],
                    "gdna_to_gdna": counts[("gdna", "gdna")],
                    "gdna_to_unresolved": counts[("gdna", "unresolved")],
                    "rna_to_rna_rate": rate(counts[("rna", "rna")], rna_total),
                    "rna_to_gdna_rate": rate(counts[("rna", "gdna")], rna_total),
                    "gdna_to_rna_rate": rate(counts[("gdna", "rna")], gdna_total),
                    "gdna_to_gdna_rate": rate(counts[("gdna", "gdna")], gdna_total),
                    "guarded_total": sum(
                        count
                        for (guard_kind, guard_threshold, *_rest), count in result["guarded"].items()
                        if guard_kind == kind and guard_threshold == threshold
                    ),
                }
            )
    guarded_rows = []
    for (kind, threshold, true_source, original_detail, zc), count in result["guarded"].items():
        guarded_rows.append(
            {
                "guard_kind": kind,
                "support_threshold": threshold,
                "true_source": true_source,
                "original_detail": original_detail,
                "zc": zc,
                "count": count,
            }
        )
    guarded_rows.sort(key=lambda row: (row["guard_kind"], row["support_threshold"], -row["count"]))
    return {"metrics": metric_rows, "guarded": guarded_rows}


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(
    args: argparse.Namespace,
    support_records: Counter[str],
    result: dict[str, Any],
    rows: dict[str, list[dict[str, Any]]],
) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "local_splice_guard_metrics.tsv",
        rows["metrics"],
        [
            "guard_kind",
            "support_threshold",
            "rna_total",
            "gdna_total",
            "rna_to_rna",
            "rna_to_gdna",
            "rna_to_unresolved",
            "gdna_to_rna",
            "gdna_to_gdna",
            "gdna_to_unresolved",
            "rna_to_rna_rate",
            "rna_to_gdna_rate",
            "gdna_to_rna_rate",
            "gdna_to_gdna_rate",
            "guarded_total",
        ],
    )
    write_tsv(
        args.out_dir / "local_splice_guarded_breakdown.tsv",
        rows["guarded"],
        ["guard_kind", "support_threshold", "true_source", "original_detail", "zc", "count"],
    )
    summary = {
        "bam": str(args.bam),
        "source_map": result["source_map"],
        "window_size": args.window_size,
        "support_record_counts": dict(support_records),
        "scan_record_counts": dict(result["records"]),
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def main() -> None:
    args = parse_args()
    support, support_records = collect_support(args)
    result = scan_with_guards(args, support)
    rows = build_rows(result)
    write_outputs(args, support_records, result, rows)
    print(f"counted {result['records']['fragments_counted']:,} fragments")
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
