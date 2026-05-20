#!/usr/bin/env python3
"""Profile alignment/tag features for VCaP true-gDNA false-RNA calls."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--threads", type=int, default=8)
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


def cigar_class(read: pysam.AlignedSegment) -> str:
    cigar = read.cigartuples or []
    if any(op == 3 for op, _length in cigar):
        return "contains_N"
    if any(op in {4, 5} for op, _length in cigar):
        return "soft_or_hard_clipped"
    if cigar and all(op == 0 for op, _length in cigar):
        return "all_M"
    return "other"


def template_bin(read: pysam.AlignedSegment) -> str:
    length = abs(int(read.template_length or 0))
    if length == 0:
        return "0_or_missing"
    if length <= 200:
        return "1-200"
    if length <= 350:
        return "201-350"
    if length <= 500:
        return "351-500"
    if length <= 1000:
        return "501-1000"
    return ">1000"


def mapq_bin(mapq: int) -> str:
    if mapq >= 60:
        return ">=60"
    if mapq >= 30:
        return "30-59"
    if mapq >= 10:
        return "10-29"
    return "0-9"


def zw_bin(value: Any) -> str:
    try:
        zw = float(value)
    except (TypeError, ValueError):
        return "missing"
    if zw < 0.25:
        return "<0.25"
    if zw < 0.50:
        return "0.25-0.50"
    if zw < 0.70:
        return "0.50-0.70"
    if zw < 0.90:
        return "0.70-0.90"
    if zw < 0.99:
        return "0.90-0.99"
    return ">=0.99"


def scan(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    records = Counter()
    feature_counts = Counter()
    combo_counts = Counter()
    by_truth_pred = Counter()
    numeric_sums = defaultdict(float)
    numeric_counts = Counter()

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen"] += 1
            if not is_representative(read):
                records["records_skipped"] += 1
                continue
            records["fragments_counted"] += 1
            true_source = source_map.get(flowcell_from_qname(read.query_name), "unknown")
            if true_source not in {"rna", "gdna"}:
                continue
            detail = predicted_detail_from_zf(get_tag(read, "ZF", None))
            predicted_binary = binary_pool(detail)
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))
            nh = int(get_tag(read, "NH", 1) or 1)
            nm = int(get_tag(read, "NM", 0) or 0)
            key = (true_source, predicted_binary, detail)
            by_truth_pred[key] += 1
            features = {
                "cigar_class": cigar_class(read),
                "nh": "1" if nh == 1 else "2" if nh == 2 else "3-10" if nh <= 10 else ">10",
                "nm": "0" if nm == 0 else "1" if nm == 1 else "2" if nm == 2 else ">=3",
                "mapq": mapq_bin(read.mapping_quality),
                "template": template_bin(read),
                "zc": zc,
                "zs": zs,
                "zw": zw_bin(get_tag(read, "ZW", None)),
            }
            for feature_name, feature_value in features.items():
                feature_counts[(true_source, predicted_binary, detail, feature_name, feature_value)] += 1
            combo_counts[(true_source, predicted_binary, detail, zc, zs, features["cigar_class"], features["nh"], features["nm"])] += 1
            numeric_sums[(true_source, predicted_binary, detail, "template_abs")] += abs(int(read.template_length or 0))
            numeric_counts[(true_source, predicted_binary, detail, "template_abs")] += 1
            numeric_sums[(true_source, predicted_binary, detail, "mapq")] += read.mapping_quality
            numeric_counts[(true_source, predicted_binary, detail, "mapq")] += 1

    return {
        "source_map": source_map,
        "records": records,
        "feature_counts": feature_counts,
        "combo_counts": combo_counts,
        "by_truth_pred": by_truth_pred,
        "numeric_sums": numeric_sums,
        "numeric_counts": numeric_counts,
    }


def build_rows(result: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    feature_rows = []
    totals = result["by_truth_pred"]
    for (true_source, predicted_binary, detail, feature_name, feature_value), count in result[
        "feature_counts"
    ].items():
        total = totals[(true_source, predicted_binary, detail)]
        feature_rows.append(
            {
                "true_source": true_source,
                "predicted_binary": predicted_binary,
                "predicted_detail": detail,
                "feature": feature_name,
                "value": feature_value,
                "count": count,
                "rate_within_group": count / total if total else 0.0,
            }
        )
    feature_rows.sort(key=lambda row: (row["true_source"], row["predicted_binary"], row["predicted_detail"], row["feature"], -row["count"]))

    combo_rows = []
    for (true_source, predicted_binary, detail, zc, zs, cigar, nh, nm), count in result[
        "combo_counts"
    ].most_common():
        total = totals[(true_source, predicted_binary, detail)]
        combo_rows.append(
            {
                "true_source": true_source,
                "predicted_binary": predicted_binary,
                "predicted_detail": detail,
                "zc": zc,
                "zs": zs,
                "cigar_class": cigar,
                "nh": nh,
                "nm": nm,
                "count": count,
                "rate_within_group": count / total if total else 0.0,
            }
        )

    group_rows = []
    for key, total in totals.items():
        true_source, predicted_binary, detail = key
        row = {
            "true_source": true_source,
            "predicted_binary": predicted_binary,
            "predicted_detail": detail,
            "count": total,
        }
        for metric in ["template_abs", "mapq"]:
            n = result["numeric_counts"][(true_source, predicted_binary, detail, metric)]
            row[f"mean_{metric}"] = result["numeric_sums"][(true_source, predicted_binary, detail, metric)] / n if n else 0.0
        group_rows.append(row)
    group_rows.sort(key=lambda row: (row["true_source"], row["predicted_binary"], row["predicted_detail"]))
    return {"features": feature_rows, "combos": combo_rows, "groups": group_rows}


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(args: argparse.Namespace, result: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "feature_counts.tsv",
        rows["features"],
        ["true_source", "predicted_binary", "predicted_detail", "feature", "value", "count", "rate_within_group"],
    )
    write_tsv(
        args.out_dir / "combo_counts.tsv",
        rows["combos"],
        ["true_source", "predicted_binary", "predicted_detail", "zc", "zs", "cigar_class", "nh", "nm", "count", "rate_within_group"],
    )
    write_tsv(
        args.out_dir / "group_summary.tsv",
        rows["groups"],
        ["true_source", "predicted_binary", "predicted_detail", "count", "mean_template_abs", "mean_mapq"],
    )
    summary = {"bam": str(args.bam), "source_map": result["source_map"], "record_counts": dict(result["records"])}
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def main() -> None:
    args = parse_args()
    result = scan(args)
    rows = build_rows(result)
    write_outputs(args, result, rows)
    print(f"counted {result['records']['fragments_counted']:,} fragments")
    print(f"wrote {args.out_dir}")


if __name__ == "__main__":
    main()
