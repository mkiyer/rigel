#!/usr/bin/env python3
"""Evaluate post-hoc guards for ambiguous unspliced RNA calls in VCaP mix."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08

DETAIL_ORDER = ["mrna", "nrna", "gdna", "unresolved"]
BINARY_ORDER = ["rna", "gdna", "unresolved"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--zw-threshold",
        type=float,
        action="append",
        default=[0.25, 0.50, 0.70, 0.90, 0.99, 1.01],
        help="Reassign guarded RNA winners to gDNA when ZW is below this threshold. 1.01 means all.",
    )
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


def zw_float(value: Any) -> float | None:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def is_guard_candidate(detail: str, zc: str, zs: str, guard_kind: str) -> bool:
    if detail not in {"mrna", "nrna"} or zs != "unspliced":
        return False
    if guard_kind == "ambig_unspliced":
        return zc in {"ambig_same_strand", "ambig_opp_strand"}
    if guard_kind == "ambig_opp_unspliced":
        return zc == "ambig_opp_strand"
    if guard_kind == "ambig_same_unspliced":
        return zc == "ambig_same_strand"
    if guard_kind == "nrna_ambig_unspliced":
        return detail == "nrna" and zc in {"ambig_same_strand", "ambig_opp_strand"}
    if guard_kind == "mrna_ambig_unspliced":
        return detail == "mrna" and zc in {"ambig_same_strand", "ambig_opp_strand"}
    return False


def empty_confusion() -> Counter[tuple[str, str]]:
    return Counter()


def scan(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    thresholds = sorted(set(args.zw_threshold))
    guard_kinds = [
        "baseline",
        "ambig_unspliced",
        "ambig_opp_unspliced",
        "ambig_same_unspliced",
        "nrna_ambig_unspliced",
        "mrna_ambig_unspliced",
    ]
    confusion = {(kind, threshold): empty_confusion() for kind in guard_kinds for threshold in thresholds}
    detail_confusion = {
        (kind, threshold): Counter() for kind in guard_kinds for threshold in thresholds
    }
    guarded_counts = Counter()
    guarded_by_truth = Counter()
    candidate_zw_bins = Counter()
    records = Counter()

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
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))
            zw = zw_float(get_tag(read, "ZW", None))
            if detail in {"mrna", "nrna"} and zc in {"ambig_same_strand", "ambig_opp_strand"} and zs == "unspliced":
                bin_label = "missing" if zw is None else (
                    "<0.25" if zw < 0.25 else "0.25-0.50" if zw < 0.50 else "0.50-0.70" if zw < 0.70 else "0.70-0.90" if zw < 0.90 else "0.90-0.99" if zw < 0.99 else ">=0.99"
                )
                candidate_zw_bins[(true_source, detail, zc, bin_label)] += 1

            for kind in guard_kinds:
                for threshold in thresholds:
                    adjusted_detail = detail
                    if kind != "baseline" and is_guard_candidate(detail, zc, zs, kind):
                        if zw is None or zw < threshold:
                            adjusted_detail = "gdna"
                            guarded_counts[(kind, threshold)] += 1
                            guarded_by_truth[(kind, threshold, true_source, detail, zc)] += 1
                    predicted_binary = binary_pool(adjusted_detail)
                    confusion[(kind, threshold)][(true_source, predicted_binary)] += 1
                    detail_confusion[(kind, threshold)][(true_source, adjusted_detail)] += 1

    return {
        "source_map": source_map,
        "thresholds": thresholds,
        "guard_kinds": guard_kinds,
        "confusion": confusion,
        "detail_confusion": detail_confusion,
        "guarded_counts": guarded_counts,
        "guarded_by_truth": guarded_by_truth,
        "candidate_zw_bins": candidate_zw_bins,
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
                    "zw_threshold": threshold,
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
                    "guarded_total": result["guarded_counts"][(kind, threshold)],
                }
            )

    detail_rows = []
    for (kind, threshold), counts in result["detail_confusion"].items():
        for true_source in ["rna", "gdna"]:
            total = sum(counts[(true_source, detail)] for detail in DETAIL_ORDER)
            for detail in DETAIL_ORDER:
                count = counts[(true_source, detail)]
                detail_rows.append(
                    {
                        "guard_kind": kind,
                        "zw_threshold": threshold,
                        "true_source": true_source,
                        "predicted_detail": detail,
                        "count": count,
                        "row_rate": rate(count, total),
                    }
                )

    guarded_rows = []
    for (kind, threshold, true_source, original_detail, zc), count in result[
        "guarded_by_truth"
    ].items():
        guarded_rows.append(
            {
                "guard_kind": kind,
                "zw_threshold": threshold,
                "true_source": true_source,
                "original_detail": original_detail,
                "zc": zc,
                "count": count,
            }
        )
    guarded_rows.sort(key=lambda row: (row["guard_kind"], row["zw_threshold"], -row["count"]))

    zw_rows = []
    for (true_source, detail, zc, bin_label), count in result["candidate_zw_bins"].items():
        zw_rows.append(
            {
                "true_source": true_source,
                "predicted_detail": detail,
                "zc": zc,
                "zw_bin": bin_label,
                "count": count,
            }
        )
    zw_rows.sort(key=lambda row: (row["true_source"], row["predicted_detail"], row["zc"], row["zw_bin"]))

    return {"metrics": metric_rows, "details": detail_rows, "guarded": guarded_rows, "zw_bins": zw_rows}


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(args: argparse.Namespace, result: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "guard_metrics.tsv",
        rows["metrics"],
        [
            "guard_kind",
            "zw_threshold",
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
        args.out_dir / "guard_detail_confusion.tsv",
        rows["details"],
        ["guard_kind", "zw_threshold", "true_source", "predicted_detail", "count", "row_rate"],
    )
    write_tsv(
        args.out_dir / "guarded_breakdown.tsv",
        rows["guarded"],
        ["guard_kind", "zw_threshold", "true_source", "original_detail", "zc", "count"],
    )
    write_tsv(
        args.out_dir / "candidate_zw_bins.tsv",
        rows["zw_bins"],
        ["true_source", "predicted_detail", "zc", "zw_bin", "count"],
    )
    summary = {
        "bam": str(args.bam),
        "source_map": result["source_map"],
        "record_counts": dict(result["records"]),
        "thresholds": result["thresholds"],
    }
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
