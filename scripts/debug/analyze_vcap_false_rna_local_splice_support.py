#!/usr/bin/env python3
"""Relate VCaP gDNA false-RNA calls to local splice-supported RNA evidence."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--window-size", type=int, default=10_000)
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


def fragment_start(read: pysam.AlignedSegment) -> int:
    if read.reference_start >= 0:
        return read.reference_start
    return 0


def scan(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    records = Counter()
    windows = defaultdict(Counter)

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen"] += 1
            if not is_representative(read):
                records["records_skipped"] += 1
                continue
            records["fragments_counted"] += 1
            chrom = read.reference_name or "*"
            window_start = (fragment_start(read) // args.window_size) * args.window_size
            key = (chrom, window_start)
            true_source = source_map.get(flowcell_from_qname(read.query_name), "unknown")
            detail = predicted_detail_from_zf(get_tag(read, "ZF", None))
            predicted_binary = binary_pool(detail)
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))

            windows[key]["total"] += 1
            if true_source == "rna":
                windows[key]["truth_rna_total"] += 1
            elif true_source == "gdna":
                windows[key]["truth_gdna_total"] += 1

            if detail in {"mrna", "nrna"}:
                windows[key]["pred_rna_total"] += 1
                if zs in SPLICE_SUPPORTED:
                    windows[key]["pred_rna_splice_supported"] += 1
                if zc == "unambig" and zs in SPLICE_SUPPORTED:
                    windows[key]["pred_rna_unambig_splice_supported"] += 1
                if zs == "unspliced" and zc in {"ambig_same_strand", "ambig_opp_strand"}:
                    windows[key]["pred_rna_ambig_unspliced"] += 1

            if true_source == "gdna" and predicted_binary == "rna":
                windows[key]["gdna_false_rna"] += 1
                if zs == "unspliced" and zc in {"ambig_same_strand", "ambig_opp_strand"}:
                    windows[key]["gdna_false_rna_ambig_unspliced"] += 1
            if true_source == "rna" and predicted_binary == "gdna":
                windows[key]["rna_false_gdna"] += 1

    return {"source_map": source_map, "records": records, "windows": windows, "window_size": args.window_size}


def support_bin(value: int) -> str:
    if value == 0:
        return "0"
    if value <= 10:
        return "1-10"
    if value <= 100:
        return "11-100"
    if value <= 1000:
        return "101-1000"
    return ">1000"


def build_rows(result: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    window_rows = []
    bin_counts = Counter()
    for (chrom, window_start), counts in result["windows"].items():
        false_count = counts.get("gdna_false_rna", 0)
        if false_count <= 0:
            continue
        splice = counts.get("pred_rna_splice_supported", 0)
        unambig_splice = counts.get("pred_rna_unambig_splice_supported", 0)
        truth_rna = counts.get("truth_rna_total", 0)
        row = {
            "chrom": chrom,
            "window_start_1based": window_start + 1,
            "window_end": window_start + result["window_size"],
            "gdna_source_total": counts.get("truth_gdna_total", 0),
            "gdna_false_rna": false_count,
            "gdna_false_rna_ambig_unspliced": counts.get("gdna_false_rna_ambig_unspliced", 0),
            "pred_rna_total": counts.get("pred_rna_total", 0),
            "pred_rna_splice_supported": splice,
            "pred_rna_unambig_splice_supported": unambig_splice,
            "pred_rna_ambig_unspliced": counts.get("pred_rna_ambig_unspliced", 0),
            "truth_rna_total": truth_rna,
            "rna_false_gdna": counts.get("rna_false_gdna", 0),
        }
        window_rows.append(row)
        for support_name, support_value in [
            ("splice_supported", splice),
            ("unambig_splice_supported", unambig_splice),
            ("truth_rna", truth_rna),
        ]:
            label = support_bin(int(support_value))
            bin_counts[(support_name, label, "windows")] += 1
            bin_counts[(support_name, label, "gdna_false_rna")] += false_count
            bin_counts[(support_name, label, "gdna_source_total")] += counts.get("truth_gdna_total", 0)
            bin_counts[(support_name, label, "ambig_unspliced_false")] += counts.get("gdna_false_rna_ambig_unspliced", 0)

    window_rows.sort(key=lambda row: -row["gdna_false_rna"])
    bin_rows = []
    for support_name in ["splice_supported", "unambig_splice_supported", "truth_rna"]:
        for label in ["0", "1-10", "11-100", "101-1000", ">1000"]:
            false_count = bin_counts[(support_name, label, "gdna_false_rna")]
            gdna_total = bin_counts[(support_name, label, "gdna_source_total")]
            bin_rows.append(
                {
                    "support_name": support_name,
                    "support_bin": label,
                    "windows": bin_counts[(support_name, label, "windows")],
                    "gdna_source_total": gdna_total,
                    "gdna_false_rna": false_count,
                    "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
                    "ambig_unspliced_false": bin_counts[(support_name, label, "ambig_unspliced_false")],
                }
            )
    return {"windows": window_rows, "bins": bin_rows}


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_outputs(args: argparse.Namespace, result: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(
        args.out_dir / "window_splice_support.tsv",
        rows["windows"],
        [
            "chrom",
            "window_start_1based",
            "window_end",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_ambig_unspliced",
            "pred_rna_total",
            "pred_rna_splice_supported",
            "pred_rna_unambig_splice_supported",
            "pred_rna_ambig_unspliced",
            "truth_rna_total",
            "rna_false_gdna",
        ],
    )
    write_tsv(
        args.out_dir / "splice_support_bins.tsv",
        rows["bins"],
        [
            "support_name",
            "support_bin",
            "windows",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
            "ambig_unspliced_false",
        ],
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
