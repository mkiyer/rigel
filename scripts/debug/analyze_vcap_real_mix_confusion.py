#!/usr/bin/env python3
"""Pool-confusion analysis for a real RNA/gDNA mix annotated by Rigel."""

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

POOL_ORDER = ["rna", "gdna", "unresolved"]
DETAIL_ORDER = ["mrna", "nrna", "gdna", "unresolved"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bam",
        type=Path,
        required=True,
        help="Rigel annotated BAM to scan.",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Directory for TSV/JSON outputs.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
        help="Markdown report path.",
    )
    parser.add_argument(
        "--rna-flowcell",
        action="append",
        default=[],
        help="Flowcell ID known to come from the RNA source. May be repeated.",
    )
    parser.add_argument(
        "--gdna-flowcell",
        action="append",
        default=[],
        help="Flowcell ID known to come from the gDNA source. May be repeated.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="BGZF decompression threads for pysam.",
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


def pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * value:.{digits}f}%"


def fmt_int(value: int) -> str:
    return f"{value:,}"


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
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


def scan_bam(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    record_counts = Counter()
    flowcell_counts = Counter()
    flowcell_detail_counts = Counter()
    flowcell_binary_counts = Counter()
    truth_binary_counts = Counter()
    truth_detail_counts = Counter()
    zf_counts = Counter()
    error_zc_zs_counts = Counter()
    read_length_sum = defaultdict(int)

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
            zc_value = get_tag(read, "ZC", "missing")
            zs_value = get_tag(read, "ZS", "missing")

            flowcell_counts[(flowcell, true_source)] += 1
            flowcell_detail_counts[(flowcell, true_source, predicted_detail)] += 1
            flowcell_binary_counts[(flowcell, true_source, predicted_binary)] += 1
            truth_binary_counts[(true_source, predicted_binary)] += 1
            truth_detail_counts[(true_source, predicted_detail)] += 1
            zf_counts[(flowcell, true_source, zf_value if zf_value is not None else "missing")] += 1
            if read.query_length is not None:
                read_length_sum[(flowcell, true_source)] += read.query_length

            if true_source in {"rna", "gdna"} and predicted_binary != true_source:
                error_zc_zs_counts[
                    (flowcell, true_source, predicted_binary, predicted_detail, zc_value, zs_value)
                ] += 1

    return {
        "source_map": source_map,
        "record_counts": record_counts,
        "flowcell_counts": flowcell_counts,
        "flowcell_detail_counts": flowcell_detail_counts,
        "flowcell_binary_counts": flowcell_binary_counts,
        "truth_binary_counts": truth_binary_counts,
        "truth_detail_counts": truth_detail_counts,
        "zf_counts": zf_counts,
        "error_zc_zs_counts": error_zc_zs_counts,
        "read_length_sum": read_length_sum,
    }


def build_rows(scan: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    flowcell_rows = []
    for (flowcell, source), total in sorted(scan["flowcell_counts"].items()):
        detail_counts = {
            pool: scan["flowcell_detail_counts"].get((flowcell, source, pool), 0)
            for pool in DETAIL_ORDER
        }
        binary_counts = {
            pool: scan["flowcell_binary_counts"].get((flowcell, source, pool), 0)
            for pool in POOL_ORDER
        }
        mean_read_length = scan["read_length_sum"].get((flowcell, source), 0) / total if total else 0.0
        flowcell_rows.append(
            {
                "flowcell": flowcell,
                "true_source": source,
                "fragments": total,
                "mean_read_length": f"{mean_read_length:.1f}",
                "pred_rna": binary_counts["rna"],
                "pred_gdna": binary_counts["gdna"],
                "pred_unresolved": binary_counts["unresolved"],
                "pred_mrna": detail_counts["mrna"],
                "pred_nrna": detail_counts["nrna"],
                "pred_gdna_detail": detail_counts["gdna"],
                "pred_unresolved_detail": detail_counts["unresolved"],
                "pred_rna_rate": binary_counts["rna"] / total if total else 0.0,
                "pred_gdna_rate": binary_counts["gdna"] / total if total else 0.0,
                "pred_unresolved_rate": binary_counts["unresolved"] / total if total else 0.0,
            }
        )

    confusion_count_rows = []
    confusion_matrix_rows = []
    true_sources = sorted({source for source, _ in scan["truth_binary_counts"]})
    for true_source in true_sources:
        total = sum(scan["truth_binary_counts"].get((true_source, pred), 0) for pred in POOL_ORDER)
        for predicted_pool in POOL_ORDER:
            count = scan["truth_binary_counts"].get((true_source, predicted_pool), 0)
            confusion_count_rows.append(
                {"true_source": true_source, "predicted_pool": predicted_pool, "count": count}
            )
            confusion_matrix_rows.append(
                {
                    "true_source": true_source,
                    "predicted_pool": predicted_pool,
                    "count": count,
                    "row_rate": count / total if total else 0.0,
                }
            )

    detail_rows = []
    for true_source in true_sources:
        total = sum(scan["truth_detail_counts"].get((true_source, pool), 0) for pool in DETAIL_ORDER)
        for predicted_pool in DETAIL_ORDER:
            count = scan["truth_detail_counts"].get((true_source, predicted_pool), 0)
            detail_rows.append(
                {
                    "true_source": true_source,
                    "predicted_pool": predicted_pool,
                    "count": count,
                    "row_rate": count / total if total else 0.0,
                }
            )

    zf_rows = []
    zf_items = sorted(scan["zf_counts"].items(), key=lambda item: (item[0][0], item[0][1], str(item[0][2])))
    for (flowcell, source, zf_value), count in zf_items:
        zf_rows.append(
            {"flowcell": flowcell, "true_source": source, "zf": zf_value, "count": count}
        )

    error_rows = []
    for key, count in scan["error_zc_zs_counts"].most_common():
        flowcell, true_source, predicted_binary, predicted_detail, zc_value, zs_value = key
        error_rows.append(
            {
                "flowcell": flowcell,
                "true_source": true_source,
                "predicted_pool": predicted_binary,
                "predicted_detail": predicted_detail,
                "zc": zc_value,
                "zs": zs_value,
                "count": count,
            }
        )

    return {
        "flowcell_rows": flowcell_rows,
        "confusion_count_rows": confusion_count_rows,
        "confusion_matrix_rows": confusion_matrix_rows,
        "detail_rows": detail_rows,
        "zf_rows": zf_rows,
        "error_rows": error_rows,
    }


def write_outputs(args: argparse.Namespace, scan: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)

    write_tsv(
        args.out_dir / "flowcell_summary.tsv",
        rows["flowcell_rows"],
        [
            "flowcell",
            "true_source",
            "fragments",
            "mean_read_length",
            "pred_rna",
            "pred_gdna",
            "pred_unresolved",
            "pred_mrna",
            "pred_nrna",
            "pred_gdna_detail",
            "pred_unresolved_detail",
            "pred_rna_rate",
            "pred_gdna_rate",
            "pred_unresolved_rate",
        ],
    )
    write_tsv(
        args.out_dir / "confusion_counts.tsv",
        rows["confusion_count_rows"],
        ["true_source", "predicted_pool", "count"],
    )
    write_tsv(
        args.out_dir / "confusion_matrix.tsv",
        rows["confusion_matrix_rows"],
        ["true_source", "predicted_pool", "count", "row_rate"],
    )
    write_tsv(
        args.out_dir / "detailed_predicted_pool.tsv",
        rows["detail_rows"],
        ["true_source", "predicted_pool", "count", "row_rate"],
    )
    write_tsv(
        args.out_dir / "zf_by_flowcell.tsv",
        rows["zf_rows"],
        ["flowcell", "true_source", "zf", "count"],
    )
    write_tsv(
        args.out_dir / "error_zc_zs_breakdown.tsv",
        rows["error_rows"],
        ["flowcell", "true_source", "predicted_pool", "predicted_detail", "zc", "zs", "count"],
    )

    summary = {
        "bam": str(args.bam),
        "source_map": scan["source_map"],
        "record_counts": dict(scan["record_counts"]),
    }
    with (args.out_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    write_report(args, scan, rows)


def row_value(rows: list[dict[str, Any]], true_source: str, predicted_pool: str) -> int:
    for row in rows:
        if row["true_source"] == true_source and row["predicted_pool"] == predicted_pool:
            return int(row["count"])
    return 0


def row_rate(rows: list[dict[str, Any]], true_source: str, predicted_pool: str) -> float | None:
    for row in rows:
        if row["true_source"] == true_source and row["predicted_pool"] == predicted_pool:
            return float(row["row_rate"])
    return None


def report_confusion_rows(matrix_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    report_rows = []
    for true_source in ["rna", "gdna", "unknown"]:
        source_rows = [row for row in matrix_rows if row["true_source"] == true_source]
        if not source_rows:
            continue
        total = sum(int(row["count"]) for row in source_rows)
        report_rows.append(
            {
                "true_source": true_source,
                "truth_fragments": fmt_int(total),
                "pred_rna": f"{fmt_int(row_value(matrix_rows, true_source, 'rna'))} ({pct(row_rate(matrix_rows, true_source, 'rna'))})",
                "pred_gdna": f"{fmt_int(row_value(matrix_rows, true_source, 'gdna'))} ({pct(row_rate(matrix_rows, true_source, 'gdna'))})",
                "pred_unresolved": f"{fmt_int(row_value(matrix_rows, true_source, 'unresolved'))} ({pct(row_rate(matrix_rows, true_source, 'unresolved'))})",
            }
        )
    return report_rows


def report_detail_rows(detail_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    report_rows = []
    for true_source in ["rna", "gdna", "unknown"]:
        source_rows = [row for row in detail_rows if row["true_source"] == true_source]
        if not source_rows:
            continue
        total = sum(int(row["count"]) for row in source_rows)
        row: dict[str, str] = {"true_source": true_source, "truth_fragments": fmt_int(total)}
        for pool in DETAIL_ORDER:
            row[f"pred_{pool}"] = (
                f"{fmt_int(row_value(detail_rows, true_source, pool))} "
                f"({pct(row_rate(detail_rows, true_source, pool))})"
            )
        report_rows.append(row)
    return report_rows


def report_flowcell_rows(flowcell_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    report_rows = []
    for row in sorted(flowcell_rows, key=lambda item: item["flowcell"]):
        report_rows.append(
            {
                "flowcell": row["flowcell"],
                "source": row["true_source"],
                "fragments": fmt_int(int(row["fragments"])),
                "read_length": row["mean_read_length"],
                "pred_rna": pct(float(row["pred_rna_rate"])),
                "pred_gdna": pct(float(row["pred_gdna_rate"])),
                "pred_unresolved": pct(float(row["pred_unresolved_rate"])),
            }
        )
    return report_rows


def write_report(args: argparse.Namespace, scan: dict[str, Any], rows: dict[str, list[dict[str, Any]]]) -> None:
    matrix_rows = rows["confusion_matrix_rows"]
    detail_rows = rows["detail_rows"]
    counted = scan["record_counts"]["fragments_counted"]
    pred_rna_total = sum(row_value(matrix_rows, source, "rna") for source in ["rna", "gdna", "unknown"])
    pred_gdna_total = sum(row_value(matrix_rows, source, "gdna") for source in ["rna", "gdna", "unknown"])
    pred_unresolved_total = sum(
        row_value(matrix_rows, source, "unresolved") for source in ["rna", "gdna", "unknown"]
    )
    assigned_total = pred_rna_total + pred_gdna_total
    rna_to_gdna = row_rate(matrix_rows, "rna", "gdna")
    gdna_to_rna = row_rate(matrix_rows, "gdna", "rna")
    rna_recall = row_rate(matrix_rows, "rna", "rna")
    gdna_recall = row_rate(matrix_rows, "gdna", "gdna")

    top_errors = []
    for row in rows["error_rows"][:20]:
        top_errors.append(
            {
                "flowcell": row["flowcell"],
                "true_source": row["true_source"],
                "pred": row["predicted_pool"],
                "detail": row["predicted_detail"],
                "zc": row["zc"],
                "zs": row["zs"],
                "count": fmt_int(int(row["count"])),
            }
        )

    source_lines = [f"- `{flowcell}` -> `{source}`" for flowcell, source in sorted(scan["source_map"].items())]
    if not source_lines:
        source_lines = ["- No explicit source flowcells were provided; sources are reported as `unknown`." ]

    text = f"""# VCaP RNA/gDNA Real Mix Pool Confusion - 2026-05-15

BAM: `{args.bam}`

Counting rule: one primary read1 record per fragment; secondary, supplementary, and read2 records are skipped. Truth source is derived from the query-name flowcell ID.

Source mapping:

{chr(10).join(source_lines)}

Fragments counted: {fmt_int(counted)}

## Headline

Rigel reports {pct(pred_gdna_total / counted if counted else None)} gDNA and {pct(pred_rna_total / counted if counted else None)} RNA across all counted fragments, with {pct(pred_unresolved_total / counted if counted else None)} unresolved. Among assigned RNA/gDNA fragments only, the reported composition is {pct(pred_gdna_total / assigned_total if assigned_total else None)} gDNA and {pct(pred_rna_total / assigned_total if assigned_total else None)} RNA.

Using the flowcell labels as truth, RNA-source fragments are called RNA at {pct(rna_recall)} and gDNA at {pct(rna_to_gdna)}. gDNA-source fragments are called gDNA at {pct(gdna_recall)} and RNA at {pct(gdna_to_rna)}.

## Flowcell Summary

{markdown_table(report_flowcell_rows(rows['flowcell_rows']), [
        ('Flowcell', 'flowcell'),
        ('Truth', 'source'),
        ('Fragments', 'fragments'),
        ('Mean read len', 'read_length'),
        ('Pred RNA', 'pred_rna'),
        ('Pred gDNA', 'pred_gdna'),
        ('Pred unresolved', 'pred_unresolved'),
    ])}

## RNA versus gDNA Confusion Matrix

{markdown_table(report_confusion_rows(matrix_rows), [
        ('True source', 'true_source'),
        ('Fragments', 'truth_fragments'),
        ('Pred RNA', 'pred_rna'),
        ('Pred gDNA', 'pred_gdna'),
        ('Pred unresolved', 'pred_unresolved'),
    ])}

## Detailed Rigel Pool Calls

{markdown_table(report_detail_rows(detail_rows), [
        ('True source', 'true_source'),
        ('Fragments', 'truth_fragments'),
        ('Pred mRNA', 'pred_mrna'),
        ('Pred nRNA', 'pred_nrna'),
        ('Pred gDNA', 'pred_gdna'),
        ('Pred unresolved', 'pred_unresolved'),
    ])}

## Largest RNA/gDNA Pool Errors

{markdown_table(top_errors, [
        ('Flowcell', 'flowcell'),
        ('Truth', 'true_source'),
        ('Pred', 'pred'),
        ('Detail', 'detail'),
        ('ZC', 'zc'),
        ('ZS', 'zs'),
        ('Count', 'count'),
    ])}

## Artifacts

- Flowcell summary: `{args.out_dir / 'flowcell_summary.tsv'}`
- RNA/gDNA confusion counts: `{args.out_dir / 'confusion_counts.tsv'}`
- RNA/gDNA confusion matrix: `{args.out_dir / 'confusion_matrix.tsv'}`
- Detailed predicted pools: `{args.out_dir / 'detailed_predicted_pool.tsv'}`
- ZF by flowcell: `{args.out_dir / 'zf_by_flowcell.tsv'}`
- Error ZC/ZS breakdown: `{args.out_dir / 'error_zc_zs_breakdown.tsv'}`
- Summary JSON: `{args.out_dir / 'summary.json'}`
"""
    args.report.write_text(text)


def main() -> None:
    args = parse_args()
    scan = scan_bam(args)
    rows = build_rows(scan)
    write_outputs(args, scan, rows)
    print(f"counted {scan['record_counts']['fragments_counted']:,} fragments")
    print(f"wrote {args.report}")


if __name__ == "__main__":
    main()