#!/usr/bin/env python3
"""Analyze local VCaP AR-locus Rigel annotations against flowcell truth."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08

TYPE_LABELS = {0: "INTERGENIC", 1: "INTRON", 2: "EXON"}
STRAND_LABELS = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
POOL_ORDER = ["rna", "gdna", "unresolved"]
DETAIL_ORDER = ["mrna", "nrna", "gdna", "unresolved"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True, help="Rigel annotated BAM.")
    parser.add_argument("--index-dir", type=Path, required=True, help="Rigel index directory.")
    parser.add_argument("--quant-dir", type=Path, required=True, help="Rigel quant output directory.")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for analysis outputs.")
    parser.add_argument("--report", type=Path, required=True, help="Markdown report path.")
    parser.add_argument("--chrom", default="chrX", help="Chromosome to analyze.")
    parser.add_argument("--start", type=int, default=66_929_868, help="1-based inclusive region start.")
    parser.add_argument("--end", type=int, default=69_019_700, help="1-based inclusive region end.")
    parser.add_argument("--window-size", type=int, default=10_000, help="Window size in bp.")
    parser.add_argument("--rna-flowcell", action="append", default=[])
    parser.add_argument("--gdna-flowcell", action="append", default=[])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sample-limit", type=int, default=200)
    parser.add_argument(
        "--false-rna-bam",
        type=Path,
        help="Optional BAM containing all records from true-gDNA fragments assigned to RNA.",
    )
    return parser.parse_args()


def source_map_from_args(args: argparse.Namespace) -> dict[str, str]:
    source_map: dict[str, str] = {}
    for flowcell in args.rna_flowcell:
        source_map[flowcell] = "rna"
    for flowcell in args.gdna_flowcell:
        previous = source_map.get(flowcell)
        if previous is not None and previous != "gdna":
            raise ValueError(f"Flowcell {flowcell} assigned to both RNA and gDNA")
        source_map[flowcell] = "gdna"
    return source_map


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


def binary_pool(predicted_detail: str) -> str:
    if predicted_detail in {"mrna", "nrna"}:
        return "rna"
    return predicted_detail


def is_representative(read: pysam.AlignedSegment) -> bool:
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


def cigar_has_refskip(read: pysam.AlignedSegment) -> bool:
    return any(op == 3 for op, _length in (read.cigartuples or ()))


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


def fmt_pct(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "n/a"
    return f"{100.0 * value:.{digits}f}%"


def write_tsv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, Any]], columns: list[str] | None = None) -> str:
    if not rows:
        return "_No rows._"
    selected_columns = columns or list(rows[0])
    header = "| " + " | ".join(selected_columns) + " |"
    separator = "| " + " | ".join("---" for _ in selected_columns) + " |"
    body = []
    for row in rows:
        values = []
        for column in selected_columns:
            value = row.get(column, "")
            values.append(str(value).replace("|", "\\|"))
        body.append("| " + " | ".join(values) + " |")
    return "\n".join([header, separator, *body])


def load_regions(index_dir: Path, chrom: str, start0: int, end0: int) -> dict[str, np.ndarray]:
    regions = pd.read_feather(index_dir / "regions.feather")
    subset = regions[
        (regions["ref_name"] == chrom) & (regions["end"] > start0) & (regions["start"] < end0)
    ].copy()
    subset = subset.sort_values("start")
    return {
        "start": subset["start"].to_numpy(np.int64),
        "end": subset["end"].to_numpy(np.int64),
        "type": subset["type"].to_numpy(np.int16),
        "strand": subset["strand"].to_numpy(np.int16),
    }


def classify_region_overlap(
    region_arrays: dict[str, np.ndarray],
    span_start: int,
    span_end: int,
) -> tuple[str, str, str, int]:
    starts = region_arrays["start"]
    ends = region_arrays["end"]
    types = region_arrays["type"]
    strands = region_arrays["strand"]
    if starts.size == 0:
        return "NONE", "NONE", "NONE", 0

    idx = max(0, int(np.searchsorted(starts, span_start, side="right")) - 1)
    type_bp: Counter[str] = Counter()
    strand_bp: Counter[str] = Counter()
    mask_bits = 0
    total_overlap = 0
    while idx < starts.size and starts[idx] < span_end:
        ov_start = max(span_start, int(starts[idx]))
        ov_end = min(span_end, int(ends[idx]))
        if ov_end > ov_start:
            overlap = ov_end - ov_start
            type_label = TYPE_LABELS.get(int(types[idx]), str(types[idx]))
            strand_label = STRAND_LABELS.get(int(strands[idx]), str(strands[idx]))
            type_bp[type_label] += overlap
            strand_bp[strand_label] += overlap
            mask_bits |= 1 << int(types[idx])
            total_overlap += overlap
        idx += 1

    mask_labels = [TYPE_LABELS[i] for i in range(3) if mask_bits & (1 << i)]
    exact_mask = "+".join(mask_labels) if mask_labels else "NONE"
    dominant_type = type_bp.most_common(1)[0][0] if type_bp else "NONE"
    dominant_strand = strand_bp.most_common(1)[0][0] if strand_bp else "NONE"
    return exact_mask, dominant_type, dominant_strand, total_overlap


def summarize_counter(counter: Counter[tuple[Any, ...]], columns: list[str]) -> list[dict[str, Any]]:
    rows = []
    for key, count in counter.most_common():
        row = dict(zip(columns, key, strict=True))
        row["count"] = count
        rows.append(row)
    return rows


def scan_bam(args: argparse.Namespace) -> dict[str, Any]:
    source_map = source_map_from_args(args)
    start0 = args.start - 1
    end0 = args.end
    region_arrays = load_regions(args.index_dir, args.chrom, start0, end0)

    records = Counter()
    confusion_binary = Counter()
    confusion_detail = Counter()
    flowcell_counts = Counter()
    false_target = Counter()
    false_locus = Counter()
    false_class = Counter()
    false_region = Counter()
    all_gdna_region = Counter()
    window_counts = Counter()
    window_false = Counter()
    posterior_by_group: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    false_posteriors: list[float] = []
    false_template_lengths: list[int] = []
    false_samples: list[dict[str, Any]] = []
    false_qnames: set[str] = set()

    with pysam.AlignmentFile(str(args.bam), "rb", threads=args.threads) as bam:
        for read in bam.fetch(until_eof=True):
            records["records_seen"] += 1
            if read.is_secondary:
                records["secondary_skipped"] += 1
                continue
            if read.is_supplementary:
                records["supplementary_skipped"] += 1
                continue
            if read.is_paired and not read.is_read1:
                records["read2_skipped"] += 1
                continue
            if not is_representative(read):
                records["other_skipped"] += 1
                continue

            records["fragments_counted"] += 1
            flowcell = flowcell_from_qname(read.query_name)
            true_source = source_map.get(flowcell, "unknown")
            predicted_detail = predicted_detail_from_zf(get_tag(read, "ZF", None))
            predicted_binary = binary_pool(predicted_detail)
            zc = str(get_tag(read, "ZC", "missing"))
            zs = str(get_tag(read, "ZS", "missing"))
            zl = int(get_tag(read, "ZL", -1))
            zw = get_tag(read, "ZW", None)
            if zw is not None:
                zw = float(zw)
                posterior_by_group[(true_source, predicted_binary, predicted_detail)].append(zw)

            span_start, span_end = fragment_span(read)
            exact_mask, dominant_type, dominant_strand, overlap_bp = classify_region_overlap(
                region_arrays,
                span_start,
                span_end,
            )
            window_start = ((span_start - start0) // args.window_size) * args.window_size + start0
            window_start = max(start0, int(window_start))
            window_end = min(end0, window_start + args.window_size)

            flowcell_counts[(flowcell, true_source)] += 1
            confusion_binary[(true_source, predicted_binary)] += 1
            confusion_detail[(true_source, predicted_detail)] += 1

            if true_source == "gdna":
                all_gdna_region[(exact_mask, dominant_type, dominant_strand)] += 1
                window_counts[(window_start, window_end)] += 1

            if true_source == "gdna" and predicted_binary == "rna":
                false_qnames.add(read.query_name)
                target = target_label(read, predicted_detail)
                false_target[(target, predicted_detail)] += 1
                false_locus[(zl, predicted_detail)] += 1
                false_class[(predicted_detail, zc, zs, cigar_has_refskip(read), get_tag(read, "NH", ""))] += 1
                false_region[(exact_mask, dominant_type, dominant_strand)] += 1
                window_false[(window_start, window_end)] += 1
                if zw is not None:
                    false_posteriors.append(float(zw))
                if read.template_length:
                    false_template_lengths.append(abs(int(read.template_length)))
                if len(false_samples) < args.sample_limit:
                    false_samples.append(
                        {
                            "qname": read.query_name,
                            "chrom": read.reference_name,
                            "fragment_start_1based": span_start + 1,
                            "fragment_end": span_end,
                            "predicted_detail": predicted_detail,
                            "target": target,
                            "zl": zl,
                            "zw": zw if zw is not None else "",
                            "zc": zc,
                            "zs": zs,
                            "zf": get_tag(read, "ZF", ""),
                            "zt": get_tag(read, "ZT", ""),
                            "zg": get_tag(read, "ZG", ""),
                            "zr": get_tag(read, "ZR", ""),
                            "mapq": read.mapping_quality,
                            "nh": get_tag(read, "NH", ""),
                            "nm": get_tag(read, "NM", ""),
                            "cigar": read.cigarstring or "",
                            "template_length": read.template_length,
                            "exact_region_mask": exact_mask,
                            "dominant_region_type": dominant_type,
                            "dominant_region_strand": dominant_strand,
                            "region_overlap_bp": overlap_bp,
                        }
                    )

    return {
        "source_map": source_map,
        "records": records,
        "flowcell_counts": flowcell_counts,
        "confusion_binary": confusion_binary,
        "confusion_detail": confusion_detail,
        "false_target": false_target,
        "false_locus": false_locus,
        "false_class": false_class,
        "false_region": false_region,
        "all_gdna_region": all_gdna_region,
        "window_counts": window_counts,
        "window_false": window_false,
        "posterior_by_group": posterior_by_group,
        "false_posteriors": false_posteriors,
        "false_template_lengths": false_template_lengths,
        "false_samples": false_samples,
        "false_qnames": false_qnames,
    }


def write_false_rna_bam(
    input_bam: Path,
    output_bam: Path,
    false_qnames: set[str],
    threads: int,
) -> None:
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(input_bam), "rb", threads=threads) as src:
        with pysam.AlignmentFile(str(output_bam), "wb", template=src, threads=threads) as dst:
            for read in src.fetch(until_eof=True):
                if read.query_name in false_qnames:
                    dst.write(read)


def posterior_summary(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n": 0, "mean": None, "median": None, "q05": None, "q95": None, "ge_90": None}
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "q05": float(np.quantile(arr, 0.05)),
        "q95": float(np.quantile(arr, 0.95)),
        "ge_90": float((arr >= 0.90).mean()),
    }


def build_rows(scan: dict[str, Any], args: argparse.Namespace) -> dict[str, list[dict[str, Any]]]:
    rows: dict[str, list[dict[str, Any]]] = {}

    confusion_binary_rows = []
    for true_source in sorted({key[0] for key in scan["confusion_binary"]}):
        total = sum(scan["confusion_binary"].get((true_source, pred), 0) for pred in POOL_ORDER)
        for predicted_pool in POOL_ORDER:
            count = scan["confusion_binary"].get((true_source, predicted_pool), 0)
            confusion_binary_rows.append(
                {
                    "true_source": true_source,
                    "predicted_pool": predicted_pool,
                    "count": count,
                    "row_rate": count / total if total else 0.0,
                }
            )
    rows["confusion_binary"] = confusion_binary_rows

    confusion_detail_rows = []
    for true_source in sorted({key[0] for key in scan["confusion_detail"]}):
        total = sum(scan["confusion_detail"].get((true_source, pred), 0) for pred in DETAIL_ORDER)
        for predicted_detail in DETAIL_ORDER:
            count = scan["confusion_detail"].get((true_source, predicted_detail), 0)
            confusion_detail_rows.append(
                {
                    "true_source": true_source,
                    "predicted_detail": predicted_detail,
                    "count": count,
                    "row_rate": count / total if total else 0.0,
                }
            )
    rows["confusion_detail"] = confusion_detail_rows

    target_rows = []
    false_total = sum(scan["false_target"].values())
    for (target, predicted_detail), count in scan["false_target"].most_common():
        target_rows.append(
            {
                "target": target,
                "predicted_detail": predicted_detail,
                "gdna_false_rna": count,
                "fraction_of_false_rna": count / false_total if false_total else 0.0,
            }
        )
    rows["false_targets"] = target_rows

    locus_rows = []
    locus_totals: Counter[int] = Counter()
    for (locus_id, _predicted_detail), count in scan["false_locus"].items():
        locus_totals[locus_id] += count
    for locus_id, total in locus_totals.most_common():
        locus_rows.append(
            {
                "locus_id": locus_id,
                "gdna_false_rna": total,
                "pred_mrna": scan["false_locus"].get((locus_id, "mrna"), 0),
                "pred_nrna": scan["false_locus"].get((locus_id, "nrna"), 0),
                "fraction_of_false_rna": total / false_total if false_total else 0.0,
            }
        )
    rows["false_loci"] = locus_rows

    rows["false_classes"] = summarize_counter(
        scan["false_class"],
        ["predicted_detail", "zc", "zs", "cigar_has_refskip", "nh"],
    )

    region_rows = []
    for key, false_count in scan["false_region"].most_common():
        gdna_total = scan["all_gdna_region"].get(key, 0)
        exact_mask, dominant_type, dominant_strand = key
        region_rows.append(
            {
                "exact_region_mask": exact_mask,
                "dominant_region_type": dominant_type,
                "dominant_region_strand": dominant_strand,
                "gdna_source_total": gdna_total,
                "gdna_false_rna": false_count,
                "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
                "fraction_of_false_rna": false_count / false_total if false_total else 0.0,
            }
        )
    rows["false_regions"] = region_rows

    window_rows = []
    for (window_start, window_end), gdna_total in sorted(scan["window_counts"].items()):
        false_count = scan["window_false"].get((window_start, window_end), 0)
        if gdna_total == 0 and false_count == 0:
            continue
        window_rows.append(
            {
                "chrom": args.chrom,
                "window_start_1based": window_start + 1,
                "window_end": window_end,
                "gdna_source_total": gdna_total,
                "gdna_false_rna": false_count,
                "gdna_false_rna_rate": false_count / gdna_total if gdna_total else 0.0,
            }
        )
    rows["windows"] = sorted(window_rows, key=lambda row: row["gdna_false_rna"], reverse=True)

    posterior_rows = []
    for (true_source, predicted_binary, predicted_detail), values in sorted(
        scan["posterior_by_group"].items()
    ):
        summary = posterior_summary(values)
        posterior_rows.append(
            {
                "true_source": true_source,
                "predicted_pool": predicted_binary,
                "predicted_detail": predicted_detail,
                **summary,
            }
        )
    rows["posterior_summary"] = posterior_rows
    rows["false_samples"] = scan["false_samples"]
    return rows


def load_locus_context(quant_dir: Path) -> dict[str, Any]:
    loci = pd.read_feather(quant_dir / "loci.feather")
    quant = pd.read_feather(quant_dir / "quant.feather")
    nrna = pd.read_feather(quant_dir / "nrna_quant.feather")
    summary = json.loads((quant_dir / "summary.json").read_text())
    return {
        "top_loci": loci.sort_values("n_em_fragments", ascending=False).head(12).to_dict("records"),
        "top_quant": quant.sort_values("count", ascending=False).head(20).to_dict("records"),
        "top_nrna": nrna.sort_values("count", ascending=False).head(20).to_dict("records"),
        "summary": summary,
    }


def write_report(
    report: Path,
    rows: dict[str, list[dict[str, Any]]],
    scan: dict[str, Any],
    context: dict[str, Any],
    args: argparse.Namespace,
) -> None:
    gdna_total = sum(
        row["count"] for row in rows["confusion_binary"] if row["true_source"] == "gdna"
    )
    gdna_false_rna = sum(
        row["count"]
        for row in rows["confusion_binary"]
        if row["true_source"] == "gdna" and row["predicted_pool"] == "rna"
    )
    rna_total = sum(row["count"] for row in rows["confusion_binary"] if row["true_source"] == "rna")
    rna_false_gdna = sum(
        row["count"]
        for row in rows["confusion_binary"]
        if row["true_source"] == "rna" and row["predicted_pool"] == "gdna"
    )
    false_posterior = posterior_summary(scan["false_posteriors"])
    template_summary = posterior_summary([float(v) for v in scan["false_template_lengths"]])
    calibration = context["summary"].get("calibration", {})
    regional_exposure = calibration.get("regional_exposure", {})

    lines = [
        "# VCaP AR-Locus Subset Analysis",
        "",
        f"Region: `{args.chrom}:{args.start}-{args.end}`",
        "",
        "## Summary",
        "",
        f"- Representative fragments counted: {scan['records']['fragments_counted']:,}",
        f"- True gDNA fragments called RNA: {gdna_false_rna:,}/{gdna_total:,} "
        f"({fmt_pct(gdna_false_rna / gdna_total if gdna_total else None)})",
        f"- True RNA fragments called gDNA: {rna_false_gdna:,}/{rna_total:,} "
        f"({fmt_pct(rna_false_gdna / rna_total if rna_total else None)})",
        f"- False-RNA posterior median: {false_posterior['median']:.3f}",
        f"- False-RNA posterior >= 0.90: {fmt_pct(false_posterior['ge_90'])}",
        f"- False-RNA template-length median: {template_summary['median']:.1f} bp",
        f"- Local regional exposure mode: `{regional_exposure.get('mode', 'missing')}`",
        f"- Local mean_pi_gdna: `{calibration.get('mean_pi_gdna', 'missing')}`",
        "",
        "## Binary Confusion",
        "",
        markdown_table(rows["confusion_binary"]),
        "",
        "## Detail Confusion",
        "",
        markdown_table(rows["confusion_detail"]),
        "",
        "## Top False-RNA Targets",
        "",
        markdown_table(rows["false_targets"][:20]),
        "",
        "## Top False-RNA Windows",
        "",
        markdown_table(rows["windows"][:20]),
        "",
        "## False-RNA Region Classes",
        "",
        markdown_table(rows["false_regions"]),
        "",
        "## False-RNA Fragment Classes",
        "",
        markdown_table(rows["false_classes"][:20]),
        "",
        "## Top Local Loci",
        "",
        markdown_table(context["top_loci"]),
    ]
    report.write_text("\n".join(lines) + "\n")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scan = scan_bam(args)
    if args.false_rna_bam is not None:
        write_false_rna_bam(args.bam, args.false_rna_bam, scan["false_qnames"], args.threads)
    rows = build_rows(scan, args)
    context = load_locus_context(args.quant_dir)

    write_tsv(
        args.out_dir / "confusion_binary.tsv",
        rows["confusion_binary"],
        ["true_source", "predicted_pool", "count", "row_rate"],
    )
    write_tsv(
        args.out_dir / "confusion_detail.tsv",
        rows["confusion_detail"],
        ["true_source", "predicted_detail", "count", "row_rate"],
    )
    write_tsv(
        args.out_dir / "false_targets.tsv",
        rows["false_targets"],
        ["target", "predicted_detail", "gdna_false_rna", "fraction_of_false_rna"],
    )
    write_tsv(
        args.out_dir / "false_loci.tsv",
        rows["false_loci"],
        ["locus_id", "gdna_false_rna", "pred_mrna", "pred_nrna", "fraction_of_false_rna"],
    )
    write_tsv(
        args.out_dir / "false_classes.tsv",
        rows["false_classes"],
        ["predicted_detail", "zc", "zs", "cigar_has_refskip", "nh", "count"],
    )
    write_tsv(
        args.out_dir / "false_regions.tsv",
        rows["false_regions"],
        [
            "exact_region_mask",
            "dominant_region_type",
            "dominant_region_strand",
            "gdna_source_total",
            "gdna_false_rna",
            "gdna_false_rna_rate",
            "fraction_of_false_rna",
        ],
    )
    write_tsv(
        args.out_dir / "false_windows.tsv",
        rows["windows"],
        ["chrom", "window_start_1based", "window_end", "gdna_source_total", "gdna_false_rna", "gdna_false_rna_rate"],
    )
    write_tsv(
        args.out_dir / "posterior_summary.tsv",
        rows["posterior_summary"],
        ["true_source", "predicted_pool", "predicted_detail", "n", "mean", "median", "q05", "q95", "ge_90"],
    )
    write_tsv(args.out_dir / "false_samples.tsv", rows["false_samples"], list(rows["false_samples"][0]) if rows["false_samples"] else [])

    summary = {
        "bam": str(args.bam),
        "region": {"chrom": args.chrom, "start_1based": args.start, "end": args.end},
        "source_map": scan["source_map"],
        "records": dict(scan["records"]),
        "flowcell_counts": {"|".join(map(str, key)): value for key, value in scan["flowcell_counts"].items()},
        "false_posterior": posterior_summary(scan["false_posteriors"]),
        "false_template_length": posterior_summary([float(v) for v in scan["false_template_lengths"]]),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_report(args.report, rows, scan, context, args)


if __name__ == "__main__":
    main()