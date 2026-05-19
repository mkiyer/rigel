#!/usr/bin/env python3
"""Compare VCaP RNA/gDNA truth against Rigel annotated BAM pool calls."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


BASE = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")
RNA_FLOWCELL = "C6EL5ANXX"
GDNA_FLOWCELL = "H7MFFDSXY"
POOL_COLUMNS = ("mrna", "nrna", "gdna", "unresolved")

DEFAULT_RUNS = {
    "legacy_with_mm": BASE / "with_mm" / "annotated.bam",
    "updated_uniform_off": BASE / "regional_off_v3_with_mm" / "annotated.bam",
    "updated_regional_auto": BASE / "regional_auto_v3_with_mm" / "annotated.bam",
}


def truth_source(query_name: str) -> str:
    fields = query_name.split(":")
    flowcell = fields[2] if len(fields) > 2 else ""
    if flowcell == GDNA_FLOWCELL or GDNA_FLOWCELL in query_name:
        return "gdna"
    if flowcell == RNA_FLOWCELL or RNA_FLOWCELL in query_name:
        return "rna"
    return "unknown"


def predicted_pool(read: pysam.AlignedSegment) -> str:
    try:
        zf = int(read.get_tag("ZF"))
    except KeyError:
        zf = 0

    # Priority matches earlier diagnostics and Rigel's pool bit semantics.
    if zf & 0x04:
        return "gdna"
    if zf & 0x08:
        return "nrna"
    if zf & 0x02:
        return "mrna"
    return "unresolved"


def optional_tag(read: pysam.AlignedSegment, tag: str, default: object = "") -> object:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def scan_bam(label: str, bam_path: Path) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    counts: Counter[tuple[str, str]] = Counter()
    gdna_false_rna_by_zc: Counter[str] = Counter()
    gdna_false_rna_by_zs: Counter[str] = Counter()
    gdna_false_rna_zw: list[float] = []
    true_totals: Counter[str] = Counter()
    primary_read1 = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            primary_read1 += 1
            truth = truth_source(read.query_name)
            pred = predicted_pool(read)
            counts[(truth, pred)] += 1
            true_totals[truth] += 1
            if truth == "gdna" and pred in {"mrna", "nrna"}:
                gdna_false_rna_by_zc[str(optional_tag(read, "ZC", "missing"))] += 1
                gdna_false_rna_by_zs[str(optional_tag(read, "ZS", "missing"))] += 1
                zw = optional_tag(read, "ZW", np.nan)
                try:
                    gdna_false_rna_zw.append(float(zw))
                except (TypeError, ValueError):
                    pass

    rows = []
    for truth in ("gdna", "rna", "unknown"):
        for pred in POOL_COLUMNS:
            count = counts[(truth, pred)]
            total = true_totals[truth]
            rows.append({
                "run": label,
                "true_source": truth,
                "predicted_pool": pred,
                "count": count,
                "row_fraction": count / total if total else 0.0,
            })

    details = []
    for key, value in gdna_false_rna_by_zc.most_common():
        details.append({"run": label, "dimension": "ZC", "value": key, "count": value})
    for key, value in gdna_false_rna_by_zs.most_common():
        details.append({"run": label, "dimension": "ZS", "value": key, "count": value})
    if gdna_false_rna_zw:
        zw_arr = np.asarray(gdna_false_rna_zw, dtype=np.float64)
        for name, value in {
            "n": float(len(zw_arr)),
            "mean": float(np.mean(zw_arr)),
            "median": float(np.median(zw_arr)),
            "p05": float(np.quantile(zw_arr, 0.05)),
            "p95": float(np.quantile(zw_arr, 0.95)),
        }.items():
            details.append({"run": label, "dimension": "ZW", "value": name, "count": value})
    details.append({"run": label, "dimension": "PRIMARY_READ1", "value": "total", "count": primary_read1})
    return rows, details


def summarize(confusion: pd.DataFrame) -> pd.DataFrame:
    records = []
    for run, grp in confusion.groupby("run"):
        matrix = grp.pivot_table(
            index="true_source",
            columns="predicted_pool",
            values="count",
            aggfunc="sum",
            fill_value=0,
        )
        for col in POOL_COLUMNS:
            if col not in matrix.columns:
                matrix[col] = 0
        gdna_total = float(matrix.loc["gdna", POOL_COLUMNS].sum()) if "gdna" in matrix.index else 0.0
        rna_total = float(matrix.loc["rna", POOL_COLUMNS].sum()) if "rna" in matrix.index else 0.0
        gdna_to_mrna = float(matrix.loc["gdna", "mrna"]) if gdna_total else 0.0
        gdna_to_nrna = float(matrix.loc["gdna", "nrna"]) if gdna_total else 0.0
        rna_to_gdna = float(matrix.loc["rna", "gdna"]) if rna_total else 0.0
        records.append({
            "run": run,
            "gdna_total": gdna_total,
            "rna_total": rna_total,
            "gdna_to_mrna": gdna_to_mrna,
            "gdna_to_nrna": gdna_to_nrna,
            "gdna_to_rna": gdna_to_mrna + gdna_to_nrna,
            "gdna_to_rna_rate": (gdna_to_mrna + gdna_to_nrna) / gdna_total if gdna_total else 0.0,
            "gdna_correct": float(matrix.loc["gdna", "gdna"]) if gdna_total else 0.0,
            "gdna_correct_rate": float(matrix.loc["gdna", "gdna"]) / gdna_total if gdna_total else 0.0,
            "rna_to_gdna": rna_to_gdna,
            "rna_to_gdna_rate": rna_to_gdna / rna_total if rna_total else 0.0,
            "rna_to_mrna": float(matrix.loc["rna", "mrna"]) if rna_total else 0.0,
            "rna_to_nrna": float(matrix.loc["rna", "nrna"]) if rna_total else 0.0,
        })
    return pd.DataFrame(records)


def markdown_table(frame: pd.DataFrame) -> str:
    def cell(value: object) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, float):
            return f"{value:.6g}"
        return str(value)

    cols = [str(c) for c in frame.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(cell(row[col]) for col in frame.columns) + " |")
    return "\n".join(lines)


def confusion_matrix(confusion: pd.DataFrame, run: str) -> pd.DataFrame:
    sub = confusion[confusion["run"] == run]
    matrix = sub.pivot_table(
        index="true_source",
        columns="predicted_pool",
        values="count",
        aggfunc="sum",
        fill_value=0,
    )
    for col in POOL_COLUMNS:
        if col not in matrix.columns:
            matrix[col] = 0
    return matrix[list(POOL_COLUMNS)].reset_index()


def load_exposure_summary(run_dir: Path) -> dict[str, object]:
    summary_path = run_dir / "summary.json"
    if not summary_path.exists():
        return {}
    summary = json.loads(summary_path.read_text())
    cal = summary.get("calibration", {})
    exp = cal.get("regional_exposure", {})
    app = cal.get("regional_weighting_stats", {})
    gdna = summary.get("gdna_eff_len", {}).get("value", {})
    return {
        "mode": exp.get("mode", "none"),
        "rho_ref": exp.get("rho_ref"),
        "n_at_floor": exp.get("n_at_floor"),
        "weighted_units": app.get("n_units_weighted"),
        "seen_units": app.get("n_units_seen"),
        "cross_ref_skipped": app.get("n_units_cross_ref_skipped"),
        "gdna_eff_len_median": gdna.get("median"),
        "gdna_eff_len_p95": gdna.get("p95"),
        "gdna_eff_len_max": gdna.get("max"),
    }


def write_report(
    out_dir: Path,
    confusion: pd.DataFrame,
    summary: pd.DataFrame,
    details: pd.DataFrame,
    run_paths: dict[str, Path],
) -> Path:
    ordered_runs = [run for run in run_paths if run in set(summary["run"])]
    before = "updated_uniform_off" if "updated_uniform_off" in ordered_runs else ordered_runs[0]
    after = "updated_regional_auto" if "updated_regional_auto" in ordered_runs else ordered_runs[-1]
    before_row = summary.set_index("run").loc[before]
    after_row = summary.set_index("run").loc[after]
    delta_rate = float(after_row["gdna_to_rna_rate"] - before_row["gdna_to_rna_rate"])
    delta_count = int(after_row["gdna_to_rna"] - before_row["gdna_to_rna"])
    rel = delta_count / float(before_row["gdna_to_rna"]) if before_row["gdna_to_rna"] else 0.0

    exposure_rows = []
    for run, bam_path in run_paths.items():
        exposure = load_exposure_summary(bam_path.parent)
        if exposure:
            exposure_rows.append({"run": run, **exposure})
    exposure_df = pd.DataFrame(exposure_rows)

    lines = [
        "# VCaP Regional Exposure Before/After Confusion",
        "",
        f"Input truth: `{RNA_FLOWCELL}` = RNA, `{GDNA_FLOWCELL}` = gDNA.",
        "Primary read1 only; secondary and supplementary records skipped.",
        "",
        "## Key Result",
        "",
        f"- Before/control: `{before}`",
        f"- After/regional: `{after}`",
        f"- gDNA -> RNA changed by {delta_count:+,} fragments ({delta_rate:+.6f} rate; {rel:+.2%} relative).",
        f"- gDNA -> RNA before: {int(before_row['gdna_to_rna']):,} / {int(before_row['gdna_total']):,} ({before_row['gdna_to_rna_rate']:.4%}).",
        f"- gDNA -> RNA after: {int(after_row['gdna_to_rna']):,} / {int(after_row['gdna_total']):,} ({after_row['gdna_to_rna_rate']:.4%}).",
        "",
        "## Summary Metrics",
        "",
        markdown_table(summary),
        "",
        "## Exposure Summary",
        "",
        markdown_table(exposure_df) if not exposure_df.empty else "No exposure summaries found.",
        "",
    ]

    for run in ordered_runs:
        lines.extend([
            f"## Confusion Matrix - {run}",
            "",
            markdown_table(confusion_matrix(confusion, run)),
            "",
        ])

    focus = details[details["run"].isin([before, after])].copy()
    if not focus.empty:
        lines.extend([
            "## gDNA -> RNA Detail Tags",
            "",
            markdown_table(focus),
            "",
        ])

    report_path = out_dir / "vcap_regional_v3_confusion_report.md"
    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=BASE / "regional_v3_confusion")
    parser.add_argument(
        "--runs",
        nargs="*",
        default=None,
        help="Optional label=path entries. Default compares legacy, updated off, and updated auto.",
    )
    args = parser.parse_args()

    if args.runs:
        run_paths = {}
        for item in args.runs:
            label, path = item.split("=", 1)
            run_paths[label] = Path(path)
    else:
        run_paths = DEFAULT_RUNS

    args.out_dir.mkdir(parents=True, exist_ok=True)
    all_rows = []
    all_details = []
    for label, bam_path in run_paths.items():
        print(f"scan {label}: {bam_path}", flush=True)
        rows, details = scan_bam(label, bam_path)
        all_rows.extend(rows)
        all_details.extend(details)

    confusion = pd.DataFrame(all_rows)
    details = pd.DataFrame(all_details)
    summary = summarize(confusion)
    confusion.to_csv(args.out_dir / "confusion_long.tsv", sep="\t", index=False)
    details.to_csv(args.out_dir / "gdna_false_rna_details.tsv", sep="\t", index=False)
    summary.to_csv(args.out_dir / "summary.tsv", sep="\t", index=False)
    report_path = write_report(args.out_dir, confusion, summary, details, run_paths)
    print(f"Report: {report_path}")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())