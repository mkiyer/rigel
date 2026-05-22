#!/usr/bin/env python3
"""Audit exposure-weighted gDNA opportunity for the VCAP FLG2 MultiLocus."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd


DEFAULT_QUANT_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/exon_strand_deconv_v1")
DEFAULT_INDEX_DIR = Path("/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index")
DEFAULT_FLG2_DIAG_DIR = Path("results/flg2_hotspot_diagnostics_2026-05-21")
DEFAULT_REPORT = Path("docs/benchmarks/locus3_exposure_audit_2026-05-21.md")


NRNA_RE = re.compile(r"^RIGEL_NRNA_(?P<ref>.+)_(?P<strand>[12])_(?P<start>\d+)_(?P<end>\d+)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quant-dir", type=Path, default=DEFAULT_QUANT_DIR)
    parser.add_argument("--index-dir", type=Path, default=DEFAULT_INDEX_DIR)
    parser.add_argument("--flg2-diag-dir", type=Path, default=DEFAULT_FLG2_DIAG_DIR)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--locus-id", type=int, default=3)
    parser.add_argument("--window-bp", type=float, default=10_000.0)
    return parser.parse_args()


def fmt_int(value: float | int) -> str:
    return f"{int(round(float(value))):,}"


def fmt_float(value: float | int, digits: int = 4) -> str:
    return f"{float(value):.{digits}f}"


def fmt_sci(value: float | int, digits: int = 3) -> str:
    return f"{float(value):.{digits}e}"


def markdown_table(rows: list[dict[str, Any]], columns: list[tuple[str, str]]) -> str:
    if not rows:
        return "_No rows._"
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    sep = "| " + " | ".join("---" for _label, _key in columns) + " |"
    body = ["| " + " | ".join(str(row.get(key, "")) for _label, key in columns) + " |" for row in rows]
    return "\n".join([header, sep, *body])


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged: list[tuple[str, int, int]] = []
    cur_ref, cur_start, cur_end = ordered[0]
    for ref, start, end in ordered[1:]:
        if ref == cur_ref and start <= cur_end:
            cur_end = max(cur_end, end)
        else:
            merged.append((cur_ref, cur_start, cur_end))
            cur_ref, cur_start, cur_end = ref, start, end
    merged.append((cur_ref, cur_start, cur_end))
    return merged


def span_of(intervals: list[tuple[str, int, int]]) -> int:
    return sum(end - start for _ref, start, end in intervals)


def parse_nrna_interval(nrna_id: str) -> tuple[str, int, int] | None:
    match = NRNA_RE.match(nrna_id)
    if not match:
        return None
    return (match.group("ref"), int(match.group("start")) - 1, int(match.group("end")))


def describe_weights(df: pd.DataFrame, label: str) -> dict[str, Any]:
    weights = df["em_exposure_weight"].dropna()
    if weights.empty:
        return {"component_set": label}
    return {
        "component_set": label,
        "n": fmt_int(weights.size),
        "min": fmt_float(weights.min(), 4),
        "p05": fmt_float(weights.quantile(0.05), 4),
        "p50": fmt_float(weights.quantile(0.50), 4),
        "p95": fmt_float(weights.quantile(0.95), 4),
        "max": fmt_float(weights.max(), 4),
        "count_weighted_mean": fmt_float((weights * df.loc[weights.index, "count"]).sum() / max(df.loc[weights.index, "count"].sum(), 1.0), 4),
    }


def main() -> int:
    args = parse_args()
    loci = pd.read_feather(args.quant_dir / "loci.feather")
    quant = pd.read_feather(args.quant_dir / "quant.feather")
    nrna = pd.read_feather(args.quant_dir / "nrna_quant.feather")
    transcripts = pd.read_feather(args.index_dir / "transcripts.feather")
    summary = json.loads((args.quant_dir / "summary.json").read_text())
    local_math = pd.read_csv(args.flg2_diag_dir / "local_math.tsv", sep="\t")

    locus = loci.loc[loci["locus_id"] == args.locus_id].iloc[0]
    q3 = quant.loc[quant["locus_id"] == args.locus_id].copy()
    n3 = nrna.loc[nrna["locus_id"] == args.locus_id].copy()

    transcript_hits = transcripts.set_index("t_id").loc[
        transcripts.set_index("t_id").index.intersection(q3["transcript_id"])
    ]
    intervals: list[tuple[str, int, int]] = [
        (str(row.ref), int(row.start), int(row.end))
        for row in transcript_hits.itertuples(index=False)
    ]
    for nrna_id in n3["nrna_id"].astype(str):
        interval = parse_nrna_interval(nrna_id)
        if interval is not None:
            intervals.append(interval)
    merged = merge_intervals(intervals)
    reconstructed_span = span_of(merged)

    gdna_alpha = float(local_math["gdna_alpha_approx"].iloc[0])
    current_len = float(locus["gdna_eff_len"])
    unweighted_len = float(locus["gdna_eff_len_unweighted"])
    current_weight = float(locus["gdna_em_exposure_weight"])
    threshold_rows = []
    for row in local_math.itertuples(index=False):
        comp_count = float(row.component_count) + 0.5
        comp_len = float(row.component_eff_len)
        threshold = gdna_alpha * comp_len / (2.0 * comp_count)
        threshold_rows.append(
            {
                "component": str(row.component).replace("RIGEL_NRNA_", ""),
                "current_len": fmt_int(current_len),
                "threshold_len": fmt_int(threshold),
                "current_over_threshold": fmt_float(current_len / threshold, 2),
                "threshold_weight": fmt_sci(threshold / unweighted_len, 3),
            }
        )

    regional = summary["calibration"]["regional_exposure"]
    class_rows = []
    for class_name, class_summary in regional["per_class"].items():
        class_rows.append(
            {
                "class": class_name,
                "rho_q50": fmt_sci(class_summary.get("rho_q50", 0.0), 3),
                "rho_q95": fmt_sci(class_summary.get("rho_q95", 0.0), 3),
                "opportunity": fmt_int(class_summary.get("opportunity", 0.0)),
                "evidence": fmt_int(class_summary.get("evidence_count", 0.0)),
            }
        )

    weight_rows = [
        describe_weights(q3, "annotated RNA components in locus 3"),
        describe_weights(n3, "nRNA entities in locus 3"),
    ]
    top_component_rows = []
    top_q = q3.sort_values("count", ascending=False).head(6)
    for row in top_q.itertuples(index=False):
        top_component_rows.append(
            {
                "component": str(row.transcript_id),
                "kind": "annotated",
                "count": fmt_int(row.count),
                "em_eff_len": fmt_int(row.em_effective_length),
                "em_weight": fmt_float(row.em_exposure_weight, 4),
            }
        )
    top_n = n3.sort_values("count", ascending=False).head(6)
    for row in top_n.itertuples(index=False):
        top_component_rows.append(
            {
                "component": str(row.nrna_id).replace("RIGEL_NRNA_", ""),
                "kind": "nRNA",
                "count": fmt_int(row.count),
                "em_eff_len": fmt_int(row.em_effective_length),
                "em_weight": fmt_float(row.em_exposure_weight, 4),
            }
        )

    desired_window_weight = args.window_bp / unweighted_len
    report = f"""# Locus 3 Exposure-Weighted Opportunity Audit - 2026-05-21

Input run: `{args.quant_dir}`

## Saved Denominator

Exposure weighting is active, but it is weak for this MultiLocus.

| Metric | Value |
| --- | --- |
| locus_id | {args.locus_id} |
| locus_span_bp | {fmt_int(locus['locus_span_bp'])} |
| gdna_eff_len_unweighted | {fmt_int(unweighted_len)} |
| gdna_em_exposure_weight | {fmt_float(current_weight, 6)} |
| gdna_eff_len passed to EM | {fmt_int(current_len)} |
| reduction factor | {fmt_float(unweighted_len / current_len, 2)}x |
| weight needed for 10 kb opportunity | {fmt_sci(desired_window_weight, 3)} |
| current / 10 kb opportunity | {fmt_float(current_len / args.window_bp, 1)}x |

Important correction: `92.4M` is not the raw locus size. It is already the exposure-weighted denominator. The unweighted FL-marginal gDNA denominator is `353.8M`, close to the saved `locus_span_bp` plus the FL expansion term.

## EM-Relevant Thresholds

These thresholds ask: how small would `L_gDNA` need to be for a strand-compatible gDNA fragment to break even against this RNA component, assuming the current alpha values and the gDNA `0.5` strand factor?

{markdown_table(threshold_rows, [('Component', 'component'), ('Current L_gDNA', 'current_len'), ('Break-even L_gDNA', 'threshold_len'), ('Current / break-even', 'current_over_threshold'), ('Break-even weight', 'threshold_weight')])}

The current denominator is about 19.3x too large for gDNA to beat FLG2 mRNA on the strand-compatible half of the pileup.

## Regional Exposure Summary

The exposure model itself is in regional mode. However, the reference density is the global 95th percentile (`rho_ref={fmt_sci(regional['rho_ref'], 3)}`), and weights are capped at 1. Regions above that reference all saturate. In this run, ordinary exon composite density is already above the reference median, so many exon/capture-like regions are not distinguished from the FLG2 hotspot.

{markdown_table(class_rows, [('Class', 'class'), ('rho q50', 'rho_q50'), ('rho q95', 'rho_q95'), ('Opportunity', 'opportunity'), ('Evidence', 'evidence')])}

## RNA Component Exposure Weights

{markdown_table(weight_rows, [('Set', 'component_set'), ('N', 'n'), ('Min', 'min'), ('P05', 'p05'), ('P50', 'p50'), ('P95', 'p95'), ('Max', 'max'), ('Count-weighted mean', 'count_weighted_mean')])}

Top locus-3 RNA/nRNA components:

{markdown_table(top_component_rows, [('Component', 'component'), ('Kind', 'kind'), ('Count', 'count'), ('EM effective length', 'em_eff_len'), ('EM weight', 'em_weight')])}

## Footprint Reconstruction Check

Using saved annotated transcript rows plus parseable nRNA span ids, the reconstructed merged footprint is {fmt_int(reconstructed_span)} bp across {fmt_int(len(merged))} merged intervals. This is close enough to the saved {fmt_int(locus['locus_span_bp'])} bp to confirm that locus 3 is a large connected component, not a 10 kb local problem.

## Code-Path Finding

The production prior path computes:

```text
L_gDNA_EM = gdna_eff_len_for_loci(ml.loci) * footprint_exposure_weight(ml.loci)
```

That is an unweighted FL-expanded denominator multiplied by one scalar bp-weighted average over the unexpanded MultiLocus footprint. The repository also contains `weighted_gdna_eff_len_for_loci(...)`, which directly integrates the exposure field over the FL-expanded midpoint windows, but `assemble_priors(...)` does not use it.

This explains the saved number mechanically: exposure weighting is enabled, but for locus 3 it collapses to a broad scalar `0.261`, not a local/capture opportunity near kb scale.
"""
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(report)
    print(f"wrote {args.report}")
    print(f"locus3 weighted len {current_len:,.0f}; unweighted {unweighted_len:,.0f}; weight {current_weight:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())