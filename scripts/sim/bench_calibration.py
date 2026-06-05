#!/usr/bin/env python3
"""Calibration accuracy benchmark for the synthetic hybrid-capture suite.

Evaluates how well Rigel's calibration + EM classify fragments into gDNA vs RNA,
against the simulated ground truth, across a suite of conditions
(gDNA level x strand-specificity x capture on/off).

Two complementary metrics, because they can disagree:

  * POOL (fractional)  -- the library-wide gDNA fraction from the EM abundances
    (summary.json `quantification`). This is what calibration is *for*.
  * PER-FRAGMENT (hard) -- the hard gDNA/RNA label each fragment received
    (annotated-BAM ZF tag) vs its true origin (oracle read name). Splits the
    gDNA->RNA leak (the FP-deleterious direction) from the RNA->gDNA siphon
    (the FP-safe direction), and separates IN-LOCUS gDNA (on captured exons,
    the hard case) from INTERGENIC gDNA (the easy case).

Reuses the run + fragment-decode machinery from rigel.sim.analysis, so it stays
in sync with the canonical harness. The stale analyze_calibration() in that
module assumes the pre-rebuild summary schema; this script reads the current one.

Usage:
    # evaluate existing rigel_out/ (must already be quantified):
    python scripts/sim/bench_calibration.py --sim-base <suite_dir>
    # (re)run rigel quant first, then evaluate:
    python scripts/sim/bench_calibration.py --sim-base <suite_dir> --run [--force]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
import pysam

from rigel.sim.analysis import build_index, parse_origin, run_quant

_GDNA_BIT = 0x04  # annotated-BAM ZF: gDNA-assigned bit (AF_GDNA_EM=0x05, AF_GDNA_INTERGENIC=0x25)


def _conditions(sim_base: Path) -> list[dict]:
    man = json.loads((sim_base / "manifest.json").read_text())
    return man["conditions"]


def _pool_row(sim_base: Path, cond: dict) -> dict:
    """Pool-level fractional gDNA classification from summary.json quantification."""
    name = cond["name"]
    q = json.loads((sim_base / name / "rigel_out" / "summary.json").read_text())["quantification"]
    # genomic = locus gDNA + intergenic (intergenic carries no transcript -> gDNA by nature)
    called_gdna = float(q.get("gdna_total", 0.0)) + float(q.get("intergenic_total", 0.0))
    called_rna = float(q.get("mrna_total", 0.0)) + float(q.get("nrna_total", 0.0))
    n_gdna = cond["n_gdna_observed"]
    n_rna = cond["n_mrna_observed"] + cond["n_nrna_observed"]
    truth_gf = n_gdna / (n_gdna + n_rna) if (n_gdna + n_rna) else 0.0
    called_gf = called_gdna / (called_gdna + called_rna) if (called_gdna + called_rna) else 0.0
    return {
        "truth_gdna_frac": truth_gf,
        "pool_called_gdna_frac": called_gf,
        "pool_gdna_frac_err": called_gf - truth_gf,
        "called_gdna": called_gdna,
        "called_rna": called_rna,
    }


def _fragment_confusion(sim_base: Path, cond: dict) -> dict:
    """Per-fragment hard gDNA/RNA confusion from annotated BAM tags vs oracle truth.

    Splits leaked gDNA by IN-LOCUS (ZL >= 0, on a gene -> capture-enriched, hard)
    vs INTERGENIC (ZL < 0, easy), which localises the capture failure mode.
    """
    bam_path = sim_base / cond["name"] / "annotated.bam"
    g_caught = g_leak_inlocus = g_leak_inter = 0
    r_total = r_siphon = r_tx_correct = 0
    fl_leak: list[int] = []
    fl_caught: list[int] = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            origin = parse_origin(read.query_name)
            try:
                zf = int(read.get_tag("ZF"))
            except KeyError:
                zf = 0
            try:
                zl = int(read.get_tag("ZL"))
            except KeyError:
                zl = -1
            try:
                category = read.get_tag("ZC")
            except KeyError:
                category = ""
            is_gdna = bool(zf & _GDNA_BIT) or category == "intergenic"
            fl = abs(read.template_length) or (read.query_length or 0)
            if origin.kind == "gdna":
                if is_gdna:
                    g_caught += 1
                    fl_caught.append(fl)
                else:
                    fl_leak.append(fl)
                    if zl is not None and zl >= 0:
                        g_leak_inlocus += 1
                    else:
                        g_leak_inter += 1
            else:  # mrna / nrna
                r_total += 1
                if is_gdna:
                    r_siphon += 1
                else:
                    try:
                        assigned_tx = read.get_tag("ZT")
                    except KeyError:
                        assigned_tx = ""
                    if assigned_tx and assigned_tx == (origin.transcript_id or ""):
                        r_tx_correct += 1
    g_total = g_caught + g_leak_inlocus + g_leak_inter
    g_leak = g_leak_inlocus + g_leak_inter

    def _mean(x: list[int]) -> float:
        return round(sum(x) / len(x), 1) if x else 0.0

    return {
        "n_gdna_frag": g_total,
        "n_rna_frag": r_total,
        "gdna_to_rna_leak": g_leak,
        "gdna_leak_frac": (g_leak / g_total) if g_total else 0.0,
        "gdna_leak_inlocus": g_leak_inlocus,
        "gdna_leak_intergenic": g_leak_inter,
        "rna_to_gdna_siphon": r_siphon,
        "rna_siphon_frac": (r_siphon / r_total) if r_total else 0.0,
        "rna_tx_correct_frac": (r_tx_correct / r_total) if r_total else 0.0,
        "fl_leaked_gdna_mean": _mean(fl_leak),
        "fl_caught_gdna_mean": _mean(fl_caught),
    }


def evaluate(sim_base: Path) -> pd.DataFrame:
    rows = []
    for cond in _conditions(sim_base):
        out = sim_base / cond["name"] / "rigel_out" / "summary.json"
        if not out.exists():
            print(f"  [skip] {cond['name']}: no rigel_out (run with --run)")
            continue
        row = {
            "condition": cond["name"],
            "gdna": cond["gdna_label"],
            "capture": cond["capture_label"],
            "ss": cond["strand_specificity"],
        }
        row.update(_pool_row(sim_base, cond))
        row.update(_fragment_confusion(sim_base, cond))
        rows.append(row)
    df = pd.DataFrame(rows)
    if not df.empty:
        order = {"none": 0, "high": 1}
        df = df.sort_values(
            by=["gdna", "capture", "ss"],
            key=lambda s: s.map(order) if s.name == "gdna" else s,
            ascending=[True, True, False],
        ).reset_index(drop=True)
    return df


def _print_report(df: pd.DataFrame) -> str:
    lines = ["", "=" * 100, "  RIGEL CALIBRATION BENCHMARK — gDNA vs RNA classification", "=" * 100]
    lines.append("\nPOOL-LEVEL (fractional EM gDNA fraction vs truth):")
    lines.append(f"  {'gdna':5}{'cap':4}{'ss':6} | {'truth%':>7} {'called%':>8} {'err%':>7}")
    lines.append("  " + "-" * 44)
    for _, r in df.iterrows():
        lines.append(
            f"  {r.gdna:5}{r.capture:4}{str(r.ss):6} | {r.truth_gdna_frac * 100:6.1f}% "
            f"{r.pool_called_gdna_frac * 100:7.1f}% {r.pool_gdna_frac_err * 100:+6.1f}%"
        )
    lines.append(
        "\nPER-FRAGMENT (hard label) — gDNA->RNA leak [FP-deleterious] & RNA->gDNA siphon [safe]:"
    )
    lines.append(
        f"  {'gdna':5}{'cap':4}{'ss':6} | {'leak%':>6} {'(in-locus':>10} {'interg)':>8} | "
        f"{'siphon%':>7} | {'tx-ok%':>7} | {'FL leak':>7} {'FL caught':>9}"
    )
    lines.append("  " + "-" * 86)
    for _, r in df.iterrows():
        lines.append(
            f"  {r.gdna:5}{r.capture:4}{str(r.ss):6} | {r.gdna_leak_frac * 100:5.1f}% "
            f"{r.gdna_leak_inlocus:>10d} {r.gdna_leak_intergenic:>8d} | "
            f"{r.rna_siphon_frac * 100:6.1f}% | {r.rna_tx_correct_frac * 100:6.1f}% | "
            f"{r.fl_leaked_gdna_mean:>7.0f} {r.fl_caught_gdna_mean:>9.0f}"
        )
    return "\n".join(lines)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--sim-base", type=Path, required=True, help="Suite dir (has manifest.json + conditions)"
    )
    ap.add_argument(
        "--run", action="store_true", help="Build index + run rigel quant before evaluating"
    )
    ap.add_argument(
        "--force", action="store_true", help="Re-quant even if rigel_out exists (with --run)"
    )
    args = ap.parse_args()
    sim_base = args.sim_base
    assert (sim_base / "manifest.json").exists(), f"no manifest.json in {sim_base}"

    if args.run:
        index_dir = build_index(sim_base)
        for cond in _conditions(sim_base):
            if args.force:
                qf = sim_base / cond["name"] / "rigel_out" / "quant.feather"
                if qf.exists():
                    qf.unlink()
            run_quant(sim_base, index_dir, cond["name"])

    df = evaluate(sim_base)
    report = _print_report(df)
    print(report)
    metrics_path = sim_base / "bench_calibration_metrics.tsv"
    report_path = sim_base / "bench_calibration_report.txt"
    df.to_csv(metrics_path, sep="\t", index=False)
    report_path.write_text(report + "\n")
    print(f"\n  metrics: {metrics_path}\n  report:  {report_path}")


if __name__ == "__main__":
    main()
