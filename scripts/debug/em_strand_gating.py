#!/usr/bin/env python
"""Phase-0 gating measurement (em_strand): is the co-stranded gDNA leak an ABUNDANCE bias,
or a hard-label artifact with correct soft quantification?

Joins Rigel's soft per-transcript estimate (`quant.tsv` count_em) against simulated truth
(`truth_abundances.tsv` observed_mrna_fragments), and asks:
  - Overall: is total mRNA inflated? by how much?
  - Per locus: does mRNA over-estimation track the locus gDNA load (loci.tsv `gdna`)?
    -> if mRNA_est - mRNA_true ~ proportional to gDNA, the leak biases abundance (build the fix)
    -> if mRNA_est ~= mRNA_true regardless of gDNA, the leak is a hard-label artifact (stop)
  - Compare across gDNA level x strandedness x capture to see if any bias traces to gDNA.

Usage: python scripts/debug/em_strand_gating.py [--bam-name unused]
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

DEFAULT_SUITE = "/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_capture/hyb_capture_500kb"
CONDITIONS = [
    "gdna_none_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.99_nrna_none_capture_off",
    "gdna_high_ss_0.99_nrna_none_capture_on",
    "gdna_high_ss_0.50_nrna_none_capture_off",
    "gdna_high_ss_0.50_nrna_none_capture_on",
]


def load_quant(path: Path):
    """transcript_id -> (count_em, locus_id, is_nrna)."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["transcript_id"]] = (
                float(r["count_em"]),
                int(r["locus_id"]) if r["locus_id"] not in ("", "-1") else -1,
                r.get("is_nrna", "False") in ("True", "true", "1"),
            )
    return out


def load_truth(path: Path):
    """transcript_id -> true observed mRNA fragments."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["transcript_id"]] = float(r["observed_mrna_fragments"])
    return out


def load_loci(path: Path):
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[int(r["locus_id"])] = r
    return out


def fnum(d, k):
    try:
        return float(d.get(k, 0.0) or 0.0)
    except ValueError:
        return 0.0


def analyze(cond_dir: Path):
    quant = load_quant(cond_dir / "rigel_out" / "quant.tsv")
    truth = load_truth(cond_dir / "truth_abundances.tsv")
    loci = load_loci(cond_dir / "rigel_out" / "loci.tsv")

    # per-transcript mRNA est vs true
    tot_est = tot_true = 0.0
    per_locus_est = defaultdict(float)
    per_locus_true = defaultdict(float)
    tx_rows = []
    for tid, (est, lid, is_nrna) in quant.items():
        if is_nrna:
            continue
        tru = truth.get(tid, 0.0)
        tot_est += est
        tot_true += tru
        per_locus_est[lid] += est
        per_locus_true[lid] += tru
        tx_rows.append((tid, lid, est, tru))

    return quant, truth, loci, tot_est, tot_true, per_locus_est, per_locus_true, tx_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", default=DEFAULT_SUITE)
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    suite = Path(args.suite)

    # Discover conditions from the suite (any gdna_* dir with a quant.tsv) so this works across
    # suite layouts/axes (e.g. the gDNA strand-overdispersion sweep), not a hardcoded list.
    conditions = sorted(
        p.name for p in suite.iterdir()
        if p.is_dir() and (p / "rigel_out" / "quant.tsv").exists()
    )
    for cond in conditions:
        cdir = suite / cond
        quant, truth, loci, tot_est, tot_true, ple, plt_, tx_rows = analyze(cdir)
        infl = tot_est - tot_true
        infl_pct = 100.0 * infl / tot_true if tot_true else 0.0
        print(f"\n{'='*92}\n=== {cond} ===\n{'='*92}")
        print(f"  TOTAL mRNA: est(count_em)={tot_est:9.0f}  true={tot_true:9.0f}  "
              f"inflation={infl:+8.0f} ({infl_pct:+.2f}%)")

        # Per-locus: mRNA over-estimation vs the locus gDNA load (loci.tsv 'gdna').
        # If inflation ~ proportional to gDNA, the leak biases abundance.
        print(f"  per-locus mRNA est vs true, with EM gDNA load (top {args.top} by |Δmrna|):")
        print(f"    {'locus':>6} {'mrna_est':>9} {'mrna_true':>9} {'Δmrna':>8} "
              f"{'gdna_em':>8} {'Δmrna/gdna':>11}")
        rows = []
        for lid in set(list(ple) + list(plt_)):
            est, tru = ple.get(lid, 0.0), plt_.get(lid, 0.0)
            rows.append((abs(est - tru), lid, est, tru))
        rows.sort(reverse=True)
        for _, lid, est, tru in rows[: args.top]:
            L = loci.get(lid, {})
            gd = fnum(L, "gdna")
            d = est - tru
            ratio = d / gd if gd > 0 else float("nan")
            print(f"    {lid:>6} {est:>9.0f} {tru:>9.0f} {d:>+8.0f} {gd:>8.0f} {ratio:>11.2f}")


if __name__ == "__main__":
    main()
