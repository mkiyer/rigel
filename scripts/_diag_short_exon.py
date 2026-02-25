#!/usr/bin/env python3
"""Diagnostic: Short-exon shared-region isoform estimation bias.

Scenario
--------
Two isoforms of the same gene on chr1, positive strand:

    T1+  exons (1000, 1050), (5000, 15000)   — 50 bp unique + 10 000 bp shared
    T2+  exons (2000, 2300), (5000, 15000)   — 300 bp unique + 10 000 bp shared

T1 has a tiny 50 bp first exon; T2 has a 300 bp first exon.
They share the large [5000, 15000) exon.  Fragments landing in the
shared region are ambiguous and must be resolved by the EM.

This scenario exposes the per-fragment effective length and coverage-
weight corrections: without them, the EM distributes shared fragments
uniformly, over-counting the shorter transcript (T1).

Sweep
-----
- T1 abundance: 0, 1, 10, 100, 1000, 10000, 100000
- T2 abundance: 0, 1, 10, 100, 1000, 10000, 100000
- n_fragments:  10, 100, 1000, 10000

For each (T1_abund, T2_abund, n_frags) combination:
  - build oracle BAM (perfect alignments, no gDNA/nRNA)
  - run pipeline
  - record ground truth vs pipeline mRNA counts per transcript

Output
------
Writes a tab-separated report to stdout and summary stats to stderr.
"""

import itertools
import sys
import tempfile
from pathlib import Path

import numpy as np

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import Scenario, SimConfig, run_benchmark


# =====================================================================
# Configuration
# =====================================================================

ABUNDANCES = [0, 1, 10, 100, 1000, 10_000, 100_000]
N_FRAGMENTS_LIST = [10, 100, 1000, 10_000]
GENOME_LENGTH = 20_000
SIM_SEED = 42
PIPELINE_SEED = 42


def sim_config():
    return SimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=1.0,
        seed=SIM_SEED,
    )


def run_one(t1_abund, t2_abund, n_frags, work_dir):
    """Run one scenario combo and return a dict of results."""
    # Skip if both are zero (nothing to simulate)
    if t1_abund == 0 and t2_abund == 0:
        return None

    name = f"se_t1={t1_abund}_t2={t2_abund}_n={n_frags}"
    sc = Scenario(
        name,
        genome_length=GENOME_LENGTH,
        seed=SIM_SEED,
        work_dir=work_dir / name,
    )
    sc.add_gene("g1", "+", [
        {
            "t_id": "t1",
            "exons": [(1000, 1050), (5000, 15000)],
            "abundance": t1_abund,
        },
        {
            "t_id": "t2",
            "exons": [(2000, 2300), (5000, 15000)],
            "abundance": t2_abund,
        },
    ])
    # Helper gene for strand model training (needs spliced reads)
    sc.add_gene("g_helper", "+", [
        {
            "t_id": "t_helper",
            "exons": [(16000, 16500), (17000, 17500)],
            "abundance": 50,
        },
    ])

    try:
        result = sc.build_oracle(
            n_fragments=n_frags,
            sim_config=sim_config(),
            gdna_config=None,
            nrna_abundance=0,
        )
        pr = run_pipeline(
            result.bam_path,
            result.index,
            sj_strand_tag="ts",
            seed=PIPELINE_SEED,
        )
        bench = run_benchmark(result, pr, scenario_name=name)

        # Extract per-transcript data
        t1_data = next((t for t in bench.transcripts if t.t_id == "t1"), None)
        t2_data = next((t for t in bench.transcripts if t.t_id == "t2"), None)

        row = {
            "t1_abund": t1_abund,
            "t2_abund": t2_abund,
            "n_frags": n_frags,
            "n_simulated": bench.n_simulated,
            "n_pipeline": bench.n_fragments,
            "t1_expected": t1_data.expected if t1_data else 0,
            "t1_observed": t1_data.observed if t1_data else 0.0,
            "t1_abs_diff": t1_data.abs_diff if t1_data else 0.0,
            "t1_rel_err": t1_data.rel_error if t1_data else 0.0,
            "t2_expected": t2_data.expected if t2_data else 0,
            "t2_observed": t2_data.observed if t2_data else 0.0,
            "t2_abs_diff": t2_data.abs_diff if t2_data else 0.0,
            "t2_rel_err": t2_data.rel_error if t2_data else 0.0,
            "total_expected": bench.total_expected,
            "total_observed": bench.total_observed,
            "total_abs_error": bench.total_abs_error,
        }
        return row
    finally:
        sc.cleanup()


def main():
    header = [
        "t1_abund", "t2_abund", "n_frags",
        "n_simulated", "n_pipeline",
        "t1_expected", "t1_observed", "t1_abs_diff", "t1_rel_err",
        "t2_expected", "t2_observed", "t2_abs_diff", "t2_rel_err",
        "total_expected", "total_observed", "total_abs_error",
    ]
    print("\t".join(header))

    results = []
    combos = list(itertools.product(ABUNDANCES, ABUNDANCES, N_FRAGMENTS_LIST))
    total = len(combos)
    done = 0

    with tempfile.TemporaryDirectory(prefix="diag_short_exon_") as tmp:
        work_dir = Path(tmp)
        for t1_abund, t2_abund, n_frags in combos:
            done += 1
            if t1_abund == 0 and t2_abund == 0:
                continue
            print(
                f"[{done}/{total}] t1={t1_abund} t2={t2_abund} n={n_frags}",
                file=sys.stderr,
                flush=True,
            )
            try:
                row = run_one(t1_abund, t2_abund, n_frags, work_dir)
            except Exception as e:
                print(
                    f"  ERROR: {e}",
                    file=sys.stderr,
                    flush=True,
                )
                continue
            if row is None:
                continue
            results.append(row)
            vals = [row[h] for h in header]
            line = "\t".join(
                f"{v:.4f}" if isinstance(v, float) else str(v)
                for v in vals
            )
            print(line, flush=True)

    # ---- Summary statistics to stderr ----
    if not results:
        print("No results collected.", file=sys.stderr)
        return

    print("\n" + "=" * 72, file=sys.stderr)
    print("SHORT-EXON SCENARIO — ERROR SUMMARY", file=sys.stderr)
    print("=" * 72, file=sys.stderr)

    # Geometry reminder
    print("\nGeometry:", file=sys.stderr)
    print("  T1: exons (1000,1050)+(5000,15000) = 10050bp exonic, "
          "50bp unique", file=sys.stderr)
    print("  T2: exons (2000,2300)+(5000,15000) = 10300bp exonic, "
          "300bp unique", file=sys.stderr)
    print("  Shared region: [5000,15000) = 10000bp", file=sys.stderr)

    # Filter to cases where both transcripts have non-zero expected
    both_present = [
        r for r in results
        if r["t1_expected"] > 0 and r["t2_expected"] > 0
    ]
    # Filter to cases where only one is present (purity test)
    t1_only = [r for r in results if r["t1_expected"] > 0 and r["t2_expected"] == 0]
    t2_only = [r for r in results if r["t2_expected"] > 0 and r["t1_expected"] == 0]

    for label, subset in [
        ("Both T1+T2 present", both_present),
        ("T1 only (T2=0)", t1_only),
        ("T2 only (T1=0)", t2_only),
    ]:
        if not subset:
            continue
        t1_rel = [r["t1_rel_err"] for r in subset if r["t1_expected"] > 0]
        t2_rel = [r["t2_rel_err"] for r in subset if r["t2_expected"] > 0]
        t1_abs = [r["t1_abs_diff"] for r in subset]
        t2_abs = [r["t2_abs_diff"] for r in subset]

        print(f"\n--- {label} ({len(subset)} cases) ---", file=sys.stderr)
        if t1_rel:
            arr = np.array(t1_rel)
            print(
                f"  T1 rel error:  mean={arr.mean():.3f}  "
                f"median={np.median(arr):.3f}  "
                f"max={arr.max():.3f}  "
                f"p95={np.percentile(arr, 95):.3f}",
                file=sys.stderr,
            )
        if t2_rel:
            arr = np.array(t2_rel)
            print(
                f"  T2 rel error:  mean={arr.mean():.3f}  "
                f"median={np.median(arr):.3f}  "
                f"max={arr.max():.3f}  "
                f"p95={np.percentile(arr, 95):.3f}",
                file=sys.stderr,
            )
        if t1_abs:
            arr = np.array(t1_abs)
            print(
                f"  T1 abs error:  mean={arr.mean():.1f}  "
                f"median={np.median(arr):.1f}  "
                f"max={arr.max():.1f}",
                file=sys.stderr,
            )
        if t2_abs:
            arr = np.array(t2_abs)
            print(
                f"  T2 abs error:  mean={arr.mean():.1f}  "
                f"median={np.median(arr):.1f}  "
                f"max={arr.max():.1f}",
                file=sys.stderr,
            )

    # Bias direction analysis: does the pipeline over-count T1 (short exon)?
    print(f"\n--- Bias direction (both present, n≥100) ---", file=sys.stderr)
    large_both = [
        r for r in both_present
        if r["n_frags"] >= 100
        and r["t1_expected"] > 5 and r["t2_expected"] > 5
    ]
    if large_both:
        t1_bias = [(r["t1_observed"] - r["t1_expected"]) / r["t1_expected"]
                    for r in large_both]
        t2_bias = [(r["t2_observed"] - r["t2_expected"]) / r["t2_expected"]
                    for r in large_both]
        arr1 = np.array(t1_bias)
        arr2 = np.array(t2_bias)
        n_t1_over = int(np.sum(arr1 > 0.01))
        n_t1_under = int(np.sum(arr1 < -0.01))
        n_t2_over = int(np.sum(arr2 > 0.01))
        n_t2_under = int(np.sum(arr2 < -0.01))
        print(
            f"  T1 (50bp exon):  mean bias={arr1.mean():+.3f}  "
            f"over={n_t1_over}/{len(large_both)}  "
            f"under={n_t1_under}/{len(large_both)}",
            file=sys.stderr,
        )
        print(
            f"  T2 (300bp exon): mean bias={arr2.mean():+.3f}  "
            f"over={n_t2_over}/{len(large_both)}  "
            f"under={n_t2_under}/{len(large_both)}",
            file=sys.stderr,
        )
    else:
        print("  (no qualifying cases)", file=sys.stderr)

    # Per n_frags breakdown
    print(f"\n--- Error by read depth (both present) ---", file=sys.stderr)
    for nf in N_FRAGMENTS_LIST:
        subset_nf = [r for r in both_present if r["n_frags"] == nf]
        if not subset_nf:
            continue
        t1_rel = [r["t1_rel_err"] for r in subset_nf if r["t1_expected"] > 0]
        t2_rel = [r["t2_rel_err"] for r in subset_nf if r["t2_expected"] > 0]
        all_rel = t1_rel + t2_rel
        if all_rel:
            arr = np.array(all_rel)
            print(
                f"  n_frags={nf:>6d}: "
                f"mean_rel={arr.mean():.3f}  "
                f"median_rel={np.median(arr):.3f}  "
                f"max_rel={arr.max():.3f}  "
                f"({len(subset_nf)} combos)",
                file=sys.stderr,
            )

    # Worst cases
    print(f"\n--- Top 10 worst relative errors (both present) ---",
          file=sys.stderr)
    scored = []
    for r in both_present:
        worst = max(r["t1_rel_err"], r["t2_rel_err"])
        scored.append((worst, r))
    scored.sort(key=lambda x: -x[0])
    for i, (worst_rel, r) in enumerate(scored[:10]):
        print(
            f"  {i+1}. t1_ab={r['t1_abund']:>6d}  t2_ab={r['t2_abund']:>6d}  "
            f"n={r['n_frags']:>5d}  "
            f"T1: exp={r['t1_expected']:>5d} obs={r['t1_observed']:>8.1f} "
            f"({r['t1_rel_err']:+.1%})  "
            f"T2: exp={r['t2_expected']:>5d} obs={r['t2_observed']:>8.1f} "
            f"({r['t2_rel_err']:+.1%})",
            file=sys.stderr,
        )

    print("\n" + "=" * 72, file=sys.stderr)


if __name__ == "__main__":
    main()
