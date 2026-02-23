#!/usr/bin/env python
"""
nRNA sweep: 11×11 grid of mature RNA × nascent RNA abundance.

Gene G with two isoforms:
  T1: exons (500,1000), (2000,2500), (3000,3500)  — 3-exon, includes middle exon
  T2: exons (500,1000), (3000,3500)                — 2-exon, skips middle exon

Sweep:
  mRNA abundance:  0, 100, 200, …, 1000  (11 values)
  nRNA abundance:  0, 100, 200, …, 1000  (11 values)
  → 121 runs total

Fixed: gDNA = 0, strand_specificity = 1.0, n_fragments = 10_000

Usage:
    conda run -n hulkrna python scripts/nrna_sweep.py
"""

import json
import logging
import sys
from pathlib import Path

import numpy as np

# Ensure project root is on path
project_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(project_root / "src"))

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logging.basicConfig(
    level=logging.WARNING,
    format="%(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger("nrna_sweep")
logger.setLevel(logging.INFO)

# ── Parameters ──────────────────────────────────────────────────────
MRNA_RANGE = list(range(0, 1001, 100))   # 0, 100, 200, …, 1000
NRNA_RANGE = list(range(0, 1001, 100))   # 0, 100, 200, …, 1000
N_FRAGMENTS = 10_000
STRAND_SPECIFICITY = 1.0
SEED = 42
GENOME_LENGTH = 8000


def run_single(mrna_abund: int, nrna_abund: int, work_dir: Path) -> dict:
    """Run one (mRNA, nRNA) point and return a result dict."""
    tag = f"mrna{mrna_abund}_nrna{nrna_abund}"

    # At least one pool must have abundance > 0 to simulate
    # When mRNA=0 and nRNA=0: nothing to simulate → skip
    if mrna_abund == 0 and nrna_abund == 0:
        return {
            "mrna_abundance": mrna_abund,
            "nrna_abundance": nrna_abund,
            "t1_expected": 0, "t1_observed": 0.0,
            "t2_expected": 0, "t2_observed": 0.0,
            "nrna_expected": 0, "nrna_pipeline": 0.0,
            "gdna_expected": 0, "gdna_pipeline": 0.0,
            "n_fragments": 0, "n_intergenic": 0, "n_chimeric": 0,
            "total_expected": 0, "total_observed": 0.0,
            "skipped": True,
        }

    sc = Scenario(
        tag,
        genome_length=GENOME_LENGTH,
        seed=SEED,
        work_dir=work_dir / tag,
    )

    # Per-transcript nrna_abundance via add_gene dict
    t1_spec = {
        "t_id": "t1",
        "exons": [(500, 1000), (2000, 2500), (3000, 3500)],
        "abundance": mrna_abund,
        "nrna_abundance": nrna_abund,
    }
    t2_spec = {
        "t_id": "t2",
        "exons": [(500, 1000), (3000, 3500)],
        "abundance": mrna_abund,
        "nrna_abundance": nrna_abund,
    }
    sc.add_gene("g1", "+", [t1_spec, t2_spec])

    sim_config = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100,
        strand_specificity=STRAND_SPECIFICITY,
        seed=SEED,
    )

    try:
        result = sc.build(
            n_fragments=N_FRAGMENTS,
            sim_config=sim_config,
            gdna_config=None,
            nrna_abundance=0,  # use per-transcript values from add_gene
        )
        pr = run_pipeline(
            result.bam_path,
            result.index,
            sj_strand_tag="ts",
            seed=SEED,
        )
        bench = run_benchmark(result, pr, scenario_name=tag)

        # Extract per-transcript results
        t1_acc = next((t for t in bench.transcripts if t.t_id == "t1"), None)
        t2_acc = next((t for t in bench.transcripts if t.t_id == "t2"), None)

        return {
            "mrna_abundance": mrna_abund,
            "nrna_abundance": nrna_abund,
            "t1_expected": t1_acc.expected if t1_acc else 0,
            "t1_observed": round(t1_acc.observed, 1) if t1_acc else 0.0,
            "t2_expected": t2_acc.expected if t2_acc else 0,
            "t2_observed": round(t2_acc.observed, 1) if t2_acc else 0.0,
            "nrna_expected": bench.n_nrna_expected,
            "nrna_pipeline": round(bench.n_nrna_pipeline, 1),
            "gdna_expected": bench.n_gdna_expected,
            "gdna_pipeline": round(bench.n_gdna_pipeline, 1),
            "n_fragments": bench.n_fragments,
            "n_intergenic": bench.n_intergenic,
            "n_chimeric": bench.n_chimeric,
            "total_expected": bench.total_expected,
            "total_observed": round(bench.total_observed, 1),
            "total_rna_observed": round(bench.total_rna_observed, 1),
            "skipped": False,
        }
    finally:
        sc.cleanup()


def main():
    work_dir = Path("/tmp/nrna_sweep")
    work_dir.mkdir(parents=True, exist_ok=True)

    total = len(MRNA_RANGE) * len(NRNA_RANGE)
    logger.info("Starting nRNA sweep: %d × %d = %d runs", len(MRNA_RANGE), len(NRNA_RANGE), total)

    results = []
    for i, mrna in enumerate(MRNA_RANGE):
        for j, nrna in enumerate(NRNA_RANGE):
            idx = i * len(NRNA_RANGE) + j + 1
            logger.info("[%3d/%d] mRNA=%4d  nRNA=%4d", idx, total, mrna, nrna)
            row = run_single(mrna, nrna, work_dir)
            results.append(row)

    # ── Save raw JSON ──
    out_json = work_dir / "nrna_sweep_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    logger.info("Raw results → %s", out_json)

    # ── Print summary table ──
    print("\n" + "=" * 100)
    print("nRNA SWEEP RESULTS — mRNA × nRNA grid (11 × 11)")
    print("Gene G: T1 (3-exon), T2 (2-exon, skips middle)")
    print(f"Fixed: gDNA=0, SS={STRAND_SPECIFICITY}, n_frags={N_FRAGMENTS}")
    print("=" * 100)

    # Header
    hdr = (
        f"{'mRNA':>6} {'nRNA':>6} │ "
        f"{'T1_exp':>7} {'T1_obs':>7} {'T1_err':>7} │ "
        f"{'T2_exp':>7} {'T2_obs':>7} {'T2_err':>7} │ "
        f"{'nRNA_exp':>8} {'nRNA_pip':>8} {'nRNA_err':>8} │ "
        f"{'gDNA_pip':>8}"
    )
    print(hdr)
    print("─" * len(hdr))

    for r in results:
        if r.get("skipped"):
            print(f"{r['mrna_abundance']:>6} {r['nrna_abundance']:>6} │  (skipped — both zero)")
            continue

        t1_err = r["t1_observed"] - r["t1_expected"]
        t2_err = r["t2_observed"] - r["t2_expected"]
        nrna_err = r["nrna_pipeline"] - r["nrna_expected"]

        print(
            f"{r['mrna_abundance']:>6} {r['nrna_abundance']:>6} │ "
            f"{r['t1_expected']:>7} {r['t1_observed']:>7.1f} {t1_err:>+7.1f} │ "
            f"{r['t2_expected']:>7} {r['t2_observed']:>7.1f} {t2_err:>+7.1f} │ "
            f"{r['nrna_expected']:>8} {r['nrna_pipeline']:>8.1f} {nrna_err:>+8.1f} │ "
            f"{r['gdna_pipeline']:>8.1f}"
        )

    # ── Print mRNA-error heatmap (T1 + T2 combined) ──
    print("\n\n" + "=" * 80)
    print("HEATMAP: Total mRNA error (T1+T2 observed − expected)")
    print("Rows = mRNA abundance, Cols = nRNA abundance")
    print("=" * 80)

    # Build lookup
    lookup = {(r["mrna_abundance"], r["nrna_abundance"]): r for r in results}

    # Column header
    col_hdr = f"{'':>6}"
    for nrna in NRNA_RANGE:
        col_hdr += f" {nrna:>7}"
    print(col_hdr)
    print("─" * len(col_hdr))

    for mrna in MRNA_RANGE:
        row = f"{mrna:>6}"
        for nrna in NRNA_RANGE:
            r = lookup.get((mrna, nrna))
            if r is None or r.get("skipped"):
                row += f" {'---':>7}"
            else:
                err = (r["t1_observed"] + r["t2_observed"]) - (r["t1_expected"] + r["t2_expected"])
                row += f" {err:>+7.1f}"
        print(row)

    # ── Print nRNA-error heatmap ──
    print("\n\n" + "=" * 80)
    print("HEATMAP: nRNA pipeline error (pipeline − expected)")
    print("Rows = mRNA abundance, Cols = nRNA abundance")
    print("=" * 80)

    print(col_hdr)
    print("─" * len(col_hdr))

    for mrna in MRNA_RANGE:
        row = f"{mrna:>6}"
        for nrna in NRNA_RANGE:
            r = lookup.get((mrna, nrna))
            if r is None or r.get("skipped"):
                row += f" {'---':>7}"
            else:
                err = r["nrna_pipeline"] - r["nrna_expected"]
                row += f" {err:>+7.1f}"
        print(row)

    # ── Print nRNA phantom heatmap (nRNA pipeline when nRNA truth = 0) ──
    print("\n\n" + "=" * 80)
    print("HEATMAP: nRNA pipeline count (absolute)")
    print("Rows = mRNA abundance, Cols = nRNA abundance")
    print("=" * 80)

    print(col_hdr)
    print("─" * len(col_hdr))

    for mrna in MRNA_RANGE:
        row = f"{mrna:>6}"
        for nrna in NRNA_RANGE:
            r = lookup.get((mrna, nrna))
            if r is None or r.get("skipped"):
                row += f" {'---':>7}"
            else:
                row += f" {r['nrna_pipeline']:>7.1f}"
        print(row)

    # ── Print gDNA pipeline heatmap ──
    print("\n\n" + "=" * 80)
    print("HEATMAP: gDNA pipeline count (should be ~0 everywhere)")
    print("Rows = mRNA abundance, Cols = nRNA abundance")
    print("=" * 80)

    print(col_hdr)
    print("─" * len(col_hdr))

    for mrna in MRNA_RANGE:
        row = f"{mrna:>6}"
        for nrna in NRNA_RANGE:
            r = lookup.get((mrna, nrna))
            if r is None or r.get("skipped"):
                row += f" {'---':>7}"
            else:
                row += f" {r['gdna_pipeline']:>7.1f}"
        print(row)

    print(f"\nResults saved to {out_json}")


if __name__ == "__main__":
    main()
