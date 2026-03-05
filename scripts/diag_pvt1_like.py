#!/usr/bin/env python
"""
Diagnostic: PVT1-like synthetic scenarios to expose nRNA assignment bugs.

Scenarios:
  A: No overlap, mRNA only        → baseline (should be perfect)
  B: Exon overlap, mRNA only       → strand disambiguation needed
  C: nRNA from T1, intronic overlap → nRNA reads in T2 exon region
  D: Multi-isoform T1 + overlap    → isoform ambiguity + gene ambiguity
  E: B + imperfect strand (SS=0.9) → tests robustness of strand model
  F: C + imperfect strand (SS=0.9) → nRNA + weak strand model

In all scenarios:
  T1 = POS-strand multi-exon (like PVT1)
  T2 = NEG-strand single-exon (like LINC02912), abundance=0
"""

import logging
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import Scenario, SimConfig, run_benchmark

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(message)s")
logger = logging.getLogger(__name__)

N_FRAG = 2000
SEED = 42


def run_scenario(name, sc, *, n_fragments=N_FRAG, ss=1.0,
                 nrna_abundance=0.0, verbose=True):
    """Build oracle, run pipeline, return (bench, estimator, t_ids)."""
    sim_cfg = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=ss, seed=SEED,
    )
    result = sc.build_oracle(
        n_fragments=n_fragments,
        sim_config=sim_cfg,
        nrna_abundance=nrna_abundance,
    )

    config = PipelineConfig(
        em=EMConfig(seed=SEED),
        scan=BamScanConfig(sj_strand_tag="ts"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=name)
    est = pr.estimator
    t_ids = [t.t_id for t in result.transcripts]
    t_counts_per_t = est.t_counts.sum(axis=1)

    if verbose:
        print(f"\n{'='*70}")
        print(f"Scenario: {name}")
        print(f"{'='*70}")
        print(bench.summary())
        print(f"\n  Per-transcript detail:")
        for i, tid in enumerate(t_ids):
            mrna = t_counts_per_t[i]
            nrna_em = est.nrna_em_counts[i]
            nrna_init = est.nrna_init[i]
            print(f"    {tid} (idx={i}): mRNA={mrna:.1f}, "
                  f"nRNA_em={nrna_em:.2f}, nRNA_init={nrna_init:.2f}")

        total_nrna = float(est.nrna_em_counts.sum())
        print(f"\n  Pool totals:")
        print(f"    mRNA: {bench.total_observed:.1f}")
        print(f"    nRNA: {bench.n_nrna_pipeline:.1f} "
              f"(expected: {bench.n_nrna_expected})")
        print(f"    gDNA: {bench.n_gdna_pipeline:.1f}")

    sc.cleanup()
    return bench, est, t_ids


def scenario_a():
    """A: No overlap, mRNA only."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_a_"))
    sc = Scenario("scA", genome_length=8000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0},
    ])
    return run_scenario("A: no_overlap_mRNA_only", sc)


def scenario_b():
    """B: Exon overlap — T1 exon overlaps T2 exon spatially."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_b_"))
    sc = Scenario("scB", genome_length=8000, seed=SEED, work_dir=work_dir)
    # T1 exon3 at 4000-5600 overlaps T2 exon at 4200-4800
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (4000, 5600)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(4200, 4800)], "abundance": 0},
    ])
    return run_scenario("B: exon_overlap_mRNA_only", sc)


def scenario_c():
    """C: nRNA from T1.  Intronic reads can land in T2 exon region."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_c_"))
    sc = Scenario("scC", genome_length=8000, seed=SEED, work_dir=work_dir)
    # T1 intron from 1200-2000 and 2200-5000.
    # T2 exon at 3000-3800 sits within T1's second intron.
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0},
    ])
    # nRNA_abundance=50 creates unspliced reads from T1's full span
    return run_scenario("C: nRNA_from_T1_intronic_overlap", sc,
                        nrna_abundance=50)


def scenario_d():
    """D: Multi-isoform T1 + exon overlap + nRNA."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_d_"))
    sc = Scenario("scD", genome_length=8000, seed=SEED, work_dir=work_dir)
    # Two isoforms of g1, both POS, overlapping T2
    sc.add_gene("g1", "+", [
        {"t_id": "t1a",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
         "abundance": 60},
        {"t_id": "t1b",
         "exons": [(1000, 1200), (3500, 3900), (5000, 5600)],
         "abundance": 40},
    ])
    # T2 single-exon antisense at 3500-3900 — overlaps T1b exon2
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(3500, 3900)], "abundance": 0},
    ])
    return run_scenario("D: multi_isoform_overlap_nRNA", sc,
                        nrna_abundance=30)


def scenario_e():
    """E: Exon overlap + imperfect strand specificity (SS=0.9)."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_e_"))
    sc = Scenario("scE", genome_length=8000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (4000, 5600)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(4200, 4800)], "abundance": 0},
    ])
    return run_scenario("E: exon_overlap_SS=0.9", sc, ss=0.9)


def scenario_f():
    """F: nRNA from T1 + imperfect strand (SS=0.9)."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_f_"))
    sc = Scenario("scF", genome_length=8000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5600)],
         "abundance": 100},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2", "exons": [(3000, 3800)], "abundance": 0},
    ])
    return run_scenario("F: nRNA_from_T1_SS=0.9", sc,
                        nrna_abundance=50, ss=0.9)


def scenario_g():
    """G: T2 is MULTI-EXON (nRNA prior NOT zeroed) + nRNA from T1.

    This tests whether the strand model alone can prevent wrong-strand
    nRNA assignment when T2's nRNA prior is active.
    """
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_g_"))
    sc = Scenario("scG", genome_length=10000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (7000, 7600)],
         "abundance": 100},
    ])
    # T2: multi-exon NEG, within T1's intron → nRNA prior NOT zeroed
    sc.add_gene("g2", "-", [
        {"t_id": "t2",
         "exons": [(3000, 3400), (4500, 4900)],
         "abundance": 0},
    ])
    return run_scenario("G: multi_exon_T2_nRNA_from_T1", sc,
                        nrna_abundance=50)


def scenario_h():
    """H: Many isoforms of T1 (5) + multi-exon T2 + nRNA.

    Simulates PVT1's many isoforms creating massive ambiguity.
    """
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_h_"))
    sc = Scenario("scH", genome_length=12000, seed=SEED, work_dir=work_dir)
    # 5 isoforms of g1, various exon combos, all POS
    sc.add_gene("g1", "+", [
        {"t_id": "t1a",
         "exons": [(1000, 1200), (2000, 2200), (8000, 8600)],
         "abundance": 30},
        {"t_id": "t1b",
         "exons": [(1000, 1200), (3000, 3400), (8000, 8600)],
         "abundance": 25},
        {"t_id": "t1c",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5400), (8000, 8600)],
         "abundance": 20},
        {"t_id": "t1d",
         "exons": [(1000, 1200), (3000, 3400), (5000, 5400), (8000, 8600)],
         "abundance": 15},
        {"t_id": "t1e",
         "exons": [(1000, 1200), (8000, 8600)],
         "abundance": 10},
    ])
    # T2: multi-exon NEG within intron, abundance=0
    sc.add_gene("g2", "-", [
        {"t_id": "t2",
         "exons": [(4000, 4400), (6000, 6400)],
         "abundance": 0},
    ])
    return run_scenario("H: 5_isoforms_multi_exon_T2_nRNA", sc,
                        nrna_abundance=40)


def scenario_i():
    """I: H but with SS=0.9 (stress test)."""
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_i_"))
    sc = Scenario("scI", genome_length=12000, seed=SEED, work_dir=work_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1a",
         "exons": [(1000, 1200), (2000, 2200), (8000, 8600)],
         "abundance": 30},
        {"t_id": "t1b",
         "exons": [(1000, 1200), (3000, 3400), (8000, 8600)],
         "abundance": 25},
        {"t_id": "t1c",
         "exons": [(1000, 1200), (2000, 2200), (5000, 5400), (8000, 8600)],
         "abundance": 20},
        {"t_id": "t1d",
         "exons": [(1000, 1200), (3000, 3400), (5000, 5400), (8000, 8600)],
         "abundance": 15},
        {"t_id": "t1e",
         "exons": [(1000, 1200), (8000, 8600)],
         "abundance": 10},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t2",
         "exons": [(4000, 4400), (6000, 6400)],
         "abundance": 0},
    ])
    return run_scenario("I: 5_isoforms_multi_exon_T2_nRNA_SS=0.9", sc,
                        nrna_abundance=40, ss=0.9)


def scenario_j():
    """J: Exon overlap stressor - T2 exon inside T1 exon + nRNA + SS=0.9.

    This is the worst case: T2 exon physically overlaps T1 exon AND
    T1 generates intronic nRNA reads. T2 is multi-exon so nRNA isn't
    zeroed.
    """
    work_dir = Path(tempfile.mkdtemp(prefix="hulkrna_diag_j_"))
    sc = Scenario("scJ", genome_length=12000, seed=SEED, work_dir=work_dir)
    # T1 has wide exon3 at 4000-6600
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(1000, 1200), (2000, 2200), (4000, 6600)],
         "abundance": 100},
    ])
    # T2 multi-exon NEG, exon1 OVERLAPS T1's exon3
    sc.add_gene("g2", "-", [
        {"t_id": "t2",
         "exons": [(4200, 4800), (8000, 8400)],
         "abundance": 0},
    ])
    return run_scenario("J: exon_overlap_T2_multiexon_nRNA_SS=0.9", sc,
                        nrna_abundance=50, ss=0.9)


def main():
    print("=" * 70)
    print("PVT1-like diagnostic: synthetic scenario suite")
    print("=" * 70)

    results = {}
    for label, fn in [
        ("A", scenario_a),
        ("B", scenario_b),
        ("C", scenario_c),
        ("D", scenario_d),
        ("E", scenario_e),
        ("F", scenario_f),
        ("G", scenario_g),
        ("H", scenario_h),
        ("I", scenario_i),
        ("J", scenario_j),
    ]:
        try:
            bench, est, t_ids = fn()
            t2_idx = next(i for i, t in enumerate(t_ids) if t == "t2")
            t2_mrna = est.t_counts.sum(axis=1)[t2_idx]
            t2_nrna = est.nrna_em_counts[t2_idx]
            total_nrna = float(est.nrna_em_counts.sum())
            results[label] = {
                "t2_mrna": float(t2_mrna),
                "t2_nrna": float(t2_nrna),
                "total_nrna": total_nrna,
                "n_nrna_expected": bench.n_nrna_expected,
                "n_nrna_pipeline": bench.n_nrna_pipeline,
                "total_observed": bench.total_observed,
            }
        except Exception as e:
            print(f"\n*** Scenario {label} FAILED: {e}")
            import traceback
            traceback.print_exc()
            results[label] = {"error": str(e)}

    # Summary table
    print("\n\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Scenario':<40} {'T2 mRNA':>8} {'T2 nRNA':>8} "
          f"{'Tot nRNA':>9} {'nRNA exp':>9} {'nRNA pip':>9}")
    print("-" * 78)
    for label, r in results.items():
        if "error" in r:
            print(f"  {label:<38} {'ERROR':>8}")
        else:
            flag = " ***BUG***" if r["t2_nrna"] > 0.5 or (
                r["n_nrna_expected"] == 0 and r["total_nrna"] > 5
            ) or (
                r["n_nrna_expected"] > 0
                and r["total_nrna"] > r["n_nrna_expected"] * 1.5 + 10
            ) else ""
            print(f"  {label:<38} {r['t2_mrna']:>8.1f} {r['t2_nrna']:>8.1f} "
                  f"{r['total_nrna']:>9.1f} {r['n_nrna_expected']:>9} "
                  f"{r['n_nrna_pipeline']:>9.1f}{flag}")


if __name__ == "__main__":
    main()
