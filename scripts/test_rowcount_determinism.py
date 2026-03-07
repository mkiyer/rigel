#!/usr/bin/env python3
"""
Test: full pipeline determinism focusing on quant_detail row counts.

Runs the pipeline N_RUNS times with multi-threaded scanning and EM,
and compares quant_detail row counts. This targets the exact failure
mode observed: different row counts (280690 vs 280693) between runs.
"""
import sys
import tempfile
from pathlib import Path

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42
N_FRAGS = 10_000  # enough to trigger edge cases
N_RUNS = 10


def build_large_scenario(work_dir):
    """Scenario with many overlapping transcripts → many EM loci."""
    sc = Scenario("rowcount", genome_length=100_000, seed=SEED,
                  work_dir=work_dir)

    # Gene cluster 1: 3 overlapping transcripts
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(1000, 1500), (3000, 3500), (5000, 5500)],
         "abundance": 200},
        {"t_id": "t2", "exons": [(1000, 1500), (5000, 5500)],
         "abundance": 80},
        {"t_id": "t3", "exons": [(1000, 1500), (3000, 3500)],
         "abundance": 40},
    ])

    # Gene cluster 2: negative strand
    sc.add_gene("g2", "-", [
        {"t_id": "t4", "exons": [(15000, 15500), (17000, 17500)],
         "abundance": 150},
        {"t_id": "t5", "exons": [(15000, 15500), (17000, 17500), (19000, 19500)],
         "abundance": 50},
    ])

    # Gene cluster 3: complex overlap
    sc.add_gene("g3", "+", [
        {"t_id": "t6", "exons": [(30000, 30500), (32000, 32500), (34000, 34500)],
         "abundance": 100},
        {"t_id": "t7", "exons": [(30000, 30500), (34000, 34500)],
         "abundance": 30},
        {"t_id": "t8", "exons": [(30000, 30500), (32000, 32500), (34000, 34500), (36000, 36500)],
         "abundance": 10},
    ])

    # Gene 4: single transcript (control)
    sc.add_gene("g4", "-", [
        {"t_id": "t9", "exons": [(50000, 50500), (52000, 52500)],
         "abundance": 60},
    ])

    # Gene 5: two transcripts
    sc.add_gene("g5", "+", [
        {"t_id": "t10", "exons": [(70000, 70500), (72000, 72500)],
         "abundance": 90},
        {"t_id": "t11", "exons": [(70000, 70500)],
         "abundance": 20},
    ])

    result = sc.build_oracle(
        n_fragments=N_FRAGS,
        sim_config=SimConfig(seed=SEED, strand_specificity=0.9),
        gdna_config=GDNAConfig(abundance=10.0),
    )
    return sc, result


def main():
    with tempfile.TemporaryDirectory(prefix="rigel_rowcount_") as tmp:
        work_dir = Path(tmp)
        print("Building scenario...")
        sc, result = build_large_scenario(work_dir)
        print(f"  Transcripts: {result.index.num_transcripts}")
        print(f"  BAM: {result.bam_path}")

        # Use multi-threaded config (auto threads)
        config = PipelineConfig(
            em=EMConfig(seed=SEED, n_threads=0),
            scan=BamScanConfig(sj_strand_tag="auto", n_scan_threads=0),
        )

        row_counts = []
        detail_shapes = []
        for i in range(N_RUNS):
            pr = run_pipeline(result.bam_path, result.index, config=config)
            est = pr.estimator
            detail = est.get_detail_df(result.index)
            rc = len(detail)
            row_counts.append(rc)
            detail_shapes.append(detail.shape)
            print(f"  Run {i+1}/{N_RUNS}: quant_detail rows = {rc}, "
                  f"n_frags={pr.stats.n_fragments}, "
                  f"gdna_em={est.gdna_em_count:.6f}")

        # Check if all row counts are identical
        unique_counts = set(row_counts)
        if len(unique_counts) == 1:
            print(f"\n  ALL {N_RUNS} RUNS: IDENTICAL ROW COUNT = {row_counts[0]}")
        else:
            print(f"\n  *** ROW COUNT DIVERGENCE DETECTED ***")
            print(f"  Unique row counts: {sorted(unique_counts)}")
            print(f"  All counts: {row_counts}")
            sc.cleanup()
            return 1

        # Also test: compare multi-threaded with single-threaded
        print("\n--- Comparing 1-thread vs N-thread ---")
        config_1t = PipelineConfig(
            em=EMConfig(seed=SEED, n_threads=1),
            scan=BamScanConfig(sj_strand_tag="auto", n_scan_threads=1),
        )
        pr_1t = run_pipeline(result.bam_path, result.index, config=config_1t)
        detail_1t = pr_1t.estimator.get_detail_df(result.index)
        rc_1t = len(detail_1t)
        print(f"  1-thread: rows={rc_1t}")
        print(f"  N-thread: rows={row_counts[0]}")
        if rc_1t != row_counts[0]:
            print(f"  *** MISMATCH: 1t={rc_1t} vs Nt={row_counts[0]} ***")
            sc.cleanup()
            return 1
        else:
            print(f"  MATCH")

        sc.cleanup()
        return 0


if __name__ == "__main__":
    sys.exit(main())
