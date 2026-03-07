#!/usr/bin/env python3
"""
Reproduction script: compare single-threaded vs multi-threaded BAM scan.

Creates a synthetic scenario with overlapping transcripts, runs the BAM
scanner in single-threaded mode and multi-threaded mode, and compares
the accumulator sizes (fragment counts).  Any difference proves a
structural threading bug.
"""

import sys
import tempfile
from pathlib import Path

import numpy as np

from rigel._bam_impl import BamScanner
from rigel.sim import Scenario, SimConfig, GDNAConfig


def build_scenario(work_dir: Path) -> tuple:
    """Build a scenario with multiple overlapping genes → more resolution paths."""
    sc = Scenario(
        "thread_test",
        genome_length=50_000,
        seed=42,
        work_dir=work_dir,
    )

    # Gene A: two overlapping transcripts (spliced + unspliced paths)
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(500, 1000), (2000, 2500), (3500, 4000)],
         "abundance": 200},
        {"t_id": "t2", "exons": [(500, 1000), (3500, 4000)],
         "abundance": 100},
    ])

    # Gene B: negative strand, single transcript
    sc.add_gene("g2", "-", [
        {"t_id": "t3", "exons": [(10000, 10500), (11000, 11500)],
         "abundance": 150},
    ])

    # Gene C: three transcripts with complex overlap
    sc.add_gene("g3", "+", [
        {"t_id": "t4", "exons": [(20000, 20500), (21000, 21500), (22000, 22500)],
         "abundance": 80},
        {"t_id": "t5", "exons": [(20000, 20500), (22000, 22500)],
         "abundance": 60},
        {"t_id": "t6", "exons": [(20000, 20500), (21000, 21500)],
         "abundance": 40},
    ])

    # Gene D: another negative strand gene
    sc.add_gene("g4", "-", [
        {"t_id": "t7", "exons": [(35000, 35500), (36000, 36500), (37000, 37500)],
         "abundance": 120},
    ])

    result = sc.build_oracle(
        n_fragments=5000,
        sim_config=SimConfig(seed=42),
        gdna_config=GDNAConfig(abundance=10.0),
    )
    return result


def run_scan(index, bam_path: str, threaded: bool, n_workers: int = 4):
    """Run a BAM scan and return (accumulator_size, stats_dict, finalized_raw)."""
    resolve_ctx = index._resolver
    resolve_ctx.set_gene_strands(index.g_to_strand_arr.tolist())
    resolve_ctx.set_transcript_strands(index.t_to_strand_arr.tolist())

    scanner = BamScanner(
        resolve_ctx, "none",
        skip_duplicates=True,
        include_multimap=False,
    )

    if threaded:
        result = scanner.scan(bam_path, n_workers=n_workers)
    else:
        result = scanner.scan(bam_path)

    acc_size = result["accumulator_size"]
    stats = result["stats"]

    # Finalize to get detailed fragment data
    raw = None
    if acc_size > 0:
        raw = scanner.finalize_accumulator(index.t_to_strand_arr.tolist())

    return acc_size, stats, raw


def compare_results(label_a, acc_a, stats_a, raw_a,
                    label_b, acc_b, stats_b, raw_b):
    """Compare two scan results and report differences."""
    print(f"\n{'='*60}")
    print(f"Comparing: {label_a} vs {label_b}")
    print(f"{'='*60}")

    # Compare accumulator sizes
    print(f"\n  Accumulator size: {label_a}={acc_a}, {label_b}={acc_b}", end="")
    if acc_a != acc_b:
        print(f"  *** MISMATCH (diff={acc_a - acc_b}) ***")
    else:
        print("  OK")

    # Compare key stats
    stat_keys = [
        "n_read_names", "unique", "multimapping",
        "n_fragments", "n_chimeric",
        "n_intergenic_unspliced", "n_intergenic_spliced",
        "n_with_exon", "n_with_annotated_sj", "n_with_unannotated_sj",
        "n_same_strand", "n_ambig_strand",
    ]
    print(f"\n  Key stats comparison:")
    any_stat_diff = False
    for key in stat_keys:
        va = stats_a.get(key, "N/A")
        vb = stats_b.get(key, "N/A")
        marker = "" if va == vb else "  *** MISMATCH ***"
        if va != vb:
            any_stat_diff = True
        print(f"    {key:40s}: {va:>10} vs {vb:>10}{marker}")

    # Compare finalized fragment data
    if raw_a is not None and raw_b is not None:
        size_a = raw_a["size"]
        size_b = raw_b["size"]
        print(f"\n  Finalized size: {label_a}={size_a}, {label_b}={size_b}", end="")
        if size_a != size_b:
            print(f"  *** MISMATCH (diff={size_a - size_b}) ***")
        else:
            print("  OK")

        # Compare frag_ids
        fids_a = set(np.frombuffer(raw_a["frag_id"], dtype=np.int64))
        fids_b = set(np.frombuffer(raw_b["frag_id"], dtype=np.int64))
        only_a = fids_a - fids_b
        only_b = fids_b - fids_a
        if only_a or only_b:
            print(f"\n  Frag IDs only in {label_a}: {len(only_a)}")
            if only_a and len(only_a) <= 20:
                print(f"    {sorted(only_a)}")
            print(f"  Frag IDs only in {label_b}: {len(only_b)}")
            if only_b and len(only_b) <= 20:
                print(f"    {sorted(only_b)}")
        else:
            print(f"\n  Frag ID sets: identical ({len(fids_a)} unique IDs)")

        # Compare t_indices (CSR)
        ti_a = np.frombuffer(raw_a["t_indices"], dtype=np.int32)
        ti_b = np.frombuffer(raw_b["t_indices"], dtype=np.int32)
        print(f"  Total transcript candidates: {label_a}={len(ti_a)}, "
              f"{label_b}={len(ti_b)}", end="")
        if len(ti_a) != len(ti_b):
            print(f"  *** MISMATCH ***")
        else:
            print("  OK")

    return acc_a == acc_b and not any_stat_diff


def main():
    with tempfile.TemporaryDirectory(prefix="rigel_thread_test_") as tmpdir:
        work_dir = Path(tmpdir)
        print("Building scenario...")
        result = build_scenario(work_dir)
        bam_path = str(result.bam_path)
        index = result.index
        print(f"  BAM: {bam_path}")
        print(f"  Transcripts: {index.num_transcripts}")

        # Run single-threaded
        print("\nRunning single-threaded scan...")
        acc_1t, stats_1t, raw_1t = run_scan(index, bam_path, threaded=False)
        print(f"  accumulator_size = {acc_1t}")

        # Run multi-threaded with various worker counts
        all_match = True
        for n_workers in [2, 4, 8]:
            print(f"\nRunning multi-threaded scan (n_workers={n_workers})...")
            acc_mt, stats_mt, raw_mt = run_scan(
                index, bam_path, threaded=True, n_workers=n_workers)
            print(f"  accumulator_size = {acc_mt}")

            match = compare_results(
                f"1-thread", acc_1t, stats_1t, raw_1t,
                f"{n_workers}-threads", acc_mt, stats_mt, raw_mt,
            )
            if not match:
                all_match = False

        # Also compare two multi-threaded runs against each other
        print(f"\nRunning two multi-threaded scans (n_workers=4) for self-consistency...")
        acc_mt_a, stats_mt_a, raw_mt_a = run_scan(
            index, bam_path, threaded=True, n_workers=4)
        acc_mt_b, stats_mt_b, raw_mt_b = run_scan(
            index, bam_path, threaded=True, n_workers=4)
        match = compare_results(
            "4-threads (run A)", acc_mt_a, stats_mt_a, raw_mt_a,
            "4-threads (run B)", acc_mt_b, stats_mt_b, raw_mt_b,
        )
        if not match:
            all_match = False

        print(f"\n{'='*60}")
        if all_match:
            print("ALL COMPARISONS PASSED — no fragment count divergence detected")
            return 0
        else:
            print("FAILURES DETECTED — fragment counts differ between modes")
            return 1


if __name__ == "__main__":
    sys.exit(main())
