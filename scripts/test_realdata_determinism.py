#!/usr/bin/env python3
"""
Test quant_detail row count determinism on real data.

Runs the full pipeline N times on a real BAM and compares:
- accumulator sizes
- quant_detail row counts
- key scalar outputs
"""
import sys
import time
from pathlib import Path

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import run_pipeline

N_RUNS = 3  # Real data is slow, so just 3 runs
SEED = 42


def main():
    idx_path = sys.argv[1] if len(sys.argv) > 1 else "/Users/mkiyer/Downloads/rigel_runs/rigel_index"
    bam_path = sys.argv[2] if len(sys.argv) > 2 else "/Users/mkiyer/Downloads/rigel_runs/sim_ccle_hela_cervix/gdna_high_ss_0.90_nrna_low/sim_oracle.bam"

    print(f"Index: {idx_path}")
    print(f"BAM:   {bam_path}")

    print("Loading index...")
    t0 = time.time()
    index = TranscriptIndex.load(idx_path)
    print(f"  Loaded in {time.time()-t0:.1f}s, {index.num_transcripts} transcripts")

    # Multi-threaded config
    config = PipelineConfig(
        em=EMConfig(seed=SEED, n_threads=0),
        scan=BamScanConfig(sj_strand_tag="auto", n_scan_threads=0),
    )

    results = []
    for i in range(N_RUNS):
        print(f"\n--- Run {i+1}/{N_RUNS} ---")
        t0 = time.time()
        pr = run_pipeline(bam_path, index, config=config)
        dt = time.time() - t0
        est = pr.estimator

        detail = est.get_detail_df(index)
        rc = len(detail)

        info = {
            "n_frags": pr.stats.n_fragments,
            "detail_rows": rc,
            "gdna_em": est.gdna_em_count,
            "nrna_em": est.nrna_em_count,
            "time": dt,
        }
        results.append(info)
        print(f"  detail_rows={rc}, n_frags={info['n_frags']}, "
              f"gdna_em={info['gdna_em']:.6f}, "
              f"nrna_em={info['nrna_em']:.6f}, "
              f"time={dt:.1f}s")

    # Compare
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    row_counts = [r["detail_rows"] for r in results]
    unique = set(row_counts)
    print(f"  quant_detail row counts: {row_counts}")
    if len(unique) == 1:
        print(f"  --> DETERMINISTIC (all {N_RUNS} runs = {row_counts[0]})")
    else:
        print(f"  --> NON-DETERMINISTIC: {sorted(unique)}")

    frag_counts = [r["n_frags"] for r in results]
    print(f"  n_fragments: {frag_counts}")

    gdna_vals = [r["gdna_em"] for r in results]
    print(f"  gdna_em: {gdna_vals}")

    return 0 if len(unique) == 1 else 1


if __name__ == "__main__":
    sys.exit(main())
