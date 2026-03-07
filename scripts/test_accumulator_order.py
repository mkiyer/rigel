"""Prove that nrna_init differs between runs due to FP accumulation order.

Hooks into the pipeline to extract intermediate values and compare them.
"""
import numpy as np
from pathlib import Path
import tempfile

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42
N_FRAGS = 5000
N_RUNS = 3

SIM_SS90 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.9, seed=SEED,
)


def _config(n_scan_threads):
    return PipelineConfig(
        em=EMConfig(seed=SEED, n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto", n_scan_threads=n_scan_threads),
    )


def build_scenario(tmp_dir):
    sc = Scenario("accum_test", genome_length=10000, seed=SEED,
                  work_dir=tmp_dir)
    sc.add_gene("g1", "+", [
        {"t_id": "t1",
         "exons": [(500, 800), (1200, 1400), (2000, 2300)],
         "abundance": 100},
        {"t_id": "t2",
         "exons": [(500, 800), (1600, 1800), (2000, 2300)],
         "abundance": 30},
    ])
    sc.add_gene("g2", "-", [
        {"t_id": "t3",
         "exons": [(4000, 4300), (5000, 5400)],
         "abundance": 50},
    ])
    sc.add_gene("g3", "+", [
        {"t_id": "t4",
         "exons": [(7000, 7500)],
         "abundance": 20},
    ])
    return sc


def run_and_extract(bam_path, index, config):
    """Run pipeline and extract intermediate arrays."""
    pr = run_pipeline(bam_path, index, config=config)
    est = pr.estimator

    return {
        'intronic_sense': est.transcript_intronic_sense.copy(),
        'intronic_antisense': est.transcript_intronic_antisense.copy(),
        'exonic_sense': est.transcript_exonic_sense.copy(),
        'exonic_antisense': est.transcript_exonic_antisense.copy(),
        'unspliced_sense': est.transcript_unspliced_sense.copy(),
        'unspliced_antisense': est.transcript_unspliced_antisense.copy(),
        'nrna_init': est.nrna_init.copy() if hasattr(est, 'nrna_init') and est.nrna_init is not None else None,
        'unambig_counts': est.unambig_counts.copy(),
        'gdna_em_count': est.gdna_em_count,
        'nrna_em_count': est.nrna_em_count,
        'mrna': est.get_counts_df(index).sort_values('transcript_id').reset_index(drop=True)['mrna'].to_numpy(),
        'nrna': est.get_counts_df(index).sort_values('transcript_id').reset_index(drop=True)['nrna'].to_numpy(),
    }


def compare_arrays(a, b, name):
    """Compare two arrays and report differences."""
    if a is None and b is None:
        return
    diff = np.abs(a - b)
    max_diff = np.max(diff)
    if max_diff == 0:
        print(f"  {name}: BIT-EXACT")
    else:
        max_val = max(np.max(np.abs(a)), np.max(np.abs(b)))
        n_diff = np.sum(diff > 0)
        rel = max_diff / max_val if max_val > 0 else float('inf')
        print(f"  {name}: max_abs_diff={max_diff:.6e}, rel={rel:.6e}, "
              f"n_diff={n_diff}/{len(a)}")
        # Show the actual values where they differ most
        idx = np.argmax(diff)
        print(f"    worst at index {idx}: a={a.flat[idx]:.18e}, b={b.flat[idx]:.18e}")


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = build_scenario(tmp / "sim")
        gdna_config = GDNAConfig(abundance=10.0)
        result = sc.build_oracle(
            n_fragments=N_FRAGS, sim_config=SIM_SS90,
            gdna_config=gdna_config,
        )

        # Single-thread reference
        ref = run_and_extract(result.bam_path, result.index, _config(1))

        print("=== 1-thread self-consistency ===")
        ref2 = run_and_extract(result.bam_path, result.index, _config(1))
        for key in ['intronic_sense', 'intronic_antisense',
                     'exonic_sense', 'exonic_antisense',
                     'unspliced_sense', 'unspliced_antisense',
                     'nrna_init', 'unambig_counts', 'mrna', 'nrna']:
            compare_arrays(ref[key], ref2[key], key)

        print(f"\n=== N-thread vs 1-thread (fragment order differs) ===")
        for run in range(N_RUNS):
            print(f"\n--- N-thread run {run} ---")
            cur = run_and_extract(result.bam_path, result.index, _config(0))
            for key in ['intronic_sense', 'intronic_antisense',
                         'exonic_sense', 'exonic_antisense',
                         'unspliced_sense', 'unspliced_antisense',
                         'nrna_init', 'unambig_counts', 'mrna', 'nrna']:
                compare_arrays(ref[key], cur[key], key)
            print(f"  gdna_em_count: ref={ref['gdna_em_count']:.18e}, "
                  f"cur={cur['gdna_em_count']:.18e}, "
                  f"diff={abs(ref['gdna_em_count'] - cur['gdna_em_count']):.6e}")

        sc.cleanup()


if __name__ == "__main__":
    main()
