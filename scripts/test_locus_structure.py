"""Prove whether fragment order changes locus structure (connected components).

Runs the same scenario twice: once with fragments in original order (1 thread),
once with multi-threaded scanning (non-deterministic order). Compares the
locus structure to see if fragments get assigned to different loci.
"""
import numpy as np
from pathlib import Path
import tempfile

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

SEED = 42
N_FRAGS = 5000

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
    sc = Scenario("locus_test", genome_length=10000, seed=SEED,
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


def run_and_inspect(bam_path, index, config, label):
    """Run pipeline, intercept locus structure, and return detailed info."""
    # We need to hook into the pipeline to observe pre-EM data.
    # Instead, let's use the pipeline internals directly.
    import rigel.pipeline as pl
    import rigel.locus as loc
    from rigel.scan import FragmentRouter
    from rigel.scoring import FragmentScorer

    pr = run_pipeline(bam_path, index, config=config)
    est = pr.estimator

    # Get loci info
    loci_df = est.get_loci_df()
    counts_df = est.get_counts_df(index)

    print(f"\n--- {label} ---")
    print(f"Number of loci: {len(loci_df)}")
    for _, row in loci_df.iterrows():
        lid = row['locus_id']
        print(f"  Locus {lid}: mrna={row['mrna']:.4f}, "
              f"nrna={row['nrna']:.4f}, gdna={row['gdna']:.4f}, "
              f"gdna_rate={row['gdna_rate']:.4f}")

    # Transcript-level
    tdf = counts_df.sort_values("transcript_id").reset_index(drop=True)
    for _, row in tdf.iterrows():
        if row['mrna'] > 0.01 or row['nrna'] > 0.01:
            print(f"  {row['transcript_id']}: mrna={row['mrna']:.4f}, "
                  f"nrna={row['nrna']:.4f}, locus_id={row['locus_id']}")

    return pr


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sc = build_scenario(tmp / "sim")
        gdna_config = GDNAConfig(abundance=10.0)
        result = sc.build_oracle(
            n_fragments=N_FRAGS, sim_config=SIM_SS90,
            gdna_config=gdna_config,
        )

        # Run with 1 thread (deterministic baseline)
        pr1 = run_and_inspect(
            result.bam_path, result.index,
            _config(n_scan_threads=1), "1-thread")

        # Run with N threads multiple times
        for i in range(3):
            prN = run_and_inspect(
                result.bam_path, result.index,
                _config(n_scan_threads=0), f"N-thread run {i}")

        sc.cleanup()


if __name__ == "__main__":
    main()
