#!/usr/bin/env python3
"""Test whether nRNA siphon error is explained by stochastic variance.

Experiment 1: Run the same scenario across 50 random seeds at n_rna=10,000
Experiment 2: Scale to n_rna=100,000 across 20 seeds
Experiment 3: Vary abundance ratio (nRNA/mRNA) at fixed n_rna=10,000

Reports both absolute and relative error.

Usage:
    conda activate rigel
    python scripts/debug/nrna_siphon_variance.py
"""

import gc
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "src"))

from rigel.config import BamScanConfig, EMConfig, FragmentScoringConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario, SimConfig

logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")


def build_and_run(mrna_ab, nrna_ab, n_rna, ss, seed, mode="map"):
    """Build scenario, run pipeline, return summary dict."""
    sc = Scenario(
        "siphon",
        genome_length=50000,
        seed=seed,
        work_dir=Path(tempfile.mkdtemp(prefix="rigel_sv_")),
    )
    sc.add_gene("GA", "+", [{
        "t_id": "TA1",
        "exons": [(1000, 1020), (5000, 5500), (12000, 13000)],
        "abundance": mrna_ab,
        "nrna_abundance": nrna_ab,
    }])
    sc.add_gene("GD", "+", [{
        "t_id": "TD",
        "exons": [(28000, 29000)],
        "abundance": 0,
        "nrna_abundance": 0,
    }])
    sc.add_gene("GE", "-", [{
        "t_id": "TE",
        "exons": [(40000, 41000), (45000, 46000)],
        "abundance": 0,
        "nrna_abundance": 0,
    }])
    sim_cfg = SimConfig(
        frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
        read_length=150, strand_specificity=ss, seed=seed,
    )
    result = sc.build_oracle(
        n_fragments=n_rna,
        sim_config=sim_cfg,
        nrna_abundance=0.0,
        n_rna_fragments=n_rna,
        gdna_fraction=0.0,
    )

    pipe_cfg = PipelineConfig(
        em=EMConfig(seed=seed, mode=mode),
        scan=BamScanConfig(sj_strand_tag="auto"),
        scoring=FragmentScoringConfig(),
    )

    gt_mrna = result.ground_truth_from_bam()
    gt_nrna = result.ground_truth_nrna_count_from_bam()

    pr = run_pipeline(result.bam_path, result.index, config=pipe_cfg)
    t_df = pr.estimator.get_counts_df(result.index)

    mrna_exp = sum(gt_mrna.values())
    mrna_obs = t_df["mrna"].sum()
    nrna_exp = gt_nrna
    nrna_obs = t_df["nrna"].sum()

    # Get unambig/EM breakdown
    mrna_unambig = t_df["mrna_unambig"].sum() if "mrna_unambig" in t_df.columns else 0
    mrna_em = t_df["mrna_em"].sum() if "mrna_em" in t_df.columns else 0

    sc.cleanup()
    gc.collect()

    return {
        "seed": seed,
        "n_rna": n_rna,
        "mrna_ab": mrna_ab,
        "nrna_ab": nrna_ab,
        "ratio": nrna_ab / max(mrna_ab, 1),
        "mode": mode,
        "mrna_exp": mrna_exp,
        "mrna_obs": mrna_obs,
        "mrna_unambig": mrna_unambig,
        "mrna_em": mrna_em,
        "nrna_exp": nrna_exp,
        "nrna_obs": nrna_obs,
        "mrna_err_abs": mrna_obs - mrna_exp,
        "mrna_err_rel": (mrna_obs - mrna_exp) / max(mrna_exp, 1),
    }


def print_summary(df, label):
    """Print summary stats for a set of runs."""
    print(f"\n{'=' * 80}")
    print(f"  {label}")
    print(f"{'=' * 80}")
    print(f"  {'seed':>5s} {'n_rna':>7s} {'ratio':>6s} {'mrna_exp':>9s} "
          f"{'mrna_obs':>9s} {'unambig':>8s} {'em':>8s} "
          f"{'err_abs':>8s} {'err_rel':>8s}")
    print(f"  {'-'*5:>5s} {'-'*7:>7s} {'-'*6:>6s} {'-'*9:>9s} "
          f"{'-'*9:>9s} {'-'*8:>8s} {'-'*8:>8s} "
          f"{'-'*8:>8s} {'-'*8:>8s}")
    for _, r in df.iterrows():
        print(f"  {int(r['seed']):5d} {int(r['n_rna']):7d} {r['ratio']:6.0f} "
              f"{r['mrna_exp']:9.1f} {r['mrna_obs']:9.1f} "
              f"{r['mrna_unambig']:8.1f} {r['mrna_em']:8.1f} "
              f"{r['mrna_err_abs']:+8.1f} {r['mrna_err_rel']:+8.3f}")

    print(f"\n  Aggregate statistics (N={len(df)}):")
    print(f"    mrna_exp:     mean={df['mrna_exp'].mean():.1f}, "
          f"std={df['mrna_exp'].std():.1f}, "
          f"range=[{df['mrna_exp'].min():.0f}, {df['mrna_exp'].max():.0f}]")
    print(f"    mrna_obs:     mean={df['mrna_obs'].mean():.1f}, "
          f"std={df['mrna_obs'].std():.1f}, "
          f"range=[{df['mrna_obs'].min():.1f}, {df['mrna_obs'].max():.1f}]")
    print(f"    mrna_unambig: mean={df['mrna_unambig'].mean():.1f}, "
          f"std={df['mrna_unambig'].std():.1f}")
    print(f"    err_abs:      mean={df['mrna_err_abs'].mean():+.1f}, "
          f"std={df['mrna_err_abs'].std():.1f}, "
          f"range=[{df['mrna_err_abs'].min():+.1f}, {df['mrna_err_abs'].max():+.1f}]")
    print(f"    err_rel:      mean={df['mrna_err_rel'].mean():+.3f}, "
          f"std={df['mrna_err_rel'].std():.3f}")

    # Theoretical mRNA fraction
    ab_m = df["mrna_ab"].iloc[0]
    ab_n = df["nrna_ab"].iloc[0]
    n_rna = df["n_rna"].iloc[0]
    theoretical_frac = ab_m / (ab_m + ab_n)
    theoretical_count = n_rna * theoretical_frac
    print(f"\n    Theoretical mRNA fraction: {ab_m}/({ab_m}+{ab_n}) "
          f"= {theoretical_frac:.4f}")
    print(f"    Theoretical mRNA count: {n_rna} × {theoretical_frac:.4f} "
          f"= {theoretical_count:.1f}")
    print(f"    Observed/Theoretical: "
          f"{df['mrna_obs'].mean():.1f}/{theoretical_count:.1f} "
          f"= {df['mrna_obs'].mean()/max(theoretical_count,1):.3f}")


def main():
    MRNA_AB = 20
    NRNA_AB = 500
    SS = 1.0

    # =================================================================
    # Experiment 1: Multi-seed sweep at n_rna=10,000
    # =================================================================
    print("\n" + "#" * 80)
    print("# EXPERIMENT 1: Multi-seed sweep (n_rna=10,000, 50 seeds)")
    print("#" * 80)

    rows = []
    for seed in range(50):
        r = build_and_run(MRNA_AB, NRNA_AB, n_rna=10000, ss=SS, seed=seed)
        rows.append(r)
        if (seed + 1) % 10 == 0:
            print(f"  ... completed {seed + 1}/50 seeds", flush=True)

    df1 = pd.DataFrame(rows)
    print_summary(df1, "Experiment 1: n_rna=10,000 × 50 seeds (MAP-EM)")

    # =================================================================
    # Experiment 2: Scale to n_rna=100,000 (20 seeds)
    # =================================================================
    print("\n" + "#" * 80)
    print("# EXPERIMENT 2: Scaled-up (n_rna=100,000, 20 seeds)")
    print("#" * 80)

    rows = []
    for seed in range(20):
        r = build_and_run(MRNA_AB, NRNA_AB, n_rna=100000, ss=SS, seed=seed)
        rows.append(r)
        if (seed + 1) % 5 == 0:
            print(f"  ... completed {seed + 1}/20 seeds", flush=True)

    df2 = pd.DataFrame(rows)
    print_summary(df2, "Experiment 2: n_rna=100,000 × 20 seeds (MAP-EM)")

    # =================================================================
    # Experiment 3: Vary abundance ratio (10 seeds each)
    # =================================================================
    print("\n" + "#" * 80)
    print("# EXPERIMENT 3: Abundance ratio sweep (n_rna=10,000, 10 seeds each)")
    print("#" * 80)

    ratios = [
        (20, 20),    # 1:1
        (20, 100),   # 5:1
        (20, 500),   # 25:1
        (20, 2000),  # 100:1
        (100, 500),  # 5:1, higher counts
    ]

    for mrna_ab, nrna_ab in ratios:
        rows = []
        for seed in range(10):
            r = build_and_run(mrna_ab, nrna_ab, n_rna=10000, ss=SS, seed=seed)
            rows.append(r)

        df_r = pd.DataFrame(rows)
        print_summary(df_r, f"Ratio sweep: mRNA={mrna_ab}, nRNA={nrna_ab} "
                      f"(ratio={nrna_ab/mrna_ab:.0f}:1)")

    print("\n" + "=" * 80)
    print("DONE")
    print("=" * 80)


if __name__ == "__main__":
    main()
