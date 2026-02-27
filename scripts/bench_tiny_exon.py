#!/usr/bin/env python
"""
Tiny-exon benchmark: trapezoid vs uniform bias model.

Scenario
--------
Gene g1 on "+" strand, two transcripts sharing a large second exon:

  T1: exons (1000, 1020) + (10000, 20000)   — tiny 20 bp unique exon
  T2: exons (2000, 2300) + (10000, 20000)   — 300 bp unique exon

The shared exon (10000–20000) is 10 000 bp.  Fragments landing entirely
in the shared exon are ambiguous.  Only fragments touching the unique
first exons disambiguate.

T1 has a *tiny* unique exon (20 bp) so very few fragments can land
uniquely.  A uniform bias model under-credits T1 because it assumes
reads are equally likely anywhere; the trapezoid model accounts for the
positional constraint and should recover T1 more accurately.

We sweep T1:T2 abundance ratios and compare three model configurations:

  1. trapezoid + no OVR  (bias_model="trapezoid", gamma=0)
  2. uniform   + no OVR  (bias_model="uniform",   gamma=0)
  3. uniform   + OVR     (bias_model="uniform",   gamma=1.0)

Alpha is swept over {1e-3, 1e-2, 0.1, 1.0}.
"""

import logging
import sys
from itertools import product

from hulkrna.pipeline import run_pipeline
from hulkrna.sim import Scenario, SimConfig, run_benchmark

logging.basicConfig(
    level=logging.WARNING,
    format="%(levelname)s %(message)s",
)

# ── Parameters ──────────────────────────────────────────────────────────
N_FRAGMENTS = 10_000
SEED = 42
GENOME_LENGTH = 25_000

# T1:T2 ratios (T2 abundance fixed at 128, T1 varies)
T2_ABUNDANCE = 128.0
RATIOS = [1 / 128, 1 / 64, 1 / 8, 1 / 2, 1, 2, 8, 64, 128]

ALPHAS = [1e-3, 1e-2, 0.1, 1.0]

CONFIGS = [
    ("trapezoid_noOVR", "trapezoid", 0.0),
    ("uniform_noOVR",   "uniform",   0.0),
    ("uniform_OVR",     "uniform",   1.0),
]

# ── Scenario builder ────────────────────────────────────────────────────

def make_scenario(ratio: float):
    """Build a Scenario with T1:T2 = ratio."""
    t1_abundance = T2_ABUNDANCE * ratio
    sc = Scenario(
        f"tiny_exon_r{ratio}",
        genome_length=GENOME_LENGTH,
        seed=SEED,
    )
    sc.add_gene("g1", "+", [
        {
            "t_id": "t1",
            "exons": [(1000, 1020), (10000, 20000)],
            "abundance": t1_abundance,
        },
        {
            "t_id": "t2",
            "exons": [(2000, 2300), (10000, 20000)],
            "abundance": T2_ABUNDANCE,
        },
    ])
    return sc


# ── Main ────────────────────────────────────────────────────────────────

def main():
    sim_cfg = SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=1.0, seed=SEED,
    )

    # Header
    print(
        "ratio\talpha\tconfig\t"
        "t1_expected\tt1_observed\tt1_rel_err\t"
        "t2_expected\tt2_observed\tt2_rel_err\t"
        "total_abs_err"
    )

    for ratio in RATIOS:
        # Build the scenario + oracle BAM once per ratio
        sc = make_scenario(ratio)
        with sc:
            result = sc.build_oracle(
                n_fragments=N_FRAGMENTS,
                sim_config=sim_cfg,
            )

            for alpha, (cfg_name, bias_model, gamma) in product(ALPHAS, CONFIGS):
                pr = run_pipeline(
                    result.bam_path,
                    result.index,
                    sj_strand_tag="ts",
                    seed=SEED,
                    em_prior_alpha=alpha,
                    em_prior_gamma=gamma,
                    bias_model=bias_model,
                )
                bench = run_benchmark(
                    result, pr,
                    scenario_name=f"r{ratio}_a{alpha}_{cfg_name}",
                )

                t1 = next(t for t in bench.transcripts if t.t_id == "t1")
                t2 = next(t for t in bench.transcripts if t.t_id == "t2")

                print(
                    f"{ratio:.4f}\t{alpha}\t{cfg_name}\t"
                    f"{t1.expected}\t{t1.observed:.1f}\t{t1.rel_error:.4f}\t"
                    f"{t2.expected}\t{t2.observed:.1f}\t{t2.rel_error:.4f}\t"
                    f"{bench.total_abs_error:.1f}"
                )
                sys.stdout.flush()

    print("\n# Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
