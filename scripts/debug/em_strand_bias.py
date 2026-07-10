"""Test the STRAND-SCORING bias against gDNA (the user's hypothesis). Pure EM at the density level.

3 components sharing one 2kb single-exon region: (+RNA), (−RNA), gDNA. All same length ⇒ eff-length is
irrelevant. Fragment strand-likelihoods (the SCORING model, put in log_liks; the EM adds no strand term):
  gDNA        : 0.5 for any fragment           (the log-half penalty)
  +RNA        : p_sense if frag +, 1−p_sense if frag −
  −RNA        : p_sense if frag −, 1−p_sense if frag +
Feed the PRODUCTION EM a PERFECT prior (gdna_prior=true gDNA, rna_prior=true RNA total). Sweep densities.
Contrast AMBIG (both RNA strands present) vs SINGLE-STRAND (only +RNA). Hypothesis: with BOTH strands,
EVERY gDNA fragment (either strand) is beaten 2:1 by the matching-strand RNA ⇒ gDNA under-assigned; and it
worsens with strand-specificity. Single-strand only beats the matching half of the (50/50) gDNA.
"""
import sys, math
sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests")
import numpy as np
from rigel.estimator import AbundanceEstimator
from rigel.config import EMConfig
from rigel.locus_partition import partition_and_free
from conftest import _make_locus_em_data


def run(d_pos, d_neg, d_gdna, ambig, p, iters=300):
    LP, LA, LG = math.log(p), math.log(1.0 - p), math.log(0.5)
    # counts of +strand / −strand fragments generated at these densities (deterministic split)
    n_plus = int(round(d_pos * p + d_neg * (1 - p) + d_gdna * 0.5))
    n_minus = int(round(d_pos * (1 - p) + d_neg * p + d_gdna * 0.5))
    N = 2 if ambig else 1
    t_per_unit, lik_per_unit = [], []
    for _ in range(n_plus):                       # +strand fragment
        t_per_unit.append([0, 1] if ambig else [0])
        lik_per_unit.append([LP, LA] if ambig else [LP])   # +RNA sense, −RNA antisense
    for _ in range(n_minus):                      # −strand fragment
        t_per_unit.append([0, 1] if ambig else [0])
        lik_per_unit.append([LA, LP] if ambig else [LA])   # +RNA antisense, −RNA sense
    rc = AbundanceEstimator(N, em_config=EMConfig(seed=42))
    em_data, loci, gprior, index = _make_locus_em_data(
        t_per_unit, log_liks_per_unit=lik_per_unit, num_transcripts=N,
        include_gdna=True, gdna_prior_count=float(d_gdna), gdna_log_lik=LG)
    partitions = partition_and_free(em_data, loci)
    parts = [partitions[i] for i in range(len(loci))]
    ptuples = [(pp.offsets, pp.t_indices, pp.log_liks, pp.coverage_weights, pp.count_cols,
                pp.is_spliced, pp.gdna_log_liks, pp.locus_t_indices, pp.locus_count_cols) for pp in parts]
    total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
        ptuples, [loc.transcript_indices for loc in loci], gprior, index,
        rna_prior_count=np.array([float(d_pos + d_neg)]), em_iterations=iters)
    return float(total_gdna), float(rc.em_counts.sum()), n_plus + n_minus


for p in [0.99, 0.90, 0.50]:
    print(f"\n===== strand specificity p_sense = {p}  (perfect prior; equal eff-len) =====")
    print(f"{'true: pos/neg/gdna':>22} | {'AMBIG (both strands)':>28} | {'SINGLE-STRAND (+RNA only)':>28}")
    print(f"{'':22} | {'gdna_asgn  gdna_err  rna':>28} | {'gdna_asgn  gdna_err  rna':>28}")
    R = 1000
    for G in [500, 1000, 2000, 3000]:
        gA, rA, _ = run(R // 2, R // 2, G, True, p)      # AMBIG: pos=neg=R/2
        gS, rS, _ = run(R, 0, G, False, p)               # single: pos=R, neg=0
        print(f"{R//2:>6}/{R//2:<5}/{G:<7} | {gA:>9,.0f} {gA-G:>+9,.0f} {rA:>7,.0f} "
              f"| {gS:>9,.0f} {gS-G:>+9,.0f} {rS:>7,.0f}")
