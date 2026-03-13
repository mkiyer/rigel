#!/usr/bin/env python3
"""
Trace EM dynamics for the simplest nRNA+gDNA scenario.

This script simulates the EXACT EM computation that Rigel performs for a single
locus containing TA1 (multi-exon +strand mRNA) with nRNA (NTA1) and gDNA.

Scenario:
  - TA1: 3-exon transcript, +strand, exonic_len=1520, span=12000
  - NTA1: nRNA spanning [1000, 13000], span=12000
  - gDNA: genome-wide, locus_span=12000
  - FL identical for RNA and gDNA (mean=250)
  - Various strand specificities

Fragment types generated (for a run with TA1=100, NTA1=100, SS=0.9, gf=0.3):
  - mRNA spliced sense: goes to TA1 only (not gDNA candidate)
  - nRNA unspliced sense: candidates = {TA1_mRNA, NTA1_nRNA, g_pos}
  - nRNA unspliced antisense: candidates = {TA1_mRNA, NTA1_nRNA, g_neg}
  - gDNA unspliced sense: candidates = {TA1_mRNA, NTA1_nRNA, g_pos}
  - gDNA unspliced antisense: candidates = {TA1_mRNA, NTA1_nRNA, g_neg}
  - gDNA spliced sense: candidates = {TA1_mRNA} only (gDNA gets -inf for spliced)

We'll trace through the EM iteration by iteration to understand dynamics.
"""

import numpy as np
import math


def run_em_trace(
    n_mrna_spliced_sense: int,
    n_unspliced_sense: int,      # nRNA + gDNA sense unspliced
    n_unspliced_antisense: int,  # nRNA + gDNA antisense unspliced
    SS: float,
    prior_alpha: float = 0.01,
    prior_gamma: float = 1.0,
    kappa: float = 6.0,
    max_iter: int = 200,
    use_ovr: bool = True,
    label: str = "",
):
    """
    Trace EM for a simplified locus.

    Components:
        0: TA1_mRNA
        1: NTA1_nRNA
        2: g_pos
        3: g_neg

    Fragment equivalence classes:
        EC_A: spliced sense → candidates {TA1_mRNA}  (n_mrna_spliced_sense fragments)
        EC_B: unspliced sense → candidates {TA1_mRNA, NTA1_nRNA, g_pos}
        EC_C: unspliced antisense → candidates {TA1_mRNA, NTA1_nRNA, g_neg}

    Log-likelihoods (after bias correction):
        For EC_B (unspliced sense):
            TA1_mRNA:  log(SS) + FL - log(eff_mRNA)    [strand=sense, FL on exonic len]
            NTA1_nRNA: log(SS) + FL_g - log(eff_nRNA)  [strand=sense, FL on genomic footprint]
            g_pos:     0      + FL_g - log(eff_gDNA)   [no strand term, FL on footprint]

        For EC_C (unspliced antisense):
            TA1_mRNA:  log(1-SS) + FL - log(eff_mRNA)
            NTA1_nRNA: log(1-SS) + FL_g - log(eff_nRNA)
            g_neg:     0         + FL_g - log(eff_gDNA)

        For EC_A (spliced sense):
            TA1_mRNA: log(SS) + FL - log(eff_mRNA)  [only candidate]

    Since RNA FL == gDNA FL (by design), and for unspliced fragments the footprint
    is the same regardless of whether they map to nRNA or gDNA, FL terms cancel in
    the posterior computation within each EC.

    The ONLY discriminating signal is:
        nRNA vs gDNA:  log(SS) - 0 = log(SS) for sense fragments
                       log(1-SS) - 0 = log(1-SS) for antisense fragments

    And the mRNA vs nRNA/gDNA discrimination comes from the different effective
    lengths (mRNA: 1520-250+1=1271, nRNA/gDNA: 12000-250+1=11751).
    """

    nc = 4  # TA1_mRNA, NTA1_nRNA, g_pos, g_neg
    log_SS = math.log(SS) if SS > 0 else -1e10
    log_1mSS = math.log(1 - SS) if SS < 1.0 else -1e10

    # Effective lengths after bias correction
    # mRNA exonic_len = 1520, nRNA span = 12000, gDNA locus_span = 12000
    # frag_len ~ 250
    eff_mrna = max(1520 - 250 + 1, 1)   # 1271
    eff_nrna = max(12000 - 250 + 1, 1)  # 11751
    eff_gdna = max(12000 - 250 + 1, 1)  # 11751
    log_eff = [math.log(eff_mrna), math.log(eff_nrna),
               math.log(eff_gdna), math.log(eff_gdna)]

    # Log-likelihoods per EC (after FL + bias correction)
    # We normalize relative to gDNA since FL terms cancel for nRNA vs gDNA

    # EC_A: spliced sense → {mRNA only}
    # (only mRNA is a candidate; gDNA gets -inf for spliced)
    ec_A_ll = np.array([log_SS - log_eff[0]])
    ec_A_cw = np.array([1.0])  # coverage weight
    ec_A_idx = np.array([0])   # component: TA1_mRNA
    ec_A_n = n_mrna_spliced_sense

    # EC_B: unspliced sense → {mRNA, nRNA, g_pos}
    # For unspliced fragments, mRNA uses SPLICED frag_len but nRNA/gDNA use
    # GENOMIC footprint. Actually for mRNA, the scoring uses frag_len (the
    # exonic length), not the genomic footprint. And the bias correction uses
    # eff_len = transcript_length - frag_len + 1.
    # But wait: for unspliced fragments hitting mRNA, the exonic basepairs
    # may be less than read_length (overhang). Let's simplify by assuming all
    # unspliced fragments have the same FL terms (they do, since FL distributions
    # are identical). The key difference is the bias correction.

    # Actually the FL MODEL is called with different arguments:
    #   mRNA: frag_len_log_lik(exonic_frag_len)
    #   nRNA: frag_len_log_lik(genomic_footprint)
    #   gDNA: gdna_frag_len_log_lik(genomic_footprint)
    # Since RNA FL == gDNA FL in our config, and for unspliced fragments the
    # genomic_footprint is the same for nRNA and gDNA, the FL terms for nRNA
    # and gDNA are identical.
    #
    # For mRNA, the exonic_frag_len may differ from genomic_footprint (if
    # the fragment spans an intron, the exonic bases are less). But we're
    # simplifying here to focus on the nRNA/gDNA competition.
    #
    # Let's use relative log-liks (subtract a constant so gDNA is 0):
    # The key insight: FL_rna(exonic) ≠ FL_rna(genomic) for fragments that
    # span introns. For intronic fragments, mRNA gets a very different FL
    # score. But for purely exonic unspliced fragments, they're similar.
    #
    # For simplicity, let's assume the intronic fragments are the majority
    # and focus on the nRNA vs gDNA competition:

    # Relative to gDNA baseline (0 + FL_g - log(eff_gdna)):
    # nRNA: log(SS) + FL_g - log(eff_nrna) - [0 + FL_g - log(eff_gdna)]
    #      = log(SS) + log(eff_gdna) - log(eff_nrna)
    #      = log(SS) + 0  [since eff_nrna == eff_gdna]
    #      = log(SS)
    delta_nrna_vs_gdna = log_SS  # sense: nRNA advantage over g_pos

    # For mRNA in unspliced fragments: depends heavily on fragment type.
    # Intronic fragments won't map to mRNA at all (no exonic bases).
    # Let's just track nRNA vs gDNA for now.

    ec_B_ll = np.array([log_SS, log_SS, 0.0])  # nRNA gets log(SS), gDNA gets 0
    ec_B_cw = np.array([1.0, 1.0, 1.0])
    ec_B_idx = np.array([0, 1, 2])  # mRNA, nRNA, g_pos
    # Actually let's simplify: most unspliced fragments in the intron won't
    # have mRNA as a candidate (they don't overlap exons). Let's separate:

    # For simplicity, model 2 ECs of unspliced:
    #  EC_B1: exonic unspliced sense → {mRNA, nRNA, g_pos}
    #  EC_B2: intronic unspliced sense → {nRNA, g_pos}  (mRNA not a candidate)
    # Similarly for antisense.

    # Actually the most important case is intronic fragments where it's
    # purely nRNA vs g_pos (or g_neg). Let's focus on that:

    # EC for intronic unspliced sense: candidates = {nRNA, g_pos}
    ec_intr_sense_ll = np.array([log_SS, 0.0])
    ec_intr_sense_cw = np.array([1.0, 1.0])
    ec_intr_sense_idx = np.array([1, 2])  # nRNA, g_pos

    # EC for intronic unspliced antisense: candidates = {nRNA, g_neg}
    ec_intr_anti_ll = np.array([log_1mSS, 0.0])
    ec_intr_anti_cw = np.array([1.0, 1.0])
    ec_intr_anti_idx = np.array([1, 3])  # nRNA, g_neg

    # EC for spliced sense: candidates = {mRNA}
    ec_spl_sense_ll = np.array([0.0])  # doesn't matter, only 1 candidate
    ec_spl_sense_cw = np.array([1.0])
    ec_spl_sense_idx = np.array([0])   # mRNA

    equiv_classes = [
        ("spliced_sense", ec_spl_sense_ll, ec_spl_sense_cw, ec_spl_sense_idx,
         n_mrna_spliced_sense),
        ("intronic_sense", ec_intr_sense_ll, ec_intr_sense_cw, ec_intr_sense_idx,
         n_unspliced_sense),
        ("intronic_anti", ec_intr_anti_ll, ec_intr_anti_cw, ec_intr_anti_idx,
         n_unspliced_antisense),
    ]

    # ---- OVR Warm Start ----
    # Mimics compute_ovr_prior_and_warm_start()
    unambig_totals = np.zeros(nc)
    # In Rigel, unambig = fragments with only 1 candidate.
    # Spliced sense fragments are unambig (only mRNA candidate)
    unambig_totals[0] = n_mrna_spliced_sense  # all go to mRNA

    # eligible: all components (gdna_init > 0 assumed)
    eligible = np.ones(nc)

    # theta_init starts from unambig_totals
    theta_init = unambig_totals.copy()

    # Accumulate coverage shares from ambiguous ECs
    coverage_totals = np.zeros(nc)
    n_ambiguous = 0

    if use_ovr:
        for name, ll, cw, idx, n_frags in equiv_classes:
            if len(idx) <= 1:
                continue  # skip unambig ECs
            for _ in range(n_frags):
                # Coverage-weight share (NOT likelihood-based!)
                weights = cw * eligible[idx]
                row_sum = weights.sum()
                if row_sum == 0:
                    row_sum = 1.0
                shares = weights / row_sum
                for j, c in enumerate(idx):
                    theta_init[c] += shares[j]
                    coverage_totals[c] += shares[j]
                n_ambiguous += 1

    # Prior: alpha_flat + gamma-scaled OVR
    prior = np.zeros(nc)
    if n_ambiguous > 0 and prior_gamma > 0:
        for i in range(nc):
            ovr_i = coverage_totals[i] * prior_gamma / n_ambiguous
            prior[i] = prior_alpha + ovr_i if eligible[i] > 0 else 0.0
    else:
        for i in range(nc):
            prior[i] = prior_alpha if eligible[i] > 0 else 0.0

    # Initial theta (normalized)
    state = theta_init + prior
    state /= state.sum()

    print(f"\n{'='*70}")
    print(f"  {label}")
    print(f"  SS={SS}, kappa={kappa}, gamma={'ON' if use_ovr else 'OFF'}")
    print(f"  Fragments: {n_mrna_spliced_sense} spliced sense, "
          f"{n_unspliced_sense} unspl sense, {n_unspliced_antisense} unspl anti")
    print(f"{'='*70}")
    print(f"\n  Unambig totals: mRNA={unambig_totals[0]:.0f}, nRNA={unambig_totals[1]:.0f}, "
          f"g_pos={unambig_totals[2]:.0f}, g_neg={unambig_totals[3]:.0f}")
    if use_ovr:
        print(f"  OVR coverage:   mRNA={coverage_totals[0]:.1f}, nRNA={coverage_totals[1]:.1f}, "
              f"g_pos={coverage_totals[2]:.1f}, g_neg={coverage_totals[3]:.1f}")
    print(f"  Prior:          mRNA={prior[0]:.4f}, nRNA={prior[1]:.4f}, "
          f"g_pos={prior[2]:.4f}, g_neg={prior[3]:.4f}")
    print(f"  theta_init:     mRNA={state[0]:.6f}, nRNA={state[1]:.6f}, "
          f"g_pos={state[2]:.6f}, g_neg={state[3]:.6f}")
    print()

    # ---- EM Iterations ----
    print(f"  {'Iter':>4}  {'mRNA':>10}  {'nRNA':>10}  {'g_pos':>10}  {'g_neg':>10}  "
          f"{'gDNA_tot':>10}  {'nRNA/gpos':>10}")

    for iteration in range(max_iter):
        # E-step: compute posteriors for each EC
        log_weights = np.log(state + 1e-300)

        em_totals = np.zeros(nc)

        for name, ll, cw, idx, n_frags in equiv_classes:
            if n_frags == 0:
                continue
            k = len(idx)
            # log posterior = ll + log_weights[idx]
            log_post = ll + log_weights[idx]
            # log-sum-exp normalize
            max_lp = log_post.max()
            exp_lp = np.exp(log_post - max_lp)
            post = exp_lp / exp_lp.sum()

            for j, c in enumerate(idx):
                em_totals[c] += post[j] * n_frags

        # M-step: theta_new = unambig + em + prior (unnormalized)
        theta_new = unambig_totals + em_totals + prior

        # Strand symmetry coupling
        if kappa > 2.0:
            g_pos_idx, g_neg_idx = 2, 3
            R_pos = theta_new[g_pos_idx]
            R_neg = theta_new[g_neg_idx]
            R_g = R_pos + R_neg
            if R_g > 0:
                phi = (R_pos + kappa / 2 - 1) / (R_g + kappa - 2)
                phi = max(0.0, min(1.0, phi))
                theta_new[g_pos_idx] = phi * R_g
                theta_new[g_neg_idx] = (1 - phi) * R_g

        # Normalize
        total = theta_new.sum()
        if total > 0:
            state_new = theta_new / total
        else:
            state_new = state.copy()

        # Check convergence
        delta = np.abs(state_new - state).sum()

        state = state_new.copy()

        if iteration < 10 or iteration % 10 == 0 or delta < 1e-8:
            gdna_tot = state[2] + state[3]
            ratio = state[1] / state[2] if state[2] > 1e-20 else float('inf')
            print(f"  {iteration:4d}  {state[0]:10.6f}  {state[1]:10.6f}  "
                  f"{state[2]:10.6f}  {state[3]:10.6f}  {gdna_tot:10.6f}  "
                  f"{ratio:10.4f}")

        if delta < 1e-10:
            print(f"  Converged at iteration {iteration}")
            break

    # Summary
    total_frags = n_mrna_spliced_sense + n_unspliced_sense + n_unspliced_antisense
    print(f"\n  Final theta: mRNA={state[0]:.6f}, nRNA={state[1]:.6f}, "
          f"g_pos={state[2]:.6f}, g_neg={state[3]:.6f}")
    print(f"  Final counts: mRNA={state[0]*total_frags:.1f}, "
          f"nRNA={state[1]*total_frags:.1f}, "
          f"g_pos={state[2]*total_frags:.1f}, "
          f"g_neg={state[3]*total_frags:.1f}")
    print(f"  gDNA total: {(state[2]+state[3])*total_frags:.1f}")

    return state


def main():
    # Scenario: TA1=100, NTA1=100, SS=0.9, gf=0.3
    # n_rna_fragments=10000, gdna_fraction=0.3
    #
    # With TA1=100 and NTA1=100, the sim splits 10000 RNA fragments.
    # The exact split depends on abundance × effective_length.
    # Let's compute: mrna_weight = 100 × max(0, 1520-250+1) = 100 × 1271 = 127100
    #                nrna_weight = 100 × max(0, 12000-250+1) = 100 × 11751 = 1175100
    # total = 127100 + 1175100 = 1302200
    # n_mrna = round(10000 × 127100/1302200) = round(976) = 976
    # n_nrna = 10000 - 976 = 9024
    # n_gdna = round(10000 × 0.3) = 3000

    # Of the 976 mRNA fragments, at SS=0.9:
    #   ~90% are sense, ~10% antisense
    #   Most are spliced (since TA1 is multi-exon with large introns)
    # Let's say ~90% spliced: ~790 spliced sense, ~88 spliced anti
    # ~10% unspliced: ~88 unspliced sense, ~10 unspliced anti
    # (These are rough; the exact numbers depend on fragment overlap with introns)

    # Of the 9024 nRNA fragments:
    #   SS=0.9: ~90% sense, ~10% antisense
    #   All are unspliced (nRNA spans the full pre-mRNA)
    #   sense: ~8122, antisense: ~902

    # Of the 3000 gDNA fragments:
    #   50% sense, 50% antisense (gDNA is unstranded)
    #   ~7.2% intergenic (outside locus), rest in locus
    #   Of in-locus: some spliced-looking (rare), most unspliced
    #   sense in locus: ~1392, antisense in locus: ~1392

    # Intronic unspliced sense total: 8122 (nRNA) + 1392 (gDNA) = 9514
    # Intronic unspliced antisense total: 902 (nRNA) + 1392 (gDNA) = 2294

    # The TRUE theta should be proportional to:
    #   mRNA:  976/13000 = 0.0751
    #   nRNA: 9024/13000 = 0.6942
    #   g_pos: 1500/13000 = 0.1154  (half of 3000)
    #   g_neg: 1500/13000 = 0.1154

    # But many gDNA fragments are intergenic (don't enter the locus EM).
    # Let's estimate: locus span = 12000bp out of 50000bp genome.
    # gDNA in locus: 3000 × 12000/50000 = 720
    # gDNA sense in locus: 360, antisense in locus: 360

    # TRUE theta (in-locus only):
    #   total_in_locus = 976 + 9024 + 720 = 10720
    #   mRNA:  976/10720 = 0.0910
    #   nRNA: 9024/10720 = 0.8418
    #   g_pos: 360/10720 = 0.0336
    #   g_neg: 360/10720 = 0.0336

    print("=" * 70)
    print("  EM DYNAMICS TRACER")
    print("  Scenario: TA1=100, NTA1=100, locus [1000,13000]")
    print("=" * 70)

    # Approximate fragment counts (SS=0.9, gf=0.3)
    # Spliced mRNA sense: ~790  (unambig — only maps to TA1)
    # Unspliced intronic sense: ~8122 (nRNA) + ~360 (gDNA) = ~8482
    # Unspliced intronic antisense: ~902 (nRNA) + ~360 (gDNA) = ~1262
    # (ignoring exonic unspliced for now — they add mRNA as candidate but
    #  mRNA eff_len is much shorter, so mRNA gets proportionally more)

    scenarios = [
        # (n_spl_sense, n_unspl_sense, n_unspl_anti, SS, kappa, gamma_on, label)
        (790, 8482, 1262, 0.9, 6.0, True,
         "SS=0.9, kappa=6, OVR ON (current defaults)"),
        (790, 8482, 1262, 0.9, 6.0, False,
         "SS=0.9, kappa=6, OVR OFF (gamma=0)"),
        (790, 8482, 1262, 0.9, 1000.0, True,
         "SS=0.9, kappa=1000, OVR ON"),
        (790, 8482, 1262, 0.9, 1000.0, False,
         "SS=0.9, kappa=1000, OVR OFF"),
        (790, 8482, 1262, 0.9, 1e6, False,
         "SS=0.9, kappa=1e6 (hard coupling), OVR OFF"),

        # Same but SS=1.0
        (790, 8482, 1262, 1.0, 6.0, True,
         "SS=1.0, kappa=6, OVR ON"),
        (790, 8482, 1262, 1.0, 6.0, False,
         "SS=1.0, kappa=6, OVR OFF"),
        (790, 8482, 1262, 1.0, 1000.0, False,
         "SS=1.0, kappa=1000, OVR OFF"),
        (790, 8482, 1262, 1.0, 1e6, False,
         "SS=1.0, kappa=1e6 (hard coupling), OVR OFF"),

        # SS=0.5 (worst case)
        (790, 8482, 1262, 0.5, 6.0, True,
         "SS=0.5, kappa=6, OVR ON"),
        (790, 8482, 1262, 0.5, 1e6, False,
         "SS=0.5, kappa=1e6, OVR OFF"),
    ]

    for n_spl, n_sense, n_anti, ss, kap, ovr, label in scenarios:
        run_em_trace(n_spl, n_sense, n_anti, ss,
                     kappa=kap, use_ovr=ovr, label=label)


if __name__ == "__main__":
    main()
