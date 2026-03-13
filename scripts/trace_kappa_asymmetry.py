#!/usr/bin/env python3
"""
Study the interaction between κ (strand symmetry coupling) and gDNA asymmetry.

Key insight from user: Real gDNA at a single locus follows Binomial(n, 0.5),
which has high variance at small n. With n=100 gDNA fragments, the 95% CI
for sense count is [40, 60] — a 40/60 or 60/40 split is entirely normal.
With n=15, even 5/10 or 10/5 is expected.

High κ forces g_pos ≈ g_neg, which:
  - Helps when the split is truly 50/50 (simulation artifact)
  - HURTS when there's natural binomial variance (real data)
  - Causes gDNA UNDER-estimation when sense > antisense (excess
    sense frags leak into nRNA)

This script systematically explores:
  1. Fine-grained κ sweep (6 → 1000) at symmetric gDNA
  2. Explicit asymmetry: 100 gDNA at 50/50, 40/60, 30/70, 20/80, etc.
  3. Small-count regime: 15 gDNA frags with realistic binomial draws
  4. Multiple gDNA counts (15, 50, 100, 500, 2000) × asymmetry × κ
"""

import numpy as np
import math
import sys


def run_em_quiet(
    n_mrna_spliced_sense: int,
    n_unspliced_sense: int,
    n_unspliced_antisense: int,
    SS: float,
    prior_alpha: float = 0.01,
    prior_gamma: float = 1.0,
    kappa: float = 6.0,
    max_iter: int = 500,
    use_ovr: bool = True,
):
    """Run EM and return final state + convergence info (no printing)."""
    nc = 4
    log_SS = math.log(SS) if SS > 0 else -1e10
    log_1mSS = math.log(1 - SS) if SS < 1.0 else -1e10

    # Equivalence classes (intronic only for nRNA/gDNA competition)
    ec_intr_sense_ll = np.array([log_SS, 0.0])
    ec_intr_sense_idx = np.array([1, 2])
    ec_intr_anti_ll = np.array([log_1mSS, 0.0])
    ec_intr_anti_idx = np.array([1, 3])
    ec_spl_sense_ll = np.array([0.0])
    ec_spl_sense_idx = np.array([0])

    equiv_classes = [
        (ec_spl_sense_ll, np.array([1.0]), ec_spl_sense_idx, n_mrna_spliced_sense),
        (ec_intr_sense_ll, np.array([1.0, 1.0]), ec_intr_sense_idx, n_unspliced_sense),
        (ec_intr_anti_ll, np.array([1.0, 1.0]), ec_intr_anti_idx, n_unspliced_antisense),
    ]

    # OVR warm start
    unambig_totals = np.zeros(nc)
    unambig_totals[0] = n_mrna_spliced_sense
    eligible = np.ones(nc)
    theta_init = unambig_totals.copy()
    coverage_totals = np.zeros(nc)
    n_ambiguous = 0

    if use_ovr:
        for ll, cw, idx, n_frags in equiv_classes:
            if len(idx) <= 1:
                continue
            weights = cw * eligible[idx]
            row_sum = weights.sum()
            if row_sum == 0:
                row_sum = 1.0
            shares = weights / row_sum
            for j, c in enumerate(idx):
                theta_init[c] += shares[j] * n_frags
                coverage_totals[c] += shares[j] * n_frags
            n_ambiguous += n_frags

    prior = np.zeros(nc)
    if n_ambiguous > 0 and prior_gamma > 0:
        for i in range(nc):
            ovr_i = coverage_totals[i] * prior_gamma / n_ambiguous
            prior[i] = prior_alpha + ovr_i if eligible[i] > 0 else 0.0
    else:
        for i in range(nc):
            prior[i] = prior_alpha if eligible[i] > 0 else 0.0

    state = theta_init + prior
    state /= state.sum()

    for iteration in range(max_iter):
        log_weights = np.log(state + 1e-300)
        em_totals = np.zeros(nc)

        for ll, cw, idx, n_frags in equiv_classes:
            if n_frags == 0:
                continue
            log_post = ll + log_weights[idx]
            max_lp = log_post.max()
            exp_lp = np.exp(log_post - max_lp)
            post = exp_lp / exp_lp.sum()
            for j, c in enumerate(idx):
                em_totals[c] += post[j] * n_frags

        theta_new = unambig_totals + em_totals + prior

        if kappa > 2.0:
            R_pos = theta_new[2]
            R_neg = theta_new[3]
            R_g = R_pos + R_neg
            if R_g > 0:
                phi = (R_pos + kappa / 2 - 1) / (R_g + kappa - 2)
                phi = max(0.0, min(1.0, phi))
                theta_new[2] = phi * R_g
                theta_new[3] = (1 - phi) * R_g

        total = theta_new.sum()
        state_new = theta_new / total if total > 0 else state.copy()
        delta = np.abs(state_new - state).sum()
        state = state_new.copy()
        if delta < 1e-10:
            break

    total_frags = n_mrna_spliced_sense + n_unspliced_sense + n_unspliced_antisense
    return {
        "theta": state,
        "counts": state * total_frags,
        "gdna_total": (state[2] + state[3]) * total_frags,
        "nrna_count": state[1] * total_frags,
        "gpos_count": state[2] * total_frags,
        "gneg_count": state[3] * total_frags,
        "iterations": iteration + 1,
        "total_frags": total_frags,
    }


def make_scenario(n_nrna, n_gdna_pos, n_gdna_neg, SS,
                  n_mrna_spliced=790):
    """Build fragment counts for a scenario.

    n_nrna: total nRNA fragments in locus
    n_gdna_pos: sense gDNA fragments in locus
    n_gdna_neg: antisense gDNA fragments in locus
    SS: strand specificity

    Returns (n_spliced_sense, n_unspliced_sense, n_unspliced_anti)
    """
    # nRNA splits by SS
    nrna_sense = round(n_nrna * SS)
    nrna_anti = n_nrna - nrna_sense

    # Unspliced pools (nRNA + gDNA)
    n_unspliced_sense = nrna_sense + n_gdna_pos
    n_unspliced_anti = nrna_anti + n_gdna_neg

    return n_mrna_spliced, n_unspliced_sense, n_unspliced_anti


def print_header(title):
    print(f"\n{'='*80}")
    print(f"  {title}")
    print(f"{'='*80}")


def main():
    SS = 0.9
    n_nrna = 9024   # from TA1=100, NTA1=100 scenario
    n_mrna_spl = 790

    # ================================================================
    # STUDY 1: Fine-grained κ sweep at SYMMETRIC gDNA (360/360)
    # ================================================================
    print_header("STUDY 1: Fine-grained κ sweep — SYMMETRIC gDNA (360 pos / 360 neg)")
    print(f"  n_nrna={n_nrna}, n_gdna=720 (360/360), SS={SS}")
    print(f"  TRUE: nRNA=9024, g_pos=360, g_neg=360, gDNA_total=720")
    print()
    print(f"  {'κ':>8}  {'nRNA':>8}  {'g_pos':>8}  {'g_neg':>8}  "
          f"{'gDNA_tot':>8}  {'nRNA_err%':>9}  {'gDNA_err%':>9}  "
          f"{'g+/g-':>6}  {'iters':>5}")

    n_spl, n_sense, n_anti = make_scenario(n_nrna, 360, 360, SS, n_mrna_spl)
    true_nrna = 9024
    true_gdna = 720

    kappa_values = [2.0, 4.0, 6.0, 8.0, 10, 15, 20, 30, 50, 75, 100,
                    150, 200, 300, 500, 750, 1000]
    for kap in kappa_values:
        r = run_em_quiet(n_spl, n_sense, n_anti, SS, kappa=kap)
        nrna_err = 100 * (r["nrna_count"] - true_nrna) / true_nrna
        gdna_err = 100 * (r["gdna_total"] - true_gdna) / true_gdna
        ratio = r["gpos_count"] / r["gneg_count"] if r["gneg_count"] > 0.1 else float("inf")
        print(f"  {kap:8.0f}  {r['nrna_count']:8.1f}  {r['gpos_count']:8.1f}  "
              f"{r['gneg_count']:8.1f}  {r['gdna_total']:8.1f}  "
              f"{nrna_err:+9.1f}%  {gdna_err:+9.1f}%  "
              f"{ratio:6.2f}  {r['iterations']:5d}")

    # ================================================================
    # STUDY 2: Explicit asymmetry with 720 total gDNA, various splits
    # ================================================================
    print_header("STUDY 2: Explicit gDNA asymmetry — 720 total gDNA, varying pos/neg split")
    print(f"  n_nrna={n_nrna}, SS={SS}")
    print(f"  TRUE: nRNA=9024, gDNA_total=720 (regardless of split)")
    print()

    # Test at multiple κ values
    test_kappas = [6, 20, 50, 100, 200, 500, 1000]
    # Asymmetry ratios: pos/total fraction
    asym_fracs = [0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90,
                  0.45, 0.40, 0.35, 0.30, 0.20, 0.10]
    # Sort for display
    asym_fracs = sorted(set(asym_fracs))

    print(f"  {'pos/neg':>10}  {'true_gpos':>9}  {'true_gneg':>9}  ", end="")
    for kap in test_kappas:
        print(f"  κ={kap:<5}  ", end="")
    print()
    print(f"  {'':>10}  {'':>9}  {'':>9}  ", end="")
    for kap in test_kappas:
        print(f"  {'gDNA_err%':>8}", end="")
    print()

    for frac in asym_fracs:
        gpos = round(720 * frac)
        gneg = 720 - gpos
        label = f"{gpos}/{gneg}"
        print(f"  {label:>10}  {gpos:9d}  {gneg:9d}  ", end="")

        n_spl, n_sense, n_anti = make_scenario(n_nrna, gpos, gneg, SS, n_mrna_spl)

        for kap in test_kappas:
            r = run_em_quiet(n_spl, n_sense, n_anti, SS, kappa=kap)
            gdna_err = 100 * (r["gdna_total"] - 720) / 720
            print(f"  {gdna_err:+8.1f}%", end="")
        print()

    # ================================================================
    # STUDY 3: Small gDNA count regime (15, 50 fragments)
    # ================================================================
    for n_gdna_total in [15, 50, 100, 500]:
        print_header(
            f"STUDY 3: Small-count regime — {n_gdna_total} total gDNA fragments, "
            f"varying asymmetry"
        )
        print(f"  n_nrna={n_nrna}, SS={SS}")
        print()

        # For small counts, test finer asymmetry
        if n_gdna_total <= 50:
            fracs = [0.50, 0.40, 0.60, 0.33, 0.67, 0.20, 0.80, 0.10, 0.90]
        else:
            fracs = [0.50, 0.45, 0.55, 0.40, 0.60, 0.35, 0.65, 0.30, 0.70,
                     0.20, 0.80]
        fracs = sorted(set(fracs))

        print(f"  {'pos/neg':>10}  ", end="")
        for kap in test_kappas:
            print(f"  κ={kap:<5}  ", end="")
        print()
        print(f"  {'':>10}  ", end="")
        for kap in test_kappas:
            print(f"  {'gDNA_err%':>8}", end="")
        print()

        for frac in fracs:
            gpos = round(n_gdna_total * frac)
            gneg = n_gdna_total - gpos
            label = f"{gpos}/{gneg}"

            n_spl, n_sense, n_anti = make_scenario(n_nrna, gpos, gneg, SS, n_mrna_spl)
            print(f"  {label:>10}  ", end="")

            for kap in test_kappas:
                r = run_em_quiet(n_spl, n_sense, n_anti, SS, kappa=kap)
                gdna_err = 100 * (r["gdna_total"] - n_gdna_total) / n_gdna_total
                print(f"  {gdna_err:+8.1f}%", end="")
            print()

    # ================================================================
    # STUDY 4: Binomial variance study
    # What does the EXPECTED error look like when gDNA follows Binom(n, 0.5)?
    # ================================================================
    print_header("STUDY 4: Expected error under Binomial(n, 0.5) gDNA")
    print(f"  For each (n_gdna, κ), we simulate 1000 binomial draws of the")
    print(f"  pos/neg split and compute mean |gDNA error %| and mean signed error.")
    print(f"  n_nrna={n_nrna}, SS={SS}")
    print()

    np.random.seed(42)
    n_draws = 1000
    gdna_totals = [15, 30, 50, 100, 200, 500, 720]
    test_kappas_binom = [6, 20, 50, 100, 200, 500, 1000]

    print(f"  {'n_gdna':>7}  {'metric':>10}  ", end="")
    for kap in test_kappas_binom:
        print(f"  κ={kap:<5}  ", end="")
    print()

    for n_gdna in gdna_totals:
        # Draw binomial samples
        gpos_samples = np.random.binomial(n_gdna, 0.5, size=n_draws)
        gneg_samples = n_gdna - gpos_samples

        # For each κ, compute error distribution
        errors_by_kappa = {}
        for kap in test_kappas_binom:
            errs = []
            for gpos, gneg in zip(gpos_samples, gneg_samples):
                n_spl, n_sense, n_anti = make_scenario(
                    n_nrna, int(gpos), int(gneg), SS, n_mrna_spl)
                r = run_em_quiet(n_spl, n_sense, n_anti, SS, kappa=kap)
                err_pct = 100 * (r["gdna_total"] - n_gdna) / n_gdna
                errs.append(err_pct)
            errors_by_kappa[kap] = np.array(errs)

        # Print mean signed error
        print(f"  {n_gdna:7d}  {'mean_err%':>10}  ", end="")
        for kap in test_kappas_binom:
            print(f"  {errors_by_kappa[kap].mean():+8.1f}%", end="")
        print()

        # Print mean absolute error
        print(f"  {'':>7}  {'mae%':>10}  ", end="")
        for kap in test_kappas_binom:
            print(f"  {np.abs(errors_by_kappa[kap]).mean():8.1f}%", end="")
        print()

        # Print 90th percentile absolute error
        print(f"  {'':>7}  {'p90_ae%':>10}  ", end="")
        for kap in test_kappas_binom:
            print(f"  {np.percentile(np.abs(errors_by_kappa[kap]), 90):8.1f}%", end="")
        print()

        # Print std dev
        print(f"  {'':>7}  {'std%':>10}  ", end="")
        for kap in test_kappas_binom:
            print(f"  {errors_by_kappa[kap].std():8.1f}%", end="")
        print()
        print()

    # ================================================================
    # STUDY 5: The "sense > antisense" danger zone
    # When gDNA has MORE sense than antisense, high κ binds g_pos
    # down to g_neg, causing gDNA underestimation and nRNA leakage
    # ================================================================
    print_header("STUDY 5: Sense-heavy gDNA — the danger zone for high κ")
    print(f"  100 gDNA fragments, varying sense/antisense split")
    print(f"  Shows estimated nRNA and gDNA counts (true: nRNA={n_nrna}, gDNA=100)")
    print()

    n_gdna_test = 100
    splits = [(50, 50), (55, 45), (60, 40), (65, 35), (70, 30),
              (75, 25), (80, 20), (85, 15), (90, 10)]

    print(f"  {'split':>8}  {'κ':>5}  {'est_nRNA':>9}  {'est_gDNA':>9}  "
          f"{'nRNA_err%':>9}  {'gDNA_err%':>9}  {'est_g+':>7}  {'est_g-':>7}")

    detail_kappas = [6, 50, 100, 200, 500]
    for gpos, gneg in splits:
        n_spl, n_sense, n_anti = make_scenario(n_nrna, gpos, gneg, SS, n_mrna_spl)
        for kap in detail_kappas:
            r = run_em_quiet(n_spl, n_sense, n_anti, SS, kappa=kap)
            nrna_err = 100 * (r["nrna_count"] - n_nrna) / n_nrna
            gdna_err = 100 * (r["gdna_total"] - n_gdna_test) / n_gdna_test
            label = f"{gpos}/{gneg}" if kap == detail_kappas[0] else ""
            print(f"  {label:>8}  {kap:5.0f}  {r['nrna_count']:9.1f}  "
                  f"{r['gdna_total']:9.1f}  {nrna_err:+9.1f}%  "
                  f"{gdna_err:+9.1f}%  {r['gpos_count']:7.1f}  "
                  f"{r['gneg_count']:7.1f}")
        print()  # blank line between splits

    # ================================================================
    # STUDY 6: Varying SS at different κ (asymmetric gDNA, 60/40)
    # ================================================================
    print_header("STUDY 6: κ × SS interaction — 720 gDNA at 60/40 split (432/288)")
    print(f"  n_nrna={n_nrna}")
    print()

    ss_values = [0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.99, 1.00]
    sweep_kappas = [6, 20, 50, 100, 200, 500]
    gpos_6040 = 432
    gneg_6040 = 288

    print(f"  {'SS':>5}  ", end="")
    for kap in sweep_kappas:
        print(f"  κ={kap:<5}  ", end="")
    print()
    print(f"  {'':>5}  ", end="")
    for kap in sweep_kappas:
        print(f"  {'gDNA_err%':>8}", end="")
    print()

    for ss in ss_values:
        n_spl, n_sense, n_anti = make_scenario(n_nrna, gpos_6040, gneg_6040, ss, n_mrna_spl)
        print(f"  {ss:5.2f}  ", end="")
        for kap in sweep_kappas:
            r = run_em_quiet(n_spl, n_sense, n_anti, ss, kappa=kap)
            gdna_err = 100 * (r["gdna_total"] - 720) / 720
            print(f"  {gdna_err:+8.1f}%", end="")
        print()


if __name__ == "__main__":
    main()
