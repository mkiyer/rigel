#!/usr/bin/env python
"""Standalone simulation of the locus EM for a single-exon gene.

Reproduces the exact EM logic from counter.py with known parameters
to understand why gDNA absorbs everything for single-exon genes.
"""
import numpy as np

_EM_LOG_EPSILON = 1e-300


def simulate_locus_em(
    n_units: int,
    mrna_ll: float,
    gdna_ll: float,
    nrna_ll: float,
    mrna_eff_len: float,
    gdna_eff_len: float,
    nrna_eff_len: float,
    mrna_unique: float = 0.0,
    nrna_init: float = 0.0,
    gdna_init: float = 0.0,
    prior_mrna: float = 0.5,
    prior_nrna: float = 0.0,  # zeroed for single-exon
    prior_gdna: float = 0.5,
    has_nrna_candidate: bool = True,
    max_iters: int = 100,
    verbose: bool = True,
):
    """
    n_units: number of fragments going through EM
    *_ll: log-likelihood for each component
    *_eff_len: effective length for each component
    """
    # Build per-unit candidate arrays
    # Each unit has up to 3 candidates: mRNA(0), nRNA(1), gDNA(2)
    components = []
    n_candidates_per_unit = 2 if not has_nrna_candidate else 3
    # For simplicity, treat all units identically
    offsets = [0]
    t_indices = []
    log_liks = []

    for u in range(n_units):
        t_indices.append(0)  # mRNA
        log_liks.append(mrna_ll)
        if has_nrna_candidate:
            t_indices.append(1)  # nRNA
            log_liks.append(nrna_ll)
        t_indices.append(2)  # gDNA
        log_liks.append(gdna_ll)
        offsets.append(len(t_indices))

    offsets = np.array(offsets, dtype=np.int64)
    t_indices = np.array(t_indices, dtype=np.int32)
    log_liks = np.array(log_liks, dtype=np.float64)
    seg_lengths = np.diff(offsets)

    unique_totals = np.array([mrna_unique, nrna_init, gdna_init])
    prior = np.array([prior_mrna, prior_nrna, prior_gdna])
    eff_len = np.array([mrna_eff_len, nrna_eff_len, gdna_eff_len])
    log_eff_len = np.log(eff_len)

    # Balanced initialization
    theta_init = unique_totals.copy()
    shadow_weights = np.array([theta_init[1], theta_init[2]])
    avg_shadow = float(shadow_weights.mean())
    if avg_shadow > 0 and theta_init[0] < avg_shadow:
        theta_init[0] = avg_shadow

    theta = theta_init + prior
    total = theta.sum()
    if total > 0:
        theta /= total

    if verbose:
        print(f"  Init: theta = {theta}")
        print(f"  unique_totals = {unique_totals}")
        print(f"  prior = {prior}")
        print(f"  eff_len = {eff_len}")
        print(f"  log_liks (per unit): mRNA={mrna_ll:.4f}, "
              f"nRNA={nrna_ll:.4f}, gDNA={gdna_ll:.4f}")
        print()

    for i in range(max_iters):
        # E-step
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len

        log_posteriors = log_liks + log_weights[t_indices]

        seg_max = np.maximum.reduceat(log_posteriors, offsets[:-1])
        log_posteriors -= np.repeat(seg_max, seg_lengths)
        posteriors = np.exp(log_posteriors)
        seg_sum = np.add.reduceat(posteriors, offsets[:-1])
        posteriors /= np.repeat(seg_sum, seg_lengths)

        # M-step
        em_totals = np.zeros(3, dtype=np.float64)
        np.add.at(em_totals, t_indices, posteriors)

        theta_new = unique_totals + em_totals + prior
        total = theta_new.sum()
        if total > 0:
            theta_new /= total

        delta = np.abs(theta_new - theta).sum()

        if verbose and (i < 5 or i % 10 == 0 or delta < 1e-8):
            print(f"  Iter {i:3d}: theta=[{theta_new[0]:.6f}, "
                  f"{theta_new[1]:.6f}, {theta_new[2]:.6f}] "
                  f"em_totals=[{em_totals[0]:.1f}, {em_totals[1]:.1f}, "
                  f"{em_totals[2]:.1f}] delta={delta:.2e}")

        theta = theta_new

        if delta < 1e-8:
            break

    alpha = unique_totals + em_totals + prior
    if verbose:
        print(f"\n  Converged: alpha = {alpha}")
        print(f"  mRNA count = {alpha[0]:.1f}")
        print(f"  nRNA count = {alpha[1]:.1f}")
        print(f"  gDNA count = {alpha[2]:.1f}")

        # Final posteriors for assignment
        log_weights = np.log(theta + _EM_LOG_EPSILON) - log_eff_len
        # Check one unit
        s, e = int(offsets[0]), int(offsets[1])
        unit_log_post = log_liks[s:e] + log_weights[t_indices[s:e]]
        unit_log_post -= unit_log_post.max()
        unit_post = np.exp(unit_log_post)
        unit_post /= unit_post.sum()
        print(f"\n  Per-unit posteriors: mRNA={unit_post[0]:.6f}, ", end="")
        if has_nrna_candidate:
            print(f"nRNA={unit_post[1]:.6f}, gDNA={unit_post[2]:.6f}")
        else:
            print(f"gDNA={unit_post[1]:.6f}")

        p_rna = sum(unit_post[j] for j in range(len(unit_post))
                     if t_indices[s + j] < 2)
        print(f"  p_rna_total = {p_rna:.6f}")
        print(f"  p_rna >= 0.5? {p_rna >= 0.5}")
        # How many units go to gDNA (threshold 0.5)?
        n_gdna = sum(1 for _ in range(n_units) if p_rna < 0.5)
        print(f"  Units → gDNA = {n_gdna}/{n_units}")

    return alpha


if __name__ == "__main__":
    # Scenario: single-exon gene, untrained strand models (SS=0.5)
    # Both strand likelihoods = 0.5, same insert model, same overlap
    log_strand = np.log(0.5)
    # Insert: uniform prior (no data) → log(1/(max_size+1))
    # Assume max_size=1000 → log_insert = log(1/1001) ≈ -6.909
    log_insert_uniform = np.log(1.0 / 1001)
    # Actually: with some training data, the insert model might be different
    # Let's try both scenarios

    print("="*72)
    print("Scenario A: Untrained insert models (both uniform)")
    print("="*72)
    ll_mrna = log_strand + log_insert_uniform  # + 0 (overlap=1)
    ll_gdna = log_strand + log_insert_uniform  # + 0 (splice_pen=1)
    ll_nrna = log_strand + log_insert_uniform
    print(f"log_lik: mRNA={ll_mrna:.4f}, gDNA={ll_gdna:.4f}")
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna, gdna_ll=ll_gdna, nrna_ll=ll_nrna,
        mrna_eff_len=301, gdna_eff_len=301, nrna_eff_len=301,
    )

    print("\n" + "="*72)
    print("Scenario B: Like A but gDNA eff_len slightly larger (350 vs 301)")
    print("="*72)
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna, gdna_ll=ll_gdna, nrna_ll=ll_nrna,
        mrna_eff_len=301, gdna_eff_len=350, nrna_eff_len=301,
    )

    print("\n" + "="*72)
    print("Scenario C: Like A but gDNA eff_len smaller (250 vs 301)")
    print("="*72)
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna, gdna_ll=ll_gdna, nrna_ll=ll_nrna,
        mrna_eff_len=301, gdna_eff_len=250, nrna_eff_len=301,
    )

    print("\n" + "="*72)
    print("Scenario D: What if mRNA overlap_frac is slightly < 1?")
    print("(e.g., some edge fragments with 0.95 overlap)")
    print("="*72)
    # Average overlap of 0.97 → overlap_exp=10 → 10*log(0.97) = -0.304
    ll_mrna_edge = log_strand + log_insert_uniform + 10*np.log(0.97)
    ll_gdna_edge = log_strand + log_insert_uniform  # no overlap penalty
    print(f"log_lik: mRNA={ll_mrna_edge:.4f}, gDNA={ll_gdna_edge:.4f}")
    print(f"  diff = {ll_mrna_edge - ll_gdna_edge:.4f}")
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna_edge, gdna_ll=ll_gdna_edge, nrna_ll=ll_mrna_edge,
        mrna_eff_len=301, gdna_eff_len=301, nrna_eff_len=301,
    )

    print("\n" + "="*72)
    print("Scenario E: mRNA overlap_frac=0.90 (more edge fragments)")
    print("="*72)
    ll_mrna_bad = log_strand + log_insert_uniform + 10*np.log(0.90)
    print(f"log_lik: mRNA={ll_mrna_bad:.4f}, gDNA={ll_gdna_edge:.4f}")
    print(f"  diff = {ll_mrna_bad - ll_gdna_edge:.4f}")
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna_bad, gdna_ll=ll_gdna_edge, nrna_ll=ll_mrna_bad,
        mrna_eff_len=301, gdna_eff_len=301, nrna_eff_len=301,
    )

    print("\n" + "="*72)
    print("Scenario F: Equal likelihoods, NO nRNA candidate")
    print("="*72)
    simulate_locus_em(
        n_units=500,
        mrna_ll=ll_mrna, gdna_ll=ll_gdna, nrna_ll=ll_nrna,
        mrna_eff_len=301, gdna_eff_len=301, nrna_eff_len=301,
        has_nrna_candidate=False,
    )
