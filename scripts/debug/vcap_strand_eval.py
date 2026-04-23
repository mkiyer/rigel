#!/usr/bin/env python3
"""Strand-channel-only diagnostic for gDNA calibration on VCaP mixtures.

For each sample:
  1. Run scan_and_buffer once (cached to pickle to allow reuse)
  2. Build region_stats
  3. Fit **strand-only EM**: identical to run_em but with llr_count = 0
     and llr_fl = None — the strand BB LLR is the only signal.
  4. Fit **full-joint EM**: the shipping calibration (count + strand + FL).
  5. Save a per-region feather with (k_sense, n_unspliced, n_spliced,
     mappable_bp, tx_strand, γ_full, γ_strand_only, γ_strand_alone_posterior).

Produces summary table across nominal-DNA titration.
"""
from __future__ import annotations
import argparse, json, pickle, time, sys
from dataclasses import replace
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/mkiyer/proj/rigel/src")
from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.native import detect_sj_strand_tag
from rigel.pipeline import scan_and_buffer
from rigel.calibration._stats import compute_region_stats
from rigel.calibration._em import (
    run_em, _StrandBBCache, _strand_llr_bb_rho, _strand_posterior,
    _compute_k_sense_valid, _estimate_kappa_G, _estimate_kappa_R,
    _logsumexp, _KAPPA_LO, _KAPPA_HI,
)


INDEX_DIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
SRC_ROOT = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human"
CACHE_ROOT = Path("/scratch/mkiyer_root/mkiyer0/shared_data/rigel_scratch/vcap_strand_v5/strand_eval")
CACHE_ROOT.mkdir(parents=True, exist_ok=True)


def get_scan_outputs(sample: str, index):
    """Scan BAM once, cache region_counts + fl_table + ss + mean_fl."""
    cache = CACHE_ROOT / f"scan_{sample}.pkl"
    if cache.exists():
        print(f"  [cache] loading scan from {cache}", flush=True)
        with open(cache, "rb") as f:
            return pickle.load(f)

    bam = f"{SRC_ROOT}/mctp_vcap_rna20m_{sample}/rigel/annotated.bam"
    cfg = BamScanConfig(sj_strand_tag="auto", include_multimap=True, n_scan_threads=10)
    cfg = replace(cfg, sj_strand_tag=detect_sj_strand_tag(bam))

    t0 = time.perf_counter()
    stats, sm, flm, buf, region_counts, fl_table = scan_and_buffer(bam, index, cfg)
    sm.finalize(); flm.build_scoring_models(); flm.finalize(prior_ess=200.0)
    print(f"  scan: {time.perf_counter()-t0:.1f}s  SS={sm.strand_specificity:.4f}", flush=True)
    del buf  # drop heavy fragment buffer

    payload = {
        "region_counts": region_counts,
        "fl_table": fl_table,
        "strand_specificity": float(sm.strand_specificity),
        "mean_fl": float(flm.global_model.mean),
        "n_fragments_total": int(stats.n_fragments),
    }
    with open(cache, "wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    return payload


def run_strand_only_em(
    stats_dict, strand_specificity, mean_frag_len,
    *, max_iter=50, tol=1e-4,
):
    """Strand-channel-only EM: no count LLR, no FL LLR, no spliced anchor.

    E-step posterior = sigmoid(logit(π_soft) + LLR_strand_i).
    M-step updates only π, π_soft, κ_G, κ_R; λ_G is derived from
    γ-weighted Poisson MLE (λ_G = Σ γ·k^u / Σ γ·E) for reporting only.
    """
    n_u = stats_dict["n_unspliced"]
    n_s = stats_dict["n_spliced"]
    E = stats_dict["mappable_bp"]
    n_regions = n_u.shape[0]

    mappable_floor = max(1.0, float(mean_frag_len))
    eligible = E >= mappable_floor

    # Seed gamma = π_soft for all eligible (no hard anchor — this is
    # STRAND only). Keep hard spliced anchor (k^s > 0 means RNA) as an
    # optional comparison.
    pi_soft = 0.5
    prev_pi = pi_soft
    gamma = np.full(n_regions, pi_soft, dtype=np.float64)
    cache = _StrandBBCache(stats_dict)
    cache.set_ss(strand_specificity)

    kappa_G = _KAPPA_HI
    kappa_R = _KAPPA_HI

    n_iter = 0
    converged = False
    for it in range(max_iter):
        n_iter = it + 1

        # E-step: gamma = sigmoid(logit(pi_soft) + llr_strand)
        llr_strand = _strand_llr_bb_rho(
            stats_dict, strand_specificity, kappa_G, kappa_R, cache=cache,
        )
        pi_safe = float(np.clip(pi_soft, 1e-12, 1.0 - 1e-12))
        log_odds = np.log(pi_safe / (1.0 - pi_safe)) + llr_strand
        log_odds = np.clip(log_odds, -500.0, 500.0)
        g_eligible = 1.0 / (1.0 + np.exp(-log_odds))
        gamma = np.where(eligible, g_eligible, pi_soft)

        # M-step
        w = gamma[eligible]
        pi = float(np.mean(w)) if w.size else 0.5
        soft_mask = eligible & (n_s == 0)
        if soft_mask.any():
            pi_soft = float(np.mean(gamma[soft_mask]))
        else:
            pi_soft = pi
        pi = float(np.clip(pi, 1e-12, 1.0 - 1e-12))
        pi_soft = float(np.clip(pi_soft, 1e-12, 1.0 - 1e-12))

        # κ updates (MLEs on unique (k, n) pairs)
        gamma_valid = gamma[cache.valid]
        kappa_G = _estimate_kappa_G(cache, gamma_valid)
        kappa_R = _estimate_kappa_R(cache, gamma_valid)

        delta = abs(pi_soft - prev_pi)
        if delta < tol:
            converged = True
            break
        prev_pi = pi_soft

    # Derive implied lambda_G from strand-only gamma
    w_all = gamma[eligible]
    ku_all = n_u[eligible]
    E_all = E[eligible]
    num = float(np.sum(w_all * ku_all))
    den = float(np.sum(w_all * E_all))
    lam_G = num / den if den > 1e-12 else 0.0

    return {
        "lam_G": lam_G, "pi": pi, "pi_soft": pi_soft,
        "kappa_G": kappa_G, "kappa_R": kappa_R,
        "gamma": gamma, "llr_strand": llr_strand,
        "n_iter": n_iter, "converged": converged,
        "n_eligible": int(eligible.sum()),
    }


def evaluate(sample: str, index, region_df):
    print(f"\n=== {sample} ===", flush=True)
    payload = get_scan_outputs(sample, index)
    rc = payload["region_counts"]
    fl = payload["fl_table"]
    ss = payload["strand_specificity"]
    mean_fl = payload["mean_fl"]
    total_frags = payload["n_fragments_total"]

    stats_dict = compute_region_stats(rc, region_df)

    # --- Full joint EM ---
    t0 = time.perf_counter()
    full = run_em(
        stats_dict, fl, ss,
        mean_frag_len=mean_fl, max_iterations=50, convergence_tol=1e-4,
    )
    t_full = time.perf_counter() - t0

    # --- Strand-only EM ---
    t0 = time.perf_counter()
    so = run_strand_only_em(stats_dict, ss, mean_fl)
    t_so = time.perf_counter() - t0

    # --- Strand-only posterior using the FULL-EM fitted κ (diagnostic) ---
    llr_with_full_k = _strand_llr_bb_rho(
        stats_dict, ss, full.kappa_G, full.kappa_R,
    )
    gamma_strand_full_k = _strand_posterior(llr_with_full_k, pi_0=full.pi_soft)

    # --- Key metric: implied gDNA fragment count & fraction ---
    n_u = stats_dict["n_unspliced"]
    n_s = stats_dict["n_spliced"]
    E = stats_dict["mappable_bp"]
    eligible = E >= max(1.0, mean_fl)

    # "Implied gDNA fragments" under each model = Σ γ_i · k^u_i  (spliced
    # are definitionally RNA; ineligible regions don't contribute).
    elig_mask = eligible
    def implied_frags(gamma):
        return float(np.sum(gamma[elig_mask] * n_u[elig_mask]))

    gdna_full = implied_frags(full.gamma)
    gdna_strand_only = implied_frags(so["gamma"])
    gdna_strand_k_from_full = implied_frags(gamma_strand_full_k)

    # ≈ "Truth-proxy" fraction = n_gdna / total_fragments
    # But we also need "lambda * E" total-expected-gDNA metric
    lam_full_exp = float(full.lam_G * np.sum(E[elig_mask]))
    lam_so_exp = float(so["lam_G"] * np.sum(E[elig_mask]))

    # Region-level assignment
    out = pd.DataFrame({
        "region_id": np.arange(len(n_u)),
        "ref": stats_dict["ref"],
        "n_u": n_u, "n_s": n_s, "mappable_bp": E,
        "tx_strand": stats_dict["tx_strand"], "n_pos": stats_dict["n_pos"],
        "gamma_full": full.gamma,
        "gamma_strand_only": so["gamma"],
        "gamma_strand_from_full_k": gamma_strand_full_k,
        "llr_strand_full_k": llr_with_full_k,
        "eligible": eligible,
    })
    out_path = CACHE_ROOT / f"regions_{sample}.feather"
    out.reset_index(drop=True).to_feather(out_path)

    summary = {
        "sample": sample,
        "ss": ss,
        "mean_fl": mean_fl,
        "n_fragments_total": total_frags,
        "n_eligible": int(eligible.sum()),
        # Full joint EM
        "full_lam_G": full.lam_G,
        "full_kappa_G": full.kappa_G,
        "full_kappa_R": full.kappa_R,
        "full_pi": full.pi,
        "full_pi_soft": full.pi_soft,
        "full_n_iter": full.n_iter,
        "full_converged": full.converged,
        "full_gdna_frags_implied": gdna_full,
        "full_gdna_frac_implied": gdna_full / total_frags if total_frags else 0.0,
        "full_lam_times_E": lam_full_exp,
        "full_frac_lam_times_E": lam_full_exp / total_frags if total_frags else 0.0,
        "full_time_s": t_full,
        # Strand-only EM
        "so_lam_G": so["lam_G"],
        "so_kappa_G": so["kappa_G"],
        "so_kappa_R": so["kappa_R"],
        "so_pi": so["pi"],
        "so_pi_soft": so["pi_soft"],
        "so_n_iter": so["n_iter"],
        "so_converged": so["converged"],
        "so_gdna_frags_implied": gdna_strand_only,
        "so_gdna_frac_implied": gdna_strand_only / total_frags if total_frags else 0.0,
        "so_lam_times_E": lam_so_exp,
        "so_frac_lam_times_E": lam_so_exp / total_frags if total_frags else 0.0,
        "so_time_s": t_so,
        # Full-κ strand-only posterior
        "strand_gdna_frac_with_full_k": gdna_strand_k_from_full / total_frags if total_frags else 0.0,
    }
    return summary


def nominal_fraction(sample: str) -> float:
    """Approx gDNA fraction from naming: rna20m_dna{X}m → X/(20+X)."""
    if not sample.startswith("dna"): return float("nan")
    x = int(sample[3:].rstrip("m"))
    return x / (20 + x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="+", required=True)
    ap.add_argument("--out", default=str(CACHE_ROOT / "strand_eval_summary.tsv"))
    args = ap.parse_args()

    index = TranscriptIndex.load(INDEX_DIR)
    region_df = index.region_df

    rows = []
    for s in args.samples:
        try:
            rows.append(evaluate(s, index, region_df))
        except Exception as e:
            print(f"  FAILED {s}: {e}", flush=True)
            import traceback; traceback.print_exc()
    df = pd.DataFrame(rows)
    df["nominal_gdna_frac"] = df["sample"].map(nominal_fraction)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nSaved: {args.out}", flush=True)

    cols = ["sample", "nominal_gdna_frac", "so_gdna_frac_implied",
            "full_gdna_frac_implied", "strand_gdna_frac_with_full_k",
            "so_kappa_G", "so_kappa_R", "full_kappa_G", "full_kappa_R",
            "ss", "so_n_iter", "full_n_iter"]
    print(df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
