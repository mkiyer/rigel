#!/usr/bin/env python3
"""Phase 3a: VCaP titration validation of composite-Poisson capture-class EM.

Successor to ``vcap_capture_class_phase2.py``.  The Phase 2 binary
partition (``λ_G_on`` vs ``λ_G_off`` with a hard 50% overlap cut)
could not recover ``λ_on`` under probe-subset-region dilution: a 120 bp
probe inside a 10 kb exon region contributed 100 % of the region's
counts to the off-class, collapsing the estimated enrichment.

Phase 3a replaces the binary partition with a composite-Poisson
superposition model (see ``docs/calibration/capture_class_density_plan.md``
§3a):

    rate_G_i = λ_G_on · E_on_i + λ_G_off · E_off_i,
    E_on_i   = mappable_bp_i · (overlap_bp_i / length_i),
    E_off_i  = mappable_bp_i − E_on_i.

The 2-iteration nested inner EM analytically unmixes the two Poisson
channels (strict convexity ⇒ no optimizer needed).

This script:

  1. Loads cached scan outputs per VCaP sample (same pickles as Phase 2).
  2. Annotates each region with per-bp ``(E_on, E_off)`` via the
     Agilent V4 BED.
  3. Runs the full calibration EM twice: once WITHOUT ``capture_e``
     (single-rate baseline), once WITH.
  4. Compares the four Phase 2 predictions for the composite model:
       (1) λ_on / λ_off > 1 and rising with nominal gDNA spike.
       (2) full_gdna_frac with BED tracks nominal gDNA fraction better.
       (3) Per-region γ (with BED) correlates with strand-only γ ≥
           the no-BED γ.
       (4) dna00m results are not harmed by enabling capture.

Outputs a TSV summary and per-sample regions feather (γ_composite,
γ_no_bed, γ_strand_only).
"""
from __future__ import annotations

import argparse
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, "/home/mkiyer/proj/rigel/src")
from rigel.calibration import annotate_capture_class, run_em  # noqa: E402
from rigel.calibration._em import (  # noqa: E402
    _KAPPA_HI,
    _StrandBBCache,
    _estimate_kappa_G,
    _estimate_kappa_R,
    _strand_llr_bb_rho,
)
from rigel.calibration._stats import compute_region_stats  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

# --- Paths ----------------------------------------------------------------
INDEX_DIR = "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
SCAN_CACHE_ROOT = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/rigel_scratch/vcap_strand_v5/strand_eval"
)
OUT_ROOT = Path(
    "/scratch/mkiyer_root/mkiyer0/shared_data/rigel_scratch/vcap_capture_class_phase3a"
)
OUT_ROOT.mkdir(parents=True, exist_ok=True)
BED_PATH = "/scratch/mkiyer_root/mkiyer0/shared_data/rigel/agilent-v4-targets-ucsc.bed.gz"


def _load_scan(sample: str) -> dict:
    cache = SCAN_CACHE_ROOT / f"scan_{sample}.pkl"
    if not cache.exists():
        raise FileNotFoundError(
            f"Missing scan cache {cache}. Run scripts/debug/vcap_strand_eval.py "
            f"first to populate it."
        )
    with open(cache, "rb") as fh:
        return pickle.load(fh)


def _run_strand_only_em(
    stats_dict, strand_specificity, mean_frag_len,
    *, max_iter=30, tol=1e-4, kappa_gdna_min=3.0,
):
    """Minimal strand-only EM (mirrors vcap_strand_eval.run_strand_only_em)."""
    n_u = stats_dict["n_unspliced"]
    n_s = stats_dict["n_spliced"]
    E = stats_dict["mappable_bp"]
    n_regions = n_u.shape[0]
    mappable_floor = max(1.0, float(mean_frag_len))
    eligible = E >= mappable_floor

    pi_soft = 0.5
    gamma = np.full(n_regions, pi_soft, dtype=np.float64)
    cache = _StrandBBCache(stats_dict)
    cache.set_ss(strand_specificity)
    kappa_G = kappa_R = _KAPPA_HI
    prev_pi = pi_soft

    for _ in range(max_iter):
        llr_strand = _strand_llr_bb_rho(
            stats_dict, strand_specificity, kappa_G, kappa_R, cache=cache,
        )
        pi_safe = float(np.clip(pi_soft, 1e-12, 1.0 - 1e-12))
        log_odds = np.log(pi_safe / (1.0 - pi_safe)) + llr_strand
        log_odds = np.clip(log_odds, -500.0, 500.0)
        gamma = np.where(eligible, 1.0 / (1.0 + np.exp(-log_odds)), pi_soft)

        w = gamma[eligible]
        pi = float(np.mean(w)) if w.size else 0.5
        soft_mask = eligible & (n_s == 0)
        pi_soft = float(np.mean(gamma[soft_mask])) if soft_mask.any() else pi
        pi_soft = float(np.clip(pi_soft, 1e-12, 1.0 - 1e-12))

        gamma_valid = gamma[cache.valid]
        kappa_G = _estimate_kappa_G(cache, gamma_valid, kappa_lo=float(kappa_gdna_min))
        kappa_R = _estimate_kappa_R(cache, gamma_valid)
        if abs(pi_soft - prev_pi) < tol:
            break
        prev_pi = pi_soft
    return gamma, float(pi), float(pi_soft), float(kappa_G), float(kappa_R)


def _nominal_gdna_fraction(sample: str) -> float:
    """rna20m_dna{X}m → X/(20+X) (spike mass fraction, not true gDNA frac)."""
    if not sample.startswith("dna"):
        return float("nan")
    try:
        x = int(sample[3:].rstrip("m"))
    except ValueError:
        return float("nan")
    return x / (20 + x)


def evaluate(
    sample: str,
    region_df: pd.DataFrame,
    E_on: np.ndarray,
    E_off: np.ndarray,
) -> dict:
    print(f"\n=== {sample} ===", flush=True)
    payload = _load_scan(sample)
    rc = payload["region_counts"]
    fl = payload["fl_table"]
    ss = float(payload["strand_specificity"])
    mean_fl = float(payload["mean_fl"])
    total_frags = int(payload["n_fragments_total"])

    stats_dict = compute_region_stats(rc, region_df)
    n_u = stats_dict["n_unspliced"]
    E = stats_dict["mappable_bp"]
    eligible = E >= max(1.0, mean_fl)
    tx_strand = stats_dict["tx_strand"]
    valid = eligible & (tx_strand != 0)
    valid_frags = float(n_u[valid].sum())

    # --- Strand-only γ (oracle baseline for per-region correlation check) ---
    t0 = time.perf_counter()
    gamma_strand, _, _, so_kappa_G, _ = _run_strand_only_em(
        stats_dict, ss, mean_fl, kappa_gdna_min=3.0,
    )
    t_strand = time.perf_counter() - t0

    # --- Full EM WITHOUT capture_e (single-rate baseline) ---
    t0 = time.perf_counter()
    fit_nobed = run_em(
        stats_dict, fl, ss,
        mean_frag_len=mean_fl, max_iterations=50, convergence_tol=1e-4,
        kappa_gdna_min=3.0,
    )
    t_nobed = time.perf_counter() - t0

    # --- Full EM WITH composite-Poisson capture_e = (E_on, E_off) ---
    t0 = time.perf_counter()
    fit_bed = run_em(
        stats_dict, fl, ss,
        mean_frag_len=mean_fl, max_iterations=50, convergence_tol=1e-4,
        kappa_gdna_min=3.0,
        capture_e=(E_on, E_off),
    )
    t_bed = time.perf_counter() - t0

    # --- Implied gDNA fragment counts / fractions ---------------------
    def implied_frags(gamma):
        return float(np.sum(gamma[valid] * n_u[valid]))

    gdna_nobed = implied_frags(fit_nobed.gamma)
    gdna_bed = implied_frags(fit_bed.gamma)
    gdna_strand = implied_frags(gamma_strand)

    # Expected gDNA via per-region λ·E.  No-BED: global λ · E.
    # Composite: (λ_on·E_on + λ_off·E_off).
    lam_nobed_E = float(fit_nobed.lam_G * np.sum(E[valid]))
    lam_bed_region = fit_bed.lam_G_on * E_on + fit_bed.lam_G_off * E_off
    lam_bed_E = float(np.sum(lam_bed_region[valid]))

    # --- Per-region γ correlations (on valid subset) ------------------
    v = valid & (n_u >= 2)
    if v.sum() > 50:
        rho_nobed_vs_strand = float(
            spearmanr(fit_nobed.gamma[v], gamma_strand[v]).correlation
        )
        rho_bed_vs_strand = float(
            spearmanr(fit_bed.gamma[v], gamma_strand[v]).correlation
        )
    else:
        rho_nobed_vs_strand = float("nan")
        rho_bed_vs_strand = float("nan")

    # Per-bp attribution summary.
    e_total = E_on + E_off
    on_frac_bp = float(E_on.sum() / max(e_total.sum(), 1.0))
    # Coverage-weighted split of observed counts by per-bp attribution:
    # k_on_attrib_i = k_i · (λ_on·E_on_i) / rate_G_i under the fitted
    # composite — the same unmixing the inner EM converges to.
    rate_G = fit_bed.lam_G_on * E_on + fit_bed.lam_G_off * E_off
    with np.errstate(divide="ignore", invalid="ignore"):
        on_attrib = np.where(
            rate_G > 0,
            n_u * (fit_bed.lam_G_on * E_on) / rate_G,
            0.0,
        )
    on_attrib_frac = float(on_attrib[valid].sum() / max(n_u[valid].sum(), 1.0))
    # Raw per-bp rates on on- vs off-dominant regions (diagnostic).
    on_dom = E_on > 0.5 * e_total
    off_dom = ~on_dom
    on_rate = float(
        n_u[on_dom & eligible].sum() / max(E[on_dom & eligible].sum(), 1.0)
    )
    off_rate = float(
        n_u[off_dom & eligible].sum() / max(E[off_dom & eligible].sum(), 1.0)
    )
    enrichment_ratio = (
        fit_bed.lam_G_on / fit_bed.lam_G_off
        if fit_bed.lam_G_off > 0 else float("nan")
    )

    # Persist per-region rows for plotting / deep dives.
    regions_out = pd.DataFrame({
        "region_id": np.arange(len(n_u)),
        "ref": stats_dict["ref"],
        "n_u": n_u, "mappable_bp": E,
        "tx_strand": tx_strand, "n_pos": stats_dict["n_pos"],
        "eligible": eligible,
        "E_on": E_on, "E_off": E_off,
        "gamma_strand_only": gamma_strand,
        "gamma_nobed": fit_nobed.gamma,
        "gamma_composite": fit_bed.gamma,
    })
    regions_out.reset_index(drop=True).to_feather(OUT_ROOT / f"regions_{sample}.feather")

    summary = {
        "sample": sample,
        "nominal_gdna_frac": _nominal_gdna_fraction(sample),
        "ss": ss,
        "mean_fl": mean_fl,
        "n_fragments_total": total_frags,
        "n_eligible": int(eligible.sum()),
        "n_strand_valid": int(valid.sum()),
        "valid_frags": valid_frags,
        # per-bp attribution
        "on_frac_bp": on_frac_bp,  # fraction of mappable bp attributed to on
        "on_attrib_frac_frags": on_attrib_frac,  # fraction of counts → on
        "on_dom_rate_k_over_E": on_rate,
        "off_dom_rate_k_over_E": off_rate,
        "observed_dom_rate_ratio": (
            on_rate / off_rate if off_rate > 0 else float("nan")
        ),
        # no-BED EM
        "nobed_lam_G": fit_nobed.lam_G,
        "nobed_kappa_G": fit_nobed.kappa_G,
        "nobed_kappa_R": fit_nobed.kappa_R,
        "nobed_pi": fit_nobed.pi,
        "nobed_pi_soft": fit_nobed.pi_soft,
        "nobed_n_iter": fit_nobed.n_iter,
        "nobed_converged": fit_nobed.converged,
        "nobed_gdna_frags_implied": gdna_nobed,
        "nobed_gdna_frac_implied": (
            gdna_nobed / valid_frags if valid_frags else 0.0
        ),
        "nobed_lam_times_E": lam_nobed_E,
        "nobed_time_s": t_nobed,
        # Composite-Poisson EM with capture_e
        "bed_lam_G_on": fit_bed.lam_G_on,
        "bed_lam_G_off": fit_bed.lam_G_off,
        "bed_lam_G_eff": fit_bed.lam_G,  # coverage-weighted effective rate
        "bed_enrichment_ratio_on_over_off": enrichment_ratio,
        "bed_mu_R": fit_bed.mu_R,
        "bed_sigma_R": fit_bed.sigma_R,
        "bed_kappa_G": fit_bed.kappa_G,
        "bed_kappa_R": fit_bed.kappa_R,
        "bed_pi": fit_bed.pi,
        "bed_pi_soft": fit_bed.pi_soft,
        "bed_n_iter": fit_bed.n_iter,
        "bed_converged": fit_bed.converged,
        "bed_gdna_frags_implied": gdna_bed,
        "bed_gdna_frac_implied": (
            gdna_bed / valid_frags if valid_frags else 0.0
        ),
        "bed_lam_times_E": lam_bed_E,
        "bed_time_s": t_bed,
        # strand-only baseline
        "strand_gdna_frags_implied": gdna_strand,
        "strand_gdna_frac_implied": (
            gdna_strand / valid_frags if valid_frags else 0.0
        ),
        "strand_kappa_G": so_kappa_G,
        "strand_time_s": t_strand,
        # per-region agreement with the strand-only oracle
        "rho_nobed_vs_strand": rho_nobed_vs_strand,
        "rho_bed_vs_strand": rho_bed_vs_strand,
    }
    # Compact console report
    print(
        f"  ΣE_on/ΣE={on_frac_bp:.3f}  "
        f"on-attrib count frac={on_attrib_frac:.3f}  "
        f"k/E on-dom={on_rate:.3e}  off-dom={off_rate:.3e}  "
        f"(obs ratio={summary['observed_dom_rate_ratio']:.1f}x)",
        flush=True,
    )
    print(
        f"  λ_G no-BED={fit_nobed.lam_G:.3e}  "
        f"λ_on={fit_bed.lam_G_on:.3e}  λ_off={fit_bed.lam_G_off:.3e}  "
        f"λ_eff={fit_bed.lam_G:.3e}  enrich={enrichment_ratio:.1f}x",
        flush=True,
    )
    print(
        f"  gDNA frac: nominal={summary['nominal_gdna_frac']:.3f}  "
        f"no-BED={summary['nobed_gdna_frac_implied']:.3f}  "
        f"composite={summary['bed_gdna_frac_implied']:.3f}  "
        f"strand-only={summary['strand_gdna_frac_implied']:.3f}",
        flush=True,
    )
    print(
        f"  ρ(γ_nobed, γ_strand)={rho_nobed_vs_strand:.3f}  "
        f"ρ(γ_composite, γ_strand)={rho_bed_vs_strand:.3f}",
        flush=True,
    )
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--samples", nargs="+",
        default=["dna00m", "dna01m", "dna02m", "dna05m",
                 "dna10m", "dna20m", "dna40m", "dna80m"],
    )
    ap.add_argument("--bed", default=BED_PATH)
    ap.add_argument("--pad", type=int, default=150)
    ap.add_argument("--out", default=str(OUT_ROOT / "phase3a_summary.tsv"))
    args = ap.parse_args()

    print(f"[load] index {INDEX_DIR}", flush=True)
    index = TranscriptIndex.load(INDEX_DIR)
    region_df = index.region_df

    print(f"[annotate] BED {args.bed}  pad={args.pad}", flush=True)
    t0 = time.perf_counter()
    E_on, E_off = annotate_capture_class(
        region_df, args.bed, pad=args.pad,
    )
    e_total = E_on + E_off
    print(
        f"  annotate: {time.perf_counter()-t0:.1f}s  "
        f"ΣE_on={E_on.sum():.3e}  ΣE_off={E_off.sum():.3e}  "
        f"ΣE_on/ΣE={E_on.sum()/max(e_total.sum(),1.0):.4f}  "
        f"n_regions_with_overlap={int((E_on > 0).sum()):,}",
        flush=True,
    )

    rows = []
    for s in args.samples:
        try:
            rows.append(evaluate(s, region_df, E_on, E_off))
        except Exception as e:  # noqa: BLE001
            print(f"  FAILED {s}: {e}", flush=True)
            import traceback
            traceback.print_exc()

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"\nSaved: {args.out}", flush=True)

    cols = [
        "sample", "nominal_gdna_frac",
        "nobed_gdna_frac_implied",
        "bed_gdna_frac_implied",
        "strand_gdna_frac_implied",
        "bed_lam_G_on", "bed_lam_G_off",
        "bed_enrichment_ratio_on_over_off",
        "rho_nobed_vs_strand", "rho_bed_vs_strand",
    ]
    print("\n" + df[cols].to_string(index=False))


if __name__ == "__main__":
    main()
