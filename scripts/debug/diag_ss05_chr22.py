"""Diagnose ss=0.5 gdna=2.0 chr22 calibration underestimate."""
from __future__ import annotations

import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark.chr22_calibration_sweep import (  # noqa: E402
    run_one, DEFAULT_FASTA, DEFAULT_INDEX, REF_NAME,
    load_chr22_genome, select_transcripts,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.config import BamScanConfig  # noqa: E402
from rigel.pipeline import scan_and_buffer  # noqa: E402
from rigel.calibration import run_em  # noqa: E402
from rigel.calibration._stats import compute_region_stats  # noqa: E402


def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="%(message)s")
    out_dir = Path("/tmp/chr22_v4_diag_ss05")
    out_dir.mkdir(parents=True, exist_ok=True)

    index = TranscriptIndex.load(DEFAULT_INDEX)

    bam_path = out_dir / "scenarios" / "ss0.50_gdna2.00_seed1" / "oracle.bam"
    if not bam_path.exists():
        print("Generating BAM ...")
        genome = load_chr22_genome(DEFAULT_FASTA)
        transcripts = select_transcripts(index, seed=1, max_transcripts=500)
        run_one(
            out_dir,
            genome=genome, transcripts=transcripts, index=index,
            ss=0.5, gdna_fraction=2.0, seed=1,
            n_rna_fragments=50_000, rec_idx=0, keep_bam=True,
        )

    print(f"Scanning {bam_path} ...")
    scan_cfg = BamScanConfig(sj_strand_tag="XS")
    _, strand_models, fl_models, _, region_counts, fl_table = \
        scan_and_buffer(str(bam_path), index, scan_cfg)
    strand_models.finalize()
    fl_models.build_scoring_models()
    fl_models.finalize(prior_ess=1000.0)

    mean_fl = fl_models.global_model.mean
    print(f"  mean_frag_len = {mean_fl:.1f}")
    print(f"  strand_specificity (measured) = {strand_models.strand_specificity:.3f}")

    stats = compute_region_stats(region_counts, index.region_df)

    print("Running EM with diagnostics ...")
    fit = run_em(
        stats, fl_table,
        strand_specificity=strand_models.strand_specificity,
        mean_frag_len=mean_fl,
        max_iterations=50,
        convergence_tol=1e-4,
        diagnostics=True,
    )

    print(f"\n=== Fit ===")
    print(f"  λ_pool  = {fit.lam_G:.3e}")
    print(f"  μ_R = {fit.mu_R:.3f}, σ_R = {fit.sigma_R:.3f}")
    print(f"  π = {fit.pi:.3f}, π_soft = {fit.pi_soft:.3f}")
    print(f"  strand κ_G = {fit.kappa_G:.2f}, κ_R = {fit.kappa_R:.2f}")

    ref_arr = np.asarray(stats["ref"])
    mask22 = ref_arr == REF_NAME
    n22 = int(mask22.sum())
    E = stats["mappable_bp"]
    n_u = stats["n_unspliced"]
    n_s = stats["n_spliced"]
    gamma = fit.gamma

    E22 = E[mask22]
    nu22 = n_u[mask22]
    ns22 = n_s[mask22]
    g22 = gamma[mask22]
    eligible22 = E22 >= max(1.0, mean_fl)

    print(f"\n=== chr22 region inventory ({n22} total) ===")
    print(f"  eligible:                      {eligible22.sum()}")
    print(f"    n_s>0 (hard RNA):            {((ns22>0)&eligible22).sum()}")
    print(f"    soft n_s=0 total:            {((ns22==0)&eligible22).sum()}")
    print(f"    soft with n_u>0:             {((ns22==0)&(nu22>0)&eligible22).sum()}")

    print(f"\n=== chr22 soft region γ histogram ===")
    soft22 = eligible22 & (ns22 == 0)
    bins = np.array([-1e-9, 0.001, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.01])
    print(f"  {'γ range':<20} {'nreg':>6} {'Σk_u':>10} {'ΣE':>14} {'rate':>11}")
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        m = soft22 & (g22 >= lo) & (g22 < hi)
        nr = int(m.sum())
        if nr == 0:
            continue
        sn = float(nu22[m].sum())
        se = float(E22[m].sum())
        rate = sn / max(se, 1.0)
        print(f"  [{lo:.3f},{hi:.3f}) {nr:6d} {sn:10,.0f} {se:14,.0f} {rate:.3e}")

    print(f"\n=== chr22 λ-MLE diagnostics ===")
    ksum_soft = nu22[soft22].sum()
    esum_soft = E22[soft22].sum()
    print(f"  ALL soft eligible (γ-weighted): "
          f"Σγk/Σγe = {(g22*nu22)[soft22].sum()/max((g22*E22)[soft22].sum(),1e-6):.3e}")
    print(f"  ALL soft eligible (unweighted): "
          f"Σk/ΣE = {ksum_soft/max(esum_soft,1):.3e}")
    print(f"  truth on chr22 (100k / chr22_len): {100_000 / len(index.region_df[index.region_df['ref']==REF_NAME]['length'].sum()):.3e}" if False else f"  truth on chr22: ~1.97e-3")

    # Look at regions that SHOULD be gdna-dominated (no RNA overlap)
    # These are regions with tx_pos=0 and tx_neg=0
    region_df = index.region_df
    r22 = region_df[region_df["ref"] == REF_NAME].reset_index(drop=True)
    intergenic_mask = (~r22["tx_pos"].astype(bool)) & (~r22["tx_neg"].astype(bool))
    print(f"\n=== chr22 intergenic regions (tx_pos=0, tx_neg=0): {intergenic_mask.sum()} ===")
    inter_idx = np.where(mask22)[0][intergenic_mask.values]
    if len(inter_idx) > 0:
        E_i = E[inter_idx]
        nu_i = n_u[inter_idx]
        ns_i = n_s[inter_idx]
        g_i = gamma[inter_idx]
        elig_i = E_i >= max(1.0, mean_fl)
        print(f"  eligible:                {elig_i.sum()}")
        print(f"  Σn_u (eligible):         {nu_i[elig_i].sum():,.0f}")
        print(f"  ΣE   (eligible):         {E_i[elig_i].sum():,.0f}")
        print(f"  mean γ:                  {g_i[elig_i].mean():.4f}")
        print(f"  γ>0.5 count:             {(g_i[elig_i]>0.5).sum()}")
        print(f"  λ-MLE if γ=1: Σn/ΣE =   {nu_i[elig_i].sum()/max(E_i[elig_i].sum(),1):.3e}")
        # γ histogram on intergenic
        print("  γ hist on chr22 intergenic eligible:")
        for lo, hi in zip(bins[:-1], bins[1:]):
            m = elig_i & (g_i >= lo) & (g_i < hi)
            if m.any():
                print(f"    [{lo:.3f},{hi:.3f}): {int(m.sum()):6d} regions "
                      f"Σk={nu_i[m].sum():8,.0f} ΣE={E_i[m].sum():10,.0f}")

    print(f"\n=== EM trace ===")
    print("  iter   λ_G        μ_R      σ_R      π     π_soft  mean_γ_soft")
    for h in fit.history:
        print(f"  {h['iter']:4d}  {h['lam_G']:.3e}  {h['mu_R']:+.3f}  {h['sigma_R']:.3f}"
              f"  {h['pi']:.3f}  {h['pi_soft']:.3f}  {h['mean_gamma_soft']:.4f}")


if __name__ == "__main__":
    main()
