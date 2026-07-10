"""Dissect the calibration→EM bridge (assemble_priors + the IPR eff-len) on the stranded flagship.

The stranded-flagship calibration PRIOR is good (contained error −16K) yet the post-EM leak is large. This
locates the leak between calibration and the trusted C++ EM, per-locus, by comparing three gDNA fractions:

  prior_frac   = gdna_prior_count / (gdna_prior_count + rna_prior_count)   — what assemble_priors handed the EM
  oracle_frac  = gdna_expected    / (gdna_expected   + rna_expected)       — the truth (split-BAM origins)
  em_frac      = gdna_observed    / (gdna_observed   + rna_observed)       — what the EM produced

If prior_frac ≈ oracle but em_frac << prior → the EM/eff-len leaks a GOOD prior (downstream).
If prior_frac << oracle → assemble_priors under-states gDNA (the bridge mis-assembles).
Also reports gdna_eff_len_em (the IPR contraction) on the worst loci.

    OMP_NUM_THREADS=1 python scripts/debug/fb_bridge_dissect.py [cond]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SUITE = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
DEFAULT = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
_EPS = 1e-9


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    nf = pd.read_csv(SUITE / "net_flow_per_locus.tsv", sep="\t")
    nf = nf[nf["condition"] == cond].copy()
    loci = pd.read_feather(SUITE / cond / "rigel_out" / "loci.feather")
    nf = nf.merge(loci[["locus_id", "locus_span_bp", "n_regions_touched", "gdna_eff_len_em"]],
                  on="locus_id", how="inner", suffixes=("", "_l"))
    nf = nf[(nf["gdna_prior_count"] + nf["rna_prior_count"]) > _EPS]  # drop empty-prior loci (intergenic)

    gp, rp = nf["gdna_prior_count"].to_numpy(), nf["rna_prior_count"].to_numpy()
    ge, re_ = nf["gdna_expected"].to_numpy(), nf["rna_expected"].to_numpy()
    go, ro = nf["gdna_observed"].to_numpy(), nf["rna_observed"].to_numpy()
    nf["prior_frac"] = gp / np.maximum(gp + rp, _EPS)
    nf["oracle_frac"] = ge / np.maximum(ge + re_, _EPS)
    nf["em_frac"] = go / np.maximum(go + ro, _EPS)
    nf["leak"] = nf["net_gdna_to_rna"]

    print(f"=== {cond} : calibration→EM bridge dissection ({len(nf)} loci) ===")
    print(f"  Σ gdna_expected={ge.sum():,.0f}  gdna_observed={go.sum():,.0f}  "
          f"Σ net_gdna_to_rna (leak)={nf['leak'].sum():,.0f}")

    def corr(a, b):
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 3 or a[m].std() < _EPS or b[m].std() < _EPS:
            return float("nan")
        return float(np.corrcoef(a[m], b[m])[0, 1])
    pf, of, ef = nf["prior_frac"].to_numpy(), nf["oracle_frac"].to_numpy(), nf["em_frac"].to_numpy()
    wt = (ge + re_)  # weight by locus size
    print("\n  gDNA-FRACTION agreement (mass-weighted mean |Δ|, and correlation):")
    print(f"    prior  vs oracle:  mean|Δ|={np.average(np.abs(pf-of), weights=wt):.3f}  corr={corr(pf, of):+.2f}"
          f"   ← is assemble_priors' SPLIT right?")
    print(f"    EM     vs oracle:  mean|Δ|={np.average(np.abs(ef-of), weights=wt):.3f}  corr={corr(ef, of):+.2f}")
    print(f"    EM     vs prior :  mean|Δ|={np.average(np.abs(ef-pf), weights=wt):.3f}  corr={corr(ef, pf):+.2f}"
          f"   ← does the EM track the prior, or leak below it?")
    # signed: does the EM systematically sit BELOW the prior (leak) ?
    print(f"    signed mean (EM_frac − prior_frac), mass-weighted = {np.average(ef-pf, weights=wt):+.3f}  "
          f"(negative ⇒ EM calls LESS gDNA than the prior ⇒ leak downstream of a good prior)")

    print("\n  TOP loci by |leak|  (prior→oracle→EM gDNA frac; gdna_eff_len_em = the IPR-contracted length):")
    print(f"  {'locus':>6} {'leak':>9} {'span':>7} {'nreg':>4} {'prior_f':>7} {'oracle_f':>8} {'em_f':>6} "
          f"{'gd_efflen':>9} {'gd_exp':>8} {'gd_obs':>8}")
    nf["absleak"] = nf["leak"].abs()
    for _, r in nf.sort_values("absleak", ascending=False).head(15).iterrows():
        print(f"  {int(r['locus_id']):>6} {r['leak']:>+9,.0f} {int(r['locus_span_bp']):>7} "
              f"{int(r['n_regions_touched']):>4} {r['prior_frac']:>7.2f} {r['oracle_frac']:>8.2f} "
              f"{r['em_frac']:>6.2f} {r['gdna_eff_len_em']:>9.1f} {r['gdna_expected']:>8,.0f} "
              f"{r['gdna_observed']:>8,.0f}")


if __name__ == "__main__":
    main()
