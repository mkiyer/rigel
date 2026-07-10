"""Ablation: is the flagship AMBIG-exon gDNA under-call owned by the GLOBAL prior — and is it the MEAN
(μ_global, the imputed gDNA fraction/direction) or the STRENGTH (N_global, the deference) that's broken?

Monkeypatches bp_solver._global_logprior (the count-space global) and re-runs production calibrate() on the
cached flagship payload; measures the per-region contained-gDNA error vs the by-origin oracle, by
strand_class × region_type. Calibration-prior level (no re-quant). Conditions:

  baseline            — the shipped global (seed ρ_global + exp-deference N_global).
  #1 oracle_ceiling   — global pins each region toward the ORACLE f_g at FULL strength (N=mass). Ceiling: if
                        this recovers the AMBIG exons, the global is the entire bottleneck.
  #2a oracle_mean     — ORACLE μ_global, but the CURRENT (deference) N_global. Isolates the MEAN/direction.
  #2b strong_strength — CURRENT μ_global, but FULL strength N=mass. Isolates the STRENGTH/deference.

    OMP_NUM_THREADS=1 python scripts/debug/ablation_global.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.simplex_sweep import _binom_pseudo  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}

# module-level ablation state read by the patched global
_ORACLE_MU = None  # per chain node: oracle f_g for REGION nodes, NaN for boundary nodes
_MODE = "baseline"
_orig_global = bp._global_logprior


def _patched_global(fgg, mass_global, eff_global, rho_global, sigma2_g, var_mean):
    cur_mu = np.clip(rho_global * eff_global / np.maximum(mass_global, _EPS), 0.0, 1.0)
    s2 = max(float(sigma2_g.predict(np.array([max(rho_global, _EPS)]))[0]), 0.0)
    cur_N = np.full_like(mass_global, rho_global * rho_global / max(var_mean + s2, _EPS))  # M-independent
    use_oracle = np.isfinite(_ORACLE_MU)  # only override region nodes; keep current on boundaries
    if _MODE == "oracle_ceiling":
        mu = np.where(use_oracle, np.nan_to_num(_ORACLE_MU), cur_mu)
        N = np.where(use_oracle, mass_global, cur_N)
    elif _MODE == "oracle_mean":
        mu = np.where(use_oracle, np.nan_to_num(_ORACLE_MU), cur_mu)
        N = cur_N
    elif _MODE == "strong_strength":
        mu = cur_mu
        N = np.where(use_oracle, mass_global, cur_N)
    else:
        return _orig_global(fgg, mass_global, eff_global, rho_global, sigma2_g, var_mean)
    return _binom_pseudo(fgg[None, :], mu[:, None], N[:, None])


def run(mode, blob, ra, g_or, scl, tcl):
    global _MODE
    _MODE = mode
    cal = calibrate(
        payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
        gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig(),
    )
    g_pr = np.asarray(cal.mass_gdna_contained, float)
    err = g_pr - g_or
    out = {"total_prod": float(g_pr.sum()), "net": float(err.sum())}
    for s, t in [(3, 2), (2, 2), (1, 2), (0, 0)]:  # AMBIG.exon, NEG.exon, POS.exon, NONE.intergenic
        m = (scl == s) & (tcl == t)
        out[f"{SC[s]}.{TC[t]}"] = float(err[m].sum())
    return out


def main():
    global _ORACLE_MU
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    r_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    oracle_fg = g_or / (g_or + r_or + _EPS)

    # build chain to map region index -> chain node, and set _ORACLE_MU (NaN on boundary nodes)
    from rigel.calibration.node_chain import REGION, build_node_chain
    payload = blob["payload_full"]
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    kind = np.asarray(chain.kind)
    refidx = np.asarray(chain.ref_idx)
    mu = np.full(chain.n_nodes, np.nan)
    reg_mask = kind == REGION
    mu[reg_mask] = oracle_fg[refidx[reg_mask]]
    _ORACLE_MU = mu

    sig = ra.signature
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    bp._global_logprior = _patched_global  # install the switch (baseline mode forwards to original)

    oracle_ambig = float(g_or[(scl == 3) & (tcl == 2)].sum())
    print(f"=== flagship {COND} : global-prior ablations (calibration-prior contained gDNA) ===")
    print(f"oracle: total contained gDNA={g_or.sum():,.0f}  AMBIG.exon gDNA={oracle_ambig:,.0f}\n")
    cols = ["AMBIG.exon", "NEG.exon", "POS.exon", "NONE.intergenic"]
    print(f"{'mode':>16} {'total_prod':>12} {'net_err':>10} | " + " ".join(f"{c:>14}" for c in cols))
    print("-" * 96)
    for mode in ["baseline", "oracle_ceiling", "oracle_mean", "strong_strength"]:
        o = run(mode, blob, ra, g_or, scl, tcl)
        print(f"{mode:>16} {o['total_prod']:>12,.0f} {o['net']:>+10,.0f} | "
              + " ".join(f"{o[c]:>+14,.0f}" for c in cols))
    print(f"\n(net_err columns = prod − oracle; AMBIG.exon target is 0 against oracle {oracle_ambig:,.0f}.)")


if __name__ == "__main__":
    main()
