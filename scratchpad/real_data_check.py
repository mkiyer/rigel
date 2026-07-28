"""REAL DATA — the gDNA hyperprior + the P-2 pin fix on cached cfRNA payloads.

⚠ **Neither the landscape nor the pin fix has ever been run on real data** (HANDOFF_17 §4 step 4). There is
no oracle here, so this is not an accuracy test — it is the set of checks that CAN be made without one:

  * **it runs at all**, at genome scale, at both refit settings, with no NaN / degenerate output;
  * **the fitted landscape's SHAPE** — cfRNA is capture-enriched, so a bimodal `P(log ρ_g)` is the expected
    signature and a collapse to one mode would be a red flag (memory: `cfrna_sample_characteristics`);
  * **what the prior actually does to the answer** — refit=0 vs refit=1 on the same payload, which on the
    synthetic suite is worth −0.026 mwae and here can only be reported as a shift, not scored;
  * **the two standing debts' fitted scalars on real data** — `ω_graft` is on record as fitting **10× apart
    on two real samples** and as *expected to be fragile*; this prints it so the claim stays measured.

⚠ Held-out predictive likelihood is the only oracle-free *quality* criterion available on real data and it
is **biased toward under-smoothing** (it scores the landscape as a predictive object, which wants
population ⊛ noise, while its use as a prior wants the population). Not used to select anything here; the
toy-derived `knn_scale` is carried across, per the plan.

Run: OMP_NUM_THREADS=1 python scratchpad/real_data_check.py [--samples LBX0190 ...]
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import pickle
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
_EPS = 1e-12
_LN10 = np.log(10.0)


def load(name):
    d = pickle.load(open(CF / f"{name}.pkl", "rb"))
    return (d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]))


def landscape_shape(prior):
    """(n modes at a coarse prominence, the two component locations, mass above log10 rho = 0)."""
    if prior is None:
        return None
    x = np.asarray(prior.log_rho) / _LN10          # natural log -> decades
    p = np.exp(np.asarray(prior.logP) - np.max(prior.logP))
    p = p / max(p.sum(), _EPS)
    # deepest interior valley = the two-component split (the `two_component` convention)
    lo, hi = int(0.05 * p.size), int(0.95 * p.size)
    k = lo + int(np.argmin(p[lo:hi])) if hi > lo else p.size // 2
    a, b = p[:k], p[k:]
    loc_a = float(np.average(x[:k], weights=a)) if a.sum() > 0 else np.nan
    loc_b = float(np.average(x[k:], weights=b)) if b.sum() > 0 else np.nan
    return dict(split=float(x[k]), mass_lo=float(a.sum()), mass_hi=float(b.sum()),
                loc_lo=loc_a, loc_hi=loc_b, above0=float(p[x > 0.0].sum()),
                span=float(np.ptp(prior.logP)), n_train=int(prior.n_train))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", nargs="*",
                    default=[p.stem for p in sorted(CF.glob("*.pkl"))])
    a = ap.parse_args()
    cfg = PipelineConfig()
    print(f"{'sample':10s} {'refit':>5s} {'nodes':>8s} {'sec':>6s} {'f_g mass-wt':>12s} "
          f"{'rho_g global':>13s} {'kappa':>7s} {'od_g':>7s} {'od_r':>7s}")
    keep = {}
    for name in a.samples:
        payload, ra, sm, gp, rp = load(name)
        for refit in (0, 1):
            dbg: dict = {}
            cc = dataclasses.replace(cfg.calibration, calib_refit_iters=refit)
            t0 = time.time()
            res = calibrate(payload, ra, sm, gp, rp, cc, _debug=dbg)
            dt = time.time() - t0
            cap = dbg["capture"]
            mass = np.asarray(cap["mass_global"], float)
            fg = np.asarray(cap["f_g"], float)
            ok = mass > _EPS
            wf = float(np.average(fg[ok], weights=mass[ok]))
            prs = dbg["calibration_priors"]
            print(f"{name:10s} {refit:5d} {int(ok.sum()):8d} {dt:6.1f} {wf:12.4f} "
                  f"{res.gdna_density_global:13.6f} {res.rna_sense_frac:7.4f} "
                  f"{prs.gdna_strand_overdispersion:7.4f} {prs.rna_strand_overdispersion:7.4f}")
            assert np.isfinite(fg[ok]).all(), f"{name} r{refit}: non-finite f_g"
            assert (fg[ok] >= -1e-9).all() and (fg[ok] <= 1 + 1e-9).all(), f"{name} r{refit}: f_g range"
            keep[(name, refit)] = (fg, mass, dbg)

    print("\n\n=== the fitted LANDSCAPE on real data (from the refit=1 run's own Phase-2 fit) ===")
    print("   bimodality is the expected signature of hybrid capture; one mode would be a red flag.\n")
    print(f"{'sample':10s} {'n_train':>8s} {'split':>7s} {'loc lo':>8s} {'loc hi':>8s} "
          f"{'mass hi':>8s} {'>0 dec':>8s} {'logP span':>10s}")
    for name in a.samples:
        payload, ra, sm, gp, rp = load(name)
        dbg: dict = {}
        calibrate(payload, ra, sm, gp, rp,
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
        from rigel.calibration.calibrate import _fit_gdna_hyperprior
        cap = dbg["capture"]
        pr = _fit_gdna_hyperprior(
            dbg["chain"], dbg["belief"], dbg["statics"], ra,
            np.asarray(cap["mass_global"], float), np.asarray(cap["eff_global"], float),
            strength=cfg.calibration.gdna_prior_strength)
        sh = landscape_shape(pr)
        if sh is None:
            print(f"{name:10s} {'DECLINED':>8s}")
            continue
        print(f"{name:10s} {sh['n_train']:8d} {sh['split']:7.2f} {sh['loc_lo']:8.2f} "
              f"{sh['loc_hi']:8.2f} {sh['mass_hi']:8.3f} {sh['above0']:8.4f} {sh['span']:10.1f}")

    print("\n\n=== what the PRIOR does to the answer (refit 0 -> 1), by node class ===")
    print(f"{'sample':10s} {'mean |d f_g|':>13s} {'mass-wt f_g r0':>15s} {'r1':>8s} "
          f"{'nodes moved >0.05':>18s}")
    for name in a.samples:
        f0, m0, _ = keep[(name, 0)]
        f1, _, _ = keep[(name, 1)]
        ok = m0 > _EPS
        d = np.abs(f1 - f0)[ok]
        w = m0[ok]
        print(f"{name:10s} {np.average(d, weights=w):13.4f} "
              f"{np.average(f0[ok], weights=w):15.4f} {np.average(f1[ok], weights=w):8.4f} "
              f"{np.average(d > 0.05, weights=w):17.1%}")

    print("\n\n=== the standing debts' fitted scalars on real data (ROADMAP: expected to be FRAGILE) ===")
    for name in a.samples:
        _, _, dbg = keep[(name, 0)]
        glv = dbg["capture"].get("_glv", [])
        if not glv:
            print(f"{name:10s}  no graft edges recorded")
            continue
        om = [f"{g['omega']:.4f}" for g in glv[:4]]
        npair = sum(int(g["n_pairs"]) for g in glv)
        print(f"{name:10s}  omega_graft (per strand, first sweeps) = {', '.join(om)}   "
              f"n_pairs={npair}")


if __name__ == "__main__":
    main()
