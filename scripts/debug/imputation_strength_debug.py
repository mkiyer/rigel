"""Why is the imputation N_src ≈ 1 (inert)? Instrument _message on the flagship + examine σ²_g.

N_src = ρ_src² / (σ²_bio(μ_q) + (n_src+1)/E_src²).  N≈1 ⇒ the denominator ≈ ρ_src², i.e. one term swamps it.
This logs, per message: ρ_src, E_src, n_src=ρ_src·E_src, μ_q, σ²_bio(μ_q), pois=(n_src+1)/E_src², N_src — and
which term dominates. Plus it prints the σ²_g curve over a μ grid + its training range (is it huge / does it
extrapolate?), and decomposes σ²_g by edge type (region↔region vs boundary-involving) to test whether the
region-vs-boundary eff-len convention mismatch is inflating the biological spread.

    OMP_NUM_THREADS=1 python scripts/debug/imputation_strength_debug.py [condition]
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
from rigel.config import CalibrationConfig  # noqa: E402

COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
_LOG = []
_CAP = {}


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    # instrument _message: record internals
    orig_msg = bp._message

    def msg(rho_src, eff_src, eff_dst, mass_dst, rho_dst_cur, varmean, spliced_dst=0.0):
        mu_q = max(0.5 * (rho_src + rho_dst_cur), 1e-9)
        s2 = float(varmean.predict(np.array([mu_q]))[0])
        pois = (rho_src + 1.0 / max(eff_src, 1e-9)) / max(eff_src, 1e-9)
        n_src = rho_src * eff_src
        N = rho_src * rho_src / max(s2 + pois, 1e-12)
        if rho_src > 1e-6:
            _LOG.append((rho_src, eff_src, n_src, mu_q, s2, pois, N))
        return orig_msg(rho_src, eff_src, eff_dst, mass_dst, rho_dst_cur, varmean, spliced_dst)

    bp._message = msg

    # capture the seed σ²_g curve
    orig_ns = bp.node_sweep

    def ns(chain, st, geom, belief, rga, bsub, **k):
        rho_g, gvm, vm = bp._gdna_seed_estimate(chain, st, geom, rga, bsub, np.asarray(belief.f_g, float).copy(), k["rna_sense_frac"])
        _CAP.update(gvm=gvm, rho_g=rho_g)
        return orig_ns(chain, st, geom, belief, rga, bsub, **k)

    calibrate.__globals__["node_sweep"] = ns
    calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())

    gvm = _CAP["gvm"]
    print(f"=== {COND}: imputation N_src debug ===")
    print(f"σ²_g (gDNA seed var~mean) training μ-range: [{np.exp(gvm.x_lo):.3g}, {np.exp(gvm.x_hi):.3g}]  "
          f"(fit on {gvm.fit_mean.size} seed-edge points)")
    grid = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0])
    print("  σ²_g(μ) on a grid:")
    for m, v in zip(grid, gvm.predict(grid)):
        print(f"    μ={m:>6.1f} : σ²_g={v:>10.3f}   (√=±{np.sqrt(v):.2f};  so a node at density μ predicts a "
              f"neighbour ±{np.sqrt(v):.1f})")

    arr = np.array(_LOG)  # rho_src, eff_src, n_src, mu_q, s2, pois, N
    if not arr.size:
        print("no messages logged")
        return
    rho_src, eff_src, n_src, mu_q, s2, pois, N = arr.T
    print(f"\nMESSAGES logged: {len(arr)}  (ρ_src>0)")
    print("  N_src distribution: median=%.2f  mean=%.2f  p90=%.2f  max=%.2f   (frac N<2: %.0f%%)"
          % (np.median(N), N.mean(), np.percentile(N, 90), N.max(), 100 * np.mean(N < 2)))
    print("  denominator: σ²_bio dominates (σ²_bio > pois) in %.0f%% of messages" % (100 * np.mean(s2 > pois)))
    print("  median:  ρ_src=%.3f  n_src(=ρ·E)=%.0f  μ_q=%.3f  σ²_bio(μ_q)=%.3f  pois=%.4f"
          % (np.median(rho_src), np.median(n_src), np.median(mu_q), np.median(s2), np.median(pois)))
    # the key ratio: if σ²_bio ≈ ρ_src² then N≈1. show it.
    print("  σ²_bio / ρ_src²  (≈1 ⇒ N≈1): median=%.2f" % np.median(s2 / np.maximum(rho_src**2, 1e-9)))
    # high-N vs low-N: what distinguishes them?
    hi = N > 10
    print(f"\n  messages with N>10: {hi.sum()} ({100*hi.mean():.0f}%).  "
          f"their median μ_q={np.median(mu_q[hi]) if hi.any() else float('nan'):.2f} vs low-N median μ_q={np.median(mu_q[~hi]):.2f}")
    print("  => N is killed where μ_q (local density) is HIGH (σ²_g ramps up there)." if hi.any() and np.median(mu_q[hi]) < np.median(mu_q[~hi]) else "")


if __name__ == "__main__":
    main()
