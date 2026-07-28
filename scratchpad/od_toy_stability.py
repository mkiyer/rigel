"""READINESS CHECK — does the information-currency shrinkage stay STABLE on the small toy genomes?

Owner's concern (2026-07-28): *"some of our toy scenarios don't have a large enough genome, so we need
some stability just for the sake of our development cycle and tests, not for real data."*

That is the right thing to check, and it is NOT obvious which way it goes. The change swaps the shrinkage
currency from SEED NODES to PAIRS:

    shipped:   od = (n_seed_nodes*od_mom + 30*od0) / (n_seed_nodes + 30)
    proposed:  od = (I*od_mom          + W*od0)  / (I + W),   I = SUM N(N-1)/2,  W = 909

A toy has FEW seeds but its seeds may be DEEP (a small genome sequenced to the same depth concentrates
fragments), and pairs grow as N^2. So the toy could easily end up with MORE information in the new currency
than in the old one -- i.e. the prior would bind LESS on toys than it does today, which is the opposite of
the stability the owner wants.

Reported per synthetic condition: seeds, pairs, the prior's share of the posterior under each rule, the
resulting od, and whether the `denom <= 0` fallback fires (where od0 becomes the ENTIRE answer).

Run: OMP_NUM_THREADS=1 python scratchpad/od_toy_stability.py
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from od_outliers import CF, SUITE, capture  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration import gdna_strand as GS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
W_NEW = 909.0            # derived: 1/tau^2, tau^2 = maxent variance on [0, od_max] with mean od0
W_OLD = 30.0             # shipped, in SEED-NODE units


def stats(sense, total, weight, kappa, mu_half):
    """Return (n_seed_nodes, I_pairs, od_mom, fallback) using INTEGER counts (the float trap)."""
    t = np.rint(np.asarray(total, float))
    s = np.asarray(sense, float)
    if mu_half:
        w = np.clip(np.asarray(weight, float), 0.0, 1.0)
        node_mean = 0.5 * w + kappa * (1.0 - w)
        nc = np.rint(w * t)
        cm = 0.5
    else:
        node_mean = np.full(t.shape, kappa)
        nc = t
        cm = kappa
    ok = t > 0
    ex = (s - t * node_mean) ** 2 - t * node_mean * (1.0 - node_mean)
    sc = np.maximum(nc * (nc - 1.0), 0.0) * cm * (1.0 - cm)
    n_seed = int(np.sum(ok & (nc > 0)))
    I = float(np.sum(np.maximum(nc * (nc - 1.0), 0.0) / 2.0)[()]) if ok.any() else 0.0
    num, den = float(ex[ok].sum()), float(sc[ok].sum())
    fb = den <= 0
    od_mom = num / den if den > 0 else np.nan
    return n_seed, I, od_mom, fb


def shrink(od_mom, n, W, od0, ceil):
    if not np.isfinite(od_mom):
        return od0
    return float(np.clip((n * od_mom + W * od0) / (n + W), 0.0, ceil))


def main():
    cfg = PipelineConfig()
    od0 = GS.overdispersion_for_beta(cfg.calibration.gdna_strand_prior_alpha_beta)
    ceil = GS._MAX_OVERDISPERSION
    print(f"od0 = {od0:.6f}   ceiling = {ceil:.4f}   W_old = {W_OLD} seed-nodes   W_new = {W_NEW} pairs\n")
    ix = TranscriptIndex.load(str(SUITE / "rigel_index"))
    ra = RegionArrays.from_region_df(ix.region_df, ix.ref_name_to_id)
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

    print("=== SYNTHETIC (truth od = 0 by construction) — gDNA fit ===\n")
    print(f"{'condition':46s} {'seeds':>7s} {'pairs I':>11s} | {'prior% OLD':>11s} {'prior% NEW':>11s} | "
          f"{'od OLD':>8s} {'od NEW':>8s} {'fallback':>9s}")
    rows = []
    for cond in conds:
        inp = _scan_and_truth(SUITE, cond, ix, cfg, Path("/tmp/rigel_selfsolve"),
                              SUITE / "_selfsolve_cache")
        g = capture(lambda inp=inp: calibrate(
            inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
            np.asarray(inp["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, w, k, _ = g["gdna"]
        n, I, om, fb = stats(s, t, w, k, True)
        a = shrink(om, n, W_OLD, od0, ceil)
        b = shrink(om, I, W_NEW, od0, ceil)
        rows.append((cond, n, I, a, b, fb))
        print(f"{cond[5:]:46s} {n:7d} {I:11,.0f} | {100 * W_OLD / (n + W_OLD):10.2f}% "
              f"{100 * W_NEW / (I + W_NEW):10.2f}% | {a:8.4f} {b:8.4f} {'YES' if fb else '':>9s}")

    A = np.array([r[3] for r in rows])
    B = np.array([r[4] for r in rows])
    print(f"\n  gDNA od across 32 synthetic conditions (truth 0):")
    print(f"    OLD  mean {A.mean():.4f}  max {A.max():.4f}  n>0.05: {(A > 0.05).sum()}")
    print(f"    NEW  mean {B.mean():.4f}  max {B.max():.4f}  n>0.05: {(B > 0.05).sum()}")
    print(f"    max |NEW-OLD| = {np.abs(B - A).max():.4f}   fallbacks: {sum(r[5] for r in rows)}/32")

    print("\n\n=== THE STABILITY QUESTION: does the prior bind MORE or LESS on the toy? ===\n")
    pr_old = np.array([100 * W_OLD / (r[1] + W_OLD) for r in rows])
    pr_new = np.array([100 * W_NEW / (r[2] + W_NEW) for r in rows])
    print(f"  prior share OLD (seed-node currency): median {np.median(pr_old):6.3f}%  "
          f"range {pr_old.min():.3f}-{pr_old.max():.3f}%")
    print(f"  prior share NEW (pair currency):      median {np.median(pr_new):6.3f}%  "
          f"range {pr_new.min():.3f}-{pr_new.max():.3f}%")
    print(f"\n  => the prior binds {'MORE' if np.median(pr_new) > np.median(pr_old) else 'LESS'} "
          f"under the new currency on the toy.")

    print("\n\n=== REAL DATA CONTROL (the change must not move these) ===\n")
    print(f"{'sample':10s} {'seeds':>9s} {'pairs I':>13s} | {'od OLD':>8s} {'od NEW':>8s}")
    import pickle
    for name in sorted(p.stem for p in CF.glob("*.pkl")):
        d = pickle.load(open(CF / f"{name}.pkl", "rb"))
        g = capture(lambda d=d: calibrate(
            d["payload"], d["region_arrays"], d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]), np.asarray(d["rna_fl_pmf"]),
            dataclasses.replace(cfg.calibration, calib_refit_iters=0)))
        s, t, w, k, _ = g["gdna"]
        n, I, om, fb = stats(s, t, w, k, True)
        print(f"{name:10s} {n:9,d} {I:13,.0f} | {shrink(om, n, W_OLD, od0, ceil):8.4f} "
              f"{shrink(om, I, W_NEW, od0, ceil):8.4f}")


if __name__ == "__main__":
    main()
