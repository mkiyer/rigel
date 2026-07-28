"""Part 2: ground-truth validation (synthetic, TRUE od = 0) + the robust-trim (breakdown) curve.

Run: OMP_NUM_THREADS=1 python scratchpad/od_reject_frames2.py
"""

from __future__ import annotations

import pickle
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from od_reject_frames import (  # noqa: E402
    CEIL,
    PRIOR_OD,
    bh_threshold,
    eff_draws,
    info_weighted_od,
    pmin_two_sided,
    tarone_threshold,
    terms,
    two_sided_p,
    w_opt,
)

REAL = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_seeds.pkl"
SYNTH = "/Users/mkiyer/proj/rigel/scratchpad/od_reject_synth.pkl"
_EPS = 1e-12


def prep(sense, total, weight, kappa):
    excess, scale, nc, node_mean = terms(sense, total, weight, kappa)
    ok = (total > 0) & (nc > 0)
    return excess[ok], scale[ok], nc[ok], node_mean[ok], sense[ok], total[ok]


def pvals(nc, sense, total, node_mean, a):
    rho = 1.0 / (2.0 * a + 1.0)
    ni = np.rint(nc).astype(np.int64)
    ki = np.clip(np.rint(sense * np.where(total > 0, nc / np.maximum(total, _EPS), 1.0)), 0, ni)
    ki = ki.astype(np.int64)
    key = ni * 100003 + ki
    uk, inv = np.unique(key, return_inverse=True)
    un, ukk = uk // 100003, uk % 100003
    pu = np.where(un >= 2, two_sided_p(ukk, un, 0.5, rho), 1.0)
    pm = np.where(ni >= 2, pmin_two_sided(np.maximum(ni, 2), 0.5, rho), 1.0)
    return pu[inv], pm


def trim_curve(excess, scale, nc, rho_plug, fracs):
    """(d) BREAKDOWN frame: drop the seeds with the largest od_s until eps of the WEIGHT is gone."""
    ok = scale > 0
    ods, n = excess[ok] / scale[ok], nc[ok]
    w = w_opt(n, max(rho_plug, 1e-4))
    order = np.argsort(-ods)  # most-dispersed first (contamination is one-sided: od_s -> +1)
    cw = np.cumsum(w[order]) / np.sum(w)
    out = []
    for eps in fracs:
        cut = np.searchsorted(cw, eps)
        keep = order[cut:]
        out.append(float(np.sum(w[keep] * ods[keep]) / np.sum(w[keep])))
    return out


def one(tag, sense, total, weight, kappa, truth=None):
    excess, scale, nc, node_mean, sense, total = prep(sense, total, weight, kappa)
    m = excess.size
    pooled = float(np.sum(excess) / np.sum(scale))
    rho_iw, _, _ = info_weighted_od(excess, scale, nc, PRIOR_OD)
    d = eff_draws(nc[scale > 0], max(rho_iw, 1e-4))
    print(f"\n{tag}")
    print(
        f"   seeds={m:8d} medn={np.median(nc):5.1f} maxn={nc.max():6.0f} "
        f"| pooled={pooled:8.4f} | info-weighted={rho_iw:7.4f} | eff Beta-draws={np.sum(d):9.1f}"
        + (f" | TRUTH={truth}" if truth is not None else "")
    )
    for a in (2.0, 3.0):
        p, pm = pvals(nc, sense, total, node_mean, a)
        tt, kk = tarone_threshold(pm, 1.0)
        rej = p < tt
        keep = ~rej
        mom_t = float(np.sum(excess[keep]) / np.sum(scale[keep]))
        rho_t, _, _ = info_weighted_od(excess[keep], scale[keep], nc[keep], PRIOR_OD)
        bt, bk = bh_threshold(p, 0.05)
        print(
            f"     a={a:.0f}: Bonferroni/BH rejects {int(np.sum(p < 1.0 / m)):d}/{bk:d}   |"
            f"  TARONE(m_eff={kk:5d}, t={tt:.2e}) rejects {int(rej.sum()):5d} "
            f"({rej.mean():.3%} of seeds, {np.sum(scale[rej]) / np.sum(scale):6.2%} of pairs)"
            f" -> pooled {mom_t:8.4f}, info-w {rho_t:7.4f}"
        )
    fr = [0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50]
    tc = trim_curve(excess, scale, nc, rho_iw, fr)
    print("     TRIM (info-weighted, drop top-od_s until eps of WEIGHT removed):")
    print("       eps  " + " ".join(f"{f:8.3f}" for f in fr))
    print("       od   " + " ".join(f"{v:8.4f}" for v in tc))


def main():
    print("#" * 104)
    print("# GROUND TRUTH: synthetic suite, multinomial simulator => TRUE od = 0")
    print("#" * 104)
    syn = pickle.load(open(SYNTH, "rb"))
    for k, v in syn.items():
        s, t, w, kap, _ = v["gdna"]
        one("SYNTH  " + k, s, t, w, kap, truth=0.0)
    print()
    print("#" * 104)
    print("# REAL cfRNA")
    print("#" * 104)
    real = pickle.load(open(REAL, "rb"))
    for k in ("LBX0190", "LBX0588", "MO_3021", "vcap"):
        s, t, w, kap, _ = real[k]["gdna"]
        one("REAL   " + k, s, t, w, kap)

    # ---- why BH/Bonferroni are structurally dead: the p-value ATOMS -------------------------
    print()
    print("#" * 104)
    print("# WHY THE p-FRAME IS DEAD HERE: the attainable p-value atoms at small n (null Beta(3,3))")
    print("#" * 104)
    for n in (1, 2, 3, 4, 5, 10, 20):
        ps = sorted({round(float(two_sided_p(k, n, 0.5, 1 / 7.0)), 4) for k in range(n + 1)})
        print(f"   n={n:3d}: attainable two-sided p-values = {ps}")


if __name__ == "__main__":
    main()
