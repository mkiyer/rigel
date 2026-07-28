"""What the screen actually REMOVES on real data, vs the truncation bias it introduces."""
from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import od_tb_core as T
from od_tb_agg import D, estimator_terms, seed_view, size_hist, pooled_bias
from od_tb_bias import ALPHA_BONF, min_n_rejectable

CACHE: dict = {}


def keep_flags(nc, sense_int, a_null, alpha):
    """Per-seed keep decision using the seed's own integer (n, k)."""
    n = np.rint(nc).astype(np.int64)
    k = np.rint(sense_int).astype(np.int64)
    keep = np.ones(n.size, bool)
    nmin = min_n_rejectable(a_null, alpha)
    if nmin < 0:
        return keep, 0
    cand = np.nonzero(n >= nmin)[0]
    for i in cand:
        ni = int(n[i])
        key = (ni, a_null, alpha)
        if key not in CACHE:
            CACHE[key] = T.two_sided_p(ni, a_null)
        ki = int(np.clip(k[i], 0, ni))
        keep[i] = CACHE[key][ki] >= alpha
    return keep, len(cand)


def main():
    print("=== WHAT THE SCREEN REMOVES ON REAL DATA (pooled MoM, no prior shrinkage) ===")
    print("    od_raw = pooled MoM on all seeds; od_scr = pooled MoM on survivors.")
    print("    'trunc bias' = the exact truncation-bias factor at od = od_scr (from tb_agg).\n")
    for a_null in (2.0, 3.0):
        print(f"--- a_null = {a_null:.0f} ---")
        print(f"  {'sample':10s} {'alpha':>9s} {'#rej':>6s} {'%seeds':>8s} {'%num':>9s} "
              f"{'%scale':>8s} {'od_raw':>8s} {'od_scr':>8s} {'delta':>9s} "
              f"{'bias@od_scr':>11s} {'bias/delta':>11s}")
        for name, d in D.items():
            s, t, w, k = seed_view(d)
            ex, sc, nc = estimator_terms(s, t, w, k)
            # the estimator's own gDNA-component sense count: for w<1 the gDNA share of `sense`
            # is not observed; these seeds have w ~ 1 so use the raw (n, k) as the screen would.
            od_raw = ex.sum() / sc.sum()
            u, c = size_hist(nc)
            for alpha in (ALPHA_BONF, 1e-4, 1e-3, 1e-2, 5e-2):
                keep, ncand = keep_flags(nc, s, a_null, alpha)
                nrej = int((~keep).sum())
                od_scr = ex[keep].sum() / sc[keep].sum() if sc[keep].sum() > 0 else np.nan
                delta = od_raw - od_scr
                od_eval = float(np.clip(max(od_scr, 1e-4), 1e-4, 0.2))
                bias_fac = pooled_bias(u, c, a_null, alpha, od_eval)["od_hat"] / od_eval
                bias_abs = (1.0 - bias_fac) * od_eval
                print(f"  {name:10s} {alpha:9.3g} {nrej:6d} {nrej/s.size:8.4%} "
                      f"{ex[~keep].sum()/ex.sum():9.2%} {sc[~keep].sum()/sc.sum():8.4%} "
                      f"{od_raw:8.4f} {od_scr:8.4f} {delta:+9.4f} "
                      f"{bias_abs:11.2e} {abs(bias_abs)/max(abs(delta),1e-12):11.2e}")
        print()


if __name__ == "__main__":
    main()
