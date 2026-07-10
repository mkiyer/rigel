"""WHERE does a clean boundary (u_pos=2, u_neg=191 → strand says +RNA, f_g≈0) get corrupted to f_g≈0.68?

Captures the EXACT ψ inputs the converged sweep feeds the boundary node at 5000 (detected by its counts),
then rebuilds f_g term by term — strand only → +Jeffreys → +global → +gDNA msg → +RNA msg — and prints, for
each term, the f_g median AND the log-contribution at f_g=0 vs the corrupted f_g, so the culprit term is visible.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_boundary_corruption_trace.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambig_nascent_trace import build  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice  # noqa: E402


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)

        cap = {}
        orig_ll = bp._local_loglik

        def spy(u_pos, u_neg, *a, **k):
            # boundary 5000: ~2 pos / ~191 neg crossing, single-strand + (allow_pos, not allow_neg)
            if u_pos[0] < 5 and 185 <= u_neg[0] <= 196 and bool(a[2][0]) and not bool(a[3][0]):
                cap["u"] = (u_pos.copy(), u_neg.copy())
                cap["a"] = a
                cap["k"] = dict(k)
            return orig_ll(u_pos, u_neg, *a, **k)

        bp._local_loglik = spy
        calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
        bp._local_loglik = orig_ll
        if "a" not in cap:
            print("boundary node not captured")
            return

        u_pos, u_neg = cap["u"]
        a = cap["a"]
        kw = cap["k"]
        lat = _simplex_lattice(int(cfg.calibration.sweep_n_grid))
        fpg, fng, fgg = lat
        kappa, odg, odr = a[4], a[5], a[6]
        N = float(u_pos[0] + u_neg[0])

        # term inventory
        glp = np.asarray(kw["global_logprior"])
        mu_g = float(kw["gdna_imp_frac"][0])
        n_g = float(kw["gdna_imp_count"][0])
        mu_p = float(kw["rna_imp_frac"][0][0])
        n_p = float(kw["rna_imp_count"][0][0])

        def med(psi):
            return float(_fg_median(psi, fgg)[0])

        def at(psi, fg_target):
            # log-posterior (normalised) summed over lattice points nearest the target f_g
            post = psi - logsumexp(psi)
            j = int(np.argmin(np.abs(fgg - fg_target)))
            return float(post[0, j])

        # strand only
        p_strand = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat)
        # + Jeffreys
        p_jeff = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat,
                               strand_obs=kw.get("strand_obs"))
        # + global
        p_glob = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat,
                               strand_obs=kw.get("strand_obs"), global_logprior=glp)
        # + gDNA msg
        p_gd = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat,
                             strand_obs=kw.get("strand_obs"), global_logprior=glp,
                             gdna_imp_frac=kw["gdna_imp_frac"], gdna_imp_count=kw["gdna_imp_count"])
        # + RNA msg (full)
        p_full = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat,
                               strand_obs=kw.get("strand_obs"), global_logprior=glp,
                               gdna_imp_frac=kw["gdna_imp_frac"], gdna_imp_count=kw["gdna_imp_count"],
                               rna_imp_frac=kw["rna_imp_frac"], rna_imp_count=kw["rna_imp_count"])

        print(f"\n=== boundary node corruption trace  (u_pos={u_pos[0]:.0f} u_neg={u_neg[0]:.0f} N={N:.0f}, "
              f"κ={kappa:.4f}, od_g={odg:.3f}, od_r={odr:.3f}) ===")
        print(f"  messages fed in: gDNA μ_g={mu_g:.3f} N_g={n_g:.2f} ; RNA+ μ_p={mu_p:.3f} N_p={n_p:.2f}")
        print(f"  global term: logprior[f_g=0]={glp[0,0]:.2f}  logprior[f_g=.68]={glp[0,int(np.argmin(np.abs(fgg-0.68)))]:.2f}")
        print(f"\n  {'term added':>16} {'f_g median':>11} {'logpost@0':>10} {'logpost@.68':>12}")
        for name, psi in (("strand only", p_strand), ("+ Jeffreys", p_jeff), ("+ global", p_glob),
                          ("+ gDNA msg", p_gd), ("+ RNA msg (full)", p_full)):
            print(f"  {name:>16} {med(psi):>11.3f} {at(psi,0.0):>10.2f} {at(psi,0.68):>12.2f}")
        sc.cleanup()


if __name__ == "__main__":
    main()
