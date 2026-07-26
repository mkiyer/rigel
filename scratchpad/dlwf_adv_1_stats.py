"""ADVERSARY #1 — the STATISTICAL properties of the ONE-SIDED DerSimonian-Laird deflation as specified.

The DL law's algebra is right (b_hat^2 -> b^2). This script attacks how it is USED in this solver:

  A. the b=0 floor            E[max(0, chi2_1 - 1)] and the resulting precision cap
  B. one-sided vs symmetric   the plan inflates ONLY the message; honest random-effects inflates BOTH.
                              Who pays when the OWN belief is the biased one?
  C. double use               the own belief deflates the message AND is a fusion partner -> the fused
                              posterior's CLAIMED precision vs its ACTUAL precision (calibration)
  D. two messages             left+right both compared against the same own belief -> a wrong own belief
                              isolates the node by damping BOTH neighbours
  E. per-strand v_own         v_own,R = fg^2/tau is the RNA-TOTAL law; applying it per strand omits the
                              tilt variance -> which way does the error go, and how big?

Everything is in lambda/log space with Gaussian errors (the regime the derivation assumes), so any defect
found here is a defect of the ESTIMATOR'S USE, not of the linearization.

    OMP_NUM_THREADS=1 python scratchpad/dlwf_adv_1_stats.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np
from scipy import integrate, stats

rng = np.random.default_rng(20260725)


def hdr(s):
    print(f"\n{'='*100}\n{s}\n{'='*100}")


def sec_A(N):
    hdr("A. THE b=0 FLOOR — DL is positively biased at zero mismatch; the precision is CAPPED")
    # exact E[max(0, X-1)], X ~ chi2_1
    f = lambda x: (x - 1.0) * stats.chi2.pdf(x, 1)  # noqa: E731
    c, _ = integrate.quad(f, 1.0, np.inf)
    print(f"   E[max(0, chi2_1 - 1)] = {c:.6f}   (MC {np.mean(np.maximum(0.0, rng.chisquare(1, N) - 1)):.6f})")
    print("   => at TRUE b=0:  E[excess] = c*(v_msg + v_own)")
    print("      E[sigma^2_delivered] = v_msg + c*(v_msg+v_own) = (1+c)*v_msg + c*v_own")
    print("\n   effective precision at b=0, as a function of the own belief's weakness r = v_own/v_msg:")
    print("   (p_full = 1/v_msg; p_eff = E[1/(v_msg+excess)]; 'never damped' = P(chi2_1 < 1) = 0.6827)")
    v_msg = 0.05
    print(f"   {'r=v_own/v_msg':>14} {'E[p_eff]/p_full':>16} {'median damp':>12} {'P(damp<0.1)':>12} "
          f"{'E[sig2]/v_msg':>14}")
    for r_ in (0.0, 0.25, 1.0, 4.0, 16.0, 64.0, 256.0):
        v_own = r_ * v_msg
        G = rng.normal(0.0, np.sqrt(v_msg + v_own), N)
        exc = np.maximum(0.0, G * G - v_msg - v_own)
        damp = v_msg / (v_msg + exc)
        print(f"   {r_:14.2f} {np.mean(damp):16.4f} {np.median(damp):12.4f} "
              f"{np.mean(damp < 0.1):12.4f} {np.mean(v_msg+exc)/v_msg:14.3f}")
    print("\n   READ: the mean delivered precision loses ~32% at ANY r (the 1-0.6827 tail), and the")
    print("   distribution is BIMODAL: 68% of messages untouched, ~32% annihilated. The larger the own")
    print("   belief's weakness r, the harder the annihilated tail is hit (E[sig2] grows linearly in r).")
    print("   There is no crossover at which a WEAK own belief stops throttling: v_own enters excess with")
    print("   a + sign in the threshold but a + sign in E[G^2] too, and the two cancel EXACTLY -> the")
    print("   *rate* of spurious damping (31.7%) is INVARIANT to v_own; only its SEVERITY grows.")


def _fuse(p_o, lam_o, p_m, lam_m):
    p = p_o + p_m
    return (p_o * lam_o + p_m * lam_m) / np.maximum(p, 1e-300), p


def sec_BC(N):
    hdr("B/C. ONE-SIDED DL: who pays when the OWN belief is the biased study? + fusion calibration")
    v_msg, v_own = 0.05, 0.20
    print(f"   v_msg={v_msg}  v_own={v_own}   truth lambda=0; total between-source bias b split as")
    print("   b = b_msg - b_own  (DL can only see the GAP; it cannot tell which study is displaced).")
    print(f"\n   {'b_msg':>6} {'b_own':>6} | {'RMSE noDL':>10} {'RMSE plan':>10} {'RMSE symDL':>10} "
          f"{'RMSE own-only':>13} | {'plan claimed/actual prec':>24}")
    for b_msg, b_own in ((2.0, 0.0), (1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (0.0, 2.0), (1.5, -1.5), (1.0, 1.0)):
        lam_m = b_msg + rng.normal(0.0, np.sqrt(v_msg), N)
        lam_o = b_own + rng.normal(0.0, np.sqrt(v_own), N)
        G = lam_m - lam_o
        exc = np.maximum(0.0, G * G - v_msg - v_own)
        p_m0, p_o0 = 1.0 / v_msg, 1.0 / v_own
        p_m_dl = 1.0 / (v_msg + exc)
        p_o_dl = 1.0 / (v_own + exc)
        post_no, pr_no = _fuse(p_o0, lam_o, p_m0, lam_m)
        post_pl, pr_pl = _fuse(p_o0, lam_o, p_m_dl, lam_m)
        post_sy, pr_sy = _fuse(p_o_dl, lam_o, p_m_dl, lam_m)
        rmse = lambda x: float(np.sqrt(np.mean(x**2)))  # noqa: E731
        # calibration of the PLAN's claimed precision: claimed pr_pl vs actual 1/MSE
        cal = float(np.mean(pr_pl)) * float(np.mean(post_pl**2))
        print(f"   {b_msg:6.1f} {b_own:6.1f} | {rmse(post_no):10.4f} {rmse(post_pl):10.4f} "
              f"{rmse(post_sy):10.4f} {rmse(lam_o):13.4f} | {cal:24.2f}x")
    print("\n   READ: rows with b_own>0 (the OWN belief displaced, the message correct) are where the")
    print("   one-sided rule costs: it kills the only estimator that is right. 'claimed/actual prec' > 1")
    print("   means the fused posterior is OVER-CONFIDENT (the forbidden direction) -- the plan claims")
    print("   the own arm's undiminished precision even though the gap PROVED one of the two is wrong.")


def sec_D(N):
    hdr("D. TWO MESSAGES vs ONE OWN BELIEF — a wrong own belief ISOLATES the node")
    v_msg, v_own = 0.05, 0.20
    print(f"   v_msg={v_msg} (each side)  v_own={v_own};  left and right messages are INDEPENDENT and")
    print("   both UNBIASED for the truth; only the own belief is displaced by b_own.")
    print(f"\n   {'b_own':>6} | {'P(both damped<0.5)':>19} {'P(neither)':>11} | {'RMSE noDL':>10} "
          f"{'RMSE plan':>10} | {'post prec plan':>15}")
    for b_own in (0.0, 0.5, 1.0, 2.0, 3.0):
        lam_l = rng.normal(0.0, np.sqrt(v_msg), N)
        lam_r = rng.normal(0.0, np.sqrt(v_msg), N)
        lam_o = b_own + rng.normal(0.0, np.sqrt(v_own), N)
        e_l = np.maximum(0.0, (lam_l - lam_o) ** 2 - v_msg - v_own)
        e_r = np.maximum(0.0, (lam_r - lam_o) ** 2 - v_msg - v_own)
        p_l, p_r = 1.0 / (v_msg + e_l), 1.0 / (v_msg + e_r)
        d_l, d_r = v_msg * p_l, v_msg * p_r
        p_o = 1.0 / v_own
        post = (p_o * lam_o + p_l * lam_l + p_r * lam_r) / (p_o + p_l + p_r)
        post0 = (p_o * lam_o + lam_l / v_msg + lam_r / v_msg) / (p_o + 2.0 / v_msg)
        print(f"   {b_own:6.1f} | {np.mean((d_l<0.5)&(d_r<0.5)):19.3f} {np.mean((d_l>=0.5)&(d_r>=0.5)):11.3f} | "
              f"{np.sqrt(np.mean(post0**2)):10.4f} {np.sqrt(np.mean(post**2)):10.4f} | "
              f"{np.median(p_o+p_l+p_r):15.4g}")
    print("\n   READ: DL cannot use the AGREEMENT BETWEEN the two messages as evidence against the own")
    print("   belief -- each message is judged only against the own belief. Two independent neighbours")
    print("   agreeing with each other and disagreeing with 'own' is the strongest possible evidence that")
    print("   OWN is wrong, and the plan responds by killing both. RMSE grows monotonically in b_own.")


def sec_E(N):
    hdr("E. PER-STRAND v_own — the plan applies the RNA-TOTAL law v_own,R = fg^2/tau to EACH strand")
    print("   Truth: log f_pos = log f_R + log s_pos (s_pos = the tilt share), so")
    print("          Var(log f_pos) = Var(log f_R) + Var(log s_pos) + 2Cov  >=  Var(log f_R) when the tilt")
    print("   is a free DOF. The plan omits Var(log s_pos) => v_own too SMALL => excess too LARGE =>")
    print("   OVER-damping (the SAFE direction), by the additive amount Var(log s_pos).")
    print(f"\n   {'tilt sd (log s)':>16} {'v_own plan':>11} {'v_own true':>11} {'excess infl.':>13} "
          f"{'damp plan':>10} {'damp true':>10} {'over-damp x':>12}")
    v_R, v_msg = 0.20, 0.05
    for sd_t in (0.0, 0.1, 0.3, 0.6, 1.0):
        v_true = v_R + sd_t**2
        G = rng.normal(1.0, np.sqrt(v_msg + v_true), N)  # a real b=1 mismatch
        e_plan = np.maximum(0.0, G * G - v_msg - v_R)
        e_true = np.maximum(0.0, G * G - v_msg - v_true)
        d_plan = v_msg / (v_msg + e_plan)
        d_true = v_msg / (v_msg + e_true)
        print(f"   {sd_t:16.2f} {v_R:11.3f} {v_true:11.3f} {np.mean(e_plan-e_true):13.4f} "
              f"{np.median(d_plan):10.4f} {np.median(d_true):10.4f} "
              f"{np.median(d_true)/max(np.median(d_plan),1e-12):12.2f}")
    print("\n   READ: bounded and SAFE in direction, but note it only bites where the tilt is free AND")
    print("   tau_own>0 -- i.e. AMBIG nodes with a nonzero intron-factory tau. For single-strand nodes the")
    print("   tilt is structurally locked (Var(log s)=0) and the plan's law is EXACT there.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# DL one-sided-deflation statistics   draws={N:,}")
    sec_A(N)
    sec_BC(N)
    sec_D(N)
    sec_E(N)


if __name__ == "__main__":
    main()
