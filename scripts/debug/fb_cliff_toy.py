"""Toy: a PRINCIPLED, smooth (no-gate) fix for FB over-smoothing the hybrid-capture gDNA density cliff.

A minimal Gaussian-chain micro-model of the calibration message passing (the SCHEDULE — FB vs Jacobi — not the
lattice). Each node has a local belief (mean, prec). A message src->dst is `rho_src | theta ~ N(theta, v_msg)`,
`v_msg = 1/prec_src + sigma2_bio + pois` — a confident source speaks strongly, as in `node_sweep`.

PROBLEM (proved below): with a fixed sigma2_bio the message assumes gDNA density is spatially CONTINUOUS, so FB's
whole-chain relay smooths right across the capture cliff (enriched exon abutting a depleted intron). FB == Jacobi
at its fixed point — Jacobi only looked better because the benchmark under-converged it.

PRINCIPLED FIX (robust message, smooth — no gating): the per-edge spatial variance is itself uncertain and large
at a cliff. Put sigma2_comm ~ Scaled-Inv-chi2(nu, sigma2_bio); marginalizing gives a Student-t message ==
a Gaussian whose precision carries a smooth robust weight
    w = (nu + 1) / (nu + r^2),   r^2 = (rho_src - theta_hat_dst)^2 / v_msg
where theta_hat_dst is the dst's own (message-free) estimate. Agreement (r~0) -> w~1 (normal smoothing);
disagreement / cliff (r large) -> w->0 smoothly. The source-dest difference is now a factor in the model.

    python scripts/debug/fb_cliff_toy.py
"""
from __future__ import annotations

import numpy as np

S2 = 0.01     # sigma2_bio: cross-node communication variance (the fitted overdispersion)
POIS = 0.002  # sampling floor
N_EXON, N_INTRON = 5, 5
N = N_EXON + N_INTRON
TRUE = np.array([0.9] * N_EXON + [0.1] * N_INTRON)


def _wt(rho_src, theta_hat, v_msg, nu):
    """Smooth robust weight (Student-t E-step): 1 where src~dst, ->0 smoothly as they disagree. nu=None -> 1."""
    if nu is None:
        return 1.0
    r2 = (rho_src - theta_hat) ** 2 / max(v_msg, 1e-12)
    return (nu + 1.0) / (nu + r2)


def fb(mu0, t0, nu=None):
    """Forward-backward, one pass; robust message weight uses the dst's LOCAL mean as theta_hat."""
    fm, ft = mu0.copy(), t0.copy()
    amode, aprec = np.zeros(N), np.zeros(N)
    for i in range(N):
        if i - 1 >= 0:
            vmsg = 1.0 / max(ft[i - 1], 1e-12) + S2 + POIS
            w = _wt(fm[i - 1], mu0[i], vmsg, nu)
            aprec[i], amode[i] = w / vmsg, fm[i - 1]
            p = t0[i] + aprec[i]
            fm[i], ft[i] = (t0[i] * mu0[i] + aprec[i] * amode[i]) / p, p
    bm, bt = mu0.copy(), t0.copy()
    bmode, bprec = np.zeros(N), np.zeros(N)
    for i in range(N - 1, -1, -1):
        if i + 1 < N:
            vmsg = 1.0 / max(bt[i + 1], 1e-12) + S2 + POIS
            w = _wt(bm[i + 1], mu0[i], vmsg, nu)
            bprec[i], bmode[i] = w / vmsg, bm[i + 1]
            p = t0[i] + bprec[i]
            bm[i], bt[i] = (t0[i] * mu0[i] + bprec[i] * bmode[i]) / p, p
    p = t0 + aprec + bprec
    return (t0 * mu0 + aprec * amode + bprec * bmode) / p


def jacobi(mu0, t0, passes, nu=None):
    """K synchronous Jacobi passes; robust weight vs the dst LOCAL mean (the message-free anchor)."""
    m, t = mu0.copy(), t0.copy()
    for _ in range(passes):
        nm, nt = mu0.copy(), t0.copy()
        for i in range(N):
            num, den = t0[i] * mu0[i], t0[i]
            for j in (i - 1, i + 1):
                if 0 <= j < N:
                    vmsg = 1.0 / max(t[j], 1e-12) + S2 + POIS
                    pr = _wt(m[j], mu0[i], vmsg, nu) / vmsg
                    num, den = num + pr * m[j], den + pr
            nm[i], nt[i] = num / den, den
        m, t = nm, nt
    return m


def run(label, mu0, t0, nu):
    fbm = fb(mu0, t0, nu)
    jinf = jacobi(mu0, t0, 200, nu)
    ex, intr = slice(0, N_EXON), slice(N_EXON, N)
    nu_s = "Gaussian" if nu is None else f"robust nu={nu}"
    print(f"\n=== {label}  [{nu_s}] ===   exon true 0.9 | intron true 0.1")
    print("  exon means: " + " ".join(f"{x:.2f}" for x in fbm[ex]) + "   intron: " +
          " ".join(f"{x:.2f}" for x in fbm[intr]))
    print(f"  EXON mean  FB={fbm[ex].mean():.3f}  Jac@inf={jinf[ex].mean():.3f}   "
          f"INTRON mean FB={fbm[intr].mean():.3f}   (true exon 0.900 / intron 0.100)")


def main():
    mu_gw = np.array([0.30] * N_EXON + [0.10] * N_INTRON)  # genome-wide: exon weak, pulled to low rho_global
    t_gw = np.array([1.0] * N_EXON + [15.0] * N_INTRON)
    mu_ea = np.array([0.90] * N_EXON + [0.10] * N_INTRON)  # enrichment-aware: exon pinned near 0.9
    t_ea = np.array([8.0] * N_EXON + [15.0] * N_INTRON)

    print("########## THE PROBLEM (Gaussian message over-smooths; FB == Jacobi fixed point) ##########")
    run("genome-wide prior", mu_gw, t_gw, nu=None)
    run("enrichment-aware prior", mu_ea, t_ea, nu=None)

    print("\n########## THE FIX (robust message — smooth, no gate) ##########")
    for nu in (8, 4, 2):
        run("enrichment-aware prior", mu_ea, t_ea, nu=nu)
    run("genome-wide prior", mu_gw, t_gw, nu=2)

    print("\n########## SANITY: a SMOOTH chain (no cliff) must be unharmed by the robust weight ##########")
    mu_flat = np.full(N, 0.5)
    t_flat = np.full(N, 2.0)
    glob_true = np.full(N, 0.5)  # noqa: F841  (documents intent: uniform truth)
    run("uniform chain (true 0.5 everywhere)", mu_flat, t_flat, nu=None)
    run("uniform chain (true 0.5 everywhere)", mu_flat, t_flat, nu=2)


if __name__ == "__main__":
    main()
