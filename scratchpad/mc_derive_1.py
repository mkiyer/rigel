"""Independent derivation agent #1 — per-component message-variance model, MC-validated.

DELTA METHOD on log-densities. Source per-component densities
    rho_g = f_g * M / E_g ,  rho_nu = f_R * M / E_R ,  rho_mu = S / E_spl   (measured, pure RNA)
with f_R = 1 - f_g the UNSPLICED-RNA (continue) fraction and the SINGLE composition dof
    lam = logit f_g ,  Var(lam) = 1/tau_lam .
=> d log f_g = (1-f_g) dlam ,  d log f_R = -f_g dlam  (perfectly ANTI-correlated: ONE dof).
The count M is a Poisson count SHARED by g and nu: Var(log M) = 1/n. The spliced brings its OWN 1/n_s.

Transport seed (source density log-variances):
    v_g   = (1-f_g)^2/tau + 1/n      (composition (+) sampling)
    v_nu  =    f_g^2 /tau + 1/n      (composition (+) sampling, SAME M)
    v_mu  =            1/n_s         (sampling only; composition-CERTAIN)

MESSAGE mode delivered to psi (per component c):
    f_c^dst = rho_c^src * r * E_c^dst / M_dst ,   r = rho_tot(dst)/rho_tot(src).
Because rho_tot(dst) = M_dst * B_dst, the explicit /M_dst CANCELS M_dst exactly:
    f_c^dst = rho_c^src * E_c^dst * B_dst / rho_tot(src)     (B_dst = dst composition bracket, fixed op-point)
=> the delivered mode is M_dst-independent  (divM-conversion law).

GRAFT (boundary->exon): matched set, RNA is the SUM rho_R = rho_nu + rho_mu, enrichment (r) CANCELS.
    Var(log f_g^dst) = [a_nu + a_mu(1-f_g)]^2 / tau + a_mu^2 (1/n + 1/n_s)        (graft-gDNA)
    Var(log f_R^dst) = a_g^2 [ (1 - w_mu f_g)^2 / tau + w_mu^2 (1/n + 1/n_s) ]    (graft-RNA-sum)
  shares a_c = rho_c/rho_tot ; within-RNA w_mu = rho_mu/rho_R, w_nu = rho_nu/rho_R.
  ANCHOR (f_g->1, no spliced): a_nu,a_mu->0 => graft-gDNA var -> 0 (composition CERTAIN); gDNA DENSITY
  var -> 1/n (finite). Ratio form k=rho_g/rho_R is singular here; per-component form is FINITE.

PEEL (exon->boundary): rho_nu = rho_R(x)/r - rho_mu is a DIFFERENCE; enrichment does NOT cancel.
    Var(log rho_nu) = u^2 (Var log rho_R(x) + Var log r) + (u-1)^2 (1/n_s) ,  u = T/rho_nu, T=rho_R(x)/r
  Var(log r) = sigma^2_transfer is LOAD-BEARING here (weights u>=1).

sigma^2_transfer = Var(log r) = Var(log rho_tot^dst) + Var(log rho_tot^src),
    Var(log rho_tot) = 1/n + [ (1/E_g - 1/E_r)/B ]^2 Var(f_g)     (enrichment_frame.composition_logvar)

COMBINE: psi receives the gDNA message on log f_g AND the RNA message on log f_R from ONE source as
INDEPENDENT Gaussians on the SAME lam dof. But they are a deterministic multiple of each other (rank-1,
correlation -1). Independent combine DOUBLES the lam-information => stated var = 1/2 the joint truth.
"""

from __future__ import annotations

import numpy as np

rng = np.random.default_rng(20260724)


def _report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<48} pred {pred:12.6g}   emp {emp:12.6g}   rel {rel:7.2%}")
    return rel, rel < tol


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


# ── GRAFT: boundary -> exon, RNA SUM, enrichment cancels ──────────────────────────────────────────────────
def graft(N, f_g, tau, n, n_s, E_g, E_r, E_spl, M, S, E_g_dst, E_r_dst, fgd=0.4, label=""):
    lam0 = np.log(f_g / (1.0 - f_g))
    lam = rng.normal(lam0, np.sqrt(1.0 / tau), N)  # composition dof
    Md = M * np.exp(rng.normal(-0.5 / n, np.sqrt(1.0 / n), N))  # node count (shared g & nu)
    Sd = S * np.exp(rng.normal(-0.5 / n_s, np.sqrt(1.0 / n_s), N))  # independent spliced count
    fg = _sig(lam)
    rho_g = fg * Md / E_g
    rho_nu = (1.0 - fg) * Md / E_r
    rho_mu = Sd / E_spl
    rho_tot = rho_g + rho_nu + rho_mu  # matched set, WITH spliced
    # dst composition bracket, fixed operating point (the lazy rho_tot iterate)
    B_dst = fgd / E_g_dst + (1.0 - fgd) / E_r_dst
    fgm = _sig(lam0)
    fg_dst = rho_g * E_g_dst * B_dst / rho_tot
    fR_dst = (rho_nu + rho_mu) * E_r_dst * B_dst / rho_tot
    emp_g = float(np.var(np.log(fg_dst)))
    emp_R = float(np.var(np.log(fR_dst)))

    # operating-point shares
    rg, rn, rm = fgm * M / E_g, (1.0 - fgm) * M / E_r, S / E_spl
    rtot = rg + rn + rm
    rR = rn + rm
    a_g, a_nu, a_mu = rg / rtot, rn / rtot, rm / rtot
    w_mu = rm / rR
    cnt = 1.0 / n + 1.0 / n_s
    pred_g = (a_nu + a_mu * (1.0 - fgm)) ** 2 / tau + a_mu**2 * cnt
    pred_R = a_g**2 * ((1.0 - w_mu * fgm) ** 2 / tau + w_mu**2 * cnt)
    print(f" GRAFT {label}  (f_g={fgm:.3f}, w_mu={w_mu:.3f}, a_g={a_g:.3f})")
    r1, ok1 = _report("graft-gDNA  Var(log f_g^dst)", pred_g, emp_g)
    r2, ok2 = _report("graft-RNA-sum  Var(log f_R^dst)", pred_R, emp_R)

    # divM-conversion: mode is M_dst-independent. Recompute f_g^dst the FULL way with a RANDOM M_dst and
    # r = rho_tot_dst/rho_tot_src, rho_tot_dst = M_dst*B_dst. Must equal the reduced form draw-by-draw.
    Mdst = M * np.exp(rng.normal(0.0, 1.0, N))  # arbitrary, large-spread destination count
    r_full = (Mdst * B_dst) / rho_tot
    fg_full = rho_g * r_full * E_g_dst / Mdst
    max_dev = float(np.max(np.abs(fg_full - fg_dst) / np.maximum(fg_dst, 1e-300)))
    print(f"      divM-conversion: max |Δlog f_g| over random M_dst (100x spread) = {max_dev:.2e}"
          f"   {'OK ' if max_dev < 1e-12 else '***'}")
    return [("graft-gDNA", pred_g, emp_g, r1, ok1),
            ("graft-RNA-sum", pred_R, emp_R, r2, ok2),
            ("divM-conversion", 0.0, max_dev, max_dev, max_dev < 1e-12)]


# ── ANCHOR limit: pure gDNA, no spliced ───────────────────────────────────────────────────────────────────
def anchor(N, f_g, tau, n, E_g, E_r, M):
    lam0 = np.log(f_g / (1.0 - f_g))
    lam = rng.normal(lam0, np.sqrt(1.0 / tau), N)
    Md = M * np.exp(rng.normal(-0.5 / n, np.sqrt(1.0 / n), N))
    fg = _sig(lam)
    rho_g = fg * Md / E_g
    rho_nu = (1.0 - fg) * Md / E_r
    rho_tot = rho_g + rho_nu
    # gDNA DENSITY log-variance -> should approach 1/n + (1-f_g)^2/tau  (finite; ->1/n as f_g->1)
    emp_dens = float(np.var(np.log(rho_g)))
    pred_dens = 1.0 / n + (1.0 - f_g) ** 2 / tau
    # graft-gDNA fraction variance -> a_nu^2/tau (composition), collapses toward 0 as f_g->1
    B_dst = f_g / E_g + (1.0 - f_g) / E_r
    fg_dst = rho_g * E_g * B_dst / rho_tot
    emp_frac = float(np.var(np.log(fg_dst)))
    a_nu = ((1.0 - f_g) * M / E_r) / (f_g * M / E_g + (1.0 - f_g) * M / E_r)
    pred_frac = a_nu**2 / tau
    print(f" ANCHOR  (f_g={f_g:.4f})")
    r1, ok1 = _report("anchor gDNA DENSITY Var(log rho_g)", pred_dens, emp_dens)
    r2, ok2 = _report("anchor gDNA FRACTION Var(log f_g^dst)", pred_frac, emp_frac)
    print(f"      (ratio form k=rho_g/rho_R would be SINGULAR: k={f_g/(1-f_g)*E_r/E_g:.3g};"
          f" per-component density var={pred_dens:.4g} is FINITE)")
    return [("anchor-limit-density", pred_dens, emp_dens, r1, ok1),
            ("anchor-limit-fraction", pred_frac, emp_frac, r2, ok2)]


# ── sigma^2_transfer = Var(log r) ─────────────────────────────────────────────────────────────────────────
def transfer_var(N, f_g_s, f_g_d, tau_s, tau_d, n_s_cnt, n_d_cnt, E_g, E_r, M_s, M_d):
    def draw_rhotot(f_g, tau, ncnt, Mmean):
        lam0 = np.log(f_g / (1.0 - f_g))
        lam = rng.normal(lam0, np.sqrt(1.0 / tau), N)
        Mv = Mmean * np.exp(rng.normal(-0.5 / ncnt, np.sqrt(1.0 / ncnt), N))
        fg = _sig(lam)
        return Mv * (fg / E_g + (1.0 - fg) / E_r), lam0

    rt_s, l0s = draw_rhotot(f_g_s, tau_s, n_s_cnt, M_s)
    rt_d, l0d = draw_rhotot(f_g_d, tau_d, n_d_cnt, M_d)
    emp = float(np.var(np.log(rt_d)) + np.var(np.log(rt_s)))  # independent => Var(log r) adds

    def clv(f_g, tau, ncnt):
        # composition_logvar: 1/n + [(1/E_g-1/E_r)/B]^2 Var(f_g).  Var(f_g) = (f_g(1-f_g))^2 Var(lam).
        B = f_g / E_g + (1.0 - f_g) / E_r
        var_fg = (f_g * (1.0 - f_g)) ** 2 / tau
        return 1.0 / ncnt + ((1.0 / E_g - 1.0 / E_r) / B) ** 2 * var_fg

    pred = clv(f_g_s, tau_s, n_s_cnt) + clv(f_g_d, tau_d, n_d_cnt)
    print(" TRANSFER  sigma^2_transfer = Var(log r) = Var(log rho_tot^dst)+Var(log rho_tot^src)")
    r1, ok1 = _report("transfer-var  Var(log r)", pred, emp)
    return [("transfer-var", pred, emp, r1, ok1)]


# ── PEEL: exon -> boundary, RNA DIFFERENCE, sigma^2_transfer load-bearing ──────────────────────────────────
def peel(N, rho_R_x, var_log_rhoR, r, var_log_r, rho_mu, n_s, label=""):
    Rx = rho_R_x * np.exp(rng.normal(-0.5 * var_log_rhoR, np.sqrt(var_log_rhoR), N))
    rr = r * np.exp(rng.normal(-0.5 * var_log_r, np.sqrt(var_log_r), N))
    rm = rho_mu * np.exp(rng.normal(-0.5 / n_s, np.sqrt(1.0 / n_s), N))
    T = Rx / rr
    nu = T - rm
    keep = nu > 0
    T0 = rho_R_x / r
    nu0 = T0 - rho_mu
    u = T0 / nu0
    var_logT = var_log_rhoR + var_log_r
    pred = u**2 * var_logT + (u - 1.0) ** 2 * (1.0 / n_s)
    emp = float(np.var(np.log(nu[keep])))
    print(f" PEEL {label}  (u={u:.2f}, kept {keep.mean():.1%})")
    r1, ok1 = _report("peel-diff  Var(log rho_nu)", pred, emp, tol=0.12)
    # ablate Var(log r): show it is load-bearing
    frac_r = u**2 * var_log_r / pred
    print(f"      Var(log r) carries {frac_r:.0%} of the peel variance (load-bearing)")
    return [("peel-diff", pred, emp, r1, ok1)]


# ── COMBINE: two independent messages on log f_g & log f_R from ONE source (rank-1) ────────────────────────
def combine(N, f_g, tau_T, E_g=100.0, E_r=100.0):
    """Source delivers ONE composition estimate (lam_hat ~ N(lam*, 1/tau_T)). It is expressed as two
    messages: log f_g and log f_R, each with its transported precision. psi combines them INDEPENDENTLY on
    the single lam axis. Compare psi's STATED posterior var to the JOINT-TRUTH var (1/tau_T)."""
    lam_star = np.log(f_g / (1.0 - f_g))
    # transported per-component log-fraction variances from the SINGLE dof (no count, pure composition):
    v_g = (1.0 - f_g) ** 2 / tau_T
    v_R = f_g**2 / tau_T
    p_g, p_R = 1.0 / v_g, 1.0 / v_R
    lam_grid = np.linspace(lam_star - 6.0 / np.sqrt(tau_T), lam_star + 6.0 / np.sqrt(tau_T), 4001)
    log_fg = np.log(_sig(lam_grid))
    log_fR = np.log(_sig(-lam_grid))
    map_est = np.empty(N)
    stated_var = np.empty(N)
    lam_hat = rng.normal(lam_star, np.sqrt(1.0 / tau_T), N)
    for i in range(N):
        mo_g = np.log(_sig(lam_hat[i]))
        mo_R = np.log(_sig(-lam_hat[i]))
        psi = -0.5 * p_g * (log_fg - mo_g) ** 2 - 0.5 * p_R * (log_fR - mo_R) ** 2
        post = np.exp(psi - psi.max())
        post /= post.sum()
        mean = post @ lam_grid
        map_est[i] = mean
        stated_var[i] = post @ (lam_grid**2) - mean**2
    actual_var = float(np.var(map_est))  # true scatter of the estimate = joint-truth var (~1/tau_T)
    stated = float(np.mean(stated_var))  # psi's OWN stated posterior variance
    ratio = stated / actual_var
    print(f" COMBINE  (f_g={f_g:.3f}, tau_T={tau_T})")
    print(f"      joint-truth Var(lam)      (= actual MAP scatter) = {actual_var:.6g}  (1/tau_T={1/tau_T:.6g})")
    print(f"      two-message STATED Var(lam) (psi curvature)      = {stated:.6g}")
    print(f"      RATIO two-message/joint-truth = {ratio:.4f}   (rank-1 double-count => expect ~0.5)")
    return ratio


def main():
    N = 400_000
    print("#" * 100)
    print("# Per-component message-variance model — independent MC validation (agent #1)\n")
    laws = []

    print("=" * 100)
    laws += graft(N, f_g=0.40, tau=1.0 / 0.02, n=800, n_s=600, E_g=110.0, E_r=200.0,
                  E_spl=100.0, M=900.0, S=1500.0, E_g_dst=380.0, E_r_dst=290.0, fgd=0.45, label="[A subst spliced]")
    print()
    laws += graft(N, f_g=0.40, tau=1.0 / 0.02, n=400, n_s=5000, E_g=110.0, E_r=200.0,
                  E_spl=100.0, M=200.0, S=20000.0, E_g_dst=380.0, E_r_dst=290.0, fgd=0.30, label="[C mature-dom]")
    print()
    # no-spliced degenerate: M cancels from the fraction, count term -> 0
    laws += graft(N, f_g=0.40, tau=1.0 / 0.02, n=800, n_s=10**9, E_g=110.0, E_r=200.0,
                  E_spl=100.0, M=900.0, S=1e-6, E_g_dst=380.0, E_r_dst=290.0, fgd=0.40, label="[B no-spliced]")

    print("=" * 100)
    laws += anchor(N, f_g=0.980, tau=1.0 / 0.02, n=500, E_g=110.0, E_r=200.0, M=1000.0)

    print("=" * 100)
    laws += transfer_var(N, f_g_s=0.30, f_g_d=0.55, tau_s=1.0 / 0.03, tau_d=1.0 / 0.02,
                         n_s_cnt=600, n_d_cnt=1200, E_g=110.0, E_r=200.0, M_s=800.0, M_d=1500.0)

    print("=" * 100)
    laws += peel(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.004,
                 rho_mu=0.10, n_s=1500, label="[D u~2]")
    print()
    laws += peel(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.002,
                 rho_mu=0.18, n_s=4000, label="[E u~10 tight]")

    print("=" * 100)
    r_sym = combine(2000, f_g=0.40, tau_T=50.0, E_g=100.0, E_r=100.0)
    print()
    r_asym = combine(2000, f_g=0.25, tau_T=80.0, E_g=100.0, E_r=100.0)

    print("=" * 100)
    print("\nSUMMARY (law, pred, emp, rel_err, pass):")
    for nm, pred, emp, rel, ok in laws:
        print(f"  {nm:<24} pred={pred:12.6g} emp={emp:12.6g} rel={rel:8.3%} {'PASS' if ok else 'FAIL'}")
    print(f"\n  combine ratio (symmetric)  = {r_sym:.4f}  (expect ~0.5)")
    print(f"  combine ratio (asymmetric) = {r_asym:.4f}")


if __name__ == "__main__":
    main()
