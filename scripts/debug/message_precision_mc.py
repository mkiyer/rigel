"""MESSAGE-PRECISION DERIVATION — Monte-Carlo validation of the GRAFT and PEEL variance laws.

Grounding (owner, 2026-07-23): **RNA is RNA.** There is one RNA species; "mature" and "nascent" are not two
entities to track but two ROUTES the same RNA takes at a junction — it either splices to another exon
(observed in the spliced channel) or continues contiguously into the next genomic region (observed unspliced).
Pass-0 solves the UNSPLICED pool; the spliced count is a direct measurement of the departing flux.

Notation at one exon|intron junction. Boundary ``b`` (a point, crossing flux), exon ``x`` (a region, contained):

    rho_mu(b) = S_b / E_spl        RNA that SPLICES OUT at b        MEASURED (motif-stranded, unambiguous RNA)
    rho_nu(b) = f_r(b)*M_b / E_r   RNA that CONTINUES through b     IMPUTED (needs the composition solve)
    rho_g(b)  = f_g(b)*M_b / E_g   gDNA crossing b                  IMPUTED
    r         = e(x)/e(b)          capture enrichment ratio, estimated as rho_tot(x)/rho_tot(b, WITH spliced)

THE FOUR CLAIMS THIS VALIDATES
------------------------------
T1  GRAFT is a MATCHED-SET transport and enrichment CANCELS. Including the spliced in the boundary's RNA
    makes its component set equal the exon's ({gDNA, RNA}), so the message is a pure transport of the
    content ratio ``k = rho_g/rho_R`` and ``r`` drops out of the final fraction identically.

T2  Var(log k) for the graft. k = rho_g / (rho_nu + rho_mu) with rho_g, rho_nu sharing the mass M_b:

        Var(log k) = w_mu^2 * (1/n_b + 1/n_s)  +  [1/f_g + w_nu/(1-f_g)]^2 * Var(f_g)
        w_mu = rho_mu/rho_R,  w_nu = rho_nu/rho_R,  w_mu + w_nu = 1

    Item E's share-weighting FALLS OUT of the delta method rather than being postulated. Degenerate check:
    with no spliced (w_mu = 0) this collapses to Var(logit f_g) EXACTLY, because k is then a ratio of two
    densities built from the same mass and the counting cancels.

T3/T4  PEEL is a DIFFERENCE and needs the enrichment scaling:

        rho_nu(b) = rho_R(x)/r - rho_mu(b)                                   (T3)
        Var(log rho_nu) = u^2*Var(log T) + (u-1)^2*Var(log rho_mu),  u = T/rho_nu,  T = rho_R(x)/r   (T4)
        Var(log T) = Var(log rho_R(x)) + Var(log r)

    The duality: a SUM carries convex weights in [0,1] (a minority component contributes quadratically
    LITTLE); a DIFFERENCE carries weights >= 1 (subtracting similar numbers DESTROYS precision). Same delta
    method, opposite regime. This is why ``Var(log r)`` (R2) is load-bearing on the PEEL and cancels on the
    GRAFT.

T5  The destination Jacobian: Var(log f_g(dst)) = (1 - f_g(DST))^2 * Var(log k) — evaluated at the
    DESTINATION's f_g, not the source's (the transport changes f_g because the eff-lengths differ).

    OMP_NUM_THREADS=1 python scripts/debug/message_precision_mc.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260723)


def _beta_draws(mean, var, size):
    """Beta samples with the requested mean/variance (the composition solve's posterior on f_g)."""
    m = float(mean)
    v = min(float(var), m * (1.0 - m) * 0.98)  # Beta admits var < m(1-m)
    c = m * (1.0 - m) / v - 1.0
    return rng.beta(m * c, (1.0 - m) * c, size=size)


def _lognormal(mean, var_log, size):
    """Positive draws with the requested Var(log ·) — the Poisson count in log space."""
    s = np.sqrt(max(var_log, 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def _report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(emp, 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<52} predicted {pred:12.6g}   empirical {emp:12.6g}   rel {rel:7.2%}")
    return rel < tol


def graft(n_draws, f_g, var_fg, n_b, n_s, E_g, E_r, E_spl, M_b, S_b, Egx, Erx, verbose=True):
    """T1/T2/T5 — the boundary -> exon GRAFT."""
    fg = _beta_draws(f_g, var_fg, n_draws)
    Mb = _lognormal(M_b, 1.0 / n_b, n_draws)
    Sb = _lognormal(S_b, 1.0 / n_s, n_draws)

    rho_g = fg * Mb / E_g
    rho_nu = (1.0 - fg) * Mb / E_r
    rho_mu = Sb / E_spl
    k = rho_g / (rho_nu + rho_mu)

    # --- predicted Var(log k) at the operating point ---
    rg, rn, rm = f_g * M_b / E_g, (1.0 - f_g) * M_b / E_r, S_b / E_spl
    rR = rn + rm
    w_mu, w_nu = rm / rR, rn / rR
    cnt = w_mu**2 * (1.0 / n_b + 1.0 / n_s)
    pred = cnt + (1.0 / f_g + w_nu / (1.0 - f_g)) ** 2 * var_fg          # the f_g-scale form
    # THE SHIPPING FORM: reparameterized onto the solver's own λ = logit(f_g). Algebraically identical, but
    # the singular 1/f_g and 1/(1-f_g) cancel against dλ/df_g, and it consumes `NodeDeconv.lam_var` natively.
    var_lam = var_fg / (f_g * (1.0 - f_g)) ** 2
    pred_logit = cnt + (1.0 - w_mu * f_g) ** 2 * var_lam
    emp = float(np.var(np.log(k)))
    ok2 = _report(f"T2  Var(log k) LOGIT form  [w_mu={w_mu:.3f}]", pred_logit, emp)
    assert abs(pred_logit - pred) <= 1e-9 * max(pred, 1e-300), "the two T2 forms must be identical"

    # --- T1: enrichment cancels. Scale EVERY source density by an arbitrary e; the message must not move. ---
    def msg_fg(rho_g_, rho_R_):
        kk = rho_g_ / rho_R_
        return kk * Egx / (kk * Egx + Erx)

    base = msg_fg(rg, rR)
    scaled = [msg_fg(rg * e, rR * e) for e in (0.001, 1.0, 1000.0, 1e6)]
    ok1 = max(abs(s - base) for s in scaled) < 1e-12
    if verbose:
        print(f"  {'OK ' if ok1 else '***'} T1  enrichment cancels on the matched-set graft"
              f"{'':<15} max |Δf_g| over e ∈ [1e-3, 1e6] = {max(abs(s - base) for s in scaled):.2e}")

    # --- T5: the destination Jacobian ---
    fgx = k * Egx / (k * Egx + Erx)
    pred5 = (1.0 - base) ** 2 * emp  # use the EMPIRICAL Var(log k) to isolate the Jacobian claim
    emp5 = float(np.var(np.log(fgx)))
    ok5 = _report("T5  Var(log f_g(dst)) = (1-f_g(DST))^2*Var(log k)", pred5, emp5)
    # the same thing with the SOURCE Jacobian, to show it is the wrong one
    wrong5 = (1.0 - f_g) ** 2 * emp
    print(f"      (source-Jacobian form would predict {wrong5:.6g} — "
          f"off by {abs(wrong5 - emp5) / emp5:.1%}; f_g src={f_g:.3f} dst={base:.3f})")
    return ok1 and ok2 and ok5


def peel(n_draws, rho_R_x, var_log_rhoR, r, var_log_r, rho_mu, n_s, verbose=True):
    """T3/T4 — the exon -> boundary PEEL (a difference, enrichment-scaled)."""
    Rx = _lognormal(rho_R_x, var_log_rhoR, n_draws)
    rr = _lognormal(r, var_log_r, n_draws)
    rm = _lognormal(rho_mu, 1.0 / n_s, n_draws)
    T = Rx / rr
    nu = T - rm
    keep = nu > 0
    frac = keep.mean()

    T0 = rho_R_x / r
    nu0 = T0 - rho_mu
    u = T0 / nu0
    var_logT = var_log_rhoR + var_log_r
    pred = u**2 * var_logT + (u - 1.0) ** 2 * (1.0 / n_s)
    emp = float(np.var(np.log(nu[keep])))
    ok = _report(f"T4  Var(log rho_nu)  [u={u:.2f}, kept {frac:.1%}]", pred, emp, tol=0.12)
    return ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# Monte-Carlo validation of the message-precision laws   draws={N:,}\n")

    print("GRAFT (boundary -> exon) — matched set {gDNA, RNA}, enrichment cancels")
    print("\n [A] junction with substantial spliced flux (w_mu ~ 0.5)")
    graft(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=600, E_g=110.0, E_r=200.0,
          E_spl=100.0, M_b=900.0, S_b=1500.0, Egx=380.0, Erx=290.0)

    print("\n [B] DEGENERATE CHECK — no spliced (w_mu = 0): must collapse to Var(logit f_g) exactly")
    fgm, vfg = 0.40, 0.004
    lam_var = vfg / (fgm * (1.0 - fgm)) ** 2
    print(f"      Var(logit f_g) = Var(f_g)/[f_g(1-f_g)]^2 = {lam_var:.6g}")
    graft(N, f_g=fgm, var_fg=vfg, n_b=800, n_s=10**9, E_g=110.0, E_r=200.0,
          E_spl=100.0, M_b=900.0, S_b=1e-9, Egx=380.0, Erx=290.0)

    print("\n [C] mature-DOMINATED junction (w_mu -> 1): counting no longer cancels")
    graft(N, f_g=0.40, var_fg=0.004, n_b=400, n_s=5000, E_g=110.0, E_r=200.0,
          E_spl=100.0, M_b=200.0, S_b=20000.0, Egx=380.0, Erx=290.0)

    print("\n\nPEEL (exon -> boundary) — a DIFFERENCE, enrichment-scaled; Var(log r) is load-bearing")
    print("\n [D] comfortable residual (u ~ 2)")
    peel(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.004,
         rho_mu=0.10, n_s=1500)
    print("\n [E] TIGHT residual (u ~ 10) — the difference destroys precision")
    peel(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.002,
         rho_mu=0.18, n_s=4000)
    print("\n [F] Var(log r) ablated to 0 — shows how much of the peel's variance it carries")
    peel(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.0,
         rho_mu=0.10, n_s=1500)


if __name__ == "__main__":
    main()
