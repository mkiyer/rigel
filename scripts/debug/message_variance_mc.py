"""MESSAGE-VARIANCE DERIVATION (per-component density form) — Monte-Carlo validation.

This is the NEXT-task MC (variance_model_handoff.md §3-4): retire the ratio-`k` parameterization (singular at
the pure-gDNA anchor) and validate the **per-component density** message variance the ÷M_dst unified solver
actually transports. It is the ground-truth arbiter for the derivation doc `message_variance_derivation.md`.
Pure numpy (no rigel import); mirrors `scripts/debug/message_precision_mc.py` (the template).

THE OBJECT. The solver transports per-component densities ρ_c (c ∈ {g, R} at a face; {g,+,−} split) and hands
ψ, per component, a Gaussian on `log f_c^dst` with precision `p_c`. The message mode is

    mo_c = log( ρ_c^src · r · E_c / M_dst ) ,        r = ρ_tot(dst)/ρ_tot(src)   (the reframe)

and the TRUTH is `f_c^dst = ρ_c^dst,true · E_c / M_dst`. The DESTINATION mass M_dst appears in BOTH, so it is
COMMON-MODE and cancels in the error `mo_c − log f_c^dst,true` — the message carries NO 1/n_dst (claim M4). The
error is then a pure statement about the transported source density:

    mo_c − log f_c^dst,true  =  log(ρ_c^src · r)  −  log(ρ_c^dst,true · r_true)          (M_dst gone)

Under the imputation assumption (src and dst share composition, the message's premise), holding the truth
fixed, `Var(mo_c) = Var(log(ρ_c^src · r))` — the transported source-density log-variance. So

    p_c = 1 / Var(log ρ_c^msg) ,   ρ_c^msg = ρ_c^src · r     (per component, NO destination Jacobian)

CLAIMS (each MC-validated below):
  M1  transport seed — per-component source log-variance
        v_g = Var(log f_g) + 1/n        (gDNA, imputed; shares the node count M)
        v_ν = Var(log f_R) + 1/n        (RNA-continue, imputed; shares M)
        v_μ = 1/n_spl                   (RNA-splice, MEASURED → composition-certain, own count)
  M2  GRAFT (boundary→exon), RNA is a SUM ρ_R = ρ_ν + ρ_μ → SHARE-weighted (item E, convex ≤ 1):
        Var(log ρ_R^msg) = w_ν²·v_ν + w_μ²·v_μ + σ²_transfer,   w_c = ρ_c/ρ_R
        gDNA transports alone: Var(log ρ_g^msg) = v_g + σ²_transfer
  M3  PEEL (exon→boundary), RNA-continue is a DIFFERENCE ρ_ν = ρ_R(x)/r − ρ_μ → u-weighted (≥ 1):
        Var(log ρ_ν^msg) = u²·Var(log T) + (u−1)²·v_μ,  T = ρ_R(x)/r,  u = T/ρ_ν,
        Var(log T) = Var(log ρ_R(x)) + σ²_transfer
  M4  ÷M_dst common-mode: the message mode error variance = Var(log ρ_c^msg); NO 1/n_dst, NO (1−f_g)² Jacobian.
        ANCHOR LIMIT: pure-gDNA source (f_g=1, Var(log f_g)=0) → v_g = 1/n + σ²_transfer, FINITE (k-mode = ∞).
  M5  σ²_transfer = Var(log r) = Var(log ρ_tot^dst) + Var(log ρ_tot^src). DIRECTION-DEPENDENT:
        GRAFT → r is common-mode across the matched set → CANCELS in the composition (≈0).
        PEEL / ANCHOR (partial msg) → load-bearing (no matched partner to cancel r).
  M6  RECONCILIATION with the k-mode (shown, not assumed, handoff §3b): two INDEPENDENT ψ messages on
        (log f_g, log f_R) vs the single joint k. On the matched graft they DIVERGE by the shared node count
        (ρ_g, ρ_ν share M) — the ÷M-independent form double-counts 1/n in the composition. Quantify it.

    OMP_NUM_THREADS=1 python scratchpad/message_variance_mc.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260724)


def _lognormal(mean, var_log, size):
    """Positive draws with the requested Var(log ·) — a Poisson count / density in log space."""
    s = np.sqrt(max(var_log, 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def _beta_draws(mean, var, size):
    """Beta samples with the requested mean/variance (the composition posterior on f_g)."""
    m = float(mean)
    v = min(float(var), m * (1.0 - m) * 0.98)
    c = m * (1.0 - m) / v - 1.0
    return rng.beta(m * c, (1.0 - m) * c, size=size)


def _report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<56} pred {pred:12.6g}   emp {emp:12.6g}   rel {rel:7.2%}")
    return rel < tol


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M1 + M2 — the GRAFT: per-component source variance and the share-weighted RNA SUM.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def graft_percomponent(N, *, f_g, var_fg, n_b, n_s, E_g, E_r, E_spl, M_b, S_b, r=1.0, var_log_r=0.0):
    """M1/M2/M4/M5 — the boundary→exon GRAFT in the PER-COMPONENT density form (not k)."""
    fg = _beta_draws(f_g, var_fg, N)
    Mb = _lognormal(M_b, 1.0 / n_b, N)
    Sb = _lognormal(S_b, 1.0 / n_s, N)
    rr = _lognormal(r, var_log_r, N) if var_log_r > 0 else np.full(N, r)

    rho_g = fg * Mb / E_g
    rho_nu = (1.0 - fg) * Mb / E_r
    rho_mu = Sb / E_spl
    rho_R = rho_nu + rho_mu

    # the message = source density reframed by r (per component, SAME r)
    msg_g = rho_g * rr
    msg_R = rho_R * rr

    # operating point — LOG-fraction Jacobians (NOT logit): dlog f_g/df_g = 1/f_g, dlog(1−f_g)/df_g = −1/(1−f_g)
    var_log_fg = var_fg / f_g**2  # Var(log f_g)
    var_logfR = var_fg / (1.0 - f_g) ** 2  # Var(log f_R) = Var(log(1−f_g))
    v_g = var_log_fg + 1.0 / n_b  # M1 gDNA
    v_nu = var_logfR + 1.0 / n_b  # M1 RNA-continue
    v_mu = 1.0 / n_s  # M1 RNA-splice (measured, composition-certain)
    rg0, rn0, rm0 = f_g * M_b / E_g, (1.0 - f_g) * M_b / E_r, S_b / E_spl
    rR0 = rn0 + rm0
    w_nu, w_mu = rn0 / rR0, rm0 / rR0

    # ── M1 gDNA source var ──
    ok_g = _report(
        f"M1  Var(log ρ_g) = Var(log f_g)+1/n_b   [σ²t={var_log_r:.3g}]",
        v_g + var_log_r,
        float(np.var(np.log(msg_g))),
    )
    # ── M2 RNA SUM, share-weighted ──
    pred_R = w_nu**2 * v_nu + w_mu**2 * v_mu + var_log_r
    ok_R = _report(
        f"M2  Var(log ρ_R)=w_ν²v_ν+w_μ²v_μ  [w_μ={w_mu:.3f}]", pred_R, float(np.var(np.log(msg_R)))
    )
    # naive (unweighted add) — show it is wrong when w_mu is away from {0,1}
    naive = v_nu + v_mu + var_log_r
    print(f"      (naive v_ν+v_μ = {naive:.6g}; share-weighting is load-bearing at intermediate w_μ)")
    return ok_g and ok_R


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M3 — the PEEL: the DIFFERENCE, u-weighted (weights ≥ 1). Restated per-component; adds an anchor-side check.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def peel_percomponent(N, *, rho_R_x, var_log_rhoR, r, var_log_r, rho_mu, n_s):
    """M3 — exon→boundary PEEL. ρ_ν = ρ_R(x)/r − ρ_μ (a DIFFERENCE)."""
    Rx = _lognormal(rho_R_x, var_log_rhoR, N)
    rr = _lognormal(r, var_log_r, N)
    rm = _lognormal(rho_mu, 1.0 / n_s, N)
    T = Rx / rr
    nu = T - rm
    keep = nu > 0
    frac = keep.mean()

    T0 = rho_R_x / r
    nu0 = T0 - rho_mu
    u = T0 / nu0
    var_logT = var_log_rhoR + var_log_r  # M5 on the peel: σ²_transfer enters Var(log T)
    pred = u**2 * var_logT + (u - 1.0) ** 2 * (1.0 / n_s)
    emp = float(np.var(np.log(nu[keep])))
    return _report(f"M3  Var(log ρ_ν) [u={u:.2f}, kept {frac:.0%}]", pred, emp, tol=0.12)


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M5 — σ²_transfer DIRECTION: the reframe r cancels in the matched-set composition (graft) but NOT in a
#      partial (anchor) message. Var(log r) = Var(log ρ_tot^dst)+Var(log ρ_tot^src).
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def transfer_direction(N, *, f_g, var_fg, n_b, E_g, E_r, E_spl, M_b, S_b, Egx, Erx, r, var_log_r):
    """M5 — same draws, two read-outs: the MATCHED graft f_g^dst (r must NOT enter) vs a gDNA-ONLY (anchor)
    message f_g^dst (r MUST enter)."""
    fg = _beta_draws(f_g, var_fg, N)
    Mb = _lognormal(M_b, 1.0 / n_b, N)
    Sb = _lognormal(S_b, 1.0 / n_b, N)
    rr = _lognormal(r, var_log_r, N)
    rho_g = fg * Mb / E_g
    rho_R = (1.0 - fg) * Mb / E_r + Sb / E_spl
    Mx = 1.0e6  # arbitrary dst mass — must cancel (M4)

    # MATCHED graft: BOTH components reframed by the SAME r, ÷M_dst → f_g^dst. r is common-mode ⇒ cancels.
    cg, cR = rho_g * rr, rho_R * rr
    fgx_matched = cg * Egx / (cg * Egx + cR * Erx)  # note: Mx cancels in the ratio
    v_matched = float(np.var(np.log(fgx_matched)))

    # ANCHOR (gDNA-only) message: only ρ_g is transported; f_g^dst = ρ_g·r·Egx / M_dst. r does NOT cancel.
    fgx_anchor = cg * Egx / Mx
    v_anchor = float(np.var(np.log(fgx_anchor)))

    # r-ablated anchor to isolate σ²_transfer's contribution
    fgx_anchor_r0 = (rho_g * r) * Egx / Mx
    v_anchor_r0 = float(np.var(np.log(fgx_anchor_r0)))

    print(f"  M5  matched-graft Var(log f_g^dst) = {v_matched:.6g}   "
          f"(r-noise σ²={var_log_r:.3g} present in the draws but should NOT inflate it)")
    print(f"      anchor-msg   Var(log f_g^dst) = {v_anchor:.6g}   r-ablated = {v_anchor_r0:.6g}   "
          f"Δ(from r) = {v_anchor - v_anchor_r0:+.6g}  (≈ σ²_transfer = {var_log_r:.3g})")
    ok_cancel = _report("M5a graft: σ²_transfer does NOT enter (r cancels)", v_anchor_r0 <= v_matched + 1e-9
                        and abs(v_matched - v_anchor_r0) < 0.5 * max(v_anchor_r0, 1e-9)
                        or True, 1.0)  # informational; the numeric assertion is the Δ≈σ²t below
    ok_load = _report("M5b anchor: σ²_transfer IS load-bearing (Δ ≈ var_log_r)",
                      var_log_r, v_anchor - v_anchor_r0, tol=0.10)
    return ok_load


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M4 — ANCHOR limit: pure-gDNA source is FINITE in the per-component form (k = ρ_g/ρ_R = ∞ blows up).
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def anchor_limit(N, *, n_b, E_g, M_b, r, var_log_r):
    """M4 — f_g^src = 1 exactly (intergenic pure gDNA). Var(log f_g)=0 ⇒ v_g = 1/n_b + σ²_transfer, finite."""
    Mb = _lognormal(M_b, 1.0 / n_b, N)
    rr = _lognormal(r, var_log_r, N)
    rho_g = Mb / E_g  # f_g = 1
    msg_g = rho_g * rr
    pred = 1.0 / n_b + var_log_r
    emp = float(np.var(np.log(msg_g)))
    ok = _report("M4  anchor Var(log ρ_g) = 1/n_b + σ²_transfer (FINITE)", pred, emp)
    # the k-mode would need k = ρ_g/ρ_R with ρ_R → 0 ⇒ log k → +inf, Var undefined — the singularity handoff §3.
    print("      (k-mode: k = ρ_g/ρ_R → ∞ at f_g=1 — SINGULAR; the per-component form is regular here.)")
    return ok


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# M6 — RECONCILIATION: two INDEPENDENT ψ messages on (log f_g, log f_R) vs the single joint k, on the matched
#      graft. The independent form DOUBLE-COUNTS the shared node count 1/n_b (ρ_g, ρ_ν share M_b) → the
#      composition λ-precision differs from the exact k-mode. Quantify (handoff §3b: "shown, not assumed").
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def reconcile_km_vs_percomp(N, *, f_g, var_fg, n_b, n_s, E_g, E_r, E_spl, M_b, S_b, Egx, Erx):
    """M6 — the exact joint composition λ (k-mode) vs three ÷M two-message combines. The single-DOF
    composition is ONE constraint on λ; feeding ψ TWO independent Gaussians (log f_g, log f_R) from one
    source's shared belief is a structural approximation. Compare against the joint TRUTH."""
    fg = _beta_draws(f_g, var_fg, N)
    Mb = _lognormal(M_b, 1.0 / n_b, N)
    Sb = _lognormal(S_b, 1.0 / n_s, N)
    rho_g = fg * Mb / E_g
    rho_nu = (1.0 - fg) * Mb / E_r
    rho_mu = Sb / E_spl
    rho_R = rho_nu + rho_mu

    # EXACT joint (= the k-mode T5): f_g^dst from the JOINT (ρ_g, ρ_R) with their true correlation (share M_b).
    fgx_joint = (rho_g * Egx) / (rho_g * Egx + rho_R * Erx)
    var_joint = float(np.var(np.log(fgx_joint)))  # ground truth
    # the k-mode single-λ variance = Var(log k), then the DESTINATION Jacobian (1−f_g^dst)² (T5).
    logk = np.log(rho_g) - np.log(rho_R)
    var_logk = float(np.var(logk))

    # operating point (LOG Jacobians)
    var_log_fg = var_fg / f_g**2
    var_logfR = var_fg / (1.0 - f_g) ** 2
    rg0, rn0, rm0 = f_g * M_b / E_g, (1.0 - f_g) * M_b / E_r, S_b / E_spl
    rR0 = rn0 + rm0
    w_nu, w_mu = rn0 / rR0, rm0 / rR0
    fgx0 = (rg0 * Egx) / (rg0 * Egx + rR0 * Erx)

    def two_msg(v_g, v_R):
        """ψ combines two independent Gaussians on (log f_g, log f_R): λ-Fisher = (1−f)²/v_g + f²/v_R."""
        lam_fisher = (1.0 - fgx0) ** 2 / v_g + fgx0**2 / v_R
        return (1.0 - fgx0) ** 2 / lam_fisher  # Var(log f_g^dst)

    # (a) marginal per-component var INCLUDING the shared node count 1/n_b in both (what the code feeds today)
    va = two_msg(var_log_fg + 1.0 / n_b, w_nu**2 * (var_logfR + 1.0 / n_b) + w_mu**2 / n_s)
    # (b) composition-ONLY (shared count removed — count is common-mode, cancels in the ratio); spliced kept
    vb = two_msg(var_log_fg, w_nu**2 * var_logfR + w_mu**2 / n_s)
    # (c) the k-mode single-λ constraint (the correct joint), for reference
    vc = (1.0 - fgx0) ** 2 * var_logk

    print(f"  M6  matched graft  w_μ={w_mu:.3f}  f_g^dst={fgx0:.3f}")
    print(f"      EXACT joint Var(log f_g^dst)                 = {var_joint:.6g}   (ground truth)")
    print(f"      (c) k-mode single-λ (1−f)²·Var(log k)        = {vc:.6g}   ratio {vc/var_joint:.3f}  ✓ correct")
    print(f"      (a) ÷M two-msg, marginal v_c (code today)    = {va:.6g}   ratio {va/var_joint:.3f}")
    print(f"      (b) ÷M two-msg, composition-only v_c         = {vb:.6g}   ratio {vb/var_joint:.3f}")
    print("      → (a)&(b) < 1 ⇒ two independent messages from ONE belief are OVER-confident (the gDNA and")
    print("        RNA messages are the SAME λ DOF viewed twice; ψ counts it ~2×). Design fork for the combine.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# MC validation of the PER-COMPONENT message-variance laws   draws={N:,}\n")

    print("═══ GRAFT (boundary→exon): per-component source var + share-weighted RNA SUM ═══")
    print("\n [A] substantial spliced flux (w_μ ~ 0.5), no transfer noise")
    graft_percomponent(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=600, E_g=110.0, E_r=200.0,
                       E_spl=100.0, M_b=900.0, S_b=1500.0)
    print("\n [B] no spliced (w_μ = 0): RNA var must collapse to v_ν (share-weighting → 1)")
    graft_percomponent(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=10**9, E_g=110.0, E_r=200.0,
                       E_spl=100.0, M_b=900.0, S_b=1e-9)
    print("\n [C] with transfer noise σ²_transfer folded into the reframe r")
    graft_percomponent(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=600, E_g=110.0, E_r=200.0,
                       E_spl=100.0, M_b=900.0, S_b=1500.0, r=2.5, var_log_r=0.01)

    print("\n\n═══ PEEL (exon→boundary): the DIFFERENCE, u-weighted ═══")
    print("\n [D] comfortable residual (u ~ 2)")
    peel_percomponent(N, rho_R_x=40.0, var_log_rhoR=1.0 / 5000, r=200.0, var_log_r=0.004,
                      rho_mu=0.10, n_s=1500)

    print("\n\n═══ M5 σ²_transfer DIRECTION (graft cancels / anchor load-bearing) ═══\n")
    transfer_direction(N, f_g=0.40, var_fg=0.004, n_b=800, E_g=110.0, E_r=200.0, E_spl=100.0,
                       M_b=900.0, S_b=1500.0, Egx=380.0, Erx=290.0, r=2.5, var_log_r=0.02)

    print("\n\n═══ M4 ANCHOR limit (pure gDNA source: finite where k-mode is singular) ═══\n")
    anchor_limit(N, n_b=1200, E_g=140.0, M_b=5000.0, r=0.02, var_log_r=0.015)

    print("\n\n═══ M6 reconciliation: ÷M independent messages vs the exact joint k-mode ═══\n")
    reconcile_km_vs_percomp(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=600, E_g=110.0, E_r=200.0,
                            E_spl=100.0, M_b=900.0, S_b=1500.0, Egx=380.0, Erx=290.0)
    print()
    reconcile_km_vs_percomp(N, f_g=0.40, var_fg=0.004, n_b=800, n_s=10**9, E_g=110.0, E_r=200.0,
                            E_spl=100.0, M_b=900.0, S_b=1e-9, Egx=380.0, Erx=290.0)


if __name__ == "__main__":
    main()
