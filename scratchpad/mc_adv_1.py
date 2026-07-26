"""ADVERSARIAL verification of the message-variance laws M1-M6 (verifier #1).

Goal: REFUTE. For each law find a regime where it is WRONG or silently OVER-confident (linearization
under-states variance). Every attack is MC'd. Then the M6 combine fork is tested with a FAITHFUL replica
of what ψ actually does (two Gaussians on log f_g and log(1-f_g)), across regimes, plus the two proposed
resolutions (single-λ via Var(log k); the anchor where single-λ-via-k is singular).

    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel \
        && OMP_NUM_THREADS=1 python /Users/mkiyer/proj/rigel/scratchpad/mc_adv_1.py
"""

from __future__ import annotations

import numpy as np

rng = np.random.default_rng(20260725)


def _ln(mean, var_log, size):
    s = np.sqrt(max(var_log, 0.0))
    return float(mean) * np.exp(rng.normal(-0.5 * s * s, s, size=size))


def _beta(mean, var, size):
    m = float(mean)
    v = min(float(var), m * (1.0 - m) * 0.999)
    c = m * (1.0 - m) / v - 1.0
    return rng.beta(m * c, (1.0 - m) * c, size=size)


def rep(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<58} pred {pred:12.6g}  emp {emp:12.6g}  rel {rel:7.2%}  "
          f"{'OVERCONF' if pred < emp and rel > tol else ''}")
    return rel < tol


N = 800_000
print(f"# ADVERSARIAL MC  draws={N:,}\n")

# ══════════════════════════════════════════════════════════════════════════════════════════════
# ATTACK M1 — the delta-method Var(log f_g)=var_fg/f_g^2 UNDER-states near the vertices / at large var.
# The transport seed's composition term is a linearization; if it under-states, every downstream law
# inherits over-confidence. Push f_g→0, f_g→1, and large var_fg.
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("═══ ATTACK M1 — delta-method Var(log f_g) under-statement (over-confidence) ═══")
for f0, v0, n in [(0.40, 0.004, 800), (0.08, 0.004, 800), (0.04, 0.010, 400),
                  (0.95, 0.002, 800), (0.015, 0.00020, 2000)]:
    fg = _beta(f0, v0, N)
    Mb = _ln(900.0, 1.0 / n, N)
    rho_g = fg * Mb / 110.0
    pred = v0 / f0**2 + 1.0 / n
    emp = float(np.var(np.log(rho_g)))
    rep(f"M1 v_g  f_g={f0:.3f} var={v0:.4f} n={n}", pred, emp)
    # RNA arm at the same operating point
    rho_r = (1.0 - fg) * Mb / 200.0
    predR = v0 / (1.0 - f0) ** 2 + 1.0 / n
    empR = float(np.var(np.log(rho_r)))
    rep(f"M1 v_nu f_g={f0:.3f}", predR, empR)

# ══════════════════════════════════════════════════════════════════════════════════════════════
# ATTACK M2 — the graft SUM share-weighting at intermediate w_mu with a THIN spliced count.
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ ATTACK M2 — graft SUM share-weight, intermediate w_mu, thin n_s ═══")
for f0, ns, Sb in [(0.40, 30, 300.0), (0.40, 8, 120.0), (0.60, 50, 400.0)]:
    nb = 800
    fg = _beta(f0, 0.004, N)
    Mb = _ln(900.0, 1.0 / nb, N)
    Sb_d = _ln(Sb, 1.0 / ns, N)
    rho_nu = (1.0 - fg) * Mb / 200.0
    rho_mu = Sb_d / 100.0
    rho_R = rho_nu + rho_mu
    rn0 = (1.0 - f0) * 900.0 / 200.0
    rm0 = Sb / 100.0
    rR0 = rn0 + rm0
    w_nu, w_mu = rn0 / rR0, rm0 / rR0
    v_nu = 0.004 / (1.0 - f0) ** 2 + 1.0 / nb
    v_mu = 1.0 / ns
    pred = w_nu**2 * v_nu + w_mu**2 * v_mu
    emp = float(np.var(np.log(rho_R)))
    rep(f"M2 Var(log rho_R)  w_mu={w_mu:.3f} n_s={ns}", pred, emp)

# ══════════════════════════════════════════════════════════════════════════════════════════════
# ATTACK M3 — the PEEL over a u-sweep. This is the KNOWN breakdown; map the over-confidence and
# check where it crosses the 8% line. Real junctions: median u=2.3, p90=9.0.
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ ATTACK M3 — peel DIFFERENCE over a u-sweep (over-confidence map) ═══")
rho_R_x, var_log_rhoR, r, var_log_r, n_s = 40.0, 1.0 / 5000, 200.0, 0.004, 1500
for rho_mu in [0.05, 0.10, 0.14, 0.16, 0.18, 0.19]:
    Rx = _ln(rho_R_x, var_log_rhoR, N)
    rr = _ln(r, var_log_r, N)
    rm = _ln(rho_mu, 1.0 / n_s, N)
    nu = Rx / rr - rm
    keep = nu > 0
    T0 = rho_R_x / r
    nu0 = T0 - rho_mu
    u = T0 / nu0
    var_logT = var_log_rhoR + var_log_r
    pred = u**2 * var_logT + (u - 1.0) ** 2 * (1.0 / n_s)
    emp = float(np.var(np.log(nu[keep])))
    rep(f"M3 peel  u={u:5.2f}  kept={keep.mean():.1%}", pred, emp, tol=0.08)

# ══════════════════════════════════════════════════════════════════════════════════════════════
# ATTACK M4 — anchor near (not at) the vertex: f_g=0.99, 0.999. Is the finite per-component form
# still faithful when the Beta sits a hair off the wall?
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ ATTACK M4 — anchor near the f_g=1 wall ═══")
for f0, v0, n in [(1.0, 0.0, 1200), (0.99, 1e-4, 1200), (0.995, 2e-5, 3000)]:
    if v0 > 0:
        fg = _beta(f0, v0, N)
    else:
        fg = np.ones(N)
    Mb = _ln(5000.0, 1.0 / n, N)
    rr = _ln(0.02, 0.015, N)
    msg_g = fg * Mb / 140.0 * rr
    pred = (v0 / f0**2 if v0 > 0 else 0.0) + 1.0 / n + 0.015
    emp = float(np.var(np.log(msg_g)))
    rep(f"M4 anchor  f_g={f0:.3f} var={v0:.1e}", pred, emp)

# ══════════════════════════════════════════════════════════════════════════════════════════════
# ATTACK M5 — (a) is the graft cancellation EXACT, or does a component-set MISMATCH re-introduce r?
#            (b) does Var(log r)=Var(dst)+Var(src) OVER-state when src/dst SHARE fragments (correlation)?
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ ATTACK M5 — graft cancellation & the additive Var(log r) independence assumption ═══")
# (a) matched set: r must cancel. Introduce r-noise; matched-graft f_g must not move.
f0, v0, nb, ns = 0.40, 0.004, 800, 600
fg = _beta(f0, v0, N)
Mb = _ln(900.0, 1.0 / nb, N)
Sb = _ln(1500.0, 1.0 / ns, N)
rr = _ln(2.5, 0.05, N)   # BIG transfer noise
rho_g = fg * Mb / 110.0
rho_R = (1.0 - fg) * Mb / 200.0 + Sb / 100.0
Egx, Erx = 380.0, 290.0
# matched: both reframed by SAME rr
fgx_m = (rho_g * rr * Egx) / (rho_g * rr * Egx + rho_R * rr * Erx)
v_m = float(np.var(np.log(fgx_m)))
# r-ablated matched
fgx_m0 = (rho_g * 2.5 * Egx) / (rho_g * 2.5 * Egx + rho_R * 2.5 * Erx)
v_m0 = float(np.var(np.log(fgx_m0)))
print(f"  M5a matched-graft Var(log f_g): with r-noise σ²=0.05 = {v_m:.6g}  r-ablated = {v_m0:.6g}  "
      f"Δ={v_m - v_m0:+.2e}  (should be ~0 ⇒ r cancels)")
rep("M5a graft: σ²_transfer does NOT enter", v_m0, v_m, tol=0.03)
# (b) shared-fragment correlation: build dst and src totals that SHARE a common fragment pool.
# If they share, Var(log r) < Var(dst)+Var(src). Quantify the over-statement (safe direction).
common = _ln(1.0, 1.0 / 400, N)      # shared pool (both nodes see it)
priv_d = _ln(1.0, 1.0 / 400, N)
priv_s = _ln(1.0, 1.0 / 400, N)
rho_tot_d = 0.5 * common + 0.5 * priv_d   # 50% shared
rho_tot_s = 0.5 * common + 0.5 * priv_s
r_shared = rho_tot_d / rho_tot_s
emp_shared = float(np.var(np.log(r_shared)))
pred_indep = float(np.var(np.log(rho_tot_d))) + float(np.var(np.log(rho_tot_s)))
print(f"  M5b Var(log r): additive(indep) pred={pred_indep:.6g}  emp(50%-shared)={emp_shared:.6g}  "
      f"ratio={pred_indep / emp_shared:.2f}  ⇒ additive {'OVER' if pred_indep > emp_shared else 'UNDER'}-states"
      f" (OVER=under-confident=SAFE)")

# ══════════════════════════════════════════════════════════════════════════════════════════════
# THE M6 COMBINE FORK — FAITHFUL replica of ψ. ψ receives from ONE source belief two Gaussians:
#   gDNA on log f_g  (prec p_g),  RNA on log(1-f_g)=log f_R (prec p_R).
# The λ-precision ψ forms is  p_g·(1-f)^2 + p_R·f^2  (at the mode) ⇒ Var(log f_g^dst)=(1-f)^2/that.
# Compare to the JOINT TRUTH (transport the correlated source belief), across regimes, and test the
# fixes. Key adversarial question: is single-λ-via-Var(log k) enough (does it keep the spliced?), and
# does it break at the anchor (only gDNA live)?
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ M6 COMBINE FORK — faithful ψ two-Gaussian combine vs joint truth ═══")


def m6(tag, *, f_g, var_fg, n_b, n_s, E_g, E_r, E_spl, M_b, S_b, Egx, Erx,
       rna_indep_only=False):
    """rna_indep_only: the RNA arm carries ONLY the independent spliced measurement (NO shared τ) —
    tests whether summing genuinely-independent arms is correct (ratio→1, no double-count)."""
    fg = _beta(f_g, var_fg, N)
    Mb = _ln(M_b, 1.0 / n_b, N)
    Sb = _ln(S_b, 1.0 / n_s, N)
    rho_g = fg * Mb / E_g
    rho_nu = (1.0 - fg) * Mb / E_r
    rho_mu = Sb / E_spl
    rho_R = rho_nu + rho_mu
    # JOINT TRUTH — transport the correlated (ρ_g, ρ_R) to the destination fraction.
    fgx_joint = (rho_g * Egx) / (rho_g * Egx + rho_R * Erx)
    var_joint = float(np.var(np.log(fgx_joint)))
    logk = np.log(rho_g) - np.log(rho_R)
    var_logk = float(np.var(logk))

    # operating point
    rg0, rn0, rm0 = f_g * M_b / E_g, (1.0 - f_g) * M_b / E_r, S_b / E_spl
    rR0 = rn0 + rm0
    w_nu, w_mu = rn0 / rR0, rm0 / rR0
    fgx0 = (rg0 * Egx) / (rg0 * Egx + rR0 * Erx)
    var_log_fg = var_fg / f_g**2
    var_logfR = var_fg / (1.0 - f_g) ** 2

    # (a) CODE TODAY: two-message, composition τ in BOTH arms + shared count in both.
    v_g_a = var_log_fg + 1.0 / n_b
    v_R_a = w_nu**2 * (var_logfR + 1.0 / n_b) + w_mu**2 / n_s
    lam_fish_a = (1.0 - fgx0) ** 2 / v_g_a + fgx0**2 / v_R_a
    va = (1.0 - fgx0) ** 2 / lam_fish_a

    # (c) single-λ via Var(log k) — the joint constraint, one message, spliced INCLUDED via w_mu.
    vc = (1.0 - fgx0) ** 2 * var_logk

    # (d) the CONSENSUS FIX faithful: gDNA arm carries the FULL composition λ (via Var(log k)); RNA
    #     arm carries ONLY the independent extra. Implemented as single-λ = (c). (kept for clarity.)
    if rna_indep_only:
        # RNA arm carries ONLY spliced measurement precision (independent), gDNA arm ONLY its own
        # composition — genuinely independent ⇒ summing is CORRECT.
        v_g_e = var_log_fg
        v_R_e = w_mu**2 / n_s if w_mu > 0 else 1e18
        lam_fish_e = (1.0 - fgx0) ** 2 / v_g_e + fgx0**2 / v_R_e
        ve = (1.0 - fgx0) ** 2 / lam_fish_e
    else:
        ve = np.nan

    print(f"  [{tag}]  w_mu={w_mu:.3f}  f_g^src={f_g:.3f}  f_g^dst={fgx0:.3f}")
    print(f"      joint TRUTH Var(log f_g^dst)     = {var_joint:.6g}")
    print(f"      (a) code-today two-msg           = {va:.6g}   ratio {va / var_joint:.3f}  "
          f"{'OVERCONF' if va < var_joint else 'ok'}")
    print(f"      (c) single-λ via Var(log k)      = {vc:.6g}   ratio {vc / var_joint:.3f}")
    if rna_indep_only:
        print(f"      (e) indep-arms (no shared τ)     = {ve:.6g}   ratio {ve / var_joint:.3f}")
    return va / var_joint, vc / var_joint


# regime sweep: vary w_mu (spliced share) — from pure composition (no spliced) to mature-dominated.
r_a, r_c = [], []
for tag, ns, Sb in [("no-spliced", 10**12, 1e-12), ("light-spl", 4000, 300.0),
                    ("mid-spl", 600, 1500.0), ("heavy-spl", 200, 8000.0)]:
    ra, rc = m6(tag, f_g=0.40, var_fg=0.004, n_b=800, n_s=ns, E_g=110.0, E_r=200.0,
                E_spl=100.0, M_b=900.0, S_b=Sb, Egx=380.0, Erx=290.0)
    r_a.append(ra)
    r_c.append(rc)
print(f"\n  → code-today ratio range across w_mu: [{min(r_a):.3f}, {max(r_a):.3f}]  "
      f"(over-conf factor {1 / max(r_a):.2f}×–{1 / min(r_a):.2f}×)")
print(f"  → single-λ (Var log k) ratio range:    [{min(r_c):.3f}, {max(r_c):.3f}]  "
      f"(≈1 ⇒ correct, spliced KEPT)")

# symmetric-frame control: source f_g == dst f_g ⇒ pure factor should be exactly 0.5 (2×).
print("\n  [symmetric-frame control: Egx/Erx tuned so f_g^dst≈f_g^src ⇒ isolate the pure factor]")
# choose Egx=E_g, Erx=E_r so dst frame == src frame
m6("symmetric", f_g=0.40, var_fg=0.004, n_b=800, n_s=10**12, E_g=110.0, E_r=200.0,
   E_spl=100.0, M_b=900.0, S_b=1e-12, Egx=110.0, Erx=200.0)

# ══════════════════════════════════════════════════════════════════════════════════════════════
# M6 ANCHOR TEST — only gDNA is live (pure-gDNA seam). single-λ-via-k is SINGULAR (k→∞). Show the
# per-component gDNA ÷M mode is finite AND that a two-message combine with p_R=0 is NOT over-confident
# here (only one arm active ⇒ no double-count). This is the load-bearing unstranded source.
# ══════════════════════════════════════════════════════════════════════════════════════════════
print("\n═══ M6 ANCHOR — only gDNA live: single-λ-via-k singular; per-component ÷M finite ═══")
n_b = 1200
Mb = _ln(5000.0, 1.0 / n_b, N)
rr = _ln(0.02, 0.015, N)      # enrichment into the anchor's frame (load-bearing)
rho_g = Mb / 140.0            # f_g = 1 exactly
msg_g = rho_g * rr
Mx = 1e6
fgx_anchor = msg_g * 380.0 / Mx     # the ÷M gDNA density mode (a FRACTION, <1)
v_anchor = float(np.var(np.log(fgx_anchor)))
pred_anchor = 1.0 / n_b + 0.015
print(f"  per-component gDNA ÷M mode Var(log f_g^dst) = {v_anchor:.6g}  pred(1/n+σ²t)={pred_anchor:.6g}  "
      f"rel {abs(pred_anchor - v_anchor) / v_anchor:.2%}  (FINITE)")
# k-mode: k = rho_g/rho_R with rho_R→0 ⇒ log k → +inf, Var(log k) undefined/blows up.
print("  k-mode: rho_R→0 ⇒ k→∞, Var(log k) SINGULAR ⇒ single-λ-via-k CANNOT be used at the anchor.")
print("  ⇒ the fix MUST retain the per-component ÷M gDNA mode for pure-gDNA partial messages;")
print("    'single-λ' is single-λ-via-k for both-live PLUS ÷M-gDNA for partial (NOT one uniform rule).")
