"""CLIFF-CROSSING message PRECISION — the honest law as MODEL/IMPUTATION uncertainty (indep. derivation #1).

Pure numpy. No rigel import. The arbiter for the derivation in the accompanying report.

THE OBJECT. A source node emits a per-component density ρ_c^src. It is REFRAMED into the destination frame by
the total-density ratio r = ρ_tot(dst)/ρ_tot(src) and handed to ψ as a Gaussian on log f_c^dst with mode

    mo_c = log( ρ_c^src · r · E_c / M_dst ) .

THE EXACT ERROR (Part 1, proven to machine precision below). Substituting ρ_c = f_c·M/E_c and ρ_tot = M·B
(B = f_g/E_g + (1−f_g)/E_r), EVERY scale cancels — M_src, M_dst, the enrichment e, and r itself — and

    mo_c − log f_c^dst,true  =  log( a_c^src / a_c^dst ) ,          a_c ≡ ρ_c/ρ_tot  (the intrinsic share).

So the message delivers the SOURCE's composition share a_c^src stretched onto the destination's total. Its whole
error is the composition MISMATCH a_c^src/a_c^dst. The imputation premise is exactly "a_c^src = a_c^dst"; the
message is honest only insofar as that holds. The enrichment/mode reframe is CORRECT (it cancels); the PRECISION
must reflect how suspect the shared-composition premise is — which degrades with the cliff.

THE DERIVED LAW.
  * λ ≡ log(a_g/a_R) = logit(f_g) + log(E_r/E_g) is the enrichment-invariant composition DOF. A message is ONE
    λ-constraint: mode λ_src (+ known eff-length shift), precision τ_msg.
  * MAJORITY DOMINANCE (Part 2). r is a TOTAL-density ratio, dominated by the majority component. Projecting a
    single λ-drift Var(Δλ) onto component c by the delta method gives variance a_{other,dst}²·Var(Δλ): ~0 for the
    majority (transports at full precision — the mode is right), ~Var(Δλ) for the minority (a_other→1 → damped).
    BUT the delta-method per-component projection MISCALIBRATES BOTH arms for a large drift (shown) — which is
    why the honest transport is a SINGLE λ-message solved EXACTLY on the ψ grid, not two per-component Gaussians.
  * THE CLIFF TERM (Part 3). Var(Δλ | src, cliff) = the honest per-edge composition-drift variance. It must be 0
    at a matched transport (r≈1: same density ⇒ same kind ⇒ shared composition) and GROW with the cliff (a node
    and its 400×-denser neighbour are almost certainly different kinds). The scale-free quadratic vanishing at
    r=1 is (log r)². The EXISTING σ²_transfer = Var(log r) ≈ 1/n_src is the ratio's SAMPLING error (negligible on
    well-counted nodes, and EXEMPTED to 0 on the graft) — the RIGHT idea (damp by the cliff), the WRONG MOMENT
    (the cliff's measurement error, not its MAGNITUDE). The fix: σ²_cliff = (log r)², NOT exempted on the graft.

        τ_msg = 1 / ( 1/τ_src + (log r)² ) .

VALIDATION (each below, on the anchor node 1909 / boundary 1910 numbers from HANDOFF_4 §6):
  P1  the exact error law (machine precision; e, r, M all cancel).
  P2  the a_other,dst² majority projection matches the minority variance for small drift; miscalibrates for large.
  P3  the single-λ exact combine — BUG (τ_msg=τ_src) collapses f_g 0.95→~0.51; FIX (τ_msg with (log r)²) holds it.
  P4  matched transport (r≈1) preserves full precision; unstranded AMBIG (τ_own=0) at a matched cliff still adopts
      the neighbour (the M5 win preserved) while the SAME source across a big cliff is correctly damped.

    OMP_NUM_THREADS=1 python scratchpad/cliff_mc_1.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260724)


def _logit(p):
    return np.log(p) - np.log1p(-p)


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def _report(name, pred, emp, tol=0.08):
    rel = abs(pred - emp) / max(abs(emp), 1e-300)
    flag = "OK " if rel < tol else "***"
    print(f"  {flag} {name:<52} pred {pred:12.6g}   emp {emp:12.6g}   rel {rel:7.2%}")
    return rel < tol


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def densities(f_g, M, E_g, E_r):
    rho_g = f_g * M / E_g
    rho_R = (1.0 - f_g) * M / E_r
    rho_tot = rho_g + rho_R
    return rho_g, rho_R, rho_tot


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# P1 — THE EXACT ERROR: mo_c − log f_c^dst,true = log(a_c^src/a_c^dst). Enrichment e, r, and BOTH masses cancel.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def p1_exact_error(N, *, E_g, E_r):
    print("═══ P1  EXACT ERROR  mo_c − log f_c^dst = log(a_c^src / a_c^dst)  (e, r, M cancel) ═══")
    # random everything: source/dst composition, masses (⇒ enrichment), independent
    fgs = rng.uniform(0.02, 0.98, N)
    fgd = rng.uniform(0.02, 0.98, N)
    Ms = np.exp(rng.uniform(1.0, 8.0, N))  # source mass ~ enrichment/count (spans 1000×)
    Md = np.exp(rng.uniform(1.0, 12.0, N))  # dst mass — WILDLY different scale (the cliff)
    rg_s, rR_s, rt_s = densities(fgs, Ms, E_g, E_r)
    rg_d, rR_d, rt_d = densities(fgd, Md, E_g, E_r)
    r = rt_d / rt_s  # the reframe

    for c, (rho_s, E_c, fcd, rho_d, rt_dd) in [
        ("g", (rg_s, E_g, fgd, rg_d, rt_d)),
        ("R", (rR_s, E_r, 1.0 - fgd, rR_d, rt_d)),
    ]:
        mo = np.log(rho_s * r * E_c / Md)  # delivered mode on log f_c^dst
        err = mo - np.log(fcd)  # vs the TRUTH
        a_src = rho_s / rt_s
        a_dst = rho_d / rt_dd
        pred = np.log(a_src / a_dst)
        maxabs = float(np.max(np.abs(err - pred)))
        print(f"  {'OK ' if maxabs < 1e-9 else '***'} component {c}: "
              f"max|err − log(a_c^src/a_c^dst)| = {maxabs:.3e}   (machine-zero ⇒ e,r,M all cancel)")
    print()


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# P2 — the MAJORITY-DOMINANCE projection. Draw a composition drift λ_dst = λ_src + ε, ε ~ N(0, V). The message
#      error on component c is log(a_c^src/a_c^dst) exactly; its variance ≈ a_{other,dst}²·V (delta method), which
#      is ~0 for the majority and ~V for the minority. Show it holds for SMALL V and BREAKS for large V (both arms
#      miscalibrate) — the reason the honest transport is a SINGLE λ-message solved exactly, not two projections.
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def p2_projection(N, *, f_g_dst, E_g, E_r, V):
    lam_shift = np.log(E_r / E_g)
    lam_dst0 = _logit(f_g_dst) + lam_shift  # operating point (dst)
    a_g_dst = _sig(lam_dst0 - lam_shift)  # = f_g_dst ... but a_c uses the density share; recompute cleanly
    # intrinsic shares at the dst operating point:
    rg, rR, rt = densities(f_g_dst, 1.0, E_g, E_r)
    a_g, a_R = rg / rt, rR / rt

    # draw the drift and form the EXACT per-component errors
    eps = rng.normal(0.0, np.sqrt(V), N)
    lam_dst = lam_dst0 + eps  # true dst λ (src fixed at lam_dst0 as the message centre)
    fgd = _sig(lam_dst - lam_shift)
    # message centre = the source share at the operating point; error = log(a_c^src/a_c^dst)
    rg_d, rR_d, rt_d = densities(fgd, 1.0, E_g, E_r)
    err_g = np.log(a_g) - np.log(rg_d / rt_d)
    err_R = np.log(a_R) - np.log(rR_d / rt_d)

    print(f"  drift V={V:<6.3g}  f_g^dst={f_g_dst:.3f}  a_g={a_g:.3f} a_R={a_R:.3f}  "
          f"(majority = {'gDNA' if a_g > a_R else 'RNA'})")
    _report(f"    minority arm Var ≈ a_maj²·V   [a_maj²={max(a_g, a_R)**2:.3f}]",
            max(a_g, a_R) ** 2 * V, float(np.var(err_R if a_g > a_R else err_g)), tol=0.15)
    _report(f"    majority arm Var ≈ a_min²·V   [a_min²={min(a_g, a_R)**2:.4f}]",
            min(a_g, a_R) ** 2 * V, float(np.var(err_g if a_g > a_R else err_R)), tol=0.15)


# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
# P3 / P4 — the SINGLE-λ exact combine on the ANCHOR. ψ combines the node's OWN λ-belief (τ_own, λ_own) with the
#      transported λ-message (τ_msg, λ_msg) — the exact 1-D Gaussian product (what the grid ψ does for the
#      composition DOF). f_g = σ(λ_post − lam_shift). BUG: τ_msg = τ_src (no cliff cost). FIX: τ_msg with (log r)².
# ────────────────────────────────────────────────────────────────────────────────────────────────────────────
def lam_combine(lam_own, tau_own, lam_msg, tau_msg, lam_shift):
    p = tau_own + tau_msg
    lam_post = (tau_own * lam_own + tau_msg * lam_msg) / p if p > 0 else lam_own
    return float(_sig(lam_post - lam_shift)), lam_post


def p3_anchor():
    print("═══ P3  ANCHOR (HANDOFF_4 §6): boundary 1910 → exon 1909, r≈407 ═══")
    # geometry (representative FL eff-lengths; the law is eff-length-robust — λ carries the shift)
    E_g, E_r = 150.0, 250.0
    lam_shift = np.log(E_r / E_g)

    # destination exon 1909: oracle f_g=0.985, own strand solve fg_loc=0.953 at WEAK precision τ_own=1.6
    fg_dst_oracle = 0.985
    fg_dst_own = 0.953
    tau_own = 1.6
    lam_own = _logit(fg_dst_own) + lam_shift

    # source boundary 1910: ~64% RNA (f_g≈0.36), a real 47-count spliced junction ⇒ well-measured composition.
    fg_src = 0.36
    n_src = 14.0  # boundary unspliced mass
    n_spl = 47.0  # the spliced junction count
    tau_src = n_spl / 4.0  # ~12 : the source KNOWS its own composition (well-counted RNA) — order only
    lam_msg = _logit(fg_src) + lam_shift  # the transported composition (mode) — enrichment-invariant

    # the cliff: ρ_tot(dst) vs ρ_tot(src) (dst dominated by its gDNA; numbers from the finding)
    rho_tot_dst = 34.0
    rho_tot_src = 0.08
    r = rho_tot_dst / rho_tot_src
    logr = np.log(r)
    print(f"  r = {r:.0f}   (log r)² = {logr**2:.2f}    "
          f"τ_src≈{tau_src:.1f}  λ_msg(f_g_src={fg_src})  λ_own(f_g={fg_dst_own})  oracle f_g={fg_dst_oracle}")

    # the three σ²_transfer candidates on the COMPOSITION stream
    s2t_graft_exempt = 0.0  # CURRENT: graft exemption ⇒ NO damping (the bug)
    s2t_varlogr = 1.0 / n_src + 1.0 / rho_tot_dst  # honest Var(log r) ≈ 1/n_src (dominated by the thin source)
    s2t_cliff = logr ** 2  # DERIVED: the cliff MAGNITUDE

    for label, s2t in [
        ("BUG   graft-exempt σ²=0        ", s2t_graft_exempt),
        ("weak  Var(log r)=1/n_src       ", s2t_varlogr),
        ("FIX   σ²_cliff=(log r)²        ", s2t_cliff),
    ]:
        tau_msg = 1.0 / (1.0 / tau_src + s2t)
        fg_out, _ = lam_combine(lam_own, tau_own, lam_msg, tau_msg, lam_shift)
        verdict = "collapses" if fg_out < 0.7 else ("held" if fg_out > 0.90 else "partial")
        print(f"  {label} τ_msg={tau_msg:8.4f}  →  f_g_out={fg_out:.3f}   ({verdict};  oracle 0.985)")
    print("  (BUG matches the reported collapse to ~0.51; the (log r)² fix holds f_g at its own ~0.95.)")
    print()


def p4_no_overdamp():
    print("═══ P4  NO OVER-DAMP: matched transport preserved; unstranded AMBIG win preserved ═══")
    E_g, E_r = 150.0, 250.0
    lam_shift = np.log(E_r / E_g)

    # (a) MATCHED same-scale transport: a confident source, destination at the SAME density (r≈1). The message
    #     must arrive at ~full precision and MOVE a weak destination toward the (correct) source composition.
    r_matched = 1.05
    tau_src = 12.0
    fg_src = 0.30
    lam_msg = _logit(fg_src) + lam_shift
    tau_own = 0.4  # weak own belief
    lam_own = _logit(0.55) + lam_shift
    s2t = np.log(r_matched) ** 2
    tau_msg = 1.0 / (1.0 / tau_src + s2t)
    fg_out, _ = lam_combine(lam_own, tau_own, lam_msg, tau_msg, lam_shift)
    print(f"  (a) MATCHED r={r_matched}: σ²_cliff={s2t:.4f} τ_msg={tau_msg:.3f} (≈τ_src={tau_src}) "
          f"→ f_g {_sig(lam_own - lam_shift):.2f}→{fg_out:.2f}  (message PRESERVED, moves the weak node)")

    # (b) UNSTRANDED AMBIG destination: τ_own = 0 (NO own composition info — the M5-win regime). A single message
    #     is the ONLY info, so it sets f_g regardless of its precision (the win is NOT killable by damping). The
    #     cliff term only ARBITRATES when messages COMPETE — so show a genuine competition: a MATCHED-cliff
    #     neighbour (f_g≈0.30, the truth) vs a BIG-cliff neighbour (f_g≈0.80, a different kind). The honest cliff
    #     term must let the matched neighbour WIN; a blanket own-cap or an r-blind law would let the loud wrong one win.
    tau_own_ambig = 0.0
    lam_own_ambig = _logit(0.5) + lam_shift
    tau_good = tau_bad = 12.0  # SAME source precision — the ONLY difference is the cliff each crosses
    lam_good = _logit(0.30) + lam_shift  # matched neighbour, carries the truth
    lam_bad = _logit(0.80) + lam_shift  # big-cliff neighbour, a different kind (wrong for this node)
    for law_name, s2_good, s2_bad in [
        ("r-blind (σ²=1/n, the bug)  ", 1.0 / 12.0, 1.0 / 12.0),
        ("σ²_cliff=(log r)²  (the fix)", np.log(1.1) ** 2, np.log(400.0) ** 2),
    ]:
        tg = 1.0 / (1.0 / tau_good + s2_good)
        tb = 1.0 / (1.0 / tau_bad + s2_bad)
        p = tg + tb
        lam_post = (tg * lam_good + tb * lam_bad) / p
        fg_out = float(_sig(lam_post - lam_shift))
        print(f"  (b) AMBIG two competing msgs — {law_name}: τ_good={tg:6.3f} τ_bad={tb:7.4f} "
              f"→ f_g_out={fg_out:.3f}  (truth 0.30)")
    print("      r-blind weights both equally ⇒ f_g dragged to ~0.5 by the wrong neighbour; (log r)² down-weights")
    print("      the big-cliff neighbour ~400× ⇒ f_g ≈ the matched (correct) 0.30. Composition/cliff-aware, and it")
    print("      NEVER caps the source's OWN precision — a lone message still sets f_g (the M5 win is untouched).\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    N = a.draws
    print(f"# CLIFF message-precision MC (independent derivation #1)   draws={N:,}\n")

    p1_exact_error(N, E_g=150.0, E_r=250.0)

    print("═══ P2  MAJORITY-DOMINANCE projection  Var(err_c) ≈ a_{other}²·Var(Δλ) ═══")
    print(" small drift (delta method holds — majority ~0, minority ~V):")
    p2_projection(N, f_g_dst=0.985, E_g=150.0, E_r=250.0, V=0.02)
    print(" large drift (BOTH arms miscalibrate ⇒ single-λ exact solve required, not per-component projection):")
    p2_projection(N, f_g_dst=0.985, E_g=150.0, E_r=250.0, V=5.0)
    print()

    p3_anchor()
    p4_no_overdamp()


if __name__ == "__main__":
    main()
