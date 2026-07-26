"""CRITIC #1 — head-to-head attack on the three cross-cliff precision laws.

Each derivation was validated on the scenarios its author CHOSE. This script runs all three laws through the
scenarios where they DIVERGE — the failure modes each author's own MC did not stress:

  LAWS (all expressed as a cliff PENALTY added to the message's honest at-origin log-variance 1/τ_src, then
        τ_msg = 1/(1/τ_src + penalty); fed to the SAME exact 1-D λ-combine ψ does):
    #1  (log r)²                              — cliff MAGNITUDE   (count-independent)
    #2  max(0, (max(1,r)-1)/n_src)            — coverage discount n_eff=n_src/r (count-scaled, r-linear)
    #3  max(0, G² - 1/τ_src - 1/τ_own)        — DerSimonian-Laird mismatch vs the dst's OWN belief (G=λ_msg-λ_own)

  ATTACKS:
    T1  MATCHED-COMPOSITION BIG CLIFF, message NEEDED (own belief weak+WRONG). The owner's case C generalized:
        two co-expressed exons, same composition, 100x different expression (r=100). #1/#2 damp by r → KILL a
        correct, needed message → dst stuck at its wrong own belief. #3 keeps it (G≈0). Is the cliff the right
        variable, or the mismatch?
    T2  AMBIG dst (tau_own=0), TWO competing messages: one matched+correct, one big-cliff+wrong. #1 down-weights
        the bad one; #3 has v_own=∞ ⇒ NO damping on either ⇒ the wrong loud message drags f_g. #3 leaves AMBIG
        UNPROTECTED — the exact regime the M5 win lives in.
    T3  #3 SINGLE-DRAW noise. In the solver G is ONE observed gap per edge, not a 400k-average. At a genuinely
        matched edge (b=0) a single G² is noisy+positively-biased ⇒ #3 randomly KILLS good matched messages.
        Quantify the delivered-precision distribution for a b=0 edge under single-draw DL.
    T4  The AMBIG-WIN-under-CAPTURE survival test for #1/#2: does damping the cross-capture-cliff message that
        CARRIES the unstranded win gut it? Lone message on an AMBIG node — does the MODE still land?

    OMP_NUM_THREADS=1 python scratchpad/cliff_crit_1.py [--draws 400000]
"""

from __future__ import annotations

import argparse

import numpy as np

rng = np.random.default_rng(20260724)


def _logit(p):
    return np.log(p) - np.log1p(-p)


def _sig(x):
    return 1.0 / (1.0 + np.exp(-x))


def lam_combine(lam_own, tau_own, lam_msg, tau_msg):
    """Exact 1-D Gaussian product on λ — what the ψ grid does for the composition DOF."""
    p = tau_own + tau_msg
    return _sig((tau_own * lam_own + tau_msg * lam_msg) / p) if p > 0 else _sig(lam_own)


# ── the three cliff penalties (added to 1/τ_src) ──
def pen1(r, n_src, tau_src, G):
    return np.log(r) ** 2


def pen2(r, n_src, tau_src, G):
    return max(0.0, (max(1.0, r) - 1.0) / n_src)


def pen3(r, n_src, tau_src, G, v_own):
    # DerSimonian-Laird: max(0, G² - v_src - v_own). v_own=inf ⇒ 0 (AMBIG).
    if not np.isfinite(v_own):
        return 0.0
    return max(0.0, G * G - 1.0 / tau_src - v_own)


def tau_msg_of(pen, tau_src):
    return 1.0 / (1.0 / tau_src + pen) if pen >= 0 and tau_src > 0 else (tau_src if pen == 0 else 0.0)


# ════════════════════════════════════════════════════════════════════════════════════════════════════════════
def T1_matched_comp_cliff():
    print("═══ T1  MATCHED-COMPOSITION BIG CLIFF, message NEEDED (own belief weak+WRONG) ═══")
    print("  Two co-expressed exons, SAME composition f_g=0.30, 100x different expression ⇒ r=100 by DENSITY")
    print("  alone. dst's own strand-solve is weak & WRONG (f_own=0.60, τ_own=0.30). The neighbour's correct")
    print("  f_g=0.30 message is the ONLY route to truth. Which laws keep it?\n")
    r = 100.0
    n_src = 50.0
    tau_src = 10.0          # source knows its composition well
    truth = 0.30
    lam_msg = _logit(0.30)  # CORRECT, matched composition
    f_own, tau_own = 0.60, 0.30
    lam_own = _logit(f_own)
    G = lam_msg - lam_own
    v_own = 1.0 / tau_own
    laws = [
        ("#1  (log r)²           ", pen1(r, n_src, tau_src, G)),
        ("#2  (max(1,r)-1)/n_src ", pen2(r, n_src, tau_src, G)),
        ("#3  DL max(0,G²-vs-vo) ", pen3(r, n_src, tau_src, G, v_own)),
        ("    NO-DAMP (baseline) ", 0.0),
    ]
    print(f"   r={r}  G=λ_msg-λ_own={G:+.2f}  truth f_g={truth}  own f_g={f_own} (τ_own={tau_own})\n")
    for name, pen in laws:
        tm = tau_msg_of(pen, tau_src)
        fg = lam_combine(lam_own, tau_own, lam_msg, tm)
        verdict = "RECOVERS" if abs(fg - truth) < 0.08 else ("stuck-wrong" if fg > 0.5 else "partial")
        print(f"   {name} penalty={pen:8.3f}  τ_msg={tm:8.4f}  → f_g={fg:.3f}   ({verdict})")
    print("   READ: the message is CORRECT and needed. #1/#2 damp by the DENSITY cliff (composition matched!)")
    print("   and strand the dst at its wrong 0.60; #3 sees G≈matched and keeps it → 0.30. The owner's invariant")
    print("   ('bigger cliff ⇒ more suspect') is a PRIOR that misfires when the big cliff is composition-MATCHED.\n")


def T2_ambig_competing():
    print("═══ T2  AMBIG dst (τ_own=0): matched-correct msg vs big-cliff-WRONG msg — who protects? ═══")
    print("  The M5-win regime. dst has NO own composition (τ_own=0). Two messages compete: a matched neighbour")
    print("  (r=1.1) carrying truth f_g=0.30, and a big-cliff neighbour (r=400) of a DIFFERENT kind, f_g=0.80.\n")
    truth = 0.30
    n_src = 20.0
    tau_src = 12.0
    lam_good, lam_bad = _logit(0.30), _logit(0.80)
    r_good, r_bad = 1.1, 400.0
    # AMBIG: no own belief. #3's G is undefined (no λ_own) ⇒ v_own=∞ ⇒ penalty 0 for BOTH.
    for name, pgood, pbad in [
        ("#1  (log r)²           ", pen1(r_good, n_src, tau_src, 0), pen1(r_bad, n_src, tau_src, 0)),
        ("#2  (max(1,r)-1)/n_src ", pen2(r_good, n_src, tau_src, 0), pen2(r_bad, n_src, tau_src, 0)),
        ("#3  DL (v_own=∞ on AMBIG)", 0.0, 0.0),
    ]:
        tg = tau_msg_of(pgood, tau_src)
        tb = tau_msg_of(pbad, tau_src)
        p = tg + tb
        fg = _sig((tg * lam_good + tb * lam_bad) / p)
        verdict = "correct" if abs(fg - truth) < 0.1 else "DRAGGED wrong"
        print(f"   {name} τ_good={tg:7.3f} τ_bad={tb:8.4f}  → f_g={fg:.3f}  (truth {truth}; {verdict})")
    print("   READ: #1/#2 down-weight the big-cliff wrong message ~400x ⇒ f_g≈0.30. #3 cannot see a mismatch")
    print("   (no own belief to form G) ⇒ both messages full ⇒ f_g dragged to ~0.5 by the wrong loud one.")
    print("   #3 offers ZERO cliff protection on exactly the AMBIG nodes where messages are the only info.\n")


def T3_dl_single_draw(N):
    print("═══ T3  #3 SINGLE-DRAW noise: at a genuinely MATCHED edge (b=0), one G² randomly kills good msgs ═══")
    print("  In the solver G is ONE observed gap per edge — NOT a 400k average (as agent-#3's MC used). Draw the")
    print("  single-edge G for a matched message (true bias b=0) and report the delivered τ_msg distribution.\n")
    tau_src = 12.0
    v_src = 1.0 / tau_src
    for tau_own in (2.0, 0.5):
        v_own = 1.0 / tau_own
        # matched: msg and own both estimate the same truth, no systematic bias
        msg = rng.normal(0.0, np.sqrt(v_src), N)
        own = rng.normal(0.0, np.sqrt(v_own), N)
        G = msg - own
        pen = np.maximum(0.0, G ** 2 - v_src - v_own)  # single-draw DL excess, PER EDGE
        tm = 1.0 / (1.0 / tau_src + pen)
        frac_killed = float(np.mean(tm < 0.5 * tau_src))
        print(f"   τ_own={tau_own} (v_own={v_own:.2f}): matched edge, TRUE penalty=0")
        print(f"      delivered τ_msg: median {np.median(tm):6.3f}  mean {np.mean(tm):6.3f}  "
              f"(honest={tau_src})   p10={np.percentile(tm,10):.3f} p90={np.percentile(tm,90):.3f}")
        print(f"      fraction of matched edges damped >2x: {frac_killed:5.1%}  "
              f"(these are GOOD messages #3 randomly kills)\n")
    print("   READ: #3's excess is a SINGLE-DRAW rectified statistic (E[max(0,χ²-shift)]>0). On matched edges it")
    print("   over-damps a large, RANDOM fraction — harmless only if the killed message AGREES with the own")
    print("   belief; on a node whose own self-solve is confidently WRONG it silences the good contradicting msg.\n")


def T4_ambig_win_capture():
    print("═══ T4  AMBIG-WIN-under-CAPTURE survival (#1/#2): does damping the cross-cliff msg gut the win? ═══")
    print("  The unstr-capON win: an AMBIG exon (τ_own=0) whose CORRECT composition arrives ONLY as a message")
    print("  from an intronic neighbour across the capture cliff (r large). A LONE message — does the MODE land?\n")
    truth = 0.40
    lam_msg = _logit(truth)   # the correct composition, carried across the cliff
    n_src, tau_src = 30.0, 8.0
    for r in (1.5, 40.0, 400.0):
        print(f"   r={r:6.1f}:")
        for name, pen in [
            ("#1 (log r)²          ", pen1(r, n_src, tau_src, 0)),
            ("#2 (max(1,r)-1)/n_src", pen2(r, n_src, tau_src, 0)),
        ]:
            tm = tau_msg_of(pen, tau_src)
            # LONE message on an AMBIG node (τ_own=0): posterior mode = message mode regardless of τ_msg>0
            fg = lam_combine(0.0, 0.0, lam_msg, tm)
            print(f"      {name} τ_msg={tm:9.5f}  → f_g={fg:.3f}  (truth {truth}; mode {'lands' if abs(fg-truth)<0.02 else 'OFF'})")
    print("   READ: a LONE message's MODE lands regardless of how damped its precision is (τ_msg>0 always) — so")
    print("   #1/#2 do NOT gut a win carried by lone modes. The win is at risk ONLY where the damped message must")
    print("   COMPETE (T2 shows #1/#2 win that; T1 shows they LOSE when the competitor is a wrong own belief).\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=400_000)
    a = ap.parse_args()
    print(f"# CRITIC #1 — head-to-head cliff-law attack   draws={a.draws:,}\n")
    T1_matched_comp_cliff()
    T2_ambig_competing()
    T3_dl_single_draw(a.draws)
    T4_ambig_win_capture()


if __name__ == "__main__":
    main()
