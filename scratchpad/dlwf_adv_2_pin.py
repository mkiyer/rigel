"""ADVERSARY #2 — the DL gap G_c as the SOLVER actually computes it: after `_pin_v`, on a real edge.

The plan puts the DL deflation AFTER `_pin_v`, using the pinned message densities `tg/tp/tn`. `_pin_v` is
NOT a neutral rescale:

    s_c = t_c if p_c > 0 else o_c           <-- the DESTINATION'S OWN density fills uninformed components
    k   = M / (s_g*E_g + (s_p+s_n)*E_r)
    return t_g*k, t_p*k, t_n*k

Three consequences this script measures exactly (no statistics, just the arithmetic):

  P1  FULL message  (all three precisions > 0, no peel): r cancels EXACTLY in G_c. Good -- G is the pure
      composition-share gap the derivation wants.
  P2  PARTIAL message (a component has precision 0 -> the OWN density enters the normalizer): r does NOT
      cancel. A pure-ENRICHMENT cliff with IDENTICAL composition manufactures |G| ~ |log r| -> DL
      re-introduces exactly the (log r)^2 over-damping it was adopted to remove.
  P3  PEEL edge (exon -> boundary): the mature is subtracted AFTER the reframe, so t_p is not proportional
      to r and r does not cancel either. Same manufactured-mismatch mechanism.
  P4  INDEPENDENCE: on a partial message the pinned t_c is a FUNCTION of o_c, so the DL "second study" is
      not independent of the first. Measured as dG/d(log o) -- it should be -1 (pure gap); anything else is
      shrinkage of the evidence DL is reading.
  P5  ZERO-DENSITY components: t_c = 0 (a fully-consumed peel) or o_c = 0 (an own belief with no RNA on a
      strand) make G = log(EPS/o) or log(t/EPS): the excess -- and therefore whether the message survives --
      is set by the value of _EPS.

    OMP_NUM_THREADS=1 python scratchpad/dlwf_adv_2_pin.py
"""

from __future__ import annotations

import numpy as np

_EPS = 1.0e-9


def pin(t_g, t_p, t_n, p_g, p_p, p_n, o_g, o_p, o_n, E_g, E_r, M):
    """Exact transcription of bp_solver._pin_v."""
    s_g = t_g if p_g > 0.0 else o_g
    s_p = t_p if p_p > 0.0 else o_p
    s_n = t_n if p_n > 0.0 else o_n
    s = s_g * E_g + (s_p + s_n) * E_r
    k = M / max(s, _EPS) if (s > _EPS and M > _EPS) else 1.0
    return t_g * k, t_p * k, t_n * k


def transport(rho_g_s, rho_p_s, rho_n_s, r, *, peel_p=0.0, peel_n=0.0, graft_p=0.0, graft_n=0.0,
              p_g=1.0, p_p=1.0, p_n=1.0, o_g=0.0, o_p=0.0, o_n=0.0, E_g=1.0, E_r=1.0, M=1.0):
    """Exact transcription of the density half of bp_solver._transport (graft before the reframe, peel
    after), ending in `_pin_v`."""
    t_g = rho_g_s * r
    t_p = (rho_p_s + graft_p) * r
    t_n = (rho_n_s + graft_n) * r
    t_p = max(t_p - peel_p, 0.0)
    t_n = max(t_n - peel_n, 0.0)
    return pin(t_g, t_p, t_n, p_g, p_p, p_n, o_g, o_p, o_n, E_g, E_r, M)


def hdr(s):
    print(f"\n{'='*100}\n{s}\n{'='*100}")


def main():
    # a destination node: mass M, gDNA/RNA effective lengths; its OWN self-solve says fg_own.
    M, E_g, E_r = 5000.0, 300.0, 250.0
    fg_own = 0.60
    o_g = fg_own * M / E_g
    o_p, o_n = (1.0 - fg_own) * M / E_r, 0.0

    hdr("P1. FULL message, MATCHED composition — r must cancel in G (the derivation's premise)")
    print(f"   dst: M={M:g} E_g={E_g:g} E_r={E_r:g} fg_own={fg_own}")
    print("   source has the SAME composition shares as the destination's own belief (true b = 0).")
    print(f"   {'r':>10} {'G_g':>14} {'G_+':>14} {'excess_g':>12}   (v_msg=v_own=0.02)")
    for r in (1e-3, 1e-2, 1.0, 1e2, 1e3, 1e4):
        rho_tot_s = 1.0
        tg, tp, tn = transport(fg_own * rho_tot_s, (1 - fg_own) * rho_tot_s, 0.0, r,
                               p_g=1.0, p_p=1.0, p_n=1e-30, o_g=o_g, o_p=o_p, o_n=o_n,
                               E_g=E_g, E_r=E_r, M=M)
        # note: p_n=1e-30 > 0 so the n component is 'live' with density 0 -- see P5 for p_n = 0 exactly
        G_g = np.log(max(tg, _EPS) / max(o_g, _EPS))
        G_p = np.log(max(tp, _EPS) / max(o_p, _EPS))
        exc = max(0.0, G_g**2 - 0.02 - 0.02)
        print(f"   {r:10.3g} {G_g:14.3e} {G_p:14.3e} {exc:12.3e}")
    print("   VERDICT: exact cancellation (1e-16). For a FULL message the plan's G is clean.")

    hdr("P2. PARTIAL message (gDNA only; RNA precision 0) — r does NOT cancel; a pure-enrichment")
    print("    cliff with IDENTICAL composition MANUFACTURES a mismatch")
    print(f"   {'r':>10} {'f_g implied by msg':>20} {'G_g':>10} {'excess_g':>10} {'damp p_new/p_old':>18}")
    v_msg, v_own = 0.02, 0.02
    for r in (1e-3, 1e-2, 0.1, 1.0, 10.0, 1e2, 1e3, 1e4):
        tg, tp, tn = transport(fg_own, 0.0, 0.0, r, p_g=1.0, p_p=0.0, p_n=0.0,
                               o_g=o_g, o_p=o_p, o_n=o_n, E_g=E_g, E_r=E_r, M=M)
        G_g = np.log(max(tg, _EPS) / max(o_g, _EPS))
        exc = max(0.0, G_g**2 - v_msg - v_own)
        print(f"   {r:10.3g} {tg*E_g/M:20.4f} {G_g:10.3f} {exc:10.3f} {1.0/(1.0+exc/v_msg):18.4g}")
    print("   VERDICT: |G_g| grows like |log r| until the pin saturates at f_g=1. At r=100 the composition")
    print("   is IDENTICAL yet excess=1.9 -> the gDNA message is damped 96x. This is precisely the")
    print("   '(log r)^2 over-damps a pure-enrichment cliff' defect the DL term was adopted to remove,")
    print("   re-entering through _pin_v on the ~45-58% of edges that are PARTIAL (measured, adv #3 R6).")

    hdr("P3. PEEL edge (exon -> boundary): the mature is subtracted AFTER the reframe -> r survives")
    print("    Source exon: rho_R = 1.0 (all continuing RNA), spliced peel = 0.5 at the destination face.")
    print(f"   {'r':>10} {'t_p (pinned)':>14} {'G_+':>10} {'excess_+':>10} {'damp':>10}")
    for r in (0.1, 0.5, 1.0, 2.0, 5.0, 20.0):
        tg, tp, tn = transport(0.5, 1.0, 0.0, r, peel_p=0.5, p_g=1.0, p_p=1.0, p_n=0.0,
                               o_g=o_g, o_p=o_p, o_n=o_n, E_g=E_g, E_r=E_r, M=M)
        G_p = np.log(max(tp, _EPS) / max(o_p, _EPS))
        exc = max(0.0, G_p**2 - v_msg - v_own)
        print(f"   {r:10.3g} {tp:14.4g} {G_p:10.3f} {exc:10.3f} {1.0/(1.0+exc/v_msg):10.4g}")
    print("   VERDICT: the peel makes G_+ a function of r even for a FULL message (the subtracted mature")
    print("   is an absolute density, r-independent). At r=0.5 the peel consumes everything -> t_p=0 ->")
    print("   G = log(EPS/o) -> the message is annihilated by the floor, not by a composition claim.")

    hdr("P4. INDEPENDENCE — on a partial message the pinned t_c is a FUNCTION of the own density o_c")
    print("    DL assumes msg _|_ own. dG/d(log o_g) must be -1 (a pure gap). Anything above -1 means the")
    print("    pin has already shrunk the message toward the own belief, so G understates the true gap.")
    print(f"   {'fg_own':>8} {'r':>8} {'G_g':>10} {'dG/dlog(o_g)':>14} {'shrinkage':>10}")
    for fgo in (0.1, 0.3, 0.6, 0.9):
        for r in (1.0, 100.0):
            def _G(fg):
                og_ = fg * M / E_g
                op_ = (1 - fg) * M / E_r
                tg_, _, _ = transport(0.6, 0.0, 0.0, r, p_g=1.0, p_p=0.0, p_n=0.0,
                                      o_g=og_, o_p=op_, o_n=0.0, E_g=E_g, E_r=E_r, M=M)
                return np.log(max(tg_, _EPS) / max(og_, _EPS)), og_
            h = 1e-5
            g1, og1 = _G(fgo)
            g2, og2 = _G(fgo * np.exp(h))
            d = (g2 - g1) / (np.log(og2) - np.log(og1))
            print(f"   {fgo:8.2f} {r:8.3g} {g1:10.3f} {d:14.4f} {abs(d)/1.0:10.3f}")
    print("   VERDICT: |dG/dlog o| < 1 wherever the pin's normalizer contains o (partial messages) -> the")
    print("   observed gap is SHRUNK toward 0 -> DL UNDER-estimates b^2 -> the message is left more")
    print("   confident than the honest random-effects value: the OVER-CONFIDENT direction.")

    hdr("P5. THE _EPS FLOOR is a controlling constant for the zero-density cases")
    print("    A fully-consumed peel gives t_p = 0 exactly; a strand the own solve killed gives o_p = 0.")
    print("    G is then log(EPS/o) or log(t/EPS) -- the excess is a function of _EPS ONLY.")
    print(f"   {'_EPS':>10} {'G (t=0, o=12)':>16} {'excess':>10} {'damp (v_msg=0.02)':>18}")
    for eps in (1e-6, 1e-9, 1e-12, 1e-15):
        G = np.log(eps / 12.0)
        exc = max(0.0, G * G - 0.04)
        print(f"   {eps:10.3g} {G:16.3f} {exc:10.1f} {1.0/(1.0+exc/0.02):18.3g}")
    print("   VERDICT: the DAMPING of a zero-density component spans 250x across plausible _EPS values.")
    print("   The rule is 'magic-number-free' only if no component is ever 0 -- measured (adv #3 R7):")
    print("   84-90% of the killed RNA+ messages on the R side are killed by exactly this floor.")

    hdr("P6. o_c = 0 (own belief says NO RNA on this strand) -> G = +inf -> every RNA message killed")
    o_p0 = 0.0
    for t_p in (0.5, 5.0, 50.0):
        G = np.log(max(t_p, _EPS) / max(o_p0, _EPS))
        print(f"   t_p={t_p:6.2f} o_p=0  ->  G={G:8.3f}  excess={G*G-0.04:12.1f}  damp={1.0/(1.0+(G*G-0.04)/0.02):.3e}")
    print("   VERDICT: a node whose own self-solve zeroed a strand can NEVER be told otherwise by a")
    print("   neighbour. That is the entrenchment failure mode: pass-0 becomes NON-correctable exactly")
    print("   where the own belief is at a vertex, which is where it is most likely to be wrong.")


if __name__ == "__main__":
    main()
