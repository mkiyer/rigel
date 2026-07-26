"""Adjudicator MC — verify the load-bearing claims of the three DL reports.

C1  v_deliv = max(v_msg, G^2 - v_own)  (the self-limiting identity)  -> pin-safety threshold sqrt(2)*sigma_own
C2  |G_lam|/2 <= max_c |G_c| <= |G_lam|   for two NORMALISED compositions (the B-lam vs A-comp 1x-4x claim)
C3  the _EPS FLOOR sensitivity: what does psi actually FEEL as _EPS moves 1e-6 .. 1e-15 ?
C4  the b=0 spurious-damping rate and its cost (R7)
C5  does v_own need +1/n (audit F1)?  -> mo_own = log(og*E_g/M) == log(f_g_own) EXACTLY (M cancels)
"""

import numpy as np

rng = np.random.default_rng(20260725)


def dl(p, gap, v_own):
    v_msg = 1.0 / p
    vo = v_own if np.isfinite(v_own) else 0.0
    b2 = max(0.0, gap * gap - v_msg - vo) if np.isfinite(v_own) else 0.0
    return 1.0 / (v_msg + b2) if b2 > 0 else p


print("=" * 78)
print("C1  self-limiting identity  v_deliv = max(v_msg, G^2 - v_own)")
print("=" * 78)
bad = 0
for _ in range(20000):
    vm = 10 ** rng.uniform(-4, 1)
    vo = 10 ** rng.uniform(-4, 1)
    g = rng.normal(0, 3)
    got = 1.0 / dl(1.0 / vm, g, vo)
    want = max(vm, g * g - vo)
    if abs(got - want) > 1e-9 * max(1.0, abs(want)):
        bad += 1
print(f"  mismatches over 20000 random draws: {bad}   -> identity {'HOLDS' if bad == 0 else 'FAILS'}")
# pin-safety: message out-weighs own  <=>  v_deliv < v_own  <=>  |G| < sqrt(2)*sigma_own  (when vm<vo)
for vo in (0.01, 0.1, 1.0):
    so = np.sqrt(vo)
    for k in (1.0, 1.40, 1.4143, 1.45, 2.0):
        vm = vo / 100.0
        vd = 1.0 / dl(1.0 / vm, k * so, vo)
        print(f"  v_own={vo:<5} |G|={k:5.3f}*sig_own -> w_msg/w_own = {vo / vd:8.4f}"
              f"  {'MSG WINS' if vd < vo else 'own wins'}")

print()
print("=" * 78)
print("C2  |G_lam|/2 <= max_c|G_c| <= |G_lam|  (A-comp under-damps tau by 1x..4x)")
print("=" * 78)
lo, hi = 9e9, -9e9
for _ in range(200000):
    fg_m, fg_o = rng.uniform(1e-3, 1 - 1e-3), rng.uniform(1e-3, 1 - 1e-3)
    Gg = np.log(fg_m / fg_o)
    GR = np.log((1 - fg_m) / (1 - fg_o))
    Gl = Gg - GR
    if abs(Gl) < 1e-6:
        continue
    ratio = Gl * Gl / max(Gg * Gg, GR * GR)  # excess_lam / excess_maxcomp  (>=1 means A-comp under-damps)
    lo, hi = min(lo, ratio), max(hi, ratio)
print(f"  G_lam^2 / max_c G_c^2   range over 2e5 draws: [{lo:.4f}, {hi:.4f}]   (theory [1,4])")
for fo, fm in ((0.95, 0.30), (0.50, 0.10), (0.50, 0.90), (0.20, 0.80), (0.99, 0.97)):
    Gg, GR = np.log(fm / fo), np.log((1 - fm) / (1 - fo))
    print(f"  f_own={fo:<5} f_msg={fm:<5} G_g={Gg:+7.3f} G_R={GR:+7.3f} G_lam={Gg - GR:+7.3f}"
          f"   ratio={((Gg - GR) ** 2) / max(Gg * Gg, GR * GR):5.2f}x")

print()
print("=" * 78)
print("C3  _EPS FLOOR sensitivity — a message with t_c == 0 (dead component)")
print("=" * 78)
print("  A dead component's message mode is mo = log(_EPS); psi feels CURVATURE p_deliv and, at a")
print("  point x on the grid, a PULL p_deliv*(x-mo).  Both are reported below.  o_c = 0.25 (own share),")
print("  x = -10 (the L=10 log-odds window edge, the furthest psi can be dragged).")
print()
o, vo, vm, x = 0.25, 0.05, 0.02, -10.0
print(f"  {'_EPS':>8} {'G':>9} {'p_deliv':>10} {'|pull| @x=-10':>15} {'p_undamped':>11}")
for e in (1e-6, 1e-9, 1e-12, 1e-15):
    mo = np.log(e)
    G = mo - np.log(o)
    p = dl(1.0 / vm, G, vo)
    print(f"  {e:8.0e} {G:9.3f} {p:10.5f} {p * abs(x - mo):15.5f} {1.0 / vm:11.2f}")
print("  -> p_deliv scales as 1/G^2 = 1/(log _EPS)^2 : LOGARITHMIC in _EPS, and the PULL as 1/|log _EPS|.")
print("  -> 1e-6 -> 1e-15 (9 ORDERS of _EPS) moves the pull by:")
p6 = dl(1.0 / vm, np.log(1e-6) - np.log(o), vo) * abs(x - np.log(1e-6))
p15 = dl(1.0 / vm, np.log(1e-15) - np.log(o), vo) * abs(x - np.log(1e-15))
print(f"     {p6:.5f} -> {p15:.5f}  = {p6 / p15:.2f}x   (vs {1.0 / vm:.0f}x if NOT damped at all)")

print()
print("  the MASK alternative (excess=0 when t_c<=0): psi gets mo=log(_EPS) at FULL precision")
for e in (1e-6, 1e-9, 1e-15):
    mo = np.log(e)
    print(f"    _EPS={e:.0e}: p={1.0 / vm:.2f}  pull @x=-10 = {(1.0 / vm) * abs(x - mo):10.2f}"
          f"   ({(1.0 / vm) * abs(x - mo) / (dl(1.0 / vm, mo - np.log(o), vo) * abs(x - mo)):.0f}x the floored-DL pull)")

print()
print("=" * 78)
print("C4  b=0 spurious damping (R7) — rate and COST")
print("=" * 78)
print(f"  P(chi2_1 > 1) = {2 * (1 - 0.8413447460685429):.4f}   E[max(0,chi2_1-1)] = 0.483941 (exact)")
for ratio in (0.1, 1.0, 10.0, 100.0):
    vm0, vo0 = 0.02, 0.02 * ratio
    G = rng.normal(0.0, np.sqrt(vm0 + vo0), 400000)
    vd = np.maximum(vm0, G * G - vo0)
    print(f"  v_own/v_msg={ratio:6.1f}: P(damped)={np.mean(vd > vm0 * (1 + 1e-12)):.4f}"
          f"  E[v_deliv]/v_msg={np.mean(vd) / vm0:8.3f}"
          f"  median v_deliv/v_msg={np.median(vd) / vm0:7.3f}"
          f"  P(v_deliv > v_own)={np.mean(vd > vo0):.4f}")
print("  -> the b=0 spurious damping NEVER pushes the message below the own belief's weight by more than")
print("     the chi2 tail allows, and at b=0 the message AGREES so the fused MODE is unmoved.")

print()
print("=" * 78)
print("C5  audit F1 — does v_own need +1/n ?   mo_own = log(og*E_g/M)")
print("=" * 78)
for M, E_g, fg in ((1234.0, 87.3, 0.42), (7.0, 5.5, 0.9), (1e5, 1e3, 0.01)):
    og = fg * M / E_g
    print(f"  M={M:<9g} E_g={E_g:<7g} f_g={fg:<5}  ->  mo_own=log(og*E_g/M)={np.log(og * E_g / M):+.12f}"
          f"   log(f_g)={np.log(fg):+.12f}   diff={np.log(og * E_g / M) - np.log(fg):+.2e}")
print("  -> M cancels EXACTLY: mo_own IS log f_g^own, so Var(mo_own) is COMPOSITION-ONLY.")
print("     Adding 1/n_dst (audit F1) would price a mass noise that is not in the object. It RAISES")
print("     v_own -> LOWERS v_deliv = max(v_msg, G^2-v_own) -> makes messages STRONGER. Reject.")
for n in (10, 100, 1000):
    vo_c, vo_f = 0.05, 0.05 + 1.0 / n
    G = 1.0
    print(f"  n={n:<5} v_own comp-only={vo_c:.4f} vs F1 {vo_f:.4f}  ->  v_deliv {max(0.02, G * G - vo_c):.4f}"
          f" vs {max(0.02, G * G - vo_f):.4f}   (F1 message is {max(0.02, G * G - vo_c) / max(0.02, G * G - vo_f):.3f}x STRONGER)")
