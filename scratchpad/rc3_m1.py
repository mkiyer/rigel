"""RC3 — exposure of the two DERIVED-BUT-UNWIRED guards (M1 near-wall, M3 peel-u) and of the
missing-frame pass-through, measured on the shipped solve.

    OMP_NUM_THREADS=1 python scratchpad/rc3_m1.py [cond ...]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from rc3_replay import build, load, oracle  # noqa: E402

_EPS = 1e-9
conds = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
]

sig_minor, f_minor, w_em = [], [], []
nframe = np.zeros(3)
uu, u_em = [], []
for cond in conds:
    ctx = load(cond)
    S = build(ctx)
    fo = oracle(ctx)
    ok = np.isfinite(fo) & (S["mass"] > _EPS)
    em = np.where(ok, np.abs(S["solved"] - fo) * S["mass"], 0.0)
    tau, fg = S["tau_own"], S["fg_loc"]
    live = (tau > _EPS) & ~S["struct"]
    # the transported arm's log-variance, per the FOUNDATION Jacobians (`own_composition_logvar`)
    with np.errstate(divide="ignore", invalid="ignore"):
        v_g = np.where(live, (1 - fg) ** 2 / np.maximum(tau, _EPS), np.nan)
        v_r = np.where(live, fg**2 / np.maximum(tau, _EPS), np.nan)
    minor_is_g = fg < 0.5
    fm = np.where(minor_is_g, fg, 1 - fg)
    vm = np.where(minor_is_g, v_g, v_r)
    m = live & np.isfinite(vm)
    sig_minor.append(np.sqrt(vm[m]))
    f_minor.append(fm[m])
    w_em.append(em[m])
    for D, val in ((S["A"], S["vl"]), (S["B"], S["vr"])):
        framed = D["r"] > 0.0
        nframe += [val.sum(), (val & (D["r"] == 1.0)).sum(), em[val & (D["r"] == 1.0)].sum()]
        pk = D["peel"] & (D["tp_post_peel"] > _EPS)
        uu.append(D["tp_pre_peel"][pk] / D["tp_post_peel"][pk])
        u_em.append(em[pk])
    print(f"  {cond}", flush=True)

s = np.concatenate(sig_minor)
f = np.concatenate(f_minor)
w = np.concatenate(w_em)
print(f"\n── M1 near-wall exposure (the log-Jacobian under-states Var(log f_c) 36–92% for a small, "
      f"high-CV minority arm) ──\nn nodes emitting a composition = {s.size:,}")
for lo, hi in ((0, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, np.inf)):
    m = (s >= lo) & (s < hi)
    print(f"  σ(log f_minor) ∈ [{lo},{hi}): n={m.sum():>7,}  err-mass={w[m].sum():>12,.0f} "
          f"({w[m].sum() / max(w.sum(), 1):>5.1%})  median f_minor={np.median(f[m]) if m.any() else 0:.4g}")
print(f"\n── missing FRAME (r forced to 1.0: the message is delivered in the SOURCE's frame at full "
      f"precision) ──\n  valid edges={nframe[0]:,.0f}  r==1 pass-through={nframe[1]:,.0f} "
      f"({nframe[1] / nframe[0]:.2%})  dst err-mass={nframe[2]:,.0f}")
U = np.concatenate(uu)
print(f"\n── M3 peel u (linearization valid only to u≈3) ── n={U.size:,} "
      f"median={np.median(U):.3g} p75={np.percentile(U, 75):.3g} p90={np.percentile(U, 90):.3g} "
      f"p99={np.percentile(U, 99):.4g} frac(u>3)={np.mean(U > 3):.2%}")
