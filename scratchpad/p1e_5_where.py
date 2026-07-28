"""P1e (2b) — WHERE the surprise actually lives: by how many components the claim supplies, and by mass.
Plus the P1b worked node (2651) message-by-message.

    OMP_NUM_THREADS=1 python scratchpad/p1e_5_where.py
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import p1e_lib as L  # noqa: E402

ORDER = [
    "gdna300_ss0.99_present_capOFF",
    "gdna300_ss0.50_present_capOFF",
    "gdna300_ss0.99_present_capON",
    "gdna100_ss0.50_present_VERYSTRONG",
    "none_ss0.50_present_capOFF",
    "gdna300_ss0.50_none_capOFF",
]
for name in ORDER:
    inp, dbg = L.solve(L.CONDS[name])
    t, nf = L.message_table(inp, dbg)
    live = t["nsup"] > 0
    print(f"\n{'=' * 118}\n{name}\n{'=' * 118}")
    print(f"  {'stratum':<34}{'n':>7}{'reads':>14}{'z2 med':>10}{'z2 p90':>11}"
          f"{'%z2>1':>9}{'%bhat2>0':>10}{'mean bhat2':>12}{'read-share':>12}")
    tot = t["M"][live].sum()
    for tag, m in (
        ("ALL live", live),
        ("  lambda-EMITTING (both comps)", live & t["lam_emit"]),
        ("  ONE-component claims", live & ~t["lam_emit"]),
        ("    ... gDNA-only", live & ~t["lam_emit"] & t["sup_g"] & ~t["sup_r"]),
        ("    ... RNA-only", live & ~t["lam_emit"] & ~t["sup_g"] & t["sup_r"]),
        ("  GRAFT x lambda-emitting", live & t["graft"] & t["lam_emit"]),
        ("  GRAFT x one-component", live & t["graft"] & ~t["lam_emit"]),
    ):
        if m.sum() == 0:
            continue
        z = t["z2"][m]
        print(f"  {tag:<34}{int(m.sum()):>7}{t['M'][m].sum():>14,.0f}{np.median(z):>10.4f}"
              f"{np.percentile(z, 90):>11.2f}{100 * np.mean(z > 1):>8.1f}%"
              f"{100 * np.mean(t['bhat2'][m] > 0):>9.1f}%{np.mean(t['bhat2'][m]):>12.4f}"
              f"{100 * t['M'][m].sum() / max(tot, 1):>11.1f}%")
    # mass carried by z2>1 messages
    hi = live & (t["z2"] > 1)
    print(f"  ==> messages with z2 > 1 carry {100 * t['M'][hi].sum() / max(tot, 1):.1f}% of the "
          f"destination mass ({int(hi.sum())} messages)")

# the P1b worked node
name = "gdna300_ss0.50_present_capOFF"
inp, dbg = L.solve(L.CONDS[name])
t, nf = L.message_table(inp, dbg)
print(f"\n{'=' * 118}\nP1b's worked node 2651 @ {name}\n{'=' * 118}")
k = np.flatnonzero(t["dst"] == 2651)
print(f"  node: class={nf['cls'][2651]}  M={nf['M'][2651]:,.0f}  E_g={nf['E_g'][2651]:,.0f}  "
      f"E_r={nf['E_r'][2651]:,.0f}  oracle f_g={nf['fo'][2651]:.4f}  solved f_g={nf['fg'][2651]:.4f}  "
      f"sd={np.sqrt(nf['var_g'][2651]):.4f}")
for i in k:
    print(f"  msg df={t['df'][i]} src={t['src'][i]} ({t['cls_src'][i]}) graft={bool(t['graft'][i])}: "
          f"S={t['S'][i]:,.0f}  M={t['M'][i]:,.0f}  claim/obs={t['S'][i] / t['M'][i]:.3f}  "
          f"delta={t['delta'][i]:+.4f}  aSa={t['aSa'][i]:.5f}  1/n_dst={1 / t['n_dst'][i]:.2e}  "
          f"z2={t['z2'][i]:,.1f}  bhat2={t['bhat2'][i]:.4f}")
    print(f"      alpha=(g {t['alpha_g'][i]:.3f}, p {t['alpha_p'][i]:.3f}, n {t['alpha_n'][i]:.3f})  "
          f"s=(g {t['s_g'][i]:.4g}, p {t['s_p'][i]:.4g}, n {t['s_n'][i]:.4g})  "
          f"sigma_cm2={t['s2cm'][i]:.4g}")
    print(f"      DERIVED damping: dv_gDNA={t['dv_g'][i]:.4f}  dv_RNA={t['dv_r'][i]:.4f}  "
          f"gDNA share={t['dv_g'][i] / max(t['dv_g'][i] + t['dv_r'][i], 1e-12):.3f}   "
          f"delivered f_g={t['fg_msg'][i]:.4f} vs oracle {t['fo'][i]:.4f}")
