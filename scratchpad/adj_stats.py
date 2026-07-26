"""ADJ — verify the brief's target-set arithmetic from /tmp/suite_nodes.npz (no solver in the loop)."""

from __future__ import annotations

import numpy as np

d = np.load("/tmp/suite_nodes.npz", allow_pickle=True)
cols = list(d.keys())
print("cols:", cols)
cond = d["cond"]
cls = d["cls"]
mass = d["mass"].astype(float)
orc = d["oracle"].astype(float)
slf = d["self"].astype(float)
sol = d["solved"].astype(float)
tau = d["tau_own"].astype(float)
ok = np.isfinite(orc) & (mass > 1e-9)
err_sol = np.abs(sol - orc)
err_self = np.abs(slf - orc)
em = err_sol * mass

tot = em[ok].sum()
print(f"\nALL live nodes           n={ok.sum():7d}  err-mass={tot:14,.0f}  mwae={em[ok].sum()/mass[ok].sum():.4f}")
full = ok & (tau > 0)
zero = ok & (tau <= 0)
for nm, m in (("tau_own>0 (full-rank)", full), ("tau_own==0", zero)):
    print(f"{nm:<24} n={m.sum():7d}  err-mass={em[m].sum():14,.0f} ({100*em[m].sum()/tot:5.1f}%)  "
          f"mwae={em[m].sum()/mass[m].sum():.4f}  self-mwae={(err_self[m]*mass[m]).sum()/mass[m].sum():.4f}")

T = full & (err_sol > err_self + 0.02)
print(f"\nTARGET (full-rank, hurt) n={T.sum():7d}  err-mass={em[T].sum():14,.0f} ({100*em[T].sum()/tot:5.1f}%)  "
      f"mwae={em[T].sum()/mass[T].sum():.4f}  self-mwae={(err_self[T]*mass[T]).sum()/mass[T].sum():.4f}")
for c in np.unique(cls[T]):
    m = T & (cls == c)
    print(f"   {c:<12} n={m.sum():6d} err-mass={em[m].sum():12,.0f} ({100*em[m].sum()/em[T].sum():5.1f}%)")
# stranded x capture split
strd = np.array(["ss_0.99" in c or "ss_0.9" in c for c in cond])
capon = np.array(["capture_on" in c for c in cond])
for nm, m in (("stranded&capON", T & strd & capon), ("stranded&capOFF", T & strd & ~capon),
              ("unstr&capON", T & ~strd & capon), ("unstr&capOFF", T & ~strd & ~capon)):
    print(f"   {nm:<16} n={m.sum():6d} err-mass={em[m].sum():12,.0f} ({100*em[m].sum()/em[T].sum():5.1f}%)")

# direction of the error on the target set: is it f_g pushed DOWN?
dn = T & (sol < orc)
print(f"\nTARGET direction: solved<oracle on {100*em[dn].sum()/em[T].sum():.1f}% of target err-mass "
      f"(n={dn.sum()}/{T.sum()})")

# channel precision shares on the target set (error-mass weighted precision)
for ch, mo, cm in (("gdna", "mo_g", "cm_g"), ("rna+", "mo_p", "cm_p"), ("rna-", "mo_n", "cm_n")):
    m = T
    print(f"  {ch}: mean mode={np.average(d[mo][m].astype(float), weights=em[m]):.3f} "
          f"mean prec={np.average(d[cm][m].astype(float), weights=em[m]):9.2f}")
print(f"  lam : mean prec={np.average(d['c_tau'][T].astype(float), weights=em[T]):9.2f}")

# the gDNA-mode bias by oracle bin, on the target set and on all full-rank
mo_g = d["mo_g"].astype(float)
print("\nmo_g − oracle by oracle bin (mass-weighted):")
bins = [(0.0, 0.3), (0.3, 0.6), (0.6, 0.9), (0.9, 0.99), (0.99, 1.01)]
for lo, hi in bins:
    m = full & (orc >= lo) & (orc < hi)
    if m.sum():
        print(f"  [{lo:.2f},{hi:.2f})  n={m.sum():6d}  mean(mo_g−oracle)={np.average(mo_g[m]-orc[m], weights=mass[m]):+.4f}"
              f"   solved={np.average(sol[m], weights=mass[m]):.4f} self={np.average(slf[m], weights=mass[m]):.4f}"
              f" oracle={np.average(orc[m], weights=mass[m]):.4f}")
