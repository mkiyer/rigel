"""RC1 — select the worst full-rank EXON nodes where the messages HURT, ranked by err*mass."""

from __future__ import annotations

import numpy as np

d = np.load("/tmp/suite_nodes.npz", allow_pickle=True)
o = {k: d[k] for k in d.files}
err = np.abs(o["solved"] - o["oracle"])
serr = np.abs(o["self"] - o["oracle"])
em = err * o["mass"]
full = o["tau_own"] > 0.0
hurt = err > serr + 0.02
tgt = full & hurt

print(f"target nodes {int(tgt.sum()):,}  err-mass {em[tgt].sum():,.0f} "
      f"({em[tgt].sum() / em.sum():.1%} of suite)")
print(f"  self mwae {np.average(serr[tgt], weights=o['mass'][tgt]):.4f} -> "
      f"solved {np.average(err[tgt], weights=o['mass'][tgt]):.4f}")

# class / cond breakdown
print("\nby class:")
for c in np.unique(o["cls"][tgt]):
    m = tgt & (o["cls"] == c)
    print(f"  {c:<12}{int(m.sum()):>7}{em[m].sum():>12,.0f}{em[m].sum() / em[tgt].sum():>8.1%}")

print("\ntop 25 by err-mass (all classes):")
idx = np.argsort(np.where(tgt, em, -1))[::-1][:25]
hdr = ("cond", "node", "cls", "dof", "mass", "or", "self", "solv", "tau", "c_tau",
       "lam", "cm_g", "mo_g", "cm_p", "mo_p", "cm_n", "mo_n", "nl", "nr")
print(f"{'#':>3} {'cond':<46}{'node':>6}{'cls':>10}{'dof':>7}{'mass':>10}"
      f"{'orac':>7}{'self':>7}{'solv':>7}{'tau':>8}{'c_tau':>9}{'lam':>7}"
      f"{'cm_g':>9}{'mo_g':>7}{'cm_p':>9}{'mo_p':>7}{'cm_n':>9}{'mo_n':>7}"
      f"  {'L':<11}{'R':<11}")
for r, i in enumerate(idx):
    print(f"{r:>3} {o['cond'][i]:<46}{o['node'][i]:>6}{o['cls'][i]:>10}{o['dof'][i]:>7}"
          f"{o['mass'][i]:>10,.0f}{o['oracle'][i]:>7.3f}{o['self'][i]:>7.3f}{o['solved'][i]:>7.3f}"
          f"{o['tau_own'][i]:>8.3g}{o['c_tau'][i]:>9.3g}{o['lam_fg'][i]:>7.3f}"
          f"{o['cm_g'][i]:>9.3g}{o['mo_g'][i]:>7.3f}{o['cm_p'][i]:>9.3g}{o['mo_p'][i]:>7.3f}"
          f"{o['cm_n'][i]:>9.3g}{o['mo_n'][i]:>7.3f}  {o['nl_cls'][i]:<11}{o['nr_cls'][i]:<11}")
del hdr

print("\ntop 20 EXON only:")
te = tgt & (o["cls"] == "exon")
idx = np.argsort(np.where(te, em, -1))[::-1][:20]
for r, i in enumerate(idx):
    print(f"{r:>3} {o['cond'][i]:<46}{o['node'][i]:>6}{o['dof'][i]:>7}"
          f"{o['mass'][i]:>10,.0f}{o['oracle'][i]:>7.3f}{o['self'][i]:>7.3f}{o['solved'][i]:>7.3f}"
          f"{o['tau_own'][i]:>8.3g}{o['c_tau'][i]:>9.3g}{o['lam_fg'][i]:>7.3f}"
          f"{o['cm_g'][i]:>9.3g}{o['mo_g'][i]:>7.3f}{o['cm_p'][i]:>9.3g}{o['mo_p'][i]:>7.3f}"
          f"{o['cm_n'][i]:>9.3g}{o['mo_n'][i]:>7.3f}  {o['nl_cls'][i]:<11}{o['nr_cls'][i]:<11}")

# which conditions dominate
print("\nby cond (top 12):")
cs = {}
for c in np.unique(o["cond"][tgt]):
    cs[c] = em[tgt & (o["cond"] == c)].sum()
for c, v in sorted(cs.items(), key=lambda kv: -kv[1])[:12]:
    print(f"  {c:<50}{v:>12,.0f}")

# the anchor
m = (o["cond"] == "gdna_gdna300_ss_0.99_nrna_none_capture_on") & (o["node"] == 2197)
i = int(np.flatnonzero(m)[0])
print("\nANCHOR 2197:", {k: o[k][i] for k in o})
