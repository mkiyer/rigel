"""RC2 (A) step 1 — reproduce the gDNA-measurement mode bias from the suite table and find its correlates.

Pure re-analysis of /tmp/suite_nodes.npz (no solve). Establishes the population, the binned bias, and which
STRUCTURAL covariates it tracks — before any transport decomposition.
"""

from __future__ import annotations

import numpy as np

d = dict(np.load("/tmp/suite_nodes.npz", allow_pickle=True))
for k in d:
    d[k] = np.asarray(d[k])

oracle, self_, solved, mass = d["oracle"], d["self"], d["solved"], d["mass"]
err, serr = np.abs(solved - oracle), np.abs(self_ - oracle)
em = err * mass
full = d["tau_own"] > 0.0
hurt = full & (err > serr + 0.02)
exon = d["cls"] == "exon"
live_g = d["cm_g"] > 0.0

print(f"all nodes {err.size:,}  err-mass {em.sum():,.0f}  mwae {np.average(err, weights=mass):.4f}")
print(f"full-rank {full.sum():,}  hurt {hurt.sum():,}  hurt err-mass {em[hurt].sum():,.0f} "
      f"({em[hurt].sum() / em.sum():.1%})")

BINS = [(0.0, 0.30), (0.30, 0.60), (0.60, 0.90), (0.90, 0.99), (0.99, 1.01)]


def binned(lab, m, col="mo_g"):
    print(f"\n── {lab}   n={int(m.sum()):,}")
    print(f"{'oracle bin':<14}{'n':>7}{'mass':>13}{'mean ' + col:>11}{'mean oracle':>13}"
          f"{'bias':>9}{'mean self':>11}{'self bias':>11}{'mean solved':>13}")
    for lo, hi in BINS:
        s = m & (oracle >= lo) & (oracle < hi)
        if not s.any():
            continue
        w = mass[s]
        mo, orc = np.average(d[col][s], weights=w), np.average(oracle[s], weights=w)
        sf = np.average(self_[s], weights=w)
        print(f"[{lo:.2f},{hi:.2f})   {int(s.sum()):>7}{w.sum():>13,.0f}{mo:>11.3f}{orc:>13.3f}"
              f"{mo - orc:>+9.3f}{sf:>11.3f}{sf - orc:>+11.3f}{np.average(solved[s], weights=w):>13.3f}")


binned("FULL-RANK EXONS w/ live gDNA msg (the reported population)", full & exon & live_g)
binned("  ... restricted to HURT", hurt & exon & live_g)
binned("FULL-RANK exons, live gDNA msg — RNA+ measurement mode", full & exon & live_g, col="mo_p")
binned("FULL-RANK exons, live gDNA msg — lambda message implied f_g", full & exon & live_g, col="lam_fg")

# ── is the bias unweighted too, or a mass artifact? ──
print("\n\n=== unweighted (plain mean) on full-rank exons w/ live gDNA msg ===")
m = full & exon & live_g
for lo, hi in BINS:
    s = m & (oracle >= lo) & (oracle < hi)
    if not s.any():
        continue
    print(f"[{lo:.2f},{hi:.2f})  n={int(s.sum()):>6}  mo_g={d['mo_g'][s].mean():.3f} "
          f"oracle={oracle[s].mean():.3f}  bias={d['mo_g'][s].mean() - oracle[s].mean():+.3f} "
          f"median bias={np.median(d['mo_g'][s] - oracle[s]):+.3f}")

# ── does the mode track the NEIGHBOUR's oracle rather than the node's own?  The messages come from the
#    left/right neighbours, so if the mode is an honest imputation it should track THEIR composition. ──
print("\n\n=== mo_g vs the node's own oracle vs its NEIGHBOURS' oracle (full-rank exons, live msg) ===")
nl, nr, ml, mr = d["nl_oracle"], d["nr_oracle"], d["nl_mass"], d["nr_mass"]
both = m & np.isfinite(nl) & np.isfinite(nr)
wn = ml + mr
nb = np.where(wn > 0, (nl * ml + nr * mr) / np.maximum(wn, 1e-9), np.nan)
for lo, hi in BINS:
    s = both & (oracle >= lo) & (oracle < hi)
    if not s.any():
        continue
    w = mass[s]
    print(f"[{lo:.2f},{hi:.2f})  n={int(s.sum()):>6}  mo_g={np.average(d['mo_g'][s], weights=w):.3f}  "
          f"own oracle={np.average(oracle[s], weights=w):.3f}  "
          f"nbr oracle(mass-wtd)={np.average(nb[s], weights=w):.3f}  "
          f"|mo_g-own|={np.average(np.abs(d['mo_g'][s] - oracle[s]), weights=w):.3f}  "
          f"|mo_g-nbr|={np.average(np.abs(d['mo_g'][s] - nb[s]), weights=w):.3f}")

# ── neighbour CLASS composition of the hurt exons ──
print("\n=== neighbour classes of hurt exons (err-mass) ===")
hm = hurt & exon
pairs: dict[str, float] = {}
for a_, b_ in zip(d["nl_cls"][hm], d["nr_cls"][hm]):
    pairs[f"{a_}|{b_}"] = 0.0
for a_, b_, e_ in zip(d["nl_cls"][hm], d["nr_cls"][hm], em[hm]):
    pairs[f"{a_}|{b_}"] += e_
for k, v in sorted(pairs.items(), key=lambda kv: -kv[1])[:10]:
    print(f"   {k:<26}{v:>12,.0f}  {v / em[hm].sum():>6.1%}")

# ── conditions ──
print("\n=== hurt err-mass by condition (top 12) ===")
cm: dict[str, float] = {}
for c_, e_ in zip(d["cond"][hurt], em[hurt]):
    cm[c_] = cm.get(c_, 0.0) + e_
for k, v in sorted(cm.items(), key=lambda kv: -kv[1])[:12]:
    print(f"   {k:<52}{v:>12,.0f}")
