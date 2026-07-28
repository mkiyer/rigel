"""M5 — the partition AFTER adding TSS/TES as partition events: how much worse
does the tiny-node problem get, and what does the anchor rule recover?"""

import numpy as np
import pandas as pd

from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"
FL = 200

reg = pd.read_feather(D + "regions.feather")
iv = pd.read_feather(D + "intervals.feather")
tx = pd.read_feather(D + "transcripts.feather")

ex = iv[iv["interval_type"] == 0]
real = set(map(int, tx.index[(~tx["is_synthetic"].astype(bool)) & (~tx["is_nrna"].astype(bool))]))
g = ex.groupby("t_index").agg(ref=("ref", "first"), start=("start", "min"), end=("end", "max"))
g = g[[int(i) in real for i in g.index]]
term = {}
for r, s, e in zip(g["ref"].to_numpy(), g["start"].to_numpy(), g["end"].to_numpy()):
    term.setdefault(str(r), set()).update((int(s), int(e)))

sig = reg["signature"].to_numpy(np.int64)
is_exon = ((sig & BIT_EXON_POS) != 0) | ((sig & BIT_EXON_NEG) != 0)

new_len, new_exon = [], []
for ref, gr in reg.groupby("ref_name", sort=False):
    starts = gr["start"].to_numpy(np.int64)
    ends = gr["end"].to_numpy(np.int64)
    exb = is_exon[gr.index.to_numpy()]
    o = np.argsort(starts)
    starts, ends, exb = starts[o], ends[o], exb[o]
    t = np.array(sorted(term.get(ref, ())), dtype=np.int64)
    for s0, e0, xb in zip(starts, ends, exb):
        if t.size:
            lo, hi = np.searchsorted(t, s0, "right"), np.searchsorted(t, e0, "left")
            cuts = t[lo:hi]
        else:
            cuts = np.empty(0, np.int64)
        pts = np.concatenate(([s0], cuts, [e0]))
        d = np.diff(pts)
        new_len.append(d)
        new_exon.append(np.full(d.size, xb))

NL = np.concatenate(new_len)
NE = np.concatenate(new_exon)
OL = reg["length"].to_numpy(np.int64)

print(f"partition TODAY   : {len(OL):,} regions, median {int(np.median(OL))} bp")
print(f"partition + TSS/TES: {len(NL):,} regions, median {int(np.median(NL))} bp"
      f"  ({100*(len(NL)/len(OL)-1):+.1f} %)")
for lbl, L, E in (("TODAY", OL, is_exon), ("+TSS/TES", NL, NE)):
    z = L < FL
    print(f"  {lbl:9s} contained-eff==0 : {100*z.mean():5.1f} % of regions,"
          f" {100*z[E].mean():5.1f} % of EXON regions"
          f" | median EXON region {int(np.median(L[E]))} bp")
print()
print("under the ANCHOR rule every region has eff-length == its own length,")
print(f"so 0.0 % of regions are information-free (vs {100*(NL<FL).mean():.1f} % above);")
print(f"smallest region {NL.min()} bp -> expected count rho*{NL.min()} , honestly imprecise, never undefined.")
