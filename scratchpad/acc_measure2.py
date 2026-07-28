"""M3b — the crossing census restricted to EXONIC start positions (the RNA-relevant
view), plus the same census on the 10 Mb synthetic suite index for comparison.

Also M4: how much of a region's own bp lies within FL of one of its interfaces
(the "edge zone"), which is where the contained/crossing split is decided.
"""


import numpy as np
import pandas as pd

from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS

FL = 200


def census(d, label, exon_only):
    reg = pd.read_feather(d + "regions.feather")
    sig = reg["signature"].to_numpy(np.int64)
    is_exon = ((sig & BIT_EXON_POS) != 0) | ((sig & BIT_EXON_NEG) != 0)
    tot = 0
    cross = np.zeros(6, np.int64)
    for ref, gr in reg.groupby("ref_name", sort=False):
        starts = gr["start"].to_numpy(np.int64)
        ends = gr["end"].to_numpy(np.int64)
        ex = is_exon[gr.index.to_numpy()]
        order = np.argsort(starts)
        starts, ends, ex = starts[order], ends[order], ex[order]
        interior = ends[:-1]
        if ends[-1] - starts[0] <= FL:
            continue
        keep = ex if exon_only else np.ones(len(starts), bool)
        for s0, e0 in zip(starts[keep], ends[keep]):
            hi = min(e0, ends[-1] - FL)
            if hi <= s0:
                continue
            s = np.arange(s0, hi, dtype=np.int64)
            k = np.searchsorted(interior, s + FL, "right") - np.searchsorted(interior, s, "right")
            tot += s.size
            cross += np.bincount(np.minimum(k, 5), minlength=6)
    if tot == 0:
        print(f"{label}: no positions")
        return
    nc = tot - cross[0]
    print(f"{label}  ({'EXONIC' if exon_only else 'ALL'} starts, n={tot:,})")
    print(f"    crosses 0 : {100*cross[0]/tot:6.2f} %   >=1 : {100*nc/tot:6.2f} %"
          f"   >=2 : {100*(nc-cross[1])/tot:6.2f} %   >=3 : {100*cross[3:].sum()/tot:6.2f} %")
    print(f"    of the CROSSERS, {100*(nc-cross[1])/max(nc,1):5.2f} % cross more than one boundary")


HUMAN = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"
TOY = "/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb/rigel_index/"

for d, name in ((HUMAN, "HUMAN v7"), (TOY, "TOY ambig_dense_10mb")):
    try:
        census(d, name, exon_only=True)
        census(d, name, exon_only=False)
    except Exception as e:  # noqa: BLE001
        print(f"{name}: {e}")
    print()

# M4 — region-size distribution on the toy, for the same tiny-node census
for d, name in ((HUMAN, "HUMAN v7"), (TOY, "TOY ambig_dense_10mb")):
    try:
        reg = pd.read_feather(d + "regions.feather")
        L = reg["length"].to_numpy(np.int64)
        sig = reg["signature"].to_numpy(np.int64)
        is_exon = ((sig & BIT_EXON_POS) != 0) | ((sig & BIT_EXON_NEG) != 0)
        z = L < FL
        print(f"{name}: {len(L):,} regions, median {int(np.median(L))} bp; "
              f"contained-eff==0 (L<{FL}): {100*z.mean():.1f} % overall, "
              f"{100*z[is_exon].mean():.1f} % of EXON regions "
              f"({int(is_exon.sum()):,} exon regions)")
    except Exception as e:  # noqa: BLE001
        print(f"{name}: {e}")
