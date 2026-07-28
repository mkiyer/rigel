"""Size the UN-MERGED partition: every distinct exon endpoint becomes a region
boundary (which automatically preserves every TSS/TES, since a terminus is an exon
endpoint). Compare against today's merged-signature partition.
"""

import numpy as np
import pandas as pd

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"

reg = pd.read_feather(D + "regions.feather")
iv = pd.read_feather(D + "intervals.feather")
tx = pd.read_feather(D + "transcripts.feather")
refl = pd.read_feather(D + "ref_lengths.feather")
print("ref_lengths columns:", refl.columns.tolist())

ex = iv[iv["interval_type"] == 0]
real = set(map(int, tx.index[(~tx["is_synthetic"].astype(bool)) & (~tx["is_nrna"].astype(bool))]))
ex = ex[[int(i) in real for i in ex["t_index"].to_numpy()]]

lens_all, lens_by_ref = [], {}
n_cut = 0
for ref, gr in ex.groupby("ref", sort=False):
    cuts = np.union1d(gr["start"].to_numpy(np.int64), gr["end"].to_numpy(np.int64))
    n_cut += cuts.size
    # region lengths of the un-merged partition on this ref (interior only)
    d = np.diff(cuts)
    lens_all.append(d)
    lens_by_ref[ref] = cuts

L_new = np.concatenate(lens_all)
L_old = reg["length"].to_numpy(np.int64)

print()
print(f"today (merged signature)      : {len(L_old):9,} regions, median {int(np.median(L_old)):6,} bp")
print(f"UN-merged (every exon endpoint): {L_new.size + len(lens_by_ref):9,} regions "
      f"(+{100*((L_new.size+len(lens_by_ref))/len(L_old)-1):.0f} %), "
      f"median {int(np.median(L_new)):6,} bp")
print()
qs = [1, 5, 10, 25, 50, 75, 90, 99]
print("  un-merged region length percentiles:",
      {q: int(np.percentile(L_new, q)) for q in qs})
for t in (1, 10, 50, 200):
    print(f"  regions shorter than {t:4d} bp: {100*(L_new < t).mean():5.1f} %"
          f"   (they hold {100*L_new[L_new<t].sum()/L_new.sum():5.2f} % of the spanned bp)")
print()
print(f"  distinct exon endpoints total : {n_cut:,}")
print(f"  1 bp regions                  : {int((L_new==1).sum()):,}")
