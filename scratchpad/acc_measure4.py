"""M1b — the contained effective length under a REALISTIC FL pmf (not a delta).

`region_eff_length(L) = E[max(0, L - l + 1)]` is exactly 0 only when the FL pmf has
no mass below L. With a real, dispersed FL it is positive but collapses. Report the
honest quantity: E and the *efficiency* E/L, versus the anchor rule's E == L exactly.
"""

import numpy as np
import pandas as pd

from rigel.calibration.effective_length import region_eff_length
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS

D = "/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index_v7/"


def gauss_pmf(mu, sd, n=2001):
    x = np.arange(n, dtype=np.float64)
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    p[0] = 0.0
    return p / p.sum()


reg = pd.read_feather(D + "regions.feather")
L = reg["length"].to_numpy(np.int64)
sig = reg["signature"].to_numpy(np.int64)
is_exon = ((sig & BIT_EXON_POS) != 0) | ((sig & BIT_EXON_NEG) != 0)

PMFS = {
    "cfRNA-like  mu=180 sd=60": gauss_pmf(180, 60),
    "gDNA-like   mu=300 sd=90": gauss_pmf(300, 90),
    "short       mu=120 sd=40": gauss_pmf(120, 40),
}

print(f"{len(L):,} regions ({int(is_exon.sum()):,} exon); median exon length "
      f"{int(np.median(L[is_exon]))} bp\n")
print(f"{'FL pmf':26s} {'pop':6s} {'med E':>9s} {'med E/L':>9s} "
      f"{'E<1':>7s} {'E/L<0.1':>8s} {'E/L<0.5':>8s}")
print("-" * 78)
for name, pmf in PMFS.items():
    eff = region_eff_length(L.astype(np.float64), pmf)
    for pop, m in (("all", np.ones(len(L), bool)), ("exon", is_exon)):
        e, ll = eff[m], L[m].astype(np.float64)
        ratio = e / ll
        print(f"{name:26s} {pop:6s} {np.median(e):9.1f} {np.median(ratio):9.3f} "
              f"{100*(e<1).mean():6.1f}% {100*(ratio<0.1).mean():7.1f}% "
              f"{100*(ratio<0.5).mean():7.1f}%")
print("\nANCHOR RULE: E == L exactly, for every region and every FL pmf")
print("             -> med E/L = 1.000, E<1 only for 1 bp regions, E/L<0.5 = 0.0 %")

# where the mass is: exon regions weighted by their own length
e = region_eff_length(L.astype(np.float64), PMFS["cfRNA-like  mu=180 sd=60"])
w = L[is_exon].astype(np.float64)
r = (e[is_exon] / w)
print(f"\ncfRNA-like FL, EXON regions, bp-weighted mean efficiency E/L = "
      f"{np.average(r, weights=w):.3f}; unweighted mean = {r.mean():.3f}")
print(f"exon regions with E < 10 (fewer than 10 start positions of opportunity): "
      f"{100*(e[is_exon]<10).mean():.1f} %")
