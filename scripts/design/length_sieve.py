"""Point 2: does the CONTAINED count carry gDNA-vs-RNA composition information?

A molecule of length L fits inside a region of length ell in (ell - L + 1) positions, or not at all if
L > ell. So the expected contained count for component c is  rho_c * A_c(ell),  with

    A_c(ell) = E_c[ max(0, ell - L + 1) ]      <- the component's own effective length

If the edge flux already gives the TOTAL density rho = rho_g + rho_r with no FL model (the Sum 1/L
result), then the contained count adds a SECOND equation, and the two are solvable for the split
exactly when A_g(ell) != A_r(ell).  Discriminability is the ratio A_g/A_r: 1.0 = no information,
0.0 = perfect (only RNA can fit).
"""

import numpy as np

L = np.arange(1, 1201)


def eff(mu, sd, ell):
    """A_c(ell) = E[max(0, ell - L + 1)] under a clipped-normal fragment length."""
    w = np.exp(-0.5 * ((L - mu) / sd) ** 2)
    w /= w.sum()
    return float((w * np.maximum(ell - L + 1, 0)).sum())


def show(name, mug, sdg, mur, sdr):
    print(f"\n{name}:  gDNA L ~ {mug}+-{sdg}   RNA L ~ {mur}+-{sdr}")
    print(f"  {'region ell':>9} {'A_gdna':>10} {'A_rna':>10} {'A_g/A_r':>9}   discriminability")
    for ell in (25, 50, 100, 150, 167, 200, 250, 300, 500, 1000, 5000):
        ag, ar = eff(mug, sdg, ell), eff(mur, sdr, ell)
        ratio = ag / ar if ar > 1e-12 else float("nan")
        if not np.isfinite(ratio):
            bar = "both components excluded"
        elif ratio < 0.05:
            bar = "**** near-pure RNA sieve"
        elif ratio < 0.5:
            bar = "***  strong"
        elif ratio < 0.85:
            bar = "**   usable"
        elif ratio < 0.97:
            bar = "*    weak"
        else:
            bar = "     none"
        print(f"  {ell:>9} {ag:>10.3f} {ar:>10.3f} {ratio:>9.4f}   {bar}")


# cfRNA: gDNA is nucleosome-protected (~167 bp mono-nucleosome, tight); RNA differs by protocol.
show("cfRNA-like, well separated", 167, 15, 90, 35)
show("cfRNA-like, modest separation", 167, 30, 200, 60)
show("worst case: identical FL", 180, 50, 180, 50)

print("\n--- how much of the human v8 region partition sits in the informative zone? ---")
print("  measured earlier this project: 8.2 % of regions < 10 bp, 56.7 % < 200 bp,")
print("  median region length 151 bp  -> the MAJORITY of regions are shorter than a fragment.")
