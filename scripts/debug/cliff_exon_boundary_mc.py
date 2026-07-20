"""Monte-Carlo validation of the EXON region <-> intron-exon BOUNDARY message arithmetic — the mature
reconciliation layered on the eff-length-frame log-odds shift (docs/calibration/cliff_message_derivation.md §8).

Physical link: mature RNA has the SAME density rho_m in two frames — contained-unspliced in the exon body
(eff-len region_eff_length(RNA_FL, L)) and spliced-crossing at the junction (eff-len ~fl_mean(RNA_FL); the code
uses the half-triangle spliced_side_eff_length for the MASS deposit — the arithmetic is identical, only the
count/mass convention differs). The exon's RNA (contained) = nascent + mature; the boundary's UNSPLICED crossing
= nascent only (mature splices out -> spliced reads). So:
  * boundary -> exon: ADD the mature density (from the boundary spliced) to the nascent, then the shift.
  * exon -> boundary: SUBTRACT the mature density (the boundary's own spliced) from the exon RNA, then the shift.

Geometry: exon E = [0, Le]; a fragment with start in E and end > Le either crosses Le unspliced (gDNA / nascent,
into the intron) or, for mature, SPLICES (spans E -> downstream exon in mature coords). Ground-truth deposits
recover rho_g, rho_n, rho_m in both frames; we check frame-consistency, the reconciliation, and the f_g
reconstruction in BOTH directions across gDNA/RNA FL distributions + enrichment invariance. Run in the rigel env.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import fl_mean, region_eff_length


def make_pmf(kind: str, m: float, s: float, L: int = 800) -> np.ndarray:
    x = np.arange(1, L)
    if kind == "gauss":
        p = np.exp(-0.5 * ((x - m) / s) ** 2)
    elif kind == "gamma":
        p = x ** ((m / s) ** 2 - 1) * np.exp(-x / (s * s / m))
    elif kind == "unif":
        p = ((x >= m - s) & (x <= m + s)).astype(float)
    else:
        raise ValueError(kind)
    p = np.maximum(p, 0.0)
    pmf = np.zeros(L)
    pmf[1:] = p / p.sum()
    return pmf


def _dep(pmf, rho, Le, Nbase, rng):
    """(contained-in-E, crossing/spliced-at-Le) counts for a component of density rho, start ~ U[0, Le]."""
    n = int(Nbase * rho)
    s = rng.uniform(0, Le, n)
    ell = rng.choice(len(pmf), size=n, p=pmf).astype(float)
    end = s + ell
    return float((end <= Le).sum()), float((end > Le).sum())


def mc(gpmf, rpmf, Le=2500, rho_g=0.4, rho_n=0.3, rho_m=0.6, e=1.0, Nbase=8_000_000, seed=1):
    rng = np.random.default_rng(seed)
    gc, gx = _dep(gpmf, rho_g, Le, Nbase, rng)     # gDNA: contained, crossing
    nc, nx = _dep(rpmf, rho_n, Le, Nbase, rng)     # nascent RNA: contained, crossing (unspliced)
    mc_, ms = _dep(rpmf, rho_m, Le, Nbase, rng)    # mature RNA: contained-in-exon, SPLICED
    gc, gx, nc, nx, mc_, ms = (v * e for v in (gc, gx, nc, nx, mc_, ms))  # enrichment scales all -> cancels
    fg_exon = gc / (gc + nc + mc_)                 # exon: gDNA / (gDNA + nascent + mature), contained unspliced
    fg_bnd = gx / (gx + nx)                         # boundary UNSPLICED: gDNA / (gDNA + nascent) -- NO mature
    return dict(gc=gc, gx=gx, nc=nc, nx=nx, mc=mc_, ms=ms, fg_exon=fg_exon, fg_bnd=fg_bnd)


def eff(gpmf, rpmf, Le):
    return dict(
        EgC=region_eff_length(np.array([float(Le)]), gpmf)[0],
        ErC=region_eff_length(np.array([float(Le)]), rpmf)[0],
        EgX=fl_mean(gpmf), ErX=fl_mean(rpmf), Espl=fl_mean(rpmf),
    )


def main():
    cases = [
        ("gDNA=RNA gauss200", make_pmf("gauss", 200, 40), make_pmf("gauss", 200, 40)),
        ("gDNA300 / RNA150", make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30)),
        ("gDNA120 / RNA350", make_pmf("gauss", 120, 25), make_pmf("gauss", 350, 60)),
        ("gamma300 / gamma180", make_pmf("gamma", 300, 90), make_pmf("gamma", 180, 60)),
        ("unif250 / unif150", make_pmf("unif", 250, 100), make_pmf("unif", 150, 60)),
    ]
    Le = 2500
    print("EXON <-> intron-exon BOUNDARY (rho_g=0.4, rho_n=0.3, rho_m=0.6). All ratios must be ~1.000:")
    print(f"{'FL dist':22s} {'ρg X/C':>7s} {'ρn X/C':>7s} {'ρm spl/C':>9s} {'recon/ρn':>9s} "
          f"{'B->E fg':>16s} {'E->B fg':>16s}")
    for name, g, r in cases:
        d = mc(g, r, Le=Le); E = eff(g, r, Le)
        # frame-consistent densities (each recovered two ways; ratio ~1 confirms the eff-length shift)
        rg_C, rg_X = d["gc"] / E["EgC"], d["gx"] / E["EgX"]
        rn_C, rn_X = d["nc"] / E["ErC"], d["nx"] / E["ErX"]
        rm_C, rm_S = d["mc"] / E["ErC"], d["ms"] / E["Espl"]
        rho_RNA_exon = (d["nc"] + d["mc"]) / E["ErC"]        # exon RNA density = nascent + mature
        recon = rho_RNA_exon - rm_S                          # subtract spliced-mature -> should equal ρ_n
        # (B->E) add mature: f_g_exon from boundary densities + boundary spliced
        Mg = rg_X * E["EgC"]; Mr = (rn_X + rm_S) * E["ErC"]; fgE = Mg / (Mg + Mr)
        # (E->B) subtract mature: f_g_bnd from exon densities - boundary spliced
        mg = rg_C * E["EgX"]; mn = max(rho_RNA_exon - rm_S, 0.0) * E["ErX"]; fgB = mg / (mg + mn)
        print(f"{name:22s} {rg_X/rg_C:7.3f} {rn_X/rn_C:7.3f} {rm_S/rm_C:9.3f} {recon/rn_X:9.3f} "
              f"{d['fg_exon']:7.4f}/{fgE:7.4f} {d['fg_bnd']:7.4f}/{fgB:7.4f}")
    print("\ncolumns: ρ_c recovered two frames (X=crossing/spliced, C=contained) ratio; recon/ρn (reconciliation);")
    print("         B->E and E->B f_g  true/predicted. Per-strand: apply the same to each RNA strand (the")
    print("         junction spliced is ONE strand); gate the ± mature by whether the boundary carries spliced")
    print("         on the exon's strand/side (WITH vs WITHOUT spliced absorption).")
    print("\nenrichment invariance (gDNA300/RNA150, e=1 vs 100):")
    for e in (1.0, 100.0):
        d = mc(make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30), e=e)
        print(f"  e={e:5.0f}: fg_exon={d['fg_exon']:.4f}  fg_bnd={d['fg_bnd']:.4f}")


if __name__ == "__main__":
    main()
