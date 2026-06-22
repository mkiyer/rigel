"""Under capture, is the true gDNA density predictable WITHIN an enrichment class (→ a stratified baseline /
imputation could solve the AMBIG exons) or locally noisy at CV~2 everywhere (→ only the genome-wide global
remains)?

Computes the ORACLE (by-origin gDNA-only) density per chain node, classifies each node, and reports the
WITHIN-CLASS spread (median, CV) vs the genome-wide ρ_global. Also reports adjacent-edge CV split by whether
the edge stays within a class vs crosses an enrichment boundary (exon-region ↔ exon-edge-boundary).

  within-exon CV ≪ adjacent CV  ⇒ the region↔boundary CV~2 is a real ENRICHMENT CLIFF; exons are mutually
                                   predictable ⇒ an enriched-exon baseline (stratified global) is the lever.
  within-exon CV ≈ 2            ⇒ gDNA is locally noisy even within a class; nothing predicts it; the
                                   genome-wide global is the only honest fallback.

    OMP_NUM_THREADS=1 python scripts/debug/class_density_spread.py [condition]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration.bp_solver import build_node_geometry  # noqa: E402
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_from_signature  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402

_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
TC = {0: "intergenic", 1: "intron", 2: "exon"}


def cv(x, w=None):
    x = np.asarray(x, float)
    if x.size == 0:
        return float("nan"), float("nan"), 0
    if w is None:
        m, s = x.mean(), x.std()
    else:
        w = np.asarray(w, float)
        m = np.sum(w * x) / max(w.sum(), _EPS)
        s = np.sqrt(np.sum(w * (x - m) ** 2) / max(w.sum(), _EPS))
    return m, s / max(m, _EPS), x.size


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs_full = CalibrationSubstrate.from_payload(payload, ra)
    bs_full = BoundarySubstrate.from_payload(payload)
    geom = build_node_geometry(ch, cs_full, bs_full, ra, blob["gdna_pmf"], blob["rna_pmf"])

    is_reg = np.asarray(ch.kind) == REGION
    is_bnd = np.asarray(ch.kind) == BOUNDARY
    refidx = np.asarray(ch.ref_idx, np.int64)
    EGl, EGr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)

    # oracle gDNA mass per node (region contained; boundary = both sides summed for density, eff = avg of sides)
    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    bsg = BoundarySubstrate.from_payload(blob["payload_gdna"])
    gR = np.asarray(csg.contained.mass_unspliced, float)
    gBl, gBr = np.asarray(bsg.left.mass_unspliced, float), np.asarray(bsg.right.mass_unspliced, float)
    Bn = gBl.shape[0]
    bidx = np.clip(refidx, 0, Bn - 1)

    dens = np.zeros(ch.n_nodes)
    eff = np.zeros(ch.n_nodes)
    dens[is_reg] = gR[np.clip(refidx[is_reg], 0, gR.shape[0] - 1)] / np.maximum(EGl[is_reg], _EPS)
    eff[is_reg] = EGl[is_reg]
    # boundary: per-side density, then take the higher-eff side as representative (the one with real crossing)
    bl = gBl[bidx] / np.maximum(EGl, _EPS)
    br = gBr[bidx] / np.maximum(EGr, _EPS)
    use_r = is_bnd & (EGr > EGl)
    dens[is_bnd] = np.where(use_r[is_bnd], br[is_bnd], bl[is_bnd])
    eff[is_bnd] = np.where(use_r[is_bnd], EGr[is_bnd], EGl[is_bnd])
    gcount = dens * eff  # oracle gDNA count on the representative face

    # node class
    sig = np.asarray(ra.signature)
    rtype = np.array([coarse_type_from_signature(int(s)) for s in sig])
    blr = np.asarray(bs_full.left_region, np.int64)
    brr = np.asarray(bs_full.right_region, np.int64)
    R = rtype.shape[0]
    lt = np.where((blr[bidx] >= 0), rtype[np.clip(blr[bidx], 0, R - 1)], -1)
    rt = np.where((brr[bidx] >= 0), rtype[np.clip(brr[bidx], 0, R - 1)], -1)
    cls = np.full(ch.n_nodes, "?", dtype=object)
    cls[is_reg] = ["reg:" + TC.get(int(rtype[r]), "?") for r in refidx[is_reg]]
    # boundary class by its flank pair
    for n in np.where(is_bnd)[0]:
        a, b = int(lt[n]), int(rt[n])
        lo, hi = sorted((a, b))
        cls[n] = "bnd:" + "-".join(TC.get(x, "ext") for x in (lo, hi))

    print(f"=== {COND}: oracle gDNA density spread by node class ===")
    g_global = gcount[(dens > 0)].sum() / np.maximum(eff[(dens > 0)].sum(), _EPS)
    print(f"genome-wide ρ_global (mass-weighted, occupied nodes) = {g_global:.3f}\n")
    print(f"{'class':>22} {'n':>5} {'n_gpos':>6} {'med_dens':>9} {'CV(dens)':>9} {'med/global':>10}")
    print("-" * 66)
    classes = sorted(set(cls.tolist()))
    for c in classes:
        m = (cls == c)
        d = dens[m]
        dp = d[d > 0]
        med = np.median(dp) if dp.size else 0.0
        _, cvd, _ = cv(dp)
        print(f"{c:>22} {m.sum():>5} {dp.size:>6} {med:>9.3f} {cvd:>9.2f} {med/max(g_global,_EPS):>10.2f}")

    # adjacent-edge CV: within-class vs cross-class (oracle densities, sampling-naive squared-diff CV)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    within_r, cross_r = [], []
    within_m, cross_m = [], []
    for nbr in (left, right):
        idx = np.where(nbr >= 0)[0]
        s = nbr[idx]
        ok = (dens[idx] > 0) & (dens[s] > 0)
        idx, s = idx[ok], s[ok]
        same = cls[idx] == cls[s]
        mu = 0.5 * (dens[idx] + dens[s])
        rd = np.abs(dens[idx] - dens[s]) / np.maximum(mu, _EPS)  # relative diff = edge CV proxy
        within_r.append(rd[same]); within_m.append(mu[same])
        cross_r.append(rd[~same]); cross_m.append(mu[~same])
    wr, cr = np.concatenate(within_r), np.concatenate(cross_r)
    print("\nadjacent-edge relative density jump |Δρ|/μ  (both endpoints gDNA-occupied):")
    print(f"  within-class edges:  n={wr.size:>5}  median |Δρ|/μ = {np.median(wr):.2f}")
    print(f"  cross-class edges:   n={cr.size:>5}  median |Δρ|/μ = {np.median(cr):.2f}")
    # enriched-exon focus: how tight is the exon-region class relative to the genome-wide global?
    ex = (cls == "reg:exon") & (dens > 0)
    if ex.any():
        _, cvex, nex = cv(dens[ex])
        print(f"\nENRICHED exon-region nodes: n={nex}  median dens={np.median(dens[ex]):.3f}  "
              f"CV={cvex:.2f}  (vs genome-wide ρ_global {g_global:.3f}: enriched {np.median(dens[ex])/max(g_global,_EPS):.1f}×)")
        print("  ⇒ a stratified enriched-exon baseline predicts an exon to ±%.0f%%, vs the genome-wide global's "
              "wrong center." % (100 * cvex))


if __name__ == "__main__":
    main()
