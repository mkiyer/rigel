"""Anatomy of the flagship error nodes: the AMBIG (overlapping opposite-strand) exon regions that hold the
~198K under-called gDNA. Are they SHORT (contained-eff-len degenerate, gDNA actually in their boundary
crossings) or LONG (real contained gDNA, error is pure strand-unidentifiability + wrong global mean)?

For each AMBIG-strand exon region (sorted by oracle contained gDNA), prints: length vs gDNA fragment length,
oracle contained gDNA mass / eff / density, the adjacent boundary crossing gDNA mass/density, oracle f_g, and
total contained mass. Summarises how much of the 198K lives in short-vs-long exons.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_exon_anatomy.py [condition]
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
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402

_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs_full = CalibrationSubstrate.from_payload(payload, ra)
    bs_full = BoundarySubstrate.from_payload(payload)
    geom = build_node_geometry(ch, cs_full, bs_full, ra, blob["gdna_pmf"], blob["rna_pmf"])

    is_reg = np.asarray(ch.kind) == REGION
    refidx = np.asarray(ch.ref_idx, np.int64)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    EGl, EGr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)

    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    csr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    bsg = BoundarySubstrate.from_payload(blob["payload_gdna"])
    gR = np.asarray(csg.contained.mass_unspliced, float)
    rR = np.asarray(csr.contained.mass_unspliced, float)
    gBl, gBr = np.asarray(bsg.left.mass_unspliced, float), np.asarray(bsg.right.mass_unspliced, float)
    Bn = gBl.shape[0]
    tot_full = np.asarray(cs_full.contained.mass_unspliced, float)

    sig = np.asarray(ra.signature)
    sclass = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tclass = np.array([coarse_type_from_signature(int(s)) for s in sig])
    rlen = np.asarray(ra.end, float) - np.asarray(ra.start, float)
    # mean gDNA fragment length for the "short" threshold
    pmf = np.asarray(blob["gdna_pmf"], float)
    flen = float(np.sum(np.arange(pmf.size) * pmf) / max(pmf.sum(), _EPS))

    # AMBIG exon REGION chain nodes
    amb = np.where(is_reg & (sclass[np.clip(refidx, 0, sclass.shape[0] - 1)] == 3)
                   & (tclass[np.clip(refidx, 0, tclass.shape[0] - 1)] == 2))[0]
    r = refidx[amb]
    g_cont = gR[np.clip(r, 0, gR.shape[0] - 1)]
    order = amb[np.argsort(-g_cont)]

    print(f"=== {COND}: AMBIG-exon anatomy (gDNA frag length ≈ {flen:.0f} bp) ===")
    print(f"AMBIG-exon region nodes: {amb.size}  total oracle contained gDNA: {g_cont.sum():,.0f}\n")
    print(f"{'node':>6} {'len':>5} {'short?':>6} {'M_g_cont':>9} {'E_cont':>7} {'ρ_g_cont':>9} "
          f"{'M_cross_g':>10} {'oracle_fg':>9} {'totmass':>8}")
    print("-" * 88)
    short_g = 0.0
    long_g = 0.0
    for n in order[:20]:
        rr = int(refidx[n])
        mc, ec = gR[min(rr, gR.shape[0] - 1)], EGl[n]
        rc = mc / max(ec, _EPS)
        mcross = 0.0
        for nb in (left[n], right[n]):
            if nb < 0:
                continue
            bi = min(int(refidx[nb]), Bn - 1)
            mcross += gBl[bi] + gBr[bi]
        ofg = mc / max(mc + rR[min(rr, rR.shape[0] - 1)], _EPS)
        sh = "SHORT" if rlen[rr] < flen else "long"
        print(f"{n:>6} {rlen[rr]:>5.0f} {sh:>6} {mc:>9.0f} {ec:>7.1f} {rc:>9.2f} "
              f"{mcross:>10.0f} {ofg:>9.3f} {tot_full[min(rr,tot_full.shape[0]-1)]:>8.0f}")
    # global short vs long split of the oracle gDNA
    is_short = rlen[r] < flen
    short_g = g_cont[is_short].sum()
    long_g = g_cont[~is_short].sum()
    print(f"\noracle contained gDNA in SHORT (<{flen:.0f}bp) AMBIG exons: {short_g:,.0f} "
          f"({100*short_g/max(g_cont.sum(),_EPS):.0f}%)")
    print(f"oracle contained gDNA in LONG  AMBIG exons:        {long_g:,.0f} "
          f"({100*long_g/max(g_cont.sum(),_EPS):.0f}%)")
    # for the contained density: is E_cont tiny (degenerate) on the mass-bearing AMBIG exons?
    ec_all = EGl[amb]
    big = g_cont > np.median(g_cont[g_cont > 0]) if (g_cont > 0).any() else np.zeros(amb.size, bool)
    print(f"\nmedian E_contained on AMBIG exons holding >median gDNA: {np.median(ec_all[big]):.1f} bp "
          f"(if ≪ region length ⇒ contained density is degenerate/inflated)")
    print(f"median region length of those: {np.median(rlen[r][big]):.0f} bp")


if __name__ == "__main__":
    main()
