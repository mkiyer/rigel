"""Is the enriched-exon ↔ boundary gDNA-density cliff REAL (sharp capture) or an eff-len artifact?

For the highest-density exon regions, print the oracle gDNA mass / eff-len / density on the exon (contained)
and on its adjacent clean-exon boundary (crossing). Two outcomes:

  M_crossing tiny vs M_contained, E comparable  ⇒ REAL: few gDNA fragments cross the exon edge (capture is
      spatially sharp; the flank is genuinely depleted) ⇒ neighbours can't inform the exon ⇒ stratified global.
  M_crossing substantial but E_crossing ≫ E_contained ⇒ ARTIFACT: the crossing eff-len overstates support and
      buries a real signal ⇒ reconciling the eff-len convention would let imputation carry the exon.

    OMP_NUM_THREADS=1 python scripts/debug/boundary_efflen_check.py [condition]
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
    bsg = BoundarySubstrate.from_payload(blob["payload_gdna"])
    gR = np.asarray(csg.contained.mass_unspliced, float)
    gBl, gBr = np.asarray(bsg.left.mass_unspliced, float), np.asarray(bsg.right.mass_unspliced, float)
    Bn = gBl.shape[0]

    sig = np.asarray(ra.signature)
    rtype = np.array([coarse_type_from_signature(int(s)) for s in sig])
    # region length (bp) to judge vs fragment length
    rstart = np.asarray(ra.start, float)
    rend = np.asarray(ra.end, float)
    rlen = rend - rstart

    # exon-region nodes with oracle contained density
    ex_nodes = np.where(is_reg & (rtype[np.clip(refidx, 0, rtype.shape[0] - 1)] == 2))[0]
    dens_ex = gR[np.clip(refidx[ex_nodes], 0, gR.shape[0] - 1)] / np.maximum(EGl[ex_nodes], _EPS)
    order = ex_nodes[np.argsort(-dens_ex)][:12]

    print(f"=== {COND}: enriched-exon ↔ boundary gDNA mass/eff/density ===")
    print("exon (contained) | its boundary neighbours (crossing, gDNA side)\n")
    print(f"{'exonNode':>8} {'len':>5} {'M_cont':>8} {'E_cont':>7} {'ρ_cont':>7} || "
          f"{'bndNode':>7} {'M_cross':>8} {'E_cross':>7} {'ρ_cross':>7}")
    print("-" * 88)
    rows = []
    for n in order:
        r = refidx[n]
        mc, ec = gR[r], EGl[n]
        rc = mc / max(ec, _EPS)
        line = f"{n:>8} {rlen[r]:>5.0f} {mc:>8.0f} {ec:>7.1f} {rc:>7.2f} || "
        for nb in (left[n], right[n]):
            if nb < 0:
                continue
            b = int(refidx[nb])
            bi = min(b, Bn - 1)
            # gDNA side of the boundary = the exon-facing side; pick the side with larger gDNA mass
            ml, mr = gBl[bi], gBr[bi]
            if mr >= ml:
                mx, ex_ = mr, EGr[nb]
            else:
                mx, ex_ = ml, EGl[nb]
            rows.append((mc, mx))
            print(line + f"{nb:>7} {mx:>8.1f} {ex_:>7.1f} {mx/max(ex_,_EPS):>7.3f}")
            line = " " * 41 + "|| "
        print()
    if rows:
        mc_all = np.array([r[0] for r in rows])
        mx_all = np.array([r[1] for r in rows])
        print(f"summary over {len(rows)} exon-boundary pairs:")
        print(f"  median M_contained(exon) = {np.median(mc_all):,.0f}")
        print(f"  median M_crossing(bnd)   = {np.median(mx_all):,.1f}   "
              f"(ratio crossing/contained = {np.median(mx_all)/max(np.median(mc_all),_EPS):.4f})")
        print("\n  ⇒ if the crossing MASS itself is ~0 (not just the density), the enrichment is genuinely sharp;")
        print("    the boundary sees no gDNA to impute, so neighbours can't carry the exon (stratified global).")


if __name__ == "__main__":
    main()
