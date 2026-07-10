"""Phase 3b(i): debug the mature-RNA-projection arithmetic in the spliced-derived ρ_g.

The unstranded zero-gDNA + no-nascent + capture condition shows a phantom gDNA across ALL exons (dissect: +152K).
The spliced-derived estimator ρ_g = clip(M_unspliced − ρ_mature·E_rna_contained, 0)/E_gdna over-attributes MATURE
RNA to gDNA. Hypothesis: under capture the junction-crossing density ρ_mature = M_spliced(B)/E_rna_crossing(B)
UNDER-states the intra-exon mature density, so the projection ρ_mature·E_rna_contained under-counts the contained
mature → residual > 0 → phantom.

For each single-strand EXON with a flanking clean boundary, compares PROJECTED mature (ρ_mature·E_rna_contained)
to the TRUE contained mature (oracle; nrna_none ⇒ payload_rna = mature only) and to the contained unspliced M.

    OMP_NUM_THREADS=1 python scripts/debug/mature_projection_check.py [cond ...]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402

_EPS = 1e-9
DEFAULT = ["gdna_none_ss_0.50_nrna_none_capture_on", "gdna_none_ss_0.50_nrna_none_capture_off"]


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs = CalibrationSubstrate.from_payload(payload, ra)
    bs = BoundarySubstrate.from_payload(payload)
    geom = bp.build_node_geometry(ch, cs, bs, ra, blob["gdna_pmf"], blob["rna_pmf"])

    is_reg = np.asarray(ch.kind) == REGION
    is_bnd = np.asarray(ch.kind) == BOUNDARY
    refidx = np.asarray(ch.ref_idx, np.int64)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    ERl, ERr = np.asarray(geom.eff_rna_left), np.asarray(geom.eff_rna_right)
    Ml = np.asarray(geom.mass_left)
    SPl, SPr = np.asarray(geom.spliced_pos_left), np.asarray(geom.spliced_pos_right)
    SNl, SNr = np.asarray(geom.spliced_neg_left), np.asarray(geom.spliced_neg_right)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    node_rtype, rtype = bp._node_region_type(ch, ra)
    Rn = rtype.shape[0]
    blr, brr = np.asarray(bs.left_region, np.int64), np.asarray(bs.right_region, np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(refidx, 0, Bn - 1)
    lt = np.where((blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, Rn - 1)], -1)
    rt = np.where((brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, Rn - 1)], -1)
    clean = is_bnd & ((((lt == 0) | (lt == 1)) & (rt == 2)) | (((rt == 0) | (rt == 1)) & (lt == 2)))

    # true contained mature (nrna_none ⇒ payload_rna = mature only)
    rR = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    r = np.clip(refidx, 0, rR.shape[0] - 1)

    def project(i, strand):
        """returns (rho_mature, projected_contained_mature, e_crossing) using i's flanking clean boundary."""
        SPL, SPR = (SPl, SPr) if strand == 1 else (SNl, SNr)
        best = None
        lb = int(left[i])
        if lb >= 0 and clean[lb]:
            best = (SPR[lb], ERr[lb])
        rb = int(right[i])
        if rb >= 0 and clean[rb]:
            cand = (SPL[rb], ERl[rb])
            if best is None or cand[0] > best[0]:
                best = cand
        if best is None:
            return None
        m_spl, e_spl = best
        rho_mature = m_spl / max(e_spl, _EPS)
        return rho_mature, rho_mature * ERl[i], e_spl

    ss = np.where(is_reg & ((scl[r] == 1) | (scl[r] == 2)) & (tcl[r] == 2) & (Ml > 0))[0]
    tot_M = tot_proj = tot_true = tot_resid = 0.0
    ecross = []
    econt = []
    s_emp = []  # empirical spliced support per region: S = M_spliced / ρ_true = m_spl·E_rna_contained/true_mature
    n_ok = 0
    for i in ss:
        pr = project(i, int(scl[r[i]]))
        if pr is None:
            continue
        rho_mature, m_proj, e_spl = pr
        true_mat = rR[r[i]]
        m_spl = rho_mature * e_spl  # the deposited spliced mass at the flanking boundary
        tot_M += Ml[i]
        tot_proj += m_proj
        tot_true += true_mat
        tot_resid += max(Ml[i] - m_proj, 0.0)
        ecross.append(e_spl)
        econt.append(ERl[i])
        if true_mat > 0:
            s_emp.append(m_spl * ERl[i] / true_mat)
        n_ok += 1
    print(f"\n=== {cond} ===  ({n_ok} single-strand exons with a clean flanking boundary)")
    print(f"  Σ contained unspliced M (≈mature here)  = {tot_M:>12,.0f}")
    print(f"  Σ TRUE contained mature (oracle)        = {tot_true:>12,.0f}")
    print(f"  Σ PROJECTED mature  ρ_mature·E_rna_cont  = {tot_proj:>12,.0f}   "
          f"(projected / true = {tot_proj/max(tot_true,1):.3f})")
    print(f"  Σ residual = Σ clip(M − projected, 0)    = {tot_resid:>12,.0f}   ← the spliced-derived gDNA phantom")
    if ecross:
        from rigel.calibration.effective_length import fl_mean
        rna_fl = blob["rna_pmf"]
        print(f"  eff-len USED for ρ_mature: E_rna_crossing(boundary) median = {np.median(ecross):.1f}   "
              f"E_rna_contained(region) median = {np.median(econt):.1f}")
        print(f"  EMPIRICAL spliced support S_emp = M_spliced/ρ_true: median = {np.median(s_emp):.1f}   "
              f"(the eff-len that WOULD give factor-1)")
        print(f"  RNA FL mean = {fl_mean(rna_fl):.1f}   →  S_emp/E_rna_crossing = {np.median(s_emp)/max(np.median(ecross),1):.2f}"
              f"   (the under-projection factor; ~read-scale vs fragment-scale)")


def main():
    print("PHASE 3b(i): mature-RNA-projection arithmetic (projected vs true contained mature)")
    for c in (sys.argv[1:] or DEFAULT):
        try:
            run(c)
        except Exception as e:  # noqa: BLE001
            print(f"\n=== {c} === FAILED: {e}")


if __name__ == "__main__":
    main()
