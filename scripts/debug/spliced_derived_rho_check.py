"""Step-0 for the UNSTRANDED enrichment idea: does spliced-derived ρ_g (strand-FREE) track the oracle gDNA?

For each single-strand EXON region R with a flanking clean (intron↔exon) boundary B, estimate ρ_g WITHOUT
using the read strand, from the motif-stranded SPLICED (pure mature RNA) signal:

    ρ_mature   = M_spliced(B, R's motif-strand, R-facing) / E_rna_crossing(B)
    M_mature(R)= ρ_mature · E_rna_contained(R)
    ρ_g(R)     = clip( M_unspliced(R) − M_mature(R), 0 ) / E_gdna_contained(R)     # residual = gDNA(+nascent)

Compares ρ_g_spliced vs the by-origin ORACLE ρ_g and vs the strand-derived ρ_g=f_g_strand·T, over single-strand
exons, on a STRANDED (ss0.99 — cross-check) and an UNSTRANDED (ss0.5 — the target) condition.

    OMP_NUM_THREADS=1 python scripts/debug/spliced_derived_rho_check.py [cond ...]
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
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
DEFAULT = ["gdna_gdna300_ss_0.99_nrna_none_capture_on", "gdna_gdna300_ss_0.50_nrna_none_capture_on"]


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs = CalibrationSubstrate.from_payload(payload, ra)
    bs = BoundarySubstrate.from_payload(payload)
    geom = bp.build_node_geometry(ch, cs, bs, ra, blob["gdna_pmf"], blob["rna_pmf"])
    stat = bp.build_node_statics(ch, cs, bs, ra)
    cap = {}
    orig = bp.node_sweep
    calibrate.__globals__["node_sweep"] = lambda c, s, g, b, rg, bsub, **k: (cap.update(k), orig(c, s, g, b, rg, bsub, **k))[1]
    calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    calibrate.__globals__["node_sweep"] = orig
    kappa = cap["rna_sense_frac"]
    odg = cap.get("gdna_strand_overdispersion", 0.0)
    odr = cap.get("rna_strand_overdispersion", 0.0)

    is_reg = np.asarray(ch.kind) == REGION
    is_bnd = np.asarray(ch.kind) == BOUNDARY
    refidx = np.asarray(ch.ref_idx, np.int64)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    EGl, EGr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    ERl, ERr = np.asarray(geom.eff_rna_left), np.asarray(geom.eff_rna_right)
    Ml = np.asarray(geom.mass_left)
    SPl, SPr = np.asarray(geom.spliced_pos_left), np.asarray(geom.spliced_pos_right)
    SNl, SNr = np.asarray(geom.spliced_neg_left), np.asarray(geom.spliced_neg_right)

    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    gR = np.asarray(csg.contained.mass_unspliced, float)
    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])
    r = np.clip(refidx, 0, gR.shape[0] - 1)
    rho_oracle = gR[r] / np.maximum(EGl, _EPS)
    T = Ml / np.maximum(EGl, _EPS)
    sc, tc = scl[r], tcl[r]

    # clean exon boundaries
    node_rtype, rtype = bp._node_region_type(ch, ra)
    Rn = rtype.shape[0]
    blr = np.asarray(bs.left_region, np.int64)
    brr = np.asarray(bs.right_region, np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(refidx, 0, Bn - 1)
    lt = np.where((blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, Rn - 1)], -1)
    rt = np.where((brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, Rn - 1)], -1)
    clean = is_bnd & ((((lt == 0) | (lt == 1)) & (rt == 2)) | (((rt == 0) | (rt == 1)) & (lt == 2)))

    lat = _simplex_lattice(60)
    fgg = lat[2]

    def strand_fg(i):
        psi = _local_loglik(stat.u_pos[i:i+1], stat.u_neg[i:i+1], stat.spliced_pos[i:i+1], stat.spliced_neg[i:i+1],
                            stat.free_pos[i:i+1], stat.free_neg[i:i+1], kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i+1])
        return float(_fg_median(psi, fgg)[0])

    def spliced_rho_g(i, strand):
        """spliced-derived ρ_g for region i, using its flanking clean boundary on the given motif-strand."""
        SPL, SPR = (SPl, SPr) if strand == 1 else (SNl, SNr)
        # region-facing flanking clean boundary: i is RIGHT region of left[i], LEFT region of right[i]
        best = None
        lb = int(left[i])
        if lb >= 0 and clean[lb]:
            best = (SPR[lb], ERr[lb])    # R faces lb's RIGHT side
        rb = int(right[i])
        if rb >= 0 and clean[rb]:
            cand = (SPL[rb], ERl[rb])    # R faces rb's LEFT side
            if best is None or cand[0] > best[0]:
                best = cand
        if best is None:
            return np.nan
        m_spliced, e_spliced = best
        rho_mature = m_spliced / max(e_spliced, _EPS)
        m_mature_R = rho_mature * ERl[i]
        return max(Ml[i] - m_mature_R, 0.0) / max(EGl[i], _EPS)

    ss = np.where(is_reg & ((sc == 1) | (sc == 2)) & (tc == 2) & (Ml > 0))[0]
    amb = np.where(is_reg & (sc == 3) & (tc == 2) & (Ml > 0))[0]
    rg_spl = np.array([spliced_rho_g(i, int(sc[i])) for i in ss])
    rg_str = np.array([strand_fg(i) * T[i] for i in ss])
    rg_or = rho_oracle[ss]
    ok = np.isfinite(rg_spl) & (rg_or >= 0)

    def corr(a, b, m):
        m = m & np.isfinite(a) & np.isfinite(b)
        return float(np.corrcoef(a[m], b[m])[0, 1]) if m.sum() > 3 else float("nan")

    print(f"\n=== {cond}  (κ={kappa:.3f}, (2κ−1)²={(2*kappa-1)**2:.3f}) ===")
    print(f"  single-strand exons with a clean flanking boundary: {int(ok.sum())}/{ss.size}")
    print(f"  ρ_g_spliced vs ORACLE      corr = {corr(rg_spl, rg_or, ok):.3f}")
    print(f"  ρ_g_strand  vs ORACLE      corr = {corr(rg_str, rg_or, ok):.3f}   (the read-strand method)")
    print(f"  ρ_g_spliced vs ρ_g_strand  corr = {corr(rg_spl, rg_str, ok):.3f}")
    en = ok & (rg_or > 1.0)
    if en.sum() > 3:
        print(f"  enriched (oracle ρ_g>1): median ρ_g spliced={np.median(rg_spl[en]):.2f} "
              f"strand={np.median(rg_str[en]):.2f} oracle={np.median(rg_or[en]):.2f}  (n={int(en.sum())})")

    # The real test: fit ê(z) on spliced-derived ρ_g and recover the withheld AMBIG exons (strand-free).
    from rigel.calibration.variance_model import MonotoneMean

    def z_of(i):
        v = []
        lb = int(left[i])
        rb = int(right[i])
        if lb >= 0 and clean[lb]:
            v.append(Mr_face(lb))
        if rb >= 0 and clean[rb]:
            v.append(Ml_face(rb))
        return float(np.mean(v)) if v else np.nan

    def Ml_face(b):
        return Ml[b] / max(EGl[b], _EPS)

    def Mr_face(b):
        return np.asarray(geom.mass_right)[b] / max(EGr[b], _EPS)

    z = np.array([z_of(i) for i in range(ch.n_nodes)])
    fitm = ok & (rg_spl > 0) & np.isfinite(z[ss])
    if fitm.sum() > 20:
        zf = z[ss][fitm]
        yf = rg_spl[fitm]
        ef = EGl[ss][fitm]
        ehat = MonotoneMean.fit(zf, yf, weight=ef, recal_weight=ef)
        amb_ok = amb[np.isfinite(z[amb])]            # AMBIG exons with a usable z (clean flanking boundary)
        pred = ehat.predict(z[amb_ok])
        Tamb = T[amb_ok]
        Mamb = Ml[amb_ok]
        fg = np.clip(pred / np.maximum(Tamb, _EPS), 0, 1)
        trueg = gR[r[amb_ok]]
        rR = np.asarray(CalibrationSubstrate.from_payload(blob['payload_rna'], ra).contained.mass_unspliced, float)
        net = (fg * Mamb).sum() - trueg.sum()
        oamb_fg = trueg / np.maximum(trueg + rR[r[amb_ok]], _EPS)
        cc = float(np.corrcoef(fg, oamb_fg)[0, 1])
        print(f"  → ê(z) fit on SPLICED-derived ρ_g, applied to {amb_ok.size} withheld AMBIG exons (true gDNA {trueg.sum():,.0f}):")
        print(f"      pred {(fg*Mamb).sum():,.0f}  net {net:+,.0f} ({100*net/max(trueg.sum(),1):+.1f}%)  corr(f_g) {cc:.2f}")


def main():
    print("STEP-0: spliced-derived (strand-free) ρ_g vs oracle — the unstranded enrichment signal")
    for c in (sys.argv[1:] or DEFAULT):
        try:
            run(c)
        except Exception as e:  # noqa: BLE001
            print(f"\n=== {c} === FAILED: {e}")


if __name__ == "__main__":
    main()
