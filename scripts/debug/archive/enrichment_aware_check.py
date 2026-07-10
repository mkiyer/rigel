"""Empirical backbone for the enrichment-aware framework.

Tests the central claim: under capture the gDNA density is ~uniform WITHIN the captured-exon class, so an
enrichment-aware global that learns ρ_g_captured from the SOLVABLE single-strand exons (no global needed —
strand alone) and applies f_g = ρ_g_captured / T_node (T = node total density) recovers the oracle f_g of the
withheld AMBIG exons. Reports:

  1. ρ_g_captured learned from single-strand LONG exons — via ORACLE f_g (ceiling) and via STRAND-SOLVE f_g
     (the realizable bootstrap estimate, no global). Within-class CV = the prior precision.
  2. For each AMBIG exon: predicted f_g = clip(ρ_g_captured / T, 0, 1) vs oracle f_g; net gDNA mass recovered.
  3. Same for the depleted classes (introns/intergenic) — does the SAME framework leave them alone?

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_aware_check.py [condition]
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
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    from rigel.calibration.substrate import BoundarySubstrate
    cs_full = CalibrationSubstrate.from_payload(payload, ra)
    bs_full = BoundarySubstrate.from_payload(payload)
    geom = bp.build_node_geometry(ch, cs_full, bs_full, ra, blob["gdna_pmf"], blob["rna_pmf"])

    is_reg = np.asarray(ch.kind) == REGION
    refidx = np.asarray(ch.ref_idx, np.int64)
    EGl = np.asarray(geom.eff_gdna_left)
    Ml = np.asarray(geom.mass_left)

    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    csr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    gR = np.asarray(csg.contained.mass_unspliced, float)
    rR = np.asarray(csr.contained.mass_unspliced, float)

    sig = np.asarray(ra.signature)
    sclass = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tclass = np.array([coarse_type_from_signature(int(s)) for s in sig])
    rlen = np.asarray(ra.end, float) - np.asarray(ra.start, float)
    pmf = np.asarray(blob["gdna_pmf"], float)
    flen = float(np.sum(np.arange(pmf.size) * pmf) / max(pmf.sum(), _EPS))

    # per-region-node observables on the contained face
    r = np.clip(refidx, 0, gR.shape[0] - 1)
    T = Ml / np.maximum(EGl, _EPS)                       # total density (observable for ALL nodes)
    rho_g_oracle = gR[r] / np.maximum(EGl, _EPS)         # oracle gDNA density
    oracle_fg = gR[r] / np.maximum(gR[r] + rR[r], _EPS)
    sc = sclass[r]
    tc = tclass[r]
    long_ = rlen[r] >= flen

    # --- 1. learn ρ_g_captured from SOLVABLE single-strand long exons ---
    ss_exon = is_reg & ((sc == 1) | (sc == 2)) & (tc == 2) & long_ & (Ml > 0)
    rho_cap_oracle = float(np.median(rho_g_oracle[ss_exon]))
    cv_oracle = float(np.std(rho_g_oracle[ss_exon]) / max(np.mean(rho_g_oracle[ss_exon]), _EPS))

    # strand-solve f_g for the single-strand exons (no global, no imputation) → realizable bootstrap estimate
    stat = bp.build_node_statics(ch, cs_full, bs_full, ra)
    lat = _simplex_lattice(60)
    fgg = lat[2]
    kappa = blob.get("kappa", 0.99)
    # recover kappa & overdispersions the way calibrate would (approx: use defaults from a quick calibrate)
    from rigel.calibration.calibrate import calibrate
    cap = {}
    orig = bp.node_sweep
    def wrap(chain, s, g, b, rga, bsub, **k):
        cap.update(k)
        return orig(chain, s, g, b, rga, bsub, **k)
    calibrate.__globals__["node_sweep"] = wrap
    calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    calibrate.__globals__["node_sweep"] = orig
    kappa = cap["rna_sense_frac"]; odg = cap.get("gdna_strand_overdispersion", 0.0); odr = cap.get("rna_strand_overdispersion", 0.0)

    def strand_fg(i):
        psi = _local_loglik(stat.u_pos[i:i+1], stat.u_neg[i:i+1], stat.spliced_pos[i:i+1],
                            stat.spliced_neg[i:i+1], stat.free_pos[i:i+1], stat.free_neg[i:i+1],
                            kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i+1])
        return float(_fg_median(psi, fgg)[0])

    ss_idx = np.where(ss_exon)[0]
    rho_cap_strand_vals = np.array([strand_fg(i) * T[i] for i in ss_idx])
    rho_cap_strand = float(np.median(rho_cap_strand_vals))
    cv_strand = float(np.std(rho_cap_strand_vals) / max(np.mean(rho_cap_strand_vals), _EPS))

    print(f"=== {COND}: enrichment-aware framework empirical check (frag len ≈ {flen:.0f}) ===")
    print(f"single-strand LONG exons (solvable bootstrap set): n={ss_exon.sum()}")
    print(f"  ρ_g_captured  ORACLE-f_g  : median={rho_cap_oracle:6.2f}  within-class CV={cv_oracle:.2f}")
    print(f"  ρ_g_captured  STRAND-solve: median={rho_cap_strand:6.2f}  within-class CV={cv_strand:.2f}  "
          f"(realizable, NO global)")
    print(f"  (genome-wide ρ_global for comparison ≈ 0.46 → 33× too low)\n")

    # --- 2. predict AMBIG exons ---
    def report(mask, label, rho_cap):
        idx = np.where(mask)[0]
        if not idx.size:
            return
        Ti = T[idx]
        pred = np.clip(rho_cap / np.maximum(Ti, _EPS), 0.0, 1.0)
        orc = oracle_fg[idx]
        # gDNA mass: contained gDNA = f_g * total_contained_mass
        tot_mass = Ml[idx]
        pred_g = pred * tot_mass
        true_g = gR[r][idx]
        net = pred_g.sum() - true_g.sum()
        finite = np.isfinite(orc)
        corr = float(np.corrcoef(pred[finite], orc[finite])[0, 1]) if finite.sum() > 2 else float("nan")
        mae = float(np.mean(np.abs(pred[finite] - orc[finite])))
        print(f"{label:>26}: n={idx.size:>4}  true_gDNA={true_g.sum():>10,.0f}  pred_gDNA={pred_g.sum():>10,.0f}  "
              f"net={net:>+10,.0f}  f_g MAE={mae:.3f}  corr={corr:.3f}")

    print("predicted f_g = clip(ρ_g_captured / T_node, 0, 1)   [enrichment-aware global, captured baseline]")
    print("baseline for comparison: production assigns these ~f_g≈0 (genome-wide global 33× too low)\n")
    report(is_reg & (sc == 3) & (tc == 2), "AMBIG exon (withheld)", rho_cap_strand)
    report(ss_exon, "SS exon (self, sanity)", rho_cap_strand)
    # depleted classes: the SAME captured baseline must NOT wreck them — but they'd use the DEPLETED baseline.
    print("\nNB: depleted classes (introns/intergenic) would use the DEPLETED baseline, not ρ_g_captured;")
    print("    shown here only to confirm the captured baseline is wrong for them (framework must stratify):")
    report(is_reg & (tc == 1) & (sc != 3), "intron (depleted class)", rho_cap_strand)
    report(is_reg & (tc == 0), "intergenic (depleted class)", rho_cap_strand)


if __name__ == "__main__":
    main()
