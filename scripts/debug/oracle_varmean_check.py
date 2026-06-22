"""Is the fitted σ²_g (which kills imputation under capture) REAL biological spread or an ARTIFACT?

The fitted gDNA var~mean has CV~2 under capture ⇒ N_src≈1/CV²≈1 ⇒ imputation inert. This computes the
ORACLE adjacent-node gDNA-density spread (from the by-origin gDNA-only payload — no RNA contamination, no
strand-deconv error) on the SAME node faces/eff-lengths the seeds use, and compares its var~mean to the
fitted one. Splits seed edges into structural↔structural (clean M/E) vs strand-involving (f_g·M/E, noisy).

  oracle σ²_g ≈ fitted  ⇒ the spread is REAL (adjacent gDNA densities genuinely differ ~2× under capture;
                          imputation legitimately can't rescue the AMBIG exons — fallback-to-global is right)
  oracle σ²_g ≪ fitted  ⇒ ARTIFACT (RNA contamination / strand-deconv noise / region-boundary eff-len
                          convention inflates the fit); fixing it unlocks imputation for the flagship

    OMP_NUM_THREADS=1 python scripts/debug/oracle_varmean_check.py [condition]
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
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.calibration.variance_model import MonotoneVarMean  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
_CAP = {}


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

    # capture the exact seed arrays _fit_seed_varmean is trained on
    orig = bp._fit_seed_varmean

    def cap(chain, dens, eff, is_seed, seed_w):
        _CAP.update(chain=chain, dens=np.asarray(dens, float), eff=np.asarray(eff, float),
                    is_seed=np.asarray(is_seed, bool), seed_w=np.asarray(seed_w, float))
        return orig(chain, dens, eff, is_seed, seed_w)

    bp._fit_seed_varmean = cap
    calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    bp._fit_seed_varmean = orig

    ch = _CAP["chain"]
    dens, eff, is_seed, seed_w = _CAP["dens"], _CAP["eff"], _CAP["is_seed"], _CAP["seed_w"]
    is_reg = np.asarray(ch.kind) == REGION
    is_bnd = np.asarray(ch.kind) == BOUNDARY
    refidx = np.asarray(ch.ref_idx, np.int64)

    # ORACLE gDNA mass per node face, matching the seed's representative face.
    payload = blob["payload_full"]
    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    bsg = BoundarySubstrate.from_payload(blob["payload_gdna"])
    gR = np.asarray(csg.contained.mass_unspliced, float)
    gBl = np.asarray(bsg.left.mass_unspliced, float)
    gBr = np.asarray(bsg.right.mass_unspliced, float)
    # re-derive the seed's exon-on-right choice exactly as _gdna_seed_estimate: representative face uses
    # left mass for regions; for boundaries the side picked depends on exon_on_right. We don't have that mask
    # here, so reconstruct gDNA mass on the SAME face by matching the seed total-mass face value.
    chain2 = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    assert chain2.n_nodes == ch.n_nodes
    # total mass per face from the production geometry (left for region; for boundary pick whichever side's
    # eff matches the seed eff — robust to the exon_on_right choice):
    from rigel.calibration.bp_solver import build_node_geometry
    cs_full = CalibrationSubstrate.from_payload(payload, ra)
    bs_full = BoundarySubstrate.from_payload(payload)
    geom = build_node_geometry(ch, cs_full, bs_full, ra, blob["gdna_pmf"], blob["rna_pmf"])
    EGl, EGr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    # which face did the seed use? the one whose eff equals `eff` (within tol).
    use_right = is_bnd & (np.abs(eff - EGr) < np.abs(eff - EGl))
    Bn = gBl.shape[0]
    bidx = np.clip(refidx, 0, Bn - 1)
    g_mass = np.where(is_reg, gR[np.clip(refidx, 0, gR.shape[0] - 1)],
                      np.where(use_right, gBr[bidx], gBl[bidx]))
    oracle_dens = g_mass / np.maximum(eff, _EPS)

    struct = is_seed & (seed_w >= 0.999)   # weight 1 ⇒ structural seed (clean M/E)
    strand = is_seed & (seed_w < 0.999) & (seed_w > 0)  # (2κ-1)² weight ⇒ strand-deconv seed (f_g·M/E)

    def edge_varmean(node_dens, mask):
        left, right = np.asarray(ch.left), np.asarray(ch.right)
        means, raws, offs, both_struct, wts = [], [], [], [], []
        for nbr in (left, right):
            idx = np.where((nbr >= 0) & mask)[0]
            s = nbr[idx]
            keep = mask[s]
            idx, s = idx[keep], s[keep]
            dr, sr, de, se = node_dens[idx], node_dens[s], eff[idx], eff[s]
            ok = (dr > 0) & (sr > 0)
            means.append(0.5 * (dr[ok] + sr[ok]))
            raws.append((dr[ok] - sr[ok]) ** 2)
            offs.append(dr[ok] / de[ok] + sr[ok] / se[ok])
            both_struct.append(struct[idx][ok] & struct[s][ok])
            wts.append(np.minimum(seed_w[idx][ok], seed_w[s][ok]))
        c = lambda p: np.concatenate(p) if p else np.zeros(0)  # noqa: E731
        return c(means), c(raws), c(offs), c(both_struct), c(wts)

    fm, fr, fo, fbs, fw = edge_varmean(dens, is_seed)          # fitted (seed) edges
    om, orr, oo, obs, ow = edge_varmean(oracle_dens, is_seed)  # oracle gDNA densities, same edges
    fit_fitted = MonotoneVarMean.fit_offset(fm, fr, fo, fw)
    fit_oracle = MonotoneVarMean.fit_offset(om, orr, oo)
    # oracle restricted to structural↔structural edges (cleanest: no strand-deconv noise, no RNA)
    ssm = obs
    fit_oracle_ss = MonotoneVarMean.fit_offset(om[ssm], orr[ssm], oo[ssm]) if ssm.any() else None

    print(f"=== {COND}: is the gDNA σ²_g real or an artifact? ===")
    print(f"seed edges: {fm.size} total, {int(obs.sum())} structural↔structural\n")
    grid = np.array([0.1, 0.5, 1.0, 2.0, 5.0, 10.0])
    hdr = f"{'μ':>6} | {'FITTED σ²_g':>12} {'(CV)':>6} | {'ORACLE σ²_g':>12} {'(CV)':>6} | {'ORACLE-ss σ²_g':>14} {'(CV)':>6}"
    print(hdr)
    print("-" * len(hdr))
    for m in grid:
        ff = float(fit_fitted.predict(np.array([m]))[0])
        oof = float(fit_oracle.predict(np.array([m]))[0])
        oss = float(fit_oracle_ss.predict(np.array([m]))[0]) if fit_oracle_ss else float("nan")
        cv = lambda v: np.sqrt(max(v, 0)) / m  # noqa: E731
        print(f"{m:>6.1f} | {ff:>12.3f} {cv(ff):>6.2f} | {oof:>12.3f} {cv(oof):>6.2f} | {oss:>14.3f} {cv(oss):>6.2f}")
    print("\nimplied N_src = μ²/σ²_g  (at μ where it matters, e.g. μ=2 enriched edge):")
    for m in (1.0, 2.0, 5.0):
        ff = float(fit_fitted.predict(np.array([m]))[0])
        oof = float(fit_oracle.predict(np.array([m]))[0])
        print(f"  μ={m}: N(fitted)={m*m/max(ff,_EPS):>6.1f}  N(oracle)={m*m/max(oof,_EPS):>7.1f}")


if __name__ == "__main__":
    main()
