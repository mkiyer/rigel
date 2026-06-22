"""Validation harness for the IMPUTATION-FIRST, continuous enrichment model (design pre-flight, no prod code).

Tests the proposed second-line-of-defense: solve an AMBIG node's gDNA density by imputing it from its flanking
gDNA-CLEAN boundary crossings, via a CONTINUOUS transfer ê(x)=E[ρ_g | z(x)] fit on the solvable nodes — with
NO global and NO strand for the withheld nodes. Everything realizable (strand-solve for f_g, never the oracle);
oracle used ONLY to score.

Per condition it:
  1. strand-solves f_g for every node (no global), forms ρ_g = f_g·T on the SOLVED set (single-strand exons),
  2. predictor z(x) = max flanking boundary total crossing density (observable, gDNA-clean by structure),
  3. fits the transfer ê = E[ρ_g | z]: log-log linear (interpretable: offset=edge→interior, slope) with
     **Duan smearing** bias-correction, AND isotonic (monotone non-parametric) for comparison,
  4. applies ê to the WITHHELD AMBIG exons → f_g = clip(ρ_g/T) → recovery vs the by-origin oracle,
  5. reports off-capture INVARIANCE (the transfer must collapse to identity: slope≈1, offset≈0, smear≈1, so
     the message reduces to the current identity-mean) and the AMBIG net/corr vs the flat baseline & genome-wide.

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_imputation_harness.py [cond1 cond2 ...]
default conditions: gdna300/ss0.99/nrna_none × capture {on, off}
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
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
DEFAULT = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
]


def strand_params(blob, ra, payload):
    """kappa + strand overdispersions, the way calibrate derives them (no global involved)."""
    cap = {}
    orig = bp.node_sweep

    def wrap(c, s, g, b, rg, bsub, **k):
        cap.update(k)
        return orig(c, s, g, b, rg, bsub, **k)

    calibrate.__globals__["node_sweep"] = wrap
    calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    calibrate.__globals__["node_sweep"] = orig
    return cap["rna_sense_frac"], cap.get("gdna_strand_overdispersion", 0.0), cap.get("rna_strand_overdispersion", 0.0)


def run(cond):
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    cs = CalibrationSubstrate.from_payload(payload, ra)
    bs = BoundarySubstrate.from_payload(payload)
    geom = bp.build_node_geometry(ch, cs, bs, ra, blob["gdna_pmf"], blob["rna_pmf"])
    stat = bp.build_node_statics(ch, cs, bs, ra)
    kappa, odg, odr = strand_params(blob, ra, payload)

    is_reg = np.asarray(ch.kind) == REGION
    refidx = np.asarray(ch.ref_idx, np.int64)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    EGl, EGr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    Ml, Mr = np.asarray(geom.mass_left), np.asarray(geom.mass_right)

    csg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    csr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    gR = np.asarray(csg.contained.mass_unspliced, float)
    rR = np.asarray(csr.contained.mass_unspliced, float)
    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    r = np.clip(refidx, 0, gR.shape[0] - 1)
    T = Ml / np.maximum(EGl, _EPS)
    ofg = gR[r] / np.maximum(gR[r] + rR[r], _EPS)           # oracle f_g (scoring only)
    sc = scl[r]
    tc = tcl[r]
    lat = _simplex_lattice(60)
    fgg = lat[2]

    def strand_fg(i):  # realizable f_g: strand likelihood only (no global, no imputation)
        psi = _local_loglik(stat.u_pos[i:i + 1], stat.u_neg[i:i + 1], stat.spliced_pos[i:i + 1],
                            stat.spliced_neg[i:i + 1], stat.free_pos[i:i + 1], stat.free_neg[i:i + 1],
                            kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i + 1])
        return float(_fg_median(psi, fgg)[0])

    # predictor z(x): max flanking boundary total (unspliced) crossing density — observable for ALL nodes
    def z_of(i):
        v = 0.0
        for nb in (left[i], right[i]):
            if nb < 0:
                continue
            v = max(v, Ml[nb] / max(EGl[nb], _EPS), Mr[nb] / max(EGr[nb], _EPS))
        return v
    z = np.array([z_of(i) if is_reg[i] else 0.0 for i in range(ch.n_nodes)])

    exon = is_reg & (tc == 2)
    ss = np.where(exon & ((sc == 1) | (sc == 2)) & (Ml > 0) & (z > 0))[0]
    amb = np.where(exon & (sc == 3) & (Ml > 0))[0]

    # realizable response on the SOLVED single-strand exons: ρ_g = f_g_strand · T
    rho_solved = np.array([strand_fg(i) for i in ss]) * T[ss]
    keep = rho_solved > 0
    xs, ys = np.log(z[ss][keep]), np.log(rho_solved[keep])

    # ---- fit ê = E[ρ_g | z] : log-log linear + Duan smearing (bias-corrected back-transform) ----
    A = np.vstack([np.ones_like(xs), xs]).T
    coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
    a, p = float(coef[0]), float(coef[1])
    resid = ys - A @ coef
    smear = float(np.mean(np.exp(resid)))                  # Duan smearing factor (corrects E[exp]≠exp E)
    R2 = 1.0 - np.sum(resid ** 2) / max(np.sum((ys - ys.mean()) ** 2), _EPS)

    def ehat_loglin(zq):
        zq = np.maximum(zq, _EPS)
        return smear * np.exp(a + p * np.log(zq))

    # ---- mean-targeting estimator: binned CONDITIONAL MEAN E[ρ_g|z] in LINEAR space (no back-transform
    #      bias), monotonized by pool-adjacent-violators. This is the unbiased-net counterpart. ----
    def pava(v):  # pool adjacent violators → monotone non-decreasing
        v = v.astype(float).copy()
        w = np.ones_like(v)
        i = 0
        while i < len(v) - 1:
            if v[i] > v[i + 1] + 1e-12:
                new = (w[i] * v[i] + w[i + 1] * v[i + 1]) / (w[i] + w[i + 1])
                v[i] = new
                w[i] += w[i + 1]
                v = np.delete(v, i + 1)
                w = np.delete(w, i + 1)
                if i > 0:
                    i -= 1
            else:
                i += 1
        return v, w  # collapsed; needs re-expansion — instead use index map below

    lz_s = np.log(np.maximum(z[ss][keep], _EPS))
    rho_lin = rho_solved[keep]
    nb = 12
    qs = np.unique(np.quantile(lz_s, np.linspace(0, 1, nb + 1)))
    centers, bmeans = [], []
    for i in range(len(qs) - 1):
        sel = (lz_s >= qs[i]) & (lz_s <= qs[i + 1]) if i == len(qs) - 2 else (lz_s >= qs[i]) & (lz_s < qs[i + 1])
        if sel.sum() >= 3:
            centers.append(float(np.mean(lz_s[sel])))
            bmeans.append(float(np.mean(rho_lin[sel])))   # LINEAR-space conditional mean (unbiased)
    centers, bmeans = np.array(centers), np.array(bmeans)
    bmeans_mono = np.maximum.accumulate(bmeans)            # enforce monotone non-decreasing in z
    ehat_binned = lambda zq: np.interp(np.log(np.maximum(zq, _EPS)), centers, bmeans_mono)  # noqa: E731

    # flat captured baseline = median ρ_g over enriched solved exons (the "global" stand-in)
    enr = rho_solved > 1.0
    flat = float(np.median(rho_solved[enr])) if enr.any() else float(np.median(rho_solved))

    trueg = gR[r[amb]]

    def score(pred_rho, lab):
        fg = np.clip(pred_rho / np.maximum(T[amb], _EPS), 0, 1)
        pg = fg * Ml[amb]
        net = pg.sum() - trueg.sum()
        fin = np.isfinite(ofg[amb])
        corr = float(np.corrcoef(fg[fin], ofg[amb][fin])[0, 1]) if fin.sum() > 2 else float("nan")
        mae = float(np.mean(np.abs(fg[fin] - ofg[amb][fin])))
        print(f"    {lab:>34}: pred {pg.sum():>9,.0f} (true {trueg.sum():>9,.0f})  net {net:>+9,.0f}  "
              f"corr {corr:>5.2f}  MAE {mae:.3f}")

    print(f"\n=== {cond} ===")
    print(f"  solved single-strand exons: {ss.size}  |  withheld AMBIG exons: {amb.size}  (true gDNA {trueg.sum():,.0f})")
    print(f"  transfer ê=E[ρ_g|z]  log-log: log ρ_g = {a:+.3f} {p:+.3f}·log z   R²={R2:.3f}  "
          f"smear={smear:.3f}   [identity ⇒ a≈0,p≈1,smear≈1]")
    print(f"    edge→interior factor exp(a)·smear ≈ {np.exp(a)*smear:.2f}×  (capture present ⇔ ≫1; ≈1 ⇒ off-capture)")
    print("  AMBIG recovery (NO global, NO strand — imputation only):")
    score(ehat_binned(z[amb]), "continuous ê (binned cond-mean, monotone)")
    score(ehat_loglin(z[amb]), "continuous ê (log-log + Duan smear)")
    score(np.full(amb.size, flat), f"flat baseline {flat:.1f}/T (third line)")
    score(np.full(amb.size, 0.46), "genome-wide 0.46/T")
    return a, p, smear, R2


def main():
    conds = sys.argv[1:] or DEFAULT
    print("IMPUTATION-FIRST continuous enrichment — validation harness")
    fits = {}
    for c in conds:
        try:
            fits[c] = run(c)
        except Exception as e:  # noqa: BLE001
            print(f"\n=== {c} === FAILED: {e}")
    print("\n--- off-capture INVARIANCE check (transfer must ≈ identity: p≈1, a≈0, smear≈1) ---")
    for c, (a, p, smear, R2) in fits.items():
        tag = "OFF" if "capture_off" in c else "ON "
        ident = abs(p - 1) < 0.25 and abs(a) < 0.5 and abs(smear - 1) < 0.3
        print(f"  [{tag}] {c}: a={a:+.3f} p={p:+.3f} smear={smear:.3f} R²={R2:.3f}  "
              f"{'≈IDENTITY ✓' if ident else 'enrichment transfer (non-identity)'}")


if __name__ == "__main__":
    main()
