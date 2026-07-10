"""Ablation toward a SIMPLE, UNIFIED node-pair imputation model that yields BOTH a density prediction and a
count-space precision from ONE fit. No production code; flagship capture-ON (off-capture invariance already
validated in enrichment_phase0.py).

The model: ρ_dst = ê(ρ_src), a monotone node-pair regression (source = flanking gDNA-clean boundary, dest =
region interior). Two estimator families compared:
  - loglog_recal : weighted log-log OLS + a count-recalibration scale (the Phase-0 stand-in).
  - poisson_glm  : a Poisson/NB GLM predicting the dest gDNA COUNT (log link, offset log E_dst, term log ρ_src).
                   Net-unbiased BY CONSTRUCTION (score eqn ⇒ Σμ=Σy, no recal hack); the dispersion φ IS the
                   biological overdispersion ⇒ mean AND variance from one parametric fit.
Three predictor variants (use as many node-pairs as possible):
  - max   : louder of the two flanking clean boundaries (Phase-0).
  - mean  : average of the (≤2) flanking clean boundary densities.
  - pairs : each flanking clean boundary is its OWN training/prediction pair (≈2× data; predict = average).

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_ablation.py
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
from rigel.calibration.simplex_sweep import _fg_median, _fg_var, _local_loglik, _simplex_lattice  # noqa: E402
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate  # noqa: E402
from rigel.calibration.variance_model import MonotoneMean  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"


def poisson_glm(logx, y, logE, w, iters=40):
    """Quasi-Poisson IRLS: log E[y] = b0 + b1*logx + logE(offset). Returns predict(ρ from logx), (b0,b1), φ, r.
    Net-unbiased: the b0 score equation gives Σ w·μ = Σ w·y on the fit set. φ = Pearson dispersion (>1 ⇒
    biological overdispersion); NB r = MoM size."""
    X = np.vstack([np.ones_like(logx), logx]).T
    b = np.array([np.log(max(np.average(y, weights=w), _EPS)) - np.average(logE, weights=w), 0.0])
    for _ in range(iters):
        eta = X @ b + logE
        mu = np.exp(np.clip(eta, -30, 30))
        W = w * mu
        working = X @ b + (y - mu) / np.maximum(mu, _EPS)   # IRLS working response (offset held fixed)
        WX = X * W[:, None]
        b_new = np.linalg.solve(X.T @ WX + 1e-9 * np.eye(2), X.T @ (W * working))
        if np.max(np.abs(b_new - b)) < 1e-8:
            b = b_new
            break
        b = b_new
    mu = np.exp(np.clip(X @ b + logE, -30, 30))
    dof = max(w.sum() - 2, 1.0)
    phi = float(np.sum(w * (y - mu) ** 2 / np.maximum(mu, _EPS)) / dof)   # quasi-Poisson dispersion
    # NB size r from MoM: Var = μ + μ²/r  ⇒ r = Σwμ² / Σ[w((y-μ)² - μ)]  (clip to positive)
    num = float(np.sum(w * mu ** 2))
    den = float(np.sum(w * ((y - mu) ** 2 - mu)))
    r = num / den if den > 0 else np.inf
    return (lambda lz: np.exp(np.clip(b[0] + b[1] * np.asarray(lz), -30, 30))), (float(b[0]), float(b[1])), phi, r


def loglog_recal(logx, y, E, w):
    """Weighted log-log OLS + count-recal scale (Phase-0)."""
    yl = np.log(np.maximum(y, _EPS))
    X = np.vstack([np.ones_like(logx), logx]).T
    W = np.sqrt(w)
    cf, *_ = np.linalg.lstsq(X * W[:, None], yl * W, rcond=None)
    raw = lambda lz: np.exp(cf[0] + cf[1] * np.asarray(lz))  # noqa: E731
    g_obs = y * E
    s = float(np.sum(g_obs) / max(np.sum(raw(logx) * E), _EPS))
    return lambda lz: s * raw(lz)  # noqa: E731


def main():
    index, blob = build_or_load_cache(COND, False)
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
    kappa, odg, odr = cap["rna_sense_frac"], cap.get("gdna_strand_overdispersion", 0.0), cap.get("rna_strand_overdispersion", 0.0)

    is_reg = np.asarray(ch.kind) == REGION
    is_bnd = np.asarray(ch.kind) == BOUNDARY
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
    ofg = gR[r] / np.maximum(gR[r] + rR[r], _EPS)
    sc, tc = scl[r], tcl[r]

    # clean_exon_bnd + region-facing per-side densities
    blr, brr = np.asarray(bs.left_region, np.int64), np.asarray(bs.right_region, np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(refidx, 0, Bn - 1)
    R = tcl.shape[0]
    lt = np.where((blr[bi_] >= 0) & is_bnd, tcl[np.clip(blr[bi_], 0, R - 1)], -1)
    rtt = np.where((brr[bi_] >= 0) & is_bnd, tcl[np.clip(brr[bi_], 0, R - 1)], -1)
    clean = is_bnd & ((((lt == 0) | (lt == 1)) & (rtt == 2)) | (((rtt == 0) | (rtt == 1)) & (lt == 2)))
    d_left, d_right = Ml / np.maximum(EGl, _EPS), Mr / np.maximum(EGr, _EPS)

    def edges(i):  # region-facing clean boundary densities (a list, ≤2)
        v = []
        lb = left[i]
        if lb >= 0 and clean[lb]:
            v.append(d_right[lb])     # i is lb's right region
        rb = right[i]
        if rb >= 0 and clean[rb]:
            v.append(d_left[rb])      # i is rb's left region
        return v

    lat = _simplex_lattice(60)
    fgg = lat[2]

    def sfg(i):
        psi = _local_loglik(stat.u_pos[i:i + 1], stat.u_neg[i:i + 1], stat.spliced_pos[i:i + 1],
                            stat.spliced_neg[i:i + 1], stat.free_pos[i:i + 1], stat.free_neg[i:i + 1],
                            kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i + 1])
        return float(_fg_median(psi, fgg)[0]), float(_fg_var(psi, fgg)[0])

    ss = [i for i in np.where(is_reg & ((sc == 1) | (sc == 2)) & (tc == 2) & (Ml > 0))[0] if edges(i)]
    amb = [i for i in np.where(is_reg & (sc == 3) & (tc == 2) & (Ml > 0))[0] if edges(i)]
    w_strand = (2 * kappa - 1) ** 2

    # precompute solved-node response (realizable) once
    fit_rho, fit_var, fit_E, fit_edges = {}, {}, {}, {}
    for i in ss:
        fg, v = sfg(i)
        fit_rho[i], fit_var[i], fit_E[i], fit_edges[i] = fg * T[i], v, EGl[i], edges(i)
    trueg = np.array([gR[r[i]] for i in amb])

    def build_xy(mode):
        """Return (logx, y=ρ, E, w) training arrays for predictor `mode`."""
        lx, y, E, w = [], [], [], []
        for i in ss:
            es = fit_edges[i]
            wt = w_strand / (1.0 + fit_var[i]) * fit_E[i]
            if mode == "pairs":
                for e in es:
                    lx.append(np.log(max(e, _EPS))); y.append(fit_rho[i]); E.append(fit_E[i]); w.append(wt / len(es))
            else:
                zz = max(es) if mode == "max" else float(np.mean(es))
                lx.append(np.log(max(zz, _EPS))); y.append(fit_rho[i]); E.append(fit_E[i]); w.append(wt)
        return map(np.array, (lx, y, E, w))

    def predict_amb(predict_rho, mode):
        out = []
        for i in amb:
            es = edges(i)
            if mode == "pairs":
                out.append(float(np.mean([predict_rho(np.log(max(e, _EPS))) for e in es])))
            else:
                zz = max(es) if mode == "max" else float(np.mean(es))
                out.append(float(predict_rho(np.log(max(zz, _EPS)))))
        return np.array(out)

    def score(pred_rho):
        Tamb = np.array([T[i] for i in amb])
        Mamb = np.array([Ml[i] for i in amb])
        oamb = np.array([ofg[i] for i in amb])
        fg = np.clip(pred_rho / np.maximum(Tamb, _EPS), 0, 1)
        net = (fg * Mamb).sum() - trueg.sum()
        corr = float(np.corrcoef(fg, oamb)[0, 1])
        return 100 * net / max(trueg.sum(), 1), corr

    print(f"=== {COND}: unified node-pair imputation ablation (true AMBIG gDNA {trueg.sum():,.0f}) ===")
    print(f"solved ss-exon pairs: {len(ss)}   withheld AMBIG: {len(amb)}\n")
    print(f"{'estimator':>14} {'predictor':>8} | {'net %':>8} {'corr':>6} | {'dispersion / σ²_bio':>22}")
    print("-" * 64)
    for mode in ("max", "mean", "pairs"):
        lx, y, E, w = build_xy(mode)
        # loglog_recal
        f_ll = loglog_recal(lx, y, E, w)
        n_ll, c_ll = score(predict_amb(lambda lz: f_ll(lz), mode))
        print(f"{'loglog_recal':>14} {mode:>8} | {n_ll:>+7.1f}% {c_ll:>6.2f} | {'(recal scale; var separate)':>22}")
        # poisson_glm (count-space; mean + dispersion in ONE fit)
        g = y * E
        f_pg, (b0, b1), phi, rnb = poisson_glm(lx, g, np.log(np.maximum(E, _EPS)), w)
        n_pg, c_pg = score(predict_amb(lambda lz: f_pg(lz), mode))
        s2 = np.log(max(phi, 1.0 + _EPS))   # biological σ²_bio ≈ log(dispersion) (CV²→log-var)
        print(f"{'poisson_glm':>14} {mode:>8} | {n_pg:>+7.1f}% {c_pg:>6.2f} | φ={phi:>5.2f} r={rnb:>6.1f} σ²_bio≈{s2:.2f} (slope {b1:.2f})")
        # MonotoneMean (the PRODUCTION primitive: monotone-spline log-log mean + count-recal scale)
        mm = MonotoneMean.fit(np.exp(lx), y, weight=w, recal_weight=E)
        n_mm, c_mm = score(predict_amb(lambda lz: float(mm.predict(np.array([np.exp(lz)]))[0]), mode))
        print(f"{'MonotoneMean':>14} {mode:>8} | {n_mm:>+7.1f}% {c_mm:>6.2f} | (production spline; var via var~mean)")
    print("\nfor reference, Phase-0 (loglog_recal/max) was net -2.3% corr 0.53; genome-wide global net -51% corr 0.28")

    # --- PRECISION: the residual around the LEARNED mean ê vs around the IDENTITY mean (ρ_src=z). The current
    # var~mean fits (ρ_dst−ρ_src)² — which under capture is inflated by the systematic edge→interior 3× GAP,
    # reading the bias as "variance" ⇒ huge σ²_g ⇒ N≈1 (the original inert-imputation problem). Centering on ê
    # removes the gap ⇒ the residual is the TRUE conditional spread ⇒ small σ²_bio ⇒ usable precision. ---
    lx, y, E, w = build_xy("mean")
    f_ll = loglog_recal(lx, y, E, w)
    z_src = np.exp(lx)
    log_resid_ehat = np.log(np.maximum(y, _EPS)) - np.log(np.maximum(f_ll(lx), _EPS))   # ρ_dst vs ê(ρ_src)
    log_resid_ident = np.log(np.maximum(y, _EPS)) - np.log(np.maximum(z_src, _EPS))      # ρ_dst vs ρ_src
    s2_ehat = float(np.average(log_resid_ehat ** 2, weights=w))    # log-space biological variance (CV² proxy)
    s2_ident = float(np.average(log_resid_ident ** 2, weights=w))
    print("\nPRECISION — σ²_bio must be CONDITIONAL on the level (this is exactly what var~mean provides), NOT a")
    print("global scalar (which depleted log-noisy exons dominate). log-space residual around ê, binned by z:")
    zz = np.exp(lx)
    qs = np.quantile(zz, [0, 0.25, 0.5, 0.75, 1.0])
    for a, b in zip(qs[:-1], qs[1:]):
        m = (zz >= a) & (zz <= b)
        if m.sum() < 4:
            continue
        s2 = float(np.average(log_resid_ehat[m] ** 2, weights=w[m]))
        rho_med = float(np.median(y[m]))
        print(f"  z∈[{a:6.2f},{b:6.2f}]  ρ_g~{rho_med:6.2f}  σ²_bio={s2:5.2f}  N_cap≈{1/max(s2,_EPS):4.1f}")
    enr = zz > np.median(zz)
    s2e = float(np.average(log_resid_ehat[enr] ** 2, weights=w[enr]))
    print(f"  ENRICHED half (the 198K-bearing nodes): σ²_bio={s2e:.2f} ⇒ N_cap≈{1/max(s2e,_EPS):.1f}  "
          "(global scalar {:.2f} is irrelevant — depleted-dominated)".format(s2_ehat))
    print("  ⇒ ê fixes the MEAN (the recovery); var~mean(level) gives per-level precision — usable where gDNA lives.")


if __name__ == "__main__":
    main()
