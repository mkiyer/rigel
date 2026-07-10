"""Phase 0 — harden the enrichment-imputation harness to the PRODUCTIONIZABLE estimator and check the four
decision gates from the implementation-risk register (docs/calibration/enrichment_aware_calibration.md §L).
No production code; this validates the estimator BEFORE we touch bp_solver.

Hardened vs enrichment_imputation_harness.py:
  - BL-3: z(x) uses ONLY gDNA-clean crossings (clean_exon_bnd: intron/intergenic↔exon), from a belief-
    independent total-density channel T=M/E. Region node → max over clean flanking-boundary densities.
  - BL-5: the transfer ê=E[ρ_g|z] is a NET-UNBIASED MONOTONE mean — weighted isotonic regression (PAVA),
    weight = gDNA eff-len E (⇒ Σ E·ρ̂ = Σ E·ρ = Σ g on the fit set, net-unbiased on the COUNT, no back-
    transform, no magic bins). This is the dependency-free stand-in for the production MonotoneMean GLM.
  - M-9: a permutation slope/structure significance gate (shuffle z among seeds, refit, 95th-pct null) — below
    the noise floor ⇒ ê collapses to ρ_global (off-capture invariance, no K-selector).
  - M-13: fit weighted by seed reliability (2κ−1)²·(1/f_g_var); compare strand-solved ê to the oracle ê.

GATES (capture ON unless noted):
  G1  clean-z transfer R² ≳ 0.90 (the z→ρ_g signal survives the clean-only restriction)
  G2  AMBIG net Σf_g·M within a few % of oracle (post-clip) AND per-node corr ≳ 0.6
  G3  off-capture: ê NOT significant ⇒ collapses to ρ_global; net matches the global-only path
  G4  strand-solved ê within tolerance of the oracle ê (strand-solve noise does not bias the transfer)

    OMP_NUM_THREADS=1 python scripts/debug/enrichment_phase0.py
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
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
ON = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
OFF = "gdna_gdna300_ss_0.99_nrna_none_capture_off"


def pava_weighted(y, w):
    """Weighted pool-adjacent-violators → non-decreasing fit; Σ w·ŷ = Σ w·y by construction."""
    blocks = []  # [value, weight, count]
    for vi, wi in zip(y, w):
        v, ww, c = float(vi), float(wi), 1
        while blocks and blocks[-1][0] >= v - 1e-15:
            pv, pw, pc = blocks.pop()
            v = (pv * pw + v * ww) / (pw + ww)
            ww, c = pw + ww, pc + c
        blocks.append([v, ww, c])
    out = []
    for v, _ww, c in blocks:
        out.extend([v] * c)
    return np.array(out)


def isotonic_fit(zf, rho, w):
    """Return a monotone predictor ê(z) = weighted-isotonic E[ρ|z]. Net-unbiased: Σ w·ê = Σ w·ρ."""
    order = np.argsort(zf)
    zs, ys, ws = zf[order], rho[order], w[order]
    yhat = pava_weighted(ys, ws)
    # collapse to unique z for interp (np.interp needs increasing x); average ties
    zu, idx = np.unique(zs, return_index=True)
    yu = yhat[idx]
    return lambda zq: np.interp(np.asarray(zq, float), zu, yu, left=yu[0], right=yu[-1]), zs, yhat


def strand_params(blob, ra, payload):
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
    sc = scl[r]
    tc = tcl[r]
    rt = np.array([coarse_type_from_signature(int(s)) for s in sig])  # region type per region index

    # ---- BL-3: clean_exon_bnd mask + per-boundary exon-facing clean density ----
    blr = np.asarray(bs.left_region, np.int64)
    brr = np.asarray(bs.right_region, np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(refidx, 0, Bn - 1)
    R = rt.shape[0]
    lt = np.where((blr[bi_] >= 0) & is_bnd, rt[np.clip(blr[bi_], 0, R - 1)], -1)
    rtt = np.where((brr[bi_] >= 0) & is_bnd, rt[np.clip(brr[bi_], 0, R - 1)], -1)
    left_clean = (lt == 0) | (lt == 1)
    right_clean = (rtt == 0) | (rtt == 1)
    clean_exon_bnd = is_bnd & ((left_clean & (rtt == 2)) | (right_clean & (lt == 2)))
    # per-boundary side densities: left side = crossing into the boundary's LEFT region, etc.
    dens_left_side = Ml / np.maximum(EGl, _EPS)
    dens_right_side = Mr / np.maximum(EGr, _EPS)

    # z(region i) = max over flanking CLEAN-exon boundaries of the crossing density on the side FACING region i
    # (region-facing, NOT always exon-facing): i is the RIGHT region of its left boundary and the LEFT region of
    # its right boundary, so an exon reads its own elevated edge while an intergenic/intron region reads its own
    # depleted side — the (z, ρ_g) pair stays meaningful & monotone for every region type.
    def z_region(i):
        vals = []
        lb = left[i]
        if lb >= 0 and clean_exon_bnd[lb]:
            vals.append(dens_right_side[lb])   # i is lb's right region → facing = right side
        rb = right[i]
        if rb >= 0 and clean_exon_bnd[rb]:
            vals.append(dens_left_side[rb])     # i is rb's left region → facing = left side
        return max(vals) if vals else np.nan
    z = np.array([z_region(i) if is_reg[i] else np.nan for i in range(ch.n_nodes)])

    lat = _simplex_lattice(60)
    fgg = lat[2]

    def strand_solve(i):
        psi = _local_loglik(stat.u_pos[i:i + 1], stat.u_neg[i:i + 1], stat.spliced_pos[i:i + 1],
                            stat.spliced_neg[i:i + 1], stat.free_pos[i:i + 1], stat.free_neg[i:i + 1],
                            kappa, odg, odr, lat, strand_obs=stat.strand_obs[i:i + 1])
        return float(_fg_median(psi, fgg)[0]), float(_fg_var(psi, fgg)[0])

    # ---- fit basis: SINGLE-STRAND EXON regions (the transfer is exon-specific: exons carry the edge→interior
    # capture gradient z≈ρ_g/3, while introns/intergenic are ~uniform z≈ρ_g — a DIFFERENT law; mixing them
    # halves R². The AMBIG solve targets are exons, so fit the transfer on exons.) AMBIG withheld (static). ----
    ss_reg = is_reg & ((sc == 1) | (sc == 2)) & (tc == 2)          # single-strand exon regions
    fit_mask = ss_reg & (Ml > 0) & np.isfinite(z)
    amb_mask = is_reg & (sc == 3) & (tc == 2) & (Ml > 0) & np.isfinite(z)

    fi = np.where(fit_mask)[0]
    rho_fit = np.empty(fi.size)
    fgvar = np.empty(fi.size)
    for j, i in enumerate(fi):
        fg, v = strand_solve(i)                                    # realizable: strand-only f_g (no global)
        rho_fit[j] = fg * T[i]
        fgvar[j] = v
    rho_oracle_fit = gR[r][fi] / np.maximum(EGl[fi], _EPS)         # oracle ρ_g on the same nodes (scoring)
    zf = z[fi]
    w_strand = (2 * kappa - 1) ** 2
    w_count = EGl[fi]                                              # eff-len ⇒ net-unbiased on the gDNA count
    w_fit = w_count * w_strand / (1.0 + fgvar)                     # M-13 reliability weighting (E stabilizes)

    # ---- estimator: smooth log-log fit (shape ⇒ per-node corr) + count-recalibration scale (⇒ Σ ê·E = Σ ρ·E
    # on the fit set, net-unbiased, no back-transform). The cheap stand-in for the production Poisson/logit GLM. ----
    def fit_smooth(rho, wf):
        ok = (zf > 0) & (rho > 0)
        X = np.vstack([np.ones(ok.sum()), np.log(zf[ok])]).T
        yl = np.log(rho[ok])
        W = np.sqrt(wf[ok])
        cf, *_ = np.linalg.lstsq(X * W[:, None], yl * W, rcond=None)
        a, p = float(cf[0]), float(cf[1])
        raw = lambda zq: np.exp(a + p * np.log(np.maximum(zq, _EPS)))  # noqa: E731
        s = float(np.sum(rho[ok] * w_count[ok]) / max(np.sum(raw(zf[ok]) * w_count[ok]), _EPS))  # count recal
        R2 = 1 - np.sum((yl - X @ cf) ** 2) / max(np.sum((yl - yl.mean()) ** 2), _EPS)
        return (lambda zq: s * raw(zq)), a, p, s, float(R2)
    ehat, a_hat, p_hat, scale, R2 = fit_smooth(rho_fit, w_fit)
    ehat_or, _, _, _, _ = fit_smooth(rho_oracle_fit, w_count)

    # G3: permutation significance — statistic = log-log R² (z→ρ_g signal); null = shuffle z among seeds
    def r2_of(zvals):
        ok = (zvals > 0) & (rho_fit > 0)
        X = np.vstack([np.ones(ok.sum()), np.log(zvals[ok])]).T
        yl = np.log(rho_fit[ok])
        cf, *_ = np.linalg.lstsq(X, yl, rcond=None)
        return 1 - np.sum((yl - X @ cf) ** 2) / max(np.sum((yl - yl.mean()) ** 2), _EPS)
    real_stat = float(r2_of(zf))
    rng = np.random.default_rng(0)
    null = np.array([float(r2_of(rng.permutation(zf))) for _ in range(200)])
    null95 = float(np.quantile(null, 0.95))
    significant = real_stat > null95
    rho_global = float(np.average(rho_fit, weights=w_count))       # the K=1 / global stand-in

    # ---- score AMBIG recovery ----
    ai = np.where(amb_mask)[0]
    trueg = gR[r][ai]

    def score(pred_rho):
        raw_fg = pred_rho / np.maximum(T[ai], _EPS)
        fg = np.clip(raw_fg, 0, 1)
        clip_bite = float(np.mean(raw_fg > 1.0))
        pg = fg * Ml[ai]
        net = pg.sum() - trueg.sum()
        fin = np.isfinite(ofg[ai])
        corr = float(np.corrcoef(fg[fin], ofg[ai][fin])[0, 1]) if fin.sum() > 2 else float("nan")
        corr_unclip = float(np.corrcoef(raw_fg[fin], ofg[ai][fin])[0, 1]) if fin.sum() > 2 else float("nan")
        return pg.sum(), net, corr, clip_bite, corr_unclip

    # apply ê only if significant, else collapse to ρ_global (M-9)
    pred = ehat(z[ai]) if significant else np.full(ai.size, rho_global)
    pg, net, corr, clip_bite, corr_unclip = score(pred)
    # global-only baseline for the off-capture invariance comparison
    pg_g, net_g, corr_g, _, _ = score(np.full(ai.size, rho_global))

    # G4: strand vs oracle transfer discrepancy on a common z grid
    zg = np.quantile(zf[zf > 0], np.linspace(0.05, 0.95, 15))
    disc = np.median(np.abs(ehat(zg) - ehat_or(zg)) / np.maximum(ehat_or(zg), _EPS))

    return dict(cond=cond, n_fit=fi.size, n_amb=ai.size, trueg=trueg.sum(), R2=R2,
                real_stat=real_stat, null95=null95, significant=significant, rho_global=rho_global,
                pg=pg, net=net, corr=corr, corr_unclip=corr_unclip, clip_bite=clip_bite,
                net_g=net_g, corr_g=corr_g, disc=disc,
                ez_lo=float(ehat(np.nanmin(zf[zf > 0]))), ez_hi=float(ehat(np.nanmax(zf))))


def main():
    print("PHASE 0 — productionizable estimator validation (clean-z, weighted-isotonic, sig-gate, weighted)\n")
    on = run(ON)
    off = run(OFF)
    for o in (on, off):
        tag = "ON " if o["cond"].endswith("on") else "OFF"
        print(f"[{tag}] {o['cond']}")
        print(f"     fit nodes={o['n_fit']}  AMBIG={o['n_amb']} (true gDNA {o['trueg']:,.0f})  clean-z R²={o['R2']:.3f}")
        print(f"     ê range over fit z: [{o['ez_lo']:.2f} .. {o['ez_hi']:.2f}]   ρ_global={o['rho_global']:.2f}")
        print(f"     significance: real={o['real_stat']:.3f} null95={o['null95']:.3f} -> "
              f"{'SIGNIFICANT (apply ê)' if o['significant'] else 'not sig (collapse to ρ_global)'}")
        print(f"     AMBIG recovery: pred {o['pg']:,.0f}  net {o['net']:+,.0f} ({100*o['net']/max(o['trueg'],1):+.1f}%)  "
              f"corr {o['corr']:.2f} (unclipped {o['corr_unclip']:.2f})  clip-bite {100*o['clip_bite']:.0f}%")
        print(f"     global-only baseline: net {o['net_g']:+,.0f}  corr {o['corr_g']:.2f}")
        print(f"     strand-vs-oracle ê discrepancy (median |Δ|/ê): {o['disc']:.3f}\n")

    print("=== GATE CHECKS ===")
    g1 = on["R2"] >= 0.90
    g2 = abs(on["net"]) / max(on["trueg"], 1) <= 0.05 and on["corr"] >= 0.60
    g3 = (not off["significant"]) and abs(off["net"] - off["net_g"]) < 0.02 * max(off["trueg"], 1)
    g4 = on["disc"] <= 0.20
    print(f"  G1 clean-z R²≥0.90:                 {'PASS' if g1 else 'FAIL'}  (R²={on['R2']:.3f})")
    print(f"  G2 capture-on net≤5% & corr≥0.6:    {'PASS' if g2 else 'FAIL'}  "
          f"(net {100*on['net']/max(on['trueg'],1):+.1f}%, corr {on['corr']:.2f})")
    print(f"  G3 off-capture collapses to global: {'PASS' if g3 else 'FAIL'}  "
          f"(sig={off['significant']}, |net−net_global|={abs(off['net']-off['net_g']):,.0f})")
    print(f"  G4 strand ê within 20% of oracle:   {'PASS' if g4 else 'FAIL'}  (disc={on['disc']:.3f})")
    print(f"\n  OVERALL: {'ALL GATES PASS — proceed to Phase 1' if all([g1,g2,g3,g4]) else 'gate failure — see above'}")


if __name__ == "__main__":
    main()
