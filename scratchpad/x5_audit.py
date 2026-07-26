"""ADVERSARIAL AUDIT of scratchpad/x4_estimator.py — independent re-derivation of the headline numbers.

Checks:
  R1  RAW (un-subtracted) Var and mean of log phi under the four claim forms.  Is "V:max = -0.000"
      a real collapse of the raw spread, or an artefact of subtracting a Poisson term that was
      computed with the WRONG (shipped-claim) weights?
  R2  a FORM-APPROPRIATE Poisson for each claim form (weights + counts recomputed per form).
  R3  the w_mu distribution on graft edges — the ledger's omega is charged to the SPLICED
      measurement, and it enters the claim's variance as w_mu^2 * omega.  om1 vs T is only
      apples-to-apples if w_mu ~ 1.
  R4  is TABLE C's `om2_corr == T(2j)` an algebraic identity?  print both sides of the algebra.
  R5  is the share observable tautological with the target?  measure corr(x, log share) with the
      target denominator REPLACED by max(mu) (the hypothesis) vs by the oracle exon density.
  R6  fragment-sharing: does the exon's contained spliced mass Rs overlap the junction's spliced?
  R7  bootstrap CI on V:max and V:this.
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
tri = lambda n: polygamma(1, np.maximum(np.asarray(n, float), 1.0))  # noqa: E731

print("=" * 130)
for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pool(k, isr=isr, idx=idx, inp=inp):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    E_r = us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    SNPf = (us["spl_n_pos_l"], us["spl_n_pos_r"])
    SNNf = (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))

    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rec = []
    for face, nbr in ((1, li), (0, ri)):
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            n_spl = float(SNPf[face][b] + SNNf[face][b])
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not np.isfinite(rho_nu_b[b]):
                continue
            claim = rho_nu_b[b] + mu
            if not (claim > _EPS):
                continue
            rec.append(dict(i=int(i), b=int(b), face=face, mu=float(mu), n_spl=n_spl,
                            rn=float(rho_nu_b[b]), n_nu=float(Ru[b]), D=float(rho_R_ex[i]),
                            n_Rex=float(Ru[i] + Rs[i]), Rs_ex=float(Rs[i]), Ru_ex=float(Ru[i]),
                            S_bnd=float(SPf[face][b] + SNf[face][b]),
                            live=bool(n_spl > 0.0 and mu > _EPS)))
    R = {k: np.asarray([r[k] for r in rec]) for k in rec[0]}
    live = R["live"]

    byex: dict[int, dict] = {}
    for k in np.flatnonzero(live):
        byex.setdefault(int(R["i"][k]), {})[int(R["face"][k])] = k
    pairs = [(v[1], v[0]) for v in byex.values() if 0 in v and 1 in v]
    kl = np.asarray([p[0] for p in pairs], int)
    kr = np.asarray([p[1] for p in pairs], int)
    two = np.concatenate([kl, kr])
    mu_t = R["mu"][two]
    mu_o = np.concatenate([R["mu"][kr], R["mu"][kl]])
    n_t = R["n_spl"][two]
    n_o = np.concatenate([R["n_spl"][kr], R["n_spl"][kl]])
    rn = R["rn"][two]
    D = R["D"][two]
    nnu = R["n_nu"][two]
    nRex = R["n_Rex"][two]

    print(f"\n### {cond[5:]}   pairs={len(pairs)}  2j-edges={two.size}  live={live.sum()}")

    # ── R1/R2: the four claim forms, RAW and with a FORM-APPROPRIATE Poisson ──────────────────
    def pois_for(cl, mu_used, n_used):
        wn = rn / cl
        wm = mu_used / cl
        t_mu = np.where(n_used > 0, wm * wm * tri(np.maximum(n_used, 1.0)), 0.0)
        return wn * wn * tri(nnu) + t_mu + tri(nRex)

    print(f"{'form':<10}{'Var raw':>10}{'mean':>9}{'med':>9}{'E[P] shipped-wts':>18}"
          f"{'E[P] form-wts':>15}{'prem shipped':>14}{'prem form':>11}")
    P_ship = pois_for(rn + mu_t, mu_t, n_t)  # x4's `pshared`: shipped weights, this junction's count
    forms = {
        "this": (rn + mu_t, mu_t, n_t),
        "max": (rn + np.maximum(mu_t, mu_o), np.maximum(mu_t, mu_o),
                np.where(mu_t >= mu_o, n_t, n_o)),
        "sum": (rn + mu_t + mu_o, mu_t + mu_o, n_t + n_o),
        "2x": (rn + 2 * mu_t, 2 * mu_t, n_t),
    }
    for nm, (cl, mu_u, n_u) in forms.items():
        lz = np.log(np.maximum(cl, 1e-30) / D)
        Pf = pois_for(cl, mu_u, n_u)
        print(f"{nm:<10}{np.var(lz):>10.4f}{np.mean(lz):>9.4f}{np.median(lz):>9.4f}"
              f"{np.mean(P_ship):>18.4f}{np.mean(Pf):>15.4f}"
              f"{np.var(lz) - np.mean(P_ship):>14.4f}{np.var(lz) - np.mean(Pf):>11.4f}")

    # R7 bootstrap on V:max premise (form weights) and V:this
    rng = np.random.default_rng(0)
    bs = {"this": [], "max": []}
    for _ in range(400):
        s = rng.integers(0, len(pairs), len(pairs))
        sel = np.concatenate([kl[s], kr[s]])
        mt = R["mu"][sel]
        mo = np.concatenate([R["mu"][kr[s]], R["mu"][kl[s]]])
        nt2 = R["n_spl"][sel]
        no2 = np.concatenate([R["n_spl"][kr[s]], R["n_spl"][kl[s]]])
        rn2, D2 = R["rn"][sel], R["D"][sel]
        nnu2, nR2 = R["n_nu"][sel], R["n_Rex"][sel]

        def _p(cl, mu_u, n_u, rn2=rn2, nnu2=nnu2, nR2=nR2):
            wn, wm = rn2 / cl, mu_u / cl
            return (wn * wn * tri(nnu2)
                    + np.where(n_u > 0, wm * wm * tri(np.maximum(n_u, 1.0)), 0.0) + tri(nR2))
        for nm, (cl, mu_u, n_u) in (("this", (rn2 + mt, mt, nt2)),
                                    ("max", (rn2 + np.maximum(mt, mo), np.maximum(mt, mo),
                                             np.where(mt >= mo, nt2, no2)))):
            lz = np.log(np.maximum(cl, 1e-30) / D2)
            bs[nm].append(float(np.var(lz) - np.mean(_p(cl, mu_u, n_u))))
    for nm in ("this", "max"):
        a = np.asarray(bs[nm])
        print(f"  bootstrap prem[{nm}]  median {np.median(a):.4f}   95% CI "
              f"[{np.percentile(a, 2.5):.4f}, {np.percentile(a, 97.5):.4f}]")

    # ── R3: w_mu ─────────────────────────────────────────────────────────────────────────────
    wm_live = R["mu"][live] / (R["rn"][live] + R["mu"][live])
    wm2 = wm_live ** 2
    print(f"  w_mu on LIVE graft edges: p25 {np.percentile(wm_live, 25):.3f}  med "
          f"{np.median(wm_live):.3f}  p75 {np.percentile(wm_live, 75):.3f}   "
          f"mean w_mu^2 = {wm2.mean():.3f}   (E1 charges the SPLICED arm; it enters the claim as w_mu^2)")

    # ── R4: the om2 == T(2j) algebra ─────────────────────────────────────────────────────────
    xl, xr = np.log((R["rn"] + R["mu"])[kl] / R["D"][kl]), np.log((R["rn"] + R["mu"])[kr] / R["D"][kr])
    Vl, Vr = float(np.var(xl)), float(np.var(xr))
    Cv = float(np.cov(xl, xr, ddof=0)[0, 1])
    sh = float(np.mean(tri(R["n_Rex"][kl])))
    np_ = len(pairs)
    side_all = P_ship - tri(nRex)
    sidel, sider = side_all[:np_], side_all[np_:]
    Pl, Pr, Pc = Vl - sidel.mean() - sh, Vr - sider.mean() - sh, Cv - sh
    om2 = (Vl + Vr - 2 * Cv - sidel.mean() - sider.mean()) / 2
    print(f"  R4 algebra: Pl {Pl:.4f} Pr {Pr:.4f} Pc {Pc:.4f} | om2=(Pl+Pr)/2-Pc = "
          f"{(Pl + Pr) / 2 - Pc:.4f} vs direct {om2:.4f} ; (Pl+Pr)/2 = {(Pl + Pr) / 2:.4f} "
          f"=> om2/T = {om2 / ((Pl + Pr) / 2):.3f} is FORCED by Pc={Pc:.4f}")

    # ── R5: is log(share) tautological with the target? ──────────────────────────────────────
    share = mu_t / np.maximum(mu_t + mu_o, _EPS)
    x = np.log((rn + mu_t) / D)
    x_hyp = np.log((rn + mu_t) / (rn + np.maximum(mu_t, mu_o)))  # target replaced by the HYPOTHESIS
    print(f"  R5 corr(x, log share) = {np.corrcoef(x, np.log(share))[0, 1]:.3f} ; "
          f"corr(x_hypothetical, log share) = {np.corrcoef(x_hyp, np.log(share))[0, 1]:.3f}"
          f"  (the 2nd uses NO oracle exon density at all)")

    # ── R6: fragment sharing between the exon's contained spliced and the junctions' spliced ─
    exs = np.asarray([byex[k] for k in byex])  # noqa: F841
    ex_ids = np.asarray(list(byex.keys()))
    Rs_ex = R["Rs_ex"][kl]
    S_pair = R["S_bnd"][kl] + R["S_bnd"][kr]
    frac = Rs_ex / np.maximum(R["n_Rex"][kl], _EPS)
    print(f"  R6 exon CONTAINED spliced mass Rs/(Ru+Rs): med {np.median(frac):.4f} "
          f"p90 {np.percentile(frac, 90):.4f} ; median Rs_ex {np.median(Rs_ex):.1f} vs median "
          f"junction spliced mass sum {np.median(S_pair):.1f}  (n_exons {ex_ids.size})")
