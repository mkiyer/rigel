"""AUDIT part 2 — the FRAME the estimator is charged in, over all 14 capture-OFF conditions.

E1 = Var(log(mu_l/mu_r))/2 is a variance of the SPLICED DENSITY.  The target T = Var(log phi) is a
variance of the WHOLE CLAIM rho_nu + mu.  A premise error d in the mu arm arrives in the claim as
w_mu*d (M2's delta method, `graft_rna_logvar`), so the like-for-like comparison is

        E1 * E[w_mu^2]   vs   T          not      E1 vs T.

Prints, per condition: E1, E[w_mu^2] on the 2j edges, the frame-corrected E1, T(2j), both ratios,
and the empirical frame factor Var(g2)/Var(g1) as an independent read of E[w_mu^2].
Also: the FORM-APPROPRIATE Poisson for the max claim form, and the median log phi.
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
CAP_OFF = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna5_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_present_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_off",
]
index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
tri = lambda n: polygamma(1, np.maximum(np.asarray(n, float), 1.0))  # noqa: E731

hdr = (f"{'condition':<32}{'E1(om1)':>9}{'E[wmu2]':>9}{'Var g2/g1':>11}{'E1*wmu2':>9}"
       f"{'om2':>8}{'T(2j)':>8}{'E1/T':>7}{'E1*w/T':>8}{'om2/T':>7}"
       f"{'V:max ship':>11}{'V:max form':>11}{'med max':>9}")
print(hdr)
print("-" * len(hdr))
acc = []
for cond in CAP_OFF:
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
    SNPf, SNNf = (us["spl_n_pos_l"], us["spl_n_pos_r"]), (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    D_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rec = []
    for face, nbr in ((1, li), (0, ri)):
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            ns = float(SNPf[face][b] + SNNf[face][b])
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not np.isfinite(D_ex[i]) or D_ex[i] <= _EPS or not np.isfinite(nu_b[b]):
                continue
            if not (nu_b[b] + mu > _EPS):
                continue
            rec.append((int(i), face, mu, ns, float(nu_b[b]), float(D_ex[i]),
                        float(Ru[b]), float(Ru[i] + Rs[i]), bool(ns > 0 and mu > _EPS)))
    A = np.asarray(rec, dtype=object)
    ii = A[:, 0].astype(int)
    fc = A[:, 1].astype(int)
    mu = A[:, 2].astype(float)
    ns = A[:, 3].astype(float)
    rn = A[:, 4].astype(float)
    D = A[:, 5].astype(float)
    nnu = A[:, 6].astype(float)
    nRex = A[:, 7].astype(float)
    lv = A[:, 8].astype(bool)

    byex: dict[int, dict] = {}
    for k in np.flatnonzero(lv):
        byex.setdefault(int(ii[k]), {})[int(fc[k])] = k
    prs = [(v[1], v[0]) for v in byex.values() if 0 in v and 1 in v]
    kl = np.asarray([p[0] for p in prs], int)
    kr = np.asarray([p[1] for p in prs], int)
    two = np.concatenate([kl, kr])
    mo = np.concatenate([mu[kr], mu[kl]])
    no = np.concatenate([ns[kr], ns[kl]])

    g1 = np.log(mu[kl] / mu[kr])
    p1 = float(np.mean(tri(ns[kl]) + tri(ns[kr])))
    om1 = max(0.0, float(np.var(g1)) - p1) / 2
    cl_t = rn[two] + mu[two]
    x = np.log(cl_t / D[two])
    g2 = x[: len(prs)] - x[len(prs):]
    wmu = mu[two] / cl_t
    wnu = rn[two] / cl_t

    def pois(cl, mu_u, n_u):
        wn, wm = rn[two] / cl, mu_u / cl
        return (wn * wn * tri(nnu[two])
                + np.where(n_u > 0, wm * wm * tri(np.maximum(n_u, 1.0)), 0.0) + tri(nRex[two]))

    P_ship = pois(cl_t, mu[two], ns[two])
    side = P_ship - tri(nRex[two])
    p2 = float(np.mean(side[: len(prs)] + side[len(prs):]))
    om2 = max(0.0, float(np.var(g2)) - p2) / 2
    T2j = float(np.var(x)) - float(np.mean(P_ship))
    Ew2 = float(np.mean(wmu ** 2))
    mx = np.maximum(mu[two], mo)
    cl_mx = rn[two] + mx
    lz = np.log(cl_mx / D[two])
    Vmax_ship = float(np.var(lz)) - float(np.mean(P_ship))
    Vmax_form = float(np.var(lz)) - float(np.mean(pois(cl_mx, mx, np.where(mu[two] >= mo, ns[two], no))))
    r = (cond[5:].replace("_nrna", "").replace("_capture", ""), om1, Ew2,
         float(np.var(g2) / np.var(g1)), om1 * Ew2, om2, T2j, om1 / T2j, om1 * Ew2 / T2j,
         om2 / T2j, Vmax_ship, Vmax_form, float(np.median(lz)))
    acc.append(r)
    print(f"{r[0]:<32}{r[1]:>9.3f}{r[2]:>9.3f}{r[3]:>11.3f}{r[4]:>9.3f}{r[5]:>8.3f}{r[6]:>8.3f}"
          f"{r[7]:>7.2f}{r[8]:>8.2f}{r[9]:>7.2f}{r[10]:>11.4f}{r[11]:>11.4f}{r[12]:>9.4f}")
m = [float(np.mean([a[j] for a in acc])) for j in range(1, 13)]
print("-" * len(hdr))
print(f"{'MEAN over 14 capture-OFF':<32}{m[0]:>9.3f}{m[1]:>9.3f}{m[2]:>11.3f}{m[3]:>9.3f}{m[4]:>8.3f}"
      f"{m[5]:>8.3f}{m[6]:>7.2f}{m[7]:>8.2f}{m[8]:>7.2f}{m[9]:>11.4f}{m[10]:>11.4f}{m[11]:>9.4f}")
