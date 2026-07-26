"""RC2 (A/B) — the GRAFT probe: how big is the boundary's measured-mature RNA claim, what is `r` doing to
it, and why does the DL mismatch term not kill it?

Same verbatim `_transport` rebuild as `rc2_a_decomp.py`, but instrumented at the graft:

    relay_rna   = rp[src]·r·E_r/M          the RNA the source RELAYED (imputed)
    graft_rna   = gp·r·E_r/M               the boundary's MEASURED mature spliced density, ×r
    graft_nor   = gp·E_r/M                 …the same measurement NOT reframed  (the counterfactual)
    true_rna    = 1 − oracle
plus the DerSimonian–Laird inputs (G, v_msg, v_own, b̂², p_pre → p_post) for the RNA arms.

    OMP_NUM_THREADS=1 python scratchpad/rc2_a_graft.py [COND]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.enrichment_frame import mismatch_deflate, mismatch_gap, transfer_logvar  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.node_init import own_composition_logvar  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9
BINS = [(0.0, 0.30), (0.30, 0.60), (0.60, 0.90), (0.90, 0.99), (0.99, 1.01)]

ap = argparse.ArgumentParser()
ap.add_argument("conds", nargs="*", default=["gdna_gdna300_ss_0.99_nrna_none_capture_on",
                                             "gdna_gdna300_ss_0.99_nrna_present_capture_on",
                                             "gdna_gdna100_ss_0.99_nrna_present_capture_on"])
ap.add_argument("--nodes", default="2197,1909,3237,2173,1523,2683")
a = ap.parse_args()

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
POOL: dict[str, list] = {}


def run(cond, want=()):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    n = np.asarray(cap["f_g"]).shape[0]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    og, op, on = us["og"], us["op"], us["on"]
    is_bnd, ex_a = us["is_bnd"].astype(bool), us["is_exon"].astype(bool)
    fp_a, fn_a = us["fp"].astype(bool), us["fn"].astype(bool)
    li, ri = us["left"], us["right"]
    logvar_tot, tau_own = us["logvar_tot"], us["tau_own"]
    struct = us["struct_lock"].astype(bool)
    ni_fg = np.asarray(cap["fg_loc"])
    v_own_g, v_own_r = own_composition_logvar(ni_fg, tau_own, struct)
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    vl, vr = li >= 0, ri >= 0
    sl, sr = np.clip(li, 0, n - 1), np.clip(ri, 0, n - 1)

    def transport(src, valid, df, sf, arrs, dst_face_v, src_face_v):
        rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
        framed = valid & (src_face_v[src] > _EPS) & (dst_face_v > _EPS)
        r = np.where(framed, dst_face_v / np.maximum(src_face_v[src], _EPS), np.where(valid, 1.0, 0.0))
        graft = ex_a & is_bnd[src] & valid
        gp = np.where(graft, spl_p_f[sf][src], 0.0)
        gn = np.where(graft, spl_n_f[sf][src], 0.0)
        tg, tp, tn = rg[src] * r, (rp[src] + gp) * r, (rn[src] + gn) * r
        relay_only = (rp[src] + rn[src]) * r
        graft_only = (gp + gn) * r
        graft_nor = gp + gn
        s2t = transfer_logvar(logvar_tot, logvar_tot[src], graft)

        def _dv(p_, s2_=s2t):
            return np.where(valid & (p_[src] > 0.0), 1.0 / (1.0 / np.maximum(p_[src], _EPS) + s2_), 0.0)

        tpg, tpp, tpn = _dv(pg), _dv(pp), _dv(pn)
        tmg, tmp, tmn = _dv(mg), _dv(mp), _dv(mn)
        _sp, _sn = np.where(graft, SP[sf][src], 0.0), np.where(graft, SN[sf][src], 0.0)
        _s2 = np.where(np.isfinite(s2t), s2t, 0.0)
        _spc = np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2), 0.0)
        _snc = np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2), 0.0)
        tpp, tpn = tpp + _spc, tpn + _snc
        tmp, tmn = tmp + _spc, tmn + _snc
        peel = is_bnd & ex_a[src] & valid
        tp = np.where(peel, np.maximum(tp - spl_p_f[df], 0.0), tp)
        tn = np.where(peel, np.maximum(tn - spl_n_f[df], 0.0), tn)
        tp, tpp, tmp = np.where(fp_a, tp, 0.0), np.where(fp_a, tpp, 0.0), np.where(fp_a, tmp, 0.0)
        tn, tpn, tmn = np.where(fn_a, tn, 0.0), np.where(fn_a, tpn, 0.0), np.where(fn_a, tmn, 0.0)
        sg = np.where(tpg > 0.0, tg, og)
        sp_ = np.where(tpp > 0.0, tp, op)
        sn_ = np.where(tpn > 0.0, tn, on)
        s = sg * E_g + (sp_ + sn_) * E_r
        k = np.where((s > _EPS) & (M > _EPS), M / np.maximum(s, _EPS), 1.0)
        tg, tp, tn = tg * k, tp * k, tn * k
        g_g, c_g = mismatch_gap(tg, og)
        g_p, c_p = mismatch_gap(tp, op)
        g_n, c_n = mismatch_gap(tn, on)
        dl = dict(G_p=g_p, G_n=g_n, G_g=g_g, c_p=c_p, c_n=c_n, pre_pp=tpp.copy(), pre_pn=tpn.copy(),
                  pre_mp=tmp.copy(), pre_mn=tmn.copy(), v_own_r=v_own_r, v_own_g=v_own_g)
        tpg = mismatch_deflate(tpg, g_g, c_g, v_own_g)
        tpp = mismatch_deflate(tpp, g_p, c_p, v_own_r)
        tpn = mismatch_deflate(tpn, g_n, c_n, v_own_r)
        tmg = mismatch_deflate(tmg, g_g, c_g, v_own_g)
        tmp = mismatch_deflate(tmp, g_p, c_p, v_own_r)
        tmn = mismatch_deflate(tmn, g_n, c_n, v_own_r)
        dl.update(post_pp=tpp, post_pn=tpn, post_mp=tmp, post_mn=tmn)
        return dict(tg=tg, tp=tp, tn=tn, tpg=tpg, tpp=tpp, tpn=tpn, tmg=tmg, tmp=tmp, tmn=tmn,
                    r=r, k=k, graft=graft, relay_only=relay_only, graft_only=graft_only,
                    graft_nor=graft_nor, dl=dl, src=src, valid=valid)

    fwd = tuple(us[f"fwd_{x}"] for x in ("g", "p", "n", "pg", "pp", "pn", "mg", "mp", "mn", "tau"))
    bwd = tuple(us[f"bwd_{x}"] for x in ("g", "p", "n", "pg", "pp", "pn", "mg", "mp", "mn", "tau"))
    A = transport(sl, vl, 0, 1, fwd, uni["rho_lf"], uni["rho_rf"])
    B = transport(sr, vr, 1, 0, bwd, uni["rho_rf"], uni["rho_lf"])
    # fidelity on the RNA-arm precisions ψ actually got
    f_mp = np.max(np.abs(A["tmp"] + B["tmp"] - uni["cm_p"]))
    f_mn = np.max(np.abs(A["tmn"] + B["tmn"] - uni["cm_n"]))
    print(f"\n### {cond}   RNA-precision fidelity |cm_p|={f_mp:.1e} |cm_n|={f_mn:.1e}")

    Gp_, Gn_, Rp_, Rn_ = _oracle_per_node(inp, chain)
    G_, R_ = Gp_ + Gn_, Rp_ + Rn_
    fo = np.where(G_ + R_ > _EPS, G_ / np.maximum(G_ + R_, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(n)])
    mass = np.asarray(cap["mass_global"])

    for T in (A, B):
        m = (np.isfinite(fo) & (mass > 0) & (tau_own > 0) & (cls == "exon") & T["valid"] & T["graft"]
             & ((T["tmp"] + T["tmn"]) > 0))
        if not m.any():
            continue
        sc = E_r[m] / np.maximum(M[m], _EPS)
        d = T["dl"]
        rna_p_live = fp_a[m]
        Gsel = np.where(rna_p_live, d["G_p"][m], d["G_n"][m])
        pre = np.where(rna_p_live, d["pre_mp"][m], d["pre_mn"][m])
        post = np.where(rna_p_live, d["post_mp"][m], d["post_mn"][m])
        for k, v in dict(
            oracle=fo[m], mass=mass[m], true_rna=1.0 - fo[m],
            relay_rna=T["relay_only"][m] * sc, graft_rna=T["graft_only"][m] * sc,
            graft_nor=T["graft_nor"][m] * sc, deliv_rna=(T["tp"][m] + T["tn"][m]) * sc,
            r=T["r"][m], G=Gsel, v_own=d["v_own_r"][m], pre_p=pre, post_p=post,
            own_rna=1.0 - ni_fg[m],
        ).items():
            POOL.setdefault(k, []).append(np.asarray(v, np.float64))

    for i in want:
        if i >= n:
            continue
        print(f"\n  node {i} [{cls[i]}] oracle={fo[i]:.4f} own={ni_fg[i]:.4f} solved={cap['f_g'][i]:.4f} "
              f"tau_own={tau_own[i]:.4g} v_own_r={v_own_r[i]:.4g}")
        for side, T in (("L", A), ("R", B)):
            if not T["valid"][i]:
                continue
            sc = E_r[i] / max(M[i], _EPS)
            d = T["dl"]
            live = "p" if fp_a[i] else "n"
            print(f"    {side} src={int(T['src'][i])}[{cls[int(T['src'][i])]}] r={T['r'][i]:.4g} "
                  f"graft={bool(T['graft'][i])} kpin={T['k'][i]:.4g}")
            print(f"       RNA frac: relay={T['relay_only'][i] * sc:.4f}  graft(×r)="
                  f"{T['graft_only'][i] * sc:.4f}  graft(NO r)={T['graft_nor'][i] * sc:.4f}  "
                  f"delivered={(T['tp'][i] + T['tn'][i]) * sc:.4f}   TRUE={1 - fo[i]:.4f}   "
                  f"own={1 - ni_fg[i]:.4f}")
            print(f"       DL[{live}]: G={d['G_' + live][i]:+.3f} G^2={d['G_' + live][i] ** 2:.3f} "
                  f"v_own_r={d['v_own_r'][i]:.4g} contra={bool(d['c_' + live][i])} "
                  f"| mode-prec {d['pre_p' + live][i]:.4g}→{d['post_p' + live][i]:.4g}"
                  f" | meas-prec {d['pre_m' + live][i]:.4g}→{d['post_m' + live][i]:.4g}")


want = [int(x) for x in a.nodes.split(",")] if a.nodes else []
for i, c in enumerate(a.conds):
    run(c, want=want if i == 0 else ())

D = {k: np.concatenate(v) for k, v in POOL.items()}
w = D["mass"]
print(f"\n{'=' * 150}\nGRAFT MAGNITUDE vs TRUTH, full-rank EXON destinations with a live RNA measurement "
      f"(n={w.size:,})\n{'=' * 150}")
print(f"{'oracle bin':<13}{'n':>6}{'TRUE rna':>10}{'own rna':>9}{'relay rna':>11}{'graft(xr)':>11}"
      f"{'graft(NO r)':>13}{'delivered':>11}{'  |  median r':>14}{'p90 r':>10}"
      f"{'  DL: mean G^2':>15}{'mean v_own':>12}{'% fired':>9}{'prec pre':>10}{'prec post':>11}")
for lo, hi in BINS:
    m = (D["oracle"] >= lo) & (D["oracle"] < hi)
    if not m.any():
        continue
    ww = w[m]

    def av(k):
        return np.average(D[k][m], weights=ww)

    fired = np.average((D["post_p"][m] < D["pre_p"][m] - 1e-9).astype(float), weights=ww)
    print(f"[{lo:.2f},{hi:.2f})  {int(m.sum()):>6}{av('true_rna'):>10.4f}{av('own_rna'):>9.4f}"
          f"{av('relay_rna'):>11.4f}{av('graft_rna'):>11.4f}{av('graft_nor'):>13.4f}"
          f"{av('deliv_rna'):>11.4f}{np.median(D['r'][m]):>17.3f}{np.quantile(D['r'][m], 0.9):>10.1f}"
          f"{av('G') ** 2:>15.2f}{av('v_own'):>12.3g}{fired:>9.1%}{av('pre_p'):>10.2f}{av('post_p'):>11.2f}")
