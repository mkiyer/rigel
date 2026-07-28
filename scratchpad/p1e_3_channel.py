"""P1e step 3 — THE CHANNEL DECISION, measured.

Two questions, both answered against the oracle with no code change:

(A) DIRECTION. On the subset where the surprise FIRES (b2 > 0), how much of the rank-1 inflation
    Delta Sigma = b2 * s s^T / sd2^2 reaches the SPLIT (lambda) versus the LEVELS?
        Delta Var(lambda) = b2*(s_g - s_R)^2/sd2^2 ,  Delta v_c = b2*s_c^2/sd2^2
    and is the correction's DIFFERENTIAL part a_g - a_R = delta*(s_g - s_R)/sd2 pointed the right way,
    i.e. against the true composition error eps_g - eps_R (oracle)?

(B) EFFECT. Offline prototype: re-fuse the two messages with P1e-damped precisions and score the
    composition the node is HANDED, against the oracle, exactly as `wp1_prototype.py` does. Three arms:
      lvl   -- damp only the per-component LEVEL precisions (cm_*/mode-fusion), tau untouched
      lam   -- damp only the tau (lambda) stream
      both  -- the full rank-1 law
    scored on the density fuse (f_g from the fused densities) AND on the lambda message the packet
    actually delivers (f_g = sigmoid(lam_msg)).

    OMP_NUM_THREADS=1 python scratchpad/p1e_3_channel.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def msg_state(ent, us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B):
    df, src, valid = ent["df"], ent["src"], ent["valid"]
    og, op, on = us["og"], us["op"], us["on"]
    face = uni["rho_lf"] if df == 0 else uni["rho_rf"]
    q = np.maximum(face / np.maximum(M * B, _EPS) - 1.0, 0.0)
    cM = q / (1.0 + q)
    tg, tp, tn = ent["tg"], ent["tp"], ent["tn"]
    p3 = np.stack([ent["tpg"], ent["tpp"], ent["tpn"]], axis=1)
    d3 = np.stack([tg, tp, tn], axis=1)
    o3 = np.stack([og, op, on], axis=1)
    E3 = np.stack([E_g, E_r, E_r], axis=1)
    sup = p3 > 0.0
    m3 = np.where(sup, d3, o3) * E3
    S = m3.sum(axis=1)
    ok = valid & (S > _EPS) & (M > _EPS)
    alpha = m3 / np.maximum(S, _EPS)[:, None]
    v3 = np.where(sup, 1.0 / np.maximum(p3, _EPS), np.inf)
    scm = np.maximum(ent["s2t"], 0.0) + 1.0 / np.maximum(ent["n_src"], _EPS)
    w3 = np.maximum(np.where(sup, v3, 0.0) - scm[:, None], 0.0)
    s_c = np.where(sup, scm[:, None] + alpha * w3, 0.0)
    sd2 = scm + (alpha * alpha * w3).sum(axis=1)
    delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
    nu = cM * cM / np.maximum(n_node, _EPS)
    b2 = np.where(ok, np.maximum(0.0, delta * delta - sd2 - nu), 0.0)
    aR = alpha[:, 1] + alpha[:, 2]
    sR = np.where(aR > _EPS, (alpha[:, 1] * s_c[:, 1] + alpha[:, 2] * s_c[:, 2])
                  / np.maximum(aR, _EPS), 0.0)
    dv = b2[:, None] * s_c * s_c / np.maximum(sd2, _EPS)[:, None] ** 2   # RANK1_s level inflation
    a1 = np.maximum(alpha.sum(axis=1), _EPS)
    dvc = np.repeat((b2 / (a1 * a1))[:, None], 3, axis=1)                # COMMON-direction inflation
    dlam = b2 * (s_c[:, 0] - sR) ** 2 / np.maximum(sd2, _EPS) ** 2       # the split inflation
    graft = ent["graft"]
    peel = is_bnd & is_exon[src] & valid
    return dict(df=df, src=src, valid=valid, ok=ok, delta=delta, sd2=sd2, nu=nu, b2=b2, alpha=alpha,
                s_c=s_c, sR=sR, dv=dv, dlam=dlam, scm=scm, p3=p3, d3=d3, sup=sup, S=S,
                graft=graft, peel=peel, plain=valid & ~graft & ~peel, dvc=dvc,
                partial=valid & ~sup.all(axis=1))


ROWS = []
DIR = {}
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
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    n = M.shape[0]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    n_node = np.where(~is_bnd, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    f_cur = np.clip(np.asarray(cap["_uni"][-2]["fg_out"], float), 0.0, 1.0)
    B = f_cur / np.maximum(E_g, _EPS) + (1.0 - f_cur) / np.maximum(E_r, _EPS)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    dcls = np.array([CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)

    A = msg_state(cap["_pin"][-2], us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B)
    Bm = msg_state(cap["_pin"][-1], us, uni, M, E_g, E_r, is_bnd, is_exon, n_node, B)

    # ── (A) the direction, b2-weighted, on the firing subset ─────────────────────────────────────────
    rg_o = fo * M / np.maximum(E_g, _EPS)
    rR_o = (1.0 - fo) * M / np.maximum(E_r, _EPS)
    for tag, X in (("L", A), ("R", Bm)):
        for name in ("plain", "graft", "peel"):
            m = X[name] & X["ok"] & (X["b2"] > 0)
            if m.sum() < 20:
                continue
            wgt = X["b2"][m]
            comf = float(np.average(X["scm"][m] / np.maximum(X["sd2"][m], _EPS), weights=wgt))
            lamratio = float(np.average(
                X["dlam"][m] / np.maximum(X["dv"][m, 0] + X["dlam"][m], _EPS), weights=wgt))
            d = DIR.setdefault(name, [0, 0.0, 0.0, 0, 0, []])
            d[0] += int(m.sum())
            d[1] += comf * m.sum()
            d[2] += lamratio * m.sum()
            tt0 = np.asarray(cap["_dl"][-2 if tag == "L" else -1]["tau_post"], float)
            mt = m & (tt0 > 0)
            if mt.any():  # how big is the SPLIT inflation against the split precision it damps?
                d[5].append(tt0[mt] * X["dlam"][mt])
            # the DIFFERENTIAL attribution test against the oracle
            gm = m & np.isfinite(fo) & (fo > 0.02) & (fo < 0.98) & (X["d3"][:, 0] > 0) & \
                ((X["d3"][:, 1] + X["d3"][:, 2]) > 0)
            if gm.any():
                eps_lam = (np.log(np.maximum(X["d3"][gm, 0], 1e-300) / np.maximum(rg_o[gm], 1e-300))
                           - np.log(np.maximum(X["d3"][gm, 1] + X["d3"][gm, 2], 1e-300)
                                    / np.maximum(rR_o[gm], 1e-300)))
                a_lam = X["delta"][gm] * (X["s_c"][gm, 0] - X["sR"][gm]) / np.maximum(X["sd2"][gm], _EPS)
                d[3] += int((np.sign(a_lam) == -np.sign(eps_lam)).sum())
                d[4] += int(gm.sum())

    # ── (B) the offline effect on the composition the node is handed ────────────────────────────────
    # The per-message lambda is INVARIANT to `_pin_v`'s common factor (it cancels from
    # log(tg*E_g) - log((tp+tn)*E_r)), so it is reconstructable exactly from the PRE-pin capture; the
    # per-message tau is `_capture["_dl"]`'s post-DL value.
    sel = ((dcls == "exon") & np.isfinite(fo) & (M > _EPS)
           & np.asarray(cap["solvable"], bool) & A["ok"] & Bm["ok"])
    row = {"cond": cond[5:], "n": int(sel.sum()), "reads": float(M[sel].sum()),
           "part": float(np.mean(np.concatenate([A["partial"][A["ok"]], Bm["partial"][Bm["ok"]]])))}
    tau = [np.asarray(cap["_dl"][-2]["tau_post"], float), np.asarray(cap["_dl"][-1]["tau_post"], float)]
    lam = []
    for X in (A, Bm):
        tg, tR = X["d3"][:, 0], X["d3"][:, 1] + X["d3"][:, 2]
        lam.append(np.where((tg > _EPS) & (tR > _EPS),
                            np.log(np.maximum(tg * E_g, _EPS)) - np.log(np.maximum(tR * E_r, _EPS)), 0.0))
    for arm in ("shipped", "lvl", "lvlC", "lam"):
        pg, dens, tt = [], [], []
        for k, X in enumerate((A, Bm)):
            p = X["p3"].copy()
            if arm in ("lvl", "lvlC"):
                v = np.where(X["sup"], 1.0 / np.maximum(p, _EPS), np.inf)
                dd = X["dv"] if arm == "lvl" else X["dvc"]
                p = np.where(X["sup"], 1.0 / np.maximum(v + dd, _EPS), 0.0)
            t = tau[k].copy()
            if arm == "lam":
                t = np.where(t > 0.0, 1.0 / np.maximum(1.0 / np.maximum(t, _EPS) + X["dlam"], _EPS), 0.0)
            pg.append(p)
            dens.append(X["d3"])
            tt.append(t)
        pa, pb = pg[0], pg[1]
        cd = np.where((pa + pb) > _EPS, (pa * dens[0] + pb * dens[1]) / np.maximum(pa + pb, _EPS), 0.0)
        gm_, rm_ = cd[:, 0] * E_g, (cd[:, 1] + cd[:, 2]) * E_r
        f = np.where((gm_ + rm_) > _EPS, gm_ / np.maximum(gm_ + rm_, _EPS), np.nan)
        m2 = sel & np.isfinite(f)
        row["D_" + arm] = float(np.average(np.abs(f[m2] - fo[m2]), weights=M[m2]))
        ct = tt[0] + tt[1]
        lm = np.where(ct > _EPS, (tt[0] * lam[0] + tt[1] * lam[1]) / np.maximum(ct, _EPS), np.nan)
        fl = 1.0 / (1.0 + np.exp(-np.clip(lm, -50, 50)))
        m3 = sel & np.isfinite(fl) & (ct > _EPS)
        row["L_" + arm] = float(np.average(np.abs(fl[m3] - fo[m3]), weights=M[m3]))
        row["L_n"] = int(m3.sum())
    ROWS.append(row)

print("\n═══ (A) the DIRECTION on the firing subset (b2-weighted) ═══")
print(f"{'edge':>7}{'n fire':>9}{'common frac':>13}{'lam/(lam+lvl_g)':>18}{'diff sign OK':>14}"
      f"{'med tau*dlam':>14}{'p90':>9}")
for name, d in DIR.items():
    td = np.concatenate(d[5]) if d[5] else np.zeros(1)
    print(f"{name:>7}{d[0]:>9}{d[1] / d[0]:>13.3f}{d[2] / d[0]:>18.3f}"
          f"{(d[3] / d[4] if d[4] else float('nan')):>13.1%}"
          f"{float(np.median(td)):>14.3g}{float(np.percentile(td, 90)):>9.3g}")

print("\n═══ (B) exons, mass-weighted |f_g - oracle| of the composition the node is HANDED ═══")
print("     DENSITY fuse (the LEVEL channels)          |  LAMBDA message (the SPLIT channel, sigmoid)")
print(f"{'condition':<40}{'n':>5}{'D:shipped':>11}{'D:rank1':>9}{'D:common':>10}{'D:lam':>8}"
      f"{'| L:n':>8}{'L:shipped':>11}{'L:lam':>8}{'part%':>7}")
for r in ROWS:
    print(f"{r['cond']:<40}{r['n']:>5}{r['D_shipped']:>11.4f}{r['D_lvl']:>9.4f}{r['D_lvlC']:>10.4f}"
          f"{r['D_lam']:>8.4f}{r['L_n']:>8}{r['L_shipped']:>11.4f}{r['L_lam']:>8.4f}{r['part']:>7.1%}")
