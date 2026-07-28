"""P1e step 2 — the AUDIT measurements (v2).

(a) c = d delta / d log M_dst per destination kind -> the CORRECTED destination-count leg of the null;
(b) where the violation LIVES: share of sum(delta^2) by edge class x destination kind;
(c) the eff-length BOUND, and the DECOMPOSITION that explains why it is exceeded:
        -delta = log(Lambda_src) + log(kappa * B_dst / B_src')
    Lambda_src = (sum_c rho_c^relay E_c^src)/M_src is the SOURCE's own conservation violation, which the
    relay's fuse-after-pin creates and which passes straight through to the destination;
(d) the CHANNEL split of the rank-1 inflation, ON THE SUBSET WHERE IT FIRES;
(e) the ATTRIBUTION test against the oracle (interior nodes only);
(f) the M7 overlap.

    OMP_NUM_THREADS=1 python scratchpad/p1e_2_audit.py [cond ...]
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
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_none_ss_0.50_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
AGG: dict = {}
LAM: list = []
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
    chain, cap = dbg["chain"], dbg["capture"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    n = M.shape[0]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    og, op, on = us["og"], us["op"], us["on"]
    n_node = np.where(~is_bnd, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    f_cur = np.clip(np.asarray(cap["_uni"][-2]["fg_out"], float), 0.0, 1.0)
    B = f_cur / np.maximum(E_g, _EPS) + (1.0 - f_cur) / np.maximum(E_r, _EPS)
    # (c) the SOURCE's own conservation ratio, for each relay direction
    lam_f = (us["fwd_g"] * E_g + (us["fwd_p"] + us["fwd_n"]) * E_r) / np.maximum(M, _EPS)
    lam_b = (us["bwd_g"] * E_g + (us["bwd_p"] + us["bwd_n"]) * E_r) / np.maximum(M, _EPS)
    live = M > _EPS
    LAM.append((cond[5:], float(np.median(np.abs(np.log(np.maximum(lam_f[live], 1e-12))))),
                float(np.median(np.abs(np.log(np.maximum(lam_b[live], 1e-12)))))))

    print(f"\n═══ {cond[5:]} ═══")
    for ent, dl, lam_src in zip(cap["_pin"][-2:], cap["_dl"][-2:], (lam_f, lam_b), strict=True):
        df, src, valid = ent["df"], ent["src"], ent["valid"]
        face = uni["rho_lf"] if df == 0 else uni["rho_rf"]
        q = np.maximum(face / np.maximum(M * B, _EPS) - 1.0, 0.0)
        cM = q / (1.0 + q)
        tg, tp, tn = ent["tg"], ent["tp"], ent["tn"]
        tpg, tpp, tpn = ent["tpg"], ent["tpp"], ent["tpn"]
        s2t, nsrc, graft = ent["s2t"], ent["n_src"], ent["graft"]
        peel = is_bnd & is_exon[src] & valid
        plain = valid & ~graft & ~peel
        p3 = np.stack([tpg, tpp, tpn], axis=1)
        d3 = np.stack([tg, tp, tn], axis=1)
        o3 = np.stack([og, op, on], axis=1)
        E3 = np.stack([E_g, E_r, E_r], axis=1)
        sup = p3 > 0.0
        m3 = np.where(sup, d3, o3) * E3
        S = m3.sum(axis=1)
        ok = valid & (S > _EPS) & (M > _EPS)
        alpha = m3 / np.maximum(S, _EPS)[:, None]
        v3 = np.where(sup, 1.0 / np.maximum(p3, _EPS), np.inf)
        scm = np.maximum(s2t, 0.0) + 1.0 / np.maximum(nsrc, _EPS)
        w3 = np.maximum(np.where(sup, v3, 0.0) - scm[:, None], 0.0)
        s_c = np.where(sup, scm[:, None] + alpha * w3, 0.0)
        sd2 = scm + (alpha * alpha * w3).sum(axis=1)
        delta = np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS))
        nu = cM * cM / np.maximum(n_node, _EPS)
        b2 = np.maximum(0.0, delta * delta - sd2 - nu)
        bnd = np.abs(np.log(np.maximum(E_g, _EPS) / np.maximum(E_r, _EPS)))
        aR = alpha[:, 1] + alpha[:, 2]
        sR = np.where(aR > _EPS, (alpha[:, 1] * s_c[:, 1] + alpha[:, 2] * s_c[:, 2])
                      / np.maximum(aR, _EPS), 0.0)
        d_lam = (s_c[:, 0] - sR) ** 2
        d_lvl = s_c[:, 0] ** 2 + sR * sR
        # residual after removing the SOURCE's own violation
        res = delta + np.log(np.maximum(lam_src[src], 1e-12))
        rg_o = fo * M / np.maximum(E_g, _EPS)
        rR_o = (1.0 - fo) * M / np.maximum(E_r, _EPS)
        eps_g = np.log(np.maximum(tg, 1e-300) / np.maximum(rg_o, 1e-300))
        eps_R = np.log(np.maximum(tp + tn, 1e-300) / np.maximum(rR_o, 1e-300))
        aa = np.where(sd2 > _EPS, delta / np.maximum(sd2, _EPS), 0.0)
        corr_g, corr_R = aa * s_c[:, 0], aa * sR
        for name, msk in (("plain", plain), ("graft", graft), ("peel", peel)):
            for dk, dm in (("reg", ~is_bnd), ("bnd", is_bnd)):
                m = msk & dm & ok
                if m.sum() < 20:
                    continue
                a = AGG.setdefault((name, dk), dict(
                    n=0, d2=0.0, b2=0.0, lam=0.0, lvl=0.0, over=0, fires=0,
                    z2m=[], cm=[], comf=[], ad=[], ares=[]))
                a["n"] += int(m.sum())
                a["d2"] += float((delta[m] ** 2).sum())
                a["b2"] += float(b2[m].sum())
                a["lam"] += float((b2[m] * d_lam[m] / np.maximum(sd2[m], _EPS) ** 2).sum())
                a["lvl"] += float((b2[m] * d_lvl[m] / np.maximum(sd2[m], _EPS) ** 2).sum())
                a["over"] += int((np.abs(delta[m]) > bnd[m] + 1e-9).sum())
                a["fires"] += int((b2[m] > 0).sum())
                a["z2m"].append(delta[m] ** 2 / np.maximum(sd2[m] + nu[m], _EPS))
                a["cm"].append(cM[m])
                a["ad"].append(np.abs(delta[m]))
                a["ares"].append(np.abs(res[m]))
                f = m & (b2 > 0)
                if f.any():
                    a["comf"].append(scm[f] / np.maximum(sd2[f], _EPS))
        gm = graft & ok & np.isfinite(fo) & (fo > 0.02) & (fo < 0.98) & (b2 > 0) & (tg > 0) & ((tp + tn) > 0)
        if gm.sum() > 20:
            worse_g = np.abs(eps_g[gm]) > np.abs(eps_R[gm])
            moves_g = np.abs(corr_g[gm]) > np.abs(corr_R[gm])
            print(f"  df={df} GRAFT b2>0 interior n={int(gm.sum()):4d}   correction targets the WRONGER arm "
                  f"{float(np.mean(worse_g == moves_g)):5.1%}   sign(corr)=-sign(err): "
                  f"gDNA {float(np.mean(np.sign(corr_g[gm]) == -np.sign(eps_g[gm]))):5.1%}  "
                  f"RNA {float(np.mean(np.sign(corr_R[gm]) == -np.sign(eps_R[gm]))):5.1%}   "
                  f"med|eps_g|={np.median(np.abs(eps_g[gm])):.2f} med|eps_R|={np.median(np.abs(eps_R[gm])):.2f}")
        mm = ok & (np.abs(dl["G_g"]) > 0) & (np.abs(dl["G_lam"]) > 0)
        if mm.sum() > 50:
            print(f"  df={df} M7 overlap  corr(delta,G_g) = {float(np.corrcoef(delta[mm], dl['G_g'][mm])[0, 1]):+.3f}"
                  f"   corr(delta,G_lam) = {float(np.corrcoef(delta[mm], dl['G_lam'][mm])[0, 1]):+.3f}"
                  f"   (n={int(mm.sum())})")

print("\n\n═══ the SOURCE's own conservation violation  median |log Lambda_src|  (relay fwd / bwd) ═══")
for c, a, b in LAM:
    print(f"  {c:<46} fwd {a:6.3f}   bwd {b:6.3f}   (ratio {np.exp(a):5.2f}x / {np.exp(b):5.2f}x)")

print("\n═══ AGGREGATE over conditions ═══")
tot = sum(a["d2"] for a in AGG.values())
print(f"{'edge':>6} {'dst':>4} {'n':>7} {'sh d2':>7} {'med|d|':>7} {'med|res|':>9} {'medz2':>7} "
      f"{'fires':>6} {'E[b2]':>7} {'>bnd':>6} {'lam/lvl':>8} {'com|fire':>9} {'med c_M':>8}")
for (name, dk), a in sorted(AGG.items()):
    z2 = np.concatenate(a["z2m"])
    comf = np.concatenate(a["comf"]) if a["comf"] else np.zeros(1)
    print(f"{name:>6} {dk:>4} {a['n']:>7} {a['d2'] / tot:>6.1%} "
          f"{float(np.median(np.concatenate(a['ad']))):>7.3f} "
          f"{float(np.median(np.concatenate(a['ares']))):>9.3f} {float(np.median(z2)):>7.3f} "
          f"{a['fires'] / a['n']:>5.1%} {a['b2'] / a['n']:>7.3f} {a['over'] / a['n']:>5.1%} "
          f"{a['lam'] / max(a['lvl'], 1e-12):>8.3f} {float(np.median(comf)):>9.3f} "
          f"{float(np.median(np.concatenate(a['cm']))):>8.3f}")
