"""BOUNDARY LEVEL step 2 — the level is ALWAYS `(1-f_g)*M/E_r`, so fuse the estimates of f_g, don't rank them.

w1 showed the MASS-identity level (`rho_nu = (M - tg*E_g)/E_r`, i.e. the message's gDNA density claim closed
against the boundary's OWN observed mass) beats HEAD's subtraction on the metric that matters — the delivered
f_g — in every regime.  It also exists at EVERY seam, which is what kills level-precedence case 3.

Every candidate level has the same shape:  rho_nu = (1 - f_hat) * M / E_r.  So there is no precedence to
choose: there are two independent estimators of the same f_hat and they FUSE by inverse variance —

    OWN   f_g from the node's own self-solve,  precision tau_own   (0 on unstranded/no-factory => inert)
    MASS  f_g = phi = tg*E_g/M,                from the imputed gDNA DENSITY + the node's own count

and the fuse is done on the LEVEL in log space with the derived variances

    v_own  = Var(log f_R) + 1/n                      (transport_seed_logvar, already shipped)
    v_mass = [ 1/n + phi^2 * v_g ] / (1 - phi)^2     (M11, derived here: rho_nu*E_r = M - tg*E_g)

This script scores the fuse against each part, on JUNCTION-BEARING edges only (where the peel acts), and
checks the calibration of Var(log w).

    OMP_NUM_THREADS=1 python scratchpad/w2_fuse.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
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
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)

print(f"{'condition':<44}{'stratum':<20}{'n':>6}{'w_true':>8}"
      f"{'  dfg: flux':>12}{'own':>8}{'mass':>8}{'fuse':>8}{'oracle_w':>10}{'  z2(fuse)':>11}")
print("-" * 146)

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
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    n = kind.shape[0]

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = (pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg"))
    fo = np.where(G + Ru > _EPS, G / np.maximum(G + Ru, _EPS), np.nan)

    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    op, on = us["op"], us["on"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    tau_own, lv = us["tau_own"], us["logvar_tot"]
    fg_own = np.asarray(cap["fg_loc"], float)
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    NSP = (np.asarray(us["spl_n_pos_l"]) + np.asarray(us["spl_n_neg_l"]),
           np.asarray(us["spl_n_pos_r"]) + np.asarray(us["spl_n_neg_r"]))
    n_node = np.where(isr, np.asarray(us["n_unspl_l"]),
                      np.asarray(us["n_unspl_l"]) + np.asarray(us["n_unspl_r"]))
    rho_lf, rho_rf = uni["rho_lf"], uni["rho_rf"]
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    solv = np.asarray(cap["solvable"], bool)

    acc: dict[str, list] = {}
    for tag, (src_i, dst_face, src_face, df, fw) in (
        ("L", (li, rho_lf, rho_rf, 0, "fwd")),
        ("R", (ri, rho_rf, rho_lf, 1, "bwd")),
    ):
        s = np.clip(src_i, 0, n - 1)
        valid = src_i >= 0
        peel = is_bnd & is_exon[s] & valid
        rg, rp, rn = us[f"{fw}_g"], us[f"{fw}_p"], us[f"{fw}_n"]
        pg, pp, pn = us[f"{fw}_pg"], us[f"{fw}_pp"], us[f"{fw}_pn"]
        framed = valid & (src_face[s] > _EPS) & (dst_face > _EPS)
        r = np.where(framed, dst_face / np.maximum(src_face[s], _EPS), np.where(valid, 1.0, 0.0))
        tg, A = rg[s] * r, (rp[s] + rn[s]) * r
        mu = (np.where(SPf[df] > _EPS, SPf[df] / np.maximum(ESP[df], _EPS), 0.0)
              + np.where(SNf[df] > _EPS, SNf[df] / np.maximum(ESP[df], _EPS), 0.0))
        s2t = lv + lv[s]
        v_g = np.where(pg[s] > 0, 1.0 / np.maximum(pg[s], _EPS) + s2t, np.inf)
        v_A = np.where((pp[s] + pn[s]) > 0, 1.0 / np.maximum(pp[s] + pn[s], _EPS) + s2t, np.inf)
        v_mu = np.where(NSP[df] > 0, 1.0 / np.maximum(NSP[df], _EPS), np.inf)
        inv_n = 1.0 / np.maximum(n_node, _EPS)

        nu_true = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), 0.0)
        w_true = np.where(nu_true + mu > _EPS, nu_true / np.maximum(nu_true + mu, _EPS), 1.0)

        # ── OWN level ──
        fgo = np.clip(fg_own, 0.0, 1.0)
        v_fr_own = np.where(tau_own > _EPS, fgo * fgo / np.maximum(tau_own, _EPS), np.inf)
        nu_own = op + on
        v_own = np.where((tau_own > _EPS) & (nu_own > _EPS), v_fr_own + inv_n, np.inf)

        # ── MASS level (M11) ──
        phi = np.where(M > _EPS, tg * E_g / np.maximum(M, _EPS), np.inf)
        resid = np.where(M > _EPS, M - tg * E_g, 0.0)
        nu_mass = np.maximum(resid, 0.0) / np.maximum(E_r, _EPS)
        v_mass = np.where(
            (phi < 1.0) & (pg[s] > 0) & (nu_mass > _EPS),
            (inv_n + phi * phi * np.where(np.isfinite(v_g), v_g, 0.0))
            / np.maximum((1.0 - np.minimum(phi, 1.0)) ** 2, _EPS), np.inf)

        # ── FUSE in log space ──
        po = np.where(np.isfinite(v_own) & (nu_own > _EPS), 1.0 / np.maximum(v_own, _EPS), 0.0)
        pm = np.where(np.isfinite(v_mass) & (nu_mass > _EPS), 1.0 / np.maximum(v_mass, _EPS), 0.0)
        pt = po + pm
        nu_f = np.where(pt > _EPS, np.exp((po * np.log(np.maximum(nu_own, _EPS))
                                           + pm * np.log(np.maximum(nu_mass, _EPS)))
                                          / np.maximum(pt, _EPS)), 0.0)
        v_f = np.where(pt > _EPS, 1.0 / np.maximum(pt, _EPS), np.inf)

        def share(nu, v):
            w = np.where(nu + mu > _EPS, nu / np.maximum(nu + mu, _EPS), 1.0)
            vw = np.where(np.isfinite(v), (1.0 - w) ** 2 * (np.where(np.isfinite(v), v, 0.0)
                                                            + np.where(np.isfinite(v_mu), v_mu, 0.0)),
                          np.inf)
            return w, vw

        w_own, vw_own = share(nu_own, v_own)
        w_mass, vw_mass = share(nu_mass, v_mass)
        w_fuse, vw_fuse = share(nu_f, v_f)
        w_flux = np.where(A > _EPS, np.maximum(A - mu, 0.0) / np.maximum(A, _EPS), 1.0)
        # case-3 semantics: no level at all -> w=0 at zero precision
        w_fuse = np.where(np.isfinite(v_f), w_fuse, 0.0)

        sel = peel & np.isfinite(fo) & (M > _EPS) & (A > _EPS) & solv
        far = ri if df == 0 else li
        cl_far = np.array([cls[j] if j >= 0 else "edge" for j in far])
        acc.setdefault("rows", []).append(dict(
            sel=sel, M=M, tg=tg, A=A, E_g=E_g, E_r=E_r, fo=fo, w_true=w_true, mu=mu,
            cl_far=cl_far, has_own=np.isfinite(v_own), has_mass=np.isfinite(v_mass),
            ws={"flux": (w_flux, None), "own": (w_own, vw_own), "mass": (w_mass, vw_mass),
                "fuse": (w_fuse, vw_fuse)},
        ))

    def dfg_of(R, m, w):
        num = R["tg"][m] * R["E_g"][m]
        den = num + R["A"][m] * w * R["E_r"][m]
        return np.abs(np.where(den > _EPS, num / np.maximum(den, _EPS), 1.0) - R["fo"][m])

    def report(label, mask_fn):
        tot = 0
        sums = {k: 0.0 for k in ("flux", "own", "mass", "fuse", "oracle")}
        wsum = 0.0
        z_num = z_den = 0.0
        wt_num = 0.0
        for R in acc["rows"]:
            m = R["sel"] & mask_fn(R)
            if not m.any():
                continue
            tot += int(m.sum())
            wgt = R["M"][m]
            wsum += float(wgt.sum())
            wt_num += float(np.sum(R["w_true"][m] * wgt))
            for k, (we, _v) in R["ws"].items():
                sums[k] += float(np.sum(dfg_of(R, m, we[m]) * wgt))
            sums["oracle"] += float(np.sum(dfg_of(R, m, R["w_true"][m]) * wgt))
            we, vw = R["ws"]["fuse"]
            good = m & (R["w_true"] > _EPS) & np.isfinite(R["ws"]["fuse"][1]) & (we > _EPS)
            if good.any():
                d = np.log(we[good]) - np.log(R["w_true"][good])
                z_num += float(np.sum(d * d)); z_den += float(np.sum(vw[good]))
        if wsum <= 0:
            return
        z2 = z_num / z_den if z_den > 0 else np.nan
        print(f"{'':<44}{label:<20}{tot:>6}{wt_num / wsum:>8.3f}"
              f"{sums['flux'] / wsum:>12.3f}{sums['own'] / wsum:>8.3f}{sums['mass'] / wsum:>8.3f}"
              f"{sums['fuse'] / wsum:>8.3f}{sums['oracle'] / wsum:>10.3f}{z2:>11.2f}")

    print(f"{cond[5:]:<44}")
    report("junction-bearing", lambda R: R["mu"] > _EPS)
    report("  far=intron", lambda R: (R["mu"] > _EPS) & (R["cl_far"] == "intron"))
    report("  far=exon", lambda R: (R["mu"] > _EPS) & (R["cl_far"] == "exon"))
    report("  no own evidence", lambda R: (R["mu"] > _EPS) & ~R["has_own"])
    report("  no mass level", lambda R: (R["mu"] > _EPS) & ~R["has_mass"])
    report("all peel edges", lambda R: np.ones(len(R["M"]), bool))
