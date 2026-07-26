"""BOUNDARY LEVEL study — score every candidate estimator of the continuing share `w` at a PEEL edge.

The peel message (exon -> boundary) must decide how much of the exon's RNA CONTINUES unspliced across the
seam.  `w = rho_nu/(rho_nu+rho_mu)`.  Three routes have been tried and all bottomed out on the LEVEL
(rho_nu): the boundary's own self-solve, the far intron, and "no evidence" (v=inf, which silences the
channel).  This script adds and scores a fourth that needs neither a factory nor the strand:

  MASS  --  the boundary's own mass identity `M = rho_g*E_g + rho_nu*E_r` closed with the message's own
            gDNA claim:   rho_nu = (M - tg*E_g)/E_r  ==>  w = (M - tg*E_g)/(A*E_r).
            This is the generic DENSITY DECONVOLUTION (the intron factory's own primitive) with the gDNA
            density prior supplied by the neighbour instead of the intergenic pool.  It exists at EVERY
            seam, including exon|exon and every seam of a low-gDNA library.

  FLUX  --  HEAD's subtraction expressed as a share: w = 1 - rho_mu/A  (M3's u-amplified variance).

Scored against the ORACLE share, and — more importantly — against the DELIVERED f_g the message would
carry, which is what mwae actually measures.

    OMP_NUM_THREADS=1 python scratchpad/w1_level.py [cond ...]
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


def mw(x, w):
    x = np.asarray(x, float)
    m = np.isfinite(x) & (np.asarray(w, float) > 0)
    return float(np.average(x[m], weights=np.asarray(w, float)[m])) if m.any() else np.nan


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
    Rup = pool("mat_uns_pos") + pool("nas_uns_pos")
    Run = pool("mat_uns_neg") + pool("nas_uns_neg")
    Ru = Rup + Run
    fo = np.where(G + Ru > _EPS, G / np.maximum(G + Ru, _EPS), np.nan)

    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    og, op, on = us["og"], us["op"], us["on"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    tau_own, lv = us["tau_own"], us["logvar_tot"]
    SPf = (us["SP_l"], us["SP_r"])
    SNf = (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    NSP = (
        np.asarray(us["spl_n_pos_l"]) + np.asarray(us["spl_n_neg_l"]),
        np.asarray(us["spl_n_pos_r"]) + np.asarray(us["spl_n_neg_r"]),
    )
    n_unspl = np.asarray(us["n_unspl_l"]) + np.asarray(us["n_unspl_r"])
    n_node = np.where(isr, np.asarray(us["n_unspl_l"]), n_unspl)
    rho_lf, rho_rf = uni["rho_lf"], uni["rho_rf"]
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array(
        [CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)]
    )

    rows = []
    for tag, (src_i, dst_face, src_face, df, fw) in (
        ("L", (li, rho_lf, rho_rf, 0, "fwd")),
        ("R", (ri, rho_rf, rho_lf, 1, "bwd")),
    ):
        s = np.clip(src_i, 0, n - 1)
        valid = src_i >= 0
        peel = is_bnd & is_exon[s] & valid
        if not peel.any():
            continue
        rg, rp, rn = us[f"{fw}_g"], us[f"{fw}_p"], us[f"{fw}_n"]
        pg, pp = us[f"{fw}_pg"], us[f"{fw}_pp"]
        pn = us[f"{fw}_pn"]
        framed = valid & (src_face[s] > _EPS) & (dst_face > _EPS)
        r = np.where(framed, dst_face / np.maximum(src_face[s], _EPS), np.where(valid, 1.0, 0.0))
        tg = rg[s] * r
        Ap, An = rp[s] * r, rn[s] * r
        A = Ap + An
        mu_p = np.where(SPf[df] > _EPS, SPf[df] / np.maximum(ESP[df], _EPS), 0.0)
        mu_n = np.where(SNf[df] > _EPS, SNf[df] / np.maximum(ESP[df], _EPS), 0.0)
        mu = mu_p + mu_n
        s2t = lv + lv[s]
        v_g = np.where(pg[s] > 0, 1.0 / np.maximum(pg[s], _EPS) + s2t, np.inf)
        v_A = np.where((pp[s] + pn[s]) > 0, 1.0 / np.maximum(pp[s] + pn[s], _EPS) + s2t, np.inf)
        v_mu = np.where(NSP[df] > 0, 1.0 / np.maximum(NSP[df], _EPS), np.inf)

        # ── truth ──
        nu_true = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), 0.0)
        w_true = np.where(nu_true + mu > _EPS, nu_true / np.maximum(nu_true + mu, _EPS), 1.0)

        # ── estimators (node-level, both strands pooled: one share per edge) ──
        w_flux = np.where(A > _EPS, np.maximum(A - mu, 0.0) / np.maximum(A, _EPS), 1.0)
        v_flux_nu = np.where(  # M3: u^2 v_A + (u-1)^2 v_mu on the DIFFERENCE
            (A - mu) > _EPS, (A / np.maximum(A - mu, _EPS)) ** 2 * v_A
            + (mu / np.maximum(A - mu, _EPS)) ** 2 * v_mu, np.inf)
        nu_flux = np.maximum(A - mu, 0.0)

        phi = np.where(M > _EPS, tg * E_g / np.maximum(M, _EPS), np.inf)
        nu_mass = np.where(M > _EPS, np.maximum(M - tg * E_g, 0.0) / np.maximum(E_r, _EPS), 0.0)
        w_mass = np.where(A > _EPS, np.minimum(nu_mass / np.maximum(A, _EPS), 1.0), 1.0)
        v_mass_nu = np.where(  # (1/n_M + phi^2 v_g)/(1-phi)^2
            phi < 1.0, (1.0 / np.maximum(n_node, _EPS) + phi * phi * v_g)
            / np.maximum((1.0 - phi) ** 2, _EPS), np.inf)

        nu_self = op + on
        w_self = np.where(nu_self + mu > _EPS, nu_self / np.maximum(nu_self + mu, _EPS), 1.0)
        has_own = tau_own > _EPS

        # far node across the seam (opposite the exon source)
        far = np.where(df == 0, ri, li)
        fx = np.clip(far, 0, n - 1)
        far_ok = (far >= 0) & (tau_own[fx] > _EPS) & (us["rho_node0"][fx] > _EPS) & (
            us["rho_node0"] > _EPS)
        r_far = np.where(far_ok, us["rho_node0"] / np.maximum(us["rho_node0"][fx], _EPS), 1.0)
        nu_far = np.where(far_ok, (op[fx] + on[fx]) * r_far, 0.0)
        w_far = np.where(nu_far + mu > _EPS, nu_far / np.maximum(nu_far + mu, _EPS), 1.0)

        # inverse-variance combine of FLUX and MASS on the LEVEL (log space)
        pf = np.where(np.isfinite(v_flux_nu) & (nu_flux > _EPS), 1.0 / np.maximum(v_flux_nu, _EPS), 0.0)
        pm = np.where(np.isfinite(v_mass_nu) & (nu_mass > _EPS), 1.0 / np.maximum(v_mass_nu, _EPS), 0.0)
        lf = np.log(np.maximum(nu_flux, _EPS))
        lm = np.log(np.maximum(nu_mass, _EPS))
        nu_comb = np.where(pf + pm > _EPS, np.exp((pf * lf + pm * lm) / np.maximum(pf + pm, _EPS)), 0.0)
        w_comb = np.where(nu_comb + mu > _EPS, nu_comb / np.maximum(nu_comb + mu, _EPS), 1.0)

        sel = peel & np.isfinite(fo) & (M > _EPS) & (A > _EPS) & np.asarray(cap["solvable"], bool)
        cl_far = np.array([cls[j] if j >= 0 else "edge" for j in far])
        rows.append(dict(
            tag=tag, sel=sel, w_true=w_true, mu=mu, A=A, tg=tg, M=M, E_g=E_g, E_r=E_r,
            fo=fo, has_own=has_own, far_ok=far_ok, cl_far=cl_far, phi=phi, n_node=n_node,
            ests={"flux": (w_flux, v_flux_nu), "mass": (w_mass, v_mass_nu),
                  "self": (w_self, np.full(n, np.nan)), "far": (w_far, np.full(n, np.nan)),
                  "comb": (w_comb, np.full(n, np.nan))},
        ))

    print(f"\n{'=' * 118}\n{cond[5:]}\n{'=' * 118}")
    # aggregate over both directions
    def agg(mask_fn, label):
        tot = 0
        line = {}
        for R in rows:
            m = R["sel"] & mask_fn(R)
            if not m.any():
                continue
            w = R["M"][m]
            tot += int(m.sum())
            for k, (we, _v) in R["ests"].items():
                dw = np.abs(we[m] - R["w_true"][m])
                # delivered f_g under this share
                fgd = np.where(
                    R["tg"][m] * R["E_g"][m] + R["A"][m] * we[m] * R["E_r"][m] > _EPS,
                    R["tg"][m] * R["E_g"][m]
                    / np.maximum(R["tg"][m] * R["E_g"][m] + R["A"][m] * we[m] * R["E_r"][m], _EPS),
                    1.0)
                dfg = np.abs(fgd - R["fo"][m])
                a = line.setdefault(k, [0.0, 0.0, 0.0])
                a[0] += float(np.sum(dw * w)); a[1] += float(np.sum(dfg * w)); a[2] += float(np.sum(w))
            a = line.setdefault("_true", [0.0, 0.0, 0.0])
            fgt = np.where(
                R["tg"][m] * R["E_g"][m] + R["A"][m] * R["w_true"][m] * R["E_r"][m] > _EPS,
                R["tg"][m] * R["E_g"][m]
                / np.maximum(R["tg"][m] * R["E_g"][m] + R["A"][m] * R["w_true"][m] * R["E_r"][m], _EPS),
                1.0)
            a[0] += 0.0; a[1] += float(np.sum(np.abs(fgt - R["fo"][m]) * w)); a[2] += float(np.sum(w))
        if not line:
            return
        wt = mw(np.concatenate([R["w_true"][R["sel"] & mask_fn(R)] for R in rows if (R["sel"] & mask_fn(R)).any()]),
                np.concatenate([R["M"][R["sel"] & mask_fn(R)] for R in rows if (R["sel"] & mask_fn(R)).any()]))
        out = f"  {label:<24} n={tot:<6} w_true={wt:.3f} |"
        for k in ("flux", "mass", "self", "far", "comb", "_true"):
            if k in line and line[k][2] > 0:
                out += f" {k}: dw={line[k][0] / line[k][2]:.3f} dfg={line[k][1] / line[k][2]:.3f} |"
        print(out)

    print("  (dw = |w_hat - w_true|, dfg = |delivered f_g - oracle f_g|; both mass-weighted)")
    agg(lambda R: np.ones(n, bool), "ALL peel edges")
    agg(lambda R: R["cl_far"] == "intron", "  far = intron")
    agg(lambda R: R["cl_far"] == "exon", "  far = exon")
    agg(lambda R: R["has_own"], "  boundary tau>0")
    agg(lambda R: ~R["has_own"] & ~R["far_ok"], "  CASE 3 (no evidence)")
    agg(lambda R: R["mu"] <= _EPS, "  no spliced flux")
    agg(lambda R: R["mu"] > _EPS, "  junction-bearing")
