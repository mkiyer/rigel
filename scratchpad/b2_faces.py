"""BOUNDARY STUDY step 2 — is the TWO-SIDED reframe actually implemented, and does it recover the capture step?

OWNER'S MODEL. A junction boundary has two sides with DIFFERENT frames:
  * the side facing the EXON  — its total density INCLUDES the spliced fragments, and that content matches the
    exon's content, so   rho_face_exon(b) / rho_tot(exon)  ==  e(b)/e(exon)   = the pure capture step;
  * the side facing the INTRON — its total density EXCLUDES the spliced (they have spliced away), and that
    content matches the intron's, so  rho_face_intron(b) / rho_tot(intron) == e(b)/e(intron).
Getting this right is what makes a rho_tot RATIO equal a CAPTURE STEP instead of capture x content-ratio —
which is exactly the transport-ratio defect measured earlier (|log r/r_true| = 0.851 on boundary->exon).

`node_total_density`'s docstring already states this intent. But `bp_solver` selects the with-spliced density
by `accept_l/accept_r = (SP+SN on that face) > 0` — the face that CARRIES spliced mass — not by the face that
FACES AN EXON. This measures (a) whether those two coincide, and (b) how each candidate selector scores
against the TRUE capture step, taken from oracle gDNA densities (gDNA content is uniform, so a ratio of true
gDNA densities IS a capture ratio).

    OMP_NUM_THREADS=1 python scratchpad/b2_faces.py
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
from rigel.calibration.node_geometry import _node_region_type, node_total_density  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), 1e-12))


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    E_g = us["E_g"]
    n = len(E_g)
    rho_g_true = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)  # == capture, up to a constant
    li, ri = us["left"], us["right"]
    rt, _ = _node_region_type(chain, ra)
    CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
    mass = np.asarray(cap["mass_global"])
    # the model's own face densities, at the SHIPPED belief
    ru, rw = node_total_density(chain, geo, np.asarray(cap["f_g"]))
    rs = rw - ru
    SPl = np.asarray(geo.spliced_pos_left, float) + np.asarray(geo.spliced_neg_left, float)
    SPr = np.asarray(geo.spliced_pos_right, float) + np.asarray(geo.spliced_neg_right, float)
    acc_l, acc_r = SPl > _EPS, SPr > _EPS
    is_bnd = ~isr
    junc = is_bnd & ((SPl + SPr) > _EPS)

    # (a) do "the face carrying spliced" and "the face facing an exon" coincide?
    ex_l = np.zeros(n, bool)
    ex_r = np.zeros(n, bool)
    for i in np.flatnonzero(is_bnd):
        ex_l[i] = li[i] >= 0 and cls[li[i]] == "exon"
        ex_r[i] = ri[i] >= 0 and cls[ri[i]] == "exon"
    both = junc & (li >= 0) & (ri >= 0)
    agree = ((acc_l == ex_l) & (acc_r == ex_r))[both]
    print(f"\n{'=' * 112}\n{cond[5:]}   junction boundaries with two flanks: {int(both.sum())}\n{'=' * 112}")
    print(f"  selector agreement — 'carries spliced' == 'faces an exon': {agree.mean():.1%}")
    print(f"    of junction boundaries: exon on BOTH flanks {np.mean((ex_l & ex_r)[both]):.1%}, "
          f"exon on ONE flank {np.mean((ex_l ^ ex_r)[both]):.1%}, NEITHER {np.mean((~ex_l & ~ex_r)[both]):.1%}")
    print(f"    spliced deposited on BOTH faces {np.mean((acc_l & acc_r)[both]):.1%}, "
          f"ONE face {np.mean((acc_l ^ acc_r)[both]):.1%}")

    # (b) score each selector: does rho_face(b)/rho_tot(nbr) recover the TRUE capture step e(b)/e(nbr)?
    print(f"\n  |log( rho_face(b)/rho_tot(nbr) ) − log( e(b)/e(nbr) )|   mass-weighted, lower is better")
    print(f"  {'selector':<34}{'-> EXON flank':>15}{'-> INTRON flank':>17}{'ALL':>10}")
    # ── EXON-SIDE VARIANTS: how should the spliced be folded into the exon-facing total density? ──
    Mn = us["M"]
    E_r = us["E_r"]
    Sl = SPl; Sr = SPr
    ESl = np.maximum(np.asarray(geo.eff_spl_left, float), _EPS)
    ESr = np.maximum(np.asarray(geo.eff_spl_right, float), _EPS)
    S_tot = Sl + Sr
    variants = {
        "(a) shipped  rho_u + S/E_spl": (ru + np.where(acc_l, rs, 0.0), ru + np.where(acc_r, rs, 0.0)),
        "(b) rho_u + S/E_r": (ru + np.where(acc_l, Sl / np.maximum(E_r, _EPS), 0.0),
                              ru + np.where(acc_r, Sr / np.maximum(E_r, _EPS), 0.0)),
        "(c) rho_u * (M+S)/M  [same blend]": (ru * np.where(acc_l, 1.0 + Sl / np.maximum(Mn, _EPS), 1.0),
                                              ru * np.where(acc_r, 1.0 + Sr / np.maximum(Mn, _EPS), 1.0)),
        "(d) unspliced only": (ru, ru),
    }
    print(f"\n  EXON-FACING side only — how to fold the spliced in:")
    print(f"  {'variant':<38}{'|log err| vs capture step':>26}")
    for lab, (fl2, fr2) in variants.items():
        acc = []
        for nbr, face in ((li, fl2), (ri, fr2)):
            s_ = np.clip(nbr, 0, n - 1)
            m = (both & (nbr >= 0) & (cls[s_] == "exon") & np.isfinite(rho_g_true) & (rho_g_true > 0)
                 & np.isfinite(rho_g_true[s_]) & (rho_g_true[s_] > 0) & (face > _EPS))
            if not m.any():
                continue
            acc.append((np.abs((lg(face[m]) - lg(ru[s_[m]])) - (lg(rho_g_true[m]) - lg(rho_g_true[s_[m]]))),
                        mass[m]))
        if acc:
            d = np.concatenate([a[0] for a in acc]); w = np.concatenate([a[1] for a in acc])
            print(f"  {lab:<38}{np.average(d, weights=w):>26.3f}")
    print(f"  median S/M at junction boundaries: {np.median((S_tot/np.maximum(Mn,_EPS))[both]):.3f}"
          f"   median E_spl/E_r: {np.median(((ESl+ESr)/np.maximum(E_r,_EPS))[both]):.3f}")

    sel = {
        "shipped (carries spliced)": (acc_l, acc_r),
        "structural (faces an exon)": (ex_l, ex_r),
        "always WITH spliced": (np.ones(n, bool), np.ones(n, bool)),
        "never (unspliced only)": (np.zeros(n, bool), np.zeros(n, bool)),
    }
    for lab, (wl, wr) in sel.items():
        fl = ru + np.where(wl, rs, 0.0)  # density on the LEFT-facing side
        fr = ru + np.where(wr, rs, 0.0)
        out = {}
        for tag, want in (("exon", "exon"), ("intron", "intron"), ("all", None)):
            acc = []
            for nbr, face in ((li, fl), (ri, fr)):
                m = (both & (nbr >= 0) & np.isfinite(rho_g_true) & (rho_g_true > 0)
                     & np.isfinite(rho_g_true[np.clip(nbr, 0, n - 1)])
                     & (rho_g_true[np.clip(nbr, 0, n - 1)] > 0) & (face > _EPS))
                if want is not None:
                    m = m & (cls[np.clip(nbr, 0, n - 1)] == want)
                if not m.any():
                    continue
                s = np.clip(nbr, 0, n - 1)[m]
                model = lg(face[m]) - lg(ru[s])  # the neighbour is a REGION ⇒ no spliced there
                true = lg(rho_g_true[m]) - lg(rho_g_true[s])
                acc.append((np.abs(model - true), mass[m]))
            if acc:
                d = np.concatenate([a[0] for a in acc])
                w = np.concatenate([a[1] for a in acc])
                out[tag] = np.average(d, weights=w)
        print(f"  {lab:<34}{out.get('exon', np.nan):>15.3f}{out.get('intron', np.nan):>17.3f}"
              f"{out.get('all', np.nan):>10.3f}")
