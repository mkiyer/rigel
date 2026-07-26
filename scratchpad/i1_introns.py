"""INTRON STUDY step 1 — is the factory's precision HONEST, and what stops an intron self-solving?

OWNER'S MODEL. The density deconvolution peels gDNA against a measured background: it delivers the LAMBDA
axis (gDNA vs RNA-total) and a precision, but NOT the tilt. So:
  * a SINGLE-STRAND intron has its tilt structurally locked ⇒ lambda (factory) + theta (structure) is a
    COMPLETE self-solve; it should need no messages at all;
  * an AMBIG intron knows its gDNA/RNA split but not which strand the RNA sits on ⇒ it needs messages for
    theta ONLY, never for lambda.

Known already (suite_dissect / pass0_error_table): introns are the ONLY class where MESSAGES MAKE THINGS
WORSE (self 0.0103 -> solved 0.0133 suite-wide), and 90.2 % of intron error sits on the most-confident
quartile. And tau_factory spans TEN ORDERS OF MAGNITUDE across scenarios (0.215 ... 1.2e10), which is the
owner's suspicion of NB over-confidence.

This measures, per scenario and split single/AMBIG:
  - the self-solve vs the solved answer (does message passing help or hurt?)
  - tau_factory's CALIBRATION: z2 = E[(lam_self - lam_true)^2] / E[1/tau_factory]   (1.0 = honest)

    OMP_NUM_THREADS=1 python scratchpad/i1_introns.py
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
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna5_ss_0.50_nrna_present_capture_on",
    "gdna_none_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def lgt(x):
    x = np.clip(np.asarray(x, np.float64), 1e-6, 1 - 1e-6)
    return np.log(x / (1 - x))


def mw(x, w):
    return np.average(x, weights=w) if np.sum(w) > 0 else np.nan


print(f"{'condition':<44}{'dof':<7}{'n':>5}{'reads':>10}{'ERR':>9}{'self':>8}{'solved':>8}"
      f"{'Δmsg':>8}|{'orcFg':>7}{'selfFg':>7}|{'med tau':>11}{'z2 (fac)':>10}{'z2 self':>9}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    mass = np.asarray(cap["mass_global"])
    ship, self_fg = np.asarray(cap["f_g"]), np.asarray(cap["fg_loc"])
    tau = np.asarray(cap["_tau0_lam"])
    vg_loc = np.asarray(cap["vg_loc"])  # the message-free local posterior Var(log f_g)
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb = fp & fn
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    is_intron = (kind == REGION) & (rt == 1)
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool) & is_intron
    first = True
    for lab, dm in (("single", ok & ~amb), ("AMBIG", ok & amb)):
        if dm.sum() < 3:
            continue
        w = mass[dm]
        e_self = np.abs(self_fg - fo)[dm]
        e_ship = np.abs(ship - fo)[dm]
        # factory calibration on the lambda axis
        dlam = lgt(self_fg[dm]) - lgt(fo[dm])
        t = tau[dm]
        good = t > _EPS
        z2 = (np.average(dlam[good] ** 2, weights=w[good])
              / np.average(1.0 / t[good], weights=w[good])) if good.any() else np.nan
        # the same, but against the node's OWN posterior variance (what it actually reports)
        vgl = vg_loc[dm]
        g2 = np.isfinite(vgl) & (vgl > _EPS)
        dlf = np.log(np.clip(self_fg[dm], 1e-9, 1)) - np.log(np.clip(fo[dm], 1e-9, 1))
        z2s = (np.average(dlf[g2] ** 2, weights=w[g2]) / np.average(vgl[g2], weights=w[g2])) if g2.any() else np.nan
        head = cond[5:] if first else ""
        first = False
        print(f"{head:<44}{lab:<7}{int(dm.sum()):>5}{w.sum():>10,.0f}{(e_ship * w).sum():>9,.0f}"
              f"{mw(e_self, w):>8.4f}{mw(e_ship, w):>8.4f}{mw(e_ship, w) - mw(e_self, w):>+8.4f}|"
              f"{mw(fo[dm], w):>7.4f}{mw(self_fg[dm], w):>7.4f}|"
              f"{np.median(t[good]) if good.any() else np.nan:>11.3g}{z2:>10.2f}{z2s:>9.2f}")
