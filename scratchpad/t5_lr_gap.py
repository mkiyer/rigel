"""TARGET step 5 — can the LEFT/RIGHT message disagreement price the mismatch where DL is inert?

Step 4: the gDNA channel's error is the TRANSPORT RATIO. `r = rho_tot(dst)/rho_tot(src)` is the capture step
ONLY when the two nodes share a composition -- capture scales both components alike, so
rho_tot ratio = (capture step) x (content ratio). intergenic->boundary (both pure gDNA) is EXACT (0.000);
boundary->exon (mixed vs mostly-gDNA) is off by e^0.851 = 2.3x. That is the imputation premise being false,
which is exactly what M7's DL b_hat^2 prices -- but DL needs the destination's own self-solve as its second
study, and at kappa = 1/2 there is none (tau_own = 0 everywhere), so DL is bit-identically INERT.

THE IDEA. A node with two neighbours receives TWO INDEPENDENT estimates of its own composition. That is a
two-study meta-analysis with NO prior required, exactly the DL setting -- just a different pair. Where the two
sides agree, the premise held; where they disagree, it did not.

This tests whether the left/right gap is actually PREDICTIVE of the node's error before anything is built.

    OMP_NUM_THREADS=1 python scratchpad/t5_lr_gap.py
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
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), 1e-12))


print("Exons with BOTH messages live. G_lr = |log(rho_g^left / rho_g^right)| — the two sides' disagreement.")
print(f"{'condition':<44}{'n':>5}{'ERRshare':>9}|{'G_lr band':>12}{'n':>6}{'ERR':>10}"
      f"{'mwae':>8}{'|Δmode|':>9}{'corr':>7}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    mass = np.asarray(cap["mass_global"])
    ship = np.asarray(cap["f_g"])
    E_g, M = us["E_g"], us["M"]
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    is_exon = (kind == REGION) & (rt == 2)
    ok = np.isfinite(fo) & (mass > _EPS) & np.asarray(cap["solvable"], bool)
    err = np.abs(ship - fo) * mass
    E = err[ok].sum()
    ag, bg = uni["ag"], uni["bg"]
    both = ok & is_exon & (ag > _EPS) & (bg > _EPS)
    if both.sum() < 10:
        continue
    G_lr = np.abs(lg(ag) - lg(bg))
    # the delivered fused gDNA mode error, for the correlation
    dmode = np.abs(lg(np.clip(np.exp(uni["mo_g"]), 1e-12, 1.0)) - lg(fo))
    w = mass[both]
    c = np.corrcoef(G_lr[both], dmode[both])[0, 1]
    first = True
    for lo, hi in ((0, 0.25), (0.25, 1.0), (1.0, 3.0), (3.0, 1e9)):
        m = both & (G_lr >= lo) & (G_lr < hi)
        if m.sum() < 3:
            continue
        wm = mass[m]
        head = (f"{cond[5:]:<44}{int(both.sum()):>5}{err[both].sum() / E:>9.1%}|" if first
                else f"{'':<44}{'':>5}{'':>9}|")
        first = False
        print(f"{head}{f'{lo:g}-{hi:g}':>12}{int(m.sum()):>6}{err[m].sum():>10,.0f}"
              f"{np.average(np.abs(ship - fo)[m], weights=wm):>8.4f}"
              f"{np.average(dmode[m], weights=wm):>9.3f}{c if lo == 0 else np.nan:>7.3f}")
