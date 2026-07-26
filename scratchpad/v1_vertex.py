"""VERTEX STUDY — why does an RNA-free intron solve to f_g = 0.982 instead of 1.000?

The chain established in HANDOFF_10 S10.3: psi's reference -> intron f_g = 0.9821 not 1.0000 -> the
continuing share reads 0.294 instead of 0 -> RNA over-claim at the seam -> intron contamination. And
`_rna_arm`'s own docstring says the `+1/2*log(1-f_g)` reference "is what bounds the f_g -> 1 vertex today,
and it is the ONLY thing doing so."

Three candidates for the 1.8 %, and they are separable:
  (a) the REFERENCE   — Beta(1/2,1/2) via +1/2 log f_g + 1/2 log(1-f_g), a proper prior with zero mass AT the vertex;
  (b) the READOUT     — `_posterior_median_fg` takes the grid MEDIAN, which for a one-sided posterior piled
                        against a vertex sits below the mass;
  (c) the GRID        — L = 10 caps f_g at sigma(10) = 0.99995, which is NOT binding at 0.982.

This reconstructs the intron lambda-posterior offline from the captured factory factor and scores each.

    OMP_NUM_THREADS=1 python scratchpad/v1_vertex.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import log_expit

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex_logodds import _logodds_grid  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
cc = cfg.calibration
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
print(f"grid: n_grid={cc.sweep_n_grid}  n_grid_ss={cc.sweep_n_grid_single_strand}  L={cc.sweep_logodds_window}"
      f"   sigma(L) = {1 / (1 + np.exp(-cc.sweep_logodds_window)):.6f}")


def readouts(psi, lam, fg):
    """(median, mean of f_g, mode) from an unnormalized log-posterior on the lambda grid."""
    p = np.exp(psi - psi.max(axis=1, keepdims=True))
    p /= np.maximum(p.sum(axis=1, keepdims=True), 1e-300)
    cw = np.cumsum(p, axis=1)
    med = fg[np.clip((cw < 0.5).sum(axis=1), 0, fg.size - 1)]
    return med, p @ fg, fg[np.argmax(p, axis=1)]


for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cc, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    is_intron = (kind == REGION) & (rt == 1)
    mass = np.asarray(cap["mass_global"])
    ip = cap["intron_prior"]
    if ip is None:
        print(f"\n{cond[5:]}: no intron prior")
        continue
    lp = np.asarray(ip, np.float64)
    lam, fg = _logodds_grid(lp.shape[1], float(cc.sweep_logodds_window))
    ok = is_intron & np.isfinite(fo) & (mass > _EPS) & (np.ptp(lp, axis=1) > _EPS)
    w = mass[ok]
    ref = 0.5 * log_expit(lam)[None, :] + 0.5 * log_expit(-lam)[None, :]
    variants = {
        "factory + reference (SHIPPED)": lp[ok] + ref,
        "factory ONLY (no reference)": lp[ok],
        "reference ONLY (no factory)": np.repeat(ref, int(ok.sum()), axis=0),
    }
    print(f"\n{'=' * 104}\n{cond[5:]}   introns n={int(ok.sum())}   oracle f_g = "
          f"{np.average(fo[ok], weights=w):.4f}   shipped self-solve = {np.average(np.asarray(cap['fg_loc'])[ok], weights=w):.4f}")
    print(f"  {'psi variant':<34}{'median':>10}{'mean f_g':>10}{'mode':>10}{'|med-orc|':>11}{'|mean-orc|':>12}")
    for lab, psi in variants.items():
        med, mean, mode = readouts(psi, lam, fg)
        print(f"  {lab:<34}{np.average(med, weights=w):>10.4f}{np.average(mean, weights=w):>10.4f}"
              f"{np.average(mode, weights=w):>10.4f}{np.average(np.abs(med - fo[ok]), weights=w):>11.4f}"
              f"{np.average(np.abs(mean - fo[ok]), weights=w):>12.4f}")
    # how much lambda-mass does the reference move? the tail decays like exp(-lambda/2)
    psi_s = lp[ok] + ref
    p = np.exp(psi_s - psi_s.max(axis=1, keepdims=True))
    p /= np.maximum(p.sum(axis=1, keepdims=True), 1e-300)
    lam_med = lam[np.clip((np.cumsum(p, axis=1) < 0.5).sum(axis=1), 0, lam.size - 1)]
    print(f"  lambda at the shipped median: {np.average(lam_med, weights=w):.2f}   "
          f"(f_g = {1 / (1 + np.exp(-np.average(lam_med, weights=w))):.4f});   "
          f"the reference contributes {0.5 * np.average(lam_med, weights=w):.2f} nats of decay there")
