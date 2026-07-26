"""SINGLE-STRAND x CAPTURE, step 4 — is the message CALIBRATED? realized error vs stated variance.

Step 3 showed no single channel owns the damage, so the question is not "which channel" but "are the
messages over-confident". For each channel this measures the RELIABILITY ratio

    z2 = E[(mode - truth)^2] / E[stated variance]        (== 1 if honestly calibrated)

in the channel's own space (log f_c for the measurement streams, lambda for the composition stream), against
the ORACLE. Split capture off/on x region class. This is exactly the DL b_hat^2 quantity, but scored against
truth instead of against the node's own self-solve.

    OMP_NUM_THREADS=1 python scratchpad/ss_cap_4_calib.py
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
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), _EPS))


print(f"{'condition':<42}{'cls':<9}{'chan':<5}{'n':>6}{'|Δmode|':>9}{'sd_stat':>9}"
      f"{'z2':>9}{'bias':>9}{'prec':>10}")
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    uni = cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    T = G + R
    fo_g = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
    fo_p = np.where(T > _EPS, Rp / np.maximum(T, _EPS), np.nan)
    fo_n = np.where(T > _EPS, Rn / np.maximum(T, _EPS), np.nan)
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    is_amb = fp & fn
    mass = np.asarray(cap["mass_global"])
    solvable = np.asarray(cap["solvable"], bool)
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(len(mass))])
    base = np.isfinite(fo_g) & (mass > _EPS) & ~is_amb & solvable

    # channel -> (delivered mode in its own space, truth in the same space, precision)
    # Score on the COMPOSITION axis (gDNA vs RNA-total), never per strand: the diagnostic's oracle pools use
    # the opposite +/- convention from `free_pos`/`free_neg`, and the tilt is a separate DOF anyway. A
    # single-strand node has exactly ONE live RNA arm, so that arm carries the whole RNA claim.
    fo_R = fo_p + fo_n
    lam_true = lg(fo_g) - lg(fo_R)
    mo_R = np.where(fp, uni["mo_p"], uni["mo_n"])
    cm_R = np.where(fp, uni["cm_p"], uni["cm_n"])
    chans = {
        "g": (uni["mo_g"], lg(fo_g), uni["cm_g"], fo_g > 1e-6),
        "R": (mo_R, lg(fo_R), cm_R, fo_R > 1e-6),
        "lam": (uni["lam_msg"], lam_true, uni["c_tau"], (fo_g > 1e-6) & (fo_R > 1e-6)),
    }
    for cl in ("exon", "boundary", "intron"):
        for nm, (mo, tr, pr, live) in chans.items():
            m = base & (cls == cl) & (np.asarray(pr) > _EPS) & live & (np.asarray(mo) > -20)
            if m.sum() < 5:
                continue
            w = mass[m]
            dm = np.asarray(mo)[m] - tr[m]
            v = 1.0 / np.asarray(pr)[m]
            z2 = np.average(dm * dm, weights=w) / np.average(v, weights=w)
            print(f"{cond[5:] if (cl == 'exon' and nm == 'g') else '':<42}{cl if nm == 'g' else '':<9}"
                  f"{nm:<5}{int(m.sum()):>6}{np.average(np.abs(dm), weights=w):>9.3f}"
                  f"{np.sqrt(np.average(v, weights=w)):>9.3f}{z2:>9.1f}"
                  f"{np.average(dm, weights=w):>+9.3f}{np.average(np.asarray(pr)[m], weights=w):>10.1f}")
    print()
