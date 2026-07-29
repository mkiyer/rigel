"""E2 step 0b — the CLEAN per-face vs node-level frame test (finding 3, absolute footing).

The first pass (`enrichment_relay_validation.py`) compared a TOTAL-density estimator to a gDNA-ONLY truth,
which introduced a systematic -0.52-nat offset (a boundary crossing is more gDNA-ish than its flanking
region -- mature cannot cross a junction, and the crossing eff-length frame differs from the contained one).
That offset is COMMON-MODE to the FACE and NODE rows, so it cancels in the FACE-vs-NODE ordering but
contaminates the ABSOLUTE centering. This pass removes it by comparing TOTAL-to-TOTAL, using the oracle's
per-SIDE gDNA AND RNA (`OracleTruth.override_masses`) to build a genuine per-face TOTAL density truth.

THE TEST. For each boundary B and each face (facing its left / right region neighbour), with node-level
oracle f_g held fixed so only the GEOMETRY frame varies:

    FACE estimator   rho = m_face * [ f_g/E_g_face + (1-f_g)/E_r_face ]        the face's own mass + eff-lens
    NODE estimator   rho = (m_l+m_r) * [ f_g/(E_gl+E_gr) + (1-f_g)/(E_rl+E_rr) ]   pooled node-level
    TRUE face        rho = g_face/E_g_face + r_face/E_r_face                   oracle gDNA + RNA on THAT face

    residual = log(estimator) - log(TRUE face)

Because the accumulator mass is an exact partition of the oracle reads (m_face = g_face + r_face), a scale
offset cannot arise; the only residual is the composition-FRAME error (node f_g applied on a face whose true
composition differs -- the bounding-lemma tolerance). GATE: FACE residual centres near 0 (the face frame is
absolutely right); NODE residual is larger and grows with capture (the pooled frame is wrong at a cliff).

Side mapping (verified): boundary B between region L=left[B] and R=right[B]. B's LEFT face faces L and IS
L's RIGHT side (override_masses.mass_*_right[L]); B's RIGHT face faces R and IS R's LEFT side
(override_masses.mass_*_left[R]).

Usage:  OMP_NUM_THREADS=1 python scripts/debug/perface_frame_validation.py
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from oracle import OracleTruth  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-12
_TINY_EFF = 1.0
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
WORK = Path("/tmp/rigel_perface")

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
CONDS = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

# (capstate, frame) -> [residual];  also a paired (capstate) -> [|FACE| - |NODE|] for the head-to-head
acc: dict[tuple, list] = defaultdict(list)
paired: dict[str, list] = defaultdict(list)
n_scen = 0

for cond in CONDS:
    if "gdna_none" in cond:
        continue  # need gDNA present for a meaningful density truth
    capst = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    inp = _scan_and_truth(SUITE, cond, index, cfg, WORK, SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, geom, cap = dbg["chain"], dbg["geometry"], dbg["capture"]

    # the per-SIDE oracle truth (re-partition the BAM by read origin; reuse the cached full payload)
    orc = OracleTruth.from_bam(
        str(SUITE / cond / "sim_oracle.bam"), index, cfg, WORK, cond,
        full_payload=inp["payload"], boundary_mass_tol=0.5,  # matches _scan_and_truth
    )
    om = orc.override_masses(ra)
    g_left_reg, r_left_reg = np.asarray(om["mass_gdna_left"]), np.asarray(om["mass_rna_left"])
    g_right_reg, r_right_reg = np.asarray(om["mass_gdna_right"]), np.asarray(om["mass_rna_right"])
    n_scen += 1

    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    ridx = np.asarray(chain.ref_idx, np.int64)
    ml, mr = np.asarray(geom.mass_left, float), np.asarray(geom.mass_right, float)
    egl, egr = np.asarray(geom.eff_gdna_left, float), np.asarray(geom.eff_gdna_right, float)
    erl, err_ = np.asarray(geom.eff_rna_left, float), np.asarray(geom.eff_rna_right, float)

    # node-level oracle f_g at the boundary (from the summed boundary pools)
    bp = inp["boundary_pools"]
    Gb = np.asarray(bp["gdna_pos"], float) + np.asarray(bp["gdna_neg"], float)
    Rb = (np.asarray(bp["mat_uns_pos"], float) + np.asarray(bp["mat_uns_neg"], float)
          + np.asarray(bp["nas_uns_pos"], float) + np.asarray(bp["nas_uns_neg"], float))

    for b in np.where(~is_reg)[0]:
        bi = ridx[b]
        fgo = Gb[bi] / (Gb[bi] + Rb[bi]) if (Gb[bi] + Rb[bi]) > _EPS else np.nan
        if not np.isfinite(fgo):
            continue
        egl_b, egr_b, erl_b, err_b = egl[b], egr[b], erl[b], err_[b]
        ml_b, mr_b = ml[b], mr[b]
        node_m, node_eg, node_er = ml_b + mr_b, egl_b + egr_b, erl_b + err_b
        if node_m <= _EPS or node_eg <= _TINY_EFF:
            continue
        node_est = node_m * (fgo / node_eg + (1.0 - fgo) / node_er)

        # LEFT face faces region L = left[b]; that face IS L's RIGHT side
        L = left[b]
        if L >= 0 and is_reg[L] and ml_b > _EPS and egl_b > _TINY_EFF:
            rL = ridx[L]
            g_f, r_f = g_right_reg[rL], r_right_reg[rL]
            true_face = g_f / egl_b + r_f / erl_b
            face_est = ml_b * (fgo / egl_b + (1.0 - fgo) / erl_b)
            if true_face > _EPS and face_est > _EPS and node_est > _EPS:
                rF = float(np.log(face_est) - np.log(true_face))
                rN = float(np.log(node_est) - np.log(true_face))
                acc[(capst, "FACE")].append(rF)
                acc[(capst, "NODE")].append(rN)
                paired[capst].append(abs(rF) - abs(rN))

        # RIGHT face faces region R = right[b]; that face IS R's LEFT side
        R = right[b]
        if R >= 0 and is_reg[R] and mr_b > _EPS and egr_b > _TINY_EFF:
            rR = ridx[R]
            g_f, r_f = g_left_reg[rR], r_left_reg[rR]
            true_face = g_f / egr_b + r_f / err_b
            face_est = mr_b * (fgo / egr_b + (1.0 - fgo) / err_b)
            if true_face > _EPS and face_est > _EPS and node_est > _EPS:
                rF = float(np.log(face_est) - np.log(true_face))
                rN = float(np.log(node_est) - np.log(true_face))
                acc[(capst, "FACE")].append(rF)
                acc[(capst, "NODE")].append(rN)
                paired[capst].append(abs(rF) - abs(rN))


def _stat(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    if not v.size:
        return 0, np.nan, np.nan, np.nan, np.nan
    return v.size, np.median(v), v.mean(), v.std(), np.median(np.abs(v))


print(f"\nPER-FACE FRAME VALIDATION (total-to-total) -- {n_scen} gDNA-present conditions.")
print("residual = log(estimator density) - log(TRUE face total density).  node f_g held fixed.")
print("GATE: FACE centres on ~0 (absolutely right); NODE deviates and grows with capture.\n")
print(f"      {'frame':<6} {'capture':<7} | {'n':>7} {'median':>9} {'mean':>9} {'sd':>8} {'med|.|':>8}")
for capst in ("OFF", "ON ", "VSTR"):
    for fr in ("FACE", "NODE"):
        n, md, mn, sd, mad = _stat(acc.get((capst, fr), []))
        if n:
            print(f"      {fr:<6} {capst:<7} | {n:>7} {md:>9.3f} {mn:>9.3f} {sd:>8.3f} {mad:>8.3f}")
    print()

print("=" * 84)
print("HEAD-TO-HEAD (paired, same edge):  |FACE residual| - |NODE residual|  (negative ==> FACE wins)")
print(f"      {'capture':<7} | {'n':>7} {'median':>9} {'mean':>9} {'% FACE better':>14}")
for capst in ("OFF", "ON ", "VSTR"):
    v = np.asarray(paired.get(capst, []))
    v = v[np.isfinite(v)]
    if v.size:
        print(f"      {capst:<7} | {v.size:>7} {np.median(v):>9.3f} {v.mean():>9.3f} {100.0*np.mean(v < 0):>13.1f}%")
print()
print("VERDICT: if FACE median|.| < NODE median|.| and FACE centres near 0 at ON/VSTR, per-face is the")
print("absolutely-correct frame and the -0.52 offset was indeed a total-vs-gDNA artifact of the first pass.")
