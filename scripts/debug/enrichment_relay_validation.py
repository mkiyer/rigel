"""E2 step 0 — VALIDATE the general concept: a message scaled by the total-density ratio lands in the
destination's frame, INCLUDING when it has relayed through data-free / tiny regions (owner, 2026-07-23).

A tiny region (length < fragment length) cannot contain a fragment -- every overlapping fragment spills into
the flanking boundaries -- so it has ~zero contained mass, zero precision, and messages RELAY through it. The
owner's claim: enrichment ratios (ratios of total densities) are NEVERTHELESS valid across such a relay,
because the ratio is formed between the two FRAMED nodes on either side of the gap, and the tiny region is
transparent to it.

TRUTH BASIS. gDNA content is genomically ~uniform, so the oracle gDNA density rho_g_oracle(x) = G(x)/E_g(x)
is a clean probe of the capture efficiency e(x): rho_g_oracle(C)/rho_g_oracle(A) = e(C)/e(A), the TRUE
enrichment ratio. (Restricted to gDNA-PRESENT scenarios, where rho_g_oracle is meaningful.) An estimator's
enrichment step delta = log rho_est(C) - log rho_est(A) is VALID iff delta ~ log[e(C)/e(A)].

THREE ESTIMATORS of the node's total density:
  BLENDED oracle f_g   M*[f_g/E_g + (1-f_g)/E_r]   with ORACLE composition   -- the framework's ideal
  BLENDED solved f_g   ... with the pass-0 SOLVED composition                -- the real-world quantity
  gDNA-frame M/E_g     belief-free                                           -- the E0 fallback

THE GATE (owner's concept):
  (1) residual = delta - true ~ 0  on the measurable edges  -> the ratio measures the enrichment;
  (2) RELAY residual ~ DIRECT residual                      -> the tiny region is transparent.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/enrichment_relay_validation.py
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-12
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
# a node is DATA-FREE (relays, no frame of its own) iff it has ~no mass OR cannot contain a fragment (eff<1)
_TINY_EFF = 1.0

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
CONDS = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

# (capstate, gaptype, estimator) -> [residual]      residual = delta_est - delta_true
acc: dict[tuple, list] = defaultdict(list)
# separately, the per-face vs node-level frame question at boundary endpoints
face_acc: dict[tuple, list] = defaultdict(list)
n_scen = 0

ESTS = ("BLENDED oracle", "BLENDED solved", "gDNA-frame M/E_g")

for cond in CONDS:
    # true enrichment probe needs gDNA present; gdna_none has ~none -> skip for the truth comparison
    if "gdna_none" in cond:
        continue
    capst = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    strand = "str" if "ss_0.99" in cond else "uns"
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, geom, cap = dbg["chain"], dbg["geometry"], dbg["capture"]
    n_scen += 1
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    order = [int(x) for x in np.asarray(chain.order)]

    # node-level pooled geometry (both faces for a boundary)
    ml, mr = np.asarray(geom.mass_left, float), np.asarray(geom.mass_right, float)
    mass = np.where(is_reg, ml, ml + mr)
    egl, egr = np.asarray(geom.eff_gdna_left, float), np.asarray(geom.eff_gdna_right, float)
    erl, err_ = np.asarray(geom.eff_rna_left, float), np.asarray(geom.eff_rna_right, float)
    E_g = np.where(is_reg, egl, egl + egr)
    E_r = np.where(is_reg, erl, erl + err_)
    f_g_solved = np.clip(np.asarray(cap["f_g"], float), 0.0, 1.0)

    # oracle per-node masses -> oracle f_g (unspliced) + the true-enrichment gDNA density probe
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    Gtot, Rtot = Gp + Gn, Rp + Rn
    with np.errstate(divide="ignore", invalid="ignore"):
        f_g_oracle = np.where(Gtot + Rtot > _EPS, Gtot / np.maximum(Gtot + Rtot, _EPS), 1.0)
    rho_g_true = np.where(E_g > _EPS, Gtot / np.maximum(E_g, _EPS), np.nan)  # ~ e(x), the truth

    # the three estimators of node total density
    est = {
        "BLENDED oracle": mass * (f_g_oracle / np.maximum(E_g, _EPS) + (1.0 - f_g_oracle) / np.maximum(E_r, _EPS)),
        "BLENDED solved": mass * (f_g_solved / np.maximum(E_g, _EPS) + (1.0 - f_g_solved) / np.maximum(E_r, _EPS)),
        "gDNA-frame M/E_g": mass / np.maximum(E_g, _EPS),
    }

    # a node has a FRAME iff it carries data: real mass AND it can contain/deposit a fragment
    framed = (mass > _EPS) & (E_g > _TINY_EFF) & np.isfinite(rho_g_true) & (rho_g_true > _EPS)
    datafree = ~framed

    # walk the chain in genomic order; for each framed node, find its nearest framed neighbour to the LEFT
    # via chain.left (skipping data-free nodes). Tag the edge DIRECT (adjacent) or RELAY (>=1 skip).
    for i in order:
        if not framed[i]:
            continue
        j = left[i]
        skipped = 0
        while j >= 0 and datafree[j]:
            j = left[j]
            skipped += 1
        if j < 0 or not framed[j]:
            continue
        gap = "DIRECT" if skipped == 0 else "RELAY "
        true_step = float(np.log(rho_g_true[i]) - np.log(rho_g_true[j]))
        for en in ESTS:
            a, c = est[en][j], est[en][i]
            if a > _EPS and c > _EPS:
                delta = float(np.log(c) - np.log(a))
                acc[(capst, gap, en)].append(delta - true_step)

    # PER-FACE vs NODE-LEVEL, at boundary endpoints only. For a boundary C with a framed region neighbour A,
    # is the FACE facing A a better frame than the pooled node-level? Compare each estimator's residual using
    # (i) node-level pooled E and (ii) the single face of C facing A.
    for i in order:
        if is_reg[i] or not framed[i]:
            continue  # boundary endpoints
        for side, nbr_arr, face in (("L", left, 0), ("R", right, 1)):
            a = nbr_arr[i]
            if a < 0 or not framed[a] or not is_reg[a]:
                continue
            # face-level total density of the boundary on the side facing the region neighbour
            m_face = ml[i] if face == 0 else mr[i]
            eg_face = egl[i] if face == 0 else egr[i]
            er_face = erl[i] if face == 0 else err_[i]
            if m_face <= _EPS or eg_face <= _TINY_EFF:
                continue
            fgo = f_g_oracle[i]
            rho_face = m_face * (fgo / max(eg_face, _EPS) + (1.0 - fgo) / max(er_face, _EPS))
            rho_node = est["BLENDED oracle"][i]
            rho_reg = est["BLENDED oracle"][a]
            # the boundary's own face gDNA density truth (face gDNA mass / face E_g) vs the region's
            g_face = (Gp[i] + Gn[i]) * (m_face / max(mass[i], _EPS))  # apportion node gDNA to the face by mass
            true_face = float(np.log(max(g_face / max(eg_face, _EPS), _EPS)) - np.log(rho_g_true[a]))
            if rho_reg > _EPS:
                if rho_face > _EPS:
                    face_acc[(capst, "FACE ")].append(float(np.log(rho_face) - np.log(rho_reg)) - true_face)
                if rho_node > _EPS:
                    face_acc[(capst, "NODE ")].append(float(np.log(rho_node) - np.log(rho_reg)) - true_face)


def _stat(v):
    v = np.asarray(v)
    v = v[np.isfinite(v)]
    return v.size, (np.median(v) if v.size else np.nan), (v.mean() if v.size else np.nan), (v.std() if v.size else np.nan)


print(f"\nENRICHMENT-RELAY VALIDATION -- {n_scen} gDNA-present conditions.")
print("residual = (estimator enrichment step) - (TRUE enrichment step from oracle gDNA density).")
print("GATE: residual ~ 0, and RELAY ~ DIRECT (the tiny region is transparent to the ratio).\n")
for capst in ("OFF", "ON ", "VSTR"):
    rows = [(g, e) for g in ("DIRECT", "RELAY ") for e in ESTS if acc.get((capst, g, e))]
    if not rows:
        continue
    print(f"  --- capture {capst.strip()} ---")
    print(f"      {'gap':<7} {'estimator':<18} | {'n':>6} {'median':>9} {'mean':>9} {'sd':>8}")
    for g, e in rows:
        n, md, mn, sd = _stat(acc[(capst, g, e)])
        if n:
            print(f"      {g:<7} {e:<18} | {n:>6} {md:>9.3f} {mn:>9.3f} {sd:>8.3f}")
    print()

print("=" * 92)
print("PER-FACE vs NODE-LEVEL boundary frame (residual vs the region neighbour, oracle f_g):")
print(f"      {'frame':<7} {'capture':<7} | {'n':>6} {'median':>9} {'mean':>9} {'sd':>8}")
for capst in ("OFF", "ON ", "VSTR"):
    for fr in ("FACE ", "NODE "):
        n, md, mn, sd = _stat(face_acc.get((capst, fr), []))
        if n:
            print(f"      {fr:<7} {capst:<7} | {n:>6} {md:>9.3f} {mn:>9.3f} {sd:>8.3f}")
print()
print("=" * 92)
print("VERDICT (read the RELAY vs DIRECT rows): if RELAY residual tracks DIRECT residual and both centre")
print("near 0 on the BLENDED estimators, the enrichment ratio is valid across tiny-region relays.")
