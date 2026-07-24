"""WHY DO STRANDED CASES REGRESS WHEN RNA MESSAGES ARE TURNED ON? A precision interrogation.

On strand-specific data the per-node strand likelihood is the sharpest signal we have, so any regression means
a message is overriding it. Two questions, answered with node-level evidence:

  Q1  Is the RNA relay's CONTINUITY premise even true? gDNA is a genomic background and is continuous.
      RNA is TRANSCRIPT-specific. `free_s` only checks that a strand is PRESENT on both flanks -- not that it
      is the SAME transcript. A boundary with no spliced channel is a pure SIGNATURE CHANGE, which is exactly
      where a transcript starts/ends or another gene's feature begins, so the flanking RNA densities need not
      match at all. Measured here as the observed RNA-density step across each boundary class.

  Q2  WHERE does the error concentrate, and does it track message precision OVERPOWERING strand evidence?
      For each node we have tau0_lam (the node's own strand Fisher information -- 0 on unstranded, large on
      stranded) and the incoming message precisions. If the regression is "messages drown a pristine strand
      signal", the damaged nodes must be those where message precision is large RELATIVE to tau0.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/stranded_regression_dissect.py
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

from rigel.calibration.bp_solver import REGION, node_global_geometry  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
STRANDED = sorted(
    d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists() and "ss_0.99" in d.name
)

q1_rows: dict[str, list] = {}
q2 = []

for cond in STRANDED:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, cap, st, g = dbg["chain"], dbg["capture"], dbg["statics"], dbg["geometry"]
    kind = np.asarray(chain.kind)
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    rtype, _ = _node_region_type(chain, ra)
    mass, effg = node_global_geometry(chain, g)
    effr = np.where(
        kind == REGION,
        np.asarray(g.eff_rna_left, float),
        np.asarray(g.eff_rna_left, float) + np.asarray(g.eff_rna_right, float),
    )
    spl = (
        np.asarray(g.spliced_pos_left, float) + np.asarray(g.spliced_neg_left, float)
        + np.asarray(g.spliced_pos_right, float) + np.asarray(g.spliced_neg_right, float)
    )
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    fg = np.asarray(cap["f_g"])
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    tau0 = np.asarray(cap["_tau0_lam"])
    pg, pp, pn = np.asarray(cap["prec_g"]), np.asarray(cap["prec_p"]), np.asarray(cap["prec_n"])

    # --- Q1: the observed RNA-density step across each boundary class ---
    rho_r = np.where((effr > _EPS) & (mass > _EPS), (1.0 - np.clip(fg, 0, 1)) * mass / np.maximum(effr, _EPS), np.nan)
    for b in np.where(kind != REGION)[0]:
        lf, rg = left[b], right[b]
        if lf < 0 or rg < 0:
            continue
        if not ((fp[lf] and fp[rg]) or (fn[lf] and fn[rg])):
            continue  # no strand is continuous -> the relay would not fire anyway
        a, c = rho_r[lf], rho_r[rg]
        if not (np.isfinite(a) and np.isfinite(c) and a > _EPS and c > _EPS):
            continue
        cls = "no-spliced boundary (N4a scope)" if spl[b] <= _EPS else "splice junction (N4b scope)"
        q1_rows.setdefault(cls, []).append(abs(np.log(c) - np.log(a)))

    # --- Q2: where is the error, and does it track message precision vs strand evidence? ---
    ok = np.isfinite(fo) & (mass > _EPS)
    for i in np.where(ok)[0]:
        q2.append((abs(fg[i] - fo[i]), mass[i], tau0[i], pg[i], pp[i] + pn[i], kind[i] != REGION, rtype[i]))

print(f"\nSTRANDED (ss_0.99) scenarios: {len(STRANDED)}\n")
print("Q1 — is RNA continuous across the boundaries the relay fires on?")
print("    |log rho_R(right) - log rho_R(left)| across boundaries with a continuous strand:")
hdr = f"    {'boundary class':<34} | {'n':>6} | {'median':>7} {'p75':>7} {'p90':>7}"
print(hdr)
print("    " + "-" * (len(hdr) - 4))
for k in sorted(q1_rows):
    v = np.array(q1_rows[k])
    print(f"    {k:<34} | {v.size:>6} | {np.median(v):>7.3f} {np.percentile(v,75):>7.3f} {np.percentile(v,90):>7.3f}")

a = np.array(q2, dtype=float)
err, m, t0, pgv, prv, isb, rt = a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4], a[:, 5] > 0.5, a[:, 6]
print("\nQ2 — does the error concentrate where MESSAGES outweigh the node's own STRAND evidence?")
tot = float(np.average(err, weights=m))
print(f"    mass-weighted |err| over all nodes: {tot:.4f}")
strand_live = t0 > _EPS
ratio = np.where(strand_live, (pgv + prv) / np.maximum(t0, _EPS), np.inf)
for lo, hi, lab in [(0, 0.1, "msg << strand   (<0.1x)"), (0.1, 1, "msg <  strand"), (1, 10, "msg >  strand"), (10, np.inf, "msg >> strand  (>10x)")]:
    sel = strand_live & (ratio >= lo) & (ratio < hi)
    if sel.sum() < 5:
        continue
    print(f"    {lab:<26} n={sel.sum():>6}  mass%={100*m[sel].sum()/m.sum():>5.1f}  |err|={np.average(err[sel], weights=m[sel]):.4f}")
sel = ~strand_live
if sel.sum():
    print(f"    {'NO strand evidence (tau0=0)':<26} n={sel.sum():>6}  mass%={100*m[sel].sum()/m.sum():>5.1f}  |err|={np.average(err[sel], weights=m[sel]):.4f}")
