"""TARGET step 4 — is the gDNA channel's mode error in the SOURCE density or in the TRANSPORT ratio?

For the gDNA component the transport is exact IN PRINCIPLE: gDNA is genomically uniform in CONTENT, so the
ratio of true gDNA densities between two nodes IS the capture step, and

    rho_g(src) * r_true * E_g(dst) / M(dst)  ==  f_g^true(dst)      exactly.

So every bit of the gDNA channel's mode error is either (a) the source's own gDNA density being wrong, or
(b) `r` failing to be the true gDNA capture step (it is a ratio of TOTAL densities, which is only a proxy and
only valid when the two nodes share a composition -- and an intron vs its exon emphatically do not).

This decomposes the delivered gDNA claim on every boundary->exon and intron->boundary edge into those two,
by substituting the ORACLE capture step and the ORACLE source density one at a time.

    OMP_NUM_THREADS=1 python scratchpad/t4_gdna_transport.py
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
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.50_nrna_none_capture_on"

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
dbg: dict = {}
calmod.calibrate(
    inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
    np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
)
chain, cap = dbg["chain"], dbg["capture"]
uni, us = cap["_uni"][-1], cap["_uni_static"]
Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
G, R = Gp + Gn, Rp + Rn
T = G + R
fo = np.where(T > _EPS, G / np.maximum(T, _EPS), np.nan)
E_g, E_r, M = us["E_g"], us["E_r"], us["M"]
li, ri = us["left"], us["right"]
n = len(M)
mass = np.asarray(cap["mass_global"])
rt, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind)
CLSN = {0: "intergenic", 1: "intron", 2: "exon"}
cls = np.array([CLSN.get(int(rt[i]), "?") if kind[i] == REGION else "boundary" for i in range(n)])
rho_g_true = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)  # gDNA is uniform in CONTENT ⇒ this IS capture
rho_l0, rho_r0, rho_node = us["rho_l0"], us["rho_r0"], us["rho_node0"]


def lg(x):
    return np.log(np.maximum(np.asarray(x, np.float64), 1e-12))


rows = []
for src_i, relay, face in ((li, "fwd_g", rho_r0), (ri, "bwd_g", rho_l0)):
    s = np.clip(src_i, 0, n - 1)
    rg = us[relay][s]  # the relay's gDNA density AT THE SOURCE
    r_mod = np.where(face[s] > _EPS, rho_node / np.maximum(face[s], _EPS), 1.0)
    r_true = rho_g_true / np.maximum(rho_g_true[s], _EPS)  # the TRUE capture step
    ok = ((src_i >= 0) & np.isfinite(fo) & (fo > 1e-6) & (mass > _EPS) & (rg > _EPS)
          & np.isfinite(r_true) & (r_true > 0) & (face[s] > _EPS) & (rho_node > _EPS)
          & np.isfinite(rho_g_true[s]) & (rho_g_true[s] > 0))
    rows.append((ok, s, rg, r_mod, r_true))

print(f"{COND}\n{'=' * 118}")
print("Delivered gDNA claim  f_g^claim = rho_g(src)*r*E_g/M , decomposed. All logs, mass-weighted.")
print(f"{'edge (src->dst)':<28}{'n':>6}{'|err| model':>12}{'|err| ORACLE r':>16}"
      f"{'|err| oracle rho_g':>20}{'|log r/r_true|':>16}")
for lab, sel in (("intron -> boundary", lambda s, d: (cls[s] == "intron") & (cls[d] == "boundary")),
                 ("boundary -> exon", lambda s, d: (cls[s] == "boundary") & (cls[d] == "exon")),
                 ("exon -> boundary", lambda s, d: (cls[s] == "exon") & (cls[d] == "boundary")),
                 ("boundary -> intron", lambda s, d: (cls[s] == "boundary") & (cls[d] == "intron")),
                 ("intergenic -> boundary", lambda s, d: (cls[s] == "intergenic") & (cls[d] == "boundary"))):
    acc = {k: [] for k in ("a", "b", "c", "d", "w")}
    for ok, s, rg, r_mod, r_true in rows:
        dst = np.arange(n)
        m = ok & sel(s, dst)
        if not m.any():
            continue
        claim_mod = lg(rg[m] * r_mod[m] * E_g[m] / np.maximum(M[m], _EPS))
        claim_orc = lg(rg[m] * r_true[m] * E_g[m] / np.maximum(M[m], _EPS))
        claim_src = lg(rho_g_true[s][m] * r_mod[m] * E_g[m] / np.maximum(M[m], _EPS))
        tgt = lg(fo[m])
        acc["a"].append(np.abs(claim_mod - tgt)); acc["b"].append(np.abs(claim_orc - tgt))
        acc["c"].append(np.abs(claim_src - tgt)); acc["d"].append(np.abs(lg(r_mod[m]) - lg(r_true[m])))
        acc["w"].append(mass[m])
    if not acc["w"]:
        continue
    a, b, c, dd, w = (np.concatenate(acc[k]) for k in ("a", "b", "c", "d", "w"))
    print(f"{lab:<28}{w.size:>6}{np.average(a, weights=w):>12.3f}{np.average(b, weights=w):>16.3f}"
          f"{np.average(c, weights=w):>20.3f}{np.average(dd, weights=w):>16.3f}")

# ── CANDIDATE prior-free frame ratios, scored against the TRUE capture step ──────────────────────────────
print(f"\n{'edge':<26}{'n':>6}" + "".join(f"{k:>14}" for k in
      ("r_tot (ship)", "r_M/E_g", "r_M/E_r", "r_geo(tot,M/Eg)")))
Mg = M / np.maximum(E_g, _EPS)
Mr = M / np.maximum(E_r, _EPS)
for lab, sel in (("intron -> boundary", lambda s, d: (cls[s] == "intron") & (cls[d] == "boundary")),
                 ("boundary -> exon", lambda s, d: (cls[s] == "boundary") & (cls[d] == "exon")),
                 ("exon -> boundary", lambda s, d: (cls[s] == "exon") & (cls[d] == "boundary")),
                 ("boundary -> intron", lambda s, d: (cls[s] == "boundary") & (cls[d] == "intron")),
                 ("ALL edges", lambda s, d: np.ones(len(s), bool))):
    acc = {k: [] for k in ("t", "g", "r", "b", "w")}
    for ok, s, rg, r_mod, r_true in rows:
        m = ok & sel(s, np.arange(n))
        if not m.any():
            continue
        rt_ = lg(r_true[m])
        acc["t"].append(np.abs(lg(r_mod[m]) - rt_))
        acc["g"].append(np.abs(lg(Mg[m] / np.maximum(Mg[s][m], _EPS)) - rt_))
        acc["r"].append(np.abs(lg(Mr[m] / np.maximum(Mr[s][m], _EPS)) - rt_))
        acc["b"].append(np.abs(0.5 * (lg(r_mod[m]) + lg(Mg[m] / np.maximum(Mg[s][m], _EPS))) - rt_))
        acc["w"].append(mass[m])
    if not acc["w"]:
        continue
    t, g, r_, b, w = (np.concatenate(acc[k]) for k in ("t", "g", "r", "b", "w"))
    print(f"{lab:<26}{w.size:>6}" + "".join(f"{np.average(x, weights=w):>14.3f}" for x in (t, g, r_, b)))
print("  (mass-weighted |log r_candidate / r_true| — lower is a better capture-step estimator)")

print("\n  'model' = as shipped;  'ORACLE r' = model source density, TRUE capture step;")
print("  'oracle rho_g' = TRUE source density, model r.  Lower is better; these isolate the two defects.")

# the true capture step distribution per edge kind, vs the model's r
print(f"\n{'edge':<28}{'median r_model':>16}{'median r_true':>16}{'ratio':>10}")
for lab, sel in (("intron -> boundary", lambda s, d: (cls[s] == "intron") & (cls[d] == "boundary")),
                 ("boundary -> exon", lambda s, d: (cls[s] == "boundary") & (cls[d] == "exon")),
                 ("exon -> boundary", lambda s, d: (cls[s] == "exon") & (cls[d] == "boundary")),
                 ("boundary -> intron", lambda s, d: (cls[s] == "boundary") & (cls[d] == "intron"))):
    rm, rtv = [], []
    for ok, s, rg, r_mod, r_true in rows:
        m = ok & sel(s, np.arange(n))
        rm.append(r_mod[m]); rtv.append(r_true[m])
    rm, rtv = np.concatenate(rm), np.concatenate(rtv)
    if rm.size < 5:
        continue
    print(f"{lab:<28}{np.median(rm):>16.3f}{np.median(rtv):>16.3f}"
          f"{np.median(rm) / max(np.median(rtv), _EPS):>10.3f}")
