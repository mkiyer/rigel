"""WHERE is the boundary's gDNA signal lost? Decompose the information chain, per unstranded boundary POSITION
across the 20 ss0.50 scenarios (hold node fixed, vary gDNA level), all correlations with the boundary's ORACLE
f_g unless noted:

  CEILING     corr(bnd_oracle, flank_oracle)   — is the answer even IN the neighbour? (arithmetic-independent)
  SRC-QUALITY corr(flank_oracle, flank_solved) — does the neighbour region SOLVE its own f_g right?
  TRANSFER    corr(flank_solved, msg_implied)  — does the message faithfully carry the neighbour's solve?
  DELIVERED   corr(bnd_oracle, msg_implied)    — what the boundary actually receives (measured: ~ −0.07)
  SOLVED      corr(bnd_oracle, bnd_solved)     — the boundary's final solve (~0.03)

`flank` = the better-correlated of the two flanks per position. If CEILING is high but DELIVERED is low, the info
is there and the arithmetic loses it (owner's thesis). The stage where corr collapses names the bug."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
_SS = sys.argv[1] if len(sys.argv) > 1 else "0.50"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_") and _SS in d.name)

# per (nid, scenario) columns for boundaries; flank values gathered from the region-node oracle/solve
cols = {k: [] for k in ("nid", "bo", "bs", "msg", "msgL", "msgR", "Lo", "Ls", "Ro", "Rs", "ltyp", "rtyp")}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    node_type, _ = _node_region_type(chain, ra)
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    fg = np.asarray(cap["f_g"]); kind = np.asarray(chain.kind); mass = np.asarray(cap["mass_global"])
    mode_g, prec_g = np.asarray(cap["mode_g"]), np.asarray(cap["prec_g"])
    msg = np.where(prec_g > _EPS, np.clip(np.exp(np.clip(mode_g, -30, 0.0)), 0, 1), np.nan)
    a_fwd, b_bwd = cap["a_fwd"], cap["b_bwd"]  # fwd=left flank msg, bwd=right flank msg
    def _im(mo, pr):
        return np.where(np.asarray(pr) > _EPS, np.clip(np.exp(np.clip(np.asarray(mo), -30, 0.0)), 0, 1), np.nan)
    msgL = _im(a_fwd[0], a_fwd[1]); msgR = _im(b_bwd[0], b_bwd[1])
    left, right = np.asarray(chain.left), np.asarray(chain.right)
    bnd = (kind != REGION) & np.isfinite(fo) & (mass > _EPS)
    for i in np.where(bnd)[0]:
        L, Rr = int(left[i]), int(right[i])
        cols["nid"].append(int(i)); cols["bo"].append(fo[i]); cols["bs"].append(fg[i]); cols["msg"].append(msg[i])
        cols["msgL"].append(msgL[i]); cols["msgR"].append(msgR[i])
        cols["Lo"].append(fo[L] if L >= 0 else np.nan); cols["Ls"].append(fg[L] if L >= 0 else np.nan)
        cols["Ro"].append(fo[Rr] if Rr >= 0 else np.nan); cols["Rs"].append(fg[Rr] if Rr >= 0 else np.nan)
        cols["ltyp"].append(int(node_type[L]) if L >= 0 else -1)
        cols["rtyp"].append(int(node_type[Rr]) if Rr >= 0 else -1)
C = {k: np.array(v) for k, v in cols.items()}


def _c(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 6 or x[ok].std() < 1e-6 or y[ok].std() < 1e-6:
        return np.nan
    return np.corrcoef(x[ok], y[ok])[0, 1]


# per position: pick the flank (L or R) whose oracle best correlates with the boundary oracle → the "ceiling" flank
uniq = np.unique(C["nid"])
_fk = ("ceiling", "srcq", "src2bnd", "transfer1", "transferC", "delivered", "solved")
res = {k: [] for k in _fk + ("ltyp", "rtyp")}
for u in uniq:
    m = C["nid"] == u
    bo, bs, msgc = C["bo"][m], C["bs"][m], C["msg"][m]
    cL, cR = _c(bo, C["Lo"][m]), _c(bo, C["Ro"][m])
    # choose the flank with the larger |ceiling| — its own oracle/solve/message
    if (np.isnan(cL) and np.isnan(cR)):
        continue
    useL = (np.nan_to_num(abs(cL), nan=-1) >= np.nan_to_num(abs(cR), nan=-1))
    fo_f = C["Lo"][m] if useL else C["Ro"][m]
    fs_f = C["Ls"][m] if useL else C["Rs"][m]
    msg_f = C["msgL"][m] if useL else C["msgR"][m]  # the CHOSEN flank's OWN message (isolates the mode)
    res["ceiling"].append(cL if useL else cR)
    res["srcq"].append(_c(fo_f, fs_f))                 # flank oracle → flank solve
    res["src2bnd"].append(_c(bo, fs_f))                # flank SOLVE → boundary oracle (the transferable signal)
    res["transfer1"].append(_c(fs_f, msg_f))           # flank solve → its OWN message (the pure mode)
    res["transferC"].append(_c(fs_f, msgc))            # flank solve → the COMBINED message
    res["delivered"].append(_c(bo, msgc))              # boundary oracle → combined message
    res["solved"].append(_c(bo, bs))
    res["ltyp"].append(int(np.round(np.median(C["ltyp"][m]))))
    res["rtyp"].append(int(np.round(np.median(C["rtyp"][m]))))
Rd = {k: np.array(v, float) if k in _fk else np.array(v) for k, v in res.items()}


def mean(a):
    a = a[np.isfinite(a)]
    return (a.mean(), len(a)) if len(a) else (float("nan"), 0)


print(f"[ss_{_SS}] {len(conds)} scenarios, {len(uniq)} boundary positions\n")
print("INFORMATION CHAIN (per-position mean corr, ceiling-flank chosen):")
for lab, key in (("CEILING    corr(bnd_or,   flank_or)", "ceiling"),
                 ("SRC-QUAL   corr(flank_or, flank_solved)", "srcq"),
                 ("SRC→BND    corr(bnd_or,   flank_solved)", "src2bnd"),
                 ("TRANSFER1  corr(flank_solved, flank_msg)", "transfer1"),
                 ("TRANSFER-C corr(flank_solved, combined_msg)", "transferC"),
                 ("DELIVERED  corr(bnd_or,   combined_msg)", "delivered"),
                 ("SOLVED     corr(bnd_or,   bnd_solved)", "solved")):
    mv, n = mean(Rd[key])
    print(f"  {lab:>44} : {mv:>7.3f}  (n={n})")

# TRACE: the highest-ceiling, worst-transfer positions — show the source-solved f_g vs the message it produces
print("\nTRACE — 4 high-ceiling positions: does the message track the source's solved f_g across scenarios?")
score = []
for u in uniq:
    m = C["nid"] == u
    cL, cR = _c(C["bo"][m], C["Lo"][m]), _c(C["bo"][m], C["Ro"][m])
    if np.isnan(cL) and np.isnan(cR):
        continue
    useL = np.nan_to_num(abs(cL), nan=-1) >= np.nan_to_num(abs(cR), nan=-1)
    ceil = cL if useL else cR
    t1 = _c((C["Ls"] if useL else C["Rs"])[m], (C["msgL"] if useL else C["msgR"])[m])
    if np.isfinite(ceil) and np.isfinite(t1):
        score.append((ceil - abs(t1), u, useL, ceil, t1))
score.sort(reverse=True)
for _, u, useL, ceil, t1 in score[:4]:
    m = C["nid"] == u
    order = np.argsort(C["bo"][m])  # sort by boundary oracle to see monotonicity
    bo = C["bo"][m][order]
    fs = (C["Ls"] if useL else C["Rs"])[m][order]
    mg = (C["msgL"] if useL else C["msgR"])[m][order]
    fo_f = (C["Lo"] if useL else C["Ro"])[m][order]
    print(f"  node {u} (flank={'L' if useL else 'R'}, ceiling={ceil:.2f}, transfer1={t1:.2f}):")
    print(f"    bnd_oracle : {np.array2string(bo, precision=2, floatmode='fixed', max_line_width=200)}")
    print(f"    flank_orac : {np.array2string(fo_f, precision=2, floatmode='fixed', max_line_width=200)}")
    print(f"    flank_SOLV : {np.array2string(fs, precision=2, floatmode='fixed', max_line_width=200)}")
    print(f"    flank_MSG  : {np.array2string(mg, precision=2, floatmode='fixed', max_line_width=200)}")

print("\nBy ceiling-flank type (is the ceiling higher when the informative flank is an exon vs intron?):")
_TN = {-1: "term", 0: "intgc", 1: "intron", 2: "exon"}
# the chosen flank type = whichever side we used; recompute per position isn't stored, so bucket by max ceiling side
for t in (2, 1, 0):
    # a position "has an exon flank" if either flank is type t
    m = (Rd["ltyp"] == t) | (Rd["rtyp"] == t)
    ce, n = mean(Rd["ceiling"][m]); de, _ = mean(Rd["delivered"][m]); so, _ = mean(Rd["solved"][m])
    print(f"  flank includes {_TN[t]:>7}: ceiling={ce:>6.3f}  delivered={de:>6.3f}  solved={so:>6.3f}  (n_pos={n})")
