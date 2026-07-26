"""RC2 (A) — how far does the ANCHOR measurement precision (`mg_own = struct_lock·prec_g`) travel, and how
much does it decay?  The measurement stream carries only PRECISION; the mode `mo_g` comes from a different
(density-fusion) stream — so this asks whether a far-travelled anchor precision is still large."""
from __future__ import annotations
import dataclasses, importlib, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node
from selfsolve_diag import _scan_and_truth
from rigel.calibration.bp_solver import REGION
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
CLS = {0: "intergenic", 1: "intron", 2: "exon", -1: "boundary"}
_EPS = 1e-9
index = TranscriptIndex.load(str(SUITE / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
for cond in ("gdna_gdna300_ss_0.99_nrna_none_capture_on", "gdna_gdna300_ss_0.99_nrna_present_capture_on"):
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg = {}
    calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                     np.asarray(inp["rna_fl_pmf"]),
                     dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    uni, us = cap["_uni"][-1], cap["_uni_static"]
    n = np.asarray(cap["f_g"]).shape[0]
    left, right = us["left"], us["right"]
    seed = us["mg_own"] > 0.0            # only struct_lock (intergenic REGION) nodes seed the stream
    # BFS chain distance to the nearest seed
    dist = np.full(n, 1 << 20, np.int64); dist[seed] = 0
    for _ in range(200):
        ch = False
        for a_, b_ in ((left, right), (right, left)):
            v = a_ >= 0
            nd = np.where(v, dist[np.clip(a_, 0, n - 1)] + 1, 1 << 20)
            upd = nd < dist
            if upd.any(): dist = np.where(upd, nd, dist); ch = True
            _ = b_
        if not ch: break
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    rt, _ = _node_region_type(chain, ra); kind = np.asarray(chain.kind)
    cls = np.array([CLS[int(rt[i])] if kind[i] == REGION else "boundary" for i in range(n)])
    mass = np.asarray(cap["mass_global"]); cm_g = np.asarray(uni["cm_g"])
    mo_g = np.exp(np.asarray(uni["mo_g"]))
    m = np.isfinite(fo) & (mass > 0) & (us["tau_own"] > 0) & (cls == "exon") & (cm_g > 0)
    print(f"\n### {cond}   seeds(struct_lock)={int(seed.sum())}  full-rank exons w/ live gDNA msg={int(m.sum())}")
    print(f"{'hops to seed':<14}{'n':>6}{'median cm_g':>13}{'p90 cm_g':>11}{'mean mo_g':>11}{'mean oracle':>13}{'bias':>8}")
    for lo, hi in ((1, 4), (4, 8), (8, 16), (16, 32), (32, 1 << 20)):
        s = m & (dist >= lo) & (dist < hi)
        if not s.any(): continue
        print(f"[{lo},{hi if hi < 1000 else 'inf'})".ljust(14) + f"{int(s.sum()):>6}{np.median(cm_g[s]):>13.1f}"
              f"{np.quantile(cm_g[s], .9):>11.1f}{np.average(mo_g[s], weights=mass[s]):>11.3f}"
              f"{np.average(fo[s], weights=mass[s]):>13.3f}"
              f"{np.average(mo_g[s], weights=mass[s]) - np.average(fo[s], weights=mass[s]):>+8.3f}")
