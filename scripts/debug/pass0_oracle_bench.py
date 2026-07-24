"""PASS-0 vs ORACLE — the full 32-scenario ambig_dense_10mb battery, per scenario, as a paired A/B.

This is the pass-0 benchmark: for every node with mass, compare the pass-0 solved f_g against the oracle f_g
(from sim_oracle.bam), mass-weighted. Reports EVERY scenario -- stranded (ss_0.99) as well as unstranded
(ss_0.50), all capture states -- because a message change can help one arm and regress the other, and an
aggregate hides exactly that.

Metrics per scenario:
  mwae  mass-weighted mean |f_g - oracle|   <- the headline (lower is better)
  corr  mass-weighted corr(f_g, oracle)     <- does the solve TRACK the truth (higher is better)
  and both split region / boundary.

Usage:
    OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --arm base   # RIGEL_P0_EXON_RNA=0
    OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --arm p0
    OMP_NUM_THREADS=1 python scripts/debug/pass0_oracle_bench.py --report
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import importlib
import os
import sys
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

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
OUT = Path(os.environ.get("P0_BENCH_OUT", "/tmp/pass0_oracle_bench.tsv"))

ap = argparse.ArgumentParser()
ap.add_argument("--arm", help="label for this run's rows, e.g. base / p0 / dircb")
ap.add_argument("--vs", default="base", help="baseline arm label for --report")
ap.add_argument("--new", default="p0", help="treatment arm label for --report")
ap.add_argument("--report", action="store_true")
args = ap.parse_args()


def wstat(fg, fo, w):
    """mass-weighted mean |err| and mass-weighted correlation."""
    if fg.size < 3 or w.sum() <= 0:
        return float("nan"), float("nan")
    mwae = float(np.average(np.abs(fg - fo), weights=w))
    ma, mb = np.average(fg, weights=w), np.average(fo, weights=w)
    va = np.average((fg - ma) ** 2, weights=w)
    vb = np.average((fo - mb) ** 2, weights=w)
    if va < 1e-12 or vb < 1e-12:
        return mwae, float("nan")
    return mwae, float(np.average((fg - ma) * (fo - mb), weights=w) / np.sqrt(va * vb))


if args.arm:
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    rows = []
    for cond in conds:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
        dbg: dict = {}
        _fac = os.environ.get('P0_FACTORY')
        cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0,  # PASS-0 ONLY
                                 **({} if _fac is None else {'intron_factory': _fac != '0'}))
        calmod.calibrate(
            inp["payload"], ra, inp["strand_model"],
            np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
        )
        chain, cap = dbg["chain"], dbg["capture"]
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        G, R = Gp + Gn, Rp + Rn
        fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
        fg = np.asarray(cap["f_g"])
        mass = np.asarray(cap["mass_global"])
        is_reg = np.asarray(chain.kind) == REGION
        ok = np.isfinite(fo) & (mass > _EPS)
        rec = {"arm": args.arm, "cond": cond}
        for tag, m in (("all", ok), ("reg", ok & is_reg), ("bnd", ok & ~is_reg)):
            e, c = wstat(fg[m], fo[m], mass[m])
            rec[f"mwae_{tag}"], rec[f"corr_{tag}"] = e, c
        rec["mass"] = float(mass[ok].sum())
        rows.append(rec)
        print(f"  {args.arm:>4} {cond:<48} mwae={rec['mwae_all']:.4f} corr={rec['corr_all']:.3f}", flush=True)
    hdr = list(rows[0].keys())
    exists = OUT.exists()
    with OUT.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
        if not exists:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows -> {OUT}")

if args.report:
    d: dict = {}
    with OUT.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[(r["arm"], r["cond"])] = r
    conds = sorted({k[1] for k in d})

    def show(title, sel):
        print(f"\n{title}")
        hdr = (
            f"{'scenario':<44} | {'mwae base':>9} {'mwae P0':>8} {'delta':>8} |"
            f" {'corr base':>9} {'corr P0':>8} {'delta':>7} |"
        )
        print(hdr)
        print("-" * len(hdr))
        agg = []
        for c in conds:
            if not sel(c):
                continue
            b, p = d.get((args.vs, c)), d.get((args.new, c))
            if not b or not p:
                continue
            eb, ep = float(b["mwae_all"]), float(p["mwae_all"])
            cb, cp = float(b["corr_all"]), float(p["corr_all"])
            mk = "  BETTER" if ep < eb - 0.002 else ("  worse" if ep > eb + 0.002 else "")
            print(
                f"{c:<44} | {eb:>9.4f} {ep:>8.4f} {ep-eb:>+8.4f} |"
                f" {cb:>9.3f} {cp:>8.3f} {cp-cb:>+7.3f} |{mk}"
            )
            agg.append((eb, ep, cb, cp, float(b["mass"])))
        if agg:
            a = np.array([x[:4] for x in agg])
            w = np.array([x[4] for x in agg])
            print(
                f"{'  MASS-WEIGHTED MEAN':<44} | {np.average(a[:,0],weights=w):>9.4f}"
                f" {np.average(a[:,1],weights=w):>8.4f} {np.average(a[:,1]-a[:,0],weights=w):>+8.4f} |"
                f" {np.nanmean(a[:,2]):>9.3f} {np.nanmean(a[:,3]):>8.3f}"
                f" {np.nanmean(a[:,3]-a[:,2]):>+7.3f} |"
            )
            print(f"{'  scenarios better / worse':<44} | "
                  f"{int((a[:,1] < a[:,0]-0.002).sum())} better, "
                  f"{int((a[:,1] > a[:,0]+0.002).sum())} worse, "
                  f"{int((np.abs(a[:,1]-a[:,0]) <= 0.002).sum())} flat")

    show("=== CAPTURE OFF ===", lambda c: "capture_off" in c)
    show("=== CAPTURE ON ===", lambda c: "capture_on" in c)
    show("=== CAPTURE VERYSTRONG ===", lambda c: "verystrong" in c)
    show("=== STRANDED (ss_0.99) — all capture states ===", lambda c: "ss_0.99" in c)
    show("=== UNSTRANDED (ss_0.50) — all capture states ===", lambda c: "ss_0.50" in c)
    show("=== ALL 32 ===", lambda c: True)
