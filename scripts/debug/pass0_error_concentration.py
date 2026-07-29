"""PASS-0 vs ORACLE — WHERE IS THE ERROR? The full ambig_dense_10mb suite, by node class and by axis.

`pass0_oracle_bench.py` answers *how much* error each scenario carries; this answers **where it is
concentrated**. Same substrate (cached scan + oracle, `calib_refit_iters=0` — PASS-0 ONLY, no refit, no
downstream EM), but it reports the **error MASS** ``Σ mass·|f_g − oracle|`` rather than only the
mass-weighted mean, because a class can have a large per-node error while holding almost none of the
library's mass (and so almost none of the damage), or the reverse.

Four node classes: intergenic / intron / exon / boundary. Three rollups: per scenario, per class, and per
scenario AXIS (gDNA level × strandedness × capture × nRNA), plus the worst (scenario × class) cells.

    OMP_NUM_THREADS=1 python scripts/debug/pass0_error_concentration.py --arm base
    RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/pass0_error_concentration.py --arm uni
    OMP_NUM_THREADS=1 python scripts/debug/pass0_error_concentration.py --report --arm base
    OMP_NUM_THREADS=1 python scripts/debug/pass0_error_concentration.py --report --arm uni --vs base
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

from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
OUT = Path(os.environ.get("P0_CONC_OUT", "/tmp/pass0_error_concentration.tsv"))
CLASSES = ("intergenic", "intron", "exon", "boundary")

ap = argparse.ArgumentParser()
ap.add_argument("--arm")
ap.add_argument("--report", action="store_true")
ap.add_argument("--vs", default=None, help="baseline arm to diff against in --report")
args = ap.parse_args()


def _axes(cond: str):
    """(gdna level, strandedness, nRNA, capture) parsed from the scenario name."""
    p = cond.split("_")
    g = p[1] if len(p) > 1 else "?"
    ss = next((p[i + 1] for i, x in enumerate(p) if x == "ss"), "?")
    nr = next((p[i + 1] for i, x in enumerate(p) if x == "nrna"), "?")
    cap = next((p[i + 1] for i, x in enumerate(p) if x == "capture"), "?")
    return g, ss, nr, cap


# ─────────────────────────────────────── measure ───────────────────────────────────────
if args.arm and not args.report:
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    rows = []
    for cond in conds:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
        dbg: dict = {}
        _fac = os.environ.get('P0_FACTORY')
        cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0,  # PASS-0 ONLY
                                 **({} if _fac is None else {'intron_factory': _fac != '0'}))
        calmod.calibrate(inp["payload"], ra, inp["strand_model"],
                         np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
        chain, cap = dbg["chain"], dbg["capture"]
        Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
        G, R = Gp + Gn, Rp + Rn
        fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
        fg = np.asarray(cap["f_g"])
        mass = np.asarray(cap["mass_global"])
        rt, _ = _node_region_type(chain, ra)
        cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
        ok = np.isfinite(fo) & (mass > _EPS)
        err = np.abs(fg - fo)
        for ci, cname in enumerate(CLASSES):
            m = ok & (cls == ci)
            if not m.sum():
                continue
            w = mass[m]
            rows.append({
                "arm": args.arm, "cond": cond, "cls": cname, "n": int(m.sum()),
                "mass": float(w.sum()),
                "errmass": float((w * err[m]).sum()),        # <- the damage
                "mwae": float(np.average(err[m], weights=w)),
                "mean_fg": float(np.average(fg[m], weights=w)),
                "mean_oracle": float(np.average(fo[m], weights=w)),
            })
        tot = float((mass[ok] * err[ok]).sum())
        print(f"  {args.arm:>6} {cond:<48} errmass={tot:12.1f} "
              f"mwae={np.average(err[ok], weights=mass[ok]):.4f}", flush=True)
    hdr = list(rows[0].keys())
    exists = OUT.exists()
    with OUT.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
        if not exists:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows -> {OUT}")

# ─────────────────────────────────────── report ───────────────────────────────────────
if args.report:
    recs = [r for r in csv.DictReader(OUT.open(), delimiter="\t")]
    arm = args.arm or "base"
    A = [r for r in recs if r["arm"] == arm]
    B = [r for r in recs if r["arm"] == args.vs] if args.vs else []
    if not A:
        sys.exit(f"no rows for arm {arm!r} in {OUT}")
    conds = sorted({r["cond"] for r in A})

    def get(rows, cond, cls):
        for r in rows:
            if r["cond"] == cond and r["cls"] == cls:
                return r
        return None

    def F(r, k):
        return float(r[k]) if r else float("nan")

    TOT = sum(F(r, "errmass") for r in A)
    print(f"\n{'='*118}\nPASS-0 vs ORACLE — ERROR CONCENTRATION   arm={arm}"
          f"{'   vs ' + args.vs if args.vs else ''}   suite=ambig_dense_10mb (32 scenarios)")
    print(f"total error mass Σ mass·|f_g−oracle| = {TOT:,.0f}\n{'='*118}")

    # ---- 1. per scenario × class: error-mass SHARE of the whole suite ----
    print("\n### 1. ERROR-MASS SHARE (% of the suite's total error), per scenario × node class")
    print(f"{'scenario':<46}{'interg':>9}{'intron':>9}{'exon':>9}{'bound':>9}{'TOTAL%':>9}{'mwae':>9}")
    print("-" * 100)
    for c in conds:
        sh = [100.0 * F(get(A, c, k), "errmass") / TOT if get(A, c, k) else 0.0 for k in CLASSES]
        em = sum(F(get(A, c, k), "errmass") for k in CLASSES if get(A, c, k))
        mm = sum(F(get(A, c, k), "mass") for k in CLASSES if get(A, c, k))
        print(f"{c:<46}" + "".join(f"{x:>9.2f}" for x in sh)
              + f"{sum(sh):>9.2f}{em/max(mm,1e-9):>9.4f}")

    # ---- 2. per class rollup ----
    print("\n### 2. BY NODE CLASS (whole suite)")
    print(f"{'class':<12}{'n_nodes':>10}{'mass':>14}{'mass%':>8}{'errmass':>14}{'err%':>8}"
          f"{'mwae':>9}{'mean_fg':>9}{'oracle':>9}" + (f"{'Δmwae':>9}" if B else ""))
    print("-" * (94 + (9 if B else 0)))
    MTOT = sum(F(r, "mass") for r in A)
    for k in CLASSES:
        rs = [r for r in A if r["cls"] == k]
        if not rs:
            continue
        n = sum(int(r["n"]) for r in rs)
        ms = sum(F(r, "mass") for r in rs)
        em = sum(F(r, "errmass") for r in rs)
        fgm = sum(F(r, "mean_fg") * F(r, "mass") for r in rs) / max(ms, 1e-9)
        orm = sum(F(r, "mean_oracle") * F(r, "mass") for r in rs) / max(ms, 1e-9)
        line = (f"{k:<12}{n:>10,}{ms:>14,.0f}{100*ms/MTOT:>8.1f}{em:>14,.0f}"
                f"{100*em/TOT:>8.1f}{em/max(ms,1e-9):>9.4f}{fgm:>9.3f}{orm:>9.3f}")
        if B:
            rb = [r for r in B if r["cls"] == k]
            mb = sum(F(r, "mass") for r in rb)
            eb = sum(F(r, "errmass") for r in rb)
            line += f"{em/max(ms,1e-9) - eb/max(mb,1e-9):>+9.4f}"
        print(line)

    # ---- 3. per scenario AXIS ----
    print("\n### 3. BY SCENARIO AXIS (error mass share + mwae)")
    for ai, aname in enumerate(("gDNA level", "strandedness", "nRNA", "capture")):
        buckets: dict[str, list] = {}
        for c in conds:
            buckets.setdefault(_axes(c)[ai], []).append(c)
        print(f"\n  {aname}")
        print(f"    {'value':<16}{'errmass':>14}{'err%':>8}{'mass%':>8}{'mwae':>9}"
              + (f"{'Δmwae':>9}" if B else ""))
        for v, cs in sorted(buckets.items()):
            em = sum(F(get(A, c, k), "errmass") for c in cs for k in CLASSES if get(A, c, k))
            ms = sum(F(get(A, c, k), "mass") for c in cs for k in CLASSES if get(A, c, k))
            line = (f"    {v:<16}{em:>14,.0f}{100*em/TOT:>8.1f}{100*ms/MTOT:>8.1f}"
                    f"{em/max(ms,1e-9):>9.4f}")
            if B:
                eb = sum(F(get(B, c, k), "errmass") for c in cs for k in CLASSES if get(B, c, k))
                mb = sum(F(get(B, c, k), "mass") for c in cs for k in CLASSES if get(B, c, k))
                line += f"{em/max(ms,1e-9) - eb/max(mb,1e-9):>+9.4f}"
            print(line)

    # ---- 4. the worst cells ----
    print("\n### 4. WORST (scenario × class) CELLS — where the damage actually is")
    cells = sorted(A, key=lambda r: -F(r, "errmass"))[:20]
    print(f"{'scenario':<46}{'class':<11}{'errmass':>12}{'err%':>7}{'cum%':>7}"
          f"{'mwae':>9}{'mean_fg':>9}{'oracle':>9}")
    print("-" * 110)
    cum = 0.0
    for r in cells:
        cum += 100 * F(r, "errmass") / TOT
        print(f"{r['cond']:<46}{r['cls']:<11}{F(r,'errmass'):>12,.0f}"
              f"{100*F(r,'errmass')/TOT:>7.2f}{cum:>7.2f}{F(r,'mwae'):>9.4f}"
              f"{F(r,'mean_fg'):>9.3f}{F(r,'mean_oracle'):>9.3f}")
