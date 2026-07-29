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
from rigel.index import TranscriptIndex, load_manifest  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
RUNS = Path("/Users/mkiyer/Downloads/rigel_runs")
OUT = Path(os.environ.get("P0_BENCH_OUT", "/tmp/pass0_oracle_bench.tsv"))

ap = argparse.ArgumentParser()
ap.add_argument("--arm", help="label for this run's rows, e.g. base / p0 / dircb")
ap.add_argument("--vs", default="base", help="baseline arm label for --report")
ap.add_argument("--new", default="p0", help="treatment arm label for --report")
ap.add_argument("--report", action="store_true")
ap.add_argument(
    "--suite",
    default=os.environ.get("P0_SUITE", "ambig_dense_10mb"),
    help="suite directory name under ~/Downloads/rigel_runs (default ambig_dense_10mb, THE bench). "
    "⚠ ambig_dense_10mb's annotation has ZERO mergeable adjacencies, so its v8 partition is "
    "byte-identical to v7 and it CANNOT show the partition effect; quick_3to1_5mb (+25.0 %) can.",
)
ap.add_argument(
    "--coarsen",
    action="store_true",
    help="permit --report to diff arms recorded on DIFFERENT partitions (v7 vs v8). Refused by "
    "default: the v8 partition refines v7, so a raw mwae diff across it confounds the partition "
    "change with whatever the arm varied. See docs/CARRY_FORWARD.md F8.",
)
args = ap.parse_args()
SUITE = RUNS / args.suite

#: Provenance columns. `partition` is the index format ("v7" = merged signature regions,
#: "v8" = the splice graph); `mass_kind` names what the mwae weight actually IS, because under
#: the accumulator v5 rework it changes from a fractional mass to an integer count; `suite` names
#: the condition set, because the partition effect is suite-dependent (measured 2026-07-29: +0.0 %
#: on ambig_dense_10mb, +25.0 % on quick_3to1_5mb, +38.7 % on the human annotation) and two suites'
#: rows are not comparable. Rows written before a column existed are read back as its default,
#: which is what they were.
_DEFAULT_PROVENANCE = {
    "partition": "v7",
    "mass_kind": "mass_frac",
    "suite": "ambig_dense_10mb",
}


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
    manifest = load_manifest(index.index_dir) or {}
    partition = f"v{int(manifest.get('format_version', 7))}"
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
    rows = []
    for cond in conds:
        inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
        dbg: dict = {}
        _fac = os.environ.get('P0_FACTORY')
        _str = os.environ.get('P0_PRIOR_STRENGTH')  # gDNA-hyperprior temperature (W4 sweep)
        cc = dataclasses.replace(cfg.calibration, calib_refit_iters=int(os.environ.get('P0_REFIT', '0')),  # 0=PASS-0
                                 **({} if _fac is None else {'intron_factory': _fac != '0'}),
                                 **({} if _str is None else {'gdna_prior_strength': float(_str)}))
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
        rec["partition"] = partition
        rec["mass_kind"] = "mass_frac"  # v5 W5a turns this into an integer count
        rec["suite"] = args.suite
        rec["refit"] = int(cc.calib_refit_iters)
        rows.append(rec)
        print(f"  {args.arm:>4} {cond:<48} mwae={rec['mwae_all']:.4f} corr={rec['corr_all']:.3f}", flush=True)
    hdr = list(rows[0].keys())
    exists = OUT.exists()
    if exists:
        # The header is written once, so appending rows with a DIFFERENT column set silently
        # misaligns every later read. Refuse instead.
        with OUT.open() as fh:
            on_disk = (fh.readline().rstrip("\n").split("\t")) if fh else []
        if on_disk and on_disk != hdr:
            raise SystemExit(
                f"{OUT} has columns {on_disk} but this run writes {hdr}. Point P0_BENCH_OUT at a "
                f"fresh file rather than appending a different schema."
            )
    with OUT.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr, delimiter="\t")
        if not exists:
            w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows ({partition}, refit={cc.calib_refit_iters}) -> {OUT}")

if args.report:
    d: dict = {}
    with OUT.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for k, v in _DEFAULT_PROVENANCE.items():
                r[k] = r.get(k) or v  # rows predating the column were recorded on v7 / fractional mass
            d[(r["arm"], r["cond"])] = r
    conds = sorted({k[1] for k in d})

    def _provenance(arm):
        return {
            (r["partition"], r["mass_kind"], r["suite"]) for (a, _), r in d.items() if a == arm
        }

    prov_b, prov_p = _provenance(args.vs), _provenance(args.new)
    for arm, prov in ((args.vs, prov_b), (args.new, prov_p)):
        if len(prov) > 1:
            raise SystemExit(f"arm {arm!r} mixes provenances {sorted(prov)}; re-record it.")
    # ⚠ A suite mismatch is NEVER waivable — --coarsen aggregates a refined partition back onto its
    # own coarse parent, which is meaningless across two different condition sets and genomes.
    suites_b = {p[2] for p in prov_b}
    suites_p = {p[2] for p in prov_p}
    if suites_b and suites_p and suites_b != suites_p:
        raise SystemExit(
            f"REFUSING to diff {args.vs!r} on {sorted(suites_b)} against {args.new!r} on "
            f"{sorted(suites_p)}: different suites are different genomes and different condition "
            f"sets. There is no aggregation that makes them comparable."
        )
    if prov_b and prov_p and prov_b != prov_p and not args.coarsen:
        raise SystemExit(
            f"REFUSING to diff {args.vs!r} {sorted(prov_b)} against {args.new!r} {sorted(prov_p)}.\n"
            "The v8 partition REFINES v7, so a raw per-object mwae diff across it confounds the "
            "partition change with whatever this arm varied — and effective lengths are "
            "superadditive, so the two populations are not comparable object-for-object "
            "(measured: sum E(children)/E(whole) = 0.765). Record a baseline on the SAME partition, "
            "or pass --coarsen if you have genuinely aggregated to a common object set."
        )

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
