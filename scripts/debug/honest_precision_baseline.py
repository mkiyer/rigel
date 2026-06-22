"""Honest-precision Step-0 validation harness — faithful, reproducible calibration-prior baseline.

The honest-precision work (docs/calibration/honest_precision_unified_design.md) changes the imputation
message + global-prior *delivery* and the σ²_bio fit. This harness is the **measurement contract** that
must be in place BEFORE any precision change, so a change is judged against a recorded baseline rather
than a remembered (and, as Step-0 found, STALE) number.

What it measures — the **calibration-prior layer** (fast: calibrate() on cached payloads, no quant):
per study condition, the per-region deconvolved gDNA mass vs a **by-origin oracle** (the BAM split into
gdna-only / rna-only by the simulator read name, each scanned through the production accumulator). It
reports, **signed**:
  * phantom = Σ_i max(prod_i − oracle_i, 0)   (gDNA OVER-call — the zero-gDNA failure)
  * siphon  = Σ_i max(oracle_i − prod_i, 0)   (gDNA UNDER-call — the capture failure)
  * net, totals, and the breakdown by (strand_class × region_type) so AMBIG-exon (the #1 error site)
    and the single-strand-exon Jeffreys-push (Bug C) are each visible.

This is a **LOCAL dev gate**, not a pytest test: the suite (~500 MB BAMs) lives under ~/Downloads, not
in the repo. The portable invariants (factor-1-uniform, the node-level defer / two-sided unit tests for
Step 2) stay in pytest.

The post-EM `net_gdna_to_rna` gate (does a prior change propagate through quant?) is the SECONDARY,
slower layer: regenerate via the `calibration-benchmark` skill (`scripts/sim/evaluate_suite.py`) and read
`net_flow_per_condition.tsv`. This harness records the current net-flow values for reference but does not
re-quant.

Usage (always OMP_NUM_THREADS=1 for bit-reproducibility):
    OMP_NUM_THREADS=1 python scripts/debug/honest_precision_baseline.py                 # measure + diff vs baseline
    OMP_NUM_THREADS=1 python scripts/debug/honest_precision_baseline.py --update-baseline  # (re)record the baseline
    OMP_NUM_THREADS=1 python scripts/debug/honest_precision_baseline.py --rebuild        # force re-scan (after a scan/deposit change)
    ... [--only stranded-0,capture]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# reuse the by-origin oracle cache builder (full / gdna-only / rna-only payloads per condition)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
BASELINE_PATH = Path(__file__).resolve().parent / "honest_precision_baseline.json"

# The study conditions (suite: ~/Downloads/rigel_runs/quick_3to1_5mb). All capture_on, nrna_none — the
# conditions the docs measured. capture-unstranded (gdna300/ss0.50) is INCLUDED (it exists in the suite;
# the docs' "missing" note was wrong — it was only absent from the dissect default set).
CONDITIONS = {
    "stranded-0": "gdna_none_ss_0.99_nrna_none_capture_on",  # truth: zero gDNA, stranded
    "unstranded-0": "gdna_none_ss_0.50_nrna_none_capture_on",  # truth: zero gDNA, unstranded (κ≈½ floor)
    "capture": "gdna_gdna300_ss_0.99_nrna_none_capture_on",  # truth: real gDNA (oracle f_g≈0.84), stranded
    "capture-unstranded": "gdna_gdna300_ss_0.50_nrna_none_capture_on",  # real gDNA + κ≈½
}

SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}

# Diff-gate tolerance: a metric is a REGRESSION if it moves by more than max(ABS_TOL, REL_TOL·|baseline|).
ABS_TOL = 200.0
REL_TOL = 0.01


def measure(cond_dir: str, rebuild: bool) -> dict:
    """Run production calibrate on the cached full payload; compare per-region gDNA to the by-origin oracle."""
    index, blob = build_or_load_cache(cond_dir, rebuild)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cal = calibrate(
        payload=blob["payload_full"],
        region_arrays=ra,
        strand_model=blob["strand_full"],
        gdna_fl_pmf=blob["gdna_pmf"],
        rna_fl_pmf=blob["rna_pmf"],
        config=CalibrationConfig(),
    )
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    g_pr = np.asarray(cal.mass_gdna_contained, float)
    err = g_pr - g_or  # + = OVER-call (phantom), − = UNDER-call (siphon)

    rdf = index.region_df.reset_index(drop=True)
    sig = rdf["signature"].to_numpy()
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    by_class = {}
    for s in (3, 1, 2, 0):
        for t in (2, 1, 0):
            m = (scl == s) & (tcl == t)
            if not m.any():
                continue
            by_class[f"{SC[s]}.{TC[t]}"] = {
                "n": int(m.sum()),
                "oracle": float(g_or[m].sum()),
                "prod": float(g_pr[m].sum()),
                "net": float(err[m].sum()),
            }
    return {
        "condition": cond_dir,
        "total_oracle": float(g_or.sum()),
        "total_prod": float(g_pr.sum()),
        "phantom": float(np.maximum(err, 0.0).sum()),  # gDNA OVER-call
        "siphon": float(np.maximum(-err, 0.0).sum()),  # gDNA UNDER-call
        "net_err": float(err.sum()),
        "gdna_density_global": float(cal.gdna_density_global),
        "rna_sense_frac": float(cal.rna_sense_frac),
        "ambig_exon_gdna": by_class.get("AMBIG.exon", {}).get("prod", 0.0),
        "by_class": by_class,
    }


def _fmt(x: float) -> str:
    return f"{x:,.0f}"


def _regression(cur: float, base: float) -> bool:
    return abs(cur - base) > max(ABS_TOL, REL_TOL * abs(base))


def print_report(results: dict, baseline: dict | None) -> int:
    """Print the per-condition signed table; diff vs baseline if present. Returns # of regressed metrics."""
    n_reg = 0
    tracked = ["total_prod", "phantom", "siphon", "net_err", "ambig_exon_gdna"]
    for name, r in results.items():
        print("\n" + "=" * 100)
        print(f"  {name}   [{r['condition']}]")
        print(f"  oracle gDNA(contained)={_fmt(r['total_oracle'])}   "
              f"gdna_density_global={r['gdna_density_global']:.4g}   rna_sense_frac={r['rna_sense_frac']:.3f}")
        print("=" * 100)
        base = (baseline or {}).get(name)
        hdr = f"  {'metric':>18} {'current':>14}"
        if base:
            hdr += f" {'baseline':>14} {'Δ':>14}  flag"
        print(hdr)
        for k in tracked:
            line = f"  {k:>18} {_fmt(r[k]):>14}"
            if base:
                d = r[k] - base.get(k, 0.0)
                flag = ""
                if _regression(r[k], base.get(k, 0.0)):
                    flag = "  <-- MOVED"
                    n_reg += 1
                line += f" {_fmt(base.get(k, 0.0)):>14} {d:>+14,.0f}{flag}"
            print(line)
        # by-class breakdown (net error per strand_class × region_type)
        print(f"\n  {'class.type':>16} {'nreg':>5} {'g_oracle':>12} {'g_prod':>12} {'net_err':>12}")
        for key, c in r["by_class"].items():
            print(f"  {key:>16} {c['n']:>5} {_fmt(c['oracle']):>12} {_fmt(c['prod']):>12} {c['net']:>+12,.0f}")
    return n_reg


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-baseline", action="store_true", help="(re)record the baseline JSON")
    ap.add_argument("--rebuild", action="store_true", help="force re-scan of the BAMs (after a scan change)")
    ap.add_argument("--only", default="", help="comma-separated subset of condition names")
    args = ap.parse_args()

    only = set(s.strip() for s in args.only.split(",") if s.strip())
    conds = {k: v for k, v in CONDITIONS.items() if not only or k in only}

    results = {}
    for name, cond_dir in conds.items():
        print(f"[measure] {name} ({cond_dir}) ...", file=sys.stderr)
        results[name] = measure(cond_dir, args.rebuild)

    baseline = None
    if BASELINE_PATH.exists() and not args.update_baseline:
        baseline = json.loads(BASELINE_PATH.read_text())

    n_reg = print_report(results, baseline)

    if args.update_baseline:
        BASELINE_PATH.write_text(json.dumps(results, indent=2))
        print(f"\n[baseline] recorded {len(results)} conditions → {BASELINE_PATH}")
    elif baseline is not None:
        print(f"\n[gate] {n_reg} metric(s) moved beyond tol (abs={ABS_TOL:.0f}, rel={REL_TOL:.0%}) vs baseline.")
        print("       (Step-0: this records movement; no judgment — Step 1 & Step 2 are judged together.)")
    else:
        print("\n[gate] no baseline on disk; re-run with --update-baseline to record one.")


if __name__ == "__main__":
    main()
