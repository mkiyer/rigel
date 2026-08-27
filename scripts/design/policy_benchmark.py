#!/usr/bin/env python3
"""HOW DOES EACH MESSAGE POLICY SCORE, PER CONDITION, AGAINST CERTIFIED TRUTH? — the standing
policy benchmark, on either substrate, never pooled.

⭐⭐⭐ **THIS IS THE INSTRUMENT A POLICY CHANGE IS JUDGED BY.** It runs each named policy over
every cached condition of a panel and reports the whole-library gDNA error in FRAGMENTS,
split by axis (REGION / BOUNDARY) and by stranded / unstranded, one row per condition.

⭐⭐ **THE BAR IS NOT "BEAT SILENT".** `SilentPolicy` is the measured floor and on
strand-specific data it is very hard to improve on: a sighted exon's own strand solve is
excellent, so a message can mostly only disturb it. The goal message propagation exists for is
UNSTRANDED data, where the strand channel is dead and the local answer is a default rather
than a measurement. So read the two halves of the table differently:

* **unstranded rows (`ss 0.50`)** — this is where a policy has to WIN.
* **stranded rows (`ss 0.99`)** — this is where a policy has to do as LITTLE HARM as possible.
  A ratio near 1.00x against silence is a pass; a large ratio is the thing to fix.

⛔ **NEVER POOL THE ROWS.** A panel total hides a sign flip between strata, and the two halves
above are being judged against different bars. The per-condition table is the result.

⭐ **SUBSTRATES.** `--panel test` is the method-development test chromosome (small, cached, a
few seconds for the whole sweep — the loop to develop in); `--panel ladder` is the 16-condition
benchmark ladder (the shipping judgement). ⚠ A toy and the panel have inverted a ranking
before (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`), so a claim must name its substrate
and a change is only real once the ladder agrees.

⛔ It scores the CALIBRATION result against the oracle's certified per-slot truth. It runs no
EM and re-scans nothing: every condition is read from its cached scan.

    python scripts/design/policy_benchmark.py --panel test
    python scripts/design/policy_benchmark.py --panel ladder --policies silent relay
    python scripts/design/policy_benchmark.py --panel test --conditions gdna_g50_ss_0.50_nrna_file_capture_off
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.fl import build_fl_models  # noqa: E402
from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION  # noqa: E402
from rigel.calibration.sj_opportunity import crossing_probability_from_index  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

RUNS = Path.home() / "Downloads" / "rigel_runs"

#: the two substrates, each: (index dir, the dir holding `oracle_cache/<condition>/`)
PANELS = {
    "test": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios"),
    "ladder": (RUNS / "suite" / "rigel_index", RUNS / "suite" / "ladder"),
}

#: policy name -> the `CalibrationConfig` fields that install it. `rna_anchor` is live only
#: under the relay (`DESIGN.md` §6b.3), so the other arms carry it False and say so.
POLICIES = {
    "silent": dict(message_propagation=False, rna_anchor=False),
    "relay": dict(message_propagation=True, message_policy="relay", rna_anchor=True),
    "message": dict(message_propagation=True, message_policy="message", rna_anchor=False),
}


def _truth_gdna(cache_dir: Path) -> dict:
    """The certified per-object gDNA fragment counts, by axis, from `slot_truth.npz`."""
    truth = dict(np.load(cache_dir / "slot_truth.npz"))
    out = {}
    for kind, axis in ((REGION, "region"), (BOUNDARY, "boundary")):
        is_k = np.asarray(truth["kind"]) == kind
        obj = np.asarray(truth["obj"], np.int64)[is_k]
        arr = np.zeros(int(obj.max()) + 1 if obj.size else 0)
        arr[obj] = np.asarray(truth["n_gdna"], np.float64)[is_k]
        out[axis] = arr
    return out


def score_condition(index, region_arrays, sj, boundary_flags, cache_dir, policies):
    """One condition, every policy: `sum |estimate - truth|` per axis, in fragments."""
    cache = read_scan_cache(cache_dir / "_main", index)
    payload = cache.payload
    fl = build_fl_models(
        payload,
        sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
        gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
    )
    kwargs = dict(
        region_arrays=region_arrays,
        strand_model=cache.strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        sj=sj,
        boundary_flags=boundary_flags,
    )
    truth = _truth_gdna(cache_dir)
    rows = {}
    for name in policies:
        started = time.perf_counter()
        result = calibrate(
            payload=payload,
            config=dataclasses.replace(CalibrationConfig(), **POLICIES[name]),
            **kwargs,
        )
        region = np.asarray(result.mass_gdna_region, np.float64)
        boundary = np.asarray(result.mass_gdna_boundary, np.float64)
        rows[name] = dict(
            region=float(np.abs(region - truth["region"]).sum()),
            boundary=float(np.abs(boundary - truth["boundary"]).sum()),
            seconds=time.perf_counter() - started,
        )
        rows[name]["total"] = rows[name]["region"] + rows[name]["boundary"]
    return rows


def _stranded(condition: str) -> bool:
    """`ss 0.99` conditions are strand-specific; `ss 0.50` are unstranded."""
    return "_ss_0.99_" in condition


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", choices=sorted(PANELS), default="test")
    ap.add_argument("--policies", nargs="+", default=["silent", "relay", "message"])
    ap.add_argument("--conditions", nargs="+", default=None, help="default: all cached")
    args = ap.parse_args()

    for name in args.policies:
        if name not in POLICIES:
            raise SystemExit(f"unknown policy {name!r} — expected any of {sorted(POLICIES)}")

    index_dir, panel_dir = PANELS[args.panel]
    oracle = panel_dir / "oracle_cache"
    if not oracle.is_dir():
        raise SystemExit(
            f"no oracle caches at {oracle} — build the panel first "
            f"(`scripts/sim/panel.py status` names the next stage)"
        )
    conditions = args.conditions or sorted(
        d.name for d in oracle.iterdir() if (d / "slot_truth.npz").exists()
    )
    if not conditions:
        raise SystemExit(f"no certified conditions under {oracle} (missing slot_truth.npz)")

    index = TranscriptIndex.load(str(index_dir))
    region_arrays = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    sj = build_sj_geometry_arrays(index)
    boundary_flags = build_boundary_flags_array(index)

    width = max(len(c) for c in conditions)
    header = f"{'condition':<{width}}  " + "  ".join(f"{p:>12}" for p in args.policies)
    if "silent" in args.policies:
        header += "   vs silent"
    print(f"⭐ {args.panel} panel — whole-library |gDNA estimate − truth|, in FRAGMENTS")
    print("   unstranded rows are where a policy must WIN; stranded rows are where it must do")
    print("   as little HARM as possible. Never pool them.\n")
    print(header)
    print("-" * len(header))

    table = {}
    for condition in conditions:
        rows = score_condition(
            index, region_arrays, sj, boundary_flags, oracle / condition, args.policies
        )
        table[condition] = rows
        cells = "  ".join(f"{rows[p]['total']:>12,.0f}" for p in args.policies)
        ratios = ""
        if "silent" in args.policies:
            floor = rows["silent"]["total"]
            ratios = "   " + " ".join(
                f"{p}={rows[p]['total'] / floor:.2f}x"
                for p in args.policies
                if p != "silent" and floor > 0
            )
        print(f"{condition:<{width}}  {cells}{ratios}", flush=True)

    print("\nper-axis (region / boundary):")
    for condition, rows in table.items():
        detail = "   ".join(
            f"{p}: {rows[p]['region']:>10,.0f} /{rows[p]['boundary']:>10,.0f}"
            for p in args.policies
        )
        print(f"{condition:<{width}}  {detail}")

    if "silent" in args.policies and len(args.policies) > 1:
        print("\nthe two bars, counted separately (never pooled):")
        for half, want in (("unstranded", "WIN"), ("stranded", "minimal harm")):
            rows = [c for c in table if _stranded(c) == (half == "stranded")]
            if not rows:
                continue
            for p in args.policies:
                if p == "silent":
                    continue
                better = sum(table[c][p]["total"] < table[c]["silent"]["total"] for c in rows)
                worst = max(
                    (table[c][p]["total"] / table[c]["silent"]["total"] for c in rows if table[c]["silent"]["total"] > 0),
                    default=float("nan"),
                )
                print(
                    f"   {half:<11} ({want:<13}) {p:>8}: beats silence on {better}/{len(rows)}"
                    f"   worst row {worst:.2f}x"
                )
    seconds = sum(r[p]["seconds"] for r in table.values() for p in args.policies)
    print(f"\n{len(conditions)} condition(s) x {len(args.policies)} policies — {seconds:.1f}s of calibrate")
    return 0


if __name__ == "__main__":
    sys.exit(main())
