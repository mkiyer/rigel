#!/usr/bin/env python3
"""How does each message policy score, per condition, against certified truth? This is the standing
policy benchmark and the instrument a message-policy change is judged by. It runs every named
policy over each cached condition of a panel and reports the whole-library gDNA error in
fragments, `sum |estimate - truth|` against the certified per-slot truth in `slot_truth.npz`,
split by axis (region / boundary), one row per condition. It runs no EM and re-scans nothing:
each condition is read from its cached scan in the drained frame the truth is certified in.
Read the two halves separately and never pool them: unstranded rows (`ss 0.50`) are where a
policy must win against `silent`, the measured floor; stranded rows (`ss 0.99`) are where it
must do minimal harm, so a ratio near 1.00x is a pass; a panel total hides a sign flip between
them. `--panel test` is the development loop (seconds); `--panel ladder` is the shipping judgement,
and a toy and the panel have inverted a ranking before, so a claim names its substrate.
`--by-class` sums the same error per node class (certified stratum, boundaries split
by terminus flag, exons by reach: licensed intron face / edge only / walled) to rank where a
policy's remaining error sits; on the test chromosome the capture-OFF rows are dominated by the
designed shadow-transcription floor, identical in every arm, so read the capture-ON rows there.
`--set SECTION.FIELD=VALUE` (repeatable) applies any config value on top of every policy's fields,
the same spelling `calibration_vs_oracle.py` takes, so a grid arm is read on the panel and on the
metric from one config value.

Usage::

    python scripts/design/policy_benchmark.py --panel test                                   # the development loop
    python scripts/design/policy_benchmark.py --panel ladder --policies silent transfer       # the shipping judgement
    python scripts/design/policy_benchmark.py --panel test --conditions gdna_g50_ss_0.50_nrna_file_capture_off
    python scripts/design/policy_benchmark.py --panel ladder --policies silent transfer --by-class
    python scripts/design/policy_benchmark.py --panel test --set calibration.sweep_logodds_step=0.1    # a config arm
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
from rigel.calibration.splice_graph import FLAG_TERMINUS  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

from _shared import set_field  # noqa: E402

RUNS = Path.home() / "Downloads" / "rigel_runs"

#: the two substrates, each: (index dir, the dir holding `oracle_cache/<condition>/`)
PANELS = {
    "test": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios"),
    # the same chromosome under the two other probe designs (junction-only and sparse probes)
    "test_junction": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_junction"),
    "test_sparse": (RUNS / "test_reference" / "idx", RUNS / "test_reference" / "scenarios_probes_sparse"),
    "ladder": (RUNS / "suite" / "rigel_index", RUNS / "suite" / "ladder"),
}

#: policy name -> the `CalibrationConfig` fields that install it
POLICIES = {
    "silent": dict(message_policy="silent"),
    "transfer": dict(message_policy="transfer"),
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


def _slot_classes(truth: dict, payload, boundary_flags) -> np.ndarray:
    """Each slot's node class: its certified stratum, a boundary tagged ``[term]`` when a transcript
    terminus sits on it, an exon tagged by its REACH — ``(licensed intron face)`` when some adjacent
    intron|exon face carries no terminus, else ``(edge only)`` when it has a gene edge, else
    ``(walled: no licensed face)`` — the classes the message rebuild's holes are named in."""
    kind = np.asarray(truth["kind"])
    obj = np.asarray(truth["obj"], np.int64)
    strata = np.asarray(truth["stratum"]).astype(str)
    chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    is_b = kind == BOUNDARY
    term = np.zeros(kind.shape[0], bool)
    term[is_b] = (np.asarray(boundary_flags, np.uint16)[obj[is_b]] & FLAG_TERMINUS) != 0
    cls = strata.astype(object)
    for i in np.flatnonzero(is_b & term):
        cls[i] = strata[i] + " [term]"
    for e in np.flatnonzero(strata == "R exon"):
        faces = [b for b in (left[e], right[e]) if b >= 0 and is_b[b]]
        if any(strata[b] == "B exon|intron" and not term[b] for b in faces):
            cls[e] = "R exon (licensed intron face)"
        elif any(strata[b] == "B gene edge" for b in faces):
            cls[e] = "R exon (edge only)"
        else:
            cls[e] = "R exon (walled: no licensed face)"
    return cls.astype(str)


def score_condition(
    index, region_arrays, sj, boundary_flags, cache_dir, policies, by_class=False, settings=()
):
    """One condition, every policy: `sum |estimate - truth|` per axis, in fragments — and, with
    ``by_class``, the same error summed per node class (``rows[name]["classes"]``, each value
    ``(slots, mass, error)``). ``settings`` are ``--set`` specs applied on top of every policy's
    fields."""
    cache = read_scan_cache(cache_dir / "_main", index)
    # the drained frame: `calibration_inputs` drains at the production seed and builds the
    # production fl models — the frame `slot_truth.npz` is certified in, so estimate and truth
    # speak one tally.
    kw = calibration_inputs(cache, index)
    payload = kw["payload"]
    kwargs = dict(
        region_arrays=region_arrays,
        strand_model=kw["strand_model"],
        gdna_fl_pmf=kw["gdna_fl_pmf"],
        rna_fl_pmf=kw["rna_fl_pmf"],
        sj=sj,
        boundary_flags=boundary_flags,
    )
    truth = _truth_gdna(cache_dir)
    if by_class:
        slots = dict(np.load(cache_dir / "slot_truth.npz", allow_pickle=True))
        classes = _slot_classes(slots, payload, boundary_flags)
        kind_s = np.asarray(slots["kind"])
        obj_s = np.asarray(slots["obj"], np.int64)
        truth_s = np.asarray(slots["n_gdna"], np.float64)
        mass_s = np.asarray(slots["count"], np.float64)
    rows = {}
    for name in policies:
        started = time.perf_counter()
        config = PipelineConfig(calibration=dataclasses.replace(CalibrationConfig(), **POLICIES[name]))
        for spec in settings:
            config = set_field(config, spec)
        result = calibrate(payload=payload, config=config.calibration, **kwargs)
        region = np.asarray(result.count_gdna_region, np.float64)
        boundary = np.asarray(result.count_gdna_boundary, np.float64)
        rows[name] = dict(
            region=float(np.abs(region - truth["region"]).sum()),
            boundary=float(np.abs(boundary - truth["boundary"]).sum()),
            seconds=time.perf_counter() - started,
        )
        rows[name]["total"] = rows[name]["region"] + rows[name]["boundary"]
        if by_class:
            is_r = kind_s == REGION
            est = np.zeros(kind_s.shape[0])
            est[is_r] = region[obj_s[is_r]]
            est[~is_r] = boundary[obj_s[~is_r]]
            err = np.abs(est - truth_s)
            rows[name]["classes"] = {
                k: (int((classes == k).sum()), float(mass_s[classes == k].sum()), float(err[classes == k].sum()))
                for k in np.unique(classes)
            }
    return rows


def _stranded(condition: str) -> bool:
    """`ss 0.99` conditions are strand-specific; `ss 0.50` are unstranded."""
    return "_ss_0.99_" in condition


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", choices=sorted(PANELS), default="test")
    ap.add_argument("--policies", nargs="+", default=["silent", "transfer"])
    ap.add_argument("--conditions", nargs="+", default=None, help="default: all cached")
    ap.add_argument(
        "--by-class",
        action="store_true",
        help="also sum each policy's error per node class (stratum, terminus flag, exon reach)",
    )
    ap.add_argument(
        "--set",
        dest="settings",
        action="append",
        default=[],
        metavar="SECTION.FIELD=VALUE",
        help="override one config field on top of every policy's fields, typed from the field; "
        "repeatable (e.g. --set calibration.sweep_logodds_step=0.1)",
    )
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
    for spec in args.settings:
        print(f"⭐ --set {spec} on every policy")
    print("   unstranded rows are where a policy must WIN; stranded rows are where it must do")
    print("   as little HARM as possible. Never pool them.\n")
    print(header)
    print("-" * len(header))

    table = {}
    for condition in conditions:
        rows = score_condition(
            index,
            region_arrays,
            sj,
            boundary_flags,
            oracle / condition,
            args.policies,
            by_class=args.by_class,
            settings=tuple(args.settings),
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

    if args.by_class:
        print("\nby NODE CLASS — slots, mass, then each policy's |err| and its share of that policy's row:")
        for condition, rows in table.items():
            print(f"\n{condition}")
            names = [k for k in sorted(rows[args.policies[0]]["classes"])]
            names.sort(key=lambda k: -rows[args.policies[-1]]["classes"][k][2])
            print(f"   {'class':<38}{'slots':>8}{'mass':>12}" + "".join(f"{p:>12}{'share':>7}" for p in args.policies))
            for k in names:
                n, mass, _ = rows[args.policies[0]]["classes"][k]
                cells = "".join(
                    f"{rows[p]['classes'][k][2]:>12,.0f}"
                    f"{(rows[p]['classes'][k][2] / rows[p]['total'] if rows[p]['total'] else 0):>7.1%}"
                    for p in args.policies
                )
                print(f"   {k:<38}{n:>8,}{mass:>12,.0f}{cells}")

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
