#!/usr/bin/env python3
"""How does each message policy score, per condition, against certified truth?

The instrument a message-policy change is judged by. It runs every named policy over each cached condition of
a panel and reports the whole-library gDNA error in fragments, `sum |estimate - truth|` against the certified
per-slot truth in `slot_truth.npz`, split by axis (region / boundary), one row per condition. It runs no EM and
re-scans nothing: each condition is read from its cached scan in the drained frame the truth is certified in.

Read the strand halves apart, never pooled: unstranded rows (`ss 0.50`) must win against `silent`, the floor;
stranded rows (`ss 0.99`) must do minimal harm, so a ratio near 1.00x passes. A panel total hides a sign flip
between them, and a row at any other specificity (the test panel's `ss 0.70`) is counted apart, against
neither bar. `--panel test` is the development loop (seconds) and `--panel ladder` the shipping judgement; a toy
and the panel can rank differently, so a claim names its substrate.

`--by-class` sums the same error per node class (certified stratum, boundaries split by terminus flag, exons by
reach: licensed intron face / edge only / walled) to show where a policy's remaining error sits. On the test
chromosome the capture-OFF rows are dominated by the designed shadow-transcription floor, identical in every
arm, so read the capture-ON rows there. `--set SECTION.FIELD=VALUE` (repeatable) applies a config value on top
of every policy's fields, spelled as `calibration_vs_oracle.py` takes it.

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

from rigel.calibration.calibrate import calibrate
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain
from rigel.calibration.splice_graph import FLAG_TERMINUS
from rigel.config import CalibrationConfig, PipelineConfig
from rigel.index import TranscriptIndex
from _scan_cache import calibration_inputs, read_scan_cache

from _shared import PANELS, add_set_option, configured, refuse_unconsumed, sibling, strandedness

#: policy name -> the `CalibrationConfig` fields that install it
POLICIES = {
    "silent": dict(message_policy="silent"),
    "transfer": dict(message_policy="transfer"),
}

#: the `PipelineConfig` sections this instrument hands to the code it measures: `calibrate` alone, on a cached scan
CONSUMED_SECTIONS = frozenset({"calibration"})
#: the keys every policy sets itself; a `--set` of one would make every row the same policy
VARIED_KEYS = frozenset(f"calibration.{k}" for fields in POLICIES.values() for k in fields)


def _slot_classes(truth: dict, payload, boundary_flags) -> np.ndarray:
    """Each slot's node class: its certified stratum, a boundary tagged ``[term]`` when a transcript
    terminus sits on it, an exon tagged by its reach: ``(licensed intron face)`` when some adjacent
    intron|exon face carries no terminus, else ``(edge only)`` when it has a gene edge, else
    ``(walled: no licensed face)``."""
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


def _class_errors(slots: dict, classes: np.ndarray, region: np.ndarray, boundary: np.ndarray) -> dict:
    """``{class: (slots, mass, error)}``: one policy's ``|estimate − truth|`` summed per node class."""
    kind_s = np.asarray(slots["kind"])
    obj_s = np.asarray(slots["obj"], np.int64)
    mass_s = np.asarray(slots["count"], np.float64)
    is_r = kind_s == REGION
    est = np.zeros(kind_s.shape[0])
    est[is_r] = region[obj_s[is_r]]
    est[~is_r] = boundary[obj_s[~is_r]]
    err = np.abs(est - np.asarray(slots["n_gdna"], np.float64))
    return {
        k: (int((classes == k).sum()), float(mass_s[classes == k].sum()), float(err[classes == k].sum()))
        for k in np.unique(classes)
    }


def score_condition(index, panel_dir, condition, policies, by_class=False, settings=()):
    """One condition, every policy: `sum |estimate - truth|` per axis, in fragments — and, with
    ``by_class``, the same error summed per node class (``rows[name]["classes"]``, each value
    ``(slots, mass, error)``). ``settings`` are ``--set`` specs applied on top of every policy's
    fields."""
    truth_dir = panel_dir / "oracle_cache" / condition
    cache = read_scan_cache(panel_dir / "scan_cache" / condition, index)
    # the drained frame: `calibration_inputs` drains at the production seed and builds the production fl
    # models, the frame `slot_truth.npz` is certified in, so estimate and truth count the same fragments
    kw = calibration_inputs(cache, index)
    truth = sibling("_oracle_arms.py").certified_gdna(truth_dir)
    if by_class:
        slots = dict(np.load(truth_dir / "slot_truth.npz", allow_pickle=True))
        classes = _slot_classes(slots, kw["payload"], kw["boundary_flags"])
    rows = {}
    for name in policies:
        started = time.perf_counter()
        config = configured(settings, PipelineConfig(calibration=dataclasses.replace(CalibrationConfig(), **POLICIES[name])))
        result = calibrate(config=config.calibration, **kw)
        region = np.asarray(result.count_gdna_region, np.float64)
        boundary = np.asarray(result.count_gdna_boundary, np.float64)
        rows[name] = dict(
            region=float(np.abs(region - truth["region"]).sum()),
            boundary=float(np.abs(boundary - truth["boundary"]).sum()),
            seconds=time.perf_counter() - started,
        )
        rows[name]["total"] = rows[name]["region"] + rows[name]["boundary"]
        if by_class:
            rows[name]["classes"] = _class_errors(slots, classes, region, boundary)
    return rows


# ── the report ───────────────────────────────────────────────────────────────────────────────────


def _print_totals(table: dict, condition: str, policies: list[str], width: int) -> None:
    """One condition's line: each policy's total, and its ratio to `silent`, the floor."""
    rows = table[condition]
    cells = "  ".join(f"{rows[p]['total']:>12,.0f}" for p in policies)
    ratios = ""
    if "silent" in policies:
        floor = rows["silent"]["total"]
        ratios = "   " + " ".join(
            f"{p}={rows[p]['total'] / floor:.2f}x"
            for p in policies
            if p != "silent" and floor > 0
        )
    print(f"{condition:<{width}}  {cells}{ratios}", flush=True)


def _print_by_class(table: dict, policies: list[str]) -> None:
    print("\nby NODE CLASS — slots, mass, then each policy's |err| and its share of that policy's row:")
    for condition, rows in table.items():
        print(f"\n{condition}")
        names = [k for k in sorted(rows[policies[0]]["classes"])]
        names.sort(key=lambda k: -rows[policies[-1]]["classes"][k][2])
        print(f"   {'class':<38}{'slots':>8}{'mass':>12}" + "".join(f"{p:>12}{'share':>7}" for p in policies))
        for k in names:
            n, mass, _ = rows[policies[0]]["classes"][k]
            cells = "".join(
                f"{rows[p]['classes'][k][2]:>12,.0f}"
                f"{(rows[p]['classes'][k][2] / rows[p]['total'] if rows[p]['total'] else 0):>7.1%}"
                for p in policies
            )
            print(f"   {k:<38}{n:>8,}{mass:>12,.0f}{cells}")


def _print_bars(table: dict, policies: list[str]) -> None:
    """Each policy against `silent` on the two strand halves, counted separately."""
    print("\nthe two bars, counted separately (never pooled):")
    # a library in neither strand half (the test panel's ss 0.70) meets neither bar: counted apart
    apart = sorted({strandedness(c) for c in table} - {"unstranded", "stranded"})
    for half, want in (("unstranded", "WIN"), ("stranded", "minimal harm"),
                       *((a, "no bar: apart") for a in apart)):
        rows = [c for c in table if strandedness(c) == half]
        if not rows:
            continue
        for p in policies:
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


def _conditions(panel_dir: Path, given) -> list[str]:
    """The named conditions, else every certified one of the panel."""
    oracle = panel_dir / "oracle_cache"
    if not oracle.is_dir():
        raise SystemExit(
            f"no oracle caches at {oracle} — build the panel first "
            f"(`scripts/sim/panel.py status` names the next stage)"
        )
    conditions = given or sorted(
        d.name for d in oracle.iterdir() if (d / "slot_truth.npz").exists()
    )
    if not conditions:
        raise SystemExit(f"no certified conditions under {oracle} (missing slot_truth.npz)")
    return conditions


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
    add_set_option(ap, "override one config field on top of every policy's fields, typed from the field; "
                       "repeatable (e.g. --set calibration.sweep_logodds_step=0.1)")
    args = ap.parse_args()
    refuse_unconsumed(args.settings, CONSUMED_SECTIONS, "policy_benchmark.py", varied=VARIED_KEYS)

    for name in args.policies:
        if name not in POLICIES:
            raise SystemExit(f"unknown policy {name!r} — expected any of {sorted(POLICIES)}")
    index_dir, panel_dir = PANELS[args.panel]
    conditions = _conditions(panel_dir, args.conditions)
    index = TranscriptIndex.load(str(index_dir))

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
        table[condition] = score_condition(
            index, panel_dir, condition, args.policies, by_class=args.by_class, settings=tuple(args.settings)
        )
        _print_totals(table, condition, args.policies, width)

    print("\nper-axis (region / boundary):")
    for condition, rows in table.items():
        detail = "   ".join(
            f"{p}: {rows[p]['region']:>10,.0f} /{rows[p]['boundary']:>10,.0f}"
            for p in args.policies
        )
        print(f"{condition:<{width}}  {detail}")
    if args.by_class:
        _print_by_class(table, args.policies)
    if "silent" in args.policies and len(args.policies) > 1:
        _print_bars(table, args.policies)
    seconds = sum(r[p]["seconds"] for r in table.values() for p in args.policies)
    print(f"\n{len(conditions)} condition(s) x {len(args.policies)} policies — {seconds:.1f}s of calibrate")
    return 0


if __name__ == "__main__":
    sys.exit(main())
