#!/usr/bin/env python
"""⭐ **THE DELIVERABLE, SCORED AGAINST TRUTH**: does the second pass reach the gDNA/RNA composition?

    Baseline it replaces: `LEDGER.md` S5.f, which measured `f_gdna` before ANY of the fragment-length
    work — "within 6 % of truth on 3 of 4 gdna100 conditions", with `gdna100 ss0.50 capture_on` 25 % low.

⛔ **THIS IS THE QUESTION THE WHOLE FRAGMENT-LENGTH TRACK EXISTS TO ANSWER.** C0–C2.6 and the second pass
(P0–P4.2) are all upstream plumbing: one definition of fragment length, then an accurate one, then an
*unbiased* one. None of that is the product. The product is the library's composition, and
`ACCUMULATOR_DESIGN.md` §7.2 measured the coupling — **a 10 % length-model error is worth 0.010–0.026 of
composition** — so a length error that went from +27 % to +0.00 % should be visible here or the coupling
is not what it was thought to be.

⭐ **ONE THING VARIED.** The same cached scan, the same index, the same config; the only difference is
whether the side buffer has been **drained** before calibration reads the tally. The undrained arm is
exactly what shipped before P4.

⛔ **DO NOT SCORE ON THE ZERO-gDNA ARM ALONE** (`CARRY_FORWARD.md` §3 trap 19). Truth there is
`f_gdna = 0` *exactly*, so any change that lowers the estimate scores better — a one-sidedness that has
already reversed a verdict in this project once. The **gdna100** arm carries the real signal: truth is
5 M mRNA against 5 M gDNA fragments, so `f_gdna = 0.5`.

Usage::

    python scripts/design/calibration_truth_ab.py [--index DIR] [--pilot DIR] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def truth_f_gdna(condition_dir: Path) -> float | None:
    """The library's TRUE gDNA fragment fraction, from the simulator's own origin counts.

    ⚠ Read from ``truth_summary.json``'s ``origin_counts`` rather than from the condition NAME. "gdna100"
    is a rate knob, not a fraction — it happens to give 5 M gDNA against 5 M mRNA, i.e. 0.5, and inferring
    that from the string would be one rename away from silently wrong.
    """
    path = condition_dir / "truth_summary.json"
    if not path.is_file():
        return None
    counts = json.loads(path.read_text()).get("origin_counts", {})
    gdna = float(counts.get("gdna", 0.0))
    total = gdna + float(counts.get("mrna", 0.0)) + float(counts.get("nrna", 0.0))
    return gdna / total if total > 0 else None


def f_gdna_of(result) -> float:
    """``f_gdna`` as `LEDGER.md` S5.f defined it, so the two baselines are comparable.

    ⚠ The sum runs over nodes AND edges. gDNA lives on both — contained in a node or crossing a line —
    and summing only one axis reports a library's gDNA as a fraction of part of itself.
    """
    g = float(np.asarray(result.mass_gdna_node).sum() + np.asarray(result.mass_gdna_edge).sum())
    r = float(np.asarray(result.mass_rna_node).sum() + np.asarray(result.mass_rna_edge).sum())
    return g / (g + r) if (g + r) > 0 else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None, help="where the truth files live")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=0, help="the drain's multinomial seed")
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    if not args.pilot.is_dir():
        print(f"no pilot scan-cache dir at {args.pilot}", file=sys.stderr)
        return 2
    suite = args.suite or args.pilot.parent

    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.fl import build_fl_models
    from rigel.config import CalibrationConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer
    from rigel.scan_cache import index_derived_inputs, read_scan_cache

    index = TranscriptIndex.load(str(args.index))
    derived = index_derived_inputs(index)
    config = CalibrationConfig()

    def run(payload, strand_model):
        fl = build_fl_models(payload)
        return calibrate(
            payload=payload,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf,
            rna_fl_pmf=fl.rna_pmf,
            config=config,
            **derived,
        )

    names = args.conditions or sorted(p.name for p in args.pilot.iterdir() if p.is_dir())
    rows = []
    for name in names:
        cache = read_scan_cache(args.pilot / name, index)
        truth = truth_f_gdna(suite / name)
        start = time.perf_counter()
        before = f_gdna_of(run(cache.payload, cache.strand_model))
        drained = _drain_side_buffer(cache.payload, index, cache.strand_model, seed=args.seed)
        after = f_gdna_of(run(drained, cache.strand_model))
        rows.append(
            {
                "condition": name,
                "truth_f_gdna": truth,
                "undrained_f_gdna": before,
                "drained_f_gdna": after,
                "held": int(cache.payload.deferred.n_fragments),
                "seconds": time.perf_counter() - start,
            }
        )
        print(f"  {name:<44} done in {rows[-1]['seconds']:.0f} s")

    print()
    print("═══ ⭐ f_gdna against TRUTH — undrained (what shipped before P4) vs drained ═══")
    print(
        f"{'condition':<44} {'truth':>7} {'undrained':>10} {'err':>8} {'⭐ drained':>11} {'err':>8} {'move':>8}"
    )
    print("-" * 106)
    for r in rows:
        t = r["truth_f_gdna"]
        if t is None:
            print(f"{r['condition']:<44} {'—':>7}")
            continue
        # ⚠ Absolute error, not relative: truth is 0 exactly on the zero-gDNA arm, where a relative error
        # is undefined and a ratio would read as a division blow-up rather than a good answer.
        eb, ea = r["undrained_f_gdna"] - t, r["drained_f_gdna"] - t
        print(
            f"{r['condition']:<44} {t:>7.4f} {r['undrained_f_gdna']:>10.4f} {eb:>+8.4f} "
            f"{r['drained_f_gdna']:>11.4f} {ea:>+8.4f} {abs(ea) - abs(eb):>+8.4f}"
        )
    print("   `move` is |err| after − |err| before: NEGATIVE is an improvement.")
    print("   ⛔ Judge on the gdna100 rows. The zero-gDNA rows are saturated at truth = 0 exactly, so any")
    print("      change that lowers the estimate 'improves' them — trap 19, which has reversed a verdict here.")

    scored = [r for r in rows if r["truth_f_gdna"]]
    if scored:
        contaminated = [r for r in scored if r["truth_f_gdna"] > 0.1]
        if contaminated:
            eb = np.mean([abs(r["undrained_f_gdna"] - r["truth_f_gdna"]) for r in contaminated])
            ea = np.mean([abs(r["drained_f_gdna"] - r["truth_f_gdna"]) for r in contaminated])
            print()
            print(
                f"⭐ mean |error| on the {len(contaminated)} CONTAMINATED conditions: "
                f"{eb:.4f} → {ea:.4f}  ({100 * (ea - eb) / eb:+.1f} %)"
            )

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
