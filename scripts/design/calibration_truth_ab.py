#!/usr/bin/env python
"""⭐ **THE DELIVERABLE, SCORED AGAINST TRUTH**: does the second pass reach the gDNA/RNA composition?

⛔ **THIS IS THE QUESTION THE WHOLE FRAGMENT-LENGTH TRACK EXISTS TO ANSWER.** The length work and the
second pass are upstream plumbing: one definition of fragment length, then an accurate one, then an
*unbiased* one. None of that is the product. The product is the library's composition, and the coupling
was measured at **a 10 % length-model error is worth 0.010–0.026 of composition** — so a length error
that went from +27 % to +0.00 % should be visible here, or the coupling is not what it was thought to be.

⚠ **Three sentences of this docstring were WRECKAGE from the numbered-label rename** — a dangling
"Baseline it replaces:," with no referent, and a rule citation that had become
"TRAPS: two-divisors-opposite-sign–TRAPS: pure-and-length-censored.6". They named documents and labels
that no longer exist and are deleted rather than guessed at (2026-08-11).

⭐ **ONE THING VARIED.** The same cached scan, the same index, the same config; the only difference is
whether the side buffer has been **drained** before calibration reads the tally. The undrained arm is
exactly what shipped before P4.

⛔ **DO NOT SCORE ON THE ZERO-gDNA ARM ALONE**. Truth there is
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


def truth_length_pmf(condition_dir: Path, kind: str, max_size: int) -> "np.ndarray | None":
    """The simulator's OWN post-capture length distribution for one origin class, as a pmf.

    ⭐⭐ **THIS IS THE CEILING INSTRUMENT** (`docs/TRAPS.md` measure-the-ceiling-first). Handing `calibrate`
    the right answer for one channel says what *perfecting* that channel is worth, before any of the work
    to perfect it is done — and it is available whenever the simulator writes truth. It is what showed
    that a perfect RNA length model was worth 0.0004 while the gDNA pool nobody was ranking was worth
    22 %: an A/B tells you whether a change helped, a ceiling tells you whether to start.

    ⚠ Read from ``truth_fragment_lengths.tsv``, which is **post-capture empirical** — the realised
    distribution, not the configured ``frag_mean``. Capture selects for length, so the configured
    parameters describe a library that was never sequenced
    (`docs/TRAPS.md` capture-selects-for-length).
    """
    path = condition_dir / "truth_fragment_lengths.tsv"
    if not path.is_file():
        return None
    pmf = np.zeros(max_size + 1, dtype=np.float64)
    with open(path) as handle:
        next(handle)
        for line in handle:
            row_kind, length_text, count_text, _fraction = line.rstrip("\n").split("\t")
            if row_kind != kind:
                continue
            length = int(length_text)
            if 0 <= length <= max_size:
                pmf[length] += float(count_text)
    total = pmf.sum()
    return pmf / total if total > 0 else None


def f_gdna_of(result) -> float:
    """``f_gdna`` in FRAGMENT units — ``CalibrationResult``'s own conserved counts, the one definition.

    ⛔⛔ **THIS USED TO SUM INCIDENCES AND SCORE THEM AGAINST A FRAGMENT COUNT.** It added
    ``mass_gdna_node + mass_gdna_edge`` over the raw banks, where an edge term books ``max(K,1)`` per
    fragment, and compared the ratio to ``truth_summary.json``'s ``origin_counts`` — which are real
    molecules. The units did not match across the subtraction, and because the two components' K
    inflations differ they did not cancel: on ladder g50 capture_off the incidence ratio reads
    **0.3851** against a truth of **0.5085**, while the conserved counts reproduce **0.5085** exactly.
    It also silently omitted the junction axis, so it disagreed with ``pipeline.py``'s version too.

    ⭐ Reading the result's own field is the point: the count is assembled once, in ``calibrate``, where
    each axis is converted by its own population's ``mass / count``. A consumer recombining the banks
    itself is how the tree came to hold three different answers to this question.
    """
    g, r = result.library_gdna_fragments, result.library_rna_fragments
    return g / (g + r) if (g + r) > 0 else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    # ⚠ `--pilot` is a HISTORICAL NAME: it is simply the scan-cache ROOT, and it has pointed at panels
    # other than the pilot for a long time. Kept rather than renamed because instruments and notes call
    # it by this name; what it means is documented here and in the module docstring.
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT,
                    help="the scan-cache ROOT — a directory of <condition>/ caches")
    # ⭐ THE LADDER'S FULL SCAN IS CACHED ONE LEVEL DEEPER. `pass0_vs_oracle` writes the undrained main
    # payload to `<oracle_cache>/<condition>/_main`, in exactly the `write_scan_cache` layout this
    # reads — so the 36-condition panel needs a subdirectory, not a second 27 GB copy of the scans.
    ap.add_argument("--cache-subdir", default="",
                    help="subdirectory under each condition holding the cache; use `_main` to read an "
                         "oracle_cache built by pass0_vs_oracle.py / prior_vs_oracle.py")
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None, help="where the truth files live")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--seed", type=int, default=0, help="the drain's multinomial seed")
    ap.add_argument(
        "--ceiling",
        action="store_true",
        help="also run the three TRUTH-pmf arms: exact gDNA, exact RNA, both. This is the scoping "
        "number for the whole length phase — what perfecting each length model is worth.",
    )
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    if not args.pilot.is_dir():
        print(f"no pilot scan-cache dir at {args.pilot}", file=sys.stderr)
        return 2
    suite = args.suite or args.pilot.parent

    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.config import CalibrationConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer
    from rigel.scan_cache import index_derived_inputs, read_scan_cache

    index = TranscriptIndex.load(str(args.index))
    derived = index_derived_inputs(index)
    config = CalibrationConfig()
    # ⚠ The same de-tilt production uses, built once off the annotation, so it is identical in both
    # arms and cannot be what moved them.
    crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = gdna_opportunity_from_index(index, 4096)

    def run(payload, strand_model, *, gdna_pmf=None, rna_pmf=None):
        """Calibrate. ``gdna_pmf`` / ``rna_pmf`` override the fitted model — that is the ceiling arm."""
        fl = build_fl_models(payload, junction_opportunity=crossing, gdna_opportunity=gdna_opp)
        return calibrate(
            payload=payload,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf if gdna_pmf is None else gdna_pmf,
            rna_fl_pmf=fl.rna_pmf if rna_pmf is None else rna_pmf,
            config=config,
            **derived,
        )

    names = args.conditions or sorted(p.name for p in args.pilot.iterdir() if p.is_dir())
    rows = []
    for name in names:
        root = args.pilot / name
        if args.cache_subdir:
            root = root / args.cache_subdir
        cache = read_scan_cache(root, index)
        truth = truth_f_gdna(suite / name)
        start = time.perf_counter()
        before = f_gdna_of(run(cache.payload, cache.strand_model))
        drained = _drain_side_buffer(cache.payload, index, cache.strand_model, seed=args.seed)
        after = f_gdna_of(run(drained, cache.strand_model))
        row = {
            "condition": name,
            "truth_f_gdna": truth,
            "undrained_f_gdna": before,
            "drained_f_gdna": after,
            "held": int(cache.payload.deferred.n_fragments),
        }
        if args.ceiling:
            max_size = int(drained.max_length)
            exact_g = truth_length_pmf(suite / name, "gdna", max_size)
            exact_r = truth_length_pmf(suite / name, "rna", max_size)
            # ⚠ A zero-gDNA condition has no gDNA truth histogram at all, so its exact-gDNA arm does not
            # exist. Reported as None rather than silently falling back to the fitted model, which would
            # read as "the ceiling arm changed nothing".
            row["exact_gdna_f_gdna"] = (
                f_gdna_of(run(drained, cache.strand_model, gdna_pmf=exact_g))
                if exact_g is not None
                else None
            )
            row["exact_rna_f_gdna"] = (
                f_gdna_of(run(drained, cache.strand_model, rna_pmf=exact_r))
                if exact_r is not None
                else None
            )
            row["exact_both_f_gdna"] = (
                f_gdna_of(run(drained, cache.strand_model, gdna_pmf=exact_g, rna_pmf=exact_r))
                if exact_g is not None and exact_r is not None
                else None
            )
        row["seconds"] = time.perf_counter() - start
        rows.append(row)
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

    if args.ceiling:
        print()
        print("═══ ⭐⭐ THE CEILING — calibrate handed the simulator's OWN length distribution ═══")
        print("   Same drained payload, same everything; the fragment-length pmf is the only thing varied.")
        print(
            f"{'condition':<44} {'truth':>7} {'⭐ drained':>11} {'exact gDNA':>11} "
            f"{'exact RNA':>11} {'both exact':>11}"
        )
        print("-" * 108)
        for r in rows:
            t = r["truth_f_gdna"]
            cells = "".join(
                f"{r.get(key):>11.4f}" if r.get(key) is not None else f"{'—':>11}"
                for key in ("drained_f_gdna", "exact_gdna_f_gdna", "exact_rna_f_gdna", "exact_both_f_gdna")
            )
            print(f"{r['condition']:<44} {t if t is None else f'{t:.4f}':>7}{cells}")

        contaminated = [r for r in rows if (r["truth_f_gdna"] or 0.0) > 0.1]
        if contaminated:
            print()
            print(f"⭐ mean |error| over the {len(contaminated)} CONTAMINATED conditions:")
            base = np.mean([abs(r["drained_f_gdna"] - r["truth_f_gdna"]) for r in contaminated])
            for label, key in (
                ("shipped (drained)", "drained_f_gdna"),
                ("the EXACT gDNA length distribution", "exact_gdna_f_gdna"),
                ("the EXACT RNA length distribution", "exact_rna_f_gdna"),
                ("⭐⭐ BOTH exact — the ceiling on the whole length phase", "exact_both_f_gdna"),
            ):
                values = [r for r in contaminated if r.get(key) is not None]
                if len(values) != len(contaminated):
                    print(f"   {label:<56} — (missing on {len(contaminated) - len(values)} row(s))")
                    continue
                mean = np.mean([abs(r[key] - r["truth_f_gdna"]) for r in values])
                delta = f"{100 * (mean - base) / base:+.1f} %" if base > 0 else "—"
                print(f"   {label:<56} {mean:.4f}   {delta}")
            print("   ⚠ A ceiling is what perfecting a channel is WORTH, not a result. Trap 31: it is")
            print("     available whenever the simulator writes truth, and it costs one afternoon.")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
