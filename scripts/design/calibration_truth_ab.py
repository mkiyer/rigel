#!/usr/bin/env python
"""How does the deliverable, the library's gDNA fraction, score against truth, and what is perfecting
each fragment-length PMF worth?

Off one cached scan per condition it runs `calibrate` twice, one thing varied: with the side buffer
undrained and with it drained by `rigel.pipeline._drain_side_buffer`, and scores each `f_gdna` (in
fragment units, from `CalibrationResult`'s own conserved counts) against the simulator's origin counts in
`truth_summary.json`. With `--ceiling` it adds three arms that hand `calibrate` the simulator's own
post-capture fragment-length distributions in place of the fitted ones (exact gDNA, exact RNA, both),
which prices what perfecting each length model is worth before any work to perfect it. That "length" is
the fl PMF inside the opportunity model, not a fragment-length composition channel, which this adds
nowhere. The ceiling reaches only what `calibrate` reads: the effective-length shrinkage is built by
`pipeline.py` outside every ceiling's patch point, so it is inside no ceiling number here. Judge on the
contaminated rows and never on the zero-gDNA rows alone, where truth is exactly 0 and any change that
lowers the estimate scores better (`TRAPS: zero-target-guards-are-one-sided`). Every row is stamped with
the shipped message policy, because the number moves with it and a row from a study configuration is
otherwise indistinguishable from a shipped one; there is no switch for it here. No EM runs and no
re-scan happens.

Usage::

    python scripts/design/calibration_truth_ab.py                       # the ladder's scan caches
    python scripts/design/calibration_truth_ab.py --ceiling --json out.json
    python scripts/design/calibration_truth_ab.py --scan-cache DIR --cache-subdir _main --index DIR
    python scripts/design/calibration_truth_ab.py --conditions <cond> ... --seed 1
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
DEFAULT_SCAN_CACHE = _RUNS / "suite" / "ladder" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"


def truth_f_gdna(condition_dir: Path) -> float | None:
    """The library's true gDNA fragment fraction, from the simulator's own origin counts.

    Read from ``truth_summary.json``'s ``origin_counts`` and never parsed from the condition name, which
    is one rename away from silently wrong even where it happens to spell the fraction.
    """
    path = condition_dir / "truth_summary.json"
    if not path.is_file():
        return None
    counts = json.loads(path.read_text()).get("origin_counts", {})
    gdna = float(counts.get("gdna", 0.0))
    total = gdna + float(counts.get("mrna", 0.0)) + float(counts.get("nrna", 0.0))
    return gdna / total if total > 0 else None


def truth_length_pmf(condition_dir: Path, kind: str, max_size: int) -> "np.ndarray | None":
    """The simulator's own post-capture length distribution for one origin class, as a pmf.

    Read from ``truth_fragment_lengths.tsv``, the realised distribution, not the configured
    ``frag_mean``: capture selects for length, so the configured parameters describe a library that was
    never sequenced. ``None`` when the class has no fragments (a zero-gDNA condition).
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
    """``f_gdna`` in fragment units, from ``CalibrationResult``'s own conserved counts.

    The result's field is read rather than the banks recombined: a boundary bank books incidences, not
    fragments, so a ratio over raw banks is in the wrong unit against the truth's molecule counts.
    """
    g, r = result.library_gdna_fragments, result.library_rna_fragments
    return g / (g + r) if (g + r) > 0 else 0.0


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--scan-cache", type=Path, default=DEFAULT_SCAN_CACHE,
                    help="the scan-cache ROOT — a directory of <condition>/ caches")
    # the oracle cache holds the undrained main payload one level deeper, at
    # `<oracle_cache>/<condition>/_main`, in the `write_scan_cache` layout this reads
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

    if not args.scan_cache.is_dir():
        print(f"no scan-cache dir at {args.scan_cache}", file=sys.stderr)
        return 2
    suite = args.suite or args.scan_cache.parent

    from rigel.calibration.calibrate import calibrate
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.config import CalibrationConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer
    from rigel.scan_cache import index_derived_inputs, read_scan_cache

    index = TranscriptIndex.load(str(args.index))
    derived = index_derived_inputs(index)
    config = CalibrationConfig()
    # the same de-tilt production uses, built once off the annotation, so it is identical in both
    # arms and cannot be what moved them
    crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = gdna_opportunity_from_index(index, 4096)

    def run(payload, strand_model, *, gdna_pmf=None, rna_pmf=None):
        """Calibrate. ``gdna_pmf`` / ``rna_pmf`` override the fitted model; that is the ceiling arm."""
        fl = build_fl_models(payload, sj_opportunity=crossing, gdna_opportunity=gdna_opp)
        return calibrate(
            payload=payload,
            strand_model=strand_model,
            gdna_fl_pmf=fl.gdna_pmf if gdna_pmf is None else gdna_pmf,
            rna_fl_pmf=fl.rna_pmf if rna_pmf is None else rna_pmf,
            config=config,
            **derived,
        )

    names = args.conditions or sorted(p.name for p in args.scan_cache.iterdir() if p.is_dir())
    rows = []
    for name in names:
        root = args.scan_cache / name
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
            # part of the measurement, not metadata: the number moves with the message policy, and a
            # saved row without it cannot be attributed to a configuration
            "messages": "on" if config.message_policy != "silent" else "off",
            "truth_f_gdna": truth,
            "undrained_f_gdna": before,
            "drained_f_gdna": after,
            "held": int(cache.payload.deferred.n_fragments),
        }
        if args.ceiling:
            max_size = int(drained.max_length)
            exact_g = truth_length_pmf(suite / name, "gdna", max_size)
            exact_r = truth_length_pmf(suite / name, "rna", max_size)
            # a zero-gDNA condition has no gDNA truth histogram, so its exact-gDNA arm does not exist;
            # reported as None rather than silently falling back to the fitted model, which would read
            # as "the ceiling arm changed nothing"
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
    messages = "on" if config.message_policy != "silent" else "off"
    print(
        f"═══ ⭐ f_gdna against TRUTH — undrained (what shipped before the drain) vs drained "
        f"· messages={messages} ═══"
    )
    print(
        f"{'condition':<44} {'truth':>7} {'undrained':>10} {'err':>8} {'⭐ drained':>11} {'err':>8} {'move':>8}"
    )
    print("-" * 106)
    for r in rows:
        t = r["truth_f_gdna"]
        if t is None:
            print(f"{r['condition']:<44} {'—':>7}")
            continue
        # absolute error, not relative: truth is 0 exactly on the zero-gDNA arm, where a relative error
        # is undefined
        eb, ea = r["undrained_f_gdna"] - t, r["drained_f_gdna"] - t
        print(
            f"{r['condition']:<44} {t:>7.4f} {r['undrained_f_gdna']:>10.4f} {eb:>+8.4f} "
            f"{r['drained_f_gdna']:>11.4f} {ea:>+8.4f} {abs(ea) - abs(eb):>+8.4f}"
        )
    print("   `move` is |err| after − |err| before: NEGATIVE is an improvement.")
    print("   ⛔ Judge on the CONTAMINATED rows (g05/g50/g98). The g00 rows are saturated at truth = 0")
    print("      exactly, so any change that lowers the estimate 'improves' them —")
    print("      TRAPS: zero-target-guards-are-one-sided, which has reversed a verdict here.")

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
            print("   ⚠ A ceiling is what perfecting a channel is WORTH, not a result.")
            print("     TRAPS: measure-the-ceiling-first — it is available whenever the simulator writes")
            print("     truth, and it costs one afternoon.")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
