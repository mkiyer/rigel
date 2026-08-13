#!/usr/bin/env python
"""Re-measure 's gap table against the CURRENT anchor, and score it on TRUTH.

⭐ **The falsification this area never had.** On a zero-gDNA condition every fragment is RNA, so the
unconditional anchor and the RNA pool describe **one population** and any gap between them is bias.
§2 measured that gap at **+11.6 % mean / +71.1 % sd** with the scanner's histogram as the anchor; TRAPS: a-purity-filter-is-a-length-filter
built the accumulator's own histogram and measured **+7.7 % / +32.0 %**, isolating the residual as the
junction-opportunity tilt (TRAPS: divide-by-a-probability's target). This script says which of those the *shipped* code is on.

⚠ It reads the anchor **through ``build_fl_models``**, not off the payload directly — the question is
what the tool is wired to use, not what is available to it.

⭐ **The truth panel (``--truth``) is TRAPS: pure-and-length-censored.6's gate**, The simulator writes
``truth_fragment_lengths.tsv`` beside every condition, so the library's realized fragment-length support
is known **exactly** and none of the targets below is chosen:

    G-tail   the anchor's mass above the library's TRUE ceiling must be 0, and `dropped_too_long`
             must collapse — a fragment longer than any molecule in the library is an uncut intron
    G-sd     the anchor's sd against the truth's sd
    G-gdna   ⛔ THE CONTROL. `DNA_INTERGENIC` is pure gDNA and gDNA has NO INTRONS TO MISS, so it was
             already exact to five decimals. If it MOVES, the fix reached fragments with no introns —
             which is impossible, and therefore a bug
    TRAPS: a-variance-cannot-fix-a-bias       how much mass leaves `RNA_SPLICED` because its `L` depends on an unsequenced intron
    TRAPS: variance-fitted-on-the-belief       the residual above the true ceiling, which is where a mate gap holding TWO introns shows
             up (`transcript_has_implicit_intron_in_gap` returns the FIRST match and stops)

⛔ **Nothing here is tuned.** Every target is read from the truth file or is a control that must not
move. If a residual will not close, it is measured and reported — never closed with a constant


Usage::

    python scripts/design/fl_anchor_gap.py [--index DIR] [--pilot DIR] [--truth/--no-truth]
    python scripts/design/fl_anchor_gap.py --json before.json     # record, then diff after a change
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

# The tail thresholds the audit and quote, so the two are comparable.
TAIL_REGION_BOUNDS = (500, 600, 700, 800)


def moments(counts: np.ndarray) -> tuple[float, float, int]:
    """Mean, sd and support ceiling of a raw 1-bp histogram."""
    c = np.asarray(counts, dtype=np.float64)
    total = c.sum()
    if total <= 0:
        return float("nan"), float("nan"), 0
    bins = np.arange(c.size, dtype=np.float64)
    mean = float((bins * c).sum() / total)
    var = float((bins * bins * c).sum() / total) - mean * mean
    nz = np.nonzero(c)[0]
    return mean, float(np.sqrt(max(var, 0.0))), int(nz[-1]) if nz.size else 0


def tail_mass(counts: np.ndarray, region_bound: int) -> float:
    """Fraction of the histogram's mass at ``length >= region_bound``."""
    c = np.asarray(counts, dtype=np.float64)
    total = c.sum()
    return float(c[region_bound:].sum() / total) if total > 0 else float("nan")


def read_truth(path: Path) -> dict[str, np.ndarray]:
    """``truth_fragment_lengths.tsv`` as ``{kind: histogram}``, indexed by length in bp.

    ⚠ Read from the simulator's own record of what it WROTE, not from anything the scanner produced.
    That is the whole point: the ceiling below is the library's, and it is a fact about the data rather
    than a threshold anybody picked.
    """
    rows: dict[str, dict[int, float]] = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        kind_i, len_i, count_i = header.index("kind"), header.index("fragment_length"), header.index("count")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= count_i:
                continue
            rows.setdefault(parts[kind_i], {})[int(parts[len_i])] = float(parts[count_i])
    out = {}
    for kind, hist in rows.items():
        size = max(hist) + 1
        arr = np.zeros(size, dtype=np.float64)
        for length, count in hist.items():
            arr[length] = count
        out[kind] = arr
    return out


def _pct(actual: float, reference: float) -> float:
    return 100.0 * (actual - reference) / reference if reference else float("nan")


def _drain(cache, index, region_types, junctions, *, seed: int):
    """The drained payload for one cached scan — production's own route, via the pipeline's helper.

    ⭐ Calls `pipeline._drain_side_buffer` rather than repeating its three steps, so this measurement
    cannot drift from what `rigel quant` actually does.
    """
    from rigel.pipeline import _drain_side_buffer

    return _drain_side_buffer(cache.payload, index, cache.strand_model, seed=seed)


def measure(
    cond_name: str, payload, truth: dict[str, np.ndarray] | None, crossing=None, gdna_opp=None
) -> dict:
    """Every number this script reports for one condition, as a plain dict (so it can be JSON-diffed).

    ⚠ Takes the PAYLOAD, not the cache, so the same measurement runs on pass one's tally and on the
    drained one — which is what makes the P4 panel a before/after with one thing varied.

    ⭐ ``crossing`` is the junction-opportunity de-tilt, and passing it is what makes the RNA panel
    report the pool the tool FITS FROM rather than the raw histogram. ``None`` reports the raw pool,
    which is the before-C3 arm.
    """
    from rigel.calibration.fl import build_fl_models, rna_fl_mass
    from rigel.calibration.junction_opportunity import detilt_pool

    fl = build_fl_models(payload, junction_opportunity=crossing, gdna_opportunity=gdna_opp)
    anchor, rna = np.asarray(fl.global_counts), rna_fl_mass(payload)
    if crossing is not None:
        rna = detilt_pool(rna, crossing)

    a_mean, a_sd, a_top = moments(anchor)
    r_mean, r_sd, r_top = moments(rna)
    row = {
        "condition": cond_name,
        "anchor": {"mean": a_mean, "sd": a_sd, "ceiling": a_top, "n": int(np.sum(anchor))},
        "rna_pool": {"mean": r_mean, "sd": r_sd, "ceiling": r_top, "n": int(np.sum(rna))},
        "gap_vs_anchor": {"mean_pct": _pct(r_mean, a_mean), "sd_pct": _pct(r_sd, a_sd)},
        "anchor_tail": {str(region_bound): tail_mass(anchor, region_bound) for region_bound in TAIL_REGION_BOUNDS},
        "qc": {
            "deposited": payload.qc.deposited,
            "dropped_too_long": payload.qc.dropped_too_long,
            "deferred_undetermined_gap": payload.qc.deferred_undetermined_gap,
            "introns_absorbed": payload.qc.introns_absorbed,
            # ⭐ The umbrella census, so the deferred total can be read as its three subclasses: is the
            # open question RNA-or-gDNA, which-structure, or both? `sj_implicit_fragments` used to sit
            # here and is GONE — TRAPS: a-variance-cannot-fix-a-bias is deleted, and "this L was partly inferred" no longer selects
            # anything: a fragment deposits when exactly ONE hypothesis survives, however it got there.
            "gap_resolved_spliced": payload.gap_resolution.gap_resolved_spliced,
            "gap_deferred_rna_or_gdna": payload.gap_resolution.gap_deferred_rna_or_gdna,
            "gap_deferred_which_introns": payload.gap_resolution.gap_deferred_which_introns,
            "gap_deferred_both": payload.gap_resolution.gap_deferred_both,
        },
    }
    row["qc"]["dropped_too_long_frac"] = (
        payload.qc.dropped_too_long / (payload.qc.deposited + payload.qc.dropped_too_long)
        if (payload.qc.deposited + payload.qc.dropped_too_long)
        else float("nan")
    )

    # ⭐ The gDNA panel reports the histogram the model is FITTED FROM — `fl.gdna_counts` — not a raw
    # pool sum. With `gdna_opp` that is the four pools each divided by its own opportunity; without it,
    # the contained pair. ⛔ Reading `gdna_fl_mass` directly here would report the four pools POOLED RAW,
    # which is a quantity nothing in the tool uses.
    gdna_pool = np.asarray(fl.gdna_counts, dtype=np.float64)
    g_mean, g_sd, g_top = moments(gdna_pool)
    row["gdna_pool"] = {
        "mean": g_mean,
        "sd": g_sd,
        "ceiling": g_top,
        "n": int(np.sum(gdna_pool)),
        "tail": {str(region_bound): tail_mass(gdna_pool, region_bound) for region_bound in TAIL_REGION_BOUNDS},
    }

    if truth is None:
        return row

    # ⭐ Which truth series is the anchor's population? On a zero-gDNA condition every deposited
    # fragment is RNA, so it is `mrna`; on gdna100 the pools are scored against their own series and the
    # anchor is a mixture, which is why the tail gate is quoted on the zero-gDNA conditions.
    for name, key, series in (
        ("anchor_vs_mrna", "mrna", anchor),
        ("rna_pool_vs_mrna", "mrna", rna),
        ("gdna_pool_vs_gdna", "gdna", gdna_pool),
    ):
        if key not in truth or truth[key].sum() <= 0:
            continue
        t_mean, t_sd, t_top = moments(truth[key])
        s_mean, s_sd, s_top = moments(series)
        row[name] = {
            "truth_mean": t_mean,
            "truth_sd": t_sd,
            "truth_ceiling": t_top,
            "mean_pct": _pct(s_mean, t_mean),
            "sd_pct": _pct(s_sd, t_sd),
            "ceiling": s_top,
            # ⭐ G-tail: mass ABOVE the library's true ceiling. Truth is 0 there by construction, so
            # this number is not a comparison — it is an absolute count of impossible molecules.
            "mass_above_truth_ceiling": tail_mass(series, t_top + 1),
            "tail": {str(region_bound): tail_mass(series, region_bound) for region_bound in TAIL_REGION_BOUNDS},
            "truth_tail": {str(region_bound): tail_mass(truth[key], region_bound) for region_bound in TAIL_REGION_BOUNDS},
        }
    return row


def _fmt(value: float, width: int = 9, places: int = 5) -> str:
    return f"{value:{width}.{places}f}" if value == value else " " * (width - 3) + "nan"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT, help="the scan-cache directory")
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument(
        "--suite",
        type=Path,
        default=None,
        help="where the conditions' truth files live (default: the scan cache's parent)",
    )
    ap.add_argument("--json", type=Path, default=None, help="also write every number to this file")
    ap.add_argument("--no-truth", action="store_true", help="skip the truth panel")
    ap.add_argument(
        "--drain",
        action="store_true",
        help="⭐ P4's gate: measure each condition BEFORE and AFTER the second pass drains the side "
        "buffer, so the tail and the gDNA control are read as a delta",
    )
    ap.add_argument("--seed", type=int, default=0, help="the drain's multinomial seed")
    ap.add_argument(
        "--raw-pool",
        action="store_true",
        help="⚠ report the RAW RNA_SPLICED histogram instead of the de-tilted pool the tool fits "
        "from — the before-C3 arm, and the only reason to want it is an A/B",
    )
    args = ap.parse_args()

    if not args.pilot.is_dir():
        print(f"no pilot scan-cache dir at {args.pilot}", file=sys.stderr)
        return 2
    suite = args.suite or args.pilot.parent

    from rigel.index import TranscriptIndex
    from rigel.scan_cache import read_scan_cache

    # ⚠ A cache is REFUSED unless it describes the index it is loaded against (graph_hash, reach_digest
    # and payload_schema_digest). That refusal is the point of the keys, so it is reported, not caught.
    index = TranscriptIndex.load(str(args.index))
    crossing = None
    if not args.raw_pool:
        from rigel.calibration.junction_opportunity import crossing_probability_from_index

        # ⚠ Built once, off the ANNOTATION — it does not depend on the condition or on the drain, so
        # the same divisor scores every arm and cannot itself be the thing that moved.
        crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = None
    if not args.raw_pool:
        from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index

        # ⭐ The four gDNA pools' own opportunities, likewise annotation-derived and condition-independent.
        # ⚠ `--raw-pool` drops BOTH divisors, so it is the before arm for the gDNA model as well as the RNA
        # one — and with both dropped the gDNA pool falls back to the CONTAINED pair, not the raw four.
        gdna_opp = gdna_opportunity_from_index(index, 4096)
    region_types = junctions = None
    if args.drain:
        from rigel.calibration.splice_graph import (
            build_junction_edge_arrays,
            build_region_partition_arrays,
        )

        _c, _o, region_types = build_region_partition_arrays(index)
        junctions = build_junction_edge_arrays(index)

    conditions = sorted(p for p in args.pilot.iterdir() if p.is_dir())
    rows = []

    print(f"{'condition':<44} {'anchor mean':>11} {'RNA mean':>9} {'d%':>7} "
          f"{'anchor sd':>10} {'RNA sd':>8} {'d%':>7} {'ceilings':>12}")
    print("-" * 116)
    for cond in conditions:
        try:
            cache = read_scan_cache(cond, index)
        except Exception as exc:  # noqa: BLE001 — a refused cache is a result, not a crash
            print(f"{cond.name:<44} REFUSED: {exc}")
            continue
        truth_path = suite / cond.name / "truth_fragment_lengths.tsv"
        truth = None if args.no_truth or not truth_path.is_file() else read_truth(truth_path)
        row = measure(cond.name, cache.payload, truth, crossing, gdna_opp)
        if args.drain:
            row["drained"] = measure(
                cond.name,
                _drain(cache, index, region_types, junctions, seed=args.seed),
                truth,
                crossing,
                gdna_opp,
            )
        rows.append(row)
        print(f"{cond.name:<44} {row['anchor']['mean']:11.1f} {row['rna_pool']['mean']:9.1f} "
              f"{row['gap_vs_anchor']['mean_pct']:+6.1f}% {row['anchor']['sd']:10.1f} "
              f"{row['rna_pool']['sd']:8.1f} {row['gap_vs_anchor']['sd_pct']:+6.1f}% "
              f"{row['anchor']['ceiling']:5d} vs {row['rna_pool']['ceiling']:<5d}")

    print()
    print("⭐ The ZERO-gDNA rows (`gdna_none_*`) are the falsification: every fragment there is RNA,")
    print("   so anchor and RNA pool describe ONE population and the gap is bias, not composition.")
    print("   §2 (scanner anchor): +11.6 % / +71.1 %.   TRAPS: a-purity-filter-is-a-length-filter (accumulator anchor): +7.7 % / +32.0 %.")

    scored = [r for r in rows if "anchor_vs_mrna" in r]
    if scored:
        print()
        print("═══ G-tail · the anchor against the library's TRUE support "
              " ═══")
        print(f"{'condition':<44} {'true ceil':>9} {'anchor ceil':>11} {'>ceiling':>9} "
              f"{'>=700':>9} {'too_long':>10} {'frac':>8}")
        print("-" * 116)
        for r in scored:
            a = r["anchor_vs_mrna"]
            print(f"{r['condition']:<44} {a['truth_ceiling']:9d} {a['ceiling']:11d} "
                  f"{_fmt(a['mass_above_truth_ceiling'])} {_fmt(a['tail']['700'])} "
                  f"{r['qc']['dropped_too_long']:10d} {_fmt(r['qc']['dropped_too_long_frac'], 8, 5)}")
        print("   ⭐ `>ceiling` and `>=700` must reach 0 on the zero-gDNA rows: the truth files say the")
        print("      library contains no such molecule. `too_long` is mass thrown away above 1000 bp.")

        print()
        print("═══ G-sd · mean and sd against truth, and G-gdna · ⛔ THE CONTROL THAT MUST NOT MOVE ═══")
        print(f"{'condition':<44} {'anchor d.mean':>13} {'anchor d.sd':>12} "
              f"{'gDNA d.mean':>12} {'gDNA d.sd':>10} {'gDNA>=600':>10} {'truth>=600':>11}")
        print("-" * 122)
        for r in scored:
            a = r["anchor_vs_mrna"]
            g = r.get("gdna_pool_vs_gdna")
            if g is None:
                print(f"{r['condition']:<44} {a['mean_pct']:+12.1f}% {a['sd_pct']:+11.1f}% "
                      f"{'—':>12} {'—':>10} {'—':>10} {'—':>11}")
                continue
            print(f"{r['condition']:<44} {a['mean_pct']:+12.1f}% {a['sd_pct']:+11.1f}% "
                  f"{g['mean_pct']:+11.1f}% {g['sd_pct']:+9.1f}% "
                  f"{_fmt(g['tail']['600'], 10)} {_fmt(g['truth_tail']['600'], 11)}")
        print("   ⛔ The two right-hand columns are G-gdna. gDNA has NO INTRONS TO MISS, so they were")
        print("      already equal to five decimals. A change there means the fix reached fragments")
        print("      with no introns — impossible, therefore a bug.")

        print()
        print("═══ the side buffer · how each gap resolved · TRAPS: variance-fitted-on-the-belief · the residual above the ceiling ═══")
        print(f"{'condition':<44} {'deposited':>10} {'resolved':>9} {'defer':>8} {'frac':>7} "
              f"{'rna|dna':>8} {'which':>7} {'both':>6} {'RNA pool n':>11} {'TRAPS: variance-fitted-on-the-belief resid':>9}")
        print("-" * 130)
        for r in scored:
            qc = r["qc"]
            offered = qc["deposited"] + qc["deferred_undetermined_gap"]
            frac = qc["deferred_undetermined_gap"] / offered if offered else float("nan")
            print(f"{r['condition']:<44} {qc['deposited']:10d} {qc['gap_resolved_spliced']:9d} "
                  f"{qc['deferred_undetermined_gap']:8d} {_fmt(frac, 7, 4)} "
                  f"{qc['gap_deferred_rna_or_gdna']:8d} {qc['gap_deferred_which_introns']:7d} "
                  f"{qc['gap_deferred_both']:6d} "
                  f"{r['rna_pool']['n']:11d} {_fmt(r['anchor_vs_mrna']['mass_above_truth_ceiling'])}")
        print("   ⛔ `defer` is HELD, not dropped: those fragments are in the side buffer, whole, and the")
        print("      identity is deposited + deferred + dropped_* == offered. ⚠ Between S1 and S3 the")
        print("      tally is deliberately THINNER — do NOT read the anchor's accuracy as a regression.")
        print("   ⚠ The deferred population is the LONG one: a longer gap admits more hypotheses. So the")
        print("      surviving anchor is biased SHORT and every junction-opportunity number in the docs")
        print("      must be re-measured AFTER the drain, not before.")
        print("   ⚠ `TRAPS: variance-fitted-on-the-belief resid` is the mass still above the true ceiling. Measure it; do NOT close it with")
        print("      a constant.")

    drained_rows = [r for r in rows if "drained" in r and "anchor_vs_mrna" in r]
    if drained_rows:
        print()
        print("═══ ⭐ P4 · THE TAIL · the anchor against the library's TRUE ceiling, BEFORE vs AFTER "
              "the drain ═══")
        print(f"{'condition':<44} {'true ceil':>9} {'ceil before':>11} {'ceil after':>10} "
              f"{'>ceil before':>12} {'>ceil after':>11} {'held':>7}")
        print("-" * 122)
        for r in drained_rows:
            b, a = r["anchor_vs_mrna"], r["drained"]["anchor_vs_mrna"]
            print(f"{r['condition']:<44} {b['truth_ceiling']:>9d} {b['ceiling']:>11d} "
                  f"{a['ceiling']:>10d} {_fmt(b['mass_above_truth_ceiling'], 12)} "
                  f"{_fmt(a['mass_above_truth_ceiling'], 11)} "
                  f"{r['drained']['qc']['deferred_undetermined_gap'] or r['qc']['deferred_undetermined_gap']:>7d}")
        print("   ⭐ THE GATE: no fragment above the library's true longest molecule. The truth files say")
        print("      the library contains none, so `>ceil after` is an absolute count of impossible")
        print("      molecules — not a comparison, and nothing about it is tunable.")
        print("   ⚠ The drain returns the LONG fragments to the tally, so `ceil after` rising is EXPECTED")
        print("      and is the point; what must not appear is mass above the TRUE ceiling.")

        print()
        print("═══ ⛔ P4 · THE gDNA CONTROL, before vs after ═══")
        print(f"{'condition':<44} {'gdna n before':>13} {'gdna n after':>12} {'d.mean before':>13} "
              f"{'d.mean after':>12} {'d.sd after':>10}")
        print("-" * 118)
        for r in drained_rows:
            gb, ga = r.get("gdna_pool_vs_gdna"), r["drained"].get("gdna_pool_vs_gdna")
            if gb is None or ga is None:
                continue
            print(f"{r['condition']:<44} {r['gdna_pool']['n']:>13d} "
                  f"{r['drained']['gdna_pool']['n']:>12d} {gb['mean_pct']:>+12.2f}% "
                  f"{ga['mean_pct']:>+11.2f}% {ga['sd_pct']:>+9.2f}%")
        print("   ⛔ On a ZERO-gDNA condition every fragment is RNA, so ANY growth in the gDNA pool is")
        print("      RNA contaminating it — a drained fragment that chose ∅ and landed contained in an")
        print("      intronic region. ⚠ On gdna100 growth is LEGITIMATE: a real gDNA fragment whose mate")
        print("      gap spans an annotated intron is genuinely ambiguous and was genuinely held.")
        print("   ⚠ So this is not TRAPS: pure-and-length-censored.6's control. There the fix could not reach a fragment with no")
        print("      introns, so any movement was a bug; here the drain reaches fragments it is supposed")
        print("      to. Read the ZERO-gDNA rows as the false-positive rate of the ∅ choice.")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
