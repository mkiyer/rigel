#!/usr/bin/env python
"""Re-measure FRAGMENT_LENGTH_AUDIT.md §2's gap table against the CURRENT anchor, and score it on TRUTH.

⭐ **The falsification this area never had.** On a zero-gDNA condition every fragment is RNA, so the
unconditional anchor and the RNA pool describe **one population** and any gap between them is bias.
§2 measured that gap at **+11.6 % mean / +71.1 % sd** with the scanner's histogram as the anchor; C1
built the accumulator's own histogram and measured **+7.7 % / +32.0 %**, isolating the residual as the
junction-opportunity tilt (C3's target). This script says which of those the *shipped* code is on.

⚠ It reads the anchor **through ``build_fl_models``**, not off the payload directly — the question is
what the tool is wired to use, not what is available to it.

⭐ **The truth panel (``--truth``) is C2.6's gate**, `docs/SPEC_GAP_INTRONS.md` §4. The simulator writes
``truth_fragment_lengths.tsv`` beside every condition, so the library's realized fragment-length support
is known **exactly** and none of the targets below is chosen:

    G-tail   the anchor's mass above the library's TRUE ceiling must be 0, and `dropped_too_long`
             must collapse — a fragment longer than any molecule in the library is an uncut intron
    G-sd     the anchor's sd against the truth's sd
    G-gdna   ⛔ THE CONTROL. `DNA_INTERGENIC` is pure gDNA and gDNA has NO INTRONS TO MISS, so it was
             already exact to five decimals. If it MOVES, the fix reached fragments with no introns —
             which is impossible, and therefore a bug
    D1       how much mass leaves `RNA_SPLICED` because its `L` depends on an unsequenced intron
    D3       the residual above the true ceiling, which is where a mate gap holding TWO introns shows
             up (`transcript_has_implicit_intron_in_gap` returns the FIRST match and stops)

⛔ **Nothing here is tuned.** Every target is read from the truth file or is a control that must not
move. If a residual will not close, it is measured and reported — never closed with a constant
(`CARRY_FORWARD.md` §3 trap 12).

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

#: The tail thresholds the audit and JUNCTION_OPPORTUNITY.md §4 quote, so the two are comparable.
TAIL_CUTS = (500, 600, 700, 800)


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


def tail_mass(counts: np.ndarray, cut: int) -> float:
    """Fraction of the histogram's mass at ``length >= cut``."""
    c = np.asarray(counts, dtype=np.float64)
    total = c.sum()
    return float(c[cut:].sum() / total) if total > 0 else float("nan")


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


def measure(cond_name: str, cache, truth: dict[str, np.ndarray] | None) -> dict:
    """Every number this script reports for one condition, as a plain dict (so it can be JSON-diffed)."""
    from rigel.calibration.fl import build_fl_models, gdna_fl_mass, rna_fl_mass

    payload = cache.payload
    fl = build_fl_models(payload)
    anchor, rna = np.asarray(fl.global_counts), rna_fl_mass(payload)

    a_mean, a_sd, a_top = moments(anchor)
    r_mean, r_sd, r_top = moments(rna)
    row = {
        "condition": cond_name,
        "anchor": {"mean": a_mean, "sd": a_sd, "ceiling": a_top, "n": int(np.sum(anchor))},
        "rna_pool": {"mean": r_mean, "sd": r_sd, "ceiling": r_top, "n": int(np.sum(rna))},
        "gap_vs_anchor": {"mean_pct": _pct(r_mean, a_mean), "sd_pct": _pct(r_sd, a_sd)},
        "anchor_tail": {str(cut): tail_mass(anchor, cut) for cut in TAIL_CUTS},
        "qc": {
            "deposited": payload.qc.deposited,
            "dropped_too_long": payload.qc.dropped_too_long,
            "dropped_ambiguous_path": payload.qc.dropped_ambiguous_path,
            "sj_implicit_fragments": payload.qc.sj_implicit_fragments,
            "introns_absorbed": payload.qc.introns_absorbed,
        },
    }
    row["qc"]["dropped_too_long_frac"] = (
        payload.qc.dropped_too_long / (payload.qc.deposited + payload.qc.dropped_too_long)
        if (payload.qc.deposited + payload.qc.dropped_too_long)
        else float("nan")
    )

    # ⛔ THE CONTROL. `DNA_INTERGENIC` is contained-in-an-intergenic-node, i.e. pure gDNA — and gDNA has
    # no introns to miss. It read 0.00023 at >=600 bp against a truth of 0.00024 BEFORE this work; if it
    # moves, the change reached fragments that have no introns, which cannot happen.
    gdna_pool = gdna_fl_mass(payload)
    g_mean, g_sd, g_top = moments(gdna_pool)
    row["gdna_pool"] = {
        "mean": g_mean,
        "sd": g_sd,
        "ceiling": g_top,
        "n": int(np.sum(gdna_pool)),
        "tail": {str(cut): tail_mass(gdna_pool, cut) for cut in TAIL_CUTS},
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
            "tail": {str(cut): tail_mass(series, cut) for cut in TAIL_CUTS},
            "truth_tail": {str(cut): tail_mass(truth[key], cut) for cut in TAIL_CUTS},
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
        row = measure(cond.name, cache, truth)
        rows.append(row)
        print(f"{cond.name:<44} {row['anchor']['mean']:11.1f} {row['rna_pool']['mean']:9.1f} "
              f"{row['gap_vs_anchor']['mean_pct']:+6.1f}% {row['anchor']['sd']:10.1f} "
              f"{row['rna_pool']['sd']:8.1f} {row['gap_vs_anchor']['sd_pct']:+6.1f}% "
              f"{row['anchor']['ceiling']:5d} vs {row['rna_pool']['ceiling']:<5d}")

    print()
    print("⭐ The ZERO-gDNA rows (`gdna_none_*`) are the falsification: every fragment there is RNA,")
    print("   so anchor and RNA pool describe ONE population and the gap is bias, not composition.")
    print("   §2 (scanner anchor): +11.6 % / +71.1 %.   C1 (accumulator anchor): +7.7 % / +32.0 %.")

    scored = [r for r in rows if "anchor_vs_mrna" in r]
    if scored:
        print()
        print("═══ G-tail · the anchor against the library's TRUE support "
              "(docs/SPEC_GAP_INTRONS.md §4) ═══")
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
        print("═══ D1 · mass barred from the pure-RNA pool · D3 · the residual above the ceiling ═══")
        print(f"{'condition':<44} {'deposited':>10} {'sj_implicit':>12} {'frac':>8} "
              f"{'ambig_path':>11} {'RNA pool n':>11} {'D3 resid':>9}")
        print("-" * 116)
        for r in scored:
            qc = r["qc"]
            frac = qc["sj_implicit_fragments"] / qc["deposited"] if qc["deposited"] else float("nan")
            print(f"{r['condition']:<44} {qc['deposited']:10d} {qc['sj_implicit_fragments']:12d} "
                  f"{_fmt(frac, 8, 4)} {qc['dropped_ambiguous_path']:11d} "
                  f"{r['rna_pool']['n']:11d} {_fmt(r['anchor_vs_mrna']['mass_above_truth_ceiling'])}")
        print("   ⚠ `sj_implicit` now counts MIXED fragments too (D1), so both the pool's mean and its")
        print("      size move, and for two reasons at once — the intron is cut AND the fragment leaves.")
        print("   ⚠ `D3 resid` is the mass still above the true ceiling. A mate gap holding TWO annotated")
        print("      introns keeps only the first cut, and this is where that shows up. Measure it; do")
        print("      NOT close it with a constant.")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
