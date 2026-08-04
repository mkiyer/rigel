"""⭐ **THE SECOND PASS, SCORED PER FRAGMENT against the simulator's own truth.**

    `docs/DESIGN.md` §4

Pass one holds every fragment whose unsequenced mate gap has more than one surviving explanation. Pass
two scores the candidates from pass-one evidence alone, draws one by multinomial sample, and re-deposits.
This asks the only question that matters about that: **did it pick the right one?**

⭐ **The oracle BAM's read name IS the per-fragment truth** — `ENST00000458178.2:1281-1331:f:281` names
the source template and the molecule's own coordinates, so a held fragment's true `L` is known exactly
and none of the numbers below is a target anyone chose. Three separate defects were found by scoring
against it, including the annihilation bug (P4.1), where an all-zero score factor annihilated the other
two and collapsed the record to a coin toss.

⛔ **THIS EXISTED ONLY AS AN AD-HOC MEASUREMENT UNTIL 2026-08-03**, which meant the second pass's
headline number — "90.5 % exact, +0.12 bp" — could not be re-derived when the panel was regenerated.
That is what this file is for.

**How a held record is linked to its truth.** The side buffer stores the fragment's clipped genomic
extent, not its name. So the oracle BAM is indexed by `(ref, leftmost, rightmost)` over each fragment's
two mates, and only keys whose fragments all agree on one true `L` are used — hence "with an
unambiguous true length". ⚠ Ambiguous keys are **counted and reported**, never silently dropped: if they
were a large share, the population being scored would not be the population being held.

**What the assigned `L` is.** `L` = genomic span minus cut introns (the one definition, C0-C2). For a
held record that is `end - start` less every observed intron (cut under every hypothesis) and every
intron implied by the CHOSEN hypothesis.

    python scripts/design/second_pass_accuracy.py [--index DIR] [--pilot DIR] [--suite DIR]
        [--conditions gdna_none_ss_0.99_nrna_none_capture_off] [--seed 0]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402
import pysam  # noqa: E402

from rigel.sim.read_name import parse_origin  # noqa: E402

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: ⭐ The zero-gDNA stranded condition is the falsification condition: every fragment is RNA, so a held
#: fragment's true `L` is its transcript-space length with no gDNA population to confuse the scoring.
DEFAULT_CONDITION = "gdna_none_ss_0.99_nrna_none_capture_off"


def truth_by_extent(bam_path: Path) -> tuple[dict[tuple[int, int, int], int], int]:
    """``(ref_id, leftmost, rightmost) -> true L`` from the oracle BAM, plus the ambiguous-key count.

    ⚠ The BAM's own reference ids are used, and the caller must translate the payload's ref index
    through ``index.ref_names``. The scanner seeds its ref-id space in ``ref_names`` order precisely so
    the two spaces agree (`index.py`'s ``set_ref_names``), but "precisely so" is not "therefore" — the
    translation is explicit here.
    """
    spans: dict[str, list] = {}
    truths: dict[tuple[int, int, int], set[int]] = defaultdict(set)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.reference_end is None:
                continue
            key = read.query_name
            entry = spans.get(key)
            if entry is None:
                origin = parse_origin(key)
                spans[key] = [
                    read.reference_id,
                    read.reference_start,
                    read.reference_end,
                    int(origin.end) - int(origin.start),
                    1,
                ]
                continue
            entry[1] = min(entry[1], read.reference_start)
            entry[2] = max(entry[2], read.reference_end)
            entry[4] += 1
            if entry[4] == 2:  # both mates seen; retire the fragment
                truths[(entry[0], entry[1], entry[2])].add(entry[3])
                spans[key] = entry
    unique = {key: next(iter(values)) for key, values in truths.items() if len(values) == 1}
    ambiguous = len(truths) - len(unique)
    return unique, ambiguous


def assigned_lengths(deferred, choices: np.ndarray) -> np.ndarray:
    """The `L` the drain assigned each held fragment: span minus observed minus chosen-hypothesis introns."""
    n = deferred.n_fragments
    span = deferred.end - deferred.start

    def intron_total(offsets: np.ndarray, flat: np.ndarray, index: np.ndarray) -> np.ndarray:
        lo, hi = offsets[index], offsets[index + 1]
        total = np.zeros(len(index), dtype=np.int64)
        for i, (a, b) in enumerate(zip(lo, hi)):
            if b > a:
                pairs = flat[2 * a : 2 * b].reshape(-1, 2)
                total[i] = int((pairs[:, 1] - pairs[:, 0]).sum())
        return total

    observed = intron_total(
        deferred.observed_intron_offsets, deferred.observed_introns, np.arange(n)
    )
    global_hypothesis = deferred.hypothesis_offsets[:-1] + choices
    implied = intron_total(
        deferred.hypothesis_intron_offsets, deferred.hypothesis_introns, global_hypothesis
    )
    return span - observed - implied


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None, help="where the oracle BAMs live")
    ap.add_argument("--conditions", nargs="*", default=[DEFAULT_CONDITION])
    ap.add_argument("--seed", type=int, default=0, help="the drain's multinomial seed")
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    suite = args.suite or args.pilot.parent
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.junction_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_junction_edge_arrays,
        build_node_partition_arrays,
    )
    from rigel.index import TranscriptIndex
    from rigel.scan_cache import read_scan_cache
    from rigel.second_pass import choose_hypotheses, score_held_fragments

    index = TranscriptIndex.load(str(args.index))
    _cuts, _offsets, node_types = build_node_partition_arrays(index)
    junctions = build_junction_edge_arrays(index)
    crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = gdna_opportunity_from_index(index, 4096)
    bam_ref_id_of_payload_ref: dict[int, int] = {}

    rows = []
    for name in args.conditions:
        cache = read_scan_cache(args.pilot / name, index)
        payload, deferred = cache.payload, cache.payload.deferred
        if deferred.n_fragments == 0:
            print(f"  {name}: nothing held; the side buffer is empty")
            continue

        # ⚠ Reproduce `pipeline._drain_side_buffer`'s call exactly — the same de-tilted RNA pool and the
        # same `rna_sense_frac`. A scorer fed anything else is not the scorer being measured.
        fl = build_fl_models(payload, junction_opportunity=crossing, gdna_opportunity=gdna_opp)
        scores = score_held_fragments(
            payload,
            fl_models=fl,
            rna_sense_frac=cache.strand_model.p_r1_sense,
            node_types=node_types,
            junctions=junctions,
        )
        choices = choose_hypotheses(scores, payload, seed=args.seed)
        assigned = assigned_lengths(deferred, choices)

        bam = suite / name / "sim_oracle.bam"
        if not bam.is_file():
            print(f"  {name}: no oracle BAM at {bam}", file=sys.stderr)
            continue
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            for payload_ref, ref_name in enumerate(index.ref_names):
                try:
                    bam_ref_id_of_payload_ref[payload_ref] = handle.get_tid(str(ref_name))
                except (KeyError, ValueError):
                    bam_ref_id_of_payload_ref[payload_ref] = -1
        truth, ambiguous = truth_by_extent(bam)

        # ⭐ Scored per candidate count as well as pooled: a record with 5 hypotheses is a harder draw
        # than one with 2, and a single percentage hides whether the score is working or the easy cases
        # simply dominate.
        runs = np.diff(deferred.hypothesis_offsets)
        errors: list[int] = []
        by_count: dict[int, list[int]] = defaultdict(list)
        unmatched = 0
        for i in range(deferred.n_fragments):
            key = (
                bam_ref_id_of_payload_ref.get(int(deferred.ref[i]), -1),
                int(deferred.start[i]),
                int(deferred.end[i]),
            )
            true_length = truth.get(key)
            if true_length is None:
                unmatched += 1
                continue
            error = int(assigned[i]) - int(true_length)
            errors.append(error)
            by_count[int(runs[i])].append(error)

        matched = len(errors)
        exact = sum(1 for error in errors if error == 0)
        error_array = np.asarray(errors, dtype=np.float64)
        row = {
            "condition": name,
            "held": int(deferred.n_fragments),
            "matched": matched,
            "unmatched": unmatched,
            "ambiguous_truth_keys": ambiguous,
            "exact": exact,
            "exact_fraction": exact / matched if matched else None,
            "mean_error_bp": float(error_array.mean()) if matched else None,
            "median_abs_error_bp": float(np.median(np.abs(error_array))) if matched else None,
            "by_candidate_count": {
                count: {"n": len(values), "exact_fraction": float(np.mean(np.asarray(values) == 0))}
                for count, values in sorted(by_count.items())
            },
        }
        rows.append(row)

        print(f"\n═══ ⭐ {name} ═══")
        print(f"  held fragments                        {row['held']:>12,}")
        print(f"  matched to an unambiguous true length {matched:>12,}  ({100 * matched / row['held']:.1f} %)")
        print(f"  unmatched (ambiguous or absent key)   {unmatched:>12,}")
        print(f"  ambiguous truth keys in the BAM       {ambiguous:>12,}")
        print(f"  ⭐ exactly the right length            {100 * row['exact_fraction']:>11.1f} %")
        print(f"  ⭐ mean error                          {row['mean_error_bp']:>+11.2f} bp")
        print(f"     median |error|                     {row['median_abs_error_bp']:>11.1f} bp")
        print("  by candidate count:")
        for count, stats in row["by_candidate_count"].items():
            print(f"     {count} candidates  n={stats['n']:>9,}  exact {100 * stats['exact_fraction']:.1f} %")

    if args.json and rows:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
