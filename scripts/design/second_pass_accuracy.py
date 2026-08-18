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
headline number — "90.5 % exact, +0.12 bp", measured on the since-deleted `pilot` panel — could not be
re-derived when the panel was regenerated. That is what this file is for, and it has now done it: on the
16-condition ladder's `gdna_g00_ss_0.99_nrna_none_capture_off`, 2026-08-17, **90.8 % exact** and a mean
error of **−0.03 bp** over 336,351 of 339,431 held fragments matched to an unambiguous true length.
⚠ Different panel, so this is a re-derivation of the QUANTITY and not a reproduction of the number.

**How a held record is linked to its truth.** The side buffer stores the fragment's clipped genomic
extent, not its name. So the oracle BAM is indexed by `(ref, leftmost, rightmost)` over each fragment's
two mates, and only keys whose fragments all agree on one true `L` are used — hence "with an
unambiguous true length". ⚠ Ambiguous keys are **counted and reported**, never silently dropped: if they
were a large share, the population being scored would not be the population being held.

**What the assigned `L` is.** `L` = genomic span minus EXCISED introns — the one definition. For a
held record that is `end - start` less every observed intron (excised under every hypothesis) and every
intron implied by the CHOSEN hypothesis.
⚠ That sentence said *"minus region_bound introns … (region_bound under every hypothesis)"* until
2026-08-17, and it was a BULK-RENAME CASUALTY: `5126cc89` mapped `cut -> REGION_BOUND` wholesale, and
here `cut` was ordinary English for *excised*, not the REGION_BOUND concept. That is
`TRAPS: two-masks-one-name` committed by a rename, and it is exactly what `rename_census.py` refuses to
do automatically. The identifiers it renamed in the same commit (`payload.region_bounds`,
`ref_region_bound_offsets`) are correct and untouched.

⛔ **ITS DEFAULT PANEL WAS DELETED AND THE FILE STILL POINTED AT IT.** `~/…/suite/pilot/scan_cache` went
with the 2026-08-13 rebuild, so a bare run died on `FileNotFoundError` at a `manifest.json` under a
directory that no longer exists — and `docs/SUCCESS.md` invokes this file bare. The default is now the
16-condition ladder's own scan cache, and the default condition is `gdna_g00_ss_0.99_…` (the ladder
spells the zero-gDNA rung `g00`; the pilot spelled it `none`). ⚠ **Four sibling instruments still
declare the same deleted path** — `held_flux_census.py`, `calibration_truth_ab.py`,
`gdna_pool_census.py`, `fl_anchor_gap.py`, each with its own `DEFAULT_PILOT` literal. One shared home is
the right final shape; five copies is how they came to rot together.

    python scripts/design/second_pass_accuracy.py [--index DIR] [--scan-cache DIR] [--suite DIR]
        [--conditions gdna_g00_ss_0.99_nrna_none_capture_off] [--seed 0]
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
#: ⛔ The 16-condition ladder, which is the ONLY panel on disk. The name is `scan_cache` and not
#: `pilot` because the flag pointing here must not be named after a panel that was deleted.
DEFAULT_SCAN_CACHE = _RUNS / "suite" / "ladder" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: ⭐ The zero-gDNA stranded condition is the falsification condition: every fragment is RNA, so a held
#: fragment's true `L` is its transcript-space length with no gDNA population to confuse the scoring.
#: ⚠ `g00`, not `none` — the rebuilt ladder spells its zero rung with the same `gNN` grammar as the rest.
DEFAULT_CONDITION = "gdna_g00_ss_0.99_nrna_none_capture_off"


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


def intron_total(offsets: np.ndarray, flat: np.ndarray, index: np.ndarray) -> np.ndarray:
    """Total intron length for each selected CSR row. ⚠ Module level so the per-hypothesis table below
    uses the SAME expression the assigned length does — two implementations of one quantity is how the
    two come to disagree."""
    lo, hi = offsets[index], offsets[index + 1]
    total = np.zeros(len(index), dtype=np.int64)
    for i, (a, b) in enumerate(zip(lo, hi)):
        if b > a:
            pairs = flat[2 * a : 2 * b].reshape(-1, 2)
            total[i] = int((pairs[:, 1] - pairs[:, 0]).sum())
    return total


def assigned_lengths(deferred, choices: np.ndarray) -> np.ndarray:
    """The `L` the drain assigned each held fragment: span minus observed minus chosen-hypothesis introns."""
    n = deferred.n_fragments
    span = deferred.end - deferred.start
    observed = intron_total(
        deferred.observed_intron_offsets, deferred.observed_introns, np.arange(n)
    )
    global_hypothesis = deferred.hypothesis_offsets[:-1] + choices
    implied = intron_total(
        deferred.hypothesis_intron_offsets, deferred.hypothesis_introns, global_hypothesis
    )
    return span - observed - implied


def _defect_table(deferred, choices, truth, ref_map, payload, index, sj, fl) -> None:
    """⭐⭐ **DOES THE SJ/BOUNDARY OPPORTUNITY MISMATCH ACTUALLY FLIP A CHOICE?**

    Per held record: which hypotheses reproduce the TRUE length, whether those are genomic or spliced,
    which one was drawn, and the BOTTLENECK exonic reach of the sj involved. ⭐ The bottleneck is
    the right summary because ``second_pass`` scores a spliced path with ``_bottleneck`` (a ``min``) over
    its sj, so the weakest sj is what the path is judged on.

    ⛔ **The control is the other direction.** If truth-spliced records are lost to the genomic
    hypothesis at short reach AND truth-genomic records are lost to the spliced one at the same rate,
    that is difficulty, not bias. Both columns are printed, always.

    ⚠ The reach bins are multiples of the library's OWN mean fragment length, because the mechanism is
    ``A_j(w)/(w-1)`` — a function of reach RELATIVE to ``w`` — and not of any absolute number of bases.
    """
    from rigel.calibration.splice_graph import build_sj_geometry_arrays
    from rigel.second_pass import _sj_id
    from rigel.types import Strand

    geom = build_sj_geometry_arrays(index)
    bottleneck_reach = np.minimum(
        np.asarray(geom.reach_lo, np.float64), np.asarray(geom.reach_hi, np.float64)
    )
    region_bounds = payload.region_bounds
    n_hyp = int(deferred.hypothesis_offsets[-1])
    implied_all = intron_total(
        deferred.hypothesis_intron_offsets, deferred.hypothesis_introns, np.arange(n_hyp)
    )
    observed_all = intron_total(
        deferred.observed_intron_offsets, deferred.observed_introns, np.arange(deferred.n_fragments)
    )
    span_all = deferred.end.astype(np.int64) - deferred.start.astype(np.int64)

    mu = float((np.arange(fl.global_pmf.shape[0]) * fl.global_pmf).sum())
    boundaries = (0.0, 0.5 * mu, 1.0 * mu, 2.0 * mu, 4.0 * mu, float("inf"))
    labels = ("<0.5 mu", "0.5-1 mu", "1-2 mu", "2-4 mu", ">=4 mu")
    #: [bin] -> [n_true_spliced, chose_genomic, n_true_genomic, chose_spliced]
    tally = np.zeros((len(labels), 4), np.int64)
    no_reach = 0

    for i in range(deferred.n_fragments):
        key = (ref_map.get(int(deferred.ref[i]), -1), int(deferred.start[i]), int(deferred.end[i]))
        true_length = truth.get(key)
        if true_length is None:
            continue
        h0, h1 = int(deferred.hypothesis_offsets[i]), int(deferred.hypothesis_offsets[i + 1])
        base = int(span_all[i]) - int(observed_all[i])
        kinds, lengths = [], []
        reaches: list[float] = []
        region_bound_lo = int(payload.ref_region_bound_offsets[int(deferred.ref[i])])
        region_bound_hi = int(payload.ref_region_bound_offsets[int(deferred.ref[i]) + 1])
        motif = int(deferred.sj_strand[i])
        for h in range(h0, h1):
            introns = [tuple(p) for p in deferred.hypothesis_introns_of(h).tolist()]
            kinds.append(bool(introns))  # True == spliced
            lengths.append(base - int(implied_all[h]))
            for a, b in introns:
                jid = _sj_id(
                    sj, region_bounds, region_bound_lo, region_bound_hi, a, b,
                    motif if motif != int(Strand.NONE) else int(deferred.hypothesis_sj_strand[h]),
                )
                if jid >= 0:
                    reaches.append(float(bottleneck_reach[jid]))
        # ⛔ A record only speaks to this defect if BOTH kinds are on the ballot; otherwise there is no
        # choice for the opportunity mismatch to tilt (`TRAPS: could-the-arm-have-fired`).
        if not (any(kinds) and not all(kinds)):
            continue
        if not reaches:
            no_reach += 1
            continue
        correct = [h - h0 for h in range(h0, h1) if lengths[h - h0] == int(true_length)]
        if not correct:
            continue
        chosen = int(choices[i])
        b = int(np.searchsorted(np.asarray(boundaries[1:], np.float64), min(reaches), side="right"))
        b = min(b, len(labels) - 1)
        true_spliced = all(kinds[c] for c in correct)
        true_genomic = all(not kinds[c] for c in correct)
        if true_spliced:
            tally[b, 0] += 1
            tally[b, 1] += int(not kinds[chosen])
        elif true_genomic:
            tally[b, 2] += 1
            tally[b, 3] += int(kinds[chosen])

    print("\n  ⭐⭐ SPLICED-vs-GENOMIC CONFUSION, by the sj's BOTTLENECK EXONIC REACH")
    print(f"     (bins are multiples of this library's own mean fragment length, mu = {mu:.1f} bp;")
    print("      only records with BOTH a genomic and a spliced candidate are counted)")
    print(f"     {'reach':<12}{'n true-SPLICED':>16}{'-> chose genomic':>18}"
          f"{'n true-GENOMIC':>16}{'-> chose spliced':>18}")
    print("     " + "-" * 80)
    for b, label in enumerate(labels):
        ns, cg, ng, cs = (int(v) for v in tally[b])
        if ns == 0 and ng == 0:
            continue
        print(f"     {label:<12}{ns:>16,}{100 * cg / ns if ns else 0:>17.1f}%"
              f"{ng:>16,}{100 * cs / ng if ng else 0:>17.1f}%")
    ts, tg = int(tally[:, 0].sum()), int(tally[:, 2].sum())
    print("     " + "-" * 80)
    print(f"     {'ALL':<12}{ts:>16,}{100 * int(tally[:, 1].sum()) / ts if ts else 0:>17.1f}%"
          f"{tg:>16,}{100 * int(tally[:, 3].sum()) / tg if tg else 0:>17.1f}%")
    print(f"     records with no resolvable sj: {no_reach:,}")
    print("     ⛔ THE DEFECT PREDICTS the 'chose genomic' column RISES as reach falls, while the")
    print("       'chose spliced' control does NOT. Both rising together is difficulty, not bias.")
    if ts == 0 and tg == 0:
        print("     ⛔ EMPTY — no record had both kinds on the ballot, so this table tests NOTHING.")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--scan-cache", type=Path, default=DEFAULT_SCAN_CACHE,
                    help="a panel's scan_cache/ directory (default: the 16-condition ladder's)")
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None, help="where the oracle BAMs live")
    ap.add_argument("--conditions", nargs="*", default=[DEFAULT_CONDITION])
    ap.add_argument("--seed", type=int, default=0, help="the drain's multinomial seed")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument(
        "--defect",
        action="store_true",
        help="⭐⭐ the SPLICED-vs-GENOMIC confusion, stratified by the sj's EXONIC REACH. "
             "The genomic hypothesis reads `boundary_unspliced_inv_length_sum`, whose `1/(w-1)` deposit "
             "cancels its `w-1` opportunity EXACTLY, so E = rho. The spliced hypothesis reads "
             "`sj_inv_length_sum`, same quantum but a TAPERED opportunity, so "
             "E = rho * E[A_j(w)/(w-1)] < rho — the shortfall growing as the flanking exons shorten. "
             "If that biases the choice, records whose truth is SPLICED are lost to the GENOMIC "
             "hypothesis at short reach and the effect must vanish at long reach.",
    )
    args = ap.parse_args()

    # ⛔ FAIL FAST AND NAME WHAT IS MISSING. A deleted panel used to surface here as a bare
    # `FileNotFoundError` on a `manifest.json` several directories below the thing that was actually
    # absent, which reads as a corrupt cache rather than an absent panel.
    missing = [n for n in args.conditions if not (args.scan_cache / n / "manifest.json").is_file()]
    if not args.scan_cache.is_dir() or missing:
        # ⚠ The absent thing is named off `is_dir()`, never off an EMPTY listing: a scan-cache
        # directory that exists and holds nothing is a different fault from one that was never built,
        # and keying the message on `have` reported the second when it was the first.
        have = sorted(p.name for p in args.scan_cache.iterdir()) if args.scan_cache.is_dir() else []
        raise SystemExit(
            f"⛔ no scan cache for "
            f"{'the directory itself' if not args.scan_cache.is_dir() else missing} under "
            f"{args.scan_cache}\n"
            + (f"   present there: {', '.join(have)}\n" if have else
               "   the directory itself holds no condition.\n" if args.scan_cache.is_dir() else "")
            + "   Build it with:  python scripts/sim/panel.py cache --config "
              "scripts/sim/configs/gdna_ladder.yaml\n"
              "   (`panel.py status` first — every stage is expensive and resumable, and it names the "
              "next one.)"
        )

    suite = args.suite or args.scan_cache.parent
    from rigel.calibration.fl import build_fl_models
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
    from rigel.calibration.sj_opportunity import crossing_probability_from_index
    from rigel.calibration.splice_graph import (
        build_sj_arrays,
        build_region_partition_arrays,
    )
    from rigel.index import TranscriptIndex
    from rigel.scan_cache import read_scan_cache
    from rigel.second_pass import choose_hypotheses, score_held_fragments

    index = TranscriptIndex.load(str(args.index))
    _region_bounds, _offsets, region_types = build_region_partition_arrays(index)
    sj = build_sj_arrays(index)
    crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = gdna_opportunity_from_index(index, 4096)
    bam_ref_id_of_payload_ref: dict[int, int] = {}

    rows = []
    for name in args.conditions:
        cache = read_scan_cache(args.scan_cache / name, index)
        payload, deferred = cache.payload, cache.payload.deferred
        if deferred.n_fragments == 0:
            print(f"  {name}: nothing held; the side buffer is empty")
            continue

        # ⚠ Reproduce `pipeline._drain_side_buffer`'s call exactly — the same de-tilted RNA pool and the
        # same `rna_sense_frac`. A scorer fed anything else is not the scorer being measured.
        fl = build_fl_models(payload, sj_opportunity=crossing, gdna_opportunity=gdna_opp)
        scores = score_held_fragments(
            payload,
            fl_models=fl,
            rna_sense_frac=cache.strand_model.p_r1_sense,
            region_types=region_types,
            sj=sj,
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

        if args.defect:
            _defect_table(
                deferred, choices, truth, bam_ref_id_of_payload_ref, payload, index, sj, fl
            )

    if args.json and rows:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
