#!/usr/bin/env python
"""⭐ **D-3, MEASURED**: how often does a held hypothesis have ZERO flux evidence, and is the zero RIGHT?

    Decision: D-3 · Phase: its P2
    Scores through the shipped `rigel.second_pass.score_held_fragments` — the question is what the
    tool does, not what could be done.

D-3 asks whether a junction with zero flux is *impossible* or *merely unobserved*::

    ⛔ A hard zero makes the hypothesis unselectable and can empty the score vector.
    ⛔ A pseudocount is a magic number and will not be invented.
    ⚠  Measure first: the held fragments are excluded from the tally they are scored against, so a
       junction used ONLY by deferred fragments reads zero — how often that happens decides this.

⭐ **A COUNT OF ZEROS IS NOT AN ANSWER.** A zero is *correct* when the transcript is not expressed and
*wrong* when its only evidence is held. So every zero here is pushed into one of four causes, and the
annotated-but-empty ones are scored against the simulator's own `truth_abundances.tsv`, where
``observed_mrna_fragments`` says whether the molecule existed at all. Nothing is tuned: the only targets
are truth columns and a structural identity.

The four causes, and why the split matters
------------------------------------------

======================  ============================================================================
``unannotated``         the implied intron resolves to no junction slot (`jid < 0`). ⚠ Not a flux
                        question at all — a lookup miss, and a different decision from D-3's
``annotated_empty``     the slot exists and `sj_inv_length_sum == 0`. ⭐ **THIS IS D-3'S POPULATION**
``no_evidence_set``     ⛔ ∅ only, and it is an ARTEFACT: the scorer's contiguous-edge set is empty.
                        See :func:`_lines_inside_inclusive` — the deposit rule says which lines
                        distinguish ∅ from a spliced path, and the shipped scorer asks for a strict
                        subset of them, empty whenever the intron spans exactly one node
``zero_edge``           ∅'s evidence set is non-empty but some edge in it carries no unspliced mass
======================  ============================================================================

⛔ **The ∅ arm is reported TWICE**, under the shipped exclusive rule and under the rule derived from the
deposit, because the two differ and the difference is not small. That is one thing varied, and it is the
only way the ∅ zero rate can be read as data rather than as the artefact.

Usage::

    python scripts/design/held_flux_census.py [--index DIR] [--pilot DIR] [--json out.json]
    python scripts/design/held_flux_census.py --conditions gdna_none_ss_0.99_nrna_none_capture_off
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

_RUNS = Path.home() / "Downloads" / "rigel_runs"
DEFAULT_PILOT = _RUNS / "suite" / "pilot" / "scan_cache"
DEFAULT_INDEX = _RUNS / "suite" / "rigel_index"

#: The fixed-point scale `inv_length_quantum` deposits in. Decoded, never chosen.
INV_LENGTH_SCALE = float(1 << 32)


def _lines_strictly_inside(cuts: np.ndarray, lo: int, hi: int, start: int, end: int) -> tuple[int, int]:
    """⛔ **THE PRE-D-6 RULE, kept only so the measurement stays reproducible.**

    ``rigel.second_pass`` asked for the lines ``start < c < end`` until 2026-08-02. That drops the donor
    and acceptor lines — the two that are guaranteed to exist and guaranteed to discriminate — and
    returns an EMPTY range when the intron spans one node. The shipped rule is now
    :func:`rigel.second_pass._distinguishing_lines`, derived from the deposit; this function exists so
    the ``artefact`` column below can still be computed. ⛔ Do not import it into ``src/``.
    """
    first = int(np.searchsorted(cuts[lo:hi], start, side="right"))
    last = int(np.searchsorted(cuts[lo:hi], end, side="left"))
    return first, last


def _read_truth_abundance(path: Path) -> dict[str, float]:
    """``transcript_id -> observed_mrna_fragments`` from the simulator's own record of what it wrote."""
    out: dict[str, float] = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        t_i = header.index("transcript_id")
        obs_i = header.index("observed_mrna_fragments")
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > obs_i:
                out[parts[t_i]] = float(parts[obs_i])
    return out


def census(payload, scored, junctions, t_ids, truth: dict[str, float] | None) -> dict:
    """Decompose every zero in ``scored.terms.density`` into its cause. Mirrors the scorer, then checks it.

    ⚠ The mirror is **verified, not asserted**: the reconstructed zero mask is compared against the
    shipped scorer's own density array and the run aborts on any disagreement. A census that has drifted
    from the code it describes is worse than no census.
    """
    from rigel.second_pass import _distinguishing_lines, _junction_id
    from rigel.types import Strand

    deferred = payload.deferred
    cuts = payload.cut_positions
    density = scored.terms.density
    n_introns = np.diff(deferred.hypothesis_intron_offsets)
    is_spliced = n_introns > 0

    cause = np.zeros(deferred.n_hypotheses, dtype=np.int8)  # 0 = positive
    CAUSE = {1: "unannotated", 2: "annotated_empty", 3: "no_evidence_set", 4: "zero_edge"}
    #: ⭐ The ∅ arm under the PRE-D-6 rule — the artefact's size, isolated.
    empty_old_zero = np.zeros(deferred.n_hypotheses, dtype=bool)
    reconstructed = np.zeros(deferred.n_hypotheses, dtype=bool)

    #: annotated-but-empty junction slot -> how many held hypotheses claim it
    empty_slot_claims: Counter[int] = Counter()
    #: and the supporting transcripts of the hypotheses that claim one
    empty_slot_supporters: dict[int, set[int]] = {}

    sj_flux = np.asarray(payload.sj_inv_length_sum, dtype=np.float64).sum(axis=1)
    edge_flux = np.asarray(payload.edge_unspliced_inv_length_sum, dtype=np.float64).sum(axis=1)

    for i in range(deferred.n_fragments):
        ref = int(deferred.ref[i])
        observed_motif = int(deferred.sj_strand[i])
        h0, h1 = int(deferred.hypothesis_offsets[i]), int(deferred.hypothesis_offsets[i + 1])
        cut_lo = int(payload.ref_cut_offsets[ref])
        cut_hi = int(payload.ref_cut_offsets[ref + 1])
        edge_base = int(payload.ref_edge_offsets[ref])

        contested = [
            tuple(pair)
            for h in range(h0, h1)
            for pair in deferred.hypothesis_introns_of(h).tolist()
        ]

        for h in range(h0, h1):
            introns = [tuple(pair) for pair in deferred.hypothesis_introns_of(h).tolist()]
            implied = int(deferred.hypothesis_sj_strand[h])
            if introns:
                motif = observed_motif if observed_motif != int(Strand.NONE) else implied
                worst = 0
                for a, b in introns:
                    jid = _junction_id(junctions, cuts, cut_lo, cut_hi, a, b, motif)
                    if jid < 0:
                        worst = max(worst, 1)
                    elif sj_flux[jid] <= 0.0:
                        worst = max(worst, 2)
                        empty_slot_claims[jid] += 1
                        empty_slot_supporters.setdefault(jid, set()).update(
                            int(t) for t in deferred.supporting_t_of(h)
                        )
                cause[h] = worst
                reconstructed[h] = worst != 0
                empty_old_zero[h] = worst != 0
            else:
                for rule, sink in (
                    (_distinguishing_lines, "shipped"),
                    (_lines_strictly_inside, "old"),
                ):
                    values = []
                    for a, b in contested:
                        first, last = rule(cuts, cut_lo, cut_hi, a, b)
                        values.extend(
                            edge_flux[edge_base + line - 1] for line in range(first, last)
                        )
                    zero = (not values) or min(values) <= 0.0
                    if sink == "shipped":
                        cause[h] = 0 if not zero else (3 if not values else 4)
                        reconstructed[h] = zero
                    else:
                        empty_old_zero[h] = zero

    # ⛔ The self-check. A cause table that disagrees with the scorer is describing something else.
    actual = density <= 0.0
    if not np.array_equal(actual, reconstructed):
        n_bad = int((actual != reconstructed).sum())
        raise SystemExit(
            f"the census diverged from `score_held_fragments` on {n_bad} hypotheses. The decomposition "
            f"mirrors the scorer's own rho branch, so a divergence means one of them has moved."
        )

    def rate(mask, sub=None):
        sub = np.ones_like(mask) if sub is None else sub
        n = int(sub.sum())
        return {"n": n, "zero": int((mask & sub).sum()), "frac": float((mask & sub).sum() / n) if n else float("nan")}

    zero_f = scored.terms.length_likelihood <= 0.0
    offsets = deferred.hypothesis_offsets
    per_record_zero = np.add.reduceat(actual.astype(np.int64), offsets[:-1]) if deferred.n_fragments else np.zeros(0, np.int64)
    per_record_n = np.diff(offsets)
    #: ⭐ The case that matters more than the empty vector: the zero is DECISIVE — it eliminates a
    #: hypothesis while the record still has a live one, so the draw never sees it.
    decisive = (per_record_zero > 0) & (per_record_zero < per_record_n)
    all_zero_rho = per_record_zero == per_record_n

    out = {
        "n_fragments": int(deferred.n_fragments),
        "n_hypotheses": int(deferred.n_hypotheses),
        "rho_zero": {
            "all": rate(actual),
            "spliced": rate(actual, is_spliced),
            "empty_shipped": rate(actual, ~is_spliced),
            "empty_pre_d6": rate(empty_old_zero, ~is_spliced),
        },
        "cause": {
            name: int(((cause == code) & is_spliced).sum() + ((cause == code) & ~is_spliced).sum())
            for code, name in CAUSE.items()
        },
        "cause_spliced": {name: int(((cause == code) & is_spliced).sum()) for code, name in CAUSE.items()},
        "cause_empty": {name: int(((cause == code) & ~is_spliced).sum()) for code, name in CAUSE.items()},
        "f_zero": {
            "all": rate(zero_f),
            "spliced": rate(zero_f, is_spliced),
            "empty": rate(zero_f, ~is_spliced),
        },
        "records": {
            "n": int(deferred.n_fragments),
            "undecided": int(scored.n_undecided),
            "undecided_frac": float(scored.n_undecided / max(deferred.n_fragments, 1)),
            "all_hypotheses_zero_rho": int(all_zero_rho.sum()),
            "decisive_zero": int(decisive.sum()),
            "decisive_zero_frac": float(decisive.sum() / max(deferred.n_fragments, 1)),
        },
        "empty_slots": {
            "n_distinct": len(empty_slot_claims),
            "n_claims": int(sum(empty_slot_claims.values())),
            "claims_per_slot_max": max(empty_slot_claims.values()) if empty_slot_claims else 0,
        },
    }

    # ⭐ THE TRUTH ARM. An annotated junction with no flux is CORRECT when nothing expressed it and WRONG
    # when its supporters did express — the second is D-3's self-exclusion, and it is the only version
    # of the zero that costs an answer.
    if truth is not None and empty_slot_claims:
        expressed_slots = 0
        expressed_claims = 0
        silent_claims = 0
        for jid, claims in empty_slot_claims.items():
            supporters = empty_slot_supporters.get(jid, set())
            observed = max(
                (truth.get(t_ids[t], 0.0) for t in supporters if 0 <= t < len(t_ids)),
                default=0.0,
            )
            if observed > 0:
                expressed_slots += 1
                expressed_claims += claims
            else:
                silent_claims += claims
        out["empty_slots"]["truth"] = {
            "slots_with_an_expressed_supporter": expressed_slots,
            "slots_all_supporters_silent": len(empty_slot_claims) - expressed_slots,
            "claims_on_expressed": expressed_claims,
            "claims_on_silent": silent_claims,
            "frac_claims_on_expressed": float(
                expressed_claims / max(expressed_claims + silent_claims, 1)
            ),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--suite", type=Path, default=None, help="where the truth files live")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--json", type=Path, default=None)
    args = ap.parse_args()

    if not args.pilot.is_dir():
        print(f"no pilot scan-cache dir at {args.pilot}", file=sys.stderr)
        return 2
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
    from rigel.second_pass import score_held_fragments

    index = TranscriptIndex.load(str(args.index))
    # ⚠ The scorer reads the DE-TILTED RNA pool in production; omitting the divisor here would score
    # the held fragments against a different length model than the one that actually decides them.
    crossing = crossing_probability_from_index(index, 4096)
    gdna_opp = gdna_opportunity_from_index(index, 4096)
    _, _, node_types = build_node_partition_arrays(index)
    junctions = build_junction_edge_arrays(index)
    t_ids = index.t_df["t_id"].to_numpy()

    names = args.conditions or sorted(p.name for p in args.pilot.iterdir() if p.is_dir())
    rows = []
    for name in names:
        cache = read_scan_cache(args.pilot / name, index)
        payload = cache.payload
        scored = score_held_fragments(
            payload,
            fl_models=build_fl_models(
                payload, junction_opportunity=crossing, gdna_opportunity=gdna_opp
            ),
            # ⭐ Pass 1's own strand model, not calibration's — the second pass runs BEFORE calibration
            # and `rna_sense_frac` is the Beta posterior mean of exactly this.
            rna_sense_frac=cache.strand_model.p_r1_sense,
            node_types=node_types,
            junctions=junctions,
        )
        truth_path = suite / name / "truth_abundances.tsv"
        truth = _read_truth_abundance(truth_path) if truth_path.is_file() else None
        row = census(payload, scored, junctions, t_ids, truth)
        row["condition"] = name
        row["p_r1_sense"] = float(cache.strand_model.p_r1_sense)
        rows.append(row)
        print(f"  {name:<44} {row['n_hypotheses']:>8d} hyps  done")

    print()
    print("═══ D-3 · rho == 0, by hypothesis kind ═══")
    print(f"{'condition':<44} {'spliced':>16} {'⭐ ∅ shipped':>16} {'∅ pre-D-6':>16}")
    print("-" * 96)
    for r in rows:
        z = r["rho_zero"]
        print(
            f"{r['condition']:<44} {z['spliced']['frac']:>15.4f} "
            f"{z['empty_shipped']['frac']:>15.4f} {z['empty_pre_d6']['frac']:>15.4f}"
        )
    print("   ⛔ ONE THING VARIED: the contiguous-edge set. `shipped` is the rule the deposit implies;")
    print("      `pre-D-6` is the strictly-between rule that shipped until 2026-08-02.")

    print()
    print("═══ why the spliced zeros are zero — ⭐ only `annotated_empty` is D-3's population ═══")
    print(f"{'condition':<44} {'unannotated':>12} {'annot_empty':>12} {'of spliced':>11}")
    print("-" * 84)
    for r in rows:
        c = r["cause_spliced"]
        n = r["rho_zero"]["spliced"]["n"]
        print(
            f"{r['condition']:<44} {c['unannotated']:>12d} {c['annotated_empty']:>12d} "
            f"{c['annotated_empty'] / n if n else float('nan'):>10.4f}"
        )

    print()
    print("═══ the consequence, per RECORD ═══")
    print(f"{'condition':<44} {'records':>9} {'undecided':>10} {'frac':>8} {'⭐ decisive':>11} {'frac':>8}")
    print("-" * 96)
    for r in rows:
        d = r["records"]
        print(
            f"{r['condition']:<44} {d['n']:>9d} {d['undecided']:>10d} "
            f"{d['undecided_frac']:>8.5f} {d['decisive_zero']:>11d} {d['decisive_zero_frac']:>8.5f}"
        )
    print("   ⚠ `undecided` is the score vector EMPTYING — the failure D-3 names. `decisive` is the")
    print("      commoner and quieter one: the zero eliminates a hypothesis while another survives, so")
    print("      the draw never sees it. A pseudocount would change `decisive`, not `undecided`.")

    if any("truth" in r["empty_slots"] for r in rows):
        print()
        print("═══ ⭐ IS THE ZERO RIGHT? annotated-empty slots, scored on truth ═══")
        print(f"{'condition':<44} {'slots':>7} {'expressed':>10} {'silent':>8} {'claims@expr':>12} {'frac':>8}")
        print("-" * 96)
        for r in rows:
            t = r["empty_slots"].get("truth")
            if t is None:
                continue
            print(
                f"{r['condition']:<44} {r['empty_slots']['n_distinct']:>7d} "
                f"{t['slots_with_an_expressed_supporter']:>10d} {t['slots_all_supporters_silent']:>8d} "
                f"{t['claims_on_expressed']:>12d} {t['frac_claims_on_expressed']:>8.4f}"
            )
        print("   ⭐ `silent` slots are zero because NOTHING EXPRESSED THEM — the zero is correct and")
        print("      killing the hypothesis is the right answer. `expressed` slots are D-3's real case:")
        print("      the molecule existed and its only evidence is held.")

    if args.json:
        args.json.write_text(json.dumps(rows, indent=2, sort_keys=True))
        print(f"\nwrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
