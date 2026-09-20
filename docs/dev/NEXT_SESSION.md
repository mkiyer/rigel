# NEXT SESSION — the twice-counted crossing is repaired in snapshots; the ladder and the report are the next reads (2026-09-20)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md`, what "done" means is `docs/SUCCESS.md`, the lessons are
`docs/TRAPS.md` cited by name.

## Where the tool is

Everything through `f77d43e1` is landed. The session of 2026-09-20 changed ONE factor in `src/`:
`assemble_priors` converts the gDNA component's crossing support by the accumulator's own `q`, as it converts
the crossing count (`EQUATIONS.md` §11). It is PREPARED AS SNAPSHOTS AND AWAITS THE OWNER'S GO
(`~/Downloads/rigel_runs/prototypes/2026-09-20_siphon_mechanism/commits/`), with the five goldens that moved
(0.2–2.1 % on two- and three-transcript scenarios), two new gates, one rewritten gate, and the doc moves.

## What was measured, in one paragraph

The capture-ON nascent siphon is the gDNA component priced at `q̄` of its density: the opportunity counted a
crossing start at every boundary its fragment crosses, the pseudocount counted the fragment once, and under
capture the crossing support is 64–84 % of a probed locus's opportunity. Proven with one E-step from the TRUE
counts (off capture a fixed point, on capture +24 % per step, 73 % of it the opportunity); the twins, the init,
the rulers, the fragment-length term, the priors and the scorer's length law were each switched off and cleared
with numbers. The repair takes `g50 ss.99 ON` from +541,216 to +32,908 nascent and −626,550 to −56,319 gDNA for
0.3 points of transcript error, with `calibration_vs_oracle.py`, `zero_controls.py` and `policy_benchmark.py`
identical. At `on_fraction` 0.10 the siphon is +522,205, so the repair is worth the full amount.

## What to do next, in order

1. Read the full ladder re-measured on the patched tree (`arms/qa_ladder_base_q.jsonl`, `_base_reseed_q`) per
   stratum, keep `g98` apart, and rebuild the report with the `ladder-report` skill once the snapshots land.
2. The transcript table's 0.3-point cost at `g50 ss.99 ON` is the isoform ruler's +6–9 % over-statement standing
   unmasked (`ISSUES: ruler-witness-geometry-on-transcript-panels`, `ruler_vs_truth.py`): the two errors had
   cancelled. That is now the largest in-scope item on the transcript number.
3. The residual +9 % per-step drift with the opportunity at its oracle value is the shadow as an unpinned
   hypothesis at the probed exons: the measured per-transcript prior on components with exclusive evidence
   (`ISSUES: per-transcript-prior-lane`), second.
4. Two owner decisions: the index's overlapping synthetic shadows (`ISSUES: overlapping-synthetic-shadows`)
   and `g98`'s RNA prior floor (`ISSUES: rna-prior-floor-at-pure-gdna-loci`).

## Where everything is

* The session's harnesses, traces with true origins, the truth E-steps and the A/Bs:
  `~/Downloads/rigel_runs/prototypes/2026-09-20_siphon_mechanism/` (`baseline/` and `ab/` are the four
  instruments before and after, from the same tree state; `patch/` the applied edit and runners).
* The 0.10 condition: `~/Downloads/rigel_runs/suite/ladder_nrna_lo/` (one condition, cached, certified).
* The teaching page: the Artifact "The Twice-Counted Crossing".
