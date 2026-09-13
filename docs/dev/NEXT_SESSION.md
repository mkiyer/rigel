# NEXT SESSION — start here (2026-09-13, after W7; W8 prepared)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C, and §7 is the design of
the next item, W8. This file is only how to begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, the commit after `e15fc0df`). Read `CLAUDE.md`, then
> `docs/dev/NEXT_SESSION.md` and `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` — the plan holds the whole
> picture; its PRE-PORT WORKLIST (W1–W9) is the session's agenda, W1–W7 are done and committed, and §7 is
> the design of W8, the vocabulary rulings and the hygiene ledger. Run `scripts/design/preflight.py` and the
> suite before touching anything; `CLAUDE.md`'s baseline line is the count to reproduce (3,435 passed /
> 2 xfail / 3,437 collected).
>
> THE STANDING RULINGS. The bar is elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel
> (`policy_benchmark.py --panel test|ladder`, the two halves apart), the suite, and timing and memory on
> back-to-back profiler pairs. A rename is a pure restructure and IS proven bit-identical: freeze fresh
> references first (`rename_identity.py --freeze`; every reference on disk predates the one-lattice landing),
> then `--check --stage <name>` after each stage. One mechanism at a time. No magic numbers. Present each
> step's plan as what will change and what gates it, and wait for my go. The owner drives commits. Never
> `ruff format scripts/`. Patch `calibrate` through `importlib.import_module("rigel.calibration.calibrate")`,
> never `import … as`. Any config value is an arm through `--set SECTION.FIELD=VALUE` on
> `calibration_vs_oracle.py`, `policy_benchmark.py` and `profiler.py` (one parser,
> `scripts/design/_shared.set_field`).
>
> W8 BEGINS WITH THREE RULINGS THAT ARE MINE, then the mechanics. Plan §7 has the census counts and the
> candidates: for `drain` the recommendation is `settle`; for `row` it is `profile`, with the question of
> whether a delivered row and a level lane's `profile` are one word; for `face` it is to keep the word and
> close the entry. Take my rulings, then: freeze the identity references; teach `rename_census.py` the three
> tokens and dump every site; rename one token at a time, docs in the same commit; then the hygiene ledger
> in §7's order, each its own commit. Then W9, the two xfails.

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## Where the worklist stands

W1–W7 DONE and committed (`2e92c6fd` W5, `c54cc5e1` W6, `e15fc0df` W7). Next: W8 (plan §7), W9 the two
xfails; and a coverage census of `src/` (the suite and every instrument's self-test under coverage, the
never-executed lines reviewed one by one) is proposed as the instrument for dead code beyond the config's
switches. The port begins only when the owner is satisfied with the Python.

## What the last session landed (all three committed)

* W5 — one λ lattice, `sweep_logodds_step` 0.2, `sweep_n_tilt` 60 explicit; the second grid, its regrid and
  the bracket scaling gone; `--set` on the three instruments (`DESIGN.md` §6b.15).
* W6 — the tunables census: 11 fields → 7, four switches gone byte-identically, the message policy one
  field, the oracle instrument's flags folded into `--set`, the dead pre-sweep AMBIG cube removed;
  `background_abundance` kept as an unruled decision (`ISSUES: background-abundance-pair-unruled`).
* W7 — the diagnostic capture as `blocks.SweepCapture` (25 fields), readers on fields, the extra ψ solve gone.

Two facts met on the way, not fixed: `tests/calibration/test_solvability_audit.py` and `test_worst_objects.py`
cannot be collected ALONE (the instruments they load import `_shared`, which an earlier test puts on the
path in a full run) — a one-line path insert in their loader closes it; and the θ quadrature's zero-gDNA
sensitivity is filed (`ISSUES: theta-quadrature-at-zero-gdna`), the next design item touching the cube.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee.
* The two lattice knobs are the model of a kept knob: dimensionless, a guarantee, a recorded ladder, not on
  the CLI. The EM's `gdna_em_llr_bias` stays (the owner's).
* The arcsine coordinate stays refused; the vertex atom is parked with its ceiling.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
