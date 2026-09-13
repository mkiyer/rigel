# NEXT SESSION — start here (2026-09-13, after W7; W8 and W9 ruled; W10–W13 are the road to the port)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C, and §8 is the design of the
remaining steps, W10–W13. This file is only how to begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, the commit after `4dc8776e`). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md`
> and `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` — the plan holds the whole picture; its PRE-PORT WORKLIST
> (W1–W13) is the session's agenda: W1–W7 are done and committed, W8 (the renames) and W9 (the two xfails)
> are ruled out of this thread by decision, and §8 is the design of what remains before the port — W10 the
> hygiene ledger's pure cleanups, W11 the coverage census, W12 the θ quadrature at zero gDNA, W13 the port's
> prerequisites — in that order. Run `scripts/design/preflight.py` and the suite before touching anything;
> `CLAUDE.md`'s baseline line is the count to reproduce (3,435 passed / 2 xfail / 3,437 collected).
>
> THE STANDING RULINGS. The bar is elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel
> (`policy_benchmark.py --panel test|ladder`, the two halves apart), the suite, and timing and memory on
> back-to-back profiler pairs. A pure restructure or a rename IS proven bit-identical: freeze fresh
> references first (`rename_identity.py --freeze`; every reference on disk predates the one-lattice landing),
> then `--check --stage <name>` after each stage; a removal is proven on the default rows of both substrates
> (re-record them first: `calibration_vs_oracle.py --json` on every condition of both panels, the pattern of
> the last session's `run_arms.py`). One mechanism at a time. No magic numbers. Present each step's plan as
> what will change and what gates it, and wait for my go. The owner drives commits. Never `ruff format
> scripts/`. Patch `calibrate` through `importlib.import_module("rigel.calibration.calibrate")`, never
> `import … as`. Any config value is an arm through `--set SECTION.FIELD=VALUE` on
> `calibration_vs_oracle.py`, `policy_benchmark.py` and `profiler.py`.
>
> W10 FIRST, each item its own commit, none moving the metric: the stale comment above the second pass's
> `build_fl_models`; `mass_*_boundary` → `count_*_boundary`; the six bank-readers to the drained frame, one
> at a time with a recorded before/after; the moment tests deleted with the length channel. Then W11, the
> coverage census — a measurement, then my decisions. Then W12, the θ quadrature — DERIVE first (the bias
> is a λ-dependence of the quadrature's error), PROTOTYPE outside `src/`, A/B on both panels with both zero
> controls. Then W13, the re-capture of the deep library's sweeps for the port's tolerance gate. Then the port.

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## Where the worklist stands

W1–W7 DONE and committed (`2e92c6fd` W5, `c54cc5e1` W6, `e15fc0df` W7). W8 RULED (no rename; the terms
stay). W9 RULED (both xfails deferred: the antisense casualty to the prior-assembly session, the two-sided
exon row to the calibration-accuracy thread; both stay as executable records). Next: W10, W11, W12, W13 (plan
§8), then the port (plan §F). The port begins only when the owner is satisfied with the Python.

## What the last session landed (all committed)

* W5 — one λ lattice, `sweep_logodds_step` 0.2, `sweep_n_tilt` 60 explicit; the second grid, its regrid and
  the bracket scaling gone; `--set` on the three instruments (`DESIGN.md` §6b.15).
* W6 — the tunables census: 11 fields → 7, four switches gone byte-identically, the message policy one
  field, the oracle instrument's flags folded into `--set`, the dead pre-sweep AMBIG cube removed;
  `background_abundance` kept as an unruled decision (`ISSUES: background-abundance-pair-unruled`).
* W7 — the diagnostic capture as `blocks.SweepCapture` (25 fields), readers on fields, the extra ψ solve gone.

Two facts met on the way, not fixed: `tests/calibration/test_solvability_audit.py` and `test_worst_objects.py`
cannot be collected ALONE (the instruments they load import `_shared`, which an earlier test puts on the path
in a full run) — a one-line path insert in their loader closes it, a W10-sized item if wanted; and the θ
quadrature's zero-gDNA sensitivity is W12.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee.
* The two lattice knobs are the model of a kept knob: dimensionless, a guarantee, a recorded ladder, not on
  the CLI. The EM's `gdna_em_llr_bias` stays (the owner's).
* `drain`, `row` and `face` stay (2026-09-13). The arcsine coordinate stays refused; the vertex atom is parked.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
