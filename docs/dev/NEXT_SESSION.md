# NEXT SESSION — start here (2026-09-13, after W5 landed; uncommitted on top of `8fa0e684`)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C. This file is only how to begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (`8fa0e684` plus the uncommitted W5 landing, or the commit that carries it). Read
> `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md` and `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` — the plan
> holds the whole picture; its PRE-PORT WORKLIST (W1–W9) is the session's agenda, W1–W5 are done, and the
> next item is W6, the tunables census. Run `scripts/design/preflight.py` and the suite before touching
> anything; `CLAUDE.md`'s baseline line is the count to reproduce (3,435 passed / 2 xfail / 3,437 collected).
>
> THE STANDING RULINGS. Bit-identity is no longer the bar for the Python work: the bar is elegant, simple,
> efficient, clear, concise, maintainable code, judged on the oracle metric (`calibration_vs_oracle.py`,
> per stratum, both zero controls), the panel (`policy_benchmark.py --panel test|ladder`, the two halves
> apart), the suite, and timing and memory on back-to-back profiler pairs. A pure restructure is still
> proven bit-identical (`rename_identity.py --check` against a reference frozen from the current tree —
> the `memory_identity_*` references predate the one-lattice landing and must be re-frozen first;
> `sweep_replay.py replay --dir … --call N --tolerance`, whose step-4 captures likewise predate it and
> now replay as a recorded design change, not a bit gate). One mechanism at a time. No magic numbers:
> a knob is derived or measured, its ladder recorded in its docstring, and neither of the two lattice
> knobs is exposed to users. Present each step's plan as what will change and what gates it, and wait
> for my go. The owner drives commits. Never `ruff format scripts/`. Patch `calibrate` through
> `importlib.import_module("rigel.calibration.calibrate")`, never `import … as`. Any config value is an
> arm through `--set SECTION.FIELD=VALUE` on `calibration_vs_oracle.py`, `policy_benchmark.py` and
> `profiler.py` (one parser, `scripts/design/_shared.set_field`), so nothing in `src/` moves to price one.
>
> W6 IS THE TUNABLES CENSUS: `CalibrationConfig`'s fields, each classed live / derived / dead with its
> derivation or its measured ladder named, and anything unearned removed; judged by `module_census.py`
> and the suite. The two lattice fields are the model for what a kept knob looks like
> (`sweep_logodds_step`, `sweep_n_tilt`: dimensionless, a stated guarantee, a recorded ladder, not on the CLI).

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## What W5 landed (2026-09-13; the ruling is `DESIGN.md` §6b.15, the refusals `ISSUES: the-second-lambda-grid-and-its-regrid`)

One λ lattice for every consumer of ψ's grid, parametrised by its step: `CalibrationConfig.sweep_logodds_step`
= 0.2 (101 points at the floor bracket, ~220 on the refits; `calibrate.lattice_points` derives K at every
bracket, so the step is the invariant the retired `_scaled_grid` held). `sweep_n_grid`,
`sweep_n_grid_single_strand`, `_regrid_global`, `CompositionPriors.regrid`, `_scaled_grid` and the CLI's
`--sweep-n-grid-single-strand` are gone; `sweep_n_tilt` is an explicit 60, decoupled from K. Judged: the
oracle metric at par or better in scope (0.993 / 0.998 / 1.000 of the retired pair, `g00` 0.994, deferred
+1.3 %); the landed tree byte-identical to the study's own `--set` arm on all 46 conditions; the panel's
two bars unchanged on the ladder and repaired on the test chromosome (15/20 → 19/20, worst row 1.22× →
1.00×); the deep library at wall 1.066×, `calibrate` 1.08×, peak +1.8 GB against a same-session baseline
(the sweeps' own ψ solves 0.96–0.97×, the passes 1.10×; the old tree's `init_beliefs` never received a tilt
count and ran its pre-sweep cube at K_t = K, which the landing corrects — that solve is 25 s at 60). 21 goldens regenerated, magnitudes read
first (counts ≤ 1.5e-3 relative; a tiny toy's `em_effective_length` ≤ 8.4 %). The scratch record of every
arm is this session's `w5/FINDINGS.md` (not in the repo).

Decisions on record from the study: K_t stays 60 — 30 and 15 break the zero control on ψ's own θ
quadrature (`ISSUES: theta-quadrature-at-zero-gdna`, the next design item touching the cube); a step
finer than 0.2 buys ≤ 1.3 % in scope for 1.33–1.71× wall because the cube scales with K × K_t; a read-out
insensitive to K and a step derived from the sharpest posterior were both refused with their numbers.

A lead found on the way, not taken (one mechanism at a time): `region_geometry.init_beliefs` solves every
slot, but `_type_belief` keeps the signature-binary `{0,0,1}` for AMBIG slots and discards their cube — that
pre-sweep cube is dead work (part of the 25 s init ψ on the deep library), removable bit-identically.

## Where the worklist stands

W1–W5 DONE. Next: W6 the tunables census, W7 the capture as a typed record, W8 the vocabulary rulings,
W9 the two xfails. The port begins only when the owner is satisfied with the Python.

## Decisions on record

* Float64 for the whole of ψ; ONE solver (a single-strand slot is the cube with a one-cell tilt grid);
  ONE λ lattice (this landing).
* The arcsine coordinate stays refused; the vertex atom is parked with its ceiling (≤ 1 % on stranded
  in-scope rows).
* The message cache's on/off switch (plan §E) is the owner's; it holds ~4 GB of the 11.4 GB peak.
* The refit count stays; the scan's thread split is the owner's (`ISSUES: scan-thread-split-starves-the-workers`).
* Parallelism waits for the port; threads are refuted for the Python passes.
