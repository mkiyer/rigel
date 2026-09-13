# NEXT SESSION — start here (2026-09-12, after commit `2ebfaea2`)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`; the PRE-PORT WORKLIST is the section before its §C, and §6 is the design
of the next item, the grid study. This file is only how to begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, `2ebfaea2`). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md` and
> `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` — the plan holds the whole picture; its PRE-PORT WORKLIST
> (W1–W9) is the session's agenda, W1–W4 are done and committed, and §6 is the design of W5, the grid
> study. Run `scripts/design/preflight.py` and the suite before touching anything; `CLAUDE.md`'s baseline
> line is the count to reproduce (3,435 passed / 2 xfail / 3,437 collected).
>
> THE STANDING RULINGS. Bit-identity is no longer the bar for the Python work: I am not concerned with
> minuscule changes; the bar is elegant, simple, efficient, clear, concise, maintainable code, judged on
> the oracle metric (`calibration_vs_oracle.py`, per stratum, both zero controls), the panel
> (`policy_benchmark.py --panel test|ladder`), the suite, and timing and memory on back-to-back profiler
> pairs. A pure restructure is still proven bit-identical (`rename_identity.py --check` against
> `~/Downloads/rigel_runs/arms/memory_identity_*.json`; `sweep_replay.py replay --dir
> ~/Downloads/rigel_runs/perf/sweeps_MO_3021_step4 --call 0..3 [--tolerance]`). One mechanism at a time.
> No magic numbers. Present each step's plan as what will change and what gates it, and wait for my go.
> The owner drives commits. Never `ruff format scripts/`. Patch `calibrate` through
> `importlib.import_module("rigel.calibration.calibrate")`, never `import … as`.
>
> W5 IS A CAREFUL SEARCH OVER THE GRID, NOT A TWO-ARM A/B. Today ψ runs on two λ grids — the coarse one
> (`sweep_n_grid` 60, scaled to 138 by the bracket) for the AMBIG cube, the message rows and the factory
> rows, and the fine one (`sweep_n_grid_single_strand` 256 → 513) for the single-strand read-out, with
> `_regrid_global` interpolating priors and rows between them per tile. Do not hastily conclude that two
> grids are needed. Understand the grid's behaviour first: accuracy on the metric as a function of K, per
> stratum and per slot class; where the median read-out's quantisation sits and how it scales with K and
> with mass; the coupling of K to the bracket (`_scaled_grid` holds dλ fixed); what the regrid costs in
> accuracy; and cost — time, and memory, since every message row and the cache scale with K. Then look
> for the simple, elegant design — one grid at a derived K, or a read-out that makes the grid's size
> irrelevant — and find it with tenacity. The first instrument change is `calibration_vs_oracle.py
> --set SECTION.FIELD=VALUE` (as `profiler.py` has), so a grid arm is a config value and nothing in
> `src/` moves until a ruling is reached. Then W6–W9.

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

The bit-identity baselines are the memory-step tree (2026-09-12): `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step4`
(four captured sweeps; `sweep_replay.py replay --dir … --call 0..3 [--tolerance]`, `--block-slots N|none`;
the refit sweeps' captures carry their message cache, so their replay exercises the cache path and ψ,
and sweep 0 the whole message layer) and the three references `~/Downloads/rigel_runs/arms/memory_identity_*.json`
(`rename_identity.py --check`; the real BAMs are under `~/Downloads/rigel_runs/cfrna/mctp_<lib>_*/bam/`).
Every older capture and reference (`sweeps_MO_3021*` before step4, `locus_identity_*`,
`onesolver_identity_*`, `cleanup_identity_*`, `perf_identity_LBX0190.json`) describes earlier code.

## Where the worklist stands

W1 (the replay's tolerance report), W2 (one ψ solver in float64, the unread strand variances deleted),
W3 (`calibrate`, `_solve_block`, `solve_chain` as named stages) and W4 (memory: peak 19.2 → 11.4 GB)
are DONE and committed (`987ce4e7`, `d5cfd771`, `2ebfaea2`). Next: W5 the grid study (plan §6), then W6
the tunables census, W7 the capture as a typed record, W8 the vocabulary rulings, W9 the two xfails. The
port begins only when the owner is satisfied with the Python.

## Decisions on record

* Float64 for the whole of ψ; ONE solver (a single-strand slot is the cube with a one-cell tilt grid).
* The arcsine coordinate stays refused; the vertex atom is parked with its ceiling (≤ 1 % on stranded
  in-scope rows).
* The message cache's on/off switch (plan §E) is the owner's; it holds ~4 GB of the 11.4 GB peak.
* The refit count stays; the scan's thread split is the owner's (`ISSUES: scan-thread-split-starves-the-workers`).
* Parallelism waits for the port; threads are refuted for the Python passes.
