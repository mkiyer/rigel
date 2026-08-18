# `scripts/` — developer tooling

Developer-facing tooling that is **not** part of the shipped `rigel` package.

⚠ **This file was rewritten on 2026-08-11 because it had rotted into fiction.** It documented a
`benchmarking/` package (deleted the same day), a `calibration/` directory and a Tier-2 `debug/`
toolkit — **neither of the last two has existed for months** — and it did not mention `design/` at all,
which is 56 files and the toolkit people actually use. ⛔ A directory listing takes one second; prefer it
to this file if the two disagree, and then fix this file.

⚠ **Every count in the table below is hand-carried prose and nothing gates it, so re-derive rather than
trust** (`TRAPS: re-record-the-baseline`): `ls scripts/<dir>/*.py | grep -v __init__ | wc -l`. What IS
gated is the membership — `tests/test_scripts_index.py` holds `design/` against `CLAUDE.md`'s table and
`profiling/` against this file, in both directions.

## The four directories

| dir | what it is |
|---|---|
| ⭐⭐ `design/` | **the instrument shelf — 59 files** (re-derived 2026-08-17, after `flgap_study_cache.py` and `prior_yardstick.py` were deleted; the row said 56, and its own rule is to prefer the directory listing and then fix the file). The debug loop, the panel harnesses, the truth-scoring instruments. ⛔ Indexed in `CLAUDE.md`'s table, and `tests/test_scripts_index.py` enforces that the index and the disk agree in BOTH directions and that every file still IMPORTS |
| ⭐ `sim/` | **8 files** — thin CLI wrappers over the simulator engine in `src/rigel/sim/`, plus the panel YAML configs in `sim/configs/`. ⭐⭐ **`panel.py` is the one entry point** — build, simulate, cache, score, report |
| ⭐ `profiling/` | **where the time and the memory go — 2 files, and this row is their index.** ⭐⭐ **`profiling/profiler.py`** — the whole pipeline, wall clock AND per-phase peak RSS across `scan` / `calibrate` / `quant`, plus `held` (what a phase hands ON rather than borrows), with optional cProfile. ⭐ **`--sweep SECTION.FIELD=V1,V2,…` puts any `PipelineConfig` field on the x-axis**, which is how compute is priced against memory: `calibration.sweep_n_grid` 30 → 120 costs **5.2 s → 16.2 s** in calibrate and **+1.1 GB** peak (10 M fragments, 4 threads). `--self-test` falsifies the memory instrumentation with no BAM and no index (3/3, and it pins the one hole it cannot cover). ⭐ **`profiling/scan_profile.py`** — the scan alone, swept across thread and chunk budgets, and the process to launch under Instruments/xctrace or py-spy. ⛔⛔ **PROFILE A HIGH-DEPTH REAL RNA-seq LIBRARY, NOT cfRNA AND NOT A PANEL CONDITION** (owner, 2026-08-17) — `docs/TESTING.md` §7 is the ruling |
| `publishing/` | release scripts (`release.sh`, `post_release.sh`, `conda_publish.sh`) |

### ⭐⭐ `profiling/` — the decision, taken 2026-08-17

**Extend the gate to a third tree; do not delete the directory.** The performance work before 0.8.0
needs these two instruments, and the reason the choice had to be made at all is that the tree had
**rotted twice** with the one defect class `tests/test_scripts_index.py` already catches next door —
`pyspy_driver.py` read `sys.argv[1]` at import time (found 2026-08-11), `scan_profile.py` imported two
names `profiler.py` did not export so even `--help` raised (found 2026-08-17). Only the gate's REACH
differed (`TRAPS: a-green-suite-hid-five-dead-instruments`). `ALL_SCRIPTS` is now
`design + sim + profiling`, the docstring gate covers all three, and **this row is the profiling tree's
index** — a file there that no row names fails the suite, and a row naming a deleted file fails it too.

**Deleted in the same pass, as legacy:** `configs/` (8 YAML files naming a `--config` flag that no
instrument has, a `scripts/profiler.py` path that does not exist, a `benchmark.py` deleted months ago,
and `em_mode` / `em_iterations` / `em_convergence_delta` keys absent from `src/rigel/config.py`); and
`pyspy_driver.py`, which was a second way to run the pipeline once — `profiler.py` does that and more,
so `py-spy record -o pipeline.svg -- python scripts/profiling/profiler.py --bam … --index …` replaces
it exactly.

⛔ **A BASENAME STILL COLLIDES, AND THE RENAME IS AN OWNER CALL.** `design/scan_profile.py` (the
accumulator's ns/fragment, regressed over several BAMs) and `profiling/scan_profile.py` (the scan's own
wall time and RSS across thread budgets) are different instruments sharing a filename. The gate no
longer *breaks* on it — case ids are tree-qualified (`design/scan_profile.py`,
`profiling/scan_profile.py`), so a failure names the tree — but a reader who greps the basename still
finds two files. ⭐ Proposed, not done: rename `design/scan_profile.py` to `accumulator_cost.py`, which
is what it measures; it has a row in `CLAUDE.md` that moves with it.

## ⭐⭐⭐ The simulation + benchmarking workflow

**One command per stage, each resumable, each gated.** `docs/TESTING.md` §2 is the long form.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1

python scripts/sim/panel.py build    --config scripts/sim/configs/gdna_ladder.yaml
python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py cache    --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py score    --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py report   --config scripts/sim/configs/gdna_ladder.yaml
```

⛔ **`score` needs BOTH caches and `cache` builds both.** The scan cache makes calibration re-runnable
without rescanning; the ORACLE cache is the origin-split truth (`gdna` / `mrna` / `nrna` partitions plus
the undrained `_main` payload) that every truth-scoring instrument reads. Until 2026-08-11 the oracle
cache had no first-class step at all — it was a side effect of running `pass0_vs_oracle.py
--oracle-cache`, which is why the documented recipe could not reach a scored table.

## Lint

```bash
ruff check src/ tests/ scripts/
```

Everything under `scripts/` must pass. ⛔ **Never `ruff format scripts/`** — the instruments' aligned
tables and banner comments are load-bearing and the formatter destroys them.
