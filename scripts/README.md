# `scripts/` — developer tooling

Developer-facing tooling that is **not** part of the shipped `rigel` package.

⚠ **This file was rewritten on 2026-08-11 because it had rotted into fiction.** It documented a
`benchmarking/` package (deleted the same day), a `calibration/` directory and a Tier-2 `debug/`
toolkit — **neither of the last two has existed for months** — and it did not mention `design/` at all,
which is 56 files and the toolkit people actually use. ⛔ A directory listing takes one second; prefer it
to this file if the two disagree, and then fix this file.

## The four directories

| dir | what it is |
|---|---|
| ⭐⭐ `design/` | **the instrument shelf — 56 files.** The debug loop, the panel harnesses, the truth-scoring instruments. ⛔ Indexed in `CLAUDE.md`'s table, and `tests/test_scripts_index.py` enforces that the index and the disk agree in BOTH directions and that every file still IMPORTS |
| ⭐ `sim/` | thin CLI wrappers over the simulator engine in `src/rigel/sim/`, plus the panel YAML configs in `sim/configs/`. ⭐⭐ **`panel.py` is the one entry point** — build, simulate, cache, score, report |
| `profiling/` | profiling drivers (`profiler.py`, `pyspy_driver.py`, `scan_profile.py`) + configs |
| `publishing/` | release scripts (`release.sh`, `post_release.sh`, `conda_publish.sh`) |

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
