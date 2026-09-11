# `scripts/` — developer tooling

Developer-facing tooling that is not part of the shipped `rigel` package. Membership is gated:
`tests/test_scripts_index.py` holds `design/` against `CLAUDE.md`'s instrument table and `profiling/`
against this file, in both directions, and checks that every script imports and has a docstring.

## The four directories

| dir | what it is |
|---|---|
| `design/` | The instrument shelf: the debug loop, the panel harnesses, the truth-scoring instruments. Indexed by question in `CLAUDE.md`'s table; `docs/SUCCESS.md` has the run order |
| `sim/` | Thin CLI wrappers over the simulator engine in `src/rigel/sim/`, plus the panel YAML configs in `sim/configs/`. `panel.py` is the one entry point: build, simulate, cache, score, report |
| `profiling/` | Where the time and the memory go, and this row is its index. `profiling/profiler.py` runs the whole pipeline on one library as a tree of named stages (the index load, the scan, the second pass, each index-derived geometry build, every calibration stage down to each sweep's builders and passes, every quant stage down to the locus EM), with calls, inclusive and self seconds, peak and held RSS per stage; `--set SECTION.FIELD=VALUE` overrides any `PipelineConfig` field, `--scan-only` runs the scan alone for tuning its knobs, `--compare A.json B.json` joins two reports stage by stage, and `--self-test` falsifies the probe machinery with no BAM. `profiling/sweep_replay.py` captures every calibration sweep's inputs from one real run and replays one sweep in isolation, timing it and comparing its result bit for bit, which is the loop for optimising the sweep. Whether a change moved any number on the whole pipeline is `design/rename_identity.py --bam`'s question. Profile the real libraries, never a panel condition (`docs/TESTING.md` §7) |
| `publishing/` | Release scripts (`release.sh`, `post_release.sh`, `conda_publish.sh`); `docs/PUBLISHING.md` is the procedure |

## The simulation + benchmarking workflow

One command per stage, each resumable, each gated. `docs/TESTING.md` §2 is the long form.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1

python scripts/sim/panel.py build    --config scripts/sim/configs/gdna_ladder.yaml
python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py cache    --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py score    --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py report   --config scripts/sim/configs/gdna_ladder.yaml
```

`score` needs both caches and `cache` builds both: the scan cache makes calibration re-runnable without
rescanning, and the oracle cache is the origin-split truth (`gdna` / `mrna` / `nrna` partitions plus the
undrained `_main` payload) that every truth-scoring instrument reads.

## Lint

```bash
ruff check src/ tests/ scripts/
```

Everything under `scripts/` must pass. Never `ruff format scripts/`: the instruments' aligned tables and
banner comments are load-bearing and the formatter destroys them.
