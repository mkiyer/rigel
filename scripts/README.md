# `scripts/` — developer tooling

Developer tooling that is not part of the shipped `rigel` package. `tests/test_scripts_index.py` holds `design/`
against `CLAUDE.md`'s instrument table, and `profiling/` and the `_`-prefixed helper modules of `design/` and `sim/`
against this file, in both directions, and checks that every script imports and has a docstring.

## The four directories

| dir | what it is |
|---|---|
| `design/` | The instruments: the debug loop, the panel harnesses, the truth-scoring instruments. Indexed by question in `CLAUDE.md`'s table; each module docstring is its manual. The `_`-prefixed helper modules answer no question of their own: `design/_shared.py` (paths and a panel's conditions, `--set`, strata and the release scope, sharding), `design/_oracle.py` (the calibration oracle: the production accumulator run on the BAM split by true fragment origin, checked sum-to-full), `design/_oracle_arms.py` (the oracle cache, the certified per-object truth and the per-object scorer) and `design/_scan_cache.py` (the scan cache every calibration instrument reads) |
| `sim/` | The whole-genome condition grid over the simulator engine in `src/rigel/sim/`, plus the panel configs in `sim/configs/`. `simulate_reads.py` simulates a panel from its YAML through the helper modules `sim/_sim_config.py` (the YAML and its parser), `sim/_sim_abundance.py` (base abundances and the nascent modes), `sim/_sim_grid.py` (the grid, one condition at a time, and the manifest) and `sim/_sim_truth.py` (the truth tables); `design_suite_probes.py` designs the capture panel. `panel.py` drives a panel through build, simulate, cache, score and report; the reference it builds on is carved by `build_suite_reference.py`, and the test chromosome is rendered from `sim/test_reference/test_chr.yaml` by `build_test_reference.py` (recipe: `docs/TESTING.md` §0a) |
| `profiling/` | Where the time and the memory go; this row is its index. `profiling/profiler.py` runs the whole pipeline on one library as a tree of named stages (the index load, the scan, the second pass, each index-derived geometry build, every calibration stage down to each sweep's builders and passes, every quant stage down to the locus EM), with calls, inclusive and self seconds, and peak and held RSS per stage. `--set SECTION.FIELD=VALUE` overrides any `PipelineConfig` field, `--scan-only` runs the scan alone for tuning its knobs, `--compare A.json B.json` joins two reports stage by stage, and `--self-test` falsifies the probe machinery with no BAM. Whether a change moved any number is `design/rename_identity.py --bam`'s question. Profile real libraries, never a panel condition |
| `publishing/` | The release scripts: `release.sh` cuts the release, `post_release.sh` finishes it once PyPI has the sdist, `conda_publish.sh` uploads a conda package |

## The simulation and benchmarking workflow

One command per stage, each resumable; `panel.py status` names the next one.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1

python scripts/sim/panel.py build    --config scripts/sim/configs/gdna_ladder.yaml
python scripts/sim/panel.py simulate --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py cache    --config scripts/sim/configs/gdna_ladder.yaml --jobs 8
python scripts/sim/panel.py score    --config scripts/sim/configs/gdna_ladder.yaml --jobs 2   # ~8.5 GB per condition
python scripts/sim/panel.py report   --config scripts/sim/configs/gdna_ladder.yaml
```

`score` needs both caches and `cache` builds both: the scan cache makes calibration re-runnable without
rescanning, and the oracle cache is the origin-split truth that every truth-scoring instrument reads (`gdna` /
`mrna` / `nrna` partitions, checked sum-to-full against the scan cache, with the RNA split again by transcript
strand into `rna_pos` / `rna_neg`).

## Lint

```bash
ruff check src/ tests/ scripts/
```

Everything under `scripts/` must pass. Never `ruff format scripts/`: the formatter destroys the instruments'
aligned tables and banner comments.
