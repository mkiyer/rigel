# `scripts/` — developer tooling

Developer-facing tooling that is **not** part of the shipped `rigel` package. Two tiers, with
different maintenance standards.

## Tier 1 — maintained tooling (lint-clean, no relaxations)

Referenced by the CI/skills/docs and expected to keep working across releases:

| dir | what it is |
|---|---|
| `sim/` | synthetic-suite drivers (`simulate_suite.py`, `evaluate_suite.py`, `simulate_reads.py`, `gen_ambig_dense.py`, …) + suite YAML configs — used by the `calibration-benchmark` skill |
| `benchmarking/` | the real-data benchmark package, run as `python -m scripts.benchmarking` |
| `publishing/` | release scripts (`release.sh`, `post_release.sh`, `conda_publish.sh`) |
| `profiling/` | profiling drivers (`profiler.py`, `pyspy_driver.py`, `scan_profile.py`) + configs |
| `calibration/` | calibration config generators |

These are held to the same lint standard as `src/` and `tests/`.

## Tier 2 — `debug/`: a working diagnostics toolkit

A **small, named** set of current calibration/EM diagnostics (oracle + net-flow + siphon tools). The
canonical tools — `oracle.py`, `oracle_reattribute.py`, `oracle_leak_trace.py`, `benchmark_ab_report.py`,
`benchmark_ab_render.py`, `dissect_loci.py`, `toy_prod.py`, `pass_trace.py`, `localize_siphon.py`,
`layer_trace.py`, `locus_component_audit.py`, `em_multiplicity*.py`, … — are documented in
`docs/TRAPS.md`.

**Policy (so this dir does not regrow into an archive):**

- Keep it a small toolkit of tools that are **current**. When a diagnostic is finished with, **delete it**
  — git history preserves it. Do **not** accumulate a `scripts/debug/archive/` (one was deleted in the
  v0.7.0 productionizing pass; ~110 one-off scripts, all recoverable from git).
- Scratch/exploration that isn't a reusable tool should not be committed here at all.
- `debug/` scripts may use quick-diagnostic idioms (one-liner `if x: y`, `a = ...; b = ...`, post-`sys.path`
  imports, compute-and-inspect locals); those specific lint rules are relaxed for `scripts/debug/**` in
  `pyproject.toml` (`[tool.ruff.lint.per-file-ignores]`). **Real errors — undefined names, redefinition,
  syntax — are still enforced**, so a script that has rotted against the current API fails lint and is a
  delete candidate.

## Lint

```bash
ruff check src/ tests/ scripts/
```

Everything under `scripts/` must pass (with the `debug/` relaxations above).
