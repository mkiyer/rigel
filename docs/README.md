# Rigel documentation

Rigel is a Bayesian RNA-seq quantifier that jointly models mRNA, nascent RNA (nRNA), and genomic
DNA (gDNA) contamination. This is the documentation map. The **code is the source of truth**; these
docs describe the current system (v0.7.0).

## Start here

| doc | what it is |
|---|---|
| [`MANUAL.md`](MANUAL.md) | user manual — install, the `index` / `quant` / `sim` / `export` CLI, options, outputs |
| [`METHODS.md`](METHODS.md) | the scientific method (the mRNA/nRNA/gDNA model + the three pipeline stages) |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | the code map and data flow (scan → calibrate → quant; Python modules + C++ extensions) |
| [`parameters.md`](parameters.md) | configuration & CLI parameter reference |
| [`BENCHMARKING.md`](BENCHMARKING.md) | the synthetic benchmark suite + net-fragment-flow evaluation |
| [`SIMULATOR.md`](SIMULATOR.md) | the synthetic read simulator (`rigel sim`) |
| [`PUBLISHING.md`](PUBLISHING.md) | building + publishing `rigel-rnaseq` to PyPI |
| [`TODO.md`](TODO.md) | open items |

## Calibration (`calibration/`)

The gDNA-vs-RNA deconvolution stage: a bipartite region↔boundary **belief-propagation sweep** (single
forward-backward pass) over the accumulator payload. Three living references:

| doc | what it is |
|---|---|
| [`CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) | **★ AUTHORITATIVE** — the information & precision model: the count-zero-information principle and the three composition sources (strand likelihood, cross-node imputation, global gDNA prior). Read this first; if any other note disagrees, this wins. |
| [`calibration_prior_production_reference.md`](calibration/calibration_prior_production_reference.md) | what actually ships on `main`: the 2-pass per-node solve, the M2 background prior + the Phase-2 gDNA-density KDE (product-of-experts), and how the per-locus Dirichlet prior is assembled |
| [`oracle_and_benchmarking.md`](calibration/oracle_and_benchmarking.md) | the one trusted debug/benchmark path — the accumulator-consistent oracle (`scripts/debug/oracle.py`) and the net-fragment-flow metric |
| [`archive/`](calibration/archive/) | the long calibration-evolution history (~120 superseded design/plan/diagnosis docs) — kept for provenance, **not** current |

## Accumulator (`accumulator/`)

The C++ fractional accumulator (the calibration sufficient statistics).

| doc | what it is |
|---|---|
| [`00_design.md`](accumulator/00_design.md) | **the spec** — per-region/boundary mass + FL pools; matched byte-for-byte by `tests/native/_accumulator_reference.py` |
| `01_implementation.md`, `01_junction_strand.md`, `accumulator_fragment_span_redesign.md`, `audit_phase1.md` | historical build/design/audit records |

## Simulator (`sim/`)

Design notes behind the `rigel sim` simulator (usage lives in [`SIMULATOR.md`](SIMULATOR.md)).

| doc | what it is |
|---|---|
| [`sim/architecture_redesign.md`](sim/architecture_redesign.md) | simulator architecture |
| [`sim/hybrid_capture_simulation.md`](sim/hybrid_capture_simulation.md) | hybrid-capture modeling |
| [`sim/suite_driver_unification_plan.md`](sim/suite_driver_unification_plan.md) | suite-driver design |

## Performance

| doc | what it is |
|---|---|
| [`performance/README.md`](performance/README.md) | memory & compute optimization — a **historical snapshot (v0.6.2)**: profiling infrastructure + the ranked roadmap |

---

**Note on history.** The `docs/` tree has been pruned repeatedly of superseded planning documents (the
long calibration-evolution journey now lives under [`calibration/archive/`](calibration/archive/)). That
history remains in git. The docs above describe the **current** system;
[`calibration/CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) is the authoritative
reference when code and any older note disagree.
