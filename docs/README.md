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

The gDNA-vs-RNA deconvolution stage. **Mid-migration between two solvers; NOT ready to ship.**
**Re-consolidated 2026-07-24: 51 more stale docs moved to [`archive/`](calibration/archive/).** The live set
is small; **do not reference `archive/`.**

**Read these — the current set:**

| doc | what it is |
|---|---|
| [`ROADMAP.md`](calibration/ROADMAP.md) | **★ START HERE** — dev status + goals: the pass-0 → gDNA-hyperprior → re-solve pipeline, what is landed vs open, the two standing DEBTS (`ω_graft` and P1e), and the ordered path to production. |
| [`CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) | **★ AUTHORITATIVE INVARIANT** — the count-zero-information principle and the three composition sources (strand + imputation + population prior). Still correct; read second. |
| [`variance_model_handoff.md`](calibration/variance_model_handoff.md) | **the variance-model derivation, to be REDONE** — why the composition solver breaks the old density-uniformity variance model; the per-component form to derive; the peel; the transfer-variance redo; the ω decision. |
| [`unified_solver_design.md`](calibration/unified_solver_design.md) | the target solver's **mode** design (the reframe + ÷M_dst). **Its precision/variance sections are superseded** by the handoff above. |
| [`gdna_intron_factory_design.md`](calibration/gdna_intron_factory_design.md) | a **shipped** feature — the intron gDNA factory (`intron_factory=True` default). |
| [`SESSION_2026_07_26_HANDOFF_14.md`](calibration/SESSION_2026_07_26_HANDOFF_14.md) | the LIVE session handoff — state, what landed, the do-not-re-run table, and how to run every gate. Earlier `SESSION_*_HANDOFF_*.md` files are superseded history. |
| [`archive/`](calibration/archive/) | **provenance, not current** — all superseded roadmaps, derivations, and per-session notes. Do not cite. |

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
