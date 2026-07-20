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

The gDNA-vs-RNA deconvolution stage: a bipartite region↔boundary **belief-propagation sweep** over the
accumulator payload. **Consolidated 2026-07-18 to a single master reference + a small set of anchors; the
per-phase roadmaps, message-layer WIP, and investigation notes are now under [`archive/`](calibration/archive/).**

**Read these — the current, converged set:**

| doc | what it is |
|---|---|
| [`CALIBRATION_MASTER.md`](calibration/CALIBRATION_MASTER.md) | **★ START HERE** — the single master reference: the goal, the pipeline in order, the two NPMLE roles (enrichment/transfer-variance vs the gDNA hyperprior), the AMBIG two-root problem, the two-solver reconciliation, what's implemented vs. what needs restructuring, and the phased plan P1–P5. |
| [`CALIBRATION_ARCHITECTURE.md`](calibration/CALIBRATION_ARCHITECTURE.md) | **★ AUTHORITATIVE INVARIANT** — the count-zero-information principle and the three composition sources. If any note disagrees about the *invariant*, this wins; `CALIBRATION_MASTER.md` is authoritative about the *design and plan*. |
| [`reference_prior_derivation.md`](calibration/reference_prior_derivation.md) | reviewer-facing — the **resolved** Beta(½,½) reference (§10) and the exact 1-D/2-D solver reconciliation (θ=arcsin(τ) cancels the tilt; class gap 0). |
| [`background_reference_derivation.md`](calibration/background_reference_derivation.md) | reviewer-facing — the aggregate background DNA level `ρ_bg` (the data-gated global hyperprior). |
| [`transfer_variance_formal_derivation.md`](calibration/transfer_variance_formal_derivation.md) | reviewer-facing — `σ²_transfer` message precision via projection onto the enrichment NPMLE (the NPMLE's Role A). |
| [`density_imputation_precision.md`](calibration/density_imputation_precision.md) | reviewer-facing (**derivation + proposal, for prototyping/review**) — the message precision for imputing a **density** across nodes that differ in total density, effective length, and composition: the 3-term variance (per-component count `1/ρE` + enrichment gap `σ²_transfer` + coverage non-uniformity). |
| [`cliff_message_derivation.md`](calibration/cliff_message_derivation.md) | reviewer-facing — the cross-node message **mode** across the enrichment cliff (the eff-length-frame log-odds shift); §9 records the rejected "hybrid enrichment-corrected density" unifiers and why the exon floor is identifiability, not the mode. |
| [`oracle_and_benchmarking.md`](calibration/oracle_and_benchmarking.md) | the one trusted debug/benchmark path — the accumulator-consistent oracle (`scripts/debug/oracle.py`) and the net-fragment-flow metric. |
| [`archive/`](calibration/archive/) | superseded roadmaps/plans, the message-layer WIP series (relay derivation, mature-gate, spliced channel), and per-session investigation notes — **provenance, not current**. Durable findings are folded into `CALIBRATION_MASTER.md` §5–§9. |

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
