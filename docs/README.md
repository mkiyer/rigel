# Rigel documentation

Rigel is a Bayesian RNA-seq quantifier that jointly models mRNA, nascent RNA (nRNA), and genomic
DNA (gDNA) contamination. This is the documentation map.

## Start here

| doc | what it is |
|---|---|
| [`MANUAL.md`](MANUAL.md) | user manual — install, `index`/`quant`/`sim` CLI, options, outputs |
| [`METHODS.md`](METHODS.md) | the scientific method (the mRNA/nRNA/gDNA model + the three stages) |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | the code map and data flow (scan → calibrate → quant; Python modules + C++ extensions) |
| [`parameters.md`](parameters.md) | configuration parameters reference |

## Calibration (`calibration/`)

The gDNA-vs-RNA deconvolution stage — the most active area. **Consolidated 2026-06-13** to a single
authoritative plan + the theory/diagnostic references; ~28 superseded planning docs moved to
[`calibration/archive/`](calibration/archive/) (in git for history).

| doc | what it is |
|---|---|
| [`CALIBRATION_PLAN.md`](calibration/CALIBRATION_PLAN.md) | **★ SINGLE SOURCE OF TRUTH** — current state → target architecture → phased path to production. The per-node inverse-variance fusion (`w=I/(I+β)`); the target = **belief propagation on a per-locus region chain with 2-simplex grid messages** (edge-factor boundaries, component-specific process noise `Q_c`); all resolved decisions (issues A–H, Q1–Q3); risks (`σ²_RNA` unsolved, `γ_ij` unbuildable→class-chains, perf chunking, the +8.7pt regression bar); code + docs consolidation actions |
| [`calibration_theory.md`](calibration/calibration_theory.md) | **authoritative theory** — the acyclic single-pass pipeline, the per-node model, the EM interface |
| [`deconvolution_generative_model.md`](calibration/deconvolution_generative_model.md) | theory: the three-species (gDNA/nascent/mature) generative model + what is identifiable |
| [`deconvolution_implementation.md`](calibration/deconvolution_implementation.md) | theory→code: the capture-aware blueprint (boundaries as the atoms) — background for the sweep |
| [`strand_deconvolution_explained.md`](calibration/strand_deconvolution_explained.md) | ground-up explainer: the strand deconvolution, `Beta(½,½)` prior, why precision = `N·(2κ−1)²` |
| [`fl_consistency_diagnostic.md`](calibration/fl_consistency_diagnostic.md) | accuracy diagnostic: the splice density-vs-count eff-length-ratio fix (quantified) |
| [`capture_effective_length_design.md`](calibration/capture_effective_length_design.md) | capture-aware effective lengths for the EM components (the gDNA IPR contraction) |
| [`accumulator_fragment_span_redesign.md`](calibration/accumulator_fragment_span_redesign.md) | gDNA FL span fix (shipped design record) |
| [`CALIBRATION_TODO.md`](calibration/CALIBRATION_TODO.md) | live tracker of open issues |
| [`archive/`](calibration/archive/) | ~29 superseded planning/design docs (count-cascade, simplex/sweep/message-passing plans, count-trust, the phase A–C plans, the redesign series, …) — history, not current |

## Accumulator (`accumulator/`)

The C++ fractional accumulator (calibration sufficient statistics).

| doc | what it is |
|---|---|
| [`00_design.md`](accumulator/00_design.md) | **the spec** (per-region/boundary mass + FL pools; matched byte-for-byte by `tests/native/_accumulator_reference.py`) |
| `01_implementation.md`, `audit_phase1.md` | historical build/audit records |

## Simulator (`sim/`) + benchmarking

| doc | what it is |
|---|---|
| [`SIMULATOR.md`](SIMULATOR.md) | the synthetic read simulator (`rigel sim`) |
| [`BENCHMARKING.md`](BENCHMARKING.md) | the synthetic benchmark suite + net-fragment-flow evaluation |
| [`sim/architecture_redesign.md`](sim/architecture_redesign.md) | simulator architecture |
| [`sim/hybrid_capture_simulation.md`](sim/hybrid_capture_simulation.md) | hybrid-capture modeling |
| [`sim/suite_driver_unification_plan.md`](sim/suite_driver_unification_plan.md) | suite driver design |

## Release

| doc | what it is |
|---|---|
| [`PUBLISHING.md`](PUBLISHING.md) | building + publishing `rigel-rnaseq` to PyPI |

---

**Note on history.** In June 2026 the `docs/` tree was pruned of ~200 superseded planning documents
(the long calibration-evolution journey: `archive/`, `benchmarks/`, `acc_caljointmodel/`,
`futureprs/`, `em_strand/`, `mappability/`). That history remains in git. The docs above describe
the **current** system; `calibration/calibration_theory.md` is the authoritative reference when code
and any older note disagree.
