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

The gDNA-vs-RNA deconvolution stage — the most active area.

| doc | what it is |
|---|---|
| [`calibration_theory.md`](calibration/calibration_theory.md) | **authoritative** theory + acyclic inference + EM interface |
| [`CALIBRATION_TODO.md`](calibration/CALIBRATION_TODO.md) | live tracker of open issues (#1 strand-clean robustness, #2 count overdispersion, #3 capture leak, …) |
| [`count_overdispersion_integration_plan.md`](calibration/count_overdispersion_integration_plan.md) | count-overdispersion design (shipped; Issue #2) |
| [`strand_cleaning_robustness_design.md`](calibration/strand_cleaning_robustness_design.md) | Issue #1 options analysis (precision-weighted shrinkage) |
| [`strand_clean_global_target_design.md`](calibration/strand_clean_global_target_design.md) | Issue #1 proposed fix: shrink to the global gDNA fraction |
| [`density_phase1_local_imputation_design.md`](calibration/density_phase1_local_imputation_design.md) | count-clue density design (shipped) |
| [`density_phase2_dna_fraction_design.md`](calibration/density_phase2_dna_fraction_design.md) | deferred DNA-fraction lever (future) |
| [`accumulator_fragment_span_redesign.md`](calibration/accumulator_fragment_span_redesign.md) | gDNA FL span fix (shipped design record) |

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
