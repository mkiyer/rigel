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
| [`decoupled_calibration_design.md`](calibration/decoupled_calibration_design.md) | **current architecture**: the decoupled strand/count handoff (design + implementation guide) |
| [`deconvolution_roadmap.md`](calibration/deconvolution_roadmap.md) | **ACTIVE ROADMAP** (three-signal complex-locus resolution): unspliced-crossing gDNA + spliced-crossing per-strand mature + strand-on-residual; mature subtraction; phased A–D plan. Operationalizes `deconvolution_generative_model.md` + `deconvolution_implementation.md` |
| [`deconvolution_generative_model.md`](calibration/deconvolution_generative_model.md) | theory: the complete three-species (gDNA/nascent/mature) generative model + what is identifiable |
| [`deconvolution_implementation.md`](calibration/deconvolution_implementation.md) | theory→code: the capture-aware S0–S5 blueprint (boundaries as the atoms) |
| [`phaseA_mature_imputation_plan.md`](calibration/phaseA_mature_imputation_plan.md) | **Phase A (shipped 83b3610)**: per-strand mature-density imputation (`mature_density.py`, RNA mirror of `density_model`); Issue D resolved |
| [`phaseB_mature_subtraction_plan.md`](calibration/phaseB_mature_subtraction_plan.md) | **Phase B (Step 0 done; merged into C)**: mature-subtraction units verified; the finding that subtraction only helps AMBIG |
| [`phaseC_ambig_resolution_plan.md`](calibration/phaseC_ambig_resolution_plan.md) | **Phase C — ACTIVE plan, rev 3**: the iterative **propagation** deconvolution — a graph of region+boundary nodes, each a 3-term {RNA+,RNA−,gDNA} partition, solved from seeds by propagation to convergence (boundaries deconvolve their unspliced crossings into all 3 terms too); consolidates density_model+mature_density+strand-blend+run_fill, retires the splice fraction / `w` blend / `ρ` |
| [`archive/joint_deconvolution.md`](calibration/archive/joint_deconvolution.md) | archived: the retired joint count×strand product (resurrection guide) |
| [`CALIBRATION_TODO.md`](calibration/CALIBRATION_TODO.md) | live tracker of open issues |
| [`remaining_phases.md`](calibration/remaining_phases.md) | **AUTHORITATIVE forward plan** — the remaining work as clean Phase 1 (count posterior) / 2 (FP quantile) / 3 (nascent removal), implementation-grade |
| [`count_channel_capture_design.md`](calibration/count_channel_capture_design.md) | background: the original count-channel roadmap (superseded by `remaining_phases.md` for sequencing) |
| [`count_mean_bias_design.md`](calibration/count_mean_bias_design.md) | **Phase 4-mean**: the splice-junction gDNA-fraction (eligibility predicate, 2-/3-term, strand-resolved sweep) |
| [`count_posterior_design.md`](calibration/count_posterior_design.md) | **Phase 4-var**: count posterior variance (`var∝mean²` from paired anchor disagreements) |
| [`phase2_design.md`](calibration/phase2_design.md) | **Phase 2** (shipped + validated): count posterior variance behavior, the blend decision (keep (2κ−1)²), the FP-rate quantile knob |
| [`sequential_calibration_redesign.md`](calibration/sequential_calibration_redesign.md) | **AUTHORITATIVE redesign** (current): sequential strand→count pipeline — strand cleans + communicates likelihood precision (`N·(2ss−1)²`), count imputes the cleaned density field + 3-term; phased plan |
| [`strand_deconvolution_explained.md`](calibration/strand_deconvolution_explained.md) | ground-up explainer: the strand deconvolution, the `Beta(½,½)` Jeffreys prior, why precision = `N·(2ss−1)²`, the `ss=½` bimodal subtlety |
| [`phase3_implementation_plan.md`](calibration/phase3_implementation_plan.md) | **SUPERSEDED** by the redesign — the direct-3-term Phase 3 (kept for 3-term/AMBIG-identifiability rationale) |
| [`fl_consistency_diagnostic.md`](calibration/fl_consistency_diagnostic.md) | **accuracy diagnostic** (Phase-3 prep): the splice `f_b·M_region` density-vs-count bias — exact eff-length-ratio fix, quantified (16–25pt at typical exons, present in the benchmark), density-path already consistent |
| [`capture_effective_length_design.md`](calibration/capture_effective_length_design.md) | **capture-aware effective lengths** for all EM components (the gDNA-only IPR contraction bug + fix) |
| [`count_local_dispersion_design.md`](calibration/count_local_dispersion_design.md) | parked: local count-dispersion estimator (Phase 4-var extension) |
| [`strand_clean_robust_deferred.md`](calibration/strand_clean_robust_deferred.md) | deferred: precision-weighted robust strand-clean concept (subsumed by the strand module) |
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
