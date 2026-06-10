# Rigel architecture

The current code map and data flow. For the scientific method see `METHODS.md`; for calibration
theory see `calibration/calibration_theory.md`; for usage see `MANUAL.md`.

Rigel is a Bayesian RNA-seq quantifier that jointly models **mRNA**, **nascent RNA (nRNA)**, and
**genomic DNA (gDNA) contamination**. A single-pass C++ BAM scanner feeds a per-region
**calibration** stage (deconvolve the library into gDNA vs RNA) and a locus-level **EM** solver.
PyPI package `rigel-rnaseq`; import and CLI are `rigel`.

## Pipeline (`pipeline.py`)

Three stages, one BAM pass:

```
BAM ──scan_and_buffer──▶ FragmentBuffer + AccumulatorPayload + trained models
                              │
                              ├── calibrate ──▶ CalibrationResult (per-region gDNA/RNA mass + scalars)
                              │
                              └── quant_from_buffer ──▶ AbundanceEstimator (per-locus EM)
```

1. **Scan** (`scan_and_buffer`): the C++ htslib reader resolves fragments against the reference
   index, trains the strand and fragment-length models from unique mappers, buffers resolved
   fragments into a columnar `FragmentBuffer`, **and** deposits per-region/per-boundary fractional
   fragment mass into the C++ **accumulator** (4 channels: unspliced ±, spliced sense/antisense)
   → an `AccumulatorPayload`.

2. **Calibrate** (`calibration.calibrate`): an **acyclic single-pass** count×strand deconvolution of
   the accumulator payload into gDNA vs RNA per node, fitting the library hyperparameters
   (`gdna_density_global`, `rna_sense_frac`, the strand and count overdispersions). Output:
   `CalibrationResult`. See `calibration/calibration_theory.md`.

3. **Quantify** (`quant_from_buffer`): bridges calibration → per-locus prior
   (`calibration.priors.assemble_priors`), scores fragments (`scan.FragmentRouter` → CSR
   `ScoredFragments`), builds loci by connected components (`locus.build_multi_loci`), partitions
   the CSR per locus, and runs the per-locus EM. Each locus is an independent subproblem with
   **`n_t + 1` components** (one per transcript row + one gDNA). The calibration prior enters as two
   per-locus Dirichlet scalars (`gdna_prior_count`, `rna_prior_count`).

## Python modules

| module | role |
|---|---|
| `cli.py` | CLI: `index`, `quant`, `sim`, `export` |
| `pipeline.py` | orchestrator: scan → calibrate → quant |
| `config.py` | frozen config dataclasses (`EMConfig`, `BamScanConfig`, `PipelineConfig`, `CalibrationConfig`, `FragmentScoringConfig`) |
| `index.py` | reference index build/load (feather from GTF+FASTA); `region_df` = the calibration region partition |
| `scoring.py` | fragment likelihood scoring (RNA/gDNA FL models, strand/coverage/splice penalties) |
| `scan.py` | `FragmentRouter` — scores the buffer → global CSR `ScoredFragments` |
| `scored_fragments.py` / `locus_partition.py` | CSR fragment container + per-locus partitioning |
| `locus.py` | locus / `MultiLocus` construction (connected components), nRNA spans |
| `estimator.py` | `AbundanceEstimator` — per-locus EM dispatch, output dataframes |
| `scan_payload.py` / `_accumulator.py` | `AccumulatorPayload` schema + Python `Accumulator` wrapper |
| `buffer.py` | memory-efficient columnar fragment buffer |
| `strand_model.py` / `frag_length_model.py` | model training from unique mappers |
| `splice.py` / `splice_blacklist.py` | splice-type encoding + artifact blacklist |
| `annotate.py` | optional annotated-BAM writer (per-fragment tags ZT/ZG/ZF/ZC/ZL/ZW) |
| `native.py` | interface to the C++ extension modules |

### Calibration package (`calibration/`)

| module | role |
|---|---|
| `calibrate.py` | the acyclic calibrator orchestrator |
| `substrate.py` | `CalibrationSubstrate`: payload → per-region 3-view (contained/left/right) sufficient stats |
| `region_arrays.py` / `regions.py` | region geometry (`RegionArrays`, partition from `index.region_df`) |
| `signature.py` | 4-bit exon/intron×strand signature + `strand_class` (POS/NEG/NONE/AMBIG) |
| `strand_balance.py` | RNA strand mean `rna_sense_frac` (κ) from the spliced channel |
| `density_model.py` | count clue: per-region gDNA density via local imputation + strand cleaning |
| `gdna_strand.py` | gDNA & RNA strand Beta-Binomial overdispersions (shared pooled-MoM core) |
| `count_dispersion.py` | gDNA count overdispersion (NB MoM); concentration `N_eff = N/(1+αN)` |
| `joint_deconv.py` | per-node joint count×strand deconvolution into gDNA/RNA |
| `derive.py` | `gdna_density_global` + geometric gDNA length from the deconvolved masses |
| `effective_length.py` / `fl.py` | FL-marginal effective lengths; gDNA/RNA/global FL pmfs |
| `priors.py` | `assemble_priors`: `CalibrationResult` → per-locus Dirichlet prior + gDNA eff-len |
| `result.py` / `errors.py` | `CalibrationResult` schema + invariants; exceptions |

## C++ extensions (`src/rigel/native/`)

Five nanobind modules (C++17, `-O3`, LTO; SIMD `fast_exp.h`; OpenMP):

| module | source | purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | single-pass htslib parser, fragment grouping, model training, **fractional accumulator** |
| `_em_impl` | `em_solver.cpp` | per-locus EM (`n_t+1` components), connected components, OpenMP, Kahan/digamma |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring, bias correction, SIMD (`-ffast-math`) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges_bind.cpp` | vendored cgranges interval overlap |

The accumulator must match the Python reference (`tests/native/_accumulator_reference.py`)
byte-for-byte. After editing `src/rigel/native/`, re-run `pip install --no-build-isolation -e .`
in the activated `rigel` conda env.

## nRNA architecture

nRNA components are not per-transcript shadows. Rigel builds a global table of unique nRNA spans
keyed by `(ref, strand, start, end)` and shares each span across transcripts with matching
coordinates. Each unique span is materialized as an ordinary **transcript row** in `index.t_df`
(flagged `is_nrna`/`is_synthetic`), so the per-locus EM treats it like any transcript component —
hence `n_t + 1` components per locus, not a separate mRNA/nRNA pair per transcript.

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim [options]
```

Input BAM must be name-sorted with the `NH` tag.
