# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

Rigel is a Bayesian RNA-seq transcript quantification tool that jointly models mRNA, nascent RNA (nRNA), and genomic DNA contamination (gDNA). It uses a single-pass C++ BAM scanner, a per-region **calibration** stage that deconvolves the library into gDNA vs RNA, and a locus-level EM solver. Python package name is `rigel-rnaseq` on PyPI; the import and CLI are `rigel`.

> **Calibration (acyclic, calibration-v6).** The calibration stage is an **acyclic single-pass** count×strand deconvolution (no iterative EM/M-step — the old `density.py`/`estep.py`/`exposure.py`/`sweep.py`/`mstep.py` modules were removed in the acyclic redesign). Theory: `docs/caljointmodel/` (`01_generative_model.md`, `03_inference.md`, `04_interface_contract.md`); fractional-accumulator spec in `docs/accumulator/`; the gDNA strand-overdispersion line of work in `docs/em_strand/`. The full pipeline (scan → calibrate → **quant**) runs end-to-end: `quant_from_buffer` + `calibration.priors.assemble_priors` are wired in `run_pipeline`.

## Build & Development

> **ALL build, test, and lint commands MUST run inside the activated `rigel`
> conda environment.** It contains the full toolchain — htslib, C/C++
> compilers, `scikit-build-core`, `nanobind` — and the C++ build locates
> htslib via `$CONDA_PREFIX` (so the env must be active, not merely on `PATH`).
> In a non-interactive shell, prefix commands with:
> ```bash
> source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
> ```

```bash
# One-time environment setup (conda-forge + bioconda channels)
mamba env create -f mamba_env.yaml

conda activate rigel    # ← required before every build / test / lint

# Build and install (editable, C++ compiled via scikit-build-core + nanobind)
pip install --no-build-isolation -e .

# Install with dev extras (pytest)
pip install --no-build-isolation -ve ".[dev]"
```

After changing any C++ code in `src/rigel/native/`, you must re-run `pip install --no-build-isolation -e .` (in the activated `rigel` env) to recompile.

## Testing

```bash
pytest tests/ -v                          # all tests
pytest tests/calibration/ -v              # calibration unit tests
pytest tests/test_em_impl.py::test_name   # single test function
pytest tests/ --update-golden             # regenerate golden output files
pytest tests/ --cov=rigel                 # with coverage
```

Golden output regression tests live in `tests/golden/` (feather, TSV, JSON). Use `--update-golden` to regenerate after intentional output changes.

## Linting

Ruff is configured in `pyproject.toml`: Python 3.12 target, 100-char line length.

```bash
ruff check src/ tests/
ruff format src/ tests/
```

## Architecture

### Pipeline (`pipeline.py`)

1. **BAM Scan** (`scan_and_buffer`): C++ htslib single-pass BAM reader. Resolves fragments against the reference index, trains strand/fragment-length models from unique mappers, buffers resolved fragments into a columnar `FragmentBuffer`, **and** deposits per-region/per-boundary fractional fragment mass into the C++ **accumulator** (4 channels: unspliced ±, spliced sense/antisense) → an `AccumulatorPayload`.

2. **Calibration** (`calibration.calibrate`): a per-region EM over the accumulator payload. It deconvolves each region's fragment mass into **gDNA vs RNA** and fits the library hyperparameters — `ρ_0` (gDNA density), `φ` (NB count dispersion), `ε_s` (gDNA splice-artifact rate), `ρ_d_bb`/`ρ_r_bb` (gDNA/RNA strand Beta-Binomial dispersions), `κ_rna` (RNA sense fraction) — plus a per-region exposure posterior `ω` and deconvolved gDNA/RNA mass. Output: `CalibrationResult`.

3. **Quantification** (`quant_from_buffer`): bridges calibration → per-locus prior (`calibration.priors.assemble_priors`), scores fragments (`scan.FragmentRouter` → CSR `ScoredFragments`), builds loci via connected components (`locus.build_multi_loci`), partitions the CSR per locus (`locus_partition.partition_and_free`), and runs per-locus EM. Each locus is an independent subproblem with **`n_t + 1` components** — one per transcript row (annotated mRNA and nRNA spans alike) plus one gDNA component. The calibration prior enters as **two per-locus Dirichlet scalars** (`alpha_rna_add`, `alpha_gdna_add`) that set the gDNA-vs-RNA split; the EM distributes RNA mass among the compatible transcripts.

### Python Module Roles

- `cli.py` — CLI entry point with subcommands: `index`, `quant`, `sim`
- `pipeline.py` — orchestrator connecting scan → calibrate → quant
- `config.py` — frozen dataclasses (`EMConfig`, `BamScanConfig`, `PipelineConfig`, `CalibrationConfig`, `FragmentScoringConfig`)
- `index.py` — reference index build/load (feather files from GTF+FASTA); `region_df` defines the calibration region partition
- `scoring.py` — fragment likelihood scoring (`FragmentScorer`, RNA/gDNA FL models + strand/coverage/splice penalties)
- `scan.py` — `FragmentRouter` (scores the buffer → global CSR `ScoredFragments`)
- `scored_fragments.py` / `locus_partition.py` — CSR fragment container + per-locus partitioning
- `locus.py` — locus / `MultiLocus` construction (connected components), nRNA spans
- `estimator.py` — `AbundanceEstimator` (per-locus EM dispatch, output dataframes)
- `scan_payload.py` / `_accumulator.py` — `AccumulatorPayload` schema + Python `Accumulator` wrapper
- `buffer.py` — memory-efficient columnar fragment buffer
- `strand_model.py` / `frag_length_model.py` — model training from unique mappers
- `splice.py` / `splice_blacklist.py` — splice-type encoding + artifact blacklist
- `native.py` — public interface to the C++ extension modules

**Calibration package (`calibration/`):**

- `calibrate.py` — the **acyclic single-pass** calibrator orchestrator (RNA strand balance → count-clue density → gDNA/RNA strand overdispersion → count overdispersion → joint count×strand deconvolution → derive)
- `substrate.py` — `CalibrationSubstrate`: payload → per-region 3-view (contained / left / right) sufficient statistics
- `region_arrays.py` / `regions.py` — region geometry (`RegionArrays`, build partition from `index.region_df`)
- `signature.py` — region 4-bit strand/type signature + `strand_class` (POS/NEG/NONE/AMBIG)
- `density_model.py` — count clue: per-region/boundary gDNA **density** via the region↔boundary sweep + count-observability masks (`region_count_observable` / `boundary_count_observable`); strand-cleaned by the strand *mean*
- `strand_balance.py` / `strand_summary.py` — **RNA** strand *mean*: `rna_sense_frac` (used by the decode). `StrandBalance.rna_strand_overdispersion` here is a QC-only thin-count power diagnostic (`1/(n_obs+3)`), distinct from the decode's RNA overdispersion (see `gdna_strand.py`)
- `gdna_strand.py` — **both** strand Beta-Binomial overdispersions (shared component-agnostic MoM core): `gdna_strand_overdispersion` (mean ½, fit from count-observable seed regions + boundary sides) and `rna_strand_overdispersion` (mean κ, fit from boundary-side spliced counts). Both applied symmetrically in `strand_likelihood` with the same default prior, so unstranded data is uninformative (see `docs/em_strand/05`)
- `count_dispersion.py` — gDNA **count** overdispersion (count-side twin of `gdna_strand`): per-type NB method-of-moments `α` (contained regions vs crossing boundary sides) shrunk toward the global pooled-seed trend (α₀=1 geometric fallback). The count-prior concentration is the overdispersion-limited effective count `N_eff = N/(1+α·N)` — the honest precision that replaced the old categorical `defer_to_strand` zeroing (see `docs/calibration/count_overdispersion_integration_plan.md`)
- `joint_deconv.py` — per-node joint count×strand deconvolution into gDNA/RNA (`strand_loglik` two-component variance); also exposes `boundary_side_seeds` for the gDNA strand fit
- `derive.py` — `gdna_density_global` + per-node exposure (+ exposure-weighted gDNA length) from the deconvolved masses
- `effective_length.py` — FL-marginal effective lengths (region / boundary)
- `fl.py` — gDNA / RNA / global fragment-length pmfs (empirical-Bayes smoothed)
- `result.py` — `CalibrationResult` schema + intrinsic invariants
- `priors.py` — `assemble_priors`: `CalibrationResult` → per-locus Dirichlet prior
- `errors.py` — calibration exceptions

### C++ Extensions (`src/rigel/native/`)

Five nanobind modules compiled via CMakeLists.txt:

| Module | Source | Purpose |
|--------|--------|---------|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | Single-pass htslib BAM parser, fragment grouping, model training, **fractional accumulator** (per-region/boundary mass + FL pools) |
| `_em_impl` | `em_solver.cpp` | Per-locus EM (`n_t + 1` components), connected components, partition scatter helpers, OpenMP |
| `_scoring_impl` | `scoring.cpp` | Fragment likelihood scoring, bias correction, SIMD optimized |
| `_resolve_impl` | `resolve.cpp` | Fragment-to-transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges_bind.cpp` | Vendored cgranges interval overlap library |

Key C++ details:
- C++17, compiled with `-O3`, LTO enabled
- SIMD: `fast_exp.h` provides AVX2/AVX-512/NEON code paths for exp()
- EM solver uses Kahan compensated summation and digamma for VBEM mode
- OpenMP for parallel BAM scanning and batch locus EM
- `_scoring_impl` uses `-ffast-math`; `_em_impl` uses `-ffp-contract=fast` (preserves Kahan)
- The accumulator must match the Python reference (`tests/native/_accumulator_reference.py`) byte-for-byte

### nRNA Architecture

nRNA components are not per-transcript shadows. Rigel builds a global table of unique nRNA spans keyed by `(ref, strand, start, end)` and shares each span across transcripts with matching genomic coordinates, reducing redundant states in isoform-rich loci. Each unique span is materialized as an ordinary **transcript row** in `index.t_df` (flagged `is_nrna` / `is_synthetic`), so the per-locus EM treats it like any other transcript component — hence `n_t + 1` components per locus (all transcript rows + one gDNA), not a separate mRNA/nRNA pair per transcript.

## Test Infrastructure

- `tests/` (top level) + `tests/calibration/` (calibration units) + `tests/native/` (accumulator spec / byte-for-byte) — ~60 test modules; fixtures in `conftest.py` (synthetic mini GTF/FASTA: 3 transcripts, 2 genes)
- Golden output tests (`test_golden_output.py`) compare against `tests/golden/` reference files
- Scenario BAM files in `tests/scenarios/` and `tests/scenarios_aligned/`
- CI runs on Ubuntu + macOS, Python 3.12 + 3.13

## CLI Subcommands

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim [options]   # generate synthetic test scenarios
```

Input BAM must be name-sorted with `NH` tag.
