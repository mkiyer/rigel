# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

Rigel is a Bayesian RNA-seq transcript quantification tool that jointly models mRNA, nascent RNA (nRNA), and genomic DNA contamination (gDNA). It uses a single-pass C++ BAM scanner, a per-region **calibration** stage that deconvolves the library into gDNA vs RNA, and a locus-level EM solver. Python package name is `rigel-rnaseq` on PyPI; the import and CLI are `rigel`.

> **Calibration (iterative odds-propagation simplex sweep — single production path).** The calibration stage deconvolves each contained region's unspliced mass into the 2-simplex **(f_rna₊, f_rna₋, f_g)** — sense-RNA / antisense-RNA / gDNA (calibration models **only RNA-vs-gDNA**; the per-locus EM separates nascent from mature downstream) — by an **odds-propagation grid sum-product** (`simplex_sweep.deconv_regions_sweep`), iterated over an **all-gDNA bootstrap**. Each node combines **three Bayesian terms**: a **strand likelihood** (the Beta-Binomial tilt of the per-strand counts — capture-invariant *direction*), a **count prior** (the local gDNA-density imputation on **RAW contained + STRAND-CLEANED boundaries**, weighted by its `var~mean` precision `τ_count` — *magnitude*), and a **global gDNA prior** (the foundation `ρ_global` at precision `τ_global`). It then **propagates** the per-strand RNA:gDNA log-odds along same-strand exon stretches so AMBIG (overlapping opposite-strand) nodes are resolved from their stranded neighbours. **Iteration:** Pass 0 = all-gDNA init (every unspliced fragment is gDNA ⇒ `ρ_global` = the count-observable total density, a deliberate over-estimate); each pass re-fits the gDNA `var~mean` (monotone P-spline, `variance_model.MonotoneVarMean`) + `ρ_global` on the previous pass's gDNA estimate, then re-solves, converging on per-node `f_g` (`config.sweep_n_grid`=60 lattice, `config.sweep_max_passes`=4). The boundary **sides** are deconvolved **once** (the per-node strand/count combine `g = w·g_strand + (1−w)·g_count`, `w = I/(I+I₀)`, `I = N·(2κ−1)²`, `I₀ = 10` = `config.gdna_strand_info_scale`, in `strand_deconv.deconv_sides`) as the fixed boundary gDNA anchors for the var~mean's boundary→region imputation, and feed the per-locus prior via the boundary-flux transport. The legacy per-node **fusion** region combine (`strand_deconv.deconv_regions`, `use_propagation` flag) was **deleted** in the Phase-4 teardown; the pre-teardown dual-path state is preserved at git tag `pre-fusion-teardown`. **The implementation is still being completed** — known open issues (count bias at AMBIG, iteration tightening, geom2 inflation) are tracked for fixing within this single framework. Theory: `docs/calibration/CALIBRATION_PLAN_v2.md` (authoritative plan) + `calibration_theory.md`; the iterative design in `docs/calibration/iterative_bootstrap_design.md` / `per_node_deconv_hierarchy_design.md` / `stage2_wiring_dryrun.md`; fractional-accumulator spec in `docs/accumulator/00_design.md`; docs index in `docs/README.md`. The full pipeline (scan → calibrate → **quant**) runs end-to-end: `quant_from_buffer` + `calibration.priors.assemble_priors` are wired in `run_pipeline`.

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

2. **Calibration** (`calibration.calibrate`): an **acyclic single-pass deconvolution** over the accumulator payload (no iterative EM/M-step). It deconvolves each region's fragment mass into **gDNA vs RNA** and fits the library hyperparameters — `ρ_0` (gDNA density), `ρ_d_bb`/`ρ_r_bb` (gDNA/RNA strand Beta-Binomial dispersions), `κ_rna` (RNA sense fraction) — plus the deconvolved gDNA/RNA mass. Output: `CalibrationResult`.

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

- `calibrate.py` — the calibrator orchestrator (RNA strand balance κ → gDNA/RNA strand overdispersion → `strand_deconvolve` cleans the boundary crossings → count prior `g_count` = imputation on raw contained + cleaned boundaries + 3-term splice → `deconv_sides` ONCE for the boundary anchors → **iterative `deconv_regions_sweep`**: per pass re-fit `ρ_global` + the gDNA `var~mean` on the running gDNA estimate, solve the per-node 3-term simplex + odds propagation, converge on `f_g` → derive). All-gDNA init; `sweep_n_grid`/`sweep_max_passes`
- `substrate.py` — `CalibrationSubstrate`: payload → per-region 3-view (contained / left / right) sufficient statistics
- `region_arrays.py` / `regions.py` — region geometry (`RegionArrays`, build partition from `index.region_df`)
- `signature.py` — region 4-bit strand/type signature + `strand_class` (POS/NEG/NONE/AMBIG)
- `density_model.py` — **count module**: per-region gDNA **density** via local boundary-anchored imputation + count-observability masks (`region_count_observable` / `boundary_count_observable`); returns `count_gdna_frac`. Operates on whatever per-view counts it is handed (`gdna_counts`: `None` ⇒ raw `pos+neg`; the calibrator passes **RAW contained + STRAND-CLEANED boundaries**, so exon/AMBIG imputation drops the nascent the count clue can't see)
- `strand_balance.py` / `strand_summary.py` — **RNA** strand *mean*: `rna_sense_frac` (used by the decode). `StrandBalance.rna_strand_overdispersion` here is a QC-only thin-count power diagnostic (`1/(n_obs+3)`), distinct from the decode's RNA overdispersion (see `gdna_strand.py`)
- `gdna_strand.py` — **both** strand Beta-Binomial overdispersions (shared component-agnostic MoM core): `gdna_strand_overdispersion` (mean ½, fit from count-observable seed regions + boundary sides) and `rna_strand_overdispersion` (mean κ, fit from boundary-side spliced counts). Both applied symmetrically in `strand_likelihood` with the same default prior, so unstranded data is uninformative (see `docs/calibration/calibration_theory.md` §4.3)
- `strand_deconv.py` — **strand module** (`strand_deconvolve`: per node the Beta-Binomial posterior gDNA fraction `g_strand` read at the FP-quantile + its information `I=N·(2κ−1)²`; weak Beta(½,½) prior) + the **count cleaning** `cleaned_gdna_count` (removes the strand-identified RNA from a count by `(w·g_strand+(1−w))`, `w=I/(I+I₀)`; applied to the boundary crossings, degrades to a no-op at low `I`) + the **per-node strand/count combine** `deconv_sides` (`g = w·g_strand + (1−w)·g_count`, `w=I/(I+I₀)`; via `_deconv_per_node`) used for the boundary anchors; strand-unobservable nodes (κ=½ / AMBIG) are count-only (`w=0`); FP-rate quantile `gdna_deconv_quantile` applied here (no-op at ½). Also exposes `boundary_side_seeds` for the gDNA strand fit. (The fusion region combine `deconv_regions` was removed in the Phase-4 teardown — regions now go through `simplex_sweep`.)
- `derive.py` — `gdna_density_global` (QC scalar) + the geometric gDNA length from the deconvolved masses
- `effective_length.py` — FL-marginal effective lengths (region / boundary)
- `capture_eff_length.py` — capture-aware **EM** effective lengths (`transcript_capture_eff_lengths`): contract each transcript's eff-len by the per-region gDNA-enrichment IPR over its region set, with an evidence-weighted shrinkage `w=G/(G+n_reg)` toward no-contraction on sparse gDNA
- `splice_junction.py` — `region_splice_gdna_frac`: convert the boundary density fraction `f_b` to the region count fraction via the gDNA/RNA region eff-length ratio (FL-consistency of the splice partition)
- `simplex.py` — per-node 2-simplex primitives (`_simplex_lattice`, `_mixture_strand_loglik`) shared by `simplex_sweep`
- `simplex_sweep.py` — **the production region deconv** (`deconv_regions_sweep`): odds-propagation grid sum-product on the per-reference region chain — per node a 3-term lattice belief (strand likelihood + count prior `τ_count` + global gDNA prior `τ_global`, with a node-class prior: Jeffreys at single-strand, global at AMBIG) propagated along same-strand exon stretches; called iteratively by `calibrate.py`
- `variance_model.py` — **gDNA `var~mean`** (`MonotoneVarMean`: monotone P-spline / SCAM, GCV-λ, IRLS-robust, power-law fallback) + `varmean_points` / `fit_direct_varmean` / `fit_imputation_varmean` builders (DIRECT own-count vs IMPUTATION boundary→region prediction-error, split by `region_count_observable`); yields the sweep's per-node `τ_count`
- `mature_density.py` / `rna_variance.py` — RNA density + splice-junction-pair var~mean model (Phase 2a)
- `run_fill.py` — run/neighbour helpers (`same_ref_left_right`) for the boundary-side and chain geometry
- `fl.py` — gDNA / RNA / global fragment-length pmfs (empirical-Bayes smoothed)
- `result.py` — `CalibrationResult` schema + intrinsic invariants
- `priors.py` — `assemble_priors`: `CalibrationResult` → per-locus Dirichlet prior (+ the gDNA-component IPR effective length)
- `errors.py` — calibration exceptions

### C++ Extensions (`src/rigel/native/`)

Five nanobind modules compiled via CMakeLists.txt:

| Module | Source | Purpose |
|--------|--------|---------|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | Single-pass htslib BAM parser, fragment grouping, model training, **fractional accumulator** (per-region/boundary mass + FL pools) |
| `_em_impl` | `em_solver.cpp` | Per-locus EM (`n_t + 1` components), connected components, partition scatter helpers, OpenMP |
| `_scoring_impl` | `scoring.cpp` | Fragment likelihood scoring, bias correction, SIMD optimized |
| `_resolve_impl` | `resolve.cpp` | Fragment-to-transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges/cgranges_bind.cpp` | Vendored cgranges interval overlap library |

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
