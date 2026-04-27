# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Summary

Rigel is a Bayesian RNA-seq transcript quantification tool that jointly models mRNA, nascent RNA (nRNA), and genomic DNA contamination (gDNA). It uses a single-pass C++ BAM scanner plus a locus-level EM solver. Python package name is `rigel-rnaseq` on PyPI; the import and CLI are `rigel`.

## Build & Development

```bash
# Environment setup (requires conda/mamba with conda-forge + bioconda channels)
mamba env create -f mamba_env.yaml
conda activate rigel

# Build and install (editable, with C++ compilation via scikit-build-core + nanobind)
pip install --no-build-isolation -e .

# Install with dev extras (pytest)
pip install --no-build-isolation -ve ".[dev]"
```

After changing any C++ code in `src/rigel/native/`, you must re-run `pip install --no-build-isolation -e .` to recompile.

## Testing

```bash
pytest tests/ -v                          # all tests
pytest tests/test_em_impl.py -v           # single test file
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

### Two-Stage Pipeline (`pipeline.py`)

1. **BAM Scan** (`scan_and_buffer`): C++ htslib-based single-pass BAM reader. Resolves fragments against the reference index, trains strand/fragment-length models from unique mappers, and buffers results into a columnar `FragmentBuffer`.

2. **Quantification** (`quant_from_buffer`): Runs Simple Regional Deconvolution (SRD) calibration on the FragmentBuffer to recover the library-wide gDNA fraction `π_pool` and the gDNA fragment-length model. Iterates the buffer to build CSR EM data (`scan.FragmentRouter`), constructs loci via connected components (`locus.build_loci`), derives per-locus Dirichlet priors from per-fragment SRD posteriors (`compute_locus_priors_from_partitions`), and runs per-locus EM with `2*n_t + 1` components (mRNA + nRNA per transcript, plus one gDNA).

### Python Module Roles

- `cli.py` — CLI entry point with subcommands: `index`, `quant`, `sim`
- `pipeline.py` — Thin orchestrator connecting the two stages
- `config.py` — Frozen dataclasses (`EMConfig`, `BamScanConfig`, `PipelineConfig`)
- `index.py` — Reference index build/load (produces feather files from GTF+FASTA)
- `scoring.py` — Fragment likelihood scoring (`FragmentScorer`)
- `locus.py` — Locus construction, nRNA initialization, per-locus Dirichlet priors
- `calibration/_simple.py` — SRD orchestrator (`calibrate_gdna`)
- `calibration/_categorize.py` — Pass 0 per-fragment 7-way categorization
- `calibration/_fl_mixture.py` — Pass 2 1-D fragment-length mixture EM
- `calibration/_fl_empirical_bayes.py` — Pass 3 empirical-Bayes Dirichlet smoothing
- `calibration/_result.py` — `CalibrationResult` schema
- `estimator.py` — `AbundanceEstimator` class, EM dispatch
- `scan.py` — `FragmentRouter` (CSR builder from scored fragments)
- `buffer.py` — Memory-efficient columnar fragment buffer with CSR layout
- `strand_model.py` / `frag_length_model.py` — Model training from unique mappers
- `priors.py` — Empirical Bayes prior computation (kappa estimation)
- `native.py` — Public interface to all five C++ extension modules

### C++ Extensions (`src/rigel/native/`)

Five nanobind modules compiled via CMakeLists.txt:

| Module | Source | Purpose |
|--------|--------|---------|
| `_bam_impl` | `bam_scanner.cpp` | Single-pass htslib BAM parser, fragment grouping, model training |
| `_em_impl` | `em_solver.cpp` | Equivalence class EM with SQUAREM acceleration, OpenMP parallelism |
| `_scoring_impl` | `scoring.cpp` | Fragment likelihood scoring, bias correction, SIMD optimized |
| `_resolve_impl` | `resolve.cpp` | Fragment-to-transcript resolution via cgranges interval tree |
| `_cgranges_impl` | `cgranges_bind.cpp` | Vendored cgranges interval overlap library |

Key C++ details:
- C++17, compiled with `-O3`, LTO enabled
- SIMD: `fast_exp.h` provides AVX2/AVX-512/NEON code paths for exp()
- EM solver uses Kahan compensated summation and digamma for VBEM mode
- OpenMP for parallel BAM scanning and batch locus EM
- `_scoring_impl` uses `-ffast-math`; `_em_impl` uses `-ffp-contract=fast` (preserves Kahan)

### nRNA Architecture

nRNA components are not per-transcript shadows. Instead, Rigel builds a global table of unique nRNA spans keyed by `(ref, strand, start, end)` and shares each span across transcripts with matching genomic coordinates. This reduces redundant states in isoform-rich loci.

## Test Infrastructure

- 28 test modules in `tests/` with fixtures in `conftest.py`
- `conftest.py` provides synthetic mini GTF/FASTA fixtures (3 transcripts, 2 genes)
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
