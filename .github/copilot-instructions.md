# Rigel Project Instructions

## Environment

- **Always** activate the conda environment before running any command: `conda activate rigel`
- Python 3.12+, C++17
- Dependencies managed via `mamba_env.yaml` (conda-forge + bioconda channels)

## Build

```bash
# After ANY C++ change in src/rigel/native/, you MUST recompile:
conda activate rigel && pip install --no-build-isolation -e .
```

Do not skip recompilation after editing `.cpp`, `.h`, or `CMakeLists.txt` files.

## Testing

```bash
conda activate rigel
pytest tests/ -v                          # all tests
pytest tests/test_em_impl.py -v           # single file
pytest tests/test_em_impl.py::test_name   # single test
pytest tests/ --update-golden             # regenerate golden outputs after intentional changes
```

- Known pre-existing failure: `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna` — unrelated to scoring/EM.
- Golden output regression files live in `tests/golden/`.
- 953+ tests should pass (excluding the known failure above).

## Linting

Ruff configured in `pyproject.toml`: Python 3.12 target, 100-char line length.

```bash
ruff check src/ tests/
ruff format src/ tests/
```

## Documentation

- Publish implementation plans to the `docs/` directory as markdown files. Use subfolders for organization. Avoid polluting the base `docs/` directory.

## Scripts

- Diagnostic exploration and debugging scripts are commonly needed. For organization, put them in `scripts/debug/` with description names. Avoid polluting the base `scripts/` directory.

## Simulations / Benchmarking / Profiling

- Create temporary scripts in `scripts/benchmark`, `scripts/profiling`, or `scripts/simulation` as needed. Try not to pollute the base `scripts/` directory.

## Benchmarking

`scripts/benchmark.py` contains a framework for running full-scale simulations, comparing multiple methods, and generating summary tables.

## Profiling

`scripts/profile.py` contains a framework for profiling

## Simulation

`scripts/synthetic_sim_sweep.py` contains a framework for running synthetic simulations with configurable parameters

## Architecture Summary

Two-stage pipeline: (1) C++ BAM scan → FragmentBuffer, (2) Python/C++ locus-level EM quantification.

- **C++ extensions** (`src/rigel/native/`): `_bam_impl`, `_em_impl`, `_scoring_impl`, `_resolve_impl`, `_cgranges_impl` — compiled via scikit-build-core + nanobind + CMake
- **Scoring model**: mRNA, nRNA, and gDNA candidates scored with strand, fragment-length, overhang, and NM mismatch penalties
- **EM solver**: Per-locus EM with `2*n_t + 1` components (mRNA + nRNA per transcript + one gDNA), SQUAREM acceleration, OpenMP parallel
- **Tripartite prior**: mRNA uses OVR prior, nRNA uses Dirichlet sparsity (α<1), gDNA uses Empirical Bayes anchor

## Key Constraints

- Mega-loci (large connected components from multimappers) are intentional and correct. Do NOT modify connected_components to exclude multimappers.
- `read_length` in resolve_context.h is the sum of aligned exon block spans (M/D/=/X), NOT the full query sequence length. Soft-clipped bases are excluded.
- Chimeric fragments are detected and buffered but **skipped** in scoring — they never enter EM.
- Input BAM must be name-sorted with `NH` tag (minimap2 doesn't set NH; the code handles this via secondary flag detection).

## C++ Compilation Details

- C++17, `-O3`, LTO enabled
- `_scoring_impl` uses `-ffast-math`; `_em_impl` uses `-ffp-contract=fast` (preserves Kahan summation)
- SIMD: `fast_exp.h` provides AVX2/AVX-512/NEON code paths
- multithreading for parallel BAM scanning and batch locus EM


