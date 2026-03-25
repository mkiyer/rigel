# Changelog

All notable changes to Rigel will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Fragment length histograms in summary.json**: the `fragment_length` section
  now includes full per-category histograms (trimmed zero bins) alongside the
  existing summary statistics. Each category (global, rna, gdna, intergenic,
  and per-SpliceType) reports `{summary: {...}, histogram: {range, values}}`.
  Enables downstream QC plotting and RNA-vs-gDNA distribution comparison.

- **Annotated BAM: transcript and gene index tags**: new integer tags `ZI`
  (transcript index) and `ZJ` (gene index) encode the zero-based index into
  the rigel reference, enabling compact downstream encodings without string
  lookups. Unassigned fragments get `-1`.

- **config.yaml output**: `rigel quant` writes a `config.yaml` to the output
  directory recording all resolved parameters. Rerun with
  `rigel quant --config results/config.yaml`.

### Changed

- **Removed `--confidence-threshold`**: the high-confidence EM count pathway
  has been fully removed from Python, C++, CLI, and config. The EM now focuses
  on primary assignment outputs only.

- **CLI cleanup**: default-first ordering for `--assignment-mode` and
  `--em-mode` choices; clearer help strings throughout; fixed duplicate
  `--no-*` flags from `BooleanOptionalAction`; quant I/O args can now be
  provided from YAML config.

## [0.3.2] - 2026-03-24

### Added

- **Run-reproducible config export**: `rigel quant` now writes `config.yaml` to
  the output directory with all resolved parameters and I/O paths, so runs can
  be reproduced via `rigel quant --config config.yaml`.

### Changed

- **CLI/config resolution cleanup**: quant I/O options (`--bam`, `--index`,
  `--output-dir`) can now be provided from YAML config when not set on CLI,
  with explicit validation and clearer error messages for missing required
  inputs.

- **EM output simplification**: removed the high-confidence EM count pathway
  (`confidence_threshold` and related accumulators) from Python/C++ EM plumbing,
  keeping transcript EM accounting focused on the primary assignment outputs.

- **Bioconda platform coverage**: recipe now opts into additional builds for
  `linux-aarch64` and `osx-arm64`.

### Fixed

- **Bioconda macOS packaging robustness**: conda recipe build flow was updated
  to use a dedicated `build.sh` that disables scikit-build stripping and forces
  classic macOS linking, addressing `llvm-otool`/Mach-O post-processing
  failures seen in cross-compiled macOS builds.

## [0.3.1] - 2026-03-24

### Fixed

- **Bioconda build failure on macOS x86_64**: `pkg_check_modules(HTSLIB IMPORTED_TARGET htslib)`
  caused CMake to fail when pkg-config could not resolve `bzip2` as a transitive htslib
  dependency inside the conda-build environment. Dropped `IMPORTED_TARGET` (transitive
  resolution is not needed — we link against the shared `libhts` directly). Also extended
  the htslib fallback discovery to check `$PREFIX` (set by conda-build) before
  `$CONDA_PREFIX` (set in an activated dev environment), and made the library-file
  probe loop over `.dylib`, `.so`, and `.a` suffixes.

---

## [0.3.0] - 2026-03-24

### Performance

This release focuses on performance improvements across both the C++ BAM scan
stage and the Python quantification stage, delivering a **2.23× throughput
increase** and **25% wall-time reduction** on a production-scale 1.6M-fragment
real BAM (`v0.2.0`: 22.7 s → `v0.3.0`: 17.0 s; throughput 261 K → 581 K
frags/s; peak RSS −3.2 %).

#### Phase 3a — C++ BAM scan stage

- **Zero-copy fragment buffer transfer**: `ResolveContext` adopts `std::move`
  semantics to transfer internal vectors to heap-allocated storage and expose
  them as numpy arrays without any data copies (`finalize_zero_copy()`).

- **Pre-allocated internal vectors**: Added `ResolveContext::reserve()` to
  pre-allocate all accumulator vectors before scanning each chunk, eliminating
  repeated `push_back` reallocation overhead on large chunks.

- **Configurable BGZF decompression threads**: `BamScanConfig.n_decomp_threads`
  (default `4`) controls how many htslib threads decompress the BGZF-compressed
  BAM. Previously hardcoded to 2. Set to `0` to disable multi-threaded
  decompression.

#### Phase 4 — Python quantification stage

- **Vectorized interval merging (Phase 4a)**: Deleted the `_merged_intervals()`
  generator that was called 69 K times (5.3 s cumulative). Replaced with a
  single vectorized numpy merge per locus. Reference chromosome strings are
  factorized to integer codes for fast sort/compare. Merged intervals are now
  computed once in `build_loci()` and cached on `Locus.merged_intervals`.
  `compute_gdna_locus_gammas()` reuses the cached intervals. Combined speedup:
  `build_loci` −80 %, `compute_gdna_locus_gammas` −87 %.

- **Bulk exon CSR construction (Phase 4b)**: Added
  `TranscriptIndex.build_exon_csr()` which converts the per-transcript exon
  dict to four flat numpy CSR arrays in two vectorized passes. Replaces a
  457 K-iteration Python loop that previously fed the C++ fragment scorer.
  The C++ `NativeFragmentScorer` constructor now accepts pre-built numpy arrays
  directly (removing 55 lines of dict-unpacking code). `fragment_scorer`
  speedup: −60 %.

- **Dead code removal and GC hygiene (Phase 4c)**: Deleted the unreachable
  `merge_accumulator_into()` function from `bam_scanner.cpp` (53 lines, never
  called since Phase 3b). Removed a spurious `gc.collect()` call after the
  scoring stage that was consuming 0.69 s; CPython's reference-count GC
  correctly handles the `del` + `release()` pattern without assistance.

### Fixed

- **Annotated BAM `frag_id` desynchronisation with `--no-include-multimap`**:
  When `include_multimap=False`, the write pass was not advancing `frag_id`
  for skipped multimapper groups. This created a deterministic offset between
  the scan-pass `frag_id` assignments and the write-pass lookups, causing
  incorrect per-fragment annotation tags (`ZT`, `ZG`, `ZP`, etc.) for all
  fragments after the first skipped multimapper group. Fixed by incrementing
  `frag_id` in the write pass for every qname group, including skipped ones.

### Changed

- **Lazy-loaded calibration regions**: `TranscriptIndex.region_df` and
  `region_cr` are now `functools.cached_property` attributes instead of
  eagerly loaded in `load()`. The calibration-region index is built on first
  access, reducing startup overhead for tools that access only transcript
  or gene data.

### Platform

- Dropped macOS 13 (x86_64) binary wheels. macOS arm64 wheels now target
  macOS 15.0 (`MACOSX_DEPLOYMENT_TARGET=15.0`).

---

## [0.2.0] - 2026-03-23

### Added

- **gDNA Calibration Framework** (`calibration.py`): New `GDNACalibration`
  class performs aggregate-first region-level EM before per-locus quantification.
  Three convergent signals drive per-region gDNA posterior estimates (γ):
  - *Density*: Gaussian model on log fragment density, fit to expressed vs.
    unexpressed regions.
  - *Strand balance*: Beta-Binomial model with shared overdispersion κ.
    Regions near 0.5 strand balance are consistent with gDNA; those near the
    library strand-specificity are consistent with RNA.  κ is fit from data;
    the LLR vanishes naturally when strand specificity is 0.5.
  - *Fragment length*: Separate RNA and gDNA FL models built iteratively.
  Hard constraint: any region with at least one spliced read has γ forced to 0.
  Per-region posteriors are fragment-count-weighted to the locus level and
  used as gDNA priors in the per-locus EM.

- **Shared-span nRNA Architecture**: nRNA components are defined by unique
  genomic spans `(ref, strand, start, end)` rather than 1:1 transcript
  shadows. Spans are shared across isoforms with the same genomic coordinates,
  reducing EM component count for isoform-rich loci and consolidating intronic
  read evidence.

- **Region-based Evidence Accumulation**: The C++ BAM scanner accumulates
  per-region strand-resolved counts and fragment-length observations in a
  single scan pass. This evidence feeds the calibration EM without a second
  pass over the BAM.

- **Discrete Fragment Assignment Modes**: New `--assignment-mode` parameter
  controls post-EM fragment assignment: `"fractional"` (traditional EM
  weights), `"map"` (assign to highest-posterior component), or `"sample"`
  (stochastic draw). Default: `"sample"`.

- **VBEM Solver Mode**: Variational Bayes EM (VBEM) is now the default solver
  (`--em-mode vbem`). Classic MAP-EM is still available with `--em-mode map`.

- **Posterior Pruning**: `--pruning-min-posterior` removes negligible candidates
  from the CSR data structure before EM, reducing state space and improving
  convergence speed.

- **`nrna_quant.feather`**: New output with nRNA-span-level abundance estimates.
  Complements the per-transcript nRNA values in `quant.feather`.

- **`rigel export` Subcommand**: Converts any `.feather` files in a results
  directory to TSV or Parquet:
  ```bash
  rigel export -o results/ --format tsv
  rigel export -o results/ --format parquet
  ```

- **Category-specific Fragment Length Models**: Separate FL histograms for
  spliced, unspliced, intergenic, and gDNA fragments. Intronic fragment
  lengths use the gDNA FL model from calibration.

- **Exponential Tail Decay**: Fragment lengths beyond `max_frag_length` receive
  exponential log-penalty decay (≈ −0.01/bp) rather than a flat overflow bin.

- **CI/CD Infrastructure**:
  - GitHub Actions workflow for automated testing on Linux and macOS
    (Python 3.12 + 3.13).
  - Automated wheel building and PyPI publishing via OIDC trusted publisher.
  - Binary wheels for Linux (x86_64, aarch64) and macOS (arm64).

- **Bioconda Recipe**: Template recipe at [conda/meta.yaml](conda/meta.yaml).

### Changed

- **Package Naming**: PyPI distribution name is `rigel-rnaseq` (install with
  `pip install rigel-rnaseq`). Import name, CLI command, GitHub repo, and
  Bioconda package remain `rigel`.

- **Default EM Mode**: Changed from `"map"` to `"vbem"`.

- **Default Assignment Mode**: Changed from `"fractional"` to `"sample"`.

- **`priors.py` removed**: Prior logic is now split between `calibration.py`
  (gDNA calibration EM) and `locus.py` (per-locus initialization).

- **Portable Wheel Builds**: `RIGEL_PORTABLE` CMake option omits `-march=native`
  for redistributable wheels. CI builds use this flag automatically.

### Fixed

- **Multimapper gDNA Double-counting**: Fixed a bug where multimapping fragments
  near gDNA-heavy regions were double-counted, inflating gDNA estimates
  particularly with minimap2-aligned BAMs.

- **Splice-junction Gap Fragment Length**: Simplified and corrected the gap
  correction applied to fragment lengths for reads spanning splice junctions.

### Performance

- Competitive with salmon/kallisto in pristine mRNA-only conditions while
  adding joint nRNA and gDNA deconvolution.
- C++ multithreading improvements: reduced synchronisation overhead in the
  BAM scanner and EM solver.
- Link-time optimisation (LTO) enabled across all native extension modules.

---

## [0.1.0] - 2026-03-01

Initial development release.

### Features

- Bayesian transcript quantification with joint mRNA, nascent RNA (nRNA), and genomic DNA (gDNA) deconvolution.
- Linked kinetic model coupling mRNA and nRNA through per-transcript nascent fraction β.
- Hierarchical empirical Bayes priors for nRNA (3-level: global → locus-strand → TSS-group → transcript) and gDNA (2-level: global → chromosome → locus).
- EM algorithm with SQUAREM acceleration and post-EM pruning.
- One Virtual Read (OVR) prior for high-isoform loci.
- C++ native extensions for BAM scanning (htslib), fragment resolution (cgranges), scoring, and EM solver.
- Single-pass pipeline: BAM scan → fragment routing → locus-level EM.
- Outputs: per-transcript, per-gene, and per-locus abundance estimates with mRNA/nRNA/gDNA decomposition.
- Support for stranded and unstranded RNA-seq libraries with auto-detected strand specificity.
- Fragment length models trained separately for RNA and gDNA.
- Coverage weight model for positional bias correction.

[0.3.2]: https://github.com/mkiyer/rigel/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/mkiyer/rigel/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/mkiyer/rigel/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/mkiyer/rigel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mkiyer/rigel/releases/tag/v0.1.0
