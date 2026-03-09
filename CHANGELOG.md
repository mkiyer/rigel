# Changelog

All notable changes to Rigel will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-03-09

### Added

- **Nascent RNA Decoupling Architecture**: nRNA components are now defined by unique genomic spans `(ref, strand, start, end)` rather than 1:1 transcript shadows. This dramatically reduces EM component count for genes with many isoforms sharing the same genomic coordinates (N ≪ T), improves convergence speed, and consolidates intronic read evidence.

- **Comprehensive Documentation**:
  - [docs/MANUAL.md](docs/MANUAL.md): Complete user manual with all commands, parameters, output formats, examples, and FAQ.
  - [docs/METHODS.md](docs/METHODS.md): Publication-quality mathematical methods document (14 sections) describing the complete Bayesian framework, generative model, EM algorithm, hierarchical priors, and implementation details.
  - [docs/PUBLISHING.md](docs/PUBLISHING.md): Step-by-step guide for publishing releases to PyPI and Bioconda.

- **CI/CD Infrastructure**:
  - GitHub Actions workflow for automated testing on Linux and macOS (Python 3.12 + 3.13).
  - Automated wheel building and PyPI publishing via trusted publisher (OIDC).
  - Binary wheels for Linux (x86_64, aarch64) and macOS (x86_64, arm64).

- **Bioconda Recipe**: Template recipe at [conda/meta.yaml](conda/meta.yaml) for publishing to the Bioconda channel.

### Changed

- **Package Naming**: PyPI distribution name is now `rigel-rnaseq` (install with `pip install rigel-rnaseq`) to avoid conflict with existing PyPI package. Import name, CLI command, GitHub repo, and Bioconda package remain `rigel`.

- **Portable Wheel Builds**: Added `RIGEL_PORTABLE` CMake option to build wheels without `-march=native`, ensuring binaries run on any machine of the same architecture. CI builds use this flag automatically.

### Performance

- **C++ Multithreading Optimizations**: Enhanced thread work-stealing and reduced synchronization overhead in EM solver and BAM scanner.

- **BAM Reading Parallelism**: Improved htslib thread pool utilization for BGZF decompression.

- **Double Precision EM**: Switched to `double` precision (64-bit float) throughout the EM kernel for improved numerical stability, particularly for loci with many transcripts.

- **Memory Efficiency**: Reduced peak memory usage through optimized buffer management and columnar fragment storage.

### Fixed

- **Documentation Accuracy**: Corrected gDNA model description in all documentation from "per-gene" to "per-locus" to match actual implementation. gDNA is modeled as a single shadow per locus at component index `[2T]`, not per-gene.

- **CMake Compiler Flags**: Reorganized optimization flags with separate handling for `-ffast-math` (applied to scoring but not EM to preserve Kahan summation accuracy).

### Technical

- **Build System**: Verified compatibility with scikit-build-core ≥0.4.3 and nanobind ≥2.0 for stable ABI wheel builds (CPython 3.12+).

- **Link-Time Optimization (LTO)**: Enabled interprocedural optimization for cross-translation-unit inlining and dead code elimination.

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

[0.2.0]: https://github.com/mkiyer/rigel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mkiyer/rigel/releases/tag/v0.1.0
