# Changelog

All notable changes to Rigel will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
  - Binary wheels for Linux (x86_64, aarch64) and macOS (x86_64, arm64).

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

[0.2.0]: https://github.com/mkiyer/rigel/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mkiyer/rigel/releases/tag/v0.1.0
