# Changelog

All notable changes to Rigel will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.5.0] - 2026-04-24

### Changed (Breaking)

- **Calibration replaced with Simple Regional Deconvolution (SRD v1)**
  (breaking): the v5 regional-density calibration is gone. Calibration
  now classifies every uniquely-aligned fragment into one of seven
  geometric categories (SPLICED, UNSPLICED_SENSE_EXONIC,
  UNSPLICED_ANTISENSE_EXONIC, UNSPLICED_EXONIC_AMBIG, EXON_INCOMPATIBLE,
  INTRONIC, INTERGENIC) using only the per-candidate `exon_bp` overlap
  computed by the C++ scanner, then fits a 1-D fragment-length mixture
  `pool ∝ π·gDNA(L) + (1−π)·RNA(L)` over the geometric pool
  `EXON_INCOMPATIBLE ∪ INTRONIC ∪ INTERGENIC` plus the scanner's
  intergenic-FL accumulator. Per-locus Dirichlet priors are derived from
  per-fragment posteriors via `compute_locus_priors_from_partitions`.
  No regional density, no per-region exposure, no SS-threshold magic
  numbers. See `docs/calibration/srd_v1_implementation.md` and
  `docs/calibration/srd_v1_results.md`.

- **`CalibrationResult` schema rewritten** (breaking): v5 fields
  (`lambda_gdna`, `density_*`, region-evidence-derived statistics) are
  removed. New schema includes `gdna_fl_quality`, `pi_pool`,
  `n_pool`, `n_intergenic_unique`, `mixture_converged`,
  `mixture_iterations`, `category_counts`, `exon_fit_tolerance_bp`,
  `fl_prior_ess`. Consumers reading `summary.json` must update.

- **Index format: `regions.feather` removed** (breaking): the offline
  region partition is no longer built or stored. The `--no-mappability`
  flag and `mappable_effective_length` per-region column are removed.
  Existing indexes built with v0.4.x can still be loaded; the missing
  table is tolerated. Re-indexing is **not** required.

- **`TranscriptIndex.region_df` and `region_cr` removed** (breaking):
  any external code that touched `index.region_df`, `index.region_cr`,
  or `rigel.index.build_region_table` must be updated. SRD calibration
  has no dependency on these tables.

- **C++ scanner: `RegionAccumulator` and `region_evidence` output
  removed** (breaking): the `bam_scanner.cpp` no longer emits
  `region_evidence` from `scanner.scan(...)`. Downstream consumers of
  the result dict should not reference this key.

- **`rigel.mappability.uniform_region_exposures` and
  `compute_mappable_effective_length` removed** (breaking): per-region
  mappability machinery is gone. Per-transcript effective length is
  unaffected — it is computed by
  `frag_length_models.rna_model.compute_all_transcript_eff_lens`
  during scoring.

### Added

- **SRD calibration modules**: `src/rigel/calibration/_simple.py`,
  `_categorize.py`, `_fl_mixture.py`, `_fl_empirical_bayes.py`,
  `_result.py`. Total ~400 LOC replacing ~1000 LOC of v5 calibration.

- **Per-locus γ → Dirichlet prior**: `compute_locus_priors_from_partitions`
  in `src/rigel/locus.py` derives per-locus α from per-fragment SRD
  posteriors. Replaces the v5 density-based per-locus prior; restores
  scenario-test pass rates.

- **C++ intergenic fragment-length accumulator**: the BAM scanner
  populates `frag_length_models.intergenic` for unique-mapper unspliced
  fragments that resolve to zero transcript candidates. SRD folds these
  into Pool B for a higher-statistics gDNA-FL signal.

- **Tests**: `tests/test_categorize.py` (15), `tests/test_fl_mixture.py`,
  `tests/test_fl_empirical_bayes.py`, `tests/test_calibration_simple.py`
  (10), `tests/test_locus_priors.py`. All previously-passing tests
  remain green; goldens regenerated.

### Fixed

- **`exon_fit` rule overhang bug**: the original SRD Pass 0 rule used
  `min(intron_bp_per_candidate) ≤ tol`. `intron_bp` excludes read bases
  that overhang transcript edges into intergenic territory, so 100+ bp
  gDNA-splash reads were mis-classified `UNSPLICED_SENSE_EXONIC` and
  Pool A was empty for every gDNA-rich library. New rule:
  `(read_length − max(exon_bp_per_candidate)) ≤ tol`. Regression test
  reproduces the production scenario at chr11:48244998 (192 bp read,
  87 bp exonic, 105 bp overhang).

- **Intergenic FL never folded into Pool B**: the C++ scanner drops
  zero-candidate unspliced fragments before they reach the
  FragmentBuffer (resolve_context.h:1043), so they never reached
  Pass 0. Their FL was being captured in `frag_length_models.intergenic`
  but not consumed. SRD now reads from this accumulator; dna80m gained
  976 k unique gDNA fragments.

- **Mixture EM iteration cap**: mid-π libraries (`gdna_fraction ∈
  [0.2, 0.5]`) genuinely converge slowly (per-step `|Δπ| ≈ 1e-4` near
  iter 200–300). Default `max_iter` raised 200 → 1000; cost <1 s/library.

### Removed

- v5 calibration modules (`_em.py`, `_stats.py`, `_fl_model.py` (old),
  `_calibrate.py` (old), `_annotate.py`).
- Region partition machinery: `build_region_table`,
  `TranscriptIndex.region_df`, `region_cr`, `regions.feather`,
  `mappable_effective_length` per-region column,
  `rigel.mappability.uniform_region_exposures`,
  `compute_mappable_effective_length`.
- C++ `RegionAccumulator`, `build_region_index`, `has_region_index`,
  `region_cr_` field, and the `region_evidence` output channel from
  `bam_scanner.cpp`.
- Tests: `tests/test_regions.py` (entire file), region portions of
  `tests/test_mappability.py` and `tests/test_calibration_integration.py`.

### Validation

VCaP mixture benchmarks (8 libraries, 20 M cell-line RNA-Seq combined
with 0–80 M reads of pure exome-capture DNA) all converge with
`gdna_fl_quality=good`. Recovered gDNA fraction tracks nominal within
±2 pp at intermediate contamination (5–50%). mRNA linear Pearson r ≥
0.97 across the entire ladder. See `docs/calibration/srd_v1_results.md`
for the full report.

---

## [0.4.0] - 2026-04-21

### Changed (Breaking)

- **Annotated BAM `ZF` tag redesigned** (breaking): `ZF` is now an 8-bit
  bitfield with explicit flags for every fragment outcome. New bits:
  `is_resolved` (0x01), `is_mrna` (0x02), `is_gdna` (0x04), `is_nrna` (0x08),
  `is_synthetic` (0x10), `is_intergenic` (0x20), `is_chimeric` (0x40),
  `is_multimapper_dropped` (0x80). Canonical values: `0x03` mRNA, `0x09` nRNA,
  `0x19` synthetic nRNA, `0x05` EM-gDNA, `0x25` intergenic gDNA, `0x40`
  chimeric, `0x80` dropped multimapper. Consumers that decoded the legacy
  bits (`is_resolved=0x1, is_gdna=0x2, is_nrna=0x4, is_synthetic=0x8`) or
  the legacy numeric values `{0, 1, 3, 5, 13}` must update. The old
  `is_resolved` bit did not distinguish mRNA from nRNA; readers should now
  branch on the explicit `is_mrna` / `is_nrna` / `is_gdna` bits.

- **Annotated BAM `ZC` tag narrowed** (breaking): `ZC` now reflects input
  ambiguity only and takes one of `{unambig, ambig_same_strand,
  ambig_opp_strand, multimapper}` on scored records, or `"."` on records
  that were not scored (chimeric, intergenic, dropped multimappers). The
  former `chimeric` and `intergenic` string values have been removed — that
  information is now exclusively carried in `ZF`.



### Added

- **Annotated BAM: gene name tag (`ZR`)**: the annotated BAM now writes a `ZR:Z`
  tag containing the gene name/symbol for the assigned fragment (`"."` for
  intergenic), enabling fast human-readable lookups without a separate index
  query.

- **Annotated BAM: locus ID tag (`ZL`)**: new `ZL:i` integer tag records the
  locus ID (zero-based) for each resolved fragment (`-1` if no locus was
  assigned). Enables grouping and filtering by locus in downstream tools.

- **Annotated BAM: filtered-read passthrough**: reads filtered during BAM
  scanning (QC-fail, unmapped, duplicates, unpaired) are now written through
  to the annotated BAM without annotation tags. Previously these records were
  silently dropped, making the output read count differ from the input. The
  summary JSON gains a corresponding `n_filtered_passthrough` counter.

- **`RigelIndex.nrna_to_transcripts()` / `nrna_to_genes()`**: new lookup
  methods on the index object that, given an nRNA entity ID, return the
  contributing transcript IDs or `(gene_id, gene_name)` pairs, respectively.
  Useful for interpreting nRNA quantification results.

- **VBEM clamp floor (`VBEM_CLAMP_FLOOR = 0.1`)**: VBEM SQUAREM now clamps
  all alpha values to a minimum of 0.1 after each extrapolation step.
  Prevents components from entering the digamma absorbing regime
  (`ψ(α) ≈ −1/α` for small `α`) from which recovery is impossible in
  double precision, eliminating the catastrophic zeroing failure mode that
  caused spurious zero-abundance estimates for lowly expressed transcripts.

- **Per-locus profiling statistics**: `AbundanceEstimator` gains an
  `emit_locus_stats` option that populates `estimator.locus_stats` with
  per-locus timing breakdowns (extract, bias, build_ec, warm_start,
  SQUAREM, assign phases in microseconds) and iteration counts. Useful for
  diagnosing bottlenecks on mega-loci.

- **AVX2 4-wide `fast_exp`**: `fast_exp.h` now provides a 4-wide double-
  precision `fast_exp_avx2()` using Cody-Waite range reduction and a
  degree-11 Horner polynomial with FMA. Paired with the existing AVX-512
  8-wide path, the E-step automatically selects the widest available SIMD
  lane width at compile time.

### Changed

- **Annotated BAM: `ZP` pool tag replaced by `ZF` assignment flags bitfield**
  *(breaking)*: the string pool tag `ZP` (`mRNA`/`nRNA`/`gDNA`/`intergenic`/
  `chimeric`) is replaced by `ZF:i`, a compact integer bitfield:
  - bit 0 (`0x1`) — fragment was scored and assigned by EM (`is_resolved`)
  - bit 1 (`0x2`) — EM gDNA component won (`is_gdna`)
  - bit 2 (`0x4`) — assigned transcript is single-exon (`is_nrna`)
  - bit 3 (`0x8`) — assigned transcript is a rigel-generated nRNA span (`is_synthetic`)

  Common values: `0` = unresolved/intergenic, `1` = multi-exon mRNA, `3` =
  gDNA, `5` = annotated single-exon nRNA, `13` = synthetic nRNA span. The old
  `ZP` string decoding is no longer supported.

- **nRNA transcript classification redesigned** *(breaking)*: the index
  transcript table columns `is_synthetic_nrna` and `is_nascent_equiv` are
  replaced by a cleaner two-column scheme:
  - `is_nrna` — `True` for any single-exon transcript used as an nRNA
    component (both annotated transcripts and generated spans)
  - `is_synthetic` — `True` only for rigel-generated nRNA spans
    (`RIGEL_NRNA_*` IDs); `False` for annotated single-exon transcripts

  Annotated single-exon transcripts now correctly contribute to fragment-
  length model training (only synthetic spans are excluded, as their
  genomic spans do not represent real insert sizes).

- **Partition-native EM dispatch**: the internal `batch_locus_em` C++ entry
  point is replaced by `batch_locus_em_partitioned`, which accepts
  pre-partitioned per-locus CSR arrays (`LocusPartition` objects) produced
  by the new `partition_and_free()` helper. Global CSR arrays are freed
  array-by-array immediately after scatter, substantially reducing peak
  memory usage on large samples without changing quantification results.

- **E-step memory layout**: the heap-allocated `ec.scratch` matrix (one row
  per fragment per candidate) is eliminated. The E-step now uses a
  stack-local row buffer (≤512 candidates) or a small heap vector for wider
  equivalence classes, combined with fused per-column Kahan accumulators that
  accumulate column sums in a single pass. Reduces memory traffic by
  ~50% for typical loci and improves cache utilization.

- **`splice_type` value rename**: the splice-type label `"ambiguous"` is
  renamed to `"unknown"` throughout the annotated BAM and output tables for
  clarity.

### Fixed

- **nRNA summary bug**: fixed an accounting error in the nRNA quantification
  summary that caused synthetic nRNA span counts to be misreported.

## [0.3.3] - 2026-03-25

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
[0.4.0]: https://github.com/mkiyer/rigel/compare/v0.3.3...v0.4.0
