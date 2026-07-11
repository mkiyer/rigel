# Changelog

All notable changes to Rigel will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.7.0] - 2026-07-10

### Changed (breaking)

- **Index format version 5 → 7.** Reference indices built with Rigel ≤ 0.6.4 are no longer
  loadable and must be rebuilt with `rigel index`. The calibration solver now reads splice-junction
  geometry directly from a slimmer index schema (vestigial fragment-length model columns were dropped
  at the same time), so `rigel quant` against an old index raises a clear format-version error rather
  than mis-reading it.

### Added

- **Global KDE-mode reference density for capture effective lengths.** The per-region capture
  eff-length contraction now derives its enrichment reference from a single mass-weighted
  kernel-density mode over the library, replacing the old per-transcript `ρ* = G_c / E_c` reference
  that could fire even with no gDNA present (a nascent-RNA over-siphon in capture-enriched loci).
- **Fragment-length reporting.** The pipeline now reports the separate gDNA and RNA empirical
  fragment-length distributions it fits, and drops the vestigial per-category FL models.
- **Validated accumulator-consistent calibration oracle** (`scripts/debug/oracle.py`, CI-checked via
  `tests/calibration/test_oracle.py`) plus a set of diagnostic tools (siphon localization, layer/leak
  tracing, controlled EM experiments) and canonical docs
  (`docs/calibration/oracle_and_benchmarking.md`, `message_precision_roadmap_v2.md`). Developer-facing.

### Changed

- **Unified gDNA / transcript effective-length node model.** gDNA and transcript effective lengths
  now share one per-node contraction path instead of two divergent code paths.
- **Message precision is a single production path** (belief-free Poisson disagreement variance
  `σ²_imp`, fit once before the belief-propagation pass — no per-pass `var~mean` refit and no outer
  fixed-point loop). Exploratory message-precision variants remain available behind the
  default-off `RIGEL_MSG_MODE` env switch for diagnostics.

### Fixed

- **Capture eff-length splice-junction seams** — junction seams are now imputed from their flanking
  exons in the capture eff-length node model, fixing a nascent < mature effective-length inversion in
  capture-enriched loci. Capture-off output is bit-identical to before.

## [0.6.4] - 2026-07-06

### Added

- **Calibration single-strand grid resolution (`--sweep-n-grid-single-strand`, default 256)** — a new
  advanced calibration parameter (`CalibrationConfig.sweep_n_grid_single_strand`, CLI + `--config` YAML).
  After the mixture-bridge fix, the dominant remaining per-node deconvolution error was **grid quantization
  of the single-strand gDNA fraction**: the solve reads `f_g` as the posterior median snapped to the
  log-odds grid, and a deep-count exon's sharp posterior concentrates at the nearest of the (coarse, shared)
  `sweep_n_grid=60` points, quantizing `f_g` to Δ≈0.085 steps — the true off-grid value differs by up to
  ±0.04, and × high exon mass this was ~61% of the residual. Single-strand nodes solve a cheap 1-D grid, so
  a fine grid (256) de-quantizes the readout there; it is **decoupled** from the AMBIG 2-D `(λ,τ)` cube
  (kept at `sweep_n_grid=60`, since a fine grid there is a genome-scale memory risk), with the coarse global
  prior regridded onto the fine grid. Keeps the **median** estimator (transform-invariant, vertex-safe — a
  parabolic sub-grid *mode* was rejected because it under-calls confident pure-gDNA nodes at the f_g→1
  vertex). Target calibration Σ|err| **−45.6%**; full-suite mature-transcript accuracy flat-to-better with no
  regression (abundant-tx Spearman 0.959→0.962, n_FP 413→403). Design:
  `docs/calibration/fix3_single_strand_grid_quantization.md`.

## [0.6.3] - 2026-07-06

### Added

- **Calibration gDNA-prior mixture bridge (`--gdna-prior-mixture-bridge`, default 0.01)** — a new advanced
  calibration parameter (`CalibrationConfig.gdna_prior_mixture_bridge`, exposed on the CLI and in the
  `--config` YAML). The Phase-2 gDNA-density KDE is estimated from clean (unimodal) region nodes, so it
  develops a deep **valley** between the depleted and enriched modes. A node whose current-belief gDNA
  density is a spatial **mixture** of enriched (in-probe) and depleted (off-target) positions — a capture
  boundary, or a sparsely-probed region — lands in that valley by construction and **collapses to `f_g≈0`**,
  emitting a pathologic RNA message that then crushes the gDNA fraction of its AMBIG-exon neighbours. The
  bridge mixes the KDE with a uniform "any-mixture" prior over the observed density support at weight ε,
  which **floors the valley** (no collapse) while leaving the KDE's real Gaussian tails outside the support
  intact (high-density false-positive suppression unchanged). On the capture-on failure conditions this fixes
  the worst per-node deconvolution errors (the #1 error node's gDNA fraction 0.62→0.93 vs a true 0.92) and
  **reduces the false-nascent hallucination ~22%** at the pool level (gDNA recovered from RNA); it is inert
  on capture-off and gDNA-free libraries. `0` disables it (bit-exact legacy KDE). Design + validation:
  `docs/calibration/boundary_kde_valley_collapse_and_simplex_precision.md`.

## [0.6.2] - 2026-07-05

### Fixed

- **Genome-scale calibration out-of-memory (production blocker)**: on a whole-genome (human) index the
  gDNA-density KDE prior evaluated `P(log ρ_g)` at every node × solve-grid point (~90M points for a human
  index) and built a single `(n_eval, n_train)` matrix — a **2.69 TiB** allocation that crashed
  `rigel quant` the first time it ran on real data. The KDE evaluation now (a) **tiles the query axis** so
  the pairwise matrix is bounded at any scale, and (b) for genome-scale queries **tabulates the exact
  kernel on a bandwidth-scaled lattice** spanning the query range and interpolates — the real quadratic
  tails are preserved (the lattice covers the full range, so nothing is clamped). Peak RSS on a
  971k-fragment human library went from an OOM crash to ~10 GB; small/golden inputs keep the exact
  per-point path (bit-identical). Regression tests in `tests/calibration/test_gdna_density_prior.py`.

### Changed

- **Faster, leaner calibration (real-data production prep)** — three optimizations to the per-node
  belief-propagation solve, all verified within the golden tolerance (rtol=1e-6): the calibration stage
  drops from ~110 s to **~66 s (−40%)** and peak RSS from ~10 GB to **~8.6 GB** on a 971k-fragment human
  library, with the gDNA/RNA deconvolution unchanged (16-condition benchmark net flow within run-to-run
  noise):
  - lean numpy log-sum-exp in `simplex_logodds.py` (replacing `scipy.special.logsumexp`; matches it to
    ~1e-15 without the scipy wrapper overhead) — byte-identical on the single-strand float64 path;
  - in-place `psi` accumulation + common-subexpression reuse in the AMBIG (λ,τ) cube — byte-identical;
  - **float32 AMBIG cube with float64 reductions** — the (m,K,Kt) log-posterior is stored/evaluated in
    float32 (≈½ the cube memory, ~1.7× the elementwise exp/log) while every reduction accumulates in
    float64 and the (m,K) marginals stay float64, so medians/means/variances keep precision. Deviation
    vs float64: `f_g` identical, RNA fractions ≤7.7e-7.
- **Deterministic calibration σ²_g (M2)**: replaced the bistable spline gDNA-density estimator with a
  deterministic one, eliminating cross-process nondeterminism in the prior. *(committed since v0.6.1)*

## [0.6.1] - 2026-07-04

### Added

- **`rigel index --collapse-duplicate-transcripts`**: collapse transcripts that
  share identical exon coordinates (same reference, strand, and every exon
  start/end) instead of aborting the index build. Keeps the lexicographically-
  smallest transcript ID in each duplicate group and drops the rest, logging a
  summary of what was collapsed. Such transcripts are mathematically
  unidentifiable in quantification, so collapsing to one representative is
  loss-free. The default is unchanged — a hard `ValueError` reporting the
  offending groups — but the error message now points at this flag. GENCODE
  ships a handful of byte-identical transcript annotations that previously
  required manual GTF pre-processing.

## [0.6.0] - 2026-07-04

### Highlights

Rigel 0.6.0 is a **ground-up rewrite of the calibration stage** — the step that
deconvolves each library into genomic-DNA contamination (gDNA) versus RNA before
the per-locus EM runs. The earlier "Simple Regional Deconvolution" calibrators
(the never-released SRD v1 / v2 iterations) are gone, replaced by a **bipartite
belief-propagation sweep** over a region↔boundary node chain:

- Each node's unspliced fragment mass is deconvolved into a three-way pie —
  **sense-RNA, antisense-RNA, gDNA** — from exactly three sources: a per-strand
  Beta-Binomial strand likelihood, cross-node density imputation between
  neighbouring nodes, and a population gDNA prior. (A fragment *count* carries no
  intrinsic gDNA/RNA information; it enters only as strand-likelihood precision.)
- **Boundaries are first-class nodes**: they own the one-sided, motif-stranded
  spliced RNA and feed the per-locus prior their per-side gDNA/RNA flux.
- gDNA flows genomically; per-strand RNA flows only where a transcript strand is
  structurally continuous. The calibration result enters the per-locus EM as two
  Dirichlet scalars that set the gDNA-vs-RNA split; the EM then distributes RNA
  mass across the compatible transcripts.

The full pipeline (scan → calibrate → quant) runs end-to-end. Theory and design:
`docs/calibration/CALIBRATION_ARCHITECTURE.md` (the count-zero-information
principle) and the docs index in `docs/README.md`.

> **Pre-release.** 0.6.0 is published to a personal channel for validation on
> real datasets. Expect further calibration tuning before a public release.

### Added

- **Joint mRNA / nRNA / gDNA quantification** with an in-place per-region /
  per-boundary fractional accumulator (four channels: unspliced ±, spliced
  sense / antisense) that runs inside the single-pass C++ BAM scan.
- **Fine-grained calibration region partition** built at `rigel index` time — a
  4-bit exon/intron × strand signature per region, persisted in the index and
  shared with the C++ scanner.
- **`--splicing-anchor-tolerance K`** (default 3 bp): minimum bp clearance on
  each side of an exon–intron boundary before a fragment counts as
  boundary-crossing — removes single-bp soft-clip / indel-slippage artefacts
  near GT-AG splice motifs.
- **gDNA EM controls**: `--gdna-em-llr-bias`, `--gdna-splice-penalty-unannot`.

### Changed (Breaking)

- **Index format bumped to `INDEX_FORMAT_VERSION = 5`.** Indexes built by older
  Rigel fail to load with an explicit rebuild message — rebuild with
  `rigel index --fasta … --gtf …`.
- **Calibration API rewritten.** `rigel.calibration.calibrate()` now takes the
  accumulator payload, region arrays, strand model, gDNA/RNA fragment-length
  pmfs, and config, and returns a `CalibrationResult` carrying per-region and
  per-boundary-side deconvolved gDNA/RNA mass. `quant_from_buffer()` and
  `FragmentScorer.from_models()` signatures changed to match.
- **`summary.json` calibration schema replaced.** Any consumer reading the old
  SRD keys (`srd.*`, `pi_pool`, `gdna_fl.*`) must update.
- **Native scored-fragment payloads are float32** (`log_liks`,
  `coverage_weights`, `gdna_log_liks`), promoted back to double inside the locus
  EM — lower peak RSS on large runs.

### Removed (Breaking)

- The entire SRD calibration surface and its `--cal-*` CLI flags, plus the legacy
  per-transcript RNA-prior machinery (`prior_weight_rna`, `build_prior_weight_rna`,
  `nrna_weight`, `rna_prior_count`). The production EM prior is now exactly two
  per-locus Dirichlet scalars (gDNA / RNA).

### Fixed

- **`SPLICED_IMPLICIT` misclassification of unspliced fragments**: true-gDNA
  fragments whose contiguous paired-end alignment merely overlapped an annotated
  intron were promoted to `SPLICED_IMPLICIT`, hard-gated out of the gDNA
  likelihood, and mis-attributed to nRNA. Replaced with a per-intron
  whole-containment test (one-sided slack `K = splicing_anchor_tolerance`).
- **Fast-math-safe multimapper gDNA accumulation**: removed an infinity sentinel
  in the native scorer's gDNA log-sum-exp that could emit `-inf` gDNA likelihoods
  for otherwise-valid unspliced multimappers under `-ffast-math`.
- **Accumulator deposit mis-routing on multi-reference genomes**: per-region /
  per-boundary mass was routed by first-seen reference order rather than the
  index ref-id, silently dropping deposits on multi-contig genomes.

### Internals

- Calibration rebuilt across the region↔boundary node chain (`node_chain`), the
  belief-propagation solver (`bp_solver`), and the per-node 2-simplex solve
  (`simplex_sweep` / `simplex`); see CLAUDE.md for the full module map.
- Full test suite green (1084 tests). `tests/golden/*` regenerated against the
  rewritten calibration output.

## [0.4.1] - 2026-04-27

### Changed (Breaking)

- **Calibration upgraded to Simple Regional Deconvolution (SRD v2)**
  (breaking): the calibration scanner now emits four strand-aware
  collapsed overlap counts per fragment (`exon_bp_pos/neg`,
  `tx_bp_pos/neg`) instead of per-candidate `exon_bp` arrays.
  Categorization is rewritten as pure column arithmetic over a
  `(N_CATEGORIES, 4_strand_sublabels)` matrix. Two new splice types
  are recognized: `SPLICED_IMPLICIT` (paired-end gap collapses an
  annotated intron) and `SPLICE_ARTIFACT` (alignment crosses a
  blacklisted SJ). Zero-candidate (intergenic) fragments now flow
  through the buffer and are categorized as INTERGENIC by the same
  geometric rule rather than via a side-channel accumulator.
  (The SRD v2 design/results notes referenced here have since been
  superseded; see the calibration history under `docs/calibration/archive/`.)

- **Headline gDNA-FL recovery**: the v1 saturation failure mode
  (gDNA-FL mode pinned at bin-0, bin-0 fraction up to 55%) is
  eliminated. v2 recovers gDNA-FL mode 256–336 bp on the VCaP spike
  series with bin-0 / bin-max edge fractions ≈ 0%.

- **Headline gDNA fraction shifts ~20% lower across the board**
  (correction, not regression): v1 was attributing edge-saturated
  RNA mis-routed mass to gDNA. The relative monotonicity across the
  spike series is preserved and tightened.

### Removed (Breaking)

- `FragmentLengthModels.intergenic` member and
  `observe_intergenic_batch()` method (replaced by buffer-routed
  INTERGENIC categorization).
- `StrandModels.intergenic` member (legacy diagnostic; no longer trained).
- `n_frag_length_intergenic` counter on `PipelineStats` and the
  C++ `intergenic_obs/truth/lengths` accumulators on `BamScanStats`,
  `StrandObservations`, and `FragLenObservations`.
- The `frag_length_models.intergenic` and
  `strand_models.diagnostics.intergenic` keys are no longer present
  in JSON summary output.

### Added

- `category_counts` shape widened from `(N_CATEGORIES,)` to
  `(N_CATEGORIES, 4)` with strand sub-labels `(none, pos, neg, ambig)`.
- New diagnostic fields on `CalibrationResult`: `n_spliced`,
  `n_pool_intronic_strand_pos`, `n_pool_intronic_strand_neg`.
- `splice_type` enum gains `SPLICED_IMPLICIT = 3` and
  `SPLICE_ARTIFACT = 4`.
- Buffer schema gains four int32 columns: `exon_bp_pos`,
  `exon_bp_neg`, `tx_bp_pos`, `tx_bp_neg` (16 bytes/fragment).

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

- **Configurable BGZF decompression threads**: `BamScanConfig.bgzf_threads`
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
[0.6.0]: https://github.com/mkiyer/rigel/compare/v0.4.0...v0.6.0
[0.6.1]: https://github.com/mkiyer/rigel/compare/v0.6.0...v0.6.1
[0.6.2]: https://github.com/mkiyer/rigel/compare/v0.6.1...v0.6.2
[0.6.3]: https://github.com/mkiyer/rigel/compare/v0.6.2...v0.6.3
[0.6.4]: https://github.com/mkiyer/rigel/compare/v0.6.3...v0.6.4
[0.7.0]: https://github.com/mkiyer/rigel/compare/v0.6.4...v0.7.0
