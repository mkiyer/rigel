# Rigel Parameters Reference

All parameters below apply to `rigel quant` unless noted.
Resolution order: **explicit CLI flag → YAML config file → built-in default**.

See also: [MANUAL.md](MANUAL.md) (usage walkthrough), [METHODS.md](METHODS.md) and
[ARCHITECTURE.md](ARCHITECTURE.md) (pipeline internals), [BENCHMARKING.md](BENCHMARKING.md),
[SIMULATOR.md](SIMULATOR.md), and the authoritative calibration theory in
[calibration/CALIBRATION_ARCHITECTURE.md](calibration/CALIBRATION_ARCHITECTURE.md).

---

## rigel index

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Genome FASTA (must have `.fai` index) |
| `--gtf` | required | Annotation GTF (GENCODE recommended) |
| `-o`, `--output-dir` | required | Output directory for index files |
| `--alignable-zarr PATH` | — | Alignable Zarr store built for the same genome+aligner. Provides per-base fractional mappability (for gDNA-aware effective length) and the splice-junction artifact blacklist (applied at BAM-scan time). Required unless `--no-mappability` is set. |
| `--no-mappability` | off | Explicitly opt out of mappability and the splice-artifact blacklist (synthetic genomes, stranded-only benchmarks). Mutually exclusive with `--alignable-zarr`. |
| `--splice-blacklist-min-count` | `2` | (Advanced) Minimum unique-fragment support per `(chrom, intron, read_length)` for a junction to enter the blacklist. Ignored under `--no-mappability`. |
| `--mappability-read-length` | `100` | (Advanced) Read-length bin queried in the alignable store (must match a built bin, commonly 50/75/100/125). Ignored under `--no-mappability`. |
| `--nrna-tolerance` | `20` | Max distance (bp) for clustering transcript start/end sites into shared synthetic nRNA spans |
| `--collapse-duplicate-transcripts` | off | Collapse transcripts sharing identical exon coordinates (keeps the lexicographically-smallest ID per group) instead of failing. Loss-free — such transcripts are unidentifiable in quantification. |
| `--gtf-parse-mode` | `strict` | `strict` fails on malformed GTF records; `warn-skip` logs and skips them |
| `--feather-compression` | `lz4` | Feather file compression: `lz4`, `zstd`, or `uncompressed` |
| `--no-tsv` | off | Skip writing TSV mirrors of index files |

Index artifacts (feather, plus `.tsv` mirrors unless `--no-tsv`): `transcripts`, `intervals`,
`regions` (the calibration region partition), `boundaries`, `ref_lengths`, `sj`,
`splice_blacklist`, and `manifest.json`. `regions.feather` and `boundaries.feather` are **core**
calibration artifacts.

---

## rigel quant

### Required

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted BAM file. Minimap2 and STAR are supported. `NH` tag used when present; otherwise multimappers are detected from secondary BAM flags. |
| `--index` | Rigel index directory (from `rigel index`) |
| `-o`, `--output-dir` | Output directory |

### Library and input

| Flag | Default | Description |
|------|---------|-------------|
| `--config FILE` | — | YAML file with parameter values (same keys as CLI options, with underscores); CLI flags override YAML. A `config.yaml` is written to the output directory after each run for reproducibility. |
| `--include-multimap` / `--no-include-multimap` | yes | Include multimapping reads (detected via `NH` tag or secondary flag) |
| `--keep-duplicates` / `--no-keep-duplicates` | no | Keep PCR/optical duplicate-marked reads |
| `--sj-strand-tag TAG [TAG …]` | `auto` | Splice-junction strand BAM tag(s). `auto` detects from the first 10,000 reads. Use `XS` (STAR), `ts` (minimap2), or list multiple tags to try in order (e.g. `XS ts`). |

### EM algorithm

| Flag | Default | Description |
|------|---------|-------------|
| `--em-iterations` | `1000` | Maximum EM iterations. `0` = skip EM entirely (unambiguous fragments only). |
| `--em-mode` | `vbem` | EM algorithm variant: `vbem` (Variational Bayes EM with digamma-based soft updates, default) or `map` (MAP-EM with hard `max(0, n+a−1)` updates). |
| `--assignment-mode` | `sample` | Post-EM fragment assignment mode: `sample` (draw from the posterior distribution, default), `fractional` (preserve EM posterior weights), `map` (assign each fragment to its argmax component). |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--tsv` | off | Also write `.tsv` mirrors of the quant tables |
| `--annotated-bam PATH` | — | Write an annotated BAM with per-fragment assignment tags (`ZT`, `ZG`, `ZR`, `ZI`, `ZJ`, `ZF`, `ZW`, `ZC`, `ZH`, `ZN`, `ZS`, `ZL`, `ZB`). Requires a second BAM pass. |
| `--emit-locus-stats` | off | Write per-locus EM convergence profiling (iteration counts, timing, equivalence-class stats) to `locus_stats.feather`. |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | `0` | Total thread budget for Rigel. During BAM scan the budget is split between scan workers and `--scan-bgzf-threads`; locus EM reuses the same budget (stages run serially). `0` = all cores, `1` = sequential. |
| `--scan-bgzf-threads N` | `4` | BGZF/BAM decompression threads reserved from `--threads` during the scan phase. `0` disables htslib threaded decompression. Capped so at least one scan worker remains. |
| `--scan-buffer-size GB` | `2` | Maximum scan buffer size in GiB before finalized chunks are spilled to disk. Increase with ample RAM to reduce spill I/O. |
| `--scan-fragments-per-chunk N` | `1000000` | Buffered fragments per native scan chunk before handoff to the Python buffer. Larger chunks reduce overhead but raise transient memory. |
| `--scan-read-name-batch-size N` | `512` | (Advanced) Read-name groups per native scanner input-queue item. Lower helps queue-boundary debugging; larger may reduce queue overhead at slightly more in-flight memory. |
| `--tmpdir DIR` | system temp | Directory for fragment-buffer spill files. Use a fast local SSD. |
| `--seed N` | timestamp | Random seed for reproducibility. Affects post-EM `sample` assignment. |

---

## Advanced options

These control scoring sensitivity, EM numerical behavior, and the calibration
deconvolution. The defaults suit standard total RNA-seq libraries.

### Scoring and EM

| Flag | Default | Description |
|------|---------|-------------|
| `--assignment-min-posterior` | `0.01` | Minimum component posterior for discrete assignment (`map`/`sample` modes only). Components below this are zeroed before assignment. |
| `--em-convergence-delta` | `1e-6` | EM convergence threshold for parameter updates (‖Δθ‖). |
| `--pruning-min-posterior` | `1e-4` | Remove CSR candidates with posterior below this threshold before EM. Lower keeps more candidates (conservative); `0` disables pruning. |
| `--overhang-alpha` | `0.1` | Per-base overhang penalty α ∈ [0, 1]. `0` = hard gate, `1` = no penalty. |
| `--mismatch-alpha` | `0.1` | Per-mismatch (NM tag) penalty α ∈ [0, 1]. `0` = hard gate, `1` = no penalty. |
| `--splicing-anchor-tolerance K` | `3` | Resolver-side splicing-anchor tolerance in bp. Used only for implicit-splice resolution slack around annotated introns; the fractional calibration accumulator does not interpret this value. |

### Calibration deconvolution

The calibration stage deconvolves each genomic node's unspliced fragment mass into the simplex
`(f_rna₊, f_rna₋, f_g)` via a single forward-backward belief-propagation sweep over the
region↔boundary chain (see [calibration/CALIBRATION_ARCHITECTURE.md](calibration/CALIBRATION_ARCHITECTURE.md)).
These knobs tune that sweep and the downstream gDNA handling in the EM.

| Flag | Default | Description |
|------|---------|-------------|
| `--gdna-prior-mixture-bridge ε` | `0.01` | Phase-2 gDNA-density KDE prior mixture-bridge weight ε ∈ [0, 1). Floors the KDE valley between the depleted and enriched modes so a mixture-density node (a capture boundary, a sparse-probe region) is read as an enriched/depleted mixture instead of collapsing to ~0 gDNA. `0` disables it (legacy bit-exact KDE). |
| `--sweep-n-grid-single-strand K` | `256` | Single-strand log-odds grid resolution. Single-strand nodes solve a cheap 1-D grid, so a fine grid de-quantizes the gDNA-fraction readout (the coarse shared grid snapped it to ~0.085 steps). Decoupled from the AMBIG 2-D grid (`sweep_n_grid`), which stays coarse for genome-scale memory. |
| `--gdna-em-llr-bias b` | `0.0` | gDNA false-positive-aversion: a log-odds (LLR) bias added to the gDNA component in the locus EM (`0.0` = neutral). Positive favors gDNA — it trades the gDNA→RNA leak for the FP-safe RNA→gDNA siphon. Units are nats of log-odds (e.g. `log(9) ≈ 2.20` requires ~9:1 RNA evidence before calling a fragment RNA). Reaches the EM directly, distinct from the calibration deconvolution. |

---

## YAML config key reference

YAML keys use underscores; hyphens are also accepted.
Unknown keys generate a warning and are ignored.

```yaml
# All rigel quant keys with their defaults

# Library
include_multimap: true
keep_duplicates: false
sj_strand_tag: [auto]          # [XS] for STAR, [ts] for minimap2
splicing_anchor_tolerance: 3   # resolver-side implicit-splice slack (bp)

# EM algorithm
em_iterations: 1000
em_mode: vbem                  # vbem | map
assignment_mode: sample        # sample | fractional | map

# Performance
threads: 0                          # total thread budget; 0 = all available cores
scan_bgzf_threads: 4                # BGZF decompression threads within that budget
scan_buffer_size: 2                 # GiB before scan-buffer spill
scan_fragments_per_chunk: 1000000   # buffered fragments per native scan chunk
scan_read_name_batch_size: 512      # read-name groups per scanner queue item
seed: null                          # null = use current timestamp
tmpdir: null                        # null = system temp directory

# Advanced scoring
overhang_alpha: 0.1
mismatch_alpha: 0.1
pruning_min_posterior: 0.0001

# Advanced EM
assignment_min_posterior: 0.01
em_convergence_delta: 0.000001
gdna_em_llr_bias: 0.0          # nats of log-odds; positive = gDNA-averse siphon

# Advanced calibration
gdna_prior_mixture_bridge: 0.01     # KDE valley bridge weight ε ∈ [0, 1)
sweep_n_grid_single_strand: 256     # single-strand log-odds grid resolution
```

---

## Config dataclass reference

All defaults are defined in `src/rigel/config.py` as frozen dataclasses.
The `_PARAM_SPECS` registry in `cli.py` maps CLI flag names to config fields.

| Dataclass | Purpose |
|-----------|---------|
| `PipelineConfig` | Top-level container: `em`, `scan`, `scoring`, `calibration`, plus `annotated_bam_path` and `emit_locus_stats` |
| `EMConfig` | EM algorithm: mode, iterations, convergence, assignment behavior, `gdna_em_llr_bias`, thread count |
| `BamScanConfig` | BAM parsing, filtering, buffering, threading, resolver splicing-anchor tolerance |
| `FragmentScoringConfig` | Overhang/mismatch log-penalties, optional per-SpliceType gDNA penalties, candidate pruning |
| `CalibrationConfig` | Belief-propagation sweep: strand-overdispersion priors, sweep grid resolutions, Phase-2 gDNA-density KDE prior knobs |

### Selected `BamScanConfig` fields

| Field | Default | CLI/YAML key | Description |
|-------|---------|--------------|-------------|
| `max_frag_length` | `1000` | — | Max fragment length (bp) for histogram models |
| `total_threads` | `0` | `threads` | Total scan-stage thread budget (`0` = all cores) |
| `bgzf_threads` | `4` | `scan_bgzf_threads` | BGZF decompression threads reserved from the budget |
| `fragments_per_chunk` | `1000000` | `scan_fragments_per_chunk` | Buffered fragments per native scan chunk |
| `read_name_batch_size` | `512` | `scan_read_name_batch_size` | Read-name groups per native scanner queue item |
| `buffer_size_bytes` | 2 GiB | `scan_buffer_size` | Scan-buffer memory limit before disk spill |
| `splicing_anchor_tolerance` | `3` | `splicing_anchor_tolerance` | Resolver implicit-splice slack (bp) |
| `log_every` | `1000000` | — | Progress log interval (read-name groups) |

### Selected `CalibrationConfig` fields

All are advanced; the defaults are validated. The sweep is a **single forward-backward pass**
(there is no iteration loop — the old `sweep_max_passes` / `sweep_convergence_delta` fields were
removed).

| Field | Default | CLI/YAML key | Description |
|-------|---------|--------------|-------------|
| `gdna_strand_prior_alpha_beta` | `14.0` | — | Symmetric `Beta(a, a)` shape for the gDNA strand-overdispersion prior (near-binomial null; `a ≥ 2`). |
| `gdna_strand_prior_weight` | `30.0` | — | Overdispersion prior strength in effective seed-node units. |
| `rna_strand_prior_alpha_beta` | `14.0` | — | Twin of the gDNA prior for the RNA strand Beta-Binomial (kept equal so unstranded data stays uninformative). |
| `rna_strand_prior_weight` | `30.0` | — | RNA overdispersion prior strength (twin of the gDNA weight). |
| `sweep_n_grid` | `60` | — | AMBIG per-node log-odds solve grid resolution `K`. |
| `sweep_n_grid_single_strand` | `256` | `sweep_n_grid_single_strand` | Single-strand 1-D log-odds grid resolution `K_ss`. |
| `sweep_logodds_window` | `10.0` | — | Log-odds grid half-window `L`: `λ ∈ [−L, L]`. |
| `sweep_n_tilt` | `None` | — | Inner tilt-grid resolution `K_t` for the AMBIG 2-D `(λ, τ)` solve; `None` reuses `sweep_n_grid`. |
| `gdna_prior_bandwidth` | `"silverman"` | — | Phase-2 KDE bandwidth selector: `"silverman"`, `"lscv"`, or a fixed float. |
| `gdna_prior_mixture_bridge` | `0.01` | `gdna_prior_mixture_bridge` | KDE valley mixture-bridge weight ε ∈ [0, 1). |
| `calib_kde_min_training_nodes` | `10` | — | Minimum gDNA-density teacher nodes required to fit the Phase-2 KDE; below it, calibration falls back to the pass-1 prior-free belief. |
| `calib_kde_bridge_trim_pct` | `0.5` | — | Percentile trim bounding the mixture-bridge to the training density support. |

### `summary.json` calibration schema

`summary.json` records the run's library scalars and QC. The `calibration` section is built
directly in `cli.py` from the `CalibrationResult` library scalars:

| Key | Type | Description |
|-----|------|-------------|
| `gdna_density_global` | float | Library-average gDNA density (a QC scalar) from the deconvolved masses. |
| `rna_sense_frac` | float | RNA sense fraction κ (`rna_sense_frac`), the fitted RNA strand mean. |
| `gdna_strand_overdispersion` | float | Fitted gDNA strand Beta-Binomial overdispersion ∈ [0, 1). |
| `rna_strand_overdispersion` | float | Fitted RNA strand Beta-Binomial overdispersion ∈ [0, 1). |
| `n_regions` | int | Number of calibration regions in the partition. |

The per-pool fragment-length models are emitted separately under the top-level
**`fragment_length`** section of `summary.json` (not nested under `calibration`): a `global` pool
(the scanner's raw histogram), plus the empirical `gdna` (intergenic + intronic structural pool)
and `rna` (spliced-annotated) FL models actually used for scoring, alongside the scanner's
per-splice-category histograms.

> The retired v5/v6 keys (`region_calibration`, `background_model`, `boundary_local`,
> `boundary_sweep`, `prior_table`, `strand_channels`, `diagnostics`) are **no longer emitted** —
> those calibrators have been superseded by the belief-propagation sweep.

---

## C++ compile-time constants

These constants are defined in `src/rigel/native/em_solver.cpp` and cannot be
changed at runtime. They control numerical behavior of the EM solver.

| Constant | Value | Description |
|----------|-------|-------------|
| `EM_LOG_EPSILON` | `1e-300` | Global floor for log-space operations and zero-valued priors |
| `SQUAREM_BUDGET_DIVISOR` | `3` | SQUAREM iteration budget = `em_iterations / 3` (each SQUAREM step costs 3 EM iterations) |
| `VBEM_CLAMP_FLOOR` | `0.1` | Minimum Dirichlet alpha after SQUAREM extrapolation or stabilization in VBEM mode. Prevents components from entering the digamma absorbing regime where `ψ(α) ≈ −1/α` makes recovery impossible. Has no effect on MAP-EM. |
| `ESTEP_TASK_WORK_TARGET` | `4096` | Target element-operations per parallel E-step task for load-balanced threading |
