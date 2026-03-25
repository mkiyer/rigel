# Rigel Parameters Reference

All parameters below apply to `rigel quant` unless noted.
Resolution order: **explicit CLI flag → YAML config file → built-in default**.

---

## rigel index

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Genome FASTA (must have `.fai` index) |
| `--gtf` | required | Annotation GTF |
| `-o`, `--output-dir` | required | Output directory |
| `--nrna-tolerance` | `20` | Max distance (bp) for clustering transcript start/end sites into shared nRNA spans |
| `--gtf-parse-mode` | `strict` | `strict` fails on malformed GTF records; `warn-skip` logs and skips them |
| `--feather-compression` | `lz4` | Feather file compression: `lz4`, `zstd`, or `uncompressed` |
| `--no-tsv` | off | Skip writing TSV mirrors of index files |

---

## rigel quant

### Required

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted BAM file. Minimap2 and STAR are supported. NH tag used when present; otherwise multimappers are detected from secondary BAM flags. |
| `--index` | Rigel index directory (from `rigel index`) |
| `-o`, `--output-dir` | Output directory |

### Library and input

| Flag | Default | Description |
|------|---------|-------------|
| `--config FILE` | — | YAML file with parameter values; CLI overrides YAML |
| `--include-multimap` / `--no-include-multimap` | yes | Include multimapping reads (detected via `NH` tag or secondary flag) |
| `--keep-duplicates` / `--no-keep-duplicates` | no | Keep PCR/optical duplicate-marked reads |
| `--sj-strand-tag TAG [TAG …]` | `auto` | Splice-junction strand BAM tag(s). `auto` detects from the BAM header. Common values: `XS` (STAR, HISAT2), `ts` (minimap2). Multiple tags tried in order. |

### EM algorithm

| Flag | Default | Description |
|------|---------|-------------|
| `--prior-pseudocount` | `1.0` | Total OVR prior budget C (virtual fragments). Distributed as γ×C to gDNA and (1−γ)×C coverage-weighted across RNA components. Increase to strengthen priors in low-count loci. |
| `--em-iterations` | `1000` | Maximum EM iterations. `0` = skip EM entirely (unambiguous fragments only). |
| `--assignment-mode` | `sample` | Post-EM fragment assignment mode: `sample` (draw from posterior distribution, default), `fractional` (preserve EM weights), `map` (assign to argmax component). |
| `--em-mode` | `vbem` | EM algorithm variant: `vbem` (Variational Bayes EM, default) or `map` (MAP-EM). |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--tsv` | off | Write TSV mirrors alongside Feather files |
| `--annotated-bam PATH` | — | Write annotated BAM with per-fragment assignment tags (`ZT`, `ZG`, `ZI`, `ZJ`, `ZP`, `ZW`, `ZC`, `ZH`, `ZN`, `ZS`, `ZL`). Requires a second BAM pass. |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | `0` | Thread count for BAM scanning and locus EM (these stages run serially). `0` = all available cores, `1` = sequential. |
| `--tmpdir DIR` | system temp | Directory for fragment-buffer spill files when memory is exceeded. Use a fast local SSD. |
| `--seed N` | timestamp | Random seed for reproducibility. Affects post-EM `sample` assignment. |

---

## Advanced options

These control scoring sensitivity and EM numerical behavior. The defaults are
suitable for standard total RNA-seq libraries.

| Flag | Default | Description |
|------|---------|-------------|
| `--assignment-min-posterior` | `0.01` | Minimum component posterior for discrete assignment (`map`/`sample` modes only). Components below this threshold are zeroed before assignment. |
| `--em-convergence-delta` | `1e-6` | EM convergence threshold for parameter updates (‖Δθ‖). |
| `--pruning-min-posterior` | `1e-4` | Remove CSR candidates with posterior below this threshold before running EM. Reduces state space for complex loci. Set to `0` to disable. |
| `--overhang-alpha` | `0.01` | Per-base overhang penalty α ∈ (0, 1]. Fragment score multiplied by α for each overhang base. `0` = hard gate, `1` = no penalty. |
| `--mismatch-alpha` | `0.1` | Per-mismatch (NM tag) penalty α ∈ (0, 1]. Score multiplied by α per mismatch. `0` = hard gate, `1` = no penalty. |
| `--gdna-splice-penalty-unannot` | `0.01` | Multiplier applied to the gDNA candidate score for fragments with unannotated splice junctions. Values close to 0 make gDNA attribution less likely for spliced fragments. |

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

# EM algorithm
prior_pseudocount: 1.0
em_iterations: 1000
assignment_mode: sample        # sample | fractional | map
em_mode: vbem                  # vbem | map

# Performance
threads: 0                     # 0 = all available cores
seed: null                     # null = use current timestamp
tmpdir: null                   # null = system temp directory

# Advanced scoring
overhang_alpha: 0.01
mismatch_alpha: 0.1
gdna_splice_penalty_unannot: 0.01
pruning_min_posterior: 0.0001

# Advanced EM
assignment_min_posterior: 0.01
em_convergence_delta: 0.000001
```

---

## Config dataclass reference

All defaults are defined in `src/rigel/config.py` as frozen dataclasses.
The `_PARAM_SPECS` registry in `cli.py` maps CLI flag names to config fields.

| Dataclass | Purpose |
|-----------|---------|
| `PipelineConfig` | Top-level container: `em`, `scan`, `scoring`, `calibration`, plus `annotated_bam_path` |
| `EMConfig` | EM algorithm: mode, iterations, priors, assignment behavior |
| `BamScanConfig` | BAM parsing, filtering, buffering, threading |
| `FragmentScoringConfig` | Overhang, mismatch, gDNA splice penalties, candidate pruning |
| `CalibrationConfig` | Aggregate-first gDNA calibration EM: convergence and data-gating thresholds |

### Internal-only fields (not exposed as CLI flags)

| Dataclass | Field | Default | Description |
|-----------|-------|---------|-------------|
| `BamScanConfig` | `max_frag_length` | `1000` | Max fragment length (bp) for histogram models |
| `BamScanConfig` | `chunk_size` | `1,000,000` | Fragments per buffer chunk |
| `BamScanConfig` | `max_memory_bytes` | 2 GiB | Fragment-buffer memory limit before disk spill |
| `BamScanConfig` | `log_every` | `1,000,000` | Progress log interval (read-name groups) |
| `BamScanConfig` | `n_decomp_threads` | `4` | htslib BGZF decompression threads. Set to `0` to disable multi-threaded decompression. Tuned independently from the worker thread count. |
| `CalibrationConfig` | `max_iterations` | `50` | Max calibration EM iterations |
| `CalibrationConfig` | `convergence_tol` | `1e-4` | Calibration EM convergence tolerance |
| `CalibrationConfig` | `density_percentile` | `0.10` | Percentile of density distribution used to seed gDNA initialization |
| `CalibrationConfig` | `min_gdna_regions` | `100` | Minimum eligible regions to run the calibration EM; below this, algebraic fallback is used |
| `CalibrationConfig` | `min_fl_ess` | `50` | Minimum effective sample size (Σγ·w) of FL observations needed to build a gDNA fragment-length model |
