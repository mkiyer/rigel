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
| `--alignable-zarr PATH` | — | Alignable Zarr store. Provides fractional mappability (for gDNA calibration) and the splice-artifact blacklist (applied at BAM-scan time). Required unless `--no-mappability` is set. |
| `--no-mappability` | off | Opt out of mappability and the splice-artifact blacklist. Mutually exclusive with `--alignable-zarr`. |
| `--splice-blacklist-min-count` | `2` | Minimum unique-fragment support per `(chrom, intron, read_length)` for a junction to enter the blacklist. Ignored under `--no-mappability`. |
| `--mappability-read-length` | `100` | Read-length bin used when querying the alignable store. |
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
| `--splicing-anchor-tolerance K` | `3` | Shared splice-anchor slack in bp. Used by implicit-splice classification for annotated introns in paired-end gaps and by calibration boundary-flux clearance. |

### EM algorithm

| Flag | Default | Description |
|------|---------|-------------|
| `--em-iterations` | `1000` | Maximum EM iterations. `0` = skip EM entirely (unambiguous fragments only). |
| `--assignment-mode` | `sample` | Post-EM fragment assignment mode: `sample` (draw from posterior distribution, default), `fractional` (preserve EM weights), `map` (assign to argmax component). |
| `--em-mode` | `vbem` | EM algorithm variant: `vbem` (Variational Bayes EM, default) or `map` (MAP-EM). |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--tsv` | off | Write TSV mirrors alongside Feather files |
| `--annotated-bam PATH` | — | Write annotated BAM with per-fragment assignment tags (`ZT`, `ZG`, `ZR`, `ZI`, `ZJ`, `ZF`, `ZW`, `ZC`, `ZH`, `ZN`, `ZS`, `ZL`, `ZB`). Rigel guarantees the output is collated and contains exactly the same records as the input. Requires a second BAM pass. |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | `0` | Total thread budget available to Rigel. During BAM scan, Rigel reserves `--scan-bgzf-threads` from this budget for BAM decompression and uses the remainder for scan workers. Locus EM uses the same budget because scan and EM run serially. `0` = all available cores, `1` = sequential. |
| `--scan-bgzf-threads N` | `4` | BGZF/BAM decompression threads reserved from `--threads` during the scan phase. Set to `0` to disable htslib threaded decompression. If the requested value would consume the full thread budget, Rigel caps it so at least one scan worker remains. |
| `--scan-buffer-size GB` | `2` | Maximum scan buffer size in GiB before finalized chunks are spilled to disk. Increase if you have ample RAM and want to reduce spill I/O. |
| `--scan-fragments-per-chunk N` | `1000000` | Buffered fragments per native scan chunk before data is handed to the Python fragment buffer. Larger chunks reduce handoff overhead but increase transient memory. |
| `--scan-read-name-batch-size N` | `512` | Advanced scanner scheduling knob: read-name groups per native worker-queue item. Lower values are useful for queue-boundary debugging; larger values may reduce queue overhead at the cost of slightly more in-flight memory. |
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
| `--overhang-alpha` | `0.1` | Per-base overhang penalty α ∈ (0, 1]. Fragment score multiplied by α for each overhang base. `0` = hard gate, `1` = no penalty. |
| `--mismatch-alpha` | `0.1` | Per-mismatch (NM tag) penalty α ∈ (0, 1]. Score multiplied by α per mismatch. `0` = hard gate, `1` = no penalty. |
| `--gdna-splice-penalty-unannot` | `0.01` | Multiplier applied to the gDNA candidate score for fragments with unannotated splice junctions. Values close to 0 make gDNA attribution less likely for spliced fragments. |
| `--cal-prior-ess` | `1000.0` | Empirical-Bayes evidence strength for the FL-Dirichlet shrinkage in the calibration orchestrator. Larger values shrink RNA/gDNA fragment-length distributions more aggressively toward the global FL. |
| `--cal-fl-scoring-prior-ess` | `200.0` | Effective sample size for joint RNA-vs-gDNA fragment-length contrast reliability during EM scoring. Larger values require more support in both FL pools before class-specific FL differences affect posterior odds; `0` disables weak-pool contrast damping while preserving fallback neutrality. |
| `--cal-quality-good` | `5000` | Minimum SPLICED-annotated count (rna pool) and gDNA count above which the per-FL distribution is flagged `"good"` and used without shrinkage. |
| `--cal-quality-weak` | `1` | Minimum SPLICED-annotated count and gDNA count above which the per-FL distribution is flagged `"weak"`. Below this, the pool is `"fallback"` and the raw pool model identity-shares the global FL. |

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
splicing_anchor_tolerance: 3   # shared implicit-splice and boundary-flux slack

# EM algorithm
em_iterations: 1000
assignment_mode: sample        # sample | fractional | map
em_mode: vbem                  # vbem | map

# Performance
threads: 0                          # total thread budget; 0 = all available cores
scan_bgzf_threads: 4                # BGZF decompression threads within that budget
scan_buffer_size: 2                 # GiB before scan-buffer spill
scan_fragments_per_chunk: 1000000   # buffered fragments per native scan chunk
scan_read_name_batch_size: 512      # read-name groups per scanner queue item
seed: null                     # null = use current timestamp
tmpdir: null                   # null = system temp directory

# Advanced scoring
overhang_alpha: 0.1
mismatch_alpha: 0.1
gdna_splice_penalty_unannot: 0.01
pruning_min_posterior: 0.0001

# Advanced EM
assignment_min_posterior: 0.01
em_convergence_delta: 0.000001

# Calibration (v6)
cal_prior_ess: 1000.0
cal_fl_scoring_prior_ess: 200.0
cal_quality_good: 5000
cal_quality_weak: 1
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
| `CalibrationConfig` | v6 calibration orchestrator: FL-pool quality thresholds, FL shrinkage/scoring reliability, regional gDNA density calibration |

### Selected `BamScanConfig` fields

| Dataclass | Field | Default | Description |
|-----------|-------|---------|-------------|
| `BamScanConfig` | `max_frag_length` | `1000` | Max fragment length (bp) for histogram models |
| `BamScanConfig` | `total_threads` | `0` | Total scan-stage thread budget; CLI/YAML key `threads` |
| `BamScanConfig` | `bgzf_threads` | `4` | BGZF decompression threads reserved from the scan budget |
| `BamScanConfig` | `fragments_per_chunk` | `1,000,000` | Buffered fragments per native scan chunk |
| `BamScanConfig` | `read_name_batch_size` | `512` | Read-name groups per native scanner queue item |
| `BamScanConfig` | `buffer_size_bytes` | 2 GiB | Scan-buffer memory limit before disk spill |
| `BamScanConfig` | `log_every` | `1,000,000` | Progress log interval (read-name groups) |

### `summary.json` calibration schema

The `calibration` section of `summary.json` is built from
`CalibrationResult.to_summary_dict()` and contains:

| Key | Type | Source | Description |
|-----|------|--------|-------------|
| `fl_models` | object | `FLModels.to_summary_dict()` | Per-pool raw fragment-length models (`rna`, `gdna`, `global`), score-side FL surfaces, quality flags, support totals, and RNA-vs-gDNA FL contrast reliability diagnostics. |
| `diagnostics` | object | `Diagnostics.to_summary_dict()` | Named per-region observation counts (intergenic / intron / exon / boundary-flux). |
| `n_multi_loci` | int | `CalibrationResult.n_multi_loci` | Number of MultiLoci that received a per-component prior. |
| `region_calibration` | object | `RegionCalibration` summary | Regional four-state posterior calibration, including state mass summaries, `rho_off`, `A_r`, pass count, and convergence status. |
| `background_model` | object | `BackgroundModel` summary | Off-target background model diagnostics used by regional calibration. |
| `boundary_local` | object | `BoundaryLocalPosterior` summary | Local boundary-excess posterior diagnostics. |
| `boundary_sweep` | object | `BoundarySweepResult` summary | Boundary-transfer sweep diagnostics. |
| `strand_channels` | object or null | `RegionGdnaChannelEstimate` summary | Strand-aware contained/boundary gDNA channel estimates when strand contrast is identifiable. |
| `prior_table` | object or null | `PriorTable.to_summary_dict()` | Adaptive grouped RNA/gDNA EM prior diagnostics populated after quantification. |

---

## C++ compile-time constants

These constants are defined in `src/rigel/native/em_solver.cpp` and cannot be
changed at runtime. They control numerical behavior of the EM solver.

| Constant | Value | Description |
|----------|-------|-------------|
| `EM_LOG_EPSILON` | `1e-300` | Global floor for log-space operations and zero-valued priors |
| `SQUAREM_BUDGET_DIVISOR` | `3` | SQUAREM iteration budget = `em_iterations / 3` (each SQUAREM step costs 3 EM iterations) |
| `VBEM_CLAMP_FLOOR` | `0.1` | Minimum Dirichlet alpha after SQUAREM extrapolation or stabilization in VBEM mode. Prevents components from entering the digamma absorbing regime where `ψ(α) ≈ −1/α` makes recovery impossible. At 0.1, `ψ(0.1) ≈ −10.4` gives relative weight ≈ 3×10⁻⁵ — small enough not to steal mass from other components, but large enough that a component with genuine read support can accumulate evidence and recover. Has no effect on MAP-EM. |
| `ESTEP_TASK_WORK_TARGET` | `4096` | Target element-operations per parallel E-step task for load-balanced threading |
