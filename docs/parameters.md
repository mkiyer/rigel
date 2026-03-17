# rigel Parameters Reference

All parameters below apply to the `rigel quant` subcommand.
Defaults shown in parentheses are the effective defaults from the config
dataclasses in `config.py`.  Every CLI flag defaults to `None` internally;
the final value is resolved as: **explicit CLI flag → YAML config file → dataclass default**.

---

## Required Arguments

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted or collated BAM file (must have NH tag). |
| `--index` | Directory containing rigel index files. |
| `-o`, `--output-dir` | Output directory for quantification results and models. |

## Configuration File

| Flag | Default | Description |
|------|---------|-------------|
| `--config` | — | YAML configuration file. Any key matching a CLI option (using underscores) may be set; explicit CLI flags override the file. |

---

## Common Options

### Read Filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--include-multimap` / `--no-include-multimap` | yes | Include multimapping reads. Detected via NH tag (STAR) or secondary BAM flag (minimap2). |
| `--keep-duplicates` / `--no-keep-duplicates` | no | Keep reads marked as PCR/optical duplicates. |

### Splice-Junction Strand

| Flag | Default | Description |
|------|---------|-------------|
| `--sj-strand-tag` | `auto` | BAM tag(s) for splice-junction strand. `auto` detects the tag from the BAM header. Use `XS` for STAR, `ts` for minimap2, or list multiple to check in order (e.g. `XS ts`). |

### EM Algorithm

| Flag | Default | Description |
|------|---------|-------------|
| `--em-iterations` | 1000 | Maximum EM iterations for ambiguous fragment resolution. Set to 0 for unambiguous-only quantification. |
| `--em-prior-alpha` | 0.01 | Flat Dirichlet pseudocount per eligible EM component. |
| `--em-prior-gamma` | 1.0 | OVR (One Virtual Read) prior scale factor. Controls how strongly coverage-weighted OVR priors influence the EM. 0.0 disables OVR. |
| `--em-convergence-delta` | 1e-6 | Convergence threshold for EM parameter updates. |
| `--confidence-threshold` | 0.95 | Minimum RNA-normalized posterior for high-confidence assignment. |
| `--prune-threshold` | 0.1 | Post-EM pruning evidence-ratio threshold. Components with zero unambiguous evidence and an evidence ratio (data\_count / alpha) below this value are zeroed out and the EM re-runs to redistribute mass. Set to -1 to disable. |
| `--em-mode` | `vbem` | EM algorithm variant. `map` uses MAP-EM with hard max(0, n+α−1) updates; `vbem` uses Variational Bayes EM with digamma-based soft updates. |
| `--nrna-sparsity-alpha` | 0.9 | Dirichlet α for nRNA EM components. Values < 1.0 create a sparsifying prior that drives weak nRNA components toward zero. 1.0 = flat (off), 0.9 = mild, 0.5 = moderate, 0.3 = strong. Does not affect mRNA accuracy. |
| `--gdna-prior-scale` | 1.0 | Scale factor for the Empirical Bayes gDNA prior anchor. The gDNA prior is α = 1 + scale × gdna\_init, where gdna\_init is the EB estimate. 0 = ignore EB (flat unit prior), 1.0 = default, higher = stronger anchor. |

### Scoring Penalties

| Flag | Default | Description |
|------|---------|-------------|
| `--overhang-alpha` | 0.01 | Per-base overhang penalty α ∈ [0, 1]. Each overhang base multiplies the score by α. 0 = hard gate, 1 = no penalty. |
| `--mismatch-alpha` | 0.1 | Per-mismatch (NM tag) penalty α ∈ [0, 1]. Each mismatch multiplies the score by α. 0 = hard gate, 1 = no penalty. |
| `--gdna-splice-penalty-unannot` | 0.01 | gDNA score penalty for `SPLICED_UNANNOT` fragments (unannotated splice junctions). |

### Output

| Flag | Default | Description |
|------|---------|-------------|
| `--no-tsv` | off | Skip writing human-readable TSV quantification files. |
| `--annotated-bam` | — | Write an annotated BAM with per-fragment assignment tags (ZT, ZG, ZP, ZW, ZC, ZH, ZN, ZS) to this path. Requires a second BAM pass. |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads` | 0 (all cores) | Number of threads for BAM scanning and locus EM (these stages run serially). 0 = all available cores, 1 = sequential. |
| `--tmpdir` | system temp | Directory for temporary buffer spill files when memory limits are exceeded. Use a fast SSD mount in production. |
| `--seed` | timestamp | Random seed for reproducibility. |

---

## Advanced Options

### gDNA Rate Hierarchical Shrinkage

The gDNA rate prior uses a two-level hierarchy:
**global → reference → locus**.

| Flag | Default | Description |
|------|---------|-------------|
| `--gdna-kappa-ref` | auto (MoM) | κ pulling reference gDNA rate toward the global estimate. |
| `--gdna-kappa-locus` | auto (MoM) | κ pulling locus gDNA rate toward the reference estimate. |
| `--gdna-mom-min-evidence-ref` | 50 | Min fragment evidence for reference gDNA MoM κ. |
| `--gdna-mom-min-evidence-locus` | 30 | Min fragment evidence for locus gDNA MoM κ. |
| `--gdna-kappa-min` | 2.0 | Lower clamp for MoM-estimated gDNA κ. |
| `--gdna-kappa-max` | 200.0 | Upper clamp for MoM-estimated gDNA κ. |
| `--gdna-kappa-fallback` | 5.0 | Fallback gDNA κ when too few features for MoM. |
| `--gdna-kappa-min-obs` | 20 | Minimum features for gDNA MoM κ estimation; fewer triggers fallback. |

### Strand Model

| Flag | Default | Description |
|------|---------|-------------|
| `--strand-prior-kappa` | 2.0 | Strand model prior pseudocount κ. The Beta prior is Beta(κ/2, κ/2), shrinking toward 0.5 (max entropy). Default 2.0 gives a uniform Beta(1, 1) prior. |

---

## Config Dataclass Reference

All defaults are defined in `src/rigel/config.py` as frozen dataclasses.
The CLI-to-config mapping in `cli.py` translates flag names to config fields
(e.g. `--em-iterations` → `EMConfig.iterations`).

| Dataclass | Purpose |
|-----------|---------|
| `PipelineConfig` | Top-level container holding `em`, `scan`, and `scoring` sub-configs plus `annotated_bam_path`. |
| `EMConfig` | EM algorithm, nRNA fraction shrinkage, and gDNA shrinkage parameters. |
| `BamScanConfig` | BAM reading, filtering, buffering, and threading. |
| `FragmentScoringConfig` | Overhang/mismatch log-penalties and gDNA splice penalties. |

### Internal-only config fields

These fields are set programmatically and not exposed as CLI flags:

| Dataclass | Field | Default | Description |
|-----------|-------|---------|-------------|
| `EMConfig` | `mode` | `"vbem"` | Algorithm variant: `"map"` or `"vbem"`. Exposed via `--em-mode`. |
| `BamScanConfig` | `max_frag_length` | 1000 | Maximum fragment length for histogram models (bp). |
| `BamScanConfig` | `log_every` | 1,000,000 | Log progress every N read-name groups. |
| `BamScanConfig` | `chunk_size` | 1,000,000 | Fragments per buffer chunk. |
| `BamScanConfig` | `max_memory_bytes` | 2 GiB | Max memory before disk spill. |
| `BamScanConfig` | `spill_dir` | None | Directory for spilled buffer chunks (set via `--tmpdir`). |

---

## Internal Constants Reference

The following constants are used internally by the C++ scoring, EM, and BAM
scanning kernels.  They are defined in C++ header/source files and exported
to Python via nanobind module attributes for transparency and parity testing.
They are **not** user-configurable — changing them requires a code edit and
rebuild.

### Enum-Mirror Constants (`constants.h`)

These mirror Python `IntEnum` definitions and must stay in sync.

#### IntervalType (`rigel.types.IntervalType` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `ITYPE_EXON` | 0 | Individual exon boundary interval. |
| `ITYPE_TRANSCRIPT` | 1 | Full transcript span `[start, end)`. |
| `ITYPE_UNAMBIG_INTRON` | 5 | Intronic region not overlapping any other transcript's exon. |

#### SpliceType (`rigel.splice.SpliceType` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `SPLICE_UNSPLICED` | 0 | Fragment has no splice junctions. |
| `SPLICE_SPLICED_UNANNOT` | 1 | Fragment has a splice junction not matching any annotated intron. |
| `SPLICE_SPLICED_ANNOT` | 2 | Fragment has a splice junction matching an annotated intron. |

#### FragmentClass (`rigel.buffer` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `FRAG_UNAMBIG` | 0 | Same-strand, 1 transcript, NH=1. |
| `FRAG_AMBIG_SAME_STRAND` | 1 | Same-strand, >1 transcript, NH=1. |
| `FRAG_AMBIG_OPP_STRAND` | 2 | Ambiguous-strand transcripts, NH=1. |
| `FRAG_MULTIMAPPER` | 3 | NH > 1 (multimapped molecule). |
| `FRAG_CHIMERIC` | 4 | Chimeric fragment (disjoint transcript sets). |

#### PoolCode (`rigel.annotate` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `POOL_CODE_MRNA` | 0 | Mature mRNA pool. |
| `POOL_CODE_NRNA` | 1 | Nascent RNA (pre-mRNA) pool. |
| `POOL_CODE_GDNA` | 2 | Genomic DNA pool. |
| `POOL_CODE_INTERGENIC` | 3 | Intergenic (no gene overlap). |
| `POOL_CODE_CHIMERIC` | 4 | Chimeric fragment pool. |

#### MergeOutcome (`rigel.types.MergeOutcome` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `MC_INTERSECTION` | 0 | Mate-pair transcript sets intersect. |
| `MC_INTERSECTION_NONEMPTY` | 1 | Non-empty intersection after filtering. |
| `MC_UNION` | 2 | Union of mate-pair transcript sets. |
| `MC_EMPTY` | 3 | No shared transcripts. |

#### ChimeraType (`rigel.types.ChimeraType` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `CHIMERA_NONE` | 0 | Not chimeric. |
| `CHIMERA_TRANS` | 1 | Trans-chromosomal chimera. |
| `CHIMERA_CIS_STRAND_SAME` | 2 | Cis-chromosomal chimera, same strand. |
| `CHIMERA_CIS_STRAND_DIFF` | 3 | Cis-chromosomal chimera, different strand. |

#### Strand (`rigel.types.Strand` → `constants.h`)

| Constant | Value | Description |
|----------|-------|-------------|
| `STRAND_NONE` | 0 | No strand information. |
| `STRAND_POS` | 1 | Positive (`pos`, `+`) strand. |
| `STRAND_NEG` | 2 | Negative (`neg`, `-`) strand. |
| `STRAND_AMBIGUOUS` | 3 | Ambiguous (POS \| NEG). |

### Scoring Constants (`constants.h` → `_scoring_impl`)

| Constant | Value | Source | Description |
|----------|-------|--------|-------------|
| `LOG_HALF` | −0.6931… | `log(0.5)` | Uninformative strand log-probability. |
| `TAIL_DECAY_LP` | −0.01005… | `log(0.99)` | Per-base tail decay for fragment length model beyond max_frag_length. |
| `SCORED_STACK_CAPACITY` | 64 | `scoring.cpp` | Stack-allocated buffer size for transcript candidates before heap fallback. |

### EM Solver Constants (`em_solver.cpp` → `_em_impl`)

| Constant | Value | Description |
|----------|-------|-------------|
| `EM_LOG_EPSILON` | 1e-300 | Floor for log-domain clamping to avoid log(0). |
| `MAX_FRAG_LEN` | 1,000,000 | Internal maximum fragment length (overflow bin). |
| `SQUAREM_BUDGET_DIVISOR` | 3 | SQUAREM acceleration: budget = max_iter / divisor. |
| `NRNA_FRAC_CLAMP_EPS` | 1e-8 | Floor/ceiling epsilon for nRNA fraction clamping. |
| `EM_PRIOR_EPSILON` | 1e-10 | Numerical-stability floor for prior components. |
| `ESTEP_TASK_WORK_TARGET` | 4,096 | Target element-ops per E-step parallel task for load balancing. |

### Python-Side Scoring Defaults (`scoring.py`)

These are used as defaults when no user override is provided:

| Constant | Value | Description |
|----------|-------|-------------|
| `LOG_SAFE_FLOOR` | 1e-10 | Floor for log-safe clamping. |
| `DEFAULT_GDNA_SPLICE_PENALTY_UNANNOT` | 0.01 | gDNA penalty for unannotated splice junctions. |
| `DEFAULT_OVERHANG_ALPHA` | 0.01 | Per-base overhang penalty α (maps to `--overhang-alpha`). |
| `DEFAULT_MISMATCH_ALPHA` | 0.1 | Per-mismatch penalty α (maps to `--mismatch-alpha`). |

### Constant Flow: Python ↔ C++

Constants defined in Python (`rigel.types`, `rigel.splice`, `rigel.buffer`,
`rigel.annotate`) are mirrored as `static constexpr` in `constants.h` for
compile-time use in C++ kernels.  The C++ values are exported back to Python
via nanobind module attributes (`_scoring_impl`, `_em_impl`, `_resolve_impl`)
enabling parity tests.

| Python Source | C++ Mirror | Nanobind Module |
|---------------|-----------|-----------------|
| `rigel.types.IntervalType` | `constants.h` | `_resolve_impl` |
| `rigel.types.Strand` | `constants.h` | — |
| `rigel.types.MergeOutcome` | `constants.h` | — |
| `rigel.types.ChimeraType` | `constants.h` | — |
| `rigel.splice.SpliceType` | `constants.h` | `_resolve_impl`, `_scoring_impl` (via `using`) |
| `rigel.buffer.FRAG_*` | `constants.h` | `_scoring_impl` |
| `rigel.annotate.POOL_CODE_*` | `constants.h` | `_scoring_impl` |
| `rigel.scoring.LOG_HALF` | `constants.h` | `_scoring_impl` |
| `rigel.frag_length_model._TAIL_DECAY_LP` | `constants.h` | `_scoring_impl` |
| — (C++ only) | `em_solver.cpp` | `_em_impl` |
