# hulkrna Parameters Reference

All parameters below apply to the `hulkrna quant` subcommand.
Defaults shown in parentheses are the effective defaults from the config
dataclasses in `config.py`.  Every CLI flag defaults to `None` internally;
the final value is resolved as: **explicit CLI flag → YAML config file → dataclass default**.

---

## Required Arguments

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted or collated BAM file (must have NH tag). |
| `--index` | Directory containing hulkrna index files. |
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

### nRNA Fraction Hierarchical Shrinkage

The nRNA fraction prior uses a three-level empirical Bayes hierarchy:
**global → locus-strand → TSS-group → transcript**.  At each level a
shrinkage pseudo-count κ pulls the estimate toward the parent level.
By default κ is auto-estimated via Method of Moments (MoM).

| Flag | Default | Description |
|------|---------|-------------|
| `--tss-window` | 200 | Fuzzy TSS grouping window (bp). Transcripts whose 5′ ends lie within this distance are grouped for the nRNA fraction prior. |
| `--nrna-frac-kappa-global` | auto (MoM) | κ pulling locus-strand nRNA fraction toward the global prior. |
| `--nrna-frac-kappa-locus` | auto (MoM) | κ pulling TSS-group nRNA fraction toward the locus-strand estimate. |
| `--nrna-frac-kappa-tss` | auto (MoM) | κ pulling transcript nRNA fraction toward the TSS-group estimate. Also sets the effective sample size of the final Beta prior passed to the EM. |
| `--nrna-frac-kappa-min` | 2.0 | Lower clamp for MoM-estimated κ. |
| `--nrna-frac-kappa-max` | 200.0 | Upper clamp for MoM-estimated κ. |
| `--nrna-frac-kappa-fallback` | 5.0 | Fallback κ when too few features pass the evidence filter. |
| `--nrna-frac-kappa-min-obs` | 20 | Minimum features required for MoM κ estimation; fewer triggers the fallback. |

#### MoM Minimum Evidence Thresholds

At each hierarchy level, groups with fewer fragments than the threshold
are excluded from the Method-of-Moments κ estimation.

| Flag | Default | Description |
|------|---------|-------------|
| `--nrna-frac-mom-min-evidence-global` | 50 | Min fragment evidence for global MoM κ. |
| `--nrna-frac-mom-min-evidence-locus` | 30 | Min fragment evidence for locus MoM κ. |
| `--nrna-frac-mom-min-evidence-tss` | 20 | Min fragment evidence for TSS-level MoM κ. |

### gDNA Rate Hierarchical Shrinkage

The gDNA rate prior uses a two-level hierarchy:
**global → chromosome → locus**.

| Flag | Default | Description |
|------|---------|-------------|
| `--gdna-kappa-chrom` | auto (MoM) | κ pulling chromosome gDNA rate toward the global estimate. |
| `--gdna-kappa-locus` | auto (MoM) | κ pulling locus gDNA rate toward the chromosome estimate. |
| `--gdna-mom-min-evidence-chrom` | 50 | Min fragment evidence for chromosome gDNA MoM κ. |
| `--gdna-mom-min-evidence-locus` | 30 | Min fragment evidence for locus gDNA MoM κ. |

### Strand Model

| Flag | Default | Description |
|------|---------|-------------|
| `--strand-prior-kappa` | 2.0 | Strand model prior pseudocount κ. The Beta prior is Beta(κ/2, κ/2), shrinking toward 0.5 (max entropy). Default 2.0 gives a uniform Beta(1, 1) prior. |

---

## Config Dataclass Reference

All defaults are defined in `src/hulkrna/config.py` as frozen dataclasses.
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
| `EMConfig` | `mode` | `"map"` | Algorithm variant: `"map"` or `"vbem"`. |
| `BamScanConfig` | `max_frag_length` | 1000 | Maximum fragment length for histogram models (bp). |
| `BamScanConfig` | `log_every` | 1,000,000 | Log progress every N read-name groups. |
| `BamScanConfig` | `chunk_size` | 1,000,000 | Fragments per buffer chunk. |
| `BamScanConfig` | `max_memory_bytes` | 2 GiB | Max memory before disk spill. |
| `BamScanConfig` | `spill_dir` | None | Directory for spilled buffer chunks (set via `--tmpdir`). |
