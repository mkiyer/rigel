# Rigel: Algorithm Code Path

> Bayesian RNA-seq transcript quantification with mRNA / nRNA / gDNA deconvolution.

This document traces the complete execution path of `rigel quant` from CLI invocation through final output. Every module, class, and function in the production pipeline is listed with its role, language (Python or C++), and the data it consumes and produces.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Module Inventory](#2-module-inventory)
3. [Phase 0 — Index Construction](#3-phase-0--index-construction-rigel-index)
4. [Phase 1 — Configuration & Dispatch](#4-phase-1--configuration--dispatch)
5. [Phase 2 — BAM Scanning & Model Training](#5-phase-2--bam-scanning--model-training)
6. [Phase 3 — Fragment Scoring & EM Data Prep](#6-phase-3--fragment-scoring--em-data-prep)
7. [Phase 4 — Empirical Bayes Prior Computation](#7-phase-4--empirical-bayes-prior-computation)
8. [Phase 5 — Locus-Level EM Quantification](#8-phase-5--locus-level-em-quantification)
9. [Phase 6 — Output Generation](#9-phase-6--output-generation)
10. [Phase 7 — Annotated BAM (Optional)](#10-phase-7--annotated-bam-optional)
11. [Data Flow Diagram](#11-data-flow-diagram)
12. [Locus Model](#12-locus-model)

---

## 1. Architecture Overview

Rigel is a Python/C++ hybrid. Python handles orchestration, configuration, hierarchical prior computation, and I/O. C++ handles all hot-path computation:

| Layer | Language | Purpose |
|-------|----------|---------|
| CLI & config | Python | Argument parsing, YAML merge, dataclass construction |
| Pipeline orchestrator | Python | Phase sequencing, memory management, GC |
| BAM parsing & resolution | **C++** (`_bam_impl`) | htslib reader → fragment resolution → accumulator |
| Fragment scoring | **C++** (`_scoring_impl`) | Strand/length/overhang/mismatch scoring, CSR construction |
| Interval overlap | **C++** (`_cgranges_impl`) | cgranges-based genomic interval queries |
| Fragment resolution | **C++** (`_resolve_impl`) | Overlap → transcript-set resolution, chimera detection |
| EM solver | **C++** (`_em_impl`) | Equivalence class building, SQUAREM, parallel E-step |
| Prior computation | Python | 3-tier EB shrinkage for nRNA fractions and gDNA densities |
| Model training | Python | Strand model (Beta posterior), fragment-length model (histogram) |
| Output formatting | Python | DataFrames, Feather, TSV, summary JSON |

### Codebase Size

| Component | Files | LOC |
|-----------|-------|-----|
| Python core (src/rigel/*.py) | 19 | 9,096 |
| Python sim (src/rigel/sim/) | 7 | 2,875 |
| C++ native (src/rigel/native/) | 10 | 8,225 |
| C cgranges (src/rigel/cgranges/) | 4 | 1,146 |
| **Total** | **40** | **21,342** |

---

## 2. Module Inventory

### Python Modules (src/rigel/)

| Module | LOC | Role | C++ Dependency |
|--------|-----|------|----------------|
| `cli.py` | 937 | CLI parsing, parameter resolution, output writing | None |
| `pipeline.py` | 762 | Pipeline orchestration (scan → score → EM → output) | `_bam_impl` |
| `estimator.py` | 1,470 | Abundance estimator, batch EM dispatch, output DataFrames | `_em_impl` |
| `index.py` | 994 | Reference index build/load, interval tiling, C++ resolver init | `_cgranges_impl`, `_resolve_impl` |
| `locus.py` | 857 | Connected components, per-locus EM data prep, gDNA EB priors | `_em_impl` (union-find) |
| `buffer.py` | 649 | Fragment buffer with CSR chunks, spill-to-disk, Arrow I/O | `_resolve_impl` (accumulator) |
| `scan.py` | 277 | Fragment routing: buffer → scored CSR arrays for EM | `_scoring_impl` (via ctx) |
| `scoring.py` | 348 | FragmentScorer construction, C++ NativeFragmentScorer init | `_scoring_impl` |
| `strand_model.py` | 563 | Beta-Binomial strand model (train + query) | None |
| `frag_length_model.py` | 549 | Histogram fragment-length model (train + query + eff. lengths) | None |
| `config.py` | 258 | Frozen dataclass configs (EMConfig, PipelineConfig, etc.) | None |
| `annotate.py` | 324 | AnnotationTable + second-pass BAM writing | `_bam_impl` |
| `resolution.py` | 178 | `resolve_fragment()` wrapper → C++, chimera detection | `_resolve_impl` |
| `stats.py` | 160 | PipelineStats counters dataclass | None |
| `types.py` | 213 | Enums: Strand, IntervalType, ChimeraType, MergeOutcome | None |
| `splice.py` | 76 | SpliceType enum, SpliceStrandCol column indexing | None |
| `bias.py` | 116 | BiasProfile (prefix-sum coverage model) | None |
| `transcript.py` | 155 | Transcript dataclass, GTF → Transcript reader | None |
| `gtf.py` | 198 | GTFRecord line parser (streaming) | None |

### C++ Extensions (src/rigel/native/)

| Module (Python name) | File | LOC | Exports |
|----------------------|------|-----|---------|
| `_bam_impl` | bam_scanner.cpp | 1,929 | `BamScanner`, `BamAnnotationWriter`, `detect_sj_strand_tag` |
| `_scoring_impl` | scoring.cpp | 1,562 | `NativeFragmentScorer`, `compute_fragment_weight` |
| `_em_impl` | em_solver.cpp | 2,661 | `batch_locus_em`, `run_locus_em_native`, `connected_components` |
| `_resolve_impl` | resolve.cpp + resolve_context.h | 1,239 | `FragmentResolver`, `FragmentAccumulator`, `ResolvedFragment` |
| `_cgranges_impl` | cgranges_bind.cpp | ~150 | `cgranges` (interval tree) |

### Shared C++ Headers

| Header | LOC | Purpose |
|--------|-----|---------|
| `constants.h` | 300 | Enum mirrors (Strand, SpliceType, etc.), helper types |
| `resolve_context.h` | 1,091 | FragmentResolver core algorithm, cgranges integration |
| `fast_exp.h` | 163 | Cody-Waite fast exp() with ARM NEON vectorization |
| `thread_pool.h` | 106 | Persistent E-step thread pool |
| `thread_queue.h` | 65 | Bounded SPMC work queue |
| `simd_detect.h` | 200 | Compile-time SIMD capability detection |

---

## 3. Phase 0 — Index Construction (`rigel index`)

**Entry:** `cli.py:index_command()` → `index.py:TranscriptIndex.build()`

This is a one-time offline step that builds the reference index from FASTA + GTF.

```
FASTA + GTF
    │
    ├─ load_reference_lengths()         → ref_lengths.feather
    │      (pysam FASTA index)
    │
    ├─ read_transcripts() / Transcript.read_gtf()
    │      (streaming GTF parse → Transcript objects)
    │
    ├─ transcripts_to_dataframe()       → transcripts.feather
    │
    ├─ compute_nrna_table()             → nrna.feather
    │      (group transcripts by unique genomic span)
    │
    ├─ build_splice_junctions()         → sj.feather
    │      (extract introns → sorted SJ table)
    │
    └─ build_genomic_intervals()        → intervals.feather
           (tile genome: EXON, TRANSCRIPT, INTERGENIC, UNAMBIG_INTRON)
```

**Output:** Directory containing 5+ Feather files.

---

## 4. Phase 1 — Configuration & Dispatch

**Entry:** `cli.py:main()` → `cli.py:quant_command()`

```
CLI args + optional YAML config
    │
    ├─ _build_quant_defaults()       ← defaults from frozen dataclasses
    ├─ _resolve_quant_args()         ← merge: CLI > YAML > defaults
    ├─ _build_pipeline_config()      ← construct PipelineConfig
    │       EMConfig + BamScanConfig + FragmentScoringConfig
    │
    ├─ TranscriptIndex.load()        ← load index from Feather files
    │       ├─ build collapsed interval index (vectorized numpy)
    │       ├─ build SJ exact-match map
    │       └─ build C++ FragmentResolver context:
    │              ctx.build_overlap_index()
    │              ctx.build_sj_map()
    │              ctx.build_sj_gap_index()
    │              ctx.set_metadata()
    │
    └─ run_pipeline(bam_path, index, config)
```

**Key:** The `FragmentResolver` C++ object is initialized here with all genomic interval data. It is passed by reference to all downstream C++ modules.

---

## 5. Phase 2 — BAM Scanning & Model Training

**Entry:** `pipeline.py:scan_and_buffer()`

This is a **single-pass** read of the BAM file. All work happens in C++.

```
BAM file
    │
    ├─ [Optional] detect_sj_strand_tag()     (C++, reads first 1000 spliced reads)
    │
    ├─ FragmentBuffer()                       (Python: manages chunks + spill)
    │       └─ FragmentAccumulator()          (C++: columnar accumulator)
    │
    └─ BamScanner(resolver_ctx, sj_tag_spec, skip_dups, include_multimap)
        │
        └─ scanner.scan(bam_path, n_workers)        ──── C++ ────
            │
            │   Producer thread:
            │     htslib sam_read1() → ParsedAlignment
            │     group by query-name → QnameGroup
            │
            │   Worker threads (N):
            │     QnameGroup → AssembledFragment
            │     AssembledFragment → FragmentResolver._resolve_core()
            │     ResolvedFragment → FragmentAccumulator.append()
            │     Training observations → strand/fraglen accumulators
            │
            └─ Returns: {stats, strand_observations, frag_length_observations}

    Post-scan (Python):
    │
    ├─ _replay_strand_observations()         → StrandModels (Beta posterior)
    │       exonic_spliced, exonic, intergenic sub-models
    │
    ├─ _replay_fraglen_observations()        → FragmentLengthModels (histogram)
    │       global + per-SpliceType + intergenic + gDNA sub-models
    │
    ├─ scanner.finalize_accumulator()        (C++ → raw numpy buffers)
    │
    ├─ Build _FinalizedChunk from raw buffers (numpy frombuffer + copy)
    │
    └─ Inject chunk into FragmentBuffer
```

**Output:**
- `FragmentBuffer` — one chunk of CSR-layout fragment data (18 numpy arrays per chunk)
- `StrandModels` — trained Beta strand model
- `FragmentLengthModels` — trained histogram fragment-length model
- `PipelineStats` — 40+ counters from BAM scanning

### CSR Fragment Layout (per chunk)

Each `_FinalizedChunk` stores N fragments in columnar + CSR format:

| Array | Dtype | Shape | Description |
|-------|-------|-------|-------------|
| `t_offsets` | int32 | [N+1] | CSR row pointers into t_indices |
| `t_indices` | int32 | [M] | Candidate transcript indices |
| `frag_lengths` | int32 | [M] | Per-candidate fragment length |
| `exon_bp` | int32 | [M] | Per-candidate exonic base pairs |
| `intron_bp` | int32 | [M] | Per-candidate intronic base pairs |
| `unambig_intron_bp` | int32 | [M] | Per-candidate unambiguous intronic bp |
| `splice_type` | uint8 | [N] | UNSPLICED / SPLICED_UNANNOT / SPLICED_ANNOT |
| `exon_strand` | uint8 | [N] | Strand of exonic overlap |
| `sj_strand` | uint8 | [N] | Strand from splice junction |
| `num_hits` | uint16 | [N] | NH tag (1 = unique, >1 = multimapper) |
| `merge_criteria` | uint8 | [N] | How transcript sets were merged |
| `chimera_type` | uint8 | [N] | NONE / TRANS / CIS_SAME / CIS_DIFF |
| `ambig_strand` | uint8 | [N] | Mixed-strand flag (0/1) |
| `frag_id` | int32 | [N] | BAM read-group serial number |
| `read_length` | int16 | [N] | Read length in bp |
| `genomic_footprint` | int32 | [N] | Genomic span of alignment |
| `genomic_start` | int32 | [N] | Leftmost genomic position |
| `nm` | int16 | [N] | NM tag (edit distance) |

---

## 6. Phase 3 — Fragment Scoring & EM Data Prep

**Entry:** `pipeline.py:quant_from_buffer()` → `scan.py:FragmentRouter.scan()`

This phase reads the buffered fragments, scores each fragment against every candidate transcript + nRNA + gDNA component, and produces the global CSR data structure for EM.

```
FragmentBuffer + trained models
    │
    ├─ Compute effective transcript lengths
    │       frag_length_model.compute_all_transcript_eff_lens()
    │       (Salmon-style eCDF convolution)
    │
    ├─ Build TranscriptGeometry
    │       (effective_lengths, exonic_lengths, t_to_g, mean_frag, spans)
    │
    ├─ Build FragmentScorer
    │       scoring.FragmentScorer.from_models()
    │       ├─ Extracts strand log-probs from StrandModels
    │       ├─ Extracts fragment-length LUT from FragmentLengthModels
    │       ├─ Builds per-transcript exon position cache
    │       └─ Constructs C++ NativeFragmentScorer ─── C++ ────
    │              (copies all scoring tables to C++ memory)
    │
    └─ FragmentRouter.scan(buffer)
        │
        └─ _scan_native()
            │
            │  For each chunk in buffer.iter_chunks():
            │    Convert chunk arrays to contiguous numpy
            │
            │    NativeFragmentScorer.fused_score_buffer()     ──── C++ ────
            │      │
            │      │  Pass 1 (count): iterate fragments, count candidates
            │      │  Pass 2 (fill): allocate exact arrays, score everything
            │      │
            │      │  Per fragment:
            │      │    ├─ Classify: unambiguous / ambig-same-strand /
            │      │    │            ambig-opp-strand / multimapper / chimeric
            │      │    │
            │      │    ├─ Unambiguous → deterministic assignment
            │      │    │   estimator.assign_unambig() [Python callback]
            │      │    │
            │      │    └─ Ambiguous → score all candidates:
            │      │        ├─ mRNA: strand_loglik + fraglen_loglik +
            │      │        │        overhang_penalty + mismatch_penalty
            │      │        ├─ nRNA: strand_loglik + fraglen_loglik (unspliced)
            │      │        └─ gDNA: LOG_HALF + fraglen_loglik + splice_penalty
            │      │
            │      └─ Returns: scored CSR arrays (22 output arrays)
            │
            └─ Accumulates → ScoredFragments (global CSR)
```

**Output:** `ScoredFragments` — global CSR of all ambiguous fragments with log-likelihoods for EM.

### ScoredFragments Layout

| Array | Description |
|-------|-------------|
| `offsets` | CSR row pointers (N_units + 1) |
| `t_indices` | Candidate component indices |
| `log_liks` | Per-candidate log-likelihoods |
| `count_cols` | Per-candidate SpliceStrandCol (6-column category) |
| `coverage_weights` | Per-candidate coverage weight (bias correction) |
| `tx_starts`, `tx_ends` | Transcript-space positions for bias model |
| `gdna_log_liks` | Pre-computed gDNA log-likelihoods per unit |
| `genomic_footprints` | Per-unit genomic span (for gDNA effective length) |
| `locus_t_indices` | Best transcript per unit (for locus assignment) |
| `locus_count_cols` | Category of best transcript per unit |
| `is_spliced` | Boolean per unit |
| `splice_type` | Per-unit splice classification |
| `frag_class` | Per-unit fragment class code |

---

## 7. Phase 4 — Empirical Bayes Prior Computation

**Entry:** `pipeline.py:quant_from_buffer()` (continued)

All prior computation is **pure Python** (vectorized numpy). This prepares hierarchical priors for the C++ EM solver.

```
Estimator accumulators + trained models
    │
    ├─ compute_nrna_init()                              [locus.py, Python]
    │       Per-nRNA initialization from intronic evidence
    │       nRNA_init = (sense_intronic - anti_intronic) / (2*SS - 1)
    │
    ├─ compute_eb_gdna_priors()                         [locus.py, Python]
    │   │
    │   │   3-tier hierarchical shrinkage:
    │   │
    │   ├─ Level 1: Global gDNA density
    │   │     compute_gdna_density_hybrid()
    │   │     (strand component + density component, inverse-variance weighted)
    │   │
    │   ├─ Level 2: Per-reference density shrunk toward global
    │   │     _compute_ref_gdna_densities()
    │   │     κ_ref via Method-of-Moments (estimate_kappa)
    │   │
    │   └─ Level 3: Per-locus density shrunk toward parent reference
    │         _compute_per_locus_gdna_densities()
    │         κ_locus via Method-of-Moments
    │         → gdna_init = shrunk_density × locus_exonic_bp
    │
    ├─ compute_global_gdna_density()                    [estimator.py, Python]
    │       Global-level gDNA reads/bp for nRNA fraction computation
    │
    └─ compute_nrna_frac_priors()                       [estimator.py, Python]
        │
        │   3-tier hierarchical shrinkage for Beta(α, β) nRNA fraction:
        │
        ├─ Level 1: Per-nRNA η (nascent fraction)
        │     _compute_hybrid_nrna_frac_vec()
        │     (density component: D_intron - gdna_density
        │      strand component: (sense - anti) / (2*SS - 1)
        │      combined by inverse-variance weighting)
        │
        ├─ Level 2: Per-locus-strand η (aggregated)
        │     _aggregate_nrna_frac_by_group()
        │     κ_locus via estimate_kappa()
        │
        └─ Level 3: Global η (empirical mean)
              κ_global via estimate_kappa()
              Hierarchical smoothing:
                locus-strand ← shrink toward global (κ_global)
                per-nRNA ← shrink toward locus-strand (κ_locus)
              → nrna_frac_alpha, nrna_frac_beta per nRNA
```

**Output:**
- `gdna_inits` — list of per-locus gDNA initialization values
- `nrna_init` — per-nRNA initialization from intronic evidence
- `nrna_frac_alpha`, `nrna_frac_beta` — per-nRNA Beta prior parameters

---

## 8. Phase 5 — Locus-Level EM Quantification

**Entry:** `pipeline.py:quant_from_buffer()` (continued)

```
ScoredFragments + priors
    │
    ├─ build_loci()                                     [locus.py → C++ _em_impl]
    │       connected_components() via union-find
    │       Groups transcripts linked by shared fragments into loci
    │       Returns: list[Locus] (each with transcript_indices, unit_indices)
    │
    ├─ build_locus_em_data()                            [locus.py, Python]
    │   │   Per-locus data extraction with local renumbering:
    │   │
    │   │   Component layout per locus:
    │   │     [0, n_t)              → mRNA (one per transcript)
    │   │     [n_t, n_t + n_nrna)  → nRNA (one per unique span)
    │   │     [n_t + n_nrna]       → gDNA (one shadow component)
    │   │
    │   ├─ Extract candidate CSR from global ScoredFragments
    │   ├─ Renumber global → local component indices
    │   ├─ Deduplicate: keep best log-lik per (unit, component)
    │   ├─ Add gDNA candidates for unspliced units
    │   └─ Build LocusEMInput (CSR + priors + init vectors)
    │
    └─ estimator.run_batch_locus_em()                   [estimator.py → C++ _em_impl]
        │
        └─ batch_locus_em()                             ──── C++ ────
            │
            │   For each locus (parallel over loci via thread pool):
            │
            │   1. build_equiv_classes()
            │      Group units with identical candidate sets
            │      Build dense n×k matrices per class
            │      Deterministic sorting for reproducibility
            │
            │   2. Initialize theta:
            │      mRNA: unambig_totals weighted by effective_length
            │      nRNA: nrna_init (from intronic evidence)
            │      gDNA: gdna_init (from EB prior)
            │      + Dirichlet prior (alpha) + OVR prior (gamma)
            │
            │   3. SQUAREM-accelerated EM loop:
            │      │
            │      │  E-step: parallel_estep()
            │      │    Per equivalence class:
            │      │      posterior = exp(log_lik + log(theta) - max)
            │      │      normalize rows → fractional assignments
            │      │      accumulate column sums (Kahan summation)
            │      │    Thread-local accumulation → deterministic reduction
            │      │
            │      │  M-step: hierarchical_map_em_step()
            │      │    mRNA: theta ∝ em_totals × effective_length × prior
            │      │    nRNA: enforces nascent fraction via Beta(α,β) prior
            │      │    gDNA: strand symmetry penalty (targeted excess)
            │      │
            │      │  SQUAREM acceleration:
            │      │    r = θ₁ - θ₀, v = θ₂ - 2θ₁ + θ₀
            │      │    α_sq = -‖r‖/‖v‖, θ_new = θ₀ - 2α·r + α²·v
            │      │
            │      └─ Converge when max|Δθ| < convergence_delta
            │
            │   4. Post-EM:
            │      Prune components < prune_threshold × max_theta
            │      Assign fractional counts to estimator accumulators
            │      Track high-confidence assignments (posterior ≥ threshold)
            │
            └─ Returns: (total_gdna_em, locus_mrna[], locus_nrna[], locus_gdna[])
```

**Output:** Populated estimator accumulators:
- `em_counts` — per-transcript EM-assigned counts (N_t × 6 columns)
- `em_high_conf_counts` — high-confidence assignments
- `nrna_em_counts` — per-nRNA EM assignments
- `locus_results` — per-locus summary (gene set, mRNA/nRNA/gDNA counts)

---

## 9. Phase 6 — Output Generation

**Entry:** `cli.py:_write_quant_outputs()`

```
AbundanceEstimator
    │
    ├─ estimator.get_counts_df()          → quant.feather (+ .tsv)
    │     Per-transcript: mRNA + nRNA split by splice×strand, TPM, eff_length
    │
    ├─ estimator.get_gene_counts_df()     → gene_quant.feather (+ .tsv)
    │     Per-gene: sum of transcript counts, gene-level TPM
    │
    ├─ estimator.get_loci_df()            → loci.feather (+ .tsv)
    │     Per-locus: gene set, n_transcripts, mRNA/nRNA/gDNA counts
    │
    ├─ estimator.get_detail_df()          → quant_detail.feather (+ .tsv)
    │     Per-transcript × category long format
    │
    └─ _write_run_config() + summary.json
          Comprehensive run metadata:
          version, command, timestamp, library protocol,
          alignment stats, strand model, frag-length model,
          quantification totals, all pipeline stats
```

---

## 10. Phase 7 — Annotated BAM (Optional)

**Entry:** `pipeline.py → annotate.py:write_annotated_bam()`

If `--annotated-bam` is specified, a second BAM pass stamps per-read tags:

```
Input BAM + AnnotationTable
    │
    └─ BamAnnotationWriter.write()     ──── C++ (_bam_impl) ────
           Re-reads BAM via htslib
           Matches reads to annotation table by frag_id
           Stamps tags: RG (pool), TP (best transcript),
                        GP (best gene), PP (posterior),
                        FC (fragment class), NC (n_candidates)
           Writes output BAM
```

---

## 11. Data Flow Diagram

```
                         ┌─────────────┐
                         │  FASTA + GTF │
                         └──────┬───────┘
                                │
                    ┌───────────▼────────────┐
                    │   Phase 0: Index Build  │  (offline, one-time)
                    │   Python: index.py      │
                    └───────────┬─────────────┘
                                │ Feather files
                                ▼
 ┌──────────┐      ┌────────────────────────┐
 │ BAM file │─────▶│  Phase 1: Config/Load   │
 └──────────┘      │  Python: cli.py         │
                   │  C++: FragmentResolver   │
                   └───────────┬──────────────┘
                               │ PipelineConfig + TranscriptIndex
                               ▼
                   ┌────────────────────────┐
                   │  Phase 2: BAM Scan      │  SINGLE PASS
                   │  C++: BamScanner        │
                   │  → FragmentResolver     │
                   │  → FragmentAccumulator  │
                   └───────────┬─────────────┘
                               │ FragmentBuffer + trained models + stats
                               ▼
                   ┌────────────────────────┐
                   │  Phase 3: Scoring       │
                   │  C++: fused_score_buf   │
                   │  Python: FragmentRouter │
                   └───────────┬─────────────┘
                               │ ScoredFragments (global CSR)
                               ▼
                   ┌────────────────────────┐
                   │  Phase 4: EB Priors     │
                   │  Python: priors.py      │
                   │  Python: locus.py       │
                   └───────────┬─────────────┘
                               │ gdna_inits + nrna_frac priors
                               ▼
                   ┌────────────────────────┐
                   │  Phase 5: Batch EM      │
                   │  C++: batch_locus_em    │
                   │  (SQUAREM + parallel)   │
                   └───────────┬─────────────┘
                               │ em_counts + nrna_em_counts
                               ▼
                   ┌────────────────────────┐
                   │  Phase 6: Output        │
                   │  Python: cli.py         │
                   │  Feather + TSV + JSON   │
                   └────────────────────────┘
```

---

## 12. Locus Model

Each locus is a connected component of transcripts linked by shared ambiguous fragments. The EM deconvolution model for a locus has three pools:

```
Locus components: [mRNA₁, mRNA₂, ..., mRNAₜ, nRNA₁, nRNA₂, ..., nRNAₙ, gDNA]
                   ├─── T transcripts ───┤ ├── N unique nRNA spans ──┤  ├─1─┤
```

- **mRNA (mature RNA):** One component per transcript. Weighted by effective length (Salmon-style eCDF). Strand-specific likelihood from trained Beta model. Fragment-length likelihood from trained histogram.
- **nRNA (nascent RNA):** One component per unique genomic span (ref, strand, start, end). Multiple transcripts may share the same nRNA span. Initialized from intronic evidence. Constrained by hierarchical Beta(α,β) nascent fraction prior.
- **gDNA (genomic DNA contamination):** One shadow component per locus. Initialized from hierarchical EB density estimation. Subject to strand symmetry penalty (sense ≈ antisense for gDNA). Splice penalty: unspliced = 1.0, unannotated = 0.01, annotated = 0.0.

The EM solver uses SQUAREM acceleration with convergence criterion max|Δθ| < 1e-6 and a budget of 1000 iterations (divided by 3 for SQUAREM's triple-step cost).
