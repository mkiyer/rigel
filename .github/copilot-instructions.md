# Copilot Instructions for hulkrna

## ⚠️ Agent Collaboration Policy

**Always pause and discuss design decisions before implementing.** This is a co-development project. When a task involves architectural choices, new module design, API surface changes, or non-trivial refactors, present options and trade-offs first. Implement only after alignment.

## Project Overview

hulkrna is a **Bayesian RNA-seq read counting** tool (Python 3.12+). It quantifies gene/transcript expression from paired-end BAM files using a two-pass Bayesian approach that learns strand specificity and gene abundance priors.

**Pipeline:** `hulkrna index` → `hulkrna count` (two-pass Bayesian) → `hulkrna gather` (multi-sample aggregation) → `hulkrna pileup` (coverage tracks)

## Architecture

### Core library (`src/hulkrna/`)

| Module | Responsibility |
|---|---|
| `core.py` | Shared types: `Strand` (IntEnum bitflags), `Interval`, `GenomicInterval`, `IntervalType`, `RefInterval` |
| `gtf.py` | Low-level GTF parser (`GTF` dataclass); converts 1-based inclusive → 0-based half-open |
| `transcript.py` | `Transcript` dataclass; GTF→Transcript parsing via `Transcript.read_gtf()` |
| `index.py` | `HulkIndex` class: `build()` creates index from FASTA+GTF → Feather files; `load()` reads index into cgranges trees + SJ map; query methods for fast overlap |
| `bam.py` | BAM read-pair iterator (`parse_bam_file`); CIGAR→exon blocks + splice junctions |
| `fragment.py` | Merges read pairs into `Fragment` (exons, introns as `GenomicInterval`, r2 strand flip, insert_size) |
| `query_bam.py` | Pass 1: scans BAM → Fragment → queries HulkIndex → writes intermediate hit Feather file |
| `strand_model.py` | `StrandModel` — Bayesian strand protocol inference (Beta posterior from spliced reads) |
| `bayes_count.py` | `BayesCounter` — probabilistic read→gene assignment via log-posterior + abundance priors |
| `count.py` | `CountCategory`, `CountStrand`, `CountType` IntEnum: 12-element (4 categories × 3 strands) flat index for count arrays |
| `pileup.py` | Genome-wide coverage tracks (Zarr + numpy); multiprocessing pipeline |

### CLI (`src/hulkrna/cli.py`)

Unified command-line interface registered as `hulkrna` console entry point via `pyproject.toml`. Uses argparse with `set_defaults(func=handler)` dispatch pattern.

| Subcommand | Status | Handler |
|---|---|---|
| `hulkrna index` | ✅ Implemented | `index_command()` → `index.build_index()` |
| `hulkrna count` | ✅ Pass 1 implemented | `count_command()` → `query_bam.query_bam()` (Pass 2 not yet implemented) |
| `hulkrna gather` | 🔲 Stub | Logic in `deprecated/hulkrna_gather.py` |
| `hulkrna pileup` | 🔲 Stub | Logic in `pileup.py` |

### Deprecated code (`deprecated/`)

Legacy standalone scripts moved here during refactoring. These contain working logic that needs to be migrated into library modules + CLI subcommands:

| Script | Purpose (reference for migration) |
|---|---|
| `hulkrna_index.py` | Build reference from FASTA + GENCODE GTF → Feather files |
| `hulkrna_count3.py` | Two-pass Bayesian counter (most recent/complete version) |
| `hulkrna_count.py` | Single-pass deterministic counter |
| `hulkrna_count2.py` | Incomplete intermediate refactor |
| `hulkrna_gather.py` | Aggregate per-sample counts (migrating from HDF5 to Zarr) |
| `hulkrna_strand.py` | Standalone strand protocol inference |
| `gencode.py` | Obsolete standalone parser; not used by any pipeline module |

All deprecated files will be removed once their logic is properly housed in library modules and wired through `cli.py`.

## Key Conventions

### Coordinates
All internal coordinates are **0-based, half-open** (BED convention). `GTF.from_str()` converts from GTF 1-based inclusive on parse.

### Strand representation
`Strand` is an `IntEnum` in `core.py` with bitwise semantics: `NONE=0, POS=1, NEG=2, AMBIGUOUS=3` (POS|NEG). Use `|=` to combine. Methods: `from_str()`, `to_str()`, `from_is_reverse()`, `opposite()`.

### Type system (`core.py`)
- **`Interval(start, end)`** — simple 0-based half-open interval. Used by `Transcript` for exon coordinates.
- **`GenomicInterval(ref, start, end, strand)`** — positioned interval on a chromosome. Used by `Fragment` for aligned exon blocks and splice junctions.
- **`IntervalType`** — `EXON=0, INTRON=1, INTERGENIC=2, SJ=3`. Classifies intervals for index lookup.
- **`RefInterval(ref, start, end, strand, interval_type, t_index, g_index)`** — annotated interval with index metadata. Used for both tiling intervals (cgranges) and splice junctions (exact-match lookup).

### R2 strand flip convention
In `Fragment.from_reads()`, read 2's genomic strand is always flipped via `Strand.opposite()` before merging. In a proper FR pair: r1=POS, r2=NEG → flip r2 → POS → concordant. Same-strand chimeric pairs → flip → AMBIGUOUS.

### Duplicate handling
`Fragment` does **not** track duplicate status. Duplicates are kept or discarded at BAM parse time via `parse_bam_file(skip_duplicates=...)`. Retained duplicates are treated identically to other reads.

### Intermediate hit file (`query_bam.py`)
Pass 1 writes one row per fragment × reference interval hit to a Feather file:

| Column | Type | Description |
|---|---|---|
| `frag_id` | int32 | Groups rows belonging to the same fragment |
| `ref` | string | Chromosome of the query interval |
| `start` | int32 | Genomic start (0-based) |
| `end` | int32 | Genomic end (half-open) |
| `strand` | int8 | Genomic strand (Strand enum) |
| `interval_type` | int8 | EXON/INTRON/INTERGENIC/SJ |
| `overlap` | int32 | Overlap bases for cgranges hits; 0 for SJ exact matches |
| `t_index` | int32 | Transcript index |
| `g_index` | int32 | Gene index |

Reconstructable: `combined_strand` (OR of strands), `is_spliced`, `is_chimeric`, `insert_size`.

### Gene table
The gene table is **not stored in the index**. It is derived on-the-fly from the transcript table via `HulkIndex._build_gene_table()` (pandas groupby aggregation).

### Count type encoding
`CountType` encodes a 4×3 matrix as a flat IntEnum: `index = category * 3 + strand_idx`. Categories: intron/unspliced/spliced_unannot/spliced_annot. Strands: ambiguous/sense/antisense. Count arrays are shape `(N, 12)`.

### Serialization
- **Feather** (via pyarrow/pandas) with LZ4 compression — primary intermediate format
- TSV mirrors alongside for human readability (controlled by `--no-tsv` flag)
- **Zarr** for multi-sample aggregated output and coverage tracks (replacing HDF5)

### Import style
All imports within the package use **relative imports** (e.g., `from .core import Strand`).

### CLI & entry point pattern
- `cli.py` uses argparse with subcommands; registered as `hulkrna` console entry point in `pyproject.toml`
- Each subcommand dispatches to a library function via `set_defaults(func=handler)`
- Logging uses `[START]`/`[DONE]` bracketed phases; `DEBUG`-level progress every N fragments

### Bayesian counting (`bayes_count.py`)
1. Query exon intervals via cgranges → determine `CountType`
2. **Unique gene**: assign 1.0 with sense/antisense from `StrandModelResult`
3. **Multi-gene**: log-posterior = `log(abundance_prior) + log(sj_likelihood)`; normalize via log-sum-exp; fractional assignment (threshold: `MIN_FRAC = 0.01`)

### Interval search (`index.py`)
`merge_sets_with_relaxation()` implements **progressive relaxation**: intersect all → intersect non-empty → union. Used by `HulkIndex.search_intervals()` and `HulkIndex.search_splice_junctions()`.

## Environment & Build

- **Python**: 3.12+ (conda-forge)
- **Build backend**: `hatchling` via `pyproject.toml` (`src/hulkrna` layout)
- **Environment**: managed via `mamba_env.yaml` (conda-forge + bioconda channels); use `conda`/`mamba` for deps (not pip) since `cgranges` and `pysam` are conda-only
- **Install**: `pip install -e . --no-deps` (editable, after conda deps are installed)
- **Testing**: `pytest` with `pytest-cov`; run via `pytest tests/ -v`; synthetic fixtures in `tests/conftest.py`
- **Performance path**: pure Python first, then Cython or nanobind for hot paths after functionality is validated

### Critical dependencies

`pysam` (BAM I/O, requires `XS`/`NH` tags from STAR), `cgranges` (interval overlap), `numpy` (count arrays, log-sum-exp), `pandas` (DataFrames, Feather I/O), `pyarrow` (Feather serialization), `scipy` (Beta distribution CIs, lazy-imported), `zarr`+`numcodecs` (coverage tracks + gather output).

### Test infrastructure

Tests live in `tests/` using pytest fixtures from `conftest.py`:
- `MINI_GTF` — synthetic GENCODE-style GTF (2 genes, 3 transcripts on chr1)
- `mini_gtf_file` — writes MINI_GTF to a temp file
- `mini_fasta_file` — creates a 2-chromosome FASTA with pysam `.fai` index (chr1=2000bp, chr2=500bp with no genes for intergenic testing)

## Things to Know

- The pipeline **requires** STAR-aligned BAMs with `XS` (splice strand) and `NH` (hit count) tags
- Only **paired-end** reads are supported (`read.is_paired` enforced in `bam.py`)
- Multimappers (`NH > 1`) are filtered out, not probabilistically assigned
- This is a **standalone CLI tool**, not a workflow/pipeline framework — no Snakemake, Nextflow, etc.
- Remaining subcommands (`count`, `gather`, `pileup`) need to be migrated from deprecated scripts into library modules and wired through `cli.py`
