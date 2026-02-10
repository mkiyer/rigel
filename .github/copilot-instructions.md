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
| `core.py` | Shared types: `Strand` (IntEnum bitflags), `Interval`, `GenomicInterval`, `IntervalType`, `RefInterval`; merge types: `MergeCriteria`, `MergeResult`, `EMPTY_MERGE`, `merge_sets_with_criteria()`; count enums: `CountCategory`, `CountStrand`, `CountType` |
| `gtf.py` | Low-level GTF parser (`GTF` dataclass); converts 1-based inclusive → 0-based half-open |
| `transcript.py` | `Transcript` dataclass; GTF→Transcript parsing via `Transcript.read_gtf()` |
| `index.py` | `HulkIndex` class: `build()` creates index from FASTA+GTF → Feather files; `load()` reads index into cgranges trees + SJ map; query methods for fast overlap |
| `bam.py` | BAM reader for **name-sorted/collated** BAM; `parse_read()` converts CIGAR→exon blocks + splice junctions; `parse_bam_file()` groups by read name via `itertools.groupby`, pairs mates via `_pair_reads()`, yields `list[(r1, r2)]` per read name |
| `fragment.py` | Merges read pairs into `Fragment` via `Fragment.from_reads()` (exons, introns as `GenomicInterval`, r2 strand flip) |
| `count.py` | Two-pass counting pipeline: `_resolve_fragment()` (shared resolution), `pass1_learn()` (model training, no file output), `pass2_count()` (count assignment via `ReadCounter`). No intermediate files — BAM is read twice. |
| `strand_model.py` | `StrandModel` — Bayesian strand model: 2×2 count table + Beta posterior over `p_r1_sense`; `strand_likelihood()` for counting; JSON output |
| `insert_model.py` | `InsertSizeModel` — insert size distribution histogram: observe/log_likelihood for Bayesian counting; `InsertSizeModels` — per-category container (global + intergenic + per-`CountCategory`); JSON output |

### CLI (`src/hulkrna/cli.py`)

Unified command-line interface registered as `hulkrna` console entry point via `pyproject.toml`. Uses argparse with `set_defaults(func=handler)` dispatch pattern.

| Subcommand | Status | Handler |
|---|---|---|
| `hulkrna index` | ✅ Implemented | `index_command()` → `HulkIndex.build()` |
| `hulkrna count` | ✅ Implemented (two-pass) | `count_command()` opens BAM twice, passes `bamfh.fetch(until_eof=True)` iterator to `pass1_learn(bam_iter, ...)` then `pass2_count(bam_iter, ...)`. Writes `strand_model.json`, `insert_size_models.json`, `transcript_counts.feather`, `gene_counts.feather`, `stats.json` to `--output-dir`. |
| `hulkrna gather` | 🔲 Stub | Logic in `deprecated/hulkrna_gather.py` |
| `hulkrna pileup` | 🔲 Stub | Logic in `deprecated/pileup.py` |

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
| `bayes_count.py` | Legacy probabilistic read→gene assignment via log-posterior + abundance priors |
| `hit_analysis.py` | Legacy Pass 2 code with intermediate Feather file streaming |
| `pileup.py` | Genome-wide coverage tracks (Zarr + numpy); for future `pileup` subcommand |

All deprecated files will be removed once their logic is properly housed in library modules and wired through `cli.py`.

## Key Conventions

### Coordinates
All internal coordinates are **0-based, half-open** (BED convention). `GTF.from_str()` converts from GTF 1-based inclusive on parse.

### Strand representation
`Strand` is an `IntEnum` in `core.py` with bitwise semantics: `NONE=0, POS=1, NEG=2, AMBIGUOUS=3` (POS|NEG). Use `|=` to combine. Methods: `from_str()`, `to_str()`, `from_is_reverse()`, `opposite()`.

### Type system (`core.py`)
- **`Interval(start, end)`** — simple 0-based half-open interval. Used by `Transcript` for exon coordinates.
- **`GenomicInterval(ref, start, end, strand)`** — positioned interval on a chromosome. Used by `Fragment` for aligned exon blocks and splice junctions.
- **`IntervalType`** — `EXON=0, INTRON=1, INTERGENIC=2, SJ=3, SJ_UNANNOT=4`. Classifies intervals for index lookup. `SJ` = annotated splice junction (exact match in index); `SJ_UNANNOT` = unannotated splice junction (CIGAR N-op with no index match, recorded with sentinel t_index/g_index = -1).
- **`RefInterval(ref, start, end, strand, interval_type, t_index, g_index)`** — annotated interval with index metadata. Used for both tiling intervals (cgranges) and splice junctions (exact-match lookup).

### R2 strand flip convention
In `Fragment.from_reads()`, read 2's genomic strand is always flipped via `Strand.opposite()` before merging. In a proper FR pair: r1=POS, r2=NEG → flip r2 → POS → concordant. Same-strand chimeric pairs → flip → AMBIGUOUS.

### BAM input requirements
The pipeline **requires name-sorted or collated BAM files** (e.g., `samtools sort -n` or `samtools collate`). `parse_bam_file()` uses `itertools.groupby` on `query_name` for O(1) memory grouping — coordinate-sorted BAMs will silently produce incorrect results.

### Mate pairing (`bam.py`)
`_pair_reads(r1_reads, r2_reads)` pairs mates using position-based matching: indexes r2 reads by `(reference_id, reference_start)`, matches each r1's `(next_reference_id, next_reference_start)`. Reciprocal verification ensures r2 points back to r1. Unpaired reads (no reciprocal match) are skipped. Returns `list[(r1, r2)]`.

`parse_bam_file()` groups by `query_name`, filters (QC, secondary, supplementary, unmapped, duplicate), checks NH tag, and yields `list[(r1, r2)]` per read name. The `include_multimap=False` parameter controls whether NH>1 groups are included.

### Duplicate handling
`Fragment` does **not** track duplicate status. Duplicates are kept or discarded at BAM parse time via `parse_bam_file(skip_duplicates=...)`. Retained duplicates are treated identically to other reads.

### Two-pass counting pipeline (`count.py`)

The `count` command reads the BAM file **twice** — no intermediate files:

1. **Pass 1** (`pass1_learn(bam_iter, index, ...)`): resolves fragments, trains `StrandModel` + `InsertSizeModels` (per-category) from unambiguous fragments. Returns `(stats, strand_model, insert_models)`. Takes a `pysam.AlignedSegment` iterator, not a file path.
2. **Pass 2** (`pass2_count(bam_iter, index, strand_model, insert_models, ...)`): re-resolves fragments identically, assigns fractional counts into `ReadCounter` arrays. Returns `(stats, counter)`. Also takes an iterator.

Both passes share `_resolve_fragment(frag, index)` which returns `_ResolveResult(merge_result, count_cat, exon_strand, sj_strand)` or `None` for intergenic fragments. Insert size computation is separate from resolution.

### ReadCounter (`count.py`)

- `t_counts`: `np.float32` array shape `(num_transcripts, 12)`
- `g_counts`: `np.float32` array shape `(num_genes, 12)`
- `assign(result, index, num_hits=1)`: determines `CountType` from `count_cat` + strand comparison, assigns fractional weights
- **Unique gene**: compares `exon_strand` to `g_to_strand_arr[g_idx]` → `SENSE` / `ANTISENSE`. `NONE` / `AMBIGUOUS` strand → `CountStrand.AMBIGUOUS`.
- **Multi-gene**: `CountStrand.AMBIGUOUS`, fractional `1/N_transcripts` and `1/N_genes`.
- **Multimappers**: each alignment pair weighted by `1/num_hits`.
- `get_t_counts_df()` / `get_g_counts_df()`: return pandas DataFrames with CountType column names.

**Intergenic fragments** produce no counts — only a stats counter. **Chimeric fragments** (multi-ref) are also skipped with a stats counter.

### Insert size computation (`count.py`)
Computed during Pass 1 per fragment after resolution. Two functions:

- **`_fragment_insert_size(frag)`** — simple footprint minus observed introns (no gap correction). Used for intergenic fragments.
- **`_compute_insert_size(frag, compatible_t_inds, index)`** — with gap correction:
  1. **Footprint** = `max(end) − min(start)` of fragment’s exon blocks.
  2. **Observed introns** (CIGAR N-ops) subtracted → **upper bound**.
  3. **Gaps** between consecutive exon blocks that are not observed introns → query `HulkIndex.query_gap_sjs()` for reference introns fully contained in each gap, filtered to compatible `t_index` set.
  4. Returns the **unambiguous insert size** if all compatible transcripts agree, or **-1** if they disagree.

### Gene table
The gene table is **not stored in the index**. It is derived on-the-fly from the transcript table via `HulkIndex._build_gene_table()` (pandas groupby aggregation).

### Count type encoding
`CountType` encodes a 4×3 matrix as a flat IntEnum: `index = category * 3 + strand_idx`. Categories: intron/unspliced/spliced_unannot/spliced_annot. Strands: ambiguous/sense/antisense. Count arrays are shape `(N, 12)`.

### Serialization
- **Feather** (via pyarrow/pandas) — count output tables (`transcript_counts.feather`, `gene_counts.feather`)
- TSV mirrors alongside for human readability (controlled by `--no-tsv` flag)
- **JSON** — model files (`strand_model.json`, `insert_size_models.json`) and stats file (`stats.json`, with `pass1` and `pass2` sub-dicts)
- **Zarr** for multi-sample aggregated output and coverage tracks (replacing HDF5)
- No intermediate hit files — the two-pass architecture processes BAM directly

### Import style
All imports within the package use **relative imports** (e.g., `from .core import Strand`).

### CLI & entry point pattern
- `cli.py` uses argparse with subcommands; registered as `hulkrna` console entry point in `pyproject.toml`
- Each subcommand dispatches to a library function via `set_defaults(func=handler)`
- Logging uses `[START]`/`[DONE]` bracketed phases; `DEBUG`-level progress every N fragments

### Bayesian counting (future)
`strand_model.py` provides `strand_likelihood()` and `insert_model.py` provides `log_likelihood()` for a future Bayesian Pass 2 upgrade. Currently unused by the active pipeline but retained for that purpose.

### Merge types (`core.py`)

- **`MergeCriteria`** IntEnum (INTERSECTION / INTERSECTION_NONEMPTY / UNION / EMPTY) tracks which relaxation level succeeded in progressive set merging.
- **`MergeResult`** (frozen dataclass) wraps merged `t_inds` / `g_inds` frozensets with the criteria used. Properties: `is_unique_gene`, `is_unique_transcript`, `is_empty`.
- **`EMPTY_MERGE`** — sentinel `MergeResult` with empty sets and EMPTY criteria.
- **`merge_sets_with_criteria(t_sets, g_sets)`** — progressive relaxation: intersect all → intersect non-empty → union. Returns `MergeResult`.

### Fragment resolution (`count.py`)

**`_resolve_fragment(frag, index)`** returns `_ResolveResult` NamedTuple or `None` (intergenic/empty).

**Resolution algorithm**:
1. Query each fragment exon block against cgranges → partition overlaps by interval type (EXON/INTRON). Build per-block transcript/gene sets. No overlap values are computed.
2. Query each fragment intron (CIGAR N-op) against the SJ exact-match map → track annotated vs unannotated status and SJ strand.
3. If **any EXON overlap** exists → merge EXON + SJ transcript sets together via `merge_sets_with_criteria()` for maximum specificity. Determine `CountCategory` (SPLICED_ANNOT / SPLICED_UNANNOT / UNSPLICED).
4. If **no EXON overlap** → try **INTRON fallback**: merge INTRON transcript sets. `CountCategory = INTRON`.
5. If neither → return `None` (intergenic).

**Chimeric check** happens in the caller loop — fragments with exon blocks on multiple references are skipped.

### Model training (`count.py` — Pass 1)

**`pass1_learn()`** trains `StrandModel` and `InsertSizeModels` inline during the Pass 1 loop. Returns `tuple[dict, StrandModel, InsertSizeModels]`.

**Strand model qualification**: `count_cat == SPLICED_ANNOT` + unique gene (from combined exon+SJ merge) + both `exon_strand` and `sj_strand` are POS or NEG (not NONE or AMBIGUOUS). Qualification logic lives in `pass1_learn()`, not in `StrandModel`.

**Insert size training**: Uses `_compute_insert_size()` (returns unambiguous scalar or -1). Unambiguous results are observed into both the per-category model and the global model via `InsertSizeModels.observe(size, count_cat)`. Intergenic fragments use `_fragment_insert_size()` with `count_cat=None`. **Per-category models**: one `InsertSizeModel` per `CountCategory` (INTRON, UNSPLICED, SPLICED_UNANNOT, SPLICED_ANNOT) + global + intergenic.

### Legacy Pass 2 code (`deprecated/hit_analysis.py`)

Contains `ResolvedFragment`, `resolve_fragment_hits()`, `analyze_hits()`, `iter_fragment_hits()`, `HitAnalysisStats`. Moved to `deprecated/` — retained for reference during future Bayesian Pass 2 redesign.

### Strand model (`strand_model.py`)

Single `StrandModel` dataclass — no enums, no separate result class.

- **2×2 count table**: `pos_pos`, `pos_neg`, `neg_pos`, `neg_neg` (exon alignment strand × SJ reference strand).
- **`observe(exon_strand, sj_strand)`** — increments one cell.
- **Derived**: `n_same` (pos_pos + neg_neg), `n_opposite` (pos_neg + neg_pos).
- **Beta posterior**: `alpha = prior_alpha + n_same`, `beta = prior_beta + n_opposite`.
- **`p_r1_sense`** — posterior mean P(read 1 aligns sense to gene). High ≈ FR, low ≈ RF, 0.5 ≈ unstranded.
- **`strand_likelihood(exon_strand, gene_strand)`** — returns `p_r1_sense` if same direction, `1 - p_r1_sense` if opposite, 0.5 if uninformative.
- **`read1_sense` / `read2_antisense`** — derived booleans for downstream protocol flags.
- **`to_dict()`** — full JSON/YAML-serializable summary (observations, posterior, probabilities, protocol).
- **`write_json(path)`** — writes strand model to JSON file (includes 95% CI when scipy available).

### Insert size model (`insert_model.py`)

**`InsertSizeModel`** dataclass — single histogram-based distribution.

- **Histogram**: `counts` array of shape `(max_size + 1,)`, float64. `counts[k]` = fragments with insert size *k*; `counts[max_size]` = overflow bin.
- **`max_size`** — default 1000; insert sizes ≥ this are clamped to overflow bin.
- **`observe(insert_size, weight=1.0)`** — increments histogram bin.
- **Summary stats**: `mean`, `std`, `median`, `mode` — all derived from histogram.
- **`log_likelihood(insert_size)`** — log-probability under the learned distribution with Laplace smoothing.
- **`to_dict()`** — JSON/YAML-serializable summary with trimmed histogram.
- **`write_json(path)`** — writes single insert size model to JSON file.

**`InsertSizeModels`** container — per-category + global + intergenic.

- Wraps `global_model`, `intergenic`, and `category_models` (dict keyed by `CountCategory`).
- **`observe(insert_size, count_cat=None, weight=1.0)`** — routes to global + appropriate sub-model. `count_cat=None` for intergenic.
- **`n_observations`** — delegates to `global_model`.
- **`to_dict()`** / **`write_json(path)`** / **`log_summary()`** — serialise all models.

## Environment & Build

- **Python**: 3.12+ (conda-forge)
- **Build backend**: `hatchling` via `pyproject.toml` (`src/hulkrna` layout)
- **Environment**: managed via `mamba_env.yaml` (conda-forge + bioconda channels); use `conda`/`mamba` for deps (not pip) since `cgranges` and `pysam` are conda-only
- **Install**: `pip install -e . --no-deps` (editable, after conda deps are installed)
- **Testing**: `pytest` with `pytest-cov`; run via `pytest tests/ -v`; synthetic fixtures in `tests/conftest.py`
- **Performance path**: pure Python first, then Cython or nanobind for hot paths after functionality is validated

### Critical dependencies

`pysam` (BAM I/O, requires `XS`/`NH` tags from STAR), `cgranges` (interval overlap), `numpy` (count arrays, log-sum-exp), `pandas` (DataFrames, Feather I/O), `pyarrow` (Feather serialization), `scipy` (Beta distribution CIs, lazy-imported). **Optional**: `zarr`+`numcodecs` (coverage tracks + gather output — listed under `[project.optional-dependencies] pileup` in `pyproject.toml`).

### Test infrastructure

Tests live in `tests/` using pytest fixtures from `conftest.py`:
- `MINI_GTF` — synthetic GENCODE-style GTF (2 genes, 3 transcripts on chr1)
- `mini_gtf_file` — writes MINI_GTF to a temp file
- `mini_fasta_file` — creates a 2-chromosome FASTA with pysam `.fai` index (chr1=2000bp, chr2=500bp with no genes for intergenic testing)

## Things to Know

- The pipeline **requires** STAR-aligned, **name-sorted or collated** BAMs with `XS` (splice strand) and `NH` (hit count) tags
- Only **paired-end** reads are supported (`read.is_paired` enforced in `bam.py`)
- Multimappers (`NH > 1`) are excluded by default; `--include-multimap` flag includes them in processing with `1/num_hits` fractional weighting in `ReadCounter.assign()`
- This is a **standalone CLI tool**, not a workflow/pipeline framework — no Snakemake, Nextflow, etc.
- Remaining subcommands (`gather`, `pileup`) need to be migrated from deprecated scripts into library modules and wired through `cli.py`
- The `count` command runs two passes: Pass 1 (model training) then Pass 2 (count assignment). BAM is read twice — no intermediate files.
- `hit_analysis.py`, `bayes_count.py`, and `pileup.py` have been moved to `deprecated/`. Not part of the active package.
- `bam.py` stats keys: `['total', 'qc_fail', 'unmapped', 'secondary', 'supplementary', 'duplicate', 'n_read_names', 'unique', 'multimapping', 'proper_pair', 'improper_pair', 'mate_unmapped', 'unpaired']`
- Pass 1 stats keys: `['n_fragments', 'n_chimeric', 'n_intergenic', 'n_with_exon', 'n_with_intron_fallback', 'n_with_annotated_sj', 'n_with_unannotated_sj', 'n_unique_gene', 'n_multi_gene', 'n_strand_skipped_no_sj', 'n_strand_skipped_multi_gene', 'n_strand_skipped_ambiguous', 'n_insert_unambiguous', 'n_insert_ambiguous', 'n_insert_intergenic']`
- Pass 2 stats keys: `['n_fragments', 'n_chimeric', 'n_intergenic', 'n_with_exon', 'n_with_intron_fallback', 'n_unique_gene', 'n_multi_gene']`
