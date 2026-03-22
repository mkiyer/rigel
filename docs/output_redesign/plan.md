# Rigel Output File Redesign Plan

**Date:** 2026-03-22
**Status:** APPROVED — implementation plan below

---

## 1. Current State Assessment

### 1.1 What currently exists (`rigel quant`)

| File | Format | Source | Notes |
|------|--------|--------|-------|
| `config.json` | JSON | `cli._write_run_config()` | CLI parameters only, no tool version |
| `quant.feather` (+`.tsv`) | Feather/TSV | `estimator.get_counts_df()` | Transcript-level counts |
| `gene_quant.feather` (+`.tsv`) | Feather/TSV | `estimator.get_gene_counts_df()` | Gene-level rollup |
| `loci.feather` (+`.tsv`) | Feather/TSV | `estimator.get_loci_df()` | Locus-level gDNA/mRNA breakdown |
| `quant_detail.feather` (+`.tsv`) | Feather/TSV | `estimator.get_detail_df()` | Long-format sparse QC breakdown |
| `summary.json` | JSON | `cli._write_quant_outputs()` | Library, alignment, quantification stats |
| Annotated BAM (optional) | BAM | `annotate.py` via `--annotated-bam` | Per-read assignment tags |

### 1.2 Issues Found

1. **`config.json` and `summary.json` are redundant and poorly separated.** `config.json` has parameters, `summary.json` has version + stats + full model dumps. The version should be in config, and config should be in summary.

2. **No calibration results in output.** The `GDNACalibration` object contains critical diagnostic info (gDNA density, mixing proportion, κ_strand, convergence metadata, per-region posteriors) that is only logged, never written.

3. **No dedicated nRNA output.** Synthetic nRNA transcript counts are embedded in `quant.feather` with `nrna > 0`, but there's no nRNA-level summary table. Users must filter `quant.feather` on `is_synthetic_nrna` (which isn't in the output) or `nrna > 0` to extract nRNA info.

4. **Transcript output lacks coordinates and nRNA linkage.** `quant.feather` has no `ref`, `start`, `end`, `strand`, `length`, or `nrna_id` column. Users need the index to map transcripts to genomic location or to their associated nRNA entity.

5. **Locus output lacks nRNA count.** The loci table has `mrna` and `gdna` but no `nrna` column. The three-pool model (mRNA + nRNA + gDNA) should be visible at locus level.

6. **Gene output lacks `gene_type` (biotype).** Available from the index but not included.

7. **Fragment length and strand distributions not easily extractable.** They're embedded in deeply nested `summary.json` objects (full histograms). Summary stats are present but the actual distributions require parsing nested JSON.

8. **`quant_detail.feather` is write-only.** It's not covered by golden tests, and it's unclear how useful the long-format sparse representation is to end users. It's more of an internal QC artifact.

9. **No `is_basic` / `is_mane` filter flags in output.** These are in the index but not propagated to the quant output. Users often want to filter to MANE Select or GENCODE Basic transcripts.

---

## 2. Proposed Output Files

### File 1: `summary.json` — Run Summary (REDESIGN)

Merge current `config.json` and `summary.json` into a single comprehensive run summary.

```
{
  "rigel_version": "0.x.y",
  "timestamp": "2026-03-22T12:34:56",

  "command": {
    "subcommand": "quant",
    "arguments": { ... all CLI args ... },
    "config_file": "path/to/config.yaml" or null
  },

  "configuration": {
    "em": { seed, prior_pseudocount, mode, iterations, ... },
    "scan": { skip_duplicates, include_multimap, max_frag_length, ... },
    "scoring": { overhang_log_penalty, mismatch_log_penalty, ... },
    "calibration": { max_iterations, convergence_tol, ... }
  },

  "input": {
    "bam_file": "...",
    "index_dir": "...",
    "index_nrna_tolerance": 20
  },

  "alignment_stats": {
    "total_reads": int,
    "mapped_reads": int,
    "unique_reads": int,
    "multimapping_reads": int,
    "proper_pairs": int,
    "duplicate_reads": int,
    "qc_fail_reads": int
  },

  "fragment_stats": {
    "total": int,
    "genic": int,
    "intergenic": int,
    "chimeric": int,
    "chimeric_trans": int,
    "chimeric_cis_same": int,
    "chimeric_cis_diff": int,
    "with_annotated_sj": int,
    "with_unannotated_sj": int
  },

  "strand_model": {
    "protocol": "R1-sense" | "R1-antisense",
    "strand_specificity": float,
    "p_r1_sense": float,
    "read1_sense": bool,
    "n_training_fragments": int,
    "posterior_variance": float,
    "ci_95": [lo, hi]
  },

  "calibration": {
    "gdna_density_global": float,
    "mixing_proportion": float,
    "expressed_density": float,
    "kappa_strand": float,
    "n_iterations": int,
    "converged": bool,
    "gdna_fl_mean": float | null,
    "gdna_fl_mode": int | null,
    "n_eligible_regions": int
  },

  "fragment_length": {
    "global": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "rna": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "gdna": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "intergenic": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "spliced_annot": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "spliced_unannot": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int },
    "unspliced": { "mean": float, "std": float, "median": float, "mode": int, "n_obs": int }
  },

  "quantification": {
    "n_transcripts": int,
    "n_genes": int,
    "n_loci": int,
    "n_nrna_entities": int,
    "n_unambig_assigned": int,
    "n_em_assigned": int,
    "mrna_total": float,
    "nrna_total": float,
    "gdna_total": float,
    "intergenic_total": int,
    "mrna_fraction": float,
    "nrna_fraction": float,
    "gdna_fraction": float
  }
}
```

**Changes from current:**
- Merges `config.json` into `summary.json` (one file, not two)
- Adds full `configuration` block (all frozen dataclass fields)
- Adds `calibration` section (currently missing entirely)
- Flattens `strand_model` to essential fields (removes diagnostic sub-models)
- Flattens `fragment_length` to summary stats per category (removes histograms)
- Removes `pipeline_stats` raw dump (replaced by structured `fragment_stats`)
- Adds `n_nrna_entities` to quantification

**Removed from `summary.json`:**
- Full histogram arrays (move to dedicated file if needed — see File 7)
- Raw diagnostic strand models (exonic, intergenic)
- Raw `pipeline_stats` dump with every internal counter

**Rationale:** The summary should be human-scannable and machine-parseable at a glance. Histograms, raw counts, and diagnostic models belong in separate diagnostic outputs. We want this file to answer: "What happened in this run?" in ~100 lines of JSON.

---

### File 2: `quant.tsv` — Transcript-Level Abundances (REDESIGN)

Primary quantification output. One row per transcript.

| Column | Type | Description |
|--------|------|-------------|
| `transcript_id` | str | Transcript ID (e.g., `ENST00000...`) |
| `gene_id` | str | Gene ID |
| `gene_name` | str | Gene name |
| `gene_type` | str | Gene biotype (e.g., `protein_coding`) |
| `ref` | str | Chromosome/reference |
| `strand` | str | `+` or `-` |
| `start` | int | Genomic start (0-based) |
| `end` | int | Genomic end (exclusive) |
| `length` | int | Spliced/exonic length |
| `effective_length` | float | Bias-corrected effective length |
| `locus_id` | int | EM locus (-1 if no locus) |
| `nrna_id` | str | Nascent RNA entity ID (synthetic transcript ID or `"."`) |
| `is_basic` | bool | GENCODE basic flag |
| `is_mane` | bool | MANE Select flag |
| `mrna` | float | Mature RNA fragment count |
| `mrna_unambig` | float | Uniquely-assigned mRNA |
| `mrna_em` | float | EM-assigned mRNA |
| `mrna_spliced` | float | Spliced mRNA fragments |
| `nrna` | float | Nascent RNA count (non-zero only for synthetic nRNAs and nascent-equiv) |
| `rna_total` | float | mrna + nrna |
| `tpm` | float | TPM (mRNA-based) |
| `posterior_mean` | float | Mean posterior probability |

**Changes from current:**
- Add `gene_type`, `ref`, `strand`, `start`, `end`, `length` from index
- Add `nrna_id` — the synthetic nRNA transcript this transcript maps to
- Add `is_basic`, `is_mane` flags for filtering
- Convert `strand` from int enum to `+`/`-` string for readability
- Drop `mrna_high_conf` — this is a QC/diagnostic metric, not a primary abundance estimate

**Removed:**
- `mrna_high_conf` — move to detail/diagnostic output if needed
- Integer strand encoding

**Format change:** Make TSV the primary output (not feather). Feather as optional `--feather` flag. TSV is the universal interchange format for bioinformatics tools (compatible with awk, R, pandas, spreadsheets). Feather is useful for large datasets but not standard.

---

### File 3: `gene_quant.tsv` — Gene-Level Abundances (REDESIGN)

One row per gene. Summary rollup from transcript level.

| Column | Type | Description |
|--------|------|-------------|
| `gene_id` | str | Gene ID |
| `gene_name` | str | Gene name |
| `gene_type` | str | Gene biotype |
| `ref` | str | Chromosome |
| `strand` | str | `+` or `-` |
| `start` | int | Gene start (min of transcript starts) |
| `end` | int | Gene end (max of transcript ends) |
| `n_transcripts` | int | Number of annotated transcripts |
| `locus_id` | int | Primary EM locus |
| `effective_length` | float | Abundance-weighted mean effective length |
| `mrna` | float | Gene-level mRNA count |
| `mrna_unambig` | float | Uniquely-assigned mRNA |
| `mrna_em` | float | EM-assigned mRNA |
| `mrna_spliced` | float | Spliced mRNA |
| `nrna` | float | Gene-level nascent RNA count |
| `rna_total` | float | mrna + nrna |
| `tpm` | float | TPM (mRNA-based) |

**Changes from current:**
- Add `gene_type`, `ref`, `strand`, `start`, `end`, `n_transcripts`
- Convert `strand` to `+`/`-`
- Drop `mrna_high_conf`
- TSV primary, feather optional

---

### File 4: `nrna_quant.tsv` — Nascent RNA Abundances (NEW)

One row per nRNA entity (synthetic nRNA transcript or annotated nascent equivalent).

| Column | Type | Description |
|--------|------|-------------|
| `nrna_id` | str | Nascent RNA entity ID (e.g., `RIGEL_NRNA_chr1_1_1000_5000`) |
| `gene_id` | str | Representative gene ID |
| `gene_name` | str | Representative gene name |
| `ref` | str | Chromosome |
| `strand` | str | `+` or `-` |
| `start` | int | nRNA span start |
| `end` | int | nRNA span end |
| `length` | int | nRNA span length |
| `effective_length` | float | Effective length |
| `locus_id` | int | EM locus |
| `is_synthetic` | bool | True for RIGEL-generated synthetics, False for annotated equivalents |
| `n_contributing_transcripts` | int | Number of annotated multi-exon transcripts merged into this nRNA |
| `count` | float | Nascent RNA fragment count |
| `tpm` | float | TPM (nRNA-based, normalized within nRNA pool) |

**Rationale:** This provides a dedicated view of nascent RNA biology. Users analyzing pre-mRNA, intron retention, or chromatin-associated RNA need this as a first-class output rather than filtering transcript tables. The `n_contributing_transcripts` gives users a sense of how many isoforms share each nRNA entity — useful for understanding locus complexity.

---

### File 5: `loci.tsv` — Locus-Level Results (REDESIGN)

One row per EM locus. Three-pool breakdown.

| Column | Type | Description |
|--------|------|-------------|
| `locus_id` | int | Locus identifier |
| `ref` | str | Primary chromosome |
| `n_transcripts` | int | Transcripts in locus (annotated + synthetic) |
| `n_annotated_transcripts` | int | Annotated transcripts only |
| `n_nrna_entities` | int | Synthetic nRNA transcripts in locus |
| `n_genes` | int | Genes in locus |
| `n_em_fragments` | int | Ambiguous fragments entering EM |
| `mrna` | float | mRNA count from locus EM |
| `nrna` | float | nRNA count from locus EM |
| `gdna` | float | gDNA count from locus EM |
| `total` | float | mrna + nrna + gdna |
| `gdna_rate` | float | gdna / total |
| `gdna_prior` | float | Calibration-derived gDNA prior (γ) |

**Changes from current:**
- Add `nrna` column (was missing — the three-pool model should be visible)
- Add `n_annotated_transcripts`, `n_nrna_entities` breakdown
- Add `total` column
- Rename `gdna_init` → `gdna_prior` (clearer)
- TSV primary, feather optional

---

### File 6: Annotated BAM (KEEP, optional)

No changes proposed. The current BAM tag schema (`ZT`, `ZG`, `ZP`, `ZW`, `ZC`, `ZH`, `ZN`, `ZS`) is well-designed for genome browser visualization and read-level QC. Activated via `--annotated-bam <path>`.

One minor addition worth considering: a `ZL` tag for locus_id, so reads can be grouped by EM locus in downstream analysis.

---

### File 7: `quant_detail.tsv` — Category Breakdown (RECONSIDER)

**Recommendation: Remove as a default output.** The sparse long-format (transcript × category × source) table is useful for deep QC but not for routine analysis. It creates large files for minimal value.

**Options:**
- (A) **Remove entirely.** The information can be reconstructed from the annotated BAM.
- (B) **Gate behind `--detail` flag.** Only write when explicitly requested.
- (C) **Keep but compress.** Only include rows with count > 0 (already sparse).

**Recommendation: Option B** — keep the capability behind a flag.

---

## 3. Additional Recommendations

### 3.1 Remove `config.json` — merge into `summary.json`

Two separate JSON files for a single run is confusing. `summary.json` should be the single authoritative record of the run — what was configured, what was observed, and what was produced.

### 3.2 Diagnostic histograms: `summary.json` or separate file?

The current `summary.json` dumps full fragment-length histograms (up to 1000 bins per category × 7 categories) and strand model raw counts. This bloats the summary.

**Recommendation:** Keep summary-level stats in `summary.json` (mean, std, mode, n_obs per category). If users need full histograms for plotting or downstream analysis, provide a `--diagnostics` flag that writes a `diagnostics/` subdirectory with:
- `frag_length_histograms.tsv` — size × category matrix
- `strand_model_detail.json` — all three strand models with CI
- `calibration_regions.tsv` — per-region posteriors, densities, strand fracs (very powerful diagnostic)

### 3.3 TSV as primary, Feather as optional

Most bioinformatics workflows consume TSV. Feather is fast for Python/R but not a universal format. Propose:
- TSV always written (default)
- `--feather` flag adds `.feather` mirrors
- `--no-tsv` flag suppresses TSV (for pipeline-only workflows)

This inverts the current default (feather primary, TSV optional).

### 3.4 Column naming conventions

Adopt consistent conventions:
- Snake_case for all column names
- No abbreviations except standard ones (`tpm`, `id`)
- `_id` suffix for identifiers
- `_count` or bare noun for fragment counts (currently inconsistent: `mrna` vs `n_em_fragments`)

### 3.5 Row ordering

- `quant.tsv`: Sort by `(ref, start, end)` — genomic order
- `gene_quant.tsv`: Sort by `(ref, start, end)` — genomic order
- `nrna_quant.tsv`: Sort by `(ref, start, end)` — genomic order
- `loci.tsv`: Sort by `locus_id`

### 3.6 Transcript nRNA linkage

Every annotated transcript should report which nRNA entity it maps to (via `nrna_id` in `quant.tsv`). This creates a join key between the transcript and nRNA tables. Computation: for each annotated multi-exon transcript, look up its merged (ref, strand, start, end) span and find the matching synthetic or nascent-equiv transcript ID.

Single-exon transcripts that are NOT nascent equivalents would have `nrna_id = "."` (no nRNA entity).

---

## 4. File Inventory Summary

| # | File | Status | Format |
|---|------|--------|--------|
| 1 | `summary.json` | REDESIGN (merge config.json) | JSON |
| 2 | `quant.tsv` | REDESIGN (add coords, nrna_id, flags) | TSV (+feather opt) |
| 3 | `gene_quant.tsv` | REDESIGN (add coords, gene_type) | TSV (+feather opt) |
| 4 | `nrna_quant.tsv` | **NEW** | TSV (+feather opt) |
| 5 | `loci.tsv` | REDESIGN (add nrna, rename) | TSV (+feather opt) |
| 6 | Annotated BAM | KEEP (optional, `--annotated-bam`) | BAM |
| 7 | `quant_detail.tsv` | GATE behind `--detail` flag | TSV (+feather opt) |
| 8 | `config.json` | **REMOVE** (merged into summary.json) | — |
| 9 | `diagnostics/` | **NEW** (optional, behind `--diagnostics`) | TSV/JSON |

---

## 5. Implementation Phases

### Phase 1: `summary.json` redesign
- Merge `config.json` into `summary.json`
- Add calibration section
- Flatten strand/fragment models to summary stats
- Remove raw histogram arrays from default output
- Remove `config.json` write

### Phase 2: Transcript output redesign
- Add coordinates, strand, length, gene_type, nrna_id, is_basic, is_mane
- Build nrna_id linkage (transcript → nRNA entity mapping)
- Invert TSV/feather default
- Drop mrna_high_conf

### Phase 3: Gene output redesign
- Add coordinates, gene_type, n_transcripts
- Invert TSV/feather default
- Drop mrna_high_conf

### Phase 4: nRNA output (new)
- Build `get_nrna_counts_df()` in estimator
- Include n_contributing_transcripts from index
- Compute nRNA-pool TPM

### Phase 5: Locus output redesign
- Add nrna column to locus results
- Add annotated/synthetic transcript breakdown
- Rename gdna_init → gdna_prior

### Phase 6: Detail and diagnostics
- Gate quant_detail behind --detail flag
- Implement --diagnostics output directory
- Calibration regions table
- Full histogram export

### Phase 7: Golden test updates
- Regenerate all golden outputs
- Add nrna_quant golden tests
- Add loci.nrna golden test column

---

## 6. Design Decisions

1. **`quant.tsv` includes only annotated transcripts.** Synthetic nRNAs go exclusively to `quant_nrna.tsv`. Annotated transcripts that are nascent equivalents (`is_nascent_equiv=True`) remain in `quant.tsv` with a flag column.
2. **nRNA ID format:** `RIGEL_NRNA_{ref}_{strand}_{start}_{end}` (already implemented).
3. **`is_nascent_equiv` transcripts** stay in primary `quant.tsv` (they are annotated transcripts). Their nascent role is reported via the `nrna_id` join column and the `is_nascent_equiv` flag.
4. **`quant_detail` output is removed entirely.**
5. **Column naming:** snake_case throughout; bare nouns for fragment counts (`mrna`, `nrna`, `gdna`); `_id` suffix for identifiers.
6. **Add `ZL` (locus_id) tag** to annotated BAM.

---

## 7. Implementation Plan

### Overview

Eight steps, ordered by dependency. Each step is self-contained and testable.

### Step 1: Index — store nRNA linkage and contributor counts

**Goal:** At index build time, record (a) which nRNA entity each annotated transcript maps to, and (b) how many annotated transcripts contribute to each nRNA entity. This data flows to the output tables in later steps.

**Files:**
- `src/rigel/transcript.py` — add two new fields to `Transcript` dataclass and `to_dict()`
- `src/rigel/index.py` — populate them in `create_nrna_transcripts()` and `build()`

**Changes:**

`transcript.py` — add to `Transcript` dataclass:
```python
nrna_t_index: int = -1          # t_index of the associated synthetic nRNA (-1 = none)
nrna_n_contributors: int = 0    # for synthetics: number of contributing annotated transcripts
```

Add both fields to `to_dict()` so they serialize to `transcripts.feather`.

`index.py` `create_nrna_transcripts()`:
- After Step 3 (dedup merged spans), build `merged_key_to_syn_index: dict[tuple, int]` mapping each `(ref, strand, start, end)` key to the synthetic's position in the `synthetics` list.
- After Step 5 (create synthetics), iterate multi-exon transcripts and set `t.nrna_t_index` to the synthetic's eventual `t_index` (computed as `len(transcripts) + syn_position`). For annotated nascent-equiv transcripts (step 4 `covered` set), set `t.nrna_t_index` to the equiv transcript's `t_index`.
- On each synthetic, set `nrna_n_contributors` to the count of annotated multi-exon transcripts that mapped to its span key.

**Return value change:** `create_nrna_transcripts()` additionally returns the `span_key → contributor_count` dict and `span_key → syn_list_index` mapping so the caller can assign `nrna_t_index` after `t_index` values are finalized.

Alternatively (simpler): have the function accept the `next_t_index` offset and assign `nrna_t_index` directly, since the synthetic indices are deterministic (offset + position in list).

**Tests:**
- Update `test_index_integrity.py` to verify `nrna_t_index` is populated on multi-exon transcripts
- Verify `nrna_n_contributors > 0` on synthetics
- Verify `nrna_t_index == -1` for single-exon transcripts that are NOT nascent equivalents
- Regenerate golden outputs (new columns in `transcripts.feather`)

---

### Step 2: `summary.json` redesign — merge config, add calibration

**Goal:** Single authoritative JSON file with version, configuration, inputs, and all summary statistics including calibration.

**Files:**
- `src/rigel/cli.py` — rewrite `_write_quant_outputs()` summary section; remove `_write_run_config()` call and `config.json` write; remove `config.json` path computation
- `src/rigel/calibration.py` — add `to_summary_dict()` method on `GDNACalibration`
- `src/rigel/config.py` — add `to_dict()` methods on all config frozen dataclasses

**`GDNACalibration.to_summary_dict()`:**
```python
def to_summary_dict(self) -> dict:
    d = {
        "gdna_density_global": self.gdna_density_global,
        "mixing_proportion": self.mixing_proportion,
        "expressed_density": self.expressed_density,
        "kappa_strand": self.kappa_strand,
        "n_iterations": self.n_iterations,
        "converged": self.converged,
    }
    if self.gdna_fl_model is not None:
        d["gdna_fl_mean"] = round(self.gdna_fl_model.mean, 2)
        d["gdna_fl_mode"] = int(self.gdna_fl_model.mode)
    if self.eligible is not None:
        d["n_eligible_regions"] = int(self.eligible.sum())
    return d
```

**Config `to_dict()` — use `dataclasses.asdict()`:**
```python
# On PipelineConfig:
def to_dict(self) -> dict:
    from dataclasses import asdict
    d = asdict(self)
    # Convert Path objects to strings
    if d.get("annotated_bam_path") is not None:
        d["annotated_bam_path"] = str(d["annotated_bam_path"])
    if d.get("scan", {}).get("spill_dir") is not None:
        d["scan"]["spill_dir"] = str(d["scan"]["spill_dir"])
    return d
```

**New `summary.json` structure** (as specified in Section 2 File 1 above), with these sections:
- `rigel_version`, `timestamp`
- `command` (subcommand + all CLI arguments + config file path)
- `configuration` (full `PipelineConfig.to_dict()`)
- `input` (bam_file, index_dir, index nrna_tolerance)
- `alignment_stats`
- `fragment_stats`
- `strand_model` (flattened summary from exonic_spliced model only)
- `calibration` (from `GDNACalibration.to_summary_dict()`)
- `fragment_length` (summary stats per category, no histograms)
- `quantification` (pool totals, fractions, counts)

**Removed:**
- `config.json` file entirely
- `_write_run_config()` function
- `strand_models` raw dump (replaced by flat summary)
- `frag_length_models` raw dump (replaced by flat summary)
- `pipeline_stats` raw dump (replaced by structured sections)

**Tests:**
- Update `test_cli.py` to not expect `config.json`
- Verify `summary.json` contains calibration section
- Verify `summary.json` contains full configuration

---

### Step 3: `quant.tsv` redesign — annotated transcripts only

**Goal:** Transcript-level abundances with coordinates, flags, and nRNA linkage. Only annotated transcripts (no synthetics).

**Files:**
- `src/rigel/estimator.py` — rewrite `get_counts_df()`
- `src/rigel/cli.py` — update write logic (TSV primary, feather optional)

**New `get_counts_df()` columns:**

```
transcript_id, gene_id, gene_name, gene_type,
ref, strand, start, end, length, effective_length,
locus_id, nrna_id,
is_basic, is_mane, is_nascent_equiv,
mrna, mrna_unambig, mrna_em, mrna_spliced,
nrna, rna_total,
tpm, posterior_mean
```

**Key changes:**
- **Filter out synthetics:** `mask = ~is_synthetic` applied to all arrays before DataFrame construction
- **Add coordinates:** from `index.t_df[["ref", "start", "end", "strand", "length"]]`
- **Convert strand:** `1 → "+"`, `2 → "-"`
- **Add `gene_type`:** from `index.t_df["g_type"]`
- **Add `nrna_id`:** lookup `index.t_df["nrna_t_index"]` → `index.t_df.loc[nrna_t_index, "t_id"]` (or `"."` if -1)
- **Add flags:** `is_basic`, `is_mane`, `is_nascent_equiv` from `index.t_df`
- **Drop `mrna_high_conf`**
- **`nrna` column for annotated transcripts:** 0.0 for all rows (their nRNA role is in `nrna_quant.tsv`; the `is_nascent_equiv` flag tells the user about the dual role)

**Sort order:** `(ref, start, end)` — genomic order

**Tests:**
- Update golden outputs — new columns, fewer rows (no synthetics)
- Verify no `RIGEL_NRNA_*` transcript IDs in output
- Verify `nrna_id` populated for multi-exon transcripts

---

### Step 4: `gene_quant.tsv` redesign

**Goal:** Gene-level rollup with coordinates, biotype, transcript count.

**Files:**
- `src/rigel/estimator.py` — rewrite `get_gene_counts_df()`

**New columns:**

```
gene_id, gene_name, gene_type,
ref, strand, start, end,
n_transcripts, locus_id, effective_length,
mrna, mrna_unambig, mrna_em, mrna_spliced,
nrna, rna_total,
tpm
```

**Key changes:**
- **Add coordinates:** from `index.g_df[["ref", "start", "end", "strand"]]`
- **Convert strand:** `1 → "+"`, `2 → "-"`
- **Add `gene_type`:** from `index.g_df["g_type"]`
- **Add `n_transcripts`:** from `index.g_df["num_transcripts"]` (note: this includes synthetics — should we compute annotated count only? Decision: use annotated-only count, computed at index time)
- **Drop `mrna_high_conf`**

**Sort order:** `(ref, start, end)` — genomic order

**Tests:** Update golden outputs

---

### Step 5: `quant_nrna.tsv` — new nascent RNA output

**Goal:** Dedicated nRNA entity table with counts, effective length, TPM.

**Files:**
- `src/rigel/estimator.py` — add `get_nrna_counts_df()`
- `src/rigel/cli.py` — write the new file

**`get_nrna_counts_df()` columns:**

```
nrna_id, gene_id, gene_name,
ref, strand, start, end, length, effective_length,
locus_id,
is_synthetic, n_contributing_transcripts,
count, tpm
```

**Implementation:**
- Filter `index.t_df` to rows where `is_synthetic_nrna == True` or `is_nascent_equiv == True`
- For synthetics: `nrna_id = t_id`, `count = t_total[syn_idx]`, `is_synthetic = True`
- For nascent equivalents: `nrna_id = t_id`, `count` = need to determine. These transcripts' counts go to `mrna` in the primary output. The nRNA role is informational — they mark an annotated transcript that *could* capture nascent RNA. Set `count = 0.0` with a note that their fragment counts are in the transcript table. Alternatively, they could be excluded from this table entirely. **Decision: include nascent equivalents with `count = 0.0` and `is_synthetic = False`** — this provides the complete nRNA entity catalog even though the EM doesn't route separately to them.
- `n_contributing_transcripts`: from `index.t_df["nrna_n_contributors"]`
- `effective_length`: from `estimator._t_eff_len[idx]`
- `tpm`: normalized within the nRNA pool only

**Sort order:** `(ref, start, end)` — genomic order

**Tests:** Add `nrna_quant` golden tests

---

### Step 6: `loci.tsv` redesign

**Goal:** Three-pool locus summary including nRNA.

**Files:**
- `src/rigel/estimator.py` — rewrite `get_loci_df()` to add `nrna` and breakdown columns
- `src/rigel/pipeline.py` — add `nrna` to `_build_locus_meta()` (or compute post-hoc)

**New columns:**

```
locus_id, ref,
n_transcripts, n_annotated_transcripts, n_nrna_entities, n_genes,
n_em_fragments,
mrna, nrna, gdna, total, gdna_rate, gdna_prior
```

**Implementation for `nrna` per locus:**

Post-hoc computation in `get_loci_df()`:
```python
# For each locus, sum t_counts for synthetic transcripts in that locus
syn_mask = self._synthetic_mask
t_total = self.t_counts.sum(axis=1)
locus_ids = self.locus_id_per_transcript

for locus in self.locus_results:
    lid = locus["locus_id"]
    mask = (locus_ids == lid) & syn_mask
    locus_nrna = float(t_total[mask].sum())
```

For `n_annotated_transcripts` and `n_nrna_entities`: count transcript_indices in locus where `is_synthetic_nrna` is False/True. This requires passing `is_synthetic_nrna` mask to `get_loci_df()` (via `index` parameter or stored on estimator).

**Also rename `gdna_init` → `gdna_prior`.**

**Tests:** Update golden loci outputs (new columns, renamed column)

---

### Step 7: Annotated BAM — add `ZL` (locus_id) tag

**Goal:** Add locus_id as an integer BAM tag on every annotated read.

**Files:**
- `src/rigel/annotate.py` — add `locus_id` array to `AnnotationTable`; pass to C++ writer
- `src/rigel/native/bam_scanner.cpp` — add `ZL` tag writing in `BamAnnotationWriter`
- Annotation table population code in `pipeline.py` / `scan.py`

**Changes:**

`AnnotationTable`:
- Add `locus_id: np.ndarray` (int32) to the dataclass
- Initialize to `-1` in `create()`
- Accept `locus_id` in `add()`
- Include in `_grow()` and `get()`

Population:
- When annotations are written (in `pipeline.py` scoring/routing), set `locus_id = locus_id_per_transcript[best_tid]` (with -1 fallback for intergenic)

C++ `BamAnnotationWriter.write()`:
- Accept additional `locus_id` array parameter
- Add `bam_aux_update_int(rec, "ZL", locus_id)` for each record

**Tests:**
- Update `test_annotate.py` to verify ZL tag presence and correct values
- Recompile C++ (`pip install --no-build-isolation -e .`)

---

### Step 8: CLI and format changes

**Goal:** TSV primary, feather optional; remove config.json; remove quant_detail; update file names.

**Files:**
- `src/rigel/cli.py` — rework `_write_quant_outputs()` and quant_command arg definitions

**CLI argument changes:**
- Remove `--no-tsv` flag from quant subparser
- Add `--feather` flag: `action="store_true", default=False` — writes `.feather` mirrors
- Remove `quant_detail.feather` / `.tsv` writes
- Add `quant_nrna.tsv` write
- Reorder writes: `summary.json`, `quant.tsv`, `gene_quant.tsv`, `quant_nrna.tsv`, `loci.tsv`

**File name changes:**
- `quant.feather` → `quant.tsv` (primary), `quant.feather` (optional)
- `gene_quant.feather` → `gene_quant.tsv` (primary), `gene_quant.feather` (optional)
- `loci.feather` → `loci.tsv` (primary), `loci.feather` (optional)
- NEW: `quant_nrna.tsv` (primary), `quant_nrna.feather` (optional)
- REMOVED: `config.json`, `quant_detail.feather`, `quant_detail.tsv`

**Output write logic:**
```python
quant_df.to_csv(output_dir / "quant.tsv", sep="\t", index=False)
gene_df.to_csv(output_dir / "gene_quant.tsv", sep="\t", index=False)
nrna_df.to_csv(output_dir / "quant_nrna.tsv", sep="\t", index=False)
loci_df.to_csv(output_dir / "loci.tsv", sep="\t", index=False)
if args.feather:
    quant_df.to_feather(str(output_dir / "quant.feather"))
    gene_df.to_feather(str(output_dir / "gene_quant.feather"))
    nrna_df.to_feather(str(output_dir / "quant_nrna.feather"))
    loci_df.to_feather(str(output_dir / "loci.feather"))
```

**Tests:**
- Update `test_cli.py` to check for new file names
- Remove `quant_detail` from golden test framework
- Regenerate all golden outputs
- Add `--feather` test variant

---

## 8. Golden Test Impact

All 20+ golden test scenarios must be regenerated. The golden framework in `test_golden_output.py` needs:

1. **Remove** `quant_detail` comparison
2. **Update** transcript DataFrame expected columns (add coords, flags, nrna_id; drop mrna_high_conf, synthetics)
3. **Update** gene DataFrame expected columns (add coords, gene_type; drop mrna_high_conf)
4. **Update** loci DataFrame expected columns (add nrna, n_annotated_transcripts, n_nrna_entities; rename gdna_init → gdna_prior)
5. **Add** nRNA DataFrame golden comparison
6. **Update** scalars JSON (no changes expected unless nrna_em_count computation changes)
7. **File format:** golden files switch from feather to TSV (or keep feather for fast comparison)

Run `pytest tests/ --update-golden` after all steps.

---

## 9. Execution Order

| Step | Depends on | Estimated scope |
|------|-----------|-----------------|
| 1. Index nRNA linkage | — | `transcript.py` + `index.py` + tests |
| 2. summary.json redesign | — | `cli.py` + `calibration.py` + `config.py` + tests |
| 3. quant.tsv redesign | Step 1 | `estimator.py` + `cli.py` + golden |
| 4. gene_quant.tsv redesign | — | `estimator.py` + golden |
| 5. quant_nrna.tsv (new) | Step 1 | `estimator.py` + `cli.py` + golden |
| 6. loci.tsv redesign | — | `estimator.py` + `pipeline.py` + golden |
| 7. ZL BAM tag | — | `annotate.py` + C++ + tests |
| 8. CLI + format changes | Steps 2-6 | `cli.py` + test_cli + golden regen |

Steps 1 and 2 can run in parallel. Steps 3-6 depend on Step 1 (for nrna_id linkage). Step 7 is independent. Step 8 is the final integration and golden regeneration.
