# Rigel User Manual

This is the user's reference for Rigel: how to install it, what it needs as input, the five CLI
subcommands (`index`, `quant`, `sim`, `export`, `report`) flag by flag, every output file and its
columns, what the calibration stage produces and its two user-facing knobs, and recipes for common
runs. It does not explain the statistical model or why the design is the way it is: the derivations
live in `docs/EQUATIONS.md` and the design and its rulings in `docs/DESIGN.md`. `rigel <subcommand>
--help` is always the authoritative flag list.

Rigel quantifies RNA-seq alignments while separating RNA from genomic-DNA contamination (gDNA). It
takes a name-sorted BAM and produces per-transcript, per-gene and per-locus abundance estimates,
with the gDNA fraction of each locus solved by a calibration stage that runs before the EM.

---

## Installation

### Bioconda (recommended)

```bash
conda install -c conda-forge -c bioconda rigel
```

### PyPI

```bash
pip install rigel-rnaseq
```

The PyPI package name is `rigel-rnaseq`. The CLI command, import name,
GitHub repo, and Bioconda package are all `rigel`.

### From source

```bash
git clone https://github.com/mkiyer/rigel.git
cd rigel
mamba env create -f mamba_env.yaml
conda activate rigel
pip install --no-build-isolation -e .
```

`rigel report` needs the optional `[report]` extra (`pip install 'rigel-rnaseq[report]'`, or
`conda install -c conda-forge vl-convert-python`); everything else works without it.

Verify the install, and note the two global flags that go before the subcommand:

```bash
rigel --version          # print the installed version and exit
rigel -v quant ...       # -v / --verbose: DEBUG-level logging for any subcommand
```

---

## Input requirements

### Reference

- Genome FASTA with `.fai` index (`samtools faidx genome.fa`)
- Gene annotation GTF (GENCODE format recommended; `exon` records required)
- Optional: a per-base mappability **Zarr store** built by the companion
  `alignable` tool for the same genome + aligner (see
  [Mappability](#mappability-and-the-splice-artifact-blacklist))

### BAM

- **Name-sorted or collated** — not coordinate-sorted
- **`NH` tag** used for multimapper handling when present; otherwise
  multimappers are detected from secondary BAM flags
- Splice-junction strand tag recommended for strand-model accuracy
  (`XS` from STAR/HISAT2, `ts` from minimap2, or `auto` for detection)

Name-sort an existing BAM:

```bash
samtools sort -n -@ 8 -o sample.namesorted.bam sample.bam
```

---

## Quick start

```bash
# Step 1: Build index (once per genome + annotation).
#   --no-mappability opts out of the alignable Zarr store; drop it and pass
#   --alignable-zarr for real genomes (see Mappability below).
rigel index \
    --fasta genome.fa \
    --gtf annotation.gtf \
    --no-mappability \
    -o index/

# Step 2: Quantify
rigel quant \
    --bam sample.namesorted.bam \
    --index index/ \
    -o results/sample/ \
    --threads 8 \
    --seed 42

# Step 3: Inspect outputs
rigel export results/sample/ --format tsv    # optional: feather -> TSV
head results/sample/quant.tsv
cat results/sample/summary.json
```

---

## Commands

Rigel has five subcommands: `index`, `quant`, `sim`, `export`, and `report`.

### rigel index

Builds the reference index from a genome FASTA and GTF. Run once per
genome/annotation combination; the index is shared across all samples.

```bash
rigel index --fasta genome.fa --gtf annotation.gtf --alignable-zarr map.zarr -o index/
```

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Genome FASTA (must have a `.fai` index) |
| `--gtf` | required | Annotation GTF (GENCODE recommended; `exon` records required) |
| `-o`, `--output-dir` | required | Output directory for index files |
| `--alignable-zarr PATH` | — | Alignable Zarr store built for the same genome + aligner. Provides per-base fractional mappability (for gDNA-aware effective length) **and** the splice-junction artifact blacklist (applied at BAM-scan time). Required unless `--no-mappability` is set. |
| `--no-mappability` | off | Explicitly opt out of mappability + splice-blacklist ingestion (synthetic genomes, stranded-only benchmarks). Mutually exclusive with `--alignable-zarr`. |
| `--splice-blacklist-min-count N` | `2` | (Advanced) Minimum unique-fragment support per `(chrom, intron, read_length)` for a sj to enter the blacklist. Lower admits more singletons; higher keeps only the most reproducible artifacts. Ignored under `--no-mappability`. |
| `--nrna-tolerance N` | `20` | Max distance (bp) for clustering transcript start/end sites into shared synthetic nRNA spans |
| `--collapse-duplicate-transcripts` | off | Collapse transcripts sharing identical exon coordinates (keeps the lexicographically-smallest ID). Default: fail with a report. Useful for GENCODE, which has a few byte-identical annotations. |
| `--gtf-parse-mode {strict,warn-skip}` | `strict` | `strict` fails on malformed GTF records; `warn-skip` logs warnings and skips them |
| `--feather-compression {lz4,zstd,uncompressed}` | `lz4` | Feather compression for index files |
| `--no-tsv` | off | Skip writing TSV mirrors of index files |

**Index artifacts** (Feather, plus `.tsv` mirrors unless `--no-tsv`):
`transcripts`, `intervals`, `ref_lengths`, `sj`, `splice_blacklist`, `manifest.json`, and the
two calibration partitions — `regions.feather` (the genome partition into intergenic / intron /
exon pieces) and `edges.feather` (its TSV mirror is `boundaries.tsv`: the contiguous and
splice-junction boundaries between those pieces).

### rigel quant

Runs the full single-pass pipeline: BAM scan → calibration → per-locus EM.

```bash
rigel quant \
    --bam sample.namesorted.bam \
    --index index/ \
    -o results/ \
    --threads 8 \
    --seed 42
```

Every flag is also documented by `rigel <subcommand> --help`.

**Input / output**

| Flag | Default | Description |
|------|---------|-------------|
| `--bam PATH` | — | Name-sorted BAM. Minimap2 and STAR are supported. `NH` used when present; otherwise multimappers come from secondary flags. |
| `--index DIR` | — | Directory of `rigel index` files |
| `-o`, `--output-dir DIR` | — | Output directory for results |
| `--config FILE` | — | YAML config (same keys as CLI options, underscored). CLI flags override it. See [YAML configuration](#yaml-configuration). |
| `--annotated-bam PATH` | — | Write an annotated BAM with per-fragment assignment tags (requires a second BAM pass) |
| `--tsv` | off | Also write `.tsv` mirrors of the quant tables |

**Alignment options**

| Flag | Default | Description |
|------|---------|-------------|
| `--include-multimap` / `--no-include-multimap` | on | Include multimapping reads |
| `--keep-duplicates` / `--no-keep-duplicates` | off | Keep reads marked PCR/optical duplicates |
| `--sj-strand-tag TAG [...]` | `auto` | Splice-junction strand tag. `auto` detects from the first 10,000 reads; use `XS` (STAR/HISAT2), `ts` (minimap2), or list several to try in order (e.g. `XS ts`). |

**Model parameters**

| Flag | Default | Description |
|------|---------|-------------|
| `--seed N` | timestamp | Random seed for reproducibility |
| `--em-iterations N` | `1000` | Maximum EM iterations. Set `0` for unambiguous-only quantification (skip EM). |
| `--em-mode {vbem,map}` | `vbem` | EM variant. `vbem` = Variational Bayes EM (digamma soft updates); `map` = MAP-EM with hard `max(0, n+a-1)` updates. |
| `--assignment-mode {sample,fractional,map}` | `sample` | Post-EM fragment assignment. `sample` draws from the posterior; `fractional` preserves posterior weights; `map` takes the argmax component. |

**Performance**

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | 0 (all cores) | Total thread budget. During scan it splits between scan workers and `--scan-bgzf-threads`; locus EM reuses the same budget (stages run serially). |
| `--scan-bgzf-threads N` | derived | BGZF decompression threads reserved from `--threads` during scan. The default divides the budget by the measured ratio — one decompression thread keeps about eight scan workers fed — so 4 threads spend none on it, 8 spend one and 16 spend two. `0` disables htslib threaded decompression. |
| `--scan-buffer-size GiB` | `2` | Max scan buffer before chunks spill to disk |
| `--tmpdir DIR` | system temp | Directory for buffer spill files |
| `--scan-fragments-per-chunk N` | `1000000` | Buffered fragments per scan chunk before the native scanner hands off to the Python buffer |
| `--scan-read-name-batch-size N` | `512` | (Advanced) Read-name groups per native scanner input-queue item |

**Advanced options**

| Flag | Default | Description |
|------|---------|-------------|
| `--assignment-min-posterior P` | `0.01` | Minimum posterior for a component to be eligible for discrete assignment (map/sample modes) |
| `--em-convergence-delta D` | `1e-6` | Convergence threshold for EM parameter updates |
| `--gdna-em-llr-bias B` | `0.0` | gDNA false-positive-aversion: a log-odds (LLR) bias in nats added to the gDNA component in the locus EM. Positive favors gDNA (trades the gDNA→RNA leak for an RNA→gDNA siphon), e.g. `2.2` ≈ require 9:1 RNA evidence before calling a fragment RNA. |
| `--calib-refit-iters N` | `3` | Number of times calibration re-solves after refitting its population gDNA prior. `0` gives the prior-free first solve only. |
| `--overhang-alpha A` | `0.1` | Per-base overhang penalty in `[0,1]`. `0` = hard gate, `1` = no penalty. |
| `--mismatch-alpha A` | `0.1` | Per-mismatch (`NM` tag) penalty in `[0,1]`. `0` = hard gate, `1` = no penalty. |
| `--pruning-min-posterior P` | `1e-4` | Minimum posterior for candidate pruning. Lower keeps more candidates; `0` disables pruning. |
| `--splicing-anchor-tolerance K` | `3` | Resolver-side splicing-anchor tolerance (bp) around annotated introns. The fractional calibration accumulator does not interpret this value. |
| `--emit-locus-stats` | off | Write per-locus EM convergence profiling (iteration counts, timing, equivalence-class stats) to `locus_stats.feather` |

### rigel sim

Generates a small synthetic scenario for benchmarking and development: a random
genome, a GTF, simulated reads aligned with minimap2 and collated with samtools,
and a ready-built index. The scenario is defined by a YAML file; the three CLI
numeric options provide defaults that the YAML overrides.

```bash
rigel sim --config scenario.yaml -o sim_out/
```

| Flag | Default | Description |
|------|---------|-------------|
| `--config FILE` | required | YAML scenario definition |
| `-o`, `--output-dir DIR` | required | Output directory for scenario artifacts |
| `--genome-length N` | `5000` | Genome length in bases (overridden by YAML) |
| `--seed N` | `42` | Random seed (overridden by YAML) |
| `--num-reads N` | `1000` | Number of fragments to simulate (overridden by YAML) |

The YAML schema is what `sim_command` in `src/rigel/cli.py` reads: optional
top-level `name`, `ref_name` (default `chr1`), `genome_length`, `seed`,
`n_fragments`, the read-model keys `frag_mean`, `frag_std`, `frag_min`,
`frag_max`, `read_length`, `error_rate`, and a `genes` list. Each gene has
`gene_id`, `strand`, optional `gene_name` / `gene_type`, and a `transcripts`
list whose entries carry `t_id`, `exons` (0-based half-open `[start, end]`
pairs), and optional `abundance` (default 100) and `nrna_abundance` (nascent
RNA laid across the gene span, default 0):

```yaml
genome_length: 5000
n_fragments: 1000
genes:
  - gene_id: g1
    strand: "+"
    transcripts:
      - {t_id: t1, exons: [[100, 300], [500, 700]], abundance: 100}
```

### rigel export

Converts every `.feather` output file in a results directory to TSV or
Parquet. Useful as a post-processing step or when Parquet is preferred
downstream.

```bash
rigel export results/sample/ --format tsv
rigel export results/sample/ --format parquet
```

| Argument / flag | Default | Description |
|-----------------|---------|-------------|
| `output_dir` (positional) | required | Directory containing `.feather` files from `rigel quant` |
| `-f`, `--format {tsv,parquet}` | `tsv` | Output format |

### rigel report

Builds a single self-contained HTML QC report from a `rigel quant` output
directory. It reads only the files `quant` already wrote (`summary.json` plus the
companion tables — see [Output files](#output-files)), so reports can be built or
rebuilt later, in bulk, without re-quantifying.

```bash
rigel report results/sample/                         # → results/sample/report.html
rigel report results/sample/ -o qc.html --title "Sample 42"
```

| Argument / flag | Default | Description |
|-----------------|---------|-------------|
| `output_dir` (positional) | required | A `rigel quant` output directory (must contain `summary.json`) |
| `-o`, `--output` | `<output_dir>/report.html` | Destination HTML path |
| `--title` | `Rigel QC · <sample>` | Document title |

The report covers alignment fates, fragment composition (including the splice-artifact
blacklist status), the strand model, per-category fragment-length distributions, the
mRNA / nRNA / gDNA split, capture on-target enrichment, a genome-wide gDNA-density track, and
a searchable gene-expression table, in one offline file.

Without the `[report]` extra (see Installation) the report still builds, with the
fragment-length charts omitted.

---

## Mappability and the splice-artifact blacklist

Real genomes have regions where reads cannot map uniquely. Rigel can ingest
a per-base **mappability** track — a Zarr store built by the companion
`alignable` tool (`pip install alignable`) for the *same* genome and aligner —
at index time via `--alignable-zarr`. It is used for two things:

1. **gDNA-aware effective length** — the mappable fraction shortens the
   effective length used in quantification, so unmappable stretches do not
   inflate or deflate abundance.
2. **Splice-junction artifact blacklist** — spurious sj (recurrent
   alignment artifacts) are recorded at index time and treated as unspliced
   during scoring, which the annotated BAM reports per record via the `ZB`
   tag.

For synthetic genomes, stranded-only benchmarks, or any setting where
running `alignable` is unnecessary, pass `--no-mappability` (mutually
exclusive with `--alignable-zarr`). The one advanced knob,
`--splice-blacklist-min-count`, sets the support a junction needs to enter the
blacklist; it is ignored under `--no-mappability`.

---

## YAML configuration

Any `rigel quant` flag can be set in a YAML file via `--config`. Keys are the
flag names with underscores (hyphens are also accepted), plus `bam_file`,
`index_dir`, `output_dir` and `tsv` for the I/O arguments; a key that is not a
`quant` flag is ignored with a warning. Explicit CLI flags always override the
YAML. A resolved `config.yaml` is written to the output directory after each
run for reproducibility.

```yaml
# rigel_params.yaml

# Library handling
include_multimap: true
keep_duplicates: false
sj_strand_tag: [auto]       # [XS] for STAR, [ts] for minimap2, [XS, ts] to try both

# EM algorithm (defaults are suitable for most libraries)
em_iterations: 1000
em_mode: vbem               # vbem | map
assignment_mode: sample     # sample | fractional | map

# Performance
threads: 16
seed: 42
```

```bash
rigel quant --bam sample.bam --index index/ -o results/ --config rigel_params.yaml --threads 4
rigel quant --config results/config.yaml     # rerun a completed analysis from its emitted config
```

---

## Output files

All outputs are written to `--output-dir`. Feather files (`.feather`) are
readable with `pandas.read_feather()` or `pyarrow.feather.read_feather()`.
Pass `--tsv` to also write `.tsv` mirrors, or convert afterward with
`rigel export`.

`rigel quant` writes:

| File | Contents |
|------|----------|
| `quant.feather` | Per-transcript abundance (mRNA + nRNA rows) |
| `gene_quant.feather` | Gene-level aggregation |
| `nrna_quant.feather` | nRNA entity estimates |
| `loci.feather` | Per-locus EM summary |
| `summary.json` | Run-level QC manifest + calibration scalars (`schema_version` 2) |
| `fragment_lengths.feather` | Raw fragment-length histograms, tidy `(category, length, count)` |
| `calibration_track.feather` | Per-region gDNA solution: `(ref, start, end, gdna_mass, rna_mass, gdna_density, gdna_frac)` |
| `calibration_track.bedgraph` | Per-region gDNA density as a genome-browser track (IGV / UCSC) |
| `gdna_density_kde.feather` | The fitted gDNA-density curve `(log_rho, log_density, density)` — the landscape prior, a capture diagnostic |
| `gdna_density_regions.feather` | The training-region rug for that curve `(log_rho, kind)` (downsampled) |
| `locus_stats.feather` | Per-locus EM convergence profiling — only with `--emit-locus-stats` |
| `config.yaml` | Resolved run configuration (reproducibility) |

The `calibration_*` and `gdna_density_*` files are written only when calibration
runs and, for the density curve, only when the gDNA-density prior was fit (enough
training regions). Build the HTML report from all of the above with `rigel report`.
A `config.yaml` is always written, recording all resolved parameters and I/O
paths; rerun the exact analysis with `rigel quant --config results/config.yaml`.

### quant.feather / quant.tsv

Transcript-level abundance estimates. Only annotated (non-synthetic)
transcripts are included; synthetic nRNA spans appear in `nrna_quant`.

| Column | Description |
|--------|-------------|
| `transcript_id` | Transcript ID |
| `gene_id` | Parent gene ID |
| `gene_name` | Gene name/symbol |
| `gene_type` | Biotype from GTF (e.g. `protein_coding`, `lncRNA`) |
| `ref` | Chromosome/contig |
| `strand` | Strand (`+` or `-`) |
| `start`, `end` | Genomic coordinates (0-based, half-open) |
| `length` | Spliced exonic length (bp) |
| `effective_length` | Length after fragment-length correction; the length `tpm` is normalised by, never contracted by capture |
| `em_effective_length` | The yield the EM reads: `effective_length` contracted by the capture efficiency of the transcript's own bases (fragments per unit abundance; equal to `effective_length` when the library carries no enriched gDNA mode). `count / em_effective_length` is the capture-corrected abundance, and a zero-count transcript's abundance is below about `3 / em_effective_length` |
| `locus_id` | Locus ID; `-1` if not placed in an EM locus |
| `nrna_id` | ID of the parent nRNA entity (`"."` if none or if the row is itself an nRNA entity) |
| `is_basic` | `1` if the transcript has the GENCODE `basic` tag |
| `is_mane` | `1` if MANE Select or MANE Plus Clinical |
| `is_nrna` | `1` if single-exon (indistinguishable from nascent RNA) |
| `count` | Fragment count assigned to this transcript |
| `count_unambig` | Deterministically assigned fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Fragment count from annotated-spliced fragments |
| `nrna_parent_count` | Parent nRNA entity's count (0 if `is_nrna=1` or no parent) |
| `tpm` | Transcripts per million (annotated transcripts only) |
| `tpm_total_rna` | TPM over all RNA (annotated + synthetic nRNA); comparable to `tpm` in `nrna_quant` |
| `posterior_mean` | Mean EM posterior over units assigned to this transcript |

### gene_quant.feather / gene_quant.tsv

Gene-level abundance aggregated from annotated (non-synthetic) transcripts.
`mature_count` / `nascent_count` split that over the gene's multi-exon vs
single-exon (`is_nrna`) annotated transcripts (`count = mature_count +
nascent_count`). Synthetic nRNA spans are gene-neutral (many-to-many) and are
reported per-entity in `nrna_quant`, not attributed to genes here.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene ID |
| `gene_name` | Gene name/symbol |
| `gene_type` | Biotype from GTF |
| `ref`, `strand`, `start`, `end` | Gene genomic coordinates |
| `n_transcripts` | Number of annotated transcripts |
| `locus_id` | Primary locus assigned to the gene |
| `effective_length` | Abundance-weighted mean effective transcript length |
| `count` | Total fragment count (sum of annotated transcript counts) |
| `mature_count` | Count over multi-exon annotated transcripts |
| `nascent_count` | Count over single-exon (`is_nrna`) annotated transcripts |
| `count_unambig` | Deterministic fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Spliced fragment count |
| `tpm` | Gene-level TPM |

### nrna_quant.feather / nrna_quant.tsv

nRNA entity estimates. Each row is one nRNA entity (a single-exon
transcript, annotated or synthetic) that has at least one multi-exon mRNA
child. Standalone single-exon genes with no multi-exon overlap appear in
`quant.feather` with `is_nrna=1` instead.

| Column | Description |
|--------|-------------|
| `nrna_id` | nRNA entity transcript ID |
| `gene_id` | Overlapping gene ID |
| `gene_name` | Gene name |
| `ref`, `strand`, `start`, `end` | Entity genomic coordinates |
| `length` | Entity length (bp) |
| `effective_length` | Effective length after fragment-length correction |
| `locus_id` | Locus containing this entity |
| `is_synthetic` | `1` if synthesized by rigel |
| `n_contributing_transcripts` | Number of multi-exon transcripts merged into this span |
| `n_mrna` | Number of multi-exon mRNA transcripts associated with the entity |
| `mrna_count` | Sum of mRNA children's fragment counts |
| `count` | Entity's own fragment count |
| `tpm` | Total-RNA TPM; comparable to `tpm_total_rna` in `quant.feather` |

### loci.feather / loci.tsv

Per-locus EM summary. One row per connected component of overlapping
transcripts. Each locus is an independent EM subproblem with `n_t + 1`
components (one per transcript row + one gDNA).

| Column | Description |
|--------|-------------|
| `locus_id` | Locus ID |
| `locus_span_bp` | Genomic span of the locus (bp) |
| `n_transcripts` | Total transcript rows (mRNA + nRNA spans) |
| `n_annotated_transcripts` | Annotated (non-synthetic) transcript count |
| `n_nrna_entities` | Number of nRNA entities (all single-exon transcripts) |
| `n_genes` | Number of genes |
| `n_em_fragments` | Ambiguous fragments entering EM |
| `count_unambig` | Deterministic (splice-confirmed) mRNA count that bypassed EM |
| `mrna` | Total mRNA count (EM-assigned annotated + unambiguous) |
| `nrna` | Total nRNA count |
| `gdna` | Total gDNA count |
| `total` | `mrna + nrna + gdna` |
| `gdna_rate` | `gdna / total` |
| `gdna_prior_count` | Calibrated per-locus gDNA Dirichlet prior scalar |
| `rna_prior_count` | Calibrated per-locus RNA Dirichlet prior scalar |
| `enable_gdna` | `1` if the gDNA component is enabled for this locus |
| `n_regions_touched` | Number of calibration regions overlapping the locus |
| `multi_locus_region_mass` | Fragment mass in regions shared with other loci |
| `partial_coverage_region_mass` | Fragment mass in partially covered regions |
| `gdna_eff_len_em` | Exposure-adjusted EM effective length of the locus gDNA component |
| `gdna_eff_len_per_bp` | `gdna_eff_len_em / locus_span_bp` diagnostic ratio |

### fragment_lengths.feather / fragment_lengths.tsv

Raw fragment-length histograms in tidy long form — one row per non-empty 1-bp
bin. This is the plotting substrate for the fragment-length distributions;
`summary.json` carries only their summary statistics.

| Column | Description |
|--------|-------------|
| `category` | FL model: `global`, `gdna`, `rna`, or a splice category (`unspliced`, `spliced_annot`, `spliced_unannot`, `spliced_implicit`, `splice_artifact`) |
| `length` | Fragment length (bp) |
| `count` | Fragment weight in this bin (fractional for the `gdna` structural pool) |

Per category, `count` sums to the matching `summary.json` →
`fragment_length.<category>.n_observations`. The `gdna` / `rna` categories are
the deconvolved structural pools used for scoring/calibration; the splice
categories are the scanner's raw per-fragment histograms.

### summary.json

Run-level QC manifest. It is a small, human-readable index — the bulky raw
fragment-length histograms live in the `fragment_lengths.feather` companion
(see above), not in the JSON. `schema_version` (integer) identifies the layout;
the current version is **2**. Top-level keys:

| Key | Contents |
|-----|----------|
| `schema_version` | Manifest schema version (currently `2`) |
| `rigel_version`, `timestamp` | Version and run time |
| `command` | Subcommand, resolved arguments, config-file path |
| `configuration` | All resolved pipeline parameters |
| `input` | Absolute BAM and index paths |
| `alignment_stats` | Two unit families: record-level counts (`total_reads`, `mapped_reads`, `secondary_reads`, `supplementary_reads`, `proper_pairs`, `duplicate_reads`, `qc_fail_reads`) and read-name-group counts (`read_groups`, `unique_reads`, `multimapping_reads`). The report's read-fate bar uses the read-group family so its segments share one denominator |
| `fragment_stats` | Total, genic, intergenic, and chimeric breakdowns; annotated/unannotated SJ counts; a `splice` sub-block with the per-fragment splice-type breakdown (`unspliced`, `spliced_annotated`, `spliced_unannotated`, `spliced_implicit`, `splice_artifact`) plus `sj_blacklisted`, and blacklist provenance `sj_blacklist_size` / `sj_blacklist_loaded` (so `splice_artifact: 0` reads as "detection off — no blacklist in the index" vs "detection on, none found") |
| `strand_model` | Protocol (`R1-sense` / `R1-antisense`), strand specificity, `p_r1_sense`, training count, posterior variance, 95% CI, and a `diagnostics` sub-block comparing the all-exonic model to the spliced-only model (`contamination_gap` — a positive value flags unstranded gDNA contamination) |
| `calibration` | Library-scalar calibration outputs (see below) |
| `gdna_eff_len` | Summary of the per-locus gDNA effective-length series (`em`, `per_bp`) |
| `fragment_length` | Per-category FL **summary statistics only** (`n_observations`, `mean`, `std`, `median`, `mode`, `max_size`, `overflow_count`, `overflow_fraction`) for `global`, `gdna`, `rna`, and each splice category. The raw per-bin histograms are in `fragment_lengths.feather` |
| `quantification` | `n_transcripts`, `n_genes`, `n_loci`, assignment counts, and mRNA/nRNA/gDNA totals + fractions |

The **`calibration`** block holds exactly the library-wide calibration
scalars (it is `null` if calibration did not run):

```jsonc
{
  "calibration": {
    "gdna_density_global":        <float>,  // library-average gDNA density (QC scalar)
    "rna_sense_frac":             <float>,  // kappa: sense-strand RNA fraction
    "gdna_strand_overdispersion": <float>,  // gDNA strand Beta-Binomial overdispersion
    "rna_strand_overdispersion":  <float>,  // RNA strand Beta-Binomial overdispersion
    "n_regions":                  <int>,    // number of calibration regions
    "capture": {                            // present when the gDNA track is informative
      "n_nodes":                 <int>,     // regions with gDNA signal used
      "background_mode_log_rho": <float>,   // dominant density peak by region COUNT (typical region)
      "enriched_mode_log_rho":   <float>,   // high-density peak by gDNA MASS (on-target shoulder)
      "separation_nats":         <float>,   // enriched - background (log peak-to-peak fold)
      "enrichment_factor":       <float>,   // exp(separation_nats) — peak-to-peak fold
      "mass_frac_ontarget":      <float>,   // fraction of gDNA mass in the on-target mode
      "enriched":                <bool>     // a distinct on-target mode was found
    }
  }
}
```

The `capture` block is descriptive only (no pass/fail verdict) and mass-weighted: under hybrid
capture the on-target regions are few but carry the captured gDNA mass, so weighting by mass
surfaces the on-target mode. `enrichment_factor` says how enriched; `mass_frac_ontarget` says how
much of the gDNA is actually on-target.

The RNA and gDNA fragment-length models used by scoring/calibration are
reported under the top-level **`fragment_length`** key (as
`fragment_length.rna` and `fragment_length.gdna`), alongside the global and
per-splice-category FL summaries. The `calibration` block above is the
complete set of calibration keys; a manifest from an older release may carry
keys this manual does not list, which the current version (0.7.1,
`pyproject.toml`) no longer writes.

### Annotated BAM

Produced by `--annotated-bam PATH`. Requires a second pass over the BAM.

The output is a collated BAM with exactly the same records as the input (no drops, no
duplications, every qname in one contiguous run); filtered records (QCFAIL / unmapped / unpaired /
skipped duplicates) are written through unchanged without annotation tags.

| Tag | Type | Description |
|-----|------|-------------|
| `ZT` | string | Assigned transcript ID, or `.` |
| `ZG` | string | Assigned gene ID, or `.` |
| `ZR` | string | Assigned gene name / symbol, or `.` |
| `ZI` | int | Transcript index into rigel reference (`-1` if unassigned) |
| `ZJ` | int | Gene index into rigel reference (`-1` if unassigned) |
| `ZF` | int | Assignment flags bitfield (see below) |
| `ZW` | float | Posterior probability of the assignment |
| `ZC` | string | Input-ambiguity class: `unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`; `.` for records not scored by the EM |
| `ZH` | int | Primary-hit flag: `1` for the winning alignment, `0` otherwise (reflects rigel's EM assignment, which may differ from the aligner's primary FLAG) |
| `ZN` | int | Number of competing candidate components |
| `ZS` | string | Splice type: `spliced_annot`, `spliced_unannot`, `unspliced`, or `unknown` |
| `ZL` | int | Locus ID (`-1` if no locus) |
| `ZB` | int | Number of SJs in **this record's CIGAR** matched to the splice-artifact blacklist and treated as unspliced during scoring. Record-local. Always `0` when no blacklist is loaded. Summed over the BAM it equals the Pass-1 stat `n_sj_blacklisted`. |

#### ZF assignment flags

`ZF` is an 8-bit bitfield describing the fate of the fragment (or read when
unpaired) in the EM pipeline:

| Bit | Mask | Flag | Meaning |
|-----|------|------|---------|
| 0 | 0x01 | `is_resolved` | Fragment was scored and assigned by the EM |
| 1 | 0x02 | `is_mrna` | Assigned to a spliced mRNA transcript |
| 2 | 0x04 | `is_gdna` | Assigned to the gDNA component (EM gDNA or intergenic) |
| 3 | 0x08 | `is_nrna` | Assigned to a nascent-RNA component |
| 4 | 0x10 | `is_synthetic` | Assigned to a rigel-generated synthetic nRNA span |
| 5 | 0x20 | `is_intergenic` | No annotation overlap (always paired with `is_gdna`) |
| 6 | 0x40 | `is_chimeric` | Chimeric fragment, skipped before scoring |
| 7 | 0x80 | `is_multimapper_dropped` | Multimapper dropped when `--include-multimap` is off |

Canonical `ZF` values (the only ones produced):

| ZF (hex) | Dec | Meaning |
|----------|-----|---------|
| `0x03` |   3 | mRNA assigned (multi-exon transcript) |
| `0x09` |   9 | nRNA assigned (annotated single-exon transcript) |
| `0x19` |  25 | nRNA assigned (rigel-synthesized span) |
| `0x05` |   5 | gDNA assigned by EM (annotation-overlapping) |
| `0x25` |  37 | gDNA assigned by intergenic fallback |
| `0x40` |  64 | Chimeric, not scored |
| `0x80` | 128 | Multimapper dropped, not scored |

Pysam usage:

```python
zf = read.get_tag("ZF")
is_resolved   = (zf & 0x01) != 0
is_mrna       = (zf & 0x02) != 0
is_gdna       = (zf & 0x04) != 0
is_nrna       = (zf & 0x08) != 0
is_synthetic  = (zf & 0x10) != 0
is_intergenic = (zf & 0x20) != 0
is_chimeric   = (zf & 0x40) != 0
is_mm_dropped = (zf & 0x80) != 0
```

---

## Calibration

Calibration is the middle stage of `rigel quant`: it runs once, in-process, between the BAM
scan and the per-locus EM, on the fragment tallies the scanner deposits during its single pass
(no extra read of the BAM). Rigel's model has three fragment origins — spliced RNA, unspliced
RNA (which includes nascent pre-mRNA), and genomic DNA contamination, which is unspliced,
unstranded and genome-wide (and enriched on probed exons under hybrid capture). A spliced
fragment is certified RNA; calibration's job is to split each genomic region's *unspliced* mass
into sense RNA, antisense RNA and gDNA, using the strand tilt of the counts, what neighbouring
regions and junctions say, and a library-wide prior on gDNA density fitted from the data. The
theory is in `docs/EQUATIONS.md` and the design in `docs/DESIGN.md`.

### What calibration produces

- **Library scalars**, in `summary.json` → `calibration`: `gdna_density_global`,
  `rna_sense_frac` (the sense fraction κ), the gDNA and RNA strand overdispersions, and the
  `capture` block when a distinct on-target gDNA-density mode is found.
- **A per-locus Dirichlet prior** — `gdna_prior_count` and `rna_prior_count` in `loci.feather`
  — which sets the gDNA-vs-RNA split each locus's EM starts from, plus the gDNA component's
  effective length. RNA is distributed among transcripts by the EM, not by calibration.
- **A genome-browser track** of the solved gDNA density, `calibration_track.feather` /
  `.bedgraph`, and the fitted density curve in `gdna_density_kde.feather`.

The per-locus priors are priors: a decisive per-locus likelihood overrides them, so a locus that
disagrees with calibration generally wins.

### The one knob

Advanced; the default suits standard libraries.

- `--calib-refit-iters N` (default `3`) — how many times calibration re-fits its population
  gDNA-density prior on the current solve and re-solves. `0` gives the prior-free first solve.

`--gdna-em-llr-bias` is not a calibration knob: it biases the gDNA component in the locus EM
downstream. The remaining `CalibrationConfig` fields (for example `sweep_logodds_step`, the log-odds
lattice's step, which is dimensionless and guarantees every slot's composition to within 1.25 % of its
mass) are not exposed on the CLI or in the YAML.

### When to suspect calibration is misfiring

- A high `gdna_fraction` in `summary.json` → `quantification` on a presumed-clean library, or an
  implausible `gdna_density_global`, usually means the strand signal is weak (near-unstranded
  data), so the solve is leaning on the prior. Inspect `strand_model.strand_specificity`.
- Heavy mass on nRNA components in `nrna_quant.feather` for a short-fragment library is usually
  an identifiability limit within that locus rather than a calibration error.

---

## Snakemake integration

A typical pipeline aligns (STAR, HISAT2 or minimap2), name-sorts the BAM, builds the Rigel index
once, and runs `rigel quant` per sample. The Snakefile below does the last three; every flag it
uses is documented in [Commands](#commands).

```python
# Snakefile — config.yaml supplies samples, genome_fasta, genome_gtf, rigel_index,
# and optionally alignable_zarr (else the index is built with --no-mappability)
# and rigel_config (a shared rigel quant YAML, see YAML configuration).

rule all:
    input:
        expand("results/{sample}/quant.feather", sample=config["samples"]),
        expand("results/{sample}/summary.json",  sample=config["samples"]),

rule rigel_index:
    input:  fasta = config["genome_fasta"], gtf = config["genome_gtf"],
    output: directory(config["rigel_index"]),
    params: zarr = config.get("alignable_zarr", ""),
    log:    "logs/rigel_index.log",
    shell:
        """
        rigel index --fasta {input.fasta} --gtf {input.gtf} \
            $([ -n "{params.zarr}" ] && echo "--alignable-zarr {params.zarr}" || echo "--no-mappability") \
            -o {output} > {log} 2>&1
        """

rule namesort_bam:
    input:  bam = "aligned/{sample}.bam",
    output: bam = temp("namesorted/{sample}.bam"),
    threads: 4
    shell:  "samtools sort -n -@ {threads} -o {output.bam} {input.bam}"

rule rigel_quant:
    input:  bam = "namesorted/{sample}.bam", index = config["rigel_index"],
    output:
        quant   = "results/{sample}/quant.feather",
        summary = "results/{sample}/summary.json",
    params: outdir = "results/{sample}", config = config.get("rigel_config", ""),
    log:    "logs/rigel_{sample}.log",
    threads: 8
    shell:
        """
        rigel quant --bam {input.bam} --index {input.index} -o {params.outdir} \
            --threads {threads} --seed 42 \
            $([ -n "{params.config}" ] && echo "--config {params.config}") > {log} 2>&1
        """
```

Load the per-sample tables into one matrix with pandas:

```python
import pandas as pd
samples = ["sample1", "sample2", "sample3"]
gene_quant = pd.concat(
    [pd.read_feather(f"results/{s}/gene_quant.feather").assign(sample=s) for s in samples],
    ignore_index=True,
)
count_matrix = gene_quant.pivot(index="gene_id", columns="sample", values="count")
tpm_matrix   = gene_quant.pivot(index="gene_id", columns="sample", values="tpm")
```

---

## Supported aligners

| Aligner | SJ strand tag | Notes |
|---------|--------------|-------|
| STAR | `XS` | Name-grouped output by default |
| HISAT2 | `XS` | Standard RNA-seq workflow |
| minimap2 | `ts` | Long-read and splice-aware mode |

Use `--sj-strand-tag auto` (default) to detect the tag automatically, or list
several tags to try in order (e.g. `--sj-strand-tag XS ts`).

---

## Recipes

### STAR-aligned BAM

STAR produces name-grouped BAMs with the `XS` tag:

```bash
rigel quant \
    --bam     STAR_out/Aligned.out.bam \
    --index   index/ \
    -o        results/ \
    --sj-strand-tag XS \
    --threads 16 \
    --seed    42
```

### Fully reproducible run, with TSV tables

```bash
rigel quant --bam sample.bam --index index/ -o results/ --seed 42 --threads 1 --tsv
rigel quant --config results/config.yaml      # rerun later from the emitted config
rigel export results/ --format tsv            # or convert existing Feather outputs afterward
```

`--seed` makes the post-EM assignment sampling deterministic; `--threads 1` makes the output
bit-reproducible.

### Inspect read assignments

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --annotated-bam results/annotated.bam

# Count fragments assigned to gDNA by the EM (ZF=5)
samtools view -F 256 results/annotated.bam \
    | awk '{ for(i=12;i<=NF;i++) if($i=="ZF:i:5") count++ } END { print count }'
```

### Exclude multimappers, skip the EM, profile the EM

```bash
rigel quant --bam sample.bam --index index/ -o results/ --no-include-multimap
rigel quant --bam sample.bam --index index/ -o results/ --em-iterations 0   # unambiguous counts only
rigel quant --bam sample.bam --index index/ -o results/ --emit-locus-stats  # -> locus_stats.feather
```

---

## FAQ

**Why does PyPI use `rigel-rnaseq`?**
The name `rigel` is taken on PyPI. The CLI, import, GitHub repo, and Bioconda
package are all `rigel`.

**What sort order does the BAM need?**
Name-sorted or collated — not coordinate-sorted.
Use `samtools sort -n -o out.bam in.bam` or `samtools collate`.

**What is the difference between `mrna` and `nrna`?**
`mrna` counts fragments consistent with mature, spliced RNA.
`nrna` counts fragments consistent with unspliced, nascent pre-mRNA.
Total RNA-seq libraries contain both. Use `mrna` for standard differential
expression; use `nrna` for transcription dynamics.

**Why are nRNA counts non-integer or shared?**
nRNA is estimated on shared genomic spans. Counts are pro-rated across all
transcripts that map to the same span. Isoforms with identical start and end
coordinates share one nRNA component.

**Why is `strand_specificity` near 0.5?**
Rigel trains the strand model from annotated spliced fragments only.
Unstranded libraries, or libraries with few informative splice reads, stay
near 0.5.

**How does calibration use strand information?**
The per-strand counts are calibration's own composition signal: gDNA is unstranded, RNA is not, so
a region's strand tilt says how much of its unspliced mass is gDNA, and the count sets how sharp
that reading is. An unstranded or low-count region reads a flat tilt and its composition is then
set by its neighbours and the library-wide gDNA-density prior. There is no strand-mode switch.

**Do I need the alignable Zarr store?**
For real genomes it is recommended — it provides gDNA-aware effective length
and the splice-artifact blacklist. For synthetic genomes or stranded-only
benchmarks, pass `--no-mappability` at index time.

**How much memory does Rigel use?**
The scan buffer defaults to 2 GiB (`--scan-buffer-size`); overflow spills to
disk under `--tmpdir`. Typical paired-end libraries (50–100M read pairs) fit
within the default.

**Can I run Rigel on a cluster?**
Yes. Each `rigel quant` is independent. Parallelise at the sample level via
your scheduler; use `--threads` for intra-job parallelism.

**When should I use `--annotated-bam`?**
For read-level inspection, debugging, or assignment validation. It requires a
second BAM pass and adds runtime overhead.

**Does Rigel support single-end reads?**
Single-end reads are handled but less thoroughly tested than paired-end.
Fragment-length estimation uses alignment length rather than insert size.

**Can a transcript's EM component die and never recover?**
In VBEM mode Rigel uses SQUAREM acceleration, which can overshoot and push a component's
Dirichlet alpha to zero. `VBEM_SQUAREM_PRIOR_FLOOR` (a compile-time constant in
`src/rigel/native/em_solver.cpp`, equal to `EM_LOG_EPSILON` ≈ 1e-300) is the minimum alpha after
each SQUAREM step. It is deliberately deep: `digamma` of it is so negative that the component
receives zero responsibility in the E-step and stays dead, while a component with genuine read
support recovers because the M-step adds real observations to its prior. It has no effect in
MAP-EM mode.
