# Rigel User Manual

This is the user's reference for Rigel: how to install it, what it needs as input, the five CLI
subcommands (`index`, `quant`, `sim`, `export`, `report`) flag by flag, every output file and its
columns, what the calibration stage produces and its one user-facing knob, and recipes for common
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

The charts in `rigel report` need the optional `[report]` extra (`pip install 'rigel-rnaseq[report]'`,
or `conda install -c conda-forge vl-convert-python`); everything else works without it.

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
- Optional: an **alignable Zarr store** built by the companion `alignable` tool for the
  same genome + aligner, the source of the splice-artifact blacklist (see
  [Mappability](#mappability-and-the-splice-artifact-blacklist))

### BAM

- **Paired-end** — the scan stops at the first unpaired record. Rigel pairs a read pair's alignments itself, for any
  aligner: each alignment of R1 is tried with each alignment of R2, and STAR's `HI` tag is not read. A pairing is used
  only if the two mates can be one molecule: on one reference, on opposite strands, their 5' ends facing each other once
  each 5' soft clip is walked back (a mate may run past the other's 5' end, as when a read was trimmed at its 5' end:
  the molecule spans both mates), with a fragment length from 1 to `max_frag_length` along the genome, or along one
  annotated transcript that overlaps both mates and holds every junction the pair observes; and, if it is spliced, only
  if one annotated transcript overlapping both mates holds every junction it observes. Rigel does not assemble
  transcripts: an alignment with a junction no transcript has, or with junctions no one transcript joins, is neither RNA
  it can place nor gDNA. An alignment that fails is dropped; a read pair left with none is withheld, counted in
  `summary.json` and written to the annotated BAM as unresolved
- **Name-sorted or collated** — a BAM whose header says `SO:coordinate` is refused, since its mates are not
  adjacent
- **Aligned to the index's genome** — references are matched by name: a BAM none of whose references is in the
  index, or one holding a reference at another length, is refused; references beyond the index (decoys,
  alternates, spike-ins) are named in a warning, and an alignment there is dropped, since Rigel can place nothing
  outside its index: a read pair with another alignment keeps it, and one with none is withheld, counted in
  `summary.json` (`off_index`) and in no pool
- **Complete** — a truncated or corrupt BAM is refused, never quantified on the records read before the damage
- **Chimeric alignments are not quantified** — a read the aligner writes in pieces (its record carries `SA`, or its
  5′ end is hard-clipped, or it has a supplementary record) is no molecule Rigel can place: such an alignment is
  dropped, a read pair with another alignment keeps it, and one with none is withheld, counted in `summary.json`
  (`chimeric`). STAR with `--chimOutType WithinBAM` writes a chimeric read in place of its linear alignment; with
  `--chimOutType Junctions` the BAM keeps the linear alignment and the chimeric junctions go to their own file. STAR's
  mate merging (`--peOverlapNbasesMin`) merges two mates on a chance 10-base match, and can write an ordinary pair as
  a chimera
- **`NH` tag** used for multimapper handling when present; otherwise
  multimappers are detected from secondary BAM flags. A read pair whose `NH` exceeds the records reported for a mate
  was truncated by the aligner (counted, and warned of: STAR reports every alignment with `--outSAMmultNmax -1`)
- **An `NM` tag on every record**, its edit distance, which the EM's mismatch penalty reads (STAR writes it with
  `--outSAMattributes NH HI AS nM NM`). A record without one reads as 0 edits, which can favor it over an alignment
  whose mismatches are known, so such records are counted and warned of; records of one alignment that disagree on
  `NM` are read with the largest
- **A splice-junction strand tag, for a stranded library** (`XS` from HISAT2, or from STAR with
  `--outSAMstrandField intronMotif`; `ts` from minimap2). Without one no junction trains the strand model
  and the library is quantified as unstranded; `--sj-strand-tag auto` reads whichever tag each spliced record
  carries, and cannot supply one the aligner did not write

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
| `--alignable-zarr PATH` | — | Alignable Zarr store built for the same genome + aligner. Provides the splice-junction artifact blacklist (applied at BAM-scan time). Required unless `--no-mappability` is set. Read before anything is written, so a store that cannot be read, or whose blacklist names none of the FASTA's references, refuses the build and leaves an existing index in place. |
| `--no-mappability` | off | Explicitly opt out of the alignable store and its splice blacklist (synthetic genomes, stranded-only benchmarks). Mutually exclusive with `--alignable-zarr`. |
| `--splice-blacklist-min-count N` | `2` | (Advanced) Minimum unique-fragment support per `(chrom, intron, read_length)` for a sj to enter the blacklist. Lower admits more singletons; higher keeps only the most reproducible artifacts. Ignored under `--no-mappability`. |
| `--nascent-tolerance N` | `20` | Max distance (bp) for clustering multi-exon transcripts' start and end sites into the shared synthetic spans that model their nascent RNA |
| `--collapse-duplicate-transcripts` | off | Collapse transcripts sharing identical exon coordinates (keeps the lexicographically-smallest ID; `transcript_aliases` maps each collapsed ID to the one kept). Default: fail with a report. Useful for GENCODE, which has a few byte-identical annotations. |
| `--gtf-parse-mode {strict,warn-skip}` | `strict` | `strict` fails on malformed GTF records; `warn-skip` logs warnings and skips them. Malformed: a line that cannot be parsed, a feature that ends before it starts, an exon without a `transcript_id`, an exon off its transcript's reference, strand or gene, an exon overlapping another of its transcript, and a transcript on a reference the FASTA lacks or past its end. Exons with no base between them are one exon. |
| `--feather-compression {lz4,zstd,uncompressed}` | `lz4` | Feather compression for index files |
| `--no-tsv` | off | Skip writing TSV mirrors of index files |

**Index artifacts** (Feather, plus `.tsv` mirrors unless `--no-tsv`):
`transcripts`, `intervals`, `ref_lengths`, `sj`, `splice_blacklist` (each junction's anchors and
the alignable count behind it), `transcript_aliases` (`dropped_t_id → kept_t_id` for every
transcript `--collapse-duplicate-transcripts` folded into an identical one; empty otherwise),
`manifest.json`, and the two calibration partitions — `regions.feather` (the genome partition into
intergenic / intron / exon pieces) and `edges.feather` (its TSV mirror is `boundaries.tsv`: the
contiguous and splice-junction boundaries between those pieces). The on-disk format is version 10;
`rigel quant` refuses an index of any other version, so an index built by an earlier rigel is
rebuilt with `rigel index`.

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
| `--tsv` / `--no-tsv` | off | Also write `.tsv` mirrors of the quant tables |

**Alignment options**

| Flag | Default | Description |
|------|---------|-------------|
| `--include-multimap` / `--no-include-multimap` | on | Include multimapping reads |
| `--keep-duplicates` / `--no-keep-duplicates` | off | Keep reads marked PCR/optical duplicates |
| `--sj-strand-tag TAG [...]` | `auto` | Splice-junction strand tag. `auto` reads `XS` and `ts` on every spliced record, and a record whose two tags disagree reads as ambiguous; name `XS` (HISAT2, or STAR with `--outSAMstrandField intronMotif`) or `ts` (minimap2) to read only that tag, both to take the first present on each record (`XS ts`, either order), or `none` to read no tag (the library is then quantified as unstranded). Any other value is refused. |

**Model parameters**

| Flag | Default | Description |
|------|---------|-------------|
| `--seed N` | `0` | Seed of the `sample` assignment's draw; another seed is another draw from the same posterior. |
| `--em-iterations N` | `1000` | EM iteration budget. The EM runs accelerated (SQUAREM) steps of three EM updates each, up to a third of this budget, until it converges. At least one step always runs, so `0` does not skip the EM: every fragment is still assigned. |
| `--em-mode {vbem,map}` | `vbem` | EM variant. `vbem` = Variational Bayes EM (digamma soft updates); `map` = MAP-EM: each step adds the calibration's gDNA and RNA pseudocounts to the components' fragment counts and normalises. |
| `--assignment-mode {sample,fractional}` | `sample` | Post-EM fragment assignment. `sample` gives every fragment to one transcript (or gDNA), count first: each transcript's fractional count is rounded within its EM locus, then each fragment is drawn from its own posterior toward the transcripts still short of their count, so whole counts equal the rounded fractional counts. `fractional` keeps each fragment's posterior weights. |

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
| `--em-convergence-delta D` | `1e-6` | Convergence threshold for EM parameter updates |
| `--calib-refit-iters N` | `3` | Number of times calibration re-solves after refitting its population gDNA prior. `0` gives the prior-free first solve only, and with it no capture correction: no gDNA landscape is fitted, so every capture efficiency is 1. |
| `--sweep-block-slots N` | `1000` | Calibration's working set: the chain is solved one locus block at a time, the pieces between intergenic regions merged up to this many slots per block. Performance only — the answer is the same for every value; smaller blocks use less memory per sweep. |
| `--overhang-alpha A` | `0.1` | Per-base overhang penalty in `[0,1]`. `0` = hard gate, `1` = no penalty. |
| `--mismatch-alpha A` | `0.1` | Per-mismatch (`NM` tag) penalty in `[0,1]`. `0` = hard gate, `1` = no penalty. |
| `--pruning-min-posterior P` | `1e-4` | Minimum posterior for candidate pruning, in [0, 1). Lower keeps more candidates; `0` disables pruning. |
| `--intron-overhang-tolerance K` | `3` | Intron-overhang tolerance (bp) at annotated introns: how many bases an aligned read end may lie inside an intron, at the edge of a gap the reads did not sequence or at a fragment's outer end, as when an aligner runs a read a base or two across a splice site. An annotated intron in an unsequenced gap may overlap up to K aligned bases at either end and still count as lying in the gap, and a transcript is contradicted only by an intron holding more than K aligned bases. Along a witness transcript it decides whether two mates are one molecule; it also shapes the candidates' lengths and the gap hypotheses calibration arbitrates, so it changes calibration's tally. Junctions are matched exactly, and it is not the splice blacklist's anchor test, which it never changes. |
| `--emit-locus-stats` / `--no-emit-locus-stats` | off | Write per-locus EM convergence profiling (iteration counts, timing, equivalence-class stats) to `locus_stats.feather` |

### rigel sim

Generates a small synthetic scenario for benchmarking and development: a random
genome, a GTF, simulated reads aligned with minimap2 and name-sorted with samtools
(both must be on `PATH`), and a ready-built index. The scenario is defined by a YAML file; the three CLI
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
`frag_max`, `read_length`, `error_rate`, and a `genes` list of at least one gene. Each gene has
`gene_id`, `strand`, optional `gene_name` / `gene_type`, and a `transcripts`
list whose entries carry `t_id`, `exons` (0-based half-open `[start, end]`
pairs), and optional `abundance` (default 100) and `nrna_abundance` (nascent
RNA laid across the gene span, default 0). Any other key, at any level, is an error
reported before anything is written:

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
blacklist status), the strand model, per-category fragment-length distributions, the split of
the fragments into annotated transcripts, nascent RNA, gDNA (EM), intergenic RNA and intergenic
gDNA, calibration's capture answer (the gDNA bulk density and the saturated
level over it), a genome-wide gDNA-density track, and
a searchable gene-expression table, in one offline file.

Without the `[report]` extra (see Installation) the report still builds, with every chart
omitted.

---

## Mappability and the splice-artifact blacklist

Real genomes have regions where reads cannot map uniquely, and aligners emit
recurrent spurious splice junctions there. Rigel can ingest a Zarr store built
by the companion `alignable` tool (not yet on PyPI or Bioconda: install it from a checkout with
`pip install -e /path/to/alignable`) for the *same*
genome and aligner, at index time via `--alignable-zarr`. From it the index
records a **splice-junction artifact blacklist**: those spurious sj are treated
as unspliced during scoring, which the annotated BAM reports per record via the
`ZB` tag. The store's per-base mappability is not read.

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
`quant` flag is an error that lists the keys this version knows. Explicit CLI flags always override the
YAML. A resolved `config.yaml` is written to the output directory after each
run for reproducibility.

```yaml
# rigel_params.yaml

# Library handling
include_multimap: true
keep_duplicates: false
sj_strand_tag: [auto]       # both tags on every record; [XS] for STAR, [ts] for minimap2, [XS, ts] first present

# EM algorithm (defaults are suitable for most libraries)
em_iterations: 1000
em_mode: vbem               # vbem | map
assignment_mode: sample     # sample | fractional

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
| `quant.feather` | Per-transcript abundance: one row per annotated transcript, single-exon included |
| `gene_quant.feather` | Per-gene abundance, summed over each gene's annotated transcripts |
| `nascent_quant.feather` | Nascent RNA: one row per synthetic span Rigel adds to model a multi-exon transcript's unspliced RNA |
| `loci.feather` | Per-locus EM summary |
| `summary.json` | Run-level QC manifest + calibration scalars (`schema_version` 4) |
| `fragment_lengths.feather` | Raw fragment-length histograms, tidy `(category, length, count)` |
| `calibration_track.feather` | Per-region gDNA solution: `(ref, start, end, count_gdna, count_rna, gdna_density, fraction_gdna)` |
| `calibration_track.bedgraph` | Per-region gDNA density as a genome-browser track (IGV / UCSC) |
| `locus_stats.feather` | Per-locus EM convergence profiling — only with `--emit-locus-stats` |
| `config.yaml` | Resolved run configuration (reproducibility) |

The `calibration_*` files are written only when calibration runs. Build the HTML
report from all of the above with `rigel report`.
A `config.yaml` is always written, recording all resolved parameters and I/O
paths; rerun the exact analysis with `rigel quant --config results/config.yaml`
on the rigel version that wrote it (its header names the version). Another
version refuses a key it does not have, naming the version that wrote the file:
every 0.7.1 `config.yaml` carries four keys removed in 0.8.0 and one renamed
(`splicing_anchor_tolerance`, now `intron_overhang_tolerance`).

**How the outputs name RNA.** Every table, `summary.json` key and report label puts RNA in one of two
categories. **Annotated** RNA is RNA on a transcript of the reference annotation, single-exon transcripts
included: a single-exon transcript does not splice, so it is reported with the annotation. **Nascent** RNA is
RNA on the synthetic spans Rigel adds over the genomic extent of multi-exon transcripts to model their
unspliced RNA, introns included. Beside them, `gdna` is genomic DNA, `rna` all RNA, and `intergenic` a place:
a fragment no transcript or span overlaps. A `count` is fragments and a `fraction` a share.

### quant.feather / quant.tsv

Transcript-level abundance estimates: one row per annotated transcript, single-exon included. The
synthetic nascent spans are rows of `nascent_quant`.

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
| `em_effective_length` | The yield the EM reads: `effective_length` scaled by the capture efficiency of every stretch the transcript's fragments land on (its exon pieces, the boundaries between them and its splice junctions), each read from the genomic DNA captured there or, for a junction, beside it. An efficiency is relative to uncaptured gDNA (the bulk), so `em_effective_length` equals `effective_length` off the probes and in a library without capture, and rises above it where capture enriched the transcript — to about `gdna_capture_saturation` times for a fully captured one. `count / em_effective_length` is the capture-corrected abundance, and a zero-count transcript's abundance is below about `3 / em_effective_length` |
| `locus_id` | Locus ID; `-1` if not placed in an EM locus |
| `nascent_id` | The transcript's synthetic nascent span, a row of `nascent_quant`; `"."` when it has none: a single-exon transcript, or a transcript whose span an annotated single-exon transcript covers, which reports that RNA as annotated on its own row |
| `is_basic` | True if the transcript has the GENCODE `basic` tag |
| `is_mane` | True if MANE Select or MANE Plus Clinical |
| `is_single_exon` | True for a single-exon transcript. It does not splice, so its RNA is annotated RNA like any other transcript's |
| `count` | Fragment count assigned to this transcript |
| `count_unambig` | Deterministically assigned fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Fragment count from spliced fragments |
| `count_nascent` | Fragment count of the transcript's nascent span (`nascent_id`); `0` when it has none. Every transcript the span models repeats it, so never sum it over isoforms |
| `tpm` | Transcripts per million, over the annotated transcripts |
| `tpm_total_rna` | TPM over all RNA, annotated transcripts and nascent spans together; comparable to `tpm_total_rna` in `nascent_quant` |
| `posterior_mean` | Mean EM posterior over units assigned to this transcript |

### gene_quant.feather / gene_quant.tsv

Gene-level abundance, summed over each gene's annotated transcripts: a gene's RNA is its annotated
`count`. The nascent spans are gene-neutral (one span can model the transcripts of more than one gene),
so they are reported per span in `nascent_quant` and never attributed to a gene here.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene ID |
| `gene_name` | Gene name/symbol |
| `gene_type` | Biotype from GTF |
| `ref`, `strand`, `start`, `end` | Gene genomic coordinates |
| `n_transcripts` | Number of annotated transcripts |
| `locus_id` | Primary locus assigned to the gene |
| `effective_length` | Abundance-weighted mean effective transcript length |
| `count` | Fragment count, the sum of the gene's annotated transcripts' counts |
| `count_unambig` | Deterministic fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Spliced fragment count |
| `tpm` | Gene-level TPM |

### nascent_quant.feather / nascent_quant.tsv

Nascent RNA, one row per synthetic span. Rigel clusters the start and end sites of the multi-exon
transcripts within `--nascent-tolerance` and adds one span over each merged extent, which holds the unspliced
RNA, introns included, of the transcripts it models. Where an annotated single-exon transcript already
covers the extent, Rigel adds no span: that transcript models it, and its RNA is annotated, a row of
`quant.feather`.

| Column | Description |
|--------|-------------|
| `nascent_id` | Span ID |
| `effective_length` | Effective length after fragment-length correction (the length TPM divides by) |
| `em_effective_length` | The effective length the EM reads (capture-aware) |
| `locus_id` | Locus containing this span |
| `count` | The span's own fragment count: its nascent RNA |
| `n_annotated` | Number of annotated transcripts the span models |
| `count_annotated` | Sum of those transcripts' fragment counts |
| `tpm_total_rna` | TPM over all RNA, annotated transcripts and nascent spans together; comparable to `tpm_total_rna` in `quant.feather` |

A span's coordinates are its row in the index's `transcripts` table; spans are gene-neutral, so they carry
no gene.

**Overlapping spans.** Rigel makes one span per merged extent, so a gene with several isoform clusters
carries several overlapping same-strand spans. The EM cannot tell overlapping spans apart, so how their
nascent RNA splits among them is arbitrary: read nascent RNA per gene, or per family of overlapping spans,
never per span. `quant.feather`'s `nascent_id` and `count_nascent` point each transcript at its own span,
so a transcript's `count_nascent` can read near 0 while an overlapping span holds the family's fragments.

### loci.feather / loci.tsv

Per-locus EM summary. One row per connected component of overlapping
transcripts. Each locus is an independent EM subproblem with one component per annotated
transcript and nascent span, plus a gDNA component where `has_gdna`.

| Column | Description |
|--------|-------------|
| `locus_id` | Locus ID |
| `locus_span_bp` | Genomic span of the locus (bp) |
| `n_transcripts` | Number of annotated transcripts |
| `n_nascent` | Number of nascent spans |
| `n_genes` | Number of genes |
| `n_em_fragments` | Ambiguous fragments entering EM |
| `count_unambig` | Deterministic (splice-confirmed) count on the annotated transcripts, which bypassed the EM |
| `count_annotated` | Fragment count on the annotated transcripts, single-exon included: EM-assigned plus `count_unambig` |
| `count_nascent` | Fragment count on the nascent spans |
| `count_gdna` | gDNA fragment count from the locus EM |
| `count` | `count_annotated + count_nascent + count_gdna` |
| `fraction_gdna` | `count_gdna / count` |
| `count_gdna_prior` | Calibrated gDNA pseudocount of the locus EM's Dirichlet prior |
| `count_rna_prior` | Calibrated RNA pseudocount of the locus EM's Dirichlet prior |
| `has_gdna` | True if the locus has a gDNA component: some unspliced fragment in it can be gDNA |
| `gdna_em_effective_length` | Exposure-adjusted EM effective length of the locus gDNA component |

### fragment_lengths.feather / fragment_lengths.tsv

Raw fragment-length histograms in tidy long form — one row per 1-bp bin from the
first non-empty bin to the last. This is the plotting substrate for the fragment-length distributions;
`summary.json` carries only their summary statistics.

| Column | Description |
|--------|-------------|
| `category` | FL model: `global`, `gdna`, `rna`, or one of the five observed pools (`pool_gdna_intergenic`, `pool_gdna_intronic`, `pool_gdna_intron_exon`, `pool_gdna_intergenic_exon`, `pool_rna_spliced`) |
| `length` | Fragment length (bp); the row at the maximum fragment length the scan models (1000 bp) holds every fragment at least that long |
| `count` | Fragment weight in this bin: whole fragments for `global` and the observed pools, fractional for the fitted `gdna` and `rna` models |

Per category, `count` sums to the matching `summary.json` →
`fragment_length.<category>.n_observations`. The `gdna` / `rna` categories are
the length models scoring and calibration use; the five pools are the observed
histograms they are built from.

### summary.json

Run-level QC manifest. It is a small, human-readable index — the bulky raw
fragment-length histograms live in the `fragment_lengths.feather` companion
(see above), not in the JSON. `schema_version` (integer) identifies the layout;
the current version is **4**. Top-level keys:

| Key | Contents |
|-----|----------|
| `schema_version` | Manifest schema version (currently `5`) |
| `rigel_version`, `timestamp` | Version and run time |
| `command` | Subcommand, resolved arguments, config-file path |
| `configuration` | All resolved pipeline parameters |
| `input` | Absolute BAM and index paths |
| `alignment_stats` | The alignment records: `total_reads`; `mapped_reads`, `unmapped_reads`, `secondary_reads`, `supplementary_reads`, `duplicate_reads`, `qc_fail_reads`, each the records carrying its flag; and `proper_pairs`, the primary records flagged as a proper pair. A `records` sub-block gives each record one status, summing to `total_reads`: `filtered` (unmapped, QC-failed, or a duplicate when duplicates are skipped), `used` (an alignment of a read pair's placement), `unused` (in a placed read pair but in none of its usable alignments: one the pair rule rejected, one that pairs with nothing, a supplementary record) and `unplaced` (in a read pair withheld or excluded by `--no-include-multimap`). An `sj_strand_tag` sub-block counts the spliced records the scan read (secondary and supplementary included) by what their junction-strand tags gave: `tagged` (a strand), `untagged` (none) and `conflicting` (XS and ts disagree, read as ambiguous), with `invalid` the tag values that were neither `+` nor `-`, each read as absent. `nm_tag` and `nh_tag` count the records without the tag (`missing`) and with it held as anything but an integer (`invalid`), each read as absent (`NM` 0, `NH` 1); `nm_tag.conflicting` counts the alignments whose records disagree on `NM`, each read with the largest |
| `read_pairs` | One outcome per read pair (a read name's records), summing to `total`, decided in this order: `excluded` (every record filtered, as a pair whose mates are both unmapped is); `withheld`, with no usable alignment, by reason: `off_index` when one of its alignments lies on a reference the index does not hold, else `chimeric` when one of its alignments has a mate aligned in pieces, else `unannotated_splice` when one of its alignments is a molecule whose junctions no one annotated transcript holds, else `unpaired`, which `unpaired_by_reason` splits by the furthest check of the pair rule any of its alignments reached (`missing_mate`, `reference`, `strand`, `outward`, `length`); `excluded_by_option` (`--no-include-multimap`); and `placed`, `unique` (one usable alignment) or `multimapping` (more). The report's read-fate bar shows every read pair by outcome |
| `fragment_stats` | The fragments (the placed read pairs): total, genic and intergenic counts; the annotated SJ count; a `splice` sub-block with the per-fragment splice-type breakdown (`unspliced`, `spliced_annotated`, `spliced_implicit`, `splice_artifact`; a spliced fragment's junctions are always annotated) plus `sj_blacklisted` (blacklisted junctions × records, multimapper records included: the annotated BAM's `ZB` summed), and blacklist provenance `sj_blacklist_size` / `sj_blacklist_loaded` (so `splice_artifact: 0` reads as "detection off — no blacklist in the index" vs "detection on, none found"); `hits_rejected`, the pairings of an R1 alignment with an R2 alignment tried and dropped, by reason (`missing_mate`, `reference`, `strand`, `outward`, `length`: the two mates are not one molecule; `contradictory`, `novel_junction`, `unsupported`: no one annotated transcript overlapping both mates holds the junctions it observes; `off_index`: both mates on one reference the index does not hold; `chimeric`: a mate aligned in pieces); and `nh_above_reported`, the read pairs whose `NH` exceeds the records reported for a mate (a truncated output, whose multiplicity is read from the alignments present) |
| `strand_model` | Protocol (`R1-sense` / `R1-antisense`, and `read1_sense`), strand specificity, `p_r1_sense`, training count, the per-junction strand table's QC (`sj`: `n_sj`, `n_observations` and the junctions' depth, `depth_median`, `depth_p90`, `depth_p99`, `depth_max`), posterior variance, 95% CI; the junction fit behind κ — `kappa` (the spliced reads' sense fraction, `pooled` over every junction and `genuine` over the junctions the three-class fit credits to RNA), `junction_classes` (each class's share — `genuine`, `artifact`, `reversed` — `by_junction` and `by_fragment`; shallow splice artifacts are most of the junctions and few of the fragments; `null` where the strand channel is dead, since the classes are not identified at κ = ½ and `genuine` then equals `pooled`, the κ calibration reads), `strand_live` (whether the protocol preserves strand, the decision calibration's strand channel reads) and `strand_overdispersion` (0, `binomial`, for gDNA and RNA alike); and a `diagnostics` sub-block comparing the all-exonic model to the spliced-only model (`contamination_gap` — a positive value flags unstranded gDNA contamination) |
| `calibration` | Library-scalar calibration outputs (see below) |
| `gdna_em_effective_length` | Summary of the loci table's column of that name: its `min`, `median`, `p95`, `max` |
| `fragment_length` | Per-category FL **summary statistics only** (`n_observations`, `mean`, `std`, `median`, `mode`, `max_size`, `overflow_count`, `overflow_fraction`) for `global`, `gdna`, `rna` and the five pools. The raw per-bin histograms are in `fragment_lengths.feather` |
| `quantification` | `n_transcripts` (the annotated transcripts) and `n_nascent` (the nascent spans), `n_genes`, `n_loci`, assignment counts (`n_unambig_assigned`, `n_em_assigned`, and `n_gated_out`: the fragments a hard scoring gate dropped — `--overhang-alpha 0` or `--mismatch-alpha 0` — which no pool holds), and three pools, each a `count_*` beside its `fraction_*` of the three together: `count_annotated` (the annotated transcripts), `count_nascent` (the nascent spans) and `count_gdna` (gDNA: the EM's plus the intergenic fragments, which reach no locus and are unspliced, since a spliced alignment is used only when a transcript holds its junctions), so `fraction_annotated + fraction_nascent + fraction_gdna = 1`. `count_gdna_em` and `fraction_gdna_em` are the EM's genic gDNA alone, so the intergenic gDNA is `count_gdna − count_gdna_em`, and `fraction_intergenic` its share; their count is `fragment_stats.intergenic` |

The **`calibration`** block holds exactly the library-wide calibration
scalars (it is `null` if calibration did not run):

```jsonc
{
  "calibration": {
    "gdna_density_global":        <float>,       // library-average gDNA density (QC scalar)
    "gdna_bulk_density":          <float|null>,  // the uncaptured gDNA level capture is read against, gDNA fragments/bp; null = no landscape
    "gdna_capture_saturation":    <float>,       // the saturated gDNA level over that bulk; 1 = nothing above it
    "fraction_rna_sense":         <float>,       // kappa: sense-strand RNA fraction
    "n_regions":                  <int>,         // number of calibration regions
    "n_boundaries":               <int>,         // number of contiguous boundaries between adjacent regions
    "n_sj":                       <int>          // number of splice junctions
  }
}
```

`gdna_bulk_density` and `gdna_capture_saturation` are calibration's capture answer. Calibration fits a
landscape of gDNA densities over the genome; its largest basin is the bulk, the gDNA off the probes, and
each region's capture efficiency is its own gDNA density read against that bulk, so `em_effective_length`
is scaled by those efficiencies and rises above `effective_length` where capture enriched the transcript.
`gdna_capture_saturation` is the landscape's most massive enriched level over the bulk, the level at which
capture saturates. A library without capture has a landscape of one level: its saturation is 1, every
efficiency is 1, and every `em_effective_length` equals `effective_length`. `gdna_bulk_density` is `null`
only when calibration fitted no landscape at all. The report shows both, descriptively: no pass/fail
verdict.

**Known limitation: a library with almost no gDNA.** There the bulk sits near the landscape's floor, and
the few regions where calibration calls gDNA read as enormous capture against it. A transcriptome-only
library with no capture panel reads `gdna_capture_saturation` of about 2 × 10⁷ and some thousands of
regions above e× the bulk, and the `em_effective_length` of the transcripts on those regions is inflated
by up to that factor. A saturation far above 1 on a library that was not captured is this limitation, not
capture. With no gDNA evidence at all (an empty or nearly empty library) the landscape's smoothing alone can read
a saturation near 10 and lengthen every `em_effective_length` by a few percent.

**Known limitation: an unstranded capture library.** Without a strand channel, gDNA on the captured exons is told
from RNA only through the neighbouring regions and the fragment counts, and some of it is read as RNA, most of
it as nascent RNA, so the gDNA fraction is underestimated. Stranded capture libraries, and unstranded libraries
without capture, do not share this.

The RNA and gDNA fragment-length models used by scoring/calibration are
reported under the top-level **`fragment_length`** key (as
`fragment_length.rna` and `fragment_length.gdna`), alongside the global model
and the five pools. The `calibration` block above is the complete set of
calibration keys; a manifest from an older release may carry keys this manual
does not list, which the current version no longer writes.

### Annotated BAM

Produced by `--annotated-bam PATH`. Requires a second pass over the BAM.

The output is a BAM with exactly the same records as the input (no drops, no duplications), collated —
every qname in one contiguous run. Filtered records (QC-failed, unmapped or a skipped duplicate) are written
through unchanged without annotation tags. Rigel's tags replace any tag of the same name the aligner wrote
(HISAT2's `ZS:i`, for one). The quantification tables are written before this pass begins.

| Tag | Type | Description |
|-----|------|-------------|
| `ZT` | string | Assigned transcript ID, or `.` |
| `ZG` | string | Assigned gene ID, or `.` |
| `ZR` | string | Assigned gene name / symbol, or `.` |
| `ZI` | int | Transcript index into rigel reference (`-1` if unassigned) |
| `ZJ` | int | Gene index into rigel reference (`-1` if unassigned) |
| `ZF` | int | Assignment flags bitfield (see below) |
| `ZW` | float | Posterior probability of the assignment |
| `ZC` | string | Input-ambiguity class: `unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`; for a read pair withheld with no usable alignment (`ZF` 0), `off_index` when one of its alignments lies on a reference the index does not hold, else `chimeric` when one of its alignments has a mate aligned in pieces, else `unannotated_splice` when one of its alignments is a molecule whose junctions no one annotated transcript holds, else `unpaired`; `.` for other records not scored by the EM |
| `ZH` | int | Primary-hit flag: `1` for the winning alignment, `0` otherwise (reflects rigel's EM assignment, which may differ from the aligner's primary FLAG) |
| `ZN` | int | Number of competing candidate components |
| `ZS` | string | Splice type: `spliced_annot`, `spliced_implicit`, `splice_artifact`, `unspliced`, or `unknown` |
| `ZL` | int | Locus ID (`-1` if no locus) |
| `ZB` | int | Number of SJs in **this record's CIGAR** matched to the splice-artifact blacklist and treated as unspliced during scoring. Record-local. Always `0` when no blacklist is loaded. Summed over the BAM it equals `summary.json`'s `fragment_stats.splice.sj_blacklisted`. |

#### ZF assignment flags

`ZF` is an 8-bit bitfield describing the fate of the fragment in the pipeline. Bits 1, 2 and 4 say which
kind of candidate won, named for the category the tables report it under:

| Bit | Mask | Flag | Meaning |
|-----|------|------|---------|
| 0 | 0x01 | `is_resolved` | Fragment was assigned to a pool, by the EM or deterministically |
| 1 | 0x02 | `is_annotated` | Assigned to an annotated transcript: any reference transcript, single-exon included |
| 2 | 0x04 | `is_gdna` | Assigned to gDNA (the EM's gDNA component, or unspliced intergenic) |
| 3 | 0x08 | — | Retired: never set, so no earlier `ZF` value is reused with a new meaning |
| 4 | 0x10 | `is_nascent` | Assigned to a synthetic nascent span |
| 5 | 0x20 | `is_intergenic` | No annotation overlap, so unspliced gDNA: always with `is_gdna` |
| 6 | 0x40 | — | Retired: never set |
| 7 | 0x80 | `is_multimapper_dropped` | Multimapper dropped when `--include-multimap` is off |

Canonical `ZF` values (the only ones produced):

| ZF (hex) | Dec | Meaning |
|----------|-----|---------|
| `0x00` |   0 | Unresolved: scored, but assigned no pool |
| `0x03` |   3 | Annotated: assigned to a reference transcript, single-exon included |
| `0x11` |  17 | Nascent: assigned to a synthetic nascent span |
| `0x05` |   5 | gDNA, assigned by the EM (annotation-overlapping) |
| `0x25` |  37 | Intergenic gDNA (unspliced, no annotation overlap) |
| `0x80` | 128 | Multimapper dropped, not scored |

Pysam usage:

```python
zf = read.get_tag("ZF")
is_resolved   = (zf & 0x01) != 0
is_annotated  = (zf & 0x02) != 0
is_gdna       = (zf & 0x04) != 0
is_nascent    = (zf & 0x10) != 0
is_intergenic = (zf & 0x20) != 0
is_mm_dropped = (zf & 0x80) != 0
```

---

## Calibration

Calibration is the middle stage of `rigel quant`: it runs once, in-process, between the BAM
scan and the per-locus EM, on the fragment tallies the scanner deposits during its single pass
(no extra read of the BAM). Rigel's model has three fragment origins — spliced RNA, unspliced
RNA (nascent RNA among it), and genomic DNA contamination, which is unspliced,
unstranded and genome-wide (and enriched on probed exons under hybrid capture). A spliced
fragment is certified RNA; calibration's job is to split each genomic region's *unspliced* mass
into sense RNA, antisense RNA and gDNA, using the strand tilt of the counts, what neighbouring
regions and junctions say, and a library-wide prior on gDNA density fitted from the data. The
theory is in `docs/EQUATIONS.md` and the design in `docs/DESIGN.md`.

### What calibration produces

- **Library scalars**, in `summary.json` → `calibration`: `gdna_density_global`,
  `fraction_rna_sense` (the sense fraction κ), the uncaptured gDNA level `gdna_bulk_density` that capture
  is read against, and the saturated level over it, `gdna_capture_saturation`. The junction fit behind
  κ and the strand overdispersion (0, binomial) are in `strand_model`.
- **A per-locus Dirichlet prior** — `count_gdna_prior` and `count_rna_prior` in `loci.feather`
  — which sets the gDNA-vs-RNA split each locus's EM starts from, plus the gDNA component's
  effective length. RNA is distributed among transcripts by the EM, not by calibration.
- **A genome-browser track** of the solved gDNA density, `calibration_track.feather` /
  `.bedgraph`.

The per-locus priors are priors: a decisive per-locus likelihood overrides them, so a locus that
disagrees with calibration generally wins.

### The one knob

Advanced; the default suits standard libraries.

- `--calib-refit-iters N` (default `3`) — how many times calibration re-fits its population
  gDNA-density prior on the current solve and re-solves. `0` gives the prior-free first solve, which fits no
  gDNA landscape and so applies no capture correction (every capture efficiency is 1).

The remaining `CalibrationConfig` fields (for example `sweep_logodds_step`, the log-odds
lattice's step, which is dimensionless and guarantees every slot's composition to within 1.25 % of its
mass) are not exposed on the CLI or in the YAML, save `--sweep-block-slots`, a memory setting that moves
no number.

Every fixed number the code reads — the numerical methods' settings, the QC thresholds and the resource
budgets — is in `rigel.config.CONSTANTS`, by component, each documented with why it has its value. They are
not user-configurable.

### When to suspect calibration is misfiring

- A high `fraction_gdna` in `summary.json` → `quantification` on a presumed-clean library, or an
  implausible `gdna_density_global`, usually means the strand signal is weak (near-unstranded
  data), so the solve is leaning on the prior. Inspect `strand_model.strand_specificity`.
- Heavy nascent RNA in `nascent_quant.feather` for a short-fragment library is usually
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
| STAR | `XS` | Written only with `--outSAMstrandField intronMotif`; name-grouped output by default |
| HISAT2 | `XS` | Standard RNA-seq workflow |
| minimap2 | `ts` | Short-read splice mode (`-ax splice:sr`) |

`--sj-strand-tag auto` (default) reads both tags on every spliced record: one present gives that record's strand,
two that agree give it too, and two that disagree read as ambiguous, so that record trains no strand. A tag holding
anything but `+` or `-` reads as absent. Name one tag to read only that one, or several to take the first present on
each record (e.g. `--sj-strand-tag XS ts`). The run logs what the tags gave and records it in `summary.json`
(`alignment_stats.sj_strand_tag`); `config.yaml` records the spec. Junction strands come only from the tag: if no
spliced read carries one, no junction trains the strand model, the library is quantified as unstranded, and Rigel
warns.

---

## Recipes

### STAR-aligned BAM

STAR produces name-grouped BAMs, and writes the `XS` tag only when run with `--outSAMstrandField intronMotif`:

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
rigel quant --bam sample.bam --index index/ -o results/ --threads 1 --tsv
rigel quant --config results/config.yaml      # rerun later from the emitted config
rigel export results/ --format tsv            # or convert existing Feather outputs afterward
```

The assignment's draw is seeded with a fixed default (`--seed` picks another draw), and `--threads 1`
makes the output bit-reproducible. With more threads the BAM scan's workers take batches in whatever
order they finish, which moves the last bits of the scan's sums from run to run and can move fragments
between isoforms that share them; use `--threads 1` where two runs must match bit for bit. The locus EM
itself gives the same answer at any thread count.

### Inspect read assignments

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --annotated-bam results/annotated.bam

# Tally the fragments by ZF value (3 annotated, 17 nascent, 5 gDNA by the EM, 37 intergenic gDNA,
# 33 intergenic RNA): first mate, primary, non-supplementary — one record per fragment
samtools view -f 64 -F 0x900 results/annotated.bam \
    | awk '{ for(i=12;i<=NF;i++) if($i ~ /^ZF:i:/) n[substr($i,6)]++ } END { for(v in n) print v, n[v] }' \
    | sort -n
```

### Exclude multimappers, profile the EM

```bash
rigel quant --bam sample.bam --index index/ -o results/ --no-include-multimap
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

**What is the difference between annotated transcripts and nascent RNA?**
An annotated transcript's `count` is the RNA assigned to that reference transcript, spliced fragments
and unspliced ones alike; a single-exon transcript's RNA is annotated too, since it does not splice.
Nascent RNA is the RNA assigned to the synthetic spans Rigel adds over the extent of multi-exon
transcripts: RNA that has not spliced, introns included. Total RNA-seq libraries carry both. Use the
annotated counts (`quant`, `gene_quant`) for standard differential expression, and nascent RNA
(`nascent_quant`, or `count_nascent` per transcript) for transcription dynamics, read per gene or per
family of overlapping spans (see Overlapping spans under `nascent_quant`).

**Why is `count_nascent` shared?**
One synthetic span models every multi-exon transcript whose start and end sites cluster within
`--nascent-tolerance`. `nascent_quant` holds its count once; `quant.feather`'s `count_nascent` repeats it on
every transcript the span models, so never sum it over isoforms. Counts are fractional only under
`--assignment-mode fractional`.

**Why is `strand_specificity` near 0.5?**
Rigel trains the strand model from annotated spliced fragments only. On a
stranded library it leaves out the junctions that are not genuine RNA: a
splice artifact (genomic DNA the aligner wrote as spliced) reads near 50/50,
and a reversed junction reads the opposite way. Unstranded libraries, or
libraries with few informative splice reads, stay near 0.5. The fit's
numerical settings are `rigel.config.CONSTANTS.junction_fit`, each documented there;
they are not exposed on the CLI.

**How does calibration use strand information?**
The per-strand counts are calibration's own composition signal: gDNA is unstranded, RNA is not, so
a region's strand tilt says how much of its unspliced mass is gDNA, and the count sets how sharp
that reading is. An unstranded or low-count region reads a flat tilt and its composition is then
set by its neighbours and the library-wide gDNA-density prior. There is no strand-mode setting: whether
the protocol preserves strand is decided from the junctions (`summary.json`'s `strand_model.strand_live`).

**Do I need the alignable Zarr store?**
For real genomes it is recommended — it provides the splice-artifact blacklist
(its per-base mappability is not read). For synthetic genomes or stranded-only
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
No: `rigel quant` requires paired-end alignments and stops at the first unpaired record.

**How are whole counts made, and why no per-fragment cutoff?**
In `sample` mode every fragment goes to exactly one transcript (or to gDNA), and the counts come first.
Each EM locus rounds its transcripts' fractional counts (the expected number of fragments each one
produced) so that they still sum to the locus's fragments; a transcript expecting under about half a
fragment gets none. Then each fragment is drawn from its own posterior, favouring the transcripts still
short of their count, and a fragment left with every candidate full is moved along a short chain of
swaps onto the counts. The whole counts therefore equal the fractional counts rounded, whatever the seed;
the seed decides only which fragment goes where, and `--annotated-bam` records that choice read by read.
There is deliberately no cutoff on one fragment's posterior: a minor isoform can hold 0.5 % of every
fragment it shares and still have produced dozens of them, and only its total says it is real.

**Can a transcript's EM component die and never recover?**
Yes, and only the EM's own step may do it. A component at the floor (`EM_LOG_EPSILON` ≈ 1e-300 in
`src/rigel/native/em_solver.cpp`) receives zero responsibility in the E-step and no share of the prior,
so it stays dead unless it holds fragments no other component can explain. Rigel accelerates the EM
with SQUAREM, which jumps ahead along the EM's path; a jump that would carry a still-live component
below the floor is shortened toward the plain EM step instead of being cut off there, in VBEM and
MAP-EM alike. Before 2026-09-25 such components were clamped to the floor, and which of several
transcripts sharing the same fragments survived then depended on where the EM started. VBEM keeps a
milder dependence of its own, accepted as a known limit: where fragments are genuinely ambiguous, which
of several near-equivalent explanations it settles on can depend on its starting point. Rigel always
starts from the same coverage-weighted point, so one input and one version give the EM one answer at any thread
count, and the alternatives measured no more accurate against simulated truth; a change that moves only the starting
point, as a new version may, can move such fragments between isoforms, unspliced RNA and gDNA.
