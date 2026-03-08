# Rigel User Manual

Complete usage guide for Rigel, a Bayesian RNA-seq transcript
quantification tool with joint mRNA, nascent RNA, and genomic DNA
deconvolution.

**Version:** 0.1.0

---

## Table of contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Input requirements](#input-requirements)
4. [Commands](#commands)
   - [rigel index](#rigel-index)
   - [rigel quant](#rigel-quant)
   - [rigel sim](#rigel-sim)
5. [Parameters reference](#parameters-reference)
   - [Required arguments](#required-arguments)
   - [Common options](#common-options)
   - [EM algorithm](#em-algorithm)
   - [Scoring penalties](#scoring-penalties)
   - [nRNA fraction hierarchical shrinkage](#nrna-fraction-hierarchical-shrinkage)
   - [gDNA rate hierarchical shrinkage](#gdna-rate-hierarchical-shrinkage)
   - [Strand model](#strand-model)
   - [Output options](#output-options)
   - [Performance](#performance)
6. [Configuration file](#configuration-file)
7. [Output files](#output-files)
   - [quant.tsv / quant.feather](#quanttsv--quantfeather)
   - [gene_quant.tsv / gene_quant.feather](#gene_quanttsv--gene_quantfeather)
   - [loci.tsv / loci.feather](#locitsv--locifeather)
   - [quant_detail.tsv / quant_detail.feather](#quant_detailtsv--quant_detailfeather)
   - [summary.json](#summaryjson)
   - [config.json](#configjson)
   - [Annotated BAM](#annotated-bam)
8. [Supported aligners](#supported-aligners)
9. [Recipes and examples](#recipes-and-examples)
10. [FAQ](#faq)

---

## Overview

Rigel quantifies transcript-level abundances from RNA-seq data while
jointly modeling three species of nucleic acid present in a typical
library:

| Species | Symbol | Description |
|---------|--------|-------------|
| Mature RNA | mRNA | Spliced, fully processed transcripts |
| Nascent RNA | nRNA | Unspliced pre-mRNA captured mid-transcription |
| Genomic DNA | gDNA | Background contamination from genomic DNA |

The pipeline operates in two stages:

1. **BAM scan** — A single-pass C++ scanner reads the name-sorted BAM,
   resolves each fragment against the reference index, trains strand and
   fragment-length models, and buffers fragments in memory.
2. **EM quantification** — Fragments are partitioned into independent loci
   and solved via per-locus Expectation-Maximization with SQUAREM
   acceleration.

---

## Installation

### Bioconda (recommended)

```bash
conda install -c conda-forge -c bioconda rigel
```

### PyPI

```bash
pip install rigel
```

A C++17-capable compiler is required for building the native extension.
On macOS, install Xcode Command Line Tools first:

```bash
xcode-select --install
```

### From source

```bash
git clone https://github.com/mkiyer/rigel.git
cd rigel

# Create and activate the environment
mamba env create -f mamba_env.yaml
conda activate rigel

# Install in editable mode
pip install --no-build-isolation -e .
```

### Verify installation

```bash
rigel --version
```

---

## Input requirements

### Reference files

- **Genome FASTA** — A reference genome FASTA file with an accompanying
  `.fai` index. Generate the index with `samtools faidx genome.fa`.
- **Gene annotation GTF** — A GTF file with gene and transcript annotations.
  GENCODE annotations are recommended. The GTF must contain `gene`,
  `transcript`, and `exon` features.

### BAM file

- Must be **name-sorted** or **collated** (not coordinate-sorted).
  Use `samtools sort -n` or `samtools collate` to re-sort if needed.
- Must contain the **NH tag** for multimapper detection. This tag is
  produced by STAR, HISAT2, minimap2, and most RNA-seq aligners.
- Paired-end reads are expected. Single-end reads should work but are
  less thoroughly tested.

---

## Commands

### rigel index

Build a reference index from genome FASTA and gene annotation GTF.

```
rigel index --fasta <genome.fa> --gtf <annotation.gtf> -o <index_dir>
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Genome FASTA file (must have `.fai` index) |
| `--gtf` | Yes | Gene annotation GTF file |
| `-o`, `--output-dir` | Yes | Output directory for index files |
| `--feather-compression` | No | Compression for Feather files: `lz4` (default), `zstd`, `uncompressed` |
| `--no-tsv` | No | Skip writing TSV mirror files |
| `--gtf-parse-mode` | No | `strict` (default, fails on errors) or `warn-skip` (skip malformed lines) |

**Example:**

```bash
rigel index \
    --fasta ~/ref/GRCh38.p14.genome.fa \
    --gtf ~/ref/gencode.v46.primary_assembly.annotation.gtf \
    -o ~/indices/gencode_v46
```

The index directory will contain interval-based lookup structures
(cgranges) and transcript metadata in Feather format.

---

### rigel quant

Quantify transcript abundances from an aligned BAM file.

```
rigel quant --bam <aligned.bam> --index <index_dir> -o <output_dir> [options]
```

This is the primary command. It performs a single-pass BAM scan followed
by locus-level EM quantification.

**Minimal example:**

```bash
rigel quant \
    --bam sample.bam \
    --index ~/indices/gencode_v46 \
    -o results/sample1/
```

**Full-featured example:**

```bash
rigel quant \
    --bam sample.bam \
    --index ~/indices/gencode_v46 \
    -o results/sample1/ \
    --config my_params.yaml \
    --threads 16 \
    --em-prior-alpha 0.01 \
    --em-prior-gamma 1.0 \
    --confidence-threshold 0.95 \
    --prune-threshold 0.1 \
    --annotated-bam results/sample1/annotated.bam \
    --seed 42
```

See the [Parameters reference](#parameters-reference) section for all
available options.

---

### rigel sim

Generate synthetic test scenarios for benchmarking and validation.

```
rigel sim --config <scenario.yaml> -o <output_dir> [options]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `--config` | Yes | YAML configuration defining the test scenario |
| `-o`, `--output-dir` | Yes | Output directory for scenario artifacts |
| `--genome-length` | No | Genome length in bp (default: 5000, overridden by YAML) |
| `--seed` | No | Random seed (default: 42, overridden by YAML) |
| `--num-reads` | No | Number of fragments to simulate (default: 1000, overridden by YAML) |

The simulation generates:
- A synthetic genome FASTA
- Gene annotation GTF
- Simulated BAM with known ground-truth assignments
- Oracle abundance files for benchmarking

---

## Parameters reference

All parameters apply to the `rigel quant` subcommand. Defaults are from
the config dataclasses. The resolution order is:
**CLI flag** > **YAML config file** > **dataclass default**.

### Required arguments

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted or collated BAM file (must have NH tag) |
| `--index` | Directory containing rigel index files |
| `-o`, `--output-dir` | Output directory for quantification results |

### Common options

#### Read filtering

| Flag | Default | Description |
|------|---------|-------------|
| `--include-multimap` / `--no-include-multimap` | yes | Include multimapping reads. Detected via NH tag (STAR) or secondary BAM flag (minimap2). |
| `--keep-duplicates` / `--no-keep-duplicates` | no | Keep reads marked as PCR/optical duplicates. |

#### Splice-junction strand

| Flag | Default | Description |
|------|---------|-------------|
| `--sj-strand-tag` | `auto` | BAM tag(s) for splice-junction strand. `auto` detects from the BAM header. Use `XS` for STAR, `ts` for minimap2, or list multiple to check in order. |

### EM algorithm

| Flag | Default | Description |
|------|---------|-------------|
| `--em-iterations` | 1000 | Maximum EM iterations. Set to 0 for unambiguous-only quantification. |
| `--em-prior-alpha` | 0.01 | Flat Dirichlet pseudocount per eligible EM component. Small values allow the coverage-weighted OVR prior to dominate. |
| `--em-prior-gamma` | 1.0 | OVR (One Virtual Read) prior scale factor. Controls coverage-weighted prior strength. Set to 0 to disable. |
| `--em-convergence-delta` | 1e-6 | Convergence threshold for parameter updates (max absolute change in θ). |
| `--confidence-threshold` | 0.95 | Minimum RNA-normalized posterior for high-confidence deterministic assignment. |
| `--prune-threshold` | 0.1 | Post-EM evidence-ratio threshold. Components with zero unambiguous evidence and data/prior ratio below this value are zeroed and the EM re-runs. Set to -1 to disable. |

### Scoring penalties

Per-fragment log-likelihood penalties applied during candidate scoring.
Values are in [0, 1] where 0 is a hard gate (zero probability) and 1
disables the penalty.

| Flag | Default | Description |
|------|---------|-------------|
| `--overhang-alpha` | 0.01 | Per-base overhang penalty. Each base of soft-clip or overhang multiplies the score by α. |
| `--mismatch-alpha` | 0.1 | Per-mismatch penalty. Each edit (NM tag) multiplies the score by α. |
| `--gdna-splice-penalty-unannot` | 0.01 | gDNA score penalty for unannotated splice junctions. |

### nRNA fraction hierarchical shrinkage

The nascent RNA fraction per transcript uses a three-level empirical Bayes
hierarchy: **global → locus-strand → TSS-group → transcript**. At each
level, a shrinkage pseudo-count κ pulls the estimate toward the parent.

| Flag | Default | Description |
|------|---------|-------------|
| `--tss-window` | 200 | Fuzzy TSS grouping window (bp). Transcripts whose 5' ends lie within this distance are grouped. |
| `--nrna-frac-kappa-global` | auto | κ pulling locus-strand estimate toward global. Auto-estimated via Method of Moments. |
| `--nrna-frac-kappa-locus` | auto | κ pulling TSS-group estimate toward locus-strand. |
| `--nrna-frac-kappa-tss` | auto | κ pulling transcript estimate toward TSS-group. |
| `--nrna-frac-kappa-min` | 2.0 | Lower clamp for MoM-estimated κ. |
| `--nrna-frac-kappa-max` | 200.0 | Upper clamp for MoM-estimated κ. |
| `--nrna-frac-kappa-fallback` | 5.0 | Fallback κ when insufficient features for MoM. |
| `--nrna-frac-kappa-min-obs` | 20 | Minimum features required for MoM estimation; fewer triggers fallback. |

#### MoM minimum evidence thresholds

Groups with fewer fragments than the threshold are excluded from
Method-of-Moments κ estimation at each hierarchy level.

| Flag | Default | Description |
|------|---------|-------------|
| `--nrna-frac-mom-min-evidence-global` | 50 | Min fragment evidence for global MoM κ. |
| `--nrna-frac-mom-min-evidence-locus` | 30 | Min fragment evidence for locus MoM κ. |
| `--nrna-frac-mom-min-evidence-tss` | 20 | Min fragment evidence for TSS-level MoM κ. |

### gDNA rate hierarchical shrinkage

The genomic DNA rate prior uses a two-level hierarchy:
**global → chromosome → locus**.

| Flag | Default | Description |
|------|---------|-------------|
| `--gdna-kappa-chrom` | auto | κ pulling chromosome gDNA rate toward global. Auto-estimated via MoM. |
| `--gdna-kappa-locus` | auto | κ pulling locus gDNA rate toward chromosome. |
| `--gdna-mom-min-evidence-chrom` | 50 | Min fragment evidence for chromosome gDNA MoM κ. |
| `--gdna-mom-min-evidence-locus` | 30 | Min fragment evidence for locus gDNA MoM κ. |

### Strand model

| Flag | Default | Description |
|------|---------|-------------|
| `--strand-prior-kappa` | 2.0 | Beta prior pseudocount κ for the strand model. The prior is Beta(κ/2, κ/2), shrinking toward 0.5 (maximum entropy). Default 2.0 gives a uniform Beta(1, 1) prior. |

### Output options

| Flag | Default | Description |
|------|---------|-------------|
| `--no-tsv` | off | Skip writing human-readable TSV files (produce Feather only). |
| `--annotated-bam` | — | Write an annotated BAM with per-fragment assignment tags to the specified path. Requires a second BAM pass. |

### Performance

| Flag | Default | Description |
|------|---------|-------------|
| `--threads` | 0 (all cores) | Number of threads for BAM scanning and parallel locus EM. 0 = all available cores, 1 = sequential. |
| `--tmpdir` | system temp | Directory for buffer spill files when memory limits are exceeded. Use a fast SSD. |
| `--seed` | timestamp | Random seed for reproducibility. Setting a fixed seed ensures deterministic results with `--threads 1`. |

> **Note — OpenMP thread suppression.**  On import, Rigel sets
> `OMP_NUM_THREADS=1` (via `os.environ.setdefault`) to prevent numpy's
> OpenMP runtime from spawning idle worker threads that compete with
> Rigel's own C++ parallelism for CPU cores and cache.  If you need
> multi-threaded numpy/BLAS operations in the same process, set
> `OMP_NUM_THREADS` to your desired value **before** `import rigel`.

---

## Configuration file

Parameters can be specified in a YAML file passed via `--config`. Keys
use underscores matching the CLI flag names. Explicit CLI flags always
override the config file.

**Example `config.yaml`:**

```yaml
# EM parameters
em_prior_alpha: 0.01
em_prior_gamma: 1.0
em_iterations: 1000
em_convergence_delta: 1.0e-6
confidence_threshold: 0.95
prune_threshold: 0.1

# Scoring
overhang_alpha: 0.01
mismatch_alpha: 0.1

# nRNA fraction shrinkage (auto = Method of Moments)
tss_window: 200
nrna_frac_kappa_global: auto
nrna_frac_kappa_locus: auto
nrna_frac_kappa_tss: auto

# gDNA shrinkage
gdna_kappa_chrom: auto
gdna_kappa_locus: auto

# Read filtering
include_multimap: true
keep_duplicates: false

# Performance
threads: 0
seed: 42
```

**Usage:**

```bash
rigel quant --bam sample.bam --index idx/ -o out/ --config config.yaml
```

Any CLI flag overrides the corresponding YAML key:

```bash
# Use config.yaml but override threads
rigel quant --bam sample.bam --index idx/ -o out/ \
    --config config.yaml --threads 8
```

---

## Output files

All output files are written to the directory specified by `--output-dir`.

### quant.tsv / quant.feather

Per-transcript abundance estimates. One row per transcript.

| Column | Type | Description |
|--------|------|-------------|
| `transcript_id` | string | Transcript identifier from the GTF |
| `gene_id` | string | Gene identifier |
| `gene_name` | string | Gene name (symbol) |
| `mrna` | float | Estimated mRNA (mature RNA) fragment count |
| `nrna` | float | Estimated nRNA (nascent RNA) fragment count |
| `length` | int | Transcript length in bp |
| `effective_length` | float | Effective length after fragment-length correction |
| `tpm` | float | Transcripts Per Million (mRNA-based) |

### gene_quant.tsv / gene_quant.feather

Per-gene aggregated abundances. Transcript-level estimates are summed
within each gene.

| Column | Type | Description |
|--------|------|-------------|
| `gene_id` | string | Gene identifier |
| `gene_name` | string | Gene name (symbol) |
| `mrna` | float | Total mRNA fragment count across all transcripts |
| `nrna` | float | Total nRNA fragment count across all transcripts |
| `tpm` | float | Gene-level TPM (sum of transcript TPMs) |

### loci.tsv / loci.feather

Per-locus summary. A locus is a connected component of transcripts
sharing at least one fragment.

| Column | Type | Description |
|--------|------|-------------|
| `locus_id` | int | Locus identifier |
| `chrom` | string | Chromosome |
| `start` | int | Locus start position |
| `end` | int | Locus end position |
| `n_transcripts` | int | Number of transcripts in the locus |
| `n_genes` | int | Number of genes in the locus |
| `n_fragments` | int | Total fragments assigned to this locus |
| `mrna` | float | Total mRNA count |
| `nrna` | float | Total nRNA count |
| `gdna` | float | Estimated gDNA count |

### quant_detail.tsv / quant_detail.feather

Fragment-level detail counts broken down by splice type and strand
orientation.

| Column | Type | Description |
|--------|------|-------------|
| `transcript_id` | string | Transcript identifier |
| `spliced_sense` | float | Annotated-spliced sense fragments |
| `spliced_antisense` | float | Annotated-spliced antisense fragments |
| `unspliced_sense` | float | Unspliced sense fragments |
| `unspliced_antisense` | float | Unspliced antisense fragments |
| `unannotated_sense` | float | Unannotated-junction sense fragments |
| `unannotated_antisense` | float | Unannotated-junction antisense fragments |

### summary.json

Run statistics and trained model parameters.

```json
{
  "rigel_version": "0.1.0",
  "library": {
    "protocol": "R1-antisense",
    "strand_specificity": 0.97,
    "read1_sense": false,
    "p_r1_sense": 0.03,
    "frag_length_mean": 198.5,
    "frag_length_median": 195,
    "frag_length_std": 52.3
  },
  "alignment": {
    "total_reads": 50000000,
    "mapped_reads": 48500000,
    "unique_reads": 42000000,
    "multimapping_reads": 6500000,
    "duplicate_reads": 2100000
  },
  "fragments": {
    "total": 24250000,
    "genic": 23800000,
    "intergenic": 450000,
    "chimeric": 38000
  },
  "quantification": {
    "n_transcripts": 150000,
    "n_genes": 50000,
    "n_loci": 45000,
    "mrna_total": 23000000.0,
    "nrna_total": 620000.0,
    "gdna_total": 180000.0,
    "mrna_fraction": 0.966,
    "nrna_fraction": 0.026,
    "gdna_fraction": 0.008
  }
}
```

Key fields:

- **`library.protocol`** — Detected library protocol: `R1-sense`
  (e.g., KAPA Stranded) or `R1-antisense` (e.g., Illumina dUTP/TruSeq).
- **`library.strand_specificity`** — Strand specificity in [0.5, 1.0].
  Values near 1.0 indicate a highly stranded library; 0.5 indicates
  unstranded.
- **`quantification.mrna_fraction`** — Fraction of genic fragments
  attributed to mature RNA, providing a global measure of library
  composition.

### config.json

The complete parameter configuration used for the run, including
defaults, YAML overrides, and CLI overrides. Useful for reproducibility.

### Annotated BAM

When `--annotated-bam` is specified, Rigel writes a copy of the input
BAM with per-fragment assignment tags. This requires a second pass
through the input BAM.

| Tag | Type | Description |
|-----|------|-------------|
| ZT | int | Assigned transcript index (-1 if unassigned) |
| ZG | int | Assigned gene index |
| ZP | float | Posterior probability of the assignment |
| ZW | float | Fragment coverage weight |
| ZC | int | Count column index (splice type × strand) |
| ZH | string | Hash of all candidate transcripts |
| ZN | int | Number of candidate transcripts (ambiguity level) |

---

## Supported aligners

Rigel has been tested with the following RNA-seq aligners:

| Aligner | SJ strand tag | Notes |
|---------|---------------|-------|
| **STAR** | `XS` | Recommended. Produces NH, XS, nM/NM tags. |
| **HISAT2** | `XS` | Well supported. |
| **minimap2** | `ts` | Use `--sj-strand-tag ts` or rely on `auto` detection. |

The `--sj-strand-tag auto` default detects the appropriate tag from the
BAM header. If your aligner uses a non-standard tag, specify it
explicitly.

---

## Recipes and examples

### Basic quantification

```bash
# Step 1: Build index (once per genome/annotation)
rigel index \
    --fasta ~/ref/GRCh38.fa \
    --gtf ~/ref/gencode.v46.gtf \
    -o ~/indices/gencode_v46

# Step 2: Quantify
rigel quant \
    --bam sample.namesorted.bam \
    --index ~/indices/gencode_v46 \
    -o results/sample1/
```

### Batch processing

```bash
INDEX=~/indices/gencode_v46

for BAM in data/*.namesorted.bam; do
    SAMPLE=$(basename "$BAM" .namesorted.bam)
    rigel quant --bam "$BAM" --index "$INDEX" -o "results/$SAMPLE/" \
        --threads 8 --seed 42
done
```

### Reproducible runs

For fully deterministic results, set both `--seed` and `--threads 1`:

```bash
rigel quant --bam sample.bam --index idx/ -o out/ \
    --seed 42 --threads 1
```

Multi-threaded runs with a fixed seed produce consistent results in
most cases, but edge cases in multi-threaded BAM parsing may cause
minor non-determinism in fragment ordering.

### Unambiguous-only quantification

Skip the EM entirely and count only uniquely mapped fragments:

```bash
rigel quant --bam sample.bam --index idx/ -o out/ \
    --em-iterations 0
```

### Writing an annotated BAM

```bash
rigel quant --bam sample.bam --index idx/ -o out/ \
    --annotated-bam out/annotated.bam
```

The annotated BAM can be loaded in IGV or processed with samtools to
inspect per-fragment assignments.

### Using a YAML config for a project

```bash
# project_config.yaml
cat > project_config.yaml << 'EOF'
em_prior_alpha: 0.01
em_prior_gamma: 1.0
confidence_threshold: 0.95
prune_threshold: 0.1
threads: 16
seed: 42
include_multimap: true
keep_duplicates: false
EOF

rigel quant --bam sample.bam --index idx/ -o out/ \
    --config project_config.yaml
```

---

## FAQ

### What BAM sort order does Rigel require?

Name-sorted or collated (not coordinate-sorted). Use
`samtools sort -n -o output.bam input.bam` to name-sort.

### Does Rigel support single-end reads?

Single-end reads should work but are less thoroughly tested than
paired-end. Fragment length estimation will use the alignment length
rather than the insert size.

### How does Rigel handle multimappers?

By default, multimapping reads (NH > 1) are included in the EM. Each
alignment is treated as a candidate, and the EM assigns posterior
probabilities. Use `--no-include-multimap` to exclude them.

### What does strand specificity mean?

Strand specificity measures how well the library preparation preserved
strand information. A value of 1.0 means perfectly stranded; 0.5 means
unstranded. Rigel automatically detects this from spliced reads and
uses it to weight strand-based gDNA and nRNA estimates. Libraries with
specificity above ~0.7 provide strong strand signal; below ~0.6, Rigel
falls back to density-based estimation.

### How much memory does Rigel use?

Memory usage is dominated by the fragment buffer. The default maximum
is 2 GB. If the buffer fills, fragments are spilled to disk in the
directory specified by `--tmpdir`. Typical RNA-seq experiments
(50–100M reads) fit comfortably in the default buffer.

### Can I run Rigel on a cluster?

Yes. Each `rigel quant` invocation is independent. Use your cluster's
job scheduler to run one sample per job. The `--threads` flag controls
per-job parallelism.

### How does Rigel compare to Salmon or kallisto?

Rigel differs from alignment-free quantifiers in several ways:
1. It works directly from BAM files rather than raw FASTQ.
2. It jointly models mRNA, nRNA, and gDNA rather than only mRNA.
3. It uses hierarchical Bayesian priors informed by strand and density
   signals.
4. It provides per-fragment assignment information via the annotated
   BAM.

The trade-off is that Rigel requires aligned reads and typically runs
slower than alignment-free methods.

### Why are there separate mRNA and nRNA columns?

Total RNA-seq libraries capture both spliced (mature) and unspliced
(nascent) transcripts. Rigel's linked kinetic model deconvolves these
two species, allowing you to study transcription dynamics. If you only
need standard expression levels, use the mRNA column.
