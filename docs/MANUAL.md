# Rigel User Manual

Rigel quantifies RNA-seq alignments while jointly modeling mature mRNA,
nascent RNA (nRNA), and genomic DNA contamination (gDNA). It takes aligned
BAM files as input and produces per-transcript, per-gene, and per-locus
abundance estimates.

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

Verify the install:

```bash
rigel --version
```

---

## Input requirements

### Reference

- Genome FASTA with `.fai` index (`samtools faidx genome.fa`)
- Gene annotation GTF (GENCODE format recommended; `exon` records required)

### BAM

- **Name-sorted or collated** — not coordinate-sorted
- **`NH` tag present** for multimapper handling
- Splice-junction strand tag recommended for strand-model accuracy
  (`XS` from STAR/HISAT2, `ts` from minimap2, or `auto` for detection)

Name-sort an existing BAM:

```bash
samtools sort -n -@ 8 -o sample.namesorted.bam sample.bam
```

---

## Quick start

```bash
# Step 1: Build index (once per genome + annotation)
rigel index \
    --fasta genome.fa \
    --gtf annotation.gtf \
    -o index/

# Step 2: Quantify
rigel quant \
    --bam sample.namesorted.bam \
    --index index/ \
    -o results/sample/ \
    --threads 8 \
    --seed 42

# Step 3: Inspect outputs
head results/sample/quant.tsv
cat results/sample/summary.json
```

---

## Commands

### rigel index

Builds the reference index from a genome FASTA and GTF. Run once per
genome/annotation combination; the index is shared across all samples.

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
```

| Flag | Default | Description |
|------|---------|-------------|
| `--fasta` | required | Genome FASTA (must have `.fai` index) |
| `--gtf` | required | Annotation GTF |
| `-o`, `--output-dir` | required | Output directory |
| `--nrna-tolerance` | `20` | Max distance (bp) for clustering transcript start/end sites into shared nRNA spans |
| `--gtf-parse-mode` | `strict` | `strict` fails on malformed GTF records; `warn-skip` logs warnings and skips them |
| `--feather-compression` | `lz4` | Feather compression: `lz4`, `zstd`, or `uncompressed` |
| `--no-tsv` | off | Skip writing TSV mirrors of index files |

### rigel quant

Runs the full scan-and-EM quantification pipeline.

```bash
rigel quant \
    --bam sample.namesorted.bam \
    --index index/ \
    -o results/ \
    --threads 8 \
    --seed 42
```

All flags are described in [Parameters Reference](parameters.md).

**Common options for pipeline use:**

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | 0 (all cores) | Threads for BAM scan and locus EM |
| `--seed N` | timestamp | Set for reproducibility |
| `--sj-strand-tag TAG` | `auto` | Use `XS` for STAR/HISAT2, `ts` for minimap2 |
| `--tmpdir DIR` | system temp | Spill directory for buffer overflow; use local SSD |
| `--tsv` | off | Write TSV mirrors alongside Feather outputs |
| `--config FILE` | — | YAML config file (see [YAML configuration](#yaml-configuration)) |
| `--annotated-bam PATH` | — | Write annotated BAM with per-fragment assignment tags |

### rigel export

Converts all `.feather` output files in a results directory to TSV or Parquet.
Useful as a post-processing step or when Parquet is preferred for downstream tools.

```bash
# Convert to TSV
rigel export results/sample/ --format tsv

# Convert to Parquet
rigel export results/sample/ --format parquet
```

| Argument | Default | Description |
|----------|---------|-------------|
| `output_dir` | required | Directory containing `.feather` files |
| `-f`, `--format` | `tsv` | Output format: `tsv` or `parquet` |

### rigel sim

Generates synthetic test scenarios for benchmarking and development.

```bash
rigel sim --config scenario.yaml -o sim_out/
```

---

## YAML configuration

Any `rigel quant` flag can be set in a YAML file via `--config`. Keys use
underscores (hyphens also accepted). Explicit CLI flags always override the
YAML. Unknown keys are ignored with a warning.

```yaml
# rigel_params.yaml

# Library handling
include_multimap: true
keep_duplicates: false
sj_strand_tag: [auto]       # [XS] for STAR, [ts] for minimap2, [XS, ts] to try both

# EM algorithm (defaults are suitable for most libraries)
prior_pseudocount: 1.0
em_iterations: 1000
assignment_mode: sample     # sample | fractional | map
em_mode: vbem               # vbem | map

# Performance
threads: 16
seed: 42
```

```bash
rigel quant \
    --bam sample.bam \
    --index index/ \
    -o results/ \
    --config rigel_params.yaml
```

CLI override still works:

```bash
rigel quant --bam sample.bam --index index/ -o results/ \
    --config rigel_params.yaml \
    --threads 4
```

---

## Output files

All outputs are written to `--output-dir`. Feather files (`.feather`) use
ZSTD compression and are readable with `pandas.read_feather()` or
`pyarrow.feather.read_feather()`.

A `config.yaml` is also written to the output directory, recording all
resolved parameters and I/O paths. Rerun the exact same analysis with:

```bash
rigel quant --config results/config.yaml
```

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
| `effective_length` | Length after fragment-length correction |
| `locus_id` | Locus ID; `-1` if not placed in an EM locus |
| `nrna_id` | ID of the parent nRNA entity for this transcript (`"."` if none or if the transcript is itself an nRNA entity) |
| `is_basic` | `1` if transcript has the GENCODE `basic` tag |
| `is_mane` | `1` if transcript is MANE Select or MANE Plus Clinical |
| `is_nrna` | `1` if transcript is single-exon (indistinguishable from nascent RNA) |
| `count` | Fragment count assigned to this transcript |
| `count_unambig` | Deterministically assigned fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Fragment count from annotated-spliced fragments |
| `nrna_parent_count` | Parent nRNA entity's count. Zero if `is_nrna=1` (progenitors have no parent) or if no parent nRNA exists |
| `tpm` | Transcripts per million (annotated transcripts only) |
| `tpm_total_rna` | TPM normalized over all RNA (annotated + synthetic nRNA). Comparable to `tpm` in `nrna_quant` |
| `posterior_mean` | Mean EM posterior over units assigned to this transcript |

### gene_quant.feather / gene_quant.tsv

Gene-level abundance aggregated from annotated (non-synthetic) transcript
estimates. No nRNA attribution is performed at the gene level.

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
| `count_unambig` | Deterministic fragment count |
| `count_em` | EM-assigned fragment count |
| `count_spliced` | Spliced fragment count |
| `tpm` | Gene-level TPM |

### nrna_quant.feather / nrna_quant.tsv

nRNA entity estimates. Each row is one nRNA entity (a single-exon
transcript, either annotated or synthetic) that has at least one
multi-exon mRNA child transcript. Standalone single-exon genes with
no multi-exon overlap are excluded; they appear in `quant.feather`
with `is_nrna=True`.

| Column | Description |
|--------|-------------|
| `nrna_id` | nRNA entity transcript ID |
| `gene_id` | Overlapping gene ID |
| `gene_name` | Gene name |
| `ref`, `strand`, `start`, `end` | Entity genomic coordinates |
| `length` | Entity length (bp) |
| `effective_length` | Effective length after fragment-length correction |
| `locus_id` | Locus containing this entity |
| `is_synthetic` | `1` if this entity was synthesized by rigel |
| `n_contributing_transcripts` | Number of multi-exon transcripts merged into this span |
| `n_mrna` | Number of multi-exon mRNA transcripts associated with this nRNA entity |
| `mrna_count` | Sum of mRNA children's fragment counts |
| `count` | nRNA entity's own fragment count |
| `tpm` | Total RNA TPM (denominator includes all annotated + synthetic transcripts). Directly comparable to `tpm_total_rna` in `quant.feather` |

### loci.feather / loci.tsv

Per-locus EM summary. One row per connected component of overlapping
transcripts.

| Column | Description |
|--------|-------------|
| `locus_id` | Locus ID |
| `ref` | Primary chromosome |
| `n_transcripts` | Total transcript count (mRNA + nRNA spans) |
| `n_annotated_transcripts` | Annotated (non-synthetic) transcript count |
| `n_nrna_entities` | Number of nRNA entities (all single-exon transcripts) |
| `n_genes` | Number of genes |
| `n_em_fragments` | Ambiguous fragments entering EM |
| `mrna` | Total mRNA count |
| `nrna` | Total nRNA count |
| `gdna` | Total gDNA count |
| `total` | `mrna + nrna + gdna` |
| `gdna_rate` | `gdna / total` |
| `gdna_init` | Calibrated per-locus gDNA prior (γ) used for EM initialization |

### summary.json

Run-level summary. Key sections:

| Section | Contents |
|---------|----------|
| `rigel_version`, `timestamp` | Version and run time |
| `command` | Subcommand, arguments, config file path |
| `configuration` | All resolved parameters after CLI + YAML merge |
| `input` | Absolute BAM and index paths |
| `alignment_stats` | Total, mapped, unique, multimapping, duplicate, QC-fail read counts |
| `fragment_stats` | Genic, intergenic, chimeric counts |
| `strand_model` | Protocol (`R1-sense` / `R1-antisense`), strand specificity, 95% CI |
| `calibration` | gDNA calibration summary: density, κ, π, per-ref stats |
| `fragment_length` | Per-category summary statistics (mean/std/median/mode/n_obs) and full histograms with trimmed zero bins for QC plotting |
| `quantification` | n_transcripts, n_genes, n_loci, mRNA/nRNA/gDNA totals and fractions |

### Annotated BAM

Produced by `--annotated-bam PATH`. Requires a second pass over the BAM.

| Tag | Type | Description |
|-----|------|-------------|
| `ZT` | string | Assigned transcript ID, or `.` |
| `ZG` | string | Assigned gene ID, or `.` |
| `ZR` | string | Assigned gene name / symbol, or `.` |
| `ZI` | int | Transcript index into rigel reference (`-1` if unassigned) |
| `ZJ` | int | Gene index into rigel reference (`-1` if unassigned) |
| `ZF` | int | Assignment flags bitfield (see below) |
| `ZW` | float | Posterior probability of the assignment |
| `ZC` | string | Fragment class: `unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`, `chimeric`, or `intergenic` |
| `ZH` | int | Primary-hit flag: `1` for the winning alignment, `0` otherwise. Note: this reflects rigel's EM assignment, which may differ from the aligner's primary/secondary FLAG. |
| `ZN` | int | Number of competing candidate components |
| `ZS` | string | Splice type: `spliced_annot`, `spliced_unannot`, `unspliced`, or `unknown` |
| `ZL` | int | Locus ID (`-1` if no locus) |

#### ZF assignment flags

The `ZF` tag is an integer bitfield encoding both the EM assignment result
and properties of the assigned transcript:

| Bit | Mask | Flag | Meaning |
|-----|------|------|---------|
| 0 | 0x1 | `is_resolved` | Fragment was scored and assigned by the EM |
| 1 | 0x2 | `is_gdna` | Assigned to the gDNA EM component |
| 2 | 0x4 | `is_nrna` | Assigned transcript is single-exon (nRNA candidate) |
| 3 | 0x8 | `is_synthetic` | Assigned transcript is a rigel-generated nRNA span |

Valid ZF values:

| ZF | Meaning |
|----|---------|
| 0 | Not resolved (intergenic, chimeric, or filtered) |
| 1 | Transcript assigned (multi-exon, annotated) |
| 3 | gDNA component assigned |
| 5 | Transcript assigned (single-exon nRNA, annotated) |
| 13 | Transcript assigned (single-exon nRNA, synthetic) |

Pysam usage:

```python
zf = read.get_tag("ZF")
is_resolved  = (zf & 0x1) != 0
is_gdna      = (zf & 0x2) != 0
is_nrna      = (zf & 0x4) != 0
is_synthetic = (zf & 0x8) != 0
is_transcript = is_resolved and not is_gdna
```

---

## Snakemake integration

### Workflow structure

A typical RNA-seq pipeline integrating Rigel:

1. Align reads (STAR, HISAT2, or minimap2)
2. Name-sort the BAM
3. Build the Rigel index once
4. Run `rigel quant` per sample

### Snakefile rules

```python
# Snakefile

rule all:
    input:
        expand("results/{sample}/quant.feather", sample=config["samples"]),
        expand("results/{sample}/summary.json",  sample=config["samples"]),


# Build the Rigel index — run once per genome + annotation
rule rigel_index:
    input:
        fasta = config["genome_fasta"],
        gtf   = config["genome_gtf"],
    output:
        directory(config["rigel_index"]),
    log:
        "logs/rigel_index.log",
    threads: 1
    shell:
        """
        rigel index \
            --fasta {input.fasta} \
            --gtf   {input.gtf}   \
            -o      {output}      \
            > {log} 2>&1
        """


# Name-sort BAM if not already sorted
rule namesort_bam:
    input:
        bam = "aligned/{sample}.bam",
    output:
        bam = temp("namesorted/{sample}.bam"),
    threads: 4
    shell:
        "samtools sort -n -@ {threads} -o {output.bam} {input.bam}"


# Quantify one sample
rule rigel_quant:
    input:
        bam   = "namesorted/{sample}.bam",
        index = config["rigel_index"],
    output:
        quant      = "results/{sample}/quant.feather",
        gene_quant = "results/{sample}/gene_quant.feather",
        nrna_quant = "results/{sample}/nrna_quant.feather",
        loci       = "results/{sample}/loci.feather",
        summary    = "results/{sample}/summary.json",
        run_config = "results/{sample}/config.yaml",
    log:
        "logs/rigel_{sample}.log",
    threads: 8
    params:
        outdir = "results/{sample}",
        config = config.get("rigel_config", ""),
    shell:
        """
        rigel quant \
            --bam    {input.bam}   \
            --index  {input.index} \
            -o       {params.outdir} \
            --threads {threads} \
            --seed   42 \
            $([ -n "{params.config}" ] && echo "--config {params.config}") \
            > {log} 2>&1
        """
```

### Snakemake config (config.yaml)

```yaml
# config.yaml

samples:
  - sample1
  - sample2
  - sample3

genome_fasta: /path/to/GRCh38.primary_assembly.fa
genome_gtf:   /path/to/gencode.v46.primary_assembly.annotation.gtf
rigel_index:  /path/to/rigel_index/

# Optional: path to a rigel quant YAML config shared across all samples
rigel_config: config/rigel_params.yaml
```

### Shared rigel parameters (rigel_params.yaml)

```yaml
# config/rigel_params.yaml
include_multimap: true
keep_duplicates: false
sj_strand_tag: [auto]
prior_pseudocount: 1.0
assignment_mode: sample
em_mode: vbem
```

### Loading results in Python

```python
import pandas as pd

samples = ["sample1", "sample2", "sample3"]

# Gene-level mRNA matrix
gene_dfs = [
    pd.read_feather(f"results/{s}/gene_quant.feather").assign(sample=s)
    for s in samples
]
gene_quant = pd.concat(gene_dfs, ignore_index=True)

mrna_matrix = gene_quant.pivot(index="gene_id", columns="sample", values="mrna")
tpm_matrix  = gene_quant.pivot(index="gene_id", columns="sample", values="tpm")
```

---

## Supported aligners

| Aligner | SJ strand tag | Notes |
|---------|--------------|-------|
| STAR | `XS` | Recommended; name-grouped output by default |
| HISAT2 | `XS` | Standard RNA-seq workflow |
| minimap2 | `ts` | Long-read and splice-aware mode |

Use `--sj-strand-tag auto` (default) to detect the tag automatically.

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

### Fully reproducible run

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --seed 42 --threads 1
```

A `config.yaml` is written to the output directory, so you can rerun later:

```bash
rigel quant --config results/config.yaml
```

### Output TSV tables

```bash
# During quantification
rigel quant --bam sample.bam --index index/ -o results/ --tsv

# Or convert existing Feather outputs afterward
rigel export results/ --format tsv
```

### Inspect read assignments

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --annotated-bam results/annotated.bam

# Count gDNA-assigned fragments (ZF=3)
samtools view -F 256 results/annotated.bam \
    | awk '{ for(i=12;i<=NF;i++) if($i=="ZF:i:3") count++ } END { print count }'
```

### Exclude multimappers

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --no-include-multimap
```

### Unambiguous counts only (skip EM)

```bash
rigel quant \
    --bam sample.bam --index index/ -o results/ \
    --em-iterations 0
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
Rigel trains the strand model from annotated spliced fragments only. Unstranded
libraries, or libraries with few informative splice reads, will stay near 0.5.

**How much memory does Rigel use?**
The fragment buffer defaults to 2 GiB. Overflow spills to disk under
`--tmpdir`. Typical paired-end libraries (50–100M read pairs) fit within
the default.

**Can I run Rigel on a cluster?**
Yes. Each `rigel quant` is independent. Parallelise at the sample level via
your scheduler; use `--threads` for intra-job parallelism.

**How reproducible are results?**
Set `--seed` for deterministic post-EM assignment sampling. For fully
bit-reproducible output, also set `--threads 1`.

**When should I use `--annotated-bam`?**
For read-level inspection, debugging, or assignment validation. It requires
a second BAM pass and adds some runtime overhead.

**Does Rigel support single-end reads?**
Single-end reads are handled but less thoroughly tested than paired-end.
Fragment length estimation uses alignment length rather than insert size.

**What is `VBEM_CLAMP_FLOOR` and why does it matter?**
In VBEM mode, Rigel uses SQUAREM acceleration to speed up EM convergence.
SQUAREM can overshoot, pushing a component's Dirichlet alpha to very small
values. Because VBEM E-step weights depend on `digamma(α)`, which diverges
as `−1/α` near zero, components pushed below α ≈ 0.01 enter an absorbing
regime where they can never recover — even if they have genuine read
support. `VBEM_CLAMP_FLOOR` (default 0.1) sets a minimum alpha after each
SQUAREM iteration, keeping components in the recoverable regime. This
constant is defined in `src/rigel/native/em_solver.cpp` and has no effect
on MAP-EM mode. See [parameters.md](parameters.md) for all compile-time
constants.
