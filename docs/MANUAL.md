# Rigel User Manual

Complete usage reference for Rigel 0.2.0.

---

## Overview

Rigel quantifies RNA-seq libraries while jointly modeling:

- mRNA from annotated transcripts
- nRNA from shared genomic spans
- gDNA contamination at the locus level

The quantifier runs in two stages:

1. A native BAM scan resolves fragments, trains the strand and
   fragment-length models, and buffers fragment metadata.
2. A locus-level EM solver assigns ambiguous signal across mRNA, nRNA,
   and gDNA components.

In the current implementation, nRNA is represented by a global table of unique
genomic spans `(ref, strand, start, end)`. Transcripts that share the same span
share the same nRNA component.

---

## Installation

### Bioconda

```bash
conda install -c conda-forge -c bioconda rigel
```

### PyPI

```bash
pip install rigel-rnaseq
```

### From source

```bash
git clone https://github.com/mkiyer/rigel.git
cd rigel

mamba env create -f mamba_env.yaml
conda activate rigel
pip install --no-build-isolation -e .
```

### Requirements

- Python 3.12+
- C++17-capable compiler for source builds
- `numpy`, `pandas`, `pyarrow`, `pysam`, `pyyaml`

On macOS:

```bash
xcode-select --install
```

Verify the install:

```bash
rigel --version
```

---

## Input requirements

### Reference inputs

- Genome FASTA with `.fai` index
- Gene annotation GTF with `gene`, `transcript`, and `exon` features

### BAM requirements

- Name-sorted or collated, not coordinate-sorted
- `NH` tag present for multimapper handling
- Splice-junction strand tag available when possible
- Paired-end is the main supported mode; single-end is less exercised

Typical resort commands:

```bash
samtools sort -n -o sample.name_sorted.bam sample.bam
samtools collate -o sample.collated.bam sample.bam
```

---

## Commands

### rigel index

Build the reference index used for fragment resolution.

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
```

Arguments:

| Flag | Required | Description |
|------|----------|-------------|
| `--fasta` | yes | Genome FASTA, indexed with `samtools faidx` |
| `--gtf` | yes | Annotation GTF |
| `-o`, `--output-dir` | yes | Output directory |
| `--feather-compression` | no | `lz4`, `zstd`, or `uncompressed`; default `lz4` |
| `--no-tsv` | no | Skip TSV mirrors |
| `--gtf-parse-mode` | no | `strict` or `warn-skip`; default `strict` |

Index outputs:

- `ref_lengths.feather` and optional `ref_lengths.tsv`
- `transcripts.feather` and optional `transcripts.tsv`
- `nrna.feather` and optional `nrna.tsv`
- `sj.feather` and optional `sj.tsv`
- `intervals.feather` and optional `intervals.tsv`

### rigel quant

Run the full scan-plus-EM quantification pipeline.

```bash
rigel quant --bam sample.bam --index index/ -o results/
```

Common example:

```bash
rigel quant \
    --bam sample.bam \
    --index index/ \
    -o results/ \
    --config params.yaml \
    --threads 16 \
    --annotated-bam results/annotated.bam \
    --seed 42
```

### rigel sim

Generate synthetic scenarios for benchmarking and regression testing.

```bash
rigel sim --config scenario.yaml -o sim_out/
```

Arguments:

| Flag | Required | Description |
|------|----------|-------------|
| `--config` | yes | Scenario YAML |
| `-o`, `--output-dir` | yes | Output directory |
| `--genome-length` | no | Default `5000`; YAML can override |
| `--seed` | no | Default `42`; YAML can override |
| `--num-reads` | no | Default `1000`; YAML can override |

---

## Quant parameters

Resolution order is:

1. Explicit CLI flag
2. YAML config file via `--config`
3. Built-in default

Unknown YAML keys are ignored with a warning. YAML keys may use either
underscores or hyphens.

### Required arguments

| Flag | Description |
|------|-------------|
| `--bam` | Name-sorted or collated BAM |
| `--index` | Rigel index directory |
| `-o`, `--output-dir` | Output directory |

### Read filtering and routing

| Flag | Default | Description |
|------|---------|-------------|
| `--include-multimap` / `--no-include-multimap` | yes | Include multimappers detected by `NH` or secondary flags |
| `--keep-duplicates` / `--no-keep-duplicates` | no | Keep duplicate-marked reads |
| `--sj-strand-tag` | `auto` | Auto-detect or specify one or more tags such as `XS` or `ts` |
| `--seed` | current timestamp | Used for reproducibility when set |

### Core EM settings

| Flag | Default | Description |
|------|---------|-------------|
| `--em-iterations` | `1000` | Max EM iterations; `0` means unambiguous-only quantification |
| `--em-prior-alpha` | `0.01` | Flat Dirichlet pseudocount per eligible component |
| `--em-prior-gamma` | `1.0` | OVR prior scale factor |
| `--em-convergence-delta` | `1e-6` | Convergence threshold |
| `--prune-threshold` | `0.1` | Post-EM evidence-ratio threshold; set negative to disable |
| `--confidence-threshold` | `0.95` | RNA-normalized posterior threshold for high-confidence assignment |
| `--em-mode` | `vbem` | EM algorithm variant: `map` (MAP-EM) or `vbem` (Variational Bayes EM) |
| `--nrna-sparsity-alpha` | `0.9` | Dirichlet α for nRNA components. Values < 1.0 sparsify weak nRNA toward zero |
| `--gdna-prior-scale` | `1.0` | Scale factor for EB gDNA prior anchor: α = 1 + scale × gdna\_init |

### Fragment scoring

| Flag | Default | Description |
|------|---------|-------------|
| `--overhang-alpha` | `0.01` | Per-base overhang penalty in probability space |
| `--mismatch-alpha` | `0.1` | Per-mismatch penalty from the `NM` tag |
| `--gdna-splice-penalty-unannot` | `0.01` | gDNA penalty for unannotated spliced fragments |

### gDNA prior settings

The `--gdna-prior-scale` parameter (in Core EM settings above) controls the
strength of the Empirical Bayes gDNA anchor in the per-locus prior.

| Flag | Default | Description |
|------|---------|-------------|
| `--gdna-kappa-ref` | auto | Shrinkage from reference toward global |
| `--gdna-kappa-locus` | auto | Shrinkage from locus toward reference |
| `--gdna-mom-min-evidence-ref` | `50.0` | Minimum evidence for reference MoM estimation |
| `--gdna-mom-min-evidence-locus` | `30.0` | Minimum evidence for locus MoM estimation |
| `--gdna-kappa-min` | `2.0` | Lower clamp for MoM-estimated kappa |
| `--gdna-kappa-max` | `200.0` | Upper clamp for MoM-estimated kappa |
| `--gdna-kappa-fallback` | `5.0` | Fallback kappa when evidence is insufficient |
| `--gdna-kappa-min-obs` | `20` | Minimum features needed for MoM kappa estimation |

### Strand model

| Flag | Default | Description |
|------|---------|-------------|
| `--strand-prior-kappa` | `2.0` | Beta prior pseudo-count for the RNA strand model |

### Output and performance

| Flag | Default | Description |
|------|---------|-------------|
| `--no-tsv` | off | Write Feather only |
| `--annotated-bam` | unset | Write a second-pass BAM with per-fragment assignment tags |
| `--tmpdir` | system temp | Spill directory for fragment-buffer overflow |
| `--threads` | `0` | `0` means all available cores; used for both scan and locus EM |

Rigel sets `OMP_NUM_THREADS=1` on import via `os.environ.setdefault(...)` to
avoid idle OpenMP thread pools competing with Rigel's native parallelism. If a
different OpenMP setting is required in the same process, set it before
importing Rigel.

---

## YAML configuration

Example:

```yaml
em_prior_alpha: 0.02
em_prior_gamma: 1.0
em_iterations: 1000
confidence_threshold: 0.95
prune_threshold: 0.1

include_multimap: true
keep_duplicates: false
sj_strand_tag: [XS, ts]

em_mode: vbem
nrna_sparsity_alpha: 0.9
gdna_prior_scale: 1.0

gdna_kappa_ref: null
gdna_kappa_locus: null

threads: 16
seed: 42
```

Usage:

```bash
rigel quant --bam sample.bam --index index/ -o out/ --config params.yaml
```

CLI still wins over YAML:

```bash
rigel quant \
    --bam sample.bam \
    --index index/ \
    -o out/ \
    --config params.yaml \
    --threads 8
```

---

## Output files

All quant outputs are written under `--output-dir`.

### quant.feather and quant.tsv

Transcript-level abundance table.

| Column | Description |
|--------|-------------|
| `transcript_id` | Transcript ID from the index |
| `gene_id` | Parent gene ID |
| `gene_name` | Gene symbol or name field from the GTF |
| `locus_id` | Locus identifier, or `-1` if no EM locus was formed |
| `effective_length` | Effective transcript length used for TPM |
| `mrna` | Total mRNA count |
| `mrna_unambig` | Deterministic mRNA count |
| `mrna_em` | EM-assigned mRNA count |
| `mrna_high_conf` | Unambiguous plus high-confidence EM assignments |
| `mrna_spliced` | Spliced mRNA count |
| `nrna` | nRNA count fanned back from shared nRNA spans |
| `rna_total` | `mrna + nrna` |
| `tpm` | mRNA-based TPM |
| `posterior_mean` | Mean posterior of EM assignments touching the transcript |

### gene_quant.feather and gene_quant.tsv

Gene-level aggregation.

| Column | Description |
|--------|-------------|
| `gene_id` | Gene ID |
| `gene_name` | Gene name |
| `locus_id` | Primary locus assigned to the gene, or `-1` |
| `effective_length` | Abundance-weighted mean transcript effective length |
| `mrna` | Total mRNA count |
| `mrna_unambig` | Deterministic mRNA count |
| `mrna_em` | EM-assigned mRNA count |
| `mrna_high_conf` | High-confidence mRNA count |
| `mrna_spliced` | Spliced mRNA count |
| `nrna` | Gene-level nRNA after transcript fan-out |
| `rna_total` | `mrna + nrna` |
| `tpm` | mRNA-based TPM |

### loci.feather and loci.tsv

Per-locus EM summary.

| Column | Description |
|--------|-------------|
| `locus_id` | Locus ID |
| `ref` | Primary reference represented in the locus |
| `n_transcripts` | Number of transcripts in the locus |
| `n_genes` | Number of genes in the locus |
| `n_em_fragments` | Number of ambiguous units entering the locus EM |
| `mrna` | Locus mRNA total |
| `nrna` | Locus nRNA total |
| `gdna` | Locus gDNA total |
| `gdna_rate` | `gdna / (mrna + nrna + gdna)` |
| `gdna_init` | Empirical Bayes gDNA initialization before EM |

### quant_detail.feather and quant_detail.tsv

Long-format QC table.

| Column | Description |
|--------|-------------|
| `transcript_id` | Transcript ID |
| `gene_id` | Gene ID |
| `category` | One of `unspliced`, `spliced_unannot`, `spliced_annot` |
| `source` | Either `unambig` or `em` |
| `count` | Count in that category/source cell |

### summary.json

`summary.json` contains:

- `rigel_version`, `command`, `timestamp`
- `input` with absolute BAM and index paths
- `library` with `protocol`, `strand_specificity`, `p_r1_sense`, `read1_sense`, and fragment-length summary stats
- `alignment` with BAM-level counters
- `fragments` with total, genic, intergenic, and chimeric counts
- `quantification` with transcript/gene/locus counts and global mRNA, nRNA, and gDNA totals and fractions
- `strand_models`, `frag_length_models`, and `pipeline_stats`

The `protocol` field is reported as `R1-sense` when `p_r1_sense >= 0.5`, else
`R1-antisense`.

### config.json

`config.json` records the resolved runtime configuration after merging:

- command metadata
- optional source config path
- absolute input and output paths
- every quant parameter after YAML and CLI resolution

### Annotated BAM

When `--annotated-bam` is set, Rigel performs a second BAM pass and writes tags
to every record.

| Tag | Type | Description |
|-----|------|-------------|
| `ZT` | string | Transcript ID, or `.` for intergenic or gDNA-only assignments |
| `ZG` | string | Gene ID, or `.` when not applicable |
| `ZP` | string | Assignment pool: `mRNA`, `nRNA`, `gDNA`, `intergenic`, or `chimeric` |
| `ZW` | float | Posterior probability of the chosen assignment |
| `ZC` | string | Fragment class: `unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`, `chimeric`, or `intergenic` |
| `ZH` | int | Primary-hit flag: `1` for the winning alignment, `0` otherwise |
| `ZN` | int | Number of competing candidate components |
| `ZS` | string | Splice type: `spliced_annot`, `spliced_unannot`, `unspliced`, or `ambiguous` |

---

## Supported aligners

Known splice-junction strand tags:

| Aligner | Tag | Notes |
|---------|-----|-------|
| STAR | `XS` | Recommended |
| HISAT2 | `XS` | Standard RNA-seq workflow |
| minimap2 | `ts` | Long-read and splice-aware mappings |

Use `--sj-strand-tag auto` unless automatic detection fails.

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

### Why does PyPI use `rigel-rnaseq`?

The name `rigel` is already occupied on PyPI. The CLI, import name, GitHub
repository, and Bioconda package are still `rigel`.

### Why can transcript `nrna` be non-integer or shared?

Because the EM estimates nRNA on shared genomic spans. Those shared counts
are fanned back to transcripts for transcript- and gene-level reporting.

### Why is `strand_specificity` close to `0.5`?

Rigel only trains the primary strand model from annotated spliced fragments. If
the library is unstranded, poorly stranded, or has too few informative splice
reads, the estimate will stay near the prior.

### When should I use `--annotated-bam`?

Use it for read-level inspection, debugging, and method validation. It requires
an extra BAM pass and adds some runtime overhead.

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
