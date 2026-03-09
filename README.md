<p align="center">
  <img src="docs/images/rigel_banner.png" alt="Rigel" width="100%"/>
</p>

<p align="center">
  <strong>Bayesian RNA-seq transcript quantification with joint mRNA, nascent RNA, and genomic DNA deconvolution</strong>
</p>

<p align="center">
  <a href="https://github.com/mkiyer/rigel/actions"><img src="https://img.shields.io/github/actions/workflow/status/mkiyer/rigel/ci.yml?branch=main&label=CI" alt="CI"></a>
  <a href="https://pypi.org/project/rigel-rnaseq/"><img src="https://img.shields.io/pypi/v/rigel-rnaseq" alt="PyPI"></a>
  <a href="https://anaconda.org/bioconda/rigel"><img src="https://img.shields.io/conda/vn/bioconda/rigel" alt="Bioconda"></a>
  <a href="https://pypi.org/project/rigel-rnaseq/"><img src="https://img.shields.io/pypi/pyversions/rigel-rnaseq" alt="Python"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License"></a>
</p>

---

## Overview

<p align="center">
  <img src="docs/images/rigel_overview.png" alt="Rigel Overview" width="85%"/>
</p>

**Rigel** is a Bayesian transcript quantification tool for RNA-seq data
that jointly estimates three RNA species present in a single library:

- **Rigel A — Mature RNA (mRNA):** Spliced, fully processed transcripts — the
  primary signal in most RNA-seq experiments.
- **Rigel B — Genomic DNA (gDNA):** Background contamination from genomic DNA,
  modeled on both forward and reverse strands.
- **Rigel C — Nascent RNA (nRNA):** Unspliced pre-mRNA transcripts captured
  mid-transcription, present in total RNA-seq and nascent RNA protocols.

Rigel employs a **linked mRNA–nRNA kinetic model** where mature and nascent
RNA share a single abundance parameter per transcript, coupled through a
per-transcript nascent fraction grounded in the steady-state kinetics of
transcription, splicing, and degradation. A **hybrid empirical Bayes framework**
with hierarchical shrinkage initializes gDNA and nRNA priors from strand
specificity and genomic density signals, enabling robust quantification even
in complex loci with overlapping genes.

### Key features

- **Joint three-species deconvolution** — Simultaneously separates mRNA, nRNA,
  and gDNA within a unified locus-level Expectation-Maximization (EM) framework
- **Linked kinetic model** — Couples mRNA and nRNA through biologically
  motivated steady-state parameters, reducing the effective parameter space
- **Automatic protocol detection** — Learns library strandedness (R1-sense or
  R1-antisense) and strand specificity from spliced reads with no user
  configuration
- **Hierarchical Bayesian priors** — Three-level empirical Bayes shrinkage for
  nRNA fractions (transcript → TSS group → locus → global) and two-level
  shrinkage for gDNA rates (locus → chromosome → global)
- **Coverage-weighted OVR prior** — One Virtual Read prior distributes prior
  mass according to geometric fragment evidence rather than uniform assumptions
- **High-performance C++ core** — Single-pass BAM scanning with htslib, native
  EM solver with SQUAREM acceleration, and memory-efficient columnar buffering
- **Strand-aware scoring** — Per-fragment log-likelihood incorporates strand
  probability, insert size distribution, alignment penalties, and positional
  coverage weights
- **Comprehensive output** — Per-transcript, per-gene, and per-locus abundance
  tables in both Feather and TSV formats, with optional annotated BAM output

---

## Installation

### From Bioconda (recommended)

```bash
conda install -c conda-forge -c bioconda rigel
```

### From PyPI

```bash
pip install rigel-rnaseq
```

Requires a C++17-capable compiler (GCC 7+, Clang 5+) and htslib headers.
On macOS, install Xcode Command Line Tools first (`xcode-select --install`).

### From source

```bash
# Clone the repository
git clone https://github.com/mkiyer/rigel.git
cd rigel

# Create conda environment with all dependencies
mamba env create -f mamba_env.yaml
conda activate rigel

# Install in development mode
pip install --no-build-isolation -e .
```

### Requirements

- **Python** ≥ 3.12
- **C++17 compiler** (GCC 7+, Clang 5+, or MSVC 19.14+)
- **htslib** (provided via pysam)
- **Runtime dependencies:** numpy, pandas, pyarrow, pysam, pyyaml

---

## Quick start

### 1. Build a reference index

```bash
rigel index \
    --fasta genome.fa \
    --gtf annotation.gtf \
    -o my_index/
```

The index step parses the GTF annotation and genome FASTA to build an
interval-based reference structure used for fragment resolution. It
requires a FASTA index (`.fai`); run `samtools faidx genome.fa` first if
one does not exist.

### 2. Quantify a BAM file

```bash
rigel quant \
    --bam aligned.bam \
    --index my_index/ \
    -o results/
```

The input BAM must be **name-sorted** (or collated) and should contain the
`NH` tag for multimapper detection (produced by STAR, HISAT2, minimap2, etc.).

### 3. Examine results

```bash
# Per-transcript abundances
head results/quant.tsv

# Per-gene aggregates
head results/gene_quant.tsv

# Run summary (protocol, strand specificity, fragment counts)
cat results/summary.json
```

---

## Output files

| File | Description |
|------|-------------|
| `quant.feather` / `quant.tsv` | Per-transcript abundances: transcript ID, gene ID, mRNA count, nRNA count, effective length, TPM |
| `gene_quant.feather` / `gene_quant.tsv` | Per-gene aggregated abundances |
| `loci.feather` / `loci.tsv` | Per-locus summary with component counts and gDNA estimates |
| `quant_detail.feather` / `quant_detail.tsv` | Fragment-level assignment counts by splice type and strand |
| `summary.json` | Run statistics: library protocol, strand specificity, fragment length distribution, alignment counts, quantification totals |
| `config.json` | Full parameter configuration used for the run |
| `annotated.bam` (optional) | Input BAM with per-fragment assignment tags (ZT, ZG, ZP, ZW, ZC, ZH, ZN) |

---

## How it works

Rigel implements a **two-stage single-pass pipeline**:

**Stage 1 — BAM Scan & Model Training.** A C++ BAM scanner reads the
name-sorted input once, resolving each fragment against the reference index
to determine exonic, intronic, and intergenic overlap. Spliced reads with
known junction strand train a **strand model** (Beta posterior on R1-sense
probability) and a **fragment length distribution** (separate histograms for
RNA and gDNA). Fragments are buffered in memory-efficient columnar format
with optional disk spill.

**Stage 2 — Locus-Level EM Quantification.** Fragments are routed into
compressed sparse row (CSR) format and partitioned into independent loci
via connected-component analysis. For each locus with *T* transcripts, the
EM solves a **2T + 1 component** mixture model (T mRNA + T nRNA + 1 gDNA).
Transcript abundances θ and nascent fractions β are updated alternately:
θ on the simplex via MAP-EM with Dirichlet + OVR priors, β via Beta
posterior MAP estimates. SQUAREM acceleration typically achieves convergence
in 10–50× fewer iterations than standard EM. Post-EM pruning removes
components with insufficient evidence.

For complete mathematical details, see the [Methods](docs/METHODS.md)
document. For all parameters and usage options, see the
[Manual](docs/MANUAL.md).

---

## Documentation

| Document | Description |
|----------|-------------|
| [Manual](docs/MANUAL.md) | Complete usage guide: commands, parameters, output formats, configuration |
| [Methods](docs/METHODS.md) | Detailed algorithmic description suitable for citation in publications |
| [Parameters](docs/parameters.md) | Quick reference for all tunable parameters and defaults |

---

## Citing Rigel

If you use Rigel in your research, please cite:

> Iyer MK. **Rigel: Bayesian RNA-seq transcript quantification with joint
> mRNA, nascent RNA, and genomic DNA deconvolution.** (2026).
> https://github.com/mkiyer/rigel

---

## License

Rigel is distributed under the [GNU General Public License v3.0](LICENSE).

---

## Contributing

Contributions are welcome. Please open an issue to discuss proposed changes
before submitting a pull request.

```bash
# Run the test suite
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=rigel --cov-report=term-missing
```

