<p align="center">
  <img src="docs/images/rigel_banner.png" alt="Rigel" width="100%"/>
</p>

<p align="center">
  <strong>Bayesian RNA-seq quantification with joint mRNA, nascent RNA, and genomic DNA modeling</strong>
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

Rigel quantifies RNA-seq alignments while explicitly modeling three sources of
signal in the same library:

- Mature RNA (mRNA)
- Nascent RNA (nRNA)
- Genomic DNA contamination (gDNA)

A single-pass native BAM scan feeds a calibration stage that separates RNA from gDNA
genome-wide, and a locus-level EM then assigns the RNA to transcripts. Nascent RNA is
represented as one component per unique genomic span `(ref, strand, start, end)`, shared by
every transcript with that span, so isoforms that start and end at the same coordinates do
not multiply nRNA states.

### Key features

- Joint mRNA, nRNA, and gDNA quantification in one locus-level model
- Single-pass C++ BAM scanner using htslib, with memory-bounded buffering and spill-to-disk support
- Automatic strand-model training from annotated spliced fragments; protocol auto-detection (`R1-sense` / `R1-antisense`)
- gDNA/RNA calibration that solves each genomic region's unspliced mass into sense RNA, antisense RNA and gDNA, with a gDNA-density prior fitted from the library itself
- Calibrated per-locus gDNA priors feeding the EM
- MAP-EM and Variational Bayes EM (VBEM, default) solver modes with SQUAREM acceleration
- Discrete fragment assignment: `fractional`, `map`, or `sample` (default) post-EM assignment modes
- Parallel BAM scanning and parallel locus EM controlled through one `--threads` setting
- Feather and TSV outputs plus optional annotated BAM output with per-fragment assignment tags

---

## Installation

### Bioconda

```bash
conda install -c conda-forge -c bioconda rigel
```

### PyPI

```bash
pip install rigel-rnaseq

# with the optional HTML QC report ('rigel report'):
pip install 'rigel-rnaseq[report]'
```

The PyPI package name is `rigel-rnaseq` because `rigel` is already taken on
PyPI. The import name and CLI stay `rigel`.

The `[report]` extra pulls in `vl-convert-python` (on conda:
`conda install -c conda-forge vl-convert-python`); without it `rigel report` still builds a
report, minus the fragment-length charts.

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
- C++17-capable compiler
- Runtime dependencies (from `pyproject.toml`): `pysam>=0.22`, `numpy>=1.26`, `pandas>=2.1`,
  `pyarrow>=14.0`, `pyyaml>=6.0`, `scipy>=1.11`

On macOS, install Xcode Command Line Tools first:

```bash
xcode-select --install
```

---

## Quick start

### 1. Build an index

```bash
samtools faidx genome.fa          # the FASTA needs a .fai index
rigel index \
    --fasta genome.fa \
    --gtf annotation.gtf \
    --no-mappability \
    -o index/
```

`rigel index` requires either `--alignable-zarr PATH` (a per-base mappability store built by
the companion `alignable` tool for the same genome and aligner, recommended for real genomes)
or `--no-mappability` to opt out explicitly.

### 2. Quantify a BAM

```bash
rigel quant \
    --bam sample.bam \
    --index index/ \
    -o results/ \
    --tsv
```

Input BAM requirements:

- Name-sorted or collated
- `NH` tag present for multimapper handling
- Splice-junction strand tag available for best strand-model training (`XS` or `ts`, or let Rigel auto-detect)

### 3. Inspect outputs

Outputs are Feather files; `--tsv` above also writes `.tsv` mirrors (or convert afterward with
`rigel export results/ --format tsv`).

```bash
head results/quant.tsv
head results/gene_quant.tsv
head results/nrna_quant.tsv
head results/loci.tsv
cat results/summary.json
```

### 4. Build a QC report (optional)

```bash
rigel report results/ -o results/report.html
```

Produces a single self-contained HTML QC file from the files `rigel quant` already wrote, so
reports can be built later and in bulk. Requires the `[report]` extra (see
[Installation](#pypi)).

The five subcommands are `index`, `quant`, `sim` (a small synthetic scenario from a YAML
file), `export` (Feather to TSV or Parquet) and `report`. `rigel --version` prints the version;
`rigel -v <subcommand>` (`--verbose`) enables DEBUG-level logging. See
[docs/MANUAL.md](docs/MANUAL.md) for every flag.

---

## Output files

| File | Description |
|------|-------------|
| `quant.feather` / `quant.tsv` | Transcript-level abundance table (annotated mRNA + nRNA rows) with `count`, `count_unambig`, `count_em`, `count_spliced`, `tpm`, `tpm_total_rna`, effective lengths, and per-transcript QC columns |
| `gene_quant.feather` / `gene_quant.tsv` | Gene-level aggregates derived from transcript estimates |
| `nrna_quant.feather` / `nrna_quant.tsv` | nRNA-span-level abundance estimates (one row per unique genomic nRNA span) |
| `loci.feather` / `loci.tsv` | Per-locus EM summary |
| `summary.json` | Library protocol, strand specificity, per-category fragment-length summary statistics, the calibration scalars, alignment counts, and global quantification totals |
| `fragment_lengths.feather` | Raw per-bin fragment-length histograms, tidy `(category, length, count)` |
| `calibration_track.feather` / `.bedgraph` | Per-region gDNA solution; the bedGraph is a genome-browser track (IGV / UCSC) |
| `gdna_density_kde.feather` / `gdna_density_regions.feather` | The fitted gDNA-density curve and its training-region rug (written when the density prior is fit) |
| `config.yaml` | Resolved run configuration (parameters, I/O paths). Rerun with `rigel quant --config config.yaml` |
| `report.html` | Optional self-contained QC report, built by `rigel report` (see step 4) |
| `locus_stats.feather` | Optional per-locus statistics, emitted only with `--emit-locus-stats` |
| `annotated.bam` | Optional annotated BAM with per-fragment assignment tags, written with `--annotated-bam` (a second BAM pass); the same records as the input, collated |

`tpm` is normalized over annotated transcripts only; `tpm_total_rna` normalizes over all RNA
(annotated + synthetic nRNA spans) and is comparable to the `nrna_quant` TPM column.

---

## How it works

Rigel runs one native BAM pass feeding three stages: **scan**, **calibrate**, and **quantify**.

```
FASTA + GTF ──▶ rigel index ──▶ regions / boundaries partition + transcript tables
BAM ──▶ 1. scan       C++ single pass: resolve fragments, train the strand and fragment-length
                      models, buffer the ambiguous fragments, tally per-region/boundary mass
    ──▶ 2. calibrate  split each region's unspliced mass into sense RNA / antisense RNA / gDNA
                      (strand tilt + neighbour messages + a gDNA-density prior fitted from the
                      library) ──▶ two Dirichlet scalars per locus
    ──▶ 3. quantify   per-locus EM (VBEM or MAP-EM, SQUAREM, OpenMP across loci) with one
                      component per transcript row plus one gDNA component ──▶ counts, TPM
```

A spliced fragment is certified RNA, so calibration's problem is the unspliced mass. Its
answer enters the EM as a prior that a decisive locus likelihood can override. The design and
its rulings are in [docs/DESIGN.md](docs/DESIGN.md); the derivations are in
[docs/EQUATIONS.md](docs/EQUATIONS.md).

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/MANUAL.md](docs/MANUAL.md) | The user manual: CLI reference, defaults, configuration rules, output schema |
| [docs/DESIGN.md](docs/DESIGN.md) | What is built, and the rulings behind it |
| [docs/EQUATIONS.md](docs/EQUATIONS.md) | The derivations the implementation depends on |
| [docs/SUCCESS.md](docs/SUCCESS.md) | How performance is measured: the accumulator, then calibration against an oracle |
| [docs/TESTING.md](docs/TESTING.md) | The benchmark panels, how to build them, and what the test suite can judge |
| [docs/ROADMAP.md](docs/ROADMAP.md) | The ranked view of what is next for the 0.8.0 release |
| [docs/ISSUES.md](docs/ISSUES.md) | The issue log: open problems and the record of what was refused |
| [docs/TRAPS.md](docs/TRAPS.md) | Mistakes already made, as named rules |
| [docs/PUBLISHING.md](docs/PUBLISHING.md) | Release workflow for PyPI and Bioconda |

---

## Citing Rigel

If you use Rigel in research, cite the repository for now:

> Iyer MK. Rigel: Bayesian RNA-seq quantification with joint mRNA, nascent RNA,
> and genomic DNA modeling. 2026. https://github.com/mkiyer/rigel

---

## License

Rigel is distributed under the [GNU General Public License v3.0](LICENSE).

---

## Development

```bash
conda activate rigel
pip install --no-build-isolation -e ".[dev]"   # rebuild after any src/rigel/native/ change
python -m pytest tests/ -q
ruff check src/ tests/ scripts/ && ruff format src/ tests/
```

