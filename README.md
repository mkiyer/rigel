# hulkrna

Bayesian RNA-seq transcript quantification with built-in gDNA and nascent
RNA deconvolution.  Runs a single-pass BAM scan followed by an EM
estimator, producing per-transcript, per-gene, and per-locus abundance
tables.

## Install

```bash
# From source (requires a C++17 compiler and conda/mamba)
mamba env create -f mamba_env.yaml
conda activate hulkrna
pip install --no-build-isolation -e .
```

## Quick start

```bash
# Build a transcript index from a GTF
hulkrna index -g annotation.gtf -o my_index

# Quantify a BAM file
hulkrna quant -i my_index -b aligned.bam -o results/
```

Output is written to `results/` as TSV files (`transcript.tsv`,
`gene.tsv`, `loci.tsv`) plus a JSON summary.

## Parameters

See [docs/parameters.md](docs/parameters.md) for the full list of
tunable parameters and their defaults.





