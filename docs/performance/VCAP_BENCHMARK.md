# VCaP Prostate Benchmark

End-to-end benchmark comparing **rigel** (VBEM, MAP-EM, oracle alignment) against **salmon** and **kallisto** on a large-scale simulated VCaP prostate cancer RNA-seq dataset.

## Overview

- **Cell line**: VCaP (prostate cancer, CCLE)
- **Simulation**: 50M RNA + 25M gDNA fragments = 75M total paired-end reads
- **Condition**: `gdna_high_ss_0.90_nrna_none` (50% gDNA, strand specificity 0.90, no nRNA)
- **Read length**: 101 bp, fragment size mean=250 (σ=50)
- **Abundances**: Derived from real VCaP salmon quantification
- **Genome**: GRCh38 + spike-in controls (`genome_controls.fasta.bgz`)
- **Annotation**: GENCODE + controls (`genes_controls.gtf.gz`)

## Prerequisites

```bash
conda activate rigel
pip install --no-build-isolation -e .   # rebuild if any C++ changes
```

### Required tools (all available in the rigel conda environment)

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | 2.30+ | Splice-aware short-read alignment |
| samtools | 1.23+ | BAM sorting (name-sorted output) |
| paftools.js + k8 | (bundled with minimap2) | GTF → BED12 junction file |
| salmon | 1.11+ | Pseudo-alignment quantification |
| kallisto | 0.52+ | Pseudo-alignment quantification |
| rigel | (this repo) | Transcript quantification with gDNA/nRNA modeling |

### Required indexes

| Index | Path | Build command |
|-------|------|---------------|
| Rigel | `.../hulkrna/refs/human/rigel_index` | `rigel index --fasta genome.fa --gtf annotation.gtf -o rigel_index/` |
| Salmon | `.../hulkrna/refs/human/salmon_index` | `salmon index -t transcripts.fa -i salmon_index/` |
| Kallisto | `.../hulkrna/refs/human/kallisto.idx` | `kallisto index -i kallisto.idx transcripts.fa` |

## Directory Layout

```
/scratch/.../rigel_benchmarks/ccle_vcap_prostate/
├── manifest.json                         # simulation parameters
├── truth_abundances_nrna_none.tsv        # ground truth abundances
├── minimap2_junctions.bed                # BED12 splice junctions (auto-generated)
└── runs/
    └── gdna_high_ss_0.90_nrna_none/      # condition directory
        ├── sim_R1.fq.gz                  # simulated FASTQ (R1)
        ├── sim_R2.fq.gz                  # simulated FASTQ (R2)
        ├── sim_oracle.bam                # oracle BAM (true alignments, 6.3 GB)
        ├── minimap2/
        │   └── aligned.bam               # minimap2 name-sorted BAM (6.6 GB)
        ├── salmon/
        │   └── quant.sf.gz               # salmon transcript quantification
        ├── kallisto/
        │   ├── abundance.tsv             # kallisto transcript quantification
        │   ├── abundance.h5
        │   └── run_info.json
        └── rigel/
            ├── vbem/                     # rigel VBEM output
            ├── map/                      # rigel MAP-EM output
            └── oracle/                   # rigel on oracle BAM
```

## Configuration

The benchmark config is at `scripts/benchmarking/configs/vcap.yaml`:

```yaml
benchmark_dir: /scratch/mkiyer_root/mkiyer0/shared_data/rigel_benchmarks/ccle_vcap_prostate
rigel_index: /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index
bam_relpath: minimap2/aligned.bam
salmon_index: /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/salmon_index
kallisto_index: /scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/kallisto.idx
threads: 8
seed: 42

rigel_configs:
  vbem:
    em_mode: vbem
  map:
    em_mode: map
  oracle:
    em_mode: vbem
    bam_relpath: sim_oracle.bam

analysis:
  output_dir: results/vcap
  log2_pseudocount: 1.0
  parse_fastq_truth: false
```

## Running the Benchmark

All commands use the benchmarking package via `python -m scripts.benchmarking`.

### Step 0: Check status

```bash
python -m scripts.benchmarking status -c scripts/benchmarking/configs/vcap.yaml
```

### Step 1: Align reads with minimap2

Generates a name-sorted BAM from simulated FASTQ using minimap2 with splice junction annotations. Takes ~70 min for 75M reads on 8 threads. Peak memory ~17 GB.

```bash
python -m scripts.benchmarking align \
  -c scripts/benchmarking/configs/vcap.yaml -v
```

This runs:
```
minimap2 -ax splice:sr --secondary=yes -N 20 -t 8 \
  -j minimap2_junctions.bed \
  genome_controls.fasta.bgz R1.fq.gz R2.fq.gz | \
  samtools sort -n -@ 4 -o aligned.bam
```

The BED12 junctions file is auto-generated on first run via `paftools.js gff2bed`.

Use `--dry-run` to preview commands without executing.

### Step 2: Run salmon and kallisto

Runs pseudo-alignment tools on the simulated FASTQ. Salmon ~37 min, kallisto ~16 min.

```bash
python -m scripts.benchmarking run-tools \
  -c scripts/benchmarking/configs/vcap.yaml -v
```

Options:
- `--tools salmon` or `--tools kallisto` to run a single tool
- `--force` to re-run even if output exists

### Step 3: Run rigel quantification

Runs rigel quant for each named config. The `vbem` and `map` configs use the minimap2 BAM; the `oracle` config uses the simulator's oracle BAM.

```bash
# Run all configs
python -m scripts.benchmarking run \
  -c scripts/benchmarking/configs/vcap.yaml -v

# Run specific configs only
python -m scripts.benchmarking run \
  -c scripts/benchmarking/configs/vcap.yaml --configs vbem map -v
```

### Step 4: Analyze results

Compares all tool outputs against ground truth.

```bash
python -m scripts.benchmarking analyze \
  -c scripts/benchmarking/configs/vcap.yaml \
  -o results/vcap
```


**Key**: The `--mem=96G` is critical. The oracle BAM buffers ~64.5M fragments with dense multimapping and requires ~60-80 GB during EM quantification. The minimap2-aligned BAM configs (`vbem`, `map`) use less memory (~25-30 GB) due to fewer secondary alignments.



## Current State

| Step | Status | Notes |
|------|--------|-------|
| minimap2 alignment | **Done** | 6.6 GB BAM, 4277s |
| salmon | **Done** | quant.sf.gz, 2241s |
| kallisto | **Done** | abundance.tsv, 949s |
| rigel vbem (minimap2) | Pending | |
| rigel map (minimap2) | Pending | |
| rigel oracle | Pending | Requires >48 GB RAM — submit via SLURM with `--mem=96G` |
| analysis | Pending | Waiting for rigel runs |

## Benchmarking Package Structure

```
scripts/benchmarking/
├── __init__.py
├── __main__.py        # CLI: align, run, run-tools, analyze, status
├── config.py          # BenchmarkConfig + YAML loader
├── runner.py          # rigel quant runner
├── aligner.py         # minimap2 alignment runner
├── tools.py           # salmon/kallisto runners
├── analysis.py        # metrics + report generation
└── configs/
    └── vcap.yaml      # VCaP benchmark config
```

## Simulation Parameters (from manifest.json)

| Parameter | Value |
|-----------|-------|
| RNA fragments | 50,000,000 |
| gDNA fragments | 25,000,000 (rate=0.5) |
| Strand specificity | 0.90 |
| nRNA | none |
| Fragment size | mean=250, σ=50, range=[50, 1000] |
| gDNA fragment size | mean=250, σ=100, range=[100, 1000] |
| Read length | 101 bp |
| Error rate | 0.0 |
| Seed | 42 |
| Abundance source | Real VCaP salmon quant (CCLE) |
