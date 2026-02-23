# hulkrna Benchmarking Runbook

This runbook documents setup and usage for the regional benchmarking workflow built around [scripts/benchmarking/benchmark_region_competition.py](scripts/benchmarking/benchmark_region_competition.py) and [scripts/benchmarking/aggregate_benchmarks.py](scripts/benchmarking/aggregate_benchmarks.py).

## What this benchmark runs

- Transcript-level tools: `hulkrna_mm` (multimap enabled), `salmon`, `kallisto`
- Gene-level tools: `hulkrna_mm`, `salmon`, `kallisto`, optional `htseq-count`
- Genome aligner for hulkrna/htseq input BAMs: `minimap2` or `hisat2`
- Condition grid: regions × gDNA rates × nRNA rates × strand specificities
- Pool-aware evaluation: mature RNA, nascent RNA, genomic DNA, plus RNA vs DNA

gDNA and nRNA abundance models:

- `gDNA_abundance = sum(transcript_abundances) * gDNA_rate`
- `nRNA_abundance = sum(transcript_abundances) * nRNA_rate`
- Per transcript, `t.nrna_abundance = t.abundance * nRNA_rate`

## Requirements

## 1) Python + hulkrna package

From repo root:

```bash
mamba env create -f mamba_env.yaml
conda activate hulkrna
pip install -e .
```

## 2) External command-line tools

Required for all runs:

- `samtools`
- `salmon`
- `kallisto`
- one aligner: `minimap2` (default) or `hisat2` + `hisat2-build`

Optional:

- `htseq-count` in a dedicated env (default env name expected by script: `htseq`)

Example install for missing tools:

```bash
conda install -n hulkrna -c bioconda -c conda-forge salmon kallisto hisat2
conda create -n htseq -c bioconda -c conda-forge htseq
```

## 3) Reference inputs

You need:

- genome FASTA (indexed by `samtools faidx` as needed)
- gene annotation GTF or GTF.GZ
- region definitions via either:
  - `region` list in YAML/CLI (`chr:start-end`), or
  - `region_file` TSV (`chrom<TAB>start1<TAB>end1<TAB>label`)

The 10-region benchmark TSV already exists at [scripts/benchmarking/regions_benchmark.tsv](scripts/benchmarking/regions_benchmark.tsv).

## Config-first workflow (recommended)

The benchmark script supports `--config` YAML with CLI overrides.

Run:

```bash
PYTHONPATH=src conda run -n hulkrna python scripts/benchmarking/benchmark_region_competition.py \
  --config scripts/benchmarking/benchmark_example.yaml
```

Aggregate:

```bash
PYTHONPATH=src conda run -n hulkrna python scripts/benchmarking/aggregate_benchmarks.py \
  --input-dir <benchmark_outdir> \
  --output-dir <benchmark_outdir>
```

## YAML schema (keys you will usually set)

Top-level:

- `genome`, `gtf`, `outdir`
- `aligner`: `minimap2` or `hisat2`
- `threads`
- one of:
  - `region: ["chr:start-end", ...]`
  - `region_file: /path/to/regions.tsv`

Nested sections:

- `simulation`
  - `n_fragments`, `sim_seed`, `pipeline_seed`
  - `frag_mean`, `frag_std`, `frag_min`, `frag_max`
  - `read_length`, `error_rate`
- `gdna`
  - `rates`, `rate_labels`
  - `frag_mean`, `frag_std`, `frag_min`, `frag_max`
- `nrna`
  - `rates`, `rate_labels`
- `abundance`
  - `mode` (`random|uniform|file`), `seed`, `min`, `max`, optional `file`
- `htseq`
  - `enabled` (bool), `conda_env`
- `runtime`
  - `keep_going` (bool)

Standalone:

- `strand_specificities: [0.95, 0.99, 1.0]`
- `verbose: true|false`

## Direct CLI usage

Minimal example:

```bash
PYTHONPATH=src conda run -n hulkrna python scripts/benchmarking/benchmark_region_competition.py \
  --genome /path/genome.fa \
  --gtf /path/genes.gtf.gz \
  --region chr11:5225000-5310000 \
  --outdir /path/out \
  --aligner hisat2 \
  --gdna-rates 0.0,0.25 \
  --nrna-rates 0.0,0.25 \
  --strand-specificities 1.0 \
  --n-fragments 5000 \
  --threads 4
```

Useful flags:

- `--region-file` for multi-region runs
- `--gdna-rate-labels` and `--nrna-rate-labels` for readable condition names
- `--no-htseq` to disable htseq
- `--htseq-conda-env` to change htseq env name
- `--keep-going` to continue after per-region failures

## Output layout

Per run output directory:

- `summary.json`, `summary.csv`, `summary.md`, `diagnostics.md`
- per-region directory containing per-condition outputs
- each condition has `per_transcript_counts.csv` and `per_gene_counts.csv`
- each condition also has `per_pool_counts.csv` with tool-vs-truth pool totals

Pool definitions:

- `mature_rna`: non-`nrna_` transcript truth/prediction
- `nascent_rna`: `nrna_` truth plus hulkrna nascent estimate (pseudoaligners report 0)
- `genomic_dna`: simulated gDNA truth plus hulkrna gDNA estimate (pseudoaligners report 0)
- `rna`: mature_rna + nascent_rna
- `dna`: genomic_dna

Condition directory naming:

- `gdna_<label>_nrna_<label>_ss_<strand_specificity>`

Aggregator outputs:

- `aggregate_summary.json`
- `aggregate_transcript_metrics.csv`
- `aggregate_gene_metrics.csv`
- `aggregate_pool_metrics.csv`
- `aggregate_per_tx.csv`
- `aggregate_summary.md`

## Included config templates

- Quick example config: [scripts/benchmarking/benchmark_example.yaml](scripts/benchmarking/benchmark_example.yaml)
- Full 10-region production config: [scripts/benchmarking/benchmark_full_10_regions.yaml](scripts/benchmarking/benchmark_full_10_regions.yaml)

## One-command launcher

Use [scripts/benchmarking/launch_benchmark.sh](scripts/benchmarking/launch_benchmark.sh) to run benchmark + aggregation in one command:

```bash
bash scripts/benchmarking/launch_benchmark.sh scripts/benchmarking/benchmark_example.yaml
```

Options:

- `-e|--env <conda_env>` (default: `hulkrna`)
- `--skip-aggregate` (run benchmark only)
- `-o|--outdir <path>` (override config outdir for this run)

## Multi-seed runner

You can also use [scripts/benchmarking/run_benchmarks.sh](scripts/benchmarking/run_benchmarks.sh) for repeated seeds:

```bash
bash scripts/benchmarking/run_benchmarks.sh /path/bench_out 50000 101 202 303 404 505
```

This wrapper runs per-seed benchmarks and then aggregates across seeds.