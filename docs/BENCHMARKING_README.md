# rigel Benchmarking Runbook

This runbook documents setup and usage for the regional benchmarking workflow built around `benchmark_region_competition.py` and `aggregate_benchmarks.py` in `scripts/benchmarking/`.

## What this benchmark runs

- **Transcript-level tools**: `rigel` (one or more parameterizations x aligners), `salmon`, `kallisto`
- **Gene-level tools**: same as above plus optional `htseq-count` (per aligner)
- **Aligners**: `oracle` (perfect truth alignment), `minimap2`, `hisat2` — each can be parameterized independently
- **Condition grid**: regions × gDNA rates × nRNA rates × strand specificities
- **Pool-aware evaluation**: mature RNA, nascent RNA, genomic DNA, plus RNA vs DNA totals
- **Per-read accuracy**: annotated BAM output with read-level truth comparison

gDNA and nRNA abundance models:

- `gDNA_abundance = sum(transcript_abundances) * gDNA_rate`
- `nRNA_abundance = sum(transcript_abundances) * nRNA_rate`
- Per transcript, `t.nrna_abundance = t.abundance * nRNA_rate`

## Requirements

### 1) Python + rigel package

From repo root:

```bash
mamba env create -f mamba_env.yaml
conda activate rigel
pip install -e .
```

### 2) External command-line tools

Required for all runs:

- `samtools`
- `salmon`
- `kallisto`
- one or more aligners: `minimap2`, `hisat2` + `hisat2-build`

Optional:

- `htseq-count` in a dedicated env (default env name: `htseq`)

Example install for missing tools:

```bash
conda install -n rigel -c bioconda -c conda-forge salmon kallisto hisat2
conda create -n htseq -c bioconda -c conda-forge htseq
```

### 3) Reference inputs

You need:

- genome FASTA (indexed by `samtools faidx` as needed; supports `.bgz`)
- gene annotation GTF or GTF.GZ
- region definitions via either:
  - `region` list in YAML/CLI (`ref:start-end`), or
  - `region_file` TSV (`ref<TAB>start1<TAB>end1<TAB>label`)

The 10-region benchmark TSV already exists at `scripts/benchmarking/regions_benchmark.tsv`.

---

## Config-first workflow (recommended)

The benchmark script supports `--config` YAML with CLI overrides.

Run:

```bash
PYTHONPATH=src conda run -n rigel python scripts/benchmarking/benchmark_region_competition.py \
  --config scripts/benchmarking/benchmark_example.yaml
```

Aggregate:

```bash
PYTHONPATH=src conda run -n rigel python scripts/benchmarking/aggregate_benchmarks.py \
  --input-dir <benchmark_outdir> \
  --output-dir <benchmark_outdir>
```

---

## YAML Configuration Schema

Below is the complete schema for benchmark YAML files. All sections are optional except `genome`, `gtf`, and at least one region.

### Top-level keys

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `genome` | string (path) | **yes** | Path to genome FASTA. |
| `gtf` | string (path) | **yes** | Path to gene annotation GTF or GTF.GZ. |
| `outdir` | string (path) | **yes** | Output directory for all benchmark results. |
| `threads` | int | no | Number of threads (default: 4). |
| `verbose` | bool | no | Emit verbose log messages (default: `false`). |
| `strand_specificities` | list of float | no | Strand specificity values to sweep (default: `[0.95, 0.99, 1.0]`). |

### `region` — Genomic regions to benchmark

A YAML list of region entries. Each entry is either a **plain coordinate string** or a **named region dict**.

```yaml
# Plain coordinate strings:
region:
  - "chr11:5225000-5310000"
  - "chr7:55019017-55211628"

# Named regions (recommended — names appear in output):
region:
  - HBB: chr11:5225000-5310000
  - EGFR: chr7:55019017-55211628
  - BRCA1: chr17:43044295-43170245
```

Each named region is a single-key YAML dict `{Name: "ref:start-end"}`. Coordinates are 1-based inclusive.

Alternatively, use `region_file` to point to a TSV:

```yaml
region_file: /path/to/regions.tsv
```

TSV format: `ref<TAB>start<TAB>end<TAB>label` (1-based start, 0-based end as per BED conventions).

### `aligners` — Aligner parameterizations

When present, this section replaces the simpler top-level `aligner` key. Each entry is a **named aligner configuration** with an aligner `type` and optional parameters.

```yaml
aligners:
  oracle: {}                    # Perfect truth alignment
  minimap2: {}                  # minimap2 with default params
  minimap2_nosec:               # minimap2 without secondary alignments
    type: minimap2
    secondary: false
  hisat2_plain:                 # hisat2 with linear index
    type: hisat2
    k: 50
    index_type: plain
    splice_sites_at_align: true
  hisat2_graph:                 # hisat2 with graph index
    type: hisat2
    k: 50
    index_type: graph
    splice_sites_at_build: true
    exons_at_build: true
```

If `type` is omitted, the key name itself is used as the type (so `oracle: {}` works).

**Supported aligner types and parameters:**

| Type | Parameter | Type | Default | Description |
|------|-----------|------|---------|-------------|
| `oracle` | *(none)* | — | — | Perfect alignment from simulator truth. |
| `minimap2` | `preset` | string | `"splice:sr"` | minimap2 preset string. |
| | `secondary` | bool | `true` | Retain secondary alignments. |
| | `bed_guided` | bool | `true` | Use BED12 annotation for splice junctions. |
| `hisat2` | `k` | int | `50` | Max alignments to report (`-k`). |
| | `secondary` | bool | `true` | Report secondary alignments. |
| | `no_unal` | bool | `true` | Suppress unaligned reads in output. |
| | `index_type` | string | `"plain"` | `"plain"` (linear) or `"graph"` (HGFM). |
| | `splice_sites_at_align` | bool | `true` | Pass known splice sites at alignment time. |
| | `splice_sites_at_build` | bool | `false` | Embed splice sites in graph index. |
| | `exons_at_build` | bool | `false` | Embed exon annotations in graph index. |

**Legacy `aligner` key** (still supported but deprecated):

```yaml
aligner: minimap2          # or: aligner: [minimap2, hisat2]
```

### `rigel_configs` — Rigel parameterizations

A mapping of named configurations, each specifying keyword arguments forwarded to `run_pipeline()`. An empty dict `{}` uses pipeline defaults. Multiple configs are compared against each other in the output.

```yaml
rigel_configs:
  default: {}
  strict_em:
    em_convergence_delta: 1.0e-5
  loose_em:
    em_convergence_delta: 1.0e-3
  custom:
    em_iterations: 2000
    overhang_alpha: 3.0
    gdna_splice_penalty_unannot: -5.0
```

**Supported parameters** (all optional — omitted values use pipeline defaults):

| Parameter | Type | Description |
|-----------|------|-------------|
| `em_convergence_delta` | float | EM convergence threshold (default: 1e-6). |
| `em_iterations` | int | Maximum EM iterations (default: 1000). |
| `em_pseudocount` | float | EM pseudo-count for regularization. |
| `overhang_alpha` | float | Overhang penalty alpha. |
| `mismatch_alpha` | float | Mismatch penalty alpha. |
| `gdna_splice_penalty_unannot` | float | gDNA splice penalty for unannotated junctions. |

### Tool naming in output

Tool names in CSV/JSON output follow these patterns:

- **Single rigel config**: `rigel_<aligner>` (e.g., `rigel_oracle`, `rigel_minimap2`)
- **Multiple rigel configs**: `rigel_<config>_<aligner>` (e.g., `rigel_default_oracle`, `rigel_strict_em_hisat2_plain`)
- **htseq**: `htseq_<aligner>` (e.g., `htseq_oracle`, `htseq_hisat2_plain`)
- **Pseudo-aligners**: `salmon`, `kallisto`

### `simulation` — Read simulation parameters

```yaml
simulation:
  n_fragments: 50000          # Number of fragments to simulate (default: 5000)
  sim_seed: 101               # RNG seed for simulation
  pipeline_seed: 101          # RNG seed for rigel pipeline
  frag_mean: 250.0            # Mean fragment length
  frag_std: 50.0              # Fragment length std dev
  frag_min: 50                # Minimum fragment length
  frag_max: 1000              # Maximum fragment length
  read_length: 150            # Read length in bp
  error_rate: 0.0             # Per-base sequencing error rate (0.0 = pristine)
```

### `gdna` — Genomic DNA contamination

```yaml
gdna:
  rates: [0.0, 0.1, 0.25]    # gDNA contamination rates to sweep
  rate_labels: [none, low, quarter]  # Human-readable labels (must match rates length)
  frag_mean: 350.0            # gDNA fragment length mean
  frag_std: 100.0             # gDNA fragment length std dev
  frag_min: 100               # gDNA minimum fragment length
  frag_max: 1000              # gDNA maximum fragment length
```

### `nrna` — Nascent RNA contamination

```yaml
nrna:
  rates: [0.0, 0.25]         # nRNA contamination rates to sweep
  rate_labels: [none, quarter]  # Human-readable labels
```

### `abundance` — Transcript abundance assignment

```yaml
abundance:
  mode: random                # "random", "uniform", or "file"
  seed: 123                   # RNG seed for random mode
  min: 1.0                    # Minimum abundance (random mode)
  max: 100000.0               # Maximum abundance (random mode)
  file: /path/to/abundances.tsv  # TSV with transcript_id, abundance columns (file mode)
```

### `htseq` — HTSeq-count configuration

```yaml
htseq:
  enabled: true               # Enable/disable htseq-count (default: true)
  conda_env: htseq            # Conda environment with htseq installed
```

### `runtime` — Execution control

```yaml
runtime:
  keep_going: true            # Continue after per-region failures (default: false)
```

---

## Complete YAML example

```yaml
genome: /path/to/genome.fasta.bgz
gtf: /path/to/genes.gtf.gz
outdir: /path/to/benchmark_output
threads: 8

region:
  - HBB: chr11:5225000-5310000
  - EGFR: chr7:55019017-55211628
  - BRCA1: chr17:43044295-43170245

aligners:
  oracle: {}
  minimap2: {}
  hisat2_plain:
    type: hisat2
    k: 50
    index_type: plain
    splice_sites_at_align: true

rigel_configs:
  default: {}

simulation:
  n_fragments: 50000
  sim_seed: 101
  pipeline_seed: 101
  frag_mean: 250.0
  frag_std: 50.0
  frag_min: 50
  frag_max: 1000
  read_length: 150
  error_rate: 0.0

gdna:
  rates: [0.0]
  rate_labels: [none]
  frag_mean: 350.0
  frag_std: 100.0
  frag_min: 100
  frag_max: 1000

nrna:
  rates: [0.0]
  rate_labels: [none]

strand_specificities: [1.0]

abundance:
  mode: random
  seed: 101
  min: 0.01
  max: 100000.0

htseq:
  enabled: true
  conda_env: htseq

runtime:
  keep_going: true

verbose: true
```

---

## Direct CLI usage

Minimal example:

```bash
PYTHONPATH=src conda run -n rigel python scripts/benchmarking/benchmark_region_competition.py \
  --genome /path/genome.fa \
  --gtf /path/genes.gtf.gz \
  --region chr11:5225000-5310000 \
  --outdir /path/out \
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

---

## Output layout

### Per-run top-level outputs

| File | Description |
|------|-------------|
| `summary.json` | Full results as JSON (all regions, conditions, tools). |
| `summary.csv` | Tidy long-form CSV with per-tool metrics. |
| `summary.md` | Markdown summary table. |
| `diagnostics.md` | Detailed diagnostic report. |

### Per-condition outputs

Each condition directory is named `gdna_<label>_nrna_<label>_ss_<spec>` and contains:

| File | Description |
|------|-------------|
| `per_transcript_counts.csv` | Truth vs predicted counts per transcript per tool. |
| `per_gene_counts.csv` | Truth vs predicted counts per gene per tool. |
| `per_pool_counts.csv` | Pool-level (mRNA/nRNA/gDNA) truth vs observed per tool. |

### Per-aligner outputs

Each aligner directory is named `align_<aligner_name>` and contains:

| File | Description |
|------|-------------|
| `reads_namesort.bam` | Name-sorted input BAM (from aligner or oracle). |
| `annotated_<config>.bam` | Annotated BAM with per-fragment assignment tags. |
| `accuracy_<config>/per_read_accuracy.csv` | Per-read truth comparison detail. |
| `accuracy_<config>/read_accuracy_summary.csv` | Aggregate read accuracy stats. |

### Annotated BAM tags

Every read in the annotated BAM carries these tags:

| Tag | Type | Description |
|-----|------|-------------|
| `ZT` | Z (string) | Assigned transcript ID (`"."` for intergenic/gDNA-only). |
| `ZG` | Z (string) | Assigned gene ID (`"."` for intergenic). |
| `ZP` | Z (string) | Pool: `mRNA`, `nRNA`, `gDNA`, `intergenic`, `chimeric`. |
| `ZW` | f (float) | Posterior probability of assignment. |
| `ZC` | Z (string) | Fragment class: `unique`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`, `chimeric`, `intergenic`. |
| `ZH` | i (int) | Primary hit flag: 1 = winning alignment, 0 = secondary. |
| `ZN` | i (int) | Number of EM candidate components. |
| `ZS` | Z (string) | Splice type: `spliced_annot`, `spliced_unannot`, `unspliced`, `ambiguous`. |

### Per-read accuracy verdicts

Each simulated read is compared to truth and assigned a verdict:

| Verdict | Description |
|---------|-------------|
| `correct_transcript` | Correct pool AND correct transcript (or gDNA pool match). |
| `correct_gene` | Correct pool, correct gene, but wrong isoform. |
| `correct_pool` | Correct pool but wrong gene/transcript. |
| `wrong_pool` | Assigned to a different pool (e.g., mRNA assigned as gDNA). |
| `intergenic` | Assigned to intergenic when truth was RNA or gDNA. |

### Pool definitions

| Pool | Truth source | Prediction source |
|------|-------------|-------------------|
| `mature_rna` | Non-`nrna_` transcript reads | rigel mRNA counts |
| `nascent_rna` | `nrna_` prefixed reads | rigel nascent estimate |
| `genomic_dna` | `gdna:` prefixed reads | rigel gDNA count |
| `rna` | mature_rna + nascent_rna | combined RNA |
| `dna` | genomic_dna | gDNA total |

### Aggregator outputs

When `aggregate_benchmarks.py` is run, it produces:

| File | Description |
|------|-------------|
| `aggregate_summary.json` | Combined JSON across regions. |
| `aggregate_transcript_metrics.csv` | Per-tool transcript-level metrics. |
| `aggregate_gene_metrics.csv` | Per-tool gene-level metrics. |
| `aggregate_pool_metrics.csv` | Per-tool pool-level metrics. |
| `aggregate_per_tx.csv` | Per-transcript counts across all regions. |
| `aggregate_summary.md` | Markdown summary. |

---

## One-command launcher

Use `scripts/benchmarking/launch_benchmark.sh` to run benchmark + aggregation in one command:

```bash
bash scripts/benchmarking/launch_benchmark.sh scripts/benchmarking/benchmark_example.yaml
```

Options:

- `-e|--env <conda_env>` (default: `rigel`)
- `--skip-aggregate` (run benchmark only)
- `-o|--outdir <path>` (override config outdir for this run)

## Included config templates

| File | Description |
|------|-------------|
| `benchmark_example.yaml` | Full example with multiple aligners and rigel configs. |
| `benchmark_pristine_10_regions.yaml` | Pristine (no error/contamination) 10-region benchmark. |
| `benchmark_full_10_regions.yaml` | Full 10-region with gDNA/nRNA sweeps. |
| `benchmark_small.yaml` | Quick smoke test config. |
