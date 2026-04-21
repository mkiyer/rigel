# Benchmarking Guide

This guide explains how to set up, run, and analyze full-scale simulation benchmarks for rigel using the `scripts/benchmarking/` package.

## Prerequisites

1. **Conda environment** — all commands assume you have activated the rigel environment:
   ```bash
   conda activate rigel
   ```

2. **Rigel installed** — rigel must be installed in the environment so the `rigel` CLI is available:
   ```bash
   pip install --no-build-isolation -e .
   ```

3. **Benchmark directory** — a directory containing simulation data, aligned BAMs, and ground truth. See [Benchmark Directory Layout](#benchmark-directory-layout) below for the required structure.

---

## Quick Start

```bash
conda activate rigel
cd /path/to/rigel

# 1. Check what's available and what's been run
python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml

# 2. Run rigel quant on all conditions (skips those already done)
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml

# 3. Analyze results vs ground truth
python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml \
  -o results/benchmark_report
```

---

## Benchmark Directory Layout

The benchmarking system expects a specific directory structure. The benchmark directory is the single source of truth for all inputs, outputs, and ground truth.

```
<benchmark_dir>/
├── manifest.json                         # simulation parameters & condition definitions
├── truth_abundances_nrna_none.tsv        # ground truth (nrna_none conditions)
├── truth_abundances_nrna_rand.tsv        # ground truth (nrna_rand conditions)
├── sample_sheet.tsv                      # sample metadata (optional)
├── rigel_index/                          # pre-built rigel index
│   ├── transcripts.feather
│   ├── regions.feather
│   ├── sj.feather
│   ├── intervals.feather
│   └── ref_lengths.feather
├── salmon_index/                         # pre-built salmon index (optional)
└── runs/
    └── <condition_name>/
        ├── sim_oracle.bam                # oracle BAM with true read origins
        ├── sim_R1.fq.gz                  # simulated R1 FASTQ
        ├── sim_R2.fq.gz                  # simulated R2 FASTQ
        ├── rigel_star/
        │   └── annotated.bam            # name-sorted STAR-aligned BAM (input)
        ├── rigel/
        │   └── <config_name>/           # rigel output (one per named config)
        │       ├── quant.feather
        │       ├── gene_quant.feather
        │       ├── nrna_quant.feather
        │       ├── loci.feather
        │       ├── summary.json
        │       └── config.yaml
        └── salmon/                       # pre-computed salmon output (optional)
            ├── quant.sf.gz
            └── quant.genes.sf.gz
```

### Required files

To run benchmarks (the `run` subcommand), you need:

- `rigel_index/` — a rigel index built from the same genome + GTF used for simulation
- `runs/<condition>/rigel_star/annotated.bam` — name-sorted aligned BAM for each condition

To analyze results (the `analyze` subcommand), you additionally need:

- `manifest.json` — contains simulation parameters and per-condition metadata
- `truth_abundances_nrna_none.tsv` and/or `truth_abundances_nrna_rand.tsv` — ground truth transcript abundances

### How to create a benchmark directory from scratch

If you're setting up a new benchmark (not using the existing armis2 data), follow these steps:

1. **Build the rigel index:**
   ```bash
   rigel index --fasta genome.fasta --gtf genes.gtf -o <benchmark_dir>/rigel_index
   ```

2. **Run the simulation** to generate FASTQ and truth files. The simulator must produce:
   - Paired-end FASTQ files (`sim_R1.fq.gz`, `sim_R2.fq.gz`)
   - Truth abundance TSV files
   - A `manifest.json` describing conditions and simulation parameters

3. **Align reads** to produce name-sorted BAMs:
   ```bash
   minimap2 -ax splice:sr -j junctions.bed --secondary=yes -N 20 -t 8 \
     genome.fasta sim_R1.fq.gz sim_R2.fq.gz | \
     samtools sort -n -o annotated.bam
   ```
   Place the BAM at `runs/<condition>/rigel_star/annotated.bam`.

4. **(Optional) Run salmon** for comparison:
   ```bash
   salmon quant -i <benchmark_dir>/salmon_index \
     -l A -1 sim_R1.fq.gz -2 sim_R2.fq.gz -o salmon_out
   ```
   Copy or move the output to `runs/<condition>/salmon/`.

### manifest.json format

The manifest must contain a `simulation` block and a `conditions` list:

```json
{
  "simulation": {
    "n_rna_fragments": 10000000,
    "sim_seed": 42,
    "frag_mean": 250.0,
    "frag_std": 50.0,
    "read_length": 101,
    "error_rate": 0.0
  },
  "conditions": [
    {
      "name": "gdna_none_ss_1.00_nrna_none",
      "gdna_rate": 0.0,
      "strand_specificity": 1.0,
      "nrna_label": "none",
      "n_rna": 10000000,
      "n_gdna": 0,
      "truth_abundances": "truth_abundances_nrna_none.tsv"
    }
  ]
}
```

### Truth abundance file format

Tab-separated with these columns:

| Column | Description |
|--------|-------------|
| `transcript_id` | Transcript identifier (e.g. ENST00000456328.2) |
| `gene_id` | Gene identifier (e.g. ENSG00000290825.1) |
| `gene_name` | Gene symbol (e.g. DDX11L2) |
| `ref` | Chromosome / reference name |
| `strand` | + or - |
| `mrna_abundance` | mRNA abundance (fractional, sums to ~1M) |
| `nrna_abundance` | nRNA abundance (0 for nrna_none conditions) |
| `total_rna` | mrna_abundance + nrna_abundance |
| `n_exons` | Number of exons |
| `spliced_length` | Transcript spliced length in bp |
| `genomic_span` | Genomic span in bp |

The analysis script automatically scales these fractional abundances to match the simulated fragment count from the manifest.

### Condition naming convention

Conditions follow the pattern: `gdna_{none|low|high}_ss_{0.50|0.90|1.00}_nrna_{none|rand}`

| Dimension | Values | Meaning |
|-----------|--------|---------|
| `gdna` | `none`, `low`, `high` | gDNA contamination level (0%, 20%, 50%) |
| `ss` | `0.50`, `0.90`, `1.00` | Strand specificity (unstranded → perfectly stranded) |
| `nrna` | `none`, `rand` | Nascent RNA contamination (none or random per-transcript) |

---

## Config File Reference

The YAML config file controls all three subcommands. Create a new config by copying `scripts/benchmarking/configs/default.yaml`.

### Full config example

```yaml
# ── Paths ──────────────────────────────────────────────────────────
benchmark_dir: /path/to/hulkrna_benchmarks
rigel_index: /path/to/hulkrna_benchmarks/rigel_index

# BAM input path relative to each condition directory
bam_relpath: rigel_star/annotated.bam

# ── Global defaults ───────────────────────────────────────────────
threads: 8
seed: 42

# ── Named rigel configurations ────────────────────────────────────
#
# Each key is a config name; output goes to runs/<cond>/rigel/<name>/
# Values map directly to `rigel quant` CLI params.
# Omitted params use rigel defaults.
rigel_configs:
  default:
    em_mode: vbem
    assignment_mode: sample

  map_high_prior:
    em_mode: map
    prior_pseudocount: 5.0

  fractional:
    assignment_mode: fractional

  aggressive_pruning:
    pruning_min_posterior: 0.01
    overhang_alpha: 0.001

# ── Conditions ─────────────────────────────────────────────────────
# Omit or set to null to auto-discover all conditions.
conditions:
  - gdna_none_ss_1.00_nrna_none
  - gdna_high_ss_0.90_nrna_none

# ── Analysis options ──────────────────────────────────────────────
analysis:
  output_dir: results/benchmark_report
  log2_pseudocount: 1.0
  parse_fastq_truth: false
```

### Config fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `benchmark_dir` | path | (required) | Root of the benchmark directory |
| `rigel_index` | path | `<benchmark_dir>/rigel_index` | Path to the rigel index |
| `bam_relpath` | string | `rigel_star/annotated.bam` | BAM input path relative to each condition directory |
| `threads` | int | `4` | Default thread count for `rigel quant` |
| `seed` | int | `42` | Random seed |
| `rigel_configs` | dict | `{}` | Named configurations (see below) |
| `conditions` | list or null | `null` | Condition filter; null = auto-discover all |
| `analysis.output_dir` | string | `results/benchmark_report` | Analysis output directory |
| `analysis.log2_pseudocount` | float | `1.0` | Pseudocount for log-space correlation |
| `analysis.parse_fastq_truth` | bool | `false` | Parse FASTQ for per-fragment truth (slow) |

### Available rigel config parameters

Each named config under `rigel_configs` accepts any `rigel quant` CLI parameter using underscores. The runner translates these to CLI flags automatically.

| Config key | CLI flag | Default | Description |
|------------|----------|---------|-------------|
| `em_mode` | `--em-mode` | `vbem` | EM variant: `vbem` or `map` |
| `prior_pseudocount` | `--prior-pseudocount` | `1.0` | Prior budget in virtual fragments |
| `em_iterations` | `--em-iterations` | `1000` | Max EM iterations |
| `em_convergence_delta` | `--em-convergence-delta` | `1e-6` | Convergence threshold |
| `assignment_mode` | `--assignment-mode` | `sample` | Post-EM assignment: `sample`, `fractional`, `map` |
| `assignment_min_posterior` | `--assignment-min-posterior` | `0.01` | Min posterior for assignment eligibility |
| `overhang_alpha` | `--overhang-alpha` | `0.1` | Per-base overhang penalty (0=hard gate, 1=no penalty) |
| `mismatch_alpha` | `--mismatch-alpha` | `0.1` | Per-mismatch penalty (0=hard gate, 1=no penalty) |
| `gdna_splice_penalty_unannot` | `--gdna-splice-penalty-unannot` | `0.01` | gDNA splice penalty for unannotated junctions |
| `pruning_min_posterior` | `--pruning-min-posterior` | `1e-4` | Min posterior for candidate pruning |
| `threads` | `--threads` | from global | Thread count (overrides global) |
| `seed` | `--seed` | from global | Random seed (overrides global) |
| `buffer_size` | `--buffer-size` | `4` | Buffer size in GiB |
| `include_multimap` | `--include-multimap` | `true` | Include multimappers |
| `sj_strand_tag` | `--sj-strand-tag` | `auto` | BAM tag for splice junction strand |
| `annotated_bam` | `--annotated-bam` | none | Write annotated BAM |
| `tsv` | `--tsv` | `false` | Also write TSV outputs |

You can also override the BAM input per config with the special `bam_relpath` key:

```yaml
rigel_configs:
  oracle:
    bam_relpath: sim_oracle.bam
    em_mode: vbem
```

---

## Subcommands

### `status` — Check what's been run

```bash
python -m scripts.benchmarking status -c <config.yaml>
```

Prints a matrix showing each condition, whether its BAM exists, the completion status of each named rigel config, and what other tools (salmon, kallisto) are available.

Example output:

```
Rigel configs (2):
  default: em_mode=vbem, assignment_mode=sample
  map_high_prior: em_mode=map, prior_pseudocount=5.0

Conditions: 13 selected / 13 available

condition                    bam      rigel/default  rigel/map_high_prior  other_tools
───────────────────────────  ───────  ─────────────  ────────────────────  ───────────
gdna_none_ss_1.00_nrna_none  ok       done           —                     salmon
gdna_high_ss_0.90_nrna_none  ok       —              —                     salmon
...

Done: 1  Pending: 25  Missing BAM: 2
```

### `run` — Execute rigel quant

```bash
# Run all conditions × all configs (skip if output exists)
python -m scripts.benchmarking run -c <config.yaml>

# Specific conditions only
python -m scripts.benchmarking run -c <config.yaml> \
  --conditions gdna_none_ss_1.00_nrna_none gdna_high_ss_0.90_nrna_none

# Specific configs only
python -m scripts.benchmarking run -c <config.yaml> \
  --configs default map_high_prior

# Preview commands without executing
python -m scripts.benchmarking run -c <config.yaml> --dry-run

# Force re-run even if output already exists
python -m scripts.benchmarking run -c <config.yaml> --force
```

The runner constructs and executes `rigel quant` subprocess commands. By default it skips conditions where `quant.feather` already exists in the output directory. Output goes to `runs/<condition>/rigel/<config_name>/`.

A typical run on 10M-fragment simulated data takes 3–4 minutes per condition with 8 threads.

### `analyze` — Compare outputs against ground truth

```bash
# Analyze all conditions and tools
python -m scripts.benchmarking analyze -c <config.yaml> -o results/benchmark_report

# Analyze specific conditions
python -m scripts.benchmarking analyze -c <config.yaml> \
  -o results/benchmark_report \
  --conditions gdna_none_ss_1.00_nrna_none
```

The analyzer discovers all available tool outputs for each condition:
- **rigel outputs** — any `rigel/<config>/` directory containing `quant.feather`
- **Legacy rigel outputs** — any `rigel_*` directory (e.g. `rigel_star/`) containing `quant.feather`
- **salmon** — `salmon/quant.sf.gz`
- **kallisto** — `kallisto/abundance.tsv`

#### Analysis outputs

| File | Description |
|------|-------------|
| `transcript_metrics.csv` | Per-tool, per-condition transcript-level accuracy metrics |
| `gene_metrics.csv` | Per-tool, per-condition gene-level accuracy metrics |
| `pool_summary.csv` | Pool-level (mRNA/nRNA/gDNA) truth vs predicted (rigel only) |
| `rigel_details.csv` | Rigel calibration details (strand model, gDNA density, etc.) |
| `stratified_metrics.csv` | Metrics stratified by expression level bins |
| `per_transcript_detail.parquet` | Full per-transcript truth vs predicted (all tools, all conditions) |
| `condition_summary.csv` | Condition metadata and fragment counts |
| `report.md` | Human-readable markdown summary |

#### Metrics computed

- **Correlation**: Pearson R, Spearman R, log2 Pearson R, log2 Spearman R
- **Error**: MAE, RMSE, MAPE, WARE (weighted absolute relative error)
- **Classification**: precision, recall, F1 (expressed vs not expressed)
- **Pool-level** (rigel): mRNA/nRNA/gDNA fragment count predictions vs truth, with relative errors and calibration estimates

---

## Common Workflows

### Test a parameter change

1. Create a new named config in your YAML:
   ```yaml
   rigel_configs:
     default:
       em_mode: vbem
       assignment_mode: sample
     my_experiment:
       em_mode: vbem
       prior_pseudocount: 3.0
   ```

2. Run only the new config:
   ```bash
   python -m scripts.benchmarking run -c <config.yaml> --configs my_experiment
   ```

3. Analyze to compare against baseline and salmon:
   ```bash
   python -m scripts.benchmarking analyze -c <config.yaml> -o results/my_experiment
   ```
   The transcript_metrics.csv will contain rows for `rigel/default`, `rigel/my_experiment`, and `salmon` for each condition, making it easy to compare.

### Run on a subset of conditions

Use `--conditions` on any subcommand to restrict to specific conditions:

```bash
# Quick test on the simplest condition
python -m scripts.benchmarking run -c <config.yaml> \
  --conditions gdna_none_ss_1.00_nrna_none

# Or set conditions permanently in the config:
# conditions:
#   - gdna_none_ss_1.00_nrna_none
#   - gdna_high_ss_0.90_nrna_none
```

### Compare rigel to salmon

Salmon outputs are pre-computed in the benchmark directory. The analysis automatically discovers and includes them. The gene-level analysis handles salmon's broken `quant.genes.sf.gz` by aggregating transcript-level counts to genes using the truth file's transcript-to-gene mapping.

### Investigate a specific transcript

Load the per-transcript detail parquet file:

```python
import pandas as pd

detail = pd.read_parquet("results/benchmark_report/per_transcript_detail.parquet")

# Find worst-predicted transcripts for rigel
rigel = detail[detail["tool"] == "rigel/default"]
worst = rigel.nlargest(20, "abs_error")[
    ["transcript_id", "gene_name", "truth_scaled", "predicted", "abs_error", "rel_error", "condition"]
]
print(worst.to_string())
```

---

## Existing Benchmark Data (Armis2)

Pre-built benchmark data is available on the armis2 cluster at:

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/
```

This contains 13 simulation conditions (the 14th, `gdna_low_ss_0.90_nrna_rand`, failed during simulation and is excluded). The simulation uses CAPAN-1 pancreatic cancer cell line expression profiles as the basis for transcript abundances.

**Simulation parameters:**
- 10M RNA fragments per condition
- Fragment length: 250 ± 50 bp (range 50–1000)
- Read length: 101 bp
- Error rate: 0.0 (perfect reads)
- Seed: 42

**Conditions sweep:**

| gDNA | gDNA rate | nRNA | Strand specificity | Description |
|------|-----------|------|--------------------|-------------|
| none | 0% | none | 0.50, 0.90, 1.00 | Pure mRNA, varying strandedness |
| low | 20% | none | 0.50, 0.90, 1.00 | Low gDNA contamination |
| high | 50% | none | 0.50, 0.90, 1.00 | High gDNA contamination |
| none | 0% | rand | 0.50, 0.90, 1.00 | Nascent RNA contamination |
| low | 20% | rand | 0.50 | gDNA + nRNA (1 condition) |

The default config file (`scripts/benchmarking/configs/default.yaml`) is pre-configured to use this data.

---

## Package Structure

```
scripts/benchmarking/
├── __init__.py            # package init
├── __main__.py            # CLI dispatcher (run, analyze, status)
├── config.py              # BenchmarkConfig dataclass + YAML loader
├── runner.py              # rigel quant CLI runner
├── analysis.py            # metrics + report generation
└── configs/
    └── default.yaml       # default config for armis2 benchmark data
```

## Troubleshooting

**`'rigel' command not found`** — The runner uses subprocess to call `rigel quant`. Make sure the rigel package is installed (`pip install --no-build-isolation -e .`) and the conda environment is activated before running benchmarks.

**`BAM not found` / `MISSING` in status** — The BAM file at `<bam_relpath>` doesn't exist for that condition. Check that alignment was completed and files are in the correct location. Two conditions (`gdna_low_ss_0.50_nrna_none` and `gdna_none_ss_0.50_nrna_rand`) may be missing BAMs if their simulations failed.

**Analysis finds no rigel outputs** — Run `status` first to confirm outputs exist. The analysis looks for `quant.feather` inside `rigel/<config>/` and `rigel_*/` directories. If you ran rigel outside the benchmarking system, ensure the output is in the expected location.

**Slow analysis** — The majority of analysis time is spent loading feather files (254K transcripts × n conditions). Restrict to specific conditions with `--conditions` to speed up iteration. Avoid `parse_fastq_truth: true` unless you need per-fragment ground truth — it takes ~10 minutes per condition to parse 10M FASTQ reads.
