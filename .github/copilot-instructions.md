# Rigel Project Instructions

## Theory

### The problem with 'genes'

The concept of a 'gene' is a convenient abstraction for organizing the transcriptome, but it is not a fundamental biological entity. The 'rigel' tool *must* be designed to handle the complexities and ambiguities that arise from this abstraction, especially when dealing with complex loci with many overlapping transcripts.

- The tool must operate at the transcript level, not the gene level. Grouping of transcripts into genes is sometimes arbitrary and can lead to issues when transcripts overlap multiple genes or when genes have complex splicing patterns.
- We output gene-level summaries for user convenience, but the core quantification and EM modeling should be transcript-centric. This allows for more accurate handling of multimappers and complex loci.


## Environment

- **Always** activate the conda environment before running any command: `conda activate rigel`
- Python 3.12+, C++17
- Dependencies managed via `mamba_env.yaml` (conda-forge + bioconda channels)

## Build

```bash
# After ANY C++ change in src/rigel/native/, you MUST recompile:
conda activate rigel && pip install --no-build-isolation -e .
```

Do not skip recompilation after editing `.cpp`, `.h`, or `CMakeLists.txt` files.

## Testing

```bash
conda activate rigel
pytest tests/ -v                          # all tests
pytest tests/test_em_impl.py -v           # single file
pytest tests/test_em_impl.py::test_name   # single test
pytest tests/ --update-golden             # regenerate golden outputs after intentional changes
```

- Known pre-existing failure: `tests/test_calibration.py::TestStrandLLR::test_biased_toward_ss_favors_rna` — unrelated to scoring/EM.
- Golden output regression files live in `tests/golden/`.

### Test Failure Investigation

- When a test fails, DO NOT reconfigure the test or change tolerances to make the test pass.
- Consider a test failure as a sentinel event with critical information that alludes to a larger problem with the code or methodology. 
- Conduct analysis to understand WHY the test failed, determine the precise root cause of the failure, and evaluate possible solutions
- If the explanation and fix for the failed test is not clear, you should propose a plan for further investigation.

## Linting

Ruff configured in `pyproject.toml`: Python 3.12 target, 100-char line length.

```bash
ruff check src/ tests/
ruff format src/ tests/
```

## Documentation

- Publish implementation plans to the `docs/` directory as markdown files. Use subfolders for organization. Avoid polluting the base `docs/` directory.

## Scripts

- Diagnostic exploration and debugging scripts are commonly needed. For organization, put them in `scripts/debug/` with description names. Avoid polluting the base `scripts/` directory.

## Simulations / Benchmarking / Profiling

- Create temporary scripts in `scripts/benchmark/`, `scripts/profiling/`, or `scripts/simulation/` as needed. Try not to pollute the base `scripts/` directory.

## Benchmarking

`scripts/benchmark.py` contains a framework for running full-scale simulations, comparing multiple methods, and generating summary tables.

### Minimap2 short read RNA-Seq alignment

Minimap2 accepts an annotated junctions file (BED12 format) to improve spliced alignment. See https://lh3.github.io/2025/04/18/short-rna-seq-read-alignment-with-minimap2.

We want to generate a comprehensive splice junction annotation file (BED12) for minimap2. More junctions likely improves spliced alignment sensitivity. There is a built-in script with minimap2 for converting GTF to BED12:

```
paftools.js gff2bed gencode.gtf.gz > gencode.bed  # if not done already
```

How to run minimap2:

```bash
minimap2 -ax splice:sr \
    -j junctions.bed \
    --secondary=yes -N 20 -t 8 \
    genome.fasta R1.fq.gz R2.fq.gz | \
  samtools sort -n -o output.bam
```


## Profiling

`scripts/profiler.py` contains a framework for profiling.

## Simulation

`scripts/synthetic_sim_sweep.py` contains a framework for running synthetic simulations with configurable parameters.

## Synthetic Benchmark Framework

### Directory Layout

```
scripts/benchmark/
├── configs/              # YAML benchmark definitions (sweep grids)
│   ├── locus_simple_baseline.yaml
│   ├── locus_simple_em_prior.yaml
│   ├── locus_simple_em_mode.yaml
│   ├── locus_simple_strand.yaml
│   └── locus_simple_scoring.yaml
├── golden/               # Gold-standard results (TSV + JSON per sweep)
│   ├── locus_simple_baseline/
│   ├── locus_simple_em_prior/
│   ├── locus_simple_em_mode/
│   ├── locus_simple_strand/
│   └── locus_simple_scoring/
├── analyze_golden.py     # Summary statistics across all sweeps
└── analyze_deep.py       # Deep drill-down: nRNA siphon, gDNA leakage, param sensitivity
```

### How to Run Benchmarks

```bash
conda activate rigel

# Run a single sweep (output goes to golden/ for comparison)
python scripts/synthetic_sim_sweep.py \
  -c scripts/benchmark/configs/<name>.yaml \
  -o scripts/benchmark/golden/<name>

# Run all sweeps
for cfg in scripts/benchmark/configs/*.yaml; do
  name=$(basename "$cfg" .yaml)
  python scripts/synthetic_sim_sweep.py -c "$cfg" -o "scripts/benchmark/golden/$name"
done
```

### How to Analyze Results

```bash
# Summary statistics (aggregate error metrics per sweep)
python scripts/benchmark/analyze_golden.py

# Deep analysis (nRNA siphon breakdown, gDNA leakage, per-parameter sensitivity)
python scripts/benchmark/analyze_deep.py
```

### How to Add a New Benchmark Config

Create a YAML file in `scripts/benchmark/configs/`. The sweep framework (`synthetic_sim_sweep.py`) supports these sweepable parameter dimensions:

- **Simulation params**: `strand_specificity`, `gdna_fraction`, per-entity `n_rna_fragments`
- **EM params** (`em_config`): `prior_pseudocount`, `mode` (map/vbem), `convergence_delta`, `max_iterations`
- **BAM scan params** (`scan_config`): (none currently sweepable)
- **Scoring params** (`scoring_config`): `overhang_log_penalty`, `mismatch_log_penalty`

Each dimension in `sweep_dims` creates a Cartesian product grid. Use `n_rna_fragments` + `gdna_fraction` mode (not `n_fragments`) to keep RNA counts fixed while varying gDNA contamination.

### Benchmark Workflow (for Copilot)

When asked to run benchmarks, compare to gold standard, and analyze results:

1. **Recompile** if any C++ changes: `pip install --no-build-isolation -e .`
2. **Run** all benchmark configs (or a specific one) using the commands above
3. **Analyze** with `analyze_golden.py` (summary) and `analyze_deep.py` (deep dive)
4. **Compare** new results against golden TSVs — look for regressions in mRNA/nRNA/gDNA relative error
5. **Report** findings with focus on: (a) nRNA siphon magnitude, (b) gDNA leakage at low contamination, (c) any hyperparameter sensitivity changes
6. **Root cause** any regressions by examining which NTA1/TA1 ratio or ss/gdna_fraction combinations degraded

### Known Benchmark Findings (baseline)

- **nRNA siphon** is the dominant failure mode: when NTA1 >> TA1 (ratio ≥ 5), the nRNA component absorbs mRNA fragments causing up to ~37% mRNA error. This is a fundamental identifiability issue, not a hyperparameter problem.
- **gDNA overestimation** at low contamination (gdna_fraction=0.3): relative error ~1.5–2.0×, improving at higher contamination.
- **Scoring penalties** (overhang, mismatch) show zero sensitivity in synthetic data (no mismatches in simulated reads).
- **EM hyperparameters** (prior_pseudocount 0.1–5.0, MAP vs VBEM) have negligible effect.

## Architecture Summary

Two-stage pipeline: (1) C++ BAM scan → FragmentBuffer, (2) Python/C++ locus-level EM quantification.

- **C++ extensions** (`src/rigel/native/`): `_bam_impl`, `_em_impl`, `_scoring_impl`, `_resolve_impl`, `_cgranges_impl` — compiled via scikit-build-core + nanobind + CMake
- **Scoring model**: mRNA, nRNA, and gDNA candidates scored with strand, fragment-length, overhang, and NM mismatch penalties
- **EM solver**: Per-locus EM with `2*n_t + 1` components (mRNA + nRNA per transcript + one gDNA), SQUAREM acceleration, OpenMP parallel
- **Calibration**: gDNA density, strand balance, and fragment length deconvolution via `calibrate_gdna()` in `calibration.py` — EM mixture model on genomic regions with hard constraints and adaptive shrinkage

## Key Constraints

- Mega-loci (large connected components from multimappers) are difficult to avoid with real data. The EM solver must be robust to loci with thousands of transcripts and millions of fragments.
- `read_length` in resolve_context.h is the sum of aligned exon block spans (M/D/=/X), NOT the full query sequence length. Soft-clipped bases are excluded.
- Chimeric fragments are detected and buffered but **skipped** in scoring — they never enter EM.
- Input BAM must be name-sorted with `NH` tag (minimap2 doesn't set NH; the code handles this via secondary flag detection).

## C++ Compilation Details

- C++17, `-O3`, LTO enabled
- `_scoring_impl` uses `-ffast-math`; `_em_impl` uses `-ffp-contract=fast` (preserves Kahan summation)
- SIMD: `fast_exp.h` provides AVX2/AVX-512/NEON code paths

## Armis2 SLURM Cluster — Full-Scale Simulation Benchmarking

### Unified Benchmark Directory

All simulation data, tool outputs, and indexes are consolidated in:

```
/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_benchmarks/
├── manifest.json                       # simulation parameters & condition definitions
├── truth_abundances_nrna_none.tsv      # ground truth (nrna_none conditions)
├── truth_abundances_nrna_rand.tsv      # ground truth (nrna_rand conditions)
├── sample_sheet.tsv                    # sample metadata
├── rigel_index/                        # pre-built rigel index
├── salmon_index/                       # pre-built salmon index
└── runs/<condition>/
    ├── sim_oracle.bam                  # oracle BAM with true read origins
    ├── sim_R1.fq.gz, sim_R2.fq.gz     # simulated paired-end FASTQ
    ├── rigel_star/
    │   └── annotated.bam               # name-sorted STAR-aligned BAM (input)
    ├── rigel/<config_name>/            # rigel output per named config
    │   ├── quant.feather               # transcript-level quantification
    │   ├── gene_quant.feather          # gene-level quantification
    │   ├── nrna_quant.feather          # nRNA quantification
    │   ├── loci.feather                # locus-level summary
    │   ├── summary.json                # run stats (calibration, strand model, etc.)
    │   └── config.yaml                 # rigel config used
    └── salmon/
        ├── quant.sf.gz                 # salmon transcript-level quantification
        └── quant.genes.sf.gz           # salmon gene-level quantification
```

New rigel runs go into `runs/<condition>/rigel/<config_name>/`.

### Reference Data (original locations)

All references are under `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/`:

| Asset | Path |
|-------|------|
| Rigel index | `.../refs/human/rigel_index` |
| Salmon index | `.../refs/human/salmon_index` |
| STAR index | `.../refs/human/star_index` |
| Genome FASTA | `.../refs/human/genome.fasta.gz` |
| Genome+controls FASTA | `.../refs/human/genome_controls.fasta.bgz` |
| GTF annotation | `.../refs/human/genes.gtf.gz` |
| GTF+controls annotation | `.../refs/human/genes_controls.gtf.gz` |

### Condition Naming Convention

`gdna_{none|low|high}_ss_{0.50|0.90|1.00}_nrna_{none|rand}`

- `gdna_none/low/high` — gDNA contamination level (0%, low%, high%)
- `ss_0.50/0.90/1.00` — strand specificity (0.5=unstranded, 1.0=perfectly stranded)
- `nrna_none/rand` — nascent RNA contamination (none or random)

13 simulation conditions available (14 total; `gdna_low_ss_0.90_nrna_rand` failed).

### Truth Columns

`transcript_id, gene_id, gene_name, ref, strand, mrna_abundance, nrna_abundance, total_rna, n_exons, spliced_length, genomic_span`

### Running Benchmarks

The benchmarking system is in `scripts/benchmarking/` with three subcommands:

```bash
conda activate rigel

# Check what's been run and what's pending
python -m scripts.benchmarking status -c scripts/benchmarking/configs/default.yaml

# Run rigel quant across conditions with named configs
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml

# Run on specific conditions only
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml \
  --conditions gdna_none_ss_1.00_nrna_none gdna_high_ss_0.90_nrna_none

# Dry-run (print commands without executing)
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --dry-run

# Force re-run even if output exists
python -m scripts.benchmarking run -c scripts/benchmarking/configs/default.yaml --force
```

### Benchmark Analysis

```bash
# Analyze all tool outputs vs truth
python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml \
  -o results/benchmark_report

# Analyze specific conditions only
python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/default.yaml \
  -o results/benchmark_report --conditions gdna_none_ss_1.00_nrna_none
```

### Benchmark Config Format

YAML config at `scripts/benchmarking/configs/default.yaml` defines:
- `benchmark_dir` — path to unified benchmark directory
- `rigel_configs` — named configs mapping to `rigel quant` CLI params
- `conditions` — optional condition filter (null = discover all)
- `analysis` — analysis options (output_dir, log2_pseudocount, parse_fastq_truth)

Output goes to `runs/<condition>/rigel/<config_name>/`.

### Benchmarking Package Structure

```
scripts/benchmarking/
├── __init__.py
├── __main__.py        # CLI dispatcher (run, analyze, status)
├── config.py          # BenchmarkConfig dataclass + YAML loader
├── runner.py          # rigel quant CLI runner
├── analysis.py        # metrics + report generation
└── configs/
    └── default.yaml   # default benchmark config
```


