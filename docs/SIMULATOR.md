# Rigel Simulator Manual

Rigel includes two simulator workflows:

- `scripts/sim/simulate_suite.py` builds a random synthetic mini-genome, writes a GTF, assigns random transcript abundances, and generates a grid of FASTQ/oracle-BAM simulation conditions.
- `scripts/sim/simulate_reads.py` simulates reads from an existing FASTA/GTF pair using a YAML config.

Use the mini-genome suite when you want hard synthetic cases with random gene locations, strands, exon structures, isoforms, antisense overlaps, gDNA contamination, nRNA spike-ins, and optional hybrid capture. Use the existing-reference workflow when you already have a FASTA, GTF, and optional probe file.

## Environment

Run simulator commands from the repository root:

```bash
conda activate rigel
```

For ordinary development runs, keep outputs outside the repository or under a scratch directory. Simulation outputs can become large quickly.

## Random Mini-Genome Suite

The suite workflow generates a synthetic reference and then simulates all selected conditions.

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/synthetic_smoke \
  --profile smoke \
  --n-rna 100000
```

Suite configs can provide the same defaults:

```bash
python scripts/sim/simulate_suite.py \
  --config scripts/sim/configs/hybrid_capture_500kb.yaml \
  --outdir /scratch/$USER/rigel_sim/hybrid_capture_500kb
```

Explicit CLI flags override values from `--config`.

Main outputs:

```text
<outdir>/
├── reference/
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genes.gtf
│   └── capture_probes_<label>.tsv  # only when generated capture configs are enabled
├── truth_abundances_nrna_<label>.tsv
├── manifest.json
└── gdna_<label>_ss_<value>_nrna_<label>[_capture_<label>]/
    ├── sim_R1.fq.gz
    ├── sim_R2.fq.gz
    └── sim_oracle.bam
```

The `_capture_<label>` suffix is added only when a config uses a capture sweep.
Legacy single-capture configs keep the older condition names.

The default condition grid is controlled by `--profile`:

- `smoke`: no nRNA spike-in, 5 gDNA rates, 2 strand-specificity settings.
- `full`: 4 nRNA ratios, 5 gDNA rates, 2 strand-specificity settings.

Useful options:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/synthetic_full \
  --profile full \
  --genome-length 10000000 \
  --n-genes 50 \
  --n-rna 1000000 \
  --seed 42 \
  --n-workers 4
```

Reference-generation options include `--genome-length`, `--n-genes`,
`--min-isoforms`, `--max-isoforms`, `--target-transcripts`, and
`--antisense-fraction`. Use `--target-transcripts 0` for uniform isoform counts
between the min/max bounds.

To generate only the random reference files:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/reference_only \
  --reference-only
```

To rerun only missing condition outputs:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/synthetic_smoke \
  --profile smoke \
  --skip-existing
```

To generate only selected conditions, pass condition directory names:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/synthetic_smoke \
  --conditions gdna_none_ss_0.99_nrna_none gdna_high_ss_0.50_nrna_none
```

## Generated Hybrid Capture Probes

The mini-genome suite can generate transcript-coordinate capture probes automatically. This lets random transcriptome simulations stress capture enrichment without requiring a hand-authored probe file.

```bash
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/capture_smoke \
  --profile smoke \
  --n-rna 100000 \
  --capture-fraction 0.25 \
  --probe-length 120 \
  --probe-density 1.0
```

Generated probes are written to:

```text
<outdir>/reference/capture_probes.tsv
```

When you configure more than one capture setting, each enabled setting gets a
label-specific probe file:

```text
<outdir>/reference/capture_probes_on.tsv
<outdir>/reference/capture_probes_dense.tsv
```

The file is a tab-separated transcript-coordinate table:

```text
transcript_id    start    end
GENE0001.1       11       131
GENE0001.1       142      262
```

Coordinates are 0-based, half-open transcript coordinates. `simulate_suite.py` passes this file into the same hybrid-capture model used by the whole-genome simulator.

Capture probe design parameters:

- `--capture-fraction`: fraction of eligible transcripts selected for targeting. `0.0` disables generated capture. Eligible transcripts must be at least `--probe-length` bases long.
- `--probe-length`: probe length in bases. Default: `120`.
- `--probe-density`: fraction of the maximum non-overlapping tiling to emit for each captured transcript. Must be between `0` and `1`. Default: `1.0`.
- `--capture-off-target-weight`: baseline weight for every legal fragment start. Default: `1.0`.
- `--capture-binding-per-base`: additional weight per fragment/probe overlap base. Default: `10.0`.
- `--capture-gdna-split-penalty`: multiplier for probes split by introns when projected to gDNA or pre-mRNA. Default: `0.2`.
- `--capture-min-overlap`: minimum fragment/probe overlap required before capture weight is added. Default: `1`.

Probe tiling rules:

- Transcripts are partitioned into captured and uncaptured groups using the suite seed.
- Uncaptured transcripts receive no probes.
- Captured transcripts shorter than `probe_length` cannot receive probes and are treated as ineligible.
- For each captured transcript, the maximum number of non-overlapping probes is `floor(transcript_length / probe_length)`.
- `probe_density` scales that maximum and rounds down. If density is positive, at least one probe is emitted for each captured eligible transcript.
- Probes are non-overlapping, equally spaced, and centered by distributing slack as evenly as possible around and between probes.
- Different isoforms can still generate probes that overlap after genomic projection, especially when exons are shared. The simulator handles this by using the best single probe overlap for each fragment; overlapping or duplicate probes do not stack extra capture strength.

Examples:

- Transcript length `1200`, probe length `120`, density `1.0`: 10 probes.
- Transcript length `1319`, probe length `120`, density `1.0`: 10 probes, centered with slack distributed as gaps.
- Transcript length `1200`, probe length `120`, density `0.5`: 5 probes.

To simulate both uncaptured and captured reads for the same synthetic reference,
use `capture.configs`. Each entry becomes another condition-grid dimension:

```yaml
capture:
  probe_length: 120
  probe_density: 1.0
  off_target_weight: 1.0
  binding_per_base: 10.0
  gdna_split_penalty: 0.2
  min_overlap: 1
  configs:
    - label: off
      enabled: false
    - label: on
      fraction: 0.5
```

This creates condition directories such as:

```text
gdna_none_ss_0.99_nrna_none_capture_off/
gdna_none_ss_0.99_nrna_none_capture_on/
```

Capture configs are paired for head-to-head comparisons: the transcript
abundance draw and expressed/unexpressed transcript set are shared across all
capture labels in the run. For the random mini-genome suite, capture labels also
share the same simulation seed within each fixed nRNA/gDNA/strand condition;
only the capture weighting changes between `capture_off` and `capture_on`.

You can add more capture configs by adding more entries. Values inside an entry
override the defaults above it:

```yaml
capture:
  probe_length: 120
  off_target_weight: 1.0
  configs:
    - label: off
      enabled: false
    - label: weak
      fraction: 0.5
      binding_per_base: 3.0
    - label: strong
      fraction: 0.5
      binding_per_base: 15.0
    - label: dense
      fraction: 0.5
      probe_density: 1.0
      binding_per_base: 10.0
```

Use short, path-safe labels because they become part of the condition directory
name and the manifest `capture_label` field.

Checked-in generated-capture suite configs:

- `scripts/sim/configs/hybrid_capture_500kb.yaml`: 500 kb genome, 10 genes, 1-5 isoforms per gene, capture off/on with 50% transcript capture in the `on` config.
- `scripts/sim/configs/hybrid_capture_5mb.yaml`: 5 Mb genome, 100 genes, 1-5 isoforms per gene, capture off/on with 50% transcript capture in the `on` config.

## Existing FASTA/GTF Workflow

Use `simulate_reads.py` when you already have a reference and annotation. The command reads a YAML config:

```bash
python scripts/sim/simulate_reads.py \
  --config scripts/sim/configs/sim_example.yaml
```

Minimal config shape:

```yaml
genome: /path/to/genome.fa
gtf: /path/to/genes.gtf
outdir: /scratch/$USER/rigel_sim/example
transcript_filter: all

simulation:
  n_rna_fragments: 1000000
  sim_seed: 42
  frag_mean: 250
  frag_std: 50
  frag_min: 80
  frag_max: 800
  read_length: 150
  error_rate: 0.0
  n_workers: 1

abundance:
  mode: random
  seed: 42
  min: 1.0
  max: 10000.0
  frac_expressed: 0.5

gdna:
  rates: [0.0, 0.1, 0.5]
  rate_labels: [none, low, med]
  frag_mean: 350
  frag_std: 100
  frag_min: 100
  frag_max: 1000

nrna:
  mode: additive_ratio
  ratios: [0.0]
  ratio_labels: [none]

strand_specificities: [0.99, 0.50]
oracle_bam: true
verbose: true
```

To use a single custom hybrid-capture probe file in this workflow, add a
`capture` block:

```yaml
capture:
  probes: /path/to/capture_probes.tsv
  probe_format: transcript
  off_target_weight: 1.0
  binding_per_base: 10.0
  gdna_split_penalty: 0.2
  min_overlap: 1
```

To compare capture off versus one or more capture probe sets in the same run,
use `capture.configs`:

```yaml
capture:
  probe_format: transcript
  off_target_weight: 1.0
  binding_per_base: 10.0
  gdna_split_penalty: 0.2
  min_overlap: 1
  configs:
    - label: off
      enabled: false
    - label: panel_a
      probes: /path/to/panel_a_probes.tsv
    - label: panel_b
      probes: /path/to/panel_b_probes.bed
      probe_format: bed12
      binding_per_base: 15.0
```

Each capture config is crossed with every nRNA, gDNA, and strand-specificity
setting. The labels become condition suffixes, for example
`gdna_high_ss_0.99_nrna_none_capture_panel_a`. Defaults under `capture:` apply to
all entries unless an entry overrides them.

Supported probe formats:

- `transcript`: tab-separated `transcript_id`, `start`, `end` columns in 0-based half-open transcript coordinates.
- `bed12`: genomic BED12 probe intervals. Split probes are projected onto transcript, pre-mRNA, and gDNA coordinate spaces.
- `auto`: detect transcript TSV headers or BED12-looking rows.

## Truth Files And Read Names

The simulator writes truth abundance tables per nRNA condition:

```text
truth_abundances_nrna_<label>.tsv
```

FASTQ read names encode origin, transcript, fragment interval, strand, and index. Oracle BAM records preserve the same origin names and add perfect alignments with `NH=1` tags. These truth names are consumed by simulator analysis utilities and tests.

## Analysis Helpers

After generating a suite, use the analysis script to evaluate outputs:

```bash
python scripts/sim/evaluate_suite.py \
  --sim-dir /scratch/$USER/rigel_sim/synthetic_smoke
```

For large benchmark runs that include Rigel, Salmon, STAR, or cluster orchestration, use the dedicated benchmarking framework under `scripts/benchmarking/`.
