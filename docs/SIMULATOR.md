# Rigel Simulator Manual

Rigel ships a synthetic RNA-seq read simulator used to generate the ground-truth test scenarios
and calibration benchmarks the tool is validated against. It is built around **one** read-generation
engine — `WholeGenomeSimulator` (`rigel/sim/wgs_engine.py`): vectorized (numpy), parallel
(fork-based sharding), FASTA-backed (pysam), generating paired-end FASTQ plus an optional
name-sorted **oracle BAM** (perfect alignments), with mRNA / nascent-RNA (nRNA) / gDNA pools,
strand-specificity, hybrid-capture weighting, and gDNA strand overdispersion.

Three user-facing entry points sit on that one engine:

| Task | Tool | Reference |
|---|---|---|
| Hand-author a **small scenario** (explicit genes/exons) on a synthetic genome | `rigel sim` (CLI) or the `Scenario` API (Python) | §1 |
| Generate a **random mini-genome suite** sweeping conditions (the benchmark) | `scripts/sim/simulate_suite.py` | §2 |
| Simulate from an **existing FASTA + GTF** reference | `scripts/sim/simulate_reads.py` | §3 |

Use a synthetic genome (`rigel sim` / `simulate_suite`) for controlled, hard test cases with known
truth; use the existing-reference workflow (`simulate_reads`) when you already have a FASTA, GTF,
and optional probe panel.

> **FASTQ vs oracle BAM.** Every workflow can emit two artifacts. **FASTQ** is meant to be aligned
> with a real aligner (minimap2) — it exercises Rigel end-to-end including aligner noise. The
> **oracle BAM** is a perfect, name-sorted alignment written directly from the known truth
> (`NH=1`, correct CIGAR/strand/`XS`), which isolates calibration/quant from the aligner. Read
> names encode ground-truth origin in both (§5).

## Environment

All commands run inside the activated `rigel` conda environment, from the repo root:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
```

Write outputs to a scratch dir outside the repo — simulations grow large quickly.

---

## 1. Single scenario — `rigel sim` / the `Scenario` API

For a targeted scenario with **explicitly defined** genes/transcripts on a small synthetic genome
(the workhorse for unit/integration tests).

### `rigel sim` (CLI)

```bash
rigel sim --config scenario.yaml --output-dir /scratch/$USER/rigel_sim/scn
```

The YAML defines the genome + genes:

```yaml
name: two_genes
genome_length: 10000
seed: 42
ref_name: chr1
# fragment model (optional; ReadSimConfig defaults shown)
frag_mean: 250
frag_std: 50
frag_min: 50
frag_max: 800
read_length: 150
error_rate: 0.0
n_fragments: 5000
genes:
  - gene_id: g1
    strand: "+"
    transcripts:
      - {t_id: t1, exons: [[100, 300], [500, 700]], abundance: 100}
  - gene_id: g2
    strand: "-"
    transcripts:
      - {t_id: t2, exons: [[2000, 2400]], abundance: 50}
```

It writes the genome FASTA (+ `.fai`), a GTF, paired FASTQ, an aligned BAM (FASTQ → minimap2), and
a built `TranscriptIndex` under `--output-dir`.

### `Scenario` API (Python)

The same capability programmatically — used throughout the test suite. `build()` aligns with
minimap2; `build_oracle()` writes a perfect oracle BAM (no aligner):

```python
from rigel.sim import Scenario, GDNAConfig

with Scenario("demo", genome_length=10_000, seed=42) as sc:
    sc.add_gene("g1", "+", [
        {"t_id": "t1", "exons": [(100, 300), (500, 700)], "abundance": 100},
    ])
    # realistic: FASTQ → minimap2 → name-sorted BAM + index
    result = sc.build(n_fragments=5_000)
    # or perfect alignments (isolate quant from the aligner), with gDNA contamination:
    result = sc.build_oracle(n_fragments=5_000, gdna_config=GDNAConfig(abundance=10))

# result.{fasta_path, gtf_path, bam_path, index, transcripts, ...}
```

`Scenario` builds the genome (`MutableGenome`) and annotation (`GeneBuilder`), then drives the one
engine. `GDNAConfig` / `ReadSimConfig` (`rigel/sim/reads.py`) are the single-condition config
surface; `Scenario` translates them to the engine internally.

---

## 2. Random mini-genome suite — `simulate_suite.py`

Generates a random synthetic reference (genes, isoforms, antisense overlaps) **and** simulates a
grid of conditions: nRNA × gDNA-rate × gDNA-strand-overdispersion × strand-specificity × capture.
This is what the calibration benchmark uses.

```bash
# Quick smoke suite
python scripts/sim/simulate_suite.py \
  --outdir /scratch/$USER/rigel_sim/smoke --profile smoke --n-rna 100000

# Config-driven (CLI flags override the config)
python scripts/sim/simulate_suite.py \
  --config scripts/sim/configs/hybrid_capture_500kb.yaml \
  --outdir /scratch/$USER/rigel_sim/hybrid_capture_500kb
```

### Output layout

```text
<outdir>/
├── reference/
│   ├── genome.fa  +  genome.fa.fai
│   ├── genes.gtf
│   └── capture_probes_<label>.tsv|.bed   # only when a capture sweep generates probes
├── truth_abundances_nrna_<label>.tsv     # molecular truth per nRNA condition
├── manifest.json                          # the condition index (params + truth/output paths + observed counts)
└── gdna_<label>[_od_<od>]_ss_<value>_nrna_<label>[_capture_<label>]/
    ├── sim_R1.fq.gz  +  sim_R2.fq.gz
    └── sim_oracle.bam                      # when oracle_bam: true
```

The `_od_<value>` slug appears only when the overdispersion sweep is used (and never for `od00`),
and `_capture_<label>` only when a capture sweep is configured — so baseline configs keep their
original condition names.

### Condition grid

`--profile` sets the nRNA/gDNA/strand defaults:

- `smoke`: nRNA `none` only; gDNA rates `[0, 0.1, 0.5, 1, 2]`; strand-specificity `[0.99, 0.50]`.
- `full`: nRNA `[none, low, med, high]`; same gDNA/strand axes.

Override any axis from the CLI or config:

```bash
python scripts/sim/simulate_suite.py --outdir … --profile smoke \
  --gdna-rates 0.0,1.0 --gdna-labels none,high \
  --strand-specificities 0.99,0.50 \
  --gdna-strand-overdispersions 0.0,0.14 --gdna-strand-overdispersion-labels od00,od14
```

**gDNA strand overdispersion** (`--gdna-strand-overdispersions`) is the intra-class correlation of
the per-region sense/antisense split around ½: `0.0` = exact 50/50 (Binomial); larger = more
region-to-region skew (`0.14 ≈ 1/7`, the calibration prior level). Each value becomes its own
condition, so the suite can stress the gDNA strand model.

### Reference + run options

- Reference generation: `--genome-length`, `--n-genes`, `--min-isoforms`, `--max-isoforms`,
  `--target-transcripts` (`0` ⇒ uniform isoform counts in `[min, max]`), `--antisense-fraction`.
- Fragments / perf: `--n-rna` (RNA fragments per condition), `--n-workers` (parallel sharding),
  `--seed`.
- `--reference-only`: write just `reference/` (genome + GTF + probes), no conditions.
- `--skip-existing`: only fill in missing condition outputs.
- `--conditions <name> …`: generate only the named condition directories.

---

## 3. Existing FASTA/GTF — `simulate_reads.py`

When you already have a reference + annotation, simulate a condition grid from it via a YAML config:

```bash
python scripts/sim/simulate_reads.py --config scripts/sim/configs/sim_example.yaml
```

```yaml
genome: /path/to/genome.fa
gtf: /path/to/genes.gtf
outdir: /scratch/$USER/rigel_sim/example
transcript_filter: all          # all | basic | mane | ccds

simulation:
  n_rna_fragments: 1000000
  sim_seed: 42
  frag_mean: 250
  frag_std: 50
  frag_min: 80
  frag_max: 800
  read_length: 150
  error_rate: 0.0
  n_workers: 4

abundance:
  mode: random                  # random | file
  seed: 42
  min: 1.0
  max: 10000.0
  frac_expressed: 0.5
  # file: /path/to/abundances.tsv   # for mode: file (salmon/kallisto/sim TSV)

gdna:
  rates: [0.0, 0.1, 0.5]
  rate_labels: [none, low, med]
  frag_mean: 350
  frag_std: 100
  frag_min: 100
  frag_max: 1000
  strand_overdispersions: [0.0, 0.14]      # optional sweep axis
  strand_overdispersion_labels: [od00, od14]

nrna:
  mode: additive_ratio          # additive_ratio | random_fraction | file
  ratios: [0.0]
  ratio_labels: [none]

strand_specificities: [0.99, 0.50]
oracle_bam: true
verbose: true
```

Every `gdna` rate × `gdna` overdispersion × `nrna` setting × `strand_specificities` × capture
config becomes one condition (same output layout + manifest as §2).

---

## 4. Hybrid capture

Both the suite and the existing-reference workflow model hybrid-capture enrichment. A fragment
start's weight is `off_target_weight + binding_per_base × best_probe_overlap` (a fragment binds at
most one probe — overlapping/duplicate probes do not stack).

### Generated probes (suite only)

The mini-genome suite can design transcript-coordinate probes from the synthetic transcriptome —
no hand-authored panel needed:

```bash
python scripts/sim/simulate_suite.py --outdir … --profile smoke \
  --capture-fraction 0.25 --probe-length 120 --probe-density 1.0
```

Design parameters: `--capture-fraction` (fraction of eligible transcripts targeted; `0` disables),
`--probe-length` (default 120), `--probe-density` (fraction of the max non-overlapping tiling,
`0..1`), `--capture-off-target-weight` (1.0), `--capture-binding-per-base` (10.0),
`--capture-gdna-split-penalty` (0.2, for probes split by introns when projected to gDNA/pre-mRNA),
`--capture-min-overlap` (1). Probes are non-overlapping, equally spaced, and centered; captured
transcripts shorter than `probe_length` are ineligible. Probe **design** lives in
`rigel/sim/capture/design.py`; the runtime sampler in `rigel/sim/capture/sampler.py`.

### Capture sweep (off vs on, head-to-head)

Add a `capture.configs` list — each entry is another condition-grid dimension, and capture labels
are **paired** (they share the abundance draw, expressed set, and per-condition seed, so only the
capture weighting differs):

```yaml
capture:
  probe_length: 120
  probe_density: 1.0
  off_target_weight: 1.0
  binding_per_base: 10.0
  gdna_split_penalty: 0.2
  min_overlap: 1
  configs:
    - {label: off, enabled: false}
    - {label: on,  fraction: 0.5}        # suite: generate probes for 50% of transcripts
    # existing-reference: give each entry a `probes:` file instead of `fraction:`
    # - {label: panel_a, probes: /path/panel_a.tsv}
    # - {label: panel_b, probes: /path/panel_b.bed, probe_format: bed12, binding_per_base: 15.0}
```

Values inside an entry override the block defaults above it. Use short, path-safe labels (they
become the `_capture_<label>` condition suffix + the manifest `capture_label`).

Supported probe formats (existing-reference workflow): `transcript` (TSV `transcript_id start end`,
0-based half-open transcript coords), `bed12` (genomic BED12, projected onto transcript/pre-mRNA/
gDNA spaces), `auto` (detect).

Checked-in suite configs: `scripts/sim/configs/hybrid_capture_500kb.yaml` (500 kb, 10 genes),
`hybrid_capture_500kb_od.yaml` (+ the overdispersion axis), `hybrid_capture_5mb.yaml` (5 Mb).

---

## 5. Truth files and read names

- **`truth_abundances_nrna_<label>.tsv`** — molecular truth per nRNA condition (transcript
  abundances, `observed_mrna_fragments`, etc.).
- **Per-condition truth** — each condition dir gets `truth_abundances.tsv`,
  `truth_fragment_lengths.tsv`, `truth_summary.json` (post-capture empirical truth + origin counts).
- **`manifest.json`** — the machine-readable index of all conditions: parameters
  (gdna/od/ss/nrna/capture, seed, fragment counts) + output/truth paths + observed origin counts.
- **Read names encode ground truth** (parsed by `rigel/sim/read_name.py`):
  - RNA: `{t_id}:{frag_start}-{frag_end}:{strand}:{index}` (nRNA: `nrna_{t_id}:…`)
  - gDNA: `gdna:[{ref}:]{start}-{end}:{strand}:{index}`

  The oracle BAM preserves these names with `NH=1` and correct CIGAR/strand tags, so truth is
  recoverable directly from any BAM/FASTQ the simulator produced.

---

## 6. Running and evaluating a benchmark

The canonical evaluator builds the index, runs `rigel quant` on each condition's oracle BAM, and
reports pool-level + per-fragment gDNA-vs-RNA classification accuracy against truth:

```bash
SUITE=/scratch/$USER/rigel_sim/hybrid_capture_500kb
python scripts/sim/bench_calibration.py --sim-base "$SUITE" --run --force
# → $SUITE/bench_calibration_metrics.tsv + bench_calibration_report.txt
```

`scripts/sim/evaluate_suite.py --sim-base "$SUITE"` produces the richer per-condition report
(calibration scalars, abundance accuracy, locus gDNA, fragment-assignment confusion). The
`calibration-benchmark` skill documents interpretation. For large multi-tool benchmarks (Salmon,
STAR, cluster orchestration) use `scripts/benchmarking/`.

---

## Package layout (for contributors)

```
rigel/sim/
  genome.py / synthetic_genome.py / annotation.py   genome + transcriptome generation
  splice_motif.py / intervals.py / bam.py / sampling.py / read_name.py   shared primitives
  reads.py                  Scenario-facing config dataclasses (ReadSimConfig / GDNAConfig)
  wgs_config.py             engine/suite config dataclasses
  wgs_engine.py             WholeGenomeSimulator — the one read-generation engine
  capture/{config,sampler,design}.py   hybrid-capture config / runtime / probe design
  truth.py / manifest.py    ground-truth tables + suite manifest
  orchestrator.py           the shared condition-grid loop (one copy; both drivers call it)
  scenario.py               Scenario API (single-condition)
  suite.py / whole_genome.py   the simulate_suite / simulate_reads frontends
  benchmark.py / analysis.py   accuracy benchmarking + per-condition analysis
```

See `docs/sim/architecture_redesign.md` for the design rationale.
