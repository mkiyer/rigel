# Synthetic 24-Condition Simulation Suite

Status: ready to run, 2026-05-12.

This suite generates a 10 Mb synthetic genome, approximately 250 transcripts,
1 million mature RNA read pairs per condition, four gDNA levels, three sparse
nRNA levels, and two strand modes. Output is written under:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24
```

## Configuration

The reusable read-generation config is:

```text
scripts/sim/configs/synthetic_24_random_nrna.yaml
```

It points to the synthetic reference files that are generated in the suite
directory:

```yaml
genome: /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24/reference/genome.fa
gtf: /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24/reference/genes.gtf
outdir: /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24
```

The mature RNA depth is fixed at 1 million read pairs in every condition:

```yaml
simulation:
  n_rna_fragments: 1_000_000
  frag_mean: 250
  frag_std: 50
  frag_min: 80
  frag_max: 800
  read_length: 150
```

Transcript expression is sampled once from a log-uniform model with 50%
expression probability:

```yaml
abundance:
  mode: random
  min: 1.0
  max: 10000.0
  frac_expressed: 0.5
```

The condition grid is 24 scenarios:

```text
4 gDNA levels x 3 nRNA levels x 2 strand modes = 24 conditions
```

gDNA is additive relative to mature RNA depth:

```yaml
gdna:
  rates: [0.0, 0.25, 1.0, 2.0]
  rate_labels: [zero, low, med, high]
```

nRNA uses sparse random per-transcript fractions. For each nonzero nRNA
scenario, 50% of expressed multi-exon transcripts are selected for nRNA, and
each selected transcript gets a random nRNA:mRNA molecular abundance fraction:

```yaml
nrna:
  mode: random_fraction
  ratio_ranges:
    - [0.0, 0.0]      # zero
    - [0.05, 0.25]    # low
    - [0.50, 1.00]    # high
  ratio_labels: [zero, low, high]
  eligible_fraction: 0.5
```

The realized total nRNA:mRNA abundance ratio is computed from the generated
truth table and recorded in `manifest.json` as `nrna_ratio`. The generated nRNA
fragment count is `round(1_000_000 * realized_ratio)`, so fragment truth follows
the sampled molecular truth rather than a fixed closed-total fraction.

## Generate The Suite

Run from the repository root:

```bash
conda activate rigel
cd /Users/mkiyer/proj/rigel
```

Create the synthetic 10 Mb reference with about 50 genes and about 250
transcripts:

```bash
python scripts/sim/simulate_suite.py \
  --outdir /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24 \
  --reference-only \
  --genome-length 10000000 \
  --n-genes 50 \
  --seed 20260512
```

Generate all 24 read/oracle-BAM conditions:

```bash
python scripts/sim/simulate_reads.py \
  --config scripts/sim/configs/synthetic_24_random_nrna.yaml \
  -j 4
```

The condition names have this form:

```text
gdna_<zero|low|med|high>_ss_<0.50|0.99>_nrna_<zero|low|high>
```

Examples:

```text
gdna_zero_ss_0.50_nrna_zero
gdna_low_ss_0.99_nrna_high
gdna_high_ss_0.99_nrna_low
```

## Run Rigel On The Suite

Run quantification and write the full evaluation report:

```bash
python scripts/sim/evaluate_suite.py \
  --sim-base /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24
```

The evaluator builds a Rigel index under the suite directory, runs `rigel quant`
for each manifest condition, and writes:

```text
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24/analysis_report.txt
```

For a quick metadata/report wiring check without running quantification or slow
fragment-level analysis:

```bash
python scripts/sim/evaluate_suite.py \
  --sim-base /Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24 \
  --skip-quant \
  --skip-frag-analysis
```

## Expected Outputs

The suite directory contains:

```text
reference/genome.fa
reference/genes.gtf
manifest.json
truth_abundances_nrna_zero.tsv
truth_abundances_nrna_low.tsv
truth_abundances_nrna_high.tsv
gdna_*_ss_*_nrna_*/sim_R1.fq.gz
gdna_*_ss_*_nrna_*/sim_R2.fq.gz
gdna_*_ss_*_nrna_*/sim_oracle.bam
```

Important manifest fields per condition:

```text
n_mrna   = 1,000,000 mature RNA read pairs
n_nrna   = realized sparse nRNA read pairs
n_gdna   = round(gdna_rate * 1,000,000)
n_total  = n_mrna + n_nrna + n_gdna
```