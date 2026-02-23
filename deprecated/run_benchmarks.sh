#!/usr/bin/env bash
# run_benchmarks.sh — Multi-seed combinatorial benchmark runner
#
# Runs benchmark_region_competition.py across multiple seeds with a
# full combinatorial grid of gDNA levels and strand specificities.
#
# Tools benchmarked:
#   hulkrna_mm       (include_multimap=True)
#   salmon
#   kallisto
#   htseq            (gene-level, enabled by default)
#
# Usage:
#   cd /Users/mkiyer/proj/hulkrna
#   bash scripts/benchmarking/run_benchmarks.sh [outdir] [n_fragments] [seeds...]
#
# Defaults:
#   outdir:       /Users/mkiyer/Downloads/hulkrna_runs/bench_combinatorial
#   n_fragments:  50000
#   seeds:        101 202 303 404 505
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Reference files
GENOME="/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genome_controls.fasta.bgz"
GTF="/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz"
REGIONS="$SCRIPT_DIR/regions_benchmark.tsv"

# Parse args
OUTDIR="${1:-/Users/mkiyer/Downloads/hulkrna_runs/bench_combinatorial}"
N_FRAGMENTS="${2:-50000}"
shift 2 2>/dev/null || true
SEEDS=("${@:-101 202 303 404 505}")
if [[ ${#SEEDS[@]} -eq 0 ]] || [[ "${SEEDS[0]}" == "" ]]; then
    SEEDS=(101 202 303 404 505)
fi

# ── Combinatorial grid ──────────────────────────────────────────────
# gDNA rates: fraction of total transcript abundance assigned to gDNA
GDNA_RATES="0.0,0.1,0.25,0.5"
GDNA_RATE_LABELS="none,low,moderate,high"

# nRNA rates: fraction of transcript abundance used as nascent abundance
NRNA_RATES="0.0,0.25"
NRNA_RATE_LABELS="none,quarter"

# Strand specificity values
STRAND_SPECS="0.95,0.99,1.0"

# HTSeq integration (set to 1 to enable; requires 'htseq' conda env)
INCLUDE_HTSEQ="${INCLUDE_HTSEQ:-1}"
HTSEQ_CONDA_ENV="${HTSEQ_CONDA_ENV:-htseq}"

# Other parameters
ABUNDANCE_MODE="random"
ABUNDANCE_MIN="0.01"
ABUNDANCE_MAX="100000"
THREADS="4"
ALIGNER="${ALIGNER:-minimap2}"

echo "======================================================"
echo "hulkrna Combinatorial Benchmark Runner"
echo "======================================================"
echo "Output:       $OUTDIR"
echo "Fragments:    $N_FRAGMENTS"
echo "Seeds:        ${SEEDS[*]}"
echo "gDNA rates:   $GDNA_RATES"
echo "nRNA rates:   $NRNA_RATES"
echo "Strand specs: $STRAND_SPECS"
echo "Aligner:      $ALIGNER"
echo "HTSeq:        $([ "$INCLUDE_HTSEQ" = "1" ] && echo "enabled (env: $HTSEQ_CONDA_ENV)" || echo "disabled")"
echo "Regions file: $REGIONS"
echo "======================================================"

cd "$PROJECT_DIR"

HTSEQ_FLAG="--htseq-conda-env $HTSEQ_CONDA_ENV"
if [[ "$INCLUDE_HTSEQ" != "1" ]]; then
    HTSEQ_FLAG="$HTSEQ_FLAG --no-htseq"
fi

for seed in "${SEEDS[@]}"; do
    OUT_SEED="$OUTDIR/seed_${seed}"
    echo ""
    echo "--- Seed $seed → $OUT_SEED ---"

    PYTHONPATH=src conda run -n hulkrna python scripts/benchmarking/benchmark_region_competition.py \
        --genome "$GENOME" \
        --gtf "$GTF" \
        --region-file "$REGIONS" \
        --outdir "$OUT_SEED" \
        --n-fragments "$N_FRAGMENTS" \
        --sim-seed "$seed" \
        --pipeline-seed "$seed" \
        --abundance-seed "$seed" \
        --gdna-rates "$GDNA_RATES" \
        --gdna-rate-labels "$GDNA_RATE_LABELS" \
        --nrna-rates "$NRNA_RATES" \
        --nrna-rate-labels "$NRNA_RATE_LABELS" \
        --strand-specificities "$STRAND_SPECS" \
        --abundance-mode "$ABUNDANCE_MODE" \
        --abundance-min "$ABUNDANCE_MIN" \
        --abundance-max "$ABUNDANCE_MAX" \
        --threads "$THREADS" \
        --aligner "$ALIGNER" \
        $HTSEQ_FLAG \
        --keep-going \
        --verbose

    echo "--- Seed $seed complete ---"
done

echo ""
echo "======================================================"
echo "All seeds complete. Aggregating results..."
echo "======================================================"

# Aggregate across seeds
PYTHONPATH=src conda run -n hulkrna python scripts/benchmarking/aggregate_benchmarks.py \
    --input-dir "$OUTDIR" \
    --output-dir "$OUTDIR"

echo "Done! Results in $OUTDIR"
