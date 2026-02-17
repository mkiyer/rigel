#!/usr/bin/env bash
# run_benchmarks.sh — Multi-seed regional benchmark runner
#
# Runs benchmark_region_competition.py across multiple seeds for
# reproducibility and variance estimation.
#
# Usage:
#   cd /Users/mkiyer/proj/hulkrna
#   bash scripts/run_benchmarks.sh [outdir] [n_fragments] [seeds...]
#
# Defaults:
#   outdir:       /Users/mkiyer/Downloads/hulkrna_runs/bench_zero_gdna
#   n_fragments:  50000
#   seeds:        101 202 303 404 505
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

# Reference files
GENOME="/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genome_controls.fasta.bgz"
GTF="/Users/mkiyer/Downloads/hulkrna_runs/refs/human/genes_controls.gtf.gz"
REGIONS="$SCRIPT_DIR/regions_benchmark.tsv"

# Parse args
OUTDIR="${1:-/Users/mkiyer/Downloads/hulkrna_runs/bench_zero_gdna}"
N_FRAGMENTS="${2:-50000}"
shift 2 2>/dev/null || true
SEEDS=("${@:-101 202 303 404 505}")
if [[ ${#SEEDS[@]} -eq 0 ]] || [[ "${SEEDS[0]}" == "" ]]; then
    SEEDS=(101 202 303 404 505)
fi

# Zero-gDNA, perfect strand specificity configuration
STRAND_SPEC="1.0"
GDNA_ABUNDANCE="0"
GDNA_THRESHOLD="0"
ABUNDANCE_MODE="random"
ABUNDANCE_MIN="0.01"
ABUNDANCE_MAX="100000"
THREADS="4"

echo "======================================================"
echo "hulkrna Regional Benchmark Runner"
echo "======================================================"
echo "Output:       $OUTDIR"
echo "Fragments:    $N_FRAGMENTS"
echo "Seeds:        ${SEEDS[*]}"
echo "Strand spec:  $STRAND_SPEC"
echo "gDNA:         $GDNA_ABUNDANCE"
echo "Regions file: $REGIONS"
echo "======================================================"

cd "$PROJECT_DIR"

for seed in "${SEEDS[@]}"; do
    OUT_SEED="$OUTDIR/seed_${seed}"
    echo ""
    echo "--- Seed $seed → $OUT_SEED ---"
    
    PYTHONPATH=src conda run -n hulkrna python scripts/benchmark_region_competition.py \
        --genome "$GENOME" \
        --gtf "$GTF" \
        --region-file "$REGIONS" \
        --outdir "$OUT_SEED" \
        --n-fragments "$N_FRAGMENTS" \
        --sim-seed "$seed" \
        --pipeline-seed "$seed" \
        --abundance-seed "$seed" \
        --strand-specificity "$STRAND_SPEC" \
        --gdna-abundance "$GDNA_ABUNDANCE" \
        --gdna-threshold "$GDNA_THRESHOLD" \
        --abundance-mode "$ABUNDANCE_MODE" \
        --abundance-min "$ABUNDANCE_MIN" \
        --abundance-max "$ABUNDANCE_MAX" \
        --threads "$THREADS" \
        --keep-going \
        --verbose
    
    echo "--- Seed $seed complete ---"
done

echo ""
echo "======================================================"
echo "All seeds complete. Aggregating results..."
echo "======================================================"

# Aggregate across seeds
PYTHONPATH=src conda run -n hulkrna python scripts/aggregate_benchmarks.py \
    --input-dir "$OUTDIR" \
    --output-dir "$OUTDIR"

echo "Done! Results in $OUTDIR"
