#!/bin/bash
# Full pristine benchmark pipeline — run with nohup
# Usage: nohup bash scripts/benchmark/run_pristine_pipeline.sh > results/vcap_pristine_pipeline.log 2>&1 &
set -euo pipefail
cd /home/mkiyer/proj/rigel

# Activate conda
eval "$(conda shell.bash hook)"
conda activate rigel

CFG=scripts/benchmarking/configs/vcap_pristine.yaml
OUTDIR=results/vcap_pristine

echo "============================================"
echo "PRISTINE BENCHMARK PIPELINE"
echo "Started: $(date)"
echo "Config: $CFG"
echo "============================================"

# Step 1: Align with minimap2
echo ""
echo "[Step 1/5] Aligning with minimap2..."
echo "  Started: $(date)"
python -m scripts.benchmarking align -c $CFG -v 2>&1
echo "  Completed: $(date)"

# Step 2: Run salmon + kallisto
echo ""
echo "[Step 2/5] Running salmon and kallisto..."
echo "  Started: $(date)"
python -m scripts.benchmarking run-tools -c $CFG -v 2>&1
echo "  Completed: $(date)"

# Step 3: Run rigel (VBEM + MAP)
echo ""
echo "[Step 3/5] Running rigel quant (VBEM + MAP)..."
echo "  Started: $(date)"
python -m scripts.benchmarking run -c $CFG -v 2>&1
echo "  Completed: $(date)"

# Step 4: Analyze
echo ""
echo "[Step 4/5] Analyzing results..."
echo "  Started: $(date)"
python -m scripts.benchmarking analyze -c $CFG -o $OUTDIR -v 2>&1
echo "  Completed: $(date)"

# Step 5: List outputs
echo ""
echo "[Step 5/5] Listing outputs..."
ls -lhR $OUTDIR/ 2>/dev/null || echo "No output directory found"
echo ""
echo "============================================"
echo "PIPELINE COMPLETE"
echo "Finished: $(date)"
echo "============================================"
