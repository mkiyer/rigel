#!/usr/bin/env bash
# Run all CCLE salmon simulations in parallel (8 jobs).
# Usage: conda activate rigel && bash scripts/yaml_sim/run_ccle_salmon_parallel.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/ccle_salmon_libs"
SIM_SCRIPT="${SCRIPT_DIR}/../sim.py"
JOBS=8

CONFIGS=("${CONFIG_DIR}"/sim_ccle_*.yaml)
TOTAL=${#CONFIGS[@]}

echo "=== CCLE Salmon Simulation Batch (parallel=$JOBS) ==="
echo "    Configs: ${TOTAL}"
echo "    Start:   $(date)"
echo ""

run_one() {
    cfg="$1"
    lib=$(basename "$cfg" .yaml | sed 's/^sim_//')
    echo "[START] $lib ($(date +%H:%M:%S))"
    python "$SIM_SCRIPT" --config "$cfg" > "/tmp/sim_${lib}.log" 2>&1
    status=$?
    if [ $status -eq 0 ]; then
        echo "[DONE]  $lib ($(date +%H:%M:%S))"
    else
        echo "[FAIL]  $lib (exit=$status) — see /tmp/sim_${lib}.log"
    fi
    return $status
}
export -f run_one
export SIM_SCRIPT

parallel --jobs "$JOBS" --line-buffer \
    run_one {} ::: "${CONFIGS[@]}"

echo ""
echo "=== All simulations complete: $(date) ==="
