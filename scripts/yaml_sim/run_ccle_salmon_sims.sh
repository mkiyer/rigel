#!/usr/bin/env bash
# Run all CCLE salmon simulations sequentially.
# Usage: bash scripts/yaml_sim/run_ccle_salmon_sims.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_DIR="${SCRIPT_DIR}/ccle_salmon_libs"
SIM_SCRIPT="${SCRIPT_DIR}/../sim.py"

TOTAL=$(ls "${CONFIG_DIR}"/sim_ccle_*.yaml 2>/dev/null | wc -l)
echo "=== CCLE Salmon Simulation Batch ==="
echo "    Configs: ${TOTAL}"
echo "    Start:   $(date)"
echo ""

COUNT=0
for cfg in "${CONFIG_DIR}"/sim_ccle_*.yaml; do
    COUNT=$((COUNT + 1))
    lib=$(basename "$cfg" .yaml | sed 's/^sim_//')
    echo "[$COUNT/$TOTAL] Running: $lib ($(date +%H:%M:%S))"
    python "$SIM_SCRIPT" --config "$cfg" 2>&1 | tail -5
    echo ""
done

echo "=== All simulations complete: $(date) ==="
