#!/bin/bash
# Re-run rigel quant on vcap mixture titration with Option-B gDNA length correction.
#
# For each library, input = existing rigel/annotated.bam (name-sorted, NH-tagged).
# Output goes to /scratch/.../runs/human_optionb/<lib>/.
#
# The base rigel/quant.feather from the old pre-Option-B run remains intact for
# comparison.

set -euo pipefail

BASE=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human
OUT=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human_optionb
INDEX=/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index
THREADS=10

LIBS=(
  mctp_vcap_rna20m_dna00m
  mctp_vcap_rna20m_dna01m
  mctp_vcap_rna20m_dna02m
  mctp_vcap_rna20m_dna05m
  mctp_vcap_rna20m_dna10m
  mctp_vcap_rna20m_dna20m
  mctp_vcap_rna20m_dna40m
  mctp_vcap_rna20m_dna80m
)

source /home/mkiyer/sw/miniforge3/etc/profile.d/conda.sh
conda activate rigel

mkdir -p "$OUT"

for lib in "${LIBS[@]}"; do
  bam="$BASE/$lib/rigel/annotated.bam"
  outdir="$OUT/$lib"
  mkdir -p "$outdir"
  if [[ -f "$outdir/quant.feather" && -f "$outdir/summary.json" ]]; then
    echo "[skip] $lib (already has output)"
    continue
  fi
  echo "=== $lib ==="
  /usr/bin/time -v rigel quant \
    --bam "$bam" \
    --index "$INDEX" \
    --output-dir "$outdir" \
    --threads "$THREADS" \
    --buffer-size 8 \
    2>&1 | tee "$outdir/run.log"
  echo "[done] $lib"
done
