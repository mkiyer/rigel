#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

usage() {
  cat <<'EOF'
Usage:
  bash scripts/benchmarking/launch_benchmark.sh <config.yaml> [options]

Runs benchmark_region_competition.py and then aggregate_benchmarks.py
against the benchmark outdir (from config, or overridden with --outdir).

Options:
  -e, --env <name>       Conda environment to run in (default: hulkrna)
  -o, --outdir <path>    Override output directory for this run
      --skip-aggregate   Run benchmark only (do not aggregate)
  -h, --help             Show this help message
EOF
}

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

CONFIG=""
CONDA_ENV="hulkrna"
OUTDIR_OVERRIDE=""
SKIP_AGGREGATE=0

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  usage
  exit 0
fi

CONFIG="$1"
shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    -e|--env)
      CONDA_ENV="$2"
      shift 2
      ;;
    -o|--outdir)
      OUTDIR_OVERRIDE="$2"
      shift 2
      ;;
    --skip-aggregate)
      SKIP_AGGREGATE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ ! -f "$CONFIG" ]]; then
  echo "Config file not found: $CONFIG" >&2
  exit 1
fi

cd "$PROJECT_DIR"

OUTDIR="$OUTDIR_OVERRIDE"
if [[ -z "$OUTDIR" ]]; then
  OUTDIR=$(PYTHONPATH=src conda run -n "$CONDA_ENV" python -c "import yaml,sys; cfg=yaml.safe_load(open(sys.argv[1])) or {}; out=cfg.get('outdir'); print(out if out else '')" "$CONFIG")
fi

if [[ -z "$OUTDIR" ]]; then
  echo "Could not determine outdir. Set it in config or pass --outdir." >&2
  exit 1
fi

echo "======================================================"
echo "hulkrna benchmark launcher"
echo "======================================================"
echo "Config:      $CONFIG"
echo "Conda env:   $CONDA_ENV"
echo "Output dir:  $OUTDIR"
echo "Aggregate:   $([[ $SKIP_AGGREGATE -eq 1 ]] && echo no || echo yes)"
echo "======================================================"

BENCH_CMD=(
  python scripts/benchmarking/benchmark_region_competition.py
  --config "$CONFIG"
)
if [[ -n "$OUTDIR_OVERRIDE" ]]; then
  BENCH_CMD+=(--outdir "$OUTDIR_OVERRIDE")
fi

PYTHONPATH=src conda run -n "$CONDA_ENV" "${BENCH_CMD[@]}"

if [[ $SKIP_AGGREGATE -eq 0 ]]; then
  PYTHONPATH=src conda run -n "$CONDA_ENV" python scripts/benchmarking/aggregate_benchmarks.py \
    --input-dir "$OUTDIR" \
    --output-dir "$OUTDIR"
fi

echo "Done. Results in: $OUTDIR"
