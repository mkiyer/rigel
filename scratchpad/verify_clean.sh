#!/bin/bash
# Behaviour-identity gate for a cleanup step.
#   Usage: bash scratchpad/verify_clean.sh <armname>            (compares against $BASE, default head0)
#          BASE=myhead bash scratchpad/verify_clean.sh <armname>
#
# ⚠⚠ TWO WAYS THIS GATE USED TO LIE, both fixed here and both worth knowing about:
#
#  1. **An empty arm scored 32/32.** The old comparison looped over the NEW arm's conditions, so a bench
#     that died (it did — OOM, under a concurrent genome-scale job) produced no rows, the loop found no
#     differences, and it printed "32/32 IDENTICAL". Now both arms must carry 32 rows or the gate fails.
#
#  2. **The baseline went stale.** `clean0` was recorded in an earlier session and by 2026-07-27 UNMODIFIED
#     HEAD no longer reproduced it — including the `mass` column, which comes from the accumulator payload
#     and no solver change can touch. The `_selfsolve_cache` had been regenerated underneath it. A stale
#     baseline reports 0/32 for a change that is in fact bit-identical, which is the expensive direction of
#     wrong. **Re-record the baseline from an actual `git stash` of HEAD in the SAME session, at BOTH
#     refits, before trusting a single comparison** — and if a HEAD-vs-baseline run is not 32/32, the
#     baseline is the thing that is broken, not your change.
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
cd /Users/mkiyer/proj/rigel
BASE="${BASE:-head0}"
ruff check src/ tests/ scripts/
python -m pytest tests/ -q 2>&1 | tail -1
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm "$1_r0" >/tmp/verify_clean_r0.log 2>&1
P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm "$1_r1" >/tmp/verify_clean_r1.log 2>&1
python - "$1" "$BASE" <<'PY'
import csv, collections, sys
a, base_arm = sys.argv[1], sys.argv[2]
r = collections.defaultdict(dict)
for x in csv.DictReader(open('/tmp/pass0_oracle_bench.tsv'), delimiter='\t'):
    r[x['arm']][x['cond']] = {k: v for k, v in x.items() if k != 'arm'}
bad = False
for suf in ('_r0', '_r1'):
    new, base = r[a + suf], r[base_arm + suf]
    if len(new) != 32 or len(base) != 32:
        print(f"  ⛔ {a+suf}: {len(new)} rows vs {base_arm+suf}: {len(base)} rows — EXPECTED 32/32. GATE INVALID.")
        bad = True
        continue
    d = [k for k in new if new[k] != base[k]]
    print(f"  {a+suf} vs {base_arm+suf}: {32-len(d)}/32 IDENTICAL" + ("" if not d else f"  ⚠ DIFFERS on {d[:3]}"))
    bad = bad or bool(d)
sys.exit(1 if bad else 0)
PY
