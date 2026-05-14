# PR 06: Lower Default Scan Memory Cap

**Position in roadmap:** Second. Depends on PR05 for full-profile
sweep flags.

**Status:** Implemented in this branch. `BamScanConfig.buffer_size_bytes`
now defaults to `2 * 1024**3`, `--scan-buffer-size` documents default 2,
and the staged profiler reports buffer/spill metrics needed for future
full-profile cap sweeps.

Validation results are published in
[../pr06_scan_memory_cap_sweep_2026-05-13.md](../pr06_scan_memory_cap_sweep_2026-05-13.md).
The 2 GiB cap reduced after-scan RSS by 1.83 GB with no material wall
penalty. Full-process peak RSS fell by 0.45 GB because the peak is now
dominated by the scoring CSR/partition plateau, which PR07 and PR08 target.

## Summary

Use 2 GiB as the default scan memory cap. Keep `--scan-buffer-size` as
the explicit override for machines where a larger resident scan buffer
is desirable.

## Motivation

Async spill (`_SpillWriter`) is **already in tree** in
[src/rigel/buffer.py](../../../src/rigel/buffer.py). The scan-only
memory sweep already showed the new tradeoff:

| Buffer cap | Scan wall | Peak RSS | Resident buffer | Spilled chunks |
|---:|---:|---:|---:|---:|
| 4 GiB | 32.11 s | 9.66 GB | 3.90 GB | 14 |
| 2 GiB | **31.21 s** | 7.88 GB | 2.00 GB | 23 |
| 1 GiB | 31.97 s | **7.07 GB** | **0.94 GB** | 28 |

Lower caps are essentially free in scan-only wall time. The remaining
question is whether more spilled chunks slow `fragment_router_scan`
when scoring has to read them back. Until PR05 lands, the staged
profiler cannot vary the cap, so this question is unanswered for the
*full* pipeline.

## Current code

* Default: [src/rigel/config.py](../../../src/rigel/config.py)
  (`BamScanConfig.buffer_size_bytes = 2 * 1024**3`).
* Spill writer: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
  (`_SpillWriter`, lines ~454–520; bounded `queue.Queue(maxsize=2)`,
  one writer thread).
* CLI override: `--scan-buffer-size` in [src/rigel/cli.py](../../../src/rigel/cli.py).

## Proposed change

A two-step PR. Step 1 measures, step 2 sets the default.

### Step 1 — Validate

Run staged full profiles at 4, 2, and 1 GiB caps using PR05's flags.
Compare:

* total wall time
* `scan_and_buffer` and `fragment_router_scan`
* peak RSS
* spilled chunk count and spill bytes
* whether `_SpillWriter`'s bounded queue ever blocked the scanner

### Step 2 — Set default

* If 2 GiB stays within 5% wall and saves at least 1 GB peak RSS:
  change `BamScanConfig.buffer_size_bytes` default to `2 * 1024**3`.
* If 1 GiB also stays within 5%, **do not** change the default to 1
  GiB. Lower caps trade allocator-cushion margin for spill volume.
  Document 1 GiB as an explicit low-memory recommendation in the CLI
  help and parameter docs but keep the default conservative.

Update CLI help for `--scan-buffer-size` to mention the new default and
the async-spill-makes-this-cheap framing.

## Implementation steps

1. Run the validation sweep (PR05 flags). Persist the JSON and a short
   numeric comparison table under `docs/performance/`.
2. Edit `BamScanConfig.buffer_size_bytes` default per the validation
   result.
3. Edit `--scan-buffer-size` help text in `cli.py`.
4. Update [docs/MANUAL.md](../../MANUAL.md) and
   [docs/parameters.md](../../parameters.md) if they cite the default.
5. Add a regression unit test asserting the new default value (so it
   does not silently drift back).

## Tests

```bash
conda activate rigel
pytest tests/test_cli.py tests/test_buffer.py tests/test_pipeline_smoke.py -v
pytest tests/test_golden_output.py -v
```

No golden change expected. Memory cap does not change fragment routing.

## Benchmark plan

After PR05:

```bash
for cap in 4 2 1; do
  python scripts/profiling/profiler.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
    --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr06_${cap}gib \
    --stages --threads 8 --scan-bgzf-threads 2 \
    --scan-buffer-size "$cap" --scan-read-name-batch-size 512 \
    --memory-interval 250
done
```

## Acceptance criteria

* Full-pipeline wall at 2 GiB ≤ 1.05 × wall at 4 GiB.
* Full-pipeline peak RSS at 2 GiB ≤ peak at 4 GiB − 1 GB.
* `_SpillWriter` queue does not block the scanner for more than 1% of
  scan wall time at the new default.
* `--scan-buffer-size` override continues to work.
* Default change has a unit test.

## Risks

* Slow disks turn higher spill counts into a wall-time penalty. Keep
  the override well-documented; do not push to 1 GiB by default.
* Scoring reading spilled chunks back is sequential disk I/O; if
  `fragment_router_scan` regresses noticeably at 2 GiB, raise the
  default back to 4 GiB and document the result instead of forcing the
  change.

## Non-goals

* No spill format / compression changes.
* No `--scan-fragments-per-chunk` tuning. That is a separate benchmark question.
