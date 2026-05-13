# PR 06: Scan Memory Budget Policy

**Position in roadmap:** Second. Depends on PR05 profiler/config
visibility so the full-profile validation can set and report scan memory
caps directly.

## Summary

Use the async spill writer from PR04 to reduce scan memory pressure. Run
full-profile validation at 4 GiB, 2 GiB, and 1 GiB scan caps, then either
lower the default scan buffer cap to 2 GiB or document it as the preferred
large-run setting.

## Motivation

The scan-only memory sweep showed that lower buffer caps are now nearly
free:

| Buffer cap | Scan wall | Peak RSS | Resident buffer | Spilled chunks |
|---:|---:|---:|---:|---:|
| 4 GiB | 32.11 s | 9.66 GB | 3.90 GB | 14 |
| 2 GiB | 31.21 s | 7.88 GB | 2.00 GB | 23 |
| 1 GiB | 31.97 s | 7.07 GB | 0.94 GB | 28 |

The old synchronous spill path made lower caps expensive. The current
async writer changes the tradeoff. This is the simplest immediate RSS
reduction available.

## Current code

- Default scan memory cap: [src/rigel/config.py](../../../src/rigel/config.py)
- Buffer and async spill lifecycle: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
- CLI `--buffer-size`: [src/rigel/cli.py](../../../src/rigel/cli.py)
- Full profiler, after PR05: [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py)

## Proposed change

Make the scan memory policy explicit and evidence-based:

- Preferred large-run default: 2 GiB if full-profile validation confirms
  no material wall-time penalty.
- Keep `--buffer-size` as an override.
- Add documentation explaining the tradeoff: lower caps spill more chunks
  but can reduce RSS substantially because spill writes are asynchronous.
- Add a regression benchmark target that compares 4/2/1 GiB caps.

## Implementation steps

1. Use PR05 profiler flags to run full quant profiles at 4 GiB, 2 GiB,
   and 1 GiB caps with identical scan/EM thread settings.
2. Compare full-run wall time, peak RSS, `scan_and_buffer`,
   `fragment_router_scan`, and spill counts. The full profile matters
   because scoring must load spilled chunks back from disk.
3. If 2 GiB is within 5% wall time of 4 GiB and saves at least 1 GB RSS,
   change `BamScanConfig.max_memory_bytes` default from `4 * 1024**3` to
   `2 * 1024**3`.
4. Update CLI help for `--buffer-size` to describe the new default and
   mention async spill.
5. Update docs and profile report references that say default 4 GiB.
6. Add a small unit test that `BamScanConfig()` default matches the
   documented default.
7. Add a profiling helper note or script command in
   [docs/performance](../) for the 4/2/1 GiB comparison.

## Tests

```bash
conda activate rigel
pytest tests/test_cli.py tests/test_buffer.py tests/test_pipeline_smoke.py -v
pytest tests/test_golden_output.py -v
```

No golden changes are expected. This PR changes memory policy, not
fragment routing or EM semantics.

## Benchmark plan

After PR05:

```bash
conda activate rigel
for cap in 4 2 1; do
  python scripts/profiling/profiler.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
    --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr06_${cap}gib \
    --stages --threads 8 --n-decomp-threads 2 \
    --max-memory-gib "$cap" --qname-batch-size 512 \
    --memory-interval 250
done
```

Use no-cProfile runs for headline timing.

## Acceptance criteria

- Full quant at 2 GiB cap is no more than 5% slower than 4 GiB.
- Full quant at 2 GiB cap reduces peak RSS by at least 1 GB.
- If 1 GiB is also neutral, document it as an aggressive low-memory
  setting but do not make it default without broader workload coverage.
- `--buffer-size` override continues to work.
- No correctness or golden-output changes.

## Risks

- Scan-only results may not predict full quant if reading more spilled
  chunks slows scoring. That is why full-profile validation is required.
- More spills can stress slower disks. Keep the override easy to find.

## Non-goals

- Do not rewrite spill serialization or compression in this PR.
- Do not change chunk size. Chunk-size tuning should be a separate
  benchmark if needed.