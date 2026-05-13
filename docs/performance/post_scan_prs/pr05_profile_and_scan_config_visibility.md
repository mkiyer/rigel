# PR 05: Profiling And Scan Config Visibility

**Position in roadmap:** First. This PR makes later measurements and
configuration changes reproducible.

## Summary

Expose the remaining scan performance controls through config, CLI, and
the full profiler. Add profile JSON fields that report the array sizes
and byte footprints responsible for peak RSS.

This PR should not change default quantification results. It is a
measurement and visibility PR.

## Motivation

The post-scan profile found two immediate visibility gaps:

- `BamScanConfig.n_decomp_threads` exists in Python, and the native
  scanner accepts it, but `rigel quant` does not expose it. The Python
  default is 4 while the native default remains 2.
- The staged profiler cannot set scan memory cap, decompression threads,
  qname batch size, or chunk size directly, even though those settings
  materially affect RSS and scan timing.

The same profile also estimates large array costs indirectly. Future PRs
need direct numbers for `em_data.n_candidates`, candidate bytes, unit
bytes, partition bytes, and scan config values.

## Current code

- Config field: [src/rigel/config.py](../../../src/rigel/config.py)
- CLI registry and parser: [src/rigel/cli.py](../../../src/rigel/cli.py)
- Full profiler config builder: [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py)
- Scan-only profiler with the desired scan flags: [scripts/profiling/scan_profile.py](../../../scripts/profiling/scan_profile.py)

## Proposed change

Add `n_decomp_threads`, `chunk_size`, and scan memory cap support to the
same declarative parameter path used by `qname_batch_size`.

Add profile output fields:

- effective scan config: `n_scan_threads`, `n_decomp_threads`,
  `chunk_size`, `qname_batch_size`, `max_memory_bytes`, `spill_dir`
- global scoring CSR shape: `n_units`, `n_candidates`, mean candidates
  per unit
- global scoring CSR byte estimate, split into candidate arrays and unit
  arrays
- partition byte estimate after `partition_and_free`
- buffer summary at end of scan: memory bytes, chunk count, spilled
  chunks, pending spills

## Implementation steps

1. Add `_ParamSpec("n_decomp_threads", "scan.n_decomp_threads")` in
   [src/rigel/cli.py](../../../src/rigel/cli.py).
2. Add a performance flag:

   ```text
   --decomp-threads N
   ```

   Use `dest="n_decomp_threads"`; document it as BGZF decompression
   threads, independent of scan worker threads.
3. Decide whether to keep `--buffer-size` as the user-facing scan memory
   knob or add an alias `--scan-max-memory-gib`. If adding an alias, map
   both to `scan.max_memory_bytes` and keep one canonical key in written
   config.
4. Add full-profiler CLI flags mirroring scan_profile:
   `--n-decomp-threads`, `--chunk-size`, `--qname-batch-size`, and
   `--max-memory-gib`.
5. Extend `_build_pipeline_config(...)` in
   [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py)
   to consume scan aliases from `RigelParams.params`:
   `n_decomp_threads`, `chunk_size`, `qname_batch_size`,
   `max_memory_bytes`, and `max_memory_gib`.
6. Add a small helper in the profiler to estimate ndarray bytes from a
   `ScoredFragments` object before it is consumed by partitioning.
7. Add a similar helper for `LocusPartition` dictionaries after scatter.
8. Write the new metrics into `profile_summary.json` and the text report.
9. Update docs that list performance parameters, including
   [docs/MANUAL.md](../../../docs/MANUAL.md) or the closest CLI parameter
   reference if present.

## Tests

```bash
conda activate rigel
pytest tests/test_cli.py tests/test_pipeline_smoke.py -v
pytest tests/test_golden_output.py -v
```

Add focused tests:

- CLI `--decomp-threads 2` writes `scan.n_decomp_threads: 2` to run
  config.
- Full profiler accepts `--n-decomp-threads 2 --max-memory-gib 2` and
  the generated JSON records those values.
- `qname_batch_size` remains visible and validated.

## Benchmark plan

Re-run the clean staged profile at the known baseline settings:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr05_visibility \
  --stages --threads 8 --n-decomp-threads 4 \
  --max-memory-gib 4 --qname-batch-size 512 --memory-interval 250
```

The stage timings should match the post-scan baseline within ordinary
run-to-run variation.

## Acceptance criteria

- No golden-output changes.
- `rigel quant --decomp-threads 2` is accepted and reflected in the run
  config.
- Full profiler can run 4 GiB, 2 GiB, and 1 GiB scan caps without YAML
  edits.
- Profile JSON reports enough data to compute bytes per candidate, bytes
  per unit, and total candidate count.
- Python and native decompression defaults are either aligned or the
  Python default is explicitly justified in docs.

## Risks

- Adding aliases can create config ambiguity. Keep written config
  canonical and test precedence.
- The profiler must not retain large objects just to measure bytes.
  Record sizes before deletion, then let existing cleanup proceed.

## Non-goals

- Do not change performance defaults in this PR except to align an
  obvious default drift if the team agrees.
- Do not implement any memory-layout optimization here.