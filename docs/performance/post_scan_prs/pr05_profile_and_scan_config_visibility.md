# PR 05: Scan Config & Profiler Visibility

**Position in roadmap:** First. Every later PR's measurement claims
depend on this landing.

**Status:** Implemented in this branch. The profiler now records
`scan_config`, `buffer_summary`, `scoring_csr`, `partition_bytes_total`,
and `partition_bytes` in `profile_summary.json`.

## Summary

Finish scan/profiler visibility on top of the pre-PR05 parameter rename,
and extend the profile JSON so memory-related claims can be quantified
without source-diving. This PR does not change quantification semantics.

## Motivation

The parameter cleanup before this PR intentionally removed the old scan
names. The canonical user-facing surface is now:

* `--threads` / `threads`: total thread budget.
* `--scan-bgzf-threads` / `scan_bgzf_threads`: BGZF decompression
  threads reserved from that budget during scan.
* `--scan-buffer-size` / `scan_buffer_size`: scan-buffer memory cap in GiB.
* `--scan-fragments-per-chunk` / `scan_fragments_per_chunk`: native scan
  chunk handoff size.
* `--scan-read-name-batch-size` / `scan_read_name_batch_size`: read-name
  groups per native scanner queue item.

The remaining gap is measurement depth:

* `rigel quant`, config YAML, `scripts/profiling/profiler.py`, and
  `scripts/profiling/scan_profile.py` all use the canonical names above.
* `profile_summary.json` does not record array cardinalities. PR06–PR10
  all want to claim "this many bytes saved" and the current JSON cannot
  back the claim.

## Current code

* CLI `_ParamSpec` registry: [src/rigel/cli.py](../../../src/rigel/cli.py)
  (line ~541 onward). Pattern: `_ParamSpec(cli_dest, config_path,
  transform)`.
* Config: [src/rigel/config.py](../../../src/rigel/config.py)
  (`BamScanConfig` lines ~107–154).
* Staged profiler: [scripts/profiling/profiler.py](../../../scripts/profiling/profiler.py)
  (CLI parser ~line 990, summary writer ~line 948).
* Scan-only profiler with the desired flag set: [scripts/profiling/scan_profile.py](../../../scripts/profiling/scan_profile.py).

## Proposed change

### Parameter visibility

Preserve the canonical names from the pre-PR05 cleanup. Do not re-add
legacy spellings such as `--buffer-size`, `--qname-batch-size`,
`--n-scan-threads`, `--n-decomp-threads`, `--chunk-size`, or
`--max-memory-gib`. The scan worker count must remain derived from
`threads - scan_bgzf_threads` with at least one scan worker.

### Profile JSON

Add fields under each profile entry:

```jsonc
{
  "scan_config": {
    "threads": 8,
    "resolved_total_threads": 8,
    "scan_worker_threads": 4,
    "scan_bgzf_threads": 4,
    "requested_scan_bgzf_threads": 4,
    "scan_fragments_per_chunk": 1000000,
    "scan_read_name_batch_size": 512,
    "scan_buffer_size_bytes": 4294967296
  },
  "scoring_csr": {
    "n_units": 30391824,
    "n_candidates": 263847291,
    "mean_candidates_per_unit": 8.68,
    "candidate_bytes": {
      "log_liks": 2110778328,
      "coverage_weights": 2110778328,
      "t_indices": 1055389164,
      "count_cols": 263847291
    },
    "unit_bytes": { /* same shape, per-unit arrays */ }
  },
  "partition_bytes_total": 6800000000,
  "buffer_summary": {
    "memory_bytes_peak": 4181590000,
    "chunks_finalized": 33,
    "chunks_spilled": 14,
    "chunks_pending_spill_peak": 2
  }
}
```

Compute byte estimates from `array.nbytes`. Capture them **before**
`partition_and_free` nulls the global arrays. Do not retain the arrays
themselves.

## Implementation steps

1. Confirm the canonical scan flags are present in both the quant CLI and
  staged profiler, and that the serialized/profile JSON uses the new
  key names only.
2. Capture `scoring_csr` byte stats inside `fragment_router_scan` after
   the scorer's `finish()` returns, before partition.
3. Capture `partition_bytes_total` after `partition_and_free` by
   summing `array.nbytes` across the returned `LocusPartition` dict.
4. Capture `buffer_summary` from the existing `FragmentBuffer` accessors
   (`memory_bytes`, spill counters in `_SpillWriter`).
5. Write fields into `profile_summary.json` and the text report.
6. Update [docs/parameters.md](../../parameters.md) and any CLI param
  reference if the profile JSON schema changes.

## Tests

```bash
conda activate rigel
pytest tests/test_cli.py tests/test_pipeline_smoke.py -v
pytest tests/test_golden_output.py -v
```

Add focused tests:

* `--scan-bgzf-threads 2` writes `scan_bgzf_threads: 2` into the
  serialized run config and `scan.bgzf_threads == 2` in `BamScanConfig`.
* The staged profiler accepts each new flag and the resulting
  `profile_summary.json` records the value.
* The byte fields exist and are positive on a small fixture run.

## Benchmark plan

Re-run the clean staged profile at the established baseline to confirm
no regression from the new instrumentation:

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr05_visibility \
  --stages --threads 8 --scan-bgzf-threads 4 \
  --scan-buffer-size 4 --scan-read-name-batch-size 512 \
  --memory-interval 250
```

## Acceptance criteria

* No golden-output changes.
* `rigel quant --scan-bgzf-threads 2` is accepted and reflected in the
  serialized run config.
* The staged profiler can run scan caps of 4 / 2 / 1 GiB without YAML
  edits.
* `profile_summary.json` contains `scan_config`, `scoring_csr` (with
  `n_candidates` and per-array `nbytes`), `partition_bytes_total`, and
  `buffer_summary`.
* Profiler wall-time overhead from the new instrumentation is under
  100 ms on the VCAP run.

## Risks

* Don't retain large arrays just to measure them. Read `nbytes`, drop
  the reference, then proceed with the existing release path.
* Do not re-add legacy aliases. Old config names should fail with a clear
  replacement message; old CLI flags should remain unregistered.

## Non-goals

* No default changes (PR06 owns the default-cap question).
* No layout changes (PR07–PR10).
