# BAM Scan Performance PR Roadmap

Date: 2026-05-13

Scope: Convert the VCAP `scan_and_buffer` profiling results into a sequence of
reviewable performance and memory pull requests.

Workload used for evidence:

- BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`
- Index: `/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index`
- Size: 66.1M BAM records, 31.99M read names, 32.46M buffered fragments
- Platform: macOS arm64, Apple Silicon M3-class machine
- Build for native stacks: `RelWithDebInfo`, `RIGEL_PROFILE_NATIVE=ON`

## Polars Migration Assessment

Polars is a strong dataframe engine, and it is often much faster than pandas for
large table transforms. For this specific bottleneck, a broad migration from
pandas and pyarrow to Polars is not expected to be a major performance win.

Reasons:

1. The dominant stage is native C++ BAM scanning, not pandas dataframe work.
   The hot frames are `BamScanner::process_qname_group_threaded`,
   `FragmentResolver::_resolve_core`, htslib tag/CIGAR parsing, queue
   synchronization, and temporary Arrow spill writes.
2. The scan spill path is not pandas-based. It converts NumPy arrays and Arrow
   list arrays to Feather/IPC via pyarrow. Polars is Arrow-backed and can read
   and write IPC/Feather, but it would still serialize the same columnar data and
   still pay compression/I/O cost unless we change the dataflow.
3. The measured spill problem is synchronous serialization on the scanner's
   Python callback path. Replacing `pyarrow.feather.write_feather` with a Polars
   write does not address the core issue: the scanner is blocked while temporary
   chunk storage is written.
4. Several Python bottlenecks in older profiles are NumPy-loop or algorithmic
   issues, not dataframe expression issues. These should be handled with NumPy,
   C++, or more direct columnar algorithms.

Recommended stance:

- Do not start with a broad Polars migration for performance.
- Keep Polars as a targeted experiment for offline benchmark/report analysis or
  isolated dataframe-heavy preprocessing if a profile shows pandas as the hot
  path.
- For BAM scan performance, focus on removing synchronous spill, reducing queue
  operations, and cutting resolver allocation churn.

## Key Measurements

| Run | Wall time | Read names/s | Spills | Peak RSS |
|---|---:|---:|---:|---:|
| s8 d4, 4GiB spill baseline | 226.6s | 141k | 14 | 9.2GB |
| s8 d2, 4GiB spill | 200.9s | 159k | 14 | 9.3GB |
| s8 d8, 4GiB spill | 198.7s | 161k | 14 | 9.4GB |
| s4 d2, 12GiB no-spill | 135.7s | 236k | 0 | 10.6GB |
| s8 d2, 12GiB no-spill | 204.2s | 157k | 0 | 11.1GB |
| s12 d2, 12GiB no-spill | 226.9s | 141k | 0 | 11.5GB |

Interpretation:

- Spilling is a wall-time-class bottleneck.
- More htslib decompression threads are not the main lever on this workload.
- More scan workers currently make the scanner slower once spill is removed,
  because per-read-name queue synchronization dominates.
- Resolver work remains CPU-heavy and allocation-heavy in every sample.

## Proposed PR Series

The PRs below are intended to be developed independently over several turns.
Each plan includes implementation scope, tests, benchmarks, acceptance criteria,
and known risks.

| PR | Title | Main target |
|---:|---|---|
| 01 | [Remove synchronous scan spill](pr01_async_spill_and_streaming.md) | Recover wall time lost to Arrow/LZ4 spill blocking the scan callback |
| 02 | [Batch qname queue work](pr02_batch_qname_queue.md) | Reduce queue mutex/condition-variable traffic and make parallelism useful |
| 03 | [Rewrite resolver set logic with scratch vectors](pr03_resolver_scratch_vectors.md) | Remove hot `unordered_set`/`unordered_map` allocation churn in `_resolve_core` |
| 04 | [Use t-aligned fragment-length vectors](pr04_t_aligned_frag_lengths.md) | Remove `compute_frag_lengths` hash-map allocation and lookup |
| 05 | [Parse SJ strand tags lazily](pr05_lazy_sj_strand_tags.md) | Avoid expensive `bam_aux_get`/aux scanning for unspliced records |
| 06 | [Store `ambig_strand` in the accumulator](pr06_accumulator_ambig_strand.md) | Remove redundant candidate scans during chunk finalization |
| 07 | [Add a common-case fragment assembly fast path](pr07_fragment_assembly_fast_path.md) | Avoid maps/sets for ordinary paired-end fragments |

Suggested implementation order:

1. PR 06, PR 05, and PR 04 are local, low-risk changes with focused tests.
2. PR 01 attacks the largest measured wall-time issue and may affect memory
   behavior, so it should land early but after the easiest invariants are fixed.
3. PR 02 should follow PR 01 so queue scaling is measured without spill noise.
4. PR 07 can land before or after PR 02; it is a common-case speed path with a
   fallback to current behavior.
5. PR 03 is the broadest native rewrite and should be developed after the
   smaller changes establish clean regression and benchmark baselines.

## Shared Validation Protocol

Every PR should preserve transcript-level semantics and exact fragment routing
unless the PR explicitly documents an intentional representation-only change.

Minimum tests for native scanner changes:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_pipeline_smoke.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
```

Recommended performance smoke after each PR:

```bash
conda activate rigel
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr_check \
  --name-prefix pr_check \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

Primary benchmark comparison points:

- `nospill_s4_d2`: 135.7s, 235.6k read names/s, 10.6GB peak RSS
- `spill_s8_d2`: 200.9s, 159.2k read names/s, 9.3GB peak RSS, 14 spills
- `nospill_s8_d2`: 204.2s, 156.7k read names/s, 11.1GB peak RSS
