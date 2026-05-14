# BAM Scan Profiling On macOS

Status: ready to run, 2026-05-13.

The VCAP RNA20M + gDNA20M profile shows that `scan_and_buffer` is the dominant
runtime component: about 200 s of a 247 s profiled run. Because this stage is
mostly native C++ plus htslib, Python `cProfile` is useful for timing boundaries
but not for explaining the work inside the scanner. On macOS arm64, the best
first tool is Apple's Instruments Time Profiler, launched from `xcrun xctrace`.

## What We Need To Separate

The native scan has three concurrent pieces:

```text
reader thread:   htslib/BGZF BAM read -> qname groups -> input queue
worker threads:  qname groups -> FragmentResolver::_resolve_core -> accumulator
main thread:     output queue -> finalize_zero_copy -> Python chunk callback/spill
```

The previous profile points to four likely explanations for the 80% wall-time
share:

- BGZF decompression or BAM I/O is starving the worker queue.
- `FragmentResolver::_resolve_core` / interval overlap is the worker hot path.
- workers are blocking on the output queue while the Python chunk callback spills
  Arrow chunks synchronously.
- the input/output queue capacities are too small for this workload.

The profiling workflow below is designed to distinguish those cases.

## Build With Native Symbols

For useful C++ call stacks, rebuild once with profiler-friendly symbols and LTO
disabled:

```bash
conda activate rigel
cd /Users/mkiyer/proj/rigel
pip install --no-build-isolation -e . \
  -Ccmake.build-type=RelWithDebInfo \
  -Ccmake.define.RIGEL_PROFILE_NATIVE=ON
```

This keeps optimization on, adds DWARF symbols and frame pointers, and disables
link-time optimization so Instruments can attribute samples to functions like
`BamScanner::scan`, `process_qname_group_threaded`, `parse_bam_record`, and
`FragmentResolver::_resolve_core`.

Return to the normal fast build afterward with:

```bash
pip install --no-build-isolation -e .
```

## Scan-Only Baseline

Run the scan harness without Instruments first. It loads the index once and runs
only `scan_and_buffer`:

```bash
conda activate rigel
cd /Users/mkiyer/proj/rigel

python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap \
  --threads 8 \
  --scan-bgzf-threads 4 \
  --scan-fragments-per-chunk 1000000 \
  --scan-buffer-size 4
```

Outputs:

```text
scan_profile_summary.json
scan_profile_summary.csv
scan_t8_bgzf4/scan_result.json
scan_t8_bgzf4/memory_timeline.csv
```

The important columns are wall time, read-name throughput, buffered-fragment
throughput, peak RSS, in-memory buffer MB, and spilled chunk count.

Measured baseline from the profiler-friendly build on 2026-05-13:

```text
baseline_t8_bgzf4 wall_time_sec=226.64
read_names_per_sec=141,131
bam_records_per_sec=291,713
peak_rss_mb=9,231
buffer_memory_mb=3,900.6
buffer_spilled_chunks=14
n_bam_records=66,113,219
n_read_names=31,985,707
n_buffered_fragments=32,460,378
```

The profiling build disables LTO, so this is expected to be slower than the
normal optimized build. Treat it as the baseline for Instruments traces, not as
the release-performance number.

## Thread And Decompression Sweep

Start with a small grid so the run finishes in a reasonable time:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap_sweep \
  --threads 4 8 12 \
  --scan-bgzf-threads 2 4 8 \
  --scan-fragments-per-chunk 1000000 \
  --scan-buffer-size 4
```

Interpretation:

- improves with more `scan_bgzf_threads`: decompression/I/O was limiting.
- improves with more `threads`: resolver/worker CPU was limiting.
- flat or worse above 8 workers: queue contention, memory bandwidth, or callback
  stall is likely.
- wall time correlates with spilled chunks: synchronous spill is part of the
  bottleneck.

## Instruments Time Profiler

Launch a single scan run under Time Profiler:

```bash
xcrun xctrace record \
  --template 'Time Profiler' \
  --output /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap/time_t8_bgzf4.trace \
  --launch -- \
  /Users/mkiyer/sw/miniforge3/envs/rigel/bin/python scripts/profiling/scan_profile.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
    --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap_xctrace \
    --name-prefix xctrace \
    --threads 8 \
    --scan-bgzf-threads 4 \
    --scan-fragments-per-chunk 1000000 \
    --scan-buffer-size 4
```

  Use the absolute conda-environment Python executable for `xctrace --launch`.
  In testing, `xctrace` did not resolve `python` from the activated shell `PATH`.

Open the `.trace` in Instruments and inspect:

- Call Tree: enable `Invert Call Tree`.
- Threads: use `Separate by Thread`.
- Libraries: do not hide system libraries initially; htslib/zlib/libdeflate are
  part of the question.
- Search for `rigel._bam_impl`, `BamScanner::scan`, `process_qname_group_threaded`,
  `_resolve_core`, `sam_read1`, `bgzf`, `libdeflate`, `condition_variable`, and
  `pyarrow`.

Expected signatures:

- Reader thread dominated by `sam_read1`, `bgzf_mt_reader`, `libdeflate`, or I/O:
  tune `scan_bgzf_threads` and consider htslib/read-ahead changes.
- Worker threads dominated by `_resolve_core`, cgranges interval queries, or
  `pair_multimapper_reads`: optimize resolver internals.
- Worker threads mostly in `std::condition_variable::wait`: queue starvation or
  output backpressure, not resolver CPU.
- Main thread with Python/pyarrow samples while workers wait: move spill I/O off
  the scan path or increase queue buffering.

## Optional Follow-Up Instruments Runs

File activity is useful if spill or BAM I/O looks suspicious:

```bash
xcrun xctrace record \
  --template 'File Activity' \
  --output /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap/file_activity_t8_bgzf4.trace \
  --launch -- \
  /Users/mkiyer/sw/miniforge3/envs/rigel/bin/python scripts/profiling/scan_profile.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
    --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap_file_activity \
    --threads 8 \
    --scan-bgzf-threads 4
```

System Trace is useful if Time Profiler shows lock/wait-heavy stacks:

```bash
xcrun xctrace record \
  --template 'System Trace' \
  --output /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap/system_t8_bgzf4.trace \
  --launch -- \
  /Users/mkiyer/sw/miniforge3/envs/rigel/bin/python scripts/profiling/scan_profile.py \
    --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
    --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
    --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_vcap_system \
    --threads 8 \
    --scan-bgzf-threads 4
```

Use System Trace to answer whether workers are runnable but CPU-starved, sleeping
on queues, or blocked behind the Python callback.

## First Questions To Answer

1. Does scan wall time improve when `scan_bgzf_threads` increases from 2 to 4 or
   8?
2. Are worker threads spending CPU in `_resolve_core`, or are they sleeping on
   queues?
3. During Arrow spills, are workers blocked on output queue push?
4. Does reducing `scan_buffer_size` increase wall time strongly? If yes, spill is
   blocking the pipeline.
5. Is one worker thread much hotter than the others? If yes, qname-group work is
   imbalanced, likely due to large multimapper groups.