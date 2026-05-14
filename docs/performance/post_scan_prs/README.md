# Post-Scan Performance PR Roadmap (revised)

Date: 2026-05-13

Context: PR01–PR04 (scan and BAM I/O) landed and reshaped the profile.
This roadmap turns the resulting full profile in
[../profile_2026-05-13_vcap_rna20m_gdna20m_post_scan_prs.md](../profile_2026-05-13_vcap_rna20m_gdna20m_post_scan_prs.md)
into an implementation-ready PR series.

The pipeline is no longer dominated by a single hotspot. Total wall time
on the VCAP RNA20M + gDNA20M run is now **73.36 s**, down from 246.8 s.
The remaining work is data-movement: large intermediate arrays, a global
CSR with a long lifetime, repeated scatter passes, and a single-thread
scoring pass.

## Current baseline

Clean staged profile from before the naming cleanup: 8 native scan workers,
4 BGZF threads, 8 EM threads, 4 GiB scan buffer cap, and read-name batch
size 512. Under the new total-budget convention, the same scan worker
split is `--threads 12 --scan-bgzf-threads 4`; future profiles should
report both the requested total budget and the resolved scan worker count.

| Stage | Wall | Share |
|---|---:|---:|
| `scan_and_buffer` | 32.78 s | 44.7% |
| `fragment_router_scan` | 14.72 s | 20.1% |
| `locus_em` | 14.30 s | 19.5% |
| `partition_and_free` | 6.15 s | 8.4% |
| `compute_eb_gdna_priors` | 4.23 s | 5.8% |
| `build_loci` | 1.14 s | 1.6% |

Peak RSS is 16.20 GB, reached after `fragment_router_scan` builds the
global CSR on top of the scan buffer high-water.

## Verified ground truth before planning

These are facts in the tree at the time of this revision; every PR below
is calibrated against them.

- The scan parameter vocabulary is canonicalized before this PR series:
  `threads` is the total budget, `scan_bgzf_threads` reserves BGZF
  decompression threads from it, `scan_buffer_size` is GiB,
  `scan_fragments_per_chunk` controls native chunk handoff size, and
  `scan_read_name_batch_size` controls read-name queue batches.
- Async spill is already implemented as `_SpillWriter` in
  [src/rigel/buffer.py](../../../src/rigel/buffer.py). It is a single
  background thread with a bounded queue. PR04 is *streaming* scan→score,
  not spill.
- `partition_and_free` is Python in
  [src/rigel/locus_partition.py](../../../src/rigel/locus_partition.py);
  it sequences 1 offsets call + 11 scatter calls into the templated
  `scatter_candidates_impl<T>` / `scatter_units_impl<T>` in
  [src/rigel/native/em_solver.cpp](../../../src/rigel/native/em_solver.cpp).
- `StreamingScorer` writes `unambig_counts` directly into a single
  estimator-owned `f64_2d_mut` array (shared mutable state). Any
  parallel scorer must give every worker its own copy and merge.
- `score_chunk()` consumes 12 buffer columns after PR10. `intron_bp`,
  `exon_bp_pos`, `exon_bp_neg`, `tx_bp_pos`, and `tx_bp_neg` are not part
  of the scan-buffer/scorer ABI.
- `enable_gdna_for_multilocus` reads `is_spliced` and `gdna_log_liks`
  from the global CSR; it must run before `partition_and_free` consumes
  those arrays.

## PR series

Status in this branch: PR05, PR06, PR07, and PR10 are implemented. The table
below keeps the original ordering for the remaining PR series.

| Order | PR | Doc | Primary target | Risk |
|---:|---|---|---|---|
| 05 | Scan config & profiler visibility | [pr05_profile_and_scan_config_visibility.md](pr05_profile_and_scan_config_visibility.md) | Implemented: reproducible scan config and memory-shape metrics | Low |
| 06 | Lower default scan memory cap | [pr06_scan_memory_budget_policy.md](pr06_scan_memory_budget_policy.md) | Implemented: default scan cap is now 2 GiB | Low |
| 07 | Float32 scored-fragment payloads | [pr07_float32_scored_fragments.md](pr07_float32_scored_fragments.md) | Implemented: 3.51 GB lower peak RSS on VCAP | Medium |
| 08 | Fused partition scatter | [pr08_fused_partition_scatter.md](pr08_fused_partition_scatter.md) | `partition_and_free` wall + memory traffic | Medium |
| 09 | Parallel streaming scorer | [pr09_parallel_streaming_scorer.md](pr09_parallel_streaming_scorer.md) | `fragment_router_scan` wall | Medium-high |
| 10 | Buffer column diet | [pr10_buffer_dtype_diet.md](pr10_buffer_dtype_diet.md) | Implemented: dead buffer/scorer fields removed; 0.46 GB lower peak RSS vs PR07 | Medium |
| 11 | Prior fast path | [pr11_prior_fast_path.md](pr11_prior_fast_path.md) | `compute_eb_gdna_priors` Python tail | Medium |
| 12 | EM high-iteration workset | [pr12_em_high_iteration_workset.md](pr12_em_high_iteration_workset.md) | Long-tail SQUAREM iterations | Research |

## Why this order

**PR05 first** because every later wall-time and RSS claim depends on
the profiler being able to set the scan parameters it reports on, and
on the JSON exposing array cardinalities (`n_candidates`, candidate /
unit byte estimates). It is small, low-risk, and instrumentation-only.

**PR06 second** because the scan-only memory sweep already shows lower
caps are nearly free. Async spill is already in tree; the only thing
left is to validate end-to-end and change a default.

**PR07 before PR09.** Parallel scoring multiplies in-flight CSR
pressure; halving payload precision first prevents the parallel scorer
from saturating memory bandwidth.

**PR08 before PR09** for the same reason: the scatter step is also
bandwidth-bound, and fusing it removes one of the two largest copy
events in the post-scan pipeline.

**PR10 in parallel with PR07/PR08** if reviewers are available. It is
isolated to the buffer layer and does not touch the scoring CSR.
Schedule independently.

**PR11 lower priority** because prior assembly already dropped from
9.5 s to 4.23 s. The remaining cost is real but small.

**PR12 is investigative.** EM is already threaded and
correctness-sensitive. Build a workset and a microbenchmark harness
before changing solver behaviour.

## Shared validation protocol

After any C++ change under `src/rigel/native/`, recompile:

```bash
conda activate rigel
pip install --no-build-isolation -e .
```

Minimum correctness checks for this PR series:

```bash
conda activate rigel
pytest tests/test_golden_output.py -v
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py -v
pytest tests/test_em_impl.py tests/test_estimator.py -v
```

For representation changes, run the full suite before merge:

```bash
pytest tests/ -q
ruff check src/ tests/
```

Performance comparison workload for the post-PR06 baseline:

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_post_scan_pr_check \
  --stages --threads 8 --scan-buffer-size 2 --memory-interval 250
```

Use cProfile only for attribution, not for headline timing.

## Non-goals for the series

- No transcript-centric modeling changes.
- No EM convergence-tolerance loosening.
- No silent removal of diagnostic outputs. PRs that change diagnostics
  must be config-gated and test both modes.
- No bundling of multiple representation changes into one PR. Each
  representation change needs its own golden-output and profile diff.
