# Post-Scan Performance PR Roadmap

Date: 2026-05-13

This roadmap turns the post-PR01-PR04 full profile in
[../profile_2026-05-13_vcap_rna20m_gdna20m_post_scan_prs.md](../profile_2026-05-13_vcap_rna20m_gdna20m_post_scan_prs.md)
into an implementation-ready PR series.

The scan PRs worked: total wall time on the VCAP RNA20M + gDNA20M run is
now 73.36 s, down from 246.8 s, and `scan_and_buffer` is now 32.78 s,
down from 200.6 s. The remaining problem is no longer one pathological
hotspot. It is data movement: large intermediate arrays, global CSR
lifetime, repeated scatter, and a single-thread scoring pass.

## Current baseline

Clean staged profile, 8 scan/EM threads, 4 GiB scan buffer cap,
`qname_batch_size=512`, current profiler defaults:

| Stage | Wall | Share |
|---|---:|---:|
| `scan_and_buffer` | 32.78 s | 44.7% |
| `fragment_router_scan` | 14.72 s | 20.1% |
| `locus_em` | 14.30 s | 19.5% |
| `partition_and_free` | 6.15 s | 8.4% |
| `compute_eb_gdna_priors` | 4.23 s | 5.8% |
| `build_loci` | 1.14 s | 1.6% |

Peak RSS is 16.20 GB. RSS rises to 9.81 GB after scan, then to 16.20 GB
after scoring/router scan, and stays there through cleanup.

## PR series

| Order | PR | Doc | Primary target | Risk |
|---:|---|---|---|---|
| 05 | Profiling and scan config visibility | [pr05_profile_and_scan_config_visibility.md](pr05_profile_and_scan_config_visibility.md) | Measurement correctness, visible `n_decomp_threads` | Low |
| 06 | Scan memory budget policy | [pr06_scan_memory_budget_policy.md](pr06_scan_memory_budget_policy.md) | Lower scan RSS using async spill | Low |
| 07 | Float32 scored-fragment payloads | [pr07_float32_scored_fragments.md](pr07_float32_scored_fragments.md) | Peak RSS and bandwidth in scoring/scatter/EM | Medium |
| 08 | Fused partition scatter | [pr08_fused_partition_scatter.md](pr08_fused_partition_scatter.md) | `partition_and_free` wall time and memory traffic | Medium |
| 09 | Parallel streaming scorer | [pr09_parallel_streaming_scorer.md](pr09_parallel_streaming_scorer.md) | `fragment_router_scan` wall time | Medium-high |
| 10 | Buffer dtype diet | [pr10_buffer_dtype_diet.md](pr10_buffer_dtype_diet.md) | Scan buffer logical size and spill volume | Medium |
| 11 | Prior fast path and eligibility fusion | [pr11_prior_fast_path.md](pr11_prior_fast_path.md) | `compute_eb_gdna_priors` and Python gathers | Medium |
| 12 | EM high-iteration workset | [pr12_em_high_iteration_workset.md](pr12_em_high_iteration_workset.md) | Remaining EM long tail | Research PR |

## Why this order

PR05 comes first because the profiler currently cannot set every scan
parameter it reports on, and the full profile JSON does not expose enough
array-cardinality information to prove memory improvements. It is a small
instrumentation PR that makes every later benchmark more trustworthy.

PR06 is next because the scan-only memory sweep already shows a nearly
free RSS reduction: 4 GiB cap was 32.11 s / 9.66 GB peak, 2 GiB cap was
31.21 s / 7.88 GB peak, and 1 GiB cap was 31.97 s / 7.07 GB peak.

PR07 and PR08 should land before PR09. Parallel scoring multiplies
in-flight output pressure; first shrink the payloads and reduce scatter
bandwidth so the parallel scorer does not simply move the bottleneck to
memory bandwidth.

PR10 can run in parallel with PR07/PR08, but it touches the scan buffer
ABI and spill serialization, so it should not be bundled with float32 EM
payload work.

PR11 is useful but lower priority: prior assembly improved from 9.5 s to
4.23 s after cleanup, so it is no longer the largest target. Keep it
separate because it changes output diagnostics policy.

PR12 is deliberately a research/measurement PR. EM is already threaded
and correctness-sensitive; do not change convergence tolerances or solver
semantics until the high-iteration loci have been characterized.

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
conda activate rigel
pytest tests/ -q
ruff check src/ tests/
```

Performance comparison workload:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_post_scan_pr_check \
  --stages --threads 8 --memory-interval 250
```

Use cProfile only for attribution, not for headline timing:

```bash
conda activate rigel
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_post_scan_pr_check_cprofile \
  --stages --cprofile --threads 8 --memory-interval 250
```

## Non-goals for the series

- Do not change transcript-centric modeling semantics.
- Do not weaken EM convergence tolerances to make profiles look better.
- Do not remove diagnostic outputs unless a PR explicitly adds a config
  switch and tests both modes.
- Do not fold several memory-layout changes into one large PR. Each
  representation change needs its own golden-output and profile diff.