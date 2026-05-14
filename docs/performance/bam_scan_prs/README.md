# BAM Scan Performance PR Roadmap (revised)

Date: 2026-05-13

This is a rewrite of the original 7-PR roadmap. The original measurements
are sound, but the order and shape of the PRs needed restructuring:

* The s4 → s8 → s12 thread sweep (135.7 s → 204.2 s → 226.9 s) shows that
  scanning gets *slower* with more workers. Until the queue pathology is
  fixed, every other measurement is contaminated.
* Direct scan-to-score streaming is not currently a valid single-pass
  design: scoring depends on finalized strand/FL models and calibration,
  and those are produced after the scan. The implementable fix is to
  keep the single-pass architecture but make chunk spill asynchronous so
  Arrow/LZ4 serialization is removed from the scanner callback path.
* The resolver rewrite should use standard sorted-vector algorithms and
  avoid building a new helper layer. The only new primitive worth adding
  is a tiny `sort_unique` helper; existing allocation-returning hot
  helpers must be rewritten in place or bypassed.
* The original PR-07 forks `build_fragment` into a fast/slow path. That
  doubles the test surface and creates drift risk. Wait for the resolver
  scratch-vector pattern to land, then apply the same pattern in place.

Workload used for evidence (unchanged from the original roadmap):

* BAM: `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`
* Index: `/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index`
* 66.1M BAM records, 31.99M read names, 32.46M buffered fragments
* macOS arm64, Apple Silicon M3-class machine
* Native build: `RelWithDebInfo`, `RIGEL_PROFILE_NATIVE=ON`

## Polars Migration Assessment

Unchanged from the original roadmap. Polars is not the right tool for
this bottleneck — the hot frames are native C++ scanning, queue
synchronisation, htslib tag/CIGAR parsing, resolver allocation, and
synchronous spill. None of these are pandas-bound. Keep Polars as a
targeted experiment for offline analysis only.

## Key Measurements (baseline before any of the new PRs)

| Run | Wall time | Read names/s | Spills | Peak RSS |
|---|---:|---:|---:|---:|
| s8 d4, 4GiB spill baseline | 226.6s | 141k | 14 | 9.2GB |
| s8 d2, 4GiB spill | 200.9s | 159k | 14 | 9.3GB |
| s8 d8, 4GiB spill | 198.7s | 161k | 14 | 9.4GB |
| s4 d2, 12GiB no-spill | **135.7s** | **236k** | 0 | 10.6GB |
| s8 d2, 12GiB no-spill | 204.2s | 157k | 0 | 11.1GB |
| s12 d2, 12GiB no-spill | 226.9s | 141k | 0 | 11.5GB |

**Three things stand out, and the PRs below attack them without
conflating unrelated changes:**

1. **Negative scan-thread scaling.** s8 is 50% slower than s4 at
   no-spill. Cause: per-read-name queue traffic. → PR 01.
2. **Allocator churn in the resolver.** Live samples show
   `unordered_set` / `unordered_map` allocation dominating
   `_resolve_core`. → PR 02.
3. **Synchronous Arrow/LZ4 spill blocks the scan callback.** Spill costs
  ≈ 65 s on the 4 GiB / s8 run. Direct scan-to-score streaming is
  blocked by calibration/model dependencies, so the practical fix is an
  asynchronous chunk store. → PR 04.

A pile of small per-record inefficiencies (lazy SJ strand parsing,
duplicated `ambig_strand` computation, FL hash map) accounts for a few
percent each and ships together as one local-hygiene PR. → PR 03.

## Revised PR Series

| Order | PR | Doc | Status |
|---:|---|---|---|
| 01 | Batch qname queue work | [pr01_batch_qname_queue.md](pr01_batch_qname_queue.md) | **must land first** — gates downstream measurements |
| 02 | Resolver scratch vectors | [pr02_resolver_scratch_vectors.md](pr02_resolver_scratch_vectors.md) | broad but representation-only rewrite after PR 01 lands |
| 03 | Small-wins resolver hygiene | [pr03_smallwins_resolver_hygiene.md](pr03_smallwins_resolver_hygiene.md) | three local micro-cleanups; easiest after PR 02 because (c) touches resolver state |
| 04 | Async chunk store for spill | [pr04_async_chunk_store.md](pr04_async_chunk_store.md) | keeps one-pass semantics; removes synchronous spill from scanner callback |
| — | Fragment assembly allocator cleanup (deferred) | [pr_deferred_fragment_assembly.md](pr_deferred_fragment_assembly.md) | revisit only if `build_fragment` is still hot after PR 02 |

### Why this order

* **PR 01 first.** Every wall-time number on the table above is a function
  of how queue traffic interacts with the rest of the pipeline. The
  s4→s8 regression is a 50% scaling cliff; if we leave it in place, every
  PR that ships after will be measured against a baseline that lies. PR
  01 is also small (one file change, one new struct).
* **PR 02 second.** Once queues are not the bottleneck, the resolver is.
  PR 02 is a focused C++ rewrite of one function with a clear acceptance
  criterion (no behavior change). Land it before structural changes so
  PR 04 doesn't carry resolver risk.
* **PR 03 third.** The changes are small, but sub-change (c) edits
  `RawResolveResult` and the resolver/accumulator handoff, so it is
  cleanest after PR 02 has settled the scratch-vector representation.
* **PR 04 fourth.** The scanner already emits chunks through a Python
  callback. PR 04 keeps that API and makes `FragmentBuffer` spill
  asynchronous. This avoids the false promise of direct scan-to-score
  concurrency while preserving the single-pass training/calibration
  model.

## Shared Validation Protocol

Every PR must preserve transcript-level semantics and exact fragment
routing unless it explicitly documents an intentional representation-only
change.

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
  --threads 4 8 12 \
  --scan-bgzf-threads 2 \
  --scan-fragments-per-chunk 1000000 \
  --scan-buffer-size 12
```

The thread sweep is mandatory. Single-thread-count benchmarks hide
scaling regressions like the one PR 01 fixes.

## Baseline comparison points

* `nospill_s4_d2`: 135.7s, 235.6k read names/s, 10.6GB peak RSS
* `spill_s8_d2`: 200.9s, 159.2k read names/s, 9.3GB peak RSS, 14 spills
* `nospill_s8_d2`: 204.2s, 156.7k read names/s, 11.1GB peak RSS
