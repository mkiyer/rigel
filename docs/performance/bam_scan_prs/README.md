# BAM Scan Performance PR Roadmap (revised)

Date: 2026-05-13

This is a rewrite of the original 7-PR roadmap. The original measurements
are sound, but the order and shape of the PRs needed restructuring:

* The s4 → s8 → s12 thread sweep (135.7 s → 204.2 s → 226.9 s) shows that
  scanning gets *slower* with more workers. Until the queue pathology is
  fixed, every other measurement is contaminated.
* Async spill is engineering complexity to work around a buffer that
  shouldn't exist on the hot path. The right fix is to stream chunks
  directly from scan to score and demote spill to a back-pressure
  safety valve.
* The original PR-03 wraps STL primitives in rigel-specific helpers
  (`sort_unique`, `intersect_sorted_into`, …). Use STL directly.
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

**Three things stand out, and all four PRs below are designed to attack
exactly one of them each:**

1. **Negative scan-thread scaling.** s8 is 50% slower than s4 at
   no-spill. Cause: per-read-name queue traffic. → PR 01.
2. **Allocator churn in the resolver.** Live samples show
   `unordered_set` / `unordered_map` allocation dominating
   `_resolve_core`. → PR 02.
3. **Synchronous Arrow/LZ4 spill blocks the scan callback.** Spill costs
   ≈ 65 s on the 4 GiB / s8 run. → PR 04.

A pile of small per-record inefficiencies (lazy SJ strand parsing,
duplicated `ambig_strand` computation, FL hash map) accounts for a few
percent each and ships together as one cosmetic-grade PR. → PR 03.

## Revised PR Series

| Order | PR | Doc | Status |
|---:|---|---|---|
| 01 | Batch qname queue work | [pr01_batch_qname_queue.md](pr01_batch_qname_queue.md) | **must land first** — gates downstream measurements |
| 02 | Resolver scratch vectors | [pr02_resolver_scratch_vectors.md](pr02_resolver_scratch_vectors.md) | broad rewrite, low risk after PR 01 lands |
| 03 | Small-wins resolver hygiene | [pr03_smallwins_resolver_hygiene.md](pr03_smallwins_resolver_hygiene.md) | three independent micro-cleanups, can land any time |
| 04 | Streaming scan→score handoff | [pr04_streaming_scan_to_score.md](pr04_streaming_scan_to_score.md) | structural, lands after PR 01–03 stabilise the inner loop |
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
  the streaming PR (04) doesn't carry resolver risk.
* **PR 03 anytime.** Three independent micro-changes packaged as one
  reviewable PR because each is too small to justify its own. Tests are
  per-change.
* **PR 04 last.** Streaming touches the public `FragmentBuffer` shape
  and the scan↔score boundary. Land it on top of a clean inner loop so
  any wall-time regression is unambiguously attributable to the
  streaming change.

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
  --n-scan-threads 4 8 12 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

The thread sweep is mandatory. Single-thread-count benchmarks hide
scaling regressions like the one PR 01 fixes.

## Baseline comparison points

* `nospill_s4_d2`: 135.7s, 235.6k read names/s, 10.6GB peak RSS
* `spill_s8_d2`: 200.9s, 159.2k read names/s, 9.3GB peak RSS, 14 spills
* `nospill_s8_d2`: 204.2s, 156.7k read names/s, 11.1GB peak RSS
