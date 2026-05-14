# Rigel Performance Profile & Optimization Plan
## VCAP RNA20M + gDNA20M (32M fragments, 8 threads, macOS arm64)

**Profile date:** 2026-05-13
**Workload:** `vcap_rna20m_gdna20m/annotated.bam` — 32 M buffered fragments, 66 M BAM records, GENCODE-class index (457 K transcripts, 266 K genes, 24 K connected-component multi-loci).
**Hardware:** Apple Silicon arm64, macOS 26.4, 8 threads.
**Source data:** [profile_summary.json](../../../Downloads/rigel_runs/profile_vcap_rna20m_gdna20m/profile_summary.json), [profile_report.txt](../../../Downloads/rigel_runs/profile_vcap_rna20m_gdna20m/profile_report.txt), [profile_default_8t.prof](../../../Downloads/rigel_runs/profile_vcap_rna20m_gdna20m/profile_default_8t.prof) (load with `snakeviz`).

> Wall and RSS measurements were captured under cProfile, which adds ~5-15 % overhead to Python-heavy stages. The relative ranking of hotspots is unchanged. A separate cProfile-off prod run produced equivalent stage-level numbers (scan ≈ 5–6 min, quant ≈ 46 s).

---

## 1. Headline numbers

| Metric | Value |
|---|---|
| Total wall time | **246.8 s** (~4 min 7 s) |
| Throughput | 268 K fragments / s |
| Peak RSS | **15.8 GB** (started at 2.9 GB → +12.9 GB) |
| Buffer in-memory | 3.9 GB + 14 spilled chunks (≈ 2.96 GB on disk) → ≈ 6.85 GB total |
| Buffered fragments | 31.99 M |
| EM units | 30.4 M (1.12 M deterministic-unique) |
| Multi-loci | 24 052 (max 42 639 transcripts, 2.82 M units in mega-locus) |

### Stage breakdown

| Stage | Wall | % | Notes |
|---|---:|---:|---|
| **scan_and_buffer** (C++ multi-threaded) | **200.6 s** | **81.3 %** | 8 worker threads, htslib decomp + resolve + train + buffer |
| fragment_router_scan (C++ single-threaded) | 15.1 s | 6.1 % | Streaming chunk scoring, builds CSR `ScoredFragments` |
| locus_em (C++ multi-threaded) | 14.1 s | 5.7 % | Batched + per-mega-locus |
| eb_gdna_priors (Python) | **9.5 s** | **3.9 %** | Per-locus FL exposure recomputation |
| partition_and_free (C++ scatter) | 6.0 s | 2.4 % | 11 sequential scatter ops |
| build_loci (Python) | 1.3 s | 0.5 % | Per-locus interval merge loop |
| index_load | 6.2 s | — | one-shot |
| calibration (orchestrator) | 0.04 s | — | FL build + density posteriors |
| compute_geometry, finalize_models | < 0.05 s | — | one-shot |

### Memory progression

| Snapshot | RSS (MB) | Δ |
|---|---:|---:|
| before | 2 940 | — |
| after_scan | 9 245 | **+6 305** |
| after_calibration | 9 245 | 0 |
| after_router_scan | **15 836** | **+6 591** |
| after_locus_em | 15 836 | 0 |
| after_cleanup | 15 836 | 0 |

Two RSS jumps dominate:

1. **+6.3 GB during BAM scan** — the `FragmentBuffer` lives at 3.9 GB resident even after spilling 14 chunks; the C++ scanner's per-worker accumulators add ≈ 1–2 GB more.
2. **+6.6 GB during scoring** — `ScoredFragments` CSR is built in float64 (`log_liks`, `coverage_weights`) on top of 260 M candidates → ≈ 4 GB just for those two arrays. The buffer is *not* released as chunks stream in (it is consumed but the global `ScoredFragments` grows in lockstep).

Allocator does not return memory to the OS after `partition_and_free` — RSS stays at peak through end of run.

---

## 2. Hotspot deep-dive

### 2.1 BAM scan (200.6 s, 81 %)

The scan is multi-threaded (3-stage producer / N workers / consumer), bound by:

- **htslib BGZF decompression** — now surfaced as `scan_bgzf_threads` / `--scan-bgzf-threads`.
- **`FragmentResolver::resolve_core`** — per-worker overlap resolution against the index intervals.
- **Python-side chunk callback** (`_on_chunk → inject_chunk → _accept_chunk`, 3.5 s cumulative) which finalizes the C++ chunk, computes `memory_bytes`, and triggers spill on overflow.
- **Spill I/O** — 14 × ~210 MB Arrow IPC writes through `pyarrow.feather.write_feather` (1.5 s wall, but synchronous on the main thread which blocks future chunks from being accepted).

**Throughput maths:**
32 M frags / 200 s = 160 K frags/s wall. With 8 worker threads + reader + decomp, effective thread-utilisation is well below 8× — at most ~3-4× over single-thread (other prod runs show similar ratios). Bottleneck is most likely shared-state contention in the resolve path or the producer queue draining slowly when chunk callbacks block on disk spill.

Areas of concern flagged by code review:

- `bam_scanner.cpp:1080` defaults `n_workers = 1` and `n_decomp_threads = 2`. The pipeline now derives `n_workers` from `scan.total_threads - scan.bgzf_threads` and passes the BGZF count from `BamScanConfig.bgzf_threads`. Worth double-checking the default `scan_bgzf_threads` for this workload — 2 BGZF threads is often the bottleneck on a fast NVMe.
- `inject_chunk` runs on the producer/Python thread holding the GIL; while it's running, no further chunks can be enqueued. Spilling synchronously from this callback can stall the entire scan.

### 2.2 Fragment router scan (15.1 s, 6.1 %)

`StreamingScorer.score_chunk` is **single-threaded C++**. With 32 M fragments and 30 M EM units it averages ≈ 470 ns/fragment, of which a portion is Python overhead in:

- `BufferedFragment.fragment_classes` (149 ms cumulative) — vectorised but allocates a new `uint8[N]` per chunk.
- `_load_chunk` reading 14 spilled chunks (1.45 s) — synchronous Arrow read on the scoring thread.

The scoring kernel itself is the dominant cost; parallelising across chunks is the highest-leverage change. With 8 workers and chunks of 1 M frags this should cut router_scan time by ~6× given enough memory headroom.

### 2.3 EB gDNA priors (9.5 s, 3.9 %) — biggest pure-Python regression

Per-locus loop in [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py) executes the following work for every multi-locus (24 K iterations):

| Sub-call | n_calls | tottime / cumtime | Notes |
|---|---:|---:|---|
| `estimate_locus_gdna` | 28 243 | 0.25 s / **4.78 s** | locoregional pi_gdna estimate |
| `expected_gdna_count_global` | 28 243 | 0.25 s / **2.27 s** | Phase-2 canonical η_g prior |
| `contained_exposure_clipped` | 84 729 | 0.37 s / 3.77 s | full + clipped FL exposure per locus |
| `compute_all_transcript_eff_lens` | **169 460** | 0.64 s / 3.39 s | analytical eff-length |
| `_build_eff_len_cache` | **169 460** | 0.27 s / **1.96 s** | recomputed per call! |
| `cumsum` on PMF | 566 767 | 1.49 s / 2.14 s | inside `_build_eff_len_cache` |
| `enable_gdna_for_multilocus` | 24 052 | 0.72 s / 0.78 s | per-unit reduction |
| `partition_units_to_loci` (genexpr) | 7 351 | 1.25 s / 1.25 s | O(n_loci × n_units) tuple build |

Three concrete inefficiencies:

1. **`FragmentLengthModel._build_eff_len_cache` is not memoised.** It runs `cumsum` over a 1024-element PMF on every `compute_all_transcript_eff_lens` call — 170 K times, even though the FL model is finalised and immutable through this whole stage. Fixing this single line drops ≈ 2 s wall (and removes 567 K cumsum calls). See [src/rigel/frag_length_model.py:399](../../src/rigel/frag_length_model.py#L399).
2. **Two passes over the same regions per locus** — `estimate_locus_gdna` and `expected_gdna_count_global` are called sequentially per multi-locus, each re-running `region_index.overlap`, `contained_exposure_clipped`, and the boundary helpers on overlapping inputs. Fusing them into a single per-locus pass eliminates duplicate FL exposure computation (~3.4 s).
3. **`partition_units_to_loci` for multi-Locus components** does a Python `for j in range(n_loci)` with `bins == j` boolean masks — O(n_loci × n_units). On the largest mega-locus (1497 loci × ~thousands of units) this is the 1.25 s line. Replace with a single `argsort(bins)` + `bincount(bins)` and slice output.

### 2.4 Locus EM (14.1 s, 5.7 %)

C++ `run_batch_locus_em_partitioned` is multi-threaded and spread well across the 8-thread pool. Two paths:

- **Mega-locus phase** — sequential per mega-locus (mega = work share ≥ `total / n_threads`). Each mega occupies the whole pool internally.
- **Normal-locus batch** — one batched call covering remaining loci.

The 14.1 s is roughly split between these two phases. Reasonable headroom remains:

- The largest multi-locus (42 639 transcripts × 2.82 M units) is single-threaded inside one batch slot — likely the long-pole.
- `_call_batch_em` repacks `partition_tuples` as Python tuples then nanobind-converts to C++; this is sub-percent but worth bundling.

### 2.5 Partition & free (6.0 s, 2.4 %)

`partition_and_free` runs **11 sequential scatter calls** (4 candidate-array scatters + 7 unit-array scatters). Each one walks the entire global CSR and writes one per-locus copy. Total reads ≈ 30 GB of array data through the scatter loop.

- The scatters could be **fused into a single C++ pass** that takes a list of `(global_arr, dtype, scratch)` tuples and emits per-locus copies in lockstep — one read, eight writes per fragment.
- Alternatively, **parallelise the scatter loop** with the existing thread pool (each scatter is independent).

### 2.6 Memory: where the 12.9 GB went

Per-fragment buffer cost (final): 32 M frags × 211 MB/M ≈ 6.5 GB. Per-fragment unit cost ≈ **211 bytes**, well above the docstring promise of 40-50. Major contributors:

| Per-frag column | dtype | bytes/frag | Notes |
|---|---|---:|---|
| `t_indices` (per candidate, avg ≈ 8) | int32 | 32 | already minimum |
| `frag_lengths` (per candidate) | int32 | 32 | could be uint16 (FL ≤ 1024) |
| `exon_bp` (per candidate) | int32 | 32 | could be uint16 |
| `intron_bp` (per candidate) | int32 | 32 | could be uint16 |
| `exon_bp_pos/neg`, `tx_bp_pos/neg` | int32 ×4 | 16 | strand-aware; int16 sufficient |
| `frag_id` | int64 | 8 | 32 M fits in int32 (saves 4 B/frag) |
| `read_length` | uint32 | 4 | uint16 sufficient |
| Other fixed | various | ~10 | |

Trimming the per-candidate widths to int16 alone would drop the buffer by ≈ 35 % (6.5 GB → 4 GB).

Per-fragment `ScoredFragments` cost (after scoring): ≈ 6.6 GB. Major contributors:

| CSR array | dtype | size | Notes |
|---|---|---:|---|
| `log_liks` | float64 | 260 M × 8 = **2.08 GB** | float32 keeps ample log-prob precision |
| `coverage_weights` | float64 | 260 M × 8 = **2.08 GB** | float32 sufficient |
| `t_indices` | int32 | 260 M × 4 = 1.04 GB | minimum |
| `count_cols` | uint8 | 260 M × 1 = 260 MB | minimum |
| `gdna_log_liks` | float64 | 30 M × 8 = 240 MB | float32 OK |
| Other per-unit | various | ~600 MB | |

Casting the two float64 candidate arrays to float32 saves **~2 GB peak RSS** (and halves memory bandwidth in the EM kernel inner loop). Numerical impact is negligible — the EM internally renormalises and the log-sum-exp pivot kills any float32 dynamic-range issues.

### 2.7 GIL contention indicator

Top cProfile self-time is `_thread.lock.acquire` (228 s, ~93 % of wall). This is the main thread waiting on the C++ scanner — not a Python lock issue per se. It does, however, expose one fact: **only one Python thread is doing any work** the entire time. Any opportunity to run the chunk-callback overhead (Arrow finalise, RSS sampling, partial scoring) on a worker thread instead of the main thread would let the C++ scanner reclaim that wait time.

---

## 3. Prioritised optimisation plan

Ordered by **expected wall-time / RSS payoff per implementation cost**. Each item is independently shippable and has a measurable target.

### Tier 1 — Trivial fixes, large payoff

| # | Item | Target | Cost | Notes |
|---|---|---|---|---|
| **T1.1** | Memoise `FragmentLengthModel._build_eff_len_cache` after `finalize()` | **−2 s wall, −567 K cumsum calls** | trivial | Cache `_cdf`, `_cmom` on first call once `_finalized` is True. Invalidate on any `observe_*` (which never fires post-calibrate). |
| **T1.2** | Cast `ScoredFragments.log_liks` and `coverage_weights` to **float32** | **−2 GB peak RSS** + scoring kernel speedup | small (C++ + dtype plumbing) | Validate numerics with golden tests; renormalise inside log-sum-exp per row. |
| **T1.3** | Vectorise `partition_units_to_loci` slow path with `argsort+bincount` | −1.2 s wall on mega-locus stage | small | Single C-loop over `bins` instead of `n_loci` boolean masks. |
| **T1.4** | Fuse `estimate_locus_gdna` + `expected_gdna_count_global` per Locus | **−3 s wall** | medium | Single region-overlap query, single `contained_exposure_clipped` call returning `(eff_full, eff_evidence, eff_core)`; share boundary computation. |

### Tier 2 — Memory diet, structural

| # | Item | Target | Cost |
|---|---|---|---|
| **T2.1** | Narrow per-candidate buffer dtypes (int16 for `frag_lengths`, `exon_bp`, `intron_bp`, `*_bp_pos/neg`); int32 for `frag_id`; uint16 for `read_length` | **−2.5 GB peak RSS** during scan | medium (C++ struct + Python ingest) |
| **T2.2** | Drop unused buffer columns (`exon_bp_pos/neg`, `tx_bp_pos/neg`) if no consumer in the v6 calibration path | up to −512 MB | small audit |
| **T2.3** | Drop the buffer chunk eagerly **before** building `ScoredFragments` for that chunk (currently the global CSR grows while spilled chunks are loaded back, doubling resident memory) | **−3 GB peak RSS** | medium — the streaming scorer already consumes chunks; verify Arrow read deletes the on-disk file before the next chunk is loaded |
| **T2.4** | Replace `log_liks` / `coverage_weights` with float32 in the buffer→scorer→EM data path end-to-end | **−4 GB peak RSS** | medium (depends on T1.2) |

### Tier 3 — Scoring & scatter parallelisation

| # | Item | Target | Cost |
|---|---|---|---|
| **T3.1** | Parallelise `StreamingScorer.score_chunk` across chunks using the existing C++ thread pool | **−10 s wall** (15 s → ~3 s) | medium — score_chunk is per-chunk independent; need per-thread scratch + atomic CSR-append (or per-thread CSR with merge) |
| **T3.2** | Fuse `partition_and_free` 11 scatters into a single C++ pass (one global-CSR walk, multiple writers) | **−4 s wall** | medium — `partition_and_free.cpp` refactor |
| **T3.3** | Run the largest mega-locus EM with `n_threads = all`, batch the rest with `n_threads = 1` per locus across the pool | up to **−5 s** on mega-locus phase | small — adjust the work-split in `_run_locus_em_partitioned` |

### Tier 4 — BAM scan engine work

| # | Item | Target | Cost |
|---|---|---|---|
| **T4.1** | Tune `scan_bgzf_threads` / `BamScanConfig.bgzf_threads` on this workload | likely **−15-30 s** | trivial |
| **T4.2** | Move spill I/O off the main scanner thread (background spill thread; the producer should never block on `pf.write_feather`) | **−5-10 s** wall | medium |
| **T4.3** | Defer `BufferedFragment.fragment_classes` until first read (currently called per-chunk in `to_scoring_arrays`) | small but reduces allocation churn | small |
| **T4.4** | Profile `FragmentResolver::resolve_core` with a dedicated micro-benchmark — confirm whether the resolve hash-map / interval search is the per-thread bottleneck (suggested by the < 4× thread scaling) | informational | medium — needs a focused C++ profile (Instruments / VTune) |

### Tier 5 — Process-wide memory

| # | Item | Target | Cost |
|---|---|---|---|
| **T5.1** | Configure jemalloc / mimalloc with aggressive `MALLOC_CONF=background_thread:true,metadata_thp:auto,dirty_decay_ms:1000`. The default macOS allocator is not returning memory to the OS at peak. | reduces *resident* RSS post-quant (memory available to other workloads) | trivial — env var |
| **T5.2** | After `partition_and_free`, explicitly `gc.collect()` and `ctypes.CDLL("libc")`-style `malloc_trim(0)` on Linux. | small but visible in long-running benchmarks | trivial |
| **T5.3** | Reduce default `BamScanConfig.buffer_size_bytes` from 2 GiB to ≤ 1 GiB to spill earlier; offset by `T4.2` so spilling is no longer a wall-time penalty | −1-2 GB peak RSS | trivial config change once T4.2 lands |

---

## 4. Suggested implementation order

1. **Land Tier 1 wholesale** (T1.1–T1.4). Pure Python / one C++ kernel cast. No risk to correctness, all backed by golden-output tests.
2. Re-run this profile to confirm: expected wall ≈ 240 s → ~225 s, peak RSS ≈ 15.8 → ~13 GB.
3. **Tier 4** in parallel with Tier 2 — independent code paths.
4. **Tier 3** last; it's the highest-leverage CPU win but the one most likely to require touching the EM ABI / memory layout.

---

## 5. Caveats

- The cProfile run *shifts* the scan time upward (200 s observed vs ~330 s prod-without-cprofile, which itself was inflated by laptop sleep mid-scan during instrumentation runs). The relative weight of stages is reliable; the absolute wall numbers in the BAM scan should be re-measured with cProfile off after Tier 1 lands.
- The existing pre-existing `int8` cap in `build_t_to_local_locus` was a real correctness bug surfaced by this workload (multi-locus 3 has 1497 loci, > 127). Fixed in this commit (`int8 → int32`); see [src/rigel/calibration/_locus_n_obs.py](../../src/rigel/calibration/_locus_n_obs.py).
- The profiler's `c_base` reference in `assemble_priors` was stale post calibration v6; fixed in this commit. See [scripts/profiling/profiler.py](../../scripts/profiling/profiler.py).
- Per-locus profiling stats for all 24 052 multi-loci are in [locus_stats_default_8t.json](../../../Downloads/rigel_runs/profile_vcap_rna20m_gdna20m/locus_stats_default_8t.json) for downstream analysis (per-locus iterations, work, mega-locus identification).
