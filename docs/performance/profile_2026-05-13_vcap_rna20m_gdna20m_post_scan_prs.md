# Rigel Post-Scan-PR Performance Profile
## VCAP RNA20M + gDNA20M after PR01-PR04 and cleanup

**Profile date:** 2026-05-13  
**Workload:** `/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam`  
**Index:** `/Users/mkiyer/Downloads/rigel_runs/refs/rigel_index`  
**Platform:** macOS 26.4.1, Apple Silicon arm64, Python 3.12.13, conda env `rigel`  
**Run shape:** 66.1 M BAM records, 31.99 M read names, 32.46 M buffered scan records, 30.39 M ambiguous EM units, 24,052 multi-loci, 457,371 transcripts.

Source profile bundles:

- Clean staged profile, no cProfile overhead: `/Users/mkiyer/Downloads/rigel_runs/profile_vcap_rna20m_gdna20m_post_pr04_nocprofile`
- Staged profile with cProfile: `/Users/mkiyer/Downloads/rigel_runs/profile_vcap_rna20m_gdna20m_post_pr04`
- Scan decompression sweep: `/Users/mkiyer/Downloads/rigel_runs/scan_decomp_sweep_post_pr04`
- Scan memory-cap sweeps: `/Users/mkiyer/Downloads/rigel_runs/scan_memory_1gib_post_pr04`, `/Users/mkiyer/Downloads/rigel_runs/scan_memory_2gib_post_pr04`

Implementation roadmap: [post_scan_prs/README.md](post_scan_prs/README.md)

The clean staged run is used for stage times and RSS. The cProfile run is used only for Python attribution; it added about 4.5 s total and especially inflated `fragment_router_scan` and prior assembly.

---

## 1. Executive Summary

PR01-PR04 achieved the intended dramatic scan improvement. The old profile was dominated by `scan_and_buffer`; the new profile is a balanced native-array pipeline where scan, scoring, EM, scatter, and priors are all visible.

| Metric | Initial profile | Current clean profile | Change |
|---|---:|---:|---:|
| Total wall time | 246.8 s | **73.36 s** | **3.36x faster** |
| Throughput | 268 K frags/s | **901 K frags/s** | **3.36x higher** |
| `scan_and_buffer` | 200.6 s | **32.78 s** | **6.12x faster** |
| Quant after scan | ~46.2 s | **40.55 s** | 1.14x faster |
| Peak RSS | 15.84 GB | **16.20 GB** | **+0.36 GB** |

The key conclusion is a little bracing but useful: wall time improved by more than 3x, but memory did not improve. Peak RSS is still about 16 GB, and the largest remaining opportunity is to reduce the amount of data materialized between the buffer, global scoring CSR, per-locus partitions, and EM subproblems.

The next highest-yield work is:

1. Lower the scan buffer memory cap now that async spill is cheap, using the canonical `--scan-buffer-size` / `scan_buffer_size` parameter.
2. Convert scoring/EM likelihood payloads from float64 to float32 or mixed precision where safe.
3. Fuse or parallelize `partition_and_free` scatter.
4. Parallelize `StreamingScorer.score_chunk` across chunks.
5. Narrow/drop buffer columns that are not consumed after calibration.

---

## 2. Current Stage Breakdown

Clean staged run from before the naming cleanup: 8 native scan workers, 4 BGZF threads, 8 EM threads, 4 GiB scan buffer cap, and read-name batch size 512. Under the new total-budget convention, the same scan worker split is `--threads 12 --scan-bgzf-threads 4`; future profiles should report both the requested total budget and the resolved scan worker count.

| Stage | Current wall | Current % | Initial wall | Speedup | Notes |
|---|---:|---:|---:|---:|---|
| `scan_and_buffer` | **32.776 s** | **44.7%** | 200.6 s | **6.12x** | PR01-PR04 landed; still largest single stage |
| `fragment_router_scan` | **14.717 s** | **20.1%** | 15.1 s | 1.03x | Single-thread native scoring/CSR builder |
| `locus_em` | **14.303 s** | **19.5%** | 14.1 s | 0.99x | C++ SQUAREM EM, already threaded |
| `partition_and_free` | **6.154 s** | **8.4%** | 6.0 s | 0.97x | 11 sequential scatter passes over global CSR |
| `compute_eb_gdna_priors` | **4.230 s** | **5.8%** | 9.5 s | **2.25x** | Previous cleanup/prior work paid off |
| `build_loci` | **1.140 s** | **1.6%** | 1.3 s | 1.14x | Native connected components plus Python assembly |
| Calibration/finalize/geometry | **0.039 s** | ~0% | ~0.1 s | n/a | Not performance relevant |

The post-scan quantification section is now 40.55 s. Its main costs are scoring, EM, scatter, and priors; none is individually as dominant as the old scan bottleneck.

---

## 3. Memory Profile

RSS progression in the clean full profile:

| Snapshot | RSS | Delta | Interpretation |
|---|---:|---:|---|
| before | 2.92 GB | - | index loaded, baseline process state |
| after scan | **9.81 GB** | **+6.89 GB** | resident `FragmentBuffer` plus native/Python allocator retention |
| after calibration | 9.81 GB | 0 | calibration itself is negligible |
| after router scan | **16.20 GB** | **+6.38 GB** | global `ScoredFragments` CSR built on top of prior allocation high-water |
| after buffer release | 16.20 GB | 0 | live objects change, macOS allocator does not return RSS |
| after locus EM | 16.20 GB | 0 | EM does not exceed scoring high-water |
| after cleanup | 16.20 GB | 0 | process RSS remains at peak |

Peak RSS occurs at about 46.6 s, shortly after scoring/router scan finishes. This is the same qualitative pattern as the initial report: scan creates the first large resident plateau; scoring creates the second.

### Scan-only memory-cap sweep

With 8 native scan workers and 2 BGZF threads (equivalent to `--threads 10 --scan-bgzf-threads 2` under the new total-budget convention), async spill makes lower memory caps nearly free in wall time:

| Buffer cap | Scan wall | Peak RSS | Resident buffer | Spilled chunks | Read names/s |
|---:|---:|---:|---:|---:|---:|
| 4 GiB | 32.11 s | 9.66 GB | 3.90 GB | 14 | 996 K/s |
| 2 GiB | **31.21 s** | 7.88 GB | 2.00 GB | 23 | **1.025 M/s** |
| 1 GiB | 31.97 s | **7.07 GB** | **0.94 GB** | 28 | 1.001 M/s |

This is the most straightforward memory lever now available. Lowering the default scan buffer cap from 4 GiB to 2 GiB should reduce scan RSS by roughly 1.8 GB on this workload with no observed scan penalty. A 1 GiB cap saves another ~0.8 GB and still shows no meaningful wall penalty in scan-only mode. PR05 now exposes scan memory/decompression parameters in the staged profiler, so future full-profile cap sweeps can be run directly from CLI flags.

### Why memory remains high

The current layout materializes several large representations in sequence:

- `_FinalizedChunk` arrays in [src/rigel/buffer.py](../../src/rigel/buffer.py): per-fragment columns plus per-candidate CSR columns. Logical buffer size is about 6.8 GB for this run; 3.9 GB remains resident at the 4 GiB cap.
- `StreamingScorer` in [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp): growing `std::vector<double>` payloads for `log_liks`, `coverage_weights`, and `gdna_log_liks` plus int arrays.
- `ScoredFragments` in [src/rigel/scored_fragments.py](../../src/rigel/scored_fragments.py): global CSR arrays are held until prior assembly and partitioning.
- `partition_and_free` in [src/rigel/locus_partition.py](../../src/rigel/locus_partition.py): scatters global arrays into per-locus copies, then nulls global arrays. This lowers live Python references but does not reduce RSS once the allocator has reached the high-water mark.

The two candidate float64 arrays are especially expensive. At the initial report's roughly 260 M RNA candidates, `log_liks` and `coverage_weights` alone account for about 4.2 GB. Moving them to float32 would save about 2.1 GB before considering partition copies and memory-bandwidth wins.

---

## 4. Focused Scan Findings

The scan PRs worked. Scan-only performance is now about 32 s for the full 66 M-record BAM, or about 1.0 M read-name groups/s.

### Decompression sweep

At 8 scan workers and a 4 GiB buffer cap:

| `scan_bgzf_threads` | Scan wall | Peak RSS | Read names/s |
|---:|---:|---:|---:|
| 2 | **32.11 s** | 9.66 GB | **996 K/s** |
| 4 | 33.18 s | 9.86 GB | 964 K/s |
| 8 | 33.69 s | 9.86 GB | 949 K/s |

More BGZF decompression threads do not help on this M3 run. The public parameter should be `--scan-bgzf-threads` / `scan_bgzf_threads`, and `--threads` should remain the total budget from which scan workers are derived. The default should be chosen intentionally rather than drifting between Python and C++.

### Remaining scan interpretation

The scan is no longer obviously decompression-bound. Increasing decompression threads slows it down, while lowering spill memory does not hurt. That points toward a remaining mix of resolver CPU, chunk finalization, and memory bandwidth. The old queue-traffic bottleneck is gone; async spill is not the limiting factor at these caps.

Recommended scan follow-up:

- Keep `scan_bgzf_threads` visible in CLI/YAML/profiler outputs and derive scan workers from the total `threads` budget.
- Consider changing the default from 4 to 2 on macOS/arm64, or use an explicit auto policy only after broader data.
- Keep `scan_buffer_size` / `--scan-buffer-size` visible in profile/CLI docs.
- Run future C++ scan sampling only after the structural memory/scoring items, because scan is no longer the only long pole.

---

## 5. Fragment Router / Scoring

`fragment_router_scan` is now the largest non-scan wall-time target: 14.72 s, about 20% of the full run. It is almost unchanged from the initial profile because the scan PRs did not target this path.

The current flow in [src/rigel/scan.py](../../src/rigel/scan.py) streams chunks from the buffer and calls `StreamingScorer.score_chunk` in [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp). The scoring kernel itself is native and not usefully attributed by cProfile; the stage timer is the reliable number.

Important details:

- Scoring is single-threaded across chunks.
- Each chunk is independent from Python's perspective: `iter_chunks_consuming()` loads/yields one chunk and releases it after scoring.
- `StreamingScorer` appends to monolithic output vectors. This is simple but serializes all chunk work and creates the large global CSR high-water.
- Output payloads use double precision for likelihoods and coverage weights even though most downstream operations are log-sum-exp normalized.

Best opportunities:

1. **Parallel chunk scoring with per-worker CSR outputs.** Each worker scores a chunk into local vectors; a final merge computes global offsets and concatenates arrays. This should reduce the 14.7 s router stage substantially, likely into the 3-6 s range depending on memory bandwidth. The design must preserve multimapper grouping assumptions and deterministic output order.
2. **Float32 likelihood/coverage payloads.** Use float32 for `log_liks`, `coverage_weights`, and probably `gdna_log_liks` until the EM inner loop needs double accumulation. This saves memory immediately and may speed scoring, scatter, and extraction through lower bandwidth pressure.
3. **Emit per-locus partitions directly or in batches.** A more structural option is to avoid building one global CSR at all: score chunks into component/locus-local buffers after connected components are known, or run a two-pass scoring where the first pass only builds connectivity. This is higher risk but attacks both memory and scatter.

---

## 6. Locus EM

`locus_em` is 14.30 s, essentially unchanged from the initial profile. The locus stats show that this is not only the giant mega-locus.

| Locus set | CPU-like summed time | Share of summed locus time |
|---|---:|---:|
| all loci with stats | 83.94 s | 100% |
| mega-locus only | 3.80 s | 4.5% |
| top 10 loci | 8.30 s | 9.9% |
| top 100 loci | 21.61 s | 25.7% |

The largest locus has 42,639 transcripts and 2,823,632 units. It uses 8 E-step threads and takes 3.80 s total, with 2.31 s in SQUAREM, 0.82 s building equivalence classes, and 0.48 s assigning posteriors. Several much smaller loci take 0.3-0.7 s because they require many SQUAREM iterations.

Summed locus substage instrumentation:

| Substage | Summed time | Interpretation |
|---|---:|---|
| `squarem_us` | 72.76 s | dominant CPU work, parallelized across loci/threads |
| `assign_us` | 4.66 s | posterior assignment cost after convergence |
| `build_ec_us` | 4.37 s | equivalence-class construction/extraction |
| `extract_us` | 1.87 s | partition to local subproblem extraction |
| `warm_start_us` | 0.29 s | negligible |

The EM solver is not the first place I would spend effort because it is already parallel and correctness-sensitive. Still, useful follow-ups exist:

- Investigate loci that hit high SQUAREM iteration counts; many smaller loci reach 333 reported iterations. Do not change tolerances casually, but inspect whether active-set pruning, better warm starts, or early detection of effectively deterministic ECs can reduce iterations without altering model semantics.
- Consider splitting the execution policy more carefully: mega-locus uses all threads internally, while the normal batch uses one thread per locus. The current split is sensible, but the top-100 distribution suggests scheduling and active-set improvements may matter more than special-casing the single mega-locus.
- Avoid repeated local extraction/sorting where possible. The extraction path in [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) remaps and sorts candidates by local component for each locus; fused partition/extraction could remove one intermediate representation.

---

## 7. Partition Scatter

`partition_and_free` remains about 6.15 s. It is a clean, isolated target.

The current implementation in [src/rigel/locus_partition.py](../../src/rigel/locus_partition.py) performs:

- one `build_partition_offsets` pass,
- four candidate-array scatters (`log_liks`, `coverage_weights`, `t_indices`, `count_cols`),
- seven unit-array scatters (`gdna_log_liks`, `locus_t_indices`, `locus_count_cols`, `is_spliced`, `frag_ids`, `frag_class`, `splice_type`).

The native scatter helpers in [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) are sequential and walk the same locus/unit structure repeatedly. This is memory-bandwidth heavy and creates per-locus copies for every array.

Recommended work:

1. **Fuse scatter into one native pass per locus.** For each locus, allocate all destination arrays, walk its units once, and copy all candidate/unit fields together. This reduces repeated traversal and should improve cache behavior.
2. **Parallelize scatter across loci.** Loci are independent after offsets are known. A thread pool over loci should be straightforward, with special handling for the mega-locus.
3. **Combine with float32 payload work.** If `log_liks` and `coverage_weights` become float32, scatter bandwidth and partition memory drop immediately.

Expected payoff: 2-4 s wall and lower memory-bandwidth pressure, with relatively low model risk.

---

## 8. gDNA Prior Assembly

`compute_eb_gdna_priors` improved from 9.5 s to 4.23 s. This is a good result from the cleanup work: effective-length caches are memoized in [src/rigel/frag_length_model.py](../../src/rigel/frag_length_model.py), unit-to-locus binning is vectorized in [src/rigel/calibration/_locus_n_obs.py](../../src/rigel/calibration/_locus_n_obs.py), and the prior loop shares per-locus scratch in [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py).

Remaining cProfile-visible costs in the prior path:

| Function/path | cProfile cumulative | Notes |
|---|---:|---|
| `_compute_locus_scratch` | 2.57 s | overlap plus FL-aware exposure setup |
| `estimate_locus_gdna` | 1.30 s | diagnostic locoregional estimate |
| `contained_exposure_clipped` | 1.25 s | still called many times, now with cached FL CDF/moment arrays |
| `partition_units_to_loci` | 1.09 s | vectorized now, still nonzero on complex multi-locus cases |
| `compute_all_transcript_eff_lens` | 0.99 s | no longer cumsum-bound, but still many vector calls |
| `enable_gdna_for_multilocus` | 0.66 s | per-locus unit gather/reduction |

Recommended work:

- Add a production fast path that computes only the EM-consumed global `eta_g` and eligibility when per-locus diagnostic decomposition is not requested. The current diagnostics are valuable, but `estimate_locus_gdna` is not consumed by EM.
- Compute `enable_gdna` during native partition/scatter or scoring rather than as a separate Python gather over `em_data` per multi-locus.
- Batch region overlaps/exposure computation where possible. This is lower priority than scoring/scatter because the whole stage is now only 4.23 s.

---

## 9. Prioritized Intervention Roadmap

This ranking balances expected yield, implementation risk, and straightforwardness.

### P0: Small, high-confidence configuration fixes

| Item | Expected payoff | Risk | Why now |
|---|---:|---|---|
| Lower default scan cap from 4 GiB to 2 GiB, or at least profile full quant at 2 GiB | 1.5-2.5 GB lower RSS on this workload | Low | Scan-only sweep shows no wall penalty; async spill changed the tradeoff |
| Keep `scan_bgzf_threads` visible in CLI/YAML/docs and align Python/C++ defaults | ~1 s scan on this run; better reproducibility | Low | Current default drift is confusing; d2 beats d4/d8 on M3 |
| Teach the full profiler to accept/report all scan params and `em_data.n_candidates` | measurement quality | Low | Needed to validate memory-cap changes end-to-end |

### P1: Highest memory payoff

| Item | Expected payoff | Risk | Notes |
|---|---:|---|---|
| Float32 scoring payloads (`log_liks`, `coverage_weights`, `gdna_log_liks`) | ~2+ GB peak RSS, faster bandwidth-bound stages | Medium | Requires native scoring/scatter/EM ABI changes and golden validation |
| Narrow buffer dtypes (`frag_lengths`, `exon_bp`, `intron_bp`, `read_length`) | 1-2+ GB logical buffer reduction | Medium | Respect `max_frag_length`; keep overflow checks/fallbacks |
| Drop or late-materialize unused buffer columns (`exon_bp_pos/neg`, `tx_bp_pos/neg`) | hundreds of MB to ~0.5 GB | Low/Medium | Audit suggests these are stored/spilled but not consumed by scoring |

### P2: Highest wall-time payoff outside scan

| Item | Expected payoff | Risk | Notes |
|---|---:|---|---|
| Parallel chunk scoring | 8-10 s wall if near-linear enough | Medium/High | Per-worker CSR plus deterministic merge; preserve multimapper semantics |
| Fused/parallel partition scatter | 2-4 s wall | Medium | Isolated native refactor, pairs well with float32 |
| Optional fast prior mode | 1-2 s wall | Low/Medium | Skip diagnostic locoregional prior details when not requested |

### P3: EM algorithm refinement

| Item | Expected payoff | Risk | Notes |
|---|---:|---|---|
| Analyze high-iteration non-mega loci | 1-5 s possible | Medium | Correctness-sensitive; no tolerance weakening without root cause |
| Active-set / deterministic-EC pruning | unknown but plausible | Medium/High | Could reduce SQUAREM work and posterior assignment |
| Fuse partition extraction with EM subproblem build | 1-3 s and less memory churn | Medium | More structural but aligns with scatter refactor |

---

## 10. Suggested Next PR Sequence

1. **PR05: Profiling/config visibility.** Use the canonical scan parameter names in the profiling harness and write `n_candidates`, candidate bytes, partition bytes, and scan config into profile JSON. This is small but prevents blind spots.
2. **PR06: Memory cap default and spill validation.** Run full profiles at 4/2/1 GiB once the profiler can set scan params. If 2 GiB remains neutral, change the default or document a recommended production profile.
3. **PR07: Float32 scoring payload.** Convert scoring output and partition payloads to float32 while preserving double accumulation inside EM where needed. Validate with full golden suite and a VCAP numerical diff.
4. **PR08: Fused partition scatter.** Replace 11 sequential scatter calls with one native fused scatter, then optionally parallelize across loci.
5. **PR09: Parallel scoring.** Add per-worker chunk scoring with deterministic CSR merge. This is the biggest remaining wall-time move, but it should come after the float32/memory work so the parallel scorer does not simply move the bottleneck to memory bandwidth.
6. **PR10: Buffer dtype diet.** Narrow safe columns and drop/late-materialize unused columns. This can run before PR09 if memory pressure blocks parallel scoring.
7. **PR11: EM convergence/active-set investigation.** Use per-locus stats to target high-iteration loci and build specific acceptance tests before changing solver behavior.

---

## 11. Bottom Line

The old problem was a scan bottleneck. The new problem is data movement.

The pipeline now scans quickly enough that the memory layout and global CSR lifecycle dominate the next round of performance work. The fastest low-risk improvement is to exploit async spill by lowering scan memory, then fix parameter visibility. The biggest structural gains are float32 scoring payloads, fused scatter, and parallel scoring. EM remains important but should be treated carefully because it is already threaded and correctness-sensitive.