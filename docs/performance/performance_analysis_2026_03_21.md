# Rigel Performance Analysis — March 21, 2026

## Executive Summary

Profiling of the rigel pipeline reveals three dominant bottlenecks that together account for **94% of wall-clock time**:

1. **gDNA Calibration** (38% large / 54% small) — Python-level `math.lgamma` via `np.frompyfunc` is 50–100× slower than `scipy.special.gammaln`; a Python for-loop in `_compute_fl_llr` processes 11M fragments at O(1) per iteration
2. **Locus EM** (28% large / 10% small) — One mega-locus with 226K transcripts and 12.7M units dominates; C++ solver is well-optimized but the sheer problem size is the bottleneck
3. **BAM Scan** (28% large / 15% small) — Already C++-optimized; limited room for improvement

Memory usage is also concerning: **21 GB peak RSS** for the large BAM (14.7M fragments), with a substantial jump during locus EM.

---

## Phase 1 Results: Before / After

Phase 1 optimizations (items 1–4) were implemented in `calibration.py` and validated against the full test suite (1,037 tests pass). Changes:

1. **`_vec_lgamma` / `_betaln_vec`** → `scipy.special.betaln` (C-compiled vectorized, replaces `np.frompyfunc(math.lgamma)`)
2. **`_compute_fl_llr` loop** → vectorized `np.minimum` + `np.add.at` (eliminates 145M `builtins.min` calls)
3. **`build_gdna_fl_model` loop** → direct array assignment to `model.counts`
4. **Removed unused** `_lgamma_vec`, `_vec_lgamma`, old `_betaln_vec` function

### Simulated BAM (large, 14.7M fragments)

| Stage | Before (s) | After (s) | Speedup | Notes |
|-------|-----------|----------|---------|-------|
| scan_and_buffer | 108.1 | 109.4 | — | Unchanged (C++) |
| **calibration** | **149.6** | **17.8** | **8.4×** | 12 iterations, identical convergence |
| quant_from_buffer | 135.1 | 134.3 | — | Unchanged |
| locus_em | 111.0 | 110.6 | — | Unchanged (C++) |
| **Total** | **392.8** | **261.5** | **1.50×** | **33% wall-time reduction** |
| Peak RSS | 21,339 MB | 20,956 MB | — | ~400 MB less (fewer temp arrays) |
| Throughput | 98,500 frags/s | 147,900 frags/s | 1.50× | |

### Real BAM (small, 309K fragments)

| Stage | Before (s) | After (s) | Speedup | Notes |
|-------|-----------|----------|---------|-------|
| scan_and_buffer | 2.5 | 2.5 | — | Unchanged |
| **calibration** | **8.9** | **1.4** | **6.4×** | 50 iterations |
| quant_from_buffer | 5.0 | 4.8 | — | Unchanged |
| locus_em | 1.7 | 1.4 | — | Unchanged |
| **Total** | **16.4** | **8.7** | **1.89×** | **47% wall-time reduction** |
| Peak RSS | 3,943 MB | 3,962 MB | — | Unchanged |
| Throughput | 149,000 frags/s | 281,000 frags/s | 1.89× | |

### Analysis

- **Calibration is no longer the dominant bottleneck.** For the large BAM, it dropped from 38% of wall time to 7%. For the small BAM, from 54% to 16%.
- **New bottleneck ranking (large BAM):** locus_em (42%) > scan_and_buffer (42%) > calibration (7%)
- **New bottleneck ranking (small BAM):** quant_from_buffer (55%) > scan_and_buffer (29%) > calibration (16%)
- Calibration convergence values are **numerically identical** before and after — the optimization is pure speedup with no accuracy change.
- The remaining calibration time (~17.8s on large BAM) is now dominated by `_compute_strand_llr_betabinom` (0.35s/iter × 12 = 4.2s), `_compute_fl_llr` vectorized (0.7s/iter × 12 = 8.4s), and `_golden_section_max` κ search.

## Phase 2 Results: Before / After

Phase 2 optimizations (items 5–7) applied on top of Phase 1:

1. **`_golden_section_max` → Brent's method** via `scipy.optimize.minimize_scalar(method='bounded')` — 34% fewer evaluations (251 vs ~455)
2. **Exon tuple construction** → `tuple(.tolist())` + `np.cumsum()` replacing `tuple(int(x) for x in ...)` generator expressions
3. **Cache `fragment_classes`** → lazy compute + `object.__setattr__` cache in `_FinalizedChunk`

### Simulated BAM (large, 14.7M fragments)

| Stage | Baseline | Phase 1 | Phase 2 | Speedup vs Baseline |
|-------|----------|---------|---------|---------------------|
| scan_and_buffer | 108.1s | 109.4s | 109.8s | — |
| **calibration** | **149.6s** | **17.8s** | **11.7s** | **12.8×** |
| fragment_scorer | 1.5s | 1.5s | 1.4s | 1.1× |
| fragment_router_scan | 18.6s | 18.3s | 18.1s | — |
| locus_em | 111.0s | 110.6s | 109.3s | — |
| **Total** | **392.8s** | **261.5s** | **254.3s** | **1.54×** |
| Peak RSS | 21,339 MB | 20,956 MB | 20,788 MB | -551 MB |
| Throughput | 98,500 | 147,900 | 152,100 | 1.54× |

### Real BAM (small, 309K fragments)

| Stage | Baseline | Phase 1 | Phase 2 | Speedup vs Baseline |
|-------|----------|---------|---------|---------------------|
| scan_and_buffer | 2.5s | 2.5s | 2.5s | — |
| **calibration** | **8.9s** | **1.4s** | **1.2s** | **7.4×** |
| fragment_scorer | 1.6s | 1.5s | 1.4s | 1.1× |
| locus_em | 1.7s | 1.4s | 1.5s | — |
| **Total** | **16.4s** | **8.7s** | **8.4s** | **1.95×** |

### Analysis

- **Calibration further reduced** by 34% vs Phase 1 (17.8s → 11.7s) from Brent's method.
- Brent's method converges in ~19 evals/call vs ~35 for golden-section (251 total vs ~455).
- Remaining calibration time dominated by `_marginal_loglik` (6.5s, 251 evals × 26ms each).
- Fragment scorer improved modestly (1.5s → 1.4s) from vectorized exon construction.
- Total function calls dropped from ~160M (baseline) to 9.25M (Phase 2) — **17× reduction**.
- **New bottleneck:** scan_and_buffer (43%) and locus_em (43%) are both C++ stages — Python-level optimization has reached diminishing returns.

---

## Test Data

| Dataset | BAM Size | Fragments | Buffered | Transcripts | Genes |
|---------|----------|-----------|----------|-------------|-------|
| **Simulated (large)** | 1.7 GB | 38.7M reads → 14.7M frags | 18.6M | 254,461 | 63,472 |
| **Real (small)** | 187 MB | 2.4M reads → 309K frags | 453K | 254,461 | 63,472 |

Both use the same human reference index (GRCh38, hg38 chromosome naming).

---

## Stage-Level Timing Breakdown

### Simulated BAM (large, 14.7M fragments)

| Stage | Time (s) | % | Notes |
|-------|----------|---|-------|
| scan_and_buffer | 108.1 | 27.5% | C++ htslib scan + resolve + buffer |
| calibration | 149.6 | 38.1% | gDNA calibration EM (12 iterations) |
| fragment_scorer | 1.5 | 0.4% | One-time setup (exon data) |
| fragment_router_scan | 18.6 | 4.7% | C++ fused scoring + routing |
| build_loci | 2.4 | 0.6% | C++ union-find connected components |
| eb_gdna_priors | 1.5 | 0.4% | cgranges overlap queries |
| **locus_em** | **111.0** | **28.3%** | Batch C++ EM solver |
| **Total** | **392.8** | 100% | **98,500 frags/s throughput** |

### Real BAM (small, 309K fragments)

| Stage | Time (s) | % | Notes |
|-------|----------|---|-------|
| scan_and_buffer | 2.5 | 15.1% | |
| calibration | 8.9 | 54.1% | 50 iterations (slow convergence) |
| fragment_scorer | 1.6 | 9.5% | Same one-time cost regardless of BAM size |
| fragment_router_scan | 0.2 | 1.2% | |
| build_loci | 0.8 | 4.9% | |
| eb_gdna_priors | 0.8 | 4.6% | |
| **locus_em** | **1.7** | **10.4%** | |
| **Total** | **16.4** | 100% | **149,000 frags/s throughput** |

### Key Observations

- **Calibration dominates** in both datasets (38–54% of wall time)
- **Locus EM** scales super-linearly with problem size (1.7s for 130K units vs 111s for 12.7M units — 97× more units → 65× more time)
- **fragment_scorer** has a **fixed ~1.5s cost** independent of BAM size (one-time per-transcript exon precomputation for 254K transcripts)
- Index loading takes **8s** (not included in profiling timers but significant for small jobs)

---

## Memory Profile

### Simulated BAM (large)

| Checkpoint | RSS (MB) | Delta | What happened |
|------------|----------|-------|---------------|
| Before pipeline | 3,960 | — | Index loaded (254K tx, 684K regions) |
| After scan | 17,513 | +13,553 | Buffer spilled 1 chunk (18.6M frags, 3.1 GB arrow), C++ accumulator |
| After finalize | 17,513 | 0 | Model finalization (negligible) |
| After calibration | 17,513 | 0 | calibrate_gdna (works on region-level data, not fragments) |
| After router_scan | 17,513 | 0 | C++ scoring converts buffer → CSR em_data |
| After buffer_release | 17,513 | 0 | Buffer released but RSS not returned to OS |
| **After locus_em** | **21,339** | **+3,826** | Locus EM allocations (one mega-locus) |
| After cleanup | 21,339 | 0 | em_data freed but RSS not returned |

**Peak RSS: 21.3 GB** for 14.7M fragments.

### Real BAM (small)

- **Flat at 3,943 MB** throughout — dominated by index loading
- The index itself consumes ~2–3 GB of baseline RSS

### Memory Hotspots

1. **FragmentBuffer spill** (3.1 GB arrow file for 18.6M fragments) — 1 chunk spilled to disk, but the in-memory accumulator holds all data before spill
2. **C++ BAM scanner accumulator** — holds all resolved fragments before finalization (~13.5 GB delta)
3. **Locus EM working set** — the mega-locus with 226K transcripts and 12.7M units requires large CSR arrays (~3.8 GB delta)
4. **Index baseline** — 254K transcripts × interval trees + cgranges occupies ~2–3 GB

---

## Hotspot Deep Dive

### Hotspot #1: `_vec_lgamma` / `_betaln_vec` in Calibration (52.7s self-time)

**Root cause**: `_vec_lgamma` wraps `math.lgamma` via `np.frompyfunc`, which:
1. Calls Python `math.lgamma` element-by-element (no SIMD, no C vectorization)
2. Returns an **object dtype array**, then casts to float64 (double allocation)
3. Called 2,808 times with arrays of ~100K–500K elements

**Call chain**: `calibrate_gdna` → `_e_step` (13 iters) → `_compute_strand_llr_betabinom` → `_betaln_vec` → `_vec_lgamma` × 3. Also called from `estimate_kappa_marginal` (golden-section search, ~35 evaluations × 13 iters × 2 calls each).

**Fix**: Replace with `scipy.special.gammaln` (C-compiled, truly vectorized):
```python
from scipy.special import gammaln
# Before: _lgamma_vec = np.frompyfunc(math.lgamma, 1, 1)
# After:  _vec_lgamma = gammaln  (or betaln directly)
```

**Expected speedup**: 50–100× for `_vec_lgamma` → **calibration drops from 150s to ~5–10s**.

### Hotspot #2: `_compute_fl_llr` Python Loop (12.6s via `builtins.min`)

**Root cause**: A Python for-loop iterates over all fragment-length observations (11.15M) per calibration iteration:
```python
for i in range(len(fl_region_ids)):
    rid = int(fl_region_ids[i])
    fl = int(fl_frag_lens[i])
    if rid < n_regions and fl > 0:
        idx = min(fl, max_size)  # ← 145M calls to builtins.min
        llr[rid] += log_ratio[idx]
```

**Fix**: Vectorize with NumPy scatter-add:
```python
valid = (fl_region_ids < n_regions) & (fl_frag_lens > 0)
fl_clamped = np.minimum(fl_frag_lens[valid], max_size)
np.add.at(llr, fl_region_ids[valid], log_ratio[fl_clamped])
```

**Expected speedup**: 100–1000× for this function → **saves ~12s per calibration run**.

### Hotspot #3: `estimate_kappa_marginal` Golden-Section Search (15s cumtime × 13 iters)

**Root cause**: Golden-section search over κ ∈ [0.01, 500] with tol=1e-4 requires ~50–100 function evaluations per call. Each evaluation computes full Beta-Binomial log-likelihood over all regions (2× `_betaln_vec`). Called 13 times (once per calibration EM iteration).

**Fix (combined with Hotspot #1)**:
1. Replace `_vec_lgamma` → `scipy.special.gammaln` (fixes the inner bottleneck)
2. Optionally switch to `scipy.optimize.minimize_scalar(method='bounded')` (Brent's method: 10–20 evaluations instead of 50–100)

**Expected speedup**: 5–10× from better optimizer + 50× from vectorized lgamma.

### Hotspot #4: `build_gdna_fl_model` Python Loop (1.3s, 27 calls)

**Root cause**: Loop over 1,000 histogram bins per call:
```python
for fl_val in range(1, max_fl + 1):
    if weighted_hist[fl_val] > 0:
        model.observe(fl_val, weight=float(weighted_hist[fl_val]))
```

**Fix**: Add batch initialization to `FragmentLengthModel`:
```python
model.set_counts(weighted_hist)  # direct array assignment
```

**Expected speedup**: 10–50× (eliminate Python loop).

### Hotspot #5: Mega-Locus in EM (111s for one locus)

**Root cause**: Connected components from multimapping reads create one mega-locus with 226,028 transcripts and 12,670,762 units. The C++ EM solver runs SQUAREM iterations on this massive problem:
- Per-iteration: CSR matrix-vector multiply over 12.7M entries
- ~1,000 max iterations (likely converges earlier, but problem is enormous)

**This is fundamentally an algorithmic challenge**, not a micro-optimization opportunity. The C++ solver is already well-optimized (SQUAREM, OpenMP, SIMD exp).

**Mitigation strategies**:
- **Locus decomposition**: Break mega-loci into sub-loci using graph partitioning (e.g., spectral clustering or METIS) to reduce per-locus size
- **Pruning**: Remove low-probability candidates (transcripts with no unique evidence) before EM
- **Early termination**: More aggressive convergence criteria for large loci
- **Hierarchical EM**: Two-level EM where the outer level handles gene-to-gene competition and inner level handles within-gene isoform competition

### Hotspot #6: `FragmentScorer.from_models` (1.5s fixed cost)

**Root cause**: Iterates over all 254K transcripts to precompute exon coordinate tuples:
```python
for t_idx in range(index.num_transcripts):
    exon_ivs = index.get_exon_intervals(t_idx)
    starts = tuple(int(x) for x in exon_ivs[:, 0])  # 1.9M generator calls
    ends = tuple(int(x) for x in exon_ivs[:, 1])     # 1.9M generator calls
```

**Fix**: Vectorize the conversion — use `.astype(int).tolist()` instead of `tuple(int(x) for x in ...)`, or better, pass raw NumPy arrays to C++ and avoid Python tuples entirely.

**Expected speedup**: 2–3× (1.5s → 0.5s).

---

## Performance Improvement Plan

### Phase 1: Quick Wins (estimated: major impact, low effort)

| # | Change | Target | Current | Expected | Priority |
|---|--------|--------|---------|----------|----------|
| 1 | Replace `_vec_lgamma` with `scipy.special.gammaln` | calibration.py | 52.7s self | <1s | **Critical** |
| 2 | Replace `_betaln_vec` with `scipy.special.betaln` | calibration.py | 42s cumtime | <1s | **Critical** |
| 3 | Vectorize `_compute_fl_llr` loop | calibration.py | 12.6s (min) | <0.1s | **Critical** |
| 4 | Vectorize `build_gdna_fl_model` loop | calibration.py | 1.3s | <0.1s | High |

**Combined Phase 1 impact: Calibration drops from ~150s to ~5–10s (simulated BAM), from ~9s to ~1–2s (real BAM).**

### Phase 2: Medium-Impact Optimizations

| # | Change | Target | Current | Expected | Priority |
|---|--------|--------|---------|----------|----------|
| 5 | Use `scipy.optimize.minimize_scalar` (Brent) instead of golden-section | calibration.py | ~35 evals/call | ~10–15 evals | Medium |
| 6 | Vectorize `FragmentScorer.from_models` exon tuple construction | scoring.py | 1.5s | 0.5s | Medium |
| 7 | Cache `fragment_classes` in `_FinalizedChunk` | buffer.py | Recomputed each time | Compute once | Low |

### Phase 3: Algorithmic — Mega-Locus Handling

| # | Change | Current | Expected | Priority |
|---|--------|---------|----------|----------|
| 8 | Graph-based locus decomposition (break mega-loci) | Single 226K-tx locus takes 111s | Reduce to sub-loci of O(1K) transcripts | **High** |
| 9 | Transcript pruning before EM (remove zero-evidence transcripts from locus) | All transcripts in connected component enter EM | Only transcripts with nonzero evidence | High |
| 10 | Adaptive convergence for large loci | Fixed convergence_delta for all loci | Relaxed convergence for mega-loci (1e-4 instead of 1e-6) | Medium |

### Phase 4: Infrastructure / Architecture

| # | Change | Notes |
|---|--------|-------|
| 11 | Index loading optimization (lazy-load intervals) | 8s index load is significant for small jobs |
| 12 | Reduce scanner accumulator memory | 13.5 GB delta during scan for 14.7M fragments |
| 13 | Stream-process mega-loci in tiled fashion | Avoid 3.8 GB locus EM working set |

---

## Memory Improvement Plan

### Current Memory Budget (Simulated BAM, 14.7M fragments)

| Component | Size (MB) | Notes |
|-----------|-----------|-------|
| Index (baseline) | ~3,960 | Interval trees, cgranges, DataFrames |
| C++ scanner accumulator | ~13,550 | 18.6M buffered fragments before spill |
| Calibration working set | ~0 | Operates on region summaries (684K regions) |
| ScoredFragments CSR | included above | Shares memory with scanner output |
| Locus EM working set | ~3,830 | Mega-locus: 226K tx × 12.7M units |
| **Peak RSS** | **~21,340** | |

### Improvement Opportunities

| # | Change | Savings | Priority |
|---|--------|---------|----------|
| M1 | **Streaming scanner output** — process fragments in chunks instead of accumulating all in memory | ~10 GB | **Critical** |
| M2 | **Spill earlier / more aggressively** — reduce max_memory_bytes to trigger spill before accumulator grows to 13 GB | ~8 GB | Critical |
| M3 | **Lazy index loading** — don't load interval arrays for transcripts not in the BAM | Up to 50% of index (~2 GB) | Medium |
| M4 | **Mega-locus decomposition** — partition into smaller loci to reduce EM working set | ~2–3 GB | High |
| M5 | **Arrow IPC for intermediate CSR** — spill ScoredFragments to disk with memory-mapped access | ~4 GB | Medium |
| M6 | **16-bit log-likelihoods** — quantize log-liks from float64 to float16/bfloat16 in CSR storage | 50% of CSR arrays | Low |

### Priority Timeline

1. **Immediate**: Phase 1 calibration fixes (items 1–4) — straightforward scipy substitutions
2. **Short-term**: Phase 3 mega-locus handling (items 8–9) — algorithmic improvement for worst-case loci
3. **Medium-term**: Memory improvements M1–M2 (streaming/spilling) — requires scanner architecture changes
4. **Long-term**: Phase 4 infrastructure changes — larger refactoring effort

---

## Appendix: cProfile Top Functions (Simulated BAM)

| Function | Self Time (s) | Calls | Cumtime (s) |
|----------|---------------|-------|-------------|
| `_thread.lock.acquire` | 306.5 | 9,764 | 339.9 |
| `calibration._vec_lgamma` | 52.7 | 2,808 | 59.3 |
| `builtins.min` | 12.6 | 145M | 12.6 |
| `ndarray.astype` | 6.9 | 2,914 | 6.9 |
| `profiler._snap_rss_current` | 2.5 | 2,446 | 65.3 |
| `calibration._marginal_loglik` | 1.6 | 455 | 23.5 |
| `threading.Condition.wait` | 1.3 | 2,441 | 441.3 |
| `numpy.asarray` | 1.3 | 19,014 | 1.3 |
| `calibration.build_gdna_fl_model` | 1.2 | 27 | 1.3 |
| `pyarrow.compute.take` | 1.0 | 18,972 | 1.0 |

Note: `_thread.lock.acquire` (306s) and `threading.Condition.wait` (441s cumtime) are cProfile artifacts from the background memory-sampling thread and C++ OpenMP threads — **not actual CPU bottlenecks**.

## Appendix: Locus Size Distribution (Simulated BAM)

- **9,486 loci** total
- **Max transcripts per locus**: 226,028 (one mega-locus with ~89% of all transcripts)
- **Max units per locus**: 12,670,762
- This mega-locus drives the 111s EM time — most other loci complete in <1ms

The mega-locus is a known consequence of multimapping reads linking distant genes into connected components. This is the primary target for algorithmic optimization.

## Appendix: Platform

- macOS-26.3-arm64 (Apple Silicon)
- Python 3.12.13
- Clang 19.1.7 (C++ compiler)
- conda-forge + bioconda dependencies
