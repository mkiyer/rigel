# EM Profiling Instrumentation

**Added:** March 28, 2026

This document describes the per-locus profiling instrumentation added to the C++ EM solver and how to completely remove it to restore the original code.

## What Was Added

### Purpose
Per-locus timing, equivalence class statistics, iteration counts, and sub-phase breakdowns for the batch locus EM. Enabled via `emit_locus_stats=True` flag — zero overhead when disabled (no timing calls, no stats collection).

### Features
When `emit_locus_stats=True`:
- **Per-locus wall time** for the total `process_locus` call
- **Sub-phase timing** (microseconds): `extract`, `bias`, `build_ec`, `warm_start`, `squarem`, `assign`
- **EC statistics**: `n_equiv_classes`, `ec_total_elements` (sum of n×k), `max_ec_width` (max k), `max_ec_depth` (max n)
- **SQUAREM iteration count** per locus
- **Thread count** used for E-step
- **Mega-locus flag** (whether the locus was Phase 1)
- **Digamma calls per E-step** (VBEM-specific: = n_components; MAP: = 0)
- Results sorted by descending wall time for easy identification of bottleneck loci

### Files Modified

#### 1. `src/rigel/native/em_solver.cpp`

**Changes:**
- Added `#include <chrono>`, `#include <nanobind/stl/string.h>`, `#include <nanobind/stl/vector.h>`
- Added `LocusProfile` struct and `hrclock` alias after the constants section
- Modified `EMResult` struct: added `int squarem_iterations` field
- Modified `run_squarem()`: tracks `completed_iterations` counter in both VBEM and MAP branches
- Modified `batch_locus_em()`:
  - Added `bool emit_locus_stats` parameter (default `false`)
  - Changed return type to include `nb::list` as 4th tuple element
  - Added `std::vector<LocusProfile> locus_profiles` allocation (only when `emit_locus_stats=true`)
  - Modified `process_locus` lambda: added `bool mega` parameter, added `std::chrono` timing around each sub-phase, collects EC stats
  - Updated Phase 1 and Phase 2 calls to pass `mega` flag (`true`/`false`)
  - Added stats serialization to `nb::list` of `nb::dict` after GIL re-acquisition
- Updated nanobind binding for `batch_locus_em`: added `emit_locus_stats` kwarg with default `false`

#### 2. `src/rigel/estimator.py`

**Changes:**
- `AbundanceEstimator.__init__()`: Added `self.locus_stats: list[dict] | None = None`
- `run_batch_locus_em()`:
  - Added `emit_locus_stats: bool = False` parameter
  - Changed unpacking from 3 to 4 return values (`total_gdna_em, locus_mrna, locus_gdna, locus_stats_raw`)
  - Passes `emit_locus_stats` to C++ `_batch_locus_em()`
  - Stores `self.locus_stats = list(locus_stats_raw)` when enabled

#### 3. `src/rigel/pipeline.py`

**Changes:**
- `_run_locus_em()`: Added `emit_locus_stats: bool = False` keyword parameter, passed through to `estimator.run_batch_locus_em()`
- `quant_from_buffer()`: Added `emit_locus_stats: bool = False` parameter, passed through to `_run_locus_em()`
- `run_pipeline()`: Passes `config.emit_locus_stats` to `quant_from_buffer()`

#### 4. `src/rigel/config.py`

**Changes:**
- `PipelineConfig`: Added `emit_locus_stats: bool = False` field

#### 5. `scripts/profiler.py`

**Changes:**
- `StageTimings`: Added `locus_stats: list | None = field(default=None, repr=False)`
- `profile_stages()`: Passes `emit_locus_stats=True` to `estimator.run_batch_locus_em()`, captures `estimator.locus_stats` into `timings.locus_stats`
- `run_profiles()`: Saves locus stats to separate `locus_stats_{config}.json` file, excludes `locus_stats` from main `profile_summary.json`

---

## How to Undo All Changes

To restore the code to its pre-instrumentation state, apply the following reversals:

### Step 1: `src/rigel/native/em_solver.cpp`

1. **Remove added includes**: Remove `#include <chrono>`, `#include <nanobind/stl/string.h>`, `#include <nanobind/stl/vector.h>` (keep the original includes)

2. **Remove `LocusProfile` struct and `hrclock`**: Delete the entire block from `using hrclock = std::chrono::steady_clock;` through the end of the `LocusProfile` struct definition (after `static constexpr int ASSIGN_SAMPLE = 2;`)

3. **Restore `EMResult`**: Remove the `int squarem_iterations = 0;` field

4. **Restore `run_squarem()`**: 
   - Remove `int completed_iterations = 0;` declaration
   - In VBEM branch: replace the `if (delta < convergence_delta) { completed_iterations = iter + 1; break; } completed_iterations = iter + 1;` block with just `if (delta < convergence_delta) break;`
   - In MAP branch: same replacement
   - Change return to `return { std::move(theta), std::move(alpha_out), std::move(em_totals) };`

5. **Restore `batch_locus_em()` return type**: Change back to 3-tuple (remove `nb::list`)

6. **Remove `emit_locus_stats` parameter** from function signature

7. **Remove `locus_profiles` vector** allocation

8. **Restore `process_locus` lambda**:
   - Remove `bool mega` parameter
   - Remove `auto locus_t0 = hrclock::now();` at start
   - Remove all `auto t1..t7 = hrclock::now();` timing points
   - Remove the entire `if (emit_locus_stats)` block at the end
   
9. **Restore Phase 1/Phase 2 calls**: Remove the `true`/`false` `mega` arguments
   - Phase 1: `process_locus(li, actual_threads, sub, local_map, pool.get())`
   - Phase 2: `process_locus(li, 1, sub, local_map)`

10. **Remove stats serialization**: Delete the `nb::list stats_list;` block and restore the 3-tuple return

11. **Restore nanobind binding**: Remove `emit_locus_stats` kwarg, restore docstring to mention 3-tuple return

### Step 2: `src/rigel/estimator.py`

1. Remove `self.locus_stats: list[dict] | None = None` from `__init__`
2. Remove `emit_locus_stats: bool = False` parameter from `run_batch_locus_em`
3. Change unpacking back: `total_gdna_em, locus_mrna, locus_gdna = _batch_locus_em(...)`
4. Remove `emit_locus_stats` from the `_batch_locus_em()` call arguments
5. Remove the `if emit_locus_stats:` / `else:` block after the call

### Step 3: `src/rigel/pipeline.py`

1. Remove `emit_locus_stats: bool = False` from `_run_locus_em` signature
2. Remove `emit_locus_stats=emit_locus_stats` from the `estimator.run_batch_locus_em()` call
3. Remove `emit_locus_stats: bool = False` from `quant_from_buffer` signature
4. Remove `emit_locus_stats=emit_locus_stats` from the `_run_locus_em()` call
5. Remove `emit_locus_stats=config.emit_locus_stats` from the `quant_from_buffer()` call in `run_pipeline`

### Step 4: `src/rigel/config.py`

1. Remove `emit_locus_stats: bool = False` from `PipelineConfig`

### Step 5: `scripts/profiler.py`

1. Remove `locus_stats: list | None = field(default=None, repr=False)` from `StageTimings`
2. Remove `timings.locus_stats = getattr(estimator, "locus_stats", None)` line
3. Remove `emit_locus_stats=True` from the `estimator.run_batch_locus_em()` call
4. Restore `"stages": asdict(r.stages),` (remove the dict comprehension filter)
5. Remove the locus stats JSON saving block

### Step 6: Recompile

```bash
conda activate rigel && pip install --no-build-isolation -e .
```

### Step 7: Delete this file

```bash
rm docs/performance/PROFILING_INSTRUMENTATION.md
```
