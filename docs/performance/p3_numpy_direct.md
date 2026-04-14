# P3: Unified Zero-Copy C++ → numpy Transfer

## Problem Statement

Three independent C++ → Python data transfer sites in the rigel native layer
use the same capsule-backed `nb::ndarray` pattern, each with its own copy of
the boilerplate (~60 lines of lambda definitions per site). A fourth transfer
path (`BamScanner::build_result()`) still uses nanobind's default `list_caster`,
which converts `std::vector<T>` → Python `list` → `np.asarray()` — a double
conversion costing ~6.4s (7.1% of samples in py-spy profiling on VCaP).

Additionally, several C++ observation vectors use `int32_t` for data that fits
in `int8_t` (strand enums {0,1,2,3}, splice type enums {0,1,2}), wasting 4×
memory and forcing dtype narrowing in Python.

### Current state: four transfer sites, three patterns

| Site | File | Pattern | Lines of Boilerplate |
|------|------|---------|---------------------|
| `FragmentAccumulator::finalize_zero_copy()` | `resolve_context.h:330` | 5 lambdas (`mk_u8`…`mk_u32`) + 17 heap-moves | ~80 |
| `StreamingScorer::finish()` | `scoring.cpp:1112` | 5 lambdas (`mk_i64`…`mk_i8`) + 19 heap-ptr transfers | ~65 |
| EM solver return paths | `em_solver.cpp:1181` | Inline capsule (raw `new[]`) | ~12 |
| `BamScanner::build_result()` | `bam_scanner.cpp:1380` | **None** — bare `std::move()` → `list_caster` | 0 (bottleneck) |

**Goal**: One template. Four sites. Zero boilerplate.

### Vectors affected in `build_result()`

| Dict Key | Current C++ Type | Target C++ Type | Size (VCaP) | Python Consumer |
|----------|-----------------|----------------|-------------|-----------------|
| `strand_observations/exonic_spliced_obs` | `vector<int32_t>` | `vector<int8_t>` | ~5M | `_replay_strand_observations` |
| `strand_observations/exonic_spliced_truth` | `vector<int32_t>` | `vector<int8_t>` | ~5M | `_replay_strand_observations` |
| `strand_observations/exonic_obs` | `vector<int32_t>` | `vector<int8_t>` | ~14M | `_replay_strand_observations` |
| `strand_observations/exonic_truth` | `vector<int32_t>` | `vector<int8_t>` | ~14M | `_replay_strand_observations` |
| `strand_observations/intergenic_obs` | `vector<int32_t>` | `vector<int8_t>` | ~3.5M | `_replay_strand_observations` |
| `strand_observations/intergenic_truth` | `vector<int32_t>` | `vector<int8_t>` | ~3.5M | `_replay_strand_observations` |
| `frag_length_observations/lengths` | `vector<int32_t>` | `vector<int32_t>` | ~15M | `_replay_fraglen_observations` |
| `frag_length_observations/splice_types` | `vector<int32_t>` | `vector<int8_t>` | ~15M | `_replay_fraglen_observations` |
| `frag_length_observations/intergenic_lengths` | `vector<int32_t>` | `vector<int32_t>` | ~3.5M | `_replay_fraglen_observations` |
| `region_evidence/counts` | `vector<double>` | `vector<double>` | ~2.7M | `np.asarray(...).reshape(n,4)` |
| `region_evidence/fl_region_ids` | `vector<int32_t>` | `vector<int32_t>` | ~13M | `pd.DataFrame` |
| `region_evidence/fl_frag_lens` | `vector<int32_t>` | `vector<int32_t>` | ~13M | `pd.DataFrame` |
| `region_evidence/fl_frag_strands` | `vector<int32_t>` | `vector<int8_t>` | ~13M | `pd.DataFrame` |

## Design

### Core abstraction: `vec_to_ndarray<T>`

A single C++17 function template in its own header replaces all per-type
lambda definitions across the entire codebase:

```cpp
// ndarray_util.h
template <typename T>
nb::object vec_to_ndarray(std::vector<T>&& v) {
    auto* heap = new std::vector<T>(std::move(v));
    size_t n = heap->size();
    nb::capsule del(heap, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        heap->data(), {n}, del).cast();
}
```

This is zero-copy: the numpy array's data pointer points directly at the
vector's internal buffer, and the capsule destructor frees it when the
Python object is garbage-collected. The captureless lambda works because
each template instantiation generates its own function pointer with the
correct `static_cast<std::vector<T>*>`.

For `scoring.cpp` which stores vectors as raw pointers (`std::vector<T>*`),
a pointer-taking overload transfers ownership from a raw heap pointer:

```cpp
template <typename T>
nb::object vec_to_ndarray(std::vector<T>* v) {
    size_t n = v->size();
    nb::capsule del(v, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        v->data(), {n}, del).cast();
}
```

### Why a dedicated header?

`resolve_context.h` is a 1400-line header concerned exclusively with
fragment resolution (resolver, accumulator, region accumulator). Adding
generic numpy utilities there forces any future consumer to include a
heavy, project-specific header for a 10-line utility. A clean
`ndarray_util.h` depends only on `<vector>`, `<cstddef>`,
`nanobind/nanobind.h`, and `nanobind/ndarray.h` — zero project-specific
dependencies.

### Why narrow types at the source?

Strand observation values are {0, 1, 2, 3} (`STRAND_NONE`, `STRAND_POS`,
`STRAND_NEG`, `STRAND_AMBIGUOUS`). Splice type values are {0, 1, 2}
(`SPLICE_UNSPLICED`, `SPLICE_SPLICED_UNANNOT`, `SPLICE_SPLICED_ANNOT`).
Region fragment strands are {1, 2} (`STRAND_POS`, `STRAND_NEG`).

Storing these in `int32_t`:
- Wastes 4× memory (e.g. strand obs: ~180 MB → ~45 MB after narrowing)
- Forces a dtype narrowing copy in Python (`np.asarray(obs, dtype=np.int8)`)
- Misrepresents the data contract — `int32` implies values up to 2 billion

With `int8_t` at the source, the Python consumer receives an `int8` ndarray
that it can use directly — no `np.asarray()`, no `.astype()`, no copy. The
C++ type documents the truth.

### Why not `nb::ndarray` return type annotations?

Nanobind can auto-convert `std::vector<T>` → `nb::ndarray` via type casters,
but only when the function's return type is declared as `nb::ndarray`. Since
`build_result()` returns `nb::dict` containing heterogeneous values, we must
manually wrap each vector. The template makes this a single-expression
operation per vector.

## Implementation Plan

### Step 1: Create `ndarray_util.h`

**New file**: `src/rigel/native/ndarray_util.h`

```cpp
/**
 * ndarray_util.h — Zero-copy std::vector<T> → numpy ndarray transfer.
 *
 * Moves a vector to heap-allocated storage and returns a capsule-backed
 * 1-D numpy array.  The capsule destructor owns the memory; the vector
 * is consumed (left moved-from) after the call.
 *
 * Two overloads:
 *   vec_to_ndarray(std::vector<T>&& v)  — move from rvalue
 *   vec_to_ndarray(std::vector<T>* v)   — take ownership of heap pointer
 */

#pragma once

#include <cstddef>
#include <vector>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

namespace nb = nanobind;

namespace rigel {

/// Move a std::vector<T> to the heap and return a capsule-backed 1-D ndarray.
/// The source vector is consumed (moved-from) after this call.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>&& v) {
    auto* heap = new std::vector<T>(std::move(v));
    size_t n = heap->size();
    nb::capsule del(heap, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        heap->data(), {n}, del).cast();
}

/// Take ownership of a heap-allocated std::vector<T>* and return a
/// capsule-backed 1-D ndarray.  The pointer must not be used after this call.
template <typename T>
nb::object vec_to_ndarray(std::vector<T>* v) {
    size_t n = v->size();
    nb::capsule del(v, [](void* p) noexcept {
        delete static_cast<std::vector<T>*>(p);
    });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(
        v->data(), {n}, del).cast();
}

}  // namespace rigel
```

**Total**: ~45 lines including header comment and include guard. Zero
project-specific dependencies.

### Step 2: Narrow C++ observation types to match data range

**File**: `src/rigel/native/bam_scanner.cpp`

#### 2a. `StrandObservations` (line 181)

**Before**:
```cpp
struct StrandObservations {
    std::vector<int32_t> exonic_spliced_obs;
    std::vector<int32_t> exonic_spliced_truth;
    std::vector<int32_t> exonic_obs;
    std::vector<int32_t> exonic_truth;
    std::vector<int32_t> intergenic_obs;
    std::vector<int32_t> intergenic_truth;
};
```

**After**:
```cpp
struct StrandObservations {
    std::vector<int8_t> exonic_spliced_obs;
    std::vector<int8_t> exonic_spliced_truth;
    std::vector<int8_t> exonic_obs;
    std::vector<int8_t> exonic_truth;
    std::vector<int8_t> intergenic_obs;
    std::vector<int8_t> intergenic_truth;
};
```

Values stored are `STRAND_POS` (1), `STRAND_NEG` (2), and occasionally
`STRAND_NONE` (0) or `STRAND_AMBIGUOUS` (3), all from `constants.h`
(`int32_t` constants). The `push_back` sites need `static_cast<int8_t>()`:

- Line 1242: `strand_obs.intergenic_obs.push_back(static_cast<int8_t>(ig_strand));`
- Line 1243: `strand_obs.intergenic_truth.push_back(static_cast<int8_t>(STRAND_POS));`
- Line 1292: `strand_obs.exonic_spliced_obs.push_back(static_cast<int8_t>(result.exon_strand));`
- Line 1293: `strand_obs.exonic_spliced_truth.push_back(static_cast<int8_t>(result.sj_strand));`
- Line 1312: `strand_obs.exonic_obs.push_back(static_cast<int8_t>(result.exon_strand));`
- Line 1313: `strand_obs.exonic_truth.push_back(static_cast<int8_t>(t_strand));`

The `merge_strand_obs` function (line 330) uses `insert()` with iterators
and requires no changes — iterator types adapt automatically to `int8_t`.

#### 2b. `FragLenObservations::splice_types` (line 196)

**Before**:
```cpp
struct FragLenObservations {
    std::vector<int32_t> lengths;
    std::vector<int32_t> splice_types;
    std::vector<int32_t> intergenic_lengths;
};
```

**After**:
```cpp
struct FragLenObservations {
    std::vector<int32_t> lengths;          // genuine int32 range (1–1000+)
    std::vector<int8_t>  splice_types;     // SpliceType enum: {0, 1, 2}
    std::vector<int32_t> intergenic_lengths;
};
```

Push site (line 1322): `fraglen_obs.splice_types.push_back(static_cast<int8_t>(result.splice_type));`

The `merge_fraglen_obs` function (line 349) adapts automatically.

#### 2c. `RegionAccumulator::fl_frag_strands` (resolve_context.h, line 1270)

**Before**: `std::vector<int32_t> fl_frag_strands;`
**After**: `std::vector<int8_t> fl_frag_strands;`

Values are always `STRAND_POS` (1) or `STRAND_NEG` (2). The `push_back`
site in `RegionAccumulator::accumulate()` needs a cast. The merge function
(`merge_from`) uses iterator-based concatenation and adapts automatically.

### Step 3: Refactor `finalize_zero_copy()` in `resolve_context.h`

Add `#include "ndarray_util.h"` to `resolve_context.h`.

Replace the 80-line lambda + heap-allocation block with direct
`vec_to_ndarray()` calls. The `ambig_strand` local vector is stack-allocated,
filled, then moved:

**Before** (~80 lines):
```cpp
auto* v_splice_type = new std::vector<uint8_t>(std::move(splice_type_));
// ... 15 more heap allocations ...
auto mk_u8 = [](std::vector<uint8_t>* v) -> nb::object { ... };
// ... 4 more lambda definitions ...
result["splice_type"] = mk_u8(v_splice_type);
// ... 16 more dict assignments ...
```

**After** (~25 lines):
```cpp
nb::dict result;
result["splice_type"]       = vec_to_ndarray(std::move(splice_type_));
result["exon_strand"]       = vec_to_ndarray(std::move(exon_strand_));
result["sj_strand"]         = vec_to_ndarray(std::move(sj_strand_));
result["num_hits"]          = vec_to_ndarray(std::move(num_hits_));
result["merge_criteria"]    = vec_to_ndarray(std::move(merge_criteria_));
result["chimera_type"]      = vec_to_ndarray(std::move(chimera_type_));
result["t_indices"]         = vec_to_ndarray(std::move(t_indices_));
result["t_offsets"]         = vec_to_ndarray(std::move(t_offsets_));
result["frag_lengths"]      = vec_to_ndarray(std::move(frag_lengths_));
result["exon_bp"]           = vec_to_ndarray(std::move(exon_bp_));
result["intron_bp"]         = vec_to_ndarray(std::move(intron_bp_));
result["ambig_strand"]      = vec_to_ndarray(std::move(ambig_strand));
result["frag_id"]           = vec_to_ndarray(std::move(frag_id_));
result["read_length"]       = vec_to_ndarray(std::move(read_length_));
result["genomic_footprint"] = vec_to_ndarray(std::move(genomic_footprint_));
result["genomic_start"]     = vec_to_ndarray(std::move(genomic_start_));
result["nm"]                = vec_to_ndarray(std::move(nm_));
result["size"]              = nb::cast(n);
size_ = 0;
```

Note: `ambig_strand` is now a stack-allocated `std::vector<uint8_t>` (no
longer `new`-allocated). `vec_to_ndarray()` moves it to the heap internally.

### Step 4: Apply `vec_to_ndarray()` to `build_result()` in `bam_scanner.cpp`

Add `#include "ndarray_util.h"` to `bam_scanner.cpp`.

**Strand observations** (now `int8_t` vectors from Step 2):
```cpp
strand_dict["exonic_spliced_obs"]   = vec_to_ndarray(std::move(strand_obs_.exonic_spliced_obs));
strand_dict["exonic_spliced_truth"] = vec_to_ndarray(std::move(strand_obs_.exonic_spliced_truth));
strand_dict["exonic_obs"]           = vec_to_ndarray(std::move(strand_obs_.exonic_obs));
strand_dict["exonic_truth"]         = vec_to_ndarray(std::move(strand_obs_.exonic_truth));
strand_dict["intergenic_obs"]       = vec_to_ndarray(std::move(strand_obs_.intergenic_obs));
strand_dict["intergenic_truth"]     = vec_to_ndarray(std::move(strand_obs_.intergenic_truth));
```

**Fragment length observations** (lengths remain `int32_t`, splice_types now `int8_t`):
```cpp
fraglen_dict["lengths"]            = vec_to_ndarray(std::move(fraglen_obs_.lengths));
fraglen_dict["splice_types"]       = vec_to_ndarray(std::move(fraglen_obs_.splice_types));
fraglen_dict["intergenic_lengths"] = vec_to_ndarray(std::move(fraglen_obs_.intergenic_lengths));
```

**Region evidence** (counts remain `double`, region_ids and frag_lens remain `int32_t`, frag_strands now `int8_t`):
```cpp
region_dict["counts"]         = vec_to_ndarray(std::move(region_acc_.counts));
region_dict["fl_region_ids"]  = vec_to_ndarray(std::move(region_acc_.fl_region_ids));
region_dict["fl_frag_lens"]   = vec_to_ndarray(std::move(region_acc_.fl_frag_lens));
region_dict["fl_frag_strands"] = vec_to_ndarray(std::move(region_acc_.fl_frag_strands));
```

### Step 5: Refactor `StreamingScorer::finish()` in `scoring.cpp`

Add `#include "ndarray_util.h"` to `scoring.cpp`.

Replace the 65-line lambda block with direct `vec_to_ndarray(pointer)` calls.
The pointer-taking overload consumes raw `std::vector<T>*` pointers:

**Before** (~65 lines):
```cpp
auto mk_i64 = [](std::vector<int64_t>* v) -> nb::object { ... };
// ... 4 more lambdas ...
auto result = nb::make_tuple(
    mk_i64(v_offsets_),
    mk_i32(v_ti_),
    // ... 17 more ...
);
```

**After** (~25 lines):
```cpp
auto result = nb::make_tuple(
    vec_to_ndarray(v_offsets_),
    vec_to_ndarray(v_ti_),
    vec_to_ndarray(v_ll_),
    vec_to_ndarray(v_ct_),
    vec_to_ndarray(v_cw_),
    vec_to_ndarray(v_ts_),
    vec_to_ndarray(v_te_),
    vec_to_ndarray(v_lt_),
    vec_to_ndarray(v_lct_),
    vec_to_ndarray(v_isp_),
    vec_to_ndarray(v_gll_),
    vec_to_ndarray(v_gfp_),
    vec_to_ndarray(v_fid_),
    vec_to_ndarray(v_fc_),
    vec_to_ndarray(v_st_),
    vec_to_ndarray(v_dti_),
    vec_to_ndarray(v_dfid_),
    vec_to_ndarray(v_chim_fid_),
    vec_to_ndarray(v_chim_stype_),
    stat_det_, stat_em_u_, stat_em_as_, stat_em_ao_,
    stat_gated_, stat_chim_, stat_mm_
);
```

Null out pointers after transfer (capsule now owns memory):
```cpp
v_offsets_ = nullptr;  v_ti_ = nullptr;  v_ll_ = nullptr;
// ... etc ...
```

### Step 6: Update Python consumers

**File**: `src/rigel/pipeline.py`

#### 6a. `_replay_strand_observations()` (line 103)

With `int8_t` vectors from C++, the data arrives as `dtype=int8` ndarray.
`StrandModel.observe_batch()` internally calls `_np.asarray()` which is a
no-op on an existing ndarray. Remove the now-unnecessary dtype conversion:

**Before**:
```python
model.observe_batch(
    np.asarray(obs, dtype=np.int8),
    np.asarray(truth, dtype=np.int8),
)
```

**After**:
```python
model.observe_batch(obs, truth)
```

The `obs` and `truth` are already `int8` ndarray — the
`_np.asarray()` inside `observe_batch` is a no-op.

#### 6b. `_replay_fraglen_observations()` (line 121)

`lengths` arrives as `int32` ndarray, `splice_types` as `int8` ndarray.
Inside `FragmentLengthModels.observe_batch()`, both are run through
`np.asarray(..., dtype=np.intp)` which handles the upcast. No change
needed.

**However**: the `observe_intergenic_batch()` call also uses
`np.asarray(intergenic_lengths, dtype=np.intp)`. This already works
correctly on an `int32` ndarray (widening copy). No change needed.

#### 6c. Region evidence (line ~308)

`np.asarray()` on a capsule-backed ndarray with matching dtype is a no-op.
`fl_frag_strands` changes from `int32` to `int8`, but the Python consumer
constructs a `pd.DataFrame` which handles mixed dtypes. No change needed.

### Step 7: Verify the `_FinalizedChunk.from_raw()` path is unaffected

`_FinalizedChunk.from_raw()` in `buffer.py` receives streaming chunk data
from `finalize_zero_copy()`. It uses `_arr()`, which checks
`isinstance(val, np.ndarray)` and returns a view or copies for dtype
narrowing. Since `finalize_zero_copy()` already returns `np.ndarray`
objects (and will continue to do so via `vec_to_ndarray`), this path is
unchanged in behavior. The only difference is that the ndarray construction
now happens through the template instead of through inline lambdas — same
memory layout, same dtypes.

## Files Modified

| File | Change | Lines Removed | Lines Added |
|------|--------|:------------:|:-----------:|
| `src/rigel/native/ndarray_util.h` | **New file**: template `vec_to_ndarray<T>` | — | ~45 |
| `src/rigel/native/resolve_context.h` | Refactor `finalize_zero_copy()`, narrow `fl_frag_strands` | ~85 | ~25 |
| `src/rigel/native/bam_scanner.cpp` | Narrow `StrandObservations`/`FragLenObservations`, use `vec_to_ndarray` in `build_result()` | ~20 | ~25 |
| `src/rigel/native/scoring.cpp` | Replace lambdas in `finish()` with `vec_to_ndarray` | ~60 | ~25 |
| `src/rigel/pipeline.py` | Remove `np.asarray(obs, dtype=np.int8)` wrappers | ~6 | ~3 |
| `tests/test_ndarray_util.py` | **New file**: unit tests for `vec_to_ndarray` | — | ~80 |

**Net**: ~170 lines removed, ~200 lines added (including new header + tests).

## Testing Strategy

### Tier 1: Unit test the template itself (`test_ndarray_util.py`)

New test file that imports the C++ modules and directly validates the
`vec_to_ndarray` contract:

```python
"""Tests for the vec_to_ndarray C++ → numpy zero-copy transfer.

Exercises the ndarray_util.h template through the public C++ module
interfaces that use it (finalize_zero_copy, build_result, scoring finish).
"""
```

**Test cases**:
1. **dtype preservation**: Verify that strand observation arrays arrive as
   `int8`, splice_types as `int8`, lengths as `int32`, counts as `float64`,
   chunk arrays preserve their original dtypes (uint8, uint16, int32, int64,
   uint32).

2. **zero-copy verification**: Confirm that capsule-backed arrays are not
   copied by checking `arr.flags.owndata is False` (numpy flag indicating
   the array does not own its data — it's backed by another object).

3. **empty vector handling**: Verify that zero-length vectors produce
   valid 0-element numpy arrays (not errors or None).

4. **value integrity**: Round-trip a known set of strand/fraglen values
   through the BAM scanner's build_result() path and verify the numpy
   arrays contain exactly the expected values.

5. **memory safety**: Ensure that the numpy array remains valid after
   the source dict is deleted (capsule prevents premature deallocation).

### Tier 2: Existing integration tests (no changes needed)

These tests exercise the full data flow and will catch any regressions:

- `test_pipeline_smoke.py` — End-to-end: oracle scenario → scan → EM → output DataFrames
- `test_golden_output.py` — Bit-exact regression against golden feather files
- `test_calibration_integration.py` — Exercises region evidence → calibration
- `test_calibration.py` — Region evidence fractional counting
- `test_pipeline_routing.py` — Scoring → routing → locus construction

### Tier 3: Existing component tests (no changes needed)

- `test_strand_model.py` — `StrandModel.observe()` / posterior computation
  (note: does not currently test `observe_batch` — Tier 1 fills this gap)
- `test_frag_length_model.py` — `FragmentLengthModel.observe()` / statistics
  (note: does not currently test `observe_batch` — Tier 1 fills this gap)
- `test_buffer.py` — `_FinalizedChunk.from_raw()` / `BufferedFragment`

### Tier 4: Golden output regeneration

After confirming all tests pass, regenerate and inspect golden outputs:

```bash
pytest tests/test_golden_output.py -v --update-golden
git diff tests/golden/
```

Expected: zero diff (same data, different transfer mechanism).

## Build and Test

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/ -v
```

## Expected Impact

### Performance
- **Eliminates double conversion**: ~120M elements no longer traverse
  C++ vector → nanobind list_caster → Python list → `np.asarray()`
- **Eliminates dtype narrowing**: strand int32→int8, splice_type int32→int8
  no longer require a copy in Python
- **Estimated savings**: 3–5s wall time on VCaP BAM (19.4M fragments)

### Memory
- **Strand observations**: ~180 MB → ~45 MB (4× reduction from int32→int8)
- **Splice types**: ~60 MB → ~15 MB (4× reduction from int32→int8)
- **Region frag_strands**: ~52 MB → ~13 MB (4× reduction from int32→int8)
- **Zero-copy ndarray**: no additional memory — numpy array backs directly
  onto the C++ vector's heap allocation

### Code quality
- **~170 lines of boilerplate removed** across three C++ files
- **One template** replaces four independent implementations of the same pattern
- **C++ types match data contracts**: `int8_t` for enums, `int32_t` for counts/lengths
- **Zero-dependency utility header**: reusable by any future C++ → numpy transfer

### Correctness
- Semantically identical: same data, same values, same Python-side dtypes
  (after `np.asarray()` promotion where needed)
- Golden output tests confirm bit-for-bit identical results
- The template is exercised ~16× per BAM scan (chunks) + 1× (build_result)
  + 1× (scoring finish) — no new code paths, just unified existing ones

## Risk Assessment

**Low risk**. The capsule-backed ndarray pattern is already battle-tested
across three sites in the codebase. This change:

1. Unifies them under one template — same machine code, less source duplication
2. Applies the pattern to one remaining site (`build_result`) — straightforward
3. Narrows integer types to match data range — the `push_back` sites cast from
   `int32_t` constants to `int8_t`, which are guaranteed to fit

The only subtle risk is the `static_cast<int8_t>()` narrowing at push_back
sites. The strand constants are `{0, 1, 2, 3}` and splice types are `{0, 1, 2}`,
all well within `int8_t` range `[-128, 127]`. This is verified by the
existing integration tests that exercise every strand/splice code path.
