# O3: Index Memory Optimization — Implementation Plan

**Date**: 2026-04-15  
**Scope**: Reduce `TranscriptIndex` memory footprint through categorical encoding, dead structure removal, and cleanup  
**Baseline**: Index loads at 2836 MB RSS before any quantification begins  

## Summary

The `TranscriptIndex` carries three classes of memory waste:

1. **String columns in DataFrames** — `ref`, `t_id`, `g_id`, `g_name`, `g_type` stored as full Python `str` objects (56+ bytes overhead each), duplicated across 457K rows even though most values repeat.
2. **Dead Python structures** — `_iv_t_set` (frozensets), `sj_map` (dict of frozensets), and `_t_exon_intervals` (dict of ndarrays) are duplicated into the C++ `FragmentResolver` during `load()`, then retained on the Python instance despite never being accessed in production code.
3. **Integer column bloat** — numeric columns in `t_df` default to int64 (8 bytes) even when int32 or int16 would suffice.

Together these account for an estimated 300–500 MB of recoverable RSS.

---

## Investigation Findings

### String columns in `t_df` (457K rows)

| Column | Unique values | Per-row overhead | Total waste | Used in hot path? |
|--------|-------------|-----------------|-------------|-------------------|
| `ref` | ~25 | 56+ bytes × 457K | ~25 MB | Once (locus.py, already re-encodes to int codes) |
| `t_id` | ~457K (1:1) | 56+ bytes × 457K | ~50 MB | Output only (estimator.py, annotate.py) |
| `g_id` | ~63K | 56+ bytes × 457K | ~40 MB | Output only (estimator.py) |
| `g_name` | ~63K | 56+ bytes × 457K | ~40 MB | Output only (estimator.py, annotate.py) |
| `g_type` | ~50 | 56+ bytes × 457K | ~25 MB | Never (only used deriving g_df) |

The `g_df` gene table (63K rows) inherits `ref`, `g_id`, `g_name`, `g_type` strings via groupby `first()` aggregation, adding another ~30 MB.

**Total string waste: ~210 MB**

### Dead Python structures retained after `load()`

| Structure | Type | Approx size | Used post-load? |
|-----------|------|-------------|-----------------|
| `_iv_t_set` | `list[frozenset[int]]` (~200K entries) | 50–100 MB | Only `query()`, which is test-only |
| `sj_map` | `dict[tuple(str,int,int,int), frozenset[int]]` | 30–50 MB | **Never** in production code |
| `_t_exon_intervals` | `dict[int, ndarray]` (~250K entries) | 30–50 MB | Once via `build_exon_csr()` in scoring.py, then dead |

All three are projected into C++ `FragmentResolver` during `load()` via CSR arrays. The Python originals persist only because they serve test-only methods (`query()`, `get_exon_intervals()`) and the one-time `build_exon_csr()` call.

**Total dead structure waste: ~110–200 MB**

### Integer column bloat in `t_df`

| Column | Current dtype | Data range | Target dtype | Savings per row |
|--------|-------------|------------|--------------|----------------|
| `start`, `end` | int64 | 0–3B | int64 | 0 (needed) |
| `strand` | int64 | 0–3 | int8 | 7 bytes |
| `t_index`, `g_index` | int64 | 0–457K | int32 | 4 bytes |
| `n_exons` | int64 | 1–~200 | int16 | 6 bytes |
| `nrna_t_index` | int64 | -1–457K | int32 | 4 bytes |
| `nrna_n_contributors` | int64 | 0–~50 | int16 | 6 bytes |
| `length` | int64 | 0–~2M | int32 | 4 bytes |

**Potential savings: ~31 bytes/row × 457K = ~14 MB** (minor, but "free" cleanup)

---

## Implementation Plan

### Phase 1: Categorical encoding of string columns

Convert the five string columns in `t_df` to `pd.Categorical` immediately after loading from feather. Categoricals store one small integer code per row plus a compact dictionary of unique values.

#### Step 1a: Add categorical conversion in `load()`

In `index.py` `load()`, immediately after `pd.read_feather()`:

```python
self.t_df = pd.read_feather(os.path.join(index_dir, TRANSCRIPTS_FEATHER))
for col in ("ref", "t_id", "g_id", "g_name", "g_type"):
    if col in self.t_df.columns:
        self.t_df[col] = self.t_df[col].astype("category")
```

This is safe because:
- `.iloc[i]` on a categorical returns the original string value
- Equality comparisons (`==`, `!=`, `isin`) work identically
- `groupby` on categoricals is faster (operates on integer codes)
- Feather/Arrow natively supports dictionary encoding, so round-trips work

#### Step 1b: Verify downstream consumers are categorical-compatible

All downstream access patterns identified in the investigation use these operations:

| Operation | Categorical-safe? |
|-----------|--------------------|
| `t_df["ref"].to_numpy(dtype=object)` (locus.py L70) | Yes — converts to strings |
| `t_df["t_id"].values` (estimator.py) | Yes — returns CategoricalArray, usable in DataFrame construction |
| `t_df["g_id"].values` / `t_df["g_name"].values` | Yes — same as above |
| `g_df = t_df.groupby("g_index").agg(ref=("ref", "first"), ...)` | Yes — `first()` on categorical returns the string value |
| Direct string comparison or indexing | Yes — categorical `.iloc[i]` returns string |

**One attention point**: `t_df["ref"].values` returns a `pd.CategoricalArray`, not a numpy object array. Code that calls `np.unique()` on it (like locus.py) already converts via `.to_numpy(dtype=object)`, so this is safe.

#### Step 1c: Apply to `g_df` as well

Since `g_df` is derived from `t_df` via `_build_gene_table()`, and groupby `first()` on categoricals returns plain strings, `g_df` will have string columns. Convert them to categorical in `_build_gene_table()`:

```python
for col in ("ref", "g_id", "g_name", "g_type"):
    if col in g_df.columns:
        g_df[col] = g_df[col].astype("category")
```

**Expected savings: ~210 MB from t_df + ~30 MB from g_df = ~240 MB**

---

### Phase 2: Free dead Python structures after C++ projection

#### Step 2a: Free `sj_map` after resolver construction

`sj_map` is never accessed after `load()` in production code. After the `ctx.build_sj_map()` call completes, set it to `None`:

```python
ctx.build_sj_map(sj_refs_l, sj_starts_l, sj_ends_l, sj_strands_l, sj_t_flat, sj_t_offsets)
self.sj_map = None  # C++ resolver owns the data now
```

Tests that access `self.sj_map` will need to be updated to either:
- Use a dedicated test fixture that retains the dict, or
- Skip the assertion (the C++ resolver behavior is tested via integration tests anyway)

#### Step 2b: Free `_iv_t_set` after resolver construction

`_iv_t_set` is only used by `query()`, which is a test-only Python method. After `ctx.build_overlap_index()` completes:

```python
ctx.build_overlap_index(...)
self._iv_t_set = None  # C++ resolver owns the data now
```

The `query()` method should be updated to raise a clear error if called after structures are freed:

```python
def query(self, exon):
    if self._iv_t_set is None:
        raise RuntimeError("query() unavailable: interval sets freed after load()")
    ...
```

Tests using `query()` should load the index into a test-specific variant that retains the structures, or test via the C++ resolver path instead.

#### Step 2c: Free `_t_exon_intervals` after `build_exon_csr()` 

`_t_exon_intervals` is needed once by `build_exon_csr()` (called from `scoring.py`). After that call, the dict is dead.

Add cleanup at the end of `build_exon_csr()`:

```python
def build_exon_csr(self):
    # ... existing CSR construction ...
    self._t_exon_intervals = None  # no longer needed
    return offsets, starts, ends
```

**Expected savings: ~110–200 MB**

---

### Phase 3: Downcast numeric columns in `t_df`

After reading the feather file, downcast integer columns to their minimum safe width:

```python
_DOWNCAST = {
    "strand": np.int8,
    "t_index": np.int32,
    "g_index": np.int32,
    "n_exons": np.int16,
    "nrna_t_index": np.int32,
    "nrna_n_contributors": np.int16,
    "length": np.int32,
}
for col, dtype in _DOWNCAST.items():
    if col in self.t_df.columns:
        self.t_df[col] = self.t_df[col].astype(dtype)
```

This is safe because:
- These columns are used as array indices (int32 is sufficient for 457K transcripts)
- Numpy implicit widening handles arithmetic correctly
- The numpy arrays extracted from them (`t_to_g_arr`, `t_to_strand_arr`) should be created with matching dtypes

Update the numpy array extraction to use consistent types:

```python
self.t_to_g_arr = self.t_df["g_index"].values  # now int32
self.t_to_strand_arr = self.t_df["strand"].values  # now int8
```

**Attention**: `t_to_g_arr` is passed to C++ code. Ensure the C++ bindings accept int32. The `set_metadata()` call and scoring C++ both currently cast via `nb::cast<i32_1d>`, so int32 is already the expected type — this change actually *removes* a hidden widening bug where int64 was passed and silently accepted.

**Expected savings: ~14 MB** (minor but consistent with the project's memory discipline)

---

### Phase 4: Test updates

#### 4a: Tests accessing `sj_map`

In `tests/test_index_integrity.py`, tests that read `index.sj_map` need adjustment. Options:

- **Best**: Add a `_load_for_testing()` classmethod that retains the structures for test use. This is a thin wrapper around `load()` that skips the freeing steps.
- **Alternative**: Use a `retain_python_structures=True` parameter on `load()`.

#### 4b: Tests accessing `_iv_t_set` / `query()`

Same approach — the test fixture loads with retained structures.

#### 4c: Tests accessing `_t_exon_intervals` / `get_exon_intervals()`

These tests should call `get_exon_intervals()` *before* `build_exon_csr()` frees the dict. Or the test fixture retains the structures via the testing classmethod.

#### 4d: Golden output tests

Categorical encoding and integer downcasting do not change the *values* in any output — only the in-memory representation. Golden output tests should pass without regeneration. Verify by running the full test suite.

---

### Phase 5: Cleanup opportunities

#### 5a: Remove `int64` → `int32` conversion in `locus.py`

Currently, `locus.py` L70 does:
```python
t_ref_strs = index.t_df["ref"].to_numpy(dtype=object, na_value="")
_ref_names, _ref_codes = np.unique(t_ref_strs, return_inverse=True)
```

After Phase 1, `t_df["ref"]` is categorical. This can use the category codes directly:

```python
ref_cat = index.t_df["ref"].cat
_ref_codes = ref_cat.codes.values  # int8 or int16
_ref_names = ref_cat.categories.values  # string array
```

This eliminates the `np.unique()` O(N log N) sort on 457K string objects.

#### 5b: Remove redundant `_build_gene_table()` groupby columns

`_build_gene_table()` computes `g_type` for every gene, but `g_type` is never used anywhere after the gene table is built. Consider dropping it:

```python
# In _build_gene_table: remove g_type from agg()
```

This simplifies the schema and removes the last `g_type` consumer.

#### 5c: Consolidate `t_to_g_arr` and `t_to_strand_arr` types

Currently these are extracted from `t_df` columns and stored as separate numpy arrays. After downcasting, they naturally become int32 and int8 respectively. Ensure all C++ consumers accept these types (they already do — C++ casts via `nb::cast<i32_1d>` and `nb::cast<i8_1d>`).

#### 5d: `intervals.feather` and `sj.feather` string `ref` columns

These feather files store `ref` as plain strings. During `load()`, the `ref` column is immediately converted to integer codes via `np.unique()`. The string column is never used further. Future index versions could store `ref` as an integer code with a separate lookup table, but this is a backwards-compatibility breaking change and not worth doing now (the loading code already handles this efficiently).

---

## Execution Order

| Step | Description | Files changed | Risk |
|------|-------------|---------------|------|
| 1a | Categorical conversion of `t_df` string columns | `index.py` | Low — no functional change |
| 1c | Categorical conversion of `g_df` string columns | `index.py` | Low |
| 3 | Downcast `t_df` integer columns | `index.py` | Low — verify C++ type expectations |
| 5a | Optimize `locus.py` ref encoding via category codes | `locus.py` | Low |
| 5c | Verify numpy array type consistency | `index.py`, `pipeline.py` | Low |
| 2a | Free `sj_map` post-projection | `index.py` | Low — test-only impact |
| 2b | Free `_iv_t_set` post-projection | `index.py` | Low — test-only impact |
| 2c | Free `_t_exon_intervals` after `build_exon_csr()` | `index.py` | Low — single use post-load |
| 4 | Update tests for freed structures | `test_index_integrity.py` | Medium — test refactoring |
| 5b | Drop `g_type` from gene table | `index.py` | Low — unused column |

## Verification

1. **Unit tests**: `pytest tests/ -v` — all 1028 tests must pass
2. **Golden output tests**: `pytest tests/test_golden_output.py -v` — no regeneration needed
3. **Profiling**: Re-run profiler to measure RSS reduction at the `before` snapshot (index load) and at peak

## Expected Impact

| Component | Current | After O3 | Savings |
|-----------|---------|----------|---------|
| String columns (t_df + g_df) | ~240 MB | ~5 MB (categorical codes) | ~235 MB |
| Dead structures (sj_map + _iv_t_set + _t_exon) | ~140–200 MB | 0 | ~140–200 MB |
| Integer bloat | ~14 MB | 0 | ~14 MB |
| **Total** | | | **~390–450 MB** |

Some of this will manifest as reduced `before` RSS (index load), some as reduced fragmentation that lowers peak RSS. Allocator fragmentation means actual RSS reduction may be 50–70% of theoretical savings.

**Conservative estimate: 200–300 MB peak RSS reduction.**
