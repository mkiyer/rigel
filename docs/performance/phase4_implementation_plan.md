# Phase 4 — Quantification-Stage Performance & Cleanup

**Date**: 2026-03-24
**Status**: Planned
**Prerequisite**: Phase 3c (streaming chunk architecture) complete

---

## Motivation

Phases 3a–3c focused on the C++ BAM scan pipeline (streaming chunks, zero-copy
finalization, merge elimination). Those phases reduced `scan_and_buffer` from
the dominant bottleneck to ~40% of wall time on real data. The remaining ~60%
is now distributed across Python-dominated quantification stages that have
clear, measurable inefficiencies.

### Baseline Profile (Real BAM, 1.6M fragments, 457K transcripts)

| Stage | Time | % | Bottleneck |
|-------|------|---|------------|
| scan_and_buffer | 9.12s | 40.2% | C++ streaming scan |
| locus_em | 3.72s | 16.4% | C++ EM solver |
| **calibration** | **2.51s** | **11.1%** | Scipy optimizer loop |
| **fragment_scorer** | **2.44s** | **10.7%** | Python 457K-iteration loop |
| **eb_gdna_priors** | **2.39s** | **10.5%** | _merged_intervals × 12K loci |
| **build_loci** | **2.21s** | **9.7%** | _merged_intervals × 12K loci |
| fragment_router_scan | 0.29s | 1.3% | C++ scoring |
| **TOTAL** | **22.70s** | | |

**Bolded stages** are the Phase 4 targets: **9.55s (42.1%)** of wall time in
pure Python code with known optimization paths. The C++ stages (scan, EM,
router) are out of scope for this phase.

### cProfile Evidence

| Function | Self Time | Calls | Root Cause |
|----------|-----------|-------|------------|
| `_merged_intervals` | 5.32s cumul. | 69,022 | Python dict+sort+merge per locus |
| `pyarrow.compute.take` | 2.21s | 24,432 | Pandas string column indexing |
| `calibration._marginal_loglik` | 1.37s | 363 | Scipy line search inner loop |
| `gc.collect` | 0.69s | 1 | Explicit collection after scoring |
| `numpy.cumsum` | 0.47s | 457,516 | Per-transcript exon data build |

---

## Design Principles

1. **Correctness first.** Every change must pass the full test suite (1018 tests + 21 golden outputs) with zero regressions.
2. **Elegant and clean.** No hacks, no workarounds. If a function is rewritten, the new version must be clearer than the old one. Dead code is removed, not commented out.
3. **Vectorize, don't loop.** Replace Python for-loops over per-transcript or per-locus data with bulk numpy/C++ operations.
4. **Compute once, use everywhere.** Data that's needed by multiple stages should be computed once and passed through, not recomputed.
5. **Measure every change.** Profile before and after each sub-phase to confirm improvements.

---

## Phase 4a — Vectorized Interval Merging & Locus Footprints

**Target**: `build_loci` (2.21s) + `eb_gdna_priors` (2.39s) = **4.60s → ~0.3s**
**Estimated improvement**: ~4.3s (19% of wall time)

### Problem

`_merged_intervals()` is called **twice per locus** — once in `build_loci()`
to compute `gdna_span`, and again in `compute_gdna_locus_gammas()` to find
overlapping calibration regions. With 12K loci and 457K transcripts, this is
69K+ calls to a Python function that builds a `defaultdict`, sorts intervals,
and merges them in a manual loop.

### Solution: Vectorized bulk computation

Replace the per-locus Python generator with a **single vectorized pass** that
computes all locus footprints and merged intervals in bulk.

#### 4a.1 — Vectorized `gdna_span` computation

Replace the `_merged_intervals` → `sum(e - s)` pattern in `build_loci` with
a vectorized numpy computation:

```python
def _compute_locus_gdna_spans(
    comp_t_offsets: np.ndarray,
    comp_t_indices: np.ndarray,
    t_ref_ids: np.ndarray,    # int32 ref_id (not string)
    t_starts: np.ndarray,
    t_ends: np.ndarray,
) -> np.ndarray:
    """Compute merged genomic footprint (bp) for each connected component.

    Uses a vectorized sort-and-scan approach:
    1. Expand all transcript intervals with their locus assignment
    2. Sort by (locus_id, ref_id, start)
    3. Detect merge boundaries with a cummax scan
    4. Sum non-overlapping spans per locus
    """
```

**Key insight**: We already have `comp_t_offsets` and `comp_t_indices` in CSR
form from the C++ union-find. We can assign each transcript a locus_id via
`np.repeat`, then do a single global sort + cummax merge.

**Requires**: `t_ref_ids` as `int32` (not pandas string column). The index
already has `ref_to_id` — we just need to expose a pre-built int32 array.

#### 4a.2 — Cache merged intervals for `compute_gdna_locus_gammas`

The `compute_gdna_locus_gammas` function needs merged `(ref, start, end)`
intervals per locus (not just the total span) for cgranges overlap queries.
Two options:

**Option A — Pre-compute merged intervals during `build_loci`**, store on
the `Locus` dataclass as a list of `(ref_id, start, end)` tuples:

```python
@dataclass(slots=True)
class Locus:
    locus_id: int
    transcript_indices: np.ndarray
    unit_indices: np.ndarray
    gdna_span: int
    merged_intervals: list[tuple[str, int, int]]  # NEW: cached
```

Then `compute_gdna_locus_gammas` simply iterates `locus.merged_intervals`
instead of recomputing them.

**Option B — Vectorized batch cgranges query.** Build flat arrays of all
merged intervals across all loci, batch-query cgranges once, then aggregate
results per locus. More complex but eliminates the per-locus Python loop
entirely.

**Decision**: **Option A.** It's simpler, eliminates the recomputation, and
the remaining per-locus cgranges loop in `compute_gdna_locus_gammas` is
already fast (0.1s for the overlap queries themselves — the 2.4s was dominated
by `_merged_intervals` recomputation).

#### 4a.3 — Cleanup

- **Delete** the `_merged_intervals()` generator function entirely — replaced
  by the vectorized `_compute_locus_gdna_spans()` and the cached
  `merged_intervals` on `Locus`.
- **Remove** the `defaultdict` import if no longer used.

### Files Modified

- [src/rigel/locus.py](../../src/rigel/locus.py) — Rewrite `build_loci()` with vectorized span computation, add `merged_intervals` to `Locus` construction, simplify `compute_gdna_locus_gammas()` to use cached intervals, delete `_merged_intervals()`.
- [src/rigel/scored_fragments.py](../../src/rigel/scored_fragments.py) — Add `merged_intervals` field to `Locus` dataclass.
- [src/rigel/index.py](../../src/rigel/index.py) — Expose `t_ref_id_arr` (int32 array mapping transcript → integer ref_id) as a cached property if not already available.

### Verification

```bash
pytest tests/ -v                    # 1018 tests
pytest tests/test_golden_output.py  # 21 golden outputs
python scripts/profiler.py \
  --bam <real_bam> --index <index> \
  --outdir results/profile_phase4a \
  --stages --threads 8
```

**Expected**: `build_loci` drops from 2.21s to ~0.15s, `eb_gdna_priors`
drops from 2.39s to ~0.15s.

---

## Phase 4b — Bulk Exon Data Construction

**Target**: `fragment_scorer` (2.44s) → **~0.05s**
**Estimated improvement**: ~2.4s (10.6% of wall time)

### Problem

`FragmentScorer.from_models()` contains a Python for-loop over all 457,513
transcripts:

```python
t_exon_data: dict = {}
for t_idx in range(index.num_transcripts):
    exon_ivs = index.get_exon_intervals(t_idx)
    if exon_ivs is not None and len(exon_ivs) > 0:
        starts = tuple(exon_ivs[:, 0].tolist())
        ends = tuple(exon_ivs[:, 1].tolist())
        lengths = exon_ivs[:, 1] - exon_ivs[:, 0]
        cumsum_before = tuple(np.concatenate(([0], np.cumsum(lengths[:-1]))).tolist())
        t_exon_data[t_idx] = (starts, ends, cumsum_before)
```

Each iteration does: dict lookup → numpy slice → `.tolist()` → `tuple()` →
`np.cumsum` → `np.concatenate` → `.tolist()` → `tuple()`. The result is a
Python dict of tuples, which the C++ `NativeFragmentScorer` constructor then
converts *back* to flat CSR arrays via another two-pass loop.

This is a round-trip: **flat arrays → dict of tuples → flat arrays**.

### Solution: Build CSR arrays directly on the Python side

The index already stores exon intervals in a per-transcript dict
(`_t_exon_intervals`). We can build the CSR arrays directly from the
underlying interval data without the intermediate dict-of-tuples:

#### 4b.1 — Add `build_exon_csr()` to `TranscriptIndex`

```python
def build_exon_csr(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build CSR arrays for per-transcript exon positions.

    Returns (offsets, starts, ends, cumsum_before) as contiguous int32
    arrays ready for direct C++ consumption.
    """
```

This method iterates the existing `_t_exon_intervals` dict **once** and
builds four flat numpy arrays: `offsets[n_transcripts + 1]`, `starts[n_exons]`,
`ends[n_exons]`, `cumsum_before[n_exons]`.

#### 4b.2 — Update `NativeFragmentScorer` to accept CSR arrays

Change the nanobind constructor to accept four numpy arrays instead of a
Python dict:

```cpp
NativeFragmentScorer(
    // ... existing params ...
    nb::ndarray<const int32_t, nb::ndim<1>> exon_offsets,
    nb::ndarray<const int32_t, nb::ndim<1>> exon_starts,
    nb::ndarray<const int32_t, nb::ndim<1>> exon_ends,
    nb::ndarray<const int32_t, nb::ndim<1>> exon_cumsum,
    // ... remaining params ...
)
```

The constructor simply copies the four arrays into its member vectors —
**no dict iteration, no tuple unpacking**.

#### 4b.3 — Update `FragmentScorer.from_models()`

Replace the 457K-iteration loop with a single `index.build_exon_csr()` call
and pass the arrays through to `NativeFragmentScorer`.

The Python-side `_t_exon_data` dict on `FragmentScorer` can be **removed**
if no other code path uses it. Check for remaining references first.

#### 4b.4 — Cleanup

- **Remove** the `_t_exon_data` field from `FragmentScorer` if unused
  outside of `NativeFragmentScorer` construction.
- **Remove** the `t_exon_data=` parameter from `NativeFragmentScorer`'s
  Python binding.
- **Simplify** the `NativeFragmentScorer` constructor in scoring.cpp — the
  two-pass dict-unpacking code (~30 lines) becomes a simple array copy
  (~5 lines).

### Files Modified

- [src/rigel/index.py](../../src/rigel/index.py) — Add `build_exon_csr()` method.
- [src/rigel/scoring.py](../../src/rigel/scoring.py) — Replace 457K loop with `build_exon_csr()`, remove `_t_exon_data` field.
- [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) — Simplify `NativeFragmentScorer` constructor to accept CSR arrays.

### Verification

Same as 4a — full test suite + golden outputs + profiler comparison.

**Expected**: `fragment_scorer` drops from 2.44s to ~0.05s.

---

## Phase 4c — Dead Code Removal & Hygiene

**Target**: Code cleanliness and minor performance wins.
**Estimated improvement**: ~0.7s (gc.collect) + code hygiene

### 4c.1 — Remove dead `merge_accumulator_into` function

[bam_scanner.cpp lines 344–397](../../src/rigel/native/bam_scanner.cpp) contains
a 53-line `merge_accumulator_into()` static function that is **never called**.
It was part of an older single-accumulator merge design superseded by Phase 3b's
per-worker accumulators and Phase 3c's streaming chunks.

**Action**: Delete the function.

### 4c.2 — Audit `gc.collect()` calls

[pipeline.py line 398](../../src/rigel/pipeline.py):

```python
del builder, ctx
buffer.release()
gc.collect()
```

This `gc.collect()` takes **0.69s** (3% of wall time) on the real BAM. The
`del` statements already drop the reference counts to zero for these objects —
CPython's reference-counting GC will free them immediately without needing a
full collection cycle.

**Action**: Remove the `gc.collect()` call after the scoring phase. Keep the
`del` statements and `buffer.release()` — those are sufficient for deterministic
cleanup. Keep the `gc.collect()` after EM (pipeline.py line 603) as it precedes
the final output stage where memory pressure matters most.

Also remove the `import gc` if no longer needed.

### 4c.3 — Remove `_merged_intervals` import remnants

After Phase 4a, the `defaultdict` import in `locus.py` may be unused. Clean
up any orphaned imports.

### 4c.4 — Update `Locus` docstring

The `gdna_span` field is undocumented in the `Locus` dataclass. Phase 4a adds
`merged_intervals`. Update the docstring to cover both fields.

### Files Modified

- [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) — Delete `merge_accumulator_into()`.
- [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — Remove `gc.collect()` after scoring, clean imports.
- [src/rigel/locus.py](../../src/rigel/locus.py) — Clean imports after 4a.
- [src/rigel/scored_fragments.py](../../src/rigel/scored_fragments.py) — Update `Locus` docstring.

### Verification

Full test suite + golden outputs. No profiler run needed (changes are
mechanical).

---

## Phase 4d — Profile & Validate

**Target**: Confirm all improvements, document results, update repo memory.

### 4d.1 — Profile real BAM

```bash
python scripts/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/real_bams/mctp_LBX0345_SI_43019_start.str.rmdup.collate.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v5_phase3b/rigel_index \
  --outdir results/profile_phase4_final \
  --stages --cprofile --threads 8
```

### 4d.2 — Profile simulated BAM (regression check)

```bash
python scripts/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v5_phase3b/gdna_none_ss_0.95_nrna_none/align_minimap2/reads_namesort.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine_v5_phase3b/rigel_index \
  --outdir results/profile_phase4_sim \
  --stages --threads 8
```

### 4d.3 — Run pristine benchmark (accuracy regression check)

```bash
python scripts/benchmark.py \
  --config /Users/mkiyer/Downloads/rigel_runs/benchmark_pristine.yaml
```

Compare MAE, Spearman, and pool counts against Phase 3c baseline.

### 4d.4 — Document results

Publish performance comparison to `docs/performance/`.

---

## Summary

| Phase | Target | Est. Savings | Risk |
|-------|--------|-------------|------|
| **4a** | Vectorized interval merge + cache | ~4.3s (19%) | Low — pure Python refactor |
| **4b** | Bulk exon CSR construction | ~2.4s (11%) | Low — eliminates round-trip |
| **4c** | Dead code + gc.collect removal | ~0.7s (3%) | None — mechanical cleanup |
| **4d** | Profile + validate | 0s | None — measurement only |
| **Total** | | **~7.4s (33%)** | |

**Projected wall time**: 22.7s → ~15.3s on the real BAM.

### Execution Order

```
4a → 4b → 4c → 4d
```

Phases 4a and 4b are independent and could be parallelized, but sequential
execution is safer for incremental verification. Phase 4c is a cleanup pass
that depends on 4a completing. Phase 4d is always last.

### What's NOT in Phase 4

- **C++ scan optimization** — requires py-spy profiling to decompose; separate
  investigation.
- **C++ EM solver optimization** — 3.7s on real data, already well-optimized
  with SQUAREM + OpenMP. Diminishing returns.
- **Calibration optimizer tuning** — 2.5s, dominated by 363 scipy minimize
  calls. Requires algorithmic change (analytical gradient, different optimizer)
  rather than mechanical speedup.
- **Multithreading for Python stages** — The GIL prevents true parallelism in
  Python. The C++ stages already use OpenMP. No easy wins here without
  multiprocessing (which adds complexity).
