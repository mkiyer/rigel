# Gene-to-Transcript Refactor Plan

## Problem Statement

The codebase has three categories of incorrect gene-level logic:

1. **Strand lookups use gene indirection** — code does `t_to_g[t_idx]` →
   `g_to_strand[g_idx]` when `t_strand_arr[t_idx]` gives the answer directly.
   The strand model itself is trained from gene-level "truth" (`gene_strand`)
   when it should use transcript strand.

2. **Fragment classification is gene-based** — fragments are classified as
   `FRAG_UNIQUE` / `FRAG_AMBIG_SAME_STRAND` / `FRAG_AMBIG_OPP_STRAND` using `n_genes`.
   This is logically wrong: classification should be strand-based and
   gene-agnostic.

3. **Estimation accumulators are gene-indexed** — `gene_sense_all`,
   `gene_antisense_all`, `gdna_gene_summary`, `gene_spans`, `intronic_spans`
   are residual gene-level arrays from an obsolete architecture.

## Design Principles

- **All estimation at transcript level.** Gene aggregation happens only at
  output time in `get_gene_counts_df()`.
- **No gene indirection for strand.** Use `t_strand_arr` (=`t_to_strand_arr`)
  everywhere. The `t_to_g` → `g_to_strand` roundtrip is never needed for
  strand.
- **Fragment classification is strand-based.** Classes:
  - `FRAG_UNIQUE` (1 transcript, NH=1) — unchanged.
  - `FRAG_AMBIG_SAME_STRAND` (>1 transcript, all same strand, NH=1) —
    replaces `FRAG_AMBIG_SAME_STRAND`.
  - `FRAG_AMBIG_OPP_STRAND` (>1 transcript, mixed strands, NH=1) —
    replaces `FRAG_AMBIG_OPP_STRAND`.
  - `FRAG_MULTIMAPPER` (NH > 1) — unchanged.
  - `FRAG_CHIMERIC` — unchanged.

---

## Phase 1 — Strand lookups: remove gene indirection ✅ DONE

**Goal:** Every place that looks up strand via `t_to_g → g_to_strand` switches
to using `t_strand_arr` (which equals `index.t_to_strand_arr`, already
available on `HulkIndex` and `ScoringContext`).

### Python changes

| File | What | Change |
|------|------|--------|
| `estimator.py` `assign_unique()` | `g_idx = t_to_g[t_idx]` → `g_to_strand[g_idx]` | Use `index.t_to_strand_arr[t_idx]` directly. Remove `g_idx` intermediate. |
| `estimator.py` `is_antisense()` | Takes `gene_strand: int` | Rename param to `ref_strand` (strand of the reference transcript). No logic change — callers change. |
| `scan.py` lines 700-714 | `g_idx = t_to_g[t_inds[0]]` → `gene_strand = g_to_strand[g_idx]` | Use `t_strand = t_strand_arr[t_inds[0]]` directly. |
| `scan.py` lines 700-714 | `counter.gene_sense_all[g_idx]` / `gene_antisense_all[g_idx]` accumulation | Keep for Phase 4 (will convert to transcript-level). For now, use g_idx derivation only for these accumulators. |
| `scoring.py` `ScoringContext` | Has `g_strand_arr` and `t_to_g` fields | Keep for now (still needed by scan accumulators until Phase 4). |

### C++ changes

| File | What | Change |
|------|------|--------|
| `resolve_context.h` `ResolveContext` | `set_gene_strands()` stores `g_to_strand_arr_` | Add `set_transcript_strands()` storing `t_strand_arr_`. |
| `resolve_context.h` `_resolve_core()` | Computes `n_genes` from `t_to_g_arr_` | Keep for now — Phase 2 replaces with strand classification. |
| `resolve_context.h` `ResolvedResult` | `get_is_unique_gene()` uses `n_genes == 1` | Keep for now — Phase 2 changes. |
| `resolve_context.h` `get_is_strand_qualified()` | Uses `n_genes == 1` | Keep for now — Phase 2 changes to use `is_unique_transcript()` or `is_same_strand()`. |
| `bam_scanner.cpp` exonic model training | `gene_idx = t_to_g[t_idx]` → `gene_strand = g_to_strand[gene_idx]` | Use `t_strand_arr_[t_idx]` directly. |
| `pipeline.py` | `resolve_ctx.set_gene_strands(...)` | Add `resolve_ctx.set_transcript_strands(index.t_to_strand_arr.tolist())`. |

### Tests

- All 689+ existing tests must pass.
- Semantic equivalence: transcript strand == gene strand for all transcripts.

---

## Phase 2 — Fragment classification: strand-based ✅ DONE

**Goal:** Replace `n_genes`-based classification with strand-based
classification. The `n_genes` field is replaced by `mixed_strand` (bool/uint8).

### New constants (buffer.py)

```python
FRAG_UNIQUE: int = 0             # 1 transcript, NH=1
FRAG_AMBIG_SAME_STRAND: int = 1  # >1 transcript, all same strand, NH=1
FRAG_AMBIG_OPP_STRAND: int = 2   # >1 transcript, both strands, NH=1
FRAG_MULTIMAPPER: int = 3        # NH > 1
FRAG_CHIMERIC: int = 4           # chimeric
```

### Data changes

| Location | Old | New |
|----------|-----|-----|
| `CoreResult` (constants.h) | `n_genes` | `mixed_strand` (bool: true if transcripts span both strands) |
| `ResolvedResult` (resolve_context.h) | `n_genes` | `mixed_strand` |
| `_FinalizedChunk` (buffer.py) | `n_genes: np.ndarray (uint8)` | `mixed_strand: np.ndarray (uint8, 0/1)` |
| `BufferedFragment` (buffer.py) | `n_genes: int` | `mixed_strand: int` |
| `ResolvedFragment` (resolution.py) | `n_genes: int` | `mixed_strand: int` |
| `NativeAccumulator.finalize()` | Computes `n_genes` per fragment from `t_to_g` | Computes `mixed_strand` per fragment from `t_strand_arr` |

### Classification logic (buffer.py `fragment_classes`)

```python
n_transcripts = np.diff(self.t_offsets).astype(np.intp)
classes = np.full(self.size, FRAG_UNIQUE, dtype=np.uint8)
# >1 transcript, all same strand → AMBIG_SAME_STRAND
classes[
    (~self.mixed_strand.astype(bool)) & (n_transcripts > 1) & (self.num_hits == 1)
] = FRAG_AMBIG_SAME_STRAND
# >1 transcript, mixed strands → AMBIG_OPP_STRAND
classes[
    (self.mixed_strand.astype(bool)) & (n_transcripts > 1) & (self.num_hits == 1)
] = FRAG_AMBIG_OPP_STRAND
# NH > 1 overrides
classes[self.num_hits > 1] = FRAG_MULTIMAPPER
# Chimeric overrides all
classes[self.chimera_type > 0] = FRAG_CHIMERIC
```

### C++ changes

- `_resolve_core()`: replace `n_genes` computation with `mixed_strand`
  computation using `t_strand_arr_` (set via `set_transcript_strands()`).
- `NativeAccumulator::finalize()`: replace `n_genes` with `mixed_strand`
  computed by checking whether all transcripts share the same strand
  (`t_strand_arr`).
- `ResolvedResult`: `get_is_unique_gene()` → `is_same_strand()` (returns
  `!mixed_strand`). `get_is_strand_qualified()` check also updated.

### Property renames

| Old | New |
|-----|-----|
| `is_unique_gene` | `is_same_strand` |
| `FRAG_AMBIG_SAME_STRAND` | `FRAG_AMBIG_SAME_STRAND` |
| `FRAG_AMBIG_OPP_STRAND` | `FRAG_AMBIG_OPP_STRAND` |

### Consumer updates

All code that checks `fc == FRAG_AMBIG_SAME_STRAND` or `fc == FRAG_AMBIG_OPP_STRAND`
or `bf.n_genes` must be updated:

- `scan.py`: `fc == FRAG_UNIQUE or fc == FRAG_AMBIG_SAME_STRAND` →
  `fc == FRAG_UNIQUE or fc == FRAG_AMBIG_SAME_STRAND`
- `scan.py`: strand accumulation condition (unique OR isoform_ambig = all
  single-gene unique fragments) — review whether AMBIG_OPP_STRAND should be
  included.
- `bam_scanner.cpp`: stat counters `n_unique_gene`, `n_multi_gene` → rename
  to `n_same_strand`, `n_mixed_strand`.
- `buffer.py` `FragmentBuffer.__init__`: `t_to_g_arr` parameter — still needed
  until Phase 4 for gDNA prior accumulators, but `finalize()` switches to
  `t_strand_arr`.
- All tests referencing `n_genes`, `FRAG_AMBIG_SAME_STRAND`, `FRAG_AMBIG_OPP_STRAND`.

---

## Phase 3 — Estimation accumulators: gene-level → transcript-level

**Goal:** Convert the remaining gene-indexed arrays to transcript-level.

### `gene_sense_all` / `gene_antisense_all` → transcript-level

- **New arrays:** `transcript_unspliced_sense[n_transcripts]`,
  `transcript_unspliced_antisense[n_transcripts]`.
- **Accumulation (scan.py):** For each unique/isoform-ambig unspliced
  fragment, distribute `1/n_cand` to each compatible transcript's sense or
  antisense accumulator (matching how `transcript_intronic_sense/antisense` is
  already done).
- **gDNA prior (locus.py `compute_eb_gdna_priors`):** Sum transcript-level
  accumulators per locus (using `locus.transcript_indices`), per chromosome
  (using transcript ref), and globally. Remove gene indirection entirely.
- **Output (`get_gene_counts_df`):** Aggregate `n_sense_all` /
  `n_antisense_all` from transcript-level using `np.add.at(g_arr, t_to_g,
  t_arr)`.

### `gdna_gene_summary` → derive from `gdna_locus_counts`

- **Current:** `gdna_gene_summary[g_idx] += per_gene` in pipeline.py.
- **New:** Remove `gdna_gene_summary`. In `get_gene_counts_df()`, compute
  gene-level gDNA by aggregating `gdna_locus_counts` (already per-transcript):
  `np.add.at(n_gdna, t_to_g, self.gdna_locus_counts.sum(axis=1))`.
- **Remove** the per-gene distribution code in pipeline.py L494-497.

### Dead code removal

| What | Where |
|------|-------|
| `gene_spans` / `intronic_spans` in `TranscriptGeometry` | config.py, pipeline.py |
| `self._gene_spans` / `self._intronic_spans` in `AbundanceEstimator` | estimator.py |
| `index.intron_span_per_gene()` | index.py |
| `index.exon_lengths_per_gene()` | index.py |
| `AbundanceEstimator.num_genes` | estimator.py (remove from __init__ params) |
| `Locus.gene_indices` | estimator.py, locus.py |

### Locus.gene_indices removal

`Locus.gene_indices` is used in two places:
1. `compute_eb_gdna_priors()` — will use `locus.transcript_indices` instead.
2. `pipeline.py` gDNA gene attribution — removed (see above).

After Phase 3, `gene_indices` can be removed from the `Locus` dataclass.

---

## Phase 4 — Cleanup

**Goal:** Remove all vestigial gene-level plumbing that is no longer needed.

- Remove `g_strand_arr` from `ScoringContext` (if no remaining consumers).
- Remove `set_gene_strands()` from C++ `ResolveContext`.
- Remove `g_to_strand_arr_` from C++ `ResolveContext`.
- Remove `t_to_g_arr` from `NativeAccumulator.finalize()` signature
  (replace with `t_strand_arr`).
- Consider removing `t_to_g` from `ScoringContext` if no remaining consumers
  in scoring path.
- Remove `t_to_g_arr` parameter from `FragmentBuffer.__init__()` if no longer
  needed.
- Update all test files.
- Run full test suite, benchmarks.

---

## Execution Order

| Phase | Scope | Risk | Effort |
|-------|-------|------|--------|
| 1 | Strand lookups | Low (semantic NOP: t_strand == g_strand) | Medium |
| 2 | Fragment classification | Medium (constant renames, C++ struct changes) | High |
| 3 | Estimation accumulators | Medium (changes gDNA prior computation) | High |
| 4 | Cleanup | Low (dead code removal) | Low |

**Phase 1 is the foundation** — it establishes `t_strand_arr` as the
universal strand source and proves semantic equivalence. Phases 2-4 build
on it.
