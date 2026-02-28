# Collapsed Interval Index — Implementation Plan

**Date:** 2026-02-27  
**Status:** Awaiting approval  

---

## Summary

Replace the per-transcript interval index with a **collapsed interval index**
that stores two types of genic intervals:

1. **EXON** — individual exon boundaries (as today)
2. **TRANSCRIPT** — full transcript genomic span `[tx_start, tx_end)`

Each unique `(ref, start, end, type)` key is stored **once** in cgranges,
mapping to a `frozenset[int]` of transcript indices. The `g_index` field is
dropped entirely — gene indices are derived from `t_to_g_arr` when needed.

### New interval model

For transcript T with exons `(100,500), (2000,2500), (4000,5000)`:

| Interval | Type | Meaning |
|---|---|---|
| `(100, 500)` | EXON | Mature RNA exon boundary |
| `(2000, 2500)` | EXON | Mature RNA exon boundary |
| `(4000, 5000)` | EXON | Mature RNA exon boundary |
| `(100, 5000)` | TRANSCRIPT | Full genic body (exons + introns) |

Fragment resolution examples:

| Fragment | Hits | Interpretation |
|---|---|---|
| F1 `(20, 99)` | INTERGENIC | Intergenic only |
| F2 `(50, 200)` | EXON + TRANSCRIPT + INTERGENIC | mRNA/nRNA/gDNA candidate |
| F3 `(200, 400)` | EXON + TRANSCRIPT | mRNA or nRNA (fully exonic) |
| F4 `(400, 700)` | EXON + TRANSCRIPT | mRNA (overhang) or nRNA (no penalty) |
| F5 `(1000, 1300)` | TRANSCRIPT only | nRNA unique intron evidence |

### Expected performance impact

- cgranges entries: **3,238 → ~950** (71% fewer)
- Per-query Python iterations: **73 → ~14 mean** (81% fewer)
- `query_exon_with_coords` time: **5.9s → ~1.2s** (estimated)
- `resolve_fragment` time: **12.1s → ~8s** (estimated)
- Pipeline total: **43.4s → ~38s** (estimated 5s savings)

---

## Implementation Steps

### Step 1: Modify `IntervalType` enum

**File:** `src/hulkrna/types.py`

```python
class IntervalType(IntEnum):
    EXON = 0
    TRANSCRIPT = 1    # was INTRON — now represents full transcript span
    INTERGENIC = 2
    SJ = 3
    SJ_UNANNOT = 4
```

**Rationale:** Renaming `INTRON → TRANSCRIPT` keeps the same integer value
(1), minimizing impact on persisted Feather files during development. The
semantic change reflects that the interval now represents the full genic
body, not just intronic gaps.

The `RefInterval` NamedTuple loses `g_index`:

```python
class RefInterval(NamedTuple):
    ref: str
    start: int
    end: int
    strand: int = Strand.NONE
    interval_type: int = IntervalType.INTERGENIC
    t_index: int = -1
```

### Step 2: Modify interval generation

**File:** `src/hulkrna/index.py`

Replace `_gen_transcript_intervals()`:

```python
def _gen_transcript_intervals(t: Transcript) -> Iterator[RefInterval]:
    """Yield exon intervals and one transcript span for a single transcript."""
    # Individual exon intervals
    for e in t.exons:
        yield RefInterval(t.ref, e.start, e.end, t.strand,
                          IntervalType.EXON, t.t_index)
    # Single transcript span (covers exons + introns)
    yield RefInterval(t.ref, t.start, t.end, t.strand,
                      IntervalType.TRANSCRIPT, t.t_index)
```

This replaces N-1 intron intervals with 1 transcript span per transcript.
For the PVT1 benchmark: 1,423 intron intervals → 379 transcript spans.

The `_gen_genomic_intervals()` function stays structurally the same — it
still tiles intergenic gaps between transcript clusters and yields
transcript intervals via `_gen_transcript_intervals()`.

The `build_genomic_intervals()` and `build()` remain unchanged
(they just iterate the generator).

`RefInterval` fields drop `g_index`. The `build_splice_junctions()` helper
also drops `g_index` from its `RefInterval` construction.

### Step 3: Build collapsed index in `load()`

**File:** `src/hulkrna/index.py`

Replace the current per-row cgranges insertion with collapsed loading.
The key data structures change from:

```python
# OLD: per-row arrays
self._iv_t_index: np.ndarray    # int[N]  — one t_index per row
self._iv_g_index: np.ndarray    # int[N]  — one g_index per row (DROPPED)
self._iv_type: np.ndarray       # int[N]  — IntervalType per row
```

To:

```python
# NEW: collapsed arrays (native Python lists for zero-conversion overhead)
self._iv_type: list[int]                 # IntervalType per collapsed label
self._iv_t_set: list[frozenset[int]]     # transcript set per collapsed label
```

Loading algorithm:

```python
# Group intervals by (ref, start, end, interval_type) → set of t_indices
collapsed: dict[tuple, set[int]] = {}
for row in iv_df.itertuples(index=False):
    key = (row.ref, row.start, row.end, int(row.interval_type))
    if key not in collapsed:
        collapsed[key] = set()
    t_idx = int(row.t_index)
    if t_idx >= 0:
        collapsed[key].add(t_idx)

# Build cgranges with collapsed entries
cr = cgranges.cgranges()
iv_type_list: list[int] = []
iv_t_set_list: list[frozenset[int]] = []

for label, ((ref, start, end, itype), t_set) in enumerate(collapsed.items()):
    cr.add(ref, start, end, label)
    iv_type_list.append(itype)
    iv_t_set_list.append(frozenset(t_set))

cr.index()
self.cr = cr
self._iv_type = iv_type_list
self._iv_t_set = iv_t_set_list
```

The `_iv_g_index` array is removed entirely. The `_t_exon_intervals` dict
remains — it's built from the raw `iv_df` EXON rows before collapsing.

The `sj_map` changes from `(ref, start, end, strand) → (frozenset[t], frozenset[g])`
to `(ref, start, end, strand) → frozenset[t]`. The `g_set` is dropped.

### Step 4: New query method

**File:** `src/hulkrna/index.py`

Replace `query_exon_with_coords()` with `query_collapsed()`:

```python
def query_collapsed(
    self, exon: GenomicInterval,
) -> list[tuple[int, int, int, frozenset[int]]]:
    """Query collapsed interval index.

    Returns
    -------
    list[tuple[int, int, int, frozenset[int]]]
        List of (hit_start, hit_end, interval_type, transcript_set).
        Interval type is IntervalType.EXON, TRANSCRIPT, or INTERGENIC.
    """
    hits: list[tuple[int, int, int, frozenset[int]]] = []
    for h_start, h_end, label in self.cr.overlap(exon.ref, exon.start, exon.end):
        hits.append((
            h_start,
            h_end,
            self._iv_type[label],       # native Python int
            self._iv_t_set[label],       # frozenset[int]
        ))
    return hits
```

Delete `query_exon()` (already unused in hot path). Delete
`query_exon_with_coords()` — replaced entirely by `query_collapsed()`.

### Step 5: Rewrite `resolve_fragment` inner loop

**File:** `src/hulkrna/resolution.py`

The inner loop changes from iterating per-transcript hits to iterating
per-boundary collapsed hits. Key semantic changes:

**Old model:**
- EXON hit → `block_exon_t.add(t_idx)`, accumulate `_t_exon_bp[t_idx]`
- INTRON hit → `block_intron_t.add(t_idx)`, accumulate `_t_intron_bp[t_idx]`
- Unambig intron BP: subtract global exon union from per-transcript intron intervals

**New model:**
- EXON hit → `block_exon_t |= t_set`, accumulate `_t_exon_bp[t_idx] += bp` for each `t_idx` in set
- TRANSCRIPT hit → `block_transcript_t |= t_set`, accumulate `_t_transcript_bp[t_idx] += bp` for each `t_idx` in set
- Derived: `_t_intron_bp[t] = _t_transcript_bp[t] - _t_exon_bp[t]` (per-transcript, clamped ≥ 0)
- Unambig intron BP: for each TRANSCRIPT hit boundary, call `_subtract_intervals((clipped_lo, clipped_hi), merged_exon)` **once** and distribute to all `t_idx` in the set

The loop structure:

```python
for exon_block in frag.exons:
    exon_strand |= exon_block.strand
    block_start = exon_block.start
    block_end = exon_block.end
    read_length += block_end - block_start

    block_exon_t: frozenset[int] = frozenset()
    block_transcript_t: frozenset[int] = frozenset()

    # Per-block: collect clipped exon intervals for unambig computation
    global_exon_intervals: list[tuple[int, int]] = []
    # Per-block: collect TRANSCRIPT hit boundaries for unambig computation
    transcript_hit_intervals: list[tuple[int, int, frozenset[int]]] = []

    for h_start, h_end, itype, t_set in index.query_collapsed(exon_block):
        clipped_lo = max(block_start, h_start)
        clipped_hi = min(block_end, h_end)
        if clipped_hi <= clipped_lo:
            continue
        bp = clipped_hi - clipped_lo

        if itype == IntervalType.EXON:
            block_exon_t |= t_set
            global_exon_intervals.append((clipped_lo, clipped_hi))
            for t_idx in t_set:
                _t_exon_bp[t_idx] += bp
        elif itype == IntervalType.TRANSCRIPT:
            block_transcript_t |= t_set
            for t_idx in t_set:
                _t_transcript_bp[t_idx] += bp
            transcript_hit_intervals.append(
                (clipped_lo, clipped_hi, t_set)
            )

    exon_t_sets.append(block_exon_t)
    transcript_t_sets.append(block_transcript_t)

    # Unambig intron BP: subtract global exon union from TRANSCRIPT hits
    if transcript_hit_intervals:
        merged_exon = _merge_intervals(global_exon_intervals)
        for clipped_lo, clipped_hi, t_set in transcript_hit_intervals:
            unambig_bp = _subtract_intervals(
                (clipped_lo, clipped_hi), merged_exon
            )
            if unambig_bp > 0:
                for t_idx in t_set:
                    _t_unambig_intron_bp[t_idx] += unambig_bp
```

**Resolution logic changes:**

The `any_intron` check becomes `any_transcript`:

```python
any_exon = any(s for s in exon_t_sets)
any_transcript = any(s for s in transcript_t_sets)
```

The merge logic uses `transcript_t_sets` where it previously used
`intron_t_sets`. The three-state categorization (SPLICED_ANNOT vs genic
vs intergenic) is otherwise identical.

**overlap_bp construction:**

The `overlap_bp` dict values change from `(exon_bp, intron_bp, unambig_intron_bp)`
to the same 3-tuple, but `intron_bp` is now derived:

```python
for t_idx in all_overlap_t:
    exon = _t_exon_bp.get(t_idx, 0)
    transcript = _t_transcript_bp.get(t_idx, 0)
    intron = max(transcript - exon, 0)
    unambig = _t_unambig_intron_bp.get(t_idx, 0)
    overlap_bp[t_idx] = (exon, intron, unambig)
```

This preserves the downstream `(exon_bp, intron_bp, unambig_intron_bp)`
tuple contract — **no changes needed in `buffer.py`, `scan.py`, or
`scoring.py`**.

### Step 6: Update `compute_overlap_profile` (standalone)

**File:** `src/hulkrna/resolution.py`

The standalone `compute_overlap_profile()` function (used by tests) is
updated to use `query_collapsed()` with the same logic as step 5. This
function is no longer called from the hot path but must remain correct
for test compatibility.

### Step 7: Update `sj_map` to drop `g_index`

**File:** `src/hulkrna/index.py`, `src/hulkrna/resolution.py`

Change `sj_map` from `(ref, start, end, strand) → (frozenset[t], frozenset[g])`
to `(ref, start, end, strand) → frozenset[t]`.

In `resolve_fragment`, update the SJ matching block:

```python
# Old:
t_set, g_set = match
sj_t_sets.append(t_set)

# New:
sj_t_sets.append(match)  # match is directly frozenset[int]
```

The strand-agnostic fallback also simplifies (union of t_sets only).

### Step 8: Remove `_iv_g_index` from HulkIndex

**File:** `src/hulkrna/index.py`

Delete `self._iv_g_index` entirely. The `g_index` column is still present
in the Feather interval table (for external tools / debugging) but is
never loaded into HulkIndex at runtime.

### Step 9: Update `exon_lengths_per_gene` and `intron_span_per_gene`

**File:** `src/hulkrna/index.py`

These methods read from the Feather file directly and use
`interval_type == EXON` filtering. Since EXON keeps its value (0), these
methods continue to work. Minor update: the `intron_span_per_gene` logic
(gene_span - exon_length) is unaffected since it uses the gene table, not
interval types.

### Step 10: Update tests

**Files:** `tests/test_types.py`, `tests/test_resolution.py`, `tests/test_index.py`,
`tests/test_scenarios.py`

**test_types.py:**
- `IntervalType.INTRON` → `IntervalType.TRANSCRIPT`
- Update integer value assertions

**test_resolution.py:**
- Mock index returns change from `(t_idx, g_idx, itype, start, end)` tuples
  to `(start, end, itype, frozenset)` tuples
- All `IntervalType.INTRON` references → `IntervalType.TRANSCRIPT`
- Mock index needs `query_collapsed()` instead of `query_exon_with_coords()`
- The `sj_map` mock returns `frozenset[int]` instead of `(frozenset, frozenset)`

**test_index.py:**
- Interval generation tests: verify TRANSCRIPT span intervals instead of
  individual intron intervals
- Verify collapsed index: each label maps to `(itype, frozenset)` not
  `(t_idx, g_idx, itype)`

**test_scenarios.py:**
- End-to-end tests should pass without changes if the overlap_bp contract
  is maintained. The `(exon_bp, intron_bp, unambig_intron_bp)` tuple
  semantics are preserved.

### Step 11: Rebuild test indexes

Any pre-built test index fixtures need to be rebuilt with the new interval
generation (TRANSCRIPT spans instead of INTRON intervals). If tests build
indexes from GTF on-the-fly, they'll automatically get the new format.

---

## Files Changed

| File | Changes |
|---|---|
| `src/hulkrna/types.py` | `INTRON → TRANSCRIPT` in IntervalType; drop `g_index` from RefInterval |
| `src/hulkrna/index.py` | Interval generation (TRANSCRIPT spans); collapsed loading; new `query_collapsed()`; delete `query_exon()`, `query_exon_with_coords()`; sj_map drops g_set; remove `_iv_g_index` |
| `src/hulkrna/resolution.py` | `resolve_fragment()` rewrite for collapsed hits; `compute_overlap_profile()` update; sj_map consumption update |
| `tests/test_types.py` | IntervalType assertions |
| `tests/test_resolution.py` | Mock index API, IntervalType references, sj_map mock |
| `tests/test_index.py` | Interval generation tests, collapsed index tests |
| `tests/test_scenarios.py` | Likely no changes (overlap_bp contract preserved) |

**Files NOT changed:**

| File | Why unchanged |
|---|---|
| `src/hulkrna/buffer.py` | overlap_bp tuple contract preserved: `(exon_bp, intron_bp, unambig_intron_bp)` |
| `src/hulkrna/scan.py` | Consumes buffer's `exon_bp`, `intron_bp`, `unambig_intron_bp` — unchanged |
| `src/hulkrna/scoring.py` | No direct interval type references |
| `src/hulkrna/estimator.py` | No direct interval type references |
| `src/hulkrna/bias.py` | No interval references |
| `src/hulkrna/fragment.py` | No interval references |

---

## Verification Plan

1. **Unit tests:** All 654 tests pass
2. **Benchmark:** PVT1 100k oracle BAM — MAE unchanged (15.86)
3. **Performance:** Profile before/after to measure actual savings
4. **Spot checks:**
   - Fragment landing purely in intron → gets TRANSCRIPT hit only → nRNA pool
   - Fragment spanning exon-intron boundary → gets EXON + TRANSCRIPT → mRNA (with penalty) or nRNA
   - Fragment fully in exon → gets EXON + TRANSCRIPT → mRNA (no penalty)
   - Fragment in intergenic region → INTERGENIC only → intergenic
   - `unambig_intron_bp` correctly excludes bases covered by any EXON hit

---

## Risk Assessment

**Low risk.** The key invariant is that the downstream `overlap_bp` tuple
`(exon_bp, intron_bp, unambig_intron_bp)` is preserved exactly. All changes
are upstream of this contract:
- `exon_bp` comes from EXON hits (unchanged)
- `intron_bp` was from INTRON hits; now derived as `transcript_bp - exon_bp` (same value)
- `unambig_intron_bp` was `_subtract_intervals(intron_iv, exon_union)`; now
  `_subtract_intervals(transcript_iv, exon_union)` — produces the same result
  because the transcript span covers the same intronic bases

The only semantic change is how hits are grouped (per-boundary vs
per-transcript), but the final per-transcript BP values are mathematically
identical.

---

## Execution order

1. **types.py** — IntervalType rename, RefInterval simplification
2. **index.py** — Interval generation, collapsed loading, new query, sj_map
3. **resolution.py** — resolve_fragment, compute_overlap_profile, sj_map
4. **Tests** — Update all test files
5. **Verify** — Run full test suite + benchmark
