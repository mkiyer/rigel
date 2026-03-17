# Fix: Multi-Chromosome gDNA Length Normalization in Mega-Loci

## Problem Statement

When minimap2 multimappers create mega-loci spanning multiple chromosomes (e.g., 217K transcripts across 43K genes on 41 chromosomes), the gDNA component's bias profile length (`locus_span`) is computed as `max(t_ends) - min(t_starts)` **across all chromosomes**. Since transcript coordinates are chromosome-relative, this produces a nonsensical value.

### Current Behavior

For the benchmark mega-locus (217,246 transcripts on 41 chromosomes):

| Metric | Current Value | Correct Value |
|--------|--------------|---------------|
| `locus_span` | 248,936,396 (≈249 Mbp) | 1,741,698,323 (≈1.74 Gbp) |
| Per-fragment gDNA penalty | −log(249M) = −19.3 nats | −log(1.74G) = −21.3 nats |
| Penalty deficit vs correct | **+2.0 nats** (under-penalized) | baseline |

The per-fragment bias correction (`apply_bias_correction_uniform`) subtracts `log(max(L − frag + 1, 1))` from each gDNA candidate's log-likelihood. The 2-nat under-penalty is enough to tip the EM equilibrium:

- **Current equilibrium**: gDNA wins by +0.74 nats over average RNA candidate → gDNA absorbs 44.5% of mega-locus
- **Corrected penalty**: gDNA loses by −1.2 nats → gDNA share drops substantially
- **Truth**: gDNA should be ~20% of total (2M gDNA / 12M fragments)

### Root Cause Chain

```
1. Multimappers connect transcripts on different chromosomes via union-find
2. Mega-locus formed: 217K transcripts, 43K genes, 41 chromosomes
3. locus_span = max(t_ends) - min(t_starts) = 249M (chromosome-blind)
4. Correct value: union of per-chromosome genomic footprints = 1.74G
5. gDNA bias correction under-penalizes by 2 nats per fragment
6. θ_gDNA (universal candidate for all unspliced fragments) accumulates
7. EM converges to 44.5% gDNA rate in mega-locus (should be ~20%)
8. 4.72M gDNA predicted vs 2M truth — 99.2% of the error in ONE locus
```

## How Length Normalization Currently Works

### Per-Fragment Bias Correction (the active mechanism)

`apply_bias_correction_uniform()` in `em_solver.cpp` mutates log-likelihoods in-place before the EM loop:

```cpp
for (size_t i = 0; i < n_candidates; ++i) {
    int64_t frag_len = tx_ends[i] - tx_starts[i];
    int64_t prof_len = profile_lengths[t_indices[i]];
    int64_t eff_len = prof_len - frag_len + 1;
    if (eff_len < 1) eff_len = 1;
    log_liks[i] -= std::log(static_cast<double>(eff_len));
}
```

For each candidate alignment:
- **mRNA**: `prof_len` = exonic transcript length → `eff_len ≈ tx_len − frag_proj + 1`
- **nRNA**: `prof_len` = genomic span of nRNA → `eff_len ≈ nrna_span − frag_proj + 1`
- **gDNA**: `prof_len` = `locus_span` → `eff_len ≈ locus_span − footprint + 1`

For gDNA, `tx_start=0` and `tx_end=genomic_footprint`, so `frag_len = footprint`.

### Per-Component E-step (disabled)

```cpp
// In map_em_step:
log_weights[i] = log(theta[i] + eps) - log_eff_len[i];
// But log_eff_len is always 0.0 (all effective lengths = 1.0)
```

This is intentionally disabled because the per-fragment correction supersedes it. Both approaches are mathematically equivalent for uniform positional models, but the per-fragment approach is more precise when fragment lengths vary.

### Where `locus_span` Feeds Into Bias Profiles

**Python path** (`build_locus_em_data` in `locus.py`):
```python
locus_start = int(t_starts.min())    # min across ALL chromosomes!
locus_end = int(t_ends.max())        # max across ALL chromosomes!
locus_span = float(locus_end - locus_start)
...
bias_profiles[gdna_idx] = int(locus_span)
```

**C++ batch path** (`batch_locus_em` in `em_solver.cpp`):
```cpp
int64_t locus_start = all_t_starts[t_arr[0]];
int64_t locus_end = all_t_ends[t_arr[0]];
for (int i = 1; i < n_t; ++i) {
    if (ts < locus_start) locus_start = ts;
    if (te > locus_end) locus_end = te;
}
int64_t locus_span = locus_end - locus_start;
...
sub.bias_profiles[sub.gdna_idx] = locus_span;
```

Both paths: chromosome-blind `max - min` across all transcripts.

## Implementation Plan

### Step 1: Precompute Per-Locus Union Genomic Footprints (Python)

Add a function in `locus.py` that computes the correct gDNA genomic footprint for each locus by computing the union of transcript intervals per chromosome.

```python
def compute_locus_gdna_spans(
    loci: list[Locus],
    index: TranscriptIndex,
) -> np.ndarray:
    """Compute gDNA genomic footprint per locus.
    
    For each locus, groups transcripts by reference (chromosome),
    merges overlapping intervals within each reference, and sums
    the merged spans across references.
    
    Returns int64[n_loci] array of gDNA span values.
    """
    t_starts = index.t_df["start"].values
    t_ends = index.t_df["end"].values
    t_refs = index.t_df["ref"].values
    
    gdna_spans = np.empty(len(loci), dtype=np.int64)
    
    for i, locus in enumerate(loci):
        t_arr = locus.transcript_indices
        if len(t_arr) == 0:
            gdna_spans[i] = 0
            continue
            
        # Group by chromosome
        refs = t_refs[t_arr]
        starts = t_starts[t_arr]
        ends = t_ends[t_arr]
        
        total_span = 0
        for ref in np.unique(refs):
            mask = refs == ref
            ref_starts = starts[mask]
            ref_ends = ends[mask]
            
            # Sort by start, merge overlapping intervals
            order = np.argsort(ref_starts)
            s_starts = ref_starts[order]
            s_ends = ref_ends[order]
            
            cur_start = s_starts[0]
            cur_end = s_ends[0]
            for j in range(1, len(s_starts)):
                if s_starts[j] <= cur_end:
                    cur_end = max(cur_end, s_ends[j])
                else:
                    total_span += cur_end - cur_start
                    cur_start = s_starts[j]
                    cur_end = s_ends[j]
            total_span += cur_end - cur_start
        
        gdna_spans[i] = max(total_span, 1)
    
    return gdna_spans
```

**Optimization**: For normal (single-chromosome, <100 tx) loci, the current `max - min` is correct and fast. Only need the full union-merge logic for multi-chromosome loci.

### Step 2: Pass Pre-Computed gDNA Spans to C++ Batch EM

Currently the C++ path computes `locus_span` inline. We'll pre-compute in Python and pass as an array.

**Changes to `estimator.py`** → `run_batch_locus_em()`:
```python
# Add gdna_spans parameter
gdna_spans = compute_locus_gdna_spans(loci, index)

# Pass to C++ alongside gdna_inits
total_gdna_em, locus_mrna, locus_nrna, locus_gdna = _batch_locus_em(
    ...,
    gdna_spans,  # NEW: int64[n_loci]
    ...,
)
```

**Changes to `em_solver.cpp`** → `batch_locus_em()`:
- Add `i64_1d gdna_spans` parameter
- In the per-locus loop, replace the inline `locus_span` computation:

```cpp
// BEFORE:
int64_t locus_start = all_t_starts[t_arr[0]];
int64_t locus_end = all_t_ends[t_arr[0]];
for (int i = 1; i < n_t; ++i) { ... }
int64_t locus_span = locus_end - locus_start;

// AFTER:
int64_t locus_span = all_gdna_spans[li];  // pre-computed, chromosome-aware
```

### Step 3: Update Python Path (`build_locus_em_data`)

In `locus.py` → `build_locus_em_data()`, accept the pre-computed gDNA span:

```python
def build_locus_em_data(
    locus, em_data, estimator, index, mean_frag, gdna_init,
    gdna_span,  # NEW: pre-computed union genomic footprint
    *, _cache=None,
):
    ...
    # BEFORE:
    # locus_span = float(locus_end - locus_start)
    
    # AFTER:
    locus_span = float(gdna_span)
    ...
```

### Step 4: Update Nanobind Bindings

The C++ function signature needs to accept the new `gdna_spans` array. Update the nanobind binding in `em_solver.cpp`.

### Step 5: Verify nRNA Bias Profiles

The nRNA spans are computed per-transcript in the C++ batch path:
```cpp
sub.bias_profiles[n_t + n] = all_t_ends[gt] - all_t_starts[gt];
```

This uses a single transcript's `end - start`, which is correct because each transcript lives on one chromosome. **No fix needed for nRNA spans**, but we should add a comment noting this is single-chromosome safe.

### Step 6: Update Tests and Run Benchmark

1. Run the full test suite to ensure no regressions
2. Run the minimap2 benchmark to verify gDNA siphon reduction
3. Update golden outputs if EM results change (they should)

## Files to Modify

| File | Change |
|------|--------|
| `src/rigel/locus.py` | Add `compute_locus_gdna_spans()`, update `build_locus_em_data()` signature |
| `src/rigel/estimator.py` | Compute and pass `gdna_spans` to C++ |
| `src/rigel/native/em_solver.cpp` | Accept `gdna_spans` array, replace inline locus_span computation |
| `src/rigel/pipeline.py` | Thread `gdna_spans` from locus builder to EM |

## Expected Impact

The per-fragment gDNA penalty increases from −19.3 to −21.3 nats (for the mega-locus). This 2-nat shift converts gDNA from a +0.74 nat winner to a −1.2 nat loser in the E-step equilibrium. The EM will converge to a substantially lower gDNA fraction.

Precise prediction requires running the EM, but the directional effect is clear: fewer fragments will be assigned to gDNA, with the freed mass flowing primarily to nRNA (which has a 3.3M deficit in the mega-locus).

## Risk Assessment

- **Normal loci (single chromosome)**: Unaffected — union of intervals on one chromosome equals `max(ends) - min(starts)` for dense loci
- **Small multi-chromosome loci**: Improved accuracy, but these are rare and small
- **Mega-loci**: Major improvement — the primary target of this fix
- **Performance**: `compute_locus_gdna_spans` iterates loci once with interval merge. Cost is O(Σ n_tx × log(n_tx)) per locus for sorting, negligible vs EM time
