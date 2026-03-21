# nRNA Redesign: Implementation Plan

## Core Insight

Nascent RNAs are just transcripts. They happen to overlap mature RNA transcripts and compete with them in the EM. Some are already present in the gene annotation (single-exon transcripts that span multi-exon genes). Others must be created as synthetic shadows because they are biologically real but not annotated.

There should be no special treatment for nascent RNA beyond:
- Recognizing which transcripts serve as nascent RNA equivalents (annotated single-exon transcripts)
- Creating synthetic transcripts where no annotated equivalent exists
- Applying an appropriate prior (sparsity) to synthetic nRNA transcripts

### The near-duplicate span problem

Transcripts with TSS/TES differing by a few bases produce distinct nRNA spans under exact-match dedup. Example: `chr1:+:1000-5000` and `chr1:+:1003-5002` become two separate entities that are 99.9% identical. Under the unified architecture these would become different synthetic nRNA transcripts, splitting fragments between near-identical components.

**Solution**: Tolerance-based merging at index build time. Cluster nearby TSS and TES coordinates independently (within a configurable tolerance in bp), then define the merged nRNA span as the **outer envelope** of each cluster — ensuring the synthetic transcript is never smaller than any contributing transcript's span. Default tolerance: **20 bp**. This is gene-independent, operating purely on coordinates per (ref, strand).


## Why the current architecture has a separate nRNA section

The current EM layout is `[mRNA(0..n_t), nRNA(n_t..n_t+n_nrna), gDNA]` with a `nrna_base` offset. This exists for one reason: **span-level deduplication**. Multiple transcripts sharing the same genomic span (ref, strand, start, end) map to a single nRNA entity, reducing EM components. In a locus with 100 transcripts and 20 unique nRNA spans, this means 120+1 components instead of 200+1.

However, the separate section creates significant architectural complexity:

1. **Union-find fan-out** (`em_solver.cpp:2160-2185`): When a fragment scores against an nRNA component, the union-find must fan out to ALL transcripts sharing that nRNA entity, linking them into the same locus. This is a primary driver of mega-loci — transcripts that happen to share the same genomic span boundaries get forcibly linked even if they share zero fragments.

2. **Post-EM fan-out** (`estimator.py:386-392`): After EM converges, nRNA counts must be redistributed from per-span to per-transcript using equal shares (maximum entropy). This is a lossy approximation that discards within-span information.

3. **Dual indexing** throughout the pipeline: Every module must track both transcript indices and nRNA indices, with `nrna_base` offset arithmetic for component addressing.

4. **Special prior logic**: The all-single-exon prior zeroing in `em_solver.cpp:1413-1424` exists only because the architecture creates nRNA entities for transcripts that don't need them.

### What the EM actually does with the two sections

Examining the EM internals, the sections receive **identical treatment** in almost all respects:
- **OVR prior**: Both mRNA and nRNA get coverage-weighted OVR (`compute_ovr_prior_and_warm_start`, lines 640-735). No section-specific prior.
- **SQUAREM / E-step / M-step**: Completely generic — no component-type awareness.
- **Only difference**: Output scattering — mRNA goes directly to `em_counts[global_t, col]` while nRNA goes to `nrna_em_counts[global_nrna]` then requires Python fan-out.

The nrna_sparsity_alpha (previously 0.9) that was part of the tripartite prior has been removed — the current code uses uniform OVR for all RNA components. There is no remaining EM-level reason to separate mRNA from nRNA.

## Proposed unified architecture

**Promote synthetic nRNA shadows to full transcripts.** Rather than maintaining a parallel nRNA index space, add synthetic transcripts to the transcript list at index build time. They get their own `t_index`, participate in the EM as regular components, and need no special indexing.

### What changes

1. **Index build**: Synthetic nRNA shadows are created as `Transcript` objects with:
   - A single exon spanning the gene's genomic footprint
   - A `is_synthetic_nrna = True` flag
   - Their own `t_index` in the transcript array
   - No separate nRNA table, no `t_to_nrna_arr`, no `nrna_base`

2. **Annotated equivalents**: Single-exon transcripts that fully contain multi-exon transcripts are flagged as `is_nascent_equiv = True`. No synthetic shadow is created for them — they already serve the role.

3. **Fragment scoring**: Synthetic nRNA transcripts score identically to current nRNA scoring — they ARE single-exon, so `ebp = span` for intronic fragments, same fragment-length model, same overhang. The existing `t_span_v <= t_exonic` guard means they never emit nRNA candidates (no nRNA section exists).

4. **Connected components**: No fan-out needed. Synthetic nRNA transcripts connect naturally through fragment co-occurrence in cgranges. When a fragment in T's intron is resolved, cgranges returns both T (via TRANSCRIPT interval) and the synthetic S (via EXON interval). They share the equivalence class → same locus. No forced linking of unrelated transcripts through metadata.

5. **EM output**: Synthetic transcripts get direct per-transcript counts, just like any other transcript. No post-EM fan-out approximation needed.

### What is eliminated

| Removed | Rationale |
|---------|-----------|
| `nrna_base` / `nrna_base_index` | No separate section |
| `t_to_nrna_arr` | No transcript→nRNA mapping |
| `nrna_df` / `nrna.feather` | No nRNA entity table |
| nRNA scoring path in `scoring.cpp` (lines 620-760, 1092-1212) | Synthetic transcripts score as mRNA |
| nRNA WTA / dedup hash map in `scoring.cpp` | Not needed |
| nRNA fan-out in union-find (`em_solver.cpp:2160-2185`) | Natural connectivity via cgranges |
| nRNA→transcript CSR in `locus.py` and `em_solver.cpp` | Not needed |
| `nrna_em_counts_out` accumulator | Direct per-transcript output |
| `_nrna_per_transcript()` fan-out in `estimator.py` | Not needed |
| All-single-exon prior zeroing (`em_solver.cpp:1413-1424`) | Not needed — no dummy entities |
| `nrna_span_`, `nrna_start_` arrays in scorer | Not needed |
| Per-nRNA strand accumulators (`acc_ua`, `acc_us`, `acc_ia`, `acc_is`) | Retain but index by transcript |

### What is gained

1. **No mega-loci from nRNA fan-out**: Connectivity is driven by actual fragment co-occurrence, not metadata. Transcripts that share no fragments won't be linked just because they have the same genomic span.

2. **One index space**: Everything is a transcript. No dual-indexing, no offset arithmetic.

3. **Direct per-transcript output**: No lossy equal-share fan-out.

4. **Simpler code**: The entire nRNA scoring path (~200 lines in `scoring.cpp`), the nRNA CSR construction, the fan-out logic, and the union-find fan-out can be removed.

5. **Same scoring fidelity**: Verified that mRNA scoring against a single-exon transcript produces identical log-likelihoods to current nRNA scoring (same fragment-length model, same overhang computation, same coverage weight when exonic_length = genomic_span).

### Tradeoffs

1. **More EM components per locus**: Without span-level dedup, a locus with 100 multi-exon transcripts and 20 unique nRNA spans has 200+1 components instead of 120+1. However:
   - Many transcripts sharing the same span is the **exact scenario where annotated single-exon equivalents already exist** (many isoforms of the same gene share the same footprint → same nRNA span → likely an annotated single-exon variant exists).
   - For loci WITHOUT annotated equivalents: one synthetic shadow per **unique span** (not per transcript) is still possible — create one synthetic transcript per unique (ref, strand, start, end) and let multiple multi-exon transcripts share it through fragment co-occurrence. This preserves dedup without requiring a separate index space.
   - The increase in components is bounded: in the worst case, it doubles n_t for a locus, but the EM's actual computation is dominated by equivalence classes (not raw component count).

2. **transcript count increases**: The transcript index grows by the number of synthetic nRNA transcripts. For GENCODE this might add ~60-80K entries (one per unique nRNA span that has no annotated equivalent). This is a moderate increase in memory but trivial compared to fragment buffer sizes.

3. **Strand accumulation**: Currently per-nRNA; would need to be per-transcript for synthetic nRNAs. Straightforward — just index by `t_index`.

### Deduplication strategy for synthetic nRNAs

To limit component count growth, create **one synthetic transcript per unique merged nRNA span**, not one per multi-exon transcript. Near-duplicate spans (TSS/TES within tolerance) are clustered into a single merged span using the outer envelope (minimum start, maximum end), ensuring the synthetic transcript fully contains every contributing transcript.

Multiple multi-exon transcripts mapping to the same merged span connect to the same synthetic transcript through fragment co-occurrence (cgranges resolves the same intronic fragment to the same synthetic transcript, regardless of which multi-exon transcript generated the nRNA signal).

This gives us span-level dedup WITHOUT a separate index space:
- Current: `t_to_nrna_arr[t0] = t_to_nrna_arr[t1] = nrna_3` (shared nRNA entity, separate index space)
- Proposed: Both t0 and t1 share fragments with synthetic transcript s0 (same genomic span, same index space)

The union-find naturally connects them through shared equivalence classes rather than forced fan-out.

## Tolerance-based nRNA span merging

### Problem

Exact-match dedup on `(ref, strand, start, end)` produces too many near-duplicate nRNA spans. Transcripts from the same biological gene locus often have TSS/TES differing by a few bases due to annotation ambiguity, alternative promoters, or polyA site variation. Without merging, each becomes a separate synthetic transcript competing for fragments in the EM.

### Algorithm: coordinate-based TSS/TES clustering

The algorithm is gene-independent, operating purely on genomic coordinates per (ref, strand) group. It clusters transcript start sites (TSS) and end sites (TES) separately, then combines the outer-envelope boundaries to define merged nRNA spans.

**Parameters:**
- `nrna_merge_tolerance`: int, default 20 bp. Maximum distance between coordinates in the same cluster.

**Step 1 — Collect multi-exon transcript coordinates by (ref, strand):**

```
For each multi-exon transcript t:
  Group by (t.ref, t.strand)
  Record (t.start, t.end, t.t_index)
```

**Step 2 — Cluster TSS (start sites) within each (ref, strand) group:**

Sort starts ascending. Sweep left-to-right:
```
cluster_min = starts[0]
for each start in sorted order:
  if start - cluster_min > tolerance:
    Close current cluster → representative_start = cluster_min  (outer envelope)
    Open new cluster: cluster_min = start
  else:
    Continue in current cluster (cluster_min stays at minimum)
Close final cluster → representative_start = cluster_min
```

For + strand transcripts: the outer envelope TSS is the **minimum** start in the cluster (leftmost). For - strand: the outer envelope TSS is the **maximum** start (rightmost, which is the outer boundary). However, since we always want the span that is >= every contributing transcript's span, we use: `representative_start = min(cluster starts)` regardless of strand. Same logic applies to TES: `representative_end = max(cluster ends)`.

This guarantees the merged span fully contains every contributing transcript.

**Step 3 — Cluster TES (end sites) within each (ref, strand) group:**

Same sweep algorithm on ends, sorted ascending:
```
cluster_max = ends[0]
for each end in sorted order:
  if end - cluster_max_start > tolerance:
    Close current cluster → representative_end = cluster_max
    Open new cluster
  else:
    cluster_max = max(cluster_max, end)
Close final cluster → representative_end = cluster_max
```

Where `cluster_max_start` tracks the first (smallest) end in the current cluster for gap detection.

**Step 4 — Define merged nRNA spans:**

Each multi-exon transcript maps to a merged span:
```
merged_span = (ref, strand, representative_start[t], representative_end[t])
```

Dedup on this merged key: one synthetic transcript per unique merged span.

### Why cluster TSS and TES independently?

TSS and TES variation are biologically independent (promoter choice vs polyA site choice). Two transcripts might share a TSS cluster but differ in TES by >tolerance, or vice versa. Independent clustering handles this correctly:

```
t1:  chr1:+:1000-5000  → TSS cluster A (min=1000), TES cluster X (max=5000)
t2:  chr1:+:1003-5002  → TSS cluster A (min=1000), TES cluster X (max=5002)
t3:  chr1:+:1001-8000  → TSS cluster A (min=1000), TES cluster Y (max=8000)

Merged spans:
  (chr1, +, 1000, 5002) for t1, t2  → one synthetic
  (chr1, +, 1000, 8000) for t3      → one synthetic
```

### Complexity

- Sorting: O(N log N) per (ref, strand) group
- Clustering sweep: O(N)
- Total: O(N log N) dominated by sorting, negligible compared to index build

## Detection algorithm: annotated equivalents + tolerance merging

```
Step 0 — Tolerance-based TSS/TES clustering:
  For each (ref, strand) group of multi-exon transcripts:
    Cluster TSS within tolerance → map each start to representative_start (cluster min)
    Cluster TES within tolerance → map each end to representative_end (cluster max)
    Compute merged_span = (ref, strand, representative_start, representative_end)
  Dedup merged spans → set of unique nRNA spans needing coverage

Step 1 — Exact match (annotated equiv detection):
  For each single-exon transcript s:
    If s.span matches any merged_span exactly:
      Mark s as is_nascent_equiv = True
      Mark that merged_span as covered

Step 2 — Full containment:
  Group single-exon transcripts by (ref, strand), sorted by start
  For each uncovered merged_span:
    Binary search for single-exon tx with start ≤ merged_span.start
    Among those, check if any has end ≥ merged_span.end
    If found: mark that tx as is_nascent_equiv, mark span as covered

Step 3 — Create synthetics:
  For each uncovered merged_span:
    Create synthetic Transcript with one exon = [merged_span.start, merged_span.end]
    Assign new t_index, flag is_synthetic_nrna = True
```

Containment criterion: `S.start ≤ merged_span.start AND merged_span.end ≤ S.end`, same ref and strand.

## Component count examples

**Locus with multi-exon T and single-exon S (annotated equiv):**

| | Current | Proposed |
|---|---------|----------|
| Components | mRNA(T), mRNA(S), nRNA(shadow), gDNA = **4** | mRNA(T), mRNA(S), gDNA = **3** |
| Intronic fragments | Score against mRNA(S) AND nRNA(shadow) | Score against mRNA(S) only |
| Siphon risk | nRNA siphons from mRNA(S) | Eliminated |

**Pure single-exon gene (S only):**

| | Current | Proposed |
|---|---------|----------|
| Components | mRNA(S), nRNA(shadow), gDNA = **3** | mRNA(S), gDNA = **2** |
| nRNA shadow | Created then zeroed | Not created |

**Multi-exon T without annotated equiv:**

| | Current | Proposed |
|---|---------|----------|
| Components | mRNA(T), nRNA(shadow), gDNA = **3** | mRNA(T), mRNA(synth_S), gDNA = **3** |
| nRNA handling | Separate index space | Regular transcript |

**100 multi-exon transcripts, 20 unique spans, no annotated equiv:**

| | Current | Proposed |
|---|---------|----------|
| Components | 100 mRNA + 20 nRNA + 1 gDNA = **121** | 100 mRNA + 20 synth + 1 gDNA = **121** |
| Connectivity | Forced fan-out from shared nRNA | Natural fragment co-occurrence |

## Phases

### Phase A: Diagnostic — COMPLETED

**Script**: `scripts/debug/nrna_diagnostic.py`

**Dataset**: GENCODE-based index — 254,461 transcripts (26,484 single-exon / 227,977 multi-exon). Current architecture: 246,058 nRNA entities.

#### Summary table

| Tol (bp) | Exact spans | Merged spans | Reduction % | Annot. equiv | Synthetic needed |
|----------|-------------|-------------|-------------|-------------|-----------------|
| 0 | 219,610 | 219,610 | 0.0% | 875 | 218,735 |
| 5 | 219,610 | 213,636 | 2.7% | 767 | 212,869 |
| 10 | 219,610 | 209,589 | 4.6% | 709 | 208,880 |
| **20** | **219,610** | **203,696** | **7.2%** | **644** | **203,052** |
| 50 | 219,610 | 192,903 | 12.2% | 570 | 192,333 |
| 100 | 219,610 | 183,291 | 16.5% | 492 | 182,799 |

#### Component count comparison

| Architecture | Total components |
|-------------|-----------------|
| Current (246K nRNA entities) | 500,519 |
| Proposed tol=0bp (218K synth) | 473,196 (−27K) |
| **Proposed tol=20bp (203K synth)** | **457,513 (−43K)** |
| Proposed tol=100bp (183K synth) | 437,260 (−63K) |

#### Key observations

1. **Annotated equivalents are rare**: Only 644–875 spans (0.3–0.4%) are covered by existing single-exon transcripts. The vast majority of nRNA spans need synthetic transcripts.
2. **Tolerance-based merging is effective**: At tol=20bp, 11,029 merged spans consolidate multiple exact spans (7.2% reduction). Diminishing returns above ~50bp.
3. **Net component reduction**: Even at tol=0, the proposed architecture eliminates 27K components (single-exon transcripts no longer need dummy nRNA). At tol=20bp: 43K fewer components.
4. **Near-duplicate examples**: LINC00635 has 25 exact spans collapsed to 1 (64 transcripts with TSS/TES differing by <10bp). SNHG5, FANCL, POLD2 show similar patterns.
5. **20bp tolerance recommended**: Good balance — 7.2% span reduction, well within biological TSS/TES variation, no risk of merging genuinely distinct regulatory elements.

### Phase B: Index changes

1. Add `nrna_merge_tolerance` parameter (default 20 bp) to index build config
2. Add `is_synthetic_nrna` and `is_nascent_equiv` fields to `Transcript` dataclass
3. Modify `compute_nrna_table()` → replace with `create_nrna_transcripts()`:
   - TSS/TES clustering with configurable tolerance
   - Detect annotated equivalents (exact match + containment against merged spans)
   - Create synthetic `Transcript` objects for remaining uncovered spans
   - Return augmented transcript list (no separate nRNA table)
3. Remove `nrna.feather` serialization; synthetic transcripts are part of `transcripts.feather`
4. Remove `t_to_nrna_arr` from index
5. Update `_gen_transcript_intervals()` — synthetic transcripts generate intervals like any single-exon transcript

### Phase C: Remove nRNA-specific machinery

1. **`scoring.cpp`**: Remove entire nRNA scoring path (two template paths: ~200 lines). Remove nRNA dedup hash map. Remove `nrna_base_`, `t_to_nrna_`, `nrna_span_`, `nrna_start_` members. Synthetic transcripts score as mRNA automatically.
2. **`em_solver.cpp`**: Remove nRNA section from `extract_locus_sub_problem` (unique nRNA collection, nRNA→transcript CSR, all-single-exon zeroing). Remove `nrna_base` from batch_locus_em. Remove nRNA fan-out from union-find.
3. **`locus.py`**: Remove nRNA→transcript CSR construction in `build_loci()`. Remove nRNA-specific logic in `build_locus_em_data()`.
4. **`estimator.py`**: Remove `nrna_em_counts`, `_nrna_per_transcript()` fan-out. Synthetic transcripts have direct per-transcript counts.
5. **`scoring.py`**: Remove nRNA array parameters from scorer construction.
6. **Strand accumulators**: Rework per-nRNA accumulators to per-transcript (or remove if not needed for synthetic transcripts).

### Phase D: Prior for synthetic nRNAs

Synthetic nRNA transcripts may benefit from a sparsity prior to prevent them from absorbing too many fragments in the absence of real nascent RNA signal. Options:
- Apply existing OVR as-is (coverage-weighted prior treats synthetics like any transcript)
- Add a per-transcript prior flag: `is_synthetic_nrna → lower prior`
- Evaluate in benchmarks whether any special treatment is needed

### Phase E: Validation

1. Recompile: `pip install --no-build-isolation -e .`
2. Run all tests: `pytest tests/ -v`
3. Update golden outputs: `pytest tests/ --update-golden`
4. Run benchmark suite: compare nRNA siphon metrics, mRNA/gDNA relative error
5. Verify locus sizes decreased (no forced fan-out)
6. Spot-check specific genes (MALAT1-like) for correct behavior

## Files affected

| File | Role | Change |
|------|------|--------|
| `src/rigel/index.py` | nRNA table → synthetic transcripts | Major rewrite of `compute_nrna_table()` |
| `src/rigel/transcript.py` | Transcript dataclass | Add `is_synthetic_nrna`, `is_nascent_equiv` fields |
| `src/rigel/native/scoring.cpp` | Fragment scoring | Major: remove ~200 lines of nRNA scoring |
| `src/rigel/native/em_solver.cpp` | EM solver | Major: remove nRNA section, fan-out, CSR |
| `src/rigel/locus.py` | Locus construction | Remove nRNA CSR, simplify `build_loci` |
| `src/rigel/estimator.py` | Abundance estimation | Remove nRNA fan-out, simplify output |
| `src/rigel/scoring.py` | Python scorer interface | Remove nRNA parameters |
| `src/rigel/scan.py` | Scored fragments | Remove `nrna_base_index` |
| `src/rigel/pipeline.py` | Orchestrator | Simplify — fewer parameters to thread |
| `src/rigel/native.py` | C++ binding interface | Remove nRNA parameters |
