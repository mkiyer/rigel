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


## UNAMBIG_INTRON and evidence gating — investigation results

### Original intent

The *intended* design was to require unambiguous intronic evidence before emitting nRNA candidates. The UNAMBIG_INTRON interval type marks genomic regions that are intronic for a transcript and do NOT overlap any exon from any transcript on either strand. If a fragment overlaps such a region, it constitutes evidence for nascent RNA (or gDNA) — exonic unspliced fragments alone are ambiguous.

### What was actually implemented

A thorough code audit (March 2026) revealed the following data flow:

**Index → resolve chain (FUNCTIONAL)**:
- `_gen_cluster_unambig_intron_intervals()` correctly identifies unambiguous intronic regions by subtracting the global exon union from each transcript's introns
- UNAMBIG_INTRON intervals are stored in cgranges with `nrna_idx` as the label (in the `t_index` field of `AnnotatedInterval`)
- The resolve step in `resolve_context.h` accumulates `t_unambig_intron_bp` from cgranges overlaps

**Index mismatch bug**: The resolve scratch buffer is written at `scratch.t_unambig_intron_bp[nrna_idx]` (during overlap accumulation) but read at `scratch.t_unambig_intron_bp[t_index]` (during output assembly). Since `nrna_idx != t_index` in general, the values reaching downstream are unreliable.

**Scoring consumption (WAS DEAD CODE — NOW REMOVED)**:
- The ONLY consumer of `unambig_intron_bp` in scoring was `has_ui = (ui_bp[k] > 0)`, which gated writes to `acc_ia`/`acc_is` (intronic sense/antisense accumulators)
- The `transcript_intronic_sense/antisense` arrays were allocated in `estimator.py`, populated by C++ scoring, but **never read by any downstream production code** — no EM, no priors, no output consumed them
- These were dead code, confirmed by the existing TODO.md note: "intronic_sense/antisense accumulators are populated by the C++ scoring pass but are not used"

**nRNA candidate emission has NO evidence gate**:
- The nRNA scoring path (`scoring.cpp` lines 1092-1210 singlemapper, 620-760 multimapper) does NOT check `unambig_intron_bp` at all
- The actual gates for nRNA candidate emission are:
  1. `stype != SPLICE_SPLICED_ANNOT` — skip splice-annotated fragments
  2. `t_span_v <= t_exonic` — skip single-exon transcripts (no intronic space)
  3. `span_bp <= 0` — skip no-overlap

### Dead code cleaned up (March 2026)

The following dead code has been removed:
- `intronic_sense`/`intronic_antisense` parameters from `fused_score_buffer()` (scoring.cpp)
- `acc_is`/`acc_ia` accumulators and `has_ui` gating logic in `score_chunk_impl()`
- `ui_bp` field from `ChunkPtrs` struct and local variable in scoring
- `transcript_intronic_sense/antisense` arrays in `estimator.py`
- Corresponding parameters in `scan.py` call
- `unambig_intron_bp` removed from scoring tuple (`to_scoring_arrays()`)
- Diagnostic reads of dead accumulators in test files

### What remains

- **UNAMBIG_INTRON interval generation** in index.py — preserved (correct infrastructure)
- **Resolution of `unambig_intron_bp`** in resolve_context.h — preserved (data still flows into buffer)
- **Buffer transport** of `unambig_intron_bp` in buffer.py — preserved (serialized in feather files)
- The nrna_idx/t_index mismatch bug exists but is harmless since no production code consumes the values
- Future evidence gating can be implemented by adding a guard in scoring.cpp that checks `unambig_intron_bp > 0` before emitting nRNA candidates

### Design decision for unified architecture

The unified architecture does NOT add evidence gating. This matches the current behavior: nRNA candidates are emitted based on transcript geometry (multi-exon) and fragment overlap, not on unambiguous intronic evidence. Evidence gating can be added as a future enhancement when needed.


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
| nRNA hard overhang gate / dedup hash map in `scoring.cpp` | Not needed |
| nRNA fan-out in union-find (`em_solver.cpp:2160-2185`) | Natural connectivity via cgranges |
| nRNA→transcript CSR in `locus.py` and `em_solver.cpp` | Not needed |
| `nrna_em_counts_out` accumulator | Direct per-transcript output |
| `_nrna_per_transcript()` fan-out in `estimator.py` | Not needed |
| All-single-exon prior zeroing (`em_solver.cpp:1413-1424`) | Not needed — no dummy entities |
| `nrna_span_`, `nrna_start_` arrays in scorer | Not needed |
| Per-nRNA unspliced strand accumulators (`acc_ua`, `acc_us`) | Rework to per-transcript indexing |

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

2. **Transcript count increases**: The transcript index grows by the number of synthetic nRNA transcripts. For GENCODE this adds ~203K entries with tol=20bp. This is a moderate increase in memory but trivial compared to fragment buffer sizes.

3. **Splice-annotated fragments**: In the unified architecture, splice-annotated fragments will resolve to the synthetic (since its EXON interval covers the region). They will score as mRNA candidates against the synthetic at competitive likelihoods. The current architecture explicitly skips nRNA scoring for splice-annotated fragments. This is a behavioral change — but the EM should naturally assign these fragments to the real multi-exon transcripts (which have higher likelihoods for spliced reads due to proper exon structure). Benchmark validation is needed to confirm.


### Deduplication strategy for synthetic nRNAs

To limit component count growth, create **one synthetic transcript per unique merged nRNA span**, not one per multi-exon transcript. Near-duplicate spans (TSS/TES within tolerance) are clustered into a single merged span using the outer envelope (minimum start, maximum end), ensuring the synthetic transcript fully contains every contributing transcript.

Multiple multi-exon transcripts mapping to the same merged span connect to the same synthetic transcript through fragment co-occurrence (cgranges resolves the same intronic fragment to the same synthetic transcript, regardless of which multi-exon transcript generated the nRNA signal).

This gives us span-level dedup WITHOUT a separate index space:
- Current: `t_to_nrna_arr[t0] = t_to_nrna_arr[t1] = nrna_3` (shared nRNA entity, separate index space)
- Proposed: Both t0 and t1 share fragments with synthetic transcript s0 (same genomic span, same index space)

The union-find naturally connects them through shared equivalence classes rather than forced fan-out.


## Tolerance-based nRNA span merging

### Algorithm: coordinate-based TSS/TES clustering

The algorithm is gene-independent, operating purely on genomic coordinates per (ref, strand) group. It clusters transcript start sites (TSS) and end sites (TES) separately, then combines the outer-envelope boundaries to define merged nRNA spans.

**Parameters:**
- `nrna_merge_tolerance`: int, default 20 bp. Maximum distance between coordinates in the same cluster.

The algorithm clusters starts and ends independently, using the outer envelope (min start, max end) to ensure the merged span fully contains every contributing transcript. TSS and TES variation are biologically independent (promoter choice vs polyA site choice), so independent clustering is correct.

### Complexity

- Sorting: O(N log N) per (ref, strand) group
- Clustering sweep: O(N)
- Total: O(N log N) dominated by sorting, negligible compared to index build


## Phase A: Diagnostic — COMPLETED

**Script**: `scripts/debug/nrna_diagnostic.py`

**Dataset**: GENCODE-based index — 254,461 transcripts (26,484 single-exon / 227,977 multi-exon). Current architecture: 246,058 nRNA entities.

### Summary table

| Tol (bp) | Exact spans | Merged spans | Reduction % | Annot. equiv | Synthetic needed |
|----------|-------------|-------------|-------------|-------------|-----------------|
| 0 | 219,610 | 219,610 | 0.0% | 875 | 218,735 |
| 5 | 219,610 | 213,636 | 2.7% | 767 | 212,869 |
| 10 | 219,610 | 209,589 | 4.6% | 709 | 208,880 |
| **20** | **219,610** | **203,696** | **7.2%** | **644** | **203,052** |
| 50 | 219,610 | 192,903 | 12.2% | 570 | 192,333 |
| 100 | 219,610 | 183,291 | 16.5% | 492 | 182,799 |

### Component count comparison

| Architecture | Total components |
|-------------|-----------------|
| Current (246K nRNA entities) | 500,519 |
| Proposed tol=0bp (218K synth) | 473,196 (−27K) |
| **Proposed tol=20bp (203K synth)** | **457,513 (−43K)** |
| Proposed tol=100bp (183K synth) | 437,260 (−63K) |


## Phased Implementation Plan

### Phase B: Index changes

**Status**: COMPLETED (2026-03-22). All 1011 tests pass.

**What was done**:
- B.1: `TranscriptIndex.build()` calls `create_nrna_transcripts()`, assigns `t_index`/`g_index` to synthetics, appends to transcript list BEFORE `compute_nrna_table()`
- B.2: `compute_nrna_table()` operates on augmented transcript list (no changes needed — automatic)
- B.3: `TranscriptIndex.load()` has backward compat for `is_synthetic_nrna`/`is_nascent_equiv` columns
- B.4: UNAMBIG_INTRON was fully removed in dead code cleanup (moot)
- B.5: Fixed 139 test failures:
  - Benchmark accountability: Added `n_synthetic_observed` to `BenchmarkResult` to account for EM mass assigned to synthetic transcript mRNA components
  - Updated `total_rna_observed` property and `assert_accountability` helpers
  - Region partition: Synthetic nRNA exons cover full gene spans, merging intronic regions (11 → 5 regions in mini_index)
  - Fragment resolution: Synthetic transcripts appear in `t_inds` for unspliced reads
  - Count assertions: Updated hardcoded transcript/gene counts across all test files
  - Golden outputs regenerated

**Code already implemented** (in `src/rigel/index.py`):
- `NRNA_MERGE_TOLERANCE = 20` constant
- `_cluster_coordinates(coords, tolerance)` helper
- `create_nrna_transcripts(transcripts, tolerance=20)` function

**Code already implemented** (in `src/rigel/transcript.py`):
- `is_synthetic_nrna: bool = False` field on `Transcript`
- `is_nascent_equiv: bool = False` field on `Transcript`
- Both fields included in `to_dict()` for serialization

**Remaining work for Phase B**:

1. **B.1**: Modify `TranscriptIndex.build()` to call `create_nrna_transcripts()`, assign `t_index`/`g_index` to synthetics, append to transcript list
2. **B.2**: Keep `compute_nrna_table()` operating on the augmented transcript list for backward compatibility — the existing scoring/EM pipeline needs `nrna_idx`, `t_to_nrna_arr`, `nrna.feather` until Phase C removes them
3. **B.3**: Add backward compat in `TranscriptIndex.load()` for the new boolean columns (old indices won't have them)
4. **B.4**: Update UNAMBIG_INTRON interval generation — currently uses `nrna_idx` as the label; need to ensure synthetics get proper nrna_idx from `compute_nrna_table()` (they will, automatically, since they're in the transcript list)
5. **B.5**: Update tests that assert hardcoded transcript counts / IDs — these will increase by the number of synthetics

**Ordering constraint**: Synthetics must be appended BEFORE `compute_nrna_table()` is called, BEFORE intervals are generated, and BEFORE transcripts.feather is written.

**Key observation**: Phase B can be deployed independently. The synthetics will participate in cgranges (getting EXON + TRANSCRIPT intervals), get their own `nrna_idx`, and enter the existing scoring/EM pipeline as regular transcripts. Their nRNA candidates will also be emitted via the existing nRNA scoring path (since `compute_nrna_table()` gives them nrna_idx entries). This is correct — both mRNA and nRNA candidates will be emitted for synthetics, and the EM will handle them.

### Phase C: Remove nRNA-specific machinery

**Prerequisite**: Phase B deployed and validated.

This is the major simplification phase. Remove the separate nRNA index space entirely.

**C.1 — scoring.cpp**: Remove entire nRNA scoring path.
- Delete nRNA dedup hash map, nRNA hard overhang gate logic in `score_chunk_impl` (singlemapper path, lines ~1092-1210)
- Delete nRNA handling in `flush_mm_group` (multimapper path, lines ~620-760)
- Remove `nrna_base_`, `t_to_nrna_`, `nrna_span_`, `nrna_start_` member variables from `NativeFragmentScorer`
- Remove corresponding constructor parameters and nanobind bindings
- **Net**: ~200 lines of C++ removed

**C.2 — em_solver.cpp**: Remove nRNA section from EM.
- Delete nRNA collection/CSR construction in `extract_locus_sub_problem()`
- Delete all-single-exon prior zeroing (lines ~1413-1424)
- Delete `nrna_base` parameter from `batch_locus_em()` signature
- Delete `nrna_em_counts_out` accumulator and output
- Delete nRNA fan-out in union-find connected components (lines ~2160-2185)
- Simplify EM layout: `[mRNA(0..n_t), gDNA]` with `n_t + 1` components
- **Net**: ~150 lines of C++ removed

**C.3 — scoring.py**: Remove nRNA parameters from Python scorer interface.
- Remove `nrna_base`, `t_to_nrna_arr`, `nrna_span_arr`, `nrna_start_arr` from `FragmentScorer` construction

**C.4 — locus.py**: Remove nRNA CSR construction.
- Remove nRNA→transcript CSR building in `build_loci()`
- Remove nRNA-specific prior logic in `build_locus_em_data()`
- Simplify `_compute_gamma()` — no nRNA totals needed

**C.5 — estimator.py**: Remove nRNA fan-out.
- Delete `_nrna_per_transcript()` method
- Delete `nrna_em_counts` array
- Synthetic transcripts have direct per-transcript counts in `em_counts`

**C.6 — scan.py**: Remove `nrna_base_index` computation.

**C.7 — pipeline.py**: Simplify — fewer nRNA-specific parameters.

**C.8 — native.py**: Update public interface — remove nRNA parameters.

**C.9 — Unspliced strand accumulators**: These currently use `nrna_idx` for indexing (`acc_us[nrna_idx]`, `acc_ua[nrna_idx]`). Must be reworked to index by `t_index` or removed. Options:
  - **Rework**: Change to `acc_us[t_idx]`, `acc_ua[t_idx]` — array sized to `n_transcripts`
  - **Remove**: They're currently not consumed by production code (only by tests for diagnostics). Could be removed as dead code. Decision: keep and rework for now (diagnostic value).

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
7. **Splice-annotated regression**: Verify that splice-annotated fragments are correctly assigned to multi-exon transcripts (not siphoned by synthetics)


## Files affected

| File | Role | Change |
|------|------|--------|
| `src/rigel/index.py` | nRNA table → synthetic transcripts | Integrate `create_nrna_transcripts()` into build (Phase B), remove `compute_nrna_table()` (Phase C) |
| `src/rigel/transcript.py` | Transcript dataclass | `is_synthetic_nrna`, `is_nascent_equiv` fields added (done) |
| `src/rigel/native/scoring.cpp` | Fragment scoring | Phase C: remove ~200 lines of nRNA scoring |
| `src/rigel/native/em_solver.cpp` | EM solver | Phase C: remove nRNA section, fan-out, CSR |
| `src/rigel/locus.py` | Locus construction | Phase C: remove nRNA CSR, simplify `build_loci` |
| `src/rigel/estimator.py` | Abundance estimation | Phase C: remove nRNA fan-out, simplify output |
| `src/rigel/scoring.py` | Python scorer interface | Phase C: remove nRNA parameters |
| `src/rigel/scan.py` | Scored fragments | Phase C: remove `nrna_base_index` |
| `src/rigel/pipeline.py` | Orchestrator | Phase C: simplify — fewer parameters to thread |
| `src/rigel/native.py` | C++ binding interface | Phase C: remove nRNA parameters |


## Changelog

- **2026-03-21**: Dead code cleanup — removed `intronic_sense/antisense` accumulators (scoring.cpp, estimator.py, scan.py), `ui_bp` consumption, `has_ui` gating logic, and `unambig_intron_bp` from scoring tuple. All 1,037 tests pass. Added UNAMBIG_INTRON investigation results section documenting the nrna_idx/t_index mismatch bug, the evidence-gating gap, and the design decision to maintain current behavior.
- **2026-03-20**: Phase A diagnostic completed on GENCODE index. 20bp tolerance recommended. Implementation infrastructure (`create_nrna_transcripts()`, `_cluster_coordinates()`, `Transcript` fields) written.
- **2026-03-19**: Initial plan — unified architecture proposal, tolerance-based merging algorithm.
