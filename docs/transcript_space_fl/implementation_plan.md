# Transcript-Space Fragment Length Computation

## Implementation Plan — Replace Gap SJ Correction

**Status:** Proposal — awaiting review before implementation

---

## 1. Evaluation: Current Approach vs. Transcript-Space Projection

### Current approach (gap SJ correction)

The current `compute_frag_lengths()` in `resolve_context.h` works in **genomic coordinate space**:

1. Sort alignment blocks (exon blocks from CIGAR) by genomic position
2. Compute genomic `footprint = last_block_end - first_block_start`
3. Subtract observed intron (CIGAR N-op) lengths → `upper`
4. Identify inter-block "gaps" not explained by observed introns
5. For each gap, query the `sj_cr_` cgranges index to find annotated introns (per-transcript) that overlap the gap
6. Per-transcript correction: `FL(t) = upper - gap_correction(t)`

**Weaknesses:**

- **Conceptual indirection.** We compute FL by *subtracting* what we don't want (introns in gaps), rather than *measuring* what we do want (spliced fragment length). This introduces edge cases whenever the subtraction doesn't perfectly reconstruct the true FL.
- **Partial overlap sensitivity.** When a read aligns a few bases into intronic space (common with minimap2), the "gap" boundaries don't exactly match the annotated intron boundaries. The overlap fix mitigates this but is still a heuristic.
- **Dedicated SJ gap cgranges index.** The `sj_cr_` index stores one entry per (intron, transcript) pair — a separate data structure that exists solely for this computation.
- **Complex code path.** The function has ~100 lines handling sorted blocks, observed intron sets, gap detection, cgranges queries, and per-transcript accumulation.

### Proposed approach (transcript-space projection)

At resolution time, we already know the candidate transcript set (`t_inds`) and have the fragment's aligned exon blocks (genomic coordinates). The proposal: for each candidate transcript, project the fragment's outermost genomic positions into transcript coordinate space and compute FL as the difference.

**Algorithm:**

```
for each candidate transcript t:
    tx_start = genomic_to_tx_pos(fragment_genomic_start, t)
    tx_end   = genomic_to_tx_pos(fragment_genomic_end,   t)
    FL(t)    = |tx_end - tx_start|
```

**Strengths:**

- **Conceptually direct.** FL is the distance between the fragment's endpoints in the coordinate system where FL is defined (spliced transcript space). No gap detection, no subtraction.
- **Inherently robust to overhangs.** Fragment endpoint overhang (into introns or beyond transcript boundaries) is handled by adding the overhang distance rather than snapping. Internal intronic overhang (at non-endpoint block boundaries) is irrelevant because only gstart/gend are projected.
- **Eliminates the SJ gap cgranges index.** The `sj_cr_` index, `sj_t_index_[]`, and `sj_strand_[]` become unnecessary for FL computation.
- **Simpler code.** The entire `compute_frag_lengths()` function reduces to: look up two transcript positions, subtract.
- **Already proven.** `genomic_to_tx_pos()` is already production code in `scoring.cpp` used for coverage weight computation. It handles strand flip and clamps correctly.

**Considerations:**

| Concern | Assessment |
|---------|-----------|
| **Performance** | Each call is O(log E) where E = number of exons per transcript (binary search). Current approach does O(G × R) where G = gaps and R = cgranges results per gap. For typical fragments (≤ 2 gaps, ~5 exons/transcript), the projection is faster. For mega-loci with 1000+ candidate transcripts, both approaches scale linearly with candidates. |
| **Data availability at resolution time** | `FragmentResolver` currently does NOT have per-transcript exon arrays. These must be added — same CSR structure already in `NativeFragmentScorer`. |
| **nRNA (nascent RNA) candidates** | nRNA spans are single synthetic exons covering the full genomic interval. Projection gives `FL = genomic_end - genomic_start` (the unspliced length), which is the correct FL for an nRNA interpretation. This is the same behavior as the current code. |
| **Strand handling** | `genomic_to_tx_pos()` already handles negative-strand transcripts via `t_len - offset`. However, for FL computation we want the *unsigned* distance in transcript space (strand-agnostic FL), so we should either skip the strand flip or take `abs(tx_end - tx_start)`. |
| **Multi-contig / chimeric fragments** | Current code already returns empty for chimeric fragments (`chimera_type != CHIMERA_NONE`). No change needed. |
| **Fragment endpoint overhang** | Fragment endpoints (gstart, gend) can fall outside exons — before the first exon, past the last exon, or in an internal intron. All three cases must add the overhang distance to the transcript-space coordinate (no clamping, no snapping). Internal block-level intronic overhang (not at gstart/gend) is safe because `compute_frag_lengths()` only projects the outermost endpoints. See Step 2 for detailed examples and rationale. |

### Verdict

**Yes — the transcript-space approach is clearly better.** It is more correct by construction, more robust to alignment artifacts, eliminates a dedicated data structure, and simplifies the code substantially. The only cost is adding the exon CSR arrays to `FragmentResolver`, which is straightforward.

---

## 2. Implementation Plan

### Step 1: Add exon CSR arrays to `FragmentResolver`

**File:** `src/rigel/native/resolve_context.h`

Add four new member vectors to `FragmentResolver`:

```cpp
// Per-transcript exon structure (CSR layout, same as NativeFragmentScorer)
std::vector<int32_t> exon_offsets_;    // [n_transcripts + 1]
std::vector<int32_t> exon_starts_;     // [total_exons] — genomic start coords
std::vector<int32_t> exon_ends_;       // [total_exons] — genomic end coords
std::vector<int32_t> exon_cumsum_;     // [total_exons] — cumulative spliced bp before each exon
std::vector<int32_t> t_length_;        // [n_transcripts] — spliced transcript lengths
```

Add a build method:

```cpp
void build_exon_index(
    const std::vector<int32_t>& offsets,
    const std::vector<int32_t>& starts,
    const std::vector<int32_t>& ends,
    const std::vector<int32_t>& cumsum,
    const std::vector<int32_t>& lengths);
```

### Step 2: Add `genomic_to_tx_pos()` to `FragmentResolver`

**File:** `src/rigel/native/resolve_context.h`

Port the existing `genomic_to_tx_pos()` from `scoring.cpp` with two critical changes:

1. **Omit the strand flip** — FL must be strand-agnostic (physical length of the molecule).
2. **Do NOT clamp to `[0, t_len]`** — fragments may overhang beyond transcript boundaries, and those bases must be counted in the FL.

**Fragment endpoint overhang:** `compute_frag_lengths()` calls `genomic_to_tx_pos()` exactly twice per candidate transcript — once for `gstart` (min of all block starts) and once for `gend` (max of all block ends). These are always the **outermost physical endpoints** of the fragment. When an endpoint falls outside an exon, those bases are physically part of the fragment and must be counted.

There are two types of overhang, with different desired behavior:

1. **Endpoint overhang** — `gstart` or `gend` falls outside an exon (before the first exon, past the last exon, or in an internal intron). These bases **must** be counted because they are the physical ends of the fragment.

2. **Internal overhang** — An alignment block extends a few bases into intronic space, but the overhang is at an inner block boundary (not at gstart or gend). This is **safe to ignore** because `compute_frag_lengths()` never queries internal block boundaries — only the outermost endpoints.

Since the function is only called for gstart/gend, the correct behavior is: **always add the overhang** when a position is outside an exon. There is no need to distinguish intron-between-exons from past-last-exon.

**Example A — Endpoint overhang into intron:**

- Transcript exons: `(1000, 2000)` and `(5000, 6000)`, transcript length = 2000
- Fragment blocks: `(1600, 1750)` and `(1855, 2005)` — gend extends 5 bp past exon 1 into the intron
- `genomic_to_tx_pos(1600) = 600` (inside exon 0)
- `genomic_to_tx_pos(2005) = 1000 + 5 = 1005` (end of exon 0 in tx space + 5 bp overhang)
- `FL = |1005 - 600| = 405` ✓ (matches genomic footprint: 2005 − 1600 = 405)

**Example B — Internal intronic overhang (naturally safe):**

- Same transcript exons: `(1000, 2000)` and `(5000, 6000)`
- Fragment blocks: `(1855, 2005)` and `(4999, 5149)` — both blocks overhang into the intron at their inner edges
- gstart = 1855 (inside exon 0), gend = 5149 (inside exon 1) — neither endpoint is in an intron
- `genomic_to_tx_pos(1855) = 855`, `genomic_to_tx_pos(5149) = 1149`
- `FL = |1149 - 855| = 294` ✓ — the internal overhangs at 2005 and 4999 are irrelevant because they are not the fragment's outer endpoints

**Example C — Endpoint overhang before first exon:**

- Transcript exons: `(3000, 3100)` and `(14400, 17000)`, transcript length = 2700
- Fragment blocks: `(2995, 3100)` and `(14400, 14445)` — gstart extends 5 bp before first exon
- `genomic_to_tx_pos(2995) = 2995 - 3000 = -5`
- `genomic_to_tx_pos(14445) = 100 + 45 = 145`
- `FL = |145 - (-5)| = 150` ✓

All three cases (before-first, intron-between, past-last) use the same logic: transcript-space position at nearest exon boundary **plus** the genomic distance into non-exonic space.

```cpp
/// Map a genomic position to transcript-space offset for FL computation.
///
/// Designed for projecting fragment ENDPOINTS (gstart, gend) only.
/// When a position falls outside an exon, the overhang distance is added
/// rather than snapped/clipped, because fragment endpoints represent the
/// physical extent of the molecule and those bases must be counted in FL.
///
/// - Strand-agnostic (no strand flip): always returns forward-strand
///   spliced coordinate.
/// - NOT clamped: returns negative values (before first exon), values
///   > t_len (past last exon), or tx-position + intronic-overhang
///   (endpoint in internal intron).
inline int32_t genomic_to_tx_pos(int32_t genomic_pos, int32_t t_idx) const {
    int32_t begin   = exon_offsets_[t_idx];
    int32_t end     = exon_offsets_[t_idx + 1];
    int32_t n_exons = end - begin;
    if (n_exons <= 0) return 0;

    const int32_t* starts = exon_starts_.data() + begin;
    const int32_t* ends   = exon_ends_.data()   + begin;
    const int32_t* cumsum = exon_cumsum_.data()  + begin;

    // bisect_right(starts, genomic_pos) - 1
    int ei = static_cast<int>(
        std::upper_bound(starts, starts + n_exons, genomic_pos) - starts) - 1;

    if (ei < 0) {
        // Before first exon: negative offset
        return genomic_pos - starts[0];
    }
    if (genomic_pos >= ends[ei]) {
        // Past end of exon[ei] — either in intron or past last exon.
        // Always add the overhang (no snapping): this is correct because
        // compute_frag_lengths() only calls this for fragment endpoints,
        // and endpoint bases must be counted in the FL.
        return cumsum[ei] + (ends[ei] - starts[ei]) + (genomic_pos - ends[ei]);
    }
    // Inside exon ei
    return cumsum[ei] + (genomic_pos - starts[ei]);
}
```

Note: The returned value is always in forward-strand spliced coordinates. Negative return values are valid (pre-transcript overhang). Values exceeding `t_len` are valid (post-transcript or intron-end overhang). Internal intronic overhangs at non-endpoint block boundaries are safe because they are never passed to this function.

### Step 3: Replace `compute_frag_lengths()`

**File:** `src/rigel/native/resolve_context.h`

Replace the current implementation with the transcript-space projection:

```cpp
std::unordered_map<int32_t, int32_t> compute_frag_lengths(
    const std::vector<ExonBlock>& exons,
    const std::vector<IntronBlock>& introns,   // unused, kept for signature compat
    const std::vector<int32_t>& t_inds,
    ResolverScratch& scratch) const             // scratch no longer needed for SJ buf
{
    std::unordered_map<int32_t, int32_t> result;
    if (exons.empty() || t_inds.empty()) return result;

    // Single exon block: FL is trivially block length (no introns possible)
    if (exons.size() == 1) {
        int32_t fl = exons[0].end - exons[0].start;
        if (fl > 0)
            for (int32_t t : t_inds) result[t] = fl;
        return result;
    }

    // Find outermost genomic positions across all alignment blocks
    int32_t gstart = exons[0].start;
    int32_t gend = exons[0].end;
    for (const auto& e : exons) {
        gstart = std::min(gstart, e.start);
        gend = std::max(gend, e.end);
    }

    for (int32_t t : t_inds) {
        int32_t tx_s = genomic_to_tx_pos(gstart, t);
        int32_t tx_e = genomic_to_tx_pos(gend, t);
        int32_t fl = std::abs(tx_e - tx_s);
        if (fl > 0) result[t] = fl;
    }
    return result;
}
```

### Step 4: Wire up exon index construction from Python

**File:** `src/rigel/native/resolve.cpp` (nanobind bindings)

Expose `build_exon_index()` to Python.

**File:** `src/rigel/pipeline.py` (or `index.py`)

Build the exon CSR arrays (same logic as `scoring.py` `_build_exon_data()`) and pass to the resolver:

```python
# Build exon CSR for resolver (same data as scorer uses)
exon_offsets, exon_starts, exon_ends, exon_cumsum, t_lengths = (
    build_exon_csr(index)
)
resolve_ctx.build_exon_index(
    exon_offsets, exon_starts, exon_ends, exon_cumsum, t_lengths
)
```

The `build_exon_csr()` helper should be factored out from the existing scorer construction code in `scoring.py` so both the resolver and scorer share the same function.

### Step 5: Remove SJ gap index (cleanup)

Once the new FL computation is verified:

- **Remove** `sj_cr_`, `sj_t_index_[]`, `sj_strand_[]` from `FragmentResolver`
- **Remove** `build_sj_gap_index()` method
- **Remove** calls to `build_sj_gap_index()` in Python (`index.py` / `pipeline.py`)
- **Remove** `sj_df` gap index construction from `index.py` if it is only used for this purpose (check if `sj_df` is used elsewhere first — it serves the SJ exact-match map too)
- **Remove** `scratch.sj_buf` / `scratch.sj_buf_cap` from `ResolverScratch` if no longer needed

### Step 6: Tests

1. **Existing tests should pass unchanged.** The new FL computation should produce identical results for fragments where the gap perfectly aligns with annotated introns (the majority case).
2. **Known edge cases to verify:**
   - Single-exon fragments (FL = block length) — unchanged
   - Spliced fragments with no gap — FL from observed intron subtraction; projection gives same answer
   - Gap fragments with exact intron match — projection gives same answer
   - Gap fragments with small overhang into intronic space — projection is *more correct* (snaps to exon boundary)
   - nRNA candidates — single giant exon, FL = genomic span
   - Chimeric fragments — skipped (no FL computed), unchanged
3. **New unit tests for `genomic_to_tx_pos()`:** Test the coordinate projection function directly with a known exon structure:
   - Position inside an exon → correct spliced offset
   - Position in an intron between exons → snaps to preceding exon end
   - Position before first exon (transcript boundary overhang) → negative value
   - Position past last exon (transcript boundary overhang) → value > t_len
   - Single-exon transcript (nRNA) — position before, inside, and after
4. **New unit tests for `compute_frag_lengths()` with overhang:**
   - **Start overhang (before first exon):** Fragment block starts 5 bp before the transcript's first exon. Verify FL counts the overhang (e.g., Example C from this plan: blocks at `(2995, 3100)` + `(14400, 14445)` with exons `(3000, 3100)` + `(14400, 17000)` → FL = 150).
   - **End overhang (past last exon):** Fragment block extends 5 bp past the transcript's last exon. Verify FL counts the overhang.
   - **Both-end overhang:** Fragment overhangs both the start and end of the transcript boundaries.
   - **Endpoint overhang into internal intron:** Fragment endpoint (`gend`) extends past an exon into the intron between two exons (Example A: blocks at `(1600, 1750)` + `(1855, 2005)` with exons `(1000, 2000)` + `(5000, 6000)` → FL must be 405, not 400). This is the critical case that distinguishes endpoint overhang from internal overhang.
   - **Internal intronic overhang (NOT at endpoints, safe):** Fragment blocks overhang into intronic space at inner block boundaries but gstart/gend are inside exons (Example B: blocks at `(1855, 2005)` + `(4999, 5149)` with same exons → FL = 294). Verify the internal overhangs don't inflate FL.
5. **Golden output:** Run `pytest tests/ --update-golden` if any golden outputs change due to improved FL accuracy.

---

## 3. Risk Assessment

| Risk | Mitigation |
|------|------------|
| FL values change for some fragments | Expected only for edge cases (overhang into introns). These changes are *improvements*. Run full test suite + benchmark to confirm. |
| Memory increase from exon CSR in resolver | The same data is already stored in the scorer. For a typical index (~200K transcripts, ~1.5M exons), the CSR is ~18 MB — negligible. |
| Regression in mega-loci performance | Unlikely — binary search is O(log E) vs. cgranges overlap O(log N + R). Profile if concerned. |
| SJ gap index removal breaks something | Check all uses of `sj_cr_` before removing. It is only used in `compute_frag_lengths()`. |

---

## 4. Implementation Order

1. Steps 1–3: Add exon CSR + `genomic_to_tx_pos()` + new `compute_frag_lengths()` (C++)
2. Step 4: Wire up Python-side exon index construction
3. Compile + run tests
4. Run benchmarks to compare FL distributions
5. Step 5: Remove SJ gap index (cleanup, separate commit)
6. Step 6: Add edge-case tests, update golden outputs

---

## 5. Why This Is Better

The fundamental insight: **fragment length is naturally defined in transcript coordinate space, not genomic space.** The current approach works backward — starting with the genomic footprint and trying to subtract introns it can't see. The proposed approach works forward — projecting endpoints into the space where FL is meaningful and measuring directly.

The `genomic_to_tx_pos()` function is already battle-tested in the scoring hot path. Reusing it in the resolver makes FL computation a first-class coordinate transformation rather than a gap-detection heuristic.
