# Fragment-to-Region Counting: Revised Implementation Plan (v2)

## 1. Overview

This document covers the design and implementation of **fragment-to-region
mapping and per-region fractional fragment counting** — the second step in the
revised gDNA calibration pipeline. The goal is to tabulate per-region fragment
evidence (strand-stratified, splice-stratified fractional counts) that will
feed downstream gDNA model training (strand symmetry κ, fragment-length
distribution).

**Prerequisite:** The genome region partition (Workstream A) is already
implemented — `build_region_table()` produces a non-overlapping, reference-
complete partition stored as `regions.feather` and loaded as `region_df` /
`region_cr` (cgranges index) on `TranscriptIndex`.

---

## 2. Critique of v1 Plan and Design Corrections

The v1 plan had five fundamental flaws that are corrected here.

### 2.1 Fractional Counting Is Required (Not Integer)

**v1 decision (wrong):** Integer counting, assign each fragment to a single
region (midpoint or max-overlap).

**Correction:** Fractional counting is mandatory. A fragment spanning two
regions at a boundary contributes proportional fractions to each. This is
mathematically correct and essential for unbiased per-region density estimation.
Integer winner-takes-all assignment creates systematic bias at region
boundaries, especially for the many small exonic regions (median 163 bp) that
are comparable to fragment size (~200–500 bp).

**Denominator:** The fragment's **total aligned length** (sum of exon block
sizes), NOT `genomic_footprint` (which includes intronic gaps for spliced
fragments).

### 2.2 Midpoint Assignment for Spliced Fragments Is Wrong

**v1 decision (wrong):** Use midpoint of genomic footprint for spliced
fragments, since aligned blocks are unavailable in the buffer.

**Correction:** The midpoint of the genomic footprint for a spliced fragment
falls in the **intron**, not on any aligned base. For a fragment with exon
blocks [100–200] and [5000–5100], the midpoint is ~2550 — deep in the intron.
Assigning this fragment to whatever region contains position 2550 is
fundamentally wrong. Each aligned block must be individually overlapped against
the region partition, and fractional counts accumulated per region from the
actual aligned bases.

### 2.3 Intergenic Uplift Is Unacceptable

**v1 decision (wrong):** For intergenic fragments (missing from the buffer),
distribute global intergenic counts proportionally across intergenic regions
by region length.

**Correction:** Intergenic fragments are the purest gDNA signal and the most
important calibration targets. Uniform-distribution uplift destroys per-region
granularity — it assumes gDNA is uniformly distributed, which is the very
hypothesis we're trying to test with per-region evidence. Every individual
intergenic fragment must be counted in its actual region(s) using the same
fractional overlap algorithm as resolved fragments.

### 2.4 The Buffer Is Too Late — Aligned Blocks Are Discarded

**v1 decision (wrong):** Operate on the `FragmentBuffer` after
`scan_and_buffer()`.

**Correction:** Region counting must happen at the `AssembledFragment` stage
in the C++ BAM scanner — AFTER `build_fragment()` constructs the fragment with
full aligned block geometry (`vector<ExonBlock> exons`), but BEFORE
`_resolve_core()` discards that geometry. This is the only point in the
pipeline where ALL fragments (including intergenic) have their complete
aligned block coordinates available.

The critical code path in `bam_scanner.cpp :: process_qname_group_threaded()`:

```cpp
AssembledFragment frag = build_fragment(r1_reads, r2_reads);
// ← HERE: frag.exons has full {ref_id, start, end, strand} per block
// ← Region counting must happen HERE, for ALL fragments

RawResolveResult cr;
bool resolved = ctx._resolve_core(frag.exons, frag.introns, ...);
if (!resolved) {
    // Intergenic → stats only, then continue (lost forever)
    continue;
}
// resolved → accumulator.append() (blocks already discarded)
```

### 2.5 Staged Development Is Required

**v1 decision (missing):** No development staging; single monolithic plan.

**Correction:** Two-stage roadmap:

- **Stage 1 (Prototype):** Pure Python proof-of-concept using pysam to read
  BAM records, group by query name, extract CIGAR blocks, and count against
  the region partition. This enables rapid iteration, testing, and validation
  without C++ compilation cycles.

- **Stage 2 (Production):** C++ integration into the BAM scanner hot path,
  adding a `RegionAccumulator` to `WorkerState` that counts fragments at the
  `AssembledFragment` stage. Single-pass, multi-threaded, zero extra I/O.

---

## 3. Architecture: Where Region Counting Fits

### 3.1 Fragment Lifecycle in the C++ Scanner

```
BAM records
    ↓  parse_bam_record()
ParsedAlignment {ref_id, exons[(s,e)], sjs[(s,e,strand)], is_reverse, nm}
    ↓  group_records_by_hit()
QnameGroup → AlignmentGroup {hits: [(r1_reads, r2_reads), ...]}
    ↓  build_fragment()                              ┐
AssembledFragment {                                   │ ALIGNED BLOCKS
    vector<ExonBlock> exons   ← {ref_id,start,end,strand} │ ARE AVAILABLE
    vector<IntronBlock> introns                       │ HERE
    int32_t nm                                        │
}                                                     │
    ↓  ← ── REGION COUNTING INSERTION POINT ── ── ── ┘
    ↓  _resolve_core()
RawResolveResult {t_inds, frag_lengths, exon_bp, ...}  ← blocks GONE
    ↓  if (!resolved) → continue  ← INTERGENIC LOST
    ↓  ResolvedFragment::from_core()
ResolvedFragment {genomic_start, genomic_footprint, ...} ← blocks GONE
    ↓  accumulator.append()
FragmentAccumulator → FragmentBuffer → Python
```

Key observations:
- `AssembledFragment.exons` is a `vector<ExonBlock>` where each block has
  `{ref_id, start, end, strand}` — exactly what we need for region overlap
- After `_resolve_core()`, only `genomic_start` and `genomic_footprint` survive
- Intergenic fragments (unresolved) are counted in global stats but their
  coordinates are never stored

### 3.2 `ExonBlock` Structure

```cpp
struct ExonBlock {
    int32_t ref_id;   // htslib tid mapped to internal ref numbering
    int32_t start;    // 0-based, inclusive
    int32_t end;      // 0-based, exclusive
    int32_t strand;   // STRAND_POS=1, STRAND_NEG=2
};
```

**Strand semantics:** `ExonBlock.strand` stores the **genomic alignment strand**
(not transcript strand). This is consistent with the region partition, which
uses genomic strand. For PE fragments, `build_fragment()` applies the Illumina
R2 strand flip — since Illumina PE sequencing produces converging reads on
opposite strands, R2's strand is flipped (POS↔NEG) so both reads express the
same genomic alignment strand of the original fragment. After the flip, R1 and
R2 exon blocks are merged into a sorted, non-overlapping set grouped by
`(ref_id, strand)`. The strand flip is the standard convention and is required
for consistent strand assignment.

### 3.3 Multimapping Fragments

**Verified in existing code:** When a read name has NH > 1 or secondary
alignments, `group_records_by_hit()` separates the alignments into distinct
hits (by HI tag, or via secondary pairing heuristics). Each hit produces a
separate `AssembledFragment` via `build_fragment()`. All AssembledFragments
from the same read name share the same `frag_id` (assigned once per
`QnameGroup` in the BAM reader loop). The variable `num_hits = max(nh,
all_hits.size())` tracks the total hit count, and `is_unique_mapper =
(num_hits == 1)`.

**For region counting:** Multi-mapping fragments are **excluded** in both
stages. A fragment group with `num_hits > 1` skips region counting entirely.
At the `AssembledFragment` level, this is straightforward — the enclosing
loop already knows `is_unique_mapper` before calling `build_fragment()`.

Note: if `include_multimap` is false in the scanner config, multi-mapping
groups return early before any AssembledFragment is built (line ~1150:
`if (is_multimap) { stats.multimapping++; if (!include_multimap) return; }`).
When `include_multimap` is true, each hit creates its own AssembledFragment
but the region accumulator skips them via the `is_unique_mapper` check.

### 3.4 `AssembledFragment` Properties

```cpp
struct AssembledFragment {
    std::vector<ExonBlock> exons;
    std::vector<IntronBlock> introns;
    int32_t nm;
    int32_t genomic_footprint() const;  // max(end) - min(start) across exons
    bool has_introns() const;           // !introns.empty()
};
```

The total aligned length is `sum(eb.end - eb.start for eb in exons)`. This is
the correct denominator for fractional counting.

### 3.5 Threading Model

The scanner uses a producer–consumer architecture:

- **Main thread:** Reads BAM, groups records by qname, pushes `QnameGroup`
  objects into a `BoundedQueue`
- **Worker threads:** Pop groups from queue, call `build_fragment()` +
  `_resolve_core()`, accumulate results into per-thread `WorkerState`

Each `WorkerState` has its own `FragmentAccumulator`, `BamScanStats`,
`StrandObservations`, and `FragLenObservations` — all merged after scan
completes. A new `RegionAccumulator` would follow the same pattern: per-worker
accumulation, post-scan merge.

---

## 4. Fractional Counting Algorithm

### 4.1 Core Algorithm

For each fragment with aligned blocks `B = [b_1, ..., b_k]`:

1. Compute **fragment aligned length**: $L = \sum_{i=1}^{k} (b_i.\text{end} - b_i.\text{start})$

2. For each block $b_i$, query the region partition for overlapping regions:
   `region_cr.overlap(ref_name(b_i.ref_id), b_i.start, b_i.end)`

3. For each overlapping region $r$, compute base-pair overlap:
   $\text{overlap}(b_i, r) = \min(b_i.\text{end}, r.\text{end}) - \max(b_i.\text{start}, r.\text{start})$

4. Accumulate fractional counts per region:
   $\text{frac}(r) = \frac{1}{L} \sum_{i : b_i \cap r \neq \emptyset} \text{overlap}(b_i, r)$

5. Add `frac(r)` to the appropriate strand × splice cell for region $r$.

**Invariant:** $\sum_r \text{frac}(r) = 1.0$ for every counted fragment (since
the region partition is reference-complete and non-overlapping, the sum of all
overlaps equals $L$).

### 4.2 Strand Determination

Strand is derived from the fragment's exon blocks, consistent with the existing
scanner logic:

```python
strand = STRAND_NONE
for block in fragment.exons:
    strand |= block.strand

# strand is now STRAND_POS, STRAND_NEG, or STRAND_AMBIGUOUS (POS|NEG)
```

- `STRAND_POS (1)`: All blocks map to positive genomic strand
- `STRAND_NEG (2)`: All blocks map to negative genomic strand
- `STRAND_AMBIGUOUS (3)`: Blocks span both strands (chimeric indicator)

Fragments with `STRAND_AMBIGUOUS` are excluded from region counting (same
treatment as chimeric fragments in the existing scanner).

**Important caveat:** This strand OR-ing procedure is only valid when
multi-mapping fragments have already been excluded. A multi-mapper with hits
on different strands would OR to `STRAND_AMBIGUOUS` even though each individual
hit has a definite strand. Since we exclude multi-mappers before region
counting (Section 4.4), the OR-ing correctly reflects the strand of each
unique-mapper's single hit.

### 4.3 Splice Determination

A fragment is considered **spliced** if it has intron blocks
(`frag.has_introns()` returns true). Splice annotation status:

- `SPLICE_UNSPLICED (0)`: No CIGAR N operations
- `SPLICE_SPLICED_UNANNOT (1)`: CIGAR N present, SJs don't match annotation
- `SPLICE_SPLICED_ANNOT (2)`: CIGAR N present, ≥1 SJ matches annotation

For the prototype stage (before transcript resolution), we can determine
splice status from `frag.has_introns()`. For annotated vs. unannotated, we
need the splice junction annotation check from `_resolve_core()`. In the
prototype, we classify as spliced/unspliced only; the full three-way
classification is available at the production stage where region counting
happens alongside resolution.

**Future improvement area:** The current binary spliced/unspliced classification
groups annotated and unannotated spliced reads together. For downstream gDNA
calibration, annotated spliced reads are strong RNA evidence while unannotated
spliced reads are weaker (they could be alignment artifacts or novel junctions).
The production stage should provide the three-way split
(`unspliced / spliced_annot / spliced_unannot`) as separate count columns,
enabling the calibration model to weight them differently. This is a
straightforward extension once `splice_type` is available from `_resolve_core()`
at the counting point.

### 4.4 Fragment Admissibility

| Fragment type | Counted? | Rationale |
|---------------|----------|-----------|
| Unique mapper, unspliced | Yes | Primary gDNA calibration signal |
| Unique mapper, spliced | Yes | RNA evidence (anti-gDNA signal) |
| Intergenic (unresolved) | Yes | Purest gDNA signal — critical |
| Multi-mapper (NH > 1) | No | Ambiguous genomic locus; excluded (Section 3.3) |
| Chimeric (multi-ref/strand) | No | Unreliable coordinates |
| Empty exon blocks | No | Guard against degenerate alignments (see below) |

**Empty exon blocks note:** The `if (frag.exons.empty()) continue` guard in the
scanner handles a rare edge case — it can occur when `build_fragment()` receives
reads whose CIGAR operations produce no aligned bases (e.g., all soft/hard
clips), or when secondary pairing produces empty read vectors. This should not
occur in well-formed BAM data but the guard is warranted as defensive coding
against malformed input.

### 4.5 Worked Example

Fragment with two exon blocks on chr1, positive strand:
- Block A: [100, 250) → 150 bp
- Block B: [500, 600) → 100 bp
- Total aligned length $L$ = 250 bp

Region partition on chr1:
- Region 0: [0, 200) — intergenic
- Region 1: [200, 500) — intronic
- Region 2: [500, 700) — exonic

Overlap computation:
- Block A ∩ Region 0: min(250, 200) - max(100, 0) = 100 bp
- Block A ∩ Region 1: min(250, 500) - max(200, 200) = 50 bp
- Block B ∩ Region 2: min(600, 700) - max(500, 500) = 100 bp

Fractional counts (unspliced_pos):
- Region 0: 100 / 250 = 0.40
- Region 1: 50 / 250 = 0.20
- Region 2: 100 / 250 = 0.40

Sum = 1.0 ✓

---

## 5. Per-Region Evidence Schema

The output is a DataFrame (and Feather file) with one row per region:

| Column | Type | Description |
|--------|------|-------------|
| `region_id` | int32 | Foreign key to `region_df` |
| `n_unspliced_pos` | float32 | Fractional unspliced count, POS strand |
| `n_unspliced_neg` | float32 | Fractional unspliced count, NEG strand |
| `n_spliced_pos` | float32 | Fractional spliced count, POS strand |
| `n_spliced_neg` | float32 | Fractional spliced count, NEG strand |

**Note:** Counts are `float32` (not int32) because of fractional counting.
A fragment spanning two regions contributes e.g. 0.6 to one and 0.4 to the
other — both columns must accept non-integer values. We use `float32` for
memory efficiency: for ~684K regions × 4 columns this saves 11 MB vs float64.
float32 provides ~7 decimal digits of precision, which is more than sufficient
for fragment counts (even a region with 1M fragments has only 7 significant
digits). Internal accumulation can use float64 to avoid rounding during
summation, then downcast to float32 for storage.

**Derived columns** (computed downstream, not stored):

- `n_unspliced = n_unspliced_pos + n_unspliced_neg`
- `n_spliced = n_spliced_pos + n_spliced_neg`
- `n_total = n_unspliced + n_spliced`
- `density_unspliced = n_unspliced / region_length`
- `strand_ratio_unspliced = n_unspliced_pos / n_unspliced` (for symmetry test)

**Design choice: 4 count columns, not 6.** The v1 plan had separate
`spliced_annot_pos/neg` columns. This is deferred; the prototype only needs
spliced vs. unspliced for gDNA calibration. Annotated splice distinction can be
added later since the production stage has access to splice type from
`_resolve_core()`.

---

## 6. Staged Implementation Roadmap

### Stage 1: Python Prototype

**Goal:** Functional proof-of-concept for development, testing, and validation.
Runs as a standalone pass over the BAM using pysam. Validates the fractional
counting algorithm and output schema before committing to C++ changes.

**Module:** `src/rigel/region_evidence.py`

**Input:** BAM file path + `TranscriptIndex` (with `region_cr`, `region_df`)

**Approach:**
1. Open BAM with pysam, iterate name-sorted records
2. Group records by query name
3. For each group, assemble fragment:
   - Extract CIGAR-derived aligned blocks `(ref_name, start, end)` per read
   - Apply R2 strand flip
   - Merge overlapping/adjacent blocks within `(ref, strand)` groups
   - Determine splice status from presence of N (skip) CIGAR ops
4. Apply admissibility filter (skip NH > 1, skip chimeric)
5. Run fractional counting algorithm (Section 4.1) against `region_cr`
6. Accumulate into a `(n_regions, 4)` float64 array (downcast to float32 for output)
7. Return DataFrame + write Feather output

**Why pysam (not the buffer):**
- The buffer lacks aligned block geometry — only `genomic_start` and
  `genomic_footprint` survive after C++ resolution
- The buffer excludes intergenic fragments entirely
- pysam provides full CIGAR access for every record in the BAM
- Development iteration is fast (no C++ recompilation)

**Performance expectation:** ~10–30 minutes for a 50M-fragment human BAM
(Python per-fragment loop with pysam). Acceptable for development and
validation; not for production.

**Fragment assembly in Python (pseudocode):**

```python
import pysam
from collections import defaultdict

def assemble_fragment(records):
    """Build aligned blocks from a group of BAM records with same qname."""
    # Key: (ref_name, strand) → list of (start, end) intervals
    exon_dict = defaultdict(list)
    intron_set = set()

    for rec in records:
        if rec.is_unmapped or rec.is_qcfail or rec.is_duplicate:
            continue

        ref = rec.reference_name
        # Strand: R1 uses is_reverse directly; R2 gets flipped
        if rec.is_read2:
            strand = STRAND_POS if rec.is_reverse else STRAND_NEG
        else:
            strand = STRAND_NEG if rec.is_reverse else STRAND_POS

        # Extract aligned blocks from CIGAR
        pos = rec.reference_start
        for op, length in rec.cigartuples:
            if op in (0, 7, 8):  # M, =, X → aligned
                exon_dict[(ref, strand)].append((pos, pos + length))
                pos += length
            elif op == 2:  # D → deletion (consumes ref)
                pos += length
            elif op == 3:  # N → skip (intron)
                intron_set.add((ref, pos, pos + length, strand))
                pos += length
            elif op == 1:  # I → insertion (consumes query only)
                pass
            elif op == 4:  # S → soft clip (consumes query only)
                pass
            elif op == 5:  # H → hard clip (consumes neither)
                pass

    # Merge overlapping/adjacent intervals per (ref, strand) group
    blocks = []
    for (ref, strand), intervals in exon_dict.items():
        intervals.sort()
        merged = [intervals[0]]
        for s, e in intervals[1:]:
            if s <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else:
                merged.append((s, e))
        for s, e in merged:
            blocks.append((ref, s, e, strand))

    is_spliced = len(intron_set) > 0
    return blocks, is_spliced
```

**Fractional counting in Python (pseudocode):**

```python
def count_fragment(blocks, is_spliced, region_cr, counts):
    """Accumulate fractional region counts for one fragment."""
    if not blocks:
        return

    # Determine fragment strand (all blocks should agree for non-chimeric)
    strands = set(strand for _, _, _, strand in blocks)
    if len(strands) != 1:
        return  # chimeric — skip
    frag_strand = strands.pop()
    if frag_strand not in (STRAND_POS, STRAND_NEG):
        return

    # Compute total aligned length (denominator)
    total_aligned = sum(end - start for _, start, end, _ in blocks)
    if total_aligned <= 0:
        return

    # Determine count column offset: unspliced=0, spliced=2; +1 for NEG
    col_base = 2 if is_spliced else 0
    col_offset = 0 if frag_strand == STRAND_POS else 1
    col = col_base + col_offset

    # Accumulate fractional overlap per region
    for ref, start, end, _ in blocks:
        for rgn_start, rgn_end, region_id in region_cr.overlap(ref, start, end):
            overlap_bp = min(end, rgn_end) - max(start, rgn_start)
            counts[region_id, col] += overlap_bp / total_aligned
```

**Main loop:**

```python
def count_region_evidence(bam_path, index):
    """Standalone BAM → region evidence using pysam."""
    n_regions = len(index.region_df)
    counts = np.zeros((n_regions, 4), dtype=np.float64)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        current_qname = None
        current_group = []

        for rec in bam:
            if rec.query_name != current_qname:
                if current_group:
                    _process_group(current_group, index.region_cr, counts)
                current_qname = rec.query_name
                current_group = [rec]
            else:
                current_group.append(rec)

        # Process last group
        if current_group:
            _process_group(current_group, index.region_cr, counts)

    return _build_dataframe(counts, index.region_df)
```

### Stage 2: C++ Production Integration

**Goal:** Single-pass, multi-threaded region counting inside the existing BAM
scanner. Zero extra BAM I/O. Counts ALL fragments (including intergenic) at
the `AssembledFragment` stage where aligned block geometry is available.

**Insertion point:** Inside `process_qname_group_threaded()` in
`bam_scanner.cpp`, immediately after `build_fragment()` and before
`_resolve_core()`:

```cpp
// Current code (bam_scanner.cpp, ~line 1196):
AssembledFragment frag = build_fragment(r1_reads, r2_reads);
stats.n_fragments++;
if (frag.exons.empty()) continue;

// ──── NEW: Region counting + FL tabulation ────
if (is_unique_mapper && !is_chimeric(frag)) {
    region_accum.count_fragment(frag, region_index, is_unique_mapper);
    // count_fragment() internally:
    //   1. Computes fractional region overlap for all blocks
    //   2. If unspliced and fully contained in one region,
    //      emits (region_id, genomic_footprint) FL observation
}
// ──── End new code ────

RawResolveResult cr;
bool resolved = ctx._resolve_core(frag.exons, frag.introns, ...);
```

**New C++ structures:**

```cpp
// Per-region counts: accumulated in double, exported as float32
// Columns: [unspliced_pos, unspliced_neg, spliced_pos, spliced_neg]
struct RegionAccumulator {
    int32_t n_regions;
    std::vector<double> counts;  // n_regions × 4, row-major (float64 for accumulation)

    // Fragment-length observations: (region_id, frag_len) pairs
    // Only for unspliced unique-mapper fragments fully contained in one region
    std::vector<int32_t> fl_region_ids;
    std::vector<int32_t> fl_frag_lens;

    explicit RegionAccumulator(int32_t n_regions)
        : n_regions(n_regions), counts(n_regions * 4, 0.0) {}

    void count_fragment(const AssembledFragment& frag,
                        const RegionIndex& region_index,
                        bool is_unique_mapper);

    void merge_into(RegionAccumulator& dst) const {
        for (size_t i = 0; i < counts.size(); i++) {
            dst.counts[i] += counts[i];
        }
        dst.fl_region_ids.insert(dst.fl_region_ids.end(),
                                  fl_region_ids.begin(), fl_region_ids.end());
        dst.fl_frag_lens.insert(dst.fl_frag_lens.end(),
                                 fl_frag_lens.begin(), fl_frag_lens.end());
    }
};
```

**`RegionIndex` wrapper:** A C++ wrapper around cgranges that accepts integer
`ref_id` (not string ref names) for fast lookup. Built from the Python
`region_cr` and passed to the scanner before scan starts.

```cpp
struct RegionIndex {
    // cgranges index per ref_id (or single flattened index)
    cr_intv_t* cr;
    // ref_id → cgranges contig name mapping
    std::vector<std::string> ref_names;

    // Query overlap for a single block
    std::vector<std::tuple<int32_t, int32_t, int32_t>>
    overlap(int32_t ref_id, int32_t start, int32_t end) const;
};
```

**Worker state extension:**

```cpp
struct WorkerState {
    ResolverScratch scratch;
    FragmentAccumulator accumulator;
    BamScanStats stats;
    StrandObservations strand_obs;
    FragLenObservations fraglen_obs;
    RegionAccumulator region_accum;  // ← NEW

    explicit WorkerState(int32_t n_transcripts, int32_t n_regions)
        : scratch(n_transcripts), region_accum(n_regions) {}
};
```

**Post-scan merge:** After all worker threads complete, merge per-worker
`RegionAccumulator` objects (same pattern as stats/strand_obs/fraglen_obs):

```cpp
// In BamScanner::scan(), after worker threads join:
RegionAccumulator merged_regions(n_regions);
for (auto& ws : worker_states) {
    ws->region_accum.merge_into(merged_regions);
}
// Export counts to Python as numpy array (downcast to float32)
// Export FL observations as two int32 numpy arrays (region_ids, frag_lens)
```

**Chimera detection at `AssembledFragment` stage:** A fragment is chimeric if
its exon blocks span multiple `(ref_id, strand)` combinations. This is
detectable directly from `frag.exons`:

```cpp
bool is_chimeric(const AssembledFragment& frag) {
    if (frag.exons.size() <= 1) return false;
    int32_t ref0 = frag.exons[0].ref_id;
    int32_t str0 = frag.exons[0].strand;
    for (size_t i = 1; i < frag.exons.size(); i++) {
        if (frag.exons[i].ref_id != ref0 ||
            frag.exons[i].strand != str0)
            return true;
    }
    return false;
}
```

### Stage Transition Criteria

Move from Stage 1 → Stage 2 when:

1. The Python prototype produces validated region counts and FL observations
   on test BAMs and human data
2. The per-region evidence schema is finalized (no further column changes)
3. The fractional counting algorithm has been verified against manual examples
4. FL containment criterion is validated (single-region fragments sufficient)
5. Performance profiling confirms the Python prototype is too slow for routine
   use (expected for production BAMs with >20M fragments)

At transition, the Python prototype becomes the **test oracle** — the C++
implementation must produce numerically identical region counts (within float32
tolerance) and identical FL observation tables on the same BAM.

---

## 7. Pipeline Integration

### 7.1 Stage 1 (Prototype) Integration

Region counting runs as a separate step, parallel to (not inside) the main
BAM scan:

```
rigel quant --bam sample.bam --index index/ -o results/
    ↓
  ┌─── scan_and_buffer(bam, index) ──→ FragmentBuffer ──→ quant_from_buffer()
  │
  └─── count_region_evidence(bam, index) ──→ region_evidence.feather
```

Both paths read the BAM independently. The prototype's pysam pass is a
separate BAM traversal. This is intentionally decoupled to avoid modifying the
existing scan-and-buffer pipeline during prototyping.

In `pipeline.py`:

```python
from rigel.region_evidence import count_region_evidence

# After scan_and_buffer completes:
region_evidence, fl_table = count_region_evidence(bam_path, index)
region_evidence.to_feather(output_dir / "region_evidence.feather")
fl_table.to_feather(output_dir / "region_frag_lengths.feather")
```

### 7.2 Stage 2 (Production) Integration

Region counting happens inside `scan_and_buffer()` during the single-pass BAM
scan. The `RegionAccumulator` is part of `WorkerState` and merged alongside
other accumulators:

```
rigel quant --bam sample.bam --index index/ -o results/
    ↓
  scan_and_buffer(bam, index)
    ├─ FragmentAccumulator  → FragmentBuffer → quant_from_buffer()
    └─ RegionAccumulator    → region_evidence.feather
```

The `BamScanner.scan()` returns the merged region counts as an additional
output array. Python side converts to DataFrame and writes Feather.

### 7.3 Synergy with Transcript Resolution

The region partition and transcript resolution indices are complementary
genomic interval structures that share the same reference coordinate system.
At Stage 2, both queries execute on the same `AssembledFragment` in the same
C++ worker thread. This opens opportunities:

- **Shared cgranges infrastructure:** The region cgranges and transcript
  cgranges use the same underlying library. A unified query interface could
  batch both lookups per exon block, amortizing memory access costs.
- **Region-aware resolution hints:** If a fragment falls entirely within an
  intergenic region (region count query returns only intergenic-flagged
  regions), transcript resolution can be skipped entirely — saving the cost
  of `_resolve_core()` for the ~40–60% of fragments that are intergenic.
- **Shared block iteration:** Both region counting and transcript resolution
  iterate over the same `frag.exons` vector. A combined pass could compute
  both results in a single exon-block loop.

These optimizations are deferred but should be explored once both systems are
production-stable.

---

## 8. Test Plan

### 8.1 Unit Tests (`tests/test_region_evidence.py`)

**Fragment assembly tests:**

| Test | Description |
|------|-------------|
| `test_assemble_single_read_unspliced` | Single alignment → one block |
| `test_assemble_paired_end_merge` | R1+R2 overlapping → merged block |
| `test_assemble_spliced_fragment` | CIGAR with N → multiple blocks + introns |
| `test_assemble_r2_strand_flip` | R2 strand is flipped (POS↔NEG) |
| `test_assemble_chimeric_detection` | Multi-ref blocks → chimeric |

**Fractional counting tests:**

| Test | Description |
|------|-------------|
| `test_single_block_single_region` | Fragment fully within one region → frac = 1.0 |
| `test_single_block_two_regions` | 200bp fragment spans boundary → correct fractions |
| `test_spliced_two_blocks_three_regions` | Spliced fragment, blocks in different regions |
| `test_fractions_sum_to_one` | Verify Σ frac(r) = 1.0 for every fragment |
| `test_intron_not_counted` | Intronic gap between blocks contributes zero overlap |

**Strand and splice stratification:**

| Test | Description |
|------|-------------|
| `test_pos_neg_separation` | Same-region POS+NEG fragments → separate columns |
| `test_spliced_unspliced_separation` | Same-region spliced+unspliced → separate columns |
| `test_ambiguous_strand_excluded` | STRAND_AMBIGUOUS → not counted |
| `test_chimeric_excluded` | Multi-ref fragment → not counted |
| `test_multimapper_excluded` | NH > 1 → not counted |

**Edge cases:**

| Test | Description |
|------|-------------|
| `test_empty_bam` | No fragments → all-zero counts |
| `test_all_regions_present` | Output has one row per region, even if count = 0 |
| `test_multi_chromosome` | Fragments on different references → correct regions |
| `test_intergenic_fragment_counted` | Fragment with no transcript overlap → still counted |

### 8.2 Integration Tests

| Test | Description |
|------|-------------|
| `test_pipeline_produces_evidence` | Full pipeline produces `region_evidence.feather` |
| `test_pipeline_produces_fl_table` | Full pipeline produces `region_frag_lengths.feather` |
| `test_evidence_golden_output` | Golden comparison for test BAMs |
| `test_fl_table_golden_output` | Golden comparison for FL observation table |
| `test_prototype_matches_production` | Stage 1 and Stage 2 produce identical output |

### 8.3 Human Genome Validation

After implementation, run on human dataset and verify:

- Total fractional counts sum to `n_counted_fragments` (within float tolerance)
- Intergenic regions have low density
- Exonic regions have high density on the annotated strand
- Strand ratio in intergenic regions is near 0.5 (gDNA symmetry expected)
- Per-region densities are plausible (no outlier regions with implausible counts)
- FL observation table has expected row count (~80% of unspliced unique mappers)
- Intergenic FL distribution matches global intergenic FL from scanner stats
- Exonic FL distribution is consistent with library prep size-selection curve

---

## 9. Performance Estimates

### Stage 1 (Python Prototype)

| Component | Estimate |
|-----------|----------|
| pysam BAM iteration | ~15 ns/record |
| Fragment assembly (Python) | ~2 μs/fragment |
| cgranges overlap query | ~200 ns/query (C library) |
| Python per-block loop | ~500 ns/block |
| **Total (50M fragments, ~1.5 blocks avg)** | **~15–30 minutes** |

Acceptable for development. Not for routine production use.

### Stage 2 (C++ Production)

| Component | Estimate |
|-----------|----------|
| Region counting per fragment | ~300 ns (1–2 cgranges queries) |
| Overhead vs. existing scan | ~5–10% additional |
| **Total (50M fragments)** | **~15 seconds (additive to existing scan)** |

Negligible compared to BAM I/O and transcript resolution.

### Memory

| Component | Human genome |
|-----------|-------------|
| Region counts array | 684K × 4 × 4 bytes = ~11 MB (float32 output) |
| Per-worker accumulation (float64, 4 workers) | ~88 MB total |
| cgranges region index | ~30 MB (already loaded) |
| FL observation table (~24M rows × 8 bytes) | ~192 MB |
| Per-worker FL vectors (4 workers) | ~768 MB peak (before merge) |

Well within budget.

---

## 10. Fragment-Length Tabulation per Region

Fragment-length distributions per region are **foundational to gDNA
classification**. The gDNA model needs to learn region-specific fragment-length
characteristics — gDNA fragments tend to follow a uniform-like distribution
while mRNA fragments follow the library prep size-selection distribution. This
section designs fragment-length tabulation as an integral part of the region
evidence framework.

### 10.1 Design Rationale

The region partition provides a natural framework for fragment-length analysis:
each region has a defined annotation type (exonic, intronic, intergenic) and a
strand. By collecting fragment-length observations per region, we can:

1. Build per-region fragment-length distributions for gDNA/RNA mixing models
2. Identify regions where the fragment-length distribution is gDNA-like
3. Use intergenic regions as a pure gDNA reference for the global FL model
4. Detect library prep artifacts that affect specific genomic contexts

### 10.2 Fragment-Region-Length Observation Table

**Approach: Store a per-fragment observation table rather than per-region
histograms.**

The output is a table of individual fragment-to-region observations for
unspliced, uniquely-mapped fragments:

| Column | Type | Description |
|--------|------|-------------|
| `region_id` | int32 | Region the fragment maps to |
| `frag_len` | int32 | Fragment length (`genomic_footprint` for unspliced) |

**Fragment eligibility for FL observations:**

| Criterion | Requirement | Rationale |
|-----------|-------------|-----------|
| Splice status | Unspliced only | Spliced fragments have ambiguous genomic length |
| Mapping | Unique mapper only | Multi-mappers have ambiguous locus |
| Region containment | See Section 10.3 | Boundary-spanning fragments need special handling |
| Strand | POS or NEG (not ambiguous) | Consistent with region counting exclusions |

**Why a table, not histograms?**

- **Memory:** A per-region histogram with ~2000 bins × 684K regions = ~5.5 GB.
  Even with sparse storage, most regions have enough fragments to populate
  many bins. A flat observation table with one (region_id, frag_len) row per
  eligible fragment is far more compact — for 30M unspliced unique fragments,
  it's ~240 MB (two int32 columns).
- **Flexibility:** Downstream models can bin the observations however they
  want (variable-width bins, KDE, parametric fits) without being locked into
  a fixed histogram resolution.
- **Aggregation:** Per-region histograms can be computed from the table on
  demand: `table.groupby('region_id')['frag_len'].value_counts()`.
- **Sparse regions:** Regions with few fragments contribute few rows (no
  wasted empty bins). Regions with zero fragments contribute zero rows.

### 10.3 Region Containment for FL Observations

**Question:** Does a fragment need to be entirely contained within a single
region to be usable for fragment-length tabulation?

**Answer: Yes, for the first implementation.** The fragment-length observation
is meaningful only when we can attribute it to a single region. A fragment
spanning two regions has an ambiguous region assignment for FL purposes — its
length characterizes neither region individually.

**Containment test:** For an unspliced fragment with interval
`[genomic_start, genomic_start + genomic_footprint)`:

```python
hits = region_cr.overlap(ref, genomic_start, genomic_start + genomic_footprint)
if len(hits) == 1:
    # Fragment entirely within one region → emit FL observation
    region_id = hits[0].label
    frag_len = genomic_footprint
    emit(region_id, frag_len)
else:
    # Fragment spans region boundary → skip for FL tabulation
    # (still counted in fractional region evidence)
    pass
```

**Impact of this restriction:** From the human genome region statistics, the
median intronic region is ~1,500 bp and the median intergenic region is ~11.5
kbp. Typical fragments are 200–500 bp. Most unspliced fragments in
intronic/intergenic regions will be fully contained within a single region.
The restriction primarily excludes fragments at exonic region boundaries
(median exonic region: 163 bp), where many fragments will straddle boundaries.
This is acceptable because:

1. Exonic regions are RNA-dominated — they're less informative for gDNA FL
   modeling anyway
2. The restriction affects only FL observations, not the fractional region
   counts (which use the full overlap algorithm)
3. We retain the vast majority of intronic and intergenic FL observations

**Future refinement:** If needed, boundary-spanning fragments could contribute
FL observations to the region with majority overlap, or could be fractionally
assigned (e.g., a 300bp fragment with 200bp in region A and 100bp in region B
could contribute to region A's FL histogram with weight 200/300).

### 10.4 Integration with Region Evidence

Fragment-length tabulation runs alongside the fractional region counting,
using the same fragment iteration loop. In the prototype (Stage 1):

```python
def count_region_evidence(bam_path, index):
    n_regions = len(index.region_df)
    counts = np.zeros((n_regions, 4), dtype=np.float64)

    # FL observation lists (will be concatenated into arrays)
    fl_region_ids = []
    fl_frag_lens = []

    # ... for each fragment:
    if is_unspliced and is_unique and len(region_hits) == 1:
        fl_region_ids.append(region_hits[0].label)
        fl_frag_lens.append(frag_len)

    # Build FL observation table
    fl_table = pd.DataFrame({
        'region_id': np.array(fl_region_ids, dtype=np.int32),
        'frag_len': np.array(fl_frag_lens, dtype=np.int32),
    })

    return region_counts_df, fl_table
```

Both outputs are written to Feather files:
- `region_evidence.feather` — per-region fractional counts (Section 5)
- `region_frag_lengths.feather` — per-fragment FL observations

### 10.5 Downstream Usage

The FL observation table enables several downstream analyses:

1. **Global intergenic FL distribution:** Filter to intergenic regions
   → pure gDNA reference FL distribution for the mixing model
2. **Per-region FL characterization:** Group by region_id → per-region
   empirical FL distributions (or parametric fits) for gDNA fraction estimation
3. **Region-type aggregate FL distributions:** Group by region type
   (exonic/intronic/intergenic) → aggregate FL distributions showing
   gDNA enrichment in intronic/intergenic vs. exonic
4. **Fragment-length filtering:** Identify fragments with lengths outside the
   expected library-prep range → candidate gDNA or degradation artifacts

### 10.6 Memory and Performance

| Component | Estimate |
|-----------|----------|
| FL table rows (30M unspliced unique frags, ~80% contained) | ~24M rows |
| FL table memory (2 × int32 × 24M) | ~192 MB |
| FL table Feather file | ~100 MB (with compression) |
| Additional per-fragment work | ~1 cgranges query (already done for counting) |

The containment check is essentially free since the region overlap query is
already performed for fractional counting. We just check whether `len(hits)
== 1` and, if so, emit the observation.

### 10.7 FL Tabulation Tests

| Test | Description |
|------|-------------|
| `test_fl_single_region_fragment` | Contained fragment → FL observation emitted |
| `test_fl_boundary_spanning_excluded` | Fragment spanning 2 regions → no FL observation |
| `test_fl_spliced_excluded` | Spliced fragment → no FL observation |
| `test_fl_multimapper_excluded` | NH > 1 → no FL observation |
| `test_fl_correct_frag_len` | `frag_len` matches `genomic_footprint` |
| `test_fl_empty_bam` | No fragments → empty FL table |
| `test_fl_intergenic_fragments` | Intergenic fragments produce FL observations |

---

## 11. Open Questions and Future Work

### 11.1 Multi-Mapper Handling

Multi-mapped fragments are excluded from region counting in both stages.
Future work could allocate multi-mapper counts probabilistically after EM,
using posterior transcript weights to distribute fractional counts across
regions.

### 11.2 Vectorized Python Fallback

If Stage 1 performance is insufficient for iterative development, a vectorized
`np.searchsorted` approach on sorted region boundaries can replace per-fragment
cgranges queries for unspliced fragments. This trades generality for speed:

```python
# Per-reference sorted boundaries
boundaries = region_df.loc[mask, 'start'].values
region_ids = region_df.loc[mask, 'region_id'].values
# Vectorized midpoint lookup
mids = starts + footprints // 2
idx = np.searchsorted(boundaries, mids, side='right') - 1
```

This is only valid for single-block (unspliced) fragments and doesn't handle
boundary-spanning fractional counting. It could serve as a fast approximation
for development iteration before the full fractional counting is needed.

### 11.3 Spliced Annotation Detail

The prototype classifies fragments as spliced/unspliced only. The production
stage has access to `splice_type` from `_resolve_core()` and can distinguish
annotated vs. unannotated splicing. Additional columns (`n_spliced_annot_pos`,
`n_spliced_annot_neg`) can be added at Stage 2 when the information is
available at zero extra cost.

### 11.4 Per-Region FL Histogram Sparse Storage

If the flat observation table proves too large for very deep sequencing
experiments (>100M fragments), a sparse histogram representation could be
used instead:

```python
# Sparse COO format: (region_id, frag_len_bin) → count
# Only non-zero bins stored
sparse_hist = scipy.sparse.coo_matrix(
    (counts, (region_ids, frag_len_bins)),
    shape=(n_regions, max_frag_len_bin)
)
```

For typical experiments the flat table is preferred for its simplicity and
flexibility. The sparse histogram option is noted for future scaling needs.

---

## 12. Summary

| Aspect | Decision |
|--------|----------|
| **Counting method** | Fractional overlap (overlap_bp / total_aligned_length) |
| **Spliced handling** | Per-block overlap, NOT midpoint — intron gaps excluded |
| **Intergenic fragments** | Counted individually at AssembledFragment stage |
| **Count type** | float32 (fractional; accumulated in float64, stored as float32) |
| **Strand source** | Genomic alignment strand from ExonBlock.strand (R2 flip applied) |
| **Splice source** | `has_introns()` (prototype); `splice_type` (production) |
| **Exclusions** | Chimeric, multi-mapped (verified per Section 3.3), ambiguous-strand |
| **Output** | Per-region DataFrame with 4 float32 count columns + FL observation table |
| **FL tabulation** | Per-fragment (region_id, frag_len) table for contained unspliced frags |
| **Module** | `src/rigel/region_evidence.py` (both stages) |
| **Stage 1** | Python prototype with pysam BAM reader |
| **Stage 2** | C++ RegionAccumulator in WorkerState, single-pass |
| **Insertion point** | After `build_fragment()`, before `_resolve_core()` |
| **Transition criterion** | Validated output + performance need |
| **Synergy** | Region + transcript queries share coordinate system; combined pass possible |
