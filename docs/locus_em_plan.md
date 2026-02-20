# Locus-Level EM: Implementation Plan

**Date:** 2026-02-18
**Status:** Draft

## 1. Problem Statement

The current EM runs a **single global instance** over ALL ambiguous
fragments in the BAM (`2*N_t + 2*N_g` components, all units in one CSR).
This has two correctness issues:

1. **Missing interdependencies.** A fragment overlapping the intron of
   Gene A and the exon of Gene B creates an implicit coupling between A
   and B. The EM must consider BOTH genes simultaneously to resolve the
   fragment. Currently, the EM does solve globally, but it confuses
   **statistical independence** with **fragment independence** — loci
   that share no fragments should be independent EM problems, while loci
   that share fragments must be solved jointly.

2. **Gene-centric gDNA model is wrong.** The current model creates one
   gDNA shadow per gene, sized to the gene's genomic span. But gDNA
   contamination is a property of the **genome**, not of genes. gDNA
   fragments originate from chromosomal DNA and should be modeled as
   such. The correct model builds a gDNA shadow from the actual
   **genomic footprint** of unspliced fragments in the locus —
   clusters of aligned read intervals on the reference — and reports
   gDNA in genomic coordinates, not attributed to genes.

3. **Scalability.** Running one global EM over the entire transcriptome
   (200k+ components for human) is wasteful when most loci are
   independent. Per-locus EM scales with the largest locus, not the
   entire genome.

## 2. Key Insight: Fragment-Connected Loci

Overlapping transcripts are grouped into **loci** via fragment
connectivity:

```
Fragments:  F1{A,B,C}  F2{A,C}  F3{D,G}  F4{C,D}  F5{M,N}

Union-find on transcript sets:
  F1 → merge(A,B,C)
  F2 → merge(A,C)         [already same component]
  F3 → merge(D,G)
  F4 → merge(C,D)         → merges {A,B,C} ∪ {D,G} = {A,B,C,D,G}
  F5 → merge(M,N)

Result:
  Locus L1 = {A, B, C, D, G}   (fragments F1,F2,F3,F4)
  Locus L2 = {M, N}            (fragment F5)
```

Loci L1 and L2 are **truly independent** — they share no fragments.
Each locus gets its own EM instance.

## 3. Architectural Overview

### Current Flow
```
BAM → fragments → buffer
                     ↓
          _scan_and_build_em_data  →  Global EMData
                                        ↓
                              ReadCounter.run_em (single EM)
                                        ↓
                              ReadCounter.assign_ambiguous
```

### Proposed Flow
```
BAM → fragments → buffer
                     ↓
          _scan_and_build_em_data  →  Global EMData  (unchanged)
                     ↓
          _build_loci(em_data)     →  list[Locus]    (NEW)
                     ↓
          for locus in loci:                          (NEW)
              _build_locus_em(locus, global_em_data)
              ReadCounter.run_locus_em(locus_em_data)
              ReadCounter.assign_locus_ambiguous(locus_em_data)
```

The first pass (`_scan_and_build_em_data`) is **unchanged** — it still
builds one global CSR of all ambiguous units/candidates. A new
**post-hoc partitioning step** groups the global CSR into independent
per-locus sub-problems.

## 4. Data Structures

### 4.1 `Locus` (new dataclass)

```python
@dataclass
class Locus:
    locus_id: int                       # sequential label
    transcript_indices: np.ndarray      # int32 — global t_indices in this locus
    gene_indices: np.ndarray            # int32 — global g_indices in this locus
    unit_indices: np.ndarray            # int32 — EM unit indices (rows) in global CSR
    gdna_intervals: list[tuple[str, int, int]]  # merged genomic intervals (ref, start, end)
    gdna_footprint_bp: int              # sum of interval sizes = gDNA effective length
    gdna_init_count: float              # unspliced fragment count initializing the gDNA shadow
```

### 4.2 `LocusEMData` (new dataclass)

A locus-local EM sub-problem with **locally-renumbered** component
indices:

```python
@dataclass
class LocusEMData:
    locus: Locus
    # CSR arrays (sliced + reindexed from global EMData)
    offsets: np.ndarray          # int64[n_local_units + 1]
    t_indices: np.ndarray        # int32[n_local_candidates] — LOCAL indices
    log_liks: np.ndarray         # float64[n_local_candidates]
    count_cols: np.ndarray       # uint8[n_local_candidates]
    locus_t_indices: np.ndarray  # int32[n_local_units]
    locus_count_cols: np.ndarray # uint8[n_local_units]
    # Component layout (local)
    n_transcripts: int           # |transcript_indices|
    n_components: int            # 2*n_t + 1 (mRNA + nRNA + 1 gDNA)
    # Mapping local ↔ global
    local_to_global_t: np.ndarray  # int32[n_transcripts]
    # Init vectors (local)
    unique_totals: np.ndarray    # float64[n_components]
    nrna_init: np.ndarray        # float64[n_transcripts]
    gdna_init: float             # scalar: unspliced fragment count for locus gDNA
    effective_lengths: np.ndarray # float64[n_components]
    prior: np.ndarray            # float64[n_components]
    # gDNA output metadata
    gdna_intervals: list[tuple[str, int, int]]  # genomic intervals for reporting
```

**Component layout per locus:**

| Local index range | Component | Description |
|---|---|---|
| `[0, n_t)` | mRNA | One per transcript in the locus |
| `[n_t, 2*n_t)` | nRNA | One per transcript (nascent shadow) |
| `[2*n_t]` | gDNA | **One** merged gDNA shadow for the locus |

Total components = `2 * n_t + 1`

### 4.3 Fragment-Cluster gDNA Shadow

The gDNA shadow is **not** derived from gene annotations. It is built
from the actual genomic footprint of fragments in the locus.

**Algorithm: build gDNA footprint**

1. Collect **all** fragment exon blocks (alignment blocks) from every
   fragment in the locus. Both spliced and unspliced fragments
   contribute their aligned blocks to define the covered region.
2. Merge overlapping/adjacent intervals into a sorted, non-overlapping
   set of genomic intervals (standard interval merge).
3. The gDNA **effective length** = sum of merged interval sizes:
   `gdna_footprint_bp = Σ (end - start) for each merged interval`

**Example:**
```
F1 spliced mRNA:  (100,200), (1000,1100), (2000,2100)
F2 unspliced:     (500,700)
F3 unspliced:     (900,1300)

All exon blocks:  (100,200), (500,700), (900,1300), (1000,1100), (2000,2100)
Merge overlaps:   (100,200), (500,700), (900,1300), (2000,2100)
   Note: (900,1300) absorbs (1000,1100) since 900 < 1000 < 1100 < 1300

gdna_footprint_bp = (200-100) + (700-500) + (1300-900) + (2100-2000) = 800
```

**gDNA initialization:** The gDNA shadow is initialized by the count
of **unspliced** fragments in the locus. Spliced fragments are
strongly indicative of RNA; unspliced fragments are ambiguous between
RNA and gDNA. The strand-corrected init formula uses antisense
unspliced fragment counts (as the current per-gene init does), but
aggregated across the entire locus rather than per-gene.

**gDNA reporting:** Fragments classified as gDNA are reported in
**genomic coordinates** (chromosome, start, end) — they are NOT
attributed to genes. gDNA is a property of the reference genome, not
of gene annotations.

## 5. Implementation Steps

### Step 1: Locus Builder (`_build_loci`)

**File:** `pipeline.py` (new function)

Build connected components from the global EMData's unit→transcript
mapping:

```
Input:  EMData (global), HulkIndex
Output: list[Locus]

Algorithm:
  1. Initialize Union-Find over transcript indices [0, N_t).
  2. For each EM unit u:
       candidates = em_data.t_indices[offsets[u]:offsets[u+1]]
       # Extract TRANSCRIPT indices (not nRNA/gDNA shadows)
       t_set = {c for c in candidates if c < nrna_base}
       union all pairs in t_set
  3. Also union deterministic unique fragments (from counter.unique_counts):
       For each transcript t with unique_counts[t].sum() > 0:
         → add to its own 1-transcript locus (or merge if overlapping).
       Actually: deterministic uniques don't create cross-transcript
       connections (they have exactly 1 transcript), so they are
       implicitly singletons unless connected by ambiguous fragments.
  4. Extract connected components → list of transcript sets.
  5. For each component, derive gene_indices from t_to_g_arr.
  6. Compute gDNA footprint:
       a. Collect all fragment exon blocks for every unit in the locus.
       b. Merge overlapping intervals.
       c. gdna_footprint_bp = Σ (end - start) for merged intervals.
       d. gdna_init_count = count of unspliced fragments in locus
          (from strand-corrected antisense accumulator).
  7. Assign unit indices: for each EM unit u, find locus of its
     first transcript index → map u to that locus.
```

**Singleton loci** (transcripts with only unique counts, no ambiguous
units) are skipped — they don't need EM.

### Step 2: Locus EM Data Builder (`_build_locus_em_data`)

**File:** `pipeline.py` or `counter.py` (new function)

For each `Locus`, extract and renumber the global EMData into a
local `LocusEMData`:

```
Input:  Locus, EMData (global), ReadCounter, HulkIndex
Output: LocusEMData

Algorithm:
  1. Build local ↔ global mappings:
       local_t[0..n_t) ↔ locus.transcript_indices
       (gene_indices retained for nRNA init, but NOT for gDNA)

  2. Build global→local index map:
       global_to_local_component[global_t] = local_t       (mRNA)
       global_to_local_component[nrna_base + global_t] = n_t + local_t  (nRNA)
       All gDNA candidates (any global gdna index) → 2*n_t  (single gDNA shadow)

  3. Slice and renumber CSR:
       For each unit u in locus.unit_indices:
         candidates = global t_indices[offsets[u]:offsets[u+1]]
         remap each candidate through global_to_local_component
         append to local offsets/t_indices/log_liks/count_cols

  4. Build local init vectors:
       unique_totals[0:n_t] = counter.unique_counts[global_t].sum(axis=1)
       unique_totals[n_t:2*n_t] = counter.nrna_init[global_t]
       unique_totals[2*n_t] = locus.gdna_init_count

  5. Build local effective lengths:
       eff_len[0:n_t] = mRNA effective lengths for local transcripts
       eff_len[n_t:2*n_t] = nRNA effective lengths (transcript spans)
       eff_len[2*n_t] = locus.gdna_footprint_bp - mean_frag + 1

  6. Build local prior:
       prior[0:n_t] = alpha  (Dirichlet pseudocount)
       prior[n_t:2*n_t] = alpha (zeroed for single-exon transcripts)
       prior[2*n_t] = 0 if gdna_init_count == 0 else alpha
```

### Step 3: Per-Locus EM Solver

**File:** `counter.py` (new method `run_locus_em`)

Run the standard EM algorithm on a `LocusEMData`:

```python
def run_locus_em(self, locus_em: LocusEMData, em_iterations=1000):
    """Run EM for a single locus.

    Returns converged theta and alpha for this locus.
    Identical algorithm to current run_em but on local arrays.
    """
```

This is a refactored version of the current `run_em()` that operates
on local vectors instead of the global `(2*N_t + 2*N_g)`-sized theta.
The algorithm (E-step, M-step, convergence) is identical.

### Step 4: Per-Locus Posterior Assignment

**File:** `counter.py` (new method `assign_locus_ambiguous`)

Run posterior assignment for a single locus, then scatter results
back into the global ReadCounter arrays:

```python
def assign_locus_ambiguous(self, locus_em: LocusEMData, ...):
    """Assign ambiguous fragments within one locus.

    Uses converged locus theta to compute posteriors.
    Scatters fractional counts back to global arrays using
    locus_em.local_to_global_t mapping.
    """
```

The classification is:
- `local_idx < n_t` → mRNA → scatter to `em_counts[global_t, col]`
- `local_idx ∈ [n_t, 2*n_t)` → nRNA → scatter to `nrna_em_counts[global_t]`
- `local_idx == 2*n_t` → gDNA → report in **genomic coordinates**
  (chromosome, fragment start, fragment end). gDNA is NOT attributed
  to any gene — it is a property of the reference genome.

### Step 5: Orchestrate in `count_from_buffer`

**File:** `pipeline.py` (modify `count_from_buffer`)

Replace the single `counter.run_em(em_data)` + `counter.assign_ambiguous(em_data)` call with:

```python
# Build global EM data (unchanged)
em_data = _scan_and_build_em_data(...)

# Compute gDNA/nRNA inits (unchanged)
gdna_init = _compute_gdna_init(...)
nrna_init = _compute_nrna_init(...)

# NEW: Build loci via fragment connectivity
loci = _build_loci(em_data, index)

logger.info(f"[LOCI] {len(loci)} loci from {em_data.n_units} units "
            f"(largest: {max(len(l.transcript_indices) for l in loci)} transcripts)")

# NEW: Per-locus EM + assignment
for locus in loci:
    locus_em = _build_locus_em_data(locus, em_data, counter, index)
    counter.run_locus_em(locus_em, em_iterations=em_iterations)
    counter.assign_locus_ambiguous(locus_em, ...)
```

### Step 6: gDNA Reporting

gDNA is **not attributed to genes**. When a fragment is classified as
gDNA by the locus EM, it is reported using its genomic coordinates:

- **Chromosome** (reference name)
- **Start** (leftmost aligned position)
- **End** (rightmost aligned position)

The output includes a per-locus gDNA summary:
- `locus_id`, `reference`, `gdna_intervals` (the merged footprint
  intervals), `gdna_footprint_bp`, `gdna_fragment_count`

This is a fundamental shift from the current model where gDNA counts
are accumulated per-gene. The new model treats gDNA as a genomic
background signal, not a gene-level quantity.

## 6. What Stays the Same

- **BAM scanning (`scan_and_buffer`)**: Unchanged. Fragments are buffered
  in arrival order.
- **Buffer data structures**: Unchanged. No locus-level buffering.
- **Fragment resolution**: Unchanged. `resolve_fragment()` still returns
  per-fragment transcript sets.
- **`_scan_and_build_em_data`**: Unchanged. Still builds one global CSR.
  The locus grouping is a post-hoc partitioning of this global CSR.
- **Deterministic assignment**: SPLICED_ANNOT + unique fragments still
  get deterministic counts. They do NOT enter the EM. Their transcript
  identity is accounted for in `unique_totals` when building locus init
  vectors.
- **Pre-EM accumulators for nRNA**: `transcript_intronic_sense/antisense`
  — unchanged. These feed nRNA init per-transcript.
- **Note**: The current per-gene gDNA accumulators (`gene_sense_all`,
  `gene_antisense_all`) will be replaced by per-locus unspliced
  fragment counting derived from the fragment cluster footprint.
- **Strand/insert models**: Shared across all loci (they are library-level
  properties, not locus-specific).

## 7. What Changes

| Component | Current | Proposed |
|---|---|---|
| EM scope | Global (single EM over all units) | Per-locus (independent EM per connected component) |
| gDNA shadow | Per-gene (N_g shadows) | Per-locus: 1 shadow built from fragment cluster footprint |
| gDNA effective length | Gene span | Sum of merged fragment-interval sizes in locus |
| gDNA init | Per-gene antisense strand correction | Per-locus unspliced fragment count |
| gDNA reporting | Attributed to genes | Reported in **genomic coordinates** (not gene-attributed) |
| EM components | `2*N_t + 2*N_g` | Per-locus: `2*n_t_local + 1` |
| `ReadCounter.run_em` | One call, global theta | One call per locus, local theta |
| `ReadCounter.assign_ambiguous` | One call, global posteriors | One call per locus, local posteriors |

## 8. Edge Cases

### 8.1 Singleton Loci (One Transcript, No Ambiguity)

If a locus contains only one transcript and all its fragments are
unique-deterministic, it has no EM units. Skip EM entirely for this
locus (already handled by "no units → skip" logic).

### 8.2 Very Large Loci

Complex genomic regions (e.g., HLA, immunoglobulin clusters) may
produce loci with hundreds of transcripts and thousands of EM units.
The per-locus EM is still efficient because the theta vector is small
(`2*n_t + 1`). The E-step vectorization (`reduceat`) handles large
numbers of units well.

Log a warning for loci exceeding a transcript threshold (e.g., >100
transcripts) for debugging.

### 8.3 Cross-Chromosome Loci

Multimapper reads can connect transcripts on different chromosomes.
The locus grouping algorithm handles this naturally — the union-find
doesn't care about genomic coordinates. The gDNA footprint is already
stored as a list of `(ref, start, end)` intervals, so cross-chromosome
loci are handled without special cases — fragment exon blocks from
different chromosomes simply produce separate intervals in the list.

### 8.4 Transcripts with No Fragments

Transcripts absent from ALL fragments (zero unique + zero ambiguous)
form singleton loci with no EM units. These are trivially resolved
(count = 0). They should not be connected to any other locus.

## 9. Testing Strategy

### 9.1 Unit Tests

- **Union-find correctness**: Given fragment→transcript sets, verify
  connected components match expected loci.
- **CSR partitioning**: Verify that slicing global EMData into per-locus
  sub-problems preserves all units/candidates.
- **Component renumbering**: Verify local→global index mapping is
  bijective and correct.
- **gDNA footprint calculation**: Given fragment exon blocks, verify
  interval merging produces correct `gdna_intervals` and
  `gdna_footprint_bp` (e.g., the 800bp example from the spec).
- **gDNA init**: Verify unspliced fragment counting for locus gDNA init.

### 9.2 Integration Tests

- **Existing test_scenarios.py**: All 8 scenarios should continue to
  pass. Most scenarios have 1-2 genes and will form 1-2 loci — behavior
  should be equivalent to current global EM.
- **New scenario: OverlappingMultiGene**: Two genes with overlapping
  exonic regions on different strands, fragments spanning both. Must
  form a single locus. gDNA shadow footprint = merged exon blocks of
  all fragments. gDNA reported in genomic coords, not genes.
- **New scenario: IndependentGenes**: Two non-overlapping genes with no
  shared fragments. Must form two independent loci, each solved
  separately. Results should be identical to solving them globally.
- **Regression: paralog multimapping**: Paralogs connected by multimapper
  reads must form a single locus.
- **gDNA footprint test**: Verify that the fragment-cluster gDNA shadow
  effective length matches the sum of merged interval sizes, not gene
  spans.

### 9.3 Diagnostic

- Add locus statistics to `PipelineStats`: `n_loci`, `max_locus_genes`,
  `max_locus_transcripts`, `max_locus_units`.
- Log per-locus breakdown for debugging.

## 10. Implementation Order

1. **`_build_loci`** — Union-find + connected components (pure function,
   independent of other changes, easy to unit test).
2. **`LocusEMData`** + **`_build_locus_em_data`** — CSR partitioning +
   renumbering (pure function, unit-testable).
3. **`run_locus_em`** + **`assign_locus_ambiguous`** — Per-locus EM
   solver (refactor from existing `run_em`/`assign_ambiguous`).
4. **Orchestration in `count_from_buffer`** — Wire it all together.
5. **Tests** — Extend test_scenarios with cross-gene tests.
6. **Cleanup** — Remove global `run_em`/`assign_ambiguous` if no longer
   needed, or keep for backward compatibility.

## 11. Performance Expectations

- **Correctness**: Identical results for non-overlapping genes. Improved
  resolution for overlapping genes (proper joint EM).
- **Speed**: Faster for typical genomes. Instead of one EM with 200k+
  components, we run many small EMs (most loci have <10 transcripts).
  The largest locus (e.g., HLA) dominates runtime, but it's still
  much smaller than the full transcriptome.
- **Memory**: Lower peak memory. Each per-locus theta is tiny compared
  to the global vector. The global CSR is still built in full, but
  the per-locus slices are views or small copies.
