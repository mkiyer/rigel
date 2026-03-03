# Plan: gDNA and nRNA Initialization Overhaul

**Date:** 2026-02-18
**Status:** Draft — Revised
**Scope:** hulkrna three-pool EM initialization pipeline

---

## 1. Overview

hulkrna is a three-pool EM-based RNA-seq quantification pipeline that
separates mature RNA (mRNA), nascent RNA (nRNA), and genomic DNA
(gDNA) contamination. The EM vector layout:

```
[0, N_t)                       — mRNA (per transcript)
[N_t, 2*N_t)                   — nRNA (per transcript)
[2*N_t, 2*N_t + 2*N_g)         — gDNA (per gene)
```

The current gDNA initialization (`_compute_simple_gdna_init`) and
nRNA initialization (`nrna_init = np.zeros(...)`) contain critical
bugs that cause systematic misclassification across all strand-
specificity regimes. This document describes the theoretical
framework, root cause analysis, and implementation plan for fixing
both initializations.

### Core Principles

Two conceptual principles govern the design:

1. **nRNA is a transcript-level phenomenon.** Nascent RNA is
   pre-mRNA being actively transcribed from a specific transcript.
   Each transcript isoform produces its own pre-mRNA with its own
   intron structure. Intronic reads are transcript-level evidence —
   a fragment overlapping intron(A) but exon(B) is nRNA evidence
   for transcript A specifically, not for transcript B. nRNA must
   be modeled, initialized, and assigned per-transcript throughout.

2. **gDNA is a locus/chromosomal/global phenomenon.** Genomic DNA
   contamination is not a property of genes or transcription — it
   is simply another molecule with overlapping sequence that must
   be separated from RNA. gDNA rates vary locally (sequence content,
   GC bias, copy number alterations — amplifications and deletions
   modify the baseline 2 copies per chromosome) and globally
   (overall gDNA contamination rate in the library). When multiple
   genes overlap on different strands, they share the same genomic
   DNA at that locus. gDNA shadows should therefore be modeled at
   the locus level (currently approximated per-gene), connected to
   chromosomal and global estimates.

---

## 2. Theoretical Framework

### 2.1 The Three-Pool Model

Every fragment in a stranded RNA-seq library originates from one of
three biological sources:

| Pool | Source | Strand behavior | Splice profile | Modeling level |
|------|--------|-----------------|----------------|----------------|
| **mRNA** | Mature messenger RNA | Stranded (gene-sense) | SPLICED_ANNOT, UNSPLICED (terminal exons) | Per transcript |
| **nRNA** | Nascent / pre-mRNA | Stranded (gene-sense) | UNSPLICED (intronic + exonic), SPLICED_UNANNOT | Per transcript |
| **gDNA** | Genomic DNA contamination | Unstranded (both strands equally) | UNSPLICED only | Per locus → gene → global |

### 2.2 The Antisense Principle

**At perfect strand specificity (SS = 1.0)**, every fragment that
aligns antisense to a gene must originate from gDNA:

- mRNA is stranded → sense only
- nRNA is stranded → sense only
- gDNA is unstranded → equal sense and antisense

Therefore:

$$\text{antisense}_g = \text{gDNA}_{\text{antisense}, g}$$

Since gDNA is unstranded, the sense half mirrors the antisense half:

$$\text{gDNA}_{\text{total}, g} = 2 \times \text{antisense}_g$$

This is the **sense projection**: we observe gDNA on the antisense
strand and project an equal contribution onto the sense strand.

**Antisense is the primary gDNA estimator.** Antisense reads arise
from exonic, intronic, AND intergenic regions — providing a dense
signal across the genome. Intergenic fragments (outside all gene
bodies) are unambiguous gDNA markers but are sparse in the human
genome due to high gene density, pseudogenes, and repeat elements.
Antisense is therefore a more robust primary signal.

Note: antisense counting for gDNA is performed at the gene level
because the EM currently has per-gene gDNA shadows (`[2*N_t,
2*N_t + 2*N_g)`). This is a practical approximation — gDNA is
really a locus-level phenomenon (see Future Work §7.1).

### 2.3 Strand-Specificity Correction

At imperfect SS (< 1.0), some sense RNA fragments are flipped to
the antisense strand by library prep noise. The observed antisense
count is a mixture:

$$\text{observed\_antisense}_g = \text{gDNA}_{\text{antisense},g} + \text{RNA\_flipped}_g$$

The flip rate is:

$$\text{flip\_rate} = 1 - \text{SS} = \min(p_{r1,\text{sense}},\; p_{r1,\text{antisense}})$$

where SS = max(p_r1_sense, 1 − p_r1_sense) from the trained
`exonic_spliced` strand model.

The corrected antisense count:

$$\text{corrected\_antisense}_g = \max\!\Big(0,\;\text{antisense}_g - \text{sense}_g \times \text{flip\_rate}\Big)$$

$$\text{gdna\_init}_g = 2 \times \text{corrected\_antisense}_g$$

**Behavior at extremes:**
- SS = 1.0: `flip_rate = 0` → `corrected = antisense` → perfect
- SS = 0.5: `flip_rate = 0.5`, pure RNA → `corrected ≈ 0` → correct

**Known limitation:** At intermediate SS, `sense` includes
gDNA_sense (not just RNA_sense), creating a small residual. For
real stranded libraries (SS ≥ 0.85), this is negligible. At SS =
0.5 (unstranded), strand-based gDNA estimation fundamentally cannot
work — this is a physics limitation, not a bug.

### 2.4 nRNA Initialization from Intronic Reads

Nascent RNA (pre-mRNA) is the only stranded source that produces
fragments overlapping introns. Both nRNA and gDNA produce intronic
fragments, but gDNA is unstranded.

**nRNA is a transcript-level phenomenon.** Each transcript isoform
has its own intron structure. A fragment overlapping intron(A) is
nRNA evidence specifically for transcript A. Therefore intronic
evidence must be tracked per-transcript, not per-gene.

For each transcript $t$:

$$\text{nrna\_init}_t = \max\!\Big(0,\;\text{intronic\_sense}_t - \text{intronic\_antisense}_t\Big)$$

The subtraction removes the gDNA contribution (symmetric across
strands), isolating the stranded nRNA signal. This is exact at
SS = 1.0 and a useful approximation at high SS.

Single-exon transcripts have no introns, so `nrna_init = 0` — nRNA
is physically identical to mRNA without introns.

### 2.5 The Regional Competition Framework

In each genomic region, the competition between pools:

**Exonic regions (unspliced fragments):**
```
total_unspliced = gDNA_antisense + gDNA_sense + nRNA_sense + mRNA_sense
               ≈ 2 × gDNA_antisense + nRNA_sense + mRNA_sense
```

**Intronic regions:**
```
total_intronic = gDNA_antisense + gDNA_sense + nRNA_sense
               ≈ 2 × gDNA_antisense + nRNA_sense
```

**mRNA anchor:** Spliced fragments are purely mRNA (gDNA/nRNA
cannot produce annotated splice junctions). The existing
deterministic routing of SPLICED_ANNOT + FRAG_UNIQUE to mRNA is
correct and unchanged.

### 2.6 Isoform Structural Ambiguity

For FRAG_AMBIG_SAME_STRAND fragments (single gene, multiple
transcripts): a fragment may be EXONIC relative to isoform A but
INTRONIC relative to isoform B. This means an unspliced sense
fragment overlapping exon(A) and intron(B) could be:

- mRNA from isoform A (exonic for A)
- nRNA from isoform A or B (genic for both)
- gDNA (unstranded)

**For gDNA estimation** (gene-level antisense counting): the
fragment's antisense status is unambiguous regardless of which
isoform it came from — all isoforms share the same gene strand.
Antisense = gDNA.

**For nRNA estimation** (transcript-level intronic counting): the
fragment carries per-candidate `intron_bp[k]`. Transcript k with
`intron_bp[k] > 0` receives a fractional intronic count
(`1/n_candidates`), reflecting the uncertainty about which
transcript produced the fragment. Transcript k with
`intron_bp[k] = 0` receives no intronic count.

---

## 3. Root Cause Analysis

### 3.1 Bug 1: `unique_counts` Only Contains SPLICED_ANNOT Entries

**Location:** `_compute_simple_gdna_init` at
[pipeline.py](src/hulkrna/pipeline.py#L862-L893);
`assign_unique` at [counter.py](src/hulkrna/counter.py#L360-L373);
routing at [pipeline.py](src/hulkrna/pipeline.py#L768-L770)

**The bug:** `_compute_simple_gdna_init` reads antisense counts from
`counter.unique_counts`. But `assign_unique()` is only called on
the deterministic path: `FRAG_UNIQUE and is_spliced_annot`
(line 768). Non-SPLICED_ANNOT unique fragments (UNSPLICED,
SPLICED_UNANNOT) are routed to the EM and never appear in
`unique_counts`.

gDNA fragments are almost always UNSPLICED → they never reach
`unique_counts` → `gdna_init` is blind to actual gDNA.

**Evidence from diagnostic sweep:**
- SS=1.0, gDNA=50: `gdna_init = 0` (invisible despite ~212
  antisense gDNA reads — all UNSPLICED, never reach unique_counts)
- SS=0.5, no gDNA: `gdna_init = 394` (flip-stranded SPLICED_ANNOT
  RNA falsely counted as antisense)

### 3.2 Bug 2: nRNA Initialized to Zero

**Location:** [pipeline.py](src/hulkrna/pipeline.py#L1046-L1047)

**The bug:** `nrna_init = np.zeros(num_transcripts)` unconditionally.
The EM must discover nRNA from likelihood alone, starting at massive
disadvantage vs gDNA (which gets `2 × antisense` init) and mRNA
(which gets `unique_counts` from SPLICED_ANNOT).

**Evidence:** At SS ≤ 0.75, up to 96% of nRNA misclassified as gDNA.

### 3.3 Bug 3: Prior Zeroing Creates Dead Zones

**Location:** [counter.py](src/hulkrna/counter.py#L399-L402)

```python
for i in range(ng):
    if self.gdna_exonic_init[i] == 0.0:
        prior[self.gdna_exonic_base_index + i] = 0.0
```

**The bug:** When `gdna_exonic_init[g] == 0.0`, the Dirichlet prior
for gDNA is zeroed. Combined with Bug 1: `init=0 → prior=0 →
theta=0 → posterior=0` forever. The EM structurally cannot recover
gDNA once the prior is zeroed.

### 3.4 Interaction Matrix

| Scenario | Bug 1 | Bug 2 | Bug 3 | Net result |
|----------|-------|-------|-------|------------|
| SS=1.0 + gDNA | UNSPLICED antisense → not in `unique_counts` → init=0 | nRNA init=0 | prior=0 → dead zone | gDNA **invisible** |
| SS=0.5, no gDNA | flip-stranded SPLICED_ANNOT → init=394 | — | prior active | False gDNA hallucination |
| SS=0.9 + nRNA | — | nRNA init=0 | — | nRNA → gDNA mislabeling ≈ 96% |

---

## 4. Implementation Steps

### Step 1: Add Pre-EM Counting Arrays

**Goal:** Accumulate (a) gene-level sense/antisense counts for gDNA
init and (b) transcript-level intronic sense/antisense counts for
nRNA init, from ALL single-gene fragments during the scan.

**Files to modify:**
- [counter.py](src/hulkrna/counter.py#L155-L220) `ReadCounter.__init__`: add arrays
- [pipeline.py](src/hulkrna/pipeline.py#L544-L850) `_scan_and_build_em_data`: add accumulation

**New arrays on `ReadCounter`:**

```python
# --- Pre-EM strand accumulators (for gDNA init) ---
# Gene-level: all single-gene fragments, any splice type.
# Used by _compute_gdna_init (strand-corrected antisense).
self.gene_sense_all = np.zeros(num_genes, dtype=np.float64)
self.gene_antisense_all = np.zeros(num_genes, dtype=np.float64)

# --- Pre-EM intronic accumulators (for nRNA init) ---
# Transcript-level: nRNA is per-transcript. Each transcript has
# its own intron structure; intron_bp[k] > 0 means this fragment
# overlaps an intron of transcript k specifically.
self.transcript_intronic_sense = np.zeros(num_transcripts, dtype=np.float64)
self.transcript_intronic_antisense = np.zeros(num_transcripts, dtype=np.float64)
```

**Accumulation logic in `_scan_and_build_em_data`:**

Inside the `for i in range(chunk.size)` loop at
[pipeline.py L759](src/hulkrna/pipeline.py#L759), for fragments
with `fc == FRAG_UNIQUE` or `fc == FRAG_AMBIG_SAME_STRAND`:

```python
# --- Pre-EM strand/intronic accumulation ---
# Independent of EM routing; runs alongside existing logic.
if fc == FRAG_UNIQUE or fc == FRAG_AMBIG_SAME_STRAND:
    g_idx = int(t_to_g[int(bf.t_inds[0])])
    gene_strand = int(index.g_to_strand_arr[g_idx])
    sm = strand_models.exonic_spliced
    is_anti = counter.is_antisense(bf.exon_strand, gene_strand, sm)

    # Gene-level sense/antisense (for gDNA init)
    if is_anti:
        counter.gene_antisense_all[g_idx] += 1.0
    else:
        counter.gene_sense_all[g_idx] += 1.0

    # Transcript-level intronic (for nRNA init)
    # Each candidate transcript k has its own intron_bp[k].
    # FRAG_UNIQUE: 1 candidate, full weight.
    # FRAG_AMBIG_SAME_STRAND: n candidates, fractional weight.
    n_cand = len(bf.t_inds)
    weight = 1.0 / n_cand
    for k, t_idx in enumerate(bf.t_inds):
        t_idx_int = int(t_idx)
        has_intron = (bf.intron_bp is not None
                      and bf.intron_bp[k] > 0)
        if has_intron:
            if is_anti:
                counter.transcript_intronic_antisense[t_idx_int] += weight
            else:
                counter.transcript_intronic_sense[t_idx_int] += weight
```

**Fragment classes included:**

| Class | Included | Rationale |
|-------|----------|-----------|
| `FRAG_UNIQUE` | Yes | Single gene, single transcript — unambiguous |
| `FRAG_AMBIG_SAME_STRAND` | Yes | Single gene, strand unambiguous; fractional per-transcript intronic |
| `FRAG_AMBIG_OPP_STRAND` | No | Cannot attribute to one gene |
| `FRAG_MULTIMAPPER` | No | Alignment uncertain |
| `FRAG_CHIMERIC` | No | Already filtered |

**Key:** This accumulation is independent of EM routing. A fragment
counted here still goes to deterministic assignment or EM as before.
The accumulation runs BEFORE the routing decision in the loop body.

### Step 2: Replace `_compute_simple_gdna_init` with `_compute_gdna_init`

**Goal:** Strand-corrected gDNA init using all-category antisense.

**File:** [pipeline.py L862-893](src/hulkrna/pipeline.py#L862-L893):
replace entire function.

**New function:**

```python
def _compute_gdna_init(
    gene_sense_all: np.ndarray,       # float64[num_genes]
    gene_antisense_all: np.ndarray,   # float64[num_genes]
    strand_models: StrandModels,
) -> np.ndarray:
    """Compute strand-corrected per-gene gDNA initialization.

    Uses all single-gene antisense fragments (any splice type),
    corrected for the expected flip-strand RNA contribution.

    gDNA is a locus-level phenomenon (unstranded genomic
    contamination). Per-gene accounting is a practical approximation
    matching the current per-gene EM shadow architecture.

    Formula:
        corrected[g] = max(0, antisense[g] - sense[g] × flip_rate)
        gdna_init[g] = 2 × corrected[g]

    At SS=1.0: flip_rate=0, corrected = antisense (perfect).
    At SS=0.5: flip_rate=0.5, corrected ≈ 0 (correct for pure RNA).
    """
    ss = strand_models.exonic_spliced.strand_specificity
    flip_rate = 1.0 - ss

    corrected = np.maximum(
        0.0,
        gene_antisense_all - gene_sense_all * flip_rate,
    )
    return 2.0 * corrected
```

**Call site update at
[pipeline.py L1036-1043](src/hulkrna/pipeline.py#L1036-L1043):**

```python
if enable_gdna_shadows:
    gdna_exonic_init = _compute_gdna_init(
        counter.gene_sense_all,
        counter.gene_antisense_all,
        strand_models,
    )
else:
    gdna_exonic_init = np.zeros(index.num_genes, dtype=np.float64)
```

### Step 3: Add `_compute_nrna_init`

**Goal:** Initialize nRNA per-transcript from transcript-level
intronic sense excess.

**File:** [pipeline.py L1046-1047](src/hulkrna/pipeline.py#L1046-L1047):
replace `nrna_init = np.zeros(...)`

**New function:**

```python
def _compute_nrna_init(
    transcript_intronic_sense: np.ndarray,       # float64[num_transcripts]
    transcript_intronic_antisense: np.ndarray,   # float64[num_transcripts]
    transcript_spans: np.ndarray,                # float64[num_transcripts]
    effective_lengths: np.ndarray,               # float64[num_transcripts]
    mean_frag: float,
) -> np.ndarray:
    """Compute per-transcript nRNA initialization from intronic evidence.

    nRNA is a transcript-level phenomenon: each transcript isoform
    produces its own pre-mRNA with its own intron structure.

    Formula:
        nrna_init[t] = max(0, intronic_sense[t] - intronic_antisense[t])

    The antisense subtraction removes the estimated gDNA contribution
    (gDNA is unstranded → equal sense/antisense in introns).

    Single-exon transcripts (no introns) get nrna_init = 0 because
    nRNA is physically identical to mRNA without introns.
    """
    nrna_init = np.maximum(
        0.0,
        transcript_intronic_sense - transcript_intronic_antisense,
    )

    # Zero for single-exon transcripts (no introns → no intronic evidence)
    exonic_approx = effective_lengths + mean_frag - 1.0
    intronic_span = np.maximum(transcript_spans - exonic_approx, 0.0)
    nrna_init[intronic_span <= 0] = 0.0

    return nrna_init
```

**Call site:**

```python
nrna_init = _compute_nrna_init(
    counter.transcript_intronic_sense,
    counter.transcript_intronic_antisense,
    transcript_spans,
    effective_lengths,
    mean_frag,
)
```

### Step 4: Fix EM Prior Zeroing

**Goal:** Remove dead-zone creation when gDNA init = 0.

**File:** [counter.py L399-402](src/hulkrna/counter.py#L399-L402):
remove the gDNA prior zeroing loop.

**Current code (remove):**

```python
for i in range(ng):
    if self.gdna_exonic_init[i] == 0.0:
        prior[self.gdna_exonic_base_index + i] = 0.0
```

**Rationale:** The old zeroing was a workaround for Bug 1 — since
`gdna_init` was unreliable, zeroing prevented hallucination. With
the corrected strand-aware init from Step 2, this workaround is
counterproductive. The default Dirichlet pseudocount (0.5) provides
a small floor so the EM can discover gDNA from likelihood even when
init is zero (e.g., gene-ambiguous antisense not counted in Step 1).

**Keep:** nRNA prior zeroing for single-exon transcripts
([counter.py L406-419](src/hulkrna/counter.py#L406-L419)) — still
correct and necessary. Single-exon nRNA is physically identical
to mRNA; allowing nRNA mass would create identifiability problems.

### Step 5: Expose Pre-EM Counts for Diagnostics

**Goal:** Make the new per-gene and per-transcript accumulators
accessible in pipeline output.

**Files:**
- [counter.py](src/hulkrna/counter.py#L762-L825) `get_gene_counts_df`:
  add gene-level columns `n_sense_all`, `n_antisense_all`
- [counter.py](src/hulkrna/counter.py#L876-L894) `gdna_summary`:
  add totals for all new accumulators
- Per-transcript output (if applicable): add
  `intronic_sense`, `intronic_antisense` from the transcript-level
  arrays

### Step 6: Update Tests and Diagnostic Sweep

**6a. Re-run diagnostic sweep:**

```bash
conda run -n hulkrna python scripts/diagnostic_two_exon_ctrl.py
```

Compare against baseline
[two_exon_ctrl_diagnostic.json](docs/two_exon_ctrl_diagnostic.json).

**6b. Update `TestTwoExonWithControl` tolerances**
([test_scenarios.py L1152+](tests/test_scenarios.py#L1152)):

- Tighten gDNA accuracy at SS=1.0 (no longer invisible)
- Tighten nRNA accuracy at SS≥0.9 (no longer zero)
- Tighten t1 overcounting bounds
- Tighten t2 false positive bounds

**6c. Add regression tests:**

| Test | Bug | Assertion |
|------|-----|-----------|
| `test_unspliced_antisense_counted` | 1 | `gene_antisense_all > 0` with UNSPLICED antisense fragments |
| `test_gdna_init_ss10` | 1+3 | `gdna_exonic_init.sum() > 0` at SS=1.0 with gDNA |
| `test_no_gdna_hallucination_ss05` | 1 | `gdna_exonic_init.sum() < 5` at SS=0.5, no gDNA |
| `test_nrna_init_from_introns` | 2 | `nrna_init.sum() > 0` when intronic sense > antisense |

---

## 5. Expected Outcomes

| Scenario | Current | After fix |
|----------|---------|-----------|
| SS=1.0, gDNA=50 | `gdna_init=0`, gDNA invisible, t1 +20 | `gdna_init ≈ gDNA_total`, gDNA detected, t1 error ≈ 0 |
| SS=1.0, no gDNA | `gdna_init=0` (correct by accident) | `gdna_init=0` (correct by design: no antisense from RNA) |
| SS=0.5, no gDNA | `gdna_init=394` (false positive) | `gdna_init ≈ 0` after strand correction |
| SS=0.5, nRNA=30% | nRNA=0%, 96% mislabeled as gDNA | `nrna_init > 0` from intronic sense excess |
| SS=0.9, nRNA=30% | nRNA detection ≈ 0% | nRNA detection dramatically improved |

---

## 6. Verification

```bash
# Run full test suite
conda run -n hulkrna python -m pytest tests/ -v

# Run 125-parameter diagnostic sweep
conda run -n hulkrna python scripts/diagnostic_two_exon_ctrl.py
```

**Key metrics per grid point:**

- `gdna_init > 0` when `gdna_abundance > 0` and `SS > 0.6`
- `gdna_init ≈ 0` when `gdna_abundance == 0` (±5 fragments)
- `nrna_init > 0` when `nrna_fraction > 0` and `SS > 0.6`
- Overall MAE for mRNA, gDNA, nRNA decreases vs baseline

---

## 7. Future Work (TODO)

### 7.1 Locus-Level gDNA Shadows

The current EM has one gDNA shadow per gene. When genes overlap on
different strands at the same genomic locus, each gene gets its own
independent gDNA shadow — but gDNA at that locus is shared. This
can lead to double-counting or inconsistent gDNA estimates.

**Goal:** Replace per-gene gDNA shadows with per-locus shadows that
span the union of overlapping genes. All transcripts at a locus
compete against a single shared gDNA component. This is the
architecturally correct model: gDNA doesn't belong to genes — it
belongs to the genomic region.

### 7.2 Empirical Bayes Shrinkage (gDNA)

Per-gene gDNA estimates are noisy at low coverage. Local gDNA rates
vary due to:

- **Sequence content:** GC bias in library prep
- **Copy number alterations:** Amplifications and deletions change
  the baseline 2 copies per chromosome
- **Chromatin accessibility:** Varies along the chromosome

**Approach:** Shrink per-gene gDNA toward a global estimate:

$$\hat{\text{gDNA}}_g^{EB} = w_g \times \hat{\text{gDNA}}_g^{\text{local}} + (1 - w_g) \times \hat{\text{gDNA}}^{\text{global}}$$

where $w_g$ depends on per-gene evidence strength (fragment count).
The global estimate can use genome-wide antisense rate and/or
intergenic density.

This connects the local (per-gene/locus), chromosomal, and global
levels of gDNA estimation into a coherent hierarchical framework.

### 7.3 Intergenic gDNA Floor

Use intergenic fragment density as a secondary gDNA estimator:
`gDNA_floor_g = intergenic_rate × gene_span_g`. Acts as safety net
when antisense signal is absent. Sparse in human genome but
unambiguous.

### 7.4 Post-EM Expression Filter

Flag genes with zero spliced/unique counts as likely gDNA false
positives. Reclassify after EM convergence.

### 7.5 Overlap Exponent Re-Tuning

With better initialization, optimal overlap exponents may change.
Re-tune after Steps 1–4 are merged.

---

## 8. Key Design Decisions

| # | Decision | Rationale |
|---|----------|-----------|
| D1 | Include FRAG_AMBIG_SAME_STRAND in pre-EM counting | Single gene, shared strand; antisense is unambiguous; excluding them discards significant signal at multi-isoform genes |
| D2 | **Transcript-level nRNA init** (not gene-level) | nRNA is pre-mRNA from a specific transcript. Each isoform has its own intron structure. Intronic reads are transcript-level evidence: `intron_bp[k] > 0` means intronic for transcript k specifically. Per-transcript counting is the only correct model. |
| D3 | **Gene-level gDNA init** (locus approximation) | gDNA is a locus/chromosomal/global phenomenon. Per-gene shadows are a practical approximation matching the current EM architecture. Future: locus-level shadows for overlapping genes (§7.1) |
| D4 | Remove gDNA prior zeroing; keep nRNA zeroing for single-exon | Prior zeroing was a workaround for unreliable init. With corrected init, the Dirichlet pseudocount provides adequate safety. Single-exon nRNA zeroing is physically correct (no introns → nRNA ≡ mRNA) |
| D5 | Use `exonic_spliced` model for flip rate | Trained from SPLICED_ANNOT (highest-confidence RNA markers); directly measures library prep strand accuracy |
| D6 | Fractional intronic weight for ISOFORM_AMBIG | For n candidates, each transcript k with `intron_bp[k] > 0` gets `1/n` weight. Reflects uncertainty about which transcript the fragment came from while correctly tracking per-transcript intronic evidence |
| D7 | Clamp `corrected_antisense` at zero | Prevents negative init from counting noise at low SS |
| D8 | Start simple, add empirical Bayes later | Get the basic framework correct first, then add hierarchical shrinkage for robustness |
