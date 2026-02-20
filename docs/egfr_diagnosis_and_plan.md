# EGFR Diagnosis & Overlap Scoring Improvement Plan

**Date:** 2026-02-20
**Status:** Analysis complete, ready for implementation

---

## 1. EGFR Diagnosis: Per-Isoform Error Profile

### The Two Problematic Transcripts

The EGFR locus has 14 annotated transcripts. The dominant error is a
**systematic misallocation between two main isoforms**:

| Transcript | Exons | Length | Unique BP | Truth | hulkrna | Bias | salmon | Bias |
|------------|------:|-------:|----------:|------:|--------:|-----:|-------:|-----:|
| ENST00000275493.7 (T1) | 28 | 9,905 | 441 (4.5%) | 17,163 | 15,132 | **-2,031** | 15,975 | -1,188 |
| ENST00000450046.2 (T2) | 28 | 9,700 | 236 (2.4%) | 8,929 | 10,585 | **+1,656** | 9,743 | +814 |

T1 is underestimated by ~2,031 counts. T2 is overestimated by ~1,656 counts.
The net ~375 leaks to other isoforms.

### Structural Analysis

T1 and T2 **share 26 of 27 splice junctions** and **9,464 exonic bases** (95.5% of T1, 97.6% of T2):

- **T1-only regions**: exon 1 start [0, 349) and last exon tail [192520, 192612) = **441 bp**
- **T2-only region**: exon 1 [90706, 90942) = **236 bp**
- **T1-only splice junction**: one (exon1 → exon2 for T1 uses intron [349, 123269))
- **T2-only splice junction**: one (exon1 → exon2 for T2 uses intron [90942, 123269))
- **Shared splice junctions**: 26

Multimap mode makes **zero difference** (average |MAE difference| = 4.25 counts), confirming this is purely an **isoform resolution problem**, not a multimapping problem.

### Why salmon does better

Salmon's MAE for T1 is 1,188 vs hulkrna's 2,031 — salmon is **42% better**. Salmon uses k-mer pseudo-alignment that creates fine-grained equivalence classes. Fragment from a k-mer position unique to T1 maps only to T1's equivalence class. Fragment from a shared position maps to the {T1, T2} equivalence class. Salmon's EM then distributes the shared-class fragments proportional to abundance × effective length, while the unique-class fragments anchor each transcript's abundance estimate.

hulkrna achieves similar anchoring for **spliced** fragments that span the unique splice junctions, and for unspliced fragments in the unique exonic regions. But the critical difference is what happens for **unspliced fragments in the shared 26 exons**: hulkrna treats them as equally compatible with both T1 and T2, while salmon's equivalence classes may still differentiate based on position within the shared exon (k-mer uniqueness within shared exons).

---

## 2. Root Cause: The Overlap Scoring Model

### Current Implementation

The log-likelihood for an mRNA candidate is computed in `_score_candidate()` ([pipeline.py](src/hulkrna/pipeline.py#L335)):

```python
log_lik = (
    log(P_strand)
    + log(P_insert)
    + overlap_exponent * log(overlap_frac)
)
```

where `overlap_frac = exon_bp / frag_length`.

### The Problem

For a 150bp unspliced fragment landing in one of the 26 shared exons of T1/T2:

| Signal | T1 value | T2 value | Discriminating? |
|--------|----------|----------|-----------------|
| `exon_bp` | 150 | 150 | **No** — both have full exonic overlap |
| `overlap_frac` | 1.0 | 1.0 | **No** — identical |
| `log(overlap_frac)` | 0.0 | 0.0 | **No** — zero penalty for both |
| `P_strand` | ~same | ~same | **No** — same gene, same strand |
| `P_insert` | ~same | ~same | **No** — same intron structure near fragment |

The EM has **zero discrimination** for ~95% of T1's fragments. It can only distinguish T1/T2 through:
1. The ~4.5% of fragments in T1-only exonic regions
2. The ~1/27 of spliced fragments spanning the T1-unique splice junction

This is insufficient for accurate abundance estimation.

### Your Example Scenario

> T1 (+) exons (500,1000) and (2000,2500), T2 (+) exons (400,1000), (1200,1700), (2000,2800).
> Fragment "A" at (499, 649).

Fragment A overlaps:
- **T1 exon1 [500,1000)**: overlap = min(649,1000) - max(499,500) = 149bp. *Fragment starts 1bp before T1's exon.*
- **T2 exon1 [400,1000)**: overlap = min(649,1000) - max(499,400) = 150bp. *Fragment fully contained.*

Current scoring:
- T1: `exon_bp=149`, `frag_length=150`, `overlap_frac=0.993`
- T2: `exon_bp=150`, `frag_length=150`, `overlap_frac=1.0`
- log-lik difference: `overlap_exponent * (log(0.993) - log(1.0))` = `1.0 * (-0.007)` = **-0.007**

With the default `overlap_exponent=1.0` for SPLICED_ANNOT/SPLICED_UNANNOT, this 1bp difference creates a **negligible** log-likelihood difference of -0.007. The `filter_by_overlap` with `min_frac_of_best=0.99` would keep both (0.993 > 0.99 × 1.0 = 0.99).

For unspliced fragments, `overlap_exponent=10.0`, giving -0.07, still small. And `filter_by_overlap` is **skipped entirely for unspliced fragments** (line 738 of resolution.py: `if splice_type != SpliceType.UNSPLICED`).

---

## 3. The Two Approaches

### Approach A: Hard Overlap Filter (Improved)

Restore `filter_by_overlap` to remove transcripts with inferior overlap, but fix the issues that caused problems in prior versions.

**Previous issues with hard filtering:**
- Used exon-only overlap, which incorrectly removed transcripts in antisense/sense overlap scenarios (gene A's exon overlaps gene B's intron)
- Premature candidate pruning in complex loci collapsed the solution space

**What would need to change:**
- Apply `filter_by_overlap` to ALL fragments (not just spliced), but use **genic overlap** (exon + intron) as the metric for unspliced fragments
- For multi-gene fragments, apply filtering **within each gene group** rather than globally across genes
- Lower the threshold to be more aggressive (e.g., `min_frac_of_best=0.95`)

**Risk:** Hard filtering is inherently brittle. In the T1/T2 case, both have `overlap_frac ≈ 1.0` for shared-region fragments, so no threshold short of 1.0 would discriminate them. Hard filtering **cannot solve this problem**.

### Approach B: Soft Scoring Enhancement (Recommended)

Keep all transcripts as candidates but improve the scoring to provide better discrimination for overlapping isoforms.

**Key Insight:** The current overlap scoring only asks "what fraction of the fragment overlaps the transcript's exons?" It never asks the reverse question: **"what fraction of the transcript is consistent with this fragment's position?"**

Consider: a fragment landing in the middle of a 9,905bp transcript is equally expected anywhere along its length. But the information content differs — a fragment at a unique position (T1-only exon) is maximally informative, while a fragment at a shared position provides no discrimination.

What we need is not a different EM — it's a better **likelihood function** that captures positional specificity.

---

## 4. Recommended Solution: Transcript Compatibility Score

### Concept

For each (fragment, transcript) pair, compute a **compatibility score** that captures how specifically this fragment pinpoints this transcript, relative to the full candidate set.

The score has two components:

1. **Exon containment score**: How well does the fragment's genomic footprint match this transcript's exon model? (Current: `exon_bp / frag_length`. This stays.)

2. **NEW: Positional uniqueness weight**: Given this fragment's position and structure, how much of the fragment's evidence uniquely supports this transcript vs others in the candidate set?

### Implementation: Fragment-Transcript Specificity Score

For each candidate transcript `t` and fragment `f`:

```
specificity(f, t) = unique_bp(f, t) / frag_length(f)
```

where `unique_bp(f, t)` is the number of bases in the fragment that overlap exons of transcript `t` but do NOT overlap exons of any other candidate transcript in the set.

For the EGFR case:
- Fragment in T1-only region [0, 349): `unique_bp(f, T1) ≈ 150`, `unique_bp(f, T2) = 0`
  → specificity(T1)=1.0, specificity(T2)=0.0 → strongly favors T1
- Fragment in shared exon [123269, 123421): `unique_bp(f, T1) = 0`, `unique_bp(f, T2) = 0`
  → specificity(T1)=0.0, specificity(T2)=0.0 → no discrimination (correct)
- Fragment in T2-only region [90706, 90942): `unique_bp(f, T2) ≈ 150`, `unique_bp(f, T1) = 0`
  → specificity(T2)=1.0, specificity(T1)=0.0 → strongly favors T2

### How to Integrate into the Log-Likelihood

The specificity score acts as an **additional log-likelihood term**:

```python
log_lik = (
    log(P_strand)
    + log(P_insert)
    + overlap_exponent * log(max(overlap_frac, eps))
    + specificity_exponent * log(max(specificity + base_specificity, eps))
)
```

where:
- `specificity_exponent` controls the strength of the positional signal (default: 2.0)
- `base_specificity` is a small floor (e.g., 0.01) that prevents zero specificity from completely killing a candidate — shared-region fragments should still be distributed by the EM

### Why This Works for All Scenarios

**Scenario 1: Single gene, multi-isoform (EGFR-like)**
- Fragments in unique exonic regions get high specificity for the correct isoform
- Fragments in shared regions get equal (low) specificity for all → EM distributes by abundance/eff_len
- Net effect: unique-region fragments anchor the abundance; shared-region fragments are distributed proportionally

**Scenario 2: Single gene, intronic overlap (nRNA/pre-mRNA)**
- Unspliced fragment in intron: exon overlap may be partial for some isoforms
- Specificity computed on exon overlap only → intronic fragments with zero exonic overlap to all transcripts get specificity=0 for all → no distortion
- nRNA scoring path handles these separately (uses genic_frac)

**Scenario 3: Multi-gene, same strand overlap**
- Fragment in gene A exon that also overlaps gene B exon
- If genes have different exon boundaries: specificity naturally discriminates
- If exons perfectly coincide: specificity = 0 for both, EM decides by abundance/eff_len

**Scenario 4: Multi-gene, antisense overlap**
- Fragment in gene A(+) exon overlapping gene B(-) intron
- Strand model already strongly discriminates sense vs antisense
- Specificity adds no harm (gene B might have zero exon overlap → low specificity anyway)

**Scenario 5: Chimeric transcripts**
- Already handled by chimera detection in resolution (separate code path)
- Non-chimeric fragments: specificity computation is orthogonal

---

## 5. Implementation Plan

### Phase 1: Compute `unique_bp` per candidate (resolution.py)

**New function**: `compute_specificity_profile()`

After `compute_overlap_profile()` returns per-candidate `(exon_bp, intron_bp)`, compute the specificity for each candidate:

```python
def compute_specificity_profile(
    frag: Fragment,
    t_inds: frozenset[int],
    index,
) -> dict[int, float]:
    """Compute per-candidate specificity score.

    For each fragment exon block, identify bases that overlap exactly
    one candidate's exons (unique bases). The specificity for transcript
    t is sum(unique_bases_for_t) / frag_length.

    Returns dict mapping t_idx → specificity (0.0 to 1.0).
    """
```

This requires querying the index for each base position's transcript coverage, which can be done efficiently by:
1. For each fragment exon block, query all overlapping intervals
2. Build a per-base transcript set (or work with interval intersections)
3. Count bases where only one candidate transcript has exonic coverage

**Optimization**: Instead of per-base counting, work with interval arithmetic. Split each exon block at all candidate transcript exon boundaries, then for each sub-interval, count how many candidate transcripts have exonic overlap. This is O(n_candidates × n_exon_boundaries) per fragment, not O(frag_length).

### Phase 2: Store specificity in buffer (buffer.py)

Add a `specificity_flat` column alongside `exon_bp_flat` and `intron_bp_flat`. One float per (fragment, candidate) pair.

### Phase 3: Integrate into scoring (pipeline.py)

In `_score_candidate()`, add the specificity term:

```python
log_lik = (
    log(P_strand)
    + log(P_insert)
    + overlap_exponent * log(max(overlap_frac, eps))
    + specificity_exponent * log(max(specificity + base_specificity, eps))
)
```

Default parameters:
- `specificity_exponent = 2.0` (moderate weight)
- `base_specificity = 0.01` (floor to keep all candidates viable)

### Phase 4: Test scenarios

Write tests covering:
1. Single-gene multi-isoform (overlapping exons, unique exon regions)
2. Fully overlapping transcripts (all shared → specificity=0 for all → no distortion)
3. Multi-gene same-strand with partial overlap
4. Multi-gene antisense overlap
5. Chimeric transcript regions

### Phase 5: Re-run EGFR benchmark

Re-run the EGFR region benchmark and compare T1/T2 bias before/after.

---

## 6. Expected Impact

### Quantitative Estimate for EGFR

T1 has 441 unique bases out of 9,905 total. With 50,000 fragments at abundance-weighted sampling, T1 gets ~17,163 fragments. Of those:
- ~4.5% (~772) land in T1-unique regions → specificity=1.0
- ~95.5% (~16,391) land in shared regions → specificity=0.0

With the specificity term, those 772 unique fragments will have much higher log-likelihood for T1, creating a stronger anchor. The EM will then distribute the 16,391 shared fragments proportionally based on the anchored abundance estimates.

Current error: 2,031 (11.8% of truth). Expected improvement: reduce to ~500-800 (3-5%), comparable to salmon's 1,188 error which benefits from k-mer level positional resolution.

### Broader Impact

This fix addresses the same fundamental issue across ALL regions:
- **EGFR**: T1/T2 with 95% shared exons → direct fix
- **FGFR2**: Multiple isoforms with high exon sharing → similar benefit
- **HOXA_cluster**: Adjacent genes with overlapping UTRs → specificity on UTR-unique bases
- **HBB_cluster**: Homologous genes (helps multimap mode further disambiguate)

### Risk Assessment

- **Low risk**: Specificity only adds signal; it never removes candidates. Worst case: specificity=0 for all candidates (fully shared region), which adds `specificity_exponent * log(0.01) ≈ -9.2` equally to all → no effect on relative rankings.
- **Medium risk**: Over-weighting unique regions could cause excessive anchoring on noisy unique-region estimates for low-abundance transcripts. Mitigated by the `base_specificity` floor and moderate exponent.

---

## 7. Alternative: Equivalence Class Compression

A complementary optimization (not required for the specificity fix, but valuable for speed and numerical stability):

Group fragments that map to the **same set of candidate transcripts** into equivalence classes. Instead of processing each fragment individually in the EM, process each equivalence class once with its count.

Benefits:
- Faster EM convergence (fewer units)
- Better numerical stability
- Natural compatibility with the specificity score (equivalence class = same transcript set + same specificity profile → same posterior)

This is what salmon/kallisto do internally. Could be a Phase 2 enhancement.
