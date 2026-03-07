# Strand Specificity Deep Dive & Implementation Plan

## Executive Summary

**The pipeline is already protocol-agnostic by design.** Both R1-sense
(e.g. KAPA Stranded) and R1-antisense (e.g. Illumina TruSeq dUTP) libraries
are handled correctly through the same code path. The concern that
`strand_specificity` might approach 0.0 for R1-sense libraries is unfounded —
the property is defined as `max(p_r1_sense, 1 − p_r1_sense)` which is
guaranteed to lie in **[0.5, 1.0]**.

No critical bugs were found. However, there are several improvements worth
making for clarity, test coverage, and consistency.

---

## 1. Architecture Overview: How Strand Handling Works

### 1.1 R2 Strand Flip (bam_scanner.cpp)

The C++ BAM scanner flips R2's alignment strand before emitting exon blocks:

```cpp
// bam_scanner.cpp L404-408
// Process R2 reads — flip strand (R2 strand flip convention)
ref_strand = rec->is_reverse() ? STRAND_NEG : STRAND_POS;
ref_strand = (ref_strand == STRAND_POS) ? STRAND_NEG : STRAND_POS;
```

After this flip, all exon blocks carry **R1's alignment strand**, regardless of
whether R1 or R2 contributed the alignment. This is the foundational
normalization step.

### 1.2 StrandModel Training (strand_model.py)

The `StrandModel` learns the library protocol from spliced reads:

| Quantity | Definition | R1-sense (KAPA) | R1-antisense (dUTP) |
|----------|-----------|-----------------|---------------------|
| `n_same` | exon_strand == SJ_strand | High | Low |
| `n_opposite` | exon_strand ≠ SJ_strand | Low | High |
| `p_r1_sense` | P(R1 aligns with gene) | ≈ 0.95 | ≈ 0.05 |
| `p_r1_antisense` | 1 − p_r1_sense | ≈ 0.05 | ≈ 0.95 |
| **`strand_specificity`** | max(p_r1_sense, 1−p_r1_sense) | **≈ 0.95** | **≈ 0.95** |
| `read1_sense` | p_r1_sense ≥ 0.5 | True | False |

**Key insight**: `strand_specificity` is always ∈ [0.5, 1.0]. It measures
*how stranded* the library is, independent of *which direction*
(R1-sense vs R1-antisense).

### 1.3 Fragment Classification (estimator.py → is_antisense)

```python
# estimator.py L1028-1035
def is_antisense(exon_strand, ref_strand, strand_model):
    p = strand_model.strand_likelihood_int(exon_strand, ref_strand)
    return p < 0.5
```

`strand_likelihood(exon, gene)` returns `p_r1_sense` when exon==gene, else
`p_r1_antisense`. For both protocols:

| Protocol | exon==gene | p returned | is_anti | Correct? |
|----------|-----------|-----------|---------|----------|
| R1-sense | True (R1 with gene) | 0.95 | False (sense) | ✓ |
| R1-sense | False (R1 opp gene) | 0.05 | True (antisense) | ✓ |
| R1-antisense | True (R1 with gene) | 0.05 | True (antisense) | ✓ |
| R1-antisense | False (R1 opp gene) | 0.95 | False (sense) | ✓ |

**The sense/antisense classification is always gene-relative, not read-relative.**

### 1.4 Scoring Context (scoring.py → FragmentScorer)

```python
# scoring.py L190-192
log_p_sense = math.log(max(p_sense, LOG_SAFE_FLOOR))       # = log(p_r1_sense)
log_p_antisense = math.log(max(p_antisense, LOG_SAFE_FLOOR))
r1_antisense = p_sense < 0.5   # True for R1-antisense, False for R1-sense
```

In the scoring loop (scan.py L186-194):
```python
same = (exon_strand == t_strand)
log_strand = log_p_sense if same else log_p_antisense
is_anti = r1_antisense if same else (not r1_antisense)
```

This correctly computes both the strand likelihood (for scoring) and the
gene-relative antisense flag (for counting).

### 1.5 nRNA/gDNA Estimation Formulas

Both the strand-based and hybrid estimators use:

- **Inputs**: Gene-relative `sense`/`antisense` counts (from §1.3)
- **Parameter**: `strand_specificity` ∈ [0.5, 1.0] (from §1.2)

The mathematical model assumes:
```
sense_count   = nRNA × SS + gDNA/2
anti_count    = nRNA × (1−SS) + gDNA/2
```

Leading to:
```
nRNA = (sense − anti) / (2SS − 1)
gDNA = 2(anti×SS − sense×(1−SS)) / (2SS − 1)
```

Since both `sense`/`anti` are gene-relative and `SS ∈ [0.5, 1.0]`, these
formulas are **completely protocol-independent**.

### 1.6 Weight Function

Inverse-variance weighting between strand and density estimators:
```
W_strand = (2SS − 1)²
```

| SS | W_strand | W_density | Behavior |
|----|----------|-----------|----------|
| 1.0 | 1.0 | 0.0 | Pure strand |
| 0.95 | 0.81 | 0.19 | Strand-dominated |
| 0.7 | 0.16 | 0.84 | Density-dominated |
| 0.5 | 0.0 | 1.0 | Pure density |

This provides smooth, continuous degradation — no hard cutoffs needed.

---

## 2. Issues Found

### 2.1 No Bugs — Protocol Handling is Correct ✓

The full data flow from BAM → R2 flip → StrandModel training → is_antisense →
sense/anti accumulation → nRNA/gDNA formulas is protocol-agnostic at every
step. The pipeline correctly handles R1-sense, R1-antisense, and
weakly-stranded libraries through the same code path.

### 2.2 Hard Cutoffs (Already Fixed)

Two hard cutoffs were found and fixed in the previous session:

1. **`pipeline.py` L506**: `if ss < 0.7: warning` → Replaced with informational
   log showing W_strand/W_density weights.
2. **`cli.py` L307**: `> 0.7` protocol label → Removed the `> 0.5`
   unstranded gate entirely; protocol always reports the more likely direction
   as `R1-sense` or `R1-antisense`.

### 2.3 Remaining Hard Cutoff: `STRAND_DENOM_MIN = 0.2`

**Location**: `locus.py` L40

Used by `compute_nrna_init()` and `compute_gdna_rate_from_strand()` to gate
the strand-only estimator when `2SS − 1 ≤ 0.2` (i.e., SS ≤ 0.6).

**Assessment**: This is a **reasonable numerical guard** for the non-hybrid
functions. The hybrid functions (used for the EB prior system) already use the
smooth `W = (2s−1)²` weighting with `_STRAND_DENOM_EPS = 0.01`, which handles
low SS gracefully. The 0.2 threshold only affects `compute_nrna_init()` (used
for per-transcript nRNA initialization) and `compute_gdna_rate_from_strand()`
(used as a fallback). Both fall back to 0.0, which is conservative.

**Recommendation**: Keep as-is. The cost of removing this guard (numerical
instability near SS=0.5) outweighs the benefit. The hybrid estimator already
provides the smooth fallback.

### 2.4 Inconsistent Threshold Constants

| Constant | Value | File | Functions |
|----------|-------|------|-----------|
| `STRAND_DENOM_MIN` | 0.2 | locus.py | `compute_nrna_init`, `compute_gdna_rate_from_strand` |
| `_STRAND_DENOM_EPS` | 0.01 | estimator.py | `_compute_hybrid_nrna_frac_vec` |
| `_STRAND_DENOM_EPS` | 0.01 | locus.py | `compute_gdna_rate_hybrid` |

The non-hybrid functions use a stricter threshold (0.2) while the hybrid
functions use a lenient one (0.01). This is intentional: the hybrid functions
have a density fallback (`1−W`) while the non-hybrid ones don't.

**Recommendation**: Add a comment explaining the difference. No code change
needed.

### 2.5 Missing R1-antisense Test Coverage

The test suite exercises R1-sense libraries thoroughly but has **no end-to-end
test with an R1-antisense library** (p_r1_sense < 0.5). While
`test_strand_model.py` has a `test_strong_rf_library` unit test, no integration
or pipeline test exercises the R1-antisense path through `is_antisense()`,
count accumulation, and EM estimation.

**Recommendation**: Add simulator support for R1-sense mode and create
integration tests. (See §3.2 below.)

### 2.6 Simulation Code: R1-antisense Only

The simulator (`sim/reads.py`) hardcodes the R1-antisense library convention:
```python
# R1 is reverse-complement of 3′ end, R2 is sense from 5′ end.
r1_seq = reverse_complement(frag_seq[-read_len:])
r2_seq = frag_seq[:read_len]
```

The strand-flip mask (for degraded SS) swaps R1↔R2, which correctly models
degradation for any protocol, but the baseline is always R1-antisense.

**Recommendation**: Add a `read1_sense` flag to `SimConfig` that swaps
the R1/R2 assignment to model R1-sense libraries. (See §3.1.)

### 2.7 Naming Clarity in FragmentScorer ✅ (Resolved)

The `anti_flag` field in `FragmentScorer` was renamed to `r1_antisense` for
clarity. The updated field names:

| Field | Meaning |
|-------|---------|
| `log_p_sense` | `log(p_r1_sense)` — log P(R1 aligns with gene) |
| `log_p_antisense` | `log(1 − p_r1_sense)` |
| `r1_antisense` | `p_r1_sense < 0.5` — True for R1-antisense protocol |

The corresponding C++ member (`r1_antisense_`) and nanobind binding were
updated in `scoring.cpp` as well.

---

## 3. Implementation Plan

### 3.1 Add R1-sense Library Support to Simulator

**Priority**: Medium — enables R1-sense testing but no production impact.

**Files**: `src/rigel/sim/reads.py`, `src/rigel/sim/oracle_bam.py`

**Changes**:
1. Add `read1_sense: bool = False` field to `SimConfig`.
   - `False` → R1-antisense convention (current behavior, R1 = RC of 3′ end)
   - `True` → R1-sense convention (R1 = sense from 5′ end, R2 = RC of 3′ end)
2. In `_gen_reads_from_mrna()` and `_gen_reads_from_premrna()`:
   ```python
   if self.config.read1_sense:
       # R1-sense: R1 is sense from 5' end, R2 is RC of 3' end
       r1_seq = frag_seq[:read_len]
       r2_seq = reverse_complement(frag_seq[-read_len:])
   else:
       # R1-antisense: R1 is RC of 3' end, R2 is sense from 5' end
       r1_seq = reverse_complement(frag_seq[-read_len:])
       r2_seq = frag_seq[:read_len]
   ```
3. Update `oracle_bam.py` similarly for BAM generation.
4. Update docstrings.

**Estimated effort**: Small (< 1 hour).

### 3.2 Add R1-antisense Integration Tests

**Priority**: Medium — validates the full pipeline with R1-antisense data.

**Files**: `tests/test_strand_model.py` (or new `tests/test_r1_antisense.py`)

**Tests to add**:
1. **Unit test**: `is_antisense()` with R1-antisense StrandModel
   (p_r1_sense < 0.5). Verify sense/antisense classification is correct.
2. **Integration test**: Run the simulator with `read1_sense=True`, scan the
   BAM, verify:
   - StrandModel learns `p_r1_sense > 0.9`, `read1_sense = True`
   - `strand_specificity > 0.9`
   - Fragment sense/antisense counts match expected distribution
   - nRNA/gDNA estimates match R1-antisense results (same underlying biology)
3. **Parametrized test**: Run the same scenario with R1-sense and R1-antisense,
   verify abundance estimates are consistent.

**Estimated effort**: Medium (1-2 hours).

### 3.3 Rename Scoring Fields for Clarity ✅ (Completed)

`anti_flag` has been renamed to `r1_antisense` across Python (`scoring.py`,
`scan.py`) and C++ (`scoring.cpp` member variable, constructor parameter, and
nanobind binding argument).

### 3.4 Document the Threshold Constants

**Priority**: Low — documentation only.

**Files**: `src/rigel/locus.py`, `src/rigel/estimator.py`

**Change**: Add inline documentation explaining why `STRAND_DENOM_MIN = 0.2`
(non-hybrid gate) differs from `_STRAND_DENOM_EPS = 0.01` (hybrid W threshold),
and why both are correct for their respective contexts.

### 3.5 Fix Misleading Simulation Docstring ✅ (Completed)

The simulation module docstring has been updated to use clear R1-sense /
R1-antisense nomenclature instead of the ambiguous FR/RF/dUTP terminology.

---

## 4. Recommended Implementation Order

| Step | Task | Priority | Risk | Status |
|------|------|----------|------|--------|
| 1 | Document threshold constants (§3.4) | Low | None | Open |
| 2 | Fix sim docstring (§3.5) | Low | None | ✅ Done |
| 3 | Add R1-sense sim support (§3.1) | Medium | Low | Open |
| 4 | Add R1-antisense integration tests (§3.2) | Medium | Low | Open |
| 5 | Rename scoring fields (§3.3) | Low | Medium | ✅ Done |

Steps 1-2 are zero-risk documentation fixes. Steps 3-4 add valuable test
coverage. Step 5 renamed `anti_flag` → `r1_antisense` across Python and C++.

---

## 5. Conclusion

The rigel pipeline correctly handles **R1-sense, R1-antisense, and
weakly-stranded** libraries through a principled three-layer design:

1. **R2 strand flip** normalizes all alignments to R1's perspective.
2. **Bayesian StrandModel** learns the protocol direction from data and
   computes `strand_specificity ∈ [0.5, 1.0]` (always protocol-independent).
3. **Gene-relative sense/antisense counting** (`is_antisense()`) uses the
   learned model to classify fragments correctly for any protocol.
4. **Downstream formulas** use only gene-relative counts and
   `strand_specificity` → protocol-independent math.

The `strand_specificity` value **cannot approach 0.0** for any library type —
it is mathematically bounded to [0.5, 1.0] by the `max()` normalization. The
user's concern was reasonable but already addressed by the existing design.

The main gap is **test coverage for R1-antisense libraries**, which can be
addressed with modest effort by extending the simulator and adding integration
tests.
