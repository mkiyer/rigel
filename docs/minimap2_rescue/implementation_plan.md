# Minimap2 Fragment Rescue: Implementation Plan

**Date**: 2026-03-22
**Status**: Planning
**Priority**: Critical — these issues cause MAE 7.58 vs oracle 1.85

## Executive Summary

The pristine benchmark (10M mRNA fragments, zero gDNA, zero nRNA) reveals
that rigel with minimap2 alignment has MAE=7.58 vs oracle MAE=1.85. Three
independent problems contribute to this gap. Each must be solved separately,
in the order presented, because later problems depend on earlier fixes.

| Problem | Type | Impact | Effort |
|---------|------|--------|--------|
| P1: Alignment-error overhang rescue | Enhancement | Moderate | Medium |
| P2: nRNA wins hard overhang gate (synthetic transcript bug) | **BUG** | **Severe** | Medium |
| P3: Robust fragment length estimation | Enhancement | High | Low-Medium |

---

## Problem 1: Alignment-Error Overhang Rescue

### Description

Short-read aligners (minimap2 `splice:sr`) sometimes fail to detect splice
junctions. Instead of splitting a read across two exons with an N-skip
(`148M500N2M`), the aligner produces a contiguous alignment (`150M`) that
extends 1–5 bp past the exon boundary into intronic space. This
creates **intronic overhang** on the mRNA transcript.

These fragments are not inherently wrong — they come from genuine mRNA
molecules. The aligner simply failed to detect the splice. We should
**rescue** them by recognizing that the tiny overhang is consistent with
an alignment error rather than true intronic origin.

### Evidence

From the pristine benchmark (500K read pairs traced):

| Overhang (bp) | Read count | Fraction |
|---------------|-----------|----------|
| 1 | 16,882 | 45.1% |
| 2 | 12,916 | 34.5% |
| 3 | 4,105 | 11.0% |
| 4 | 1,207 | 3.2% |
| 5 | 1,543 | 4.1% |
| 150 (full read) | 598 | 1.6% |
| Other | 220 | 0.6% |
| **Total** | **37,471** | **7.5% of all reads** |

The distribution sharply peaks at 1–3 bp, consistent with alignment boundary
errors. The oh=150 cases are reads aligned entirely within an intron (true
nascent/gDNA origin or severely misaligned reads).

### Root Cause (code path)

1. Minimap2 aligns a read as `150M` instead of `148M500N2M`
2. The read block spans `[exon_end - 148, exon_end + 2)`
3. In `resolve_context.h`, the cgranges query returns the EXON interval
   `[exon_start, exon_end)` and computes `e_bp = overlap = 148`
4. Overhang: `oh = rl - e_bp = 150 - 148 = 2`
5. The overhang penalty `oh * oh_log_pen_` adds a log-likelihood penalty

### Related Code

| File | Function / Lines | Role |
|------|-----------------|------|
| [resolve_context.h](../src/rigel/native/resolve_context.h) L786-810 | ITYPE_EXON overlap → `t_exon_bp` | Exon base-pair accumulation |
| [resolve_context.h](../src/rigel/native/resolve_context.h) L571-658 | `compute_frag_lengths()` | Gap SJ correction (overlap-based) |
| [scoring.cpp](../src/rigel/native/scoring.cpp) L800-860 | Unique mapper scoring loop | `oh = rl - ebp`, `log_lik += oh * oh_log_pen_` |

### Interaction with Gap SJ Correction

The gap SJ correction in `compute_frag_lengths()` handles a related but
different issue: when mates span an intron without CIGAR N-skips, the
inter-mate gap needs intron subtraction for correct fragment length.

An overlap-based correction (replacing the original strict containment
check) was prototyped and shows partial improvement:

| Metric | Strict containment | Overlap-based | Oracle |
|--------|-------------------|---------------|--------|
| SPLICED_ANNOT FL mode | 1000 | **294** | 299 |
| MAE | 8.64 | 7.58 | 1.85 |
| nRNA siphon | 215,943 | 171,189 | 3,625 |
| gDNA false pos | 3,005 | 39,896 | 1,840 |

The overlap fix corrects the FL model but causes a gDNA explosion because
fragments still hard-gated to nRNA get penalized by the now-correct FL model
and the EM pushes them to gDNA. This demonstrates that the gap SJ fix
alone is insufficient — the hard overhang gate bug (Problem 2) must also be fixed.

### Possible Approaches

**Option A: Soft overhang tolerance at the exon boundary**

During fragment resolution, when computing `e_bp`, extend each exon boundary
by a small tolerance (e.g., 3 bp). Reads that overhang by ≤ tolerance
would get `e_bp = rl` and `oh = 0`.

- Pro: Simple, targeted
- Con: Masking instead of fixing; could incorrectly rescue reads that
  genuinely extend into intronic space; tight coupling to alignment
  error profile

**Option B: Overhang-aware scoring (soft penalty instead of hard gate)**

Replace the hard hard overhang gate (`oh == m_min_oh`) with a likelihood-based
selection that considers all scoring signals (strand, FL, overhang, NM)
holistically. Candidates with small overhang would get a modest penalty
but would not be excluded if their other likelihoods are strong.

- Pro: Principled; works for any overhang magnitude; avoids arbitrary
  threshold
- Con: May increase equivalence class sizes; higher engineering complexity;
  interacts with Problem 2

**Option C: Gap SJ correction as rescue mechanism**

When a read's alignment block extends past an exon boundary, check whether
the overhang region overlaps an annotated splice junction. If so, "credit"
the overhang as an alignment error (splice junction that the aligner
missed) and set `oh = 0` for that transcript.

- Pro: Uses annotation knowledge; targeted at the specific failure mode
- Con: Requires SJ index lookup during resolution (perf cost); only works
  for annotated junctions

**Recommendation**: Evaluate Option A first (simplest) with benchmark
validation. If insufficient, pursue Option B (comprehensive) or C
(annotation-aware). Option B may be best combined with the Problem 2
fix since both modify the hard overhang gate mechanism.

### Testing Strategy

- Existing 998 unit tests (regression gate)
- Pristine benchmark: minimap2 MAE, nRNA siphon, gDNA false positives
- New unit test: synthetic fragment with 1–3 bp overhang, verify
  mRNA transcript is retained in candidate list
- Diagnostic script: `scripts/debug/prove_nrna_overhang_gate_v2.py` — counts
  overhang reads and pruning outcomes

---

## Problem 2: nRNA Wins Hard Overhang Gate (Synthetic Transcript Bug)

### Description

After introducing synthetic nRNA transcripts into the index, they
compete with annotated mRNA transcripts in the hard overhang gate
overhang gate. Since nRNA transcripts are single-exon spans covering
the entire genomic locus (introns included), any read that falls
within the locus has **zero overhang** against the nRNA transcript.
The hard overhang gate keeps only min-overhang candidates, so **nRNA always
wins** and mRNA is discarded — even when the fragment is genuinely
from mRNA with only 1 bp of alignment error.

**This is a critical bug** introduced when synthetic nRNA transcripts
were added to the transcript pool.

### Evidence (SMOKING GUN)

From `scripts/debug/prove_nrna_overhang_gate_v2.py` (500K pairs):

- **37,471 reads** (6.9%) have intronic overhang with no CIGAR N-skip
- **100%** of these result in nRNA winning the hard overhang gate (nrna_oh=0, mrna_oh>0)
- Most common overhang: 1 bp (45%), 2 bp (35%), 3 bp (11%)

Concrete example:
```
Read: ENST00000000233.10:182-552:f:573  (TRUTH: mRNA)
  Alignment: chr7:127590961, CIGAR=15S127M8S
  mRNA exon: [127590962, 127591088)
  Read block: [127590961, 127591088)  ← extends 1bp past exon start

  mRNA: e_bp=126, oh=1
  nRNA: e_bp=127, oh=0  ← single giant "exon" covers everything

  hard overhang gate: oh==0 only → nRNA WINS, mRNA DISCARDED
```

### Root Cause (code path)

1. Index build: `create_nrna_transcripts()` creates synthetic single-exon
   transcripts with `is_synthetic_nrna=True`
2. Interval build: `_gen_transcript_intervals()` generates ITYPE_EXON for
   the nRNA's single giant exon and ITYPE_TRANSCRIPT for its full span
3. Resolution: cgranges query returns both mRNA exon intervals and nRNA
   exon intervals for a read near an exon boundary
4. `e_bp` computation: for mRNA, overlap with real exon = `rl - overhang`;
   for nRNA, overlap with giant "exon" = `rl` → `oh = 0`
5. hard overhang gate (scoring.cpp L855-861): `m_min_oh = 0` (from nRNA), only
   candidates with `oh == 0` survive → mRNA is pruned

### Cascade Effect

```
Fragment with 1bp overhang
  → Hard overhang gate prunes mRNA (oh=1), keeps nRNA (oh=0)
    → FL computed against nRNA (single-exon spanning locus)
      → FL = genomic span (~5000 bp) instead of ~300 bp
        → FL model trains on inflated values (SPLICED_ANNOT mode=1000)
          → All scoring uses wrong FL model
            → Widespread quantification errors
```

### Related Code

| File | Function / Lines | Role |
|------|-----------------|------|
| [index.py](../src/rigel/index.py) L138-295 | `create_nrna_transcripts()` | Creates synthetic single-exon nRNA |
| [index.py](../src/rigel/index.py) L393-407 | `_gen_transcript_intervals()` | Generates EXON + TRANSCRIPT intervals for all transcripts |
| [scoring.cpp](../src/rigel/native/scoring.cpp) L551-564 | MM hard overhang gate | `mrna_min_oh` selection across all candidates |
| [scoring.cpp](../src/rigel/native/scoring.cpp) L799-861 | Unique hard overhang gate | `m_min_oh` selection across all candidates |

### Possible Approaches

**Option A: Exclude nRNA from hard overhang competition**

Add a `t_is_nrna` boolean array to `NativeFragmentScorer`. During the hard overhang gate
gate, compute `m_min_oh` only from mRNA candidates (where
`t_is_nrna[t_idx] == false`). Keep all mRNA candidates with `oh == mrna_min_oh`
AND all nRNA candidates unconditionally (their FL/strand penalties will
naturally suppress them in the EM if they're poor fits).

- Pro: Minimal code change; preserves pruning benefit for mRNA; nRNA competes
  fairly in the EM via likelihoods
- Con: Equivalence classes may grow slightly; nRNA always present

**Option B: Pool-separated hard overhang gates**

Run two completely independent hard overhang gates: one for mRNA candidates and one
for nRNA candidates. Each pool produces its own set of winners. Both sets
enter the EM.

- Pro: Clean separation; each pool optimized independently
- Con: Slightly more complex; unclear how to merge the two winner sets

**Option C: Replace hard overhang gate with likelihood-based pruning**

Instead of keeping only `oh == min_oh`, keep all candidates whose
total `log_lik` (strand + FL + overhang + NM) is within some threshold
of the best. This makes the gate soft and holistic.

- Pro: Principled; naturally handles nRNA (its FL penalty will be huge);
  also addresses Problem 1
- Con: Threshold selection is tricky; performance impact from larger
  equivalence classes; may need careful tuning

**Recommendation**: Start with Option A (targeted fix for the bug).
Consider Option C as a future enhancement that unifies Problems 1 and 2.

### Testing Strategy

- All 998 unit tests (regression)
- Pristine benchmark: mRNA recovery (target ≥99%), nRNA siphon (target <5K)
- New unit test: fragment with 1–3 bp overhang resolves to mRNA candidates
  (not nRNA-only) after hard overhang gating
- Verify FL model: SPLICED_ANNOT mode should be ~294 with oracle AND minimap2

---

## Problem 3: Robust Fragment Length Estimation

### Description

The RNA fragment length model is trained from SPLICED_ANNOT observations —
fragments classified as annotated-spliced unique mappers. This is the
"gold standard" because spliced reads definitively originate from processed
mRNA. However, the FL training pipeline has NO robustness filtering:

1. All FL observations are accepted without outlier rejection
2. If the FL computation is wrong (e.g., due to gap SJ correction failure),
   the corrupted value enters the training set
3. No median-based or trimmed-mean estimation
4. No validation against expected FL distribution shape (unimodal, ~150–500)

When Problems 1 and 2 are present, corrupted FL values (from nRNA-assigned
fragments) contaminate the SPLICED_ANNOT pool. This shifts the mode from
~294 to ~1000, making ALL subsequent scoring wrong.

### Evidence

| Configuration | SPLICED_ANNOT FL mode | SPLICED_ANNOT n_obs |
|---------------|----------------------|---------------------|
| Oracle alignment | 299 | 5,421,307 |
| Minimap2 (pre-fix) | 1000 | 1,764,892 |
| Minimap2 (overlap fix) | 294 | 2,892,014 |

The overlap fix recovered the mode but the observation count is much lower
(2.89M vs 5.42M for oracle), indicating many spliced fragments are still
being misclassified or lost. Additionally, UNSPLICED observations (3.23M)
have a reasonable mode (295), suggesting they could serve as a validation
reference.

### Root Cause

Fragment length observations flow through this pipeline:

```
BAM read pair
  → build_fragment() in bam_scanner.cpp
    → resolve + compute_frag_lengths() in resolve_context.h
      → get_unique_frag_length() returns FL if all candidates agree
        → observation stored with splice_type
```

The vulnerability is in `get_unique_frag_length()`: it returns the FL
value if all candidate transcripts agree on the same value. But if the
candidate set has been corrupted by the hard overhang gate bug (Problem 2), the
"agreed-upon" FL is the nRNA's inflated value.

Furthermore, `compute_frag_lengths()` may return inflated values when:
- Gap SJ correction fails (strict containment issue, now partially fixed
  with overlap approach)
- Fragment spans multiple introns but only one is corrected
- Fragment resolves to nRNA transcript with no introns to subtract

### Related Code

| File | Function / Lines | Role |
|------|-----------------|------|
| [resolve_context.h](../src/rigel/native/resolve_context.h) L92-103 | `get_unique_frag_length()` | Returns FL only if all candidates agree |
| [bam_scanner.cpp](../src/rigel/native/bam_scanner.cpp) L1305-1320 | FL training in scan loop | Collects FL + splice_type for unique mappers |
| [frag_length_model.py](../src/rigel/frag_length_model.py) L134-160 | `observe()` / `observe_batch()` | Accumulates FL histogram (no filtering) |
| [frag_length_model.py](../src/rigel/frag_length_model.py) L461-491 | `build_scoring_models()` | Copies SPLICED_ANNOT histogram to `rna_model` |

### Possible Approaches

**Option A: Post-hoc outlier rejection in FL training**

After collecting all FL observations during the BAM scan, apply a robustness
filter before building the model. For example:
- Compute median and MAD (median absolute deviation) from SPLICED_ANNOT
  observations
- Reject observations outside `median ± k * MAD` (e.g., k=5)
- Use the filtered histogram for the `rna_model`

- Pro: Simple; doesn't require changing the C++ scan; robust to corruption
- Con: Post-hoc; needs to happen before `build_scoring_models()`

**Option B: FL sanity check during training**

In `get_unique_frag_length()` or in the FL observation callback, add a
sanity check: reject FL values that are implausibly large (e.g., >1000 bp
for short-read PE data). This prevents corrupted values from ever entering
the model.

- Pro: Proactive; prevents corruption at source
- Con: Hard-coded threshold; may miss edge cases; may reject legitimate
  long-fragment libraries

**Option C: Cross-validate SPLICED_ANNOT against UNSPLICED**

Since both come from the same library prep, the true FL distribution must
be the same. If the SPLICED_ANNOT mode diverges significantly from the
UNSPLICED mode, flag the SPLICED_ANNOT model as corrupted and fall back
to UNSPLICED (or global).

- Pro: Self-validating; no hard thresholds
- Con: Delayed detection (only after training is complete); UNSPLICED may
  also be corrupted if the hard overhang gate bug affects them; adds complexity

**Option D: Weighted training with confidence scores**

Weight FL observations by the confidence of the assignment. Fragments with
high overhang or low log_lik get lower weight. This smoothly downweights
corrupted observations.

- Pro: Principled; no hard cutoffs
- Con: More complex; requires propagating weights through the training
  pipeline

**Recommendation**: Start with Option A (simple, effective). Add Option C
as a diagnostic check that triggers a warning. Option B is a good
belt-and-suspenders addition. Option D is overkill for now.

### Testing Strategy

- Pristine benchmark: SPLICED_ANNOT mode must be within ±10 of oracle
  (target ~294–304)
- New unit test: inject outlier FL observations, verify they are rejected
- Diagnostic: compare SPLICED_ANNOT vs UNSPLICED mode; warn if divergent

---

## Execution Order

The three problems must be solved in this order because of dependencies:

### Phase 1: Problem 2 (nRNA Hard Overhang Gate Bug) — HIGHEST PRIORITY

**Why first**: This is the most severe issue. The hard overhang gate bug causes mRNA
fragments to be assigned to nRNA, corrupting the FL model and cascading
into all downstream quantification. Fixing this alone should dramatically
reduce nRNA siphon and improve MAE.

**Steps**:
1. Research: review hard overhang gate code paths, confirm fix approach
2. Implement `t_is_nrna` array in scorer, modify hard overhang gate to exclude nRNA
3. Run unit tests (998 must pass)
4. Run pristine benchmark, verify:
   - nRNA siphon < 10K (oracle has 3.6K)
   - FL model SPLICED_ANNOT mode ~294
   - MAE significantly improved
5. Analyze results, iterate if needed

### Phase 2: Problem 3 (Robust FL Estimation)

**Why second**: Even after the hard overhang gate fix, some corrupted FL observations
may leak through (e.g., from genuinely misaligned reads). Adding
robustness ensures the FL model is always reliable.

**Steps**:
1. Implement outlier rejection in `build_scoring_models()`
2. Add cross-validation check (SPLICED_ANNOT vs UNSPLICED mode)
3. Run benchmark, verify FL model is stable
4. Test with intentionally corrupted inputs

### Phase 3: Problem 1 (Alignment-Error Overhang Rescue)

**Why third**: This is an enhancement, not a bug fix. With Problems 2
and 3 solved, fragments with small overhang will still be in the
candidate list (not hard-gated away) and will have reasonable FL. The
overhang penalty will slightly down-weight them, but the EM can
still assign them correctly. Problem 1 is about maximizing rescue
of these fragments.

**Steps**:
1. Evaluate whether the Problem 2 fix alone is sufficient
2. If MAE gap remains significant, implement soft overhang tolerance
3. Benchmark and validate
4. Consider the overlap-based gap SJ fix (already prototyped) as part
   of this phase — evaluate whether it helps or is now unnecessary

---

## Current State of Prototype Fixes

### Overlap-Based Gap SJ Correction (in resolve_context.h)

An overlap-based fix is currently applied to `compute_frag_lengths()`.
This replaces the original strict containment check with:

```cpp
int32_t overlap = std::min(he, ge) - std::max(hs, gs);
if (overlap > 0)
    t_gap_size[ti] += overlap;
```

**Status**: Applied and tested (998 tests pass). The fix corrects the
FL model (mode 1000 → 294) but causes gDNA explosion (3K → 40K) because
the hard overhang gate bug (Problem 2) is still active.

**Decision**: Keep or revert? This fix is directionally correct but its
full benefit can only be evaluated after Problem 2 is solved. Recommend
**keeping it** and re-evaluating during Phase 3.

---

## Benchmark Reference (Pristine: 10M mRNA, zero gDNA/nRNA, SS=0.95)

| Configuration | MAE | RMSE | Spearman | nRNA | gDNA | mRNA |
|---------------|-----|------|----------|------|------|------|
| Oracle | 1.85 | 30.39 | 0.880 | 3,625 | 1,840 | 9,994,535 |
| Minimap2 (pre-fix) | 8.64 | 124.51 | 0.677 | 215,943 | 3,005 | 9,779,947 |
| Minimap2 (overlap fix) | 7.58 | 121.76 | 0.719 | 171,189 | 39,896 | 9,787,810 |
| Salmon | 2.91 | 38.61 | 0.894 | — | — | — |
| Kallisto | 3.36 | 59.28 | 0.848 | — | — | — |

**Target after all fixes**: MAE < 3.0, nRNA < 5K, gDNA < 3K
(approaching oracle performance with genome alignment)

---

## Diagnostic Scripts

All diagnostic scripts live in `scripts/debug/`:

| Script | Purpose |
|--------|---------|
| `prove_nrna_overhang_gate_v2.py` | Proves nRNA wins hard overhang gate for 37K reads (500K pairs) |
| `trace_minimap2_fl.py` | Traces fragments with inflated FL |
| `trace_sj_gap_correction.py` | Quantifies gap SJ containment tolerance impact |
| `trace_gap_classification.py` | Classifies all gap types (near-miss, far-miss, no-SJ) |
| `compare_fl_distributions.py` | Compares FL stats from annotated BAMs |
| `analyze_pristine_benchmark.py` | Full benchmark analysis |
| `analyze_pristine_followup.py` | nRNA siphon, salmon bias, isoform redistribution |
