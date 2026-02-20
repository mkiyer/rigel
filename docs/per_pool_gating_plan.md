# Per-Pool Compatibility Gating Plan

## Status: Planned (not yet implemented)

This document describes the planned per-pool binary compatibility gating
that will be implemented **after** the exponential overhang penalty is in
place. The overhang penalty makes explicit gating less critical (large
overhang → near-zero likelihood), but gating still provides value as a
hard floor that prevents hopeless candidates from entering the EM at all,
reducing EM problem size and eliminating numerical noise.

---

## Background

### Current behaviour (pre-overhang-penalty)

`resolve_fragment()` returns **all** transcripts with any genic overlap
(exon or intron). For spliced fragments, `filter_by_overlap()` prunes
candidates whose exon overlap fraction is < 99 % of the best. For
unspliced fragments, the filter is skipped entirely — the original intent
was to avoid removing valid nRNA/gDNA candidates whose overlap is
predominantly intronic.

The problem: this means every unspliced fragment carries the full set of
overlapping transcripts into the EM, including transcripts where only 1–2
exonic bases overlap. The overlap-exponent scoring was supposed to
down-weight these, but it produces near-zero discrimination when
candidates have similar overlap fractions.

### After the overhang penalty

With the exponential overhang penalty (`b_out × log(α)`), a candidate
with even 5 bases of overhang gets a penalty of `5 × log(0.01) ≈ −23`.
This naturally drives poor candidates to near-zero posterior probability.
Gating adds a clean, deterministic floor on top of this.

---

## Design

### Principle: gate per pool, not per fragment

The key insight is that a fragment's compatibility must be evaluated
**separately** for each pool (mRNA, nRNA) because the target boundary
differs:

| Pool | Target boundary | Compatible when |
|------|----------------|-----------------|
| mRNA | Transcript exons | Fragment falls within exon footprint |
| nRNA | Gene body (exons + introns) | Fragment falls within gene span |
| gDNA | Genome | Always compatible (no gating) |

A fragment that overhangs a transcript's exon boundary into its intron is
**incompatible** as mRNA but **compatible** as nRNA. The current system
conflates these two checks.

### Gating rule

During EM data construction (`_add_transcript_candidates`,
`_add_nrna_candidates`, and `_flush_mm_group`), before scoring a
candidate, apply:

**mRNA candidate T for fragment f:**
```
overhang = frag_length − exon_bp(f, T)
if overhang > max_overhang_bp:
    skip this mRNA candidate
```

**nRNA candidate T for fragment f:**
```
overhang = frag_length − (exon_bp(f, T) + intron_bp(f, T))
if overhang > max_overhang_bp:
    skip this nRNA candidate
```

### Threshold derivation (no new parameters)

Rather than introducing separate mRNA/nRNA thresholds, derive the gate
from the existing overhang penalty parameter `α`:

```
max_overhang_bp = floor(−effective_floor / log(α))
```

Where `effective_floor` is the log-likelihood floor below which a
candidate cannot meaningfully contribute (e.g., −50, meaning probability
< 2 × 10⁻²²). With `α = 0.01`:

```
max_overhang_bp = floor(50 / 4.605) = 10
```

This means: any candidate whose overhang penalty alone would push it below
10⁻²² is gated out before entering the EM. The gate is purely an
optimization — it removes candidates that the overhang penalty would
already make negligible.

**No new user-facing parameters.** The gate threshold is derived from `α`.

---

## Implementation steps

### 1. Add gating to `_add_transcript_candidates()`

```python
def _add_transcript_candidates(bf) -> tuple[float, int]:
    fl = bf.frag_length if bf.frag_length > 0 else 1
    for k, t_idx in enumerate(bf.t_inds):
        ebp = int(bf.exon_bp[k]) if bf.exon_bp is not None else fl
        b_out = max(fl - ebp, 0)
        if b_out > max_overhang_bp:      # ← NEW: mRNA gate
            continue
        # ... score and append
```

### 2. Add gating to `_add_nrna_candidates()`

```python
def _add_nrna_candidates(bf) -> None:
    fl = bf.frag_length if bf.frag_length > 0 else 1
    for k, t_idx in enumerate(bf.t_inds):
        ebp = int(bf.exon_bp[k]) if bf.exon_bp is not None else fl
        ibp = int(bf.intron_bp[k]) if bf.intron_bp is not None else 0
        b_out_nrna = max(fl - ebp - ibp, 0)
        if b_out_nrna > max_overhang_bp:  # ← NEW: nRNA gate
            continue
        # ... score and append
```

### 3. Same changes in `_flush_mm_group()`

Both the mRNA and nRNA branches of the multimapper flush function get
identical gates.

### 4. Safety: ensure at least one candidate survives

If gating removes ALL mRNA candidates for a fragment, the fragment
effectively has no mRNA hypothesis — it competes only in nRNA/gDNA pools.
This is correct behaviour (a fragment with massive overhang on all
transcripts is likely nRNA or gDNA).

If gating removes ALL candidates across all pools, the fragment is
unassignable. Track this in stats as a new counter.

---

## Impact on the five robustness scenarios

| Scenario | mRNA gate effect | nRNA gate effect |
|----------|-----------------|-----------------|
| Single-gene multi-isoform (EGFR) | Removes isoforms whose exons don't cover fragment position | All isoforms pass (same gene body) |
| Intronic overlap | Removes transcripts from mRNA (correct) | Retains for nRNA (correct) |
| Multi-gene same-strand | Each gene's transcripts gated independently | Each gene's transcripts gated independently |
| Antisense overlap | Gene A exon / Gene B intron: B removed from mRNA, kept in nRNA | Strand model discriminates further |
| Chimeric transcripts | Gating operates per-candidate; chimera classification unchanged | Same |

---

## Testing plan

1. Unit test: verify mRNA gate removes candidate with overhang > threshold
2. Unit test: verify nRNA gate uses genic (exon+intron) boundary
3. Unit test: verify fragment with all candidates gated gets no mRNA entries
4. Integration: EGFR benchmark shows reduced candidate set sizes
5. Regression: all existing scenario tests still pass

---

## Dependencies

- **Prerequisite:** exponential overhang penalty must be implemented first
  (the gate threshold is derived from `α`)
- **No changes to:** `resolve_fragment()`, `filter_by_overlap()`,
  `compute_overlap_profile()`, buffer storage, or the EM solver itself
