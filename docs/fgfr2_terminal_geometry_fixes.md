# FGFR2 Terminal Geometry Problem — Analysis & Fix

**Date:** 2026-02-22
**Status:** Implemented (warm start); evaluation pending benchmark re-run
**Related benchmark:** `bench_pristine_10_regions`, chr10:121476640–121601584

---

## Problem Summary

The FGFR2 locus (chr10) is a massive outlier: hulkrna MAE = 382.9 vs Salmon 98.4
vs Kallisto 42.3.  The #1 offender is **ENST00000684153.1** (truth = 13,
hulkrna estimate = 3,524, error = +3,511).  Fragments are siphoned from
**ENST00000682550.1** (truth = 26,428, error = −2,923).

### Root Cause: The Terminal Geometry Problem

The two transcripts are nearly identical (shared 14 of 15 splice junctions,
all exons overlap) and differ at two ends:

| End | ENST00000682550.1 (donor) | ENST00000684153.1 (siphon) |
|-----|---------------------------|----------------------------|
| 5′ genomic (exon 1) | Exon boundary at 3382 | Exon boundary at 3474 (+92 bp extension) |
| 3′ genomic (terminal exon) | SJ 117328→119837, terminal exon 112 bp | SJ 117328→121322, terminal exon 405 bp |

**The geometric reality:** A terminal exon of width W can only generate
uniquely-mapping fragments from W start positions (fragments must start
within the exon).  An identical internal exon generates W + frag_len − 1
start positions because fragments can also begin upstream.  For the
donor's 112 bp terminal exon with 300 bp fragments, the "runway" deficit
is ~3.5×.

**How this breaks the EM:**

1. The donor's unique evidence comes exclusively from its 112 bp terminal
   exon — only **85 unique reads** out of 26,428 truth fragments.

2. The siphon's 92 bp exon extension (3382–3474) generates **402 reads**
   that correctly cross annotated SJ 3474→7058 with small anchors (1–17 bp).
   These are valid alignments to an annotated SJ — not misalignments.

3. The EM initialized from `theta = unique_counts + prior`, giving:
   siphon = 402, donor = 85 — a **4.73× ratio** that should be ~0.0005×.
   This seeds a catastrophic local optimum.

---

## Fix: Warm-Start EM Initialization (Implemented)

### Concept

Replace the unique-count-only initialization with a **marginal
compatibility warm start**: each ambiguous fragment distributes
1/k weight across its k candidate transcripts before the first
E-step.  This adds a geometry-independent baseline proportional
to `eff_len × abundance`, eliminating the fragile dependence on
unique counts from tiny terminal exons.

### Why This Directly Fixes the Runway Problem

The "runway" issue means unique counts are a noisy, geometry-biased
estimate of abundance.  A 112 bp terminal exon yields only 85 uniques
despite 26,428 true fragments — the unique signal represents 0.3% of
the transcript.  This sparse signal is trivially overwhelmed by the
siphon's 402 reads.

The warm start adds ~2,000 marginal-compatibility counts to each
transcript.  This baseline is **independent of unique-region geometry**:
it depends only on how many total fragments are compatible with the
transcript, which is proportional to effective_length × abundance.
The 85-vs-402 unique disparity becomes 3.8% noise on a stable ~2,000
baseline.

### Quantified Impact

| Metric | Before | After Warm Start |
|--------|--------|-----------------|
| Donor init (682550) | 85 | ~2,221 |
| Siphon init (684153) | 402 | ~2,419 |
| Eff. abundance ratio (siphon/donor) | **4.29×** | **0.99×** |

With near-parity initialization, the EM converges based on
likelihoods (strand, insert size, SJ matching) rather than on
the fragile unique-count seed.

### Implementation

In `counter.py`, `run_locus_em()`, the initialization was changed from:

```python
# OLD: fragile unique-count initialization
theta_init = unique_totals.copy()
# + shadow floor logic
```

To:

```python
# NEW: warm start with marginal compatibility
theta_init = unique_totals.copy()
if len(seg_lengths) > 0 and len(t_indices) > 0:
    cand_weights = 1.0 / np.repeat(
        seg_lengths.astype(np.float64), seg_lengths,
    )
    np.add.at(theta_init, t_indices, cand_weights)
```

The old "balanced initialization" (shadow floor) was removed — the
warm start subsumes it by providing a robust baseline for all
components.

### Test Results

All 675 tests pass including all 51 counter tests and all 213
scenario tests.

---

## Alternative Fixes Evaluated

### Rejected: SJ Anchor Filter

**Proposal:** Filter splice junctions with CIGAR anchors below a
minimum threshold (e.g., ≤8 bp) in `parse_read()`.

**Why rejected:** The 402 reads at the siphon's SJ are correct
alignments to annotated splice junctions, not misalignments.
Genome browser inspection confirms valid alignments.  Filtering
would reduce sensitivity globally — throwing away reads that
correctly cross annotated SJs — for a local gain at one locus.
For annotated splice junctions, tiny anchors are expected and
valid because the splice site is pre-validated.

### Rejected: Naive Runway Correction (Inflate Unique Counts)

**Proposal:** Compute each transcript's unique region fraction
and inflate unique counts by 1/fraction.

**Why rejected:** Both transcripts have small unique regions.
Inflating by the inverse fraction amplifies the siphon's false
signal more than the donor's true signal:

| | Unique | Unique fraction | Corrected |
|---|---|---|---|
| Donor | 85 | 0.031 | 2,723 |
| Siphon | 402 | 0.023 | 17,269 |

The corrected ratio worsens from 4.7× to **6.3×**.

### Future: Positional Bias Model

**Concept:** Model position-dependent fragment sampling along
each transcript to give terminal-SJ-crossing fragments higher
discriminative weight.

**Status:** Deferred.  Most principled long-term improvement but
high implementation complexity (~200+ LOC, 4-5 files).  The warm
start eliminates the grossest initialization artifacts; the
positional model would further refine per-fragment likelihoods.

---

## Appendix: Key Data Points

### Fragment Resolution Statistics

- Total fragments: 50,000
- Fragments with siphon as candidate: 40,343 (80.7%)
- Fragments *uniquely* resolved to siphon: 402
  - ALL from SJ 3474→7058 (small anchors 1–17 bp at annotated SJ)
  - ALL truth origin: other transcripts (235 from donor, 129 from 369061)
  - These are valid alignments to an annotated SJ, not misalignments
- Fragments uniquely resolved to donor: 85
  - ALL from terminal SJ 117328→119837 (112 bp terminal exon)

### EM Initialization Dynamics

| Scenario | Siphon theta_init | Donor theta_init | Eff. ratio |
|----------|-------------------|------------------|------------|
| Before (unique only) | 402 | 85 | 4.29× siphon |
| After (warm start) | ~2,419 | ~2,221 | 0.99× (near parity) |
