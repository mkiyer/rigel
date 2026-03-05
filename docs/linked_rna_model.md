# Linked mRNA–nRNA Model

## Summary

This document describes the **linked RNA model** for jointly estimating
mRNA (mature) and nRNA (nascent/pre-mRNA) transcript abundances within
hulkrna's locus-level EM framework.  The model couples mRNA and nRNA
through a shared per-transcript abundance parameter θ\_t with a
per-transcript nascent fraction nrna_frac\_t, grounded in a kinetic model of
transcription, splicing, and degradation.

## Motivation

In eukaryotic cells, most genes produce overlapping mRNA and nRNA
(pre-mRNA) molecules.  A total RNA-seq library captures both species.
The current hulkrna architecture already models mRNA and nRNA as
separate EM components per transcript, plus a single gDNA component per
locus.  However, the mRNA and nRNA components are treated as fully
independent, which:

1. **Doubles the effective parameter count** — every transcript
   contributes two independent abundance parameters.
2. **Creates massive loci** — overlapping nRNA transcripts from
   neighbouring genes merge into giant EM sub-problems.
3. **Ignores biological coupling** — observing mRNA from a transcript
   *increases* the probability of observing nRNA from the same
   transcript, and vice versa.

The linked model addresses all three issues by tying mRNA and nRNA
through a single abundance parameter per transcript.

## Kinetic Model

The biological basis is the RNA lifecycle:

```
  transcription         splicing         degradation
       ↓                  ↓                  ↓
  DNA ──→ nRNA (pre-mRNA) ──→ mRNA (mature) ──→ ∅
          rate = r_txn     rate = r_splice    rate = r_deg
```

At steady state, the ratio of nRNA to total RNA per transcript is:

```
  nrna_frac_t = nRNA_t / (mRNA_t + nRNA_t) = 1 / (1 + ρ_t)

  where  ρ_t = r_splice / r_deg  (splicing-to-degradation rate ratio)
```

Key properties:
- **Fast splicing** (ρ ≫ 1): nrna_frac → 0, most RNA is mature mRNA
- **Slow splicing** (ρ ≪ 1): nrna_frac → 1, most RNA is nascent nRNA
- nrna_frac is **per-transcript**, allowing different isoforms of the same gene
  to have different splicing kinetics

## Parameters

Per locus with T transcripts:

| Parameter | Count | Domain | Description |
|-----------|-------|--------|-------------|
| θ\_t | T | simplex (Σ θ + γ = 1) | Total (mRNA + nRNA) abundance fraction |
| nrna_frac\_t | T | \[ε, 1−ε\] | Nascent RNA fraction per transcript |
| γ\_ℓ | 1 | \[0, 1\] | gDNA fraction for the locus |

**Total parameters per locus: 2T + 1** (same as current architecture).

Derived quantities:
- `mRNA_t = θ_t × (1 − nrna_frac_t)` — mature mRNA fraction
- `nRNA_t = θ_t × nrna_frac_t` — nascent nRNA fraction

## Fragment Likelihood

For fragment f assigned to transcript t:

```
  w(f, t, mRNA) = θ_t × (1 − nrna_frac_t) × L_exonic(f, t)
  w(f, t, nRNA) = θ_t × nrna_frac_t       × L_genomic(f, t)
  w(f, gDNA)    = γ_ℓ              × L_uniform(f)
```

where:
- `L_exonic(f, t)` = exonic alignment likelihood (spliced/unspliced to exons)
- `L_genomic(f, t)` = genomic alignment likelihood (pre-mRNA span)
- `L_uniform(f)` = uniform gDNA likelihood across locus span

In log-space (implemented in the C++ EM kernel):

```
  log_weight[mRNA component i] = log(θ_t) + log(1 − nrna_frac_t) − log(eff_len_mRNA_i)
  log_weight[nRNA component i] = log(θ_t) + log(nrna_frac_t)     − log(eff_len_nRNA_i)
  log_weight[gDNA component]   = log(γ_ℓ)                 − log(eff_len_gDNA)
```

## E-Step and M-Step

### E-Step (unchanged structure)

For each fragment f, compute posterior responsibility:

```
  r(f, c) = w(f, c) / Σ_c' w(f, c')
```

where c ranges over all candidate components (mRNA, nRNA, gDNA).

### M-Step

**θ update** — total abundance per transcript:

```
  θ_t ∝ unambig_totals_t + Σ_f [r(f, t, mRNA) + r(f, t, nRNA)] + prior_t
```

Normalized so that `Σ_t θ_t + γ_ℓ = 1`.

**nrna_frac update** — nascent fraction with Beta prior (MAP estimate):

```
  nrna_assigned_t = Σ_f r(f, t, nRNA)
  total_assigned_t = Σ_f [r(f, t, mRNA) + r(f, t, nRNA)]

  nrna_frac_t = (nrna_assigned_t + α_t − 1) / (total_assigned_t + α_t + β_t − 2)
```

Clamped to `[ε, 1−ε]` where `ε = 1e-8`.

**γ update** — gDNA fraction:

```
  γ_ℓ ∝ gdna_init + Σ_f r(f, gDNA)
```

Normalized jointly with θ\_t.

## SQUAREM Integration

The SQUAREM acceleration operates on the state vector `θ` (the simplex
of total abundances for all components including gDNA).  **nrna_frac is NOT part
of the SQUAREM state vector** — it is derived deterministically from the
M-step ratio after each EM iteration.  This means:

- No boundary violations for nrna_frac from SQUAREM extrapolation
- The existing SQUAREM infrastructure requires no structural changes
- nrna_frac is simply clamped to `[ε, 1−ε]` in the M-step

## Prior Hierarchy for nrna_frac

nrna_frac\_t is estimated per-transcript with a Beta(α\_t, β\_t) prior.  The
prior parameters are computed **before EM** from the observed
spliced/unspliced exonic fragment counts, using a cascading hierarchy:

### Hierarchy Levels

1. **Transcript level** — Use the transcript's own unambig counts
2. **TSS group** — Transcripts sharing the same transcription start
   site (fuzzy 5' clustering, see below)
3. **Locus-strand group** — All transcripts on the same strand within
   the same locus
4. **Weak prior** — nrna_frac = 0.5, κ = 2 (uninformative)

### Cascade Logic

```python
MIN_EVIDENCE = 10  # minimum fragments to trust an estimate

for each transcript t:
    # Level 1: transcript's own data
    spliced_t = unique_spliced_exonic[t]
    unspliced_t = unique_unspliced_exonic[t]
    total_t = spliced_t + unspliced_t

    if total_t >= MIN_EVIDENCE:
        nrna_frac_est = unspliced_t / total_t
        κ = total_t
    # Level 2: TSS group
    elif tss_group_total[t] >= MIN_EVIDENCE:
        nrna_frac_est = tss_group_unspliced[t] / tss_group_total[t]
        κ = tss_group_total[t]
    # Level 3: locus-strand
    elif locus_strand_total[t] >= MIN_EVIDENCE:
        nrna_frac_est = locus_strand_unspliced[t] / locus_strand_total[t]
        κ = locus_strand_total[t]
    # Level 4: weak uninformative prior
    else:
        nrna_frac_est = 0.5
        κ = 2.0

    α_t = nrna_frac_est × κ
    β_t = (1 − nrna_frac_est) × κ
```

### Fuzzy TSS Grouping

Transcripts are grouped by approximate transcription start site using
**single-linkage clustering** of 5' positions within a configurable
window:

- **5' position**: `start` for + strand, `end` for − strand
- **Algorithm**: Per (chromosome, strand), sort transcripts by 5'
  position, then merge adjacent transcripts whose 5' positions differ
  by ≤ `tss_window` (default 200 bp)
- **Rationale**: Annotation fragmentation and alternative TSS usage
  creates multiple transcripts with nearly identical 5' ends.  Exact
  coordinate matching would fail to group these biologically equivalent
  start sites.

The TSS window is an advanced user-configurable parameter
(`EMConfig.tss_window`, default 200).

## Implementation Plan

### Phase 1: Data Structures and Priors (Python)

**1a. Fuzzy TSS grouping**
- Add `tss_window` parameter to `EMConfig` (default 200)
- Implement `compute_tss_groups()` function in `index.py`
  - Per (chrom, strand): sort by 5' position, single-linkage cluster
  - Returns `int32[num_transcripts]` array of group IDs
- Call from `TranscriptIndex.load()` or pipeline init

**1b. nrna_frac prior cascade**
- Add `eta_alpha` and `eta_beta` arrays to `AbundanceEstimator`
- After the scan phase (when unambig counts are available), compute
  per-transcript spliced/unspliced exonic counts
- Apply the 4-level cascade: transcript → TSS group → locus-strand →
  weak prior
- Pass nrna_frac priors into `build_locus_em_data` → `LocusEMInput`

### Phase 2: Linked Likelihood (C++ EM Solver)

- Modify `em_step_kernel` log\_weights to use linked formulation
- M-step: θ\_t = sum of mRNA+nRNA responsibilities; nrna_frac\_t from MAP ratio
- nrna_frac stored as a per-locus array alongside θ, updated each iteration
- No new CSR columns — mRNA/nRNA are already separate edges

### Phase 3: SQUAREM Integration

- θ remains in SQUAREM state vector (no changes)
- nrna_frac derived from M-step after each fixed-point map (not extrapolated)
- Clamp nrna_frac to \[ε, 1−ε\] in M-step

### Phase 4: Output Integration

- Report nrna_frac\_t in `quant.feather` and `quant_detail.feather`
- Update `loci.feather` with nrna_frac statistics
- Update `summary.json` with nrna_frac distribution metrics

## Design Decisions

### Why per-transcript nrna_frac (not per-gene)?

Different isoforms of the same gene can have different splicing rates
due to different intron structures, regulatory elements, or alternative
polyadenylation.  Per-transcript nrna_frac captures this biological variation.

### Why a hierarchy instead of a global estimate?

Global estimates (e.g., genome-wide average nrna_frac) are biologically
irrelevant because splicing kinetics vary enormously across genes and
cell types.  The hierarchy provides local estimates that fall back
gracefully when data is sparse.

### Why is nrna_frac NOT in the SQUAREM state vector?

SQUAREM accelerates convergence by extrapolating the state vector.
Extrapolating nrna_frac could violate the \[0, 1\] boundary.  Since nrna_frac is
deterministically derived from the θ M-step, it converges naturally
as θ converges.  This avoids the need for projection operators or
constrained SQUAREM variants.

### Why fuzzy TSS grouping?

Exact coordinate matching for TSS grouping fails when:
- Multiple annotations exist for the same biological TSS with slightly
  different 5' coordinates
- Alternative TSS usage produces starts within a few hundred base pairs
- Annotation databases disagree on exact 5' boundaries

Single-linkage clustering with a 200 bp window captures these cases
while remaining conservative enough to avoid merging truly distinct
start sites.
