# Capacity-Weighted EM with Per-Fragment Effective Lengths

## TL;DR

The EM currently ignores where a fragment falls on a transcript, causing
systematic bias when isoforms share a large exon but differ at a small
terminal exon.  We fix this in two complementary ways:

1. **Per-fragment effective lengths** in the E-step (Bayesian-correct
   replacement for the current global effective length).
2. **Capacity-weighted warm start** that distributes ALL fragments
   (unique and ambiguous) according to geometric difficulty, eliminating
   reliance on unique reads.

The trapezoid coverage model defines "difficulty": fragments near
transcript edges have low coverage capacity, so observing them is
stronger evidence for that transcript.  This resolves ambiguity even
when zero uniquely-mapping fragments exist.

---

## Background

### The problem

Benchmarking revealed that the largest estimation errors come from
transcripts with ~99% identical structures — same exons and introns
except for a different 5' start.  When the 5' start exon is smaller
than the fragment size, the transcript with the tiny first exon is
systematically undercounted and a competing transcript with a larger
"runway" is systematically overcounted.

Example: T1 exons `[(1000,1050), (2000,22000)]`, T2 exons
`[(1200,1500), (2000,22000)]`, abundance ratio 10:1.  T2 is
overcounted by **53%** because only ~3/1000 fragments span the splice
junction — the EM has almost no discriminating signal for the 997
fragments landing entirely in the shared 20 kb exon.

### Why the standard EM fails

The EM's per-fragment posterior is:

$$\pi(t \mid f) \propto \frac{\theta_t}{\text{eff\_len}_t} \times L_{ft}$$

where $\text{eff\_len}_t$ is a **global** per-transcript constant.
When two transcripts share a massive exon their effective lengths
are nearly identical (e.g. 19,701 vs 19,951 — only 1.3% difference).
The MLE converges to a ratio driven almost entirely by the few unique
reads, ignoring that T1 squeezed its unique reads from a 50 bp
bottleneck while T2 had 300 bp of runway.

### Why k-mer tools do better

Salmon and kallisto operate at the k-mer level (k ≈ 31).  More k-mers
can fit on a tiny first exon, giving them finer positional resolution.
hulkrna operates at the fragment level and must compensate via its
likelihood model.

---

## Phase 1: Per-Fragment Effective Length Correction

Replace the global `eff_len_t` in the E-step with a per-fragment value
`(L_t - l_{f,t} + 1)` for each candidate.  This is the Bayesian-correct
denominator — see Roberts & Pachter (2011).

### Design rationale

Per-fragment effective length is an additive term in the log-likelihood:

    log P(f | t) = log_strand + log_insert + oh × penalty + log_nm
                   - log(max(L_t - l_f + 1, 1))

For mRNA and nRNA, all inputs are available during the scan phase, so
the correction is baked directly into `log_lik` at scoring time.  For
gDNA, the locus span is not known until `build_locus_em_data()`, so the
correction is applied when gDNA candidates are constructed.  This avoids
adding per-candidate arrays to `ScanData` and keeps each correction in
the code path that naturally has the data.

The three effective-length coordinate systems:

- **mRNA**: `L_t` = transcript spliced exonic length, `l_f` = per-
  candidate SJ-corrected fragment length (varies per transcript).
- **nRNA**: `L_t` = transcript genomic span (start to end including
  introns), `l_f` = fragment genomic footprint.
- **gDNA**: `L` = true locus span (see below), `l_f` = fragment
  genomic footprint.

### Locus span for gDNA

The locus span is the full DNA region "in play":

    locus_start = min(min(transcript_starts), min(fragment_starts))
    locus_end   = max(max(transcript_ends),   max(fragment_ends))
    locus_span  = locus_end - locus_start

Using transcript bounds (not gene bounds) and extending to cover any
overhanging fragments ensures the span reflects the actual observed
footprint.  This requires `genomic_start` and `genomic_footprints`
per unit in `ScanData`.

### Implementation

**Step 1a.** Add `t_length_arr` (spliced exonic lengths, int32) and
`t_span_arr` (genomic spans, int32) to `ScoringContext`.  Set both in
`from_models()`.

**Step 1b.** In `_score_wta_mrna()`: after computing `log_lik`, subtract
`log(max(L_t - frag_len + 1, 1))` where `L_t = ctx.t_length_arr[t_idx]`
and `frag_len = frag_lengths[k]`.

**Step 1c.** In `_score_wta_nrna()`: after computing `nrna_ll`, subtract
`log(max(span_t - genomic_footprint + 1, 1))` where
`span_t = ctx.t_span_arr[t_idx]`.

**Step 1d.** Add `genomic_starts` (int32[n_units]) and `genomic_footprints`
(int32[n_units]) to `ScanData`.  Populate during the scan from
`bf.genomic_start` and `bf.genomic_footprint`.

**Step 1e.** In `build_locus_em_data()`: compute the true locus span
from transcript boundaries and fragment boundaries.  When constructing
gDNA candidates, subtract `log(max(locus_span - footprint + 1, 1))` from
each gDNA log-likelihood.

**Step 1f.** Set all effective lengths to 1.0 in `build_locus_em_data`:
`eff_len[:] = 1.0`.  The per-fragment correction already handles length
normalization.

**Step 1g.** Add `effective_length` column to the transcript-level output
DataFrame (`get_counts_df`).  Add `effective_length` to the gene-level
output DataFrame (`get_gene_counts_df`).

### Output

The reported `effective_length` in the output is the **global
analytical value** (same formula salmon uses), NOT the per-fragment
value.  The per-fragment correction changes the EM's internal
weighting but downstream TPM should use the classical formula:

    TPM_t = (count_t / eff_len_t) / Σ_k(count_k / eff_len_k) × 10^6

For gene-level output, the effective length is the analytical
effective length of the gene's merged exonic footprint.

---

## Phase 2: Capacity-Weighted Warm Start

Instead of splitting ambiguous fragments uniformly (1/k) during EM
initialization, distribute them proportional to the **inverse coverage
capacity** at each fragment's position on each candidate transcript.

### The trapezoid model

For a transcript of spliced length L and mean insert size μ:

- `w_max = min(μ, L/2)` — maximum coverage depth in the plateau
- Left ramp: positions `[0, w_max)` — coverage rises linearly
- Plateau: positions `[w_max, L - w_max)` — coverage = w_max
- Right ramp: positions `[L - w_max, L)` — coverage falls linearly

The **fragment interval weight** is `w_max / mean_capacity` where
`mean_capacity` is the integral of the trapezoid over the fragment's
span divided by the fragment length.  Fragments deep in the plateau
get weight ≈ 1.0; fragments at the edges get weight > 1.0.

### Why this works without unique reads

Consider a gene with three isoforms where no fragment maps uniquely:
- T1: exons (1000,2000), (3000,4000), (5000,6000) — L=3000
- T2: exons (1000,2000), (3000,3500) — L=1500
- T3: exons (3500,4000), (5000,6000) — L=1500

A fragment at [3000,3400] is compatible with T1 and T2.  For T1 it
lands in the middle (weight ≈ 1.0).  For T2 it's near the 3' edge
(weight ≈ 4.0).  The capacity-weighted split gives T2 80% of this
fragment's initialization weight — correctly reflecting that T2 "had
to work harder" to generate it.

### Implementation

1. Cache per-transcript exon intervals on `HulkIndex`.
2. Store `genomic_start` per fragment in the buffer.
3. Implement `genomic_to_transcript_pos()`.
4. Implement `compute_fragment_interval_weight()`.
5. Compute capacity weights during the scan phase.
6. Thread weights through `ScanData` and `LocusEMInput`.
7. Replace uniform 1/k warm start with capacity-weighted distribution.
8. Optional: capacity-weighted persistent Dirichlet prior (gamma parameter).

---

## Phase 3: Regression Tests

- **Scenario A**: T1/T2 with tiny vs large first exon, 10:1 abundance.
  Assert T2 relative error < 20% (currently 53%).
- **Scenario B**: T1/T2/T3 with zero unique reads (all-ambiguous).
  Assert reasonable accuracy for edge-bounded transcripts.
- **Scenario C**: Single-exon control — verify no regression.

---

## Decisions

- **Capacity weight as initialization, not E-step penalty**: Adding a
  position penalty to the log-likelihood penalizes the transcript you
  are trying to rescue.  The weight enters via warm start and optional
  prior, keeping the E-step generatively pure.
- **Per-fragment effective length is separate and complementary**:
  Addresses the mathematical correctness of the denominator.  Small
  impact alone (~1.3% for the pathological case) but correct in general.
- **No unique-read dependency**: The capacity weight is computed for
  every fragment-candidate pair, so it works for complex genes with
  15+ isoforms and zero uniquely-mapping fragments.
- **Output effective_length is the classical analytical value**: The
  per-fragment correction changes internal EM weighting only.
  Downstream TPM uses the standard salmon-style effective length.