# gDNA Effective-Length Normalization Audit

Date: 2026-05-15

## Summary

Rigel currently normalizes RNA and gDNA lengths in two different places:

- mRNA and synthetic nRNA are normalized inside the EM using a transcript-level
  FL-marginal effective length `Ltilde_t`.
- gDNA is normalized before EM inside the scorer by subtracting a per-fragment,
  per-anchor sampling-window term from `gdna_log_lik`; the gDNA EM component
  then receives `log_eff_len = 0`.

So gDNA is not unnormalized. The confusing part is that it is normalized in a
different parameterization from RNA. That discrepancy is a code-health problem,
and it may also be a theoretical problem when comparing gDNA against nRNA in
ambiguous unstranded loci.

## Current Implementation

### RNA Candidates

The scorer writes RNA candidate likelihoods with fragment-length probability,
strand likelihood, overhang penalty, and mismatch penalty:

```text
log_lik_rna(f, t) = log h_R(flen_f) + log strand(f, t) + penalties
```

The native EM then applies transcript effective length every E-step and during
post-EM assignment:

```text
posterior(f, t) proportional to theta_t * exp(log_lik_rna(f, t)) / Ltilde_t
```

where `Ltilde_t = sum_l h_R(l) * max(length_t - l + 1, 0)`, floored at 1.0 for
log safety. For synthetic nRNA, `length_t` is the synthetic genomic span because
synthetic nRNA rows are single-exon transcript-like components.

Implementation points:

- `src/rigel/scoring.py` passes per-transcript effective lengths to the
  estimator.
- `src/rigel/native/em_solver.cpp` populates `sub.log_eff_len` from those
  effective lengths for transcript-like components.
- `src/rigel/native/em_solver.cpp` uses `log(theta) - log_eff_len` in both the
  EM E-step and assignment.

### gDNA Candidates

The scorer writes one gDNA likelihood per fragment unit. For unspliced units,
it computes a local sampling window using an anchor transcript span plus a gDNA
fragment-length flank:

```text
L_h = transcript_span(anchor_t) + gdna_flank
e_h(f) = max(L_h - genomic_footprint_f + 1, 1)

log_lik_gdna(f) = log h_G(genomic_footprint_f)
                + log gDNA_splice_penalty
                + log(0.5)
                + log mismatch_penalty
                - log e_h(f)
```

For multimapping units, the scorer accumulates these per-hit gDNA terms by
log-sum-exp and divides by the number of gDNA-eligible hits.

The EM then appends gDNA as a locus-level component but assigns:

```text
log_eff_len_gdna = 0
```

because the scorer has already baked in the gDNA length correction.

Implementation points:

- `src/rigel/native/scoring.cpp` computes `e_h` and subtracts `log(e_h)` from
  gDNA likelihoods.
- `src/rigel/native/em_solver.cpp` explicitly sets the gDNA component's
  `log_eff_len` to zero.

## Theoretical Model

There are two coherent ways to write fragment likelihoods.

### Conditional-On-Observed-Length Model

For a component with support length `L`, condition on the observed fragment
length `l_f`:

```text
P(f | component, l_f) proportional to q(f, component) / E_component(l_f)
```

where `E_component(l_f)` is the number of valid start positions for a fragment
with that specific length. This produces per-fragment denominators such as
`L - l_f + 1`.

This is a valid model if used consistently for all components. Rigel used to be
closer to this for RNA, but it was changed because the observed-length floor and
TPM mismatch caused short-transcript artifacts.

### FL-Marginal Effective-Length Model

For a component with support `S`, integrate possible start positions over the
fragment-length distribution:

```text
Ltilde_component = sum_l h(l) * N_valid_starts(S, l)

P(f | component) proportional to h(l_f) * q(f, component) / Ltilde_component
```

This is the salmon/kallisto-style contract now used by Rigel for RNA. It makes
the component denominator independent of the observed fragment and keeps the EM
denominator consistent with TPM/exposure calculations.

Because Rigel now uses the FL-marginal model for mRNA and nRNA, the theoretically
clean gDNA treatment is to use the same contract for gDNA:

```text
posterior(f, gDNA_locus) proportional to
    theta_g * h_G(length_f) * q_g(f) / Ltilde_gdna_locus
```

## What Should `Ltilde_gdna_locus` Be?

The gDNA component is not a transcript. Its support should be the accessible
genomic/capture support for the locus, not one arbitrary RNA transcript used as
an anchor for the current fragment.

In discrete form:

```text
Ltilde_gdna(locus) = sum_l h_G(l) * sum_start W_locus(start, l)
```

where `W_locus(start, l)` is an indicator or weight that a gDNA fragment of
length `l` starting at `start`:

1. is mappable/alignment-observable,
2. is capture-accessible or library-accessible,
3. overlaps the genomic support assigned to this locus in a way that would make
   it enter the locus EM, and
4. passes the same region/routing semantics used by scoring and calibration.

For ordinary RNA-seq without capture information, `W` can start as a simple
unweighted locus-footprint exposure. For hybrid-capture data, `W` should include
probe/bait accessibility or an empirical capture profile. For repetitive
intergenic sequence, `W` should include mappability. This is where mappability-
corrected effective length belongs.

The practical first approximation is:

```text
Ltilde_gdna(locus) = sum_l h_G(l) * N_valid_starts(locus_footprint, l)
```

with a clearly documented choice of `N_valid_starts`:

- contained support: starts where the whole fragment is inside the locus
  footprint;
- overlap support: starts where the fragment overlaps the locus footprint;
- capture/accessibility support: starts weighted by mappability/capture.

The current `transcript_span(anchor) + mean_gDNA_FL` heuristic is neither a pure
contained support nor a pure overlap support. It approximates a local overlap
window for typical fragment lengths, but it is fragment- and anchor-dependent.

## Is The Current Implementation Theoretically Discrepant?

Yes, relative to Rigel's current RNA EM contract.

The current gDNA likelihood is not missing length normalization, but it uses a
conditional, per-fragment denominator:

```text
h_G(l_f) / e_h(f)
```

while RNA uses an FL-marginal component denominator:

```text
h_R(l_f) / Ltilde_t
```

Those can both be valid in isolation, but mixing them makes gDNA and nRNA
posteriors depend on different abundance parameterizations. It also hides gDNA
normalization in the scorer while RNA normalization is visible in the EM.

The discrepancy is usually small for long loci because `e_h(f)` and
`Ltilde_gdna` are both roughly the locus span. It can matter for short loci,
capture windows, repeat/mappability correction, and edge cases where nRNA and
gDNA are otherwise nearly tied.

For the `gdna_zero_ss_0.50_nrna_high` failure, this discrepancy is unlikely to be
the sole driver: scaling the gDNA prior to zero still left large gDNA siphoning.
But the inconsistency makes the behavior harder to reason about and prevents us
from cleanly adding mappability/capture-corrected gDNA exposure.

## Recommended Fix

Move gDNA length normalization into the same EM interface used for RNA.

1. Compute a per-locus `Ltilde_gdna` in Python when building locus EM inputs.
2. Remove the `-log(e_h(f))` term from gDNA scoring, leaving only fragment-level
   likelihood terms: `log h_G(length_f)`, strand symmetry, splice penalties,
   mismatch penalties, and mapping penalties.
3. Pass `log(Ltilde_gdna)` into the native locus subproblem and set
   `sub.log_eff_len[gdna_idx] = log(Ltilde_gdna)`.
4. Keep gDNA eligibility separate from gDNA prior strength and effective length.
5. Add tests that assert RNA and gDNA both use exactly one length normalizer and
   that the normalizer changes gDNA/RNA posterior odds in the expected direction.

This is mostly an interface/semantics cleanup, but it also enables the more
important biological extensions: mappability-corrected and capture-corrected
gDNA effective lengths.

## Suggested Rollout

1. Add a `gdna_effective_lengths` array parallel to `gdna_prior_count` and
   `enable_gdna`.
2. Initially compute it from the current locus footprint with a simple documented
   approximation. Reproduce current behavior as closely as possible in aggregate.
3. Change scoring to emit unnormalized gDNA log likelihoods and change EM to
   consume the new gDNA effective length.
4. Compare against synthetic 24 and the focused unstranded high-nRNA condition.
5. Only after parity/behavior is understood, replace the approximation with
   mappability/capture-aware exposure.

## Code-Health Recommendation

Even if the numerical behavior barely changes, the code should not leave gDNA
normalization hidden in `gdna_log_lik` while RNA normalization lives in
`log_eff_len`. The current layout is too easy to misread as "gDNA has no
effective-length correction." A unified `log_eff_len` path for every component
would make the EM contract explicit:

```text
component posterior = exp(fragment_log_likelihood) * theta / effective_length
```

with one exception-free rule for mRNA, nRNA, and gDNA.