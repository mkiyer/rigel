# gDNA Regional Exposure Revised Plan v3a

**Status**: design plan, not implemented.  
**Date**: 2026-05-19  
**Supersedes**: `gdna_exposure_revised_plan_v3.md`,
`gdna_exposure_revised_plan_v2.md`, and `gdna_exposure_plan_v3.md`.

## 0. Decision

Regional exposure is an assay sampling-opportunity correction, not
per-fragment evidence for molecule identity. Production EM should not mutate
fragment log-likelihood numerators with `log A(x_u)`. Instead, every component
that competes in EM receives an exposure-weighted effective length internally:

```text
L_k -> L_k^A
```

This applies to:

- the per-multilocus gDNA component;
- synthetic nRNA span components;
- annotated mRNA transcript components.

The prior used by the gDNA component must be rescaled with the same gDNA
effective-length ratio:

```text
eta_g_for_em = eta_g_unweighted * L_g^A / L_g
```

The current v3-production numerator path is the thing to undo. The numerator
was harmful when applied only to gDNA, and a symmetric numerator implementation
would add complexity without helping the production EM.

These weighted effective lengths are **EM-only denominators**. Public abundance
outputs keep raw, unweighted effective lengths and raw TPM semantics so samples
remain directly comparable across different capture designs and exposure
fields.

## 1. Model

For component `k`, let `p_k(u)` be the existing Rigel fragment likelihood for
observed unit `u`, excluding exposure. Let `A(x)` be the relative regional
capture exposure field. The capture-aware density can be written as:

```text
f_k(u) = p_k(u) * A(u) / L_k^A
```

where:

```text
L_k^A = integral p_k(s) A(s) ds
```

If the exposure numerator is common to all candidates for one observed unit,
it cancels from the EM posterior:

```text
P(k | u) = theta_k * p_k(u) / L_k^A
         / sum_j theta_j * p_j(u) / L_j^A
```

That cancellation is the important practical rule. In real data, candidate
footprints are not always identical: multimappers can imply different genomic
positions, and spliced mRNA candidates project through exons rather than one
contiguous genomic interval. Rigel's current regional field is not a trusted
candidate-level capture-likelihood model for those cases. Therefore v3a uses
only the identifiable denominator term in production.

## 2. Undo From Prior Implementation

Remove exposure numerator weighting from the production path:

- Do not call `_apply_unit_gdna_weights()` from `quant_from_buffer()` when
  `--regional-exposure auto` is active.
- Do not mutate `em_data.gdna_log_liks` with `log A(midpoint)`.
- Do not add v2-style per-candidate nRNA or mRNA numerator appliers.
- Remove the production need for `safe_unit`, cross-reference skips, and
  multimapper numerator skips.
- Keep `ScoredFragments.genomic_midpoint` only as a diagnostic/debug field.

Delete the old production path cleanly. Do not keep hidden compatibility modes
unless a separate benchmarking script needs to reproduce old results outside
the main pipeline.

## 3. Exposure Field

Reuse the existing `RegionalGdnaExposure` concept, but fix normalization before
using the field for any effective-length computation.

The v3-production bug was global cross-class normalization:

```text
rho_ref = max(rho_ref_per_class.values())
```

This made classes with different opportunity definitions, especially intronic
contained fragments versus exon-intron boundary crossings, look artificially
depleted. v3a uses a class-specific reference quantile:

```python
for cname in CLASSES:
    rho_ref_c = weighted_quantile(rho_hat[class == cname], opportunity[class == cname])
    raw = log(rho_hat_r) - log(rho_ref_c)
    log_A_r = signal_c * min(raw, 0.0)
```

Keep the old global `rho_ref` only as a summary diagnostic. It must never be
used to compute `A_r`.

Important limitation: class-specific normalization makes each class internally
consistent, but it does not prove that exon, intron, boundary, and intergenic
scales are mutually calibrated. If v3a fails the VCaP gates, the first suspect
is the exposure field construction, not the EM denominator machinery.

## 4. Weighted Effective Lengths

### 4.1 gDNA

Reuse `weighted_gdna_eff_len_for_loci()` from
`src/rigel/calibration/_exposure.py`.

For each multilocus, retain both values:

```text
L_g              = gdna_eff_len_unweighted
L_g^A            = gdna_eff_len
gdna_weight_ratio = L_g^A / L_g
```

`--regional-exposure off` and uniform exposure must delegate to the existing
unweighted path and remain bit-exact.

### 4.2 Synthetic nRNA

Add `weighted_nrna_eff_lens()` in `src/rigel/calibration/_exposure.py`.

Inputs:

- unique synthetic nRNA spans `(ref_id, strand, start, end, span_id)`;
- the RNA fragment-length model;
- the regional exposure field.

Algorithm:

- derive unique nRNA span IDs at `TranscriptIndex.load()` from existing
  synthetic transcript rows, with no index format bump;
- for each span and each positive fragment length, compute valid start
  windows;
- convert those start windows to the same coordinate convention used by the
  existing gDNA weighted integration;
- integrate `A` with `exposure.weighted_length_on_ref(...)`;
- clamp with the same minimum effective-length policy used elsewhere.

Only synthetic nRNA rows receive this span-based length. Annotated single-exon
or nascent-equivalent rows stay on the annotated transcript path below.

### 4.3 Annotated mRNA

Use a transcript-level exon exposure average. This is simple, concrete, and
matches the biological interpretation: mature mRNA opportunity lives on exons,
not on the genomic span between them.

Add `transcript_exposure_weights()` in `src/rigel/calibration/_exposure.py`:

```python
def transcript_exposure_weights(
    exon_intervals_by_t: Mapping[int, np.ndarray],
    t_to_ref_id: np.ndarray,
    exposure: RegionalGdnaExposure,
    *,
    min_weight: float = 1.0e-4,
) -> np.ndarray:
    """Return Abar_t for annotated transcript rows."""
```

For transcript `t`:

```text
Abar_t = sum_e sum_r overlap_bp(e, r) * A_r / sum_e exon_bp(e)
L_t^A  = max(Abar_t * L_t, min_effective_length)
```

Implementation details:

- use `index._t_exon_intervals` or an equivalent loaded exon-block table;
- map transcript reference with `index.t_to_ref_arr`;
- for each exon, call
  `regional_exposure.weighted_length_on_ref(ref_id, exon_start, exon_end)`;
- do not assume an exon maps to exactly one region, because region boundaries
  can split exon sequence when annotation contexts change;
- `weighted_length_on_ref()` already handles split regions, missing-reference
  identity behavior, and exact overlap lengths;
- divide by transcript exonic length, not genomic span;
- clamp `Abar_t` to `[min_weight, 1.0]`.

Then compute mRNA EM effective lengths by scaling the existing raw effective
lengths:

```text
mrna_eff_len_for_em[t] = max(raw_eff_len[t] * Abar_t, min_effective_length)
```

This deliberately avoids fragment-placement projection in v3a. If later data
show that exon-position effects within long transcripts matter, an exact
fragment-length-aware mRNA integration can be added as a refinement. It is not
needed for the first clean implementation.

### 4.4 EM Hand-Off

Build both EM and output transcript effective-length arrays before constructing
`AbundanceEstimator`:

```python
t_eff_len_output = current_effective_lengths(...)
t_eff_len_for_em = t_eff_len_output.copy()

t_eff_len_for_em[synthetic] = nrna_weighted_by_span[t_to_nrna_span_id[synthetic]]
t_eff_len_for_em[annotated] = t_eff_len_output[annotated] * Abar_t[annotated]
```

Use `t_eff_len_for_em` for native EM. Retain `t_eff_len_output` for public
`effective_length` and TPM outputs. Add explicit diagnostic fields for EM-only
weighted lengths rather than changing the meaning of existing output columns.

No native EM ABI change is required. The C++ solver already accepts
per-transcript effective lengths and per-locus gDNA effective lengths.

## 5. gDNA Prior

Current v6 prior construction computes an unweighted expected gDNA count:

```text
eta_g_unweighted = rho_g * L_g
```

If EM uses `L_g^A` but the prior stays at `eta_g_unweighted`, the prior density
becomes too strong in low-exposure loci:

```text
eta_g_unweighted / L_g^A = rho_g / mean_A_g
```

v3a therefore passes this pseudocount to EM:

```text
eta_g_for_em = eta_g_unweighted * L_g^A / L_g
```

Implementation details:

- `gdna_prior_count` should mean the actual prior consumed by EM.
- Add `gdna_prior_count_unweighted` as a diagnostic column.
- Keep `gdna_prior_count_regional` only if it remains useful as a diagnostic;
  it must not be confused with the EM prior.
- If `enable_gdna` is false, preserve the existing disable behavior.
- Guard the ratio with the same `min_value` policy used for `L_g` and `L_g^A`.

This rescaling preserves the calibrated prior rate per unit observed gDNA
opportunity. It does not imply that arbitrary global rescaling of `A` should be
treated as biological signal; `A` remains a relative field normalized by the
exposure builder.

## 6. Output Semantics

Public abundance outputs must stay on raw, unweighted lengths.

This is a hard requirement for v3a: exposure weights are sample- and
assay-specific, while TPM and `effective_length` should remain comparable
across samples. The weighted lengths are internal EM denominators used to route
ambiguous fragments under a capture-aware model. They are not the public
normalization lengths.

Therefore:

- `quant.feather effective_length`: raw transcript effective length;
- `quant.feather tpm` and `tpm_total_rna`: raw TPM using raw effective length;
- `nrna_quant.feather effective_length`: raw nRNA effective length;
- `nrna_quant.feather tpm`: raw total-RNA TPM using raw effective length;
- `gene_quant.feather effective_length`: current raw gene aggregation;
- no replacement of public TPM columns with exposure-weighted TPM.

Add diagnostics, not replacement outputs:

- transcript-level `exposure_weight` or `em_exposure_weight` where practical;
- transcript-level `effective_length_em` only if needed for debugging;
- locus-level aggregate weight ratios for mRNA, nRNA, and gDNA;
- `gdna_prior_count_unweighted` alongside EM-consumed `gdna_prior_count`.

If exposure-weighted TPM is ever useful, add it later as a clearly named new
field such as `tpm_exposure_weighted`; do not overload existing TPM semantics.

## 7. Pipeline Order

`CalibrationResult.regional_exposure` is already available before
`quant_from_buffer()` starts.

Inside `quant_from_buffer()`:

1. Read `calibration.regional_exposure`.
2. Compute unweighted transcript effective lengths using the current path.
3. If exposure is uniform/off, keep the unweighted array and preserve bit-exact
   behavior.
4. If exposure is regional, compute weighted mRNA and synthetic nRNA effective
   lengths and build `t_eff_len_for_em`.
5. Construct EM geometry with `t_eff_len_for_em`, while retaining raw lengths
    for public output.
6. Score fragments normally. Do not apply numerator exposure weights.
7. Build multiloci normally.
8. In `assemble_priors()`, compute `L_g`, `L_g^A`, and `eta_g_for_em`.
9. Partition and run native EM with weighted transcript lengths and weighted
   gDNA lengths.
10. Emit raw TPM/effective-length outputs plus EM-weight diagnostics.

## 8. Validation Contract

### 8.1 Unit Tests

- Class-normalization sentinel: when a region equals its own class reference
  quantile, `A_r` is 1 even if another class has a reference density at least
  100-fold larger.
- Uniform exposure parity: `--regional-exposure off` and auto-uniform mode are
  bit-exact against current `off` outputs.
- Weighted gDNA length parity: uniform exposure delegates to
  `gdna_eff_len_for_loci()` exactly.
- Weighted nRNA length parity: uniform exposure is bit-exact to the existing
  synthetic nRNA effective lengths.
- Transcript exposure parity: uniform exposure gives `Abar_t == 1` for every
  annotated transcript and EM mRNA lengths match raw effective lengths.
- Prior ratio check: `gdna_prior_count / gdna_eff_len` equals
  `gdna_prior_count_unweighted / gdna_eff_len_unweighted` within tolerance.
- Public output parity: under auto-uniform exposure, public `effective_length`,
  TPM, and gene/nRNA TPM columns are bit-exact to `--regional-exposure off`.

Do not require full MAP/VB EM equality under a hand-built constant
non-identity exposure field while also changing the gDNA pseudocount; changing
the absolute pseudocount strength can legitimately alter a prior-influenced
fixed point. Use uniform `A=1` for end-to-end bit-exactness, and use ratio
checks for the prior invariant.

### 8.2 Synthetic Regional Locus

Add a deterministic two-region locus test:

- high-exposure exon, low-exposure intron;
- true gDNA, synthetic nRNA, and mRNA fragments with known origins;
- compare denominator-all v3a against the current gDNA-only behavior from a
  standalone reproduction script, if needed for validation.

Expected behavior:

- denominator-all does not route low-exposure intronic gDNA preferentially into
  nRNA;
- pool counts recover truth within a narrow tolerance chosen before looking at
  the result;
- the standalone gDNA-only reproduction shows the old regression pattern if we
  need a positive-control failure mode.

### 8.3 VCaP Production Gates

Run the VCaP 50/50 RNA/gDNA sample with same-version uniform and v3a auto.
Primary read1, flowcell-truth confusion remains the gate.

| Gate | Threshold |
|---|---|
| gDNA -> RNA leak | <= uniform baseline + 1 percentage point |
| gDNA -> nRNA leak | <= uniform baseline + 0.5 percentage point |
| gDNA -> mRNA leak | <= uniform baseline + 0.5 percentage point |
| RNA -> gDNA leak | within +/- 0.5 percentage point of uniform baseline |
| mRNA -> gDNA leak | <= uniform baseline + 0.5 percentage point |
| Synthetic 24 mRNA error | <= uniform baseline + 0.5% relative error |
| Synthetic 24 nRNA error | <= uniform baseline + 0.5% relative error |
| `--regional-exposure off` | bit-exact regression suite |

If these gates fail, do not tune EM tolerances or priors first. Inspect the
exposure field and weighted effective-length ratios by class and locus.

## 9. Implementation Order

Each step should end with passing focused tests before moving on.

1. Fix class-specific `rho_ref` normalization in `_regional_exposure.py`.
2. Add the class-normalization sentinel test.
3. Add `weighted_nrna_eff_lens()` with uniform parity tests.
4. Add `transcript_exposure_weights()` and mRNA EM-length scaling with uniform
    parity tests.
5. Derive `index.nrna_span_df` and `index.t_to_nrna_span_id` at
   `TranscriptIndex.load()`.
6. Build weighted and unweighted transcript effective-length arrays in
   `pipeline.py` before constructing `AbundanceEstimator`.
7. Split estimator geometry into EM lengths and public output lengths, or add
    equivalent plumbing so output TPM uses raw lengths.
8. Rescale `gdna_prior_count` in `assemble_priors()` and add the prior-ratio
   test.
9. Delete the production call to `_apply_unit_gdna_weights()` and remove stale
   application stats from `summary.json`.
10. Keep quant, nRNA, and gene TPM/effective-length outputs raw; add only
    clearly named EM/exposure diagnostics.
11. Add the deterministic synthetic regional-locus test.
12. Run the focused test set, then full tests.
13. Rebuild only if native code changes; v3a should not require native changes.
14. Run synthetic 24 and VCaP gates before enabling `auto` as production-ready.

## 10. Low-Exposure Denominators

The counterintuitive part of denominator weighting is real: when `L_k^A` is
small, one observed fragment gives more evidence for component `k` than the
same fragment would under a large exposed opportunity. That is not a bug. It is
the same logic as ordinary effective-length correction: seeing reads from a
short transcript implies a higher molecule density than seeing the same number
from a long transcript.

The reconciliation is:

- unexposed regions should not contribute much opportunity, so they shrink
  `L_k^A`;
- a fragment observed from a low-exposure opportunity is rare under low
  abundance, so it legitimately carries more abundance information;
- if every competing component for that fragment has the same local exposure,
  the exposure numerator cancels and only relative component opportunity
  remains;
- low exposure cannot create mass without observed fragments, because EM still
  assigns only observed fragment counts plus priors.

The risk is not mathematical circularity; the risk is noisy or overly small
`A` values creating excessive density inflation for ambiguous fragments. v3a
uses explicit safeguards:

- `A_r` is bounded below by `LOG_A_FLOOR` and above by 1;
- transcript exposure weights `Abar_t` are clamped to `[min_weight, 1]`;
- every EM effective length is clamped to `min_effective_length`;
- gDNA prior density is preserved by `eta_g_for_em = eta_g * L_g^A / L_g`;
- diagnostics report low-exposure components and weight-ratio distributions;
- VCaP gates explicitly catch gDNA-to-RNA and RNA-to-gDNA overcorrection.

Do not add numerator weighting to fix low-exposure discomfort. If v3a appears
to over-prioritize low-exposure components, investigate the exposure floor,
class normalization, and component weight-ratio diagnostics.

## 11. Remaining Concerns

1. **The exposure field may still be the limiting factor.** Class-specific
   normalization fixes the known v3-production bug, but cross-class calibration
   remains approximate. Gate failures should send us back to field
   construction.
2. **mRNA transcript exposure is an approximation.** Exon-length weighted
    `Abar_t` is the right first implementation, but it ignores fragment-length
    edge effects within transcripts. If validation exposes transcript-specific
  artefacts, refine this one layer rather than changing the EM model.
3. **Prior scaling is a density-preservation rule, not a gauge-invariance
   rule.** The exposure builder should normalize away arbitrary global scale;
   the prior test should check `eta_g / L_g`, not demand identical MAP/VB
   solutions under arbitrary constant non-identity `A`.
4. **Output diagnostics need naming discipline.** Public TPM and
  `effective_length` stay raw. Any EM-only weighted length or exposure-weight
    column must be clearly named so users do not compare weighted denominators
  across samples by accident.

## 12. Bottom Line

v3a removes the numerator path and makes regional exposure a consistent
denominator correction for every EM component. The implementation should be
smaller than v2 and safer than v3-production:

- no production `log A(x_u)` mutation;
- no per-candidate exposure scoring;
- weighted effective lengths for gDNA, synthetic nRNA, and mRNA;
- rescaled gDNA prior count;
- raw public TPM/effective-length outputs, with clearly named EM-weight
  diagnostics only;
- no native ABI change.

If v3a passes the synthetic and VCaP gates, ship it behind
`--regional-exposure auto`. If it fails, the next investigation should focus
on exposure-field construction and weighted effective-length geometry, not on
restoring numerator weighting.