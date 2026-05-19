# gDNA Regional Exposure Revised Plan v3

**Status**: design plan, not implemented.
**Date**: 2026-05-19
**Supersedes**: `gdna_exposure_revised_plan_v2.md`, `gdna_exposure_plan_v3.md`.
**Implements lessons from**: the VCaP regression on v3 production, the v2
design review, and a clean rederivation of the capture-aware EM equations.

## 0. Why v3 (this document) exists

v2 was correct but over-engineered. Working from first principles, the
per-fragment exposure numerator factor `A(x_u)` is the **same constant
across every competing component at a given fragment** and therefore
**cancels exactly** in the EM E-step posterior. The numerator
manipulations in v3-production (gDNA only) and v2 (gDNA + nRNA) were not
the load-bearing mechanism — the load-bearing mechanism is the
component effective length `L_k → L_k^A`. The numerator was either
harmful (v3, asymmetric) or harmless-but-unnecessary (v2, symmetric).

The corrected design is denominator-only, applied to every component:

> Capture is molecule-agnostic. Every fragment from genomic position
> `x` was retained with probability `A(x)`, regardless of whether it
> originated as cDNA from mRNA, cDNA from nRNA, or cDNA from
> chromosomal DNA. The only EM-visible consequence is that each
> component's effective length shrinks to its
> exposure-weighted form, `L_k^A = ∫ p_k(s) A(s) ds`.

This document replaces both prior plans with that single rule and the
machinery to compute it.

---

## 1. The generative model and why denominator-only is sufficient

### 1.1 Proper capture-aware likelihood

For component `k` with intrinsic per-position fragment density
`p_k(x, ℓ)` (uniform over its footprint, weighted by its FL pmf), the
observed density of a fragment at `(x, ℓ)` after capture is:

$$
f_k(x, \ell) \;=\; \frac{p_k(x, \ell)\,A(x)}
                          {\displaystyle\int p_k(s, \ell')\,A(s)\,ds\,d\ell'}
            \;=\; \frac{p_k(x, \ell)\,A(x)}{L_k^A}
$$

where

$$
L_k^A \;=\; \sum_\ell h_k(\ell)\int A(s)\,\mathbb{1}[\text{valid for }k]\,ds
\;=\; \bar A_k\,L_k.
$$

`A(x) ∈ (0, 1]` is the assay accessibility at genomic position `x`,
independent of which molecule was captured.

### 1.2 The cancellation

The EM E-step posterior at unit `u`:

$$
P(k\mid u) \;=\; \frac{\theta_k\,f_k(u)}{\sum_j \theta_j\,f_j(u)}
           \;=\; \frac{\theta_k\,p_k(u)\,A(x_u)/L_k^A}
                     {\sum_j \theta_j\,p_j(u)\,A(x_u)/L_j^A}
           \;=\; \frac{\theta_k\,p_k(u)/L_k^A}
                     {\sum_j \theta_j\,p_j(u)/L_j^A}.
$$

`A(x_u)` is identical in every term and cancels. The EM trajectory and
fixed point depend only on the per-component effective lengths
`L_k^A` (and the unweighted per-position densities `p_k(u)`, which are
unchanged).

### 1.3 What this means in plain language

- **The numerator manipulation on `gdna_log_liks` (v3 production) and on
  per-candidate nRNA `log_liks` (v2) does nothing useful for the EM**
  *if* it is applied symmetrically to every competing component.
- **If it is applied asymmetrically** (v3 production: gDNA only), it
  injects a fictitious `log A(x_u) − log Ā_k` log-odds shift against
  the weighted component at every non-locus-mean position. That is the
  observed VCaP regression in one line.
- **The whole effect of exposure on EM is encoded in `L_k^A`.** Apply
  consistently across components, and the EM will silently and
  correctly reweight per-fragment density.

### 1.4 What `θ_k` means after the change

EM fixed point: `θ_k = (1/N) Σ_u P(k|u)` = observed fraction of
fragments routed to component `k`. The true molecule abundance `n_k`
relates to `θ_k` via:

$$
\theta_k \;\propto\; n_k \, \bar A_k,
\qquad\text{so}\qquad
\frac{n_k}{L_k} \;\propto\; \frac{\theta_k}{L_k^A}.
$$

Downstream TPM-style normalization (`θ_k / L_k`) currently produces
"observed effective abundance per unit unweighted length". After the
change, the **mathematically meaningful per-length abundance** is
`θ_k / L_k^A`. v3 ships TPM with the weighted denominator and emits the
unweighted denominator as a diagnostic for one release of backward
auditability.

---

## 2. What we remove from v3-production and v2

| Component | v3-production behaviour | v2 behaviour | v3 (this plan) |
|---|---|---|---|
| gDNA numerator (`gdna_log_liks += log A(x_u)`) | **applied** | applied | **removed** |
| nRNA per-candidate numerator (`log_liks[j] += log A(x_u)`) | not present | added | **not added** |
| mRNA per-candidate numerator | not present | guard-rail (Mode 3b) | **not added** |
| `genomic_midpoint` in `ScoredFragments` | required | required | **kept for diagnostics; not used by EM** |
| `_apply_unit_gdna_weights` in `pipeline.py` | active | active | **removed from pipeline; kept as `--regional-mode-internal=gdna_only_legacy` for benchmarking only** |
| `_apply_candidate_nrna_weights` in `pipeline.py` | n/a | added | **not added** |
| `safe_unit` mask, cross-ref skip, multimapper skip for numerator | active | active | **not needed; deleted** |
| Weighted `L_g^A` per locus | applied | applied | **applied** |
| Weighted `L_n^A` per nRNA span | not applied | applied | **applied** |
| Weighted `L_m^A` per mRNA transcript | not applied | not applied | **applied** |
| gDNA prior count `η_g` rescaling | not applied (Strategy A) | optional (Strategy A default) | **applied (Strategy B mandatory)** |

The pipeline becomes structurally simpler. The native scoring ABI does
not need to change at all for v3 (the `genomic_midpoint` field stays as
a diagnostic; if a future plan needs it, the path is already wired).

---

## 3. The exposure field `A(x)`

### 3.1 Construction (unchanged in mechanism, fixed in normalization)

Reuse the v3-production `_regional_exposure.py` skeleton:
per-region conservative-count tables `Y_r`, opportunity `E_r`,
EB-shrunk `ρ̂_r`, deterministic null-spread auto-uniform signal.

Mandatory fix vs v3-production: **class-specific reference quantiles**.

```python
for cname in CLASSES:
    log_rho_ref_k = log(_weighted_quantile(rho_hat_k, E_k, REFERENCE_QUANTILE))
    raw_log_ratio_r = min(log(rho_hat_r) - log_rho_ref_k, 0.0)
    log_A_r = signal_k * raw_log_ratio_r
    log_A_r = max(log_A_r, LOG_A_FLOOR)
```

The global-max `rho_ref` is retained in `summary.json` only as a
diagnostic; it is never used to compute `A_r`.

### 3.2 Acknowledged residual limitation

Cross-class comparability is still imperfect because INTRON's `ρ̂`
measures contained-fragment density and EXON-INTRON's `ρ̂` measures
boundary-crossing density. After class-specific normalization each
class is internally consistent. A per-locus integral `∫_M A(s) ds` that
crosses class boundaries integrates approximately-comparable
quantities, but the residual error is bounded and is the same for every
component, so it cancels in the gDNA-vs-nRNA-vs-mRNA log-odds for
fragments observed at the same position (which is precisely what EM
looks at).

A fully unified smoothed `coverage(x)` estimator is deferred to a
follow-up plan.

### 3.3 Sentinel test

Construct three synthetic classes with `rho_ref_EXON / rho_ref_INTRON
≥ 100`. Assert that a region whose `ρ̂` equals its own class reference
quantile receives `A_r` within `1e-9` of `1.0`. This locks the v3
class-normalization bug closed.

---

## 4. Weighted effective lengths

### 4.1 gDNA: per multi-locus

Already implemented as `weighted_gdna_eff_len_for_loci` in
`src/rigel/calibration/_exposure.py`. Reuse unchanged. Uniform fast
path remains bit-exact.

### 4.2 nRNA: per unique span

Add `weighted_nrna_eff_lens` to `src/rigel/calibration/_exposure.py`:

```python
def weighted_nrna_eff_lens(
    nrna_spans: NrnaSpanTable,        # (ref, strand, start, end, span_id)
    rna_fl: FragmentLengthModel,
    exposure: RegionalGdnaExposure,
    *,
    min_value: float = 1.0,
) -> np.ndarray:
    """Per-unique-span exposure-weighted RNA-FL effective length."""
```

Algorithm mirrors the gDNA path:

- one `(ref, start, end)` span at a time;
- iterate `positive_ell = nonzero(rna_fl.pmf[1:]) + 1`;
- merge per-ref valid start windows into midpoint windows;
- accumulate `rna_fl.pmf[ell] * exposure.weighted_length_on_ref(...)`;
- `min_value=1.0` clamp;
- uniform fast path is bit-exact to existing `nrna_eff_lens`.

### 4.3 mRNA: per annotated transcript

Add `weighted_mrna_eff_lens` to the same module:

```python
def weighted_mrna_eff_lens(
    t_exon_blocks: ExonBlockTable,    # per transcript, per exon: (ref, start, end)
    rna_fl: FragmentLengthModel,
    exposure: RegionalGdnaExposure,
    *,
    min_value: float = 1.0,
) -> np.ndarray:
    """Per-annotated-transcript exposure-weighted RNA-FL effective length."""
```

The integration domain is **not** a single genomic interval; it is the
union of exonic blocks for the transcript. Fragments may span multiple
exons (across spliced junctions), so the standard
`l_eff_for_transcript` algorithm — which works in transcript coordinates
and projects valid windows back to genomic intervals — is the right
template. Each projected genomic window contributes
`exposure.weighted_length_on_ref(...)` instead of its raw length.

Uniform fast path is bit-exact to existing mRNA effective lengths.

### 4.4 EM hand-off

`AbundanceEstimator` consumes one `t_eff_lens` for all transcript-like
components plus a per-locus gDNA effective length. v3 replaces both:

```python
t_eff_lens_weighted = t_eff_lens.copy()
synth = index.is_synthetic_t
annotated = ~synth
t_eff_lens_weighted[synth]     = weighted_nrna_eff_lens(...)[t_to_nrna_span_id[synth]]
t_eff_lens_weighted[annotated] = weighted_mrna_eff_lens(...)
```

The per-locus gDNA effective length passed to EM is `L_g^A` from
`weighted_gdna_eff_len_for_loci`. Native ABI unchanged.

### 4.5 nRNA span bookkeeping

`index.t_to_nrna_span_id` and `index.nrna_span_df` are derived at
`TranscriptIndex.load()` from existing transcript columns (no index
format bump). The mapping is `-1` for annotated rows.

---

## 5. gDNA prior count

The v6 Dirichlet pseudocount `η_g(M)` was calibrated as the expected
gDNA fragment count per unit **unweighted** opportunity:
`η_g = ρ_g · L_g`. The EM places prior mass on the gDNA component
proportional to `η_g / L_g`.

Under v3 the EM denominator is `L_g^A = ā_g L_g`, so the un-rescaled
prior mass becomes `η_g / L_g^A = η_g / (ā_g L_g) = ρ_g / ā_g`, which is
inflated by `1/ā_g`. To preserve the calibrated prior rate per unit
**observed** opportunity, rescale:

$$
\eta_g^{(\mathrm{v3})}(M) \;=\; \eta_g(M)\,\cdot\,\frac{L_g^A(M)}{L_g(M)}.
$$

This is mandatory in v3, not optional. The pre-rescaled `η_g` is
emitted as `gdna_prior_count_legacy` in `loci.feather` for one release
of backward auditability, then removed.

---

## 6. TPM and downstream abundance interpretation

The mathematically correct per-length abundance becomes
`θ_k / L_k^A`. Two practical choices:

- **Ship TPM = θ_k / L_k^A** (recommended): TPM denominator is now the
  weighted effective length, so reported TPMs reflect estimated true
  molecule abundance per unit unweighted footprint.
- **Continue TPM = θ_k / L_k**: would mean TPM stays at "observed
  per-effective-length", which is what we have today but is inconsistent
  with the EM having corrected for capture.

v3 ships option 1. `quant.feather` gains a `tpm_legacy_unweighted`
column for one release to allow side-by-side audit, then it is removed.

---

## 7. Pipeline order

`CalibrationResult.regional_exposure` is already built by
`rigel.calibration._orchestrator.calibrate(...)` before
`quant_from_buffer(...)` starts.

Inside `quant_from_buffer`:

1. Read `calibration.regional_exposure` (class-normalized per §3.1).
2. Build `t_eff_lens_weighted` (annotated rows via
   `weighted_mrna_eff_lens`, synthetic rows via
   `weighted_nrna_eff_lens`). Uniform-mode fast path is bit-exact.
3. Construct `TranscriptGeometry` / `AbundanceEstimator` with
   `t_eff_lens_weighted`.
4. Score fragments into global `ScoredFragments` (unchanged native ABI).
5. Build `multi_loci`.
6. `assemble_priors(... regional_exposure=...)` — computes `L_g^A` per
   locus and rescaled `η_g^{(v3)}` per locus (§5).
7. **No** per-unit or per-candidate numerator weighting on `log_liks`
   or `gdna_log_liks`. (Both removed.)
8. `partition_and_free`.
9. Run EM with `t_eff_lens_weighted` and `L_g^A`.
10. Construct TPM with weighted denominator (§6).

There is no `safe_unit` mask, no cross-ref skip, no multimapper
exposure handling, no genomic-midpoint dependency in the EM path. All
of that lived only in service of the (now removed) numerator.

---

## 8. Validation

### 8.1 Cancellation invariant (sentinel)

CI test: construct a small locus, set `exposure` to a hand-built
constant non-identity field `A(x) = A_0 ≠ 1` on every position the
locus touches. Assertions:

- `weighted_gdna_eff_len_for_loci == A_0 * gdna_eff_len_for_loci`
  (within fp tolerance).
- `weighted_nrna_eff_lens == A_0 * nrna_eff_lens` (within fp tolerance).
- `weighted_mrna_eff_lens == A_0 * mrna_eff_lens` (within fp tolerance).
- EM posteriors and assignment counts are numerically equivalent to the
  unweighted-baseline run (within fp tolerance).
- `η_g^{(v3)} = η_g · A_0`.

Uniform mode (`A_0 = 1` everywhere) must be **bit-exact**.

### 8.2 Synthetic two-region locus

`tests/test_regional_exposure_locus.py`: one deterministic locus,
5 kb high-exposure exon (`A_r = 1.0`), 95 kb low-exposure intron
(`A_r = 0.01`). 100 true gDNA, 100 true nRNA, 50 true mRNA fragments
with known positions.

Assertions:

- v3 (denominator-only on all components) recovers gDNA, nRNA, mRNA
  counts within ±10% of truth and does not over-call any pool by more
  than 15%.
- Legacy `gdna_only_legacy` mode reproduces the v3-production
  misrouting of intronic gDNA to nRNA (positive assertion that the
  diagnostic mode preserves the historical bug).
- Runtime budget: < 2 s.

### 8.3 VCaP gates (production criteria for `auto`)

| Gate | Threshold |
|---|---|
| G1: gDNA→RNA leak | ≤ same-version uniform baseline + 1 pp |
| G2: gDNA→nRNA leak | ≤ uniform baseline + 0.5 pp |
| G3: RNA→gDNA leak | within ±0.5 pp of uniform baseline |
| G4: mRNA→gDNA leak | ≤ uniform baseline + 0.5 pp |
| G5: gDNA→mRNA leak | ≤ uniform baseline + 0.5 pp |
| G6: synthetic 24-condition mRNA relative error | ≤ uniform baseline + 0.5% |
| G7: synthetic 24-condition nRNA relative error | ≤ uniform baseline + 0.5% |
| G8: `--regional-exposure off` regression suite | bit-exact vs current `off` |

G5 is new in v3: with all three components weighted, a mis-calibrated
field could push mRNA mass into gDNA at exonic positions. We gate on it
explicitly.

### 8.4 Diagnostic emissions in `loci.feather`

- `gdna_eff_len_unweighted`, `gdna_eff_len`, `gdna_eff_len_weight_ratio`
  (already in v3-production).
- `mrna_eff_len_unweighted`, `mrna_eff_len_weight_ratio` (new aggregate;
  summed over annotated transcripts in the locus).
- `nrna_eff_len_unweighted`, `nrna_eff_len_weight_ratio` (new aggregate;
  summed over nRNA spans in the locus).
- `gdna_prior_count_legacy` (pre-rescale; one release only).

---

## 9. CLI and config surface

User-facing flag unchanged:

```text
--regional-exposure {off, auto}
```

`auto` resolves to denominator-only-on-all-components at runtime. `off`
is bit-exact to current uniform behaviour.

Internal-only enum for benchmarking and forensics
(`--regional-mode-internal`, hidden, env-var configurable; not
documented in `--help`):

```text
off
gdna_only_legacy            # reproduce v3-production regression
gdna_denominator_only       # weight only L_g^A; leave nRNA, mRNA unweighted
denominator_all             # v3 default (== auto)
denominator_all_unrescaled  # v3 default but skip §5 prior rescaling
```

Production config files should never reference the internal flag.

---

## 10. Implementation checklist

Ordered to minimise blast radius per step. Each step ends at a
passing-tests state.

1. **Class-specific `rho_ref` fix** in `_regional_exposure.py` + §3.3
   sentinel test. Land first; this alone removes the worst part of the
   v3-production field bug.
2. **`weighted_nrna_eff_lens`** in `_exposure.py` + bit-exact uniform
   fast path test.
3. **`weighted_mrna_eff_lens`** in `_exposure.py` + bit-exact uniform
   fast path test.
4. **`index.t_to_nrna_span_id` / `index.nrna_span_df`** derived at
   `TranscriptIndex.load()`.
5. **Build `t_eff_lens_weighted`** in `pipeline.py` and pass to
   `AbundanceEstimator`. Bit-exact parity test when
   `exposure.mode == "uniform"`.
6. **Rescale `η_g`** per §5 in `assemble_priors` (or in the prior
   construction path that already returns `gdna_prior_count`).
7. **Remove `_apply_unit_gdna_weights`** from the production pipeline
   (move behind `--regional-mode-internal=gdna_only_legacy`).
8. **TPM denominator switch** in `AbundanceEstimator.get_quant_df`
   (or wherever TPM is computed). Add `tpm_legacy_unweighted` column.
9. **Locus diagnostic columns** per §8.4.
10. **§8.1 cancellation invariant** CI test.
11. **§8.2 synthetic two-region locus** CI test.
12. **VCaP gate sweep** (§8.3) and synthetic 24-condition regression.
13. **Decide.** If G1–G8 pass: ship behind `auto`. If any gate fails:
    do not ship; the field-construction layer (§3.1) is the next
    suspect, not the EM machinery.
14. **Document.** Mark v3-production and v2 plans superseded; link this
    document from the calibration README.

---

## 11. Risks and open questions

### 11.1 Field-quality risk

Class-specific normalization makes each class internally consistent
but does not unify cross-class scales. If `Ā_m` ends up
systematically larger than `Ā_g` because of a calibration artefact
(boundary-crossing density vs contained density), mRNA gets a
relative penalty per fragment vs gDNA. This is the failure mode G5
catches. If G5 fails, the fix is in §3, not §4–§7.

### 11.2 mRNA effective length blast radius

Every annotated transcript's effective length changes. Synthetic
24-condition is the safety net (G6, G7). Off-target transcripts (low
`Ā_m`) get smaller `L_m^A`, which inflates per-fragment mRNA density;
combined with the rescaled `η_g`, this should be self-consistent, but
we should watch for inflation of mRNA counts on completely-uncaptured
transcripts.

### 11.3 What about `nrna_quant.feather` interpretation?

Same logic as TPM: per-length abundance becomes `θ_n / L_n^A`. The
column should be added; the legacy column kept for one release.

### 11.4 Footprint vs midpoint for `A`

Denominator-only does **not** sample `A` at fragment midpoints; it
integrates `A(s)` over candidate windows. So the footprint-vs-midpoint
ambiguity that v2 inherited from v3-production is gone.

### 11.5 Numerator can still help downstream

A per-fragment `log A(x_u)` is useful for fragment-level QC scoring
(e.g., "this fragment lives in a low-exposure region; flag it") even
though it does not influence EM. Keep `genomic_midpoint` in
`ScoredFragments` as a diagnostic; do not feed it into EM.

---

## 12. Bottom line

Hybrid capture acts on cDNA molecules without knowing their biological
origin. The mathematically correct way to reflect that in EM is to
replace every component's effective length with its exposure-weighted
form `L_k^A`. The per-fragment numerator `log A(x_u)` is mathematically
redundant when applied symmetrically and actively harmful when applied
asymmetrically. v3 is denominator-only on every component, with the
gDNA prior count rescaled to keep its calibration meaningful.

The implementation is structurally smaller than both v3-production and
v2: no per-unit or per-candidate numerator mutation, no `safe_unit`
mask, no cross-reference skip, no multimapper-skip bookkeeping, no
native ABI change.

The validation contract is §8.3 gates G1–G8 on VCaP plus the synthetic
24-condition regression. If those pass, ship behind `--regional-exposure
auto`. If they fail, the next investigation is the exposure field
construction (§3), not the EM machinery (§4–§7).
