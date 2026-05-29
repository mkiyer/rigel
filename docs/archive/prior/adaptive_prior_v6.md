# Adaptive Prior v6 — Operational Sensitivity / Specificity Dial

Date: 2026-05-26
Status: implementation-ready design
Builds on: `docs/prior/adaptive_prior_v5.md`
Supersedes: nothing in v5; v6 is a strict additive extension.

## 1. Motivation

v5 is parameter-free by design: the prior is whatever the calibration model says
it is, with information entropy as the universal trust weight. That is the
correct *statistical* default.

It does not, however, cover an *operational* choice that some users legitimately
need to make:

> "I would rather lose true low-abundance RNA than admit a false-positive RNA
> call from contaminating gDNA."

This is a utility/risk statement, not a tuning statement. Calibration cannot
infer it from data because it is a property of the downstream application
(clinical diagnostics, panel screening, biomarker QC), not of the sequencing
sample.

v6 exposes exactly one operational dial that maps this risk preference onto the
v5 prior, leaves the statistical default unchanged, and reserves a separate
reporting-layer mechanism for hard clinical audit floors.

## 2. Scope

In scope:

- One scalar dial `rna_call_bias ∈ (0, 1)` (default `0.5`) that asymmetrically
  splits the v5 grouped prior between gDNA and aggregate RNA without changing
  the total prior ESS.
- Documentation and CLI/YAML wiring for that single dial.
- A reserved interface for a future clinical reporting overlay (Plan 3 from the
  evaluation), not implemented in v6.

Out of scope:

- Likelihood-level multipliers (rejected; see Section 6.1).
- Tuning ESS, entropy weight, or anything else inside the v5 core.
- Per-transcript or per-region knobs.

## 3. The Dial

### 3.1 Definition

Let `ω = rna_call_bias` with `ω ∈ (0, 1)` and default `ω = 0.5`.

Let `(α_g, α_r)` be the v5 grouped prior counts produced by
`compute_adaptive_prior()` (before the cap and the structural gate, or
equivalently after — the transformation commutes with both).

Apply the asymmetric split:

```text
α_g' = 2 (1 - ω) · α_g
α_r' = 2 ω       · α_r
```

Properties (proved trivially):

- `ω = 0.5  ⇒  α_g' = α_g, α_r' = α_r`. v5 default is preserved exactly.
- `α_g' + α_r' ≠ α_g + α_r` in general. That is intentional: the dial is a
  log-odds shift on the prior share, not a mass redistribution. If you want
  mass conservation, see Section 3.2.
- `ω → 0  ⇒  α_r' → 0`, gDNA prior doubles. Ambiguous unspliced mass goes to
  gDNA unless the data strongly contradict it.
- `ω → 1  ⇒  α_g' → 0`, RNA prior doubles. Ambiguous unspliced mass goes to
  RNA unless the data strongly contradict it.
- The dial is dimensionless and bounded. It has no calibration-noise
  sensitivity: small calibration perturbations produce small prior shifts.

### 3.2 Mass-Conserving Variant (chosen)

To keep the v5 ESS cap meaningful, v6 uses the mass-conserving variant:

```text
T   = α_g + α_r                    # v5 prior total at locus l
s_v = α_r / max(T, ε)              # v5 RNA share
s_w = 2 ω · s_v / (2 ω · s_v + 2 (1 - ω) · (1 - s_v))

α_g' = T · (1 - s_w)
α_r' = T · s_w
```

This is the logit-space shift form: it preserves `α_g' + α_r' = T`, so the v5
`MAX_ESS` cap and `U_l` cap remain exact. At `ω = 0.5`, `s_w = s_v`,
recovering v5 exactly.

Both forms are mathematically natural. v6 uses the mass-conserving form so
that v6 priors never exceed v5 priors in total ESS.

### 3.3 Where It Applies

The dial applies once, at the end of `compute_adaptive_prior()`, after the cap
and the structural gate, to the final `(α_gdna_add, α_rna_add)` pair per locus.

It does not apply to:

- Calibration internals (`p_states`, `U_r`, entropy weights).
- The structural `has_gdna_candidate_l` gate (still hard-zeros both alphas).
- The native EM solver (no native changes; the dial is a Python-side rescale).
- Any non-prior code path.

### 3.4 Why Not a Likelihood Multiplier

Plan 2 from the evaluation (`Λ` multiplier on gDNA likelihood inside the
E-step) is rejected for three reasons:

1. It conflates the likelihood (what the data say) with the loss function
   (what the user prefers). Once mixed, model log-likelihoods are no longer
   comparable across runs.
2. It is mathematically equivalent to a per-fragment log-prior shift, which is
   the role the v6 dial already fills more cleanly at the prior level.
3. It requires native EM changes for no additional expressivity.

### 3.5 Why Not Post-EM Veto (Yet)

Plan 3 (a post-EM per-fragment confidence floor `τ`) is the right mechanism
for a *clinical reporting overlay*, but it is a different abstraction:

- It changes the report, not the model. Estimator θ remains unconstrained.
- It requires per-fragment (or per-EC) RNA responsibility from converged θ,
  which rigel does not currently materialize.
- It is most useful *on top of* a biased prior (set `ω < 0.5`, then apply
  `τ` for audit), not as a replacement.

v6 reserves `--min-fragment-rna-confidence` for a future v7 reporting overlay
and documents it as a non-goal here. See Section 9.

## 4. User Interface

### 4.1 CLI

Single new flag on `rigel quant`:

```text
--rna-call-bias FLOAT   Operational sensitivity/specificity dial for the
                        unspliced gDNA-vs-RNA prior split. Range (0, 1).
                        0.5 (default) is the unbiased v5 prior.
                        < 0.5 favors gDNA assignment (RNA-specificity mode).
                        > 0.5 favors RNA assignment (RNA-sensitivity mode).
```

### 4.2 YAML

Top-level under `quant:`:

```yaml
quant:
  rna_call_bias: 0.5
```

### 4.3 EMConfig

Add one field, only:

```python
@dataclass(frozen=True, slots=True)
class EMConfig:
    ...
    rna_call_bias: float = 0.5

    def __post_init__(self):
        if not (0.0 < self.rna_call_bias < 1.0):
            raise ValueError("rna_call_bias must be in (0, 1)")
```

Nothing else in `EMConfig` changes relative to v5.

### 4.4 Removed Keys (still removed, no change from v5)

All v3/v4 keys remain rejected with the v5 migration error. v6 does not
reintroduce any of them, and `rna_call_bias` is not an alias for any removed
key.

## 5. Implementation

### 5.1 Helper Change

Extend `compute_adaptive_prior()` with a `rna_call_bias` argument:

```python
def compute_adaptive_prior(
    *,
    region_arrays: RegionArrays,
    multi_loci: list[MultiLocus],
    p_states: np.ndarray,
    unspliced_total: np.ndarray,
    has_gdna_candidate: np.ndarray,
    rna_call_bias: float = 0.5,
    max_ess: float = MAX_ESS,
) -> AdaptivePriorResult:
    ...
```

Final step inside the helper (after cap and gate):

```python
if rna_call_bias != 0.5:
    T   = alpha_gdna_add + alpha_rna_add
    s_v = np.divide(alpha_rna_add, T, out=np.zeros_like(T), where=T > 0)
    two_w   = 2.0 * rna_call_bias
    two_omw = 2.0 * (1.0 - rna_call_bias)
    num = two_w * s_v
    den = num + two_omw * (1.0 - s_v)
    s_w = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
    alpha_rna_add  = T * s_w
    alpha_gdna_add = T * (1.0 - s_w)
```

That is the entire algorithmic change.

### 5.2 Pipeline Wiring

`assemble_priors()` reads `em_config.rna_call_bias` and passes it through:

```python
result = compute_adaptive_prior(
    ...,
    rna_call_bias=em_config.rna_call_bias,
)
```

No other code path changes.

### 5.3 Output Schema

Add to `loci.feather`:

```text
prior_rna_share_v5     # α_r / (α_g + α_r) before the dial
prior_rna_share_final  # α_r' / (α_g' + α_r') after the dial
```

Add to `summary.json.prior_policy`:

```json
"rna_call_bias": 0.5,
"rna_share_shift_p50": 0.0,
"rna_share_shift_p90": 0.0
```

`rna_share_shift_p*` are the per-locus distribution of
`prior_rna_share_final - prior_rna_share_v5`, exposed so a user can confirm
how much their dial actually moved the prior.

### 5.4 Flags

Add one optional bit:

```text
0x8 PRIOR_BIAS_APPLIED          rna_call_bias != 0.5 and locus prior was shifted
```

## 6. Tests

### 6.1 Unit (`tests/test_adaptive_prior.py`, new cases)

1. `ω = 0.5` leaves every per-locus `(α_g, α_r)` identical to v5.
2. `ω → 0` drives `α_r → 0` and `α_g → T` for every locus with `T > 0`.
3. `ω → 1` drives `α_g → 0` and `α_r → T` for every locus with `T > 0`.
4. Mass conservation: `α_g' + α_r' == α_g + α_r` to float tolerance.
5. Structural gate still wins: `has_gdna_candidate_l == False` ⇒ both alphas
   exactly 0 regardless of `ω`.
6. Cap still wins: `α_g' + α_r' ≤ min(U_l, MAX_ESS)`.
7. Monotonicity: for fixed locus, `α_r'(ω)` is strictly increasing in `ω`
   on `(0, 1)` whenever `T > 0`.

### 6.2 Integration

- Default `ω = 0.5`: every existing v5 sentinel passes bit-identical to v5.
- `ω = 0.1` on a zero-true-gDNA sample: RNA recovery degrades gracefully
  (RNA loss bounded, no catastrophic siphon), gDNA component absorbs the
  shifted prior mass.
- `ω = 0.9` on a true high-gDNA sample: gDNA recovery degrades; the user is
  warned via `summary.json.rna_share_shift_p90`.

### 6.3 CLI

- `--rna-call-bias 0.5` is the default and is omitted from CLI help unless
  `--show-defaults`.
- `--rna-call-bias 0.0` and `--rna-call-bias 1.0` are rejected with a clear
  error pointing at this doc.
- YAML round-trip preserves `rna_call_bias`.

## 7. Backwards Compatibility

- v5 has not shipped as code yet; v6 lands as part of the same cutover.
- Defaulting `ω = 0.5` means a user who never sets the dial gets exactly v5
  behavior. No goldens shift.
- No previously removed key is reintroduced.

## 8. Documentation

- This file is the source of truth for the dial.
- `docs/prior/adaptive_prior_v5.md` is updated to add a one-line cross-reference
  to v6 in its "Non-Goals" section: "Operational risk dial: see v6."
- `docs/MANUAL.md` adds a short "RNA call bias" section under quantification
  options, framed as a *clinical / operational* choice, not a tuning knob,
  with example settings:
  - `0.5` — default, statistically unbiased (recommended).
  - `0.2` — high RNA specificity (recommended for diagnostic screens).
  - `0.8` — high RNA sensitivity (low-input or single-cell where missing a
    true transcript is worse than picking up gDNA).

## 9. Future Work (Reserved, Not v6)

### 9.1 Plan 3: Clinical Reporting Overlay

Reserved CLI flag: `--min-fragment-rna-confidence FLOAT` (default off).

When enabled, after EM converges:

1. For each fragment (or equivalence class) compute aggregate RNA
   responsibility from converged θ and likelihoods.
2. Fragments with `P(RNA) < τ` have their mass redirected entirely to the
   gDNA component for a *separate* reporting table.
3. Estimator θ is unchanged; a new `clinical_quant.feather` is emitted
   alongside `quant.feather`.

This requires:

- Materializing per-EC RNA responsibility (new pass in `estimator.py`).
- A second abundance accumulation that respects the veto.
- A clear schema split between "model" and "clinical" outputs so neither is
  silently corrupted.

It is intentionally deferred so the v6 cutover stays small. The reserved flag
name is documented here to prevent a future naming collision.

### 9.2 Calibration-Side ω

A future v7 could let calibration suggest a per-sample default ω from the
estimated false-positive RNA risk (e.g., from the four-state confusion
profile). v6 deliberately does not do this: user intent is the only signal
that should drive an operational dial.

## 10. Summary

v6 adds exactly one user-facing parameter to v5: `rna_call_bias ∈ (0, 1)`,
default `0.5`, which asymmetrically splits the v5 grouped prior between gDNA
and aggregate RNA while preserving total ESS. It is a mass-conserving
logit-space shift applied once at the end of `compute_adaptive_prior()`. It
requires no native changes, leaves the v5 default behavior bit-identical, and
gives clinical users a single, auditable, dimensionless dial for the
sensitivity/specificity tradeoff. A clinical reporting overlay (`τ`) is
reserved for v7.
