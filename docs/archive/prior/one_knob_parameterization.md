# One-Knob Parameterization — Principled Redesign

Date: 2026-05-26
Status: design proposal
Companion to: [adaptive_prior_v1.md](adaptive_prior_v1.md), [../newcal/new_cal_plan_v5.md](../newcal/new_cal_plan_v5.md)

## TL;DR

There is exactly **one** legitimate user-facing tuning knob in this system,
and we already half-have it: a single scalar `rna_confidence ∈ (0, 1)` that
expresses the user's asymmetric loss between calling a fragment RNA when it
was gDNA versus calling it gDNA when it was RNA. Everything else —
`aggregate_prior_strength`, `aggregate_prior_edge_count`,
`aggregate_prior_max_count`, `gdna_prior_logit_bias`,
`gdna_density_confidence`, sweep transfer prior, damping coefficients — is a
symptom of trying to compensate for the absence of that one knob being
honored as a posterior quantile. They should be **deleted** as user surface
and **learned** internally.

The remainder of this document is the derivation.

## 1. Root cause of the knob explosion

A Bayesian pipeline produces a posterior. A prior is a posterior summary
fed into the next inference stage. With four hyperparameters controlling
the same prior, the design is implicitly admitting that none of them is the
right summary. Each tries to patch a failure mode the others created:

| v3 knob | Concept it was trying to encode |
| --- | --- |
| `aggregate_prior_strength` | "How much should the prior fight the data?" — the *ESS* of the prior. |
| `aggregate_prior_edge_count` | "How do we handle one-sided posteriors?" — the *quantile* of the prior. |
| `aggregate_prior_max_count` | "Don't let the prior dominate." — a *safety cap*. |
| `gdna_prior_logit_bias` | "Lean toward RNA when ambiguous." — the *risk preference*. |

Three of these (ESS, edge handling, safety cap) are properties of the
calibration posterior and have no business being user-set. The fourth
(risk preference) is the only one the user can answer.

## 2. The single legitimate user question

> "When an unspliced fragment is ambiguous between gDNA and RNA, how much
> evidence do you want before I call it RNA?"

Answering "a lot" → conservative RNA reporting. Answering "a little" →
sensitive RNA reporting. This is a one-dimensional asymmetric loss
specification.

Under proper Bayesian decision theory with asymmetric loss, the optimal
point summary of a posterior is a **quantile**. The mapping is exact:

```
loss(report, truth):
    if truth = gDNA but report = RNA: cost = q
    if truth = RNA but report = gDNA: cost = 1 - q
optimal report = quantile_q(posterior)
```

So the right user knob is literally the quantile level.

## 3. Definition of `rna_confidence`

```
rna_confidence ∈ (0, 1)
default = 0.5
meaning = posterior quantile at which the gDNA share of unspliced mass
          is summarized when constructing the EM prior
```

Semantics:

```
gdna_share_for_prior(l) = quantile_{rna_confidence}( φ_l | calibration )
```

where `φ_l = D_l / U_l` is the gDNA share of unspliced mass at locus `l`.

- `rna_confidence = 0.5` → median → neutral, Bayes-optimal under symmetric
  loss.
- `rna_confidence = 0.9` → 90th percentile of gDNA share → "I want the
  prior to push EM toward gDNA unless RNA evidence is strong." Conservative
  RNA.
- `rna_confidence = 0.1` → 10th percentile → "I want the prior to push EM
  toward RNA unless gDNA evidence is strong." Sensitive RNA.

This subsumes:

- the existing `rna_lower_confidence` (which is `1 - rna_confidence` in
  spirit — both describe the same axis from opposite sides);
- the existing `gdna_density_confidence` (same axis, evaluated on the gDNA
  density posterior instead of the gDNA share posterior — but those
  posteriors should share their quantile cut by design);
- the v3 `gdna_prior_logit_bias` (a hand-tuned logit-space offset that was
  emulating "use the q-th quantile" through additive bias).

**Recommendation:** rename to one canonical name, e.g.
`--rna-confidence`. Deprecate the others as aliases for one release.

## 4. What gets learned (not set)

Everything below is a *property of the data and the posterior*. None of
these should appear on the CLI.

### 4.1 Prior ESS per locus

Already correctly derived in [adaptive_prior_v1.md](adaptive_prior_v1.md):
moment-match the calibration posterior's (mean, variance) to a Beta and
read off `τ_l = p(1-p)/v - 1`. The data tells us how much it knows; we
honor that.

```
share_l       = quantile_q(φ_l | calibration)           # ← user q
ess_local_l   = moment_match_ess(mean_l, var_l)         # ← learned
alpha_gdna(l) = ess_local_l · share_l
alpha_rna(l)  = ess_local_l · (1 - share_l)
```

Note: `share_l` is the **quantile**, not the mean. This is the key change
from `adaptive_prior_v1.md`, which still uses the posterior mean and tries
to encode user risk preference through `psi` shrinkage. With this design,
risk preference is in the quantile; shrinkage handles only weak-evidence
loci.

### 4.2 Sample-level shrinkage target `ψ` and its precision `τ_global`

Already correctly described in §"Sample-Level Empirical Bayes Shrinkage"
of the adaptive plan. Estimate `ψ` from inverse-variance-weighted regional
gDNA shares; estimate `τ_global` from the between-region variance after
subtracting measurement noise (method-of-moments). Both are pure empirical
Bayes — no user knobs.

When `rna_confidence ≠ 0.5`, the same quantile is applied to the
*shrunken* logit-space posterior, so the user's risk preference propagates
to ambiguous loci that rely on global shrinkage.

### 4.3 All other "knobs"

| Current knob | Replacement |
| --- | --- |
| `aggregate_prior_strength` | Deleted. ESS from moment matching. |
| `aggregate_prior_edge_count` | Deleted. Quantile naturally handles one-sided posteriors. |
| `aggregate_prior_max_count` | Internal safety cap (`max_ess`), not user-facing. |
| `gdna_prior_logit_bias` | Deleted. Folded into `rna_confidence` quantile. |
| `gdna_density_confidence` | Folded into `rna_confidence`. |
| `rna_lower_confidence` | Folded into `rna_confidence` (= `1 − rna_confidence`). |
| `density_min_eff_length` | Internal floor. Derive from `gdna_fl_mean`. |
| `density_max_exposure` | Internal floor. |
| Sweep transfer prior gain, damping η, max calibration passes | Internal numerical defaults; learned damping schedule (e.g. line-search on ELBO) preferred. |
| `fl_prior_ess`, `pool_quality_good/weak` | Internal empirical-Bayes hyperparameters; derive from pooled data. |

### 4.4 What stays exposed (and why)

Only three categories survive on the CLI:

1. **One epistemic dial**: `--rna-confidence` (default 0.5).
2. **Sample-protocol facts that cannot be learned from a single BAM**:
   - `--strand-mode auto|none|fr|rf` (auto first; manual override only
     when the auto detector flags low confidence)
   - `--paired/--single` (already inferred but allow override)
   - `--capture-mode auto|none|hybridization` (auto from coverage profile)
3. **Resource controls**: `--threads`, `--max-memory-gb`.

Everything else is a `--debug-*` flag that defaults off.

## 5. Why this is biologically and statistically natural

### 5.1 Bayesian sparsity

A well-posed Bayesian model has hyperparameters either *learned* (empirical
Bayes) or *elicited from genuine prior information*. Tuning constants that
sit between those two roles are a sign that:

- the model is misspecified (the prior cannot represent the posterior we
  want), or
- the loss function is implicit and inconsistent (every constant is doing
  decision-theory work in disguise).

Collapsing four constants into one quantile fixes the second problem
directly. The first problem (model misspecification) is what the v5
calibration plan is meant to address structurally; the two redesigns are
complementary, not competing.

### 5.2 Biological sparsity

There is one biological question the user can answer that the data cannot:
"In your sample, are you more afraid of false RNA or false gDNA?" Capture
protocol, tissue, library prep history — none of those tell us this
preference. It is a scientific stance about how the user will use the
output. One scalar is the right amount of expressiveness.

Everything else (gDNA contamination level, capture efficiency, strand
specificity, fragment-length distribution, nRNA fraction) is a property of
the library that the data itself can and should reveal.

### 5.3 Why a quantile beats a "strength" + "share" pair

A strength + share parameterization (ESS and `p`) is mathematically
equivalent to a (`α_gdna`, `α_rna`) pair, but it requires the user to set
two numbers when they really only have one preference. Mapping `p` to a
quantile of an estimated posterior lets the data choose the share *and* lets
the user choose the riskiness in a single, interpretable, dimensionally
correct way.

## 6. The full data-flow under one knob

```
BAM scan
  ↓
strand model, FL models                        (learned)
  ↓
four-state region calibration                  (learned, v5 plan)
  ↓
per-region (mean, variance) of gDNA share      (learned)
  ↓
project to locus: (μ_φ_l, σ²_φ_l)              (learned, allocation-weighted)
  ↓
empirical-Bayes shrinkage toward ψ             (learned)
  ↓
posterior over φ_l                             (learned)
  ↓
share_l = quantile_q(φ_l)   ←  q = rna_confidence   (← ONE user knob)
ess_l   = moment_match(μ, σ²)                  (learned)
  ↓
α_gdna(l) = ess_l · share_l
α_rna(l)  = ess_l · (1 - share_l)
  ↓
locus EM (transcript-level)
```

No magic numbers between BAM and EM except `q`.

## 7. Default behavior contract

With `rna_confidence = 0.5` and no other flags:

- A no-gDNA sample: `ψ → 0`, `τ_global` large → ambiguous-locus shares
  shrink near zero → priors RNA-conservative without any logit bias.
- A high-gDNA sample: `ψ` rises, `τ_global` reflects between-locus
  heterogeneity → gDNA is not globally suppressed.
- A locus with strong local strand-deconvolved gDNA evidence: small
  `σ²_φ_l` → high `ess_l` → local evidence dominates.
- A locus with weak local evidence in either sample: large `σ²_φ_l` → low
  `ess_l` → global `ψ` shrinks it; the EM is allowed to decide.
- A locus with no structural gDNA candidate: `α_gdna = α_rna = 0` as today.

No constants were tuned to make these behaviors true. They fall out of the
math.

## 8. Migration path

1. **Phase A (preserves outputs)**: Add `rna_confidence` as a no-op flag
   that maps to existing v3 constants via the documented translation:

   ```
   rna_confidence = 0.5  →  current v3 defaults
   ```

   Wire the four v3 constants behind a single conversion function so the
   CLI now exposes only `--rna-confidence` while the internals still read
   the old fields. Test goldens unchanged.

2. **Phase B (adaptive prior)**: Land the adaptive_prior_v1 plan with one
   modification: replace its "shrunken posterior mean" with the
   "shrunken posterior quantile at `q = rna_confidence`". Default
   `q = 0.5` should reproduce the adaptive-mean behavior in the limit of
   symmetric posteriors. Add golden tests for `q = 0.5`, `q = 0.1`,
   `q = 0.9`.

3. **Phase C (deletion)**: Remove the v3 constants from the public config.
   Keep them as internal `_max_ess`, `_max_global_ess` safety caps with
   hard-coded defaults. Update docs.

4. **Phase D (audit)**: Run the v3 sentinel suite with `rna_confidence` at
   {0.1, 0.5, 0.9}. Confirm:
   - 0.5 reproduces v3 sentinels within tolerance,
   - 0.9 produces strictly higher gDNA assignments (without re-tuning),
   - 0.1 produces strictly higher RNA assignments (without re-tuning).
   This is the test of correctness for the redesign: monotone response to
   the one knob, no tuning required at any operating point.

## 9. What this redesign does not solve

- **Model misspecification at the calibration layer.** If the four-state
  posterior is wrong (e.g., nRNA confused with gDNA), no choice of `q`
  fixes it. The v5 calibration plan addresses this and is the necessary
  complement.
- **Truly unidentifiable samples.** Unstranded, no-capture, single-exon
  loci with no global evidence remain irreducible. The system should
  report large `σ²_φ_l`, low `ess_l`, and a high `prior_conflict_score`.
  The right user-facing action is to inspect those loci, not to set a
  different knob.
- **Computational ergonomics.** Threads, memory, output verbosity are
  separate concerns. Keep their flags.

## 10. Summary in one sentence

**The user gets one dial — a posterior quantile that encodes their
RNA-vs-gDNA risk preference — and every other current "tuning knob" is
either deleted, demoted to an internal safety cap, or learned by empirical
Bayes from the sample itself.**

