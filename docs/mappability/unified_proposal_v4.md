# Unified gDNA Calibration — Proposal v4

_Date: 2026-04-18 · Supersedes `unified_proposal_v3.md`_

## Changes from v3

Two substantive issues raised on v3 are addressed here:

1. **`alignable` Zarr is now a hard requirement.** Bigwig and BED
   fallbacks are removed. We own `alignable`, we trust it, and locking
   the input format simplifies the index, the CLI, and future
   expansion of mappability use throughout rigel.
2. **Overdispersion is handled head-on.** The v3 estimator assumes
   Poisson counts; real gDNA is overdispersed (GC, replication
   timing, chromatin, library prep). v4 keeps the same one-page
   iteration but reinterprets it as a **quasi-Poisson estimator with
   an estimated dispersion**, and uses the dispersion to (a) inflate
   the Fisher information used for inverse-variance pooling and (b)
   define an honest membership rule for $\mathcal A$ that is
   principled rather than asymptotic.
3. **The `f_min` "magic number" is replaced** by a principled
   minimum effective length: a region is admitted to the density
   estimator only if its effective length $E_i = f_i S_i$ is at least
   one expected fragment length $\bar L_F$ (or some integer multiple
   thereof, defaulting to 1). This is the natural lower bound — we
   need at least one fragment's worth of mappable sequence to even
   *observe* a fragment-start in the region.

The estimator stays a one-page function. No new tuning knobs are
introduced; the parameters that exist are the (data-driven) dispersion
and the (data-driven) mean fragment length.

---

## 1. Mappability — `alignable` Zarr only

### 1.1. Decision

The continuous per-base mappability score from the `alignable` Zarr
store is the **single supported source** for $f_i$. No bigwig, no BED
fallback. If the Zarr is missing or unreadable, index build fails
with a clear error.

Rationale:
- We author `alignable` and control its API.
- One input path means one code path, one schema migration story, one
  set of tests.
- Future use of mappability in scoring, multimapper credit, splice
  blacklisting, etc. all benefit from a guaranteed-present continuous
  signal rather than a "BED-or-Zarr-or-none" tri-state.

The cost is operational (Zarr is larger than a BED, and is being
optimised separately). The benefit is permanent simplicity downstream.

### 1.2. Per-region exposure

For a region of span $S_i$ on reference $r$, with `alignable` per-base
score $m(r, b)$ at base $b$ and chosen read length $L$:

$$f_i = \frac{1}{S_i}\sum_{b \in \text{region}_i} m(r, b),\qquad E_i = f_i \cdot S_i$$

The integer-aggregated form $E_i = \sum_b m(r,b)$ is what gets stored
on `regions.feather` (`mappable_effective_length`, float32). $f_i$
itself is recoverable as $E_i / S_i$ when needed for diagnostics.

### 1.3. Read length

CLI flag at index build:
`--mappability-read-length INT` (default 100). Stored in the index
manifest. At quant time, rigel checks observed median read length and
warns if the deviation exceeds 30%. No interpolation, no multi-table
storage. A single value is sufficient because exposure mis-calibration
of ~10% propagates ~linearly to $\hat\lambda_G$, and the existing P₁₀
estimator's failure modes are 1–3 orders of magnitude worse.

---

## 2. The estimator under overdispersion — the critical question

### 2.1. The concern

gDNA fragmentation is not Poisson. Per-region count distributions
exhibit substantial extra-Poisson variance from:
- GC content (well-documented; Benjamini–Speed)
- Replication timing / chromatin state (heterochromatin
  under-represented in standard Illumina libraries)
- Library prep biases (PCR duplication, fragmentation chemistry)
- Probe efficiency (in capture libraries — but that's §3)

Empirically, RNA-seq and DNA-seq counts at the gene/region level
typically follow a negative binomial with dispersion $\phi$ such that
$\text{Var}(k) = \mu + \mu^2 / \phi$ with $\phi$ in the range
$5{-}50$ for gDNA-like signals. Pure Poisson is $\phi = \infty$.

A naive Poisson MLE on overdispersed data gives the right point
estimate (Poisson MLE is consistent for the rate even under
overdispersion — it is a quasi-likelihood estimator) but
**under-estimates standard errors**, which would corrupt the
inverse-variance pooling with the strand pathway.

### 2.2. The good news — point estimate stays valid

The profile-Poisson iteration in §2 of v3 is, mathematically,

$$\hat\lambda = \frac{\sum_{\mathcal A} k_i}{\sum_{\mathcal A} E_i}$$

This is a **method-of-moments / quasi-likelihood** estimator: it
equates observed mean to expected mean over $\mathcal A$. It is
consistent for $\lambda_G$ under any count distribution with
$\mathbb{E}[k_i \mid \text{region } i \text{ is null}] = \lambda_G E_i$,
including all standard overdispersion models (negative binomial,
Poisson-lognormal, beta-binomial sampling). The Poisson assumption
buys us only the variance formula and the membership rule, not the
estimate itself.

So: **the point estimate is robust to overdispersion at the level we
care about.**

### 2.3. The fix — quasi-Poisson with estimated dispersion

We add one estimated nuisance, $\hat\phi$, the **dispersion ratio**:

$$\hat\phi = \frac{1}{|\mathcal A| - 1} \sum_{i \in \mathcal A} \frac{(k_i - \hat\lambda E_i)^2}{\hat\lambda E_i}$$

This is the standard Pearson dispersion statistic from
quasi-Poisson regression. $\hat\phi = 1$ for true Poisson;
$\hat\phi > 1$ for overdispersion; $\hat\phi < 1$ for
underdispersion.

Two consequences:

**A. Honest Fisher information for pooling.** The Fisher information
of $\hat\lambda$ becomes

$$I_\text{density}(\hat\lambda) = \frac{\sum_{\mathcal A} E_i}{\hat\phi \cdot \hat\lambda}$$

This down-weights the density pathway in inverse-variance pooling
when the data are overdispersed (which is exactly when the density
pathway is less reliable). The strand pathway gets the same treatment
— estimate its own dispersion on its own residuals. No magic
constants; the relative weights emerge from the data.

**B. Honest membership rule for $\mathcal A$.** Under pure Poisson
the rule "$k_i < \lambda E_i$" is asymptotically correct: a region is
"null" if its observed count is below its expected count. Under
overdispersion, we should allow more slack. A natural rule using
$\hat\phi$:

$$i \in \mathcal A(\lambda) \iff k_i \le \lambda E_i + z \cdot \sqrt{\hat\phi \cdot \lambda E_i}$$

i.e. "the count is consistent with the null at the upper $z$-tail of
a Poisson with dispersion-inflated variance." Default $z = 2$
(roughly the upper 97.5% of the null distribution). This is the only
new "knob," and it is interpretable as a tail probability rather than
a magic percentile or fraction.

In practice we iterate: estimate $\hat\lambda$ with the simple rule
("$k_i < \lambda E_i$"), compute $\hat\phi$, then re-estimate
$\hat\lambda$ with the dispersion-aware rule, then re-compute
$\hat\phi$. Two outer iterations are sufficient (checked in tests).

### 2.4. What we do not do

- **We do not switch to a full negative-binomial likelihood with
  region-specific $\phi_i$.** That would introduce per-region
  dispersion estimation, regularisation, and EM bookkeeping — the
  exact complexity v2/v3 deliberately avoided. A single global
  $\hat\phi$ captures the bulk of the issue; per-region dispersion is
  not identifiable from a single sample anyway.
- **We do not transform counts.** No log, no Anscombe, no VST. Counts
  enter raw; the dispersion adjustment happens in the variance
  formula.
- **We do not add a winsorisation step** for high-count outliers.
  The $\mathcal A$ membership rule already excludes them naturally
  (any region with $k_i \gg \lambda E_i$ falls outside $\mathcal A$
  and contributes nothing). This is a property the profile-Poisson
  formulation gives us for free, and it is exactly the right
  behaviour for the gDNA-rate estimation problem.

### 2.5. Will it actually be okay?

Two arguments that the answer is yes:

1. **The estimator only uses regions in $\mathcal A$**, which by
   construction excludes the highest-leverage outliers (RNA-rich
   regions). The dispersion of the *remaining* regions — pure-gDNA
   regions, mostly intergenic, mostly low-count — is much closer to
   Poisson than the dispersion of all regions taken together. We
   measure $\hat\phi$ on this filtered set, where it is genuinely
   modest (probably 1–3 in most libraries, not the 10–50 you'd see if
   you tried to fit gene-level RNA counts).
2. **The point estimate is a quasi-likelihood mean equality and is
   consistent for $\lambda_G$ under any null model with the right
   first moment.** The Poisson assumption is only used for the
   uncertainty calculation, which is now explicitly corrected.

We will validate this on simulated data with controllable dispersion
($\phi \in \{1, 5, 10, 50\}$) in Phase 3 testing.

### 2.6. The full algorithm

```
# Inputs
#   counts: k_i (length N)
#   exposures: E_i = f_i * S_i (length N)
#   L_bar: mean fragment length (already estimated by FL model)
# Constants (interpretable, not tuned):
#   min_eff = L_bar              # one fragment's worth of mappable bp
#   z = 2.0                      # ~97.5% upper tail under null
#   tol = 1e-6
#   max_outer = 3, max_inner = 50

# Eligibility — only regions large enough to host a fragment
elig = (exposures >= min_eff)

# Initial λ from the eligible pool
λ ← sum(counts[elig]) / sum(exposures[elig])
φ ← 1.0

for outer in range(max_outer):
    for inner in range(max_inner):
        upper = λ * exposures + z * sqrt(φ * λ * exposures)
        A = elig & (counts <= upper)
        if A.sum() < min_regions:           # trivial-collapse guard
            break  # fall back to previous λ
        λ_new = sum(counts[A]) / sum(exposures[A])
        if abs(λ_new - λ) < tol * λ: break
        λ = λ_new
    # Update dispersion from current A
    resid_sq = (counts[A] - λ * exposures[A])**2 / (λ * exposures[A])
    φ = max(1.0, resid_sq.sum() / (A.sum() - 1))   # never below 1.0

I_density = A.sum() ⋅ ⟨E⟩_A / (φ ⋅ λ)              # Fisher info
return λ, φ, I_density, A
```

About 30 lines of Python including comments. No tuning knobs were
introduced (the `z = 2.0` is an interpretable tail threshold that we
will keep as a documented constant rather than expose as a CLI flag).

---

## 3. Replacing the `f_min` guard with an effective-length lower bound

### 3.1. The right floor

A region cannot meaningfully contribute a fragment-start density
estimate if it cannot even host a single fragment. The natural floor
is therefore expressed in **effective base pairs**, not in
mappability fraction:

$$E_i \ge \bar L_F$$

where $\bar L_F$ is the library's mean fragment length, already
estimated by `FragmentLengthModel` during the BAM scan. With a
typical $\bar L_F \approx 300$ bp, this:

- Subsumes the old "min region length 500 bp" heuristic for binary
  masks (the binary case had $f_i \in \{0, 1\}$, so $E_i = S_i$, so
  the $\bar L_F$ floor is approximately the same magnitude).
- Adapts automatically to the actual library: short-fragment
  libraries (e.g. cell-free DNA, RNA from FFPE) get a smaller floor,
  long-fragment libraries (long-read or large-insert) get a larger
  floor.
- Has a clean physical meaning: "can this region accommodate at least
  one full fragment of the typical length?"

### 3.2. Why this is better than `f_min`

`f_min = 0.1` was a fraction with no physical referent. A 5,000 bp
region with $f_i = 0.05$ has $E_i = 250$ bp — too short to host a
fragment, correctly excluded by the new rule. A 50,000 bp region
with $f_i = 0.05$ has $E_i = 2{,}500$ bp — plenty of effective
length, correctly admitted by the new rule. The old `f_min` rule
would have rejected both, throwing away the second region's signal.

### 3.3. Stronger version (optional, default off)

If we want to require a region to host more than one fragment for
inclusion (e.g. for power), the floor becomes $E_i \ge k \bar L_F$
with default $k = 1$. We expose `--min-effective-length-frags k` for
power users; the default is $k = 1$ which equates to "at least one
fragment fits."

This is a single principled parameter with units (fragment lengths),
rather than a unitless fraction. It is reportable as a count
("admitted X regions covering Y bp of effective genome"), and
auditable.

---

## 4. Hybrid capture — unchanged from v3

Section 3 of v3 (target-BED-aware two-rate estimator) carries over
verbatim. The two-rate split happens **before** the §2 estimator runs;
each subset gets its own $\hat\lambda$, $\hat\phi$, and Fisher info.
Inverse-variance pooling with the strand pathway is performed within
each subset. Without `--target-bed`, behaviour collapses to the
single-pool case.

One v3 subtlety becomes relevant here: when computing $\hat\phi$ on
the on-target subset, we likely see lower dispersion than on the
off-target subset (probes generally enrich uniformly within targets),
so the on-target density pathway gets a relatively higher Fisher info
weight in pooling. This is correct behaviour and emerges from the
data.

---

## 5. Pathway pooling — refined for quasi-likelihood

The pooled estimate becomes

$$\hat\lambda_G = \frac{I_\text{strand}\,\hat\lambda^{(\text{strand})} + I_\text{density}\,\hat\lambda^{(\text{density})}}{I_\text{strand} + I_\text{density}}$$

with $I_\text{density}$ from §2.6 and $I_\text{strand}$ similarly
computed with the strand pathway's own dispersion estimator (Pearson
chi-square / df on opposite-strand exonic counts under the strand
model). At SS = 0.5 the strand pathway's information collapses
naturally; at SS = 1.0 it dominates. No blend weight required.

A consistency diagnostic — the chi-square statistic
$(\hat\lambda^{(\text{strand})} - \hat\lambda^{(\text{density})})^2 \cdot
(I_\text{strand}^{-1} + I_\text{density}^{-1})^{-1}$ — is reported
in `summary.json`. Values $\gg 1$ flag a model-mismatch and prompt
manual review.

---

## 6. Implementation plan (revised)

### Phase 0 — Alignable Zarr probe (unchanged)

`scripts/debug/alignable_landscape_probe.py`: load the Zarr at
$L = 100$, compute $f_i$ for the existing index regions, validate the
distribution and a few chromosome-level sanity checks (centromeres,
exons, repeats).

### Phase 1 — Index changes

1. Hard requirement: `--mappability-zarr PATH` is required at index
   build. Index build aborts with a clear error if missing.
2. Compute and store `mappable_effective_length` (float32) per region.
3. Store mappability read length and Zarr provenance in the index
   manifest.
4. Bump index format version. Update affected tests.
5. Remove `compute_containment_mask` and the old binary-mask code
   paths.

### Phase 2 — Calibration

1. Implement `density_pathway_lambda(counts, exposures, mean_frag_len,
   ...)` per §2.6 — quasi-Poisson with dispersion estimate, $E_i \ge
   \bar L_F$ admission, $z$-tail membership rule.
2. Implement strand-pathway dispersion estimation symmetrically.
3. Replace the v3 inverse-variance pool with the dispersion-corrected
   version (§5).
4. Wire optional `--target-bed` for the two-rate split (carry over
   from v3 §3 / v4 §4).
5. Augment `CalibrationResult` and `summary.json` with $\hat\phi$,
   $|\mathcal A|$, $\sum_\mathcal A E_i$, and the consistency
   chi-square.

### Phase 3 — Tests

1. Synthetic Poisson, recovery within Cramér-Rao.
2. **Synthetic negative-binomial** with $\phi \in \{1, 5, 10, 50\}$.
   Verify (a) point estimate of $\hat\lambda$ stays consistent and
   (b) $\hat\phi$ is recovered within tolerance.
3. Synthetic adversarial RNA (5+ orders of magnitude dynamic range).
4. Effective-length floor: synthetic regions of varying $S_i, f_i$,
   verify only those with $E_i \ge \bar L_F$ are admitted.
5. Hybrid capture two-rate recovery.
6. No-Zarr at index time → index build fails cleanly.
7. Read-length-mismatch warning fires.

### Phase 4 — Validation

1. 24-scenario `stress_calibration.py` sweep.
2. STAR pipeline-validation BAM (unstranded + 50% gDNA, the case the
   current calibration breaks on).
3. Hybrid-capture-like real or simulated dataset.
4. Full hulkrna benchmark sweep.

### Phase 5 — Unmapped-read research probe (separate turn)

Per v2 §5.

---

## 7. Outstanding questions

1. **Default $z$ for the membership rule.** $z = 2.0$ is the obvious
   round number. A sensitivity analysis on simulated data with
   varying $\phi$ should confirm. Worth checking $z \in \{1.5, 2.0,
   2.5, 3.0\}$.
2. **Iterating dispersion — convergence guarantees.** Two outer
   iterations should suffice in practice but we should verify
   numerically and document.
3. **Strand-pathway dispersion estimator.** The strand pathway's
   error model is conditional binomial-on-Poisson; the right Pearson
   statistic needs a brief derivation rather than copy-pasting the
   density formula. Do as part of Phase 2.
4. **Alignable Zarr loader — public API.** Confirm with the
   `alignable` author (us) that the per-base accessor we want exists
   and is reasonably fast for the ~700k region queries we do at index
   build time. If not, add it to `alignable` first.

---

## 8. Summary

| Topic | v4 decision |
|---|---|
| Mappability source | `alignable` Zarr **required**; bigwig/BED removed |
| Read length | CLI flag, default 100 bp; mismatch warning at quant time |
| Estimator | Profile quasi-Poisson with estimated dispersion $\hat\phi$ (§2) |
| Overdispersion handling | $\hat\phi$ inflates Fisher info; dispersion-aware $\mathcal A$ membership rule using upper $z$-tail; point estimate stays consistent |
| Lower bound on regions | $E_i \ge \bar L_F$ (one fragment's worth of effective length) — replaces `f_min` magic number |
| Hybrid capture | Optional `--target-bed`; two-rate split before §2 estimator runs |
| Pathway pooling | Inverse-variance with dispersion-corrected Fisher info per pathway |
| Unmapped reads | BAM-scanner counting only; production integration deferred |

The estimator core remains a one-page function. The only new
parameters are data-estimated nuisances ($\hat\phi$, $\bar L_F$) and
one interpretable tail threshold ($z$). The "magic number" $f_min$
is gone, replaced by a physical floor on effective length.
Overdispersion is handled with a quasi-likelihood treatment that is
both consistent for the point estimate and honest about uncertainty.
