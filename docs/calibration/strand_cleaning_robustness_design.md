# Robust strand-cleaning — design doc (CALIBRATION_TODO Issue #1)

**Status:** design / options analysis. 2026-06-09. No implementation yet. For review.

## 1. The problem

The calibrator separates each node's unspliced fragment mass into gDNA vs RNA *by strand* using a
closed-form linear unmix (`density_model.strand_clean_gdna_frac`):

```
gdna_frac  ĝ = clip( (s − κ) / (½ − κ),  0, 1 )
```

where `s = sense/total` is the node's observed sense fraction, gDNA has sense rate ½ (unstranded,
double-stranded), and RNA has sense rate κ = `rna_sense_frac`. This is the method-of-moments unmix
of `s` between the two components.

**It is unbiased but its variance explodes as κ→½ (or N→0).** With `N` fragments and mixture sense
rate `p(g) = κ + g(½−κ)`:

```
Var(ĝ) = Var(s) / (½−κ)² = p(1−p) / [ N (½−κ)² ]
```

so the precision (Fisher information for g from the strand channel) is

```
τ_strand(g) = N (½−κ)² / [p(1−p)]   ∝   N (2κ−1)²
```

As κ→½ the strand carries **zero information** (gDNA at ½ and RNA at κ≈½ are indistinguishable), yet
the closed form returns a finite, wildly amplified number (it divides by ½−κ→0). As N→0 the same.
The estimator is *correct in expectation* but *infinitely uncertain* exactly where we lean on it
hardest. SS=½ is not an edge case to special-case — it is the maximum-uncertainty limit of a
continuum, and the estimator must degrade gracefully across that continuum.

**Empirical confirmation.** The de-bias prototype (which summed the *unclipped* clean) blew up at
ss=0.50 (κ=0.499) → a pure-RNA library was called **42.9 % gDNA**. The shipped code clips ĝ to
`[0,1]` and special-cases `|½−κ| < 1e-6 → 1.0`, which *bounds* the per-node value but (a) is a sharp,
discontinuous cliff that only triggers at κ *exactly* ½ — at κ=0.499 it still amplifies noise — and
(b) clipping does not down-weight a noisy estimate when it is *consumed*: a low-SS node's garbage
ĝ is used downstream as if it were reliable.

## 2. Where the clean is used, and what is NOT affected

The closed-form clean is a **pre-deconvolution** point estimate, used in two places:

1. **`node_gdna_density`** — `clean_count = ĝ · total` per node → the per-region gDNA *density* (the
   count clue, the global density, and the count-prior mean `count_gdna_frac`).
2. **`joint_deconv._compute_side`** — the boundary-side count-prior density *and* the gDNA strand-fit
   **seed weight** (`gdna_weight`).

**Important scoping fact:** the *final* per-node gDNA fraction (the joint deconvolution,
`_joint_per_node`) is **already robust** — it combines a count-prior Beta × a Beta-Binomial strand
*likelihood*, and when κ→½ the strand likelihood goes flat, so the posterior correctly reverts to
the count prior. The uncertainty IS modelled there. The fragility is **only** in the closed-form
point-estimate pre-uses (1) and (2), which feed the count density / seed weights *before* the
Bayesian combination. So this issue is narrow and self-contained: **make the pre-deconv clean
degrade gracefully**, without touching the (already-correct) joint posterior.

## 3. Requirements

A clean fix must:
- **Degrade gracefully** as κ→½ and N→0 — no explosion, no discontinuous cliff.
- **Preserve the unbiased estimate where strand is informative** (don't regress the stranded win).
- **Be cheap** — this is a per-node vectorized pre-step over all regions/boundaries; it cannot run a
  heavy per-node optimizer.
- **Minimal magic numbers** (project rule) — prefer quantities that are *the information content*,
  not tunables.
- **Robust across regimes** — stranded/unstranded, sparse/dense, capture on/off, with/without nRNA.

## 4. Options

### O1 — Status quo: clip + sharp `|½−κ|<1e-6` cliff
The current code. **Reject.** Discontinuous; the explosion lives at κ=0.499 (just outside the
tolerance), not at exactly ½; clipping bounds but doesn't down-weight when consumed.

### O2 — Precision-weighted shrinkage toward a fallback g₀  *(recommended)*
Treat the strand clean as a measurement of g with precision `τ = (2κ−1)²·N` (the Fisher info, up to
the slowly-varying p(1−p)≈¼), and shrink it toward a fallback mean `g₀` with prior precision `τ₀`:

```
g_robust = [ τ·ĝ + τ₀·g₀ ] / [ τ + τ₀ ]
```

Substituting ĝ=(s−κ)/(½−κ) and τ=4(½−κ)²N (using 2κ−1 = −2(½−κ)), the 1/(½−κ) **cancels**:

```
g_robust = [ 4N(½−κ)(s−κ) + τ₀·g₀ ] / [ 4N(½−κ)² + τ₀ ]
```

This is smooth and well-defined at κ=½ (→ g₀), with **no division by (½−κ)**. It is exactly the
Gaussian/conjugate posterior mean for g under a strand-count likelihood + a prior centred at g₀.
- κ→½ or N→0 ⇒ τ→0 ⇒ g_robust → g₀ (fall back, no explosion).
- strong SS + many fragments ⇒ τ≫τ₀ ⇒ g_robust → ĝ (the unbiased data estimate; win preserved).
- **No clipping cliff, no special case.** The `_UNSTRANDED_TOL` branch disappears.

Two sub-decisions (see §6): the **target g₀** and the **prior strength τ₀**.

### O3 — Full per-node Bayesian posterior (Beta-Binomial)
Model `K ~ BetaBinomial(N, p(g), od)` with a prior on g; report the posterior mean. The uncertainty
is modelled exactly (wide posterior at low SS/N ⇒ reverts to prior). **This is what the joint
deconv already does** for the final estimate. Using it *again* as the pre-deconv clean duplicates
that machinery per node and is heavier (grid/optimizer per node). O2 is the cheap Gaussian
approximation to this; the exact version is reserved for the deconv. **Defer** unless O2's Gaussian
approximation proves inadequate at very small N.

### O4 — Eliminate the pre-clean; fold all strand into the joint deconv
Don't strand-clean the count density at all. Let the count clue use the *raw* unspliced count and
let the deconv's strand likelihood do 100 % of the gDNA/RNA separation. **Conceptually cleanest**
(one Bayesian combination, strand used once). But: the count clue's *density* and the *global gDNA
density hyperparameter* would then be gDNA+RNA-contaminated, and the acyclic architecture introduced
the pre-clean precisely to get a clean global density without a density→deconv→density loop.
Removing it is a substantial restructure with its own circularity to solve. **Larger scope;** worth
noting because it also resolves a latent concern (strand is currently used *twice* — in the
count-prior mean via the clean, and again in the strand likelihood — a mild double-count at high SS).

### O5 — Reliability-gate the *influence*, not the estimate
Keep ĝ (clipped), but scale its *downstream weight* by reliability `(2κ−1)²` — e.g. the count-prior
concentration or the seed weight is multiplied by `(2κ−1)²`. Cheap, and it stops a low-SS clean from
*dominating*. But the noisy estimate still enters (just down-weighted), it does not give a usable
point value at low SS, and it conflates strand reliability with the count concentration (which the
count-overdispersion work, Issue #2, already governs). A partial measure, subsumed by O2 (O2's
denominator `τ+τ₀` *is* a reliability gate, applied to the estimate itself).

## 5. Trade-off summary

| option | graceful @ κ→½ | preserves stranded win | cost | new knobs | scope |
|---|---|---|---|---|---|
| O1 cliff (now) | no (cliff at exactly ½) | yes | trivial | tol (1e-6) | shipped |
| **O2 shrinkage** | **yes (smooth)** | **yes** | **vectorized, trivial** | **g₀, τ₀** | **local, drop-in** |
| O3 Bayesian | yes | yes | heavy (per-node) | prior | duplicates deconv |
| O4 restructure | yes (deconv handles) | yes | medium | — | architectural |
| O5 gate influence | partial | yes | trivial | weight scale | partial |

## 6. Recommendation

**Adopt O2 — precision-weighted shrinkage** (`g_robust = [4N(½−κ)(s−κ) + τ₀g₀]/[4N(½−κ)²+τ₀]`).
It is the minimal, smooth, principled fix: it models the strand's *reliability* via its own Fisher
information `(2κ−1)²·N` (not a tunable — it is the information content), reverts to a sensible
fallback exactly where strand is uninformative, preserves the unbiased estimate where strand is
strong, costs nothing (vectorized closed form), and deletes the ad-hoc `_UNSTRANDED_TOL` cliff.

The two sub-decisions, to settle together:

**(a) The target g₀.** Make it node-type-aware — this falls out of the existing structure:
- **Count-observable regions** (intron / intergenic): gDNA-pure by signature, so g≈1 — shrink
  toward **g₀ = 1** (at low SS, default to "all gDNA", which is the truth there). This also matches
  the current `TS_NONE`/`AMBIG` override.
- **Non-observable (exon) regions & boundary sides**: a genuine gDNA/RNA mix whose split is
  *unidentifiable* at the node when strand is weak — shrink toward **g₀ = a pre-deconv global gDNA
  fraction**. That global is estimable *without strand* from the count-observable regions
  (gDNA-pure), which **breaks the circularity** (the post-deconv global density is not needed as the
  target). Candidate: `Σ gDNA_count(observable) / Σ unspliced_count(touching the locus)`, or simply
  the library gDNA fraction implied by the observable density.

**(b) The prior strength τ₀.** This is the one true knob and needs the no-magic-numbers discipline.
Options, in order of preference:
- **τ₀ = a small fixed pseudo-count in the same units as τ** (e.g. the information of ~1 fully-
  stranded fragment, τ₀ = (2·1−1)²·1 = 1, i.e. "one stranded observation"). Principled, scale-free,
  and tiny — it only bites when the node's own strand info is comparably tiny (low SS *and* low N).
- τ₀ tied to the strand-overdispersion prior already in `CalibrationConfig` (reuse the existing
  prior_weight machinery for consistency).
- Cross-validate τ₀ on the SS sweep (ss 0.99 → 0.50) so the shrinkage onset matches where strand
  actually loses power.

**Validation plan.** Re-run the SS sweep (ss ∈ {0.99, 0.90, 0.75, 0.60, 0.50}) × {capture on/off} ×
gDNA levels and confirm: (i) ss=0.50 no longer produces a false-gDNA blowup (the 42.9 % → ~global-
rate); (ii) ss=0.99 is unchanged (win preserved); (iii) smooth monotone behaviour in between (no
cliff). Add a unit test asserting `g_robust` is finite and → g₀ as κ→½, and = ĝ as τ≫τ₀.

**Related concerns to track (not in scope here):**
- *Strand double-use* (O4's motivation): strand enters both the count-prior mean (via the clean) and
  the strand likelihood. O2 doesn't fix this; flag for a separate look — at high SS it mildly
  over-weights strand. If it matters, O4 (fold strand into the deconv only) is the principled
  resolution, but it's a larger architectural change.
- The shrinkage interacts with the count-overdispersion concentration (Issue #2): both gate how much
  a node's count is trusted. Keep them orthogonal — O2 governs the *mean's* reliability (strand),
  count overdispersion governs the *concentration's* honesty (count). Verify they compose cleanly.
