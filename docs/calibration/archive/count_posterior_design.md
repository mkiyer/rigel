# Phase 4 — the count-module posterior (detailed design & critique)

**Status:** design proposal for critique. 2026-06-10. Expands the Phase-4 sketch in
`count_channel_capture_design.md`. The count module currently emits a point `g_count =
clip(ρ_local·eff/M)`; this gives it an honest **posterior** (mean + variance) over the gDNA fraction,
to feed (a) the precision-weighted strand→count blend (Phase 1's weight) and (b) the FP-rate quantile
(Phase 5). **Read §1 first — it reframes what this phase can and cannot do.**

## 1. The critical distinction: variance ≠ bias (read first)

The motivation we keep returning to is the flagship +1pp regression (stranded+capture): the blend
mixes 4% of a *capture-biased* count estimate into an excellent strand estimate. It is tempting to
think "give the count an honest variance and the inverse-variance weight will down-weight it." **That
is only true if the count's error shows up as variance.** Under hybrid capture the count's error is a
**systematic bias** (the exon gDNA density imputed from depleted off-target boundaries → ~2× low), and
a biased estimate can be *low-variance* (precise but wrong): both flanking exon–intron boundaries sit
at the same enriched exon edge, so they **agree** on the same biased-low density. A variance model
built from **boundary disagreement** would see agreement → low variance → **high count weight** → the
bias is trusted. So:

> **Phase 4 (a variance model) does NOT, by itself, fix the flagship bias regression.** That regression
> needs a *less biased mean* — the point-5 unspliced-fraction / transport-style attribution that uses
> the exon's own enriched count (a separate work item, call it Phase 4-mean). Phase 4 here is the
> *variance* (Phase 4-var): it is what the **FP-rate quantile (Phase 5)** requires, and what lets the
> weight defer to strand where the imputation is *genuinely uncertain* (disagreeing/sparse anchors) —
> not where it is confidently biased.

This corrects my earlier claim that "Phase 4's bias-aware weight fixes the flagship." It would only if
we model **bias**, which is harder (we don't observe truth at runtime). So the honest plan is:
- **Phase 4-var (this doc):** the count posterior's *variance* — counting noise + imputation
  uncertainty. Enables Phase 5 and honest weighting under genuine uncertainty.
- **Phase 4-mean (separate, = point-5):** debias the exon mean under capture. This is what fixes the
  flagship and the unstranded+capture residual. The two compose: debiased mean + honest variance.

The remainder designs Phase 4-var, and §6 critiques whether/how a *bias* signal can be folded in.

## 2. What the count module emits

For every node, a posterior over its gDNA fraction `g`:
- **mean** `μ_g = g_count = clip(ρ_local·eff_len / M, 0, 1)` (unchanged).
- **variance** `σ_g²`, from two independent sources:

```
σ_g²  =  (eff/M)² · [ σ²_count(N_g)            # counting/dispersion noise of the gDNA count
                     + σ²_impute(ρ_local) · eff² ]   # imputation uncertainty of the density
```

For an **observable** node `ρ_local` is its own count (no imputation term; `σ²_impute=0`) and `μ_g≈1`
with small variance — and these nodes are usually strand-routed anyway, so their variance barely
matters. For an **imputed** (exon/AMBIG) node the **imputation term dominates** — that's the
high-value piece (§4). The counting term (§3) is the baseline/floor.

## 3. Counting noise — Poisson baseline + (optional) local dispersion

**Baseline (recommended for v1):** Poisson, `σ²_count(N_g) = N_g`. Principled, parameter-free, no
fitting. Its only flaw is being *too* precise (no overdispersion → under-states uncertainty), but the
imputation term dominates anyway, so a Poisson floor is a fine v1.

**Refinement — local dispersion as a function of count (your kernel/kNN idea).** RNA-seq counts are
NB-overdispersed (`Var = μ + α·μ²`); a *global* `α` is exactly what exploded under capture (it booked
on/off-target *mean* heterogeneity as dispersion — the teardown). Your proposal localizes it: estimate
`α(x)` from observable seed nodes **near count x**, so on-target (high count) and off-target (low
count) get their own dispersion. Two estimators:

- **kNN:** for query count `x`, take the `k` nearest seed counts; `α̂(x)` from their excess variance.
  Knob: `k` (a magic number).
- **Gaussian kernel over sorted counts (recommended):** convolve a Gaussian (bandwidth `h` in
  **log-count** space — counts span orders of magnitude) across the sorted seeds; `α̂(x)` is the
  proximity-weighted local excess variance. Smooth, no hard bin edges, weights near neighbours more.
  Knob: bandwidth `h`.

**The estimator (per kernel window around expected count μ):**
```
α̂(μ) = [ Σᵢ wᵢ ((Nᵢ − μ)² − Nᵢ) ] / [ Σᵢ wᵢ μ² ]     # NB excess variance, kernel-weighted
```
with `wᵢ = Gaussian(log Nᵢ − log μ; h)`, over observable gDNA-pure seeds.

**The caveat to critique (§6):** within a kernel window the seeds do **not** share one true rate
(capture gives each its own), so `(Nᵢ−μ)²` still conflates *true-rate variation within the window*
with *sampling dispersion* — the same conflation as the global `α`, only **localized and reduced**
(similar-count seeds have more similar rates). Narrower `h` → less conflation but noisier. So local `α`
is an *upper bound* on the true dispersion. Given that, and that imputation variance dominates, **v1
uses the Poisson baseline only**; the local-dispersion estimator (kernel/kNN) is parked in its own
design doc so the idea isn't lost: `count_local_dispersion_design.md`. (Note: the §4 variance~mean fit
gets the conflation-free benefit the count-dispersion lacks, because it pools *paired* same-node
disagreements rather than *unpaired* cross-node count spreads — so §4 is the better-founded path.)

## 4. Imputation uncertainty — boundary disagreement (the high-value piece)

For a non-observable node, `ρ_local` is imputed from its **anchoring observable boundary sides**. Per
your point, at most **two** boundaries are informative:
- **Adjacent case** (e.g. `intron│exon│intron`): the exon's own left/right observable boundary sides
  give density estimates `d_L`, `d_R`.
- **Run case** (`intron│exon│exon│exon│intron`): the *run* is anchored by the boundary at its left
  edge (`d_L`) and its right edge (`d_R`); interior nodes are carried/interpolated between them.

**Disagreement → variance.** The imputation's uncertainty is how much the anchors disagree, plus each
anchor's own count noise:
```
σ²_impute(ρ)  =  ¼·(d_L − d_R)²            # disagreement: the variance of the 2-point mean
               + ½·(σ²_{d_L} + σ²_{d_R})    # each anchor's own counting noise (§3) propagated
               + σ²_floor                    # prior floor (so agreeing high-count anchors aren't 0)
```
- **2 agreeing anchors** → small disagreement → low variance (precise — *but possibly biased*, §1).
- **2 disagreeing anchors** (capture enriches one side) → large variance → the blend defers to strand.
- **1 anchor** (reference/gene edge): no cross-check → variance = the single anchor's count noise +
  an inflated floor (we can't validate it).
- **0 anchors** (no-anchor region, global fallback): maximal variance (the count is uninformative).

**Run interiors.** A node `j` positions into the run carries from the run-edge anchors; its variance
grows with distance from the nearer anchor. Two model options:
- **Linear-interpolation variance:** treat the carried density as a linear interpolation between
  `d_L` and `d_R`; the interpolation variance is `(1−t)²σ²_{d_L} + t²σ²_{d_R} + t(1−t)(d_L−d_R)²`
  where `t` is the fractional position — maximal mid-run, where we're least anchored. Clean, parameter-
  free, and it *is* the variance of a linear estimator between two noisy endpoints. (Recommended.)
- **Distance-penalty:** the edge disagreement × a monotone function of run-depth. Heuristic; needs a
  knob. Prefer the interpolation form.

**The variance~mean fit (recommended — denoises the per-node 2-point estimate).** A single node's
`¼(d_L−d_R)²` is a 1-dof, very noisy variance. But **pool** the `(mean, disagreement)` pairs across
*all* 2-anchor imputations genome-wide and fit a smooth **`σ²_impute ~ f(mean)`** trend; then read each
node's variance off the curve at its own mean. This is the count-channel analog of DESeq2's
mean-variance trend — and it has a property the count-dispersion model (§3) lacks:

> **It is conflation-free.** `d_L` and `d_R` are two estimates of the *same* node's density — a
> **paired** sample at one true rate — so `(d_L−d_R)²` is genuine variance, **not** confounded by
> rate-variation across different nodes (the conflation that sank the global/local count `α`). The
> 2 anchors share the node's true rate; their disagreement is pure measurement + within-node-gradient
> uncertainty, both legitimate components of imputation error.

**Validated on the benchmark** (`scripts/debug/diag_variance_mean.py`, gdna1000): under **capture-on**
the pooled pairs follow **`var ∝ mean²`** (slope 2.00 in log-log) — exactly the NB law `Var = α·μ²`,
cleanly fittable with a single `α`. (Capture-off is noisier — slope ~6 — but the densities are tight
and the variance is tiny there, so it barely matters.) So the fit is: `σ²_impute(μ) = α_impute·μ²`,
`α_impute` = the pooled slope from the 2-anchor disagreements (a robust, conflation-free estimate),
with the per-node form as the fallback where the pooled fit is degenerate. The 1-anchor / 0-anchor
cases read the curve at their mean with an inflated floor.

Two fit forms to choose between (open question §6): a **parametric NB** `α_impute·μ²` (one number,
interpretable, matches the data) vs a **non-parametric kernel** `f(μ)` (no functional assumption).
The data says `μ²` fits well under capture → lean parametric NB.

## 5. Posterior parameterization & how it's consumed — NB vs Beta vs Gaussian

These aren't competing choices; they live at different layers, and your NB intuition is correct:
- **NB is the COUNT noise model.** RNA-seq/gDNA *counts* are negative-binomial (`Var = μ + α·μ²`). That
  governs the *count*-level uncertainty — and the §4 fit confirmed the imputation variance is NB-shaped
  (`var ∝ mean²`). So NB is right *for the counts/density*.
- **Beta is the FRACTION posterior.** The count module's output is a *fraction* `g ∈ [0,1]` (the gDNA
  share of the node's unspliced mass), not a count. The natural distribution over a bounded fraction is
  the **Beta**. We take the NB count/density uncertainty and *propagate* it to a variance `σ_g²` on the
  fraction, then represent that as `Beta(g)` with mean `μ_g` and `concentration = μ_g(1−μ_g)/σ_g² − 1`.
- **Gaussian is just a computational shortcut.** For the Phase-5 quantile we need `F⁻¹(q)`; a truncated
  Gaussian `(μ_g, σ_g²)` gives `μ_g + Φ⁻¹(q)·σ_g` in closed form. Use it as an approximation to the Beta
  where convenient (they agree away from 0/1).

So: **NB models the counts → propagated to a Beta over the fraction → Gaussian-approximated for the
quantile.** No conflict. Then:
- **Weight (Phase 1 refinement):** `w = I_strand / (I_strand + I_count)`, `I_count = 1/σ_g²`. Where the
  imputation is uncertain (disagreeing/sparse anchors) → `I_count` small → `w→1` → strand governs.
  Where it's confident → `I_count` large → the count keeps its share. (Subject to the §1 caveat: a
  *confident-but-biased* count still gets weight — the flagship needs Phase 4-mean.)
- **FP quantile (Phase 5):** the combined per-node posterior (strand BB ⊗ count Beta, weighted) has a
  variance; the user's quantile `q` reads `μ + Φ⁻¹(q)·σ`. The count posterior is what makes `q`
  meaningful on count-routed (unstranded/AMBIG) nodes.

## 6. Open questions & critique (let's decide these together)

1. **Variance vs bias (the big one, §1).** Do we accept that Phase 4-var does *not* fix the flagship
   (that's Phase 4-mean/point-5), and scope Phase 4-var to enabling Phase 5 + honest uncertainty? Or do
   we attempt a **bias signal**? One candidate: under capture the exon *interior* count `C` vastly
   exceeds `ρ_boundary·eff` (the boundary-imputed expectation) — that *ratio* (`C / (ρ_boundary·eff)`)
   is a per-node **bias indicator** (it's ≫1 exactly where the boundary imputation under-calls). We
   could inflate `σ_g²` (or shift `μ_g`) by it. But that's really the point-5 *mean* fix wearing a
   variance costume — cleaner to do it as the mean. **Recommend: Phase 4-var scoped to variance; bias
   handled by Phase 4-mean.**
2. **Local dispersion (§3): build it or Poisson-only?** Given the within-window conflation and that
   imputation dominates, is the kernel `α(x)` worth its bandwidth knob in v1, or do we ship Poisson and
   revisit? (Lean: Poisson v1.)
3. **2-point disagreement is noisy.** A per-node variance from 2 anchors is itself high-variance. Is
   the floor + linear-interpolation form enough, or do we want the global calibration scale `c` (§4.2)?
4. **Run-interior model:** linear-interpolation variance (parameter-free) vs distance-penalty
   (knob). (Lean: interpolation.)
5. **Beta vs truncated-Gaussian** for the posterior, and how to combine with the strand BB for the
   Phase-5 quantile (Gaussian approx vs explicit mixture).
6. **Validation:** the nascent + capture benchmark now lets us check the variance *calibrates* — does
   high `σ_g²` actually predict high per-node imputation error (against oracle truth)? That's the
   acceptance test for Phase 4-var, distinct from net-leak.

## 7. Implementation sketch

- `density_model.node_gdna_density` → also return `σ_g²` per node: the imputation variance from the
  anchoring boundary sides (track each node's `d_L`, `d_R`, their count noise, and run-position `t`),
  plus the Poisson counting term. Extend `NodeDensity` with `gdna_frac_var` (or `count_precision`).
- `strand_deconv._deconv_per_node` → the weight `w` uses `I_count = 1/σ_g²` (Phase 1 currently uses the
  fixed `(2κ−1)²`; this makes it depth/uncertainty-aware) — gated on the §1 decision.
- The FP quantile (Phase 5) consumes the combined posterior.
- Acceptance: nascent+capture benchmark net-leak **and** the per-node `σ_g²`-vs-error calibration check.

## 8. Recommendation

Ship **Phase 4-var** as: Poisson counting floor + **boundary-disagreement imputation variance**
(linear-interpolation form for runs, floor for sparse anchors), as a Beta posterior. Defer the local
dispersion model (§3) and the global calibration scale (§4.2) unless the validation demands them. Pair
it with **Phase 4-mean (point-5)** for the flagship bias — they are complementary, and the benchmark
(now with nascent + capture) is what tells us how much of the residual each one closes.
