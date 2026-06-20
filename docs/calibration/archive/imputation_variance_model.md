# Imputation variance model — separating sampling noise from biological dispersion

**Status:** the first-principles specification of the **imputation reliability** (the precision of a neighbour→
node propagation), for the Step-2 imputation layer. It refines the precision model of
`CALIBRATION_ARCHITECTURE.md` (§1.2/§2): the imputation is the **second** of the two legitimate count→precision
channels. Read §0 of the architecture first (count = statistical power, never a composition vote).

The mechanism (which nodes impute which, transcript-structure gating) is in `CALIBRATION_PLAN_v6.md` +
`rna_imputation_transcript_structure.md`; this document specifies only **how precise the propagation is**.

---

## 1. THE PROBLEM — density buries the count

Imputation propagates a neighbour's density estimate onto a node; its precision is a **modeled reliability** — a
`var~mean` curve giving the variance of the prediction error as a function of the density level. The pre-Step-2
builder fit this curve in **density space**: for each adjacency `(src → dst)`, a point `(μ = ρ_dst,
raw_var = (ρ_dst − ρ_src)²)`, fit `log σ² ~ log μ`.

Density `ρ = C/L` (count over effective length) **buries two things** that drive the real uncertainty:

| | count `C` | density `ρ` | sampling certainty |
|---|---|---|---|
| node A | 10 / 2 bp | 5 | **wild** (ρ = 5 ± 1.6) |
| node B | 1000 / 200 bp | 5 | **stable** (ρ = 5 ± 0.16) |

Same density, vastly different certainty. A density-only `var~mean` cannot tell them apart, so it fails **twice**:

- **At fit time it is contaminated.** A residual `(ρ_dst − ρ_src)²` between two low-count nodes is large because
  the *counts* are noisy, not because the *biology* is variable. The curve reads high variance at low density and
  mistakes sampling noise for biological dispersion.
- **At apply time it is blind.** Looking the curve up by density alone, it cannot see that *this* predictor has
  10 counts vs 1000 — so a noisy neighbour and a reliable one get the same imputation precision.

Both failures have one cause: **density discards the discrete-count nature of the data.** The fix is to model the
count explicitly — as statistical power, not as a composition vote (the §0 invariant).

---

## 2. THE DECOMPOSITION (first principles)

A node has a true rate (density) `λ`. We observe `C` fragments over effective length `L`; the estimate is
`ρ = C/L`. Fragments arrive by (to first order) **Poisson** sampling, so `Var(C) = C`, and by the delta method
the density estimate's sampling variance is

```
Var_samp(ρ) = Var(C) / L²  =  C / L²  =  ρ / L .
```

This is the whole point: it depends on **both** `C` and `L`. Node A: `10/2² = 2.5`. Node B: `1000/200² = 0.025`
— a 100× difference at identical density, exactly matching the table above.

Now the imputation residual `R² = (ρ_dst − ρ_src)²`. Writing each observation as truth + sampling error
(`ρ = λ + ε`, `Var(ε) = Var_samp`):

```
ρ_dst − ρ_src  =  (λ_dst − λ_src)  +  (ε_dst − ε_src) ,
```

so, taking expectations (the cross terms vanish — sampling errors are independent of the rates and of each
other):

```
E[R²]  =  (λ_dst − λ_src)²   +   Var_samp(ρ_src)  +  Var_samp(ρ_dst)
       =   σ²_bio(μ)         +   ρ_src/L_src      +  ρ_dst/L_dst .
       └── LEARN ──┘             └──────── COMPUTE ────────┘
```

The observed squared residual is the sum of:

- **`σ²_bio(μ)`** — the **biological / structural dispersion**: how much the true rate genuinely differs between
  adjacent nodes (the irreducible imputation error even with infinite counts). A function of the density level
  `μ`; this is what we want to **learn**. (gDNA is smooth ⇒ small; RNA is spiky ⇒ larger and rising with μ.)
- **`ρ_src/L_src + ρ_dst/L_dst`** — the **Poisson sampling noise** of the two density estimates. **Computed**, not
  learned, from the observed counts and lengths.

The pre-Step-2 curve fit `E[R²]` against `μ` — i.e. the *sum* — and so blended the two. The decomposition lets us
separate them.

---

## 3. THE FIT — learn only the excess

Because the Poisson term is **computed**, we learn only the *excess over it*:

```
σ²_bio(μ)  =  E[ R² − ρ_src/L_src − ρ_dst/L_dst  │  μ ] .
```

The curve now represents the genuine biological dispersion, de-confounded from count noise. On a uniform field
(`λ` constant, `σ²_bio = 0` truly), the expected residual is `ρ_src/L_src + ρ_dst/L_dst`, the subtraction gives
**0**, and the curve correctly learns no dispersion — the contamination is removed exactly (worked example §9).

**The fit nuance — a Poisson offset, not a per-point subtraction.** Subtracting per point makes some residuals
negative (sampling-dominated pairs), and `log(negative)` breaks the log-space monotone spline; flooring at zero
then logging would bias the curve *up* (keeping only the positive-excess points). The clean estimator is the
**offset model**: the response stays the raw `R² ≥ 0`, and the model is

```
E[R²]  =  σ²_bio(μ)  +  V_p ,        V_p = ρ_src/L_src + ρ_dst/L_dst   (a KNOWN per-point offset)
```

so the monotone P-spline learns `σ²_bio(μ)` such that `curve(μ) + V_p` matches `R²`. This is the standard
Poisson/variance-offset GLM idea realized inside the existing `MonotoneVarMean`; it composes with the Jensen
df-offset already there (which corrects the log-transform bias of a dof=1 variance estimate). Negatives never
arise because the offset is added to the model, not subtracted from the data.

---

## 4. THE APPLY — the count re-enters as statistical power

When we impute `dst` from a predictor `src`, the prediction of `dst`'s **true** rate is `μ_dst = ρ_src` (§5), and
its error vs `λ_dst` is `(λ_dst − λ_src) − ε_src`, with variance

```
Var(λ_dst │ ρ_src)  =  σ²_bio(μ_dst)  +  Var_samp(ρ_src)  =  σ²_bio(μ_dst)  +  ρ_src/L_src ,

  τ_impute(dst │ src)  =  1 / ( σ²_bio(μ_dst)  +  ρ_src/L_src ) .
```

The destination's own observation enters separately (its strand likelihood in the sweep, §1.1); the imputation is
a **prior**, so its precision is the dispersion plus the **predictor's** sampling noise.

This is where the count comes back — and it is the **second count→precision channel** of §0:

```
   tight prior  ←  high-count, low-dispersion predictor   (ρ_src/L_src small,  σ²_bio small)
   diffuse prior ←  low-count or high-dispersion predictor (ρ_src/L_src large OR σ²_bio large)
```

The asymmetry is automatic. Imputing the noisy node A from the stable node B: `τ = 1/(0 + 0.025) = 40` — trust
it. The reverse, B from A: `τ = 1/(0 + 2.5) = 0.4` — don't. **The count sets reliability, never composition.** A
predictor's 1000 fragments make its density estimate trustworthy; they say nothing about whether the target is
gDNA or RNA. That is the §0 invariant honored, not broken — the same way the strand BB Fisher info is.

**Density → fraction.** The sweep solves fractions, not densities. The imputation prior on the component fraction
uses the **definitional density→fraction jacobian** `(eff/mass)²` (the same normalizer kept for the global prior
in Step 1, `CALIBRATION_ARCHITECTURE.md` §2): `μ_f = μ_dst·L_dst / M_total_dst` and `τ_f = τ_impute · (M_total/eff)²`.

---

## 5. DESIGN DECISIONS

**No bootstrap.** Drawing `C_k ~ Poisson(Ĉ)` and taking `Var(C_k/L)` *converges to* `Ĉ/L²` — the closed form is
the bootstrap's exact expectation. The analytic `ρ/L` is deterministic, exact, has no `K`-sample knob, needs no
RNG, and is simpler. Resampling would only add Monte-Carlo noise to a quantity we can write down.

**No direct fit; the mean is the identity.** Modeling `ρ_dst ~ ρ_src` with a fitted slope would let an eff-len
mis-calibration hide as a "learned relationship." Adjacent nodes of the same field (gDNA genome-wide, nascent
within one transcript) share **one** underlying rate, and the eff-len redesign makes that **factor-1-under-
uniform** by construction. So the imputation mean is `μ_dst = ρ_src` (identity), with no slope parameter. A slope
that wasn't ≈1 would be an eff-len **bug to fix**, never biology to absorb — the parsimonious, honest choice.

**RNA: pool the strands for the curve, apply per-strand.** Biological dispersion of nascent RNA is a property of
*transcript structure*, not of which strand carries it — the `+` and `−` processes are symmetric. So pool both
strands' points into **one** RNA `σ²_bio(μ)` curve (more data, no asymmetry bug), and apply it per-strand using
each strand's own count and density. gDNA has its own curve (one component). This also makes the historical
"neg-fit reads the pos density" class of typo structurally impossible.

---

## 6. COMPONENT DEFINITIONS (what `C`, `L`, `ρ` are)

Per node, per component (gDNA, RNA⁺, RNA⁻):

- **count** `C = f_component · N_unspliced` — the component's expected fragment count (the Poisson sampling unit
  is the *fragment*). For a **boundary** node's RNA component, add the one-sided motif-stranded spliced count:
  `C_rna± = f_pos/neg · N_unspliced + N_spliced_±` (the spliced mass is mature RNA the boundary owns).
- **effective length** `L` — the component's FL-marginal effective support (region `gdna_region_eff_len` /
  per-side density length; the RNA twin for RNA; `effective_length.py`).
- **density** `ρ = C / L`. The sampling variance `Var_samp(ρ) = C/L² = ρ/L` uses the **count** `C` (the discrete
  sampling unit), *not* the fractional accumulator mass — mass ≈ count for contained fragments but the Poisson
  power is in the fragment count.

---

## 7. IMPLEMENTATION (Step-2 deltas)

Reuses the kept machinery (`variance_model.pair_imputation_points`, `MonotoneVarMean`). The deltas:

1. **`pair_imputation_points` carries `(C, L)` per node**, not just densities, so it can form the per-point
   Poisson offset `V_p = C_src/L_src² + C_dst/L_dst²` (= `ρ_src/L_src + ρ_dst/L_dst`). Point: `(μ = ρ_dst,
   R² = (ρ_dst − ρ_src)², V_p)`.
2. **`MonotoneVarMean.fit` accepts a per-point offset** `V_p`: it fits `σ²_bio(μ)` with the model
   `curve(μ) + V_p` against `R²` (the Poisson-offset variant), composing with the existing Jensen dof correction.
3. **Apply**: `τ_impute = 1/(σ²_bio(μ_dst) + ρ_src/L_src)`, then convert to the fraction prior via the
   `(eff/mass)²` jacobian. The mean is `μ_dst = ρ_src` (identity).
4. **gDNA + RNA** share the builder/fitter; RNA pools `±` for the curve (§5) and adds the boundary spliced count
   (§6); gDNA propagates on the strand-deconvolved / signature-known density.

---

## 8. MAGIC-NUMBER LEDGER

**Zero new constants.** `Var_samp = C/L²` is exact (Poisson). The imputation mean is the identity (slope = 1, by
the factor-1 eff-len construction). The only knobs are the existing GCV-λ and lattice resolution of
`MonotoneVarMean`. The biological-dispersion floor is structural (`σ²_bio ≥ 0`).

**NB-later (overdispersion) is a graceful extension.** If counts are negative-binomial, the true sampling
variance is `C + C²/r > C`; the Poisson floor `C` then under-states it, and the learned excess `σ²_bio` absorbs
the difference (it is still "variance beyond the computed floor"). So the model is correct-in-direction today and
sharpens when the apply-time floor is upgraded from `C` to `C + C²/r` — no structural change.

---

## 9. WORKED EXAMPLE (the contamination, removed)

Two adjacent nodes, true `λ = 5` both (uniform field ⇒ `σ²_bio = 0`). Node A: `10/2 bp`, `Var_samp = 2.5`. Node
B: `1000/200 bp`, `Var_samp = 0.025`.

- **Density-only (pre-Step-2):** the typical residual `(ρ_A − ρ_B)² ≈ Var_samp,A + Var_samp,B = 2.525`. The curve
  reads **2.5 at density 5** — high variance on a uniform field. *Wrong:* it has mistaken node A's count noise
  for biological dispersion, and will under-trust every imputation near that density.
- **Poisson-offset (Step-2):** the model fits `σ²_bio(5) + 2.525 ≈ 2.525` ⇒ `σ²_bio(5) ≈ 0`. *Correct:* the field
  is uniform; the dispersion is zero; the 2.525 is fully explained by the computed sampling floor.
- **Apply:** impute A from B ⇒ `τ = 1/(0 + 0.025) = 40` (B is a reliable predictor); impute B from A ⇒
  `τ = 1/(0 + 2.5) = 0.4` (A is a noisy predictor). The count drives reliability, exactly as it should.
