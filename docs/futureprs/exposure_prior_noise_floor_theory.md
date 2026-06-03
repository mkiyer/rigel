# Setting the Exposure Prior — A Noise-Floor / Empirical-Null Theory

**Companion to** `exposure_prior_findings_report.md`.
**Question:** how do we set the regularizing prior on the exposure dispersion `φ` so that
it **adapts to data sparsity automatically**, with *no hand-tuned constant*, gracefully
across the full continuum from a 30-region toy to a 10⁶-region human capture panel?

**Conclusion up front.** There is a clean closed-form *foundation* — the dispersion
estimate has a sampling noise floor that **unifies both kinds of sparsity (region count
and coverage) into a single quantity** `Σ_r (ρ₀ L_r)²`. But the closed form is built on a
Poisson ideal and **provably under-shrinks**, because in the sparse regime the dominant
noise is not Poisson sampling but **deconvolution mis-allocation**, which has no tractable
closed form. The robust route therefore *seemed* to be an **empirical null**: measure what the
pipeline hallucinates as dispersion from genuinely-uniform input, as a function of
`(R, coverage)`, and treat *that* as the floor.

**Then we ran the empirical null (§7), and it refuted the prior-centric framing entirely.**
Region count does *not* reduce the floor (so `Σμ²` is the wrong statistic); the
deconvolution noise dominates Poisson by 30–1000×; and — decisively — the dispersion
magnitude does not predict the failure (the worst-recovering scenario has the *lowest*
hallucinated dispersion). A prior on `φ` is the wrong lever. The real blocker is the
`(ρ₀, ω)` identifiability bifurcation (findings report §6): **stabilize the gDNA density
`ρ₀` first**, and the prior reduces to a secondary coverage-based shrinkage. The sections
below record the a-priori theory (§2–6) and then the data that overturned it (§7), because
the negative result is itself the main finding.

---

## 1. What we are estimating, and what corrupts it

Estimand: `φ = Var(ω)`, the true region-to-region variance of the gDNA exposure. Under a
uniform library `φ = 0`; under hybrid capture `φ` is large and *real*.

Observable: `M_g,r`, the (soft, fractionally-allocated) gDNA mass in region `r`, with
working mean `μ_r = ρ₀ L_r` (the expected gDNA at uniform exposure `ω = 1`).

The empirical dispersion `φ̂` computed from `{M_g,r}` does **not** measure `φ`. It measures
`φ` **plus a noise floor** with three contributions:

- **(a) Poisson sampling** — finite gDNA counts fluctuate; `Var(M_g,r) ≥ μ_r` even at
  `ω ≡ 1`.
- **(b) Deconvolution mis-allocation** — `M_g,r` is a *soft* gDNA-vs-RNA call. When
  coverage is thin and the strand/length signal is weak, the E-step spreads gDNA
  unevenly across regions in a way that *looks* like exposure variation but is an
  artifact of the allocation, not the library.
- **(c) Empty-region collapse** — regions with `M_g,r ≈ 0` acquire a collapsed posterior
  `ω̂_r → 0`; in the Gamma-shape MLE these inflate the dispersion signal (a numerical
  artifact, already mitigated by evidence-weighting in the current code).

A principled prior must subtract **(a)+(b)+(c)**, keeping only the genuine `φ`.

---

## 2. The closed-form foundation: a noise floor that unifies both sparsities

Consider the ideal case where `M_g,r` is a clean Poisson count (ignore (b),(c) for now).
Test `H₀: φ = 0` against `φ > 0`. The locally-most-powerful (score / Dean–Lawless) statistic
for Poisson overdispersion is built from

```
G = Σ_r [ (M_g,r − μ_r)² − M_g,r ].
```

Under `H₀` each term has mean 0; under the alternative
`E[(M−μ)² − M] = Var(M) − μ = φ μ²`, so

```
E[G] = φ · Σ_r μ_r².
```

The null variance, using Poisson central moments
(`Var[(Y² − Y)] = 2μ²` for `Y = M − μ`), is

```
Var₀[G] = 2 · Σ_r μ_r².
```

Hence the method-of-moments estimator and its sampling error:

```
φ̂ = G / Σ_r μ_r² ,        SE(φ̂) = √( 2 / Σ_r μ_r² ) ,        T = φ̂ / SE(φ̂).
```

**This is the foundation, and it answers the user's question about the two sparsities
directly.** The reliability of `φ̂` is governed by the *single* quantity

```
Σ_r μ_r²  =  ρ₀² · Σ_r L_r²
```

which simultaneously encodes:

- **region sparsity** — the number of terms in the sum (few regions ⇒ small `Σμ²`);
- **coverage sparsity** — the magnitude of each `μ_r = ρ₀ L_r` (shallow gDNA ⇒ small `μ_r`
  ⇒ small `Σμ²`).

Both axes the user identified collapse into one scalar. No separate handling, no two
knobs — `Σμ²` *is* the data's information content about `φ`.

### 2.1 A parameter-free shrinkage from the foundation

Treat `φ̂ ~ N(φ, SE²)` and shrink toward 0 by empirical Bayes: estimate the prior variance
of the *true* `φ` from `φ̂`'s own excess over its noise, `τ̂² = max(0, φ̂² − SE²)`. The
posterior-mean shrinkage `τ̂²/(τ̂² + SE²)` gives the **positive-part estimator**

```
φ_shrunk = φ̂ · max( 0, 1 − SE²/φ̂² )  =  φ̂ · max( 0, 1 − 1/T² ).
```

It is **parameter-free**: the only constant is `1`, the null standard deviation of the
score statistic `T` (a unit-information / score-test threshold, not a tuning knob). Its
behavior is exactly the target:

| regime | `Σμ²` | `T` | `φ_shrunk` |
|---|---|---|---|
| dense (many regions, deep) | large | large | `≈ φ̂`  (trust the data; capture preserved) |
| sparse (few regions *or* shallow) | small | `≤ 1` | `0`  (`ω → 1`, uniform) |

The crossover is at `T = 1` — when the apparent dispersion is one sampling-SD above zero —
and is **continuous** (`φ_shrunk → 0` smoothly as `T → 1⁺`), so it introduces no
discontinuity of its own.

This replaces both heuristics in the current code: the Kish region-count `N_eff` (which
saw region sparsity but was blind to coverage) becomes `Σμ²` (which sees both), and the
hand-set `κ` disappears.

---

## 3. Why the closed form is not enough (the honest part)

The foundation assumes `M_g,r` is a clean Poisson count. It is not — it is a *soft
deconvolution output*, and contribution **(b)** dominates the sparse regime. Worked
example, sparse `nrna_dc` toy with a healthy `ρ₀ ≈ 0.04` and `L_r ≈ 500–1000`:

```
μ_r ≈ 20–40 ,   Σμ² ≈ 3·10⁴ ,   SE(φ̂) = √(2/3·10⁴) ≈ 0.008.
Observed φ̂ ≈ 5  (inflated by mis-allocation)  ⇒  T ≈ 600.
```

The Poisson floor declares `φ̂ = 5` overwhelmingly significant (`T ≈ 600`) and **shrinks it
not at all** — yet we *know* the true `φ = 0` and the `5` is a deconvolution artifact. The
mis-allocation creates genuine excess variance in `M_g,r` *beyond* Poisson, so the Poisson
floor correctly says "beyond Poisson" but wrongly attributes it to exposure rather than to
the allocation. (Worse, the floor's value depends on `ρ₀`, which itself collapses in the
sparse regime — §5 of the findings report — so the closed form is doubly fragile here.)

Contribution (b) has no clean closed form: the spurious dispersion comes from *systematic*
mis-allocation correlated with region features (a region whose RNA looks gDNA-like
over-calls gDNA), which is a model-dependent bias, not an i.i.d. noise term.

---

## 4. The robust route: an empirical null

Since the floor resists closed form, **measure it**. Run the *entire* pipeline on
synthetic libraries that are **genuinely uniform** (`ω ≡ 1`, no capture, no enrichment)
and record the dispersion `φ̂_null` it reports. By construction the true `φ = 0`, so
**whatever `φ̂_null` it returns is the total hallucinated floor** — Poisson + deconvolution
+ collapse, all sources at once, with their real interactions. This is exactly Efron's
*empirical null* (the null distribution of a statistic estimated from the data-generating
process itself, rather than assumed), here obtained by parametric bootstrap over the
simulator.

Characterize it as a law over the two sparsity axes:

```
φ̂_null = Φ₀(R, c)        R = region count,  c = a coverage measure (e.g. Σμ² or median gDNA/region).
```

Then the principled, parameter-free estimator is a floor subtraction / shrinkage against
the *measured* null:

```
φ_real = max( 0, φ̂_obs − Φ₀(R, c) )           (subtractive), or
φ_real = φ̂_obs · max( 0, 1 − Φ₀(R,c)/φ̂_obs )   (multiplicative, matching §2.1 with the
                                                 empirical floor in place of SE).
```

`Φ₀` is a derived function of measurable data properties — **not a constant, not a knob**.
It will reproduce the qualitative law of §2 (decreasing in both `R` and coverage) but with
the correct *magnitude*, because it includes the deconvolution term the closed form
misses. The expectation (to be confirmed in §6) is `Φ₀` decaying like
`a/√(Σμ²) + b·(allocation-ambiguity term)`, the first matching the Poisson SE and the
second capturing (b).

This is the solution we believe is both principled and robust: the closed-form theory
tells us the **right variables** (`Σμ²`, the unification of both sparsities) and the
**right shape** (shrink by signal-to-floor); the empirical null supplies the **calibrated
magnitude** that the deconvolution coupling puts out of analytic reach.

---

## 5. The remaining coupling (out of scope here, but must be tracked)

§6 of the findings report showed that driving `ω → 1` (which *any* correct prior does in
the sparse limit) trips a separate instability: the global gDNA density `ρ₀` collapses
because gDNA inside high-RNA regions becomes undetectable once exposure can no longer
"place" it. **No exposure prior — however principled — fixes this**, because it is a
property of the joint `(ρ₀, ω)` E/M, not of the dispersion estimate. A complete solution
needs the §2/§4 prior **and** an independent stabilization of `ρ₀` (e.g. estimating the
density from detectable low-RNA regions and refusing the circular collapse). The two are
separable and should be designed separately; this note covers only the prior.

---

## 6. Validation / derivation plan (the diagnostic sweep)

One experiment both *derives* `Φ₀` and *validates* the estimator:

1. **Null sweep — derive `Φ₀(R, c)`.** Uniform-gDNA scenarios over a grid of region count
   `R` (replicate gene blocks: 1, 2, 5, 10, 20, 50, 100) × coverage (gDNA fragment budget
   per region). Record `φ̂` (the empirical null), `Σμ²`, `ρ₀`. Fit `Φ₀` and check whether
   its sparse tail exceeds the Poisson `√(2/Σμ²)` (quantifying the deconvolution term).
2. **Recovery check — uniform.** With `Φ₀` in place, `φ_real` must be ≈0 and RNA recovery
   graceful across the whole grid (no sparse collapse, no dense over-shrink).
3. **Discrimination check — capture.** On capture scenarios `φ_real` must stay large and
   the captured `ω ≫ 1` must survive (the floor must not eat real enrichment). This is the
   crucial test that the empirical null subtracts *noise* without touching *signal*.
4. **Protocol robustness.** Repeat at different gDNA fractions / strand specificities to
   confirm `Φ₀(R, c)` generalizes (the whole point — no per-dataset tuning).

If `Φ₀` is stable and low-dimensional (a smooth function of `Σμ²`, perhaps with one
coverage-ambiguity covariate), we will have an adaptive prior with **zero free constants**
that degrades gracefully across every sparsity regime — the elegant solution the project
requires.

---

## 7. Empirical null results — the a-priori theory is refuted

The §6 sweep was run (uniform-gDNA, `κ=0` so the raw floor is visible, single seed;
`scripts/debug/exposure_null_sweep.py`). The data overturn the closed-form picture.

**Region sweep** (coverage-per-block held constant, `nfrag = K·2000`):

| K | R | gDNA/reg | φ_null | Σμ² | SE_pois | φ_null/SE | RNA rec |
|---|---|---|---|---|---|---|---|
| 1 | 33 | 5.97 | **0.95** | 2.0e3 | 0.031 | 30 | **47%** |
| 2 | 65 | 15.4 | 1.82 | 4.8e4 | 0.0065 | 283 | 104% |
| 5 | 161 | 7.96 | 1.68 | 2.2e4 | 0.0095 | 176 | 91% |
| 10 | 321 | 7.67 | 1.56 | 4.6e4 | 0.0066 | 236 | 92% |
| 20 | 641 | 8.88 | **1.94** | 1.6e5 | 0.0036 | 541 | **97%** |
| 50 | 1601 | 11.1 | 1.77 | 6.2e5 | 0.0018 | 990 | 116% |

**Coverage sweep** (`R` fixed at K=5; scale the fragment budget):

| nfrag | gDNA/reg | φ_null | RNA rec |
|---|---|---|---|
| 2500 | 2.21 | 1.55 | 106% |
| 5000 | 7.31 | 1.89 | 104% |
| 10000 | 7.96 | 1.68 | 92% |
| 20000 | 14.5 | 1.08 | 61% |
| 40000 | 43.9 | 0.78 | 113% |
| 80000 | 64.8 | **0.19** | 95% |

Four results, each falsifying part of §2:

1. **Region count does NOT reduce the floor.** `φ_null` stays ≈1.7 across R = 65→1601 — it
   does not average toward 0. So `Σμ²` is the **wrong sufficient statistic**: `Σμ²` grows
   ~300× across the region sweep while `φ_null` is flat. The two sparsity axes are *not*
   unified, and they do not act alike.
2. **The deconvolution floor dominates Poisson by 30–1000×** (`φ_null/SE_pois` grows with
   R). Contribution (b) is not a correction to the Poisson floor — it *is* the floor; the
   closed form of §2 is irrelevant at every scale tested.
3. **Coverage IS the axis that shrinks the floor.** `φ_null` falls 1.55→0.19 as gDNA/region
   rises 2→65. Per-region gDNA evidence — not region count — controls the hallucinated
   dispersion.
4. **The dispersion magnitude does not predict the failure.** K=1 has the *lowest* `φ_null`
   (0.95) yet the *worst* recovery (47%); K=20 has a *higher* `φ_null` (1.94) yet recovers
   97%. Downstream robustness is **region-count averaging** (the law of large numbers in
   the per-locus sums), which is independent of, and invisible to, the dispersion estimate.

### 7.1 What this means for the prior

A prior on `φ` — however principled — is the **wrong lever**. The quantity it would
regularize (`φ_null`) neither predicts the failure nor responds to the axis (R) that
actually drives downstream robustness. And shrinking `φ → 0` (`ω → 1`) does not rescue the
sparse case: it trips the `ρ₀` collapse (findings report §6). The sparse failure is the
**(ρ₀, ω) bifurcation**, which no exposure prior can escape from either side:

- `ω` free → a few regions over-attract gDNA (no averaging at small R) → gDNA over-called;
- `ω → 1` → `ρ₀` collapses → gDNA under-called.

**Revised conclusion / recommended pivot.** Stabilize `ρ₀` *first* (make `ω → 1` safe).
With a robust density, the sparse limit `ω → 1` gives the correct answer (uniform gDNA,
density held), and the prior question collapses to a mild **coverage-based** regularization
(region count self-handles via averaging). The elegant parameter-free prior we sought is
real but **secondary**; the `(ρ₀, ω)` identifiability is **primary**. The diagnostic the
user proposed did its job — it told us we were tuning the wrong knob.

*(Caveat: single seed per cell; `φ_null` carries scenario-to-scenario noise, but the four
qualitative results above hold across all 12 cells and both axes.)*

---

### Open questions for review

1. Is the **empirical-null** route (parametric-bootstrap calibration of the dispersion
   floor) sound here, or does estimating `Φ₀` from the same simulator that drives the tool
   risk circularity we should guard against?
2. Can contribution (b) (deconvolution mis-allocation dispersion) be given a closed form
   after all — e.g. via the per-fragment allocation variances `Σ a_i(1−a_i)` the E-step
   already computes — so the floor becomes analytic again?
3. Is `Σμ² = ρ₀² Σ L_r²` the right sufficient statistic for "information about `φ`," or
   should the per-region NB Fisher information `Σ μ_r²/(1+μ_r φ)²` (which saturates at high
   coverage) be used instead?
