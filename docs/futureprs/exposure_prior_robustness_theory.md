# Theory — Robust exposure-weight calibration: the Gamma–Poisson model, shrinkage, and adaptive regularization

**Purpose.** Formalize the per-region gDNA **exposure** model, show *why* it is brittle
on small data, evaluate the proposed remedy (a stronger, data-adaptive prior that
shrinks exposures toward uniformity), and recommend a path forward. Motivated by the
count/EM deep dive: the `nrna_dc g20` failures (Mechanism A — the gDNA effective-length
collapse → over-attractive EM gDNA → catastrophic RNA loss; Mechanism B — the
neg-control region glitch) both trace to **brittle exposure weights**, and that
brittleness is a **sparse-data artifact** (verified below).

## 1. The statistical model

The accumulator gives, per calibration region `r`, a soft-allocated gDNA mass `M_{g,r}`
over a physical length `L_r`. The generative model (doc 01 §4.2, doc 03 §4):

```
M_{g,r}  ~  Poisson( ρ₀ · ω_r · L_r )        # gDNA count, given local exposure ω_r
ω_r      ~  Gamma( s, s )                     # exposure prior: E[ω]=1, Var[ω]=1/s = φ
```

`ρ₀` is the library-wide gDNA density; `ω_r` is the region's multiplicative **exposure**
(1 = neutral/uniform, >1 enriched, <1 depleted — e.g. hybrid capture). The prior mean is
**1 (uniformity)**; its strength is the shape `s = 1/φ`, where `φ = exposure_dispersion`
is the variance of `ω`.

By Gamma–Poisson conjugacy the posterior is closed-form ([exposure.py](../../src/rigel/calibration/exposure.py)):

```
ω_r | M_{g,r}  ~  Gamma( s + M_{g,r},  s + ρ₀·L_r )
ω̂_r = E[ω_r | data] = (s + M_{g,r}) / (s + ρ₀·L_r)
```

## 2. Shrinkage — exposure is a convex pull toward uniformity

Writing the raw per-region estimate `m_r = M_{g,r}/(ρ₀·L_r)`, the posterior mean is a
**convex combination of the prior mean (1) and the data (m_r)**:

```
ω̂_r = [ s·1  +  (ρ₀·L_r)·m_r ] / [ s + ρ₀·L_r ]
```

- **prior weight** `s` (in pseudo-gDNA-count units);
- **data weight** `ρ₀·L_r` (the *expected* gDNA count in the region under uniformity).

So the shrinkage is already **per-region adaptive in depth**: a region with little
expected gDNA (`ρ₀·L_r ≪ s`) is pulled to `ω̂≈1`; a deep region (`ρ₀·L_r ≫ s`) follows
its data. The crossover is set entirely by the **prior strength `s`**. Everything below
is about how to choose `s` (equivalently `φ`).

## 3. How `s` is set today, and why it is brittle

`φ = 1/s` is re-estimated every M-step by the **Gamma-shape MLE** over the per-region
exposure posteriors (`update_exposure_dispersion`, [mstep.py](../../src/rigel/calibration/mstep.py)) — an empirical-Bayes fit of the dispersion from the spread of the `ω̂_r`.
Asymptotically (many regions) this is correct. On **sparse** data it fails in a
self-reinforcing way:

1. With few regions and/or few gDNA counts per region, the observed spread of `ω̂_r` is
   dominated by **sampling noise** (Poisson fluctuation in `M_{g,r}`) and by
   **deconvolution mis-allocation** (at low strand specificity the E-step splits a
   locus's reads bimodally — some regions get `M_g≈0`, others `M_g>0`), *not* by real
   exposure variation.
2. That inflated spread inflates the estimated `φ` ⇒ small `s` ⇒ a **weak prior**.
3. A weak prior lets a region with `M_{g,r}≈0` collapse: `ω̂ = s/(s+ρ₀L) → 0`. (Concretely
   in `g20_n70_s65`: `φ≈5.6 ⇒ s≈0.18`, while `ρ₀·L≈20`, so `ω̂ ≈ 0.18/20 ≈ 0.009`.)
4. The collapsed `ω̂` feeds `gdna_eff_len = Σ ω̂·L → 0`, which makes the EM's gDNA
   component pathologically over-attractive (Mechanism A), and feeds the neg-control
   region glitch (Mechanism B). The collapse increases the apparent dispersion → step 2.

**The crux:** the prior strength `s` is *estimated from the very data it must regularize*,
and on sparse data that estimate is biased toward "weak," which is the opposite of what
robustness requires.

## 4. Empirical verification — it is a sparse-data artifact

Scaling the `nrna_dc g20_n70_s65` scenario up (more gene-pairs + proportional fragments)
makes the failure dissolve:

| gene-pairs | regions | total-RNA recovered | gDNA |
|---|---|---|---|
| 1 (toy) | 33 | **20 %** (284/1404) | 1716 |
| 5 | 161 | **91 %** (6342/6991) | 3658 |
| 20 | 641 | **104 %** (28999/27943) | 11001 |

`exposure_dispersion` stayed ~4–5.6 at *every* scale — so the recovery is **not** the
dispersion adapting; it is the **law of large numbers**: with hundreds of regions a
single locus's ω-collapse is a harmless tail, whereas in the toy it is the *only*
expressed locus. On a human genome (10⁵–10⁶ regions) the system is in the robust regime;
the toy scenarios sit in the pathological one.

**Second verification — the exposure fix alone is sufficient (no downstream floor needed).**
Forcing a stronger prior (capping the dispersion `φ`) on the toy recovers the RNA *without
touching* `gdna_eff_len`: at `φ≈0.2` the gDNA effective length stops collapsing and RNA is
recovered. This confirms the exposure is the right lever, and that fixing it at the root
makes `gdna_eff_len` correct *as a consequence* (so the downstream eff-len floor discussed
below is unnecessary — see §6.5).

| φ cap | gdna_eff_len min | t1 (exp 576) | total RNA (exp 1404) |
|---|---|---|---|
| none (φ=5.6) | 67.6 | 251 | 284 |
| 0.2 | 3870 | 588 | 1564 |
| 0.05 | 7502 | 585 | 1787 (over) |

## 5. Evaluation of the proposed remedy

The proposal — *stronger prior when sparse, weaker when dense; shrink ω→1 under sparsity*
— is **correct, and it is the principled fix**, for two reasons:

- **Statistically:** `ω=1` (uniformity) is the right *default* precisely because, with
  too few gDNA counts, the data cannot distinguish real exposure variation from sampling
  noise. Shrinking to 1 there is not a hack — it is the honest posterior under a prior
  that says "assume uniform until the evidence says otherwise."
- **It targets the root** (the brittle `ω`), upstream of the Mechanism-A symptom: if `ω̂`
  stays ≈1 on the toy (where the gDNA *is* uniform), then `gdna_eff_len ≈ Σ L = physical
  span`, the EM gDNA component is no longer over-attractive, and RNA is preserved — the
  same effect as flooring `gdna_eff_len`, but achieved correctly at the source.
- **It does not harm dense/capture data:** with many high-count regions the prior relaxes
  (`s` small relative to `ρ₀·L`), and genuine enrichment (capture, `ω≫1`) is recovered.
  The whole point of *adaptive* strength is to be strong only where the data are weak.

The honest framing the proposal gets right: we do **not** expect good quantification on
sparse data; the goal is **graceful degradation** (shrink to uniform) instead of
catastrophic collapse.

## 6. The fix — regularize the dispersion toward uniformity, data-adaptively

The Gamma-shape MLE solves `log s − ψ(s) = c`, where
`c = mean_r( E[ω_r] − E[log ω_r] ) − 1 ≥ 0` is the **dispersion signal** (`c = 0 ⇔ no
dispersion ⇔ s = ∞ ⇔ ω≡1`). Make the estimate **shrink `c` toward 0** by `κ`
pseudo-uniform regions (each contributing `0` to `c`):

```
c_reg = c · N_eff / (N_eff + κ)
```

- `N_eff` = the **evidence-weighted** effective number of regions informing the dispersion —
  **not** the raw region count. Each region is weighted by its gDNA *information* (e.g. by the
  posterior shape `α_r = 1/φ + M_{g,r}`, which the Gamma posteriors already carry); empty /
  sparse regions contribute ~0. This is essential for **hybrid capture**: a panel may have
  10⁵ regions with only ~50 enriched, and the 10⁵−50 empty regions carry no evidence about
  the exposure distribution — using the raw count would falsely declare "dense → relax the
  prior" when the actual evidence is ~50 regions. Evidence-weighting makes the adaptivity
  honest: sparse evidence ⇒ small `N_eff` ⇒ strong shrinkage ⇒ `ω→1`; abundant evidence
  (deep uniform, *or* the 50 well-covered capture targets) ⇒ large `N_eff` ⇒ `c_reg → c`
  (real exposure recovered, including focal enrichment).
- `κ` = the **prior pseudo-count on the dispersion** — the one design parameter. It is the
  "transition scale": how much evidence we demand before trusting an apparent dispersion as
  real rather than noise. Equivalent Bayesian view: a hyperprior `π(φ)` concentrated near 0
  (uniformity), `κ` = its effective sample size; the posterior mean of `φ` is the shrunk
  estimate.

This is a **hierarchical (full-Bayes) exposure model**: level 1 the Poisson counts, level 2
the `Gamma(s,s)` exposure prior, level 3 a prior on `s` favoring uniformity. The MLE today
is the degenerate κ=0 case.

**Setting `κ` principledly (the open question).** Options, in increasing data-drivenness:
1. A **unit-information** prior (`κ` = 1 pseudo-region): weakest defensible; likely too weak
   to rescue the 33-region toy (shrinks `c` by only `N/(N+1)`).
2. `κ` tied to the **reliability of a variance estimate** — a variance needs O(tens) of
   well-determined observations, so `κ` ~ that scale makes the prior dominate until `N_eff`
   reaches it. Interpretable, but a chosen scale.
3. A **reference/Jeffreys hyperprior** on `φ` — parameter-free, but may be improper for a
   variance and needs care.

I lean toward (2) with `κ` interpreted as "regions needed to trust a dispersion," chosen and
documented explicitly (a meaningful statistical quantity, not an arbitrary cliff like the
current `0.01`/`100` bounds).

## 6.5 Converged design (after external review)

An external statistical review proposed three mechanisms: **(A)** a MAP penalty on `φ` via an
Exponential hyperprior `p(φ)=λe^{−λφ}`; **(B)** a count-dependent strength floor
`s_eff = 1/φ + γ₀/(1+Σ M_g)`; **(C)** a downstream `gdna_eff_len` floor (+ strand veto). We
converge as follows:

- **Adopt the core (A ≡ §6).** The MAP penalty and §6's pseudo-count are the *same family* — a
  weakly-informative prior on `φ` favoring uniformity, overwhelmed by data when abundant. We
  keep the **pseudo-count parametrization**: "`κ` regions of perfect uniformity" is
  interpretable (Q6) and drops into the existing 1-D Gamma-shape bisection by shrinking `c`.
  (Both are 1-D root-finds; the reviewer's solver-complexity concern does not apply.)
- **Drop (B).** `γ₀/(1+Σ M_g)` is an ad-hoc form of the same shrinkage with an arbitrary shape.
  The pseudo-count/MAP prior is its principled version, and the `N_eff`-weighting above already
  ties the prior strength to the gDNA evidence.
- **Drop (C) as a primary fix.** §4's second verification shows the exposure fix alone recovers
  the RNA — the eff-len floor is unnecessary. Worse, a flat floor at the physical span has a
  **capture tension**: under capture the gDNA *is* concentrated, so a small `gdna_eff_len` is
  *correct*, and flooring it would under-attract real gDNA. The exposure fix instead yields the
  correct `gdna_eff_len` in **both** regimes (≈physical when uniform-sparse since `ω→1`;
  concentrated when dense-capture since the prior relaxes and real `ω` is recovered) — one
  mechanism, both regimes. (A capture-*aware* eff-len guard is held only as a last resort if
  validation ever shows residual collapse.)
- **Adopt the reviewer's validation methodology** (§7) — its strongest contribution.

The result is **one** elegant mechanism — the adaptive, evidence-weighted pseudo-count prior —
not three, and it is the only single-knob option that stays correct across uniform and
focally-enriched (capture) data.

## 7. Recommendation / path forward

Sequence (agreed): finalize this design → build the test scenarios that gate it → implement.

1. **Build the missing test infrastructure first** — a real blind spot today, and it gates a
   safe fix (without it, "fixing" sparse robustness could silently break capture):
   - **sparse** (downsampled fragment counts) — exposures must degrade gracefully to `ω→1`
     with no EM collapse;
   - **dense-uniform** (the scale-up above) — no estimation bias / artificial variance;
   - **dense-capture** (focal enrichment: ~tens of targets at high `ω`, the rest depleted) —
     the prior must *relax* and recover the real enrichment, proving the `N_eff`-weighting.
2. **Implement the adaptive, evidence-weighted dispersion regularization** (§6/§6.5): shrink
   the dispersion signal `c` by `N_eff/(N_eff+κ)`, with `N_eff` weighted by per-region gDNA
   information. One interpretable parameter `κ` replaces the `floor`/`ceil` magic bounds; the
   current MLE is the degenerate `κ=0` case.
3. **Acceptance test (the reviewer's φ-vs-downsampling diagnostic).** Plot recovered `φ`
   against the global downsampling ratio: a correct adaptive prior shows a clean inflection
   where `φ` smoothly drops toward 0 as the library enters the sparse trap, and holds at the
   empirical dispersion when evidence is abundant. Pin `κ` by where that inflection should sit.
   Pass criteria across all three regimes: sparse → graceful `ω→1`; dense-uniform → robust;
   dense-capture → real `ω` recovered.
4. Re-examine Mechanism B (neg-control glitch) afterward — it should soften once `ω` is no
   longer brittle; if residual, letting clear unstranded strand-channel evidence **veto** a
   gDNA call (the reviewer's other suggestion) is the lever — not an eff-len floor.

**Bottom line.** The exposure model and its per-region shrinkage are sound; the failure is
that the prior *strength* is estimated weak-when-sparse, the opposite of robust. A
data-adaptive prior that shrinks exposures toward uniformity under sparsity — strong when
the evidence is thin, relaxing as it accumulates — is the principled remedy, fixes
Mechanism A at its root, and (with capture test coverage) preserves the dense/enriched
behavior the tool must also handle.
