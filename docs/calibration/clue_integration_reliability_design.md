# Principled clue integration (design / for review)

> **SUPERSEDED — exploration log only.** The implementation foundation is now
> `count_overdispersion_integration_plan.md` (clean, phased). This file records the reasoning trail
> (the vetoed count-reliability §1–§7, the superseded strand-reliability weight §8, and the §9
> count-overdispersion result that became the chosen solution). Kept for provenance; do not implement
> from it.

**Status:** §1–§7 below proposed a *count-reliability* metric — **VETOED** (see banner). The chosen
direction is **§8: weight the count by STRAND reliability** (the count yields to strand to the extent
strand is trustworthy), validated at both SS endpoints. 2026-06-09.

> **VETO (count reliability) — accepted.** Count reliability is **local**, not global: the probe
> design + transcript geometry determine, per region, whether the boundary-crossing density over- or
> under-estimates the interior (e.g. a 5 kb single exon with a *central* probe → the boundaries are
> *depleted* relative to the enriched centre; another probe layout → enriched). We cannot know the
> probe design, so the per-region count bias is **unestimable**. A global count-reliability factor
> (§3) is therefore wrong. "Statistically precise" (thousands of fragments) ≠ "accurate."

> **CHOSEN DIRECTION — weight the count by STRAND reliability (§8).** Unlike *count* reliability,
> *strand* reliability is global and estimable (κ's distance from ½, identifiability, n_obs). Scale the
> imputed-region count concentration by `(1 − strand_trust)`: where strand is reliable it governs and
> the count fades; where strand is weak/inapplicable the count carries. Continuous in SS (no cliff),
> and it **subsumes** the AMBIG special-case (strand not observable ⇒ strand_trust 0 ⇒ count carries)
> AND fixes the weak-SS single-strand-exon floor. Validated: ss0.50 single-strand exons 0.60→~0.92
> (count kept), ss0.99 strand governs (0.77, leak stays closed). See §8.

---

## 1. The critique (accepted)

The deconv combines two clues per region: a **count prior** `Beta(mean, concentration)` and a
**strand likelihood** (Beta-Binomial over sense/antisense). Today:
- the **strand** weight is already correct — the Beta-Binomial is *flat* at κ=½ (uninformative) and
  *sharp* as κ→0/1, so its influence scales continuously with strand specificity, automatically;
- the **count** concentration is set to the **raw observed gDNA count** (often thousands). That is a
  **statistical-only** precision — it ignores the count clue's **systematic** error (the within-exon
  enrichment gradient: the imputed crossing density is ~½ the interior under capture). So the count
  clue is handed confidence it has not earned, dominates the joint, and **degrades** the result
  (the captured-exon leak: count 0.40 overrides strand 0.77).

The current mitigation is a **binary**: zero the count concentration for single-strand exons (defer
to strand), keep it for AMBIG (strand absent). This is "okay but sloppy" because:
- **strand specificity is a spectrum** — at *weak* SS the strand clue is itself unreliable, so zeroing
  an *accurate* count (e.g. a no-capture library, where the imputed density is unbiased) throws away
  the only good quantitative signal;
- the AMBIG-vs-exon split is really a proxy for "is a better clue available," but it ignores *how
  reliable* that better clue is.

**Root issue: our judgment of each clue's accuracy is poor.** We weight the count by its *count*, not
its *accuracy*.

## 2. Principle

Treat each clue as a measurement of the region's gDNA fraction with precision

```
precision = statistical_precision × reliability
```

The strand likelihood already embodies this (statistical = sense+antisense count; reliability =
κ-contrast, via the Beta-Binomial). The count prior must too: its concentration should be the
observed gDNA count **discounted by the count clue's reliability** — the degree to which the imputed
density reflects the true in-region density. Then the **ordinary Bayesian product**
`count_prior × strand_likelihood` handles the whole spectrum *automatically*:
- high SS → strand sharp → governs (count's earned weight is small);
- low SS → strand flat → the count governs **to the extent it is reliable**;
- AMBIG → strand absent → the count carries, again weighted by its reliability.

No region categories, no hard zero. The only new quantity is the **count reliability**.

## 3. Estimating count reliability from the data (the key idea)

The count clue can be *checked* exactly where a trustworthy reference exists: **single-strand exons at
identifiable strand specificity**, where the strand clue gives a reliable gDNA fraction. Compare the
count-clue fraction to the strand-clue fraction across those regions:
- they **agree** (no capture, imputed density ≈ interior) → the count clue is reliable → it should
  carry full weight everywhere it is the only signal (AMBIG, low-SS regions);
- they **disagree** systematically (capture, imputed density < interior) → the count clue is
  unreliable → discount its concentration so strand governs where sharp, and accept the genuine floor
  where strand is absent.

This is a **per-library reliability** (capture is a library property), an empirical-Bayes
hyperparameter fit from the aggregate — the same spirit as the existing overdispersion fits. Concrete
forms to evaluate (pick from data, §6): the regression slope of `count_gf` on `strand_gf`; or
`reliability = 1/(1 + residual_var/expected_var)` from the count-vs-strand residuals on the reference
regions; the count concentration is then `statistical_count × reliability`.

It **generalizes the categorical fix**: that fix is the degenerate `reliability ∈ {0,1}` switched by
AMBIG-vs-exon; this *learns* reliability and lets the Bayesian product weight continuously.

> **Optional refinement (de-bias, not just down-weight):** the same strand-vs-count comparison yields
> the systematic factor `g = strand_density / count_density` on the reference regions (the
> ~uniform-per-library enrichment gradient). Correcting `count_density ← g · count_density` library-wide
> would make the count clue *accurate* (not merely down-weighted), so it carries truthfully at low SS
> and AMBIG. This is a strand-calibrated, density-only analog of the Phase-2 DNA-fraction estimator —
> and may make Phase 2 unnecessary. Higher value, higher risk; gate behind the §6 validation.

## 4. Predicted behavior (vs current categorical fix)

| regime | current (binary) | proposed (reliability) |
|---|---|---|
| high SS, single-strand exon | strand ✓ | strand ✓ (count's small earned weight is harmless) |
| **low SS, single-strand exon, no capture** | **defers to weak strand → floor** | **accurate count carries (reliability≈1) ✓** |
| low SS, single-strand exon, capture | floor | floor (count unreliable, strand weak — genuine) |
| AMBIG, no capture | count carries ✓ | count carries ✓ (reliability≈1) |
| AMBIG, capture | full (biased-low) count | down-weighted → ~0.5 floor (count unreliable; genuine) |

The clear win is **weak-SS + low-capture** (a real regime — some unstranded protocols), where the
count is accurate and should carry but the binary throws it away. Under capture the AMBIG case is a
genuine floor either way (no clue is reliable); the de-bias refinement (§3) is what would lift it.

**Measured confirmation (antisense stress scenario, no capture; `build_antisense_stress_scenario.py`):**

| region | ss0.99 cal/oracle | ss0.50 cal/oracle | ss0.50 abs err |
|---|---|---|---|
| AMBIG overlap (count kept) | 0.91 / 0.93 | 0.81–0.90 / 0.93–0.94 | 0.04–0.12 |
| **single-strand exon (count zeroed)** | 0.96 / 0.96 ✓ | **0.60 / 0.96** | **~0.36** |

At weak SS the categorical fix is **backwards**: the AMBIG regions (count kept) recover well, while
the single-strand exons **floor at ~0.6 (err 0.36)** because their *accurate* count was zeroed and
strand is blind. And **no binary resolves it** — keeping the count would re-break high-SS + capture
(the original leak), zeroing it breaks this. Only a reliability-weighted count (≈1 here since
no-capture/unbiased, ≪1 under capture) handles both. This is the concrete motivation to implement §3.

## 5. Alternative considered — one-sided count prior

The imputed crossing density is a **lower bound** on the interior gDNA density (capture only enriches
the interior; no-capture is exact). So a count prior of the form "gDNA fraction ≥ count_gf" (one-sided)
never pulls *down* against a higher strand estimate, and provides a floor when strand is weak. Simpler
(no reliability fit), but it trades the systematic *under*-estimate for a potential *over*-estimate
when the truth sits just above the bound, and the upper value is undetermined without strand. Weaker
than the reliability/de-bias approach; noted for completeness.

## 6. Open questions

1. **Reliability estimator form** — regression slope vs residual-variance ratio vs a likelihood-based
   precision; decide from the antisense stress scenario + the condition matrix (not a priori).
2. **Circularity** — strand is used both to estimate reliability/`g` AND in the per-region likelihood.
   It is a *global* hyperparameter from the aggregate (a region's ~1/N self-contribution is
   negligible — the same argument the existing density/overdispersion fits already rely on), and on
   AMBIG (where reliability is *used*) strand is absent (no double-use). Confirm empirically.
3. **Fully-unstranded libraries** (κ≈½ everywhere) — no reliable reference to estimate reliability →
   fall back (to what? the current density·eff_len, or a conservative low weight?). This is the
   genuine unidentifiable floor; ensure the fallback is rational, not catastrophic.
4. **Per-region vs global reliability** — start global (capture is library-wide); revisit if exons
   vary (probe density differs across exons).
5. **AMBIG + capture floor** — accept the ~0.5 floor, or pursue the §3 de-bias refinement to lift it?
6. **Interaction with the gDNA/RNA overdispersion fits** — those already consume the count clue on
   seed regions; verify the reliability discount does not destabilize them.

## 7. Recommendation / decision requested

The current categorical fix is **excellent at high SS** (the common case) and correct for AMBIG at
no/low capture; its weakness is weak-SS + low-capture and the AMBIG+capture floor. Proposed path:

1. **Build the antisense stress scenario** (this turn) with ample observable + single-strand
   deconvolvable regions for training, plus opposite-strand overlaps + gDNA, across SS × capture.
2. **Measure the current categorical fix on it** — quantify exactly where it floors.
3. If the floors matter for the target use, **implement the data-driven reliability** (and possibly
   the §3 de-bias); otherwise keep the categorical fix and document the floor.

**Decision requested:** (a) implement the reliability integration now, or validate-the-categorical-fix
first then decide; (b) reliability down-weight only, or also the strand-calibrated de-bias (§3).

---

## 8. CHOSEN DIRECTION — weight the count by strand reliability

**Principle.** Each clue is unbiased *or* biased. The strand clue is (approximately) unbiased but its
precision fades with strand specificity. The count clue is statistically precise but **locally biased
in an unknown way** under capture. So: trust the strand to the extent it is reliable, and let the
count carry only the *remainder* of the belief. Concretely, scale the count-prior concentration of
**imputed (non-count-observable) regions** by the strand model's unreliability:

```
strand_trust(region) = strand_observable(region) × strand_reliability(library)   # ∈ [0,1]
count_evidence(region) = density·eff_len × (count_observable ? 1 : (1 − strand_trust))
```

- `strand_reliability(library)` — the strand estimate's **precision** (inverse variance). The strand
  estimate is `ĝ = (s − κ)/(½ − κ)`, so `Var(ĝ) = Var(s)/(½−κ)²` and **precision ∝ n_eff · (2κ−1)²**
  — the **squared** strand contrast (the per-fragment Fisher information of a strand label about g is
  `(½−κ)²`). This is the inverse-variance / BLUE weight, **not a heuristic constant**, and it is what
  the Beta-Binomial likelihood's curvature already encodes. **Use `(2κ−1)²` (squared), not linear
  `|2κ−1|`** (an earlier error). The current workable form: `strand_trust = strand_observable ·
  (2κ−1)²` (optionally × an n_eff / identifiability factor).
- `strand_observable(region)` — already computed; 0 for AMBIG (no defined sense), so AMBIG → trust 0
  → count carries, with **no special-case branch** (it subsumes the prior AMBIG fix).
- count-observable regions (intron/intergenic) keep full concentration (direct gДНК observation).

**Why this is the right shape, not a band-aid.** The strand likelihood's *effective* sample size is
limited by its (real, fitted) overdispersion, so at high SS the count's raw statistical concentration
overpowers it even though strand is the trustworthy clue — that mismatch is the entire bug, and the
previous binary zeroing was a hard patch for it. Tying the count's *weight* to the strand's
*reliability* is the continuous, estimable expression of "defer to the unbiased clue when it is
trustworthy." It never consults the (unestimable) count reliability.

**Validated (antisense stress scenario, no capture):**

| | oracle | current (binary) | strand-reliability weighting |
|---|---|---|---|
| ss0.99 single-strand exon | 0.96 | 0.96 ✓ | 0.96 ✓ (strand governs) |
| ss0.99 AMBIG | 0.93 | 0.91 ✓ | 0.91 ✓ |
| **ss0.50 single-strand exon** | 0.96 | **0.60** ✗ | **~0.92** ✓ (count carries) |
| ss0.50 AMBIG | 0.93 | 0.81–0.90 | 0.81–0.90 |

**Remaining floor (honest):** low SS **+ capture** (strand blind AND count locally biased) and
AMBIG + capture — neither clue is trustworthy. This is the genuinely-unidentifiable corner; the
**Path-3 DNA-fraction** de-bias (a same-position ratio robust to the local probe gradient — exactly
the "probe in centre" case) is the deferred enhancement that could lift the multi-exon part of it.

**Open detail:** whether to fold an `n_eff`/identifiability factor into `strand_trust` alongside
`(2κ−1)²`. Decide from the stress scenario across the SS spectrum.

---

## 9. Going deeper — WHY the count is mis-specified (live thinking; NOT done)

**Status:** the `(2κ−1)²` strand-reliability weight (§8) is the best *workable* solution now. But it
patches a symptom. The deeper diagnosis: the count prior's **precision is wrong**, in two distinct
ways.

**Beta vs NB (orientation).** NB models a **count** (magnitude); Beta models a **proportion**. The
deconv unknown is the gДНК *fraction* `g ∈ [0,1]` → Beta. They are one story: if gДНК ~ Gamma(α₁) and
RNA ~ Gamma(α₂), the **total** is Gamma(α₁+α₂) (NB-like magnitude) and the **proportion** `g` is
`Beta(α₁,α₂)`, independent of the total. The overdispersion analogy is exact: *Beta-Binomial :
Binomial :: NB : Poisson*. So our **strand** likelihood (Beta-Binomial) already models strand
overdispersion; the **count** prior (plain Beta with concentration = raw count) does not.

**Defect A — statistical over-confidence (Poisson vs NB).** `count_concentration = density·eff_len`
= the raw expected gДНК count. That is the *Poisson* precision (precision = count). RNA-seq counts are
**NB-overdispersed**: e.g. 1234 counts carries a ~99% CI of ~823–2543, not ±105. So the count's
precision is over-stated by ~an order of magnitude. Meanwhile the strand precision IS correctly
overdispersion-limited (Beta-Binomial caps effective N at ~1/φ_strand). **This asymmetry — strand
overdispersion modeled, count overdispersion not — is why the count over-powers the strand**, and is
the real thing the §8 weight is compensating for.
- **Principled fix:** model the gДНК **count overdispersion** `φ_count` (fit symmetrically with the
  strand fit, from the count-observable regions), and set the Beta concentration to the
  overdispersion-limited effective N `N/(1+φ_count·N)` (→ ~1/φ_count for large N), not the raw N.
  This restores count↔strand symmetry → the joint model's precision-weighting becomes correct →
  much less leak, and likely shrinks or removes the §8 heuristic. (This is the count-side twin of the
  existing strand Beta-Binomial.)

**Defect B — systematic local bias (the capture gradient).** Even with an honest precision, the count
*mean* for an imputed exon is biased under capture (boundary crossing under-samples the enriched
interior; probe-placement-dependent, §veto). Crucially, `φ_count` fit from observable regions
(introns/intergenic, *uniformly* depleted under capture) does **not** see this — the gradient is the
error of imputing depleted→enriched, invisible to the observable-region variance. So overdispersion
modeling **bounds** B's harm (low precision ⇒ a biased mean pulls less) but does **not remove** it.
- **Only fix for B:** the **DNA-fraction (path 3)** — a same-position ratio that cancels the local
  probe coverage → an *unbiased* count mean. With B removed and A fixed, the joint model is fully
  coherent (correct precisions, no patch, no κ-circularity).

**Synthesis / roadmap (current thinking):**
1. **Now (workable):** §8 strand-reliability weight with `(2κ−1)²`.
2. **Principled, defect A:** model count overdispersion → honest count precision; re-measure whether
   the §8 weight is still needed. *This is the next thing to investigate — it may be the real answer.*
3. **Principled, defect B:** DNA-fraction (path 3) → unbiased count mean → joint model coherent.

Open question to resolve by experiment: with count overdispersion modeled (defect A fixed), does the
joint model self-weight correctly enough that the §8 strand-reliability heuristic becomes unnecessary
(at least off-capture)? Measure φ_count on the benchmark + stress scenarios and re-run the deconv.

### 9.1 RESULT — modeling the boundary count overdispersion is the principled fix (supersedes §8)

Measured (DESeq2-style common-dispersion MoM, `CV²(ρ) − mean(1/N)`; gDNA counts from
count-observable **regions** and **exon-adjacent observable boundaries** — the imputation source):

| condition | region-seed α | **boundary-seed α** | boundary eff-N (1/α) |
|---|---|---|---|
| gdna400 / ss0.99 / **capture** | 0.038 | **1.21** | ~1 |
| antisense stress / **no capture** | ~0 | **0.0002** | ~5000 |

The **boundary** overdispersion is huge under capture (enrichment heterogeneity) and negligible
without — a measurable, capture-sensitive reliability signal (the intron/intergenic *regions* miss
it: uniformly depleted ⇒ α≈0.04 even under capture, which is why boundaries are essential, as
predicted). Setting the count-prior concentration to the **overdispersion-limited effective count**
`count_evidence = N/(1 + α·N)` (→ ~1/α for large N) then self-weights the joint deconv correctly,
**with no strand-reliability weight, no zeroing, and no AMBIG special-case**:

| condition | oracle gf | deconv gf with `N/(1+αN)` |
|---|---|---|
| capture / ss0.99 | 0.804 | **0.795** (α huge ⇒ count fades ⇒ strand governs ⇒ leak closed) |
| no-capture / ss0.50 | 0.963 | **0.883** (α tiny ⇒ count carries ⇒ the weak-SS floor 0.60→0.88) |

**This is the chosen solution — it supersedes §8 (the strand-reliability weight) and the categorical
zeroing.** It is principled (the NB count overdispersion, the count-side twin of the strand
Beta-Binomial), data-driven (α measured from the seeds), continuous (no cliff), and it feeds the
joint model an *honest* count precision so its inverse-variance combination is automatic — dissolving
both the over-weighting and the κ-circularity. Defect B (the systematic gradient *bias*) is still not
removed — only bounded by the now-tiny effective count under capture — so the DNA-fraction (path 3)
remains the eventual enhancement; but density-only is now coherent and self-weighting.

**Implementation plan:**
1. Fit the gDNA **count overdispersion** `α` (MoM → optionally the DESeq2 mean-trend + shrinkage if
   the per-seed estimates prove noisy) from count-observable regions + **exon-adjacent observable
   boundaries**. A new calibration step alongside the strand-overdispersion fits.
2. `count_evidence = N / (1 + α·N)` for the supporting count `N` (= `density·eff_len`). Use
   `α_boundary` for imputed (exon/AMBIG) regions, `α_region` for count-observable regions.
3. **Remove** the categorical `defer_to_strand` zeroing (density_model.py) — the overdispersion-honest
   concentration replaces it; AMBIG needs no special-case.
4. (RNA twin, later) the same can model RNA count overdispersion from spliced-fragment counts at
   boundaries, if an analogous reliability is wanted for the RNA side.
5. Validate across the SS × capture matrix + the antisense stress; golden regen; end-to-end benchmark.

