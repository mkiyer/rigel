# Principled clue integration (design / for review)

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

- `strand_reliability(library)` — a continuous [0,1] measure of how identifiable the strand model is
  (κ's distance from ½ and/or the existing identifiability margin). Start: `|2κ−1|`; refine from the
  strand model's own confidence statistics (the pipeline already computes an identifiability margin).
  **Not a tunable constant** — a deterministic function of the fitted strand posterior.
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

**Open detail:** the exact `strand_reliability(library)` form — `|2κ−1|` vs an identifiability-margin
based measure (accounting for n_spliced_obs). Decide from the stress scenario across the SS spectrum;
must be monotone, ∈[0,1], and free of tunable constants.

