Below is the audit + design. I'm framing this as a planning answer; I haven't changed any code.

## 1. Audit: what we have today

**Per-region "raw" estimates** — already correct and elegant:
- density_global.py `l_eff_contained(spans, gdna_fl)` → FL-PMF-weighted contained effective length. Single source of truth, delegates to `FragmentLengthModel.compute_all_transcript_eff_lens` (vectorised C-path, integer-bin spans). This is the right primitive; it's used wherever we need "fragment-start opportunities in a contained window of size L under FL distribution f".
- _exposure.py `fractional_boundary_side_exposure(lengths_bp, gdna_fl)` → the analogue for one side of a junction (formula: $\sum_l h(l)\cdot \min((l-1)/2,\,S/2)$). Available but not yet wired into the density estimator.

**Pooled global density** — current ceiling:
- density_global.py `_contained_density_for_mask` does `rho = sum(counts) / sum(eff)` per type (INTERGENIC, INTRON). One scalar per type. *Implicit model:* `N_R ~ Poisson(ρ · L_eff(R))`, slope-1 in log-log, no per-L dispersion.
- _kappa.py adds NB MoM overdispersion `κ` per type: `N_R ~ NB(ρ·L_eff(R), κ)`. Constant κ across L.
- density_loco.py `shrink_to_loco(N, L, ρ_global, κ)` = `(N + κρ) / (L + κ)`. Closed-form NB EB posterior mean. Clean.

**Gaps relative to what you want:**
1. Pooled ρ collapses the L→count relationship to a single constant. We cannot answer "for an exon of L=500bp, what is the 95th percentile of expected gDNA counts?" We only know the mean.
2. NB κ is type-global, not L-dependent. Real overdispersion typically grows with L (GC bias, repeats, copy-number).
3. Boundary mass is not yet folded into the density numerator/denominator at all. density_global.py only consumes `contained_unspliced_total`. The boundary channels are sitting unused even though `fractional_boundary_side_exposure` exists.
4. No "predict-and-subtract" mechanism: EXON regions get no gDNA estimate today (the orchestrator returns `_empty_density("EXON-INTRON")`).

## 2. Modeling choices — critique of quantile regression on `log N ~ log L`

Pros of QR:
- Distribution-free; gives quantile bands directly.
- Easy to fit with the `quantreg` formulation or a small custom solver.

Cons that bite us here:
- **Counts are discrete and zero-heavy** in small intergenic regions. `log(N)` is undefined; you need `log(N + c)` and the choice of `c` materially shifts the lower quantiles.
- **Heteroscedasticity is L-dependent in a known way**: under Poisson/NB, `Var(N|L) = μ(L) + μ(L)²/κ`. QR ignores this and reweights uniformly across L, which underweights the high-L regions where you most need the gDNA estimate (exons sit at small-to-medium L, but bias info is in long intergenic).
- **Quantile crossing** at small N (the 0.95 prediction can fall below the 0.90 prediction). Solvable but ugly.
- **Extrapolation**: a typical exon (a few hundred bp) is well within the L range of intergenic regions, so this is fine; but very long exons / very short exons need a model with a sane functional form, not a free-form quantile fit.

## 3. Recommended backbone: **NB-GLM with L-dependent mean and dispersion**, with an empirical quantile-calibration layer

The unification is much cleaner than QR and fits naturally on top of what we already have:

$$
N_R \mid L_R \sim \mathrm{NB}\bigl(\mu(L_R),\ \kappa(L_R)\bigr),\qquad
\log \mu(L) = \beta_0 + \beta_1 \log L,\qquad
\log \kappa(L) = \gamma_0 + \gamma_1 \log L.
$$

- **β₁ = 1, γ₁ = 0** recovers exactly today's behavior (pooled ρ, constant κ). So this is a strict generalization.
- **β₁ free** captures sub/superlinear scaling (mappability and GC effects shorten effective genome at large L).
- **κ(L)** captures the empirical fact that long regions are more variable.
- **Quantile interface**: `q_α(L) = NB_CDF⁻¹(α; μ(L), κ(L))`. Monotone in α by construction, no crossing.
- **Discrete-aware**: handles zeros without pseudocounts.
- **Few parameters**: 4 scalars per region class. Robustly fittable from a few thousand regions.

Fit by MLE (IRLS for the mean, MoM or profile-likelihood for the dispersion-vs-L slope). All operations vectorize over the intergenic mask we already build.

Add a non-parametric **empirical quantile calibration** on top: bin held-out intergenic regions by `log L` (10–20 bins), compute the empirical α-quantile of `N_R / q̂_α(L_R)`, and rescale. This catches systematic model misspecification without overfitting. Two-line implementation.

## 4. Region-class strategy

| Class | Role in model fit | Role in inference |
|---|---|---|
| INTERGENIC | training data (pure gDNA, large support) | reference baseline |
| INTRON-ONLY | training data (gDNA + sparse nRNA → still mostly gDNA) | optional second baseline; compare ρ̂ to intergenic to detect bias |
| EXON-only / EXON-AMBIG | **not** training data | query: predict `q_α(L_exon)` gDNA expectation, subtract from observed → RNA upper bound at `1−α` confidence |
| EXON-INTRON boundary regions | partial-mass training data (their boundary-side counts come from the same gDNA process as intron interior) | covered by the boundary handling below |

Two separate models are reasonable: a `gDNA_intergenic` model and a `gDNA_intronic` model. Compare their predictions on small EXON-INTRON regions to detect compartment-specific bias. For unstranded data both collapse pos+neg before fitting.

## 5. Boundary fragments — first-class observations, not a special case

The cleanest unification: for every region R, build a **single (count, exposure) pair**:

$$
N_R \;=\; n^{\mathrm{contained}}_R \;+\; m^{\mathrm{boundary\_left}}_R + m^{\mathrm{boundary\_right}}_R
$$
$$
L^{\mathrm{eff}}_R \;=\; \ell^{\mathrm{contained}}_R(\mathrm{span}_R) \;+\; \ell^{\mathrm{boundary\_left}}_R \;+\; \ell^{\mathrm{boundary\_right}}_R
$$

where:
- The boundary `m` values come from `PayloadArrays.boundary_*_unspliced_total` (sum of fractional bp-shares deposited from junction-crossing fragments).
- The boundary `ℓ` values come from `fractional_boundary_side_exposure(side_lengths, gdna_fl)` evaluated per side — i.e. the FL-weighted bp-share contribution a single fragment is *expected* to deposit in this side, given the side reaches some distance into R.

Then the NB-GLM is fit on `(N_R, L^eff_R)` exactly as if all evidence were contained. This is a real improvement over today because:
- Boundary mass is bp-weighted (matches our fractional accumulator).
- Tiny EXON-INTRON regions, which today contribute negligibly because they have small `contained` mass, now contribute proportionally to their actual junction-flux opportunity.
- The "asymmetry" between intergenic (mostly contained) and exon-intron (mostly boundary) is handled by one consistent (count, exposure) formula. No special EXON projection logic needed.

Constraint: the boundary side length used in `fractional_boundary_side_exposure` is the bp-length of the neighbor region (or the distance to the next edge, capped at the FL support). That's a per-junction property and is exactly what the edge-centric view from the previous turn naturally exposes; for now compute it from adjacent `RegionArrays.end/start` per ref-CSR slice.

## 6. Querying "how much of this exon is gDNA?"

For an EXON region with `L^eff` (as above) and observed total `N`:

1. Pick a confidence level α (e.g. 0.95).
2. Compute `μ̂_g(L^eff)` from the fitted gDNA NB-GLM and `κ̂_g(L^eff)`.
3. `gDNA_upper(α) = NB_CDF⁻¹(α; μ̂_g, κ̂_g)`.
4. `RNA_lower = max(0, N − gDNA_upper(α))`.

This gives a *conservative* RNA floor. For the Bayesian prior, pass `μ̂_g(L^eff)` (the predictive mean) as a pseudocount and `κ̂_g(L^eff)` as the prior strength — this finally lets the gDNA prior carry **uncertainty that scales with L**, matching what you wrote in the "gdna prior uncertainty" TODO.

## 7. Concrete implementation phases

Each phase is a self-contained, testable slice. None of them disturbs the current pooled estimator until the last one swaps it in.

**Phase D1 — unified (count, exposure) per region, still pooled.**
- Add a helper in density_global.py: `compute_unified_counts_exposure(region_arrays, payload_arrays, *, gdna_fl)` returning `(N_R, L_eff_R)` arrays as defined in §5.
- Currently this is only "what we'd feed to the model"; pooled ρ' = `sum(N)/sum(L_eff)` should equal today's pooled ρ within rounding when boundary mass is small. Diff against the existing INTERGENIC/INTRON ρ in a test.
- Side-effect win: today's boundary mass starts flowing into ρ.

**Phase D2 — NB-GLM fit.**
- New module `src/rigel/calibration/density_model.py` with `NBDensityModel` (frozen dataclass: β₀, β₁, γ₀, γ₁ + diagnostics) and `fit_nb_density(N, L_eff)` (IRLS for mean, profile likelihood for dispersion slope).
- Constrain β₁ ∈ [0.5, 1.5] and γ₁ ∈ [−0.5, 1.5] to keep the fit stable on small samples; fall back to today's pooled ρ + per-type κ when n_regions < ~50 or when the constrained MLE hits a bound.
- Implement `model.predict_mean(L)`, `model.predict_dispersion(L)`, `model.quantile(L, alpha)`.

**Phase D3 — empirical quantile calibration.**
- 10–20 log-L bins; held-out intergenic; rescale predicted quantiles by observed/predicted ratio. Persist as a small lookup on the model.

**Phase D4 — wire into orchestrator.**
- `compute_global_densities` returns `GlobalDensityTable` extended with `intergenic_model: NBDensityModel`, `intron_model: NBDensityModel`, and a thin `predict_gdna(L, alpha) -> (mean, quantile)` accessor.
- Keep the existing pooled fields for backward-compatible diagnostics; mark them as derived/legacy in the summary dict.

**Phase D5 — consumers.**
- `locus_prior` and the EXON projection in `ExonCompositeDensity` query `predict_gdna(L)` per EXON region. Subtract from contained EXON counts → initial RNA density.
- Replace the constant `gdna_prior_count` in the Bayesian EM with `μ̂_g(L^eff)` as pseudocount and `κ̂_g(L^eff)` as concentration — directly addresses the "gdna prior uncertainty" item.

**Phase D6 (later, optional) — covariate extension.**
- The same NB-GLM accepts extra log-linear covariates: GC content, mappability fraction, signature class. Drop-in additions without re-architecting. Mappability covariate also addresses the "mappability-corrected effective length" TODO without changing `L_eff` semantics — the model absorbs it as a multiplicative term.

## 8. Alternatives I considered and why I'm not recommending them

- **Pure non-parametric quantile lookup with isotonic smoothing** (bin by L, take empirical q_α per bin, monotonise). Simplest possible; respects data; no model assumptions. *Worth keeping as a fallback path* for small samples or as a sanity-check overlay on the NB-GLM. Not a replacement: it doesn't carry dispersion through to the EM prior cleanly.
- **GAM with NB family** (smooth `s(log L)`). Strictly more flexible than the log-linear NB-GLM, but adds a dependency and is hard to make deterministic and fast. Defer until a real motivating dataset appears.
- **Bayesian hierarchical** (per-locus shrinkage to a global NB prior). We already have this via `shrink_to_loco`; the NB-GLM is the per-L generalization. The hierarchical layer can sit unchanged on top.
- **Quantile regression on `log(N+0.5)`** (your original proposal). Workable but has the discreteness/crossing/heteroscedasticity issues above and doesn't compose with the NB-EM prior. If you really want quantile-regression semantics, use the NB-GLM's CDF inverse — same interface, better statistics.

## 9. Open questions to lock down before D2

1. **Two models or one with a class indicator?** I'd start with two (intergenic, intron-only) so we can *compare* and detect intron contamination. Merge later if they agree.
2. **Spliced fragments in the gDNA model?** Today's `contained_unspliced_total` excludes spliced. Boundary unspliced totals are included. Confirm we want the same exclusion for boundary contribution (yes, I think; spliced fragments are RNA evidence not gDNA evidence).
3. **Strand collapse for the unstranded path?** Sum pos+neg into `N_R`. Strand-balance kappa from §current code stays as a separate diagnostic / regularizer.
4. **How to handle the EXON-INTRON intermediate class?** Two viable routes: (a) include it in the intron training set with a soft down-weight by exon-bit fraction; (b) treat it as a third class with its own model. I'd start with (a) (simpler, more data).

Want me to go ahead and stub Phase D1 (the unified counts+exposure helper + a regression test against today's pooled ρ)? That's a small, safe first step that immediately gets boundary mass flowing into the density numerator without touching the model.