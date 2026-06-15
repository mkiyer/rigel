# Calibration — the iterative Bayesian deconvolution (PLAN v2)

**Status:** fresh authoritative plan (2026-06-15). **Supersedes** `CALIBRATION_PLAN.md` (2026-06-13),
which predates the per-node Bayesian-hierarchy redesign, the monotone-P-spline variance model, and the
iterative bootstrap. Detailed companion docs (kept, referenced here, not duplicated):
`per_node_deconv_hierarchy_design.md` (the 3-tier per-node solve), `iterative_bootstrap_design.md` (the
circularity + bootstrap), `variance_model.py` docstring (the SCAM fitter). Where this and an older note
disagree, this wins.

---

## 0. TL;DR

Calibration deconvolves each genomic region's **unspliced** fragment mass into **gDNA vs RNA** (a per-node
pie `(f₊, f₋, f_g)`; `f_g` → the per-locus EM prior). The redesign is **one uniform iterative loop**:

```
init: assume ALL unspliced mass is gDNA
repeat:
   1. FIT the var~mean variance models (DIRECT + IMPUTATION) on the current per-node gDNA estimates
   2. SOLVE every node — a 3-tier Bayesian posterior: strand likelihood + odds-propagation (priority),
      count prior (fallback, precision = confidence·1/var), global gDNA prior (foundation)
   until the per-node masses stabilize (prediction error stops falling)
```

The init "all gDNA" is the conservative worst case; the variance models *report high uncertainty where
that's wrong* (RNA-rich regions), so the count yields and the strand/global do the work; iteration then
removes RNA and the gDNA residual becomes spatially smooth and boundary-predictable. No special first pass;
convergence and per-node confidence fall out of the boundary→region prediction error.

---

## 1. Goal & the problem

- **Goal:** per-node `f_g` (gDNA fraction of unspliced mass) feeding the per-locus EM, accurate across
  hybrid-capture (whole-transcriptome → tiny targeted panels), stranded → unstranded, zero-DNA → heavy-gDNA,
  cancer CNVs. Unsupervised (we don't know the capture/CNV landscape a priori).
- **The leak it must fix:** gDNA↔RNA misattribution (esp. the nascent sink) and the complex
  overlapping-opposite-strand (AMBIG) loci.
- **The circularity it resolves:** a trustworthy *imputation* variance model needs enriched-region training
  triplets; under capture those only exist after the enriched regions are solved; solving needs the model.
  Resolved by §5 (iterate, with a self-aware variance model + all-gDNA init).

## 2. Architecture — one uniform iterative loop

Every pass is identical (no strand-only Pass 0). A pass: **(1)** fit the variance models on the current
gDNA estimates; **(2)** solve every node's posterior with those variances; repeat until converged. Pass 0's
"current estimate" is the **all-gDNA** assumption (gDNA count = total unspliced mass). This commits
calibration to **iterative** (it was acyclic single-pass) — bounded (a few passes), seed/strand/spliced-
anchored, not open EM.

## 3. The per-node solve — the 3-tier Bayesian hierarchy

`log posterior(f₊,f₋,f_g) = L_strand(tilt) + L_odds(neighbours) − ½τ_count(f_g−μ_count)² − ½τ_global(f_g−μ_global)²`,
normalized over the 2-simplex lattice (one normalization; the count+global terms compose to one
precision-weighted prior). Details: `per_node_deconv_hierarchy_design.md` §3-4.

- **Tier 1 — strand likelihood + propagation (priority, unbiased).** The 3-component mixture on the
  observed plus-fraction (the *tilt* `p=½f_g+κf₊+(1−κ)f₋`, first-class at AMBIG) + the odds-propagation
  sum-product coupling per-strand `log(f_c/f_g)` along same-strand-exon chains. **Node-class prior:**
  Beta(½,½) Jeffreys reference at single-strand nodes (the fusion's proven prior; concentrates gDNA against
  the mixture's overdispersion spread), the **global gDNA prior** at AMBIG (where the Jeffreys U-shape would
  vertex-push a phantom).
- **Tier 2 — count prior (fallback).** Gaussian pull of `f_g` toward `count_gdna_frac`, with **per-node
  precision `τ_count = confidence(μ)·1/var(μ)`** from the variance model (§4). It *yields where the strand
  is informative* (it's a prior, not co-equal) and *governs where the strand is silent* — and yields where
  the variance model is *extrapolating* (confidence→0). Replaces the magic `β` (`count_trust_beta`).
- **Tier 3 — global gDNA prior (foundation).** `μ_global = clip(ρ_global·eff_len/mass, 0, 1)`,
  `σ²_global = var~mean(ρ_global)·(eff_len/mass)²` (§4). Catches the truly-unanchored AMBIG and is the
  zero-DNA fix (seeds agree ρ_global≈0 ⇒ tight prior at f_g≈0 ⇒ no phantom). `ρ_global` is the pooled
  observable-node density, PRE-solve (distinct from `derive.gdna_density_global`, POST-solve).

The strand>count>global hierarchy is **emergent from the precisions** — no hand-tuned weight.

## 4. The variance models — `variance_model.py`

One consolidated module: the SCAM fitter + the DIRECT and IMPUTATION builders + the confidence guard.

- **The fitter — `MonotoneVarMean` (monotone-increasing P-spline / SCAM).** `log var ~ monotone spline in
  log mean` via the reparameterization `β=cumsum([α₀,exp(α₁..)])` (monotone for any α), penalized
  (GCV-λ) weighted least squares with bisquare-IRLS robustness. Smooth + monotone-by-construction +
  fringe-flexible — validated to beat LOESS (non-monotone fringe) / isotonic (staircase) / power-law
  (rigid). Pure numpy+scipy (no new deps); the fit is once-per-pass over ~10³–10⁴ points, sub-second — **no
  C++ port needed.**
- **Confidence (the bootstrap lynchpin).** `MonotoneVarMean.confidence(μ)` = local training-data density
  (`n_eff/(n_eff+1)`, kernel-weighted at the knot-spacing bandwidth) → 0 on extrapolation. The count
  precision `confidence·1/var` ⇒ an untrained mean regime *cannot* give false certainty (the count yields
  to strand/global there). This is what makes the all-gDNA Pass-0 safe automatically.
- **Two models, by estimation method:**
  - **DIRECT** `var~mean` — the variance of a node estimated from its **own** counts (Poisson/NB
    sampling + dispersion). Feeds count precision at count-observable nodes and `σ²_global`.
  - **IMPUTATION** `var~mean` — the variance of estimating a node by **boundary→region imputation** (the
    boundary↔region↔boundary prediction error). Feeds count precision at imputed/AMBIG nodes. Properly
    humbler than DIRECT (it measures real prediction error, not the anchor's own confidence — the §15
    over-trust fix).
- **All-gDNA init + the smooth/spiky insight.** At init both models train on **all** nodes using total
  unspliced mass (the all-gDNA assumption). The IMPUTATION error is then **high at enriched exons**
  (depleted-intron boundaries badly under-predict an enriched exon's total) and low on smooth regions — so
  the model correctly reports high variance where RNA is present, and the error **falls as the
  deconvolution removes the spiky RNA**, leaving a smooth, boundary-predictable gDNA residual. *gDNA is
  spatially smooth (boundary-predictable); mature RNA is spiky.* Minimizing boundary→region error therefore
  *separates smooth-gDNA from spiky-RNA* — a second lever beyond the strand. **Caveat:** it does **not**
  separate gDNA from **nascent** (both smooth) — that is the strand's job and is *invisible* to the
  prediction error; in unstranded data the smooth gDNA/nascent split is the honest identifiability limit
  (the global baseline governs).
- **Locality / span (open):** a single LOESS span over-smooths small targeted panels; the SCAM penalty
  (one GCV-λ) adapts better, but stress it on a targeted-panel sim. The span/locality calibration is the one
  remaining fitter research task (`new_var_fit_plan.md`).

## 5. The iteration — convergence & per-node confidence

- **Loop:** all-gDNA init → fit var~mean → solve → re-fit on new estimates → … → converge.
- **Convergence signal:** the IMPUTATION boundary→region **prediction error** (it falls as RNA is correctly
  removed). It tracks the *smooth/spiky (gDNA-vs-mature) separation*; it is a **monitor**, not the
  objective (the objective is the full strand-anchored posterior — minimizing prediction error alone would
  call smooth nascent "gDNA").
- **Stopping rule:** the per-pass **mass delta** stabilizes (and/or the prediction error plateaus); cap at a
  few passes. Guard against **oscillation** with the global anchor + monotone var~mean + optional damping
  (blend with the prior pass); track the delta.
- **Per-node confidence** (combine): (1) the predicted variance at the node's mean; (2) the standardized
  residual `|observed−predicted|/predicted` (does the node conform); (3) the local-data-density
  `confidence(μ)` (extrapolation guard); (4) **cross-iteration stability** (a node that stops moving has
  converged — the most reliable, iteration-only signal).

## 6. Diagnostics feature (decided: dataframes, not plots)

Optional calibration outputs behind `CalibrationConfig.emit_diagnostics` / CLI `--calibration-diagnostics`,
written as dataframes (feather+TSV) to `<out>/calibration_diagnostics/`; a companion
`scripts/plot_calibration_diagnostics.py` renders them. **No matplotlib in `rigel quant`.** Dataframes:
`varmean_{gdna,rna}` (the fitter's points+curve+confidence), `per_region` (f_g, count_frac, g_strand,
μ_global, strand_class, mass, per-node confidence, per-pass deltas), `global_scalars` (ρ_global, σ²_global,
κ, overdispersions, convergence trace). Essential for real-data QC. (`per_node_deconv_hierarchy_design.md`
§18.)

## 7. Current state (2026-06-15)

**Committed (`main`):** eff-len incidence fix (`20d43f2`), evidence-weighted eff-len shrinkage (`54e7110`),
spliced-mass withhold (`17d0209`), dead-code cleanup (`78a1e62`), CLAUDE.md sync (`7a4ad2d`).

**Uncommitted (working tree — the redesign):**
- `variance_model.py` — SCAM `MonotoneVarMean` (fit/predict/confidence/to_dataframe) + DIRECT/IMPUTATION
  builders + `varmean_points`. `test_variance_model.py` (8 tests). ✅ validated.
- `simplex_sweep.py` — node-class prior (Jeffreys single-strand + global AMBIG), the global gDNA prior, a
  per-node `count_precision` param. β=10 is the **current placeholder** for count precision.
- `calibrate.py` — `use_propagation` branch is now the **unified** sweep (hybrid retired): global prior +
  β-placeholder count + propagation; `need_count_variance` on for the sweep; passes `ρ_global` + eff-len.
- `density_model.py` — exposes `global_density` (ρ_global) + `count_rel_var` (v_rel).
- `config.py` — `sweep_n_grid`.
- `scripts/debug/scam_var_mean.py` — the fitter iteration surface (uses production `variance_model`).
- Docs: this plan, `per_node_deconv_hierarchy_design.md`, `iterative_bootstrap_design.md`.

**Current sweep result** (uncommitted, opt-in `use_propagation=True`): **wins the complex battery** (2D/1D
0.95) and is ≈parity-slightly-worse on the simple suite (flagship +9k EM-bound; ss0.5 +29k; zero-DNA ~2k) —
those three residuals are the targets of the var~mean + iteration below. **Not yet default; not committed.**
1106 full-suite + 239 calibration tests green.

## 8. Implementation plan (ordered)

1. **Variance-model fitter — finalize** (mostly done): lock the SCAM `MonotoneVarMean`; settle the locality
   /span (targeted-panel sim stress). [`variance_model.py`; `scam_var_mean.py`]
2. **DIRECT / IMPUTATION builders — by the right cases**: DIRECT = own-count variance; IMPUTATION =
   boundary→region prediction error (the all-three-observable case is the trainable one). Emit per-node
   prediction error for the convergence signal + per-node confidence.
3. **Wire count precision** = `confidence(μ)·1/var(μ)` (DIRECT at observable, IMPUTATION at imputed),
   replacing the β placeholder. **Wire `σ²_global`** = DIRECT var~mean at ρ_global, replacing flat `τ=1`.
4. **The iterative control flow** in `calibrate`: all-gDNA init → fit → solve → re-fit → converge
   (mass-delta stop; oscillation guard). Per-node confidence (4 signals).
5. **Diagnostics feature**: config/CLI flag → dataframes; companion plotter.
6. **Validate** (§9); then **flip `use_propagation` default + delete the fusion path** (separate commit);
   **consolidate** the scattered var~mean code into `variance_model.py` (rna_variance / density_variance
   become thin callers).

## 9. Validation gates (every change)

- Scoreboard fusion vs new (`/tmp/compare_pipelines.py`): flagship ss0.99 / ss0.50 / zero-DNA
  (`nrna_none` & `nrna_rnd`) / capture-off — **≥ fusion, no regression**; zero-DNA phantom ≈ fusion's ~30.
- Complex battery (`scripts/debug/complex_loci_benchmark.py`): TOTAL 2D/1D ≤ 0.95; anchor-starved families
  improve (global fallback, not vertex-push).
- Convergence: the per-pass mass delta → 0 in a few passes; prediction error monotone-ish down; no
  oscillation.
- Full suite (`pytest tests/`) green; goldens unaffected until the default flip; conservation
  `gdna+rna=total`.

## 10. Open questions

- Locality-adaptive span / monotone-spline knot count (targeted-panel stress sim).
- Stopping threshold + max passes; damping yes/no (observe oscillation first).
- Per-node confidence weighting of the 4 signals.
- Unstranded gDNA-vs-nascent: accept the global-baseline fallback (likely) vs any further signal.
- Stage-3 cancer-CNV-aware global baseline (the global prior currently assumes a smooth baseline;
  amplifications/deletions are local — a later refinement, the global prior's variance partly absorbs it).
