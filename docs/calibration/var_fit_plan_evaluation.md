# Evaluation of `new_var_fit_plan.md` against the production code

**Status:** review / decision input — 2026-06-15. Companion to `new_var_fit_plan.md` (the proposal)
and `CALIBRATION_PLAN_v2.md` (the authoritative plan). Written after the Phase-4 teardown
(`main@2f76875`), against the single production path (the iterative odds-propagation simplex sweep).

**Provenance.** This evaluation came from a code audit: every load-bearing claim below was read from the
source (with `file:line`) and the quantitative claims (the tilt identity, the two-sweep exactness, the
Jensen bias, errors-in-variables attenuation, the divergence of the proposed update rule) were checked by
direct numerical simulation against the actual functions. Where a number is an estimate rather than a
measurement it is marked.

---

## TL;DR

The proposal's strongest instinct — that **node confidence is the crux** — is correct. But the confidence
machine it describes **already exists and is already exact**, the `{f+, f-, fg}` node state **already
exists**, and the three models it proposes to *build* are (respectively) dead code that needs wiring, a
quantity nothing downstream consumes, and a refinement of a calibration stage the repo's own measurements
show is **not the bottleneck**. The dominant end-to-end (flagship) error lives in the **two-component
eff-len construction** (`priors.py` / `capture_eff_length.py`) — a stage this plan does not touch.

The var~mean machinery the proposal would build on does have **real defects** — but *different* ones than
the proposal describes (an uncorrected Jensen/log bias, a fit/predict domain mismatch, and a `geom2`
nullifier), which the proposed models would **inherit rather than fix**.

---

## 1. Verdict table

| Idea | Verdict | Single most important reason |
|---|---|---|
| **Model 1** — locus-level global gDNA var~mean | **Reject (salvage one slice)** | "gDNA uniform within a locus" is false under exactly capture + focal CNV (the cancer use-case); locus pooling starves the spline below its point floor → power-law fallback. Salvage only a pooled **intergenic-floor** point. |
| **Model 2** — gDNA imputation var~mean | **Already exists / modify** | Built as `raw_var ~ mean` (`variance_model.fit_imputation_varmean`), not the proposed `region ~ boundary` regression — and the regression form is *worse* (errors-in-variables attenuates the slope, simulated 1.0→0.61). Keep the object; fix its real defects (§3D). |
| **Model 3** — per-strand RNA imputation var~mean | **Reject / wire-only** | Components exist as **dead code** (`mature_density.py`, `rna_variance.py`, imported nowhere in production); per-strand RNA is **non-identifiable** off junctions and **never consumed** downstream (`priors.py` reads RNA total + scalar κ only). |
| **RNA-not-inverse / `{f+,f-,fg}` node state** | **Already exists** | The production lattice *is* the `(f₊, f₋, f_g)` 2-simplex (`simplex.py:117-130`). "RNA is not the inverse of gDNA" is false *as a constraint*: `f₊+f₋+f_g = 1` ⇒ RNA total ≡ `1−f_g`. The valid reading (RNA carries its own splice evidence) is already the sided spliced lower bound. |
| **Total-imputation-error minimization** | **Reject (current superior)** | Replacing a precision-weighted log-posterior with an unweighted SSE discards heteroscedasticity and double-counts the sum-to-one constraint. The sweep already *is* the principled combine; imputation error already enters as the `τ_count` it sets. |
| **Propagation + node confidence** | **Reject the rewrite (current superior)** | `_sweep_chain` is already exact forward+backward sum-product (matched brute-force joint marginalization to ~9e-16). The proposed iterative confidence-blend, implemented faithfully, **diverges** (posterior variance → 4–6% of truth, mean drifts 0.25, *worse* with more iteration). |

---

## 2. Corrections — where the proposal mis-describes the current code

Each changes what is actually left to build.

- **The node state is already a per-strand 3-simplex.** `_simplex_lattice` returns `(f_pos, f_neg,
  1−f_pos−f_neg)` (`simplex.py:117-130`); the tilt distinguishes strands via
  `p = ½·f_g + κ·f₊ + (1−κ)·f₋` (`simplex.py:141`); spliced bounds, allow-masks, and propagated log-odds
  are all per-strand (`simplex_sweep.py:78-80, 188-189, 38-43`). The only gap is that RNA is **collapsed on
  output** (`rna_mass = (1−f_g)·mass_unspl + mass_spliced`, `simplex_sweep.py:256`).

- **The propagation is already message passing, not a "static grid sum-product."** `_sweep_chain`
  (`simplex_sweep.py:141-155`) runs a forward `alpha` and backward `beta` recursion with
  `belief = psi + alpha + beta`. A per-reference region partition is a linear chain (a tree), so two sweeps
  give the **exact** posterior — verified node-by-node against brute-force joint marginalization to ~9e-16,
  including correctly decoupling at `Q=∞` AMBIG gaps.

- **Per-node confidence already exists.** `_fg_var` (`simplex_sweep.py:120-125`) is the posterior variance
  of `f_g` after propagation, surfaced as `NodeDeconv.gdna_frac_var` (`:257`). The real gap is downstream:
  nothing **consumes** it.

- **Model 2's imputation variance is built — as `raw_var ~ mean`, not `region ~ boundary`.**
  `varmean_points` (`variance_model.py`) forms per-region density measurements and a k-sample disagreement
  variance; `fit_imputation_varmean` is that fit restricted to non-region-observable nodes. The "boundary→
  region prediction error" docstring is loose: it is left-vs-right boundary disagreement regressed on the
  mean. This is the *statistically correct* object (a variance, for a precision), not a conditional mean.

- **Model 3's RNA components already exist — and are dead.** `mature_density.mature_density` (per-strand RNA
  boundary→region imputation) and `rna_variance.rna_spliced_variance` (the splice-pair var~mean) are
  imported only by each other and tests — not wired into production. `deconv_regions_sweep` is called
  **without** `q_rna`, so the RNA coupling variance is the hardcoded scalar `q_rna=0.25`
  (`simplex_sweep.py:161`), not a fit.

- **"RNA is not the inverse of gDNA" is false as a constraint.** On the simplex RNA total ≡ `1−f_g` by
  construction. The defensible reading (RNA carries independent splice evidence) is already encoded (the
  sided spliced lower bound + per-strand odds propagation). The genuine independence — `f₊` vs `f₋` — is
  degenerate off junctions (§3B).

- **(Fixed) stale docstrings** that likely seeded the "it's missing" belief: the per-node simplex was
  mis-described as "mature/nascent" (it is sense/antisense RNA + gDNA); `simplex_sweep.py` claimed gDNA is
  "never propagated" (it couples implicitly via the shared `log(f_c/f_g)` denominator); and two comments
  attributed `τ_global` to "the DIRECT var~mean at ρ_global" when it is the robust MAD of per-node densities
  (`calibrate.py`). Corrected in `main@68ec986`.

---

## 3. The deep issues (ranked)

**(A) ROI — the flagship is EM-bound, and the lever is not in this plan.** Calibration prior/truth ≈ 0.97
at ss0.99; the leak is the EM degrading an accurate prior (observed/prior ≈ 0.93). Root-caused to the
**relative eff-len construction**: the gDNA component uses the raw IPR (no FL haircut, `priors.py`
~251-263) while RNA/nascent use the FL-marginal contracted over the full span including capture-depleted
introns (`capture_eff_length.py`). The decisive sweep `gdna_eff ×0.5` recovers flagship truth (phantom
nascent 88k→15k, gDNA 2.75M→2.94M). **Every model in this plan refines an already-~97%-correct `f_g` and
cannot move the dominant error by construction.** (Corroborates `flagship_em_bound_not_calibration` /
`gdna_leak_root_cause_efflen` / `em_gdna_strand_likelihood_fix.md` §9.)

**(B) Identifiability ceiling on per-strand RNA (Model 3 / total-error).** The unspliced strand likelihood
depends *only* on the tilt `t = f₊−f₋` via the exact identity `p = ½ + (κ−½)(f₊−f₋)`
(`simplex.py:16-18`, verified to machine precision). N unspliced fragments give **one** scalar constraint
on three unknowns; the overdispersion second channel is negligible and misleadingly peaks at `f_g≈0.3`
regardless of truth. Per-strand RNA is separable *only* at junction-flanked exons (from stranded spliced
flux) — i.e. nowhere it is needed (AMBIG single-exon overlaps, introns) — and is never consumed downstream.
So Model 3 sharpens an unidentifiable, unconsumed quantity.

**(C) Double-counting in the proposed propagation (fatal correctness bug).** "Store {value, confidence},
weigh intrinsic vs incoming, propagate the delta, iterate" — implemented faithfully and compared to the
closed-form posterior — re-absorbs each node's own evidence circulating back through neighbours every pass:
after 20 iterations the posterior variance collapses to **4–6% of truth** and the mean drifts **0.25** off,
*worsening* with iteration. This is the worst failure mode (confident-and-wrong) for a prior the EM trusts.
The mathematically correct version of "forward then reverse, weighing confidence" *is* the α/β sum-product
already implemented. On scalar-vs-3-vector confidence: **neither** — `cov(f₊,f₋) ≠ 0` (simulated −0.032),
so three independent scalars are incoherent; the honest object is a 2×2 tangent covariance, which the full
lattice belief already is exactly.

**(D) Errors-in-variables + Jensen bias in the var~mean fit (real defects the models would inherit).**
  - **Jensen/log bias (serious, unfixed):** the fit is `log var = spline(log mean)` by least squares
    (`variance_model.py:107`), but `raw_var` is a k-sample estimate with k=2 (χ²₁) nodes included
    (`:282-283`) and **no df correction** (no digamma offset anywhere). `exp(E[log raw_var])`
    *under*-estimates the true variance — ~3.5× at the dominant k=2 nodes (simulated). Since
    `τ_count = 1/var`, the count prior is systematically **too confident**. Fix: add the df-aware offset
    `c_k = −(ψ((k−1)/2) − log((k−1)/2))` before fitting, or use a gamma GLM with log link.
  - **Fit/predict domain mismatch (serious):** IMPUTATION points are fit on the *boundary* mean (contained
    masked) but the curve is *queried* at the region's own density `gdna_c/eff_len` (`calibrate.py:259-260`).
    Under capture these are different physical quantities at different magnitudes → flat extrapolation. The
    "all means are seen" docstring (`variance_model.py:173-177`) is false for this subset.
  - The literal Model 2 `region ~ boundary` regression is *worse* than the current object: both axes are
    noisy estimates of the same latent → OLS attenuates the slope (simulated 1.0→0.61); and "region is
    truth" is unsound precisely in the IMPUTATION regime (capture-enriched, RNA-contaminated exons).

**(E) The `geom2` nullifier.** `geom2 = (eff_len/mass_u)²` (`calibrate.py`) is typically ≫ 1 (≈14 for a
moderate exon: eff_len=1500, mass_u=400; ≈5625 at mass_u=20). It deflates `τ_count` toward inert and
produces nonsensical fraction-variances >0.25 (a fraction's variance is bounded by 0.25). **Before
refining any imputation variance, confirm the count prior is even influential** — if `geom2` has already
rendered it near-inert at low-coverage nodes, the entire Model 2/3 line is wasted effort until the
density→fraction conversion is reworked.

---

## 4. What to actually build (prioritized)

### FIRST — outside this plan, highest ROI
**#0. Fix the two-component eff-len construction.** Make the gDNA and RNA/nascent component effective
lengths consistent (`em_gdna_strand_likelihood_fix.md` Option B: bind calibration `f_g` as the
gDNA-vs-unspliced-RNA split as an M-step regularizer). Touches `capture_eff_length.py` + `priors.py`.
Expected effect: recovers the bulk of the −7% flagship leak; collapses the phantom-nascent sink. Orthogonal
to every model below. **Defer the var-fit plan until this is honest and re-measured.**

### HIGH ROI / LOW RISK
1. **Instrument the actual `τ_count` / `geom2` / `var_d` distribution** on the gDNA suite (debug script
   only). Settles whether `geom2` has already neutered the count prior — which *gates* whether Models 2/3
   are worth anything. Cheapest decision-making information available.
2. **Clamp the `geom2` fraction-variance at the Beta bound** `count_frac·(1−count_frac)` (`calibrate.py`,
   ~1 line). Addresses the tracked **geom2-inflation** open issue.
3. **Fix the Jensen/log bias** in the var~mean (df-aware offset or gamma GLM; `variance_model.py:80-148`).
   Removes the systematic ~3.5× over-confidence in `τ_count` at k=2 nodes — the largest *statistical*
   defect in the current machinery.
4. **Upgrade the convergence stop to a prediction-error signal** (track the per-pass IMPUTATION `raw_var`
   aggregate instead of `mean|Δf_g|<1e-3`, `calibrate.py`). Addresses the **iteration-tightening** issue.
   This is the one genuinely good idea in the prose, and needs no engine rewrite.
5. **Pool an intergenic-floor var~mean point** (one/two extra DIRECT points from all `sig==0` contained
   fragments — intergenic only, *never* intron-inclusive, which amplifies a circular ρ_global bias). The
   only salvageable slice of Model 1; no locus partition or uniform-locus assumption needed.
6. **Consume `gdna_frac_var` downstream** (surface to `CalibrationResult`; precision-weight the per-locus
   Dirichlet in `assemble_priors`). Realizes the "confidence matters" instinct. **Gate on a measurement**
   that prior mis-confidence is actually costing accuracy.

### SPECULATIVE (only after a measurement justifies)
7. **Data-drive `q_rna` per-edge** from the existing `rna_variance` fit (`deconv_regions_sweep` already
   accepts per-edge `qp/qn`; only the scalar is plumbed). Fix the RNA-strand double-count first; handle the
   all-gDNA-init blowup (pass-0 spliced imputation → huge RNA error → decoupled propagation); A/B on the
   complex-loci battery. Per the flagship evidence it will not move the dominant error.
8. **Clean-triplet imputation floor** (the salvageable bit of Model 2): where a count-observable region is
   flanked by two observable boundaries, emit `err_var = (contained − mean(boundaries))²` as a floor — done
   only where region-as-truth is valid (clean intergenic/intron), not exons. First count how many clean
   triplets exist under capture (likely few).
9. **Per-strand RNA output** (Model 3 proper) — only when a downstream consumer of `f₊/f₋` exists. If built,
   decode the existing belief's `f₊/f₋` marginals + the 2×2 tangent covariance, not three scalars.

---

## 5. What NOT to build

- **The message-passing / confidence engine rewrite.** `_sweep_chain` is already exact (≈9e-16); the
  proposed iterative blend diverges into confident-wrongness. The "CRUX" is already solved.
- **A 3-scalar confidence `{conf+, conf−, conf_gdna}`.** Incoherent: `cov(f₊,f₋) ≠ 0`. The full lattice
  belief already carries the exact anisotropic uncertainty; `_fg_var` is the correct marginal for the only
  output (`f_g`).
- **Total-imputation-error minimization as the objective.** Net-negative: discards precision weighting and
  double-counts the sum-to-one constraint. The sweep already minimizes the right thing.
- **Model 1 as a locus-level estimate.** Uniform-gDNA-in-locus is false under capture/CNV; starves the
  spline; re-includes RNA-dominated exon nodes that the count-observable-only `ρ_global` deliberately
  excludes; and there is no fragment-linked locus partition at calibration time anyway (`build_multi_loci`
  runs post-calibration). Take only the intergenic-floor point (#5).
- **The literal `region_density ~ boundary_density` regression.** Wrong object (a conditional mean, not a
  precision) and worse than the current `raw_var ~ mean` (slope attenuation 1.0→0.61; "region is truth"
  unsound under capture). Keep the monotone-spline `raw_var ~ mean`; fix its Jensen bias + domain mismatch.
- **An intron-inclusive gDNA floor.** Balanced-nascent introns inflate `ρ_global`; pooling more intronic
  mass amplifies a circular bias, worst at low strandedness. Intergenic-only.

---

## 6. Note on the tracked open issues

This plan touches them only partially. **geom2 inflation** → #2 (clamp). **Iteration tightening** → #4
(prediction-error stop). **Count-bias-at-AMBIG** and the **weak-SS (ss=0.65) zero-gDNA phantom (~7.6%)**
live in the per-node count/strand **combine** — the raw contained intron count is never strand-cleaned;
only the boundary crossings are (`calibrate.py` ~170-171) — **not** in the var~mean fit. If those are the
real targets, investigate the combine, not the variance model.

**Bottom line:** the propagation crux is already solved correctly; the per-strand lattice already exists;
the var~mean has real but *different* defects than the proposal describes; and the flagship lever is in the
eff-len construction, a stage this plan does not touch. Do #0 first, then the high-ROI/low-risk set (#1–6),
measure, and only then consider #7–9.
