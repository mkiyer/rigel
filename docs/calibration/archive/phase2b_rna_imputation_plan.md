# Phase 2b — the RNA imputation model (execution-ready plan)

> **⚠️ SUPERSEDED by `CALIBRATION_PLAN_v5.md` (2026-06-16, same day).** v5 folds this RNA imputation into the
> **unified full-pie imputation** (all-pairs, per-strand, bidirectional, node-agnostic) rather than a separate
> RNA-prior bolt-on, and the §3 coherence guard becomes structural. The RNA `var~mean` reliability fit + the
> spliced-floor + FL-density details below are carried into v5 (§2-3); the separate-bolt-on wiring is replaced.

**Status:** SUPERSEDED by v5 (2026-06-16). Scope of `CALIBRATION_PLAN_v4.md` Phase 2 ("2b"). Grounded against
`main@8d9df6bd` (post-2a) + an adversarial design review. Companion: `CALIBRATION_PLAN_v4.md` (Amendment A1,
the var~mean machinery), the SHIPPED Phase-2a gDNA imputation (`variance_model.fit_imputation_varmean_current`).

## 0. What this is (and is not)

The **RNA imputation model** is the RNA analog of the shipped gDNA imputation. For each boundary↔region pair,
per strand `s`, fit a `var~mean` reliability for predicting the region's **same-strand unspliced RNA** from the
boundary's **same-strand (unspliced-RNA + spliced) crossings**, and feed it as a per-strand RNA **prior**
`(μ_rna, τ_rna)` on `f₊`/`f₋` in the sweep's per-node `ψ`. It is the symmetric twin of the gDNA prior:

| | predictor (source) | target | reliability | feeds |
|---|---|---|---|---|
| **gDNA** (2a, shipped) | boundary unspliced **gDNA** density | region gDNA density | `var~mean( (region gDNA − impute)² )` | `τ_count` → `f_g` prior |
| **RNA** (2b, this plan) | boundary (unspliced-RNA **+** spliced) density, strand `s` | region unspliced-RNA density, strand `s` | `var~mean( (region RNA − impute)² )` | `τ_rna±` → `f₊/f₋` prior |

Same machinery (the `MonotoneVarMean` SCAM fitter, current-density axis, Jensen df-offset, Bernoulli clamp),
different data. **Including the boundary *unspliced* RNA — not just spliced — is the key refinement:** the
boundary unspliced carries the same nascent+mature mix as the region target, removing the one-directional
nascent under-prediction that a spliced-only predictor would have, so the prior can be **two-sided** rather
than a lower bound.

**NOT in scope:** the region↔region propagation coupling `q_rna` (still `0.25`). That is a *different*
mechanism (transitive continuity vs this one-hop magnitude) and is **sequenced after Phase 5, not subsumed**
(see §6). The Phase-1 A/B (RNA-odds helps +1.1–5.0%) remains its live arbiter.

## 1. The model — precise

Per boundary side and strand `s` (oriented to a genomic strand exactly as the sweep orients spliced-sense,
`simplex_sweep.py:256-259` — only at single-strand-adjacent sides; see §4):

- **Predictor density** `ρ̂_rna(side,s) = (boundary_unspliced_RNA_s + boundary_spliced_s) / rna_fl_mean`, where
  `boundary_unspliced_RNA_s = crossing_unspliced_s − boundary_gDNA_s`. **`boundary_gDNA_s` MUST be the SAME
  cleaned gDNA count that feeds `μ_count`** (`cleaned_gdna_count`, `calibrate.py:185-186`) — see the coherence
  invariant §3. `rna_fl_mean` = the RNA crossing eff-length (`calibrate.py:106`).
- **Target density** `ρ_rna(region r,s) = region_unspliced_count_s · (1 − f_g_region[r]) / region_eff_len_rna[r]`,
  with `f_g_region` = the **previous pass's** `regions.gdna_frac` and `region_eff_len_rna` the RNA region
  eff-length (`calibrate.py:112`) — source and target share the **RNA** density axis (FL-consistency).
- **Fit point** (per eligible region, per strand): `mean = ρ_rna(r,s)` (current density — fit & query on the
  same axis, the 2a no-extrapolation contract), `raw_var = ¼·(ρ̂_L − ρ̂_R)²` (the two-side boundary disagreement,
  `dof=1` → Jensen `+[log(½)−ψ(½)] = +1.2704`). **Pool `+` and `−` strands into one fit** (the reliability is a
  property of the boundary→region geometry, not the strand — maximizes the point count).
- **The prior** in `_local_loglik`: `μ_rna±` = the predicted region RNA *fraction of total unspliced* on each
  strand; `τ_rna± = 1 / σ²_rna±`, `σ²_rna± = imputation_rna.predict(μ) · jac²`, `jac = region_eff_len_rna /
  mass_unspliced` (the RNA density→fraction Jacobian, the RNA analog of `geom2`), **Bernoulli-clamped** to
  `μ_rna±·(1−μ_rna±)` and capped at the node's RNA-count Poisson ceiling — exactly the `τ_count` construction
  (`calibrate.py:294-308`). Added as `ψ −= ½·τ_rna₊·(f₊ − μ_rna₊)² + ½·τ_rna₋·(f₋ − μ_rna₋)²`.

## 2. Hard guard #1 — KEEP the sided spliced lower bound

The existing sided spliced lower bound (`simplex_sweep.py:121-124`) is **retained**, not replaced. Spliced is a
**direct, gDNA-free, strand-resolved observation** of mature RNA (a spliced read cannot be gDNA); `μ_rna` is an
**imputed, gDNA-deconvolution-dependent, smoothed** estimate. A two-sided Gaussian at `μ_rna` with finite
`τ_rna` could pull `f±` **below** the directly-observed spliced density — asserting less RNA than was literally
seen, which is incoherent. Architecture: the lower bound is a **hard likelihood-side floor** (`f± ≥ observed
spliced density`); the RNA prior is a **soft Gaussian on top** — they compose on the lattice exactly as the
gDNA count prior composes with the strand mixture. (No double-count: the floor uses the region's **own
contained** spliced `c.n_spliced_sense`; the prior uses the **neighbour boundary** spliced crossing flux.)

## 3. Hard guard #2 — ONE consistent boundary deconvolution (the load-bearing correctness condition)

On the shared 2-simplex (`f₊+f₋+f_g=1`, normalized once), the gDNA prior pulls `f_g→μ_count` and the RNA prior
pulls `f±→μ_rna`. The lattice normalization guarantees the posterior sums to 1 — there is no arithmetic
violation. **But** if `μ_count` and `μ_rna` come from **two independent boundary deconvolutions**, they can
encode contradictory beliefs about the *same* boundary's gDNA-vs-RNA split (e.g. the gDNA split says `f_g=0.7`
while the RNA split independently implies `f_g=0.4`), and the lattice will silently average two priors pulling
in opposite directions — a double-count of the boundary evidence.

**Invariant (enforce as a hard assertion, not an aspiration):**
`boundary_gDNA_s + boundary_unspliced_RNA_s == crossing_unspliced_s` per strand-pooled side, with
`boundary_gDNA_s` the **same** `cleaned_gdna_count` quantity feeding `μ_count`. This is the natural construction
(`splice_junction.boundary_gdna_fraction:157-159` already builds `ρ_gDNA` and `ρ_RNA` from one split) — make it
a unit assertion so the two priors partition **complementary** (not redundant) information.

## 4. Per-strand orientation (a real wiring constraint)

`substrate.left/right.n_spliced_sense/antisense` are **transcript-sense-relative**, not genomic-strand-relative;
and `left_split/right_split.gdna_frac` is a single per-side fraction (NaN where strand-uninformative). Both are
cleanly defined **only at single-strand-adjacent sides**. So the RNA-imputation training set is restricted to
**single-strand exon regions (`TS_POS`→`+`, `TS_NEG`→`−`) with junction-bearing observable flanks** — orient
via `_side_strand_orientation` (`strand_deconv.py:255-283`), mirroring `simplex_sweep.py:256-259`. AMBIG-side
per-strand spliced is undefined and is **deferred** (consistent with deferring AMBIG/`q_rna` to Phase 5). This
restriction strictly shrinks the point count vs the gDNA fit — see the sparsity fallback (§5) and the binding
constraint (§7).

## 5. Mandatory sparsity fallback (makes 2b a strict, data-gated upgrade)

The same-strand-mature-reference requirement cuts the point count below 2a's gDNA fit; `MonotoneVarMean.fit`'s
power-law fallback fires below `max(k,8)=18` points, and the RNA pair count will often be lower. **When the
eligible same-strand RNA pairs are below the fit threshold, do NOT emit a two-sided RNA prior — fall back to
the existing sided spliced lower bound alone** (which §2 keeps regardless). Phase-2b's two-sided prior then
*activates only when the data supports it* — a strict upgrade, never a regression. (Also: the RNA prior is
naturally **inert at pass 0** — all-gDNA ⇒ boundary RNA = 0 ⇒ predictor collapses to spliced-only ⇒ the
two-sided term vanishes; it activates pass 1+. `sweep_max_passes=6` gives ≥2 post-init passes for it to act.)

## 6. File-level implementation (execution-ready)

1. **`variance_model.py` — `fit_imputation_rna_varmean(...)`** (clone `fit_imputation_varmean_current:351-400`):
   reuse `count_observable_masks` + `same_ref_left_right` + the `la/rb` side masks; eligibility = single-strand
   region with a junction-bearing observable flank carrying `>0` mature mass on `s` and both adjacent sides
   usable; `mean = ρ_rna(r,s)` (current), `raw_var = ¼(ρ̂_L − ρ̂_R)²`, `dof=1`, pooled over strands;
   `MonotoneVarMean.fit(mean, raw_var, dof=ones)`. Add to `__all__`. (Reuse the `ρ_RNA` numerator form from
   `splice_junction.boundary_gdna_fraction:157-159` — compute the density, don't call it.)
2. **`simplex_sweep.py` — `_local_loglik` + `deconv_regions_sweep`:** add params `rna_frac_pos/neg`,
   `rna_precision_pos/neg` (default `None` ⇒ no-op, exact back-compat). Insert **after** the count term
   (`:125-129`), **before** the node-class prior: `psi −= ½·τ_rna₊·(f_pos−μ_rna₊)² (+ f_neg twin)`. **KEEP**
   `:121-124` (the spliced floor). Thread the four kwargs through `deconv_regions_sweep` (`:230`, call site
   `:297`); update the `ψ` docstring.
3. **`calibrate.py` — bracketed loop wiring** (mirror the gDNA var~mean, `:275-309`): freeze the boundary RNA
   estimates once before the loop (the `cleaned_left/right` + per-side spliced + `left_split/right_split.gdna_frac`
   from `:171-186`); per pass, compute the previous-pass per-strand region RNA density from `regions.gdna_frac`,
   fit `fit_imputation_rna_varmean`, predict `μ_rna±` at each node's current RNA density, derive `τ_rna±` via the
   `(region_eff_len_rna/mass_u)²` Jacobian + Bernoulli clamp + `min(·, mass_u)` ceiling, and pass the four kwargs
   into `deconv_regions_sweep`. `region_eff_len_rna` (`:112`), `rna_fl_mean` (`:106`) already in scope.

## 7. Validation gates

- **Suite** (`python -m pytest tests/`) green; goldens regenerated only after the accuracy gate (as in 2a).
- **Coherence assertion (§3):** a unit test `boundary_gDNA + boundary_unspliced_RNA == crossing_unspliced`
  per strand-pooled side; and `μ_rna₊ + μ_rna₋ + μ_count ≈ 1` on a clean single-strand region (the
  Jacobian-targets-total-unspliced check).
- **Jensen unit test** for the new `k=2` RNA points (Δ=+1.2704); **fallback unit test** (sparse pairs ⇒
  no two-sided prior, spliced floor only).
- **Monotone-convergence test (the new cross-pass risk, §Open-1):** on a capture-on + nascent toy, assert the
  per-pass `μ_rna` / mass-delta converges monotonically (no gDNA↔RNA oscillation). If it oscillates, under-relax
  the RNA estimate fed to the next pass's fit (the bracketing provides the hook).
- **No-regression:** complex battery TOTAL ≤ 4502 (the 2a baseline); net-flow non-regressing on capture-on
  ss0.99/ss0.5 + zero-DNA (before/after on identical BAMs, as in 2a) — especially the **mature-RNA false-positive
  rate** must not rise (the RNA prior must not manufacture RNA).
- **Activation-rate measurement (gates ROI — do FIRST, §Open-2):** instrument the fraction of exon nodes where
  the two-sided RNA prior activates (eligible + above the fit threshold) on representative geometry. **If it
  activates on <~10% of exon nodes, Phase 2b is marginal and should be deprioritized below the gDNA-side
  Phases 5/6.**

## 8. Residual uncertainties / open issues (ranked)

1. **gDNA↔RNA cross-pass oscillation (MEDIUM).** A *second* iterating prior on the *same* simplex is a coupling
   2a never had: gDNA prior raises `f_g` → less `f±` → next pass's RNA predictor sees less boundary RNA → pulls
   `f±` lower → frees more `f_g`. Bracketing kills the within-pass loop; frozen κ + overdispersions make the
   direction stable so it *likely* contracts — **but this is not proven.** Resolve via the monotone-convergence
   test (§7) before trusting; damp by under-relaxing if needed.
2. **Data sparsity (HIGH — the binding constraint on whether 2b does anything).** The same-strand-mature-reference
   requirement (§4) strictly cuts below 2a's gDNA pair count; Phase-1 flagged same-strand-exon adjacencies as
   sparse. The fallback (§5) makes it non-regressing, but **measure the activation rate first** — if marginal,
   deprioritize. This is the single most important gate for sequencing.
3. **Coherence invariant (HIGH if violated, LOW if enforced).** The §3 one-deconvolution condition is load-bearing
   — enforce as a hard assertion, do not assume.
4. **Residual predictor bias (LOW, accepted).** `boundary_unspliced_RNA = crossing_unspliced − boundary_gDNA`
   inherits the gDNA-split error, so the RNA var~mean absorbs gDNA-split error into the RNA reliability — a
   coupling 2a's raw-gDNA predictor did not have. Accepted (same bracketing 2a accepts); the two-sided prior
   activates only pass 1+ when a gDNA split exists (pass 0 degrades gracefully to the spliced floor).
5. **TSS/TES nascent-smoothness break (LOW, self-mitigating).** The boundary-nascent ≈ region-nascent assumption
   breaks at promoters/terminators/IR isoforms — but the `¼(d_L−d_R)²` disagreement is large there → `τ_rna`
   small → the prior yields to the strand likelihood + the spliced floor. The var~mean is designed to learn this;
   the only requirement is the 2-side disagreement (enforced by the builder).
6. **Cost (LOW).** A second var~mean fit per pass ~doubles the refit (already the wall-time bottleneck at
   90–246 ms vs 25–33 ms sweep). No new mechanism — absorbed by the eventual C++ kernel / fit caching.
7. **`q_rna` framing (LOW, doc-only).** Reframe v4 Amendment A1's "subsumes" language: local one-hop RNA
   imputation (magnitude) and transitive region↔region propagation (continuity) are **complementary**, not
   redundant; `q_rna` is **sequenced after Phase 5, not subsumed**. Phase-1 A/B stays the live arbiter.

## 9. Recommended sequencing

Per §7/§8.2, **measure the activation rate before the full build** — it decides whether 2b is worth doing now
or should yield to Phases 5/6 (the gDNA-side unified solver + IMPUTATION re-cut, where the bigger error lives).
If activation is healthy (≥~10% of exon nodes), build per §6 behind the two hard guards (§2, §3) and the
fallback (§5), and validate per §7 (the monotone-convergence test is the new must-pass). If marginal, land the
fallback-only path (the spliced floor stays) and defer the two-sided prior.
