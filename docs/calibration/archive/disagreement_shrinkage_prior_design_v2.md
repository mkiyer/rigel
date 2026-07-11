# Message-precision fix: Poisson disagreement-variance model — design & implementation plan

**Status:** design, ready to implement (2026-07-07). Clean rewrite; supersedes the message-precision half of
`dispersion_aware_message_precision.md`. Grounded on the AMBIG-dense benchmark
(`~/Downloads/rigel_runs/ambig_dense_10mb`, `gdna300 / ss0.99 / capture-on`). Related:
`mature_rna_channel_design.md` (the mature motivation that surfaced this bug), `ambig_dense_benchmark.md`.

> **TL;DR.** Message precision is currently built from the *source node's own belief certainty* (`vbg`), which
> collapses to ~0 as the sweep sharpens beliefs — so precision runs to ~10⁹ and a single message can override
> the intrinsic strand signal (region 503: a mature message at precision 1027 flips `f_g` 0.63→0.83, phantom
> gDNA, the #1 benchmark error). The fix replaces that with a **derived** variance: the irreducible
> disagreement between adjacent nodes (estimated from the data) plus the source's honest Poisson sampling.
> No source can be more certain than adjacent nodes actually agree.

---

## 1. The defect

`bp_solver._scan` sets each message's precision from
```
base_var  = vbg[src] + pois            # source-belief log-variance + 1/count sampling
s2_edge   = max(resid² − (base_var + var_loc[dst]), 0)
pr        = 1 / max(base_var + s2_edge, _EPS)     # _EPS = 1e-9  ⇒  precision up to ~1e9
```
`vbg[src]` is the source's *self*-uncertainty; it shrinks monotonically toward 0 as the source absorbs
messages, and there is **no term for the irreducible disagreement between two adjacent nodes** — how much a
source's density fails to predict its neighbour's even when both are perfectly measured. So `base_var → 0`
and precision explodes. Measured on the benchmark: message precisions up to 2250, when the empirical
adjacent-node disagreement says nothing should exceed ~0.3. Consequence (region 503, gene G0102, a
single-strand `+` exon, `true f_g = 0.625`): an 87%-mature RNA⁺ message at precision **1027** overrides the
strand-only solve (0.626 ≈ truth) and drags `f_g` to **0.833** — a phantom-gDNA over-call, the largest single
error on the benchmark.

## 2. The model — a Poisson disagreement-variance (fully derived, density space)

**Everything is log-density.** Messages carry **densities** `ρ = count / eff_len`, the only cross-node
comparable quantity. Fractions are node-internal (a node's `f_g/f_pos/f_neg` simplex) and never travel.

**Per node.** Count `nᵢ ~ Poisson(ρᵢ·Lᵢ)` with eff-length `Lᵢ`; observed log-density `xᵢ = log(nᵢ/Lᵢ)`. By
the delta method the sampling variance of the log-density is `Var(xᵢ) ≈ 1/nᵢ`.

**Per adjacent pair.** The imputation edge is always **one boundary and one region** (never boundary-boundary
or region-region — that is the message topology). The observed disagreement decomposes into the *true*
imputation gap plus the two nodes' independent Poisson noise:
```
d_ij      = xᵢ − xⱼ = (μᵢ − μⱼ) + (εᵢ − εⱼ)
Var(d_ij) = σ²_imp  +  1/nᵢ + 1/nⱼ
```
`σ²_imp` is the irreducible adjacent-node density disagreement (biology + capture + mappability) — the
quantity we want.

**Estimator (all pairs, no threshold).** Subtract each pair's *known* sampling and average, inverse-variance
weighted so noisy low-count pairs contribute less (threshold-free, no data discarded):
```
σ²_imp = clamp₊( [ w₀·σ²_chance + Σ wᵢⱼ·((d_ij − d̄)² − 1/nᵢ − 1/nⱼ) ] / [ w₀ + Σ wᵢⱼ ] )
         wᵢⱼ = 1/(1/nᵢ + 1/nⱼ) = (nᵢ·nⱼ)/(nᵢ + nⱼ)          # harmonic form: 0 at zero count, no blowup
```
* `wᵢⱼ` in harmonic form is `0` when either count is 0 (no division by a small number).
* `w₀·σ²_chance` is a **shrinkage prior toward the chance disagreement** (§5): with abundant data
  (`Σwᵢⱼ ≫ w₀`) it is inert; as a component becomes so sparse that `Σwᵢⱼ → 0` it falls back *smoothly* to
  `σ²_chance` (the "no-information" variance → weak messages), and the denominator `≥ w₀ > 0` so the estimator
  never divides by ~0. `w₀ = 1` (one virtual chance-pair). Relevant only in v2 (sparse per-component); inert in
  v1 (total density has thousands of pairs).
* **`d_ij` and the resid used below MUST be the actual message residual `log ρ_src − log ρ_dst` as `_scan`
  forms it** (the dst-frame conversion baked in), **asserted** equal to a `_scan` value on a sample edge — never
  a `mass/eff` proxy. `ρ_boundary` (crossing eff `E[min(ℓ,L)]`) and `ρ_region` (contained eff `E[max(0,L−ℓ)]`)
  must be on the same scale, or the disagreement carries a spurious frame offset. `d̄` = the residual median
  (removes any residual systematic offset). **Validate the centered residual is ~zero-mean and not
  length-dependent**; a surviving length-dependent offset is a genuine eff-length frame bug (fallback:
  stratified/directional median). A crude `mass/eff` proxy shows a length-dependent offset (median +2.6 short →
  +1.2 long) precisely because it is *not* frame-consistent — hence the assertion. `clamp₊` = max(·, 0).

**Message variance (the payoff).** A message from source `s` targets the *destination's true density*, so it
carries the source's sampling plus the imputation gap — **not** the destination's own sampling (that lives in
the destination's local belief):
```
σ²_msg = σ²_imp + 1/n_src
pr     = 1 / σ²_msg  =  n_src / (n_src·σ²_imp + 1)     # algebraic form: denom ≥ 1, pr = 0 exactly at n_src = 0
```
The right-hand algebraic form is how `_scan` computes it — it never forms `1/n_src`, so there is no
division-by-zero and no `if`-guard/clamp (a zero-count source yields `pr = 0` smoothly, i.e. sends no message).

**What this replaces / removes** (net simpler than today, zero tuned constants):
* `vbg[src]` — dropped (the collapsing source-belief term; the bug).
* `pois` — *kept, but now derived* as `1/n_src` from the same Poisson model, not a bolt-on.
* the high-count subset — not needed; all pairs used, inverse-variance weighted. (It was also **biased**:
  high-count pairs are disproportionately exon↔boundary composition mismatches, over-estimating `σ²_imp`.)
* the pseudo-count `ν` / `resid²` shrinkage — not needed for v1: using `σ²_msg` directly is the
  fully-shrunk-to-prior limit, justified because a single `resid²` is one-sample noise. (Optional
  anomaly-detection refinement, §8.)
* `var_loc` and the `max(·,0)` clamp — dropped (fraction-space; subsumed).

**Behaviour.** High-count source → `σ²_msg ≈ σ²_imp` → precision `≈ 1/σ²_imp`, **bounded** (region 503:
`σ²_imp ≈ 3.2` → precision ≈ 0.31 vs the strand's ~50 → the strand wins; fixed). Low-count source →
`1/n_src` large → low precision (directional, as it must be). The floor is structural: `σ²_imp > 0`
empirically, so no source certainty can drive `σ²_msg` to 0.

## 3. Task scope: TWO components now (gDNA, RNA); three later

Two independent tasks, tackled one at a time:
* **Task 1 (this plan) — per-component disagreement variance.** The current message architecture has exactly
  **two components: gDNA and (total) RNA.** Each RNA message already lumps nascent + mature into one RNA density
  (per strand), so there is no mature/nascent split to give separate variances *yet*. v2 = a **two-component**
  variance model: one `σ²_imp` for gDNA, one for RNA.
* **Task 2 (deferred) — separate mature-RNA channel.** Splitting each RNA message into a mature (spliced,
  junction) sub-message and a nascent (unspliced, contiguous) sub-message (`mature_rna_channel_design.md`).
  When that lands, we **extend the same estimator/apply machinery from two components to three** (gDNA, nascent,
  mature) — §4.3. Not part of this plan.

Strand is orthogonal to component: RNA⁺ and RNA⁻ are the same component "RNA" (same imputation reliability), so
the RNA variance is estimated **pooled over both strands** and applied to both RNA messages.

### 3.1 Staged estimation: total density → two components

`σ²_imp` from **total density** is fixed (the data does not change). Iteration is meaningful only because we
move from total density to **per-component** densities, which depend on the evolving belief. The existing
two-pass structure (calibrate.py: PASS-1 → fit KDE → PASS-2) gains a disagreement-refit at each fit point:

0. **Pass-0 (pre-solve).** Belief = the existing signature-binary `init_beliefs`. Estimate one **total-density**
   `σ²_imp` over all adjacent boundary↔region pairs (§2 — this is exactly the shipped v1). Floor every message
   with `σ²_msg = σ²_imp + 1/n_src`.
1. **Pass-1.** `node_sweep` with the total-density `σ²_msg`. Confident locals win; no message is tyrannical.
2. **Refit-1 (post-pass-1), beside the existing KDE fit.** From the pass-1 belief, estimate the **two**
   `σ²_imp` — one for **gDNA**, one for **RNA** — on adjacent pairs, using each component's solved density
   `ρ_c = f_c·mass/eff_c (+ spliced/eff_spl for RNA)` and its component count `n_c` (so the Poisson term is `1/n_c`):
   * **gDNA:** `ρ_g = f_g·mass/eff_gdna`, `n_g = f_g·n`; **all** adjacent edges (gDNA is genomically universal).
   * **RNA:** the per-strand RNA density (`f_s·mass/eff_rna + spliced_s/eff_spl`), `n_rna,s = f_s·n + spliced_s`;
     only **strand-live** edges (`free_s` on both endpoints — exclude intergenic and the OFF strand of a
     single-strand node), **pooled over both strands**.
   * **Exclude AMBIG** node pairs (the belief is least trustworthy there this early).
3. **Pass-2.** `node_sweep` with the two-component `σ²_msg` (gDNA message → `σ²_imp,gDNA`; RNA⁺/⁻ →
   `σ²_imp,RNA`) + the KDE prior (fills the AMBIG τ null-space).
4. **Refit-2 (post-pass-2).** Recompute the two `σ²_imp`, **now including AMBIG** pairs (they are solved).
   Adding AMBIG changes the population, so `σ²_imp` may *legitimately loosen* here vs Refit-1 (which measured
   only the unambiguous genome and may be artificially tight) — this is expected, not a failure. **Monitor
   `σ²_imp` across Refit-1→Refit-2** empirically.
5. **Pass-3 if not converged.** Stop when `Δbelief < tol` (typ. ~2–3 passes). **Monotone sharpening is a
   *within-population* convergence aid** (once the pair population is stable, Refit-2 onward, `σ²_imp` should
   only tighten) — it is NOT a cross-population constraint and must not clamp the Refit-1→Refit-2 transition
   (where AMBIG is added). Apply it only after the population stabilizes; verify empirically.

## 4. Applying the two variances

The three messages map to the two components directly — **no blend needed** (the blend is Task 2):
* **gDNA message** → `σ²_msg = σ²_imp,gDNA + 1/n_g_src`.
* **RNA⁺ / RNA⁻ messages** → `σ²_msg = σ²_imp,RNA + 1/n_rna,s_src`.
`n_*_src` is already the component's source count in `_scan` today (gDNA `fbg·sm`; RNA `n_nasc + n_mat`), so the
algebraic precision `pr = n_src/(n_src·σ²_imp,c + 1)` needs only the right `σ²_imp,c` per block.

**Why this fixes v1's EX+IN− regression.** v1's single total-density `σ²_imp` (≈3.86) is inflated by the
gDNA-vs-RNA *composition* mismatch at exon/intron boundaries, so it over-weakens **every** message uniformly —
including the gDNA messages that were legitimately correcting a local RNA over-attribution. Splitting into two
gives gDNA its deservedly **small** `σ²_imp,gDNA` (gDNA is genomically smooth → strong gDNA messages) while RNA
keeps a larger `σ²_imp,RNA`. The gDNA correction returns; the RNA over-confidence stays capped.

### 4.3 Task-2 extension (deferred): two components → three

When the mature channel is split out (Task 2), the RNA component becomes two — **nascent** and **mature** — with
separate variances (mature is junction-spliced and capture-biased → larger; nascent is intron pre-mRNA →
smaller). Since nascent and mature both feed the same `f_pos`/`f_neg`, the (then-two) RNA sub-messages combine
via a **density-share-weighted** variance `w_mat·σ²_imp,mature + (1−w_mat)·σ²_imp,nascent` (w_mat = the source's
mature fraction of its strand RNA). The estimator simply gains a third live-edge mask (mature = exon↔junction
edges carrying spliced mass). *Deferred — do not build in this plan.*

## 5. Non-adjacent (random) pairs — diagnostic + sparse-component fallback prior

Random (non-adjacent) boundary-region pairs give `σ²_chance` (≈ 13 on the benchmark vs adjacent ≈ 6). It is
**never a per-edge precision floor** (we allow precision → 0 for uninformative sources; the §2 algebraic form
delivers that smoothly). Its two roles:
* **Diagnostic:** the adjacency correlation `ρ_adj = 1 − σ²_imp/σ²_chance` (≈ 0.54) quantifies how informative
  adjacency is; if a component's `σ²_imp ≈ σ²_chance`, adjacency carries no information for it → flag it.
* **Sparse-component fallback prior (v2):** `σ²_chance` is the shrinkage prior `w₀·σ²_chance` in the §2
  estimator. It is inert when a component has abundant adjacent data and takes over *smoothly* when the
  component is too sparse to estimate (`Σwᵢⱼ → 0`), so the estimator degrades gracefully to "no information"
  rather than dividing by ~0. This is the *variance-estimate* prior, not a per-edge clamp.

We do **not** floor precision: a wildly-disagreeing message is allowed to go silent (weak messages are safe;
revisit only if messages later prove too weak).

## 6. Milestones

**v1 — total-density global `σ²_imp` — SHIPPED (behind flag, branch `calib-disagreement-shrinkage-v1`).**
One global `σ²_imp`, `σ²_msg = σ²_imp + 1/n_src` used directly; a net simplification of `_scan`. Validated on
`ambig_dense_10mb`: region 503 `f_g` 0.833→0.606, precision max 2251→0.5; gene-level mature `spearman_abund`
0.872→0.892 / `mard` 0.398→0.341 on gdna300/unstranded/capture-on, no regression on capture-off / no-gDNA /
stranded; 194 calib tests pass. `CalibrationConfig.sweep_disagreement_shrinkage` (default off).

**v2 — TWO-component `σ²_imp` {gDNA, RNA} + staged refit (§3, §4).** The accuracy lever: gives gDNA its
deservedly-strong message back (fixing v1's uniform over-weakening and the EX+IN− regression) while keeping RNA
capped. Adds the two-component estimator, the AMBIG exclusion / strand-live filtering, and the refit-2 / pass-3
loop. **This plan.**

**Task 2 (separate mature-RNA channel) — extends v2 to THREE components** {gDNA, nascent, mature} + the RNA
density-share blend (§4.3). Deferred; independent of this plan.

## 7. Implementation plan (the v2 delta over the shipped v1)

v1 already ships: `adjacent_disagreement_variance(chain, geometry)` (total-density scalar), the `_scan`
algebraic `pr = n_src/(n_src·sigma2_imp + 1)` branch, and the `sweep_disagreement_shrinkage` flag threaded
through `node_sweep`/`calibrate.py`. v2 changes:

**`bp_solver.py`**
* Generalize the estimator to **per-component** — `adjacent_component_disagreement_variance(chain, geometry,
  belief, statics) -> (σ²_gDNA, σ²_RNA)` — computing, per adjacent boundary↔region edge, the component density
  `ρ_c` and count `n_c` from the **solved belief**: gDNA `(f_g·mass/eff_gdna, f_g·n)` on all edges; RNA the
  per-strand `(f_s·mass/eff_rna + spliced_s/eff_spl, f_s·n + spliced_s)` on `free_s` edges, pooled over ±. Each
  uses the §2 Poisson estimator (subtract `1/n_i+1/n_j`, inverse-variance weight, median-center) — reuse the
  v1 core. Masks: exclude AMBIG (refit-1), include AMBIG (refit-2). Keep the `_scan` frame-consistency assertion.
* `_scan` / `node_sweep` — `disagreement_sigma2` becomes a 2-field struct `(gdna, rna)`; the gDNA block uses
  `.gdna`, the RNA⁺/⁻ blocks use `.rna`. (v1's scalar path stays for pass-0.)

**`calibrate.py`**
* Pass-0/Pass-1 unchanged (scalar total-density σ²_imp — the v1 path).
* **Refit-1** (beside the KDE fit, ~L186–195): call the per-component estimator on the pass-1 belief (excl
  AMBIG) → `(σ²_gDNA, σ²_RNA)`; run **Pass-2** (`_sweep`) with it.
* **Refit-2** (post-pass-2): recompute incl AMBIG; **Pass-3** if `Δbelief > tol` (monotone-sharpen guard,
  within-population only). Log `σ²_gDNA`, `σ²_RNA` at each refit.

**`config.py`** — reuse the existing `sweep_disagreement_shrinkage` flag (v2 is the per-component realization of
the same mode); no new numeric hyperparameter.

## 8bis. Validation (v2, extends §9)
Same A/B protocol. v2-specific checks beyond v1: the **EX+IN− class recovers** (v1: 17.5k→60.5k; v2 target:
back toward baseline); gDNA messages strengthen (σ²_gDNA ≪ σ²_RNA, verified in the refit logs); the modest
present/on total-error win of v1 (−1.7%) grows; capture-OFF / no-gDNA still flat.

## 8. Deferred / out of scope

* **Mature as its own relay channel.** The shrinkage caps mature's *precision*; mature is still folded into the
  running belief `fbp` and relayed downstream as nascent. Whether that residual relay bleed needs separate
  surgery (keep mature out of `fbp`) is decided *after* the v2 A/B. `mature_rna_channel_design.md` holds the
  motivation. The precision bug takes priority.
* **`resid²` anomaly shrinkage.** `σ²_edge = (ν·σ²_msg + resid²)/(ν+1)` re-introduces per-edge seam detection if
  the structural emission gates + `σ²_imp` floor prove insufficient. Adds `ν` (start 1, or empirical-Bayes on
  the adjacent `resid²` spread). Deferred.
* **Low-count sampling *accuracy* (numerical part already solved).** The numerical blowup is gone (algebraic
  `pr` form + harmonic weight, §2). What remains is *accuracy*: `1/n` under-states the true log-Poisson variance
  at `n ≲ 3` (log-density non-Gaussian), so those sources' messages are mildly over-precise (still floored by
  `σ²_imp`) and inverse-variance weighting keeps them out of the *fit*. A continuity correction (`1/(n+c)` or
  exact log-Poisson) is the refinement if low-count messages misbehave in the A/B — one standard constant,
  deferred.
* **Directional `σ²_imp` (centering already handled).** The frame-offset *centering* is handled in §2 (resid
  from the actual `_scan` conversion + median centering + the zero-mean/length assertion). What is deferred is
  maintaining *separate* `σ²_imp` values for boundary→region vs region→boundary if their (post-centering)
  variances prove materially different; v1 pools them into one.

## 9. Validation

In-process A/B (`total_threads=1`) on `ambig_dense_10mb`. Surface per run: `σ²_imp`, `σ²_chance`, the derived
`ρ_adj`, and the message-precision distribution (median / p99 / max) before vs after. Targets / guards:
* **Region 503 / node 1007:** mature message precision 1027 → ~0.3; `f_g` 0.833 → ~0.63.
* **Capture-ON mature-AMBIG** (EX−IN⁺ net +32k, EX+EX⁻, EX⁺IN⁻): error collapses toward the capture-OFF level.
* **Capture-OFF & pure-nascent IN⁺IN⁻:** unchanged — must not regress the already-good cases.
* **No-gDNA condition:** no new phantom gDNA on high-count nodes (the guardrail).
* **Thin / relay nodes** (short exons, near-zero mass): precision appropriately low, not spuriously high.
* Full 24-condition net-flow + mature-transcript accuracy: flat-or-better.

Ship each milestone only on a clean result across all of the above.
