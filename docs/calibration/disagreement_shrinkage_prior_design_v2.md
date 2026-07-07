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

## 3. Staged estimation: total density → per component

`σ²_imp` from **total density** is fixed (the data does not change). Iteration is meaningful only because we
move from total density to **per-component** densities, which depend on the evolving belief. The existing
two-pass structure (calibrate.py: PASS-1 → fit KDE → PASS-2) gains a disagreement-refit at each fit point:

0. **Pass-0 (pre-solve).** Belief = the existing signature-binary `init_beliefs`. Estimate one **total-density**
   `σ²_imp` over all adjacent boundary↔region pairs (§2). Floor every message with `σ²_msg = σ²_imp + 1/n_src`.
1. **Pass-1.** `node_sweep` with the total-density `σ²_msg`. Confident locals win; no message is tyrannical.
2. **Refit-1 (post-pass-1), beside the existing KDE fit.** From the pass-1 belief, estimate **per-component**
   `σ²_imp` — one each for **gDNA, nascent RNA, mature RNA** — on adjacent pairs, using each component's solved
   density and its component count (`n_c = f_c · n`, so the Poisson term becomes `1/n_c`):
   * **Exclude AMBIG** node pairs (the belief is least trustworthy there this early).
   * **Live-component only:** gDNA on all edges; nascent only where that strand's nascent is active (exclude
     intergenic and the OFF strand of single-strand nodes); mature on exon↔boundary edges carrying spliced mass.
3. **Pass-2.** `node_sweep` with the per-component `σ²_msg` + the KDE prior (fills the AMBIG τ null-space).
4. **Refit-2 (post-pass-2).** Recompute per-component `σ²_imp`, **now including AMBIG** pairs (they are solved).
   Adding AMBIG changes the population, so `σ²_imp` may *legitimately loosen* here vs Refit-1 (which measured
   only the unambiguous genome and may be artificially tight) — this is expected, not a failure. **Monitor
   `σ²_imp` across Refit-1→Refit-2** empirically.
5. **Pass-3 if not converged.** Stop when `Δbelief < tol` (typ. ~2–3 passes). **Monotone sharpening is a
   *within-population* convergence aid** (once the pair population is stable, Refit-2 onward, `σ²_imp` should
   only tighten) — it is NOT a cross-population constraint and must not clamp the Refit-1→Refit-2 transition
   (where AMBIG is added). Apply it only after the population stabilizes; verify empirically.

## 4. The RNA message (nascent and mature share one strand fraction)

Nascent and mature both feed the same `f_pos` (or `f_neg`), so they are not separate simplex axes. The RNA⁺
message keeps one mode (the total RNA⁺ source density) and a **density-share-weighted** imputation variance:
```
w_mat        = ρ_mat_src / (ρ_mat_src + ρ_nasc_src)                       # per edge
σ²_imp_RNA⁺  = w_mat·σ²_imp_mature + (1 − w_mat)·σ²_imp_nascent
σ²_msg_RNA⁺  = σ²_imp_RNA⁺ + 1/n_rna⁺_src
```
A mature-dominated message (region 503: 87% mature) inherits the large `σ²_imp_mature` → low precision → it
cannot override the strand. This realizes "separate our confidence in mature vs nascent RNA" without splitting
the simplex. (First-order approximation; exact would be a delta-method combination — refine only if needed.)

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

**v1 — total-density global `σ²_imp` (the minimal fix, and the decision gate).** One global `σ²_imp`,
`σ²_msg = σ²_imp + 1/n_src` used directly. This is a net *simplification* of `_scan`. **Hypothesis: v1 alone
collapses region 503 (precision 1027 → ~0.3, `f_g` 0.83 → ~0.63) and the whole runaway-precision tail.** If
the v1 A/B shows the capture-ON mature-AMBIG error resolving to near capture-OFF levels, v2 may be unnecessary
or reduced to a small refinement.

**v2 — per-component `σ²_imp` + staged refit (§3, §4).** Only if v1 leaves error where mature/nascent
reliability diverge. Adds the per-component estimator, the AMBIG exclusion / live-edge filtering, the RNA
density-share blend, and the refit-2 / pass-3 loop.

## 7. Implementation plan

**`bp_solver.py`**
* New `adjacent_disagreement_variance(chain, boundary_substrate, dens, counts) -> float` — the §2 estimator
  over adjacent boundary↔region pairs. `dens`/`counts` are per-node log-density and fragment count; `resid`
  computed as `log ρ_src − log ρ_dst` via the **same** density→frame path `_scan` uses (assert equality against
  a `_scan` value on a sample edge — frame consistency is silent-failure-prone). v2: a per-component variant
  taking the solved component densities + component counts, with the AMBIG / live-edge masks.
* `_scan` — replace the three `base_var/s2_edge/pr` blocks (gDNA ~L1067, RNA-pos ~L1106, RNA-neg ~L1123) with
  the algebraic form `pr = n_src / (n_src·sigma2_imp[component] + 1)` (§2 — no `1/n_src`, no `if`, `pr=0` at
  `n_src=0`). The belief-combine (`fbg/fbp/fbn`, `vbg/vbp/vbn`) is unchanged.
* `node_sweep` — accept `sigma2_imp` (scalar in v1; a `{gdna, nascent, mature}` struct in v2) and thread to `_scan`.

**`calibrate.py`**
* Compute the total-density `σ²_imp` before PASS-1 (`_sweep(None)`, ~L181); pass it into both `_sweep` calls
  (v1 reuses the one scalar). v2: insert Refit-1 beside the KDE fit (~L186–195, excl AMBIG) and Refit-2 after
  PASS-2, with the pass-3 convergence/monotone-sharpen check.

**`config.py`**
* A single toggle to gate the new precision path for the A/B (e.g. `CalibrationConfig.sweep_precision_mode`).
  No numeric hyperparameter (v1 has none).

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
