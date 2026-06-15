# Per-node gDNA/RNA deconvolution — a Bayesian hierarchy (strand likelihood × count prior × global gDNA foundation)

**Status:** design + implementation plan (2026-06-15). Approved in principle (user). The effective-length
contraction redesign is being written separately (`new_eff_len_plan.md`); it shares the per-region gDNA
density with this design but is **not** a barrier to implementing this one.

**Goal.** Replace the two competing per-node combines —
- the fusion blend `g = w·g_strand + (1−w)·g_count`, `w = I/(I+I₀)` (`strand_deconv._deconv_per_node`), and
- the sweep's `β`-penalty + w-gated Jeffreys (`simplex_sweep._local_loglik`) —

with **one** principled per-node solve so the propagation pipeline is ≥ the fusion **everywhere** and wins on
complex (overlapping opposite-strand) loci. Then flip `use_propagation` to default and **delete the fusion
region-deconv** — one pipeline, no redundant paths.

---

## 1. Why

- Two region-deconv paths behind `use_propagation` is confusing debt. The propagation path has the higher
  ceiling (it resolves AMBIG loci via per-strand odds), so it should be the survivor.
- The fusion's per-node combine works but (a) **discards the AMBIG tilt** (forces `w=0` ⇒ count-only at
  AMBIG), (b) needs a **brittle magic `β`/`I₀`**, and (c) the sweep's variant adds a U-shaped **Jeffreys
  vertex-push** that invents phantom gDNA in zero-DNA libraries and at anchor-starved AMBIG nodes.
- Deciphered discrepancy (2026-06-15 scoreboard, `/tmp/sweep_scoreboard.py`, `/tmp/sweep_beta_test.py`):
  - **zero-DNA phantom** = coarse lattice **+ the Jeffreys vertex-push** (K=20→60 cut it 16.7k→1.1k; the
    residual is the prior pushing unanchored nodes to a vertex);
  - **ss0.50 regression** = the **count combination** (the soft `β`-penalty under-weights the count vs
    fusion's direct `g_count`; β=10→1000 pulls nrna 299k→251k → fusion 243k);
  - **flagship ss0.99** = **EM-bound** (the leak is the shared eff-len, not the calibration — see
    `em_gdna_strand_likelihood_fix.md` §9 and the separate eff-len redesign).

## 2. Principles (user, 2026-06-15)

1. **Strand is unbiased and reliable** → the priority signal; use it wherever it is informative.
2. **Count is valuable but unreliable** → a *fallback*. Its reliability is **local** (it tracks hybrid-capture
   probe-to-boundary distance, which we cannot observe) and varies node-to-node: sometimes highly accurate,
   sometimes not. It **must** take over when the strand is silent, but it must not corrupt the strand where
   the strand speaks. We *down-weight* the count by its local reliability; we do **not** bias-correct it (the
   capture bias is exon-vs-intron-specific and may not transfer).
3. **The global gDNA baseline is the foundation** → numerical stability + edge cases (truly-unanchored
   nodes). It must not overtake strand or count.

## 3. The model — a hierarchy, not a blend

For each node the unspliced mass is one pie `(f₊, f₋, f_g)` on the 2-simplex. The per-node posterior is a
**likelihood updating a prior**:

```
posterior(f₊,f₋,f_g)  ∝   L_strand(f | tilt)  ·  L_odds(f | neighbours)   ·   Prior(f_g)
                          └────── unbiased, PRIORITY ──────────────────┘      └─ count → global, FALLBACK ─┘
```

The strand-information hierarchy (**strand > count > global**) is **emergent**, not imposed by a weight: a
*sharp* likelihood overwhelms the prior (strand wins where informative); a *flat* likelihood cedes to the
prior (count, then global, take over). This is ordinary Bayesian updating — the reason the fusion needed a
hand-set `β` is that it blended two *point estimates* instead of writing the count down as a prior with an
honest width.

### 3.1 The tilt makes the strand first-class at AMBIG

The observed plus-fraction is `p = ½·f_g + κ·f₊ + (1−κ)·f₋` (κ = `rna_sense_frac`). This is a **plane through
the simplex** that constrains `f_g` — it is *not* zero information at AMBIG:

- balanced node (`p ≈ ½`) ⇒ "gDNA-dominated **or** balanced f₊/f₋",
- tilted node ⇒ "RNA-dominated on one strand."

The second equation needed to locate the point on that plane is the **propagated per-strand odds**
(`L_odds`) from neighbours. So an AMBIG node is solved by strand+propagation, **not** forced to count-only
(`w=0`). Only a node that is *both* balanced *and* has no odds inflow is genuinely strand-silent → it falls
to the count, then the global baseline.

### 3.2 The exact combine (on the lattice) — one normalization, precision-weighted shrinkage

All terms are added in log-space and the posterior is normalized **once** over the 2-simplex lattice
(`_simplex_lattice`); they are **not** separately-normalized likelihoods multiplied together:

```
log post(f) = L_strand(f|tilt) + L_odds(f|neighbours)
              − ½·τ_count(f_g − μ_count)²          (count prior)
              − ½·τ_global(f_g − μ_global)²         (global prior)      ;  then softmax-normalize over the lattice
```

Because the lattice restricts `f_g ∈ [0,1]`, each quadratic acts as a **truncated Gaussian on the simplex**
(a bounded, normalizable prior — the same device `solve_node` already uses; this is the resolution of "a
raw Gaussian on a fraction"). The two quadratics **compose into a single Gaussian prior** with

```
τ_prior = τ_count + τ_global ,   μ_prior = (τ_count·μ_count + τ_global·μ_global) / (τ_count + τ_global)
```

— i.e. a precision-weighted mean: the count governs where it is precise, the global baseline where the
count is thin. *This* is the "empirical-Bayes shrinkage" of §4c, made precise (it is one prior, not two
pieces of evidence, so there is no double-counting even though both means descend from the per-region
density).

## 4. The terms

### 4a. Strand likelihood + propagated odds — the priority (unbiased)

As already implemented: `simplex._mixture_strand_loglik` (the 3-component Gaussian-moment mixture on the
tilt; reduces *exactly* to the 2-component `strand_loglik` at a single-strand node) + the edge coupling on
the per-strand log-odds `log(f_c/f_g)` (`simplex_sweep._edge_logphi` / `_sweep_chain`). Its precision is the
intrinsic Fisher information (`N·(2κ−1)²` plus the propagated odds curvature). **Never gated to `w=0`.**

### 4b. Count prior — per-node fallback with derived, local precision

A Gaussian prior on `f_g`:  `f_g ~ N(μ_count, σ²_count)`, contributing `−½·(f_g − μ_count)² / σ²_count` to
the log-posterior.

- `μ_count = count_gdna_frac` (the local density→fraction the count module already produces).
- **Precision `1/σ²_count` is the per-node count precision — this replaces `β`.** It is local because the
  count's reliability is local:
  - the number of boundary-crossing fragments is the **observable proxy for probe-to-boundary distance** —
    many crossings ⇒ probe near ⇒ reliable ⇒ tight `σ²_count`; few/none ⇒ probe far ⇒ unreliable ⇒ wide;
  - a **count-unobservable** node (no crossings, imputed) gets imputation-level (low) precision;
  - the depletion bias is *self-correlated* with thinness (far-from-probe ⇒ depleted **and** thin ⇒ already
    down-weighted), so we are protected without estimating the bias.
- **Stage 1:** use the existing Poisson/observability variance (`density_model` /
  `_compute_side`'s `count_gdna_frac_var`) as the precision placeholder.
- **Stage 2:** swap in the **nonparametric `var~mean`** model (capture-robust) for `σ²_count` as a function
  of the node's count level.
- **Down-weight only — no bias-correction.**

### 4c. Global gDNA baseline — the foundation

A weak Gaussian prior on `f_g`:  `f_g ~ N(μ_global, σ²_global)`, with, per node,
`μ_global[i] = clip(ρ_global · region_eff_len[i] / mass_unspliced[i], 0, 1)` — identical in form to the
count module's own `count_gdna_frac` (`density_model.py:347-349`) with `ρ_global` substituted for the local
density. (So `deconv_regions_sweep` must now receive `region_eff_len`; `mass_unspliced` is already inside
the sweep via `substrate.contained`.)

- `ρ_global` = a **single pooled gDNA density** from observable nodes. **This already exists** — it is
  `density_model.node_gdna_density`'s `global_density` (`density_model.py:338-339`:
  `Σ contained_gdna[observable] / Σ region_eff_len[observable]`, `0.0` when there are no observable nodes,
  which **is** the init-uniform / zero-DNA limit). Stage 1 **exposes** this value (on `NodeDensity` or a thin
  wrapper) and passes it as `ρ_global` — it does **not** build a second estimator (the §11 one-density
  rule). It is computed **PRE-sweep** from observable nodes, and is distinct from
  `derive.gdna_density_global` (which is POST-sweep over **all** nodes — feeding that back in would create
  the density→deconv→density loop the acyclic design forbids; that loop is the deferred EB enhancement §6).
- **Init = uniform** — no capture-structure assumption (panels target a fraction of transcripts; we cannot
  assume exons are captured, so we cannot build capture structure from the seeds without biasing the model).
- `σ²_global` = a **robust between-observable-node spread** of the per-node density estimate, e.g.
  `σ_global = 1.4826·MAD(ρ_node over count-observable regions)` (weak — must not overtake count/strand).
  **Zero-observable limit:** when there are no observable nodes `ρ_global=0` and the spread is undefined →
  use a **wide-but-finite** fallback width, so an unanchored node is *gently pulled toward 0*, not pinned.
  Stage-1 validation for this width: (a) zero-DNA unanchored AMBIG `f_g → ~0`; (b) a real low-gDNA
  single-strand node (`f_g~0.1`, `N~10`) is **not** biased downward vs fusion.
- **This is the zero-DNA fix.** In a pure-RNA library `ρ_global ≈ 0`, so `μ_global ≈ 0`, so an unanchored
  node settles at `f_g ≈ 0` — **no phantom gDNA**, replacing the Jeffreys vertex-push that invented `f_g=1`.
  "The algorithm depends on knowing the global baseline gDNA level" (gDNA is present everywhere at an
  estimable baseline).
- The product `count_prior × global_prior` **is** the empirical-Bayes shrinkage (a precision-weighted mean
  of `μ_count` and `μ_global`): the count governs where it is observed and reliable, the global baseline
  governs where the count is thin/absent. No separate shrinkage step.

## 5. What this removes

| removed | replaced by |
|---|---|
| `count_trust_beta` / `I₀` for the region solve (magic constant) | the **derived per-node count precision** (4b) |
| the w-gated Jeffreys vertex-push (`jeffreys_w`) | the **informative global gDNA prior** (4c) |
| the hybrid "compute fusion + sweep, select per node" | the **unified per-node solve** (strand dominates at single-strand too) |
| the separate fusion `deconv_regions` path | nothing — deleted once the unified solve validates |

## 6. The empirical-Bayes outer loop — Stage 3, an ENHANCEMENT (deferred)

After a pass we have a gDNA posterior at **every** node, not just the depleted seeds. Re-fit `{ρ_global,
var~mean}` on **all** nodes and re-run:

- **init uniform → iterate** to escape the seed-node bias (seeds are depleted under capture and only a
  subset are captured);
- **the strand teaches the count, through the loop:** strand-resolved posteriors dominate the global re-fit,
  so the global count model becomes strand-calibrated, which then helps the still-unresolved nodes next pass.
  (This is the right home for "calibrate count against strand" — the EB loop, not a per-node residual hack.)

This turns calibration from *acyclic single-pass* into an **acyclic per-node/propagation solve inside a
bounded EB outer loop** (a few seed-/strand-anchored passes, not open EM). **Deferred:** implement the
single pass first; iteration is a documented enhancement. Convergence (the global prior feeds the per-node
solve which feeds the global prior) is the thing to watch when we add it.

## 7. Staging (implementation order)

1. **Stage 1 — the β-free per-node combine, single pass.** Implement 4a × 4b × 4c with a **uniform** global
   prior and the **existing** Poisson/observability count precision. Remove `β` and the Jeffreys-push from
   the sweep. **Validate the structure works before adding learning:** the unified solve must be ≥ fusion on
   the suite and win the complex battery, with the unified solve handling single-strand nodes (so the hybrid
   compute-both is retired).
2. **Stage 2 — nonparametric `var~mean`** for the per-node count precision (4b).
3. **Stage 3 — the EB outer loop** (§6), with convergence guards.

## 8. Validation gates (every stage)

- **Scoreboard** (`/tmp/compare_pipelines.py`, fusion vs new): flagship ss0.99, ss0.50, zero-DNA
  (`gdna_none`, both `nrna_none` and `nrna_rnd`), capture-off — the new solve **≥ fusion** (no regression).
  Zero-DNA must be clean (phantom gDNA ≈ fusion's ~30, not thousands).
- **Complex battery** (`scripts/debug/complex_loci_benchmark.py`): TOTAL 2D/1D ≤ 0.92, and the
  anchor-starved families (`full_single`, `multistrand_*`, `anchor_sj_neither`) should **improve** (the
  global-baseline fallback replaces the vertex-push that made them lose).
- **Full suite** (`pytest tests/`, 1106) green; goldens unaffected until the default flip; conservation
  `gdna + rna = total` holds.

## 9. Integration points

- `simplex_sweep.py` — rewrite `_local_loglik`: drop `jeffreys_w`, drop the bare `β` count penalty; add the
  **count-prior** term (precision = per-node count precision) and the **global-prior** term (`μ_global`,
  `σ²_global`). `deconv_regions_sweep` computes `μ_global` per node and passes the count precision.
- `calibrate.py` — the `use_propagation` branch becomes **just** the unified `deconv_regions_sweep` (drop
  the hybrid compute-both/select); compute `ρ_global` and hand it + the count precision to the sweep.
- A small **global-baseline estimator** `ρ_global` (pooled gDNA density from observable nodes; init uniform)
  — new helper (in `derive.py` or `density_model.py`).
- `config.py` — `count_trust_beta` becomes unused by the region solve (keep the field deprecated for now to
  avoid churn; remove with the fusion deletion). `sweep_n_grid` stays.
- Count precision Stage 1: reuse `density_model`'s `count_gdna_frac_var` (regions) and `_compute_side`'s
  (sides).

## 10. Open questions / risks

- **Single-strand lattice resolution.** The K=60 simplex line carries 61 `f_g` points vs the fusion's 200.
  If single-strand nodes show a residual vs fusion, options: raise `sweep_n_grid`, or a 1-D fast path for
  single-strand nodes (they don't need the 2-simplex). Validate first.
- **`σ²_global` magnitude.** Must be weak enough not to overtake count/strand, strong enough to anchor the
  truly-unanchored AMBIG. Derive from the between-node spread; validate it doesn't bias real-gDNA nodes.
- **Count precision at imputed nodes.** How wide. Stage 1 uses the existing variance; Stage 2's `var~mean`
  is the principled answer.
- **EB convergence** (Stage 3) — deferred; seed+strand anchoring should bound it.

## 11. Relationship to the effective-length redesign

`ρ_global`, the count prior, and the new eff-len contraction all consume the **same per-region gDNA
density**. Keep one density definition so the gDNA baseline used here and the enrichment used by the eff-len
contraction stay consistent. Coordinate when `new_eff_len_plan.md` lands.

## 12. Review hardening (adversarial review, 2026-06-15) — ordering + watchpoints

Verdict: **ready-with-doc-fixes**. The fixes above (§3.2, §4c) are folded in. Remaining must-knows:

**Strict implementation ordering (a correctness dependency, not a preference):**
1. **Global prior in** (§4c) → validate zero-DNA. It is the new stabilizing prior.
2. **Per-node count precision in**, drop the scalar `count_trust_beta` from the sweep call (§4b).
3. **Jeffreys vertex-push out** (`_local_loglik` `jeffreys_w` term) → re-validate zero-DNA. **Do not remove
   it before step 1** — without the global prior, single-strand `f_g≈0` nodes would have *no* stabilizing
   prior (the `solve_node` gDNA Dirichlet pseudo-count is absent from `_local_loglik`).
4. Only after all gates green: **flip `use_propagation` default and delete the fusion `deconv_regions` + the
   hybrid select** — in a *separate* commit from the math change (keep the diff bisectable;
   `count_trust_beta` stays as a deprecated unused field until the fusion path is removed).

**Scope:** the sweep redesign targets **contained regions only**. Boundary **sides keep the fusion
`deconv_sides`** (the flux transport is unchanged) — they are not swept.

**Zero-DNA mechanism is node-class-specific (validate them separately):** at single-strand nodes the sided
spliced lower bound (`_local_loglik:65-67`) already pushes `f_g` down; the global prior is a secondary
anchor there. At **AMBIG** nodes the sided spliced bound is **zeroed** (`simplex_sweep.py:150-151`), so the
global prior is the **primary** new anchor for balanced + no-odds-inflow AMBIG. Attribute phantom reduction
separately for single-strand vs AMBIG.

**Watchpoints during Stage 1:**
- **Precision clipping (critical).** `precision = 1/count_gdna_frac_var` blows up at no-anchor nodes
  (`v_rel=1`, `count_gdna_frac≈0` ⇒ `var→0` ⇒ `precision→∞`), which would wrongly *pin* `f_g` to a
  degenerate count. **Cap the count precision** so no-anchor nodes get a *low* precision (wide prior) and
  defer to the global prior, as §4b intends.
- **ss0.50 (κ=½, fully unstranded):** every node has `w=0` so the count prior alone sets `f_g`. Confirm the
  per-node precision `1/count_gdna_frac_var` is the same order as the old effective `β=10`, else nrna is
  under-estimated; gate on ss0.50 nrna ≥ fusion.
- **K=60 resolution:** before retiring the hybrid, confirm single-strand on the K=60 simplex ≥ fusion at
  `n_grid=200`; at zero-DNA with a tight global prior at `f_g=0`, check the K=60 (61 `f_g` points)
  resolves small `f_g` (target median `f_g < 0.01` for truly-zero AMBIG; test K=60 vs K=100 if needed).
- **AMBIG-only / count-starved chains** (`full_single`, `multistrand_*`, `anchor_sj_neither`): the design
  predicts these *improve* (global fallback replaces the vertex-push). If they regress, `σ²_global` is too
  wide / the global prior too weak at count-starved nodes.
- **Splice 3-term variance:** the 3-term upgrade (`calibrate.py:197-208`) rewrites `region_count_frac`
  (= `μ_count`) but not its variance — acceptable (a magnitude refinement), but verify the upgraded regions
  are not over-trusted.
- **Conservation** (`gdna+rna=total`) must hold per node after the rewrite.

**Verified false positive (no action):** `_fg_median` (`simplex_sweep.py:83-92`) is correct — the
full-posterior cumsum sorted by `f_g` *is* the marginal CDF (duplicate-`f_g` points are adjacent after
argsort). No rewrite needed.

## 13. Implementation finding (2026-06-15) — the prior must be NODE-CLASS-DEPENDENT

Implementing the literal "remove the Jeffreys, replace with the global prior" (§5) **regressed every
condition** (flagship nrna 88→115k, ss0.5 +73k, capture-off 5→17k; `/tmp/compare_pipelines.py`). Cause: the
Beta(½,½) Jeffreys is the **legitimate strand reference prior at single-strand nodes** — the fusion uses it
ungated (`strand_deconv.py:127`) and *relies* on it to concentrate `f_g` against the mixture's
overdispersion-induced spread (without it, balanced pure-gDNA nodes under-call gDNA → leak). It is *not* a
"vertex-push" there (the likelihood already favours that vertex). The Jeffreys is only harmful at **AMBIG**
nodes, where the U-shape pushes a flat-likelihood node to the `f_g=1` vertex (phantom gDNA).

**Refinement (implemented):** the prior is **node-class-dependent** (`_local_loglik`, `strand_obs` mask):
- **single-strand (TS_POS/NEG):** Beta(½,½) Jeffreys reference prior (ungated), as the fusion does.
- **AMBIG / intergenic:** the **global gDNA prior** (`μ_global`, weak `τ≈1`), replacing the harmful Jeffreys.

This supersedes §5's "remove the Jeffreys wholesale" and §4c's "global prior replaces the Jeffreys
everywhere" — the global prior replaces it **only at AMBIG**.

**Scoreboard after this refinement** (fusion vs unified sweep, no hybrid):

| condition | fusion | unified sweep | gap | remaining cause |
|---|---|---|---|---|
| flagship ss0.99 nrna | 87,574 | 96,402 | +8.8k | AMBIG + K=60 (EM-bound, second-order) |
| ss0.50 nrna | 243,446 | 272,810 | +29k | **count combination (β=10, Step 2 next)** |
| zero-DNA gdna phantom | 27 | 1,155 (+816 nrna) | +~2k | AMBIG global-prior strength |
| capture-off nrna | 5,240 | 8,427 | +3k | AMBIG |

Recovered from the broken "no-Jeffreys" state; now ≈ the earlier hybrid but as a **single unified solve**
(no compute-both). **Next: Step 2 — replace the scalar β count penalty with the per-node count precision**
(should close the ss0.5 +29k), then assess the AMBIG global-prior strength for the zero-DNA +2k.

## 14. Step 2 (per-node count precision) — done; surfaces the observed-vs-imputed reliability gap

Implemented (β eliminated): the count prior now uses the per-node precision `τ_count = 1/v_rel`
(`density_model` exposes `count_rel_var`; `need_count_variance` also turns on for `use_propagation`;
`calibrate` passes `1/max(v_rel, 1e-6)` to the sweep). `v_rel` = `1/N_own` (observable), LOESS-floored
(imputed), `1.0` (no-anchor) — bounded, no `1/var` blow-up (W3 resolved by using the *relative* variance,
not the absolute σ²).

Scoreboard (fusion vs unified sweep):

| condition | fusion | sweep (β) | sweep (τ=1/v_rel) | note |
|---|---|---|---|---|
| flagship ss0.99 | 87,735 | 96,402 | 96,866 | EM-bound, ~flat |
| ss0.50 | 242,973 | 272,810 | **258,958** | improved +29k→+16k ✓ |
| capture-off | 5,026 | 8,427 | **4,234** | now ≈ fusion (slightly better) ✓ |
| zero-DNA gdna phantom | 27 | 1,155 | **3,414** | regressed ✗ |

**Finding — the count's `1/v_rel` precision over-trusts the IMPUTED count.** ss0.5 and capture-off improved
as designed, but zero-DNA regressed: at imputed AMBIG nodes `1/v_rel` exceeds the old `β=10`, so the count
is trusted *more* — and in a pure-RNA library the imputed count is *confidently wrong* (it imputes nonzero
gDNA from noisy neighbours), which the weak global prior (`τ=1`) cannot override. `v_rel` captures
*sampling* uncertainty (does the node have enough crossings) but not the *kind* of estimate: a **directly
observed** count (own crossings) is data; an **imputed** count (no crossings) is a guess from the density
field and should defer to the global baseline.

**Refinement TRIED and REVERTED — cap imputed-count precision at the global `τ`.** Hypothesis: an imputed
count is a guess, so cap it low and let the global baseline govern. **It regressed everything** (flagship
96.8→104k, ss0.5 mrna −66k as the exon-blind global baseline over-called gDNA, zero-DNA unchanged at 3.4k).
Why it was wrong: the imputed count at **exons** is the splice-junction-upgraded fraction
(`region_count_frac`) — genuinely informative and **capture-aware**, *better* than the global baseline
there. The imputed count is only noise in **zero-DNA**, which is **locally indistinguishable** from real
signal — you can only tell globally that the count over-calls (it consistently disagrees with the strand),
which is exactly what the **Stage-3 EB loop** learns. So there is no clean Stage-1/2 local fix; the zero-DNA
residual is deferred to Stage 3.

## 15. The count-precision trade-off — β-elimination needs Stage 2's var~mean (β=10 is the Stage-1 placeholder)

Three per-node-precision variants were tried against BOTH the simple scoreboard and the complex battery:

| count precision | ss0.50 | complex 2D/1D | verdict |
|---|---|---|---|
| **β = 10 (flat scalar)** | +29k | **0.95** (wins) | best complex; ss0.5 worst |
| `1/v_rel` (per-node) | **+16k** | 1.31 (loses) | best ss0.5; breaks complex |
| `1/v_rel` capped at global τ | worse (global over-calls gDNA at exons) | 1.31 | reverted |
| `(1−w)·(1/v_rel)` | +16k | 1.22 (loses) | reverted |

**Root cause (the genuine impasse).** The count's role is **node-class-dependent**: at **count-observable**
nodes (intergenic/intron — no competing strand) it is the *primary* signal and wants *high* precision
(this is the ss0.5 win); at **AMBIG** nodes it must *yield* to the strand+propagation (the user's priority
principle — else the complex win is lost). But the imputed-count precision `1/v_rel` **inherits the
anchor's confidence, not the node's own data**: a low-N AMBIG node imputed from a confident high-N anchor
gets *high* precision and over-rides the propagation. A flat modest β can't over-ride the propagation
(`10 ≪ N`) and modestly reinforces the strand at single-strand — so it is the Goldilocks placeholder, but it
leaves ss0.5 worse than fusion.

**Resolution = Stage 2, not a Stage-1 hack.** The principled per-node precision is **observability-gated**:
`1/v_rel` (full) at count-observable nodes, *modest* at AMBIG/imputed nodes (count yields to the
propagation). That is exactly what the **nonparametric var~mean** delivers when fit on the nodes that have
data and *not* over-extended to AMBIG — i.e. Stage 2. A flat cap at the global τ was too low (the
exon-blind global baseline over-called gDNA in real-gDNA libraries); the var~mean gives a *calibrated*
modest precision instead of a hard cap.

**Decision (2026-06-15):** **Stage 1 lands with β=10 as the count placeholder** (node-class prior + global
prior + β count). It wins the complex battery (0.95) and is the principled foundation (the Bayesian
hierarchy, β-free *structure*); the literal β-elimination moves to **Stage 2** (observability-gated
var~mean), which also targets the ss0.5 +29k. Per the user's steer, the simple-suite residuals
(flagship +9k EM-bound; ss0.5 +29k → Stage 2; zero-DNA ~2k → Stage 3) are accepted milestones, each with a
named later-stage fix.

## 16. What is MISSING to make the calibration better — the model specs

The 3-tier structure is complete; what's missing is that **every precision/variance in it is still a
placeholder, and they are all the same kind of object — a `var~mean(density)` curve**. The nonparametric
machinery already exists (`density_model.density_variance_curve`, a robust log-log LOESS that is explicitly
capture-bimodality-aware; `rna_variance.py` already reuses it) — it is just scattered and used only for the
FP-quantile, not made authoritative. Three models to develop:

### Model A — the gDNA-density `var~mean` curve (the KEYSTONE; serves count precision AND the global variance)

One curve `σ²_ρ(ρ)` = the variance of a gDNA *density* estimate at density level `ρ`, fit over the
count-observable/seed nodes. Variance observables already present: a node's own Poisson `ρ²/N`, and the
two-anchor disagreement `¼(d_L−d_R)²` at mean `½(d_L+d_R)`. **Local, not global** — the log-log LOESS keeps
the depleted (low ρ, low var) and enriched (high ρ, high var) capture modes distinct, so a single global
variance never has to span that range (your point exactly). It feeds **two** tiers from one model:

- **Count precision (Stage 2):** `τ_count = 1 / [σ²_ρ(ρ_node) · (eff_len/mass)²]` (the `(eff_len/mass)²`
  propagates the density variance to the gDNA-*fraction* scale). **Observability-gated:** the curve is the
  count's reliability where the count is DATA (observable own-crossing nodes); at IMPUTED/AMBIG nodes the
  precision must stay modest so it does not over-ride the propagation (the §15 failure). This is what closes
  **ss0.5 +29k** while preserving the complex win.
- **Global-prior variance (your new piece):** replace the flat `τ_global=1` with
  `σ²_global = σ²_ρ(ρ_global) · (eff_len/mass)²`. **This is very likely the real zero-DNA fix:** in a pure-RNA
  library the seeds agree the density is ≈0, so `σ²_ρ(low)` is *small* ⇒ the global prior is *tight* at
  `μ_global≈0` ⇒ unanchored/AMBIG nodes are pinned to `f_g≈0` (no phantom). Under capture the same curve
  makes the prior appropriately *wide* at enriched densities (defers to count/strand). One model fixes both
  **ss0.5 and zero-DNA**.

What's needed from you: the robust `var~mean` *form* — the LOESS span (the one smoothing knob,
`_LOESS_SPAN=0.4`, CV-selectable), how to pool the variance observables, and how it degrades with very few
seed nodes. The plumbing (expose the curve, evaluate it at node-density and at `ρ_global`, the
observability gate) is mine.

### Model B — the per-edge RNA odds-coupling variance `q_rna` (Tier-1 propagation)

The propagation couples per-strand log-odds with a **scalar `q_rna=0.25` placeholder**. It should be a
**per-edge** `var~mean`: the RNA-density disagreement *across* each same-strand-exon boundary
(`½(d_a_right + d_b_left)`, `¼(d_a_right − d_b_left)²` from `substrate.left/right` spliced flux — all
available). `rna_variance.py` (Phase-2a) already fits the per-region RNA `σ²_RNA(μ)` with the same LOESS;
the missing piece is the *per-edge* lookup. A refinement the battery's `staircase`/`gradient` families
imply: the coupling should be **depth-aware** (the correct per-strand odds varies where transcript depth
ramps along the chain). What's needed from you: whether per-edge pooling has enough samples (likely
library-wide pooling), and the depth-aware form.

### What is NOT a missing model (so you don't build for it)

- **Multi-transcript-per-strand** (the *largest* complex-battery failure class — `double_plus_skew`,
  `triple_plus_stack`, `buried_sibling`, `double_cover`, `rare_sibling_leak`): the sweep propagates ONE
  odds value per strand per edge and cannot represent 2+ transcripts' odds. But **calibration only needs the
  per-strand RNA *total*** (gDNA-vs-RNA); the **EM resolves isoforms downstream**. So these are largely
  *benchmark artifacts* of scoring per-node `f_g` against an isoform-aware oracle — not a calibration model
  gap. (If we want the battery to stop penalizing them, score the per-strand RNA total, not isoforms.)
- **Anchor-starved AMBIG** (`full_single`, `anchor_sj_neither`, `nest_anchor_starved_core`, orphan
  islands): genuinely under-determined *locally* — no strand tilt resolves `f_g`, no neighbor odds, no own
  count. The **global prior is the floor** (the best possible); Model A's tighter, data-derived variance
  makes them **degrade gracefully** (toward `ρ_global`) instead of vertex-pushing, but cannot *resolve*
  them. Not a model to build — an inherent identifiability limit.
- **flagship ss0.99 +9k:** EM-bound (the shared eff-len), calibration-invariant — your eff-len redesign is
  the lever, not the calibration.

### The Stage-3 EB loop is the *mechanism*, not a model

Re-fitting `{ρ_global, σ²_ρ}` on the full POST-sweep posterior (not just the depleted seeds) and re-running
is what lets Models A/B escape the seed-only/init-uniform bias and learn the count's *global* reliability
("the strand teaches the count" — strand-resolved nodes dominate the re-fit). It refines A/B; it is not a
separate variance model.

**Net:** develop **Model A** (the gDNA-density `var~mean`) first — one curve closes ss0.5 *and* zero-DNA
and supplies your global variance. **Model B** (per-edge `q_rna`) is the Tier-1 refinement. Everything else
is inherent or downstream.

## 17. Stage 2 — var~mean: empirical findings, the DIRECT/IMPUTATION split, and the consolidation

**Empirical (flagship gdna300/ss0.99/cap_on, `/tmp/seed_varmean_plot.py`):** the seed/observable nodes give
**1038 var~mean points** (101 intergenic + 470 intron + 467 boundary-imputed exon) — **plenty for a LOESS;
keep the existing dataset** for iteration (a *targeted-panel* sim is needed later only to stress the
locality/adaptive-span case, which this whole-transcriptome-like sim does not exhibit). The data is
**bimodal** as predicted (a depleted off-target cluster at μ≈0.01 and an enriched on-target cluster at
μ≈0.5–10). The trend is a clean **power-law β≈2.06** (var ∝ μ² ⇒ ~constant relative variance) through the
bulk; LOESS span 0.4–0.75 tracks it; isotonic (PAVA) is monotone but step-jumpy at the high-μ enriched
scatter. So on *this* dataset a single span ≈ the power-law; the adaptive-span need is real but only bites
on small targeted panels (a separate stress sim).

**The DIRECT vs IMPUTATION split (the consolidation backbone, and the §15 fix).** Two var~mean models,
distinguished by *what the variance is of*:
- **DIRECT** — variance of an OBSERVABLE node's own density estimate (own-count Poisson + any excess
  dispersion). Used for: count precision at count-observable nodes, and the **global-baseline variance**
  σ²_global (Model A-(B)). Built from intergenic/intron regions + exon-intron/exon-intergenic boundaries.
- **IMPUTATION** — variance of *predicting* a node from its boundaries, learned two ways:
  (i) boundaries observable, region not → the two boundary estimates are replicate measurements
  (¼(d_L−d_R)²; what `density_variance_curve` already fits); (ii) all three observable → the region is the
  *known answer*, so we measure the boundary→region prediction error directly (a true validation set).
  Used for: count precision at IMPUTED/AMBIG nodes. **This is the §15 fix** — imputed nodes were
  over-trusted because they borrowed the DIRECT precision (the anchor's confidence); the IMPUTATION
  variance is properly humble (larger), so imputed AMBIG defers to the propagation instead of over-riding it.

**Method (research, voom/DESeq2 precedent):** keep it data-driven and monotone. Options that fit
"monotone + locality-adaptive + simple": (a) **adaptive-bandwidth (kNN) LOESS + isotonic post-pass** (the
recommendation — locality from kNN, monotonicity from PAVA); (b) **monotone P-spline / SCAM** (single
penalty, monotone by construction); (c) the **power-law β≈2** as a robust parametric baseline with a local
fallback (DESeq2's dispersion-trend strategy). The fixed-span LOESS is exactly voom's known over-smoothing
risk on small panels. **The span/bandwidth calibration (monotonicity + locality) is the open task** the
user is researching.

**Consolidation (part of Stage 2).** The var~mean code is scattered (`density_variance_curve`,
`_count_fraction_variance`, `_loess`, `rna_variance.py`, `count_rel_var`). Consolidate into one
`variance_model.py`: a single core var~mean fitter (the chosen monotone/adaptive method) + the DIRECT and
IMPUTATION builders, consumed by (1) count precision, (2) global-baseline variance, (3) the per-edge RNA
odds-coupling `q_rna` (Model B). One authoritative module; `rna_variance.py`/the count-variance helpers
become thin callers.

**Iterative build (Stage 3 EB):** pass 1 fits DIRECT/IMPUTATION from the depleted *seeds* only; pass 2+
re-fits from ALL deconvolved nodes (far more data, at the cost of imputation/strand-deconv noise) — this is
how the curves escape the seed-only bias and learn the true capture structure.

**BUILT (2026-06-15): the SCAM fitter** `variance_model.MonotoneVarMean` (monotone-increasing P-spline via
the `β=cumsum([α₀,exp(α₁..)])` reparameterization; GCV-λ; bisquare-IRLS robustness; pure numpy+scipy, no
new deps, no C++ — the fit is once-per-calibration, sub-second). Validated head-to-head
(`scripts/debug/scam_var_mean.py`): smooth + 100% monotone + fringe-flexible, beating LOESS (non-monotone
fringe) / isotonic (staircase) / power-law (rigid). 7 unit tests (`test_variance_model.py`); 238 calibration
tests green. **OPEN — the extraction split:** the first cut splits the triplet points by
`region_observable`, which is *confounded with mean* (DIRECT = intergenic/intron = depleted low-μ;
IMPUTATION = exon = enriched high-μ) — different μ regimes, not the reliability axis. The intended split
(user) is by *measurement method at the same node*: **IMPUTATION learned from the all-three-observable
case** (boundary predicts region, region = the known answer → the genuine prediction error), applied to the
region-imputed case. To refine before wiring.

## 18. Calibration diagnostics — optional dataframe outputs + companion plotter (decided: dataframes)

**Decision (2026-06-15): emit tabular dataframes, not plots** (a `--calibration-diagnostics` flag /
`CalibrationConfig`), with a companion plotting script. Rationale: keeps **matplotlib out of `rigel quant`**
(no plotting dep / import cost in the production tool), gives users the **raw calibration data** for their
own QC, and makes plots **reproducible from versioned data**. Essential for real-data runs.

Design:
- `CalibrationConfig.emit_diagnostics: bool` (CLI `--calibration-diagnostics`) + the run output dir.
- When on, calibration returns a `CalibrationDiagnostics` bundle of named dataframes written to
  `<out>/calibration_diagnostics/` (feather + TSV): `varmean_gdna` (the `MonotoneVarMean.to_dataframe()`
  points+curve, DIRECT/IMPUTATION tagged), `varmean_rna` (Model B), `per_region` (f_g, count_frac,
  g_strand, μ_global, strand_class, mass, observability), `global_scalars` (ρ_global, σ²_global, κ,
  overdispersions).
- Companion `scripts/plot_calibration_diagnostics.py` (or a `rigel` subcommand) reads the dir → renders the
  plots; **matplotlib lives only here.** The var~mean is the first diagnostic (the fitter already emits its
  dataframe). Built alongside the Stage-2 wiring (when the var~mean enters calibration).
