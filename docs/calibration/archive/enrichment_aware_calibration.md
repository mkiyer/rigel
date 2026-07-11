# problem

## imputation

we use imputation to help solve ambiguous nodes
imputation means we solve an ambiguous node using its neighbor nodes

hybrid capture creates enriched and depleted nodes.

when we model imputation, we try to predict a destination node from its neighboring source node. we model the variance of neighboring nodes.

with hybrid capture, the variance can be very high because going from a depleted region to an enriched region is a big change which creates a huge variance.

the high variance makes imputation imprecise.

imputation is also inaccurate because of the large changes in enrichment and depletion. neighbor nodes become a poor predictor because they made be enriched at different amounts


## global gdna mean

we depend upon a global gdna estimate to solve ambiguous regions. 

- hybrid capture creates enriched nodes and depleted nodes
- gdna density is not uniform
- the global gdna density mean of all nodes is inappropriate. hybrid capture creates a bimodal gdna distribution.
- with hybrid capture data, the global gdna variance is high

so the global gdna mean is incorrect (underestimates captured regions) and the variance is very high (low precision).


# solution is enrichment-aware models

calibration does not know whether we have hybrid capture data or uniform data. it must learn this.

the current calibration algorithm does not have any model in place for learning enrichment.

we must allow calibration to learn enrichment.

we must make no assumption initially.

at initialization, we have 'seed' nodes:
- intergenic regions are gdna and can set them
- intergenic <-> exon boundaries are gdna
- intronic regions are gdna + (sparse nascent RNA)
- exon <-> intron boundaries are gdna + (sparse nascent RNA)

### KEY POiNT ###
*We cannot fit global gdna model or imputation model from seed nodes*

**We cannot fit usable models at this stage, we MUST solve some nodes initially before we can fit models**

single-stranded nodes (strand=POS or strand=NEG) are directly solvable without additional information. they do not need a global gdna model

we NEED to solve single-stranded NODES initially. without them we cannot correctly model enrichment and we cannot correctly fit a capture-aware enrichment-aware imputation model or our global gdna model.

so after initializing nodes based on their transcript structure (signature), we need to let nodes solve themselves!

for simplicity we can assume we have strand-specific data (ss>0.99), which makes this much simpler. (with unstranded data, we are doing to need to allow initial imputation from solvable neighbor nodes because we won't have a strand model).

at this point we should have:

- gdna from 'seed' nodes
- *solved single-stranded nodes*

Now we can fit an initial global gdna model and an imputation model

# fit models

We must fit our initial gdna model and imputation models based on our solvable nodes.

as described above, "solvable" means:
- either a 'seed' node based on transcript structure (intergenic, intronic, exon-intron boundary, exon-intergenic boundary)
- or a single-stranded node (strand = pos, strand = neg)

Only (strand = AMBIG) nodes are withheld. These nodes require an intact imputation model and/or global gdna model to solve.

So now, how do we account for enrichment?

## enrichment ratio computation

1) compute the global gDNA density from solvable nodes (AFTER we initialize solve them!)
- (sum of estimated gDNA mass per node) / (node effective length)
- call this the `global_gdna_density`

2) compute node enrichment ratios:
- `enrichment_ratio = node_gdna_density / global_gdna_density`

The enrichment ratio now tells us if a node is enriched (density > global gdna density) or depleted (density < global gdna density)


## fit enrichment-aware models

the premise here:
- `node_gdna_density * node_enrichment_ratio = global_gdna_density`

now, can you help me?

we need to fit an enrichment-aware global gdna model!

we need to fit an enrichment-aware imputation model!

can you help me think about how to do this? 

publish a document that encompasses this discussion, the theory, the derivations.

we need to learn enrichment/depletion and then build models accordingly.

help me to derive this, present your best framework and let's develop this idea together!

---

# Part II — The framework (imputation-first, continuous enrichment)

> Status: design + empirical validation (no production code yet). Direction set with the user 2026-06-21:
> **imputation is the second line of defense and is developed FIRST**; the global gDNA model is the third line
> (a fallback for nodes the imputation cannot reach) and is developed after. Enrichment is modeled as a
> **continuous** quantity, not discrete classes. Every number is measured on the flagship
> (`gdna 3:1, ss=0.99, nrna_none, capture ON/OFF`) against the by-origin oracle; harness +
> diagnostics in `scripts/debug/enrichment_imputation_harness.py`, `enrichment_aware_check.py`,
> `class_density_spread.py`, `boundary_efflen_check.py`, `imputation_strength_debug.py`.

## A. Reframing: imputation is the second line, and it can stand alone

The three sources that set a node's `f_g` are ordered by the user's design:

1. **Strand** (first line, intrinsic) — solves single-strand nodes outright; flat/uninformative for AMBIG.
2. **Imputation** (second line) — solve an AMBIG node from its neighbours. **Developed here, first.**
3. **Global gDNA prior** (third line) — the population fallback, only where imputation has no informative
   neighbour. Developed later.

The reason to lead with imputation: it is **per-node** (it can tell *which* AMBIG exon is more or less enriched),
whereas the global is a single population number. Measured on the withheld AMBIG exons with **no global and no
strand**, a continuous imputation from the gDNA-clean boundaries recovers them with per-node correlation
**0.63**, versus **0.32** for a flat global-style baseline. The second line, done right, beats the third line —
and largely removes the burden the (broken) single-scalar global was carrying.

## B. The core object: a continuous enrichment transfer `ê(z)`

### B.1 Why the naïve enrichment ratio is degenerate (and the fix)

The natural definition `enrichment_ratio(x) = ρ_g(x) / global_gdna_density` is fine as a *descriptor* of the
solved nodes (depleted < 1 < enriched). But it **cannot drive imputation by normalization**, because

```
ρ_g(x) / enrichment_ratio(x) = ρ_g(x) / ( ρ_g(x) / global ) = global      for every node.
```

Every node's "de-enriched" density collapses to the same number — there is nothing left to impute. (Also: for a
**withheld AMBIG** node `ρ_g(x)` is the unknown we are solving, so its enrichment ratio is not even computable.)

**The fix:** estimate enrichment from a signal that is (i) **gDNA-clean**, (ii) **decoupled** from the node's own
`f_g`, and (iii) **observable for every node including withheld AMBIG ones**. That signal is the node's flanking
**boundary-crossing density**:

```
z(x) = max over the node's flanking BOUNDARY nodes of ( crossing unspliced mass / crossing eff-len )
```

Boundary crossings are gDNA-clean by structure (no spliced mature RNA crosses into an intron; nascent is sparse),
and `z` uses the *neighbours'* density, never the node's own gDNA — so the degeneracy is gone and `z` is defined
for AMBIG nodes. (`T = M/E`, the node's own total density, is **not** used to estimate enrichment — it mixes in
RNA, the RNA-confound. `T` enters only at the very end, in the `f_g = ρ_g/T` split.)

### B.2 The transfer, and what it subsumes

Learn a **continuous, monotone** transfer on the solved nodes (where `ρ_g` is known):

```
ê(z) = E[ ρ_g | z ]          (continuous; monotone non-decreasing; fit ONLY on solvable nodes)
```

Measured: a single smooth curve fits the solved single-strand exons with **R² = 0.946**
(`log ρ_g_interior = 1.13 + 1.08·log z`). Two consequences worth stating plainly:

- **`c_face` dissolves.** What I earlier proposed as a separate "edge-to-interior" factor (a boundary crossing
  reads ~`0.44×` the interior it borders) is simply the **offset and slope of this one curve** (`exp(1.13)≈3.1×`
  at slope 1.08). It is *what the transfer learns*, not extra infrastructure. The user's worry about modeling
  partial-probe-overlap is also handled for free: partial overlap → intermediate `z` → intermediate `ê`,
  continuously, with no class boundary to fall on the wrong side of.
- **No discrete classes, no K-selection.** Enrichment is one continuous function. There is no captured/depleted
  label to misclassify — the exact failure mode the user flagged. A node halfway between depleted and enriched
  gets a halfway prediction.

### B.3 Transfer validity (the key de-risking result)

The transfer is *trained* on single-strand nodes but *applied* to withheld AMBIG nodes — so it only works if the
`z→ρ_g` relationship is the same for both. Measured (oracle `ρ_g`):

| fit on | offset a | slope p | R² |
|---|---|---|---|
| single-strand exons | +1.126 | 1.077 | 0.946 |
| AMBIG exons (self) | +1.122 | 1.062 | 0.920 |

Intercept gap **0.004**; at matched interior `ρ_g≈16.7` the boundary `z` is identical (SS 6.00 vs AMBIG 6.18);
the single-strand-trained transfer predicts the AMBIG oracle `ρ_g` with **median ratio 0.99**. The transfer
generalizes from solved to withheld nodes essentially unbiased. This is the result that makes imputation-first
viable.

## C. Two layers: enrichment trend + de-enriched residual

The imputation is two stacked pieces — and this is what makes the off-capture behaviour correct:

1. **Enrichment trend `ê(z)`** — the systematic capture structure (active under capture; flat off-capture).
2. **De-enriched residual** `r(x) = ρ_g(x) / ê(z(x))` — what remains after the trend is removed. It is
   ≈ `ρ_base` (the underlying genomic copy-number rate), spatially **smooth**, so it imputes *reliably* between
   neighbours. This residual layer **is** the current identity-mean imputation, just operating in the normalized
   (de-enriched) space where the cliff has been removed.

A withheld AMBIG node is solved by:

```
ρ_g(x) = ê(z(x)) · r̂(x)              r̂ = neighbour-imputed residual (≈ 1 where no signal)
f_g(x) = clip( ρ_g(x) / T(x), 0, 1 )
```

For the flagship (uniform genomic rate) `r̂ ≈ 1`, so `ρ_g ≈ ê(z)` and the recovery is the transfer prediction.
Where the genomic rate genuinely varies locally (copy number, mappability), the residual carries it.

## D. The global is the third line (developed later)

The global gDNA model is only invoked where the imputation has no usable `z` — an AMBIG node with no
gDNA-clean flanking boundary, or an isolated node. There, fall back to the population level (the `z→0` / mean
limit of `ê`, i.e. the depleted baseline `ρ_base`, or a population `ê`). The global's job shrinks to *exactly*
what the user wanted: a fallback for the unsolvable, not the primary mechanism. We design it after the
imputation lands.

## E. Honest, source-governed precision

Strength is count-space and **destination-mass-independent** (no `(M/E)²` Jacobian). The key realization
(ablation 2026-06-22, `enrichment_ablation.py`): **the mean and the variance are two reads of ONE node-pair
fit**, and the precision must be **conditional on the level** (`var~mean`), not a global scalar.

- **Mean** `ê(ρ_src)` — the learned monotone transfer (replaces the identity message).
- **Variance** = the scatter of `ρ_dst` around `ê(ρ_src)`, fit by the EXISTING `var~mean`/`fit_offset` (which
  already subtracts the Poisson sampling floor → count-space-honest biological σ²_bio), **as a function of the
  level**. Measured: at enriched nodes (where the gDNA is) σ²_bio≈0.04 → **N≈24**; at depleted nodes σ²_bio≈4 →
  N≈0.2 (but they hold ~0 gDNA, so it doesn't matter). This **resolves the original "inert N≈1"**: the old code
  fit `var~mean` on the *identity* residual `(ρ_dst−ρ_src)`, pooling the systematic 3× edge→interior gap AND the
  depleted-node noise into one huge σ²_g → N≈1 everywhere. Centering on `ê` (removes the gap) and keeping the
  variance conditional gives each node the precision it deserves.

The source count enters as the count-space Fisher weight; the prior **informs but cannot overrule** real strand
evidence (it dominates only AMBIG nodes, where the strand likelihood is flat — the point).

**Framework decision (supersedes §F's GLM lean): stay with the monotone-spline `var~mean`; do NOT switch to an
NB count-GLM.** The NB appeal was count-space precision, but the count *dispersion* is the wrong scale —
measured φ≈163 (driven by the count *magnitude* at large μ, not the biological CV), implying a spurious N≈0.2.
The density-space `var~mean` (offset-subtracted, conditional) gives the right count-aware precision (N≈24 at
enriched). So we reuse the monotone-spline backend for BOTH the mean `ê` and the variance `σ²_bio(level)` — one
framework, two reads; the identity message is the special case `ê = identity`. Estimator ablation: log-log+recal
and a Poisson mean-GLM tie on recovery (corr 0.52, net within a few %); the relationship is near-log-linear
(R²≈0.95) so the spline reduces to ~linear, with the backend future-proofing real-data non-linearity. Using BOTH
flanking edges (node-pairs / mean) over the louder one is a small net gain (−2.3% → −1.3%) and is folded in.

## F. The net-calibration problem (the real open issue)

Measured on the withheld AMBIG exons, capture ON, **no global / no strand**:

| estimator | net (true 198,429) | per-node corr | MAE |
|---|---|---|---|
| continuous `ê`, log-log + Duan smearing | +63,286 (+32%) | **0.63** | 0.116 |
| continuous `ê`, binned conditional-mean (monotone) | +29,743 (+15%) | 0.45 | 0.147 |
| flat baseline `17.4/T` (the "global") | +3,357 (+1.7%) | 0.32 | 0.198 |
| genome-wide `0.46/T` | −192,745 (−97%) | 0.06 | 0.624 |

The continuous transfer **wins decisively on per-node accuracy** (corr 0.63 vs 0.32, MAE 0.116 vs 0.198) — it
knows which exons are enriched — **but over-calls the aggregate net by 15–32%.** The cause is *not* transfer
bias (median pred/oracle ratio is 0.99 — §B.3); it is a **back-transform / skew artifact**: the error
distribution is right-skewed (median ratio 0.99 but mean 1.89), and a log-space fit back-transformed by Duan
smearing amplifies the tail into the sum.

**The fix is to target the conditional mean `E[ρ_g|z]` directly in count space** rather than fit in log space
and back-transform. The risk review (§L/BL-5) sharpens this to a specific, net-unbiased-by-construction
estimator: a **monotone quasi-binomial GLM with a logit link predicting the fraction `f_g` directly** (response
`g/M`, weight `M`). Two payoffs: predictions land in [0,1] (so the `clip(ρ_g/T)` becomes a no-op — removing a
second, opposite-sign tail distortion, §L/M-6), and the GLM moment condition gives **Σ predicted = Σ observed on
the fit set by construction** — exactly the net calibration the SUM needs, with no magic bin count and no
back-transform. `MonotoneVarMean.fit_offset` cannot do this (it fits a *variance* from squared residuals); we
factor out its monotone-spline core and add a `MonotoneMean` GLM mode. **This is the Phase-0/1 gate.**

## G. Off-capture invariance is emergent (no capture flag, no K)

Off-capture, `z` (boundary density) and `ρ_g` (interior) are both ≈ `ρ_base` everywhere — uniform — so there is
no `z`-variation to fit:

| condition | offset a | slope p | R² |
|---|---|---|---|
| capture ON | +0.94 | +1.18 | 0.912 |
| capture OFF | −0.58 | **−0.006** | **0.000** |

Off-capture the transfer auto-flattens to a **constant** ≈ `ρ_base`. So:

- `ê → const` ⇒ the enrichment layer reduces to a single number = the genome-wide global (the correct
  off-capture answer);
- the **residual layer** `r = ρ_g/const` carries the local spatial variation and imputes it exactly as the
  current identity-mean message does (off-capture imputation reliability is preserved — measured median `N` 3.77
  for the current message).

So **capture detection is emergent from the fit itself** (slope/R² → 0 ⇒ no enrichment), with no capture flag
and no K-selector. Off-capture AMBIG recovery stays strong (corr 0.84–0.89).

**Important caveat (§L/M-9): the flattening is *not* fully automatic.** A monotone fit whose reparameterization
can only produce non-decreasing curves will *rectify* z-noise into a spurious positive slope rather than cancel
it. So the invariance needs an explicit **ê-vs-constant significance gate** (a permutation null on shuffled z);
below the noise floor, early-return to the **literal current `_message`** so capture-off output (and the golden
files) stay bit-identical. Equivalently, nest ê as a monotone *deviation* from `ρ_global` penalized toward zero,
so the no-enrichment limit is ê ≡ ρ_global exactly.

## H. Hard cases (honest)

- **Unstranded data.** No strand ⇒ no single-strand solve ⇒ the *enriched* observers of `ê` are gone (seeds are
  all depleted). The transfer collapses toward the depleted constant and captured AMBIG exons are systematically
  **under-called** — a real identifiability floor, not a bug this layer can cross. The user's seed→single-strand
  imputation bootstrap still gives a *depleted* baseline and partial coverage; the enriched level needs an
  external signal (fragment length: gDNA fragments run longer). Build the stranded path first; leave unstranded
  FL-ready.
- **Short exons** (interior `E → 0`, `T = M/E` is noise-over-zero). Their real gDNA lives in the boundary
  crossings; exclude them from the `ê` fit and **solve their `f_g` from the boundary mass directly**
  (`E[min(frag,L)]` is well-defined), consistent with boundaries being first-class nodes.
- **Nascent-present conditions.** When nascent RNA is non-trivial, the boundary crossings are no longer pure
  gDNA, so `z` is contaminated. Needs the boundary's own deconvolution to supply a gDNA-only `z` (the nascent is
  sparse but not zero) — a dependency to validate.

## I. Validated results (summary)

- Continuous boundary→interior transfer: **R² 0.946**; transfer SS→AMBIG **unbiased** (median ratio 0.99).
- AMBIG recovery, imputation only (no global, no strand): per-node **corr 0.63** vs 0.32 flat — imputation
  resolves per-node enrichment.
- Off-capture: transfer auto-flattens (slope ≈ 0, R² ≈ 0) → reduces to the global; recovery corr 0.84–0.89.
- Open: **net over-call +15–32%** (back-transform/skew) → needs a mean-targeting count-space fit.

## J. Implementation plan (detailed, execution-ready; no C++)

**The model in one box.** One node-pair monotone model. Per-REGION enrichment estimate `ê(z)` (z = the region's
flanking gDNA-clean boundary density) gives the prior **mean** `μ_g = clip(ê(z)·E/M, 0,1)`; the **variance**
`σ²_bio(level)` (conditional, via the existing `var~mean`) gives the precision `N = ρ̂²/(σ²_bio + sampling)`.
`ê` is fit **once, pre-loop**, on single-strand exons. It **replaces** the genome-wide global mean where `z` is
usable, falls back to `ρ_global` otherwise, and (via a significance gate) collapses to `ρ_global` off-capture.
Both `ê` and `σ²_bio` come from one monotone-spline backend; the identity message is the special case `ê=id`.

**Surfaces (all Python):** `variance_model.py` (factor the spline core; add a mean mode), `bp_solver.py`
(`_gdna_seed_estimate`: fit `ê`+`σ²_bio`+`z`; `_global_logprior`: enrichment-aware mean+precision; `node_sweep`:
wire `z`, the gate). `_message` is touched only in the optional refinement phase.

### Phase 1 — monotone-spline backend + `MonotoneMean` (`variance_model.py`)
- **Extract** `_fit_monotone_spline(x, y, w, *, k, degree, robust_iters) -> (knots, coeffs, x_lo, x_hi, edf, lam)`
  from `fit_offset`'s body (the `_knots`/`_bspline_design`/GCV-λ loop/`_fit_monotone` IRLS). `fit_offset` then
  calls it with the de-offset response and heteroskedastic weights — **behaviour-identical** (regression: existing
  calibration unit tests + golden files must stay green).
- **Add** `MonotoneMean.fit(x, y, weight, *, recal_weight) -> MonotoneMean` + `predict(x)`: the monotone-increasing
  conditional mean `E[y|x]` via the same core (response `y`, no Poisson offset), storing a **count-recalibration
  scale** `s = Σ(recal_weight·y)/Σ(recal_weight·ŷ)` so `Σ ŷ·recal_weight = Σ y·recal_weight` (net-unbiased — the
  §F fix, no NB needed). `predict` returns `s·spline(clip(x))`.
- **Gate (unit tests):** monotonicity; flat/no-signal input → ~constant = weighted mean; on a fixture from
  `enrichment_ablation.py`, reproduces `ê` and the −1.3% net.

### Phase 2 — the `z` channel + the pre-loop `ê` & `σ²_bio` fit (`bp_solver._gdna_seed_estimate`)
- **Per-region `z`** (belief-independent, precompute once like `eff_global`/`mass_global` at `:761`): the
  *region-facing* density of its flanking `clean_exon_bnd` boundaries (reuse the `:564` cleanliness test;
  `d_left=Ml/EGl`, `d_right=Mr/EGr`; region `i` reads `d_right[left[i]]` and `d_left[right[i]]`), the two edges
  **averaged**. `NaN` where no clean flank ("z unusable").
- **In `_gdna_seed_estimate`** (already once, pre-loop, on the strand-only init `f_g` = the validated basis), add:
  - fit `ê = MonotoneMean` on **single-strand exon regions** (`is_reg & strand_obs & type==exon & mass>0 &
    isfinite(z)`), response `ρ_g = f_g_init·T`, predictor `log z`, weight `(2κ−1)²·E` (M-13; `/(1+var_g)` if the
    init variance is to hand), `recal_weight = E` (net-unbiased on the gDNA count);
  - fit the **conditional `σ²_bio`** via the existing `fit_offset` on the **ê-residuals**: `mean = ê(z)`,
    `raw = (ρ_g − ê(z))²`, Poisson `offset`, weights — so `σ²_bio` is the biological excess *as a function of
    level* (the result that gave N≈24 enriched / N≈0.2 depleted);
  - the **significance gate**: a permutation null (shuffle `z` among the fit nodes, refit, take the 95th-pct
    fit-R²); if real R² ≤ floor → set `ê ≡ ρ_global` (flat) and mark "no enrichment" (off-capture collapse, M-9).
  - **return** `(rho_global, gdna_vm, var_mean, ehat, z, sigma2_bio_level)` — extend the tuple; `node_sweep`
    unpacks the new fields.
- **Gate:** off-capture → not significant → `ê≡ρ_global`, the z-prior equals the current global (path-identical);
  capture-on → `ê` matches the harness (enriched ρ_g≈17, R²≈0.91); `ê` identical across passes (static basis).

### Phase 3 — the enrichment-aware prior, replacing the global mean (`_global_logprior` / `node_sweep`)
- **Mean:** `μ_g = clip(ê(z)·E/M, 0,1)` for REGION nodes with usable `z`; else `clip(ρ_global·E/M)` (current).
  Boundaries keep `ρ_global`. This **replaces** the global mean — it is not an additive third term (M-7).
- **Precision:** `N = ρ̂_g²/(σ²_bio(ê(z)) + var_mean)` where `z` usable (conditional, M-independent — enriched
  N≈24, depleted small); else the current `n_glob`.
- Delivered through the unchanged `_binom_pseudo` / full-simplex Dirichlet path.
- **Gate (the decisive one):** AMBIG recovery net within a few % + corr ≥ the harness's 0.52, **with the existing
  identity message still active** (confirms the strong ê-prior, N≈24, dominates the weak message N≈1 — if it
  doesn't, Phase 4 becomes required); AMBIG `f_g_var` does NOT collapse (triple-count check); capture-off goldens
  bit-identical (z unusable ⇒ `ρ_global` path).

### Phase 4 — de-enriched residual message (refinement; only if Phase-3 gate shows the message fighting `ê`)
- Change the gDNA message (`_message` / line 844) to carry the de-enriched residual `r_src = ρ_g_src/ê(z_src)`,
  reconstruct `ρ̂_g_dst = ê(z_dst)·μ_r`, precision from `var~mean` on the residual scale; sources gated on the
  fit-basis, directional solved→AMBIG (BL-4, M-10). Flat `ê` ⇒ identical to the current `_message` (golden-safe).

### Phase 5 — edge cases + benchmark
- Short exons (`E→0`): inherit `f_g` from the flanking boundary beliefs (`chain_boundary_side_deconv`), don't solve
  from the degenerate contained `T` (M-14). Nascent-aware `z`: prefer intergenic↔exon crossings, down-weight
  intron↔exon by `(2κ−1)²` (M-11). Run `calibration-benchmark`: capture-on net leak improved, capture-off
  unchanged; the 4-condition signed phantom/siphon gate; monotone convergence; regenerate capture-on goldens with
  documented deltas.

### Risk register coverage
BL-1 (fit once pre-loop) → P2 · BL-2 (static masks) → P2 · BL-3 (clean, region-facing `z`) → P2 ·
BL-5 (mean-targeting, net-unbiased — count-recal, NB trap avoided per §E) → P1 · M-7 (ê replaces global) → P3 ·
M-8 (source-side precision = conditional `σ²_bio`) → P2/P3 · M-9 (significance gate) → P2 ·
M-12 (golden bit-identity off-capture) → P3 fallback · M-13 (strand weighting) → P2 ·
BL-4 / M-10 (de-enriched message, directional sources) → P4 (optional) · M-14 / M-11 (short exons, nascent z) → P5.

## K. Open problems / risks

*(Initial list; the adversarial implementation-risk review augments this — see §L.)*

1. **Net calibration (§F)** — the back-transform over-call; pick the mean-targeting count-space estimator. The
   gating risk for Phase 1.
2. **Strand-solve noise propagates into `ê`** — the realizable strand-solved fit (a=0.94, p=1.18) differs from
   the oracle fit (a=1.13, p=1.08); the response axis `ρ_g = f_g_strand·T` inherits strand-solve error.
3. **The predictor `z` for non-region nodes** — how `z` is defined for boundary nodes themselves, for exon↔exon
   boundaries (not gDNA-clean), and for nodes with no usable flanking boundary.
4. **Convergence / triple-counting** — `ê` refit each pass on beliefs that depend on `ê`; ensure the residual,
   trend, and global do not triple-count the same evidence; monotone convergence.
5. **Off-capture spurious slope** — guard so the fit returns a flat `ê` at R²≈0 rather than fitting noise.
6. **Unstranded under-call floor** and **nascent-contaminated `z`** (§H).

## L. Implementation-risk register (from adversarial review)

A 5-aspect adversarial review (bootstrap/circularity, predictor & message integration, fit calibration,
precision & convergence, edge-cases/invariance) stress-tested the plan against the current source. **Verdict:
imputation-first is sound, but the plan as written contained one unsound step** (per-pass ê refit) and several
integration blockers. All line anchors below were re-verified against `bp_solver.py` (953 lines).

### Soundness correction (supersedes §J step 5)

The validated clean basis — `strand_fg` (strand-only `f_g`, no global/imputation) — exists in production at
exactly **one moment**: the working-copy `f_g` passed to `_gdna_seed_estimate` (`bp_solver.py:786`) **before**
the sweep mutates it in place (`:813`). Refitting ê each pass would read beliefs that already embed ê → the
degenerate `ê=ρ_g/(ρ_g/global)` self-reference → a moving message target and possible non-convergence. **Fit ê
ONCE, pre-loop, on the frozen strand-only init `f_g` restricted to the static `struct_seed | strand_seed` mask**
(exactly how `ρ_global`/`gdna_vm` are already computed once). Per-pass refit is deferred unless a
monotone-damping convergence proof is produced first.

### Blockers (resolve before any code)

- **BL-1 — fit ê once, pre-loop, on the frozen strand-only basis** (not per-pass on swept beliefs). Mirror the
  existing compute-once structure at `_gdna_seed_estimate`. *(strikes §J step 5.)*
- **BL-2 — the fit basis & AMBIG gate are STATIC structural masks**, never belief-derived. Fit basis =
  `struct_seed | strand_seed` (`bp_solver.py:577-578`, AMBIG-free by construction); withhold =
  `~(struct_seed | strand_seed)`; per-point weight = `seed_w` (`:581`). Forbid any `f_g_var`/current-belief
  basis selection.
- **BL-3 — `z(x)` must use ONLY gDNA-clean crossings** (`clean_exon_bnd`, `bp_solver.py:564` — intron/intergenic↔exon
  only; never exon↔exon or AMBIG-flank, which carry RNA and would make `z` endogenous to the node's own RNA).
  Per node-kind: a **boundary** node uses its own larger-side density; a **region** node uses the max over its
  clean-crossing flanks. Source `z` from a belief-independent total-density channel `T=M/E` (geometry only,
  precomputed like `eff_global`/`mass_global` at `:761-762`), **never from `f_g`**. A node with no clean
  crossing → explicit "z unusable" → short-circuit to the global. (Note: captured **edges** are boundary nodes
  at ρ_g≈7.7 — they must get a `z`, or the model fails exactly where it matters.)
- **BL-4 — split the message; precision lives on the residual.** Operate the gDNA channel in residual space:
  emit `r_src = ρ_g(src)/ê(z_src)` as the pairwise `_message` currency (so `N_src` prices the dimensionless,
  source-local residual — honest); reconstruct `ρ_g_dst = ê(z_dst)·μ_r` at the destination *before* forming
  `μ_c`/`_binom_pseudo`; **refit `fit_gdna_varmean` on residual densities `r`, not `ρ_g`** (`:651-652`). Naïvely
  substituting `ρ_src ← ê(z_dst)·r_src` would price the mean and the count on different quantities — breaking
  honest precision. Unit-test: flat ê ⇒ byte-identical to the current `_message`.
- **BL-5 — the estimator is a count-space conditional-MEAN GLM predicting the FRACTION, not log-OLS+Duan.** Log-OLS
  + a single Duan smear is structurally net-biased on a heavy-right-tailed response (median ratio 0.99 but mean
  1.89 → +15–32% over-call; no scalar smear fixes it). And `MonotoneVarMean.fit_offset` is the wrong machinery
  (it fits a *variance* from squared residuals, floored at 0). Instead: **factor out** the shared monotone
  spline core (`_fit_monotone`/`_bspline_design`/GCV, `variance_model.py:215-234`) and **add a monotone
  penalized GLM** — a **quasi-binomial logit-link predicting `f_g` directly** (response `g/M`, weight `M`).
  Predictions land in [0,1] (so the `clip(ρ_g/T)` becomes a no-op — resolves M-6), match the `_binom_pseudo`
  currency, and the GLM moment condition gives **Σ predicted = Σ observed by construction** — the net
  unbiasedness the SUM needs, with no magic bin count and no back-transform.

### Major (resolve; can follow the blocker locks)

- **M-7 — ê REPLACES the global mean under capture; it is NOT a third additive term.** ê, the gDNA message, and
  the global all derive from the same seed set; adding all three triple-counts the evidence → false confidence
  (the very AMBIG over-confidence we're curing). Implement as ONE `w_z`-gated term: scale the global by
  `(1−w_z)`, `w_z` = ê's explained-variance/fit-R² (learned). Off-capture `w_z→0` (global full); under capture
  `w_z→1` (global silent, ê drives). *(This is the "fix message + global together" directive.)*
- **M-8 — ê needs a source-side precision** = its own GLM pointwise prediction SE (`B(BᵀWB+λP)⁻¹Bᵀ`), not the
  scalar `var_mean`. `N_ehat = ê²/(Var_pred(ê)+σ²_between)`, M-independent.
- **M-9 — off-capture flattening is NOT automatic** (correcting §G): the monotone reparameterization
  `β=cumsum([α₀,exp(α₁..)])` can only make non-decreasing curves, so it *rectifies* z-noise into a spurious
  positive slope rather than cancelling it. Need an explicit **ê-vs-constant significance gate** (permutation
  null: shuffle z among seeds, refit, take the 95th-pct slope as the noise floor); below it → early-return to
  the **literal current `_message`** (preserves capture-off golden bit-identity). Equivalently, nest ê as a
  monotone *deviation* from `ρ_global` penalized toward 0, so large-λ gives ê≡ρ_global exactly.
- **M-10 — gate message SOURCES on the fit-basis membership, not `solvable`, and make ê-messages directional
  (solved→AMBIG only).** An unsolved AMBIG node keeps init `f_g=1` and would emit phantom-gDNA residual messages
  on pass 1; the existing guard (`:770-780`) blocks only unsolvable (G1) sources, not unsolved AMBIG (G3). Cap
  any imputed node's emitted `N` by its own `1/f_g_var` (`var_g`, `:816`) so it can't manufacture confidence.
  *(Interacts with the RELAY problem.)*
- **M-11 — `z` is contaminated by nascent** under nascent-present conditions (intron↔exon crossings then carry
  nascent). Prefer **intergenic↔exon** crossings for structural-seed `z`; down-weight intron↔exon by `(2κ−1)²`.
- **M-12 — golden regression**: changing `_message`/`node_densities` re-prices every node (`test_golden_output`
  at rtol=1e-12). The M-9 early-return to the literal current path keeps capture-off goldens bit-identical;
  regenerate capture-on goldens only after net-mass validation, documenting deltas.
- **M-13 — strand-solve noise biases ê** (realizable a=0.94 vs oracle a=1.13). Weight fit points by
  `(2κ−1)²·(1/f_g_var)`; the count-space GLM reduces the offset vs log-OLS; assert ê slope/offset within tol of
  oracle in the harness.
- **M-14 — short exons** (contained `E→0`, `T` degenerate, forces `f_g→0`): inherit `f_g` from the flanking
  boundary beliefs (`chain_boundary_side_deconv` already projects boundary `f_g` onto region sides) rather than
  solving from the degenerate contained `T`.

### Phased build order (retires blockers earliest)

- **Phase 0 — harden the harness (no prod code).** Restrict `z` to clean crossings (+per-node-kind); replace
  log-OLS+Duan with the count-space logit GLM; add the permutation slope-significance gate; weight by
  `(2κ−1)²·(1/var)`. **Decision gate:** R²≈0.946 holds with clean-only z; net Σf_g·M within a few % of oracle
  (post-clip); off-capture ê≈ρ_global; strand-solved slope/offset within tol of oracle. *All four pass → go;
  any fail → the estimator isn't net-unbiased and the redesign blocks here.*
- **Phase 1 — shared monotone-spline backend + `MonotoneMean` GLM** in `variance_model.py` (logit/log link,
  offset, prediction SE). Leave `fit_offset` untouched.
- **Phase 2 — pre-loop ê fit + belief-independent T/z channels** in `_gdna_seed_estimate`/`node_sweep`. Nest ê
  as a deviation from `ρ_global`. Gate: off-capture `max|ê−ρ_global|<tol`; ê identical across passes.
- **Phase 3 — residual-space gDNA message** (`_message`, `fit_gdna_varmean` on `r`, source-side `N_ehat`,
  directional/basis-gated sources). Gate (bit-exact): flat ê ⇒ identical to current `_message`; capture-off
  goldens pass at rtol=1e-12.
- **Phase 4 — global/message handoff** (`(1−w_z)` scaling; ê replaces global mean; clip becomes a no-op). Gate:
  AMBIG `f_g_var` does not collapse (triple-count test); per-node corr ≥ 0.63 capture-on.
- **Phase 5 — edge cases + net validation** (short-exon boundary inheritance, "no usable z" short-circuit,
  nascent-aware z). Final gate: `calibration-benchmark` — capture-on net leak improved, capture-off unchanged.
- **Defer:** per-pass ê refit (needs a convergence proof); unstranded interior enrichment (bounded, not fixable
  — ê learns edge enrichment from clean crossings but has no interior observers; add a regression test that
  unstranded ê never regresses below the current global, document as accepted limit).

### Phase 0 outcome (run 2026-06-22; `scripts/debug/enrichment_phase0.py`)

The productionizable estimator was built and checked against the four gates. Two design refinements were forced
and validated:

1. **`z` is region-FACING, not exon-facing.** Each region reads the clean-boundary crossing density on the side
   facing *it* (an exon reads its own elevated edge; an intron/intergenic region reads its own depleted side).
   Using the exon-facing side for *all* adjacent regions paired depleted introns with the neighbouring exon's
   high edge density → (high z, low ρ_g) outliers that collapsed R² to 0.33.
2. **The transfer is fit on single-strand EXONS only**, not `struct_seed | strand_seed`. Exons carry the
   edge→interior capture gradient (ρ_g ≈ 3× z); introns/intergenic are ~uniform (ρ_g ≈ z) — a *different* law.
   Mixing them halved R². The AMBIG solve-targets are exons, so the transfer is exon-specific. *(Refines BL-2:
   the ê-transfer basis is the single-strand exons; the structural seeds still anchor the global `ρ_base`
   separately.)*

Estimator: a smooth weighted log-log fit (shape ⇒ per-node corr) + a **count-recalibration scale** so
`Σ ê·E = Σ ρ·E` on the fit set (net-unbiased by construction — the cheap stand-in for the production Poisson/
logit GLM; the logit GLM's extra benefit is bounding, since the clip turned out to *help* corr here, not hurt).

| gate | result | status |
|---|---|---|
| G1 clean-z transfer R² ≥ 0.90 | **0.908** | ✅ |
| G2 capture-on net ≤ 5% | **−2.3%** (vs global-only −51%) | ✅ net |
| G2 capture-on per-node corr ≥ 0.60 | **0.53** (vs global-only 0.28) | ⚠️ below bar |
| G3 off-capture collapses to global | sig-gate fires "not sig"; ê flat; net = global | ✅ |
| G4 strand-solved ê within 20% of oracle | median \|Δ\|/ê = **0.04** | ✅ |

**Verdict: the estimator is net-unbiased (the prior-relevant metric, since the per-locus prior is a SUM), the
off-capture invariant holds via the significance gate, and the transfer is realizable (strand≈oracle).** The one
miss is the per-node corr (0.53) versus an arbitrary 0.60 bar — it is ~2× the flat global (0.28) and is the
honest clean-z ceiling (the unclipped corr is *lower*, 0.35, so the clip is not the limiter; the residual layer,
~1 on the uniform-genomic flagship, would add little). Decision point: accept corr 0.53 and proceed to Phase 1,
or first invest in per-node resolution (residual layer / richer multi-boundary `z`).

---

# Part III — Unstranded extension: spliced-derived ρ_g (the strand-free response)

> Status: concept developed + empirically validated 2026-06-22; first cut "raw" (no spliced→strand level
> calibration), approved to build. This lights up the case the Phase-5 strand-contrast gate deferred:
> capture under κ→½ (unstranded). Harness: `scripts/debug/spliced_derived_rho_check.py`.

## A. The enabler — spliced fragments are *motif*-stranded

The enrichment transfer ê(z) needs a per-region gDNA-density response ρ_g. Stranded data supplies it via
strand-deconvolution (ρ_g = f_g·T). Unstranded data (κ→½) cannot — `f_g_init` collapses to the ~0.5 prior
midpoint, and the Phase-5 gate therefore *defers* to the genome-wide ρ_global.

But there is a signal that survives: **spliced fragments are stranded by the GT–AG splice MOTIF, not by the
read chemistry.** A junction-crossing fragment's transcript strand is read off the motif, so the spliced channel
carries strand even when the reads are fully unstranded. And "single-strand region/boundary" is an *annotation*
property (the signature, from the GTF — always available). So unstranded data loses only the *response* ρ_g;
the fit basis (single-strand exon nodes) and the RNA signal at their junctions are intact.

## B. The arithmetic — reconstruct ρ_g from the spliced (pure-mature) signal

Spliced fragments are pure mature mRNA. For a single-strand exon region R with a flanking clean (intron↔exon)
boundary B on R's motif-strand (reusing the existing `_message` spliced-projection geometry):

```
ρ_mature   = M_spliced(B, R-facing, R's motif-strand) / E_rna_crossing(B)   # mature per-bp density at junction
M_mature(R)= ρ_mature · E_rna_contained(R)                                  # project to R's intra-exon mature
ρ_g(R)     = clip( M_unspliced(R) − M_mature(R), 0 ) / E_gdna_contained(R)  # residual unspliced → gDNA(+nascent)
```

The spliced mass at the boundary fixes the mature density (a *lower bound* on total RNA — nascent could add);
projecting it onto R via the eff-lengths accounts for the mature share of R's unspliced mass; the residual is
split between gDNA and nascent and — absent strand — attributed to gDNA (the nascent-sparse assumption the seeds
already make).

**Two honest caveats.** (1) It is an *upper bound* on gDNA: the residual is gDNA **+ nascent** (nascent is
unspliced RNA with no splice signature). Strand-deconvolution separates gDNA from *all* RNA; the spliced method
cannot — so in the nrna-present + unstranded corner it over-attributes to gDNA (an intrinsic unstranded limit).
(2) **Single-exon transcripts** have no junction → no spliced → unrecoverable (the same gDNA-identical floor that
limits the EM; FL is the only escape, deferred).

## C. Validation (gdna 3:1, capture-on)

| | ρ_g_spliced vs ORACLE | ρ_g_strand vs ORACLE | ê-fit → AMBIG recovery |
|---|---|---|---|
| **ss0.99** (stranded) | corr 0.874 | 0.982 (strand wins) | net +11.3%, corr 0.63 |
| **ss0.5** (unstranded) | corr **0.895** | 0.768 (the 0.5·T artifact) | net **+11.0%**, corr **0.52** |

At ss0.5 the spliced method *beats* the read-strand method and matches its own ss0.99 quality — strand-immune,
as predicted. There is a level bias (enriched median 21.8 vs oracle 17.4 — the projection under-counts mature →
over-attributes gDNA), surfacing as the +11% net over-call; the first cut accepts it (a spliced→strand level
calibration is the obvious follow-up).

## D. Integration — a reliability-weighted blend (replaces the hard gate)

Only the ê-fit *response* changes; the ê machinery, significance gate, and conditional precision are unchanged.
The response becomes a blend of the two estimators, weighted by their reliability (κ is global ⇒ scalar weights):

```
w_str = (2κ−1)²            w_spl = 1 − w_str            avail = (spliced present at R's flanking boundary)
denom = w_str + w_spl·avail
ρ_g(R) = ( w_str·ρ_strand(R) + w_spl·avail·ρ_spliced(R) ) / denom
fit weight w_fit(R) = denom · E(R)      recal weight = E(R)
```

- **κ→1 (stranded):** ρ_g ≈ ρ_strand (≈4% spliced bleed — negligible; the shipped stranded result is preserved).
- **κ→½ (unstranded):** ρ_g = ρ_spliced on avail nodes — the deferred case lights up.
- **In between:** smooth hand-off — "lean harder on the spliced solution as the strand weakens."
- A node enters the fit basis if w_str is meaningful **or** spliced is available; if **neither** (unstranded
  *and* single-exon / no spliced) it is excluded → falls to the global.

This **dissolves the Phase-5 `(2κ−1)²≥0.25` strand-contrast gate**: that gate was the "defer unstranded"
placeholder; the spliced term is the unstranded path, so the gate (and its flagged constant) are removed.

Crucially, the unstranded significance test now fires on a *genuine* enrichment signal (real spliced-derived
ρ_g, corr 0.89), not the `0.5·T` artifact — so the Phase-5 phantom-gDNA regression does not return.

## E. Validation gates (before commit)

1. Production `fit_enrichment_transfer`: κ=0.99 → blend≈strand (net ≈ −1%); κ=0.5 → blend=spliced (net ≈ +11%).
2. Calibration units 184 pass / 1 (`test_ambig`); goldens green (tiny scenarios degenerate → flat → byte-identical).
3. Full 16-condition re-quant: **unstranded capture-on leak drops** (deferred case recovers), **no-gDNA
   unstranded FP stays controlled** (no phantom return), **stranded + off-capture unchanged**.

## F. Accepted limitations (first cut)

- +11% level bias (no spliced→strand calibration this pass — the obvious refinement).
- gDNA-vs-(gDNA+nascent) mismatch in nrna-present + unstranded — intrinsic unstranded limit.
- Single-exon / no-spliced + unstranded — unrecoverable (FL-deferred).
