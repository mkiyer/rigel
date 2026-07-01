# Phase 2 — the gDNA-density mixture prior: model-fit design

**Status:** design proposal, 2026-06-30. Consumes the shipped Phase-1 prior-free single-strand solver
(`origin/main d6f79c53`). Parent strategy: `PRODUCTION_DESIGN_gdna_mixture_prior.md` (§3 = this doc).
Reuses machinery in `bp_solver.py` (`node_densities`, `_gdna_seed_estimate`, `_floor_estimate`,
`fit_enrichment_transfer`'s `rho_spliced`), `simplex_logodds.py` (the grid solve), `variance_model.py`
(GCV/CV). Nothing here needs C++.

---

## 0. What this fits, in one paragraph

The calibrator's deliverable is an accurate per-locus gDNA prior. The hard case is the **AMBIG** node
(overlapping opposite-strand transcripts): its strand likelihood is flat in `f_g`, so its own data cannot
say how much of its mass is gDNA. Phase 2 learns, from the *trustworthy* single-strand (SS) nodes the
shipped solver already resolves, a **nonparametric distribution `P(log ρ_g)`** of per-node gDNA densities
— unimodal off-capture, bimodal under capture (depleted off-target ⊕ enriched on-target), multimodal under
CNV, **from one code path with one knob (bandwidth) and no capture/mode detector**. Phase 3 then prices each
AMBIG node by adding `log P(log ρ_g)` as its per-node prior term in the existing grid solve — a **drop-in**
for the Gaussian `_global_logprior` already computes. The prior fills exactly the strand tilt's null space,
its strength self-scaling by the tilt's Fisher information `∝ N(2κ−1)²` — no tuned weight.

---

## 1. The training substrate (the teacher set) — exact definition

**Source:** the Phase-1 `node_sweep` output (`NodeBelief` per chain node) → `node_densities` → per-node,
per-face gDNA density `ρ_g = f_g · M / E_gdna`. We take the **log** (densities span 5–6 decades; the
mixture lives in log space).

**Which nodes teach (confidence, not strand-class — `PRODUCTION_DESIGN §2`):**
- **SS region nodes** — `strand_class ∈ {POS, NEG}`, `mass_unspliced > 0`, solved as G2. Their `ρ_g` is
  strand-pinned (high κ) → the enriched + depleted teacher across the full density range.
- **Structural region nodes** — intergenic (pure gDNA, locked) and single-strand introns. These anchor the
  **depleted** end at *any* κ (they need no strand). Intergenic is the purest depleted class.
- **Clean-exon BOUNDARY nodes** (intron/intergenic↔exon crossing, exon-facing side) — the capture-observable
  *enriched-edge* gDNA signal; the `_gdna_seed_estimate` "structural" seeds already isolate these. Include
  them (open question Q3) — they extend the enriched mode's support where expressed SS exons are sparse.
- **EXCLUDE:** AMBIG nodes (the students — non-circular teacher/target split, `PRODUCTION_DESIGN §3.5`),
  zero-mass nodes, and TSS/TES black-hole boundaries.

**Per-sample weight `w_i`** (a KDE sample is one *observation that a gDNA-density level exists*, so we do
NOT exposure-pool — that would let a few deep nodes dominate the mode locations; `§3.5`):
```
w_i = strand_discriminability_i · reliability_i
    = (2κ−1)²                    ·  1/(Var(log f_g)_i + 1/(gcount_i+1))
```
- `(2κ−1)²` fades a κ→½ node's *strand-derived* density to 0 (unstranded → only structural teachers survive,
  which is correct — an unstranded exon's ρ_g is not trustworthy).
- `Var(log f_g)_i` = the solve's own posterior log-fraction variance (`gdna_frac_var`, already returned);
  `1/(gcount+1)` = the Poisson log-density floor. Together = the node's density precision → a noisy node
  contributes a *broad, light* bump, a confident node a *sharp, heavy* one. This is the CLOSE/EBCF idea
  (precision-weighted teachers) applied to the density estimate, not a regression.
- Structural/intergenic teachers get `(2κ−1)²→` replaced by weight 1 (locked f_g=1, always count), matching
  `_gdna_seed_estimate`'s treatment.

> **Design note (no magic numbers):** every term above already exists in the solve output. `w_i` introduces
> **no new constant**. If we later find equal-weighting gives cleaner modes, that's a one-line switch — flag
> Q2.

---

## 2. Phase 2 — the nonparametric fit

### 2.1 The estimator
Weighted Gaussian KDE in `x = log ρ_g`:
```
P̂(x) = Σ_i w_i · φ((x − x_i)/h) / (h · Σ_i w_i),     φ = standard normal kernel
```
- **Log space** (not linear ρ_g): densities are multiplicative/heteroscedastic; log makes the modes
  roughly equal-width and the kernel bandwidth meaningful across the range.
- **Gaussian kernel:** smooth, infinitely differentiable (the grid solve wants a smooth prior term; a
  boxcar/Epanechnikov would put kinks in ψ). One kernel, one knob.

### 2.2 The one knob — bandwidth `h` (needs a robust estimator + a plotting framework, `PRODUCTION_DESIGN §3.3`)
Ship **three** estimators behind a config selector, plus a plotting harness; pick production empirically
(on real data when it lands — do NOT hard-code a magic `h`):
1. **Silverman's rule of thumb** — `h = 0.9 · min(σ̂, IQR/1.349) · n_eff^{−1/5}`, `n_eff = (Σw)²/Σw²`
   (Kish effective sample size for weighting). The robust, always-available default.
2. **Likelihood cross-validation (LSCV/LCV)** — leave-one-out KDE log-likelihood maximized over `h`. This is
   the same *cross-validation principle* the `variance_model` GCV machinery uses (§ reuse below), but for a
   kernel width, not a P-spline λ — so it is a small new routine, not a call into `MonotoneVarMean`.
3. **A bimodality-preserving floor** — clamp `h` from above at the value where the depleted/enriched modes
   merge (detected as the density becoming unimodal between the intergenic floor and the enriched teachers).
   This encodes the *only* structural fact we trust (there ARE ≥2 modes under capture) without counting them.

> **Reuse `variance_model` GCV where possible (per the directive):** the P-spline GCV/`_select_lambda` in
> `variance_model.py` is for spline smoothing, not KDE bandwidth — they are different objectives. We reuse
> the *CV pattern* (score a candidate by held-out fit) and the plotting/inspection discipline, and keep the
> `MonotoneVarMean` machinery for what it already does (the message/floor spreads). Do not force a P-spline
> where a KDE is the right tool.

### 2.3 Representation for the solver (efficiency)
Pre-evaluate `logP̂(x)` on a **fine 1-D grid** of `x = log ρ_g` (e.g. `K_ρ = 1024` points spanning
`[min x_i − pad, max x_i + pad]`), store `(x_grid, logP_grid)`. Phase 3 does **linear interpolation** of
`logP_grid` at each AMBIG node's per-grid-point `log ρ_g(λ_k)`. Cost: pre-eval `O(K_ρ · n_teacher)` once;
per-node `O(K_λ)` interp. This keeps genome scale tractable (n_teacher can be ~10⁶; on-demand
`logsumexp` per grid point would be `O(K_λ · n_ambig · n_teacher)` — infeasible). Outside the grid range,
`logP̂` extrapolates as the boundary value (a flat tail = "no prior information beyond the observed range",
the honest default; the tilt/messages/residual then decide).

### 2.4 The depleted floor anchor
`_floor_estimate` already returns `ρ_floor` (intergenic+intron depleted minimum) + its spread. The KDE's
depleted mode *should* land at `ρ_floor` because the same intergenic/intron nodes are in the substrate.
**Anchoring choice (Q6):** optionally seed the KDE with a small pseudo-mass virtual sample at `log ρ_floor`
(weight = the floor's pooled effective count) so the depleted mode is pinned even when few depleted SS
nodes are present. Default: rely on the substrate; expose the virtual-sample option.

### 2.5 Deliverable: the `GdnaDensityPrior` object
```
GdnaDensityPrior:
    x_grid: (K_ρ,)          # log ρ_g grid
    logP_grid: (K_ρ,)       # log P̂(log ρ_g), normalized
    bandwidth: float        # the selected h
    n_eff: float            # effective teacher count
    def logpdf(log_rho: array) -> array     # linear-interp lookup (the Phase-3 hook)
    # diagnostics for the plotting framework:
    teacher_x, teacher_w    # the raw samples + weights (rug plot)
    modes: [(x, height)]    # detected local maxima (REPORTING ONLY — never used to classify)
```
`modes` is **diagnostic only** (for QC / the plots) — the solve never uses a hard mode assignment
(`PRODUCTION_DESIGN §3.2`, §6 ruled-out hard clustering).

---

## 3. Phase 3 — consuming the prior (the drop-in)

### 3.1 The one-line swap — applied to ALL nodes (self-scaling; Q4 resolved 2026-06-30)
`_global_logprior` currently adds, per grid point `k`, a **Gaussian** on `log ρ_g`:
```
glob_k = −½ · n_node · (log σ(λ_k) − target)² ,   target = log(clip(ρ_global·E/M, ε, 1))
       = −½ · n_node · (log ρ_g(λ_k) − log ρ_global)²          # a single-mode prior at ρ_global
```
Phase 3 replaces this **for every solvable node** with the mixture — ONE code path, no AMBIG/SS fork:
```
glob_k = prior.logpdf( log ρ_g(λ_k) ) ,   log ρ_g(λ_k) = log( σ(λ_k) · M / E_gdna )
```
Same slot in `_local_loglik_logodds` / `_solve_ambig_logodds` (the `global_logprior` `(m,K)` argument).

**Why applying to all nodes is safe (the non-circularity, made rigorous).** The teacher/student split is
now *emergent from Bayes*, not enforced by a mask:
- **The prior is weak exactly where it was trained.** On an SS node the strand-tilt Fisher information in
  the `f_g` direction is `∝ N(2κ−1)²` (large), so the node barely listens to the *population* prior it
  helped build — the self-feedback is second-order in the tilt-vs-prior precision ratio. On an AMBIG node
  the tilt is flat → the prior dominates, and AMBIG nodes **contributed nothing** to the substrate (§1
  excludes them), so an AMBIG node is still taught purely by SS/structural nodes.
- **A node is 1/n_eff of a population object.** Removing any one teacher barely moves `P̂` (O(1/n_eff)),
  so an SS node's pull toward *its own* contributed bump is negligible at genome scale (and ~1% even on a
  ~100-node toy). If a scenario ever shows measurable self-reinforcement, the strict fix is a
  **leave-one-out** `logpdf` (subtract the node's own kernel before evaluating) — available, but O(n²) and
  expected unnecessary. Do not build it until a plot demands it.

Net: the *substrate* (who teaches, §1) stays SS + structural — AMBIG excluded, non-circular. The *apply*
(who is priced, here) is uniform — Bayes' Fisher-info weighting makes it a near-no-op on the teachers and
the answer on the students. One term, one path.

### 3.2 The strand-free per-node observation (localizes isolated AMBIG loci)
For each AMBIG node add its **spliced-residual** — a per-node, strand-free, gDNA-specific *observation* of
`ρ_g` (the `rho_spliced` already computed in `fit_enrichment_transfer:796`, generalized to AMBIG nodes and
their flanking clean boundaries):
```
ρ_mature = M_spliced(motif-stranded, flanking clean boundary) / E_spliced_half_triangle
ρ_resid  = clip(M_unspliced − ρ_mature · E_rna_contained, 0) / E_gdna_contained
```
Enter it as a Gaussian on `log ρ_g` at the residual's own count precision — a Bayesian shrinkage of the
observation toward the nearest mixture mode. This is what lets a **whole-AMBIG locus with no SS neighbour**
still localize (Q9). It is an *observation updated against the prior*, NOT a regression predictor
(`PRODUCTION_DESIGN §6` — the ê(z)-on-gDNA regression is dead).

### 3.3 Self-scaling (why no weight knob)
The AMBIG posterior at grid point `k` fuses: strand-tilt BB (flat in `f_g` when balanced) + gDNA messages
from solved SS neighbours (disagreement precision) + `ρ_resid` observation + `prior.logpdf`. The tilt's
Fisher info in the `f_g` direction is `∝ N(2κ−1)²` → 0 at a balanced ridge / κ→½ / off-capture, so the
prior dominates *exactly* where the node is uninformative and recedes where the tilt pins `f_g`. The
bandwidth `h` sets the prior's sharpness (a tight KDE = a strong pull to the modes; a broad KDE = weak) —
so `h` is the *implicit* prior-strength control, one more reason it is the single knob to get right.

### 3.4 Multimodal posterior → the per-locus EM Dirichlet scalar (open, §Q7)
An AMBIG node with no resolving evidence yields an honest bimodal posterior. Options: (a) posterior **mean**
of `f_g` (simple, defensible; slightly hedged); (b) posterior **median** (the current SS convention); (c)
propagate the full distribution to the EM. Recommend **(a) mean** for the point estimate feeding
`assemble_priors`, and record the posterior spread as the prior's confidence — decide on real data.

---

## 4. Orchestration — a 3-stage single pass (no iteration)

Restructure `node_sweep` into three explicit stages (the human transcriptome has enough SS nodes to train
once; verify no iteration needed — Q8):

```
STAGE 1  (Phase 1, shipped)  Solve SS + structural nodes. AMBIG nodes held at their {0,0,1} init and
                             EMIT no messages (they carry no information yet). SS FB sweep as today.
                                    │  → SS/structural final beliefs (+ their gDNA messages)
STAGE 2  (Phase 2, NEW)      Build the training substrate (§1) from the Stage-1 densities; fit
                             GdnaDensityPrior (§2). Anchor the depleted mode at ρ_floor.
                                    │  → prior.logpdf
STAGE 3  (Phase 3, NEW)      Solve with global_logprior = prior.logpdf (§3.1) + the ρ_resid observation
                             (§3.2), applied UNIFORMLY to all nodes. AMBIG nodes fuse the (flat) tilt +
                             the mixture + ρ_resid + gDNA messages from now-solved SS neighbours; SS
                             nodes re-solve as a near-no-op (self-scaling) and MAY be frozen from Stage 1
                             for efficiency (equivalent — the prior is weak where the strand pins f_g).
                                    │  → final beliefs (all nodes)
                             chain_region_deconv / chain_boundary_side_deconv  (unchanged)
```
Anchoring/non-circularity: the prior is fit **once** on Stage-1 (belief-independent of Stage-3), consistent
with the shipped anchored-global architecture. Message-ordering detail: Stage-3 AMBIG messages read the
Stage-1 *final* SS beliefs (an AMBIG node between two SS genes gets both neighbours' solved gDNA density).

---

## 5. Implementation plan (phased; each phase self-contained + tested)

- **P2.0 — substrate extractor.** New `calibration/gdna_density_prior.py::build_training_substrate(chain,
  belief, geometry, region_arrays, boundary_substrate, kappa) -> (log_rho, weight, node_kind)`. Pure
  function of the Stage-1 output. Unit test: on a controlled toy (known f_g per node), the returned
  `log_rho` matches `node_densities`; weights fade to 0 at κ=½ for strand teachers.
- **P2.1 — the KDE + estimators.** `GdnaDensityPrior.fit(log_rho, weight, *, bandwidth='silverman'|'lscv'|
  float, floor_anchor=None)`. Unit tests: recovers a known 2-Gaussian mixture's mode locations; Silverman
  vs LSCV agree within tolerance on a clean bimodal sample; `logpdf` interp matches direct KDE eval.
- **P2.2 — the plotting framework.** `scripts/debug/plot_gdna_prior.py`: run `toy_prod` (or a benchmark
  condition), extract the substrate, overlay P̂(log ρ_g) for each bandwidth estimator + a rug of the
  teacher samples + the true per-node gDNA densities (oracle) + `ρ_floor`. This is how we *see* the fit and
  choose the estimator. Deliverable before any production bandwidth decision.
- **P2.3 — Stage-1/Stage-3 split in `node_sweep`.** Refactor so AMBIG nodes are excluded from Stage 1 and
  solved in Stage 3; verify SS-only results are byte-identical to today (AMBIG were already weak). This is
  the risky refactor — do it behind the split with the mixture prior *off* first (Stage 3 = today's weak
  global) to prove no SS regression, THEN turn on the mixture.
- **P2.4 — Phase-3 prior hook.** Wire `prior.logpdf` into `_solve_ambig_logodds`'s `global_logprior` for
  AMBIG; add the `ρ_resid` observation term. Unit test: `test_gdna_sweep_factor1_uniform` flips BACK to
  recovering ρ=0.5 (restore `atol=0.05`) — the Phase-2 prior closes the documented Phase-1 gap.
- **P2.5 — config surface.** `CalibrationConfig`: `gdna_prior_bandwidth: str|float = 'silverman'`,
  `gdna_prior_floor_anchor: bool`, and the density-grid resolution. No magic bandwidth.
- **P2.6 — validation + goldens.** AMBIG toy_prod scenarios (below) → the net-flow benchmark → regen
  goldens → full pytest.

---

## 6. Validation (per `PRODUCTION_DESIGN §9`)

- Grow `toy_prod` with AMBIG scenarios that today under-call: **partial-overlap** and **fully-encompassing
  opposite-strand** gene pairs, at the full gDNA:RNA × κ × capture ladder, plus a ~100-SS-transcript
  background so the substrate spans both modes (the enriched-mode support, O3).
- Metric: converged per-node `|f_g − oracle|` on the AMBIG nodes (partial correlation controlling for true
  RNA to avoid the co-location trap) — must fall without regressing the SS nodes.
- Then the 16-condition net-flow benchmark (skill `calibration-benchmark`), goldens, full pytest.
- Real data when it lands (synthetic can't surface prep/GC/mappability/κ subtleties — the bandwidth
  decision especially waits for real data).

---

## 7. Open questions / things to sort out

1. **Bandwidth estimator (the one knob).** Silverman default now; LSCV + the bimodality floor as options;
   *decide production on real data via the plotting framework.* Is a fixed rule stable across the
   gDNA:RNA × κ × capture ladder, or must `h` be data-adaptive per library?
2. **Substrate weighting.** Reliability-weighted (proposed, §1) vs equal-weight vs light exposure. Affects
   mode *locations* and *heights*. Inspect on the plots.
3. **Boundary nodes in the substrate.** Include clean-exon boundary crossings (proposed — they extend the
   enriched-mode support) vs regions-only (simpler, avoids the crossing-geometry density subtlety)?
4. **Prior scope — RESOLVED (2026-06-30): apply to ALL nodes (self-scaling).** One code path; the
   substrate still excludes AMBIG (§1) so students stay non-circular, and the tilt Fisher-info weighting
   makes the prior a near-no-op on the teachers (§3.1). Leave-one-out `logpdf` held in reserve if a plot
   ever shows SS self-reinforcement.
5. **Density representation.** Pre-evaluated grid + interp (proposed, fast) vs on-demand KDE (exact, slow).
   Confirm the grid resolution `K_ρ` is enough not to blur the modes.
6. **Floor anchoring.** Let the KDE find the depleted mode vs seed a `ρ_floor` virtual sample. Matters when
   depleted SS nodes are sparse.
7. **Multimodal posterior → EM scalar.** Posterior mean (proposed) vs median vs propagate the distribution.
8. **Flow.** Single 3-stage pass (proposed) vs any iteration. And: is the Stage-1/Stage-3 split worth the
   `node_sweep` refactor, or can Stage 3 be a cheap AMBIG-only re-solve *after* the existing single sweep?
9. **Isolated whole-AMBIG loci.** Rely on `ρ_resid` + the marginal prior (proposed). Sufficient, or does a
   probe-overlap covariate fallback become necessary? (Prefer not to add one — decide empirically.)
10. **Teacher spectrum coverage (O3).** If a library's SS exons are mostly unexpressed, the enriched mode is
    under-sampled and enriched AMBIG nodes extrapolate. Measure the teacher vs AMBIG `ρ_g` histograms; the
    boundary-crossing teachers (Q3) are the main mitigation.
11. **Nascent contamination (O1).** `ρ_resid = gDNA + nascent`; the depleted-mode/floor assumes nascent
    sparse. Fails in `nrna_rnd` conditions. The prior is a population object so a few nascent-contaminated
    teachers only broaden a mode; but confirm on nascent scenarios.
12. **CNV / multimodality (future).** The nonparametric prior captures extra modes for free; the
    disagreement-aware, breakpoint-respecting messages localize a node to its segment. Only mine the
    copy-number-segmentation literature if intrinsic CNV becomes a target.

---

## 8b. P2.2 empirical findings (from the plotting framework, 2026-06-30)

Built `calibration/gdna_density_prior.py` (substrate + KDE) + `scripts/debug/plot_gdna_prior.py`, run on
`toy_prod` across the κ × capture ladder. What the plots show:

- **The core hypothesis holds.** Capture-on is cleanly **bimodal** — a depleted mode (intergenic/intron)
  and an enriched mode (captured exons, with a boundary-crossing shoulder between). Capture-off is
  **unimodal** (uniform gDNA). Same code, no detector. The fitted modes track the oracle density clusters.
- **Bandwidth is delicate (confirms §7.1).** LSCV under-smooths — it splits the depleted mode into spurious
  peaks that are just the `1/E` floor of different-sized empty intergenic regions. Silverman is more robust
  but its `min(std, IQR)` scale is fooled by a peaked-bulk-plus-geometry-tail shape. **Added a principled
  noise floor** (bandwidth ≥ the weighted-median per-node Poisson log-density std, `log_rho_std`) — no magic
  number; it stops the sampling-noise fracturing.
- **Residual off-capture sub-modes are a SYSTEMATIC GEOMETRY BIAS, not sampling noise:** small exons
  (tiny `E_gdna`) and boundary crossings read a *lower* apparent gDNA density than large regions even under
  uniform gDNA (few-fragment + the depleted-neighbour drag). The noise floor doesn't remove these (they're
  bias, not variance). → **Open decision D1** (below).
- **Depleted mode sits at `1/E` under strong capture** — intergenic is nearly empty, so its "density" is the
  min-observable (E-dependent). Directionally right ("~0 off-target gDNA") but the location is artificial.
  → the `ρ_floor` virtual-sample anchor (§2.4) is the intended fix; validate it lifts the depleted mode to a
  meaningful, E-independent level.
- **Unstranded (κ→½) collapses the enriched teachers.** The strand-derived exon density carries `(2κ−1)²`,
  so at κ=½ the exon-interior teachers get weight ≈0 and the enriched mode is left to the boundary crossings
  alone (biased low). → **Open decision D2:** add the strand-free **`ρ_resid`-blended** exon teacher density
  (`(2κ−1)²·strand + (1−(2κ−1)²)·ρ_resid`, the exact blend `fit_enrichment_transfer` already computes) so
  the enriched mode survives unstranded. This is the O3/O4 mitigation and the "optimize unstranded too"
  directive; it depends on the mature-eff-len quality (O2).

**Decisions to make before P2.4 (wiring the solve):**
- **D1 — geometry-bias sub-modes:** live with them (a broad-enough bandwidth washes them into the correct
  unimodal off-capture prior, and they are all "≈uniform gDNA" anyway), OR de-bias the small-exon/boundary
  teacher density, OR down-weight/exclude sub-`E` teachers. Prefer: bandwidth + weighting handle it; confirm
  on the plots.
- **D2 — unstranded enriched teacher: DONE (2026-06-30).** Added the strand ⊕ `ρ_resid` blend
  (`_exon_spliced_residual`): the exon teacher density is `(2κ−1)²·strand + (1−(2κ−1)²)·ρ_resid` on exons
  with flanking spliced. At κ=½ the enriched exon teachers now survive (exon Σw 0→6.3; a fitted enriched
  mode reappears at ≈+0.5). **Two follow-ups surfaced (for the D1 dissection turn):** (i) the `ρ_resid`
  teachers get LOW weight because `precision` uses the strand solve's `Var(log f_g)` (unreliable when
  unstranded) — it should use the SPLICED count; so the unstranded enriched mode is still boundary-dominated
  (biased low vs the oracle's ≈−0.5 spike). (ii) the `ρ_resid` level itself may be biased by the mature
  eff-len (O2) — the +0.5 unstranded mode sits above the −0.6 stranded exon mode.
- **D3 — bandwidth production rule:** Silverman + noise-floor is the robust default; keep LSCV + the plots
  for inspection; the final rule waits for real data (§7.1).

**D1 (geometry-bias sub-modes) — DEFERRED to a dedicated dissection turn** (user, 2026-06-30): plot +
understand the small-exon/boundary density bias, the `ρ_resid` weighting (D2-i), the `ρ_resid` level bias
(D2-ii / O2), and the boundary-dominance, before deciding the fix. The current default (bandwidth +
precision weighting) stands until then.

## 8. Why this is the right design (the invariants it respects)

- **Count-zero-info** (`CALIBRATION_ARCHITECTURE.md`): the mixture is over gDNA *densities*, and enters as a
  prior on `log ρ_g` — the count still only sets precision (the teacher weights + the tilt Fisher info),
  never the composition directly.
- **No magic numbers:** one knob (bandwidth), chosen by a data estimator with a plotting framework; every
  weight/precision term already exists in the solve.
- **No detector, no hard clustering:** one KDE code path yields uni/bi/multi-modal priors; the per-node
  solve marginalizes softly; nothing ever hard-assigns a mode (`PRODUCTION_DESIGN §6`).
- **Non-circular:** teacher = SS/structural nodes, student = AMBIG; the prior is fit once, belief-
  independent of the nodes it prices — same anchored architecture as the shipped global.
- **Drop-in:** Phase 3 is one term swap in the existing grid solve; no new solver machinery
  (`PRODUCTION_DESIGN §4`).
