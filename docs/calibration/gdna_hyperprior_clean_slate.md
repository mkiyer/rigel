# The gDNA hyperprior that respects enrichment & depletion — a clean-slate design study

**Status:** strategic design study (2026-07-20). A zoom-out to think about the Role-B gDNA composition prior as
a whole before honing in. Consolidates: the released **v0.7.1 KDE** (`gdna_density_prior.py` @ tag `v0.7.1`),
this branch's **NPMLE** refactor (`npmle.py`), the crush dissection (`npmle_crush_dissection.md`), the restore
plan (`gdna_kde_restore_plan.md`), the zero-handling study (`archive/gdna_prior_zero_handling.md`), and the
`kde_vs_npmle_enriched_mode` finding. Companion figures live in the session scratchpad (EM-vs-KDE overlay, the
node-1055 trace).

---

## 0. STARTER PROMPT — orient a new session

> You are picking up the **gDNA hyperprior (Role B)** work on branch `calib-ambig-init-wip`. Pass-0 is solid and
> honest (its error is ~80% identifiability, not confident-wrong — `npmle_crush_dissection.md` §0). The **Phase-2
> refit** consults a fitted population density `P(log ρ_g)` as the composition arm of ψ, projected per node by the
> **δ-pin** `f_g ← argmax P(log(f_g·M/E))`. This coupling means: *a prior over gDNA **density**, applied at a node
> of total density `ρ_tot`, dictates `f_g ≈ ρ_g^{prior-mode}/ρ_tot`.* Hybrid capture makes gDNA density
> **non-uniform** (enriched at exons); a prior that ignores that **crushes** enriched, unstranded nodes to
> `f_g≈0` (traced number-by-number on node 1055 in `npmle_crush_dissection.md`).
>
> **The released v0.7.1 KDE solved this**; this branch's EM-NPMLE reintroduced the crush. Read §1–§4 here, then
> the restore plan. The immediate, concrete work is **§8 experiments**: with the additive kernel ON
> (`config.gdna_prior_additive`), (a) calibrate the **floor strength** (the zero-gDNA FP vs captured-recovery
> frontier), (b) settle **precision handling** — unit weight vs the "one-count-spread-by-precision" kernel (§6),
> (c) settle **bandwidth** (§7). Validate on the 32 cached `ambig_dense_10mb` scenarios, refit-vs-oracle EMD by
> node type, with node 1055 and the `none_*` scenarios as the two poles. **No new magic numbers without
> discussion; no clamps/cliffs — smooth, honest Bayesian behaviour.** Tools: `scripts/debug/refit_gdna_movie.py`,
> `hyperprior_suite_vs_oracle.py`, `pass0_error_diagnostic.py`, and this session's `kde_vs_em_compare.py`,
> `trace_crush_node.py`, `dissect_neg4.py`.

---

## 1. Why the hyperprior must respect enrichment & depletion

The composition prior enters the solve through a **δ-pin** (`DensityNPMLE.logprior` / `GdnaDensityPrior.
logpdf_kernel`): for a node with total density `ρ_tot = M/E`,

```
    logprior(f_g)  =  logP( log ρ_g )  where  ρ_g = f_g · ρ_tot .
```

So a prior over gDNA **density** is, at each node, a prior over `f_g = ρ_g / ρ_tot`. Two consequences define the
whole problem:

1. **The prior's mode dictates f_g inversely with enrichment.** If `P` says "gDNA density ≈ `ρ_g^mode`", the node
   is pushed to `f_g ≈ ρ_g^mode / ρ_tot`. A **captured (high-`ρ_tot`) node** is therefore pushed **down** unless
   `P` carries mass at the *enriched* gDNA level. gDNA density is **not uniform** — capture enriches it at exons
   just as it enriches RNA — so a prior fit as one pooled, depleted-dominated distribution is *wrong at exactly
   the enriched nodes*, and (where the strand can't resist — unstranded) it crushes them. Node 1055: true f_g
   0.902 (21,860 gDNA fragments) → refit 0.0001.
2. **The prior is the only information on unstranded nodes.** Stranded nodes are pinned by the Beta-Binomial
   strand likelihood and are safe. Unstranded (`ss_0.50`) nodes have a flat strand channel, so ψ ≈ prior, and the
   prior decides f_g almost alone. Every failure mode below is an unstranded phenomenon.

**Requirement.** The hyperprior must carry mass at BOTH the depleted background level (so true-zero-gDNA nodes go
to `f_g≈0`) AND the capture-enriched level (so real enriched-gDNA nodes are not crushed) — a **bimodal, enrichment
-respecting** density — and it must not collapse a node that lands *between* those modes.

---

## 2. The released v0.7.1 KDE — how it works (and the epsilon floor that made it work)

Module: `gdna_density_prior.py` (`GdnaDensityPrior` + `build_training_substrate`). A **weighted Gaussian KDE in
log-density space**:

```
    P̂(log ρ_g)  =  Σ_i  w_i · φ_h( log ρ_g − â_i ) / (h · Σw)  ,   â_i = log(f_g,i · M_i / E_i)  floored at 1/E_i
```

- **Additive, one kernel per node.** No EM, no competition — a minority (capture-enriched) population **cannot be
  competed away**. This is the load-bearing property (`kde_vs_npmle_enriched_mode`).
- **UNIT weight `w_i = 1`.** Every solved node is one observation of a genomic location's density. **Precision is
  deliberately NOT a weight** — precision correlates with density (the depleted floor is low-count/low-precision,
  the enriched mode high-count/high-precision), so precision-weighting would systematically down-weight the whole
  floor mode. Noise is handled by the **bandwidth**, not the weight (design §8e).
- **Bandwidth `h`** — swappable: Silverman (Kish-effective-count + robust weighted-IQR scale), likelihood-CV, or
  fixed; **floored at the weighted-median per-node noise scale** `√(Var(log f_g)+1/(g+1))` so sampling noise is
  not resolved into spurious modes.
- **Real quadratic tails in the solve.** `logpdf_kernel` (the true kernel sum, `−½((x−x_i)/h)²` tails) is used in
  the per-node ψ — **not** the clamped `logpdf` interpolation. The real tails **suppress high-density false
  positives**: a node whose implied `ρ_g` is far above the population gets a genuine penalty (a clamped flat tail
  would let it drift to `f_g≈0.5` → catastrophic FP).
- **THE EPSILON FLOOR (the mixture-bridge ε).** The KDE fit on clean unimodal region nodes has a **deep valley**
  between the depleted and enriched modes (~10² nats). A node whose belief density lands *in that valley* (a
  capture boundary, a sparse-probe region) collapses to `f_g≈0` — **the crush**. The fix: `_kde_logprior` mixes
  the KDE with a **uniform "any-mixture" bridge over the observed support at weight ε**:
  `P = (1−ε)·KDE + ε·Uniform_support`. This **floors the valley** (no collapse) while leaving the KDE's real tails
  *outside* the support intact (FP suppression unchanged). `ε=0` = bit-exact legacy KDE (which crushes); any small
  ε defeats the collapse cliff. *This is the "released KDE crushed until we added the epsilon floor and it began
  to work" — the ε bridge is that floor.*

**The two crush mechanisms — keep them distinct.** They have different fixes and both must hold:

| crush | cause | fix |
|---|---|---|
| **A — valley collapse** | a node lands in a −∞ valley *between* the prior's modes | the **epsilon floor** (ε-bridge / interior floor) — floor the valley |
| **B — silenced enriched mode** | the fit has *no* enriched mode ⇒ the prior's PEAK is depleted | the **additive kernel** (+ unit weight + no τ-discount) — keep the enriched mode |

The released KDE fixes **both** (additive+unit-weight ⇒ no B; ε-bridge ⇒ no A). This branch's NPMLE fixes A
(`floor_eps`) but reintroduces **B** — hence the crush.

---

## 3. This branch's NPMLE — the Poisson-rate reframe, and why it crushes

Module: `npmle.py` (`DensityNPMLE`). The refactor's *motivation was sound* — it targets the **zero problem** the
density KDE cannot solve (`archive/gdna_prior_zero_handling.md`): a density KDE floors `ρ_g` at `1/E`, so it
**cannot represent `ρ_g=0`**, and 64–94% of confident nodes on real cfRNA are count-0. The Poisson-rate model
`k ~ Poisson(ρE)` makes **zero native, discreteness native, precision native** — the right generative object.

But the *estimator chosen* (a competitive EM mixture) reintroduces crush-B through three compounding mechanisms
(each dissected in `npmle_crush_dissection.md` §2):

1. **Per-node τ-width** (`_cell_loglik`): the enriched nodes are the low-count/unstranded/imprecise ones → their
   likelihood bump is **broad and short** before any competition.
2. **Competitive EM** (`_em_weights`): components compete for responsibility → the sharp depleted majority wins,
   the broad enriched minority is **competed toward zero**.
3. **The `n_regions` background tower** (fit, npmle.py:264–276): the pooled background is one cell weighted by the
   *whole* intergenic count at genome-scale `ΣE` — a razor-sharp depleted mode that **dominates the EM**.

Net: the fitted density is depleted-dominated, its top `hi = max(â)+2h` is set by pass-0's own under-called gDNA
so the true enriched level clamps **off-grid**, and the δ-pin crushes. The `floor_eps` interior floor fixes A but
is powerless against B (a depleted *peak* is not an empty valley).

---

## 4. KDE vs NPMLE — advantages & disadvantages

| axis | density KDE (v0.7.1) | Poisson-rate NPMLE (this branch, EM) |
|---|---|---|
| enriched-mode survival (crush-B) | ✅ additive ⇒ minority preserved | ❌ competitive EM ⇒ minority erased |
| zero / absence (`ρ_g=0`) | ❌ density-space, floors at `1/E` ⇒ FP over-call | ✅ Poisson zero-native (`k=0` supports `ρ→0`) |
| discreteness (small-`k` comb) | ❌ smears the `k/E` comb into a fake "enriched mode" | ✅ Poisson is the correct discrete likelihood |
| precision handling | unit weight + bandwidth-floor (noise ≠ weight) | ✅ likelihood *is* the precision — but see below |
| how precision is USED | ignored except as a bandwidth floor | ❌ as τ-width in a **competitive** likelihood ⇒ discounts enriched |
| background/floor | bolt-on depleted scalar + ε-bridge (a patch) | aggregate cell — but as an `n_regions` **tower** ⇒ over-crush |
| complexity / determinism | simple, additive, deterministic | EM (more machinery); additive path is deterministic |
| valley collapse (crush-A) | fixed by the ε-bridge | fixed by `floor_eps` |

**Reading:** the KDE's *representation* (additive, unit-weight, real-tail, ε-floored) is right for **crush
avoidance**; the NPMLE's *generative framing* (Poisson-rate, zero/precision-native) is right for the **zero &
discreteness** problems. The target is the **union**: an additive, occupancy-weighted, **Poisson-aware** density
that keeps the enriched mode AND represents zero — which is exactly the `additive=True` path plus a correctly
calibrated floor and a principled precision treatment. Neither shipped form is that yet.

---

## 5. The training subset — who is in, who is out, and why

The fit reads the **solve's own deconvolved gDNA** `ĝ = f_g·M` per REGION node. Non-circularity + trustworthiness
drive inclusion:

| node class | mask | in the fit? | why |
|---|---|---|---|
| **single-strand region** (`free_pos ^ free_neg`) | `single` | **YES** | the strand pins f_g ⇒ the deconvolved gDNA is trustworthy; carries the real enriched (exon) + expressed-depleted (intron) structure |
| **structural-gDNA / intergenic** (`~free_pos & ~free_neg`) | `gonly` | **depends** — the FLOOR | locked to f_g=1; the *depleted anchor* (incl. count-0: "gDNA absent here" is the strongest depletion evidence). Released KDE + EM-hybrid: fold into the **background/floor**, not the individual modes (else it floods the depleted mode into a tower). Additive path: excluded from `sel`, represented by the weak floor. |
| **AMBIG** (`free_pos & free_neg`) | `fp & fn` | **NO** | the two-root ambiguity the prior must *resolve* — training on it is circular |
| **boundaries** | `kind != REGION` | **NO** | different geometry (crossing vs contained), short near-fixed `E` ⇒ a discreteness comb (`gdna_prior_zero_handling.md` §4); the density substrate would mistake quantization for enrichment |
| **dead** (`E≤0` or `M≤0`) | `~live` | **NO** | no observation |

Current code (`calibrate.py:_fit_gdna_hyperprior`): EM path `sel = single | gonly` (intergenic→aggregate cell
under HYBRID), precision-weighted. Additive path `sel = single` only (intergenic→weak floor), unit weight
(`var_g=None`). **Open question for review:** the additive path *excludes* intergenic entirely and leans on the
floor; the released KDE *included* zero-count intergenic/intronic as the depleted anchor. Which representation of
"absence" is better — data points at `1/E`, a Poisson `k=0` likelihood, or a separate floor — is a live §8 lever.

---

## 6. Precision — unit weight vs precision weight vs "one count, spread by precision"

This is the crux the owner surfaced, and it has three genuinely distinct options. Each node is fundamentally **one
observation** of a genomic location's gDNA density; the question is how its *uncertainty* enters.

- **(a) Precision weight (current EM).** Weight `∝ 1/var`. **Rejected by evidence.** Precision correlates with
  density: the enriched exons are low-count/imprecise, so precision-weighting **down-weights the whole enriched
  mode** — a direct cause of crush-B. It changes *occupancy* (mass), which is exactly what must be preserved.

- **(b) Unit weight + fixed/global bandwidth (released KDE).** Every node = area-1 kernel of common width `h`.
  **Occupancy = height**, so the enriched minority renders at its true share — this is why it avoids crush-B. But
  it **ignores per-node precision entirely** (except as a bandwidth *floor*): an imprecise enriched node (M≈10) is
  drawn as sharp as a precise one, over-claiming the mode's location (restore-plan OI-5).

- **(c) One count, spread by precision (the owner's proposal).** Every node still contributes **area 1**
  (occupancy preserved — no crush-B), but its kernel **width = the node's own precision**: imprecise → **fat &
  short**, precise → **tall & skinny**. This is the honest Bayesian object — a node's contribution to the
  population density *is* its posterior over its own log-density (a unit-mass Gaussian of width `τ_i`), and the
  population density is their **sum** (a mixture of per-node posteriors, i.e. an additive KDE with **per-node
  bandwidth**).

**Compare/contrast.** (c) is strictly better than both (a) and (b) for our failure modes:

| property | (a) precision-weight | (b) unit + fixed-h | **(c) unit + per-node-τ width** |
|---|---|---|---|
| occupancy (mass) preserved? | ❌ shrinks enriched mass | ✅ area 1 each | ✅ area 1 each |
| crush-B (enriched mode) | ❌ silenced | ✅ preserved | ✅ preserved |
| honest per-node precision | as *weight* (wrong axis) | ignored | ✅ as *width* (right axis) |
| imprecise enriched node | down-weighted away | drawn falsely sharp | ✅ a broad, honest bump |
| relation to the EM τ | competitive likelihood | none | **additive** likelihood — the fix |

The key distinction from the EM's τ (mechanism 1 in §3): the EM used τ as the width of a **competitive
likelihood**, so the broad enriched bumps were then *competed away*. Option (c) uses τ as the width of an
**additive unit-mass kernel** — the broad bump keeps its full mass and simply says "gDNA is around here,
imprecisely," which is the truth. **(c) is the recommended default** — it unifies "additive KDE" (§2) and
"adaptive bandwidth" (§7) into one correctly-specified object, and it is the honest per-node posterior sum.

**Per-node τ, restated.** τ = `√Var(log f_g)` is the pass-0 belief width. It is *poison as a weight* (a), *lost if
ignored* (b), and *correct as a per-node kernel width* (c). The bandwidth-floor still matters (a node can't be
sharper than its sampling wall `1/√g` regardless of belief), so the per-node width is `max(τ_i, noise_floor,
h_min)`.

---

## 7. Bandwidth options

Bandwidth controls mode resolution vs over-smoothing. Options, roughly increasing sophistication:

1. **Fixed global** (current NPMLE `npmle_bandwidth=0.15`). Simple; but one width can't be right for both the
   sharp depleted mode and the broad imprecise enriched mode ⇒ over-smooths the tight modes and under-smooths the
   fat ones. The visible cause of the KDE's over-broad captured-mode bumps in the EM-vs-KDE overlay.
2. **Silverman / Kish** (released KDE default). Data-driven global width from the (weighted-effective) sample size
   and a robust IQR scale, **floored at the median per-node noise**. Better than fixed; still one global width.
3. **Likelihood-CV** (released KDE optional). Picks `h` by cross-validated fit; costlier, can over-fit sparse
   data.
4. **Per-node adaptive = §6(c).** The width IS the per-node precision (`max(τ_i, noise, h_min)`). This is the
   principled endpoint — bandwidth stops being a single hyperparameter and becomes the honest per-observation
   posterior width. **Constraint (restore-plan OI-4/OI-5):** it must be a *density-resolution* width, and it must
   preserve occupancy (unit area). Do **not** re-import τ as a *likelihood* width or a *weight* — that is the EM
   crush. A global `h_min` (mode-separation ceiling `h ≤ Δ/2` between the −1 enriched and −3 expressed modes)
   still guards against noise-driven spurious modes.

---

## 8. Levers & experiments (the immediate work)

Three levers; the additive kernel (§6c foundation) is settled — adopt it. The two open calibrations:

**Lever 1 — additive kernel (SETTLED).** Turn `gdna_prior_additive` ON. It is the only thing that avoids crush-B
(EM-vs-KDE overlay: captured-unstranded EMD 2.09→0.59 etc.). Non-negotiable foundation.

**Lever 2 — floor strength (OPEN, the priority).** The additive KDE alone regresses the **zero-gDNA** scenarios
(`none_*` EMD 0.28→1.15) because its 1-pseudo-observation floor can't pull pass-0's zero-gDNA over-call down to
the true depleted spike, while the EM's `n_regions` tower can (but over-crushes the enriched). The right object is
a floor **at `ρ_bg`** whose strength: (i) pulls true-zero nodes to `f_g≈0`, (ii) leaves enriched nodes (whose true
`ρ_g` is *above* `ρ_bg`) alone. **Experiment 2:** additive ON, sweep the floor weight `w_floor` (and whether
`h_floor` tracks `σ_bg`); plot the **frontier** of zero-gDNA-FP (`none_*`) vs captured-recovery (`gdna100/300
ss0.50 on`, node 1055) and find the knee. Watch the stranded controls stay flat. *(Flag: `w_floor` is a new
scalar — discuss its principled value; "one pseudo-observation" was the KDE's answer, but the sim under-represents
the intergenic flood vs real data.)*

**Lever 3 — precision→bandwidth (OPEN).** Implement §6(c) — per-node unit-mass kernels of width
`max(τ_i, noise, h_min)`. **Experiment 3:** A/B (b) fixed-h vs (c) per-node-τ-width on the enriched-mode
recovery + the mode-separation (does the imprecise enriched node render honestly broad without leaking FP?).

**Experiment 1 (do first, cheap):** re-run the KDE A/B gates now that pass-0 is honest — `gdna5_capON` recovery,
`gdna_none` FP, stranded no-regression — to confirm the upstream pass-0 work removed the old blocker before
tuning the floor.

**Validation harness for all:** `refit_gdna_movie.py` / `kde_vs_em_compare.py` (refit-vs-oracle EMD by node type
over the 32 scenarios), with node 1055 and `none_*` as the two poles, and `hyperprior_suite_vs_oracle.py` for the
per-scenario density overlay. Real-data (LBX0190 / MO_3005) is the final gate for `w_floor` + bandwidth (the sim
under-represents the intergenic flood).

---

## 9. Alternative approaches (if the additive-KDE-plus-floor path stalls)

- **A. Enrichment-conditioned prior (the principled generalization).** The root defect is pooling the *observed*
  density `ρ_g = e(x)·d_g` across enrichment levels. Condition on the node's enrichment `e(x)` — already estimated
  as `mu_proj` by the Role-A enrichment NPMLE — and model the **intrinsic rate** `d_g = ρ_g/e(x)`, whose
  population *is* roughly unimodal (genomic gDNA is ~uniform pre-capture). Then a captured node's prior expects
  enriched gDNA by construction and the δ-pin no longer crushes. This is "keep the enriched mode" done from first
  principles rather than by preserving a fitted bump. Cost: couples Role B to Role A; needs care that `mu_proj`
  isn't itself composition-contaminated. **The strongest long-term candidate.**
- **B. Zero-inflated / spike-and-slab KDE.** `π·δ(ρ_g→0) + (1−π)·KDE`, `π` = the measured count-0 fraction
  (`gdna_prior_zero_handling.md` §3.2). Makes absence explicit, subsumes the floor. Density-space, so it still
  inherits the discreteness comb.
- **C. Poisson-rate deconvolution done right.** Keep the NPMLE's Poisson framing but drop the three crush
  mechanisms: additive (not competitive) weights, no `n_regions` tower (a `k=0`-native background instead of a
  tower cell), and precision as width not weight. This is essentially §6(c) expressed in rate-space and is the
  cleanest reconciliation of §4's two columns.
- **D. Redescending / robust prior** (`archive/robust_redescending_prior_design.md`, `unified_prior_redescending_
  design.md`) — bound the prior's pull so no single mode can crush a node; a different tool for crush-A/B that the
  ε-floor already approximates.

---

## 10. For a third-party reviewer — code modules & where bugs hide

**Modules to read (in order):**

1. `src/rigel/calibration/npmle.py` — `DensityNPMLE.fit` (the EM mixture + the `additive` KDE path), `_cell_loglik`
   (τ-width), `_em_weights` (competition), `_kde_density` (the additive kernel), the aggregate background cell
   (fit:264–276), the interior `floor_eps`, and **`logprior`** (the δ-pin — the load-bearing coupling).
2. `src/rigel/calibration/calibrate.py` — `_fit_gdna_hyperprior` (the training subset masks) and the Phase-2 refit
   wiring (`calib_refit_iters`, `_debug["gdna_hyperprior"]`).
3. `src/rigel/calibration/simplex_logodds.py` — `_gdna_arm` / `_solve_nodes_logodds_all`: how the `(n,K)`
   `logprior` term enters ψ and combines with the strand likelihood and messages (the crush is realized here).
4. `src/rigel/calibration/background_reference.py` — `measure_background` → `ρ_bg` / `log_rho_floor` / `σ_bg` (the
   floor location & softness).
5. `gdna_density_prior.py` @ tag `v0.7.1` (`git show v0.7.1:...`) — the released KDE: `GdnaDensityPrior`,
   `build_training_substrate`, `logpdf_kernel` (real tails) vs `logpdf` (clamped), the mixture-bridge ε.

**Where a bug would live (the high-suspicion surface):**

- **The δ-pin (`logprior`)** — `log ρ_g = log f_g + log(M/E)`. This is *correct* but is the amplifier: any error
  in the prior's mode location becomes an f_g error scaled by `ρ_tot`. Check the grid top `hi=max(â)+2h` clamps
  the enriched mode *off-grid* under an under-called pass-0 (fit:217, 255).
- **The training masks** (`_fit_gdna_hyperprior`) — an off-by-one in `single`/`gonly`/`intergenic`
  (signature-index clamping at chain.ref_idx vs region_arrays.signature) silently changes which nodes flood the
  depleted mode.
- **The `n_regions` tower** (fit:273–276) — weighting the background cell by the whole intergenic count is the
  crush-B accelerant; verify it is *dropped* in the additive path (it is: additive uses `w_floor=1`).
- **Precision as weight vs width** — `var_g` is passed (EM) or `None` (additive). Confusing the two (or
  reintroducing τ as a likelihood/weight in the additive path) re-creates the crush.
- **The floor location** — `background.log_rho_floor` sign/units (natural-log vs decade) and the `Σg/ΣE` blend;
  a mislocated floor injects the false-positive tail seen on `none_*`.
- **Role-A/Role-B leakage** — `bp_solver` must use the *enrichment* NPMLE for `project`/σ²_transfer and the *gDNA*
  hyperprior only for `logprior`; a swap would let a composition prior damp messages or vice-versa
  (`gdna_kde_restore_plan.md` §2, OI-9).

---

## Appendix — the one-paragraph synthesis

The gDNA prior is a density prior applied through `f_g = ρ_g/ρ_tot`, so it must be **bimodal and
enrichment-respecting** (depleted floor for true-zero, enriched mode for captured gDNA) and must not collapse a
node **between** modes. Two crushes: **A**, a between-mode valley (fix = the epsilon floor / ε-bridge); **B**, a
silenced enriched mode (fix = an **additive**, occupancy-preserving fit). The released v0.7.1 KDE fixed both
(additive unit-weight kernels + the ε-bridge) but was density-space (blind to `ρ_g=0`). This branch's EM-NPMLE
fixed the zero problem (Poisson-native) but reintroduced **B** (competitive EM + `n_regions` tower + τ-as-weight).
The synthesis: an **additive, Poisson-aware, occupancy-weighted** density with a **calibrated floor** and
precision expressed as **per-node kernel width** (one count, spread by its own uncertainty), optionally
generalized to condition on enrichment (model the intrinsic rate `d_g=ρ_g/e(x)`). Adopt the additive kernel; the
open work is the floor strength and the precision-as-width bandwidth.
