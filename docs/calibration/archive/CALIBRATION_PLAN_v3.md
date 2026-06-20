# Calibration — the iterative Bayesian deconvolution (PLAN v3)

> **⚠️ SUPERSEDED by `CALIBRATION_PLAN_v4.md` (2026-06-16, same day).** v4 overturns v3's "RNA chains"
> framing with the **gDNA-chains / RNA-local inversion** and adds the **unified node solver** (regions and
> boundaries solved by one 3-tier posterior, eliminating `deconv_sides` + `I₀` + side-freezing). v4 is the
> execution-ready plan. v3 is kept for the design rationale (the 6-fork synthesis) it records.

**Status:** SUPERSEDED design plan (2026-06-16). **Supersedes** `CALIBRATION_PLAN_v2.md`
(2026-06-15) and folds in the user's mental-model draft (`calibration_nodes.md`, `my_notes.md`). v2
described the iterative-bootstrap *target*; v3 records that the target is **largely built** and specifies
the principled **completion** — the resolved architecture forks, the magic-number → derived-quantity
program, and the cleanup — that takes the code from "runs but incomplete" to **production**. Where this and
any older note disagree, this wins. Companion docs kept (not duplicated): `per_node_deconv_hierarchy_design.md`
(the per-node math), `effective_length_redesign_plan.md` §8 (the eff-len IPR, SHIPPED), `variance_model.py`
docstring (the SCAM fitter).

The production bar is four properties, used as the acceptance test for every decision below:
**ROBUST** (no brittle tunable — every constant is derived from data, canonical, or a documented
numerical-resolution knob with a recorded validation), **ELEGANT** (simple; behaviour *emerges* from
honest precisions, not hand weights), **EFFICIENT** (genome-scale), **ACCURATE** (optimal given the
identifiability limits).

---

## 0. TL;DR — the thesis

Calibration deconvolves each genomic region's **unspliced** fragment mass into the per-node 2-simplex pie
`(f₊, f₋, f_g)` = sense-RNA / antisense-RNA / gDNA (RNA-vs-gDNA **only** — the per-locus EM separates
mature from nascent downstream). `f_g` feeds the per-locus EM prior; the deconvolved gDNA mass feeds the
effective-length contraction.

The elegant core, and the whole point of v3:

> **One inference object, one variance machine, one combination rule.**
> The per-node belief is solved **exactly** on the 2-simplex by a grid sum-product (belief propagation).
> A **single gDNA `var~mean` machine**, re-fit each pass on the running gDNA estimate, emits **every**
> precision in the system (count, global, RNA-coupling). The combination is ordinary **precision-weighted
> Bayesian updating**. Every node-class behaviour — pure-gDNA concentration, zero-DNA pinning, AMBIG
> deference, capture smear-cap, the priority *strand > count > global* — **emerges** from that; there are no
> hand-tuned weights in the steady state.

The loop:

```
SETUP (once):  κ + the two strand Beta-Binomial overdispersions (capture-invariant) ; the boundary SIDES
INIT:          all unspliced mass is gDNA  (the conservative worst case)
REPEAT (≤ few passes, stop on the prediction-error plateau):
   1. RE-FIT  ρ_global (count-observable nodes only) + the ONE gDNA var~mean (DIRECT + IMPUTATION)
              on the current gDNA estimate
   2. DERIVE  τ_count (per node), σ²_global (the zero-DNA pin), q_rna (per edge)  — all from that one machine
   3. SOLVE   every node's 2-simplex posterior exactly (strand likelihood + odds propagation
              + count prior + node-class prior), read the MEDIAN f_g and its variance
OUTPUT:        per-node {f₊, f₋, f_g} + uncertainty → CalibrationResult → the per-locus Dirichlet + eff-len
```

This is **not a rewrite.** The loop, the 3-tier solve, the SCAM fitter, the median readout, and the
genome-scale edge-matrix caching are already built. v3 is the principled finishing: derive the three
precisions that are still placeholders, surface the uncertainty that is computed-then-discarded, delete the
dead scaffolds, and stop the docstrings from lying.

---

## 1. The starting point — code is ahead of its docs

Established by a full plan-vs-code audit (`/tmp/calib_synthesis.md`, 2026-06-16). What is **built and
correct** (`calibrate.py:238-294`): the iterative all-gDNA bootstrap; the 3-tier per-node posterior
(`simplex_sweep._local_loglik`); `τ_count` driven by the var~mean (the v2 `β=10` placeholder is genuinely
retired — `count_precision=tau_count` is always passed); the node-class prior (Jeffreys at single-strand,
global gDNA at AMBIG); posterior **median** `f_g` + variance; the `MonotoneVarMean` SCAM fitter; the
decoupled-edge `O(P)` short-circuit that fixed the ~31 GB/chain genome-scale wall (`simplex_sweep.py:143-188`).

What is **placeholder or absent** — the work of v3 (each resolved below):

| Gap | Where |
|---|---|
| `σ²_global` is the robust-MAD of node densities, **not** the documented DIRECT-curve-at-ρ_global | `calibrate.py:269-274` |
| the IMPUTATION `var~mean` is the *confounded* left-vs-right disagreement, with **no** df/Jensen correction | `variance_model.py:277-284` |
| `q_rna = 0.25` hard-coded scalar; never passed; no config field | `simplex_sweep.py:194` |
| per-node uncertainty `_fg_var` is **computed then discarded**; `CalibrationResult` has no uncertainty field | `calibrate.py:294`, `result.py` |
| the convergence stop `1e-3` is an inline literal; no prediction-error monitor, no oscillation guard | `calibrate.py:291` |
| dead scaffolds: `count_trust_beta`, the `density_model` LOESS path, `simplex.solve_node`, `rna_variance.py`/`mature_density.py` | (D) below |
| ~6 docstrings claim "single feed-forward pass / no EM loop / acyclic" while the code iterates | `config.py:226`, `calibrate.py:1,83`, `result.py`, `__init__.py`, `priors.py`, `substrate.py` |

---

## 2. First principles (the user's, formalized)

1. **Strand is unbiased and reliable → the priority/direction signal.** Use it wherever it is informative.
   Its information is the Fisher information `I = N·(2κ−1)²` (depth × specificity); a Beta-Binomial models
   the overdispersion. It is **capture-invariant** (a probe enriches both strands of a duplex equally).
2. **Count is valuable but biased-under-capture, with *local* reliability → the magnitude fallback.** Its
   reliability tracks unobservable probe-to-boundary distance, so it varies node-to-node. We **down-weight
   it by its reliability; we never bias-correct it** (the capture bias is exon-vs-intron-specific and does
   not transfer). It must take over where the strand is silent and yield where the strand speaks.
3. **The global gDNA baseline is the foundation** for truly-unanchored nodes — numerical stability + the
   zero-DNA fix. It must never overtake strand or count.
4. **The fractional (pie) representation caps capture "smear."** A *fraction* propagated/imputed into a
   low-count region is projected onto that region's *observed* counts, so the smear is bounded by the local
   evidence — the reason the latent is a simplex, not a count.
5. **Behaviour emerges from honest precisions, not hand weights.** The `strand > count > global` hierarchy
   is the consequence of a sharp likelihood overwhelming a wide prior — not of a tuned coefficient.
6. **Capture-awareness is partitioned into exactly three places**, never duplicated:
   - **direction** → Tier-1 strand (capture-invariant);
   - **local magnitude** → Tier-2 count (boundary crossing-flux is capture-aware) **and** the eff-len IPR
     contraction (`priors.py` / `capture_eff_length.py`);
   - **baseline level** → Tier-3 global foundation — deliberately **structure-free** (uniform).
7. **Identifiability honesty.** gDNA is spatially smooth (boundary-predictable); mature RNA is spiky;
   nascent RNA is smooth. The smooth-gDNA / smooth-nascent split is **only** resolvable by the strand; in
   unstranded data it is the honest identifiability limit and the global baseline governs. Calibration
   reports this floor; it does not manufacture certainty there.

---

## 3. Assimilating & critiquing the mental model

The user's draft (`calibration_nodes.md` + `my_notes.md`) is the catalyst for v3. Three proposals are
**overridden** (with reasoning); four are **adopted** because they converge with open problems already in
the code. This section is the honest critique the user asked for.

### 3.1 Overridden — and why

**(a) "Thick" belief propagation (count/gDNA imputation flows as learned-reliability messages between
nodes).** *Overruled.* The two things the thick vision actually wants **already exist** inside the thin grid
engine: the per-node `{f, unc}` object **is** the grid belief's marginal (`_fg_median` + `_fg_var`), and the
"learned communication error" **is** the IMPUTATION `var~mean`. Replacing the exact simplex belief with a
Gaussian `{f, unc}` message would be a strict *downgrade* — it fails exactly at `f_g→0` / `f_g→1`, the
boundary where the posterior-**median** readout was introduced to kill a +8.7-point regression. And gDNA
imputation is fundamentally a **one-hop** boundary→region operation (an observable boundary directly anchors
its adjacent region); it is **not** transitively continuous the way the per-strand RNA log-odds are along an
exon chain — propagating gDNA density *across* an exon↔intron boundary would be physically wrong. So we keep
the thin engine and deliver the vision's real wants additively (§9). *The user's instinct that the
two spatial mechanisms feel un-unified is correct; the resolution is to surface the marginal and pin
`σ²_global` to the one machine, not to add a coupling channel — see the deferred refinement in §15.*

**(b) The global hyperprior as `scalar density × per-node enrichment factor`.** *Overruled — it is a
tautology, which the user spotted ("isn't this 1.0 everywhere?").* Define enrichment as `local/global`; then
the prediction `global · (local/global) = local` identically, the residual is ≡0, and the "uncertainty"
is ≡1 by construction. It is *doubly* circular: the unanchored AMBIG nodes the foundation exists to serve
have **no** local density to form a factor from. Eliminate the factor entirely (§8). Capture-enrichment is
already carried correctly by Tier-2 and the eff-len IPR (principle 6); the foundation is the structure-free
off-target floor. (Same verdict retires the related "global Γ enrichment edge-scalar" idea from
`my_notes.md` — same-class structure makes it unnecessary and the iteration recovers the signal more safely.)

**(c) A typology-*staged* Pass-0 (Group-1 seeds, Group-2 solved by strand-ALONE, Group-3 → weak prior,
then fit models, then solve).** *Overruled as an algorithm; adopted as a lens (§4).* A staged AMBIG→all-gDNA
init **manufactures** the very phantom the global prior was built to remove (the Jeffreys vertex-push), and
hard-zeroing the count at high-SS exons already regressed mRNA by 66k in testing
(`per_node_deconv_hierarchy_design.md` §14). A persistent per-node geometry object is also the wrong
direction at genome scale. The single all-gDNA-init sweep **already realizes** the Group-1/2/3 behaviour
through emergent precision — which honours the user's own principle 5 better than an explicit ladder does.

### 3.2 Adopted — the convergences

**(a) The source→destination prediction-error imputation model** (`calibration_nodes.md` "NEW IDEA"): fit
*predicted-vs-actual* on node pairs where the destination is **observable** (ground truth exists) → the
residual *is* the imputation reliability; refit each pass. **This is exactly the fix** that
`per_node_deconv_hierarchy_design.md` §17 and `var_fit_plan_evaluation.md` flagged as "intended but not
built" — the current code instead fits the *confounded* left-vs-right disagreement. Adopted as the IMPUTATION
`var~mean` re-cut (§7, Stage B, measurement-gated).

**(b) First-class per-node uncertainty** `{f±, f_g} + {unc±, unc_g}` as an output. Adopted: the variance is
already computed (`_fg_var`) and merely discarded — surface it (§9, Step 1), then use it to precision-weight
the per-locus prior (§9, Step 2, gated).

**(c) `var~mean` for both gDNA and RNA, re-fit each iteration.** Adopted as the *one-machine* elegance (§7).
The RNA `var~mean` feeds `q_rna` only (conditional, §7).

**(d) "Behaviour emerges from honest precisions, not hand weights."** Adopted as the governing principle —
it is what justifies deriving `σ²_global`, deferring the `I₀`/`q_rna` over-derivation until measured, and
reading the priority hierarchy off the precisions rather than a coefficient.

---

## 4. The node typology — an explanatory lens, not a staged algorithm

The user's Group-1/2/3 partition is the right way to **understand and explain** solvability. We keep it as
documentation and prove that the *single* all-gDNA-init sweep realizes each class through emergent precision
— there is no special-cased init.

| Class | Definition | How it is solved (emergent) |
|---|---|---|
| **0-D seed** | intergenic region / intergenic-exon boundary (`strand_class = NONE`) | only the `f_g=1` vertex survives the strand mixture; the forbid mask zeroes `f₊/f₋` (`simplex_sweep.py:104-107`). A propagation sink. |
| **1-D single-strand** | one strand active (intron, exon, exon↔intron boundary; `TS_POS`/`TS_NEG`) | the strand likelihood is **sharp** (high `I = N(2κ−1)²`) and dominates; the count prior is down-weighted by `(1−w)`, `w = I/(I+I₀)` → count yields (`:235-240`). Jeffreys Beta(½,½) reference prior concentrates `f_g` (`:99-103`). At `ss→½` the strand goes flat and the count governs — the same node, no branch. |
| **2-D AMBIG** | both strands active (`(ex+∨in+) ∧ (ex−∨in−)`) | the strand tilt + propagated odds resolve it where anchored; where balanced **and** odds-starved it is genuinely under-determined and falls to the count, then the global foundation (`μ_global`, `:99-103`). |

**"Solvable" is a continuous spectrum, not a binary.** The right reframing of the user's "can we reduce
uncertainty?" is the weight `w = I/(I+I₀)` (and the posterior variance): it slides smoothly from
strand-governed (`w→1`, high SS + depth) to count/global-governed (`w→0`, unstranded or thin). The typology
names the *corners* of this spectrum; the math handles the interior.

---

## 5. The architecture — one uniform iterative loop

**SETUP (once, pre-loop; the capture-invariant / direction quantities).**
- κ (`rna_sense_frac`) and the two Beta-Binomial overdispersions: gDNA (mean ½, seeded from count-observable
  + seam nodes) and RNA (mean κ, from boundary-side spliced). Both shrunk to the **same** default prior so
  unstranded data is uninformative. **Frozen for the run** (re-fitting them would feed the deconvolution's
  own gDNA weight back into the strand model — a feedback loop with no upside).
- The boundary **SIDES** (`deconv_sides`): solved once as the fixed gDNA anchors. They are the *known-answer
  reference* the IMPUTATION `var~mean` trains against and the per-locus flux the prior consumes; moving them
  each pass would move the convergence target (§12). **Frozen.** *(The one-time `I₀` retirement for the
  sides is a separate refactor — §15, not per-pass iteration.)*

**INIT — all-gDNA.** Every unspliced fragment is gDNA (`gdna_c = u_tot`). This is the conservative worst
case; it makes `ρ_global` a deliberate over-estimate and gives the `var~mean` support across the **full**
depleted→enriched range from pass 0 (which is *why* no `confidence(μ)` extrapolation guard is needed — every
node has training support at init; see §7).

**THE LOOP (≤ `sweep_max_passes`; stop on the prediction-error plateau, §12). Each pass RE-FITS from the
current gDNA estimate:**
1. `ρ_global` = `Σ gdna_c[observable] / Σ eff_len[observable]` — **count-observable nodes only** (never all
   nodes: an all-node average is irreducibly inflated by AMBIG/unstranded RNA the strand cannot remove, so
   it never drops and the tightening global locks in a phantom).
2. **The ONE gDNA `var~mean` machine** (`MonotoneVarMean`): DIRECT (own-count variance at observable nodes)
   and IMPUTATION (boundary→region prediction error), each read at the node's current gDNA density μ.
3. **Derive the three precisions from that one machine** (§7): `τ_count` per node; `σ²_global =
   direct.predict(ρ_global)·geom2`; `q_rna` per edge (conditional).
4. **SOLVE** every node (§6) → median `f_g` + variance → update `gdna_c`.

**Frozen across passes:** the strand overdispersions, the boundary sides, κ. **Re-fit each pass:** only the
*count-magnitude* quantities (`ρ_global`, the `var~mean`). This is the precise reading of v2's "every pass is
identical": the iteration refines magnitude/baseline; direction is fit once.

**OUTPUT.** `CalibrationResult` carries per-region/-side `gdna_mass`, `rna_mass`, `gdna_frac`, **and the
surfaced `gdna_frac_var`**; `assemble_priors` turns mass + uncertainty into the per-locus Dirichlet + the
gDNA eff-len (§9).

---

## 6. The per-node solve — the 3-tier emergent-precision hierarchy

Exact grid sum-product on the 2-simplex triangular lattice (`K = sweep_n_grid` points). The per-node log
posterior, normalized **once** over the lattice (the two quadratics compose into a single precision-weighted
Gaussian prior — one prior, no double-count):

```
log post(f₊,f₋,f_g) =  L_strand(f | tilt)              # Tier 1 — likelihood, capture-invariant DIRECTION
                     +  L_odds(f | neighbours)          # Tier 1 — the only propagated quantity
                     −  ½ · τ_count · (f_g − μ_count)²   # Tier 2 — count MAGNITUDE, down-weighted by (1−w)
                     −  ½ · τ_node  · (f_g − μ_node)²    # Tier 3 — node-class prior (see below)
```

- **Tier 1 — strand likelihood + propagated odds (priority, unbiased).** The 3-component mixture on the
  observed plus-fraction `p = ½f_g + κf₊ + (1−κ)f₋` (first-class at AMBIG: it is a plane through the simplex,
  not zero information), plus the sided spliced lower bound (single-strand), plus the **only** propagated
  quantity: the per-strand RNA:gDNA log-odds `log(f_c/f_g)` coupled along **same-strand-exon** edges with
  penalty `½(Δlo)²/q_rna`. All other edges decouple (`Q=∞` → the `O(P)` short-circuit). Precision is the
  intrinsic Fisher information; **never gated to `w=0`.**
- **Tier 2 — count prior (fallback magnitude).** Gaussian pull of `f_g` toward `μ_count = count_gdna_frac`
  (the boundary-anchored local density→fraction), at per-node precision `τ_count` (§7), **down-weighted by
  `(1−w)`**, `w = I/(I+I₀)` — so it yields where the strand is informative and governs where the strand is
  silent. Down-weight only; never bias-correct.
- **Tier 3 — node-class prior (foundation).** This term is **node-class-dependent** (the hard-won
  `per_node §13` finding — a uniform replacement regressed everything):
  - single-strand (strand-observable): **Beta(½,½) Jeffreys** reference — the correct binomial-proportion
    prior, concentrates `f_g` against the mixture's overdispersion spread (`μ_node`, `τ_node` are the
    Jeffreys form);
  - AMBIG / intergenic: the **global Gaussian** — `μ_node = μ_global = clip(ρ_global·eff_len/mass, 0, 1)` at
    `τ_node = τ_global`. The U-shaped Jeffreys would vertex-push a phantom here; the global baseline replaces
    it, so a pure-RNA library (`ρ_global≈0` ⇒ tight `τ_global`) settles an unanchored node at `f_g≈0`.

**Readout:** the posterior **median** of `f_g` (avoids the overdispersion-skew that biased the MAP) and the
posterior **variance** (the per-node confidence — surfaced in §9).

---

## 7. The variance machine — ONE `var~mean`, three uses

The keystone of the elegance: a **single** gDNA density `var~mean` curve, re-fit each pass on the running
gDNA estimate, supplies every precision. The fitter is the SCAM monotone-increasing P-spline
(`MonotoneVarMean`: `log var ~ monotone spline in log mean`; GCV-λ; bisquare-IRLS; power-law fallback for
tiny data) — validated to beat LOESS / isotonic / power-law and **kept**. Pure numpy+scipy, sub-second, no
C++ port.

**Two models, by estimation method (not by node class):**
- **DIRECT** — the variance of a node estimated from its **own** counts (Poisson/NB sampling). Feeds
  `τ_count` at count-observable nodes **and** `σ²_global`.
- **IMPUTATION** — the variance of estimating a node by boundary→region imputation. Feeds `τ_count` at
  imputed/AMBIG nodes. Properly **humbler** than DIRECT (it measures real prediction error).

**The three uses (`geom2 = (eff_len/mass)²` maps density-variance → fraction-variance):**
1. `τ_count = min(1/(var·geom2), mass_u)` per node — DIRECT at observable, IMPUTATION at imputed.
   **Clamp** `σ²_frac` at the Bernoulli bound `μ_count·(1−μ_count)` (a fraction-variance cannot exceed it).
2. **`σ²_global = direct.predict(ρ_global)·geom2`** — *the headline fix.* This is the **sampling variance of
   the baseline density estimate** ("how well do we know the baseline level"), not the between-node
   heterogeneity the current MAD placeholder wrongly answers. It is the **zero-DNA pin**: when observable
   seeds agree density≈0, `σ²_global` is tight at `μ_global≈0` ⇒ unanchored/AMBIG nodes pinned to `f_g≈0`,
   no phantom. Under capture the same curve is appropriately wide ⇒ defers to count/strand. **Four
   independent analyses converged on this** — it is the strongest cross-topic agreement.
3. `q_rna` (the per-edge RNA odds-coupling variance) — **conditional** (see below); the RNA `var~mean` (a
   revived `rna_variance.py`) supplies it only if Phase-1 measurement justifies it.

**Two corrections, both unconditional and free (Stage A):**
- **The Jensen / log-bias df correction.** `varmean_points` fits `log(raw_var)` directly, but for a sample
  variance with ν = k−1 dof, `E[log s²] = log σ² + [ψ(ν/2) − log(ν/2)]` and that bracket is **negative**, so
  the fit silently **under**-estimates `σ²` ⇒ over-confidence (verified **3.56×** at k=2, **1.78×** at k=3 —
  exactly where §15's humility is demanded). **Add the per-point positive offset
  `Δ_k = log((k−1)/2) − ψ((k−1)/2)` to `log(raw_var)` before the spline fit** (`k = kcount`, already carried
  at `variance_model.py:279`). This is the canonical digamma identity, no tunable.
- the Bernoulli clamp above.

**The IMPUTATION re-cut (Stage B — the user's "NEW IDEA", measurement-gated).** Today IMPUTATION is the
*confounded* left-vs-right disagreement split by `region_observable` (DIRECT = low-μ intron/intergenic,
IMPUTATION = high-μ exon — confounded with mean, not the reliability axis). The principled form is the
user's: learn the **held-out boundary→region prediction error** from the **all-three-observable triplets**
(both boundaries *and* the region observed → the region is the known answer → the genuine residual), on the
region-density axis. **Gate (open Q §15):** ship only after instrumenting (i) that `τ_count` survives
`geom2` as a live lever, and (ii) that clean all-three-observable triplets actually **span the exon μ-range**
— otherwise the re-cut re-introduces the domain mismatch in new clothing, and `raw_var~mean` + the Jensen
fix is the honest stopping point.

**Why no `confidence(μ)` guard** (v2 §4 called it "the lynchpin"): under all-gDNA init the curve is fit on
**every** node, so predictions interpolate and the curve **learns** genuinely-high variance at RNA-rich
means directly (the prediction error is large there) — `1/var` already yields the count where it is
unreliable. The latent risk is a **targeted panel** whose enriched-exon means lie outside the
intron/intergenic training range (flat extrapolation, `variance_model.py:179`); that is the one residual
need for the guard — deferred behind the targeted-panel stress sim (§15), not blocking.

---

## 8. The global foundation — a uniform scalar baseline

(Resolves the user's open questions A and B.)

- `ρ_global = Σ gdna_c[observable] / Σ eff_len[observable]` over **count-observable nodes only**, on the
  current gDNA estimate (iterating). Per node, `μ_global = clip(ρ_global · eff_len / mass, 0, 1)` — the
  count module's own fraction form with the *global* density substituted. Capture-awareness enters **per
  node through the node's own observed `mass`**, never through an enrichment ratio.
- `σ²_global = direct.predict(ρ_global)·geom2` (§7, use 2). Replaces the MAD placeholder.
- **No enrichment factor** (§3.1b — the tautology). Capture-enrichment lives in Tier-2 + the eff-len IPR.

**Answer to (A) — init the global scalar to zero or a floor?** *Neither.* Init it at the **all-gDNA
over-estimate** (it is just `ρ_global` under `gdna_c = u_tot`), then let it **decay** toward the truth as the
loop removes RNA from the observable nodes. A zero init makes the foundation inert on pass 0 and starves the
`var~mean` of its enriched range; a fixed nonzero floor is a magic number that would lock in a phantom in a
pure-RNA library. The over-estimate-then-decay is self-correcting; the *floor on the result* is structural
(restricting `ρ_global` to observable nodes ⇒ pure-RNA library converges to ≈0), not a constant.

**Answer to (B) — the global uncertainty, resolving "1.0 everywhere."** That value is not a bug to debug; it
is the mathematical signature of the enrichment-factor tautology (residual ≡0 ⇒ uncertainty ≡1). Eliminate
the factor. The foundation's honest uncertainty is `σ²_global = direct.predict(ρ_global)·geom2` — it answers
"how well do we know the *baseline level*" (tight when observable seeds agree density≈0; wide under
heterogeneous enriched seeds), which is exactly the quantity a foundation prior needs.

**Answer to (C) — the imputation message-passing formalization.** Messages are not a new engine (§3.1a). The
desired `{f, unc}` message **is** the grid marginal; the "learned communication error" **is** the IMPUTATION
`var~mean`, refit per pass against a **frozen** gDNA reference (the bracketing that keeps it convergent — do
*not* move the fit inside a message loop). gDNA imputation is one-hop; the only transitively-propagated
quantity is the RNA log-odds (already in the engine). The single genuine gap — runs of consecutive
non-observable regions where one-hop fails and `runfill_bidirectional` handles it *outside* the BP — is a
candidate same-class-only gDNA-coupling channel, **deferred** behind a measured failure demonstration (§15).

---

## 9. Uncertainty as a first-class output

Three independently-gated steps:

- **Step 1 (unconditional).** Add `gdna_frac_var` (+ left/right) to `CalibrationResult`; stop discarding the
  `_fg_var` the sweep already computes. Pure diagnostic surfacing; no behaviour change.
- **Step 2 (gated on the net-flow suite — open Q §15).** Precision-weight the per-locus Dirichlet by a
  **split-direction confidence** `c_locus` (inverse-variance, mass-weighted pooling of `_fg_var`),
  `a_g' = c·a_g + (1−c)·a_g_global`, `a_r'` likewise — blending an uncertain split toward the
  **global-baseline split**, *not* toward zero (which would reopen the capture leak). Discount prior-pinned
  (`strand_obs = False`) nodes — their tightness is the prior's, not the data's. **Why orthogonal:** the
  pseudocount *magnitude* `a_g + a_r` already encodes depth; `c` must act on the *direction* of the split,
  or it double-counts depth. Ship only if Phase-1 shows a measurable accuracy gain (memory: flagship is
  EM-bound / isoform-swap-dominated, the prior is ~3% of ss0.99 error — this may be an honesty feature, not
  an accuracy lever).
- **Step 3.** Keep `deconv_sides` and the strand overdispersions frozen (§5); the `I₀` retirement for the
  sides is a separate one-time refactor (§15), never per-pass.

---

## 10. Robustness — the magic-number → derived-quantity program

The honest constant budget (target ≤ ~8, `calibration_no_magic_numbers`) is only truthful once dead code
stops counting. Four buckets:

| Constant | Current | Verdict |
|---|---|---|
| **`σ²_global`** | `1.4826·MAD(densities)²·geom2` (`calibrate.py:269`) | **DERIVE NOW** → `direct.predict(ρ_global)·geom2`. The curve is already fit; just evaluate it here. |
| **Jensen df offset** | absent (`variance_model.py:282`) | **DERIVE NOW** → `Δ_k = log((k−1)/2) − ψ((k−1)/2)` added to `log(raw_var)`; `k` already carried. |
| **`σ²_frac` upper bound** | uncapped (`calibrate.py:267`) | **DERIVE NOW** → clamp at Bernoulli `μ_count(1−μ_count)`. |
| **per-node enrichment factor** | proposed (`calibration_nodes.md`) | **REJECT** — circular identity (§3.1b). |
| **Γ enrichment edge-scalar** | proposed (`my_notes.md`) | **REJECT** — unnecessary; EB loop recovers it. |
| **convergence stop** | inline `1e-3` (`calibrate.py:291`) | **PROMOTE** → `CalibrationConfig.sweep_convergence_delta` (default `1e-3`, golden-compatible) + prediction-error monitor + oscillation guard (§12). |
| **`I₀` `gdna_strand_info_scale`** | `10.0` (`config.py:288`) | **KEEP + DOCUMENT** as a validated monotone-saturating regularizer; fix the false "≈1 fragment" docstring. The emergent `w = I_strand/(I_strand+I_count)` is **RESEARCH/DEFER** — it needs an unidentifiable count-bias term; shipping it without that over-trusts the count exactly where it is most wrong. Retire `I₀` *for the sides only* via the §15 refactor. |
| **`q_rna`** | `0.25` hard-coded (`simplex_sweep.py:194`) | **DERIVE-IF-MEASURED** → per-edge RNA `var~mean` (revive `rna_variance.py`) only after Phase-1 shows propagation strength affects accuracy against an isoform-aware oracle; else **KEEP** `0.25` as a documented numerical knob with a recorded sensitivity sweep, and surface it as config. |
| **`τ_global` floor `1.0`** | inline (`calibrate.py:276`) | **KEEP-canonical + DOCUMENT** — the 1-pseudo-observation weak foundation; cap at `mass_u` stays. |
| **`count_trust_beta`** | `0.0` + elif (`simplex_sweep.py:48,89,193`) | **DELETE** — superseded by `τ_count`; dead at default. |
| **dead LOESS** `_LOESS_SPAN`/`_MIN_FIT`/`_BISQUARE_C` | `density_model.py:97-99` | **DELETE** — unreachable at quantile ½. |
| **`simplex.solve_node`; `rna_variance.py`/`mature_density.py`** | dead / import-only | **DELETE** (quarantine `rna_variance` if `q_rna` derivation is imminent). |
| **`_PRIOR_EPS=1e-3`, `sweep_n_grid K=60`** | `simplex_sweep.py:35`, `config.py:296` | **KEEP-numerical-knob + DOCUMENT** the K=20/60/100 accuracy-cost sweep. |
| **SCAM `k=18`, GCV λ-grid** | `variance_model.py:40` | **KEEP-numerical-knob** — basis upper bound; GCV-λ does the smoothing. |
| **`1.4826`, Tukey `4.685`, Beta(2,2) od=0.2 ceiling** | various | **KEEP-canonical** (normal-consistent MAD; 95%-efficiency Tukey; hard model bound). |
| **`POOL_EB_PRIOR_ESS=1000`** | `fl.py:59` | **FLAG** — over-shrinks tiny targeted panels per its own docstring; scale with pool size or move to config. |

Net steady-state model constants after this program: `I₀`, the `τ_global` floor, `q_rna` (if kept), plus the
numerical-resolution knobs (K, SCAM k) and the canonical statistics constants — within the budget, and each
either derived, canonical, or a documented-and-validated knob.

---

## 11. Efficiency — genome scale

- The grid engine is exact in two sweeps per per-reference chain. Edge `(P,P)` matrices are **cached** by
  `(q_pos, q_neg)` (≤ a handful distinct per chain), and fully-decoupled edges (`Q=∞`, the majority:
  exon↔intron / silent-strand) short-circuit to an `O(P)` constant message — this is what fixed the ~31 GB
  /chain wall. **Verify (not redesign)** the caching + chunking hold at genome scale.
- The per-pass `var~mean` refit (~10³–10⁴ points × ≤4 passes × 40-λ GCV × L-BFGS-B) is *asserted* sub-second
  but **must be measured** (open Q §15) — if the monotone fit dominates, the 4-pass loop is the bottleneck
  and the GCV grid / `maxiter` become the levers.
- Per-reference chunking caps memory; the loop is embarrassingly parallel across references.

---

## 12. Convergence & stopping

- **Mechanical stop:** promote the inline `1e-3` to `CalibrationConfig.sweep_convergence_delta` (default
  `1e-3` for golden-compatibility); break when `mean|f_g − prev_fg| < δ`.
- **Honest monitor:** trace the per-pass IMPUTATION boundary→region **prediction error** (and the mean
  posterior variance). It falls as RNA is correctly removed; it is the diagnostic that *justifies*
  `sweep_max_passes`. It is a **monitor, not the objective** (minimizing it alone would call smooth nascent
  "gDNA"; the objective is the full strand-anchored posterior).
- **Oscillation guard / damping:** add **only after** instrumenting the delta trace (v2 §10: "observe
  oscillation first"); the global anchor + monotone `var~mean` + frozen sides already bound it.
- **Frozen sides are load-bearing here:** the IMPUTATION `var~mean` (and Stage B) train against the sides as
  the fixed known-answer reference; moving them each pass would move the target and break the monitor (§5).

---

## 13. Implementation sequence (ordered, dependency-aware, measurement-gated)

Each phase's gate unblocks the next. Phases 0–1 are **truth + measurement** and decide whether the later
modelling is worth doing — this is the production discipline the user asked for (no speculative complexity).

**Phase 0 — Truth-in-docs + dead-code teardown (no behaviour change).**
- Fix the ~6 "single feed-forward / acyclic" docstrings → "iterative bootstrap"; fix `calibrate.py:4`
  (simplex axes are sense-RNA/antisense-RNA/gDNA, **not** mature/nascent); repoint broken refs
  (`phase2_design.md`, etc.) to v3 / `effective_length_redesign_plan.md` §8; rewrite the `docs/README.md`
  calibration section (it still anoints the superseded `CALIBRATION_PLAN.md`).
- Delete: `count_trust_beta` + branch, the `density_model` LOESS path, `simplex.solve_node` + companions,
  `strand_summary` gate, `rna_variance.py`/`mature_density.py` (quarantine if Phase-6 imminent); `git clean`
  the stale `.pyc`; fix/delete `scripts/debug/dissect_efflen.py`.
- **Gate:** full suite green (1122 baseline); zero production-path change.

**Phase 1 — Instrument before building (measurement; settles open Qs 1, 2, 3, 5).**
- Debug script: distributions of `τ_count`, `geom2`, `var_d`, and the realized count-prior **influence** on
  the gDNA suite; per-pass mass-delta trace (convergence/oscillation); count + μ-coverage of clean
  all-three-observable triplets; per-pass `var~mean` wall-time + peak per-chain edge memory.
- **Gate:** measured answers to — does the count lever survive `geom2`? do the loop deltas converge in ≤4
  passes without oscillation? do triplets span the exon μ-range? is refit sub-second? These decide Phases
  5–7's scope.

**Phase 2 — Cheap, unconditional correctness (Stage A).**
- `σ²_global = direct.predict(ρ_global)·geom2`; the Jensen `Δ_k` offset; the Bernoulli clamp on `σ²_frac`;
  promote the stop to config + add the prediction-error trace.
- **Gate:** zero-DNA scenarios clean (`test_ambig_scenario` phantom tol should *tighten*); gDNA-release
  benchmark net-flow ≥ baseline; no new oscillation. *(Verify the DIRECT curve has support near `ρ_global`
  before flipping — the latent `confidence(μ)` need.)*

**Phase 3 — Surface uncertainty (Step 1, unconditional).**
- `gdna_frac_var` (+ sides) onto `CalibrationResult`; thread through (diagnostic only).
- **Gate:** schema test; goldens regenerated; no behaviour change.

**Phase 4 — Diagnostics feature.**
- `CalibrationConfig.emit_diagnostics` / `--calibration-diagnostics` → `varmean_{gdna,rna}` / `per_region`
  (f_g, count_frac, μ_global, strand_class, mass, observability, `gdna_frac_var`, per-pass delta) /
  `global_scalars` dataframes (feather+TSV); companion `scripts/plot_calibration_diagnostics.py`. **No
  matplotlib in `rigel quant`.** Essential for real-data QC.
- **Gate:** dataframes render; the `var~mean` diagnostic visualizes the DIRECT/IMPUTATION split.

**Phase 5 — IMPUTATION re-cut (Stage B) — only if Phase-1 gates pass.**
- Re-cut IMPUTATION as the held-out boundary→region residual on the region-density axis (the user's NEW
  IDEA); else keep `raw_var~mean` + Jensen and record the decision.
- **Gate:** ss0.50 and capture-on net-flow improve; complex battery non-regressing.

**Phase 6 — Per-edge `q_rna` — only if Phase-1 shows propagation strength is a live lever.**
- Revive the RNA `var~mean`, derive `q_rna` per edge, surface as config; fix the RNA-strand double-count
  first (`rna_strand_bb_double_count_followup`).
- **Gate:** complex multi-exon loci improve over the `0.25` baseline against an isoform-aware oracle.

**Phase 7 — Precision-weight the per-locus prior (Step 2) — only if Phase-1 shows it moves accuracy.**
- The split-direction `c_locus` blend (§9). Else ship Step 1 only and document the prior as honest-output.
- **Gate:** net-flow at ss0.5 / AMBIG-heavy conditions improves; goldens regenerated.

---

## 14. Validation gates (every change)

- **Scoreboard** (flagship ss0.99 / ss0.50 / zero-DNA `nrna_none`+`nrna_rnd` / capture-off): **≥ baseline,
  no regression**; zero-DNA phantom ≈ baseline (must not balloon).
- **Complex battery** (`scripts/debug/complex_loci_benchmark.py`): TOTAL 2D/1D ≤ baseline; anchor-starved
  families degrade *gracefully* (toward `ρ_global`), not vertex-push.
- **Calibration-benchmark skill** (net fragment-flow, three-pool gDNA↔nRNA↔mRNA): leak ≤ baseline.
- **Convergence:** per-pass mass delta → 0 in ≤ `sweep_max_passes`; prediction error monotone-ish down; no
  oscillation.
- **Full suite** (`pytest tests/`) green; goldens unaffected until an intentional output change; node
  conservation `gdna + rna = total` holds.
- **Efficiency:** measured per-pass refit + per-chain peak memory within budget at genome scale.

---

## 15. Residual genuine open questions (measurement- or user-gated)

1. **Does the count lever survive `geom2`?** (`geom2 ≈ 14` at a moderate exon, `~5625` at `mass_u=20`.) If
   `τ_count` is near-inert at low coverage, Stage B (and even the Jensen fix) are invisible to net-flow.
   *Measure in Phase 1 before building Stage B.*
2. **Do clean all-three-observable triplets span the exon μ-range?** Stage B's held-out residual is honest
   only if such triplets exist at high μ. *Count + μ-coverage in Phase 1; if scarce/low-μ-only, keep
   `raw_var~mean` + Jensen.*
3. **Is the per-locus prior strength a live accuracy lever or an honesty feature?** Memory says flagship is
   EM-bound / isoform-swap-dominated. *Measure in Phase 1; if no gain, ship §9 Step 1 only.*
4. **Does per-edge `q_rna` beat the `0.25` scalar against an isoform-aware oracle?** Multi-transcript classes
   may be benchmark artifacts the EM resolves. *Phase-1 sensitivity sweep gates Phase 6.*
5. **Per-pass `var~mean` wall-time at genome scale.** Asserted sub-second; *measure in Phase 1.*
6. **The `confidence(μ)` guard / targeted-panel locality.** The all-gDNA init removes the need on
   whole-transcriptome data; a small **targeted panel** with enriched-exon means outside the training range
   could break flat-extrapolation. *Build a targeted-panel stress sim; add the guard only if it manufactures
   false certainty.* This is also the home of the unresolved SCAM span/locality task.
7. **Deferred refinement — same-class gDNA-coupling channel.** Re-open the "thick BP" only on a *measured*
   run-interior failure of one-hop imputation; it must be same-class-only and must not defeat the
   decoupled-edge `O(P)` short-circuit.
8. **Deferred refinement — unify the boundary-side solve** (retire `I₀` for the sides) as a one-time pre-loop
   3-tier solve, never per-pass.

---

## 16. Doc consolidation & status

- **Authoritative:** this plan; `effective_length_redesign_plan.md` §8 (eff-len, SHIPPED);
  `per_node_deconv_hierarchy_design.md` §3 (the per-node math) — update its stale §16 LOESS prose / β framing.
- **Active-design, fold in:** `var_fit_plan_evaluation.md` (the ROI/Jensen source); `new_var_fit_plan.md`
  (the SCAM span task → open Q 6).
- **Superseded → archive or banner:** `CALIBRATION_PLAN.md`, `CALIBRATION_PLAN_v2.md` (this replaces it),
  `iterative_bootstrap_design.md` (strand-only Pass-0 overruled), `simplex_ambig_count_fallback_design.md`
  (w-gated Jeffreys → node-class prior), `deconvolution_implementation.md` (its "no global ρ₀" thesis
  contradicts the shipped global tier), `new_eff_len_plan.md`, `capture_effective_length_design.md`,
  `fl_consistency_diagnostic.md`, `accumulator_fragment_span_redesign.md`, `stage2_wiring_dryrun.md`.
- **Scratch (reconcile then retire):** `calibration_nodes.md`, `my_notes.md` — their durable content is
  assimilated here (§3, §4).
- **Keep:** `splice_junction_node_architecture.md` (deferred bipartite SJ), `nascent_efflen_investigation.md`,
  `em_gdna_strand_likelihood_fix.md` (operative conclusion is its §9.4 unified-enrichment eff-len — an
  EM-side front, **not** part of calibration completion; its refuted §1-7 Option A/B must not be implemented).
- Rewrite the `docs/README.md` calibration section to point at v3.
