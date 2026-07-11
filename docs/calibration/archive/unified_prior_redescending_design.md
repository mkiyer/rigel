# Unified two-model gDNA prior via a DISTANCE-redescending background + a mixture — design derivation

**Status:** design derivation, 2026-07-06, branch (new) off `main` (post the M2 ship). This is the concrete
realization of the unified two-model vision (`unified_two_model_prior_WIP.md`) + the clean-slate architecture
(`gdna_prior_clean_slate_architecture.md`), using the user's key mechanism: a background shrinkage whose
**precision fades with the node's DISTANCE from the background mean**. It supersedes the `v_node`-scaled Cauchy
floor of `robust_redescending_prior_design.md` (which failed for the reason derived in §3).

## 1. The two jobs the prior must do (why M2 conflates them)
On an unstranded / sparse node the gDNA-vs-RNA split is set by the prior on `log ρ_g`. The prior must:
- **JOB 1 — spare enriched.** A capture-enriched node (ρ_g far above background) must NOT be dragged down.
- **JOB 2 — anchor background.** A depleted / no-gDNA node (ρ_g at background) must be held there consistently,
  so no-gDNA exons don't false-positive AND their abundance ranks stay stable.

M2's genome baseline does both with ONE knob: a Gaussian at `ρ_global` whose precision `1/(var_mean+σ²_g)` is
weak under capture (job 1) and strong at no-gDNA (job 2). That knob is the σ²_g density-dependence — the thing
that had to be a bistable spline. **The two jobs are really "weak far from background, strong at background" —
i.e. a precision that depends on DISTANCE from the background, not on a global capture scalar.**

## 2. The mechanism: distance-redescending shrinkage (derivation)
Let `r = log ρ_g − log ρ_bg` be a node's log-density distance from the background mean `ρ_bg`. A plain Gaussian
background prior `−½ r²/σ_bg²` has pull `−r/σ_bg²` (Hooke spring): **linear and unbounded → hardest on the
farthest (enriched) nodes** — the inverted influence (measured: pull 0.15 at the floor vs 7.2 on enriched,
`demo_floor_influence.py`). We want the opposite: strong at `r≈0`, vanishing for large `r`.

Use a **Student-t / Cauchy** background (a scale-mixture of Gaussians — the proper redescending prior):
```
glob(ρ_g) = −½ · log(1 + r² / σ_bg²)          r = log ρ_g − log ρ_bg
```
Its pull and EFFECTIVE precision:
```
pull(r)   = −r / (σ_bg² + r²)                  (peaks at r=σ_bg, → 0 as r→∞ : redescending)
N_eff(r)  = pull/r = 1 / (σ_bg² + r²)          (the "variance grows with distance": v_eff = σ_bg² + r²)
```
This IS the user's statement made exact — *"the amount of shrinkage / the variance depends on the distance from
the background mean; far nodes are barely affected."* `N_eff(0)=1/σ_bg²` (strong anchor), `N_eff→0` far out.

**The scale is σ_bg = the BACKGROUND spread, NOT the node's own variance `v_node`.** This is the correction that
makes it work:
- **σ_bg scale (this design):** an unstranded no-gDNA exon sits at `ρ_g≈ρ_bg` ⇒ `r≈0` ⇒ FULL precision
  `1/σ_bg²` ⇒ strongly anchored to the (near-zero) background ⇒ **JOB 2 solved** (no FP, stable ranks). An
  enriched node the data pushed to `ρ_g≫ρ_bg` ⇒ large `r` ⇒ released ⇒ **JOB 1 solved**.
- **`v_node` scale (the disqualified 2026-07-04 attempt):** `v_node` is LARGE for unstranded nodes (no strand
  info), so `N_eff(0)=1/v_node` is tiny ⇒ unstranded no-gDNA exons were only WEAKLY anchored ⇒ they drifted ⇒
  the 106k no-gDNA FP. The scale was the bug; distance-from-mean is the fix.

`σ_bg` is data-derived (the pooled intergenic+intron log-density spread, `_floor_estimate`) — no magic number.

## 3. Why redescending does NOT reintroduce the no-gDNA FP (the disqualification, resolved)
The worry: "if the pull vanishes far away, a drifted no-gDNA node won't be pulled back → FP." Resolution, two
independent guards:
1. **No-gDNA nodes don't drift far.** A no-gDNA exon's gDNA count ≈ 0, so at any candidate `f_g` its implied
   `ρ_g=f_g·M/E` is only large if `f_g` is large — and nothing in its data (flat strand) or the RNA Jeffreys
   `−log(1−f_g)` favours large `f_g`. The redescending region (large `r`) is only REACHED by nodes whose data
   pushes them there = genuinely dense (enriched) nodes. So redescending engages for real enriched nodes, not
   for no-gDNA drift. (This is why the σ_bg-scaled Cauchy is safe where the v_node one was not: it anchors the
   no-gDNA node's low-f_g region at full strength.)
2. **The Cauchy tail is still monotone** (`½log(1+r²/σ_bg²)` increases without bound, just logarithmically), so
   an implausibly-high `ρ_g` is still penalised — weakly, but combined with the Jeffreys it suppresses FP. And
   in PASS 2 the mixture's outermost mode carries a real quadratic tail (§5).

## 4. The two models (clean split)
- **Background model** `B` — from **intergenic + intron REGIONS only** (`_floor_estimate`): mean `ρ_bg`, spread
  `σ_bg`. Adaptive: `ρ_bg≈0` in no-gDNA, elevated-uniform in gdna300/capoff, depleted-off-target in
  gdna300/capon. Constant across both passes. This is the redescending shrinkage centre/scale of §2.
- **Enriched model** `E` — a KDE from **exon + boundary nodes only** (the substrate NOT used for `B`, so no
  double-use). Trained on the pass-1 solved densities. If capture: a second mode forms above `ρ_bg`. If no
  capture: it collapses onto `B` (enriched ≈ background — the user's requirement). Pass-2 only.

## 5. The passes
**PASS 1 — redescending background only, on ALL nodes.** Prior = the §2 Cauchy at (`ρ_bg`, `σ_bg`). Effect:
depleted / no-gDNA / unstranded nodes anchor to `ρ_bg` (job 2, clean pass-1 densities); nodes whose data drives
them dense (STRANDED enriched) are released (job 1) and solve to their true density → the enriched teachers.
Then fit `E` = KDE on the solved exon+boundary densities.

**PASS 2 — the mixture, on ALL nodes.** Prior on `log ρ_g` = the normalized mixture
```
P(ρ_g) = w_B · N(log ρ_bg, σ_bg)  +  w_E · E(log ρ_g)          (evaluated via logpdf_kernel: real tails)
```
Each node "shrinks toward the nearest mode" — this is the mixture's native GRAVITY, realized by the responsibility
`r_B(x)=w_B N_B/(w_B N_B + w_E E)`: near `ρ_bg` the background dominates (anchor); near the enriched mode the
enriched dominates (spare); the OUTERMOST mode's real quadratic tail springs back a node above everything (FP
suppression — the guard the redescending-floor lacked). No capture ⇒ `E≈B` ⇒ unimodal ⇒ ≈ pass-1. Capture ⇒
bimodal ⇒ enriched nodes resolve to the enriched mode.
- Weights `w_B, w_E` = the background vs enriched **population evidence** (node counts / pooled exposure) —
  data-derived, magic-free. (Sensitivity to be measured; §7.)
- **Boundaries** are scored at the crossing-density scale (below the contained enriched mode) — either give them
  a crossing-scale enriched reference (KDE_boundary_prior_review Direction A) or keep them on `B` only; do not
  let the contained-scale enriched mode crush them (the 86k leak).

## 6. The honest limit (unstranded capture)
Without strand, an enriched exon's high `M/E` is gDNA-vs-RNA **unidentifiable** at the node. Pass-1 anchors such
a node to `ρ_bg` (no data to release it), so it does NOT self-identify as enriched — the enriched mode under
FULLY-unstranded (ss0.50) capture must be learned from the **density bimodality** and the **structural**
(boundary/spliced) RNA signal, not from strand. Expect: strong wins on stranded + all no-gDNA (job 2); partial
improvement on unstranded-capon (the M2 leak regime) — bounded by identifiability, not by the prior. This is the
metric to watch; matching M2 there while fixing no-gDNA anchoring is already a net win.

## 7. Implementation (env-gated A/B, deterministic)
Reuse existing machinery; one new toggle `RIGEL_PRIOR ∈ {baseline(M2), unified}`.
1. **Background** (`bp_solver`): `_floor_estimate` already gives `ρ_bg`(=ρ_floor over intergenic+intron) +
   `σ_bg`(=√(s2_floor+var_mean_floor)). Add `_redescending_background(fgg, mass, eff, ρ_bg, σ_bg)` returning the
   §2 Cauchy term on the f_g grid. Pass-1 (`gdna_prior=None`, unified): global_lp = this, on ALL nodes (drop the
   M2 genome baseline + the floor_mask special-case).
2. **Enriched KDE** (`gdna_density_prior`): `build_training_substrate` restricted to exon+boundary (hard split,
   not the soft weight); AMBIG excluded. Fit as today (deterministic bandwidth).
3. **Mixture** (`bp_solver._kde_logprior` → a mixture term): pass-2 global_lp = logpdf_kernel of
   [background kernel (μ=log ρ_bg, h=σ_bg, w=w_B) ⊕ enriched teacher kernels (w summing to w_E)] + Jeffreys.
   (Per-sample bandwidth already prototyped in the WIP patch.) NO M1 −log f_g Jacobian (measured catastrophic).
4. **Determinism:** every input (ρ_bg, σ_bg, weights, KDE bandwidth) is a continuous data function; keep the
   interpolated `_weighted_median` (shipped). No argmin. The M2 line stays as the `baseline`/fallback.
5. **Validate:** 16-scenario net-flow (primary) + abundance, `unified` vs `baseline(M2)`; the 3 decisive cells —
   no-gDNA-unstranded (job 2 must hold ranks + FP), gdna300-capon-unstranded (job 1 / the leak, vs M2), and a
   boundary-leak probe. Cross-process determinism check. Then goldens + full suite.

## 9. RESULT — REFUTED by the benchmark (2026-07-06). M2 remains best.
Implemented (env `RIGEL_PRIOR=unified`) + measured on the 4 decisive cells (baseline=M2):

| cell | M2 spear/leak | unified spear/leak |
|---|---|---|
| none ss0.50 capon (job 2) | 0.791 / 6997 | **0.578 / 24610** |
| none ss0.50 rnd capon | 0.673 / 60503 | 0.492 / 37527 |
| gdna300 ss0.50 capon (job 1) | 0.714 / 71623 | 0.663 / **19214** |
| gdna300 ss0.99 capon (stranded) | 0.710 / 26911 | 0.726 / 28014 |

**The §3 guard was wrong.** For an UNSTRANDED node the strand is flat, so the PRIOR ALONE sets f_g — and a
redescending prior is (by construction) weak at large r. At f_g=0.5 the σ_bg-Cauchy penalizes an implausible
ρ_g by only ~2.7 nats vs M2's capped-constant ~112 nats, so an unstranded no-gDNA exon DRIFTS up → FP (job 2
BROKEN, leak 7k→24.6k). Those drifted pass-1 exons then pollute the enriched KDE (spurious high mode) → pass-2
FP. This is the same failure that disqualified the v_node floor; the σ_bg scale makes it milder (24.6k vs 106k),
not gone.

**The deeper reason (identifiability):** the g300 leak improvement (71.6k→19.2k, real — enriched released) and
the no-gDNA FP regression (7k→24.6k) are ONE mechanism — the redescending releases every far node, real-enriched
and no-gDNA-drift alike. Under unstranded data there is no per-node signal to distinguish them, so you cannot
keep the enriched release without the no-gDNA FP. M2's choice — a GLOBAL density-dependent precision that stays
strong-everywhere in no-gDNA (suppress FP) and weak-everywhere under capture (spare enriched) — is the better
trade than a per-node distance-dependence, precisely because "strong-everywhere in no-gDNA" is what a redescending
(weak-far) prior cannot provide. Strand-filtering the enriched teachers would recover job 2 but then collapses
to background-only under ss0.50 (no strand → no enriched mode) → ≈ M2. Net: the redescending does not beat M2.

**Disposition:** keep the shipped M2. This design + result documented; the identifiability wall on unstranded
capture is the real limit (needs a NON-prior signal — structural boundary/spliced RNA, or the mappability/GC
covariates — not a cleverer prior). Env toggle left in place for future study.

## 8. Why this can beat M2 where the earlier attempts couldn't
- vs **M2**: replaces the global σ²_g density-dependence with a PER-NODE distance-dependence — job 1 & job 2 are
  decoupled by construction, no conflated knob.
- vs the **`dropall`/seed unified attempt**: pass-1 no longer leaves exons unanchored — the σ_bg-scaled
  redescending anchors near-background nodes (fixing the spearman regression) while still releasing far ones.
- vs the **`v_node` Cauchy floor**: the scale is the background spread, so unstranded no-gDNA is strongly (not
  weakly) anchored — the fix for the 106k FP disqualification.
