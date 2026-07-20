# Task for an independent agent: solve the gDNA-hyperprior *projection* (the "snap-to-depleted" problem)

> **⚠ SUPERSEDED (2026-07-19).** This external-handoff prompt is kept for provenance only. Two things changed
> after it was written: **(1)** the **fitted RNA arm (S4 / `logP_r`) was tested and REJECTED** — a full-solve
> A/B showed it does not recover under message anchoring and manufactures false positives; ignore every S4 /
> "RNA-side prior" ask below. **(2)** The root-cause investigation (see `kde_vs_npmle_enriched_mode` and
> `CALIBRATION_MASTER.md`) found the snap is **structural to the NPMLE representation** — a competitive EM
> mixture with a heavy in-mixture background cell (weight = `n_regions`) starves the enriched mode, whereas the
> shipped v0.7.1 KDE (additive unit-weight kernels + a separate *capped* product-of-experts background) keeps
> it. The live plan is to **debug the gDNA-hyperprior representation toward the proven KDE design** (bandwidth
> S1/S2 and occupancy-leveling S3 remain relevant; S4 is closed). Read `CALIBRATION_MASTER.md` for the current
> roadmap; the sections below are historical.

## Your mission

In Rigel's calibration, a **gDNA-rate hyperprior** (a 1-D density `π_g(log ρ_g)` estimated by a fixed-kernel
Poisson-mixture NPMLE) is **projected onto each genomic node's composition axis** `f_g` (the gDNA fraction of
that node's fragments). Under hybrid-capture the hyperprior is **bimodal** — a **tall "depleted" (background)
mode** plus **short "enriched" (capture-enriched gDNA) modes**. The projection currently **snaps
partially-enriched nodes onto the tall depleted mode**, so their `f_g` is driven far too low and real gDNA is
mis-assigned to RNA. We need a **principled projection** that recovers the weak enriched modes **without
manufacturing false-positive gDNA** where there is none. This is the single biggest remaining error in the
calibration; solving it correctly should yield a large accuracy gain.

We have already derived the correct statistical model and narrowed the levers (below). We want you to
**pressure-test that derivation, then help design + validate the concrete fix** — especially the two pieces our
own analysis says are load-bearing but under-specified: (S3) rebuilding the enriched mode's *weight*, and (S4) a
proper RNA-side prior. We work on controlled simulations now; real capture libraries come later.

---

## 1. Orientation — what the calibration does (just enough)

Rigel jointly quantifies mRNA, nascent RNA, and genomic-DNA (gDNA) contamination from RNA-seq. The genome is
partitioned into a **chain of nodes** (regions between transcript features, and boundaries = splice junctions).
For each node we solve the composition of its **unspliced** fragment mass on the axis `λ = logit(f_g)`, where
`f_g = P(a node's unspliced fragment is gDNA)`. The solve is a belief-propagation pass over the chain; each
node's posterior over `λ` combines a **strand likelihood** (a Beta-Binomial on the ± strand counts — the only
*intrinsic* composition signal), **cross-node messages**, and a **gDNA prior**.

The **hard regime** is **unstranded + hybrid-capture**: the strand tilt's Fisher information is `∝ N(2κ−1)²`,
which is **identically 0 at κ=½ (unstranded)**. So an unstranded node has *no intrinsic evidence* about its own
composition — the gDNA prior and the messages must carry it. Capture *enriches gDNA in the targeted exons*, so
the gDNA-density distribution is bimodal (un-enriched background vs enriched-in-capture).

**Pipeline order (all implemented):** initial *prior-free* solve → fit the gDNA hyperprior on the solved gDNA →
**refit** the solve with the hyperprior as the composition prior → project to per-locus priors. The snap
happens in the **refit**: the hyperprior's tall depleted mode pulls enriched nodes down.

## 2. The invariant you must respect — count-zero-information

**A fragment *count* carries zero intrinsic information about gDNA-vs-RNA composition.** A node with 1000
fragments is not "more gDNA" than one with 10. The count enters *only* as precision (how sharply we know a
rate). Concretely: **the hyperprior `π_g` is a prior over the *absolute* gDNA rate `ρ_g`, not over composition
`f_g`.** Any mechanism that lets a node's *total density* `M/E` vote its *composition* is the cardinal bug we
keep having to root out. (See `npmle_is_enrichment_not_dna_prior` in the memory notes: total-density is an
*enrichment* signal, never a DNA-composition prior.)

## 3. The problem, with concrete numbers (sim condition `gdna5_capON`, ~5% gDNA, unstranded, capture-on)

The fitted hyperprior `π_g(log₁₀ ρ_g)`, heights normalized to the tallest mode = 1:

```
depleted mode  log10 ρ_g = −3.06   height 1.000   ← "Everest"
enriched modes         −2.01        0.176         ← the partially-captured nodes live here
                       −1.64        0.148
                       −0.18        0.112
                        1.62        0.003          ← 300× shorter
```

There are **35 region nodes** whose observed total density sits near the −2.01 enriched mode and whose **oracle
`f_g = 0.63`** (they are genuinely ~⅔ gDNA). The **prior-free initial solve** gets them ≈ right (`f_g = 0.507`).
The **refit snaps them to `f_g = 0.176`** (mostly RNA). Trace of one node (2589, total density 10^−2.31, oracle
`f_g = 0.33`): the gDNA arm `logP_g(f_g)` peaks at `f_g = 0.18` (mapping `ρ_g` onto the depleted mode −2.97,
`logP = −3.12`); every enriched-mapping `f_g` is **~4 log-units lower** (−6.5…−7.6). The depleted mode
out-votes the enriched by `e⁴ ≈ 50×`. That is the snap.

**Boundary nodes are frequently *partially* enriched** (a capture probe only partly covers a region edge), so
"which mode did this come from?" genuinely has a **mixture** answer — we must integrate over the populations
(pushing toward the most likely one), not hard-clamp to one.

## 4. The correct model (our derivation — please verify or correct it)

Node observables: unspliced mass `M` (≈ count), effective length `E`, total density `ρ_tot = M/E`, strand split
`u = (u₊,u₋)`. Latent **absolute rates** `ρ_g, ρ_r` (gDNA, RNA); `a = log ρ_g`, `b = log ρ_r`,
`λ = a − b = logit f_g`. gDNA-enrichment and RNA-expression are independent ⇒ the prior is the product
`π_g(a)·π_r(b)`. The total count is `M | ρ ~ Poisson((ρ_g+ρ_r)E)`; the strand split is the Beta-Binomial
`BB(u | σ(λ))` (flat in `f_g` at κ=½). Marginalizing the total direction, the **correct composition posterior**
is

```
P(λ | M,u) ∝ BB(u | σ(λ)) · ∫da  exp(−½·M·(a − a*(λ))²) · π_g(a) · π_r(a − λ),   a*(λ) = log(f_g · M/E)
```

i.e. **`π_g` convolved along the log-rate axis by the node's own counting kernel (precision `M`, SD `1/√M`),
read along the composition ray, jointly weighted by a real RNA prior `π_r` on the RNA the composition implies.**

The **shipped** solver instead computes `ψ(λ) = BB(u|σ(λ)) + logπ_g(log(f_g·M/E)) + ½·log(1−f_g)` — which is the
above in the **degenerate `M→∞` (δ-pin) × `π_r ≡ flat` corner**: (i) the counting kernel collapses to a point
evaluation on the ray (the node's total is treated as *known exactly*), and (ii) the RNA arm is an improper
flat reference that never penalizes "explain the whole node as RNA." With strand silent, `ψ` is then just the
*height* of `π_g` along the ray, whose argmax is the **tallest sub-`ρ_tot` mode** = the depleted Everest.
(Note: this is also a count-zero-info violation — the count-pinned total selects which mode votes.)

**Four surfaces a fix can act on:** **S1** the estimation bandwidth `h` of the NPMLE fit; **S2** the per-node
counting kernel `1/M` (dropped by the δ-pin); **S3** the fitted *weight/height ratio* of the enriched vs
depleted modes; **S4** the RNA arm `π_r` (currently flat).

## 5. What we've already tried (empirical, on the sim)

- **Bandwidth sweep (fit `h`, S1), full re-solve:** raising `h` recovers the 35 nodes only `0.176 → 0.234`
  (target 0.63), and **over-smooths past `h≈0.65`** (modes merge, recovery drops). Not sufficient alone.
- **Per-node counting kernel (S2)** re-rendering `π_g` at `H² = (h·ln10)² + 1/M`: helps only `0.234 → 0.297`
  here, because these nodes have `M≈10–55` so `1/M` barely widens `h²` (it "bites" only at `M ≲ 8`).
- **Widening the *projection* kernel (render-only):** reaches ~0.39 — still far from 0.63 — and the **full
  refit loop drags recovery *below* the render-only value** (the circularity reinforces the snap).
- **Uniform interior floor** (`floor_eps=0.02`, already present): fixes `−∞` valleys but does **nothing** to
  the *relative mode heights*, which is the actual failure.

**Conclusion so far: bandwidth is necessary (it de-approximates S1/S2 exactly) but insufficient on this case.**
The residual **enriched-mode weight (S3)** still out-votes, and/or the missing **RNA discriminant (S4)** is what
would actually penalize the all-RNA escape. Slogan from the derivation: *"h buys margin; π_r buys the
discriminant."*

## 6. Candidate levers and our current verdicts (challenge these)

1. **Bandwidth (S1 data-driven `h` + S2 counting kernel) — primary, ship-track.** The only term-for-term
   de-approximation of the yardstick. Proposed data-driven default `h = median(1/√g)/ln10` (Poisson log-SD
   floor). *Under-specified:* the right selector, and whether an adaptive/variable `h` is needed.
2. **Enrichment-skeleton (S3) — conditional, we think required here.** Bandwidth cannot rebuild an enriched
   mode whose *weight* pass-0 collapsed. Idea: detect the enriched mode's *existence/location/occupancy* from
   the **independent total-density (enrichment) landscape** — which we already fit for transfer-variance — and
   inject it into `π_g` **anchored at the absolute rate `ρ_bg · e^{δ_k}`** (background rate × detected
   enrichment), with an *occupancy* (node-count) weight, never fragment mass. *Open:* detectability limit (can
   35 enriched gDNA nodes form a total-landscape mode distinct from the RNA-expression continuum?), and how to
   avoid RNA-abundance *smearing* the composition.
3. **RNA arm `π_r` (S4) — parallel track, not one of the four "levers."** Replace the flat `½·log(1−f_g)` with
   a proper (symmetric) `logP_r` so the all-RNA escape is penalized. Our flagship diagnostics independently
   blame the flat/asymmetric RNA arm. *Open:* the right form of `π_r` that stays count-zero-info-clean.
4. **Peak-height tempering `π_g^{1/T}` — disfavored** (not derived from the model; a heuristic reshaping).
5. **Mixture-integration projection (weight modes by proximity of the node's *total* to each mode) — REJECTED:**
   its proximity factor reads the total to vote composition = the count-zero-info violation.

**FP-bounding principle we hold to:** a lever may only reshape the *width/weight/smoothing* of `π_g` (which is
fit on **deconvolved** gDNA, never the total), via operations that are **monotone/support-preserving OR anchored
to the absolute background rate `ρ_bg`** — so it can *lower* an over-dominant mode but **never create** one.
(Structural guarantee: the fit's grid top is `≤ ρ_tot`, so `ρ_g = f_g·M/E ≤ M/E` can't extrapolate above the
total — a true zero-gDNA node has only the sub-total background mode on its ray, giving `f_g ≈ 0`.)

## 7. The open questions we want your help on

1. **Is the yardstick (§4) right?** Any error in the generative model, the counting-kernel derivation, or the
   claim that the shipped `ψ` is exactly its `M→∞ × flat-π_r` corner?
2. **The enrichment-skeleton (S3):** what is the *principled* way to detect + place + weight an enriched gDNA
   mode from the total-density landscape **without** letting RNA-expression levels smear it or manufacture a
   phantom mode? Is there a clean separation of *detection* ("an enriched population exists") from *composition*
   ("this node's `f_g`")? What is the **fundamental limit of detection** — when gDNA is a few % and/or the
   capture panel is small, can the enriched gDNA mode be identified at all from this data?
3. **The RNA arm `π_r` (S4):** what is the right count-zero-info-clean form? (It must not itself become a
   density-votes-composition backdoor.)
4. **Bandwidth (S1/S2):** the right *data-driven* selector for a **Poisson-rate mixture** NPMLE (Silverman /
   Sheather-Jones analogues? likelihood cross-validation? an information-driven `h(1/M)`?), and the diagnostic
   that detects over-smoothing (mode merge). Real capture chemistry will set final widths later — we want the
   *principle* + a sim-testable default now.
5. **How do S1–S4 compose**, and what single principle bounds the sensitivity/specificity trade (recover weak
   enriched modes without amplifying false positives)?
6. **The refit circularity:** the hyperprior is fit on the solve it then corrects; a pass-0 snap starves the
   enriched weight. Is there a non-circular fit/refit scheme (or a damping) that avoids self-reinforcement?

## 8. Constraints & acceptance criteria

- **Constraints:** respect count-zero-information (§2); **no unjustified magic numbers** (pause before any new
  tuned constant); do **not** perturb the enrichment/transfer-variance NPMLE (its bandwidth is signed-off at
  0.15) — changes apply to the **gDNA composition hyperprior only**; settle what you can on sims, flag what
  must wait for real data (and note tumor **CNV** breaks the `ρ_g = ρ_bg·e` uniformity the skeleton assumes).
- **Acceptance gates (sim, `ambig_dense_10mb` suite):**
  - **Recovery:** the 35 `gdna5_capON` enriched nodes → median `f_g ≈ 0.63` (from 0.176).
  - **Specificity (must pass):** `gdna_none` (true zero gDNA) → median `f_g < ~0.02` (no manufactured gDNA).
  - **No regression:** stranded (`ss0.99`) and resolvable conditions stay correct.
  - Metric for the fitted-density comparison to oracle: **Wasserstein-1 on the log-ρ axis** (L1 saturates on
    the disjoint-support zero-gDNA case).

## 9. Docs and code to read

- `docs/calibration/gdna_projection_derivation.md` — the full derivation (yardstick + 4-lever analysis + ranked
  experiment plan Exp 1–7). **Start here.**
- `docs/calibration/gdna_background_floor_derivation.md` — the just-landed background-floor fix (the depleted
  mode's *location*), a solved sibling problem and a template for a workflow-derived + live-validated fix.
- `docs/calibration/CALIBRATION_ARCHITECTURE.md` — the count-zero-information principle (authoritative).
- `src/rigel/calibration/npmle.py` — `DensityNPMLE` (the fixed-kernel Poisson-mixture NPMLE): `.fit`
  (bandwidth `h=0.15`, `floor_eps`, grid, aggregate-background cell), `.logprior` (the δ-pin projection — the
  thing to fix), `.project` (the softmax responsibility — do **not** touch; it's the transfer-variance path).
- `src/rigel/calibration/simplex_logodds.py` — the per-node ψ solve; `_gdna_arm` and line ~267
  `psi = psi + _gdna_arm(lam, global_logprior) + _rna_arm(lam)` (the flat RNA arm is `_rna_arm`).
- `src/rigel/calibration/calibrate.py` — `_fit_gdna_hyperprior` (the training-node selection) and the Phase-2
  refit loop (`calib_refit_iters`).
- `src/rigel/calibration/background_reference.py` — `measure_background` (the `ρ_bg` scalar + the derived
  resolution floor).
- Prototyping harnesses (in the scratchpad; reproduce our numbers): `bandwidth_sweep.py`, `exp1_rescore.py`,
  `refit_gdna_movie.py` (the frame-by-frame pass-0 → refit vs oracle, EMD metric), `projection_diag.py`.

Please critique the derivation, then propose and (where possible) prototype the fix for S3 and S4 against the
§8 gates. Be adversarial about false-positive amplification and about any mechanism that reads composition from
total density.
