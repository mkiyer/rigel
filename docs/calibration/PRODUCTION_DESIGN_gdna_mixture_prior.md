# Production design — single-strand bootstrap → Bayesian gDNA mixture prior → AMBIG solve

**Status:** architecture design / handoff, 2026-06-30. This is the **strategic arc** the session converged on.
The companion `NEXT_SESSION_calibration_tasks.md` is the *immediate, tactical* work (line 1: the single-strand
solver). THIS doc is line 2: the full Bayesian model that consumes the single-strand solutions and prices the
ambiguous regions. Read both. Background memory: `calib_mature_message_and_toy_prod_driver.md`,
`disagreement_precision_shipped_flagship_ehat_bound.md`.

---

## 1. The vision — a three-phase calibrator

The calibrator's one job is to hand the per-locus EM an **accurate gDNA prior** (expected gDNA fragment count
per locus). We have repeatedly shown that *given an accurate gDNA prior, the EM solves abundance beautifully* —
so calibration is the whole game, and the gDNA prior is the deliverable. The architecture:

```
Phase 1  SOLVE single-strand nodes   (strand=NONE/POS/NEG; AMBIG held out)
         using structure + strand-tilt + message-propagation + an EXTREMELY WEAK gDNA floor (stability only).
         These nodes are solvable WITHOUT a gDNA enrichment prior. → a set of trustworthy per-node gDNA densities.
                                   │
                                   ▼
Phase 2  TRAIN a Bayesian model of the gDNA-density distribution P(ρ_g) from those single-strand solutions.
         Capture makes this a MIXTURE (depleted off-target ⊕ enriched on-target; CNV adds more modes).
         Model it NONPARAMETRICALLY (bandwidth-smoothed empirical density in log space) — no hard clustering,
         no mode-counting, no capture detector. The depleted MINIMUM (intergenic) anchors the floor.
                                   │
                                   ▼
Phase 3  SOLVE ambiguous (multi-strand / AMBIG) nodes by SOFT MARGINALIZATION over that mixture prior,
         fused with the node's own (weak) strand-tilt + messages. The prior fills exactly the tilt's null space.
                                   │
                                   ▼
         Project to per-locus gDNA prior counts → the EM.
```

**Why this order:** Phase 1 needs no gDNA prior (single-strand nodes self-resolve from strand, or from messages
+ the spliced-residual when strand is weak), so it is non-circular. Phase 2's prior is bootstrapped only from
the trustworthy Phase-1 nodes (the non-circular teacher). Phase 3's hard cases (balanced-tilt AMBIG under
enrichment) are exactly where the strand likelihood is flat in f_g, so the Phase-2 prior IS the answer there.

---

## 2. Phase 1 — the single-strand solver (the foundation; detail in NEXT_SESSION doc)

PRIORITY #1, because everything downstream is bootstrapped from it. It must be pristine **stranded AND
unstranded**. It is the existing BP solver (per-node log-odds solve + the forward-backward message sweep +
disagreement-aware message precision, already shipped) restricted to single-strand nodes, with only a weak
floor. The open work (line 1): (A) mature-projection eff-length geometry, (B) the spliced→mature MEASUREMENT
message (implemented, WIP on branch), (C) the nascent-hallucination constraint (§5). When strand is strong the
tilt pins each node; when strand is weak/unstranded, the message system (the generalized "spliced residual")
carries it — which is why the message PRECISION model is the load-bearing piece.

**Teacher definition is κ-dependent — define by CONFIDENCE, not strand-class:**
- High κ → strand-solved single-strand nodes are the teacher (measured perfect on toys).
- κ→½ (unstranded) → the strand carries zero gDNA info; the teacher becomes the *structurally* solvable nodes:
  intergenic (the depleted floor) and the **spliced-residual** (motif-stranded → strand-free → works at any κ).

---

## 3. Phase 2 — the Bayesian gDNA mixture prior (the new model, the meat)

### 3.1 The object
A distribution `P(log ρ_g)` over per-node gDNA densities, **estimated nonparametrically** (a kernel-density /
smoothed-histogram, in log space) from the Phase-1 single-strand solved gDNA densities. This single object is:
- **unimodal** off-capture (uniform gDNA),
- **bimodal** under capture (depleted off-target ⊕ enriched on-target; the enriched mode is *wider* — partial
  + full probe overlap, probe-efficiency variation),
- **multimodal** under intrinsic CNV (amplification/deletion copy-number levels) —
all from the **same code**, with no "is this capture?" detector and no "how many modes?" decision. That
robustness is the entire reason for the nonparametric choice (see §6, ruled-out alternatives).

### 3.2 Why nonparametric, not a Gaussian mixture / hard clustering
The user's explicit constraint: avoid brittle unsupervised clustering and untunable hyperparameters. A GMM
needs K (component count), initialization, label-switching handling — all brittle, and a misclassified node
*snaps* to the wrong mode (a cliff). The nonparametric density has **one** knob (bandwidth) and never classifies
hard: the per-node solve marginalizes softly over the whole density, so an ambiguous node gets an honest
(possibly bimodal) posterior that its own evidence (tilt / messages / spliced-residual) then resolves smoothly.

### 3.3 The bandwidth (the one knob — needs a robust estimator + a plotting framework)
Build, as deliverables: (a) several bandwidth estimators (Silverman's rule, GCV, likelihood cross-validation —
reuse the existing `variance_model` GCV machinery where possible); (b) a **plotting framework** to *see* the
fitted P(log ρ_g) and how each estimator behaves (number/separation of modes vs bandwidth). Decide the
production estimator empirically — on **real data** when it lands; what works on synthetic may not transfer.
Until then, expose the bandwidth + the plots so we can inspect it on every scenario. Do NOT hard-code a magic
bandwidth.

### 3.4 The floor = the truly-depleted MINIMUM, not the population mean
The shrinkage target / low anchor is the **intergenic** gDNA density (the off-target depleted minimum), a
*minimum-class* baseline — NOT the exposure-pooled mean (which averages depleted + enriched and sits above the
floor). This is the "global gDNA minimum, not mean" point. Introns approximate it but carry nascent risk; the
purest depleted class is intergenic. The mixture's depleted mode should sit at this floor.

### 3.5 The transferable quantity is the MODE LOCATIONS
The single-strand→AMBIG transfer is robust because what transfers is the **set of mode locations** (the depleted
and enriched gDNA-density LEVELS), which are a property of the capture probes / CNV, *shared across strand
class* — NOT a per-node regression slope (which the red-team showed is fragile; §6). An AMBIG node's own
evidence decides *which* mode it belongs to; the model only supplies *where the modes are*.

---

## 4. Phase 3 — solving AMBIG by soft marginalization

The per-node grid solve is **already** a correct product-of-experts over `(λ=logit f_g, τ=RNA tilt)`. To use the
mixture prior, evaluate `P(log ρ_g)` (Phase 2) at each grid point and add it as the per-node prior term
(`logsumexp` over the density — trivial, no new solver machinery). Then the per-node posterior fuses:
- the **strand-tilt likelihood** (Beta-Binomial; informative when imbalanced, flat in f_g when balanced),
- the **messages** (gDNA + nascent imputation at disagreement precision; mature measurement at count precision),
- the **mixture prior** (Phase 2).

**The prior's strength self-scales — no knob, no detector.** The tilt's Fisher information in the f_g direction
is `∝ N(2κ−1)²`; it → 0 at a balanced ridge AND at κ→½ AND off-capture. So the prior dominates *exactly* where
the node is uninformative (balanced AMBIG under capture) and recedes where the tilt pins f_g (single-strand).
This is the principled resolution of the "strong vs weak prior" tension: Bayes does it, not a tuned weight.

**The spliced-residual's role (red-team-corrected):** for an AMBIG node it is a per-node, strand-free,
gDNA-specific **observation** of ρ_g (`clip(M_unsp − ρ_mature·E_rna, 0)/E_gdna`), updated against the mixture
prior — a clean Bayesian shrinkage toward the nearest mode. It localizes even an *isolated* whole-AMBIG locus
(no single-strand neighbour) within the prior. It is NOT a predictor fed to a regression (that path is dead, §6).

**Multimodal-posterior → EM:** an AMBIG node with no resolving evidence yields a genuinely bimodal posterior
(honest uncertainty). Decide how to summarize it into the per-locus EM Dirichlet scalar — usually the messages +
spliced-residual resolve the bimodality first; otherwise propagate the distribution or take a principled point
estimate. (Open question, §7.)

---

## 5. The nascent-RNA constraint (hard rule — threads through Phases 1 & 3)

**Never hallucinate nascent RNA from a boundary↔exon discrepancy.** Definitive nascent evidence comes ONLY from
a **single-stranded intron** showing either (1) fragment density **much higher than the intergenic gDNA floor**
(the excess = nascent), or (2) **strand sense-bias** inside the intron. Intronic fragments *support* nascent;
messages must not *create* it. Scant hallucinated nascent is tolerable; gross is not. This is why the floor must
be the depleted *minimum* (the comparison baseline) and why the mature message routes to exons only. In the
mixture-prior world, an intron's nascent component is gated by the intron's own density-excess / strand — not by
imputation from expressed exon neighbours.

---

## 6. RULED OUT (do not re-derive — these were explored and rejected this arc)

- **A regression `ê(z)=E[ρ_g|z]` as the AMBIG prior mean.** Red-teamed out: (1) **predictor≡response
  degeneracy** — if the predictor is itself a gDNA estimate (e.g. the spliced-residual) and the teacher response
  is also gDNA, the fit collapses to identity and cannot de-bias; (2) the boundary-flux predictor `z`
  **saturates** at high enrichment; (3) **regression dilution** + **flat-clamp extrapolation** cap the lift
  below truth. The mixture-prior + per-node-observation framing *replaces* this — no regression of gDNA on gDNA.
- **Contained density `M/E` as a gDNA proxy.** Ratio-inflated (= gDNA + RNA); `corr(contained,gDNA)` decays from
  0.998 to 0.90 as gDNA:RNA falls 6.8→0.07. Use the gDNA-specific spliced-residual, not contained density.
- **Hard clustering / explicit K-component GMM with mode classification.** Brittle, cliff-prone, untunable K.
- **Intergenic *mean* as the floor.** Use the depleted *minimum*.
- **Trusting the hand-rolled toy harness or the hand-assembled `node_sweep` dives.** Proven unfaithful (flipped
  stranded RNA→gDNA). Always validate via the real `calibrate()` (`scripts/debug/toy_prod.py`). The earlier
  flagship node-dive numbers (0.533 anchor, 77/23 split) are DIRECTIONAL only.

---

## 7. Open design questions to resolve while building Phase 2/3

1. **Bandwidth estimator** — which rule, and is it stable across the gDNA:RNA × κ × capture ladder? (Build the
   plotting framework first; decide on real data.)
2. **Multimodal-posterior → per-locus Dirichlet scalar** — propagate the distribution vs point-estimate; what to
   do for a genuinely-uncertain isolated AMBIG node.
3. **Whole-locus-AMBIG capture loci** — no local single-strand teacher, uninformative messages → rely on the
   spliced-residual + the marginal prior. Is that enough, or is a covariate (probe-overlap) fallback needed?
   (Decide empirically; prefer NOT adding a covariate if the residual suffices.)
4. **Ordering** — current hypothesis: a single multi-stage pass (solve single-strand → train prior → solve
   AMBIG), no iteration (the human transcriptome has enough single-strand nodes to train once). Verify
   convergence is unnecessary on real data.
5. **CNV / intrinsic copy-number** — the nonparametric prior captures extra modes; the spatial messages
   (disagreement-aware, breakpoint-respecting) localize a node to its segment. Worth mining the copy-number
   segmentation literature (the modes = copy-number levels) if intrinsic CNV becomes a target.

---

## 8. Literature grounding (from the session's hypothesis-generation workflow)

The framing draws on: hierarchical / empirical Bayes shrinkage & partial pooling (Fay–Herriot; James–Stein);
**CLOSE** (Chen 2022 — precision-dependent EB, the diagnosis that low-strand-precision nodes are the enriched
ones); **EBCF / covariate-powered EB** (Ignatiadis–Wager 2019 — shrink toward a covariate, robustly; here the
"covariate" is the spliced-residual *observation*, not a regression); Gaussian belief-propagation / product-of-
experts in information form (the per-node fusion); Leroux-CAR / graph-MRF (the spatial message smoothness);
panel-of-normals / tangent normalization (the single-strand nodes as the within-sample reference set). The
cross-cutting insight: the per-node solve is already the correct product-of-experts; all the leverage is in
building the *prior inputs* (the mixture model's mode locations + the strand-free per-node observation) correctly.

---

## 9. Validation strategy (non-negotiable)

- Develop on `scripts/debug/toy_prod.py` (production-faithful: full simulator → oracle BAM → real `calibrate()`),
  growing the toys: single-strand multi-exon, microexons (region < fragment → message-only), partial-overlap and
  fully-encompassing opposite-strand AMBIG, multi-gene neighbourhoods, ~100 single-strand background transcripts
  for a realistic training substrate spanning both gDNA modes.
- Metric: **converged post-sweep per-node `|f_g − oracle|`** across the full gDNA:RNA × κ × capture ladder
  (NOT a pooled aggregate; NOT the init prior). Use **partial correlation controlling for true RNA** to avoid
  the synthetic's gDNA/expression co-location trap.
- Then the 16-condition net-flow benchmark (skill `calibration-benchmark`), regen goldens, full `pytest`.
- Real-data validation when the user's new sample lands (the synthetic can't surface prep/GC/mappability/κ
  subtleties; resolve the VCaP κ question before trusting any real density).

---

## 10. Relationship to shipped work

- **Shipped (origin/main 029fb806):** disagreement-aware message precision — the message-precision foundation
  Phases 1/3 stand on.
- **WIP (branch calib-disagreement-precision, e738d5c9):** the spliced→mature MEASUREMENT message (Phase-1
  line-1 work) + the production-faithful toy driver. Follow-ups A/B pending.
- **This doc + NEXT_SESSION_calibration_tasks.md:** the forward plan. Build Phase 1 to pristine (line 1), then
  Phase 2 (this mixture model) + Phase 3 (AMBIG), validating each on toy_prod + the benchmark.
