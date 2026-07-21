# Enriched-Mode Sensitivity — Hypothesis Generation & Benchmarking Plan

**Goal.** Make the DNA prior **sensitively detect the enriched mode**, measured by the *ultimate* endpoint: the
**projection** of a node's observed DNA density onto the fitted landscape → an anchor `μ*` that a truly-enriched
node lands *at the enriched level* (vs oracle). Landscape shape (EMD) is a means; **the projection is the metric
we optimize**. This round is **sensitivity-first** — enhance the enriched mode even at a specificity cost
(fabrication on zero-DNA); specificity is a later round (§6).

---

## 0. The metric (build this first)

For every node: observed DNA density `d = log10(ĝ/E)` (pass-0), oracle density `t = log10(G/E)`. Project `d`
onto the fitted landscape `P` → `μ*`. **Enriched sensitivity** = for truly-enriched nodes (`t` high), how close
`μ*` lands to `t`. Report per node-type, per scenario, cross-suite (ambig + quick).

**Baseline (the number to beat).** Sampling-likelihood projection `μ* = Σ_j r_j·μ_j`, `r_j ∝ P(μ_j)·N(d;μ_j,h)`,
on the current unified landscape, gdna100/gdna300 capON, truly-enriched nodes:

| scenario | observed `d` | **`μ*` (h=0.15)** | oracle `t` |
|---|---|---|---|
| gdna100 ss0.50 capON | +0.47 | **+0.46** | +0.69 |
| gdna300 ss0.50 capON | +0.80 | **+0.81** | +1.05 |

The projection does **nothing** (μ*≈d) at small kernel, and *pulls down* toward depleted (μ*→+0.28) at large
kernel — because the landscape's enriched mode is diluted to near-zero mass (the reliability weighting silences
the unstranded enriched nodes; diagnosis: they land correctly at ~+0.5 but carry `w≈0.16`). **So μ* cannot rise
to the enriched level until either (a) the landscape carries a real enriched mode, or (b) the projection reads
the mode rather than a depleted-diluted mean, or both.** Every hypothesis below targets one of those.

**Projection-function validation (the projection itself is correct — the landscape is the bottleneck).** Using
the ORACLE landscape (built from true `G`): projecting the *true* density back returns it to ±0.06 (self-
consistency ✓); projecting the *under-called observed* density onto that oracle landscape **pulls it up** toward
truth (gdna100 ss0.50: +0.47 → +0.60 vs target +0.69; stranded ss0.99: +0.65 → +0.68 ✓). Conclusions: (1) the
projection mechanism works and is a trustworthy diagnostic; (2) a real enriched mode **is** what lets it recover
under-called nodes — so §II/§III (build the mode) is the dominant lever; (3) even with a perfect landscape the
recovery is only *partial* (+0.60 vs +0.69) because the mean read-out averages the under-call with the mode — so
§I (read-out / kernel) is worth ~0.09 on its own. Both axes are validated levers.

---

## I. The projection function (the endpoint itself)

**H-I.1 — Mode/peak read-out, not the mean.** The responsibility-weighted *mean* is dragged down by the tall
depleted mode. Read `μ*` as the **posterior mode** (argmax of `r_j`) or the *nearest local maximum* of
`P(μ)·N(d;μ,h)`. A node observed at +0.47 is *closer* to a +0.7 enriched hill than to the −1.8 depleted mode, so
a mode read-out can snap it enriched even when the enriched mass is small. *Test:* argmax vs mean vs median μ*.

**H-I.2 — Projection kernel width `h_proj` (adaptive).** Small `h` → μ*≈d (no correction); large `h` → depleted
dominates. There may be an intermediate `h`, or a **per-node `h`** set by the node's own belief precision (a
confident node projects narrowly, an uncertain one is allowed to reach a mode). *Test:* sweep `h_proj`; derive
per-node `h` from `var`.

**H-I.3 — Posterior/MAP projection.** Treat it as a proper Bayesian update: node likelihood `N(d; μ, σ_obs)`
(σ_obs from the node's count/belief) × landscape prior `P(μ)` → posterior; read its mode. Principled version of
H-I.1/2 where σ_obs is derived, not tuned.

**H-I.4 — Local (proximity-gated) projection.** Restrict the responsibility to a *window* around `d` (or a
distance-decaying weight steeper than Gaussian), so a far depleted mode cannot vote at all — the node is anchored
only by nearby landscape structure. *Test:* truncated / heavy-tailed proximity kernels.

**H-I.5 — Two-mode assignment.** Fit/assume a 2-component landscape (depleted + enriched) and project = soft-assign
`d` to the nearer component, returning that component's location. Sidesteps the dilution entirely.

---

## II. Enhancing the enriched mode in the landscape

**H-II.1 — Boundaries build the enriched mode (the owner's lever).** Boundaries self-solve (spliced RNA + intron
DNA), so a real-gDNA boundary is *confident* while a `none_*` boundary confirms RNA. Diagnosis: boundaries already
carry more enriched weight than exons (w 0.20 vs 0.11). *Tests:* up-weight boundaries in the enriched region;
build the enriched mode **primarily/only** from confident boundaries; a boundary-only enriched sub-landscape.

**H-II.2 — Two-population landscape.** The depleted mode is a *free win* (confident, easy). Fit it from confident
nodes, and fit the **enriched** mode *separately* from the enriched-carriers (boundaries / high-density / spliced-
resolved), then **stack without letting the depleted mass dominate the enriched normalization**. *Test:* per-region
or per-component normalization so the enriched mode isn't drowned.

**H-II.3 — Gentler / floored reliability weighting.** The weight `w=S0/(v·g/(g+1)+S0)` is too aggressive on
enriched nodes. Try `w^α` (α<1), a **weight floor** `w←max(w,w0)`, or make the temper `g/(g+1)` weaker so
high-count enriched nodes keep mass. *Test:* α and floor sweeps; measure the enriched μ* recovery vs the zero-DNA
fabrication cost.

**H-II.4 — Density up-weighting.** Counter the dilution directly: weight a node's landscape mass by its density
(or by `enrichment`), so enriched nodes contribute proportionally more to the enriched region. Principled version:
weight by the node's *information about the enriched region* rather than uniform mass.

**H-II.5 — Adaptive (balloon / k-NN) bandwidth.** Enrichment is a **spectrum** (nodes at +0.2…+0.8), so a *larger*
kernel in the *sparse enriched region* aggregates them into one detectable mode, while the dense depleted region
keeps a small kernel. *Test:* k-NN or density-adaptive bandwidth (wide where nodes are sparse). This is the
principled form of "smoother helps the enriched mode."

**H-II.6 — Spliced-resolved gDNA at exons.** Use the spliced fragments to peel mature RNA; the residual unspliced
mass at an exon is then a *confident* gDNA measurement at enriched density — a low-variance enriched carrier that
survives the weighting. *Test:* add spliced-derived gDNA nodes to the enriched substrate.

**H-II.7 — Peak-preserving fit.** Whatever smoothing we use, protect local maxima so a real (small) enriched peak
is never smoothed into the depleted shoulder (a shape-aware / edge-preserving smoother).

---

## III. External signals to break the identifiability

Real-enriched and over-called nodes are indistinguishable *by variance*; distinguish them *by evidence*.

**H-III.1 — Enrichment covariate (Role A total density).** A node's total density (`μ_proj`, already computed) is
observed and not under-called. Use it to *place / boost* the enriched mode: a high-total-density node is enriched,
so up-weight its enriched contribution or condition the landscape on enrichment level. (This is the
enrichment-conditioning idea, redeployed for *sensitivity* not composition.)

**H-III.2 — Spliced fraction.** A node with high unspliced density but few spliced reads is more likely gDNA
(enriched); one with abundant spliced reads is RNA. Use the spliced fraction to gate/weight enriched-mode
membership — the direct signal boundaries already exploit.

**H-III.3 — Intron-neighbor DNA.** Boundaries infer DNA from adjacent introns; propagate that confident intron-DNA
to place the boundary's enriched contribution.

---

## IV. Reframing the object (bigger swings)

**H-IV.1 — Explicit mixture model.** Fit `P = π_d·Depleted + π_e·Enriched` (each a parametric or KDE component);
the projection assigns to a component. Makes the enriched mode a first-class object with its own mass `π_e`,
immune to depleted dilution.

**H-IV.2 — Rank / top-density anchor.** The enriched mode's *location* is well-defined by the densest confident
nodes (order statistics); estimate only its presence/mass. Robust to the bulk depleted mass.

**H-IV.3 — Iterative enrichment.** Project enriched nodes up → re-fit the landscape with anchored densities → the
enriched mode strengthens → repeat. Self-reinforcing; needs a guard against runaway fabrication (couple to the
boundary/spliced evidence of §III).

---

## V. Benchmarking plan

1. **Build the projection eval** into `gdna_explore_lib` (µ* per node from a landscape + projection fn) and an
   **enriched-sensitivity metric**: mean/among truly-enriched nodes, `|μ* − t|` and the fraction landing in the
   enriched region; report per node-type × scenario × suite. This is the single number every hypothesis moves.
2. **Individual ablation** — each hypothesis as a one-line variant against the baseline (µ* +0.46 → +0.69 target),
   cross-suite (ambig + quick). Rank by enriched-μ* recovery.
3. **Combinations** — the promising axes are orthogonal: {mode read-out (I.1) or adaptive `h_proj` (I.2)} ×
   {boundary-built / two-population enriched mode (II.1/2) or gentler weighting (II.3) or adaptive bandwidth (II.5)}
   × {enrichment/spliced signal (III)}. Test the Cartesian of the winners.
4. **Cross-suite + spectrum** — always both suites; break out by DNA level and capture, since the enriched mode
   only exists capture-ON and grows with DNA. Keep an eye on capOFF/zero-DNA (must not fabricate a mode there even
   in the sensitivity round — that's the specificity canary).

### Workflow families to launch (creative, parallel)
- **W-A "projection read-out"** — I.1–I.5: mode vs mean vs MAP, `h_proj`, proximity gating. (Can win *without*
  changing the landscape.)
- **W-B "enriched substrate"** — II.1/2/6 + III: boundary-built, two-population, spliced-resolved, enrichment
  covariate — which node set + signal builds the strongest enriched mode.
- **W-C "weighting & bandwidth"** — II.3/4/5/7: gentler/floored/density weighting, adaptive bandwidth, peak-
  preserving — how to stop diluting/smoothing the enriched mode away.
- **W-D "reframe"** — IV: explicit mixture, rank-anchor, iterative — bigger-swing objects, higher risk/reward.
- **Synthesis** — combine the winning read-out × substrate × weighting; report enriched-μ* recovery + the
  fabrication canary, cross-suite.

---

## VI. The tradeoff (explicit, deferred)

Every sensitivity lever risks **fabricating** an enriched mode on zero-DNA (specificity ↓). We accept that this
round and *measure* it (the capOFF/none canary), but the durable answer to both is **§III — evidence, not
variance**: boundaries and spliced fractions distinguish real enrichment from over-call, so the same signals that
build the enriched mode (sensitivity) also gate out fabrication (specificity). Sensitivity now; wire the evidence
gate in the specificity round.

---

## Appendix — where things stand
Landscape recipe (`scripts/debug/gdna_landscape_recipe.py`, unified: one reliability-mass rule, zero-native
Poisson, 1 constant) — robust cross-suite (mean EMD 0.267) but **enriched-blind** (µ* +0.46 vs +0.69). Shared
substrate cache + eval lib: `gdna_cache_build.py` / `gdna_explore_lib.py` (ambig 32 + quick 16 scenarios, region +
boundary nodes). Diagnosis figures: `docs/figures/landscape_synth.png`, `bandwidth_enriched.png`.
