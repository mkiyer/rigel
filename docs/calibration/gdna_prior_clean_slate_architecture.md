# The gDNA density prior — clean-slate architecture (mixture of background + learned enriched)

**Status:** design / synthesis, 2026-07-05, branch `calib-gdna-accuracy`. The unifying architecture the
nondeterminism investigation converged on. Supersedes the piecemeal picture in
`background_gdna_density_prior_model.md` (M1), `M2_loglinear_sigma2g_design.md` (M2), and
`calibrate_cross_process_nondeterminism.md`. Read those for the evidence; read THIS for the target design.

## 1. What calibration actually needs (the estimand)

Every node has an observed total unspliced density `d = M/E`. Calibration must split it into a gDNA density
`ρ_g` and an RNA density `ρ_r = d − ρ_g`. A node's OWN evidence (the Beta-Binomial strand tilt) constrains
`ρ_g` only when the strand is informative; on unstranded / AMBIG / sparse nodes it does not. There, the split
is decided by a **prior on `ρ_g`** — and the right prior is a **mixture over the population's gDNA-density
modes**: a node's `ρ_g` resolves toward whichever mode it belongs to.

## 2. The conflation we uncovered (the key clue)

The current global Gaussian prior was doing TWO different jobs at once, and that conflation is why no single
deterministic estimator could replace its σ²_g:

1. **Background shrinkage** (what it was designed for): pull a noisy/sparse node's `ρ_g` toward the background
   level, stabilizing zero-inflated counts.
2. **Becoming imprecise under capture** (the accidental job): the fitted σ²_g **blew up to ~31 under hybrid
   capture** and shrank to ~0 without it. That "high variance under capture" is not a background property — it
   is the prior *trying to get out of the way of enriched nodes* (a weak prior can't drag an enriched node
   down). Measured: the spline's 31 is an **artifact** of its monotone-smoothing being pulled up by the
   enriched edges; the true background overdispersion is ~15. Every principled estimator gives ~15 → a *tight*
   prior → it drags enriched nodes → the +12% leak. **The spline was crudely, accidentally re-implementing the
   enriched mode of the KDE inside the background prior's variance.**

So the background prior is confounded: one knob (σ²_g) is forced to encode both "how tight is the background"
AND "is there enrichment to avoid." Those are two different things and must be two different components.

## 3. The clean decomposition — one mixture, two kinds of component

Model the prior on `log ρ_g` as a **mixture density** with two kinds of component:

- **A background mode — ALWAYS present, accurate, tight.**
  Center `ρ_bg` and scale `σ_bg` estimated from **intergenic + intron REGIONS** (the depleted/off-target or
  contamination background). This is the easy, accurate part: `ρ_bg` = the pooled background rate, `σ_bg` = the
  background overdispersion (the true ~15, a small honest scalar) — both deterministic (weighted sums, no fit
  to select). Under capture this is the depleted mode; without capture it is the whole distribution.

- **Enriched mode(s) — LEARNED, may or may not exist.**
  Trained from the node-density distribution itself (the KDE / mixture fit). Present as a second mode ONLY when
  hybrid capture creates an enriched population; absent otherwise. The calibration must LEARN this — it cannot
  be assumed.

The prior is their **sum**. A node's `ρ_g` resolves toward the **nearest mode**.

## 4. "Gravity" = the mixture's shape (this is what the user meant, made precise)

The desired behaviour — *the further a node is from the background, the less the background pulls it* — is not a
special weakening bolted onto the background component; it is the **native shape of a mixture prior**:

- **At a mode** (background or enriched): strong pull → the node's `ρ_g` snaps to that mode (shrinkage /
  stabilization).
- **In the valley between modes**: the log-density is low and flat → weak net pull → the node is drawn to the
  *nearer* mode, not forced.
- **The background's pull on an ENRICHED node is weak** — not because the background component redescends, but
  because the enriched node **belongs to the enriched mode**, which is nearer and shields it. This is precisely
  "gravity toward the nearest mass," and it dissolves the artifact: the background component can be **tight and
  accurate** (σ_bg ≈ 15, the honest value) because sparing enriched nodes is the ENRICHED mode's job, not a job
  faked by widening the background variance.
- **Beyond the outermost mode** (a node denser than any learned mode — e.g. a false-positive high-`ρ_g` node in
  a no-gDNA library where the only mode is background): the mixture's real Gaussian tail is a **spring back** to
  the mode → the false positive is suppressed. (This is the already-verified "real KDE tails vs clamped tails"
  effect: clamped tails → ~585k FP, real tails → ~2k.)

So one mixture gives BOTH behaviours that the conflated prior could not: **spare enriched** (they sit in the
enriched mode's basin) AND **suppress no-gDNA FP** (no enriched mode → the background's spring tail pulls them
back). No capture-detection gate — the mixture is unimodal or bimodal purely from the data.

## 5. Why this resolves everything we hit

- **Determinism:** `ρ_bg`, `σ_bg` are deterministic scalars; the mixture/KDE fit must be deterministic (a
  data-derived bandwidth with no `argmin`/GCV/IRLS — the only remaining requirement). The bistable σ²_g spline
  is deleted outright; there is no density-dependent variance curve to fit at all.
- **The M2 regression is dissolved, not patched:** M2 (tight σ_bg ≈ 15) regressed ONLY because the background
  prior was still solely responsible for sparing enriched nodes. With the enriched MODE present, the tight
  accurate background is correct — enriched nodes are spared by the mode, not by an inflated variance.
- **No conflation:** "how tight is the background" (σ_bg, accurate) and "is there enrichment" (a learned mode)
  are now separate, each estimated by the right data.
- **The historical shift completes:** pre-KDE, pass-1 carried everything (hence the conflated global). The
  train(pass-1)/predict(pass-2) split lets pass-1 learn the modes and pass-2 apply the mixture; the standalone
  background Gaussian becomes just the mixture's always-present background component.

## 6. Implementation path

1. **Background component (deterministic, easy):** `ρ_bg` + `σ_bg` from intergenic+intron regions (this is
   essentially the existing `_floor_estimate` — `ρ_floor` + `s2_floor` — promoted to the background MODE of the
   mixture). Seed it into the mixture as an always-present component so the prior is never mode-less.
2. **Enriched component(s) (learned):** the existing pass-2 KDE, trained on pass-1 node densities — but its fit
   must be made **deterministic** (audit the bandwidth `_weighted_median` and any discrete step; use a
   continuous data-derived bandwidth). The pass-2 KDE dissection (the deferred capon-unstranded peak-smoothing)
   is exactly "does the enriched mode form and shield correctly."
3. **Delete** the conflated global Gaussian's density-dependent σ²_g entirely (the spline / M1 / M2 line of
   attempts) — the mixture subsumes it. Keep only the deterministic background component + the mixture.
4. **Validate:** determinism (ε ≤1×, cross-process bit-identical) + the 16-scenario benchmark (no substantial
   regression; expect the capon-unstranded leak to IMPROVE once the enriched mode shields properly).

## 7. Relationship to the immediate determinism ship

This clean architecture is the correct endpoint but couples the background (M2-tight) with the enriched-mode
shielding (the deferred pass-2 KDE redesign) — they must land together, because a tight background WITHOUT a
working enriched mode drags enriched nodes (the M2 +12%). Two ways to sequence:
- **(A)** Do the full mixture redesign now (background component + deterministic KDE that shields enriched) —
  the principled fix, larger scope.
- **(B)** Ship a **deterministic-λ spline** as an interim (reproduces the current behaviour, zero regression,
  determinism fixed) to unblock the production release, then land the clean architecture (A) as the calibration
  redesign. The interim keeps the module but is purely a determinism patch.

The architecture above is the north star either way.
