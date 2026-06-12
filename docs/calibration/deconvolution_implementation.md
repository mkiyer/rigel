# Converting the generative model to an implementation — with hybrid capture

**Status:** implementation blueprint, 2026-06-11. Bridges `deconvolution_generative_model.md` (the complete
three-species model) to code, with **hybrid capture** as the organizing constraint. Capture is what makes
this *local*: there is no global gDNA density to fit, so each region's gDNA is recovered from its own
strand and its **bounding boundaries**, and imputed from the nearest informative boundary where local
strand fails. Boundaries are not a detail — they are the atoms.

## 1. What capture does, and why it forbids a global fit

Probes target exons. A captured fragment overlaps a probe and extends a fragment-length (`FL`) beyond it.
So the gDNA density `ρ(x)` is a **local profile**: high over exons and ~`FL` into the flanking introns
(captured spill), dropping to a low off-target floor in deep introns / intergenic.

Consequences that dictate the architecture:

- **No global `ρ₀`.** The exon density and the deep-intron density differ by the enrichment factor (often
   orders of magnitude). A single fitted density — or imputing an exon from a far region — is wrong. *Fit
  locally.*
- **The contained density of an exon is unusable directly** — it is gDNA + nascent + mature-contained
  (the three-species mix). Strand cleans it, but only where strand is informative.
- **Boundaries are the clean local probes.** The crossing-unspliced flux at a boundary is **mature-free**
  (mature crossing a junction is *spliced*, not unspliced) — so it is gDNA + nascent only, no mature
  subtraction needed. And the captured gDNA spanning an exon edge is at the **on-target** density. So a
  boundary's strand-deconvolved crossing flux is a clean, local, capture-correct gDNA-density probe at
  that point. **This is why boundaries are essential for imputation** — they are the only place we read
  the local gDNA density without the mature contaminant and at the right enrichment level.

## 2. The pipeline (model → stages → modules)

Library-global parameters (legitimately global — fit from all reads) stay global; **only the density is
local.**

| stage | what it computes | from the model | module (re-based) |
|---|---|---|---|
| **S0 Library fits** | `κ`, strand overdispersions, gDNA/RNA FL pmfs | the rate `κ`, the noise | `strand_balance`, `gdna_strand`, `fl` — *intact* |
| **S1 Mature subtraction** | per node, predict mature-contained `= S·(E_mu/E_ms)` and subtract from the unspliced | §3 step 2 (spliced ⇒ mature-contained) | **new**; pure geometry `E_mu/E_ms` |
| **S2 Local strand deconvolution** | per region *and* per boundary side: split the (mature-subtracted) unspliced into gDNA(½) vs nascent(κ); emit the **likelihood** for the local gDNA density | §3 steps 3–4 | `strand_deconvolve` (Phase 1), + the S1 subtraction for regions |
| **S3 Boundary probes** | each boundary's crossing flux → local gDNA-density likelihood (mature-free, capture-correct) | the boundary is gDNA+nascent only | `strand_deconvolve` on the side views |
| **S4 Per-region local density** | combine a region's *own* strand likelihood with its two **bounding-boundary** likelihoods (multiply — same local `ρ`); impute regions whose local strand fails from the nearest informative boundary; carry per-strand nascent into AMBIG | the joint likelihood for `ρ_region`; the carry-over supplies the missing AMBIG equation | `node_gdna_density` re-based; `splice_junction` **subsumed** (the boundary probe *is* the splice fraction) |
| **S5 Per-locus prior** | aggregate gDNA vs RNA over the locus → the two Dirichlet scalars | gDNA = `ρ·E_g`; RNA = the rest + spliced | `derive`, `priors` — ~intact |

**S2/S3 are one operation** (the per-node strand likelihood), applied to the contained view (regions,
after S1 mature subtraction) and to the side views (boundaries, mature-free). **S4 is the only genuinely
new spatial logic**, and it is local by construction.

## 3. S4 in detail — local, principled, capture-aware

A region's gDNA density `ρ_r` has up to three independent local measurements, each a likelihood from S2/S3:
its **contained** strand (mature-subtracted) and its **left/right boundary** crossings. Under capture the
region is signature-homogeneous (all-exon or all-intron) and small (~`FL`), so these measure the *same*
local `ρ_r` — **combine by multiplying the likelihoods** (the ML local density; precision = the summed
curvature, derived, not a weight we invent). Then:

- **Strand-informative region** (single-strand, `κ` off ½, enough counts): its own + boundary likelihoods
  pin `ρ_r` locally. Done — no imputation, no global density.
- **AMBIG region** (overlapping `+`/`−`): contained strand is flat (two free nascent rates); the *boundary*
  crossings are likewise AMBIG. So `ρ_r` is under-determined locally → **carry the per-strand nascent
  from the nearest single-strand neighbour** (the §6 carry-over = supplying the missing equation), which
  is itself a local boundary probe, then solve. The neighbour must be **capture-compatible** (on-target ↔
  on-target): an AMBIG *exon* inherits from a nearby single-strand *exon* edge (on-target), never a deep
  intron.
- **κ ≈ ½ (unstranded)**: no strand anywhere → the only gDNA signal is the *magnitude* where the signature
  forbids RNA (intergenic / off-target, where `U ~ Poisson(ρ·E_g)`). On-target exon gDNA is then
  under-determined — the documented unstranded+capture floor. The implementation should *say so* (low
  confidence) rather than manufacture a value.
- **Imputation rule (capture-aware):** a region with no usable local likelihood inherits the nearest
  informative boundary probe **within its enrichment class** — carried along the genome by the existing
  run-fill, but gated so on-target and off-target do not cross-contaminate. The signature (exon vs
  intron/intergenic) is the enrichment label.

## 4. What is principled vs. what is still open

**Principled (derived from the model, no invented weights):**
- Each measurement is a likelihood; combining is multiplication; the "precision" is the curvature `−ℓ″`.
- The strand isolates gDNA from *all* RNA (nascent+mature), so the gDNA prior is well-posed (rev-2 §4).
- The carry-over is "supply the missing equation," not a heuristic transfer.

**Open — to settle before/while coding:**
1. **The `E_mu/E_ms` geometry (S1).** Mature-contained vs junction-spanning effective lengths as a
   function of exon length and RNA FL — pure geometry, reuses `region_eff_length`/`boundary_eff_length`,
   but must be derived carefully (it is the short-exon FL-consistency axis, now correctly per-species).
2. **The local-likelihood combination (S4).** Exact (profile the free nascent/mature and ML the local
   `ρ_r` over the region+its boundaries) vs a Gaussian/Laplace approximation. Exact is cleaner and avoids
   the "made-up weight" objection; the Laplace form *is* the rejected precision, so derive it as the
   curvature if approximating.
2.5 **Capture-class gating of imputation (S4).** How strictly to forbid on/off-target mixing — by hard
   signature class, or a soft enrichment estimate. Hard class is the simple first cut.
3. **The carry-over geometry for AMBIG** — single-strand neighbour selection along the transcript extent
   vs the run; capture-compatibility. (The toy AMBIG scenario is the test bed.)

## 5. Relationship to the existing code (this is a re-base, not a rewrite)

The current modules already have the *shape* — local density (`node_gdna_density`), boundary crossings
(the side views), the splice fraction (`splice_junction`), per-region capture-contracted lengths (the
IPR). The re-base: (a) add S1 mature subtraction; (b) make S2/S3 the generative strand likelihood (Phase
1's `strand_deconvolve`, already built); (c) replace `node_gdna_density`'s heuristic raw-count
imputation + `splice_junction`'s 2/3-term with the **S4 local-likelihood combination** (boundaries as
probes, capture-gated imputation, carry-over). `derive`/`priors`/EM are ~intact.

So the phase plan reshapes to: **Phase A** S1 mature-subtraction geometry (`E_mu/E_ms`) + validate on a
short-exon scenario; **Phase B** S2/S3 local likelihoods (extend the built `strand_deconvolve` with S1);
**Phase C** S4 local-density combination + capture-gated imputation + carry-over (the heart), validated on
the AMBIG + capture toy scenarios; **Phase D** wire into `calibrate`, retire the blend, golden, suites.
Each is plan-then-implement as before.

## 6. The one-paragraph summary

Capture removes the global density, so the implementation is local: the **boundary crossing flux**,
strand-deconvolved, is the clean (mature-free) capture-correct local gDNA-density probe; each region
combines its own strand likelihood with its bounding boundaries' to pin its local `ρ` by maximum
likelihood; AMBIG and unstranded regions inherit from the nearest *capture-compatible* informative
boundary (the carry-over supplying the missing equation), and where even that fails we report the
unstranded+capture floor rather than invent a number. No global fit, derived weights, boundaries as the
atoms.
