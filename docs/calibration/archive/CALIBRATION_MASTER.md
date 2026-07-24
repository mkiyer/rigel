# Calibration — Master Reference & Roadmap

**Status:** the single source of truth for the calibration stage. Written 2026-07-18 to consolidate a scattered
set of roadmaps, derivations, and diagnosis notes into one big-picture reference. Where any other calibration
note disagrees with this document *about the design*, prefer this document (and fix the other). The **code is the
source of truth about what is currently implemented**; §8 tracks the gap between the design here and the code.

**Audience:** the calibration developer (phase-by-phase work) and external reviewers (the derivations in §5–§7
are reviewer-facing; the detailed proofs live in the linked companion docs).

**Companion anchors** (kept; everything else is archived — see §10):
- [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) — **the one principle** (count-zero-information)
  and the three information sources. Authoritative on the *invariant*; this document is authoritative on the
  *design and plan*.
- [`reference_prior_derivation.md`](reference_prior_derivation.md) §10 — the resolved Beta(½,½) reference and
  the 1-D/2-D solver reconciliation (reviewer-signed-off).
- [`background_reference_derivation.md`](background_reference_derivation.md) — the aggregate background DNA
  level ρ_bg (the data-gated global hyperprior).
- [`transfer_variance_formal_derivation.md`](transfer_variance_formal_derivation.md) — the σ²_transfer message
  precision via projection onto the enrichment NPMLE (reviewer-signed-off).
- [`message_absorption_fix.md`](message_absorption_fix.md) — a **standing correctness item** (the diagnosed
  splice-junction RNA-imputation false-positive) to land alongside the Role-B representation restore (§9 Phase A).
- [`kde_vs_npmle_enriched_mode` (memory)] — why the shipped v0.7.1 KDE keeps the enriched mode and the current
  NPMLE silences it (§5 Mistake 2 / §9 Phase A) — the current crux.
- [`oracle_and_benchmarking.md`](oracle_and_benchmarking.md) — the trusted evaluation path.

---

## 1. What calibration does (the goal)

Calibration deconvolves each genomic **node**'s observed **unspliced** fragment mass into a composition:

```
    single-strand node:  (f_rna, f_gdna)            f_rna + f_gdna = 1
    AMBIG node:           (f_pos, f_neg, f_gdna)     f_pos + f_neg + f_gdna = 1
```

A node is a **region** (a segment of the genome between transcript boundaries) or a **boundary** (an
exon/intron junction). Calibration models **only RNA-vs-gDNA** — it does **not** split nascent from mature RNA;
that is the per-locus EM's job downstream. (Nascent-vs-mature is only identifiable via intronic RNA fragments: a
single-exon transcript is simultaneously mature and nascent, so the split is a per-locus question, not a
calibration one.)

**Outputs consumed downstream:**
- the per-region and per-boundary-side **gDNA / RNA mass split** → two per-locus Dirichlet scalars
  (`alpha_rna_add`, `alpha_gdna_add`) that set the gDNA-vs-RNA split the EM starts from;
- library-level scalars: `rna_sense_frac` (κ), the gDNA / RNA strand overdispersions, `gdna_density_global`.

**Why it matters:** an accurate gDNA fraction lets the per-locus EM attribute background-DNA mass correctly
instead of inflating transcript abundances.

---

## 2. The one principle (count-zero-information)

> **A fragment count — or equivalently a fragment density — carries ZERO intrinsic information about a node's
> gDNA/RNA composition.** The count may enter ONLY as *precision* (statistical power), never as a composition
> vote. ([`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) §0.)

100 unspliced fragments split 50/50 across the two genome strands is *equally consistent* with pure gDNA
(genomic strand is symmetric) and pure RNA in an unstranded library (κ=½). The count cannot say which. This is
the invariant every past calibration bug violated, and the one this whole roadmap is organized around. The
legitimate use of the count is as **statistical power**: how *sharply* a strand tilt is measured, or how noisy a
neighbour's density *estimate* is.

**Corollary that governs the prior (this is the crux we lost and recovered):** a prior fit on the distribution
of **total density** is a model of *density*, not of *composition*. It cannot be used to vote a node's gDNA
fraction. It can only characterize **enrichment vs depletion** and set message precision (§5).

---

## 3. The information sources for the gDNA-vs-RNA split

A node's composition is set by exactly these, in priority order ([`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) §1):

1. **Structure** (the annotation, via the region **signature**). Intergenic ⇒ no RNA possible ⇒ f_gdna = 1
   (a structural certainty, locked outside the solver). TSS/TES flanks, single-strand vs both-strand
   continuity — all structural. This is the initialization.
2. **Strand** — the only *intrinsic* per-node signal. RNA tilts the genome-strand rate toward κ; gDNA stays at
   ½. Its Fisher information for f_gdna is `∝ n(½−κ)²` — **exactly zero at κ=½** (unstranded), for any count.
   A stranded node with an informative strand is *solved by the strand alone*.
3. **Imputation (messages)** — a higher-precision neighbour propagates its solution across the node chain,
   transcript-structure-aware (gDNA flows genomically; per-strand RNA flows only where that strand is
   continuous). The message *content* is a density; its *precision* is the belief-free transfer variance (§5).
   **The message-content rule (owner-confirmed):** a **region** reports only the densities it *has*
   `(ρ_RNA+, ρ_RNA−, ρ_gDNA)` — it does no splicing subtraction, because it cannot know what splices next and
   we can **only measure** nascent, never infer it from abundance. A **boundary** owns the splicing information:
   it partitions its own directly-measured spliced (mature that left) vs unspliced (nascent + gDNA) on its own
   effective length, and the subtraction lives on the boundary's **outgoing** side (the nascent continues, the
   spliced does not). The current code violates this — it subtracts a *destination boundary's* spliced from a
   *neighbour's* RNA imputation, which goes negative and is laundered into a confident false gDNA call
   ([`message_absorption_fix.md`](message_absorption_fix.md) — the immediate next fix).
4. **The gDNA hyperprior** — the expected DNA level. Needed to break the AMBIG two-root ambiguity (§4). A
   population baseline, at population precision — never overriding a node's own strand.

The count/density is **not** a source. It sets precision only.

---

## 4. The AMBIG two-root problem — why the gDNA hyperprior exists

The gDNA hyperprior is *not* a general-purpose "make the answer better" knob. It exists for one specific,
fundamental reason: **at a node where transcripts on opposite strands overlap gDNA, the composition has two
distinct valid roots, and no per-node signal distinguishes them.**

The clearest case: an unspliced density of 10 with a **50/50 strand balance** is equally explained by

- **(a)** equal-density RNA on the + and − strands (two overlapping antisense transcripts), *or*
- **(b)** 100% gDNA with no transcription at all.

Across the whole range of strand balance, two roots persist. The strand likelihood *cannot* choose (it is
consistent with both); the messages carry density, not the root selection. **Only an external statement of "how
much DNA do we expect at this position" grounds the solve** and tells us how to apportion the unspliced mass
among RNA₊, RNA₋, and gDNA. That is the gDNA hyperprior.

**Consequences that shape the design:**
- **Single-strand nodes do not need the gDNA hyperprior** — there is no opposite-strand conflict, so the two
  roots collapse to one. They are solvable from strand + structure + messages alone.
- The hyperprior is nonetheless a **useful crutch** elsewhere: hybrid-capture enrichment, low counts, and
  sampling stochasticity make a grounding prior helpful even where it is not strictly required. But that
  usefulness is secondary; the *reason* it exists is the AMBIG two-root ambiguity.
- Because it must state "the expected DNA level," it **must be fit on DNA** — i.e. on the *deconvolved* gDNA,
  not on total density (§5). This is the whole point.

---

## 5. The TWO roles of the NPMLE — the critical distinction

We fit a Fixed-Kernel Poisson-lognormal Mixture NPMLE (`DensityNPMLE`) over per-node densities. It serves
**two entirely different purposes**, and conflating them is the mistake that broke calibration. They are
different *fits* on different *data* for different *jobs*.

### Role A — the enrichment profile → transfer variance (message precision)

- **Fit on:** **total** density at all nodes, **before** any solve (belief-free).
- **What it is:** a model of the library's hybrid-capture **enrichment/depletion landscape** — where the density
  modes are (a depleted off-target mode, one or more enriched on-target modes). It is a model of *density*, and
  says nothing about composition (§2).
- **What it is for:** message precision. Belief propagation must assign a variance to each source→destination
  density message. The solution ([`transfer_variance_formal_derivation.md`](transfer_variance_formal_derivation.md),
  reviewer-signed-off): **project** the observed source and destination densities onto the NPMLE
  (`DensityNPMLE.project` → `mu_proj`, `var_proj`) — i.e. ask "is this enriched or depleted?" — and compute
  `σ²_transfer = var_proj[dst] + (mu_proj[dst] − mu_proj[src])²`. This makes message precision
  **belief-free** (independent of the density itself) and **robust**: it compares the source and destination
  through their *population-level* enrichment assignment rather than differencing two noisy raw densities that
  might agree or disagree by chance.
- **Why it is fit before solving:** we need transfer variance to run the first-pass solver at all, and the
  enrichment landscape is knowable from total density without any composition solve. This is a *good* and
  necessary use — do not remove it.

### Role B — the gDNA hyperprior (the composition prior)

- **Fit on:** the **DECONVOLVED gDNA density** `ρ_g = f_gdna · M / E`, **after** an initial solve, with **AMBIG
  nodes excluded** (they are the target, so including them would be circular) **and intergenic nodes excluded**
  (they are structurally gDNA-locked and would flood the depleted mode — they are represented by the weak
  background floor instead; see `gdna_kde_restore_plan.md`).
- **What it is:** a model of the *DNA* density distribution — the background mode *and* the enriched-DNA modes,
  learned only from mass that was already attributed to DNA.
- **What it is for:** the §4 AMBIG root selection — the expected DNA level per position.
- **Representation (§5 Mistake 2 / `gdna_kde_restore_plan.md`):** an **additive, occupancy-weighted, fixed-bandwidth
  KDE** on the deconvolved gDNA + a **weak (1-pseudo-obs) floor** at `ρ_floor` — so the enriched mode is
  guaranteed and the background cannot dominate. **Stratification by signature class is REJECTED** (owner
  directive 2026-07-19) — not pursued.

### Two mistakes: one undone, one LIVE

**Mistake 1 — total-density fit (UNDONE).** An earlier state of this branch fit the NPMLE on **total density**
(Role A's data) and wired it into ψ as the composition prior (Role B's job) — a total-density prior votes
composition = the §2 violation. **This is fixed:** `_fit_gdna_hyperprior` now fits on the **deconvolved** gDNA
`ρ_g = f_g·M/E` after a prior-free pass-0, in a 2-pass structure (`calibrate.py`; Phase 2 landed 2026-07-19).

**Mistake 2 — the NPMLE *representation* silences the enriched mode (LIVE, the current crux).** Even fit
correctly on deconvolved gDNA, the current `DensityNPMLE` **loses the capture-enriched mode** under hybrid
capture — a tall "depleted" background mode drowns the short "enriched" modes, and partially-enriched nodes snap
to f_gdna≈0. The 2026-07-19 investigation (`kde_vs_npmle_enriched_mode`) found this is **structural to how the
NPMLE represents and combines the density**, NOT the training set (the shipped KDE also trains on
intergenic+intron+exon) and NOT primarily the bandwidth:

- **Competitive EM mixture vs additive kernels.** The NPMLE fits component weights by **EM over a shared,
  normalized mixture** — so the numerous *depleted* population **wins the weight and starves the enriched
  components** ("population concentration defeats the honest width"). The shipped v0.7.1 **KDE** instead places
  **one unit-weight kernel per node**, so every enriched node **guarantees its own bump** — nothing can compete
  it away.
- **Heavy in-mixture background vs a separate capped expert.** The NPMLE injects the background as **one
  aggregate cell inside the mixture with EM weight = `n_regions`** (the whole intergenic population) at `eff=ΣE`
  (a sharp low mode) — mechanically lowering the enriched mode's *relative* height in the normalized `logP`.
  The KDE keeps the background as a **separate M2 Gaussian CAPPED at one pseudo-observation**, combined by
  **product-of-experts** (absolute height cancels across the composition grid), with floor nodes **excluded**
  from the KDE entirely. On real data (100k+ intergenic regions) the NPMLE's `n_regions`-weighted cell is
  guaranteed to dominate; the KDE's cap makes it structurally immune.
- **Precision weighting.** The KDE uses **unit weight** deliberately (precision correlates with density, so
  precision-weighting biases the distribution's *shape*); the NPMLE re-introduced belief-width weighting, which
  only *broadens* likelihoods and does not rescue the starved enriched mode.

This is why the released v0.7.1 KDE **retained** the enriched mode and the redesign lost it — a **restore of the
KDE's structural choices**, not a new invention (§9 is re-scoped around this).

---

## 6. The calibration pipeline (the correct order)

```
 1. INITIALIZE          signature-binary: intergenic / TSS sinks → f_gdna=1 (locked);
                        single-strand → strand solve seed; AMBIG → uncertain.

 2. LIBRARY SCALARS     κ = rna_sense_frac; gDNA & RNA strand overdispersions;
                        ρ_bg = the aggregate intergenic background DNA level (data-gated;
                        ρ_bg=0 for a DNA-free library ⇒ no gDNA sway ⇒ no false positive).

 3. ENRICHMENT NPMLE    fit DensityNPMLE on TOTAL density (Role A). Use it ONLY for
    (Role A)            σ²_transfer / message precision — NOT as a solver composition input.

 4. INITIAL SOLVE       ψ = strand + reference(Beta½,½) + ρ_bg hyperprior + messages.
    (prior-free of      NO deconvolved-gDNA hyperprior yet. Single-strand nodes solve here;
     the DNA prior)     AMBIG nodes remain at the honest reference (they await step 6).

 5. FIT gDNA HYPERPRIOR fit a SEPARATE NPMLE on the DECONVOLVED gDNA density ρ_g=f_g·M/E from
    (Role B)           step 4, AMBIG excluded. (Representation must keep the enriched mode —
                        §5 Mistake 2 / §9.)

 6. REFIT SOLVE         re-solve with the gDNA hyperprior added → grounds the AMBIG two-root
                        selection. Iterate 5↔6 only if it demonstrably helps (guard against
                        the self-training drift — see the refit caution below).

 7. PROJECT             per-region / per-boundary-side gDNA-RNA mass → the per-locus prior;
                        gdna_density_global (a QC scalar).
```

**The refit caution (5↔6 iteration):** re-fitting the DNA prior on the solve's own output is self-training and
can drift toward an over-smoothed fixed point. Keep any iteration cold-restarted and gated on a truth-free
disagreement instrument (the solved between-node density disagreement should approach the *population* level,
not zero). See the archived refit study for the evidence; do not turn iteration up without measuring.

---

## 7. The per-node solver (the two modes, reconciled)

The per-node objective ([`reference_prior_derivation.md`](reference_prior_derivation.md) §10, resolved &
reviewer-signed-off):

```
    ψ = strand  +  gDNA arm  +  RNA arm  +  messages
```

- Each **arm** is that component's **fitted logP if we have one, else the Jeffreys reference** `+½·log f_c`
  (`_JEFFREYS_REF = 0.5`). Writing nothing for a component is **not** neutral — the uniform-λ grid supplies a
  Haldane (improper, vertex-amplifying) reference by default. The Beta(½,½) reference is the *derived* answer
  (Poisson-rate Jeffreys = composition Jeffreys, two routes), and it is a *structural / declared* prior — §2's
  count-zero-info invariant governs the **likelihood**, not the prior, so the reference is legitimate
  (`reference_prior_derivation.md` §10.6).
- **The gDNA arm takes the deconvolved-gDNA hyperprior (Role B) once fitted, else the reference.** It must
  **never** take the total-density NPMLE (Role A) — that is §5's mistake. **No `+½·log f_g` on top of a fitted
  gDNA arm** (that double-counts).
- **The RNA arm is the `+½·log(1−f_g)` Beta(½,½) reference — permanently.** Fitting a learned `logP_r` was
  tested and **REJECTED** (full-solve A/B, 2026-07-19: it does not recover under message anchoring and
  manufactures false positives — `gdna_projection_snap_fix`). The RNA-side constraint is carried by this
  reference plus the messages; there is no fitted RNA arm on the roadmap.

**The two solver modes are the same model** (`reference_prior_derivation.md` §10.3–§10.7):

| | 1-D single-strand (`_solve_nodes_logodds`) | 2-D AMBIG (`_solve_ambig_logodds`) |
|---|---|---|
| latent | λ = logit(f_g) | λ and the tilt θ = arcsin(τ) |
| reference | ½log f_g + ½log(1−f_g) | **identical on f_g** — the tilt term −½log(1−τ²) is cancelled EXACTLY by the θ=arcsin(τ) Jacobian |
| class gap | — | **0.0, verified** (marginal consistency by construction) |
| grid | **256** λ-points (`sweep_n_grid_single_strand`) | **60** λ × 60 θ (`sweep_n_grid`) |
| precision | float64 | float32 cube (memory: the AMBIG cube is O(m·K²)) |

So the two implementations do **not** disagree on the math. They differ only **numerically**, for performance:
the 256-vs-60 grid resolution (a single-strand node is solved 4× finer than an AMBIG node on the same axis) and
the `(1−τ²)^{−½}` θ-endpoint singularity (currently a crude half-grid-step clip; Gauss–Jacobi quadrature is the
principled fix and could let the AMBIG grid be coarser). These are real engineering debt (§9 Phase 4), **not**
the source of the calibration error. The known-red `test_gdna_sweep_zero_gdna_pin_and_monotone` (0.5638) is the
*correct* behaviour of a proper reference on an unidentifiable node (`reference_prior_derivation.md` §10.7) —
re-derive its bound when initialization lands, or retire it.

---

## 8. Implemented vs. needs-restructuring (the current gap)

| piece | design (this doc) | code today (branch `calib-ambig-init-wip`) | verdict |
|---|---|---|---|
| Initialization | signature-binary | implemented (`node_geometry.init_beliefs`) | ✅ |
| Library scalars (κ, overdispersions) | §6.2 | implemented | ✅ |
| Background ρ_bg | §6.2 data-gated hyperprior | implemented (`background_reference.py`) | ✅ |
| The Beta(½,½) reference | §7 | implemented (`_JEFFREYS_REF`, both arms, both solver modes) | ✅ |
| Messages / relay | §3.3 | implemented (`bp_solver._scan`, the coherent (λ,θ) relay) | ✅ structure |
| **RNA message content (region reports / boundary subtracts)** | §3.3 | ❌ **BUG**: dst-boundary spliced subtracts from a neighbour's RNA imputation → negative → laundered to confident false gDNA | ⛔ **fix NEXT** |
| **Honest clamp (cancelled imputation ⇒ no message)** | §3.3, §5B | ❌ `max(rho, 1/erd)` manufactures a confident "no RNA" | ⛔ **fix NEXT** |
| σ²_transfer (Role A) | §5A | implemented (`DensityNPMLE.project`, F1) | ✅ |
| Two-mode solver | §7 | implemented; reconciled model, split grids | ✅ model / ⚠ grid debt |
| **Enrichment NPMLE used ONLY for Role A** | §5A | separated: `enrichment_prior` (σ²_transfer) ⟂ `gdna_prior` (arm) | ✅ **Phase 1a landed** |
| The ρ_bg clamp/floor | §5B (rejected) | **removed** — a clamp is a cliff; background enters smoothly via the pinned component | ✅ **Phase 1a landed** |
| Initial solve = inert reference alone | §6 step 4 | `gdna_prior=None` + strand + messages + σ²_transfer | ✅ **Phase 1a landed** |
| **gDNA hyperprior (Role B), deconvolved 2-pass** | §5B, §6.5 | ✅ fit on deconvolved ρ_g in a refit (`_fit_gdna_hyperprior`, Phase 2 landed) | ✅ substrate / ⛔ **representation** |
| **gDNA-hyperprior representation (keep the enriched mode)** | §5 Mistake 2 | ⛔ competitive EM mixture + heavy `n_regions` in-mixture background **silences** the enriched mode | ⛔ **the LIVE crux — §9** |
| `logP_r` (fitted RNA hyperprior) | — | **REJECTED** — full-solve A/B failed (2026-07-19); RNA arm = reference permanently | ✖ **removed** |
| Refit order (prior-free → fit-on-deconvolved → refit) | §6 | initial prior-free solve landed; refit is Phase 2 | ◐ **half** |

**Phase 1a A/B (2026-07-18, production-faithful):** flagship gdna300-unstr-capON **0.592→0.266** (crush mass
35%→0.1%); gdna100-unstr **0.377→0.244**; stranded-capON **0.058→0.034**. Cost (the mirror over-call the DNA
hyperprior fixes in Phase 2): zero-gDNA-unstr **0.001→0.393**, gdna5-verystrong **0.035→0.364**, and the
`test_antisense_intronic[ss_0.65]` scenario over-calls. 226 calibration tests pass (1 pre-existing known-red).

**The v0.7.1 baseline is the proof this is a restore, not an invention.** Released v0.7.1 fit its (KDE) gDNA
prior on the **deconvolved** gDNA density in a **2-pass** structure — Pass 1 with *no* density prior (strand +
anchored global + floor + messages) → train the prior on Pass 1's peeled belief → Pass 2 re-solve. That is
exactly §6 steps 4–6, and the WIP branch now matches it on **substrate + stage** (Mistake 1 undone). The WIP
improvements worth keeping — the resolved §10 reference; ρ_bg as the aggregate background; the σ²_transfer NPMLE
projection (better than v0.7.1's disagreement-variance σ²_imp) — stay. **What still regressed is the prior's
*representation* (§5 Mistake 2):** v0.7.1 was an **additive unit-weight KDE** with a **separate capped
product-of-experts background** (enriched mode guaranteed, background can't dominate); the WIP is a
**competitive EM mixture** with a **heavy `n_regions`-weighted in-mixture background** (enriched mode starved).
Restoring the KDE's representational choices is the remaining Role-B work — §9.

---

## 9. The phased plan (from here to a finished calibration)

Each phase is A/B-gated on the flagship (unstranded+capture) **plus** a zero-gDNA control **plus** a stranded
control — never a single aggregate (`oracle_and_benchmarking.md`). No new magic numbers without discussion.

**LANDED (Phases 1–2).** The two NPMLE roles are separated (`enrichment_prior` = Role A / σ²_transfer only;
never in ψ); the initial solve is prior-free of the DNA prior (strand + Beta(½,½) reference + ρ_bg + messages);
and the gDNA hyperprior (Role B) is fit on the **deconvolved** ρ_g in a 2-pass refit (`_fit_gdna_hyperprior`,
`calib_refit_iters`). Mistake 1 is undone. The message-propagation correctness fix
([`message_absorption_fix.md`](message_absorption_fix.md) — the splice-junction RNA-imputation
false-positive) remains a standing correctness item and should land alongside Phase A (a smeared pass-0 belief
degrades any Role-B fit).

**The strategy from here (owner directive 2026-07-19): OPTIMIZE + DEBUG the proven design toward a plateau,
THEN add something new.** The shipped v0.7.1 KDE *worked* (kept the enriched mode); the redesign regressed on
**representation** (§5 Mistake 2). So the plan is to **restore the KDE's structural choices** in Role B — a
debugging exercise against a known-good target — not to invent new machinery. New mechanisms (e.g. the
enrichment-skeleton) come only *after* the restore plateaus below target.

**Phase A — restore the gDNA-hyperprior representation (THE live crux).** In leverage order, each step
A/B-gated on {enriched-node recovery, `gdna_none` FP, stranded no-regression}:

- **A1 — constrain the background (biggest lever).** The current aggregate background enters the mixture with
  EM weight `= n_regions`; on real data (100k+ intergenic/intron regions) it *will* dominate and silence every
  enriched mode. Move to the KDE's discipline: the background is a **separate, weak, CAPPED expert** (v0.7.1
  capped it at **one pseudo-observation**) combined by **product-of-experts** (so its absolute height cancels
  across the composition grid and cannot lower an enriched node's value), with the background/floor nodes
  **excluded** from the enriched-structure model. *This is the owner's "constrain the background so it doesn't
  blow up," and it is the single highest-leverage change.*
- **A2 — stop the enriched mode being competed away (the occupancy question).** The NPMLE's normalized EM
  mixture lets the depleted majority win the shared weight and **starve** the enriched components. Move toward
  an **additive / occupancy-preserving** representation (each node deposits its own kernel, KDE-style; or an
  EM with a per-component occupancy floor) so an enriched population of even a few dozen nodes keeps a bump.
  **The Bayesian tension (owner):** class occupancy *is* legitimate prior information — genuinely rare
  enrichment *should* demand per-node evidence — but the KDE expresses that correctly, as the **strength of a
  weak, capped background prior**, never by drowning the enriched mode's **height** in a normalized height
  contest. Level the *representation*, and let the (capped) background carry the honest prior-odds.
- **A3 — drop precision weighting (unit weight).** v0.7.1 used **unit weight** deliberately: precision
  correlates with density (the depleted floor is low-precision), so precision-weighting down-weights the whole
  floor and biases the distribution's *shape*. The NPMLE's belief-width weighting only *broadens* likelihoods
  and does not rescue the starved enriched mode. Match the KDE.
- **A4 — data-driven bandwidth with a per-node Poisson-noise floor.** Replace the fixed `npmle_bandwidth=0.15`
  on the **gDNA prior only** (never `enrichment_prior`) with v0.7.1's Silverman-on-Kish-count + noise-floor
  rule (`h = max(silverman, weighted-median per-node log-density noise)`). Bandwidth is a *secondary* lever
  here (the representation A1–A3 dominates), but an adaptive `h` is the honest home of the estimation width,
  and the **final selector must be set on real data** (real capture chemistry sets the enriched-mode widths).

**Phase B — real-data validation of the imbalance regime.** The CI sim (`ambig_dense_10mb`) is
gene-*dense* — it has ~216 intergenic vs ~868 expressed nodes, so it **under-represents** the intergenic/intron
flood that dominates real libraries by orders of magnitude. The background-cap strength (A1) and the bandwidth
selector (A4) therefore **cannot be finalized on the sim**; validate on real caches (LBX0190 cfRNA, MO_3005 in
`~/Downloads/rigel_runs`). This is where "how much do we weight the background" gets its real answer.

**Phase C — reconcile the two solvers (engineering debt).**
Address the 256-vs-60 grid split and the θ-endpoint quadrature (Gauss–Jacobi). Goal: single, consistent
resolution across node classes without the memory blow-up. Re-derive/retire the known-red test.

**Phase D — validation & finalization.**
Golden regen. Fold this document's design into `CALIBRATION_ARCHITECTURE.md` where it belongs and retire the
archived roadmaps.

**Only-if-we-plateau (not before):** if the restored KDE representation (A1–A4) still leaves genuinely
information-starved enriched islands below target on real data, *then* consider a new channel — the
enrichment-skeleton (source the enriched mode's existence/occupancy from the independent total-density
landscape, anchored at `ρ_bg·e^δ`; `gdna_projection_derivation.md` §C) or fragment-length per-node signal.
Do not reach for these until the restore is measured and plateaued.

**The honest ceiling (state it, don't chase it):** a fully-unstranded (κ=½) node that is an enriched,
spliced-free island with depleted neighbours is at the *identifiability ceiling* — strand is silent, messages
have no enriched-level neighbour to borrow from, and the gDNA hyperprior can only supply the background/stratum
level. Such nodes will not fully solve without strand or a new channel (fragment length, currently only an
aggregate 6-pool statistic — not per-node). This is a property of unstranded+capture data, not a bug. Do not
bend the design to chase it.

---

## 10. Document structure (the clean set)

**Keep — the current, converged references:**

| doc | role |
|---|---|
| **`CALIBRATION_MASTER.md`** (this) | the big-picture reference + roadmap + phased plan |
| `CALIBRATION_ARCHITECTURE.md` | the count-zero-information invariant (§2 here) |
| `reference_prior_derivation.md` | the resolved reference + 1-D/2-D reconciliation (reviewer-facing) |
| `background_reference_derivation.md` | ρ_bg, the aggregate background DNA level (reviewer-facing) |
| `transfer_variance_formal_derivation.md` | σ²_transfer via NPMLE projection = Role A (reviewer-facing) |
| `oracle_and_benchmarking.md` | the trusted evaluation path |

**Archive — superseded by this document** (move to `calibration/archive/`; provenance, not current): the
per-phase roadmaps and implementation plans now folded here (`calibration_roadmap.md`, `npmle_roadmap.md`,
`prior_ramp_and_bp_roadmap.md`, `background_reference_implementation_plan.md`), the message-layer WIP series
(`message_model_derivation.md`, `dof_pie_*`, `mature_crossing_gate.md`, `nascent_rna_sourcing_regression.md`,
`boundary_spliced_channel_design.md`, `node_prior_design.md`), the earlier status notes
(`calibration_selfsolve_status.md`, `calibration_bp_theory.md`, `calibration_prior_production_reference.md`), and
this session's investigation notes (`phase4_investigation_findings.md`, `flagship_prior_asymmetry_diagnosis.md`,
`refit_loop_study_findings.md`, `sigma2_transfer_derivation.md`, `npmle_projection_variance_design.md`,
`SESSION_HANDOFF.md`, `CLEANUP_LOG.md`). Their durable findings are captured in §5–§9 above; the files remain in
git + the archive for provenance.

*(Executed 2026-07-18: 22 docs moved to `calibration/archive/`; `docs/README.md` updated to this clean set.)*
