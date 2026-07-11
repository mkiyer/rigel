# Phase-2 gDNA-Density Prior — Design Review

**Status:** Design proposal for **third-party review**. Not shipped. Supersedes the *framing* (though
not all of the machinery) of [`PHASE2_gdna_mixture_fit_design.md`](PHASE2_gdna_mixture_fit_design.md).

**How this document was produced.** The empirical numbers below come from prototypes run through the
*production* `calibrate()` on production-faithful synthetic toys (the baseline was reproduced
byte-for-byte, so the harness is faithful). The design was then put through three independent
adversarial critiques (Bayesian/statistical, adversarial-skeptic, simplicity/no-magic-numbers). Their
findings — including a **reproduced fatal failure** and two **contradictions with this repo's own prior
design docs** — are reported honestly in §9–§10.

**Reader's guide (TL;DR).** All three reviewers **endorse the direction** — replace the *product* of two
priors with a *mixture*, remove the double-counting, and treat "background" and "enrichment" as
different estimands needing different estimators. None of them endorses **Candidate A as specified**:
its exposure-weighting, its background-width choice, and its intron-pooling each need to be resolved,
and it fails the one case the redesign exists for (a strand-flat, genuinely-enriched AMBIG node under
capture — §10). Read §6 for the candidate, then §9–§10 for why it is not yet ready, and §8 for the
alternatives (a minimal-delta refactor of the shipped code; a single-KDE-with-super-kernel unification;
and an empirical-Bayes g-model that must be reconciled with an architectural invariant first).

---

## 1. Overview

Rigel is a Bayesian RNA-seq quantifier that jointly separates mature mRNA, nascent RNA, and genomic-DNA
(gDNA) contamination. A **calibration** stage runs before per-locus quantification and deconvolves each
genomic **node** — a region or a boundary in a bipartite chain over the genome — into a gDNA fraction
`f_g` versus RNA. The deliverable of calibration is a small set of per-locus priors; the accuracy of the
whole tool rests on getting `f_g` right, especially for nodes whose own data cannot determine it.

A node's `f_g` is set by exactly three sources: an intrinsic **strand likelihood**, cross-node
**belief-propagation messages**, and a **global gDNA prior**. This document concerns the third. For two
node classes — **AMBIG** nodes (two overlapping opposite-strand genes, whose per-strand counts are
balanced and thus strand-indistinguishable from gDNA) and **pure-RNA** nodes (which cannot reach
`f_g = 0` because Beta-Binomial overdispersion leaves a residual) — the strand evidence is flat and the
global prior is decisive.

**One-sentence proposal (Candidate A, under review):** replace the current *product* of a weak
exposure-pooled floor and a unit-weight kernel density estimate (KDE) with a single *mixture*
`P(log ρ_g) = π_bg · N(log ρ_bg, s_bg²) ⊕ π_enr · KDE_exon`, combined by log-sum-exp, where the
background component is an exposure count-aggregate and the enriched component is a per-node KDE over
exon nodes only — each node used exactly once. This document presents that candidate, the empirical
evidence for it, and, prominently, the unresolved defects that three independent critiques identified.

---

## 2. Background & Problem

**The count-zero-information principle.** Rigel's calibration architecture holds that a raw fragment
**count** carries *zero* intrinsic information about whether a node is gDNA or RNA — a locus with 1,000
fragments is not more "gDNA" than one with 10. The only intrinsic gDNA/RNA signal is the **strand tilt**
of the per-strand counts (gDNA is unstranded; RNA is stranded), and the count enters *only* as the
overdispersed Fisher information (statistical power) of that strand likelihood. No code path may let a
count *vote* on composition. This is an architectural invariant, not a preference; it constrains what
any prior estimator is permitted to do (see §8).

**What a node is.** Calibration builds a bipartite chain of **region** nodes (contained genomic
intervals: intergenic, intron, exon) interleaved with **boundary** nodes (junctions). Each node carries
an unspliced fragment **mass** `M`, a gDNA **effective length** `E` (the FL-marginal length over which a
gDNA fragment could be observed), and, after solving, a gDNA fraction `f_g`. The per-node gDNA
**density** is `ρ_g = f_g · M / E`.

**Why AMBIG and pure-RNA nodes need the prior.** An AMBIG node's balanced strand counts are *identical*
to gDNA's — its strand likelihood is flat, so `f_g` is determined by messages and the prior alone. A
pure-RNA node *should* solve to `f_g = 0`, but Beta-Binomial overdispersion leaves an irreducible
residual — the **"od-floor"** — pinning it near `f_g ≈ 0.01`. Both classes need the global prior to
supply a credible "gDNA is scarce here" anchor.

**The Phase-2 flow.** (1) Solve all nodes with an *extremely weak* prior. (2) **Fit** `P(log ρ_g)` on
the solved single-strand nodes. (3) **Re-solve** all nodes (including AMBIG) under the fitted prior.
Under capture-OFF, `P(ρ_g)` is ~unimodal at background; under capture-ON, it is bimodal — a depleted
background mode plus an enriched-exon mode.

**The real-data landscape.** A transcriptome has ~75,000 genes and >500,000 transcripts, only a minority
expressed/captured in any sample. The vast majority of nodes are background (intergenic + intronic). The
background mode is over-determined; the enriched mode is a fragile minority. Any estimator must be
dominated by the background where it should be, without erasing the enrichment where it exists.

---

## 3. Current Solution (shipped)

The shipped pass-2 global prior is a **sum of two log-priors** applied to each per-node solve over the
density grid:

```
global_lp = _global_logprior(...)      # "floor" (composite; see below)
          + _kde_logprior(...)         # trained KDE
```

**The floor (`_global_logprior`).** This is *not* a single weak Gaussian; the honest description (from
`bp_solver.py`) is a composite of a flat log-density baseline, an **ê(z)** enrichment-transfer blend
(`fit_enrichment_transfer`, a `MonotoneMean`/`MonotoneVarMean` apparatus with a continuous
slope-evidence weight), and a **`_floor_estimate`** depleted-region override. `_floor_estimate` is itself
an exposure-pooled background: `ρ_floor = (1 + ΣG)/ΣE_tot` over the intergenic+intron floor mask, with a
between-node log-density spread `s2_floor` (raw spread minus the Poisson floor) and a strand discount
`w_str = (2κ−1)²`. Under the current `_PHASE1_PRIOR_FREE` flag, the total effective precision is **capped
at one pseudo-observation** (`_GLOBAL_STAB_PREC = 1.0`) — the mechanism enforcing count-zero-info at the
prior level.

**The KDE (`_kde_logprior`).** A Gaussian-kernel KDE over `log ρ_g`, fit with **unit weight** on
**all** solved single-strand nodes (zero-count intergenic/intron floored at the min-observable `1/E`;
tiny nodes with `E <` gDNA mean FL excluded), evaluated on the per-node density grid via
`ρ_g = f_g · M / E`.

**Two-pass flow.** Pass 1 uses `_global_logprior` alone. Pass 2 *adds* the trained KDE on top. Constant
inventory: `_GLOBAL_STAB_PREC = 1.0`, the KDE bandwidth (`gdna_prior_bandwidth`), the `E <` FL
exclusion, and the `1/E` floor.

---

## 4. Issues with the Current Solution

**(1) "Add" is a product-of-experts.** Adding two log-priors gives `P ∝ P_floor · P_KDE`. Semantically
this is an **AND**: a node's density must be plausible under *both* the floor Gaussian and the KDE
simultaneously. For a genuinely multimodal population — where a node is at background **OR** enriched —
the coherent object is a **mixture** (a weighted **SUM**, an OR), evaluated by log-sum-exp. The product's
tail is the *sum* of two quadratic penalties, so an enriched node is penalized for being implausible
under the background, and vice versa — exactly the wrong shape.

**(2) Double-counting.** The intergenic/intron nodes inform the prior **twice**: once through the
exposure-pooled floor (`_floor_estimate`) and again as KDE training nodes. The additive combination
counts their evidence twice.

**(3) The floor is load-bearing — the CLUE.** Using the KDE **alone** in pass-2 (replacing the floor)
*fails* the AMBIG-in-gDNA=0 test: the baseline reproduction gives `f_g = 0.376` (a false gDNA read).
Adding the weak floor back fixes it (`0.0052`). That a *weak* floor is load-bearing is the diagnostic
clue. Baseline control numbers, reproduced byte-for-byte through production `calibrate()`: **S1 (AMBIG,
gDNA=0) → 0.0052; S2 (AMBIG, gDNA=200) → 0.4577; S4 (pure-RNA) → median 0.0102.** These match the
documented anchors (0.005 / 0.458).

---

## 5. The Insight

The background gDNA density is fundamentally an **exposure-aggregate**. The signal "gDNA is globally
scarce" comes from empty intergenic regions: `G ≈ 0` fragments over enormous `E` → overwhelming evidence
that `ρ ≈ 0` — and that evidence is carried by their **exposure**, not their node count. A unit-weight,
node-count KDE structurally **cannot** reconstitute it. Three measured findings support this (all through
production `calibrate()`, pass-1 weak-floor belief):

- **(a) Zero-count nodes carry only `1/E`.** On the capture-ON zero-inflated toy, every zero-count
  background node had `ρ_g = 1/E` *exactly* (`|log ρ − (−log E)| = 0`, `corr(log ρ, −log E) = +1.000`).
  Their `log ρ` spanned 1.72 nats — driven **entirely** by `E ∈ [1746, 9746]` bp, not by rate. A
  zero-count node's density is an artifact of effective length.

- **(b) Only pooling recovers the background.** Apples-to-apples on the same gDNA eff-length basis
  (gDNA=0.3 present, 85 background nodes, ΣE = 475,432): pooled-oracle `log ρ_bg = −10.87`; pooled-solved
  `= −11.44` (agree to ~0.6 nats). The **unit-weight KDE top mode = −3.13** — **+7.75 nats (~2,300×) too
  high** — because it tracks the exon per-node cluster and out-votes the sparse background. Even a
  background-only unit KDE modes at `−8.63`, ~2.2 nats above the pooled rate, because the per-node values
  *are* the `1/E` floor, not `ΣG/ΣE`.

- **(c) The od-floor makes a spurious enriched mode — but not the way we first thought.** The original
  framing ("false gDNA count ≈ 0.01·M, constant fraction") was **refuted**: on the pure-RNA sweep (E
  fixed 800 bp, M spanning 35×), `gCOUNT ~ M^+0.22` (constant-fraction model R² = −10.9; constant-count
  model R² ≈ 0). The od-floor is a **near-constant COUNT** (~3–8 fragments), so `f_g ~ c/M` — **largest
  at small-mass exons**. The downstream *consequence* the insight rests on is confirmed: pure-RNA
  per-node densities sit well above the true (~0) background and form a spurious "enriched" mode the
  unit-weight KDE faithfully tracks (KDE mode −6.26 near the exon od-floor cluster −6.07, not the
  intergenic floor −8.66).

**Conclusion:** background and enrichment are **different estimands requiring different estimators**.
Background = an exposure count-aggregate (best **pooled**). Enrichment = a per-node distribution (best
**KDE'd**). The floor is not a patch; it supplies the one quantity a node-count KDE cannot.

---

## 6. Proposed Solution — Candidate A (exposure-weighted mixture)

Fit, from the pass-1 solved single-strand nodes partitioned by node class:

- **Background component** over intergenic+intron nodes:
  `ρ_bg = (1 + ΣG_bg) / ΣE_bg`, `s_bg² = 1/(1 + ΣG_bg)`, component = `N(log ρ_bg, s_bg²)`.
- **Enriched component** over exon nodes only: `KDE_exon` = unit-weight Silverman Gaussian KDE.
- **Mixture:** `P(x) = π_bg · N(x; log ρ_bg, s_bg²)  ⊕  π_enr · KDE_exon(x)`, combined by **log-sum-exp**,
  with `π_bg = ΣE_bg / (ΣE_bg + ΣE_exon)` (by **exposure**).

The weak floor is **disabled** — the mixture is the *sole* global prior. Each node is used **once**
(background nodes only in the pooled Gaussian; exon nodes only in the KDE) → no double-count.

**How it fixes each issue.** (1) Product → mixture: log-sum-exp gives the correct OR combinator. (2)
Double-count → removed by the disjoint partition. (3) The load-bearing floor is *promoted* to an
explicit, exposure-dominant background component, so it no longer needs to be silently added on top.

**Why it self-calibrates (measured).** S1 (no gDNA): `ΣG_bg = 49.3 / ΣE_bg = 20469 → ρ_bg = 0.0025`,
`π_bg = 0.932`; the exposure-dominant background pulls the flat-strand AMBIG node (`fg_strand = 0.376`)
down to `0.0019`. S2 (real gDNA): `ΣG_bg = 4554 / 19508 → ρ_bg = 0.233`; the *same* estimator now sits
at the elevated background, *agrees* with the AMBIG strand (`fg_strand = 0.458`), so `fg_final = 0.4577`.
The identical mechanism gives ~0 without gDNA and reads gDNA when present.

**od-floor as a consequence:** pure-RNA nodes are pulled toward the exposure-dominant near-zero
background → S4 `f_g = 0.0026`.

**Critical implementation requirement (a genuine finding).** The mixture must be evaluated
**analytically** in `logpdf`, so the background Gaussian's `−½((x−μ_bg)/s_bg)²` tail persists at all `ρ`.
The shipped grid-interpolated `logpdf` flat-extrapolates outside the training grid; with flat
extrapolation the AMBIG node (which sits *above* the grid) sees a constant prior, zero downward gradient,
and **S1 fails at 0.542**. A drop-in that swaps only the fit while keeping the flat-tail `logpdf` still
fails S1. The downward anchor lives in the background Gaussian's tail, which grid discretization silently
destroys. (This finding turns out to be central — see §9.4.)

---

## 7. Empirical Validation

All via production `calibrate()` on tiny production-faithful toys. Targets: S1 < 0.08, S2 > 0.2, S4 → ~0,
S3 bimodal.

| Scenario | Baseline (floor + KDE, product) | **Candidate A (mixture)** | Candidate B (single exposure-KDE) | Target |
|---|---|---|---|---|
| **S1** AMBIG, gDNA=0 | 0.0052 | **0.0019** ✓ | 0.3756 ✗ | < 0.08 |
| **S2** AMBIG, gDNA=200 | 0.4577 | **0.4577** ✓ | 0.2339 (weak ✓) | > 0.2 |
| **S3** capture-ON bimodal | not run | **BIMODAL** ✓ (bg mode log −11.44 vs enriched −3.0…−3.9) | not run / unimodal | bimodal |
| **S4** pure-RNA `f_g` | 0.0102 | **0.0026** ✓ | 0.0073 (unchanged) | ~0 |
| double-counts | **yes** | **no** | no | — |
| resolves all | yes (via load-bearing floor) | **yes** | **no** | — |

**Honestly resolved:** S1–S4 all pass for Candidate A on these toys, with the double-count removed and a
true mixture combinator. The baseline was reproduced exactly (harness is faithful). Candidate B is a
**clean negative** — a density-only KDE cannot supply a *fraction* anchor (a small-`E` AMBIG node
matching a background *density* implies a *large* fraction), which is decisive evidence *for* separating
the estimands.

**Not resolved / not tested:** (i) The toys are tiny (`n_bg = 10, n_exon = 2` in the AMBIG toys) — a
2-sample KDE and an `s_bg` at `n_exon = 2` are not distributions. (ii) **No end-to-end benchmark** was
run; the 16-condition net-flow suite that gates releases was not exercised, so whether Candidate A moves
the incumbent's known capture-exon residual is unknown. (iii) S2 has **no oracle** `f_g` — "passes > 0.2"
and "matches baseline 0.4577" are not correctness claims; the baseline is the thing being replaced. (iv)
No toy exercises a **strand-flat, message-isolated, genuinely-high-gDNA AMBIG** node — the case the
redesign exists for (see §10).

---

## 8. Alternatives Considered

**Candidate B — exposure-weighted single KDE (no separate background).** One KDE over all nodes, kernels
weighted by `E`. Predicted and confirmed to **FAIL** (S1 = 0.3756). The KDE mode (−1.421) actually
landed near the pooled background (−1.455), yet the AMBIG node still solved high — because a *density*
target does not translate to a *fraction* target for a small-`E` node. This is the strongest evidence for
the mixture.

**Empirical-Bayes Poisson g-model / NPMLE (Kiefer–Wolfowitz).** Model `G_i ~ Poisson(ρ_i · E_i)` with
`ρ_i ~ P(ρ)` estimated by nonparametric MLE. **Pros (per the statistical reviewer):** handles count
precision natively (a zero over huge `E` pushes mass to ρ≈0 with no `1/E` hack), no bandwidth knob,
exposure enters automatically via the Poisson mean, and the background atom emerges with the correct
width. **Cons:** (a) it estimates `P(ρ)` from **raw counts** `{G_i, E_i}` — which **violates the
count-zero-information invariant** (the count votes on composition). It is legal *only* if run on the
deconvolved `ρ_g = f_g·M/E` (where `f_g` came from strand+messages), not raw `G` — a reconciliation the
doc must make explicit before the g-model can be adopted. (b) It does **not** escape the od-floor bias —
it fits the same solver-biased counts and reproduces the same spurious low-ρ mode. So the g-model is
cleaner on several axes but does not remove the need to fix the substrate.

**Minimal-delta (recommended by the simplicity reviewer).** Reuse the shipped `_floor_estimate`
*verbatim* as the background component (it already is exposure-pooled `(1+G)/E` **with** a between-node
spread `s2_floor` **and** a strand discount `w_str` that Candidate A discards), restrict the KDE to exon
nodes, and swap `add` for exposure-weighted `logsumexp`. This introduces **one** new knob (the mixture
weight) versus Candidate A's three (partition rule, π scheme, `s_bg`), and inherits the strand discount
for free.

**Single-KDE-with-a-super-kernel.** Place the pooled background as one high-exposure, narrow-bandwidth
kernel *inside* one variable-bandwidth KDE. One object, one combinator; the partition survives only as
"which nodes pool into the super-kernel." The most elegant unification if a single machine is wanted.

---

## 9. Open Questions for the Reviewer

1. **Mixture weight.** Is `π_bg` by **exposure** correct, or should it be **node-count** or
   **class-conditioned**? An AMBIG node *is* an exon (two overlapping genes); its correct
   class-conditional prior may be the *enriched* distribution, not the exposure-dominant background.
   `π_bg` swings 0.60–0.99 across panel fractions and node-count weighting flips the sign — this is *the*
   parameter, not a detail. Should it be a per-node responsibility rather than a global scalar?
2. **Contradiction with prior art.** Exposure-weighting directly contradicts this repo's own
   `PHASE2_gdna_mixture_fit_design.md` §1 ("do **NOT** exposure-pool") and §8e ("precision weighting is
   **WRONG** — use unit weight"), and pooling intergenic+intron contradicts
   `PRODUCTION_DESIGN_gdna_mixture_prior.md` §3.4 (floor = intergenic **minimum**, not the exposure-pooled
   **mean**). Which is right, and why? (One reconciliation: the earlier docs argue against exposure-weighting
   the *enriched* component, which Candidate A does *not* do — it keeps `KDE_exon` unit-weighted and only
   uses exposure for the *background* pooling and the mixing weight. This distinction must be made
   explicit and defended, not elided by the prose "weighted by exposure".)
3. **Background width.** `s_bg² = 1/(1+ΣG_bg)` is the *sampling* variance of the pooled log-rate (how
   well we know the **mean**) — not the **population spread** of background nodes. On real data `ΣG_bg` is
   in the thousands → `s_bg → 0` → a near-delta spike with effective precision ≫ 1 pseudo-obs, violating
   the count-zero-info cap. Should the width instead be the existing `s2_floor` population spread, and
   should the one-pseudo-obs cap be reapplied?
4. **The nonparametric claim.** The load-bearing S1 fix is *entirely* the analytical background-Gaussian
   **tail** (flat-tail fails at 0.542). Should we stop calling this "nonparametric," specify the tail
   explicitly (Gaussian vs Student-t vs slab), and defend its strength (`= 1/s_bg²`)?
5. **od-floor placement.** Is the od-floor a **substrate bug** to fix in the strand solver, or should the
   prior absorb it? Candidate A's S4 pass is partly an artifact of an uncapped tail that will over-pull
   genuinely-low-but-nonzero nodes once the solver is fixed. Note the correct model: near-constant COUNT,
   `f_g ~ M^−0.78`, worst at small-mass exons — not `0.01·M`.
6. **The g-model vs the invariant.** Can an NPMLE be reconciled with count-zero-information (run on
   deconvolved ρ, not raw G)? If yes, is it the target architecture? Either way, it inherits the od-floor
   bias — so the solver fix is orthogonal.
7. **Real-data survival.** Does the scheme hold at 75k genes / capture panels / extreme zero-inflation?
   None of S1–S4 exercises this, and the release benchmark was not run.

---

## 10. Risks & Failure Modes

**Enriched-AMBIG crush (reproduced, FATAL for the stated purpose).** The adversarial reviewer built the
canonical hard case in production `calibrate()` — an unexpressed, captured, opposite-strand overlap =
pure-gDNA AMBIG, truth `f_g = 1.0` (`g=20, r=0`). Result: **stranded → 0.624 (err 0.376); unstranded →
0.376 (err 0.624)**, with a prior pull toward background of **+56 nats**. This is *exactly* the case the
redesign exists for, and it fails hardest, because Candidate A's only high-density mode is the
`KDE_exon` of *expressed* exons — a pure-gDNA enriched AMBIG node has **zero representation** in the
training substrate. It is invisible in S1–S4 because those AMBIG nodes are strand-informed or
message-supported. This bites on any capture panel with real gDNA over unexpressed /
opposite-strand-overlapping regions (ubiquitous in exome/panel data).

**Background anchor inflated by introns.** Pooling introns into `ρ_bg` folds nascent-bearing regions into
the background (introns read as real unspliced density: pooled intergenic+intron = 0.23 in S2, "not near
zero"). On deep total-RNA (nascent-rich introns) this raises `ρ_bg` above the true intergenic floor,
weakening the downward pull the anchor exists to supply. The shipped `_floor_estimate` strand-discounts
introns by `(2κ−1)²`; the Candidate A prototype pools raw `G/E` with no discount.

**Background width mis-scales in both regimes.** At low count `s_bg` is so wide it barely anchors (the
+56-nat crush happened *with* `s_bg = 0.71`); at high count it collapses to a brittle delta that crushes
genuinely-elevated real background (CNV, GC bias). Sampling precision and population spread are conflated
into one symbol.

**Class partition is a hard classifier under a continuum.** Background-vs-enriched is a continuum under
CNV and partial-probe-overlap; routing by `node_kind` mis-files a CNV-amplified intergenic region into
the low background and a depleted exon into the enriched KDE — wrong exactly where copy-number or
probe-efficiency variation matters.

**Self-reinforcing od-floor loop.** KDE-ing the od-floor exon cluster verbatim (mode −6.26, tracking the
artifact at −6.07) means the trained prior *feeds the strand-solver artifact back* into pass-2. The
mitigation is to subtract the expected Beta-Binomial residual count before fitting — a substrate fix that
shrinks the prior's job.

**Bottom line.** The *direction* is endorsed by all three reviewers — mixture over product, double-count
removal, "different estimands, different estimators." Candidate A's specific parameterization is **not**
yet defensible: the exposure weight, the `1/ΣG` width, and the intergenic+intron pooling each need to be
resolved (§9), the enriched-AMBIG crush needs a training-substrate fix or a per-node responsibility
weight, the od-floor should likely be fixed in the solver, and a real-data benchmark is mandatory before
this ships.
