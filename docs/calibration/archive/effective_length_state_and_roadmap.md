# Effective-length contraction: state of the tool + roadmap

_Last updated 2026-07-04._

This document records the design, current production state, measured accuracy, remaining error, and the
concrete path to further improvement for Rigel's **gDNA effective-length contraction** — the piece of the
calibration→EM prior that stops hybrid-capture-enriched genomic-DNA (gDNA) contamination from leaking into
RNA abundance.

## 1. The problem

Under hybrid capture, gDNA is not uniform: capture probes concentrate gDNA on the targeted exons and deplete
it elsewhere. The per-locus EM competes a **gDNA component** against the transcript components for each
fragment. A component's competitiveness scales with its per-base rate `abundance / effective_length`. If the
gDNA component uses its full genomic span as its effective length, it is under-weighted exactly where it is
concentrated (the exons) — so enriched gDNA fragments leak to mature/nascent RNA. This is the
**capture-dependent gDNA→RNA leak**; it grows with capture strength and gDNA level, and it was the dominant
residual error in the deconvolution.

## 2. The architecture (production)

The gDNA component's effective length is **contracted** toward the enriched footprint:

```
eff_gdna = θ_g / ρ*        (θ_g = total gDNA mass over the locus's nodes; ρ* = a reference density)
```

Equivalently, per node `ℓ_n = m_n/ρ*` and `eff = Σ ℓ_n` — each node's gDNA mass expressed as the length it
would occupy at the reference density `ρ*`. Two structural pieces make this correct:

**2a. Node model + transcript→node mapping.** The genome is partitioned into nodes: **regions** (contiguous
same-signature segments) and **boundaries** (the seams between them). Each node has a gDNA mass `m_n` and an
**effective support** `S_n` (the FL-marginal sampling support: `E[max(0,L−ℓ)]` for a contained region;
`½(E[min(ℓ,L_r)]+E[min(ℓ,L_{r+1})])` for a pooled seam). A component's node set is exactly the nodes it
occupies **contiguously** (`calibration/capture_eff_length.py::_transcript_node_incidence`):

- **gDNA** (a contiguous genomic interval) → ALL of the locus's nodes, *including the two outer boundaries*
  that bookend the locus (exclusively-gDNA crossing fragments live there), but **not** the intergenic
  *regions* outside the locus (dropped by the per-locus projection).
- **nascent RNA** (a single unspliced span) → all regions + all interior boundaries in its span (introns
  included), no outer boundaries.
- **mature mRNA** (spliced) → its exon regions + only the boundaries it crosses *without a splice*
  (exon-interior signature-change boundaries); its introns and splice-junction boundaries are excluded.

Fixing the boundary handling (the pre-redesign code excluded *all* boundaries from transcripts) was itself a
large win: on the target case it dropped the mature leak from ~130k to ~41k. The node model is verified to
match this spec exactly (single-transcript locus: gDNA eff = nascent eff + the 2 outer-boundary supports, a
~2.5% gap).

**2b. Reference density `ρ*` (production = `contained`).** `ρ* = G_c/E_c`, the **mass-weighted CONTAINED
(exon) density** (`E_c` = the contained-node Laplace-smoothed IPR). This reads the gDNA density where mature
RNA competes for the ambiguous contained fragments. Properties (adversarially verified, 4/5 invariants; see
below): factor-1 bit-identical under a uniform (capture-off) field; monotone contraction toward the exon
footprint under capture; bounded; no magic numbers (Laplace `+1`, `w=C/(C+1)` shrinkage, `½`-averaged seam
support — all parameter-free); FP-safe on zero-gDNA loci (`G_c=0 ⇒ factor 1`).

`calibration/priors.py::assemble_priors` computes the gDNA eff-len; `capture_eff_length.py` applies the same
`ρ*` per transcript so gDNA and a full-span nRNA sharing a locus land at the same eff-len (up to gDNA's two
outer boundaries) — parity by construction, not a special case.

## 3. Current accuracy (production `contained`)

Release benchmark (`quick_3to1_5mb`, 16 conditions, net fragment-flow metric), vs the pre-redesign baseline
`d011b487`:

| metric | baseline | **contained (production)** |
|---|---|---|
| target (gdna300, ss0.99, capture-on) | 171,268 | **130,441** (−24%) |
| capture-ON total \|net\| | 784,868 | 640,547 (−18%) |
| capture-OFF total \|net\| | 152,851 | 157,082 (~flat; near-pristine) |
| mature leak (target) | ~130k | **~41k** |

Verification (6-agent adversarial workflow): factor-1 exact across 20k random layouts (relerr ~1e-15);
monotone contraction to the exon footprint; zero magic numbers; the gDNA/nascent ρ*-identity + node-set model
confirmed. One low-severity edge fixed (a `_GDNA_EFF_LEN_FLOOR` that could exceed span on sub-basepair spans
— clamped; no real-data effect). 1084 tests pass; capture-off bit-identical.

## 4. The reference-density study — where the remaining error is

`ρ*` is the **dominant lever**: sweeping it from the depleted extreme to the enriched extreme moves the
capture-on leak across an order of magnitude (~530k → ~17k on the sweep). The study (all env-gated, default
unchanged) established:

**4a. The eff-len and calibration accuracy are COUPLED.** With *perfect* per-node gDNA masses (an "oracle
calibration" that replaces the deconvolved split with the ground truth), the production `contained` reference
*over*-contracts and gDNA over-steals mature (−58k on the target). The reason: today's deconvolution
**under-calls captured-exon gDNA by ~2×** (observed peak density ~55 vs true ~103), and `contained`'s
aggressive contraction was silently *compensating* for that under-call. So `contained` is right for
**today's (observed) calibration**, but it is not the reference we'd want once calibration is accurate.

**4b. On accurate calibration the correct reference is the ENRICHED MODE, milder.** On the oracle substrate,
mature-balance (net mature↔gDNA ≈ 0) is achieved by a reference at the **enriched shoulder/mode** of the
per-locus node-density distribution — robustly across capture strength × gDNA level — not the mass-weighted
mean (over-contracts), the peak/max (over-contracts hard), or the support-weighted mean (under-contracts).

**4c. `kmeans_midpoint` — the magic-free enriched-mode estimator (optional path, the target).** A log-space,
support-weighted **1-sided k-means (k=2)**: two centroids fit to the log-densities with the split at their
self-determined geometric midpoint; `ρ*` = the support-weighted geomean of the enriched (above-split)
cluster. It has **no tunable constant**, is **robust in the few-enriched/many-depleted real-panel regime**
(where a fixed quantile like p90 collapses into the depleted cluster), and is degenerate-safe (uniform → that
density). It is available in production via `RIGEL_RHOSTAR=kmeans`.

**4d. Why `kmeans` is NOT yet the production default.** On *observed* (non-oracle) calibration it is currently
**worse than `contained` and fails the capture-off robustness test**:

| scenario | contained | kmeans |
|---|---|---|
| ALL \|net\| | **798k** | 878k |
| capture-ON | **641k** | 693k |
| capture-OFF | **157k** | 185k ⚠ |
| gdna_none (FP) | **111k** | 134k ⚠ |

The failure mode is diagnostic: capture-off density is ~unimodal (uniform), so k-means **hallucinates a
two-cluster split out of Poisson noise** and contracts when it should stay at factor 1. `contained` is robust
here by construction (the mass-weighted mean of a uniform field is that uniform value → factor 1). k-means
needs a **cleanly bimodal** node-density landscape, which only emerges once the deconvolution is accurate.

**4e. The mature/nascent tension is fundamental to a scalar reference.** On perfect calibration, mature
competition (at exon nodes) balances at a *mild* reference; the nascent siphon (at intron/crossing nodes)
wants an *aggressive* reference. No single scalar `ρ*` serves both — confirming the nascent leak is a
genuinely separate competition, not something the mature eff-len should chase.

## 5. Remaining error, ranked

1. **Deconvolution under-calls captured-exon gDNA (~2×).** This is the highest-leverage remaining error. It
   forces `contained` to over-contract as compensation and blocks the principled `kmeans_midpoint` reference.
2. **Scalar-reference ceiling (~1% mature).** Even with perfect calibration and the enriched-mode reference,
   a single per-locus `ρ*` leaves ~1% residual mature over/under-correction (the balancer is an enriched
   *shoulder*, not cleanly estimable) and cannot simultaneously serve the nascent competition.
3. **Nascent hallucination / siphon.** A separate competition (intron/crossing fragments), unresolved by the
   eff-len; wants the opposite reference; hardest at low strand-specificity where an intron gDNA fragment is
   indistinguishable from an intron nascent fragment except by strand.

## 6. Roadmap (in dependency order)

1. **Calibration exon-gDNA accuracy** _(unlocks everything below)_. Diagnose and fix the ~2× under-call of
   captured-exon gDNA in the deconvolution (candidates: the strand likelihood, the density prior, the KDE
   prior interaction under capture). Success criterion: per-node gDNA densities approach the oracle, so the
   node-density landscape is cleanly bimodal.
2. **Flip the reference to `kmeans_midpoint`.** Once (1) lands, re-run the observed-calibration A/B; expect
   `kmeans` to match/beat `contained` and stay robust across capture on/off × stranded/unstranded × gDNA.
   Then make it the default. (A magic-free bimodality guard may still be wanted for the unimodal capture-off
   case — to be derived, not tuned.)
3. **Per-node local density at scoring** _(the exact fix)_. Replace the per-component scalar eff-len with the
   per-node density `ρ_n` used locally per fragment: gDNA competes at each node's own density — peak density
   at peak exons, typical at typical exons, intron density at introns — simultaneously, with no reference to
   estimate and no scalar-collapse ceiling. Touches the scoring/EM layer.
4. **Nascent hallucination** _(separate initiative)_. Address the intron/crossing gDNA-vs-nascent
   identifiability at the prior/likelihood layer.

## 7. Reproducing the study

All study machinery is env-gated with the production default unchanged (`RIGEL_RHOSTAR=contained`):

- `RIGEL_RHOSTAR=kmeans` — the optional enriched-mode reference (§4c).
- `RIGEL_ORACLE_CALIB=<oracle_calib.npz>` — inject perfect per-node masses (calibrate.py hook), built by
  `scripts/debug/oracle_calibration.py <condition_dir>`; the substrate for §4a/4b/4e.
- `scripts/debug/oracle_gdna_density.py` — per-region true gDNA density (the observed-vs-true peak diagnostic).

The contraction-curve families (linear / power / log-linear / clamped) are plotted for reference; the current
`ℓ_n = m_n/ρ*` is log-linear with slope 1. See the memory note `efflen_node_mapping_contained_rhostar.md`.
