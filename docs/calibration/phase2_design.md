# Phase 2 — count posterior variance + the FP-rate quantile knob (design, research-backed)

**Status:** dedicated design doc, 2026-06-11. Built from diagnostics (not assumption) so the behavior of
the count variance is understood before it is wired in. Supersedes the Phase-2 sketch in
`remaining_phases.md`. Diagnostics: `scripts/debug/diag_count_variance_intuition.py`
(`/tmp/count_variance_intuition.png`) and `diag_count_posterior_calibration.py`.

## 0. The two decisions this doc settles

1. **The blend weight stays `w=(2κ−1)²`** — *not* the Fisher/inverse-variance blend (empirically: the
   Fisher blend would let a confidently-biased count dominate even at high strand specificity). §3.
2. **The FP-rate quantile acts on a *node-type-aware* posterior** — trustworthy on strand-observable
   nodes (a real Bayesian posterior), deliberately weak on count-routed nodes (the count variance is
   anti-calibrated under capture, so it must not be used to *correct* the estimate, only to widen it). §4.

## 1. Completing the Phase-1 variance model (Poisson floor everywhere)

**Every count is Poisson at baseline — regions and boundaries alike.** A boundary's flux is the number
of fragments whose span covers a single position; under the generative model (fragment starts ~ a
Poisson process, lengths ~ FL) that count is a superposition of independent Poisson start-counts ⇒
**Poisson** (variance = mean ≈ `λ·mean_FL`). So the complete model is:

| node | baseline | + imputation term? |
|---|---|---|
| **observable region** (intron/intergenic) | Poisson on its own contained count | no — direct |
| **imputed region** (exon/AMBIG) | Poisson on its anchor counts | **yes** — variance~mean LOESS over the 2 anchors |
| **boundary side** | Poisson on its crossing count | no — a boundary is a direct measurement, not imputed |

The variance~mean LOESS is **regions-only** (only regions are imputed from two anchors); boundaries get
the Poisson floor alone. **Current gap:** `node_gdna_density` gives regions the floor + imputation, but
**boundary sides have no variance** — Phase-1 must add the per-side Poisson floor (`σ²_count` relative =
`1/N_crossing`) before Phase 2 consumes it. Overdispersion (NB on top of Poisson) is the *separate*
planned model; the Poisson floor is the correct baseline.

## 2. What the count variance actually looks like (intuition)

The count module's output is a **fraction** `g ∈ [0,1]`, so its variance is bounded by the Bernoulli
maximum `μ(1−μ)`. Consequences, confirmed on the complex suite (fig. panel A):

- **Confident at the extremes.** Observable gDNA-pure regions (`μ→1`) and pure-RNA exons (`μ→0`) have
  `σ_count ≈ 0` — the *fraction* is confident even though the *count* has Poisson noise. (The cap forces
  this; it is correct — a pure-gDNA region is confidently gDNA.)
- **Variance lives in the mid-fraction imputed exons** — the genuinely ambiguous nodes, largest at low
  count (sparse anchors) and shrinking with depth.
- The Poisson floor in fraction units is `≈ μ·(1−μ)·(1/N)`-scaled — small for deep observable regions
  (`1/√N ≈ 0.16` at capture depths, but capped toward 0 as `μ→1`), larger for boundaries (`1/√N_cross ≈
  0.25`, fewer crossing fragments).

**The critical property — the count variance is *not* an accuracy signal under capture.** Binning imputed
exon nodes by `σ_count` and comparing to the oracle error `|μ_count − truth|`:
`Spearman(σ_count, |error|) = −0.19` (ss0.99 cap_on), `−0.23` (ss0.50 cap_on), `+0.02` (ss0.99 cap_off).
Under capture the **most confident (smallest σ) count nodes are the most wrong** — the capture bias is a
*confident* error (two depleted off-target anchors agree on a biased-low density). So `σ_count` measures
imputation *disagreement*, which is blind to the systematic capture bias. This is the central fact that
shapes both decisions below.

## 3. The blend weight — keep `(2κ−1)²` (the Fisher blend fails)

The question was whether to replace `w=(2κ−1)²` with the inverse-variance blend
`w_fisher = I_strand/(I_strand+I_count) = σ²_count/(σ²_strand+σ²_count)`, hoping for an elegant
strand→count handoff as SS drops. **Diagnostics say no** (fig. panel B):

- `w_fisher` piles at **~0 for every SS regime**, including ss0.99 where `(2κ−1)²=0.96`. So the Fisher
  blend would let **count dominate even at high strand specificity** — the opposite of the desired
  handoff.
- Cause: the count variance is Bernoulli-capped to ~0 at confident fractions, so `I_count=1/σ²_count` is
  spuriously enormous and swamps `I_strand` regardless of SS.
- Combined with §2 (confident count = biased count under capture), inverse-variance weighting would
  **launder the capture bias** into the estimate — exactly what the decoupling was built to prevent.

**`(2κ−1)²` is bias-robust** (it weights the unbiased strand channel by its *effect size* / discriminability,
not by the deceptive count precision) and gives the correct behavior: strand governs at high SS
(`w→1`), count governs at unstranded (`w→0`). **Decision: keep `(2κ−1)²`.** The Fisher blend is recorded
here as evaluated-and-rejected so it isn't revisited without new evidence.

## 4. The FP-rate quantile knob — node-type-aware

The knob is still the right user abstraction (one parameter trading gDNA→RNA leak vs RNA→gDNA siphon).
But §2 dictates *how* it can act:

- **Strand-observable nodes** (stranded libraries): the strand BB posterior is a genuine Bayesian
  posterior with a meaningful variance → the quantile `g(q)=center+Φ⁻¹(q)·σ_strand` is trustworthy and
  uncertainty-aware (it bites on wide-posterior nodes). This is the high-value case (the critical
  stranded-capture scenario).
- **Count-routed nodes** (unstranded/AMBIG): `σ_count` is anti-calibrated (confident=biased), so the
  quantile must **not** be used to *correct* a confidently-biased node (it can't — small σ → no shift)
  and must **not** trust σ to mean accuracy. Honest options:
  - **(a)** use `σ_count` only as a *floor-widened* uncertainty (never shrinking a node toward
    confidence it hasn't earned), so the quantile at least widens uniformly;
  - **(b)** apply a **global log-odds shift** on count-routed nodes for the FP knob (the old
    `gdna_strand_llr_bias` mechanism, but as the *count-side* FP lever) — uniform, not pretending the
    per-node σ is informative;
  - **(c)** accept that the unstranded-capture bias is the documented identifiability floor and let the
    knob trade only what the calibrated (strand) posterior supports.

**Design (default-preserving):** `g(q) = clip(center + Φ⁻¹(q)·σ)`, `center` = today's point estimate,
`σ` = the combined std (strand BB ⊗ count, with the count term widened-only per (a)). `q=0.5 ⇒ g=center`
⇒ **bit-identical default** (the variance is consumed only when `q≠0.5`). `config.gdna_deconv_quantile`
(default 0.5); retire `gdna_strand_llr_bias` (subsumed). Record `q` in `CalibrationResult`.

## 5. Validation plan

- Unit: `q=0.5` reproduces the baseline exactly; `g(q)` monotone in `q`; count-routed nodes widen but
  don't sharpen.
- Sweep `q ∈ {0.5,0.7,0.9,0.95}` on the complex suite; confirm the **3-pool net flow** trades leak↔siphon
  monotonically and that the trade concentrates on strand-observable / wide-posterior nodes (not the
  confidently-biased count nodes, which it correctly can't move). Build the FP%→`q` curve.

## 6. Open questions

1. **Count-side FP mechanism** — (a) widened-σ quantile vs (b) global log-odds shift vs (c) accept-floor.
   Lean (b)+(c): a uniform count-side shift for the knob, and document the unstranded-capture floor.
2. **Boundary-side variance** (§1) is needed first; it is a small addition (per-side Poisson floor).
3. **Does the quantile help the 3-pool picture** (gDNA↔nRNA↔mRNA) — e.g. FP-averse `q` pulling
   FL-ambiguous reads back toward gDNA — or is that better served by the FL/strand discriminator? Measure
   on the suite before committing the knob's scope.
