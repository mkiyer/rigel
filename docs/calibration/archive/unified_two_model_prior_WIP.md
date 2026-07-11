# Unified two-model gDNA prior — DESIGN IN PROGRESS (paused 2026-07-05)

**Status:** WIP / paused. NOT shipped. The production ship is **M2** (deterministic `_LogLinearVarMean` σ²_g
line replacing the bistable P-spline) — see `M2_loglinear_sigma2g_design.md`. This document records the
unified two-model prior we designed, implemented, and benchmarked, why it does not yet net-win, and exactly
where to resume. The full experimental implementation is saved as a patch:
`scratchpad/session_experiment.patch` (the working tree was reverted to the clean M2 ship after measurement).

## 1. The motivating insight (user, validated)
The Phase-2 KDE has a bandwidth floored at each node's Poisson log-density noise (~`1/E`), and a zero-count
node's density is itself floored at `1/E`. So the DEPLETED background — dominated by zeros and short/low-E
nodes — collapses to a spurious "mode" at ~`1/E`, ABOVE the true near-zero background. **The KDE is
structurally blind below its noise floor: a good instrument for the ENRICHED mode, a bad one for the DEPLETED
mode.** This is why a KDE-only exon prior (`dropall`) worked when gDNA was present but wrecked no-gDNA.

The correct depleted-floor estimator is POOLING (many zeros over huge total exposure → a confident near-zero
rate), which `_floor_estimate` already computes: `ρ_floor = (1+ΣG)/ΣE` over intergenic+intron regions, zeros
kept. Critically, `ρ_floor` is ADAPTIVE — it tracks the actual background in every regime (≈0 in no-gDNA,
elevated-uniform in gdna300/capoff, depleted-off-target in gdna300/capon).

## 2. The architecture (as implemented in the patch)
Two models on `log ρ_g`, combined as a **MIXTURE** (not a sum — summing a background Gaussian onto the KDE
drags enriched nodes, the M2 regression):
- **Background model B** = `_floor_estimate` (`ρ_floor`, spread `√(s2_floor+var_mean_floor)`). Added to the
  KDE as ONE separate analytic component (its own wide bandwidth, so it does NOT corrupt the teacher
  bandwidth). Realized via a per-SAMPLE bandwidth in `_weighted_kde_logpdf` (teachers at the fitted `h`, B at
  σ_bg); evaluating the sum via `logpdf_kernel` IS the normalized background⊕enriched mixture.
- **Enriched model K** = the KDE, trained on SOFT-weighted teachers `w = gcount/(gcount+1)` (gcount = f_g·M),
  so it focuses on measurable-gDNA nodes and never represents the sub-noise-floor floor. Its complement
  `Σ(1−w) = Σ1/(gcount+1)` is the pooled background-evidence mass → the B component's weight (P2 = P3 are one
  choice; both deterministic, continuous, magic-free).
- **Pass 1** = floor override on intergenic/intron only, exons prior-free (`drop_nonfloor`); **Pass 2** =
  drop the genome Gaussian on exons/boundaries, the seeded mixture is the sole exon prior. Pass-2 always runs
  (don't skip when few teachers). Boundaries excluded from the KDE term (crossing-scale ≠ contained-scale).

## 3. What the benchmark showed (16-condition suite, quick_3to1_5mb)
- **The background SEED fixes the no-gDNA FP** the KDE-alone couldn't: no-gDNA-unstranded-capon leak
  35226 (`dropall`) → **10015** (seeded), near baseline's 7214. The core insight WORKS.
- **Better on gDNA-present** (enriched mode shields): g300-unstranded-capon spearman 0.703 → 0.741.
- **But net mean spearman 0.803 → 0.754** — broadly WORSE on no-gDNA (6 conditions regress >0.02).

## 4. Why it does not net-win yet (the open problem)
The genome baseline (σ²_g) does **two** jobs; the seed only replaces one:
1. **FP suppression** — the accurate background seed handles this. ✓ (validated)
2. **Consistent anchoring of every unstranded exon for good abundance RANKS** — this needs precision that is
   STRONG at no-gDNA (consistent → good ranks) and WEAK under capture (don't drag enriched). That
   strong@no-gDNA / weak@capture profile IS the density-dependent σ²_g. A fixed / floor / seed anchor can't do
   both: floor-on-all-in-pass-1 drags enriched teachers (trades the no-gDNA regression for a g300 one). So
   dropping the density-dependent genome baseline loses job (2) — and **M2 keeps exactly that, deterministically.**

The **seed-weight tension** is the same conflation relocated: `w_seed` is pulled UP by no-gDNA FP-suppression
and DOWN by not-drowning-the-enriched-mode under capture (the review's P2 risk). Count-based `w_seed` doesn't
escape it under capture (background is the genomic majority).

## 5. Dead ends found (do not repeat)
- **M1 measure Jacobian** (`−log f_g` to make the KDE a proper density in the uniform-f_g measure): THEORETICALLY
  correct, EMPIRICALLY CATASTROPHIC — g300 leak 59k→394k (over-rewards low f_g → massive gDNA under-call). The
  KDE+Jeffreys balance is tuned with the gDNA term BARE; keep it bare.
- **Precision-convolved KDE** (`h_node=√(h_pop²+σ_node²)`, Direction B): consistently HURT abundance (smears
  the enriched-mode lift on exactly the uncertain enriched exons that need it).
- **M2 + seed hybrid** (keep the genome baseline, add the seed + soft teachers): best spearman (0.802 ≈
  baseline) but gDNA leak +73% (582k) — the soft-teacher/seed KDE under-calls gDNA on enriched. Bad trade for a
  deconvolution tool.

## 6. The adversarial design review (6 lenses + synthesis) — key resolved decisions
`scratchpad/` transcripts + the synthesis: P1 = pass-1 floor-only; P2/P3 = the soft complementary weights
above; the KDE must NOT be capped (defeats AMBIG resolution); boundaries need a crossing-scale reference
(KDE_boundary_prior_review Direction A) not the contained-scale seed. Determinism blockers it surfaced (both
FIXED in the M2 ship): the KDE bandwidth `_weighted_median` searchsorted step (→ interpolated), and (deferred)
the log-odds posterior-median grid snap in `simplex_logodds` (D1 — interpolate the CDF-0.5 crossing; measured
to improve determinism 0.05%→0.007% but cost a little accuracy, so NOT shipped).

## 7. How to resume
1. Re-apply `scratchpad/session_experiment.patch` (or re-implement §2) behind a `RIGEL_PRIOR` toggle.
2. **The crux to crack:** give pass-1 unstranded exons a DETERMINISTIC density-dependent anchor that replaces
   job (2) without σ²_g — candidates: a 3-pass scheme (pass-2 KDE anchors pass-3), or a strand-power-attenuated
   background pull, or accept M2's σ²_g AND add the seed only where it strictly helps (no-gDNA), gated by a
   deterministic capture/no-gDNA signal derived from ρ_global vs ρ_floor.
3. Fix the boundary crossing-scale reference (Direction A) so M5 isn't just exclusion.
4. Ship the D1 median-interpolation determinism fix alongside (or instead do the C++ scanner ordered/fixed-point
   accumulation — the true bit-identity fix — which makes both the median snap and the spline moot).

## 8. The north star
`gdna_prior_clean_slate_architecture.md` — one mixture (adaptive background mode + learned enriched mode(s)),
gravity = the mixture's shape. The unified prior here is its first implementation; it needs the job-(2)
anchoring solved to land.
