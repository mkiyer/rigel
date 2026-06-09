# CALIBRATION_TODO — open issues toward a production-ready deconvolution

**Status:** living tracker. Opened 2026-06-09 after the count-overdispersion work (Phase A/B,
main@a23e236) + a calibration code review (main@978ad9e) + a deep dive into the residual
gDNA→RNA leak on the flagship capture scenario. Several fronts are open at once; this document
enumerates them so we stop conflating them and can attack each on its own terms.

## Where we are

- **Phase A/B (shipped, main@a23e236):** the count-prior concentration is now the
  overdispersion-limited effective count `N_eff = N/(1+α·N)` (count-side twin of the strand
  Beta-Binomial). Net gDNA→RNA leak ~halved across the gDNA-present matrix (e.g. gdna1000 cap_on
  ss0.99 net 12.0%→5.4%; cap_off ss0.50 10.2%→4.9%).
- **Cleanup (shipped, main@978ad9e):** dead code, stale docstrings, naming, magic-number naming.
- **Flagship dive (this session):** gdna400 (4:1) capture-on ss0.99 net leak = **6.03%**, decomposed
  into **4.15% calibration prior** under-call + **1.89% EM** (per-locus `gdna_expected` vs
  `gdna_prior_count` vs `gdna_observed`). Count clue correctly ~OFF under capture (α_crossing≈128
  ⇒ N_eff≈0.008; α_contained≈0.95 ⇒ N_eff≈1).
- **Strand-clean de-bias prototype (REVERTED, not shipped):** confirmed a real per-node
  posterior-median **boundary bias** (pure-gDNA node, true gf=1 → reported 0.91 at N=20, 0.97 at
  N=140) that dominates the *no-capture* leak (cap_off ss0.99 → calibration 0.7956 vs 0.80, near
  perfect) — but it does **not** fix the *capture* target (6.03→6.46%, interior-gf nodes aren't
  boundary-biased) and it **explodes unstranded** (gdna_none cap_on ss0.50 → 42.9% FP). The
  prototype was reverted; its lessons feed issues #1 and #3 below.

Memory: [[calibration_capture_leak_strand_boundary_bias]], [[calibration_benchmark_findings]],
[[gdna_leak_root_cause_efflen]], [[em_strand_per_fragment_penalty]],
[[calibration_no_magic_numbers]]. Benchmark suite was wiped; only gdna400 cap_on/off ss0.99 +
gdna_none ss0.99/0.50 cap_on are currently re-quanted. Rerun via skill `calibration-benchmark`.

---

## Issue #1 — strand-cleaning is unbiased but not robust (explodes as SS→½)

**Problem.** The closed-form strand clean `gdna_frac = (s − κ)/(½ − κ)` (and its count form
`N_gdna = N·(s−κ)/(½−κ)`) is the right **unbiased** estimator — and the de-bias prototype confirmed
it should be allowed to **oscillate, including transiently negative**, so the *aggregate* is
unbiased (clipping each node first re-introduces the boundary bias). **But** its variance is
∝ `1/(2κ−1)²`, which **diverges as κ→½**: at SS=0.50 there is *zero* strand information (infinite
uncertainty), yet the closed form returns a finite (wildly amplified) number. The estimator is
**unbiased but does not represent its own uncertainty**, so naive aggregation at low SS produces
garbage (gdna_none cap_on ss0.50 → 42.9% false gDNA).

**This is not a bug in the math — it is a missing uncertainty model.** SS=½ ⇒ maximum uncertainty
⇒ the estimate should carry ~zero precision and contribute ~nothing to any aggregate, not explode.

**Clean-slate question (to design):** how should strand-cleaning work so that it is (a) unbiased,
(b) allows transient negatives, AND (c) respects uncertainty — i.e. each cleaned count comes with a
variance ∝ `1/(2κ−1)²` (and ∝ 1/N), and aggregation is **precision-weighted** so low-SS / low-N
nodes fade smoothly to zero weight rather than dominating by noise. Likely shape: a per-node
posterior over the gDNA count (mean = linear clean, variance from strand noise + the `(2κ−1)²`
contrast), combined across nodes by inverse-variance weighting; the global κ contrast `(2κ−1)²` is
the natural precision scale (it already appears as the strand Fisher information). This subsumes the
ad-hoc `gdna_strand_llr_bias` knob and the vetoed global strand-reliability metric.

**Status:** design open. Highest-value foundational fix; #3 and #5 partly depend on it.

---

## Issue #2 — count overdispersion: new, parameterized by priors/magic numbers, under-validated

**Problem.** `count_dispersion.py` (Phase B) introduced two constants — `count_overdispersion_prior`
(α₀=1, geometric) and `count_overdispersion_prior_weight` (30, mirroring the strand fit) — plus the
per-type NB-MoM + pooled-trend shrinkage. It looks like a ~5pp net win, but it has **not** been
validated in isolation: we do not know whether α_crossing≈128 under capture (⇒ count essentially
off) is *correct distrust* or an over-estimate, whether the pooled-trend shrinkage helps or hurts,
or how sensitive the result is to α₀/weight.

**To do:** (a) validate the fitted α's against known-overdispersion synthetic seeds at realistic
(low) seed counts; (b) sensitivity-sweep α₀ and prior_weight; (c) decide whether the prior should be
parameterized differently (e.g. precision in seed-fragment rather than seed-node units); (d) confirm
the contained-vs-crossing split is the right granularity; (e) re-justify every constant against the
no-magic-numbers rule. Until validated, treat the 5pp as provisional.

**Status:** implementation shipped, validation owed.

---

## Issue #3 — flagship gdna400 capture-on ss0.99 still leaks 6% gDNA→RNA

**Problem.** Net 6.03% = 4.15% prior + 1.89% EM. The prior under-call splits: **silent loci −8.8pp**
(pure-gDNA exons, explained by the #1 boundary bias) vs **expressed loci −5.6pp** (the bulk of the
mass; *interior*-gf exonic nodes, NOT boundary-biased — the de-bias did not touch them). With the
count correctly off under capture, the deconv there is strand-only; the strand clean of a mixed node
is unbiased in the mean, so the expressed-locus deficit needs another explanation.

**Leading candidate (from [[gdna_leak_root_cause_efflen]]):** the per-region gDNA **effective
length / IPR** (`gdna_eff_len`, `assemble_priors`) inflates ~2.2× under capture, so the EM gives the
gDNA component too-low per-base density and it loses the exonic competition to RNA — a lever the
de-bias left on physical mass. Also re-examine the **EM 1.89%** (per-fragment strand scoring /
hard-vs-soft assignment; note the de-bias prototype showed the EM partly *self-corrects* a worse
prior, so the per-fragment likelihood is doing real work).

**To do:** quantify the eff_len/IPR inflation at this scenario (join `loci.feather`
`gdna_eff_len_em` to truth); test exon-granular gDNA support; separately profile the EM leak.

**Status:** open; depends partly on #1 (the silent-locus piece) and on the eff_len investigation.

---

## Issue #4 — gDNA-fraction imputation opportunity

**Problem / opportunity.** The deferred "path-3 / DNA-fraction" lever: directly estimate the gDNA
fraction (e.g. from a DNA-fraction model / external information) for regions where neither count
(capture-heterogeneous) nor strand (unstranded) is sufficient. Could anchor the calibration where
both clues are weak (the #5 worst case) and improve the capture prior (#3).

**Status:** unscoped; revisit after #1–#3 to see how much residual it would actually close.

---

## Issue #5 — unstranded + capture-on is egregious (worst case, but diagnostic)

**Problem.** gdna_none cap_on ss0.50 already mislabels ~4.5% of a pure-RNA library as gDNA under
Phase B (and 42.9% under the naive de-bias). Strand is uninformative (κ≈½) AND count is
capture-biased — the two clues that separate gDNA from RNA are both blind. Accepted as the
worst case, but **studying it isolates the count-information + capture-handling issues** with the
strand confound removed: whatever still mislabels here is a pure count/capture artifact.

**To do:** use this condition as the count/capture test bed once #1 (uncertainty-respecting clean)
lands — at κ=½ the strand clean should contribute *zero* precision, so any residual gDNA call must
come from the count clue and exposes its capture bias directly.

**Status:** worst-case, deprioritized for shipping but kept as a diagnostic lens for #1/#2.

---

## Suggested order of attack

1. **#1** (uncertainty-respecting strand-cleaning) — foundational; fixes the unstranded explosion
   and the no-capture boundary bias in one principled estimator, and is a prerequisite for cleanly
   using #5 as a diagnostic.
2. **#3 eff_len/IPR** — the leading explanation for the flagship capture leak (independent of #1).
3. **#2** count-overdispersion validation/parameterization (can proceed in parallel).
4. **#5** as the diagnostic lens, then **#4** if residual remains.
