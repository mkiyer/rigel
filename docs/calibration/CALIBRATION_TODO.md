# CALIBRATION_TODO — open issues toward a production-ready deconvolution

**Status:** living tracker. Opened 2026-06-09; major redirect 2026-06-10.

> **2026-06-10 architecture change — DECOUPLED strand/count.** The Phase-0/1 benchmark (below) proved
> the *joint* count×strand product mis-weights two unequal estimators: the count overdispersion crush
> was a load-bearing crutch that made the deconvolution defer to strand under capture, and removing it
> tripled the flagship leak. Resolution: **decouple** — route each node to the strand module
> (Beta-Binomial posterior; unbiased) OR the count module (raw density ratio; biased under capture),
> never a product, on disjoint node sets. The joint approach is archived in
> `archive/joint_deconvolution.md`; the live design is `decoupled_calibration_design.md`. This
> supersedes Phase 0/1 (`phase0_phase1_implementation_plan.md`) and the count-overdispersion line; the
> count-module accuracy work (point-5 imputation) continues under `count_channel_capture_design.md`.
> Strand-cleaning is replaced by the strand module's exact Beta-Binomial posterior (the clip-free
> robust estimator the deferred doc described). Issues #1–#5 below predate this and are re-scoped by it.

> **2026-06-10 redirect (superseded by the above).** Diagnostics (`scripts/debug/diag_*.py`) showed the dominant capture leak is
> a **count-channel** failure, not strand-cleaning: a single global NB count overdispersion mistook
> hybrid-capture's on/off-target *mean* heterogeneity for dispersion (α_crossing 0.0005→86 off→on),
> crushing the count concentration to ~0 and annihilating the count clue exactly where it is the only
> signal (unstranded+capture). **Phase 0** tore down that overdispersion model (count concentration =
> the gDNA count behind the node, Poisson precision); **Phase 1** floored the count prior on the
> node's own mean instead of Jeffreys ½ (kills the no-gDNA-exon false positive). The new direction +
> staged plan live in `count_channel_capture_design.md`; the shipped change in
> `phase0_phase1_implementation_plan.md`. Issue #2 is **reverted/closed**; Issue #1 (strand-cleaning)
> is **deferred** — concept captured in `strand_clean_robust_deferred.md`. Issues #3–#5 below predate
> the redirect (kept for context; re-evaluate against the post-Phase-0/1 benchmark).

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

**Status:** design doc written → `docs/calibration/strand_cleaning_robustness_design.md` (options
analysis + recommendation). Recommended: **O2 precision-weighted shrinkage** — treat the clean as a
measurement with precision `τ=(2κ−1)²·N` (its Fisher info) and shrink toward a fallback g₀; the
precision-weighting algebraically cancels the 1/(½−κ) blow-up (`g=[4N(½−κ)(s−κ)+τ₀g₀]/[4N(½−κ)²+τ₀]`),
smooth and → g₀ at κ=½. Sub-decisions: target g₀ (node-type-aware: observable→1, exon→pre-deconv
global) + prior strength τ₀. Highest-value foundational fix.

**TRIAL implemented 2026-06-09 (uncommitted, may revert):** O2 with the *minimal* g₀=1 (fixed
"keep raw count" fallback, not yet node-type-aware) + config `strand_clean_prior_weight` τ₀ (default
1.0). `strand_clean_gdna_frac` does `g=[4N(½−κ)(s−κ)+τ₀]/[4N(½−κ)²+τ₀]`, threaded through the region
clean (`node_gdna_density`) and boundary-side clean (`deconv_sides`); τ₀=0 reproduces legacy exactly.
Smoke: gdna400 cap_on **ss0.99 net leak 4.58% UNCHANGED** (high SS: τ≫τ₀ ⇒ shrinkage inert ⇒ win
preserved, no regression). O2's effect is at **low SS**.

**FULL 20-SCENARIO RESULT (2026-06-09): g₀=1 is behaviorally INERT** — it smooths the cliff but
reproduces its numbers (flagship 4.58→4.61%; gdna_none ss0.50 FP −4.48→−4.5%). Reason: at κ≈½ a node
defaults to g₀, and g₀=1 = "all gDNA" = exactly the old cliff. The residual concentrates in the
**unstranded (ss0.50)** conditions: gdna_none ss0.50 −4.5% (RNA→gDNA FP), unstranded+capture 13.6%
(gdna100) → 17.4% (gdna400) → 18.7% (gdna1000). **To make O2 useful, g₀ must be the LIBRARY GLOBAL
gDNA fraction at exon/non-observable nodes** (so an unstranded node defaults to the actual rate: ≈0
for gdna_none → kills the FP; ≈0.9 for gdna1000 → correct) — the node-type-aware g₀ (observable→1,
exon→pre-deconv global). Next: implement g₀=global, OR revert O2 (g₀=1 adds config+code for no
behavioral change). Mechanism (smooth shrinkage, no high-SS regression) is validated; the TARGET is
the lever.

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

## Issue #3 — flagship gdna400 capture-on ss0.99 leak (6.03% → 4.58%, gene-edge bug FIXED)

**Root-caused at fragment resolution (2026-06-09)** via `scripts/debug/introspect_locus81.py` on the
silent locus 81 (all gDNA): the dominant leak was **gene-edge boundary sides** (a stranded gene exon
abutting intergenic) defaulting to the **Jeffreys ½ prior** — both clues gated off at once: the
strand likelihood was *skipped* (`joint_deconv` required BOTH flanks to share a strand, so a
NEG↔NONE edge was wrongly strand-blind) and the count clue was *nuked* (α_crossing=127.8 ⇒ N_eff≈0).
194 such sides (~2/gene) carried ~164k gross rna_mass = ~68% of the in-locus leak.

**FIXED (this session):** intergenic (TS_NONE) is a strand **wildcard** — `_compute_side` now orients
a one-sided boundary by the stranded flank (`{POS,NONE}→POS`, `{NEG,NONE}→NEG`; only `{POS,NEG}` /
`TS_AMBIG` stay undefined). Gene-edge sides recover sfrac=½ → gDNA: locus-81 deconv 0.979→0.995,
dataset strand-blind edge leak 164k→1.5k, **flagship net leak 6.03%→4.58%**, no FP regression
(gdna_none ss0.50 4.48%→4.51%, ss0.99 0.04%). Golden regenerated; full suite green.

**Remaining 4.58%** is now the *within-gene* boundary sides' median bias (= **Issue #1**, the ~281k
strand-observable residual the tally showed unchanged: each within-gene gDNA side reads ~0.4% RNA
because the bounded posterior median truncates at gf=1) **plus the EM** (1.89%; uniform-density gDNA
component under-competes for sense-exonic gDNA under capture — `em_solver` scalar `gdna_eff_len`).
Note the small-N caveat: at low coverage the now-strand-observable edges can mis-clean by strand
noise (toy golden shifted gDNA slightly *down*) — another face of Issue #1 (no uncertainty model).

**To do (remaining):** Issue #1 (uncertainty-respecting strand clean — fixes the within-gene median
bias AND the small-N edge noise); then the EM uniform-density gDNA competition / eff_len.

**Status:** gene-edge bug FIXED (−1.45pp net, clean). Residual handed to #1 + the EM.

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
