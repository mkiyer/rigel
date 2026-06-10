# Count-channel redesign — phased implementation plan

**Status:** authoritative roadmap for the calibration work. 2026-06-10. Supersedes the pre-decoupling
version of this doc. Builds on the shipped decoupled architecture
(`decoupled_calibration_design.md`, main@8d8b896). Implementation-grade: each phase is a self-contained,
benchmark-gated PR. Open design decisions are called out inline — settle them before/within each phase.

> ## Phasing reconciliation (2026-06-10) — READ FIRST
>
> The original 1–5 plan below has been reorganized as the work matured. Current state and mapping:
>
> | phase | what it is | status | detail doc |
> |---|---|---|---|
> | **Phase 1** | precision-weighted blend `g=w·g_strand+(1−w)·g_count`, `w=(2κ−1)²` | **implemented, UNCOMMITTED** (flagship +1pp, expected to clear once 4-mean debiases count) | this doc §1 |
> | **Phase 2** | nascent 8-scenario benchmark + head-to-head sim | **committed** (branch `phase2-nascent-benchmark`); suite **not yet generated/run** | this doc §2 |
> | **Phase 3** *(clean→count)* | remove nascent from the count anchors via strand-cleaning | **ABSORBED into Phase 4-mean** — the splice-junction spliced channel + the 3-term/strand-resolved sweep is the cleaner mechanism (no separate anchor-cleaning pass) | — |
> | **Phase 4-mean** | splice-junction **gDNA-fraction** — debias the count *mean* under capture | **4-mean.1 (eligible 2-term + fallback) implemented**; 4-mean.2 sweep + 4-mean.3 (3-term / per-strand) pending | `count_mean_bias_design.md` |
> | **Phase 4-var** | count posterior **variance** (variance~mean, validated `var∝mean²`) | designed, pending | `count_posterior_design.md` (+ parked `count_local_dispersion_design.md`) |
> | **Phase 5** | FP-rate **quantile** knob over the combined posterior | designed (this doc §5), pending | this doc §5 |
>
> **Validation gate:** Phase 4-mean.3 and Phase 4-var both require the **nascent benchmark** (Phase 2)
> to be generated and run. The remaining sequence is in §"Sequencing & dependencies" (revised below).

## 0. Context & the unifying idea

The decoupling shipped a **handoff**: each node routes to the strand module (Beta-Binomial posterior)
or the count module (raw density ratio), disjoint. That recovered the flagship and fixed the
pure-RNA false positive, and isolated the residual **unstranded+capture leak (~18–21%)** as a pure
count-module problem. This redesign refines the count side and the strand/count *combination* into a
single coherent object:

> **Per node, one combined posterior over the gDNA fraction `g`** = the strand posterior and the count
> posterior, blended by a **precision weight** `w` that → 1 at high strand specificity and → 0 when
> unstranded. The reported gDNA fraction is a **user-chosen quantile** of that posterior (the
> false-positive-rate knob). The count posterior's *mean* is itself improved by feeding the count
> module **strand-cleaned** counts where strand is available.

Five phases assemble this. The theory (derived in the design discussion) in brief:

- **Weight** = the strand channel's *discriminability* `w = (2κ−1)²` — the per-fragment strand Fisher
  information / squared standardized separation between gDNA (sense ½) and RNA (sense κ). It → 1 at
  high SS, → 0 at unstranded, recovers the old `ss²` heuristic *with* grounding, and needs no count
  variance. Crucially it gates on **effect size, not significance**, so a near-unstranded library
  (κ≈½) gets `w≈0` *regardless of depth* — dissolving both the hard gate and the high-depth
  false-positive risk a significance test could never handle. We deliberately do **not** use naive
  inverse-variance weighting: it assumes both estimators are unbiased (count is biased under capture →
  MSE≠variance) and it penalises strand for overdispersion variance that averages out at locus
  aggregation (the benchmark proves strand-heavy SS=99% is correct).
- **FP knob** = a quantile `q` of the combined posterior. Under an asymmetric gDNA/RNA loss the
  Bayes-optimal point estimate *is* a quantile; `q=0.5` (median) is neutral, `q→1` is FP-averse. One
  knob, both channels, uncertainty-aware (it bites where the posterior is wide).
- **Count posterior**: to take a quantile on count-routed nodes (and to optionally sharpen the
  weight's weak-SS transition), the count module must emit a *variance*, not just a point.
- **Clean→count**: the count module's density anchors should be the strand module's *cleaned* gDNA
  counts (where strand is available), not raw counts — so nascent RNA at introns/exon-intron seams no
  longer inflates the imputation of AMBIG/exon nodes. Unstranded is unchanged (nothing to clean).

## 1. Diagnostic evidence (the problem the count side must solve)

Post-decoupling benchmark, the residual lives entirely in **unstranded + capture** (~18–21% gDNA→RNA
leak), and it is a count-module accuracy problem: under hybrid capture the on-target exon gDNA density
is imputed from **depleted off-target** observable neighbours, so the count estimate under-calls exon
gDNA by ~2× (oracle audit: imputed exon gDNA fraction ≈ 0.41 vs truth ≈ 0.91, pre-transport;
`scripts/debug/diag_imputation_truth.py`). When strand is present it rescues this (flagship 3.7%);
when unstranded, the count module is the sole estimator and the bias is exposed. The redesign attacks
it from two sides: **clean→count** (remove nascent contamination from the anchors) and the
**count posterior** (honest uncertainty so weak imputations defer appropriately).

---

## Phase 1 — precision-weighted strand→count deference (dissolve the gate)  *[was #1]*

**Goal.** Replace the hard handoff + identifiability gate with a smooth blend
`g = w·g_strand + (1−w)·g_count`, `w = (2κ−1)²` for strand-observable nodes (`w=0` otherwise). This is
the cleanest, most independent step and validates on the *existing* benchmark.

**Changes.**
- `strand_deconv._deconv_per_node`: for a strand-observable node, compute the strand posterior median
  `g_strand` (as now) **and** read `count_gdna_frac` as `g_count`; return `w·g_strand + (1−w)·g_count`.
  Non-strand-observable nodes: `w=0` → `g_count` (unchanged). `w = (2κ−1)²` from `rna_sense_frac`.
- `calibrate.py`: drop the `strand_identifiable` gate (the weight subsumes it — at κ=½, `w=0`). Keep
  the `CalibrationStrandError` on *zero spliced reads* (κ undefined). Pass `w` (or κ) into the deconv.
- `config.py`/`cli.py`: **remove** `strand_identifiability_confidence` (the param added with the
  decoupling) — the effect-size weight makes a significance gate unnecessary. Keep
  `strand_summary.strand_contrast_identifiable` only for the pipeline's *QC warning* (not a gate).
- `gdna_strand_llr_bias`: leave as-is for now (retired/replaced by the quantile in Phase 5).

**Validation (existing 20-scenario benchmark).** Stranded+capture unchanged (~3.7%, w≈0.96 ≈ the old
handoff); unstranded unchanged (w≈0 → count); gdna_none FP stays ~0 (κ≈½ → w≈0); behaviour smooth
across SS (no cliff). Add a unit test: `w=(2κ−1)²`, `g→g_count` as κ→½, `g→g_strand` as κ→1.

**Open decisions.** (a) library-level `w` (one κ) vs per-node depth-aware `w` — lead with library-level
(elegant); (b) the *shape* (`(2κ−1)²` vs a depth-aware form) in the intermediate-SS band (~90%, not in
the current benchmark) — revisit if a weak-SS scenario shows it too blunt.

---

## Phase 2 — nascent-RNA benchmark (the test harness for Phase 3)  *[was #5; sequenced early]*

**Goal.** An 8-scenario suite to exercise the regimes where clean→count matters. Built before Phase 3
because nascent is the *only* thing that distinguishes cleaned from raw imputation.

**Spec.** `gDNA 2× (2:1 gDNA:RNA) × {capture on/off} × {ss 0.50/0.99} × {nascent on/off}` = 8 scenarios,
on the existing synthetic genome/transcriptome.
- **Nascent model:** 50% of transcripts carry nascent RNA, 50% do not; each nascent-bearing
  transcript's `nrna_abundance` = a **random per-transcript ratio** of its mature abundance,
  **additive** (mature fragment counts fixed; nascent layered on top). The sim already supports
  per-transcript `nrna_abundance` (`rigel.sim.annotation`).
- **Fixed RNG + fixed nascent abundances baked into the suite config**, so a nascent-**on** scenario
  generates the *identical mature fragments* as its nascent-**off** twin plus the nascent layer →
  strictly head-to-head.

**Changes.** A suite-config generator (`scripts/sim/configs/`) + the seeded per-transcript nascent
assignment; wire into `scripts/sim/simulate_suite.py` / `evaluate_suite.py`. Confirm the net-flow
analysis counts nascent fragments as RNA (so nascent→gDNA is scored as leak) — `benchmark.py` already
tracks `n_nrna_expected`.

**Validation.** Harness only. Sanity: nascent-off scenarios reproduce the corresponding cells of the
existing suite; nascent-on adds the expected extra RNA mass.

**Open decisions.** Random-ratio range (e.g. nascent/mature ∈ [0.05, 1.0]?) and the RNG seed value —
bake both into the config for reproducibility.

---

## Phase 3 — clean→count imputation  *[ABSORBED into Phase 4-mean]*

> **Superseded 2026-06-10.** This phase's goal — stop nascent RNA from inflating AMBIG/exon gDNA
> imputation — is now met by **Phase 4-mean** more directly. The splice-junction **spliced channel** is
> an unconflated mature reference (gDNA/nascent never splice), and the **3-term fraction +
> strand-resolved sweep** (`count_mean_bias_design.md` §5, §5.1) removes nascent from the gDNA numerator
> using the strand split — the same nascent-removal this phase aimed at, without a separate
> "strand-clean the density anchors" pass. The acyclicity concern below is moot (no density→clean cycle;
> the spliced reference and the carried per-strand split are inputs, not a refit of the density). Kept
> for the record; do not implement separately. The original text follows.

**Goal.** Build the count module's density field from **strand-cleaned** gDNA counts where strand is
available, so AMBIG/exon imputation is no longer inflated by nascent. Unstranded → unchanged.

**Changes.**
- Run the strand deconvolution on **observable** nodes *first* to produce per-node cleaned gDNA
  **mass**; `density_model` reads cleaned gDNA counts (not raw) for strand-observable observable nodes,
  raw for the rest; imputes AMBIG/exon nodes from the cleaned anchors as today.
- **Clip deferral** (the prior lesson): the cleaned masses feed imputation **unclipped/aggregated**;
  clipping is deferred to the locus prior. Use the strand posterior **mean** for the density field
  (aggregates linearly), not the median.

**Acyclicity risk (must break a new cycle).** Today: κ → density(raw) → strand overdispersion (seeds
weighted by `density.count_gdna_frac`) → deconv. If density now needs *cleaned* counts, and cleaning
needs the overdispersion, that's circular. **Break it** by fitting the gDNA strand overdispersion from
**signature-pure raw seeds** (intergenic = pure gDNA, weight 1; observable boundary sides = raw
crossing ratio) *before* the density build — intergenic needs no cleaning. New order: κ → strand
overdispersion (raw signature-pure seeds) → strand-deconv observable nodes (cleaned) → count density
(from cleaned) → blend. Keep it acyclic and assert no density→overdispersion→clean→density loop.

**Validation (Phase 2 benchmark).** Stranded+capture+**nascent** improves (AMBIG/exon imputation no
longer nascent-inflated); unstranded+nascent unchanged; all nascent-**off** cells bit-stable vs Phase 1.

**Open decisions.** Mean vs median for the cleaned density field (lean mean); exactly which observable
nodes contribute cleaned vs raw at low `w` (a node at κ where w is small is barely cleaned — fine).

---

## Phase 4 — the count-module posterior  *[split into 4-mean + 4-var]*

> **Reorganized 2026-06-10.** The original "count posterior" is now **two independent efforts** with
> separate detail docs, because the **mean** (bias) and the **variance** turned out to be different
> problems (variance ≠ bias):
>
> - **Phase 4-mean — `count_mean_bias_design.md`** (the splice-junction gDNA-fraction). Debias the count
>   *mean* under capture using the spliced channel as a clean mature reference. **4-mean.1** (eligible
>   exon regions, 2-term boundary fraction + absolute-density fallback) is **implemented**
>   (`calibration/splice_junction.py`, wired into `calibrate.py`; Tests A+B). Remaining: **4-mean.2** the
>   run-interior sweep (carry-over), **4-mean.3** the 3-term fraction + strand-resolved per-strand sweep
>   state (§5, §5.1 of that doc) — which also absorbs the retired Phase 3 (clean→count).
> - **Phase 4-var — `count_posterior_design.md`** (the count *variance*). The variance~mean fit from
>   paired 2-anchor boundary disagreements (validated `var∝mean²`, conflation-free). Emits a count
>   posterior (mean + variance / Beta) so Phase 5's quantile knob can act on count-routed nodes. Parked
>   extension: `count_local_dispersion_design.md`.
>
> The original sketch follows for reference.

**Goal.** The count module emits a **distribution** over `g` (mean + variance, or a Beta), not a point.
This enables the Phase-5 quantile knob on count-routed nodes and can sharpen the Phase-1 weight's
weak-SS transition. Flagged by the team as needing careful design and likely iteration/critique.

**Design sketch (to critique).** `g_count` mean = `clip(ρ_local·eff/M)` (as now). Variance from two
sources: (a) **Poisson** counting noise on the gDNA count behind the estimate; (b) **imputation
uncertainty** for non-observable nodes — the anchoring crossing count's precision (anchor strength)
and between-anchor disagreement. Parameterise as a Beta with mean `g_count` and concentration = an
effective count that shrinks with imputation uncertainty (a no-anchor node → near-flat).

**Open decisions (the crux).** (a) exact variance decomposition and how imputation uncertainty
propagates from anchor counts → `ρ` → `g`; (b) Beta vs Gaussian parameterisation; (c) whether to then
replace the Phase-1 `(2κ−1)²` weight with a **bias-aware** blend (count's known capture bias inflates
its effective variance — the MSE-optimal weight) now that a count variance exists, or keep `(2κ−1)²`
for its robustness. Expect a round of critique here.

**Validation.** Both benchmarks: no regression vs Phase 3; the unstranded+capture leak should *narrow*
as honest count uncertainty lets weak imputations defer (and the Phase-5 knob then trade FP/sensitivity).

---

## Phase 5 — the false-positive-rate quantile knob  *[was #3]*

**Goal.** One user knob controlling the gDNA↔RNA confidence (FP rate) for **both** channels: a quantile
`q` of the per-node combined posterior. `q=0.5` (median) = neutral default.

**Changes.**
- Combined posterior per node: strand BB (where observable) ⊗ count posterior (Phase 4), weighted by
  `w`. Cheapest faithful form: Gaussian — mean `= w·μ_s+(1−w)·μ_c`, var `= w²σ_s²+(1−w)²σ_c²` — and the
  q-quantile `= mean + Φ⁻¹(q)·sd` (this is also the natural **z-score** parameterisation). Read `g` at
  `q` instead of the bare median.
- `config.py`/`cli.py`: `gdna_deconv_quantile` (default 0.5). Decision-theoretic semantics: `q` ↔ the
  tolerated asymmetric gDNA/RNA loss ratio. A thin calibration layer can map a user-facing **FP-rate
  target** → `q` (derive from the loss ratio, or calibrate empirically on the benchmark).
- Retire `gdna_strand_llr_bias` in favour of `q` (the quantile is the uncertainty-aware generalisation
  of the log-odds shift), or keep it as a point-estimate fallback — decide here.

**Validation.** Sweep `q` on both benchmarks; confirm `q>0.5` monotonically trades gDNA→RNA leak for
RNA→gDNA siphon (and vice-versa), uncertainty-aware (bites on wide-posterior nodes), `q=0.5` reproduces
Phase-4 numbers. Add the FP-rate→`q` calibration check.

**Open decisions.** Gaussian-approx vs exact mixture for the combined posterior; quantile of the
combined posterior vs of each channel then blended; the FP%→`q` calibration (analytic vs empirical).

---

## Sequencing & dependencies  *(revised 2026-06-10)*

```
Phase 1   (precision weight, gate dissolved)     ── implemented, UNCOMMITTED (flagship +1pp pending 4-mean)
Phase 2   (nascent 8-scenario benchmark + sim)   ── committed (branch); SUITE NOT YET GENERATED/RUN
Phase 4-mean.1 (eligible 2-term splice fraction) ── DONE: splice_junction.py + calibrate wiring + Tests A/B
   │                                                validate on existing 20-cond benchmark (in progress)
Phase 4-mean.2/.3 (sweep + 3-term/per-strand)    ── needs the nascent benchmark; absorbs old Phase 3
Phase 4-var   (count variance ~ mean)            ── independent of 4-mean; needs nascent benchmark to validate
Phase 1 commit (blend) ───────────────────────── once 4-mean debiases count, re-check flagship, then commit
Phase 5   (quantile FP knob)                     ── needs Phase 4-var (count posterior) + Phase 1 (combined posterior)
```

**Critical path to a release-ready count side:** (1) confirm 4-mean.1 on the existing benchmark →
(2) **generate + run the Phase-2 nascent suite** (the validation harness for everything nascent-aware) →
(3) 4-mean.2/.3 (sweep + 3-term/per-strand) and 4-var, both validated on the nascent suite →
(4) commit Phase 1 once the flagship +1pp clears → (5) Phase 5 quantile knob. Each step a
benchmark-gated PR. The old 1-5-4-2-3 ordering is superseded by this.

## Magic-number audit

- **Weight `(2κ−1)²`**: not a tunable — the per-fragment strand Fisher information / squared effect
  size; κ is fit from the data.
- **Gate**: *removed* (dissolved by the weight); `strand_identifiability_confidence` deleted.
- **FP knob `q`**: a user preference (default 0.5 = neutral median), decision-theoretic, not a fitted
  constant.
- **Count posterior**: variances are derived from counts/anchor strength — no fitted shrinkage knobs
  intended (guard against re-introducing the old overdispersion-style constants).
- **Nascent benchmark**: the random-ratio range + RNG seed are *test fixtures*, not model parameters.
