# Calibration — remaining implementation phases (clean plan)

**Status:** authoritative forward plan, for review. 2026-06-10. This document **replaces the confusing
`4-mean.1 / 4-mean.2 / 4-mean.3 / 4-var / Phase-5` nomenclature** with three plainly-numbered phases in a
discrete implementation order. The older `count_channel_capture_design.md` / `count_mean_bias_design.md`
/ `count_posterior_design.md` remain as background design records; this doc is the one to implement from.

> **Name mapping (so we don't lose the thread):**
> | this doc | old name |
> |---|---|
> | Phase 0 (shipped) | decoupling + Phase 1 blend + `4-mean.1` splice 2-term |
> | **Phase 1** | `4-var` — count-module posterior (variance) |
> | **Phase 2** | `Phase 5` — false-positive-rate quantile knob |
> | **Phase 3** | `4-mean.2` + `4-mean.3` — strand-resolved 3-term nascent removal + carry-over sweep |

---

## Part I — Overview

### Where we are (Phase 0, shipped @ d0bf62b)

- **Decoupled strand/count** + **precision-weighted blend** `g = w·g_strand + (1−w)·g_count`,
  `w = (2κ−1)²`.
- **Splice-junction gDNA-fraction (2-term)**: exon regions with a direct eligible splice-junction
  boundary get `g_count` from the boundary fraction `(crossing-unspliced/E^gDNA)/(…+crossing-spliced/E^RNA)`
  instead of the absolute-density imputation. Validated: unstranded+capture leak 15.4%→12.0%
  (critical gdna1000 18.4%→12.5%), no false positives, no regression.

### The remaining problem, split into what is and isn't reducible

The residual leak lives in **capture-on**, and decomposes into two structurally different buckets
(measured on gdna1000 ss0.50 cap_on: 562 non-observable regions):

1. **Junction-bearing regions (237)** — have a spliced reference. The 2-term already debiases them; what
   remains is (a) **nascent contamination** (the 2-term counts nascent as gDNA) and (b) the FL-consistency
   approximation. *Reducible* → **Phase 3** (nascent) and the deferred FL-consistency work.
2. **Single-exon / no-junction regions (308)** — *no* spliced reference anywhere in their run. The
   absolute density is the only signal, and under capture it is measured at depleted boundaries, not the
   enriched interior. **No splice method can reach these.** This is the classic hard case (unstranded +
   capture: strand is blind *and* count is capture-biased). The bias here is **largely irreducible** — we
   cannot manufacture an enriched-interior gDNA measurement that doesn't exist. What we *can* do is
   **represent its uncertainty honestly and give the user a principled FP-rate control** to trade the
   residual gDNA→RNA leak against an RNA→gDNA siphon. → **Phase 1 + Phase 2**.

So the count side (Phases 1–2) is the high-leverage, user-requested work for the headline number; Phase 3
completes the splice story for real-world antisense + nascent.

### Proposed implementation order (and why)

```
Phase 1 — Count-module posterior (honest uncertainty: mean + variance)
   │  foundation; independent; validatable on the nascent+capture suite (variance-calibration check)
   ▼
Phase 2 — False-positive-rate quantile knob  (the ONE user knob)
   │  needs Phase 1 (a count variance) + the Phase-0 strand posterior → combined posterior → quantile
   ▼
Phase 3 — Strand-resolved nascent removal (3-term splice fraction + carry-over sweep)
      independent of 1–2; completes the splice story for antisense(AMBIG) + nascent; LOW measurable
      impact on the *current* suites (see §Phase-3 scope), so sequenced last
```

**Why this order, not splice-first:** the headline residual (unstranded+capture single-exon) is
unreachable by splice and is exactly what Phases 1–2 address; Phase 2 is your stated "one FP-rate knob"
requirement; and Phase 3's measurable benefit needs a richer (antisense+nascent) validation genome we
don't yet have, so it is correctness/completeness rather than a benchmark mover. **Hard dependency:** only
2-after-1. Phase 3 may be reordered earlier if real-world antisense/nascent validation becomes the
priority. Each phase is a self-contained, benchmark-gated commit.

### Magic-number discipline (applies to all phases)

No fitted shrinkage constants. The only scalars introduced are: (1) `α_impute`, **fit from data** (the
variance~mean slope), not chosen; (2) `gdna_deconv_quantile`, a **user preference** (default 0.5 =
neutral), decision-theoretic, not fitted. Anything else (kernels, decays, floors-with-knobs) is called out
as an open question to avoid before it sneaks in.

---

## Part II — Phase details

---

# Phase 1 — Count-module posterior (honest uncertainty)

### 1.1 Problem & goal

The count module emits a **point** `g_count = clip(ρ_local·eff/M)`. To (a) let the blend defer to strand
where the count imputation is *genuinely* uncertain and (b) enable the Phase-2 quantile on count-routed
(unstranded / single-exon / AMBIG) nodes, the count module must emit a **posterior**: mean + variance over
the gDNA fraction `g`.

> **Critical caveat (variance ≠ bias):** an honest variance does **not** fix the capture *bias* of the
> single-exon residual — a biased-low estimate can be low-variance (two depleted boundaries agree). Phase 1
> models **uncertainty**, not bias; it is the prerequisite for the Phase-2 knob, which is what actually
> lets the user move the single-exon leak. Do not over-claim Phase 1 as a leak fix.

### 1.2 Design — the variance

Per node, `μ_g = g_count` (unchanged — the 2-term splice fraction or the absolute-density fallback), and
```
σ_g²  =  (eff/M)² · [ σ²_count(N_g)  +  σ²_impute(ρ_local)·eff² ]
```
two independent sources:

**(a) Counting noise — Poisson floor (v1).** `σ²_count(N_g) = N_g` (the gDNA count behind the estimate).
Parameter-free. Under-states overdispersion, but the imputation term dominates for the nodes that matter,
so Poisson is the v1 floor. *Local NB dispersion is explicitly deferred* (`count_local_dispersion_design.md`):
it pools unpaired cross-node counts → conflates rate-variation with dispersion (the same flaw that sank the
global `α`). Do not build it in v1.

**(b) Imputation uncertainty — boundary disagreement (the high-value piece).** A non-observable node is
imputed from ≤2 anchoring observable boundary sides with densities `d_L`, `d_R`. The disagreement *is* the
uncertainty:
```
σ²_impute  =  ¼·(d_L − d_R)²                # variance of the 2-point mean (the disagreement)
            + ½·(σ²_{d_L} + σ²_{d_R})        # each anchor's own counting noise, propagated
            + σ²_floor                       # so two agreeing high-count anchors aren't exactly 0
```
- 2 agreeing anchors → low variance (precise — *possibly biased*, §1.1).
- 2 disagreeing anchors (capture enriches one) → high variance → blend defers to strand.
- 1 anchor (gene/reference edge) → no cross-check → single-anchor count noise + an inflated floor.
- 0 anchors (global fallback) → maximal variance (count uninformative).

**The variance~mean fit (denoises the noisy per-node 2-point estimate).** A single node's `¼(d_L−d_R)²` is
a 1-dof, very noisy variance. Pool the `(mean=(d_L+d_R)/2, var=¼(d_L−d_R)²)` pairs across *all* 2-anchor
nodes genome-wide and fit a smooth trend, then read each node's variance off the curve at its own mean.
**Validated** (`scripts/debug/diag_variance_mean.py`, gdna1000): under capture-on the pooled pairs follow
**`var ∝ mean²`** (log-log slope 2.00) — the NB law `Var = α·μ²`. So:
```
σ²_impute(μ) = α_impute · μ²,   α_impute = pooled slope from the 2-anchor disagreements
```
**`α_impute` is fit from data, not chosen** (no magic number). It is **conflation-free**: `d_L`, `d_R` are
two estimates of the *same* node's rate (paired), so `(d_L−d_R)²` is genuine variance — unlike the
cross-node count dispersion that conflated rate-variation. 1-anchor / 0-anchor nodes read the curve at their
mean with an inflated floor.

**Run interiors (`intron│exon│exon│exon│intron`).** A node at fractional position `t` between the run-edge
anchors `d_L`, `d_R` uses the **linear-interpolation variance** (parameter-free — it *is* the variance of a
linear estimator between two noisy endpoints):
```
σ²(t) = (1−t)²σ²_{d_L} + t²σ²_{d_R} + t(1−t)(d_L−d_R)²     # maximal mid-run, where least anchored
```
Prefer this over a distance-penalty form (which needs a decay knob).

### 1.3 Design — the posterior object (NB → Beta → Gaussian)

These are three layers, not competing choices:
- **NB** is the count noise model (`Var = μ + α·μ²`) — governs the count/density uncertainty above.
- **Beta** is the *fraction* posterior: propagate `σ_g²` to `Beta(g)` with mean `μ_g`,
  `concentration = μ_g(1−μ_g)/σ_g² − 1`. The natural law over a bounded fraction `g ∈ [0,1]`.
- **Gaussian** is a computational shortcut for the Phase-2 quantile: `(μ_g, σ_g²)` → `μ_g + Φ⁻¹(q)·σ_g`
  in closed form; agrees with the Beta away from 0/1.

### 1.4 Files & changes

- `calibration/density_model.py` — `node_gdna_density` also returns per-node `σ_g²`: track each node's
  `d_L`, `d_R`, their Poisson count noise, and run-position `t`; pool the 2-anchor `(mean, disagreement)`
  pairs, fit `α_impute` (single robust slope), evaluate `σ²_impute` per node, add the Poisson floor.
- `calibration/density_model.py::NodeDensity` — add `gdna_frac_var: np.ndarray` (or `count_concentration`).
- `calibration/result.py` — optionally surface a library-level `α_impute` QC scalar.
- No C++/accumulator change. No config knob (α_impute is fit; the floor is the Poisson/interpolation form).

### 1.5 Unit tests

- `α_impute` recovery: synthetic 2-anchor pairs drawn from `var = α·μ²` → fit recovers `α`.
- Monotonicity: disagreeing anchors → larger `σ_g²` than agreeing; 1-anchor > 2-agreeing; 0-anchor maximal.
- Linear-interpolation variance: peaks mid-run, equals endpoint variance at `t∈{0,1}`.
- Beta concentration: `μ_g(1−μ_g)/σ_g² − 1`, and the Gaussian quantile matches away from 0/1.

### 1.6 Validation (the acceptance test is new)

The **variance-calibration check** (distinct from net-leak): on the nascent+capture suite, does a high
`σ_g²` actually predict high per-node imputation error against the **oracle** gDNA fraction? Bin nodes by
predicted `σ_g`, plot vs realized |μ_g − truth| — they should track. Net-leak should not regress (Phase 1
alone doesn't move the mean; the blend weight stays `(2κ−1)²` unless §1.7-Q2 says otherwise).

### 1.7 Open questions

1. **Floor `σ²_floor`** — is the Poisson + interpolation form enough, or do we need an explicit floor (a
   knob)? Lean: derive the floor from the anchors' own count noise (no free constant).
2. **Does the blend weight become uncertainty-aware now?** With a count variance we *could* set
   `w = I_strand/(I_strand+I_count)`, `I_count = 1/σ_g²`. But the current `(2κ−1)²` is bias-robust (it
   doesn't trust a confident-but-biased count). **Lean: keep `(2κ−1)²` for the blend; use `σ_g²` only for
   the Phase-2 quantile.** Decide explicitly.
3. **Parametric NB (`α·μ²`) vs non-parametric kernel `f(μ)`** for the fit. Data says `μ²` → lean NB.

### 1.8 Magic-number audit

`α_impute`: fit. Poisson floor: parameter-free. Interpolation variance: parameter-free. No fitted constant.

---

# Phase 2 — False-positive-rate quantile knob

### 2.1 Problem & goal

Deliver the **one user knob** for the gDNA↔RNA confidence (your stated requirement: "a false-positive rate
threshold… a beautiful, intuitive user parameter… I will tolerate a 5% FP rate"). Under an asymmetric
gDNA/RNA loss, the Bayes-optimal point estimate of a posterior **is a quantile** of it. So the knob is a
quantile `q` of the per-node combined posterior: `q=0.5` (median) = neutral default; `q→1` = FP-averse
(call more gDNA → less gDNA→RNA leak, more RNA→gDNA siphon); `q→0` = the reverse.

This is the lever for the otherwise-irreducible single-exon residual: it cannot un-bias the mean, but it
lets the user **choose the operating point** on the leak↔siphon curve, uncertainty-aware (it bites hardest
where the posterior is widest — exactly the uncertain single-exon/capture nodes).

### 2.2 Design — the combined posterior and the quantile

Per node, combine the strand BB posterior (Phase 0, where strand-observable) and the count posterior
(Phase 1), weighted by the blend weight `w`. Cheapest faithful form — **Gaussian**:
```
μ = w·μ_strand + (1−w)·μ_count
σ² = w²·σ²_strand + (1−w)²·σ²_count
g(q) = clip(μ + Φ⁻¹(q)·σ, 0, 1)            # read g at quantile q instead of the bare mean/median
```
- Count-routed nodes (`w=0`, unstranded/AMBIG/single-exon): `g(q) = clip(μ_count + Φ⁻¹(q)·σ_count)` — this
  is where the Phase-1 count variance is essential (no strand posterior exists there).
- Strand-routed nodes (`w≈1`): the strand BB variance dominates; `q` still applies.
- `q=0.5` ⇒ `g(q)=μ` ⇒ reproduces Phase-0/1 numbers exactly (a clean no-op default).

### 2.3 The user-facing parameter

- `config.CalibrationConfig.gdna_deconv_quantile: float = 0.5`, plumbed via `cli.py`.
- **Decision-theoretic semantics:** `q` ↔ the tolerated asymmetric gDNA/RNA loss ratio. A thin layer can
  map a user-facing **FP-rate target** (e.g. "≤5%") → `q`, either analytically (from the loss ratio) or by
  empirical calibration on the benchmark. Decide in §2.6.
- **Retire `gdna_strand_llr_bias`:** the quantile is the uncertainty-aware generalization of the fixed
  log-odds shift (the LLR bias moves every node equally; the quantile moves wide-posterior nodes more).
  Keep `q` as the single knob; remove the LLR bias (or alias it).

### 2.4 Files & changes

- `calibration/strand_deconv.py::_deconv_per_node` — after forming `(μ, σ)` per node, read `g = g(q)`
  instead of the median; `q` from config.
- `config.py` / `cli.py` — add `gdna_deconv_quantile` (default 0.5); remove `gdna_strand_llr_bias`.
- `calibration/result.py` — record `q` in the config echo.

### 2.5 Unit tests

- `q=0.5` reproduces the Phase-1 mean exactly (bit-stable golden).
- Monotonic: `g(q)` increases with `q`; `q>0.5` raises gDNA on wide-posterior nodes more than narrow.
- Count-routed node with no strand: `g(q)` uses the count Gaussian only.

### 2.6 Validation

Sweep `q ∈ {0.5, 0.7, 0.9, 0.95}` on both benchmarks; confirm `q` monotonically trades gDNA→RNA leak for
RNA→gDNA siphon, that the trade is **uncertainty-aware** (moves the wide-posterior single-exon/capture
nodes, not the confident ones), and `q=0.5` reproduces Phase-1. Build the **FP-rate → q calibration curve**
(empirical: leak/siphon vs `q`) so the user-facing "5% FP" maps to a concrete `q`.

### 2.7 Open questions

1. **FP-rate→`q` mapping:** analytic (loss-ratio) vs empirical (benchmark-calibrated). Lean: ship the raw
   `q` knob first; add the FP% convenience layer once the curve is measured.
2. **Combined posterior:** Gaussian approx vs explicit Beta⊗BB mixture. Lean: Gaussian (closed-form, agrees
   away from 0/1).
3. **Quantile of the combined posterior vs quantile-then-blend.** Lean: combine to `(μ,σ)` then quantile.

### 2.8 Magic-number audit

`q` is a user preference (default 0.5, neutral), not a fitted constant. `Φ⁻¹` is exact. No fitted constant.

---

# Phase 3 — Strand-resolved nascent removal (3-term splice fraction + carry-over)

### 3.1 Problem & goal

The shipped 2-term boundary fraction lumps crossing-unspliced as **gDNA + nascent** → nascent is counted as
gDNA. The nascent suite shows the resulting siphon (nascent-on net leak drops / goes negative — RNA(nascent)
being eaten as gDNA; e.g. off_0.50 +1.0%→−4.9%, on_0.99 +5.8%→+1.1%). **Goal:** where a gDNA/RNA split of
the unspliced is available, use the **3-term** fraction so nascent moves to the RNA side:
```
f = (U_gDNA / E^gDNA) / ( U_gDNA/E^gDNA  +  (U_RNA + S)/E^RNA )      # pure gDNA / total
```
The split comes from **(a) strand pre-cleaning** (sense/antisense vs κ deconvolves the crossing unspliced
into gDNA at sense ½ and nascent at sense κ) or **(b) carry-over** from a neighbour in a left→right sweep.

### 3.2 Design — strand-resolved sweep state (per §5.1 of `count_mean_bias_design.md`)

The carry must be **per strand**, for *identifiability*, not bookkeeping. In an opposite-strand overlap
(AMBIG, e.g. `E+·I−`) the unspliced decomposes as `n_pos = gDNA/2 + nascent⁺`, `n_neg = gDNA/2 + nascent⁻`
— 2 equations, 3 unknowns, **locally under-determined** (why AMBIG is excluded from the strand module).
Carrying one strand's RNA forward from a single-strand neighbour where it *was* identifiable supplies the
missing constraint: `gDNA = 2·(n_pos − nascent⁺)`. So the sweep state is the vector `{gDNA, RNA⁺, RNA⁻}`
(mirrors the accumulator's existing 4 channels — no new information, just don't sum away the strand). The
spliced strand is reliable everywhere (motif-oriented); only the unspliced strand needs the carry.

**Carry-over sweep (the folded "4-mean.2"):** runs = maximal stretches between count-observable anchors;
seeds = regions with a direct 3-term split; bidirectional forward+reverse fill of `{gDNA,RNA⁺,RNA⁻}` along
the run (mirrors the existing density run-fill). Per-region prior still collapses to the scalar
`gDNA/(gDNA+RNA⁺+RNA⁻)` — strand resolution is internal to the sweep.

### 3.3 Acyclicity

Order stays: `κ → strand overdispersion → splice 3-term (uses κ + the fitted overdispersion to split the
crossing unspliced) → blend`. The 3-term consumes already-fitted κ/overdispersion and the carried state; it
does **not** refit the density. Assert no density→clean→density loop.

### 3.4 Scope & honest expectations (why this is last)

On the **current** suites the measurable impact is **small**: (i) stranded libraries route exons through
the strand module (which already sees nascent as sense=RNA), so the 3-term only refines the ~4% count share;
(ii) unstranded libraries cannot split nascent from gDNA at all (no strand, no junction) — that residual is
**irreducible** (off_0.50 −4.9% siphon is the accepted floor); (iii) the synthetic genome has **~0 AMBIG**
non-observable regions, so the per-strand machinery isn't exercised. Phase 3's real value is **real-world
antisense (AMBIG) + nascent at intermediate strand specificity** — which needs a **richer validation
genome** (add antisense overlap + nascent). So Phase 3 is correctness/completeness; sequence it last, and
build the AMBIG validation genome as part of it.

### 3.5 Files & changes

- `calibration/splice_junction.py` — extend `region_splice_gdna_frac` to: (i) strand-deconvolve the crossing
  unspliced per side into `{gDNA, RNA⁺, RNA⁻}` (using κ + overdispersion); (ii) the bidirectional carry-over
  sweep over runs; (iii) emit the 3-term fraction where a split exists, else the 2-term, else fallback.
- `calibration/calibrate.py` — pass κ / overdispersion into the splice step (already computed upstream).
- Possibly a new `scripts/sim/configs/` antisense+nascent suite for validation (§3.4).

### 3.6 Unit tests

- 3-term recovers pure gDNA when nascent is present and a split is supplied (extends Test A's
  `boundary_gdna_fraction` 3-term cases).
- AMBIG identifiability: constructed `E+·I−` overlap; carrying `nascent⁺` from a `+`-only neighbour recovers
  `gDNA` and `nascent⁻` (the 2-eq/3-unknown resolution).
- Carry-over sweep: a run with one seeded region propagates the split to its non-seeded neighbours;
  unreachable runs keep the fallback.

### 3.7 Validation

Nascent suite: the **stranded** nascent-on conditions' nascent→gDNA siphon should shrink toward the
nascent-off baseline (per-node gDNA estimate matches the oracle *gDNA-only*, not gDNA+nascent). On a new
antisense+nascent genome: AMBIG regions deconvolve correctly. Unstranded+nascent unchanged (irreducible).

### 3.8 Open questions

1. **Build the per-strand state from the start** (recommended — it *is* Phase 3's core) vs a scalar
   stopgap. Lean: per-strand from the start.
2. **The richer validation genome** (§3.4): extend the nascent suite with antisense overlap, or a dedicated
   AMBIG suite? Needed to show benchmark impact (unit tests prove correctness, not impact).
3. **Unstranded+nascent floor:** accept the −4.9% siphon as documented-irreducible, or attempt anything?
   Lean: accept + document.

### 3.9 Magic-number audit

No new constants — the 3-term uses the already-fit κ, overdispersion, and FL effective lengths; the
carry-over reuses the existing run-fill structure.

---

## Part III — Cross-cutting

- **Single-exon hard floor.** The 308/325 single-exon/no-junction regions under unstranded+capture are the
  irreducible core (strand blind, count capture-biased, no spliced reference). Phases 1–2 give honest
  uncertainty + a user FP-rate control over them; no phase claims to un-bias them. Document this as the
  accepted limit of the method.
- **Validation harnesses.** 20-condition gDNA suite (`gdna_benchmark_5mb`); 8-condition nascent suite
  (`nascent_benchmark_5mb`). Net-flow needs the analysis run **without** `--skip-frag-analysis`. Phase 1
  adds the variance-calibration check; Phase 3 needs a new antisense+nascent genome.
- **Commit discipline.** One benchmark-gated commit per phase; regenerate golden on intended output changes.
