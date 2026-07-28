# The gDNA hyperprior — implementation plan to PRODUCTION (Role B)

**Owner-directed, 2026-07-27.** Decision: retire the Role-B δ-pin `DensityNPMLE` fit and build a production
landscape + projection. This is the working plan. Evidence base: `gdna_hyperprior_resurrection.md` §4b
(all 48 scenarios, regenerated end-to-end at HEAD).

**Secondary goal, owner-stated:** leave the codebase *healthier* — clean, concise, readable, maintainable.
Cleanup is workstream **W0**, sequenced first, not bolted on afterwards.

---

## 0. Scope — replaced, preserved, untouched

| surface | today | after |
|---|---|---|
| `calibrate._fit_gdna_hyperprior` | returns `DensityNPMLE` | returns `GdnaLandscape` (new) |
| `bp_solver:258` `global_lp` | `prior.logprior(grid, M, E)` | **same mechanism**, new object |
| `simplex_logodds._gdna_arm` | fitted prior **REPLACES** `½·log f_g` | reference **+** fitted logP (additive) |
| `npmle.DensityNPMLE` | two roles | **RETIRED, NOT DELETED** — enrichment/QC only |

* **Role A is already dead in the solver** — no `.project(` call exists outside `npmle.py`;
  `enrichment_prior` survives only for `CalibrationDiagnostics.from_prior` and toy injection. This work
  cannot disturb σ²_transfer.
* **`DensityNPMLE` is preserved** (owner): a powerful construction that may resurface. Delete only at full
  production, once we are sure it is not needed.
* **No projection-and-anchor in ψ, and no `w_anchor`.** An earlier draft proposed collapsing the landscape
  to `(μ*, v*)` plus a Gaussian pull. Dropped: (i) a Gaussian anchor **destroys the multi-modality** this
  effort exists to recover, and (ii) `μ*` conditions on the node's own counts, which ψ already carries —
  the damping factor was a fudge for a formulation error. ψ takes the **full `log P` along the ray**, and
  the posterior *is* the projection. **The shipped mechanism was never the defect** — the crush dissection
  proved it ("the crush is hop-1", the fit).
* **Prior STRENGTH is a separate, legitimate knob** (owner): the landscape is fitted from *biased* pass-0
  output, so tempering it is robustness, not a fudge — it is what lets real data overcome a wrong prior.
  Default 1 (exact Bayes); optimized in W4.

---

## 1. W0 — CODE HEALTH (do first)

**Done 2026-07-27:** deleted 16 superseded debug scripts (2,225 lines; `scripts/debug/` hyperprior layer
30 → 13 files) per `scripts/README.md`'s delete-don't-archive policy; promoted the landscape recipe out of a
bare `def recipe(s)` loaded by `exec(open(...).read())` into importable, lintable
`gdna_explore_lib.{recipe, recipe_substrate, recipe_tau_prior}` — **verified byte-identical over all 48
scenarios**. Ruff clean, 1202 tests pass.

**Remaining, sequenced with the workstreams that touch them:**

| item | evidence | when |
|---|---|---|
| `\| gonly` clause in `_fit_gdna_hyperprior` is **dead** | verified 32/32: with `background_floor=True` the `~intergenic` mask removes exactly the `gonly` regions | W1 (rewrites this function) |
| `config.gdna_prior_additive` — one caller, default `False`, selects a dead branch | grep | W4 |
| the `tau_prior` **gate** — 3 bare constants, a cliff | §2 | W1 |
| the projection's `hup=0.70, hdn=0.02, cap_up=1.0, h_proj=0.15` | §4 | W3 |
| consolidate the review harness (`scratchpad/gdna_*.py`) into one place | 5 scripts, overlapping | W1 |

**Standing rule for this work:** every refactor step carries a behaviour-identity check (byte-identical over
the cached suites) *before* any tuning is layered on it, so a behaviour change is never confused with a
modelling result.

---

## 2. W1 — THE SUBSTRATE: which nodes train the prior (from scratch, no assumptions)

Owner: *"a single magic number precision cutoff is NOT the right path forward... start from nothing and
build up our understanding of the behavior from scratch."*

### 2.1 The `tau_prior` gate — what it is, and why it goes

```python
tau_prior = clip(std over intergenic REGION nodes of log10(g/E), 0.15, 0.70)   # a population spread
admit a BOUNDARY node  iff  sqrt(Var(log f_g))/ln10  <  tau_prior
```

The *intent* is an information criterion — admit a boundary only if it pins its density tighter than the
population varies. The *implementation* fails on four counts:

1. **The "derived" threshold is a constant most of the time.** Measured: it clips to the **lower bound 0.15
   on 19/32 conditions** and **never** reaches the upper bound; observed range 0.150–0.300 against a clip
   window of [0.15, 0.70].
2. **That 0.15 is the KDE bandwidth `S0`** reused for an unrelated purpose — one constant doing two jobs,
   with no argument that the right smoothing scale is also the right admission threshold.
3. **Three bare constants** (0.15, 0.70, 0.3-fallback) and a **hard cliff**, against the standing
   "no clamps/cliffs" rule.
4. ⚠ **It is NOT redundant with the continuous weight — it is consequential.** Hypothesis tested and
   **refuted**: the gate-rejected boundaries carry **median `w ≈ 0.37`** and **12–16 % of all boundary
   weight**. So it discards a sixth of the boundary evidence by constant. Removing it is a real change and
   must be measured.

**→ The gate is removed. What replaces it is the subject of this workstream, and may be "nothing".**

### 2.2 The taxonomy — four distinct reasons a node might be excluded, which must not be conflated

| axis | question | correct form | example |
|---|---|---|---|
| **circularity** | would training on it make the prior predict itself? | **structural exclusion** | **AMBIG** — the two-root ambiguity the prior exists to resolve. Already excluded, and for this reason |
| **identifiability** | is the composition fixed by construction? | **structural inclusion** — it is the anchor | intergenic / no-admissible-RNA-strand nodes = the depleted floor |
| **precision** | how well does this node know its density? | **continuous weight**, never a cutoff | `w = ref/(v+ref)` |
| **geometry** | is it the same measurement? | **separate treatment or matched reference** | boundaries: *crossing* eff-lengths, not *contained* |

The current gate conflates precision with admission. AMBIG's exclusion is *circularity* and is not a
threshold at all — it is a structural fact, which is why it is the one exclusion that is clearly right.

### ⭐ 2.3 W1 RESULTS — LANDED 2026-07-27 (`scratchpad/gdna_w1_substrate.py`)

All arms share the same kernel and the same weights unless the arm *is* about weights, so no substrate
result is confounded with a weighting change. Boundaries INCLUDED throughout. Δ vs the new default:

| arm | ambig ΔEMD | quick ΔEMD | verdict |
|---|---|---|---|
| **removing the `tau_prior` gate** | **−0.0034** | **−0.0255** | ✅ **FREE — landed.** Recall and spurious rate unchanged |
| + AMBIG admitted | −0.0126 | +0.0200 | ⛔ keep excluded — see below |
| − the zero-count structural anchor | **+0.3192** | **+0.6015** | ✅ **critically load-bearing** — keep |
| − `gonly` (intergenic) | +0.0116 | −0.0024 | ~neutral on EMD, but spurious 0.24 → 0.39 — keep |
| **FLAT weights** (precision ignored) | **+0.1452** | +0.1363 | precision weighting is load-bearing |
| **precision CUTOFF** (a hard threshold) | **+0.2487** | **+0.4297** | ⭐ **WORSE THAN IGNORING PRECISION** |

**The headline vindicates the taxonomy.** A hard precision cutoff is *worse than having no precision
information at all* (+0.249 vs +0.145 on ambig; +0.430 vs +0.136 on quick). Precision belongs in a
continuous weight; expressing it as an admission decision is strictly harmful. That is the general reason
the `tau_prior` gate had to go, not merely that its constants were unjustified.

**AMBIG: the apparent win is SUBSTRATE-MATCHING, not better estimation.** `oracle_landscape` includes AMBIG
regions while the fit excludes them, so admitting AMBIG partly just matches the target. Scored against the
non-AMBIG population the prior is actually applied to, it **loses on both suites**:

| arm | vs oracle ALL regions | vs oracle NON-AMBIG regions |
|---|---|---|
| AMBIG out | 0.205 | **0.175** |
| AMBIG in | **0.193** | 0.180 |

So the circularity exclusion now has empirical support as well as principle. ⚠ EMD cannot settle it either
way — circularity is about the *refit's* validity, which only an end-to-end test (W4) can measure.

**Landed:** `recipe_substrate(s, mk)` is gate-free; `recipe_tau_prior` and its three constants are deleted;
`recipe_weights` is exposed so precision can be ablated independently of substrate. New baseline —
ambig **EMD 0.205, sd/oracle 1.00**; quick **0.252**. Ruff clean, 1202 tests pass.

### 2.4 The remaining study

Sweep candidate substrates and score each against the oracle-fitted prior. Axes: node class
(region / boundary / AMBIG / intergenic / zero-count structural), precision treatment (continuous weight vs
cutoff vs flat), and geometry handling.

### ⭐ W1a — DONE 2026-07-27. Both scoring confounds measured (`scratchpad/gdna_w1a.py`)

Instruments added to `gdna_explore_lib`: `oracle_pointmass` (the population at minimal smoothing),
`oracle_landscape(s, sel)` (now takes any node set, so a matched reference is possible), and shared
`spread` / `modes`. `recipe(s, sel=...)` takes a substrate override so variants A/B **with the weighting
held fixed** — verified byte-identical at the default.

**Confound 1 — THE SMEAR: measured, and it is NOT a problem. Deconvolution is not warranted.**

| subset | sd pointmass | sd landscape | inflation | EMD between the two oracles |
|---|---|---|---|---|
| **G > 0** (nodes with real gDNA) | 0.800 | 0.830 | **1.04** | **0.033** |
| G = 0 (zero-gDNA nodes) | 0.728 | 0.990 | 1.36 | **1.555** |

Where gDNA exists the two references are **the same distribution** (EMD 0.033 against fit errors of
0.21–0.28). **All** the disagreement is on zero-count nodes — and there the *pointmass is the wrong
instrument*: it pins a zero-gDNA node at `log10(1/E)`, inventing a location, whereas the Poisson landscape
spreads it downward (`P(0|ρE) ∝ e^{−ρE}` is monotone decreasing) = "ρ could be anything below the
resolution wall", which is the honest statement and exactly why the zero-native framing was adopted.
→ **W2 decision settled: ACCEPT the smear, do not deconvolve.** ⚠ And use `oracle_pointmass` **only on
G > 0 nodes**; on zeros it is not a reference, it is an artefact.

**Confound 2 — BOUNDARIES: the comparison is now fair, and the answer is genuinely AMBIGUOUS.**

*B1 — boundaries are a measurably different population (ORACLE truth, no solver involved):* their true
gDNA density sits **+0.28 to +0.42 decades ABOVE** regions on gDNA-bearing conditions (+0.55…+1.04 on
zero-gDNA) and is **3–7× tighter** (sd 0.11–0.75 vs region 0.38–1.29). So including them genuinely imports
a shifted, narrow population — the benefit on offer is sample size, and the cost is that bias.

*B2/B3 — cost of dropping boundaries, weighting held fixed:*

| suite | vs REGION (the production target) | vs MATCHED | vs POINTMASS |
|---|---|---|---|
| ambig (10 Mb) | **+0.0232** (dropping is worse) | +0.0143 | +0.0170 |
| quick (5 Mb) | **−0.0453** (dropping is BETTER) | +0.0123 | −0.0113 |

**The sign flips between suites on the production-relevant target.** Boundaries always help the fit describe
its *own* training set (matched column, both positive) — which is unsurprising and not the question.
⚠ An earlier version of this experiment dropped the reliability weight along with the boundaries and
therefore reported boundaries as uniformly good; that was an artefact of the harness, not a result.

→ **Not settled by these data.** Binary include/exclude is the wrong frame: the honest reading is that
boundaries are a *different measurement* (crossing vs contained geometry) with a known offset, so W1 should
test *representativeness-aware* inclusion, not a coin flip.

⚠ **Known counter-evidence — exclusion can HURT.** `archive/solve_gate_design.md` records a stronger
exclusion rule derived, implemented and **empirically refuted** (+0.010 r0 / +0.025 r1). And a node
reporting a genuine **enriched** mode may be exactly the one an exclusion rule drops — the enriched mode is
the scarce thing.

**Gates.** Beat the current substrate on TV-to-oracle-fitted-prior, on **both** suites, with the enriched
mode's mass/location no worse. Any exclusion rule must justify itself against "include it, weighted".

---

## 3. W2 — THE KERNEL: per-node width + a population resolution

**The measured defect.** The recipe's resolution is *entirely* per-node (each node's own Poisson count), so
a 50k-fragment node contributes a ~0.002-decade near-delta. Where nodes are dense these merge; where sparse
they stand alone and read as modes. **Mode count is therefore a function of training-set size** — 3.0 → 5.6
as `n_train` 1217 → 121 against an oracle of 3 — and the spurious-mode rate is 0.21 on the 10 Mb suite vs
**0.48** on the 5 Mb one.

**What bandwidth buys (owner).** The spurious spikes are largely *fragments of the enriched mode*. Smoothing
should consolidate them, so the reported mode becomes the centroid instead of whichever fragment is tallest.
Prediction to test: outliers like `gdna1 verystrong`'s **−0.90** shift collapse toward the typical
−0.15…−0.20 band. What smoothing **cannot** do is move that centroid — the residual band is genuine bias.

```
    h_i²  =  h_within,i²   +   h_pop²
             per-node          POPULATION RESOLUTION  (entirely missing today)
```

* `h_within` — keep the **zero-native** kernel (a zero-count node must remain the exact "gDNA ≈ 0" anchor;
  a point-estimate KDE floors at `1/E` and cannot express it). `fit_poisln` already does Poisson ⊗ belief
  width and is unused — the starting point.
* `h_pop` — chosen **empirically, no assumed rule** (owner). Sweep fixed / normal-reference
  `∝ σ̂·n_eff^(−1/5)` / k-NN adaptive.

**Negative binomial (owner's question).** `NB(0 | μ, α) = (1+αμ)^(−1/α) > 0`, so it **is** zero-native, and
`sd(log ρ) ≈ √(1/k + α)` gives a width floor that never vanishes — exactly what the near-deltas need. But α
and `h_pop` are **different quantities and must stay separately named**: α is a physical claim about count
noise and **does not shrink with n**; `h_pop` is a resolution choice and **does** (`~n^(−1/5)`). Merging
them lets a smoothing choice be fitted and reported as physical overdispersion — the `ω_graft` trap.
⚠ **α cannot be validated here:** the suite is **Poisson by construction** (ω = 0), real ω ≈ 0.02, and
without replicates α is only weakly identifiable. **Toy validates `h_pop`; only real data speaks to α.**

**Scoring bandwidth on real data, where there is no oracle.** EMD/mode-recall do not exist on cfRNA. Use
**held-out predictive likelihood** — fit on a subset, score the probability assigned to held-out nodes'
counts. No oracle, no constants; penalises over-smoothing and spurious modes alike. Protocol: confirm on the
toy that it selects the same `h` the oracle metrics prefer, then use it as *the* criterion on real data.

**"Deconvolution" vs "smear", concretely.** If every intron truly has `ρ_g = 0.01`, the population is a
spike of zero width. Node A (10 fragments) estimates 0.008, node B (10,000) estimates 0.0101. Summing
per-node posteriors yields a landscape with real **width** — but the truth had none. That excess is
measurement error leaking into the population estimate: the *smear*. It matters because the landscape is a
prior and the solve *also* accounts for this node's error — counting it twice, leaving the prior too broad.
Settle deconvolve-vs-accept on the `oracle_pointmass` measurement (HANDOFF_15's independent anchor says the
achievable prior is 33 % too broad).

### ⭐ W2 RESULTS — 2026-07-27 (`scratchpad/gdna_w2_bandwidth.py`; figures `gdna_w2_ladder.png`, `..._curves.png`)

**A GLOBAL `h_pop` is REFUTED — measured and visually explained.** Every aggregate metric prefers `h = 0`
(EMD 0.205 → 0.372 over h = 0 → 0.60; held-out log-likelihood −3.82 → −4.39 monotonically). The reason:
**the oracle's depleted mode is itself sharp**, and the fit at h = 0 is already *under*-peaked there —

| condition | oracle peak | fit h=0 | fit h=0.30 |
|---|---|---|---|
| gdna300 ss0.50 capON | 0.193 | 0.141 | **0.026** |
| gdna300 ss0.50 capOFF | 0.553 | 0.341 | **0.036** |

— so a global width flattens the one region the fit gets right (7–9× too flat) to fix the sparse enriched
tail. One width cannot serve regions whose required resolution differs ~10×.

**ABRAMSON sample-point adaptivity (`h_i ∝ f(a_i)^(−1/2)`) is ~NEUTRAL** (ambig 0.216 → 0.214; quick
0.254 → 0.242). Diagnosis: a spike **is** high pilot density at its own location, so the rule keeps it
narrow — self-reinforcing exactly what needs merging.

**✅ ADOPTED — adaptive resolution by NEAREST-NEIGHBOUR SPACING.** `h_i = scale · dist(a_i, k-th nearest
neighbour)`, `k = √n`. Spacing is the quantity that actually decides whether two kernels merge.

*It fixes the overfitting by construction* — fewer nodes ⇒ farther neighbours ⇒ wider kernels, no tuning:

| kernel | n=100 % | 50 % | 25 % | 10 % | slope | (oracle = 3) |
|---|---|---|---|---|---|---|
| h = 0 (today) | 4.0 | 4.8 | 7.0 | **9.6** | **+5.6** | |
| **kNN 0.5** | **3.0** | 2.4 | 2.0 | 2.0 | **−1.0** | exact at full data |

and it costs nothing: EMD **ambig 0.205 → 0.205** (tie), **quick 0.252 → 0.241**, and it is the *best*
arm on ambig by the oracle-free criterion too. On the zero-gDNA scenario it cuts invented modes 4.0 → 2.0.
The ladder figure shows why: the kNN curves lie **on** the h=0 curve through the sharp bulk while damping
the spurious wiggles in the sparse tail — the thing global smoothing could not do.

### ⚠ The oracle-free criterion has a BIAS, and it matters for real data

Held-out predictive likelihood and EMD **agree on ambig** (both pick kNN 0.5) but **disagree on quick**
(LL picks h = 0, EMD picks kNN 2.0). The reason is structural: held-out LL scores the landscape as a
**predictive** object, which wants population ⊛ noise, whereas its use as a **prior** wants the population
itself. **So held-out LL systematically under-smooths and is not a clean proxy for prior quality.** It
remains the only criterion computable on cfRNA, so treat it as a *lower bound* on the right bandwidth and
carry the toy-derived `scale` across, rather than re-selecting on real data.

**⚠ Metric caveat found here:** the *spurious-mode RATE* is unreliable under smoothing — it is a fraction of
fitted modes, so merging two real modes into one that lands between them both raises the numerator and
shrinks the denominator. **Use the dose-response slope as the overfitting gate**, not the spurious rate.

**Still open in W2:** `scale` (0.5 adopted provisionally — best on ambig by both criteria, never worse than
h=0; quick prefers 2.0 on EMD) needs the end-to-end W4 test to finalise, and the NB `α` remains
unmeasurable here.

**Gates.** (1) mode count flat in `n_train` — *the* overfitting gate ✅; (2) `none_*` stays a monotone
decay ✅; (3) recall not below today's; (4) support covers the oracle's mass.

---

## 4. W3 — THE PROJECTION (as the scoring instrument) and PRIOR STRENGTH

Measured, against `IDEAL(d) = E[true ρ_g | observed = d]` — the best any function of `d` can do:

| projection | mean gain | median | negative on |
|---|---|---|---|
| symmetric `h=0.15` | +0.11 | +0.09 | 0/32 |
| asymmetric (`hup`,`hdn`,`cap`) | +0.14 | −0.06 | **19/32** |
| **posterior, `σ = √Var(log f_g)/ln10`** | **+0.21** | **+0.17** | **0/32** |

The IDEAL curve is largely **flat** — the observation is often uninformative and the job is **shrinkage**,
not translation. Symmetric sits on the identity (inert); asymmetric is identity **+ a constant ≈ 0.6** and
fails on 19/32 because it is an unconditional bias correction. **Adopt the posterior form; it adds no
constant** (σ is already computed). Retire `hup`/`hdn`/`cap_up`/`h_proj`.

**σ-calibration (do this).** Gain is +0.21, not 1.0 — either the landscape is wrong or `σ_obs` is
mis-calibrated. Scale σ by `c ∈ {0.5,1,2,3,5}`, find `c*` per stratum. `c* ≈ 1` ⇒ declared variances are
honest in the units the prior consumes and W2 is the lever; `c* ≫ 1` ⇒ a **direct measurement of pass-0
over-confidence on a new axis**, feeding the existing `z2` programme. Either way it is information.

**Prior strength.** Temperature on the ψ prior term, default 1. Optimized in W4, and the metric is *"does
real data still dictate the solution"* — a wrong prior must be overcome by evidence.

### ⭐ W3 RESULTS — 2026-07-27 (`scratchpad/gdna_w3_projection.py`, figure `gdna_w3_projection.png`)

On the W2 kNN landscape, the ranking holds: **POSTERIOR** gain **+0.208** mean / **+0.177** median, negative
on **2/32**; asymmetric +0.128 / **−0.073**, negative on **19/32**; symmetric +0.098 / +0.068 (inert).
Adopt the posterior form.

**⭐ THE σ-CALIBRATION IS THE HEADLINE, AND IT IS LARGE.** Scaling the declared width by `c`:

| stratum | c=1 | c=2 | c=3 | **c=5** | c=8 | c=20 | c=60 | c\* |
|---|---|---|---|---|---|---|---|---|
| ambig all | 0.208 | 0.437 | 0.637 | **0.749** | 0.584 | 0.132 | −0.259 | **5** |
| ambig unstranded | 0.210 | 0.441 | 0.680 | **0.836** | 0.685 | 0.334 | 0.285 | **5** |
| ambig capture ON | 0.156 | 0.294 | **0.393** | 0.300 | −0.218 | −1.415 | −2.389 | **3** |
| ambig capture OFF | 0.206 | 0.478 | 0.708 | 0.940 | 1.048 | 1.167 | **1.233** | **→∞** |

**Inflating the declared width ~5× takes the projection from closing 21 % of the achievable gap to 75 %.**
The median declared width is only **0.07 decades** — so narrow the posterior is nearly inert, which is
exactly what the transfer-function figure shows (the blue curve sits on the identity).

The stratification is interpretable and not noise: **capture-OFF wants `c → ∞`** (its IDEAL curve is flat —
the observation is uninformative, so full shrinkage to the population is correct) while **capture-ON wants
`c ≈ 3`** (there the observation genuinely carries signal).

⚠ **Do NOT read `c* ≈ 5` as "pass-0 is 5× over-confident".** Optimal shrinkage responds to a node's TOTAL
error, and pass-0's error is substantially **BIAS** (the measured compression: depleted over-called +0.045,
enriched under-called −0.152). A declared *variance* does not contain the bias, so `c*` conflates genuine
variance under-declaration with un-modelled bias. That is coherent with §2.3/§4: the residual is directional.

### ⭐ W3b — WHAT `c*` IS ACTUALLY COMPENSATING FOR (derivation, 2026-07-27)

Three instruments, two of which misled; the third is the reliable one. Recorded because the wrong two are
tempting.

1. ⛔ **Adjacent-node pairs** ("gDNA is uniform, so neighbours share a true density, and their observed
   difference is pure noise"). **Premise fails:** adjacent EXON pairs differ by **0.916 decades in TRUE
   density** (different capture enrichment), so the subtraction removes the signal and the estimator has no
   power (implied c ≈ 0.3). Intergenic/intron pairs are too rare in the chain to rescue it.
2. ⛔ **`sd(z)` where `z = (obs − true)/σ_declared`.** Reads 3–33, suggesting a uniform 3–33× under-
   declaration. **It is driven entirely by outliers** — see 3.
3. ✅ **Quintiles of declared σ vs realised |error|.** The reliable instrument, and it says the declared
   widths are **well calibrated in the MEDIAN**:

| declared-σ quintile | Q1 | Q2 | Q3 | Q4 | Q5 |
|---|---|---|---|---|---|
| median declared | 0.0116 | 0.0510 | 0.0956 | 0.1300 | 0.2040 |
| median \|error\| | 0.0112 | 0.0327 | 0.0686 | 0.1128 | 0.1687 |
| **ratio** | 1.0 | 0.6 | 0.7 | 0.9 | 0.8 |

Ratios 0.6–1.7 across every quintile, tracking monotonically. **So the problem is NOT scale — it is a HEAVY
TAIL** of confidently-wrong nodes (which is precisely what the project's own `z2` programme measures).

**Consequence: uniform inflation is the wrong mechanism.** `c = 5` defends against the tail by widening
*everyone*, which (a) damages the well-calibrated majority — it makes **9/32 conditions worse**, vs 2/32 at
c=1 — and (b) at σ = 0.35 dec exceeds the **minimum oracle mode gap of 0.12 dec**, so it smears across real
structure. ⚠ Note the "modes are ~3 decades apart" intuition holds only for the depleted↔enriched
*extremes*; the enriched region is a cluster with 0.09–0.12 dec gaps.

**✅ The right mechanism is a HEAVIER-TAILED LIKELIHOOD, not a wider one.** A Student-t is exactly a Gaussian
marginalised over an *uncertain* width — which is the true state of knowledge here — so it keeps the median
node's narrow, mode-resolving kernel while automatically refusing to pin an anchor from a tail node:

| likelihood | mean gain | median | conditions made WORSE |
|---|---|---|---|
| Gaussian c=1 (today) | +0.208 | +0.177 | 2/32 |
| Gaussian **c=5** | +0.749 | +0.859 | **9/32** |
| **Student-t ν=2, c=1** | **+0.527** | +0.569 | **4/32** |
| Student-t ν=4, c=1 | +0.357 | +0.351 | 4/32 |
| Student-t ν=8, c=1 | +0.261 | +0.226 | 2/32 |

**70 % of the gain, less than half the harm, and σ is left alone.**

**ν is derivable, not tunable:** for a t-distribution the excess kurtosis is `6/(ν−4)`, so ν follows from
the measured kurtosis of the residual `z` by **method of moments** — the allowed pattern (κ, both strand
overdispersions, `ω_graft`). Fit it; do not pick it. **Build this before wiring W4.**

### ⭐ N3 — ν IS NOW FITTED, NOT PICKED (2026-07-27). `scratchpad/gdna_n3_{nu,gain}.py`

W3b's prescription was *"excess kurtosis `6/(ν−4)` ⇒ ν by method of moments"*. **That estimator is
unusable here, and the reason is structural, not numerical:** it exists only for ν > 4, so it can never
return the ν ≈ 2 the data actually want. Measured on the residual `z = (obs − true)/σ_declared`:

| suite | stratum | n | z == 0 | sd(z) | MAD-sd | excess kurt | **ν(kurt)** | ν(log-var) | **ν(MLE)** | **ν(HELD-OUT)** |
|---|---|---|---|---|---|---|---|---|---|---|
| ambig | all | 47 303 | 24.0 % | 2.43 | 0.90 | 149 | **4.04** | 1.21 | **2.20** | **3.0** |
| ambig | unstranded | 32 352 | 23.2 % | 2.58 | 0.86 | 156 | 4.04 | 1.21 | 2.20 | 3.0 |
| ambig | stranded | 14 951 | 25.9 % | 1.99 | 1.05 | 59 | 4.10 | 1.22 | 2.40 | 4.0 |
| quick | all | 14 779 | 28.5 % | 6.75 | 0.92 | 143 | **4.04** | 3.83 | **2.00** | **3.0** |

* ⛔ **ν(kurt) saturates at its own boundary (4.02–4.15) on every single stratum.** It is not measuring
  anything — the excess kurtosis is 40–280 and `4 + 6/κ` is pinned. **Do not ship it.**
* ⛔ **ν(log-var)** — `Var(log z²) = ψ'(½) + ψ'(ν/2)`, which is finite for every ν > 0 and so *should* have
  rescued the case — is destroyed by a feature of the data, not a defect: **20–40 % of `z` are EXACTLY 0**
  because observation and truth both floor at the one-count wall. Those are successes, not outliers, and
  `log z²` diverges on them. Reported on `z ≠ 0` it reads 1.05–4.40, i.e. unstable.
* ✅ **ν(MLE) = 1.6–2.8, and 2.20 / 2.00 on the two full suites** — an independent estimator landing on
  W3b's hand-picked optimum. But it consumes the ORACLE, so **it cannot run on real data.**
* ✅ **ν(HELD-OUT) = 3.0 on both suites** — fit the landscape on half the training nodes, score the other
  half's observed density under `P ⊛ t_ν(σ_i)`. **No truth anywhere: this is the one estimator deployable
  on cfRNA**, and it is the same protocol W2 adopted for the bandwidth. It sits one step light-tailed of
  the MLE — **exactly the bias direction W2 already documented** for held-out likelihood (it scores a
  *predictive* object, which wants population ⊛ noise).

**And the fitted value delivers** (gain vs `IDEAL(d) = E[true | observed = d]`, W2 kNN landscape):

| arm | ambig gain | median | mean bias | worse | quick gain | worse |
|---|---|---|---|---|---|---|
| Gaussian c=1 (W3 adopted) | +0.208 | +0.177 | +0.540 | 2/32 | +0.148 | 2/16 |
| Gaussian c=5 (⛔ rejected) | +0.749 | +0.859 | −0.616 | **9/32** | +0.557 | 4/16 |
| Student-t ν=2 (hand-picked) | +0.527 | +0.569 | **+0.121** | 4/32 | +0.338 | 2/16 |
| **Student-t ν=2.2 (oracle MLE)** | +0.504 | +0.538 | +0.159 | 4/32 | +0.323 | 2/16 |
| **Student-t ν=3.0 (held-out — DEPLOYABLE)** | **+0.425** | +0.423 | +0.274 | 4/32 | **+0.277** | 2/16 |
| Student-t ν=8 | +0.261 | +0.226 | +0.478 | 2/32 | +0.189 | 2/16 |

**The deployable fitted ν = 3 keeps 81 % of the hand-picked ν = 2's gain, at 2× the adopted Gaussian's,
with identical harm counts** — and the bare 5× σ inflation stays rejected (9/32 harmed, and it is the only
arm that over-shoots into negative bias). ν is a fitted parameter of the allowed pattern; **the constants
ledger gains no new entry.**

### ⭐⭐ CORRECTION — "a perfect landscape is worth only +0.03" IS WRONG (2026-07-27)

That figure (§2.3, and the `ROADMAP` banner) was measured with the **symmetric h=0.15 projection, which is
inert** — a projection that cannot move cannot benefit from a better landscape, so it measured the
projection's deadness, not the landscape's value. Re-measured with a *working* projection, substituting the
ORACLE landscape:

| landscape | Gaussian c=1 | Gaussian c=5 | **t ν=2** | conditions made worse (t ν=2) |
|---|---|---|---|---|
| our FIT (kNN) | +0.208 | +0.749 | **+0.527** | 4/32 |
| **ORACLE** | **+0.394** | +0.936 | **+0.748** | **0/32** |

**A perfect landscape nearly DOUBLES the gain at every setting and eliminates every harmed condition.**
So fit quality is a **first-order limiter** of the projection, not a +0.03 side-issue — the owner's
instinct ("the fit actually is what needs to be robust") is the correct read, and the earlier claim is
superseded wherever it appears.

### ⛔⛔ RETRACTED — "the enriched region is a cluster of closely-spaced modes". IT IS NOT.

**The ground truth is TWO modes, by construction** (owner, 2026-07-27): gDNA is uniform; capture partitions
nodes into captured and not-captured. Measured on ORACLE counts (`G > 50`, so `G/E` is precise):

| node class | p05 | median | p95 | **IQR** |
|---|---|---|---|---|
| intergenic | −1.31 | **−1.28** | −1.26 | **0.02** |
| intron | −1.33 | **−1.28** | −1.24 | **0.04** |
| exon | −1.30 | +1.41 | +1.68 | 0.76 |

**gDNA is uniform to within 0.02–0.04 decades** — the depleted level is a point mass. The exon population is
a mixture of that same depleted level (exons that caught no probe; note p05 = −1.30) and ONE broad enriched
mode of spread **0.32 dec**.

**The "3–8 oracle modes" were an artefact of the mode detector's prominence threshold**, which was set to
2 % of the density max and therefore counted noise wiggles:

| prominence | 0.02 | 0.05 | 0.10 | 0.20 | 0.35 |
|---|---|---|---|---|---|
| modes found on the enriched-only landscape | 7 | 5 | 3 | 2 | **1** |

**⛔ What this INVALIDATES (do not cite these):**
* the **mode-recall** and **spurious-mode-rate** columns throughout W1/W2 — they counted wiggles;
* the **"minimum mode gap 0.12 dec"** and the σ-ceiling derived from it. The real separation is
  depleted (−1.28, IQR 0.02) to enriched (+1.41) ≈ **2.7 decades**, so σ = 0.35 dec is entirely safe and the
  owner's original "enriched and depleted are light-years apart" was correct.

**✅ What SURVIVES** (none of it uses mode detection): every **EMD**-based result (W1 substrate ranking, W2
kernel ranking), the **dose-response overfitting gate** (wiggle count growing as the sample shrinks IS
overfitting, whatever we call the wiggles), the **projection gains**, the **σ-calibration**, the
**Student-t** result, and the **oracle-landscape doubling**.

**→ Metric fixed:** `gdna_explore_lib.two_component` scores the landscape the shape the truth actually has.
Re-scored on gDNA-bearing capture-ON conditions:

| arm | dep loc | dep wid | enr loc | enr wid | **enr mass** | EMD |
|---|---|---|---|---|---|---|
| **ORACLE** | −2.136 | 0.902 | 1.076 | 0.169 | **0.140** | — |
| fit: no resolution term | −2.148 | 0.944 | 1.008 | 0.185 | **0.056** | 0.306 |
| fit: kNN 0.5 | −2.146 | 0.947 | 1.034 | 0.225 | **0.056** | 0.304 |
| fit: kNN 2.0 | −2.146 | 0.959 | 1.195 | 0.373 | 0.058 | 0.297 |
| fit: GLOBAL h=0.30 | −2.133 | 0.981 | 1.108 | 0.279 | 0.052 | 0.348 |

### ⭐⭐ THE DEFECT IS THE ENRICHED COMPONENT'S **MASS**, AND NO KERNEL FIXES IT

* The **depleted** component is reproduced almost exactly (loc −2.148 vs −2.136, width 0.944 vs 0.902) —
  as expected, since gDNA is uniform and that is the easy half.
* The **enriched** component's location and width are close (1.008–1.034 vs 1.076; 0.185–0.225 vs 0.169).
* Its **MASS is 0.056 against a true 0.140 — we recover 40 %** — and **every kernel arm lands at
  0.052–0.058.** Bandwidth cannot fix it; it is set by the SUBSTRATE, the WEIGHTS, and pass-0's own
  enriched under-call.

**That is the target for the rest of this work.** It also explains the W3b oracle-vs-fit result from the
other side: the landscape the projection is shrinking toward is missing 60 % of its enriched mass, so the
projection cannot send a node there however well it is specified.

⚠ `kNN 2.0` wins EMD (0.297) but over-widens the enriched component **2.2×** (0.373 vs 0.169) — a good
example of why EMD alone is not sufficient and the two-component view is needed. **kNN 0.5 stands.**

### ⭐⭐⭐ THE DERIVATION: where the enriched mass actually goes (2026-07-27)

**Step 1 — bandwidth is ruled out ARITHMETICALLY.** For a node of true count `G` the Poisson kernel's width
on log-density is `~1/√G`. Enriched nodes carry `G` ≈ 3,400–9,100 ⇒ kernel sd **0.0045–0.0075 dec**, against
an observed enriched spread of **0.23–0.32 dec** — **30–70× larger**. No kernel can create or remove it.

**Step 2 — the two-mode ground truth is CONFIRMED; the apparent breadth is two artefacts.** Splitting
enriched exons by effective length:

| eff-length bin | Q1 (E=53) | Q2 (209) | **Q3 (442)** | Q4 (1187) | Q5 (1897) |
|---|---|---|---|---|---|
| median density | 1.01 | 1.46 | **1.44** | 1.39 | 1.38 |
| **sd** | 0.531 | 0.241 | **0.091** | 1.328 | 1.364 |

At E ≈ 442 the enriched population is **tight, sd 0.091** — one sharp mode, as the owner stated. The breadth
elsewhere is (a) **probe granularity** on small nodes (E=53 vs `probe_length` 120 — a node can sit wholly
inside or outside a probe) and (b) a captured/uncaptured **mixture** in the large bins (large exons include
uncaptured ones sitting at the depleted level). Neither is a broad mode.

**Step 3 — the deficit DECOMPOSES, and pass-0's errors are not in it:**

| | enriched mass | recovery |
|---|---|---|
| ORACLE | 0.140 | 100 % |
| our fit | 0.056 | 40 % |
| + flat weights | 0.082 | 58 % |
| **+ ORACLE counts, same substrate** | **0.081** | **57 %** |
| + AMBIG admitted (weights unchanged) | 0.070 | 50 % |

**Perfect counts buy NOTHING** (0.081 vs 0.082) — the deficit is not pass-0's accuracy. It is:
* **the reliability WEIGHT**, worth ~18 points, and it bites **only on unstranded data** (mean `w` for
  truly-enriched 0.528 vs 0.721 depleted when unstranded; **0.930 vs 0.930** when stranded);
* **the SUBSTRATE**, worth ~10 points on its own.

**Step 4 — the substrate leak is ENTIRELY AMBIG.** Of the truly-enriched nodes outside the training set,
**507/507 unstranded and 461/461 stranded are AMBIG** — zero from boundaries, zero from missing mass, zero
other. 32 % of all truly-enriched nodes are AMBIG.

### ⭐⭐⭐ N1 — THE DEFICIT IS NOW CLOSED, EXACTLY. IT IS A SUBSTRATE CENSUS EFFECT (2026-07-27)

**Figure: `figures/gdna_n1_massflow.png`.** Scripts: `scratchpad/gdna_n1_{massflow,ladder,factorial,where,arms,plot}.py`
— all pure numpy on the cached substrate, no `calibrate` re-runs.

The fitted landscape is *exactly* a weighted sum of per-node kernels, so the mass it places above a
threshold is an identity, not a model:

```
    M(t) = Σ_i  wt_i · q_i(t)          q_i(t) = the share of node i's own kernel above t
```

which decomposes over any partition of the training nodes. Evaluating it (13 gDNA-bearing capture-ON/VSTRONG
conditions, at the ORACLE's own split — the convention that reproduces §W3's table to 4 decimals):

| | ORACLE | FIT | ratio |
|---|---|---|---|
| **SUBSTRATE** — training weight on truly-enriched nodes | 0.1474 | **0.0612** | **0.415** |
| **PLACEMENT** — that weight's kernel share above the split | 0.936 | 0.817 | 0.873 |
| LEAKAGE — mass above the split from truly-depleted nodes | 0.0003 | 0.0034 | — |
| **enriched mass** | **0.1405** | **0.0561** | 0.40 |

**The deficit is ~entirely the SUBSTRATE term.** Placement is 0.873 — pass-0 puts truly-enriched nodes in
roughly the right place; it is the *census* that is wrong. And the census loss attributes, one rung at a time:

| rung (one change each) | n_train | truly-enriched share | vs oracle |
|---|---|---|---|
| 0 ORACLE substrate: live REGION nodes | 1488 | 0.1474 | 1.000 |
| 1 − AMBIG regions (circularity) | 1232 | 0.1198 | 0.813 |
| 2 − regions with no mass (+ struct_zero back) | 1150 | 0.1186 | 0.805 |
| 3 **+ BOUNDARY nodes = the production substrate** | 2277 | **0.0821** | **0.557** |
| 4 + the reliability WEIGHT | 2277 | **0.0612** | **0.415** |

and the factorial confirms the rungs very nearly compose — **0.802 (AMBIG) × 0.682 (boundaries) × 0.687
(weight) × 1.014 (counts) = 0.381 against a directly measured 0.399:**

| arm (vs the REGION reference) | recovery | | vs a MATCHED oracle (same nodes) |
|---|---|---|---|
| production | **0.399** | | **0.696** |
| + flat weights | 0.582 | | **1.014** ← perfect on its own substrate |
| + region-only (drop boundaries) | 0.547 | | 0.651 |
| + AMBIG admitted | 0.498 | | 0.699 |
| **flat + region-only + AMBIG** | **0.938** | | 0.894 |

⭐ **So the ~30 "unexplained" points were BOUNDARY DILUTION.** Boundaries are ~as numerous as regions
(1127 vs 1150) and only **5.1 %** of them are truly enriched against the regions' **12.1 %**, so adding them
nearly halves the enriched share of the census. It was never a missing-mode problem.

⭐⭐ **And the second control settles pass-0's role for good: on the production substrate with FLAT weights,
pass-0's own counts reproduce the matched oracle's enriched mass at 1.014.** Perfect counts buy nothing
because there is nothing left for them to buy.

**WHERE the mass lands** (the exact destination ledger; bands = the true depleted level / the valley / above
the split):

| population | n | weight | → depleted | → **the VALLEY** | → enriched |
|---|---|---|---|---|---|
| region × truly-ENRICHED | 139 | 0.0415 | 0.0003 | 0.0041 | **0.0372** |
| region × truly-depleted | 1011 | 0.5407 | 0.5030 | 0.0373 | 0.0005 |
| boundary × truly-ENRICHED | 58 | 0.0197 | 0.0000 | 0.0041 | 0.0155 |
| **boundary × truly-depleted** | 1069 | 0.3981 | 0.2641 | **0.1311** | 0.0028 |
| TOTAL | 2277 | 1.0000 | 0.7673 | **0.1766** | 0.0561 |

**Truly-depleted BOUNDARY nodes deliver 74 % of all the valley mass.** A boundary node's unspliced mass is a
mixture over its two flanks, so it lands *between* the two true modes — precisely in the gap the ground truth
leaves empty. That is why the FIT's own valley sits at **−0.04** while the oracle's is at **+0.72**.

⚠ **This also refutes the split as an innocent instrument (N1 lead 4).** `two_component`'s valley is stable
for the oracle (it is the truth's own gap) but the fit's is ~1.3 decades lower. **Scoring each landscape at
its OWN split reports 0.127 vs 0.140 = 91 % recovery and hides the entire defect.** Score at the oracle's
split; the §W3 convention was right.

### N1 — what the "217 split-crossers" are actually worth: ~13 % of the gap, not the lead

Reproduced (truly-enriched EXON nodes in the substrate, at the oracle's split): **stranded 9/752 exactly as
recorded**; unstranded 206/797 capture-ON, 223/1047 with VSTRONG (HANDOFF_17's 217/875 sits between the two —
the population definition it used was not written down, but the substance replicates).

| stratum | crossers | median shift | med eff | med G | med w | med var |
|---|---|---|---|---|---|---|
| unstranded, all truly-enriched | 1457 | −0.091 | 199 | 2126 | 0.493 | 0.130 |
| unstranded, **CROSSERS** | 355 (24.4 %) | **−0.293** | 241 | 2348 | 0.351 | 0.223 |
| stranded, all truly-enriched | 1098 | −0.012 | 233 | 3974 | 0.993 | 0.001 |
| stranded, **CROSSERS** | 38 (3.5 %) | −0.122 | 125 | 788 | 0.817 | 0.027 |

They are **not** the small-node probe-granularity population — their median eff-length (241) and count (2348)
are *higher* than the enriched population's. They are ordinary large exons carrying twice the declared
variance. **But they hold only 17.7 % of the truly-enriched weight (2.8 % stranded)**, so the entire
placement defect is worth ≈0.011 of the 0.084 gap. HANDOFF_17 promoted them as the best lead; they are a
real effect of secondary size.

### N1 — the WEIGHT's mechanism, and why the obvious fix is REFUTED

`w = ref/(v+ref)` does **not** discriminate enriched from depleted *within* a node class — on unstranded data
exons read 0.466 (enriched) vs 0.459 (depleted). It discriminates **classes**:

| unstranded | intergenic | intron | **exon** | boundary |
|---|---|---|---|---|
| median `w` | 1.000 | 0.954 | **0.459** | 0.746 |
| (stranded) | 1.000 | 0.988 | 0.906 | 0.957 |

and **every truly-enriched node is an exon**. So the reliability weight is effectively a node-class weight
that halves exactly the class carrying the enriched mode — an *informative* weight, not a variance
reduction. The declared precisions are honest; the mis-specification is that precision-weighting estimates
a common **mean**, whereas a population's mixing proportions are a **census**, and down-weighting removes a
node from the census.

**The implied fix — move precision from the MASS to the KERNEL WIDTH** (unit mass, kernel = Poisson ⊗
N(0, τ²), with τ = 0 on the zero-mass anchors, where `f_g·0 = 0` makes the density exact for every `f_g`;
measured, all 928 non-finite `var_gdna` in the substrate are precisely those nodes) — **is REFUTED, by its
own control:**

| arm | enr mass (n=13) | EMD/region (n=13) | EMD (ALL 32) | **EMD (ZERO-gDNA, n=9)** | enr mass on ZERO-gDNA |
|---|---|---|---|---|---|
| BASE (production recipe) | 0.0557 | 0.3045 | 0.2051 | 0.2059 | 0.0296 (6.2× oracle) |
| BND-OUT (regions only) | 0.0782 | 0.4356 | 0.2280 | **0.1190** | **0.0180 (3.8×)** |
| FLAT | 0.0816 | 0.2977 | 0.3535 | 0.7392 | 0.1605 (34×) |
| WIDTH (precision → kernel) | 0.0841 | 0.3061 | 0.3708 | 0.7777 | 0.1660 (35×) |
| WIDTH + BND-OUT | **0.1067** | **0.1704** | 0.2433 | 0.5491 | 0.1272 (27×) |
| **ctrl: CONST-τ + BND-OUT** | 0.1071 | 0.1777 | 0.2407 | 0.5354 | 0.1219 (26×) |

1. ⛔ **A single constant τ performs identically to the per-node one** (0.1777 vs 0.1704; 0.2407 vs 0.2433;
   0.5354 vs 0.5491). The per-node structure contributes nothing — **it is a global `h_pop` under another
   name**, which W2 already refuted.
2. ⛔ **It manufactures enrichment.** On zero-gDNA conditions the enriched mass goes 6.2× → 35× the oracle's.
   **Broadening kernels raises "enriched mass" for free**, so that metric is one-sided and must always be
   read against the zero-gDNA fabrication guard.
3. The `WIDTH + BND-OUT` win on the headline set (EMD 0.3045 → 0.1704) is therefore **not** a landscape
   improvement — it is two changes, one of which is a refuted smoother, scored on the one stratum that
   rewards it.

**What survives as a live candidate is BND-OUT alone**, and it is genuinely two-sided: better on `quick`
everywhere (0.3633 → 0.2454 headline, 0.2406 → 0.2048 all-16), much better on the ambig zero-gDNA guard
(0.2059 → 0.1190, false enrichment 6.2× → 3.8×), and **worse on ambig gDNA-bearing (0.3045 → 0.4356)** —
the same sign flip W1a already measured for boundaries. It is also what the shipped
`_fit_gdna_hyperprior` already does. **→ owner decision; only the W4 end-to-end test can settle it.**

### ⚖ THE CONFLICT THIS EXPOSES (owner decision needed)

The **circularity** rule and the **enriched-mass** goal are in direct opposition, and both are legitimate:

* AMBIG must be excluded — it is the two-root ambiguity the prior exists to resolve, and W1 measured that
  admitting it **loses** against the non-AMBIG population the prior is applied to (0.180 vs 0.175).
* AMBIG is **where a third of the enriched gDNA lives**, and excluding it costs 10 points of the very
  component the prior most needs. This is exactly the archived warning: *"a node reporting a genuine
  ENRICHED mode may be exactly the one an exclusion rule drops."*

**Candidate resolutions** (not yet tested): admit AMBIG at reduced weight; or exploit the architecture that
already exists — pass-0 → prior → **re-solve** produces better AMBIG estimates, so a *second* prior fit on
the re-solved AMBIG is non-circular by construction (the owner's open "iterative strategy"). The second is
strictly more principled and costs one extra fit.

---

### ⭐⭐⭐ N5 — THE LANDSCAPE'S RESOLUTION WAS A *MEASUREMENT* WIDTH, AND IT COLLAPSES AS ρ^(−1/2)

**Owner-spotted, 2026-07-27:** on a log-y axis the landscape is smooth through the depleted bulk and a comb
above `log10 ρ_g = 0`, worse the further right you go — **and the oracle does it too.** The diagnosis is the
owner's: a width set by a *linear-domain* (count) quantity, rendered on a log axis.
Figure `figures/gdna_n5_bandwidth.png`; scripts `scratchpad/gdna_n5_{bandwidth_diag,roughness,reference,plot}.py`.

**The arithmetic.** The per-node kernel is the Poisson likelihood of `g` at rate `ρE`, rendered on
`x = log₁₀ ρ`. Its curvature at the mode is

```
    d²/dx² [ g·(x·ln10 + ln E) − E·10^x ]  =  −g·(ln10)²      ⇒   sd_x = 1/(√g · ln10)  decades
```

**so the width shrinks as ρ^(−1/2)** — a node ten times denser carries ten times the counts and a √10-times
narrower kernel. Measured across the suite: **0.317 → 0.0078 decades, a 41× collapse**, dropping below
`GRID_H = 0.0290` (the axis's own resolution) above +1 decade:

| band of log₁₀ρ | (−5,−2] | (−2,−1] | (−1,0] | (0,+1] | (+1,+2.5] |
|---|---|---|---|---|---|
| median count `g` | 2 | 105 | 43 | ~175–1430 | 5100 |
| kernel sd (dec) | 0.317 | 0.143 | 0.137 | 0.041 | **0.0078** |
| **sd / GRID_H** | 10.9 | 4.9 | 4.7 | 1.4 | **0.27** |
| **fraction of nodes that are sub-grid DELTAS** | 0.00 | 0.21 | 0.11 | 0.64 | **0.88** |
| … carrying this share of the band's mass | 0.09 | 0.49 | 0.35 | 0.96 | **1.00** |

Below `GRID_H` a node is not a narrow density — **it is a delta in one cell at the wrong height**, and no
grid can render it otherwise. That is the comb, and it explains every part of the observation: smooth on the
left (low counts ⇒ wide kernels), spiky on the right, worse the higher you go, and "many modes" where the
truth has one.

**⭐ THE FIT WAS ALREADY PROTECTED. THE REFERENCE WAS NOT — AND THE REFERENCE IS THE SCORING INSTRUMENT.**
Roughness = total variation of `log P` per decade, height-free (a smooth unimodal bump scores 2–4):

| curve, band | (−5,−2] | (−2,−1] | (−1,0] | (0,+1] | **(+1,+2.5]** |
|---|---|---|---|---|---|
| **reference, PRODUCTION `oracle_landscape`** | 1.0 | 6.3 | 2.9 | 8.4 | **40.7** |
| fit, h = 0 (no resolution term) | 1.7 | 6.7 | 2.8 | 10.6 | **46.9** |
| **fit, kNN 0.5 (production)** | 1.8 | 6.2 | 0.8 | 1.6 | **5.1** |
| reference, CORRECTED (kNN 0.5) | 1.0 | 5.6 | 2.1 | 1.9 | **6.1** |

W2's kNN term already fixed the fit (46.9 → 5.1). It was never applied to the reference, which renders the
**TRUTH** with a *measurement* kernel — but `G` is not an observation of `ρ`, it **is** `ρ·E`, so there is no
posterior width to draw. **Every EMD-scored decision in W1/W2/W3 was taken against that comb**, and EMD to a
comb rewards a fit that is also a comb — the mechanical reason "EMD prefers h = 0".

*(And the owner's recollection that "the previous production version didn't have this" is correct:
`DensityNPMLE` uses `npmle_bandwidth = 0.15` **decades** — a fixed log-scale width, ≈5 grid cells. The comb
arrived with the per-node measurement kernel that replaced it.)*

**✅ FIX 1 — the reference gets a POPULATION resolution** (`oracle_landscape(s, knn_scale=…)`, the identical
kNN law the fit uses, so one resolution rule governs both). Validated against ground truth taken straight
off the oracle counts with no landscape involved — depleted **−2.244 dec (IQR 0.023, a point mass)**,
enriched **+1.030 dec, sd 0.307**:

| oracle rendering | dep loc | dep wid | enr loc | **enr wid** | enr mass | roughness >+1 |
|---|---|---|---|---|---|---|
| **PRODUCTION (measurement kernel)** | −2.136 | 0.902 | +1.076 | **0.169** ✗ 45 % too narrow | 0.1405 | **40.7** |
| **population kNN 0.5** | −2.102 | 0.897 | +1.307 | **0.328** ✅ vs truth 0.307 | 0.1399 | 6.1 |
| population kNN 1.0 | −1.858 | 1.224 | **+1.900** ✗ | 0.244 | **0.0383** ✗ | 3.6 |
| population kNN 2.0 / 4.0 | −1.99 | 1.08–1.15 | +1.9 ✗ | 0.25 | 0.08 ✗ | 2.5 / 1.7 |

kNN 0.5 is the **only** setting that reproduces the truth's enriched width *and* leaves the enriched mass
where it was (0.1399 vs 0.1405) *and* keeps the depleted component put. Past 0.5 the two modes start merging,
the `two_component` split runs away to +1.9 and the mass accounting breaks.

**✅ FIX 2 — the kernel width is floored at `GRID_H`.** Forced by the axis, not chosen: nothing narrower than
one cell is representable. Previously the floor was `knn_scale·GRID_H` = 0.0145 at scale 0.5, i.e. *below*
the grid, and the Poisson kernel itself was unfloored. Cost: EMD(old, new) = 0.0048 on the fit.

**✅ kNN 0.5 SURVIVES RE-SELECTION — now on a validated instrument.** Scored by SHAPE against the corrected
reference (substrate = regions only, per the owner's decision):

| fit kernel | dep loc | dep wid | enr loc | **enr wid** | enr mass | EMD | roughness >+1 |
|---|---|---|---|---|---|---|---|
| CORRECTED REFERENCE | −2.102 | 0.897 | +1.307 | **0.328** | 0.1399 | — | 6.1 |
| h = 0 | −2.357 | 0.530 | +1.169 | 0.284 | 0.0949 | 0.4450 | **46.9** |
| kNN 0.25 | −2.357 | 0.535 | +1.224 | 0.292 | 0.0951 | 0.4418 | 11.1 |
| **kNN 0.5** | −2.360 | 0.541 | +1.260 | **0.356** ✅ | 0.0961 | 0.4370 | **5.2** |
| kNN 1.0 | −2.366 | 0.564 | +1.332 | 0.415 | 0.0963 | 0.4288 | 3.0 |
| kNN 2.0 / 4.0 | −2.37 | 0.63–0.71 | +1.39 | 0.47 / 0.55 ✗ | 0.095 | 0.42 / 0.41 | 1.8 / 1.0 |

⚠ **EMD does not discriminate here** — it moves 0.4450 → 0.4084 monotonically across the whole sweep and
prefers the most smoothing at every reference, corrected or not. **Select the bandwidth by SHAPE**
(enriched width and mass against the reference, plus roughness); EMD is a summary, not a resolution
criterion.

⚠ **Enriched MASS is flat at ~0.095 across every kernel** — independent confirmation of N1: bandwidth cannot
touch the census, and with boundaries excluded the recovery is **0.096 / 0.140 = 0.69** (up from 0.40).

### ⭐ W1 RE-SCORED AGAINST THE CORRECTED REFERENCE — every substrate verdict SURVIVES

`scratchpad/gdna_w1_rescore.py`. Base is now REGIONS ONLY (owner decision), so `+ boundaries` appears as an
arm. ΔEMD vs base, scored against both instruments:

| arm | ambig ALL: Δ/comb | **Δ/CORRECTED** | quick ALL: Δ/comb | **Δ/CORRECTED** | enr mass Δ | rough >+1 |
|---|---|---|---|---|---|---|
| BASE (regions only) | — | — | — | — | — | 8.5 |
| − **struct_zero anchor** | +0.2611 | **+0.2560** | +0.6213 | **+0.6090** | +9.98 | 11.5 |
| **precision CUTOFF** | +0.1768 | **+0.1745** | +0.3689 | **+0.3593** | +7.68 | 10.1 |
| FLAT weights | +0.0074 | **+0.0011** | +0.1123 | **+0.0979** | +6.41 | 8.5 |
| + boundaries in | −0.0227 | −0.0251 | **+0.0415** | **+0.0381** | −0.011 | **12.9** |
| − gonly (intergenic) | −0.0337 | −0.0372 | −0.0105 | −0.0164 | +0.79 | 10.8 |
| + AMBIG in ⚠ | −0.0357 | −0.0382 | −0.0040 | −0.0034 | +1.32 | 6.5 |

**The corrected instrument changes the substrate ranking by ≤0.006 in every cell.** That is itself the
useful result: **the comb corrupted the BANDWIDTH choice — a shape question, where EMD is weak — but not
the SUBSTRATE ranking**, which is a mass-and-location question EMD handles fine. So W1 stands:

* the **zero-count structural anchor is critically load-bearing** (+0.26 / +0.61, and **+1.04 / +1.19 on the
  zero-gDNA guard**);
* a **hard precision cutoff is worse than ignoring precision entirely** (+0.175 vs +0.001 ambig; +0.359 vs
  +0.098 quick) — the general reason never to reintroduce one;
* the **weight is load-bearing where it matters**: with boundaries out, FLAT weights is nearly free on
  overall EMD (+0.001) but costs **+0.370 on the zero-gDNA guard** with 23× the reference's enriched mass.
  W1's headline "+0.1452" was inflated by the boundaries that are now excluded.
* **excluding boundaries is better supported than before** — it now wins `quick` (+0.038 to include),
  the zero-gDNA guard (+0.082), low-gDNA enriched mass, and roughness (8.5 vs 12.9).

⚠ **The `+ AMBIG in` row here is NOT the circularity test** — it is scored against an oracle that itself
contains AMBIG, which is the substrate-matching confound W1 identified. The circularity test is N2's Q2
table, scored against the NON-AMBIG population.

⚠ **`− gonly` now looks mildly better on EMD (−0.037 / −0.016)** where W1 called it neutral, but it is worse
on the zero-gDNA guard (+0.035) and on roughness (8.5 → 10.8). Left as-is; flagged, not decided.

### ⭐ THE FULL GAMUT — figures `gdna_n5_gamut_{ambig,quick}.png`. Two failure directions, both named

Every condition, fit vs corrected reference, log-y. Mass above a **fixed** `log₁₀ρ_g = 0` (no adaptive
split, so the accounting cannot move):

| gDNA level | capture OFF | capture ON | capture VSTRONG |
|---|---|---|---|
| **none** | **20×** (0.0021 vs 0.0001) | **42×** (0.0045) | **60×** (0.0064) |
| gdna1 | 0.63 | 0.63 | **0.35** |
| gdna5 | 2.2× | 0.76 | **0.24** |
| gdna100 | **6.5×** | 0.65 | **0.28** |
| gdna300 | 2.4× | 0.71 | — |

1. **FABRICATION wherever no enriched mode exists** — zero-gDNA and every capture-OFF condition. Relatively
   large (2–60×), absolutely small: **0.2–2.4 % of the landscape's mass** placed above 0 where the truth has
   ~0.01–1 %. ⚠ **This is the same surface as the additive `_gdna_arm`'s `gdna_none` regression (+0.0045)** —
   treat them as one risk, and note prior STRENGTH is the natural control.
2. **UNDER-RECOVERY, worst under the strongest capture** — capture ON is a level-independent **0.63–0.76**,
   but VSTRONG is **0.24–0.35**.

**And the VSTRONG shortfall is the RELIABILITY WEIGHT, measured, not a mystery:**

| stratum | reference | production | **flat weights** | oracle counts | flat + oracle counts |
|---|---|---|---|---|---|
| capture ON | 1.00 | 0.69 | **0.84** | 0.69 | 0.83 |
| **VSTRONG** | 1.00 | **0.29** | **1.24** | **0.28** | 0.91 |

Oracle counts buy **nothing** (0.28 vs 0.29) — N1's result again, now in the hardest stratum. Flat weights
alone take VSTRONG from 0.29 to over-shooting at 1.24. Every VSTRONG condition is unstranded, which is
exactly where N1 measured `w` to act as a node-class weight that halves exons. **So one unresolved
modelling question — how per-node precision should enter a MIXTURE estimate — holds the entire remaining
census reserve (0.69 → 0.84 on capture-ON, 0.29 → ~1.0 at VSTRONG).** The kernel-width answer to it was
refuted (§N5/N1); no replacement is derived. **Carry it into W4 as a named gate rather than blocking on it.**

### ⭐ N2 — THE ITERATIVE AMBIG RESOLUTION: premise CONFIRMED, verdict BLOCKED on W4 (2026-07-27)

`scratchpad/gdna_n2_{cache,iterative}.py`; the paired pass-0 / re-solved belief cache is
`scratchpad/gdna_refit_cache.pkl` (32 conditions, both `calib_refit_iters` settings, nothing else varied).

**Q1 — does the re-solve actually improve AMBIG?** Yes, substantially. Mean |log₁₀ρ̂ − log₁₀ρ_true| over
live AMBIG region nodes:

| stratum | n | pass-0 | re-solved | Δ | better |
|---|---|---|---|---|---|
| **all** | 6 604 | 0.4238 | **0.2454** | **−0.1784** | 20/32 |
| capture OFF | 3 078 | 0.4349 | **0.1661** | −0.2687 | **14/14** |
| stranded | 2 588 | 0.3382 | **0.0797** | −0.2584 | 9/12 |
| unstranded | 4 016 | 0.4752 | 0.3447 | −0.1304 | 11/20 |
| capture VSTRONG | 607 | 0.5613 | 0.6115 | **+0.0502** | 1/4 |

**The architecture's premise holds — a re-solve cuts AMBIG density error 42 %** — which is what the whole
`pass-0 → prior → re-solve` design predicted, now measured for the first time.

**Q2 — is prior #2 better?** No, and **the test is contaminated**:

| arm (13 gDNA-bearing capture-ON/VSTRONG) | enriched recovery | EMD vs ALL | **EMD vs NON-AMBIG** |
|---|---|---|---|
| #1 pass-0, AMBIG out (today) | 0.378 | 0.3045 | 0.2269 |
| pass-0, AMBIG in (the "circular" arm) | **0.528** | 0.2567 | **0.2186** |
| #2 re-solved, AMBIG out | 0.299 | 0.3998 | 0.2966 |
| #2 re-solved, AMBIG IN (the test) | 0.379 | 0.3654 | 0.2776 |

⚠ **`#2 re-solved, AMBIG OUT` is worse than `#1 pass-0, AMBIG OUT` on every column.** The re-solve degrades
the *non-AMBIG* nodes it touches — because **prior #1 is the shipped δ-pin `DensityNPMLE`, the object this
workstream exists to retire**, and this plan already records that it rescues `gdna_none` "for the wrong
reason". So the iterative route is scored through a prior known to be wrong-shaped. **N2 cannot be settled
before W4; re-run it the moment the landscape is wired.** The pairing is coherent, though: AMBIG nodes have
no own evidence, so almost any population prior helps them, while non-AMBIG nodes do, and a wrong prior
corrupts them.

### ⭐⭐⭐ N2 RE-RUN AFTER W4 — THE VERDICT IS IN, AND IT REVERSES (2026-07-27)

`scratchpad/gdna_n2b_iterative.py`. The blocked verdict above was scored through the δ-pin. Both caches were
**rebuilt at the current tree** (`gdna_cache_build.py`, `gdna_n2_cache.py`), so prior #1 is now the
`GdnaLandscape`. Three corrections over the first run, each one the record demands: prior #1 is the
landscape; the reference renders at **population resolution** (`knn_scale=0.5` — N5, the comb); and the
substrate is **regions only** (owner decision D1), matching the shipped `_fit_gdna_hyperprior`.

**Q1 — the premise is confirmed and STRONGER than it was through the δ-pin.** Mean
`|log₁₀ρ̂ − log₁₀ρ_true|` over live AMBIG region nodes:

| stratum | n | pass-0 | re-solved | Δ | better |
|---|---|---|---|---|---|
| **all** | 6 604 | 0.4014 | **0.2327** | **−0.1687 (−42 %)** | **30/32** (was 20/32) |
| capture OFF | 3 078 | 0.4149 | **0.2184** | −0.1965 | 14/14 |
| **capture ON** | 2 919 | 0.3503 | **0.1837** | −0.1666 | **13/14** (was contaminated) |
| stranded | 2 588 | 0.3326 | **0.1262** | −0.2063 | **12/12** |
| unstranded | 4 016 | 0.4427 | 0.2965 | −0.1462 | 18/20 |
| capture VSTRONG | 607 | 0.5330 | 0.4539 | −0.0791 | 3/4 |

**Q2 — and the iterative route now WINS, on the axis that decides.** Enriched recovery at the ORACLE's
split (n=13, gDNA-bearing capture-ON/VSTRONG); fabrication = fitted mass above `log₁₀ρ_g = 0` over the 9
ZERO-gDNA conditions, against the oracle's, where **any** enriched mass is false:

| arm | enr recovery | EMD vs NON-AMBIG | **fabrication** | n_train |
|---|---|---|---|---|
| **#1 pass-0, AMBIG out — what ships** | 0.572 | 0.2930 | **14.5×** | 1150 |
| pass-0, AMBIG in (the "circular" arm) | **0.773** | **0.2425** | ⛔ **26.4×** | 1369 |
| #2 re-solved, AMBIG out | 0.562 | 0.2945 | **4.7×** | 1150 |
| ⭐ **#2 re-solved, AMBIG IN (the test)** | **0.757** | **0.2544** | ⭐ **6.7×** | 1369 |

**The contamination is gone and the picture inverts.** Through the δ-pin, `#2 AMBIG out` was worse than
`#1` on every column (0.299 vs 0.378 recovery); with the landscape it is **flat** (0.562 vs 0.572) — the
degradation of the non-AMBIG nodes was the δ-pin's, exactly as predicted. And against what ships, the
iterative arm delivers **+0.185 enriched recovery, −0.039 EMD, and 2.2× BETTER specificity.**

⭐ **The circular arm is no longer the better trade.** It buys slightly more census (0.773 vs 0.757, a 2 %
edge) for **4× the fabrication** (26.4× vs 6.7×). The iterative route keeps 98 % of its census gain and
*improves* the guard over the shipped prior. That is the resolution §"THE CONFLICT" asked for, and it is
now measured rather than argued.

⚠ **The honest counter-reading.** Part of the specificity gain is generic shrinkage — the re-solve pulls
everything toward the prior, and enriched recovery on the zero-gDNA conditions falls 2.600 → 1.625 for the
AMBIG-out arm too. What is *not* generic is the combination: prior #2 keeps 98 % of the circular arm's
enriched census while halving its fabrication. Shrinkage alone would cost census in the same proportion.

⚠ Also unchanged: EMD vs NON-AMBIG still slightly favours the circular arm (0.2425 vs 0.2544). The
circularity argument is what rules it out, and it is principled, not empirical — W1 said EMD was never going
to settle it.

**→ OWNER DECISION.** Adopting this ships a **second** hyperprior fit (`pass-0 → prior #1 → re-solve →
prior #2 → re-solve`), i.e. one extra fit and one extra sweep per library. Not implemented; the measurement
is what was asked for.

### ⚠ N2 side-result — W1's empirical case against admitting AMBIG is STRATUM-DEPENDENT

W1 rejected admitting AMBIG on "it loses against the non-AMBIG population, 0.180 vs 0.175". **That
reproduces exactly on all 32 conditions (0.1816 vs 0.1752) — and the sign FLIPS on the 13 conditions where
an enriched mode exists at all (0.2186 vs 0.2269, AMBIG-in winning), along with enriched recovery
0.378 → 0.528.** The 19 conditions carrying the verdict are zero-gDNA and capture-OFF, which have no
enriched mode to lose. **The circularity argument is unaffected — it is principled and EMD was never going
to settle it (W1 said so) — but the empirical support cited alongside it does not hold where it matters.**

---

## 5. W4 — WIRING AND THE END-TO-END GATE

### ⭐ N4 step 1 LANDED IN THE WORKING TREE (uncommitted) — `_gdna_arm` is ADDITIVE (2026-07-27)

`simplex_logodds._gdna_arm`: `return ref + global_logprior` instead of `return global_logprior`. One line,
plus the docstring that records why the old justification was wrong:

> The reference is not an information claim to be superseded — it is the **measure** ψ is written against,
> and it is the only term bounding this arm at `f_g → 0`. `_rna_arm` bounds the other vertex and was never
> replaced, **so the two arms were not even treated alike.** Bayes composes a prior with a measure by
> ADDITION in log space; there is no double-count to avoid.

**A/B re-recorded from HEAD in-session** (`P0_BENCH_OUT=/tmp/pass0_oracle_bench_n4.tsv`), so the baseline is
not one of the stale ones HANDOFF_16 §0 warns about — base_r0 0.084116 / base_r1 0.067915, within 0.0005 of
the recorded suite figures.

| gate | result |
|---|---|
| **`refit=0` bit-identical 32/32** (prior is `None` at pass-0) | ✅ **32/32, mwae 0.084116 → 0.084116** |
| `refit=1` suite mwae | **0.067915 → 0.055724 (−0.0122, −18 %)**, 11 better / 6 worse / 15 flat |
| ruff | ✅ clean |
| pytest | 1190 pass, **12 golden failures** — goldens LAST, deliberately NOT regenerated |

| stratum | refit=1 base | additive | Δ | b/w |
|---|---|---|---|---|
| stranded ss_0.99 | 0.0289 | **0.0258** | −0.0031 | 4/0 |
| **unstranded × capON** (the historical regression site) | 0.1575 | **0.1271** | **−0.0304** | 5/2 |
| capture ON | 0.1076 | **0.0882** | −0.0193 | 9/2 |
| capture VERYSTRONG | 0.1300 | **0.0872** | −0.0428 | 2/2 |
| capture OFF | 0.0153 | 0.0166 | +0.0012 | 0/2 |
| **`gdna_none` (the FP guard)** | 0.0028 | 0.0073 | **+0.0045** | **0/3** |

⚠ **`gdna_none` regresses, and the mechanism is understood, not mysterious.** `_JEFFREYS_REF · log f_g`
→ −∞ as `f_g → 0`, so restoring it *penalises* the all-RNA vertex. Under REPLACEMENT the fitted logP was the
whole arm and its tallest (low-ρ) mode acted as a de-facto zero-gDNA prior — **which is precisely the "rescues
`gdna_none` for the wrong reason" this plan already flagged.** The additive form removes that accident and
charges the honest price. The gate as written ("every stratum improves or is flat") is therefore **not met**,
and this is an owner call: the improper-ψ fix is principled and worth −0.0122 suite-wide, against +0.0045 on
a guard that was previously passing by luck. **It is also measured through the OLD prior**, so both numbers
move again at W4.



1. `_gdna_arm` → `ref + global_lp` (**additive**). One edit; both solvers already consume `global_lp`.
   This also removes the crush mechanism (today's replacement deletes ψ's only `f_g → 0` barrier).
2. `_fit_gdna_hyperprior` returns `GdnaLandscape`; drop the dead `| gonly` and `additive`.
3. Prior strength swept; `config.gdna_prior_additive` / `npmle_bandwidth` retired if dead.
4. Development flag removed before ship (house rule: no ablation flags in production).

### ⭐⭐ W4 IS IMPLEMENTED (2026-07-27, uncommitted) — and it FAILS one named gate. Owner call.

**The code.** New module `calibration/gdna_landscape.py` (~260 lines incl. docs): a frozen `GdnaLandscape`
(`log_rho`, `logP`, `n_train`, `strength`) plus `fit_gdna_landscape`. Five private helpers, one job each —
`_grid`, `_poisson_kernels`, `knn_widths`, `_reliability`, `_render`. `calibrate._fit_gdna_hyperprior`
shrank to **substrate selection only** (the part that needs the chain) and its docstring is now the
four-axis taxonomy. `npmle.py` is **untouched** — retired in this role, still fitted for Role A / QC, and
all 30 of its unit tests still pass. **11 new unit tests**, one per property the design rests on.

**Config.** `gdna_prior_additive` **DELETED** (one caller, default `False`, selected a dead branch).
`gdna_prior_strength: float = 1.0` added — the temperature, validated `>= 0`.

**Two constants removed by derivation, not by fiat:**

1. ⭐ **`GRID = linspace(−5, 2.5, 260)` is gone.** The axis is now *exactly the domain `logprior` can be
   asked about*: ψ evaluates at `ρ_g = f_g·M/E` with `f_g ≤ 1`, so the top is `max_i log10(M_i/E_i)` (a hard
   bound — no node can be placed above its own total density) and the bottom is `min_i log10(1/E_i)`, the
   deepest one-count resolution wall. Both dominate every kernel centre, so nothing is truncated.
   ⚠ **A first attempt padded the span by the widest kernel and was measurably worse** (+0.007 ambig /
   +0.056 quick EMD, enriched width 0.42 → 0.54): the pad is set by *one* isolated node, which stretches the
   span 2–3 decades, and with a fixed point count the step coarsens and the whole landscape over-smooths.
   Recorded so it is not retried.
2. **The `1e-12` density floor is gone**, replaced by **one pseudo-node of complete ignorance** spread
   uniformly (ordinary Laplace smoothing; the "one" is one node, so there is no constant). It bounds the
   prior's log-range by `log W`, which is the *weak and correctable* property this object needs — it is
   fitted from biased pass-0 output. Behaviour-neutral: **+0.0001 suite**.

`_N_GRID = 260` and `_WIDTH_BINS = 12` remain, documented as **computational budgets** (finer is strictly
more faithful and strictly slower), not modelling choices. `_KNN_SCALE = 0.5` and `_S0` remain modelling
constants — `_S0` carrying the standing warning that it is a tuning constant in disguise and the biggest
lever left on the census.

**⭐ PRIOR STRENGTH IS MEASURED AND NEEDS NO TUNING — the default is optimal.** The sweep is **monotone on
every stratum**, `gdna_none` included:

| strength | ALL 32 | stranded | unstr × capON | capture ON | capture OFF | VSTRONG | **gdna_none** |
|---|---|---|---|---|---|---|---|
| 0.15 | 0.0767 | 0.0291 | 0.1545 | 0.1054 | 0.0313 | 0.1579 | 0.0762 |
| 0.30 | 0.0708 | 0.0278 | 0.1446 | 0.0991 | 0.0282 | 0.1398 | 0.0640 |
| 0.50 | 0.0649 | 0.0268 | 0.1349 | 0.0932 | 0.0247 | 0.1220 | 0.0510 |
| 0.75 | 0.0593 | 0.0262 | 0.1250 | 0.0874 | 0.0215 | 0.1063 | 0.0384 |
| **1.00** | **0.0562** | **0.0261** | **0.1193** | **0.0843** | **0.0197** | **0.0970** | **0.0315** |

**Tempering does not buy specificity — it makes every stratum worse, including the guard it was reserved
for.** So the knob ships at its default and adds no tuned constant. (Below ~0.15 the re-solve degenerates
toward pass-0, whose `gdna_none` is 0.0931.)

**THE GATES:**

| gate | result |
|---|---|
| **`refit=0` bit-identical 32/32** | ✅ **32/32** (prior is `None` at pass-0) |
| ruff / pytest | ✅ clean / **1193 pass**, 20 golden failures (goldens LAST, not regenerated) |
| `refit=1` suite mwae | **0.0679 → 0.0562 (−0.0117, −17 %)**, 9 better / 10 worse / 13 flat |
| **`unstranded × capON`** (historical regression site) | ✅ **0.1575 → 0.1193 (−0.0382)** |
| **`gdna_none`** (the FP guard) | ⛔ **0.0028 → 0.0315 (+0.0287), 0 better / 4 worse — GATE FAILED** |

| stratum | HEAD | + additive arm | **+ landscape** |
|---|---|---|---|
| ALL 32 | 0.0679 | 0.0557 | **0.0562** |
| stranded ss_0.99 | 0.0289 | 0.0258 | **0.0261** |
| unstranded × capON | 0.1575 | 0.1271 | **0.1193** |
| capture ON | 0.1076 | 0.0882 | **0.0843** |
| capture OFF | 0.0153 | 0.0166 | **0.0197** |
| capture VERYSTRONG | 0.1300 | 0.0872 | **0.0970** |
| **gdna_none** | **0.0028** | 0.0073 | **0.0315** |

⚠ **On suite mwae the landscape is a WASH against the additive arm alone (+0.0005).** It wins
`unstranded × capON` (−0.0080) and capture-ON (−0.0039) and loses VSTRONG (+0.0098) and `gdna_none`
(+0.0242). What it buys is not mwae: it is the object every N1/N5 measurement validated, it replaces a fit
whose good `gdna_none` behaviour this plan **pre-registered as an artifact**, and it removes two asserted
constants.

### ⭐⭐⭐ THE `gdna_none` DISSECTION — IT IS CIRCULARITY, AND THE ADDRESS IS UPSTREAM (2026-07-27)

Scripts `scratchpad/gdna_d{1,2,3,4}_*.py`. On a zero-gDNA library **every fragment called gDNA is a false
positive**, which makes this the one place the error is exactly identifiable.

**⭐ The prior is FAITHFUL, not inventive — the arithmetic closes to within 0.7 nats.** On
`none_ss_0.50_nrna_none_capture_on`, pass-0 places **84 training nodes above `log10 ρ_g = 0`, carrying
1.54 % of the training weight**, in a library with no gDNA at all. Spread over the ~50 grid cells of the
enriched range against a peak cell of ~0.05 that is `log(0.0154/50 / 0.05) = −5.1` nats; **the landscape
measures −4.4.** It is reporting what pass-0 told it.

**The δ-pin's −25 nats is not a better estimate — it is a narrow fixed kernel (`h = 0.15` dec) plus EM
competition collapsing those nodes into the dominant low mode.** It wins here by *distrusting its own
training data*, which is an accident of bandwidth, not a principle that can be ported. And the trade is
structural, the same mechanism with opposite sign in the two regimes:

> **an ADDITIVE KDE preserves a minority; an EM MIXTURE competes it away.** Preserving the minority is
> exactly why the landscape recovers the enriched mode (N1: the δ-pin's failure) and exactly why it
> preserves a false-positive tail. One property, two signs.

**Both priors help; the δ-pin helps more.** False-positive gDNA mass over the 9 zero-gDNA conditions:

| | pass-0 | δ-pin re-solve | **landscape re-solve** |
|---|---|---|---|
| FP fragments | 1,857,227 | **56,030 (0.03×)** | **629,497 (0.34×)** |

**⭐ AND THE UPSTREAM DEFECT HAS A PRECISE ADDRESS.** pass-0's own false-positive rate on zero-gDNA
libraries, by stratum:

| stratum | class | nodes | total mass | FP mass | **FP %** |
|---|---|---|---|---|---|
| **unstranded × nrna_none** | **boundary** | 142 | 98,683 | 50,926 | **51.6 %** |
| **unstranded × nrna_none** | **exon** | 917 | 3,220,458 | 987,324 | **30.7 %** |
| unstranded × nrna | boundary | 2,348 | 590,027 | 94,002 | 15.9 % |
| unstranded × nrna | exon | 1,242 | 6,508,737 | 624,016 | 9.6 % |
| stranded × * | exon / boundary | — | — | — | 1.1–2.2 % |
| any | **intron** | 2,194 | 1,794,101 | 767 | **0.04 %** |

**Introns are essentially perfect (0.04 %) — the intron factory doing its job — while unstranded exons and
boundaries with no nascent RNA are 31 % and 52 % wrong.** That is the count-zero-information wall in its
purest form: unstranded ⇒ flat strand likelihood, `nrna_none` ⇒ empty introns ⇒ the neighbours have nothing
to impute with either. **This is a pass-0 defect with a named population, and it is worth its own work.**

### ⭐⭐⭐⭐ THE SOURCE OF THE FALSE gDNA, TRACED TO THE NODE (2026-07-27)

Owner: *"If there's zero DNA there's no [gDNA]. There's no seeding at the ends and the boundaries. Where is
this DNA coming from?"* — **it is seeded by the solver, at the boundaries, by initialization.**
Scripts `scratchpad/gdna_d{5,6}_*.py`, on `gdna_none_ss_0.50_nrna_none_capture_off` (1,642,048 unspliced
fragments, **zero** true gDNA, **481,734 called gDNA = 29.3 % of the library**).

**The chain, each link measured:**

**1 — 84 % of the false gDNA sits on exons that DID receive an RNA imputation message.** Not on
evidence-starved nodes. The strand channel alone puts them at `f_g = 0.506` (the honest uninformative
answer); the message pulls them only to `0.271`.

| stratum | nodes | mass | FP | FP % | share of FP | `f_g` strand-only |
|---|---|---|---|---|---|---|
| NO evidence at all | 97 | 127,408 | 66,197 | 52.0 % | 13.7 % | 0.520 |
| **RNA message only** | **433** | **1,494,745** | **404,985** | **27.1 %** | **84.1 %** | 0.506 |
| has own spliced mass | 56 | 19,895 | 10,551 | 53.0 % | 2.2 % | 0.530 |

*(Replay fidelity verified per node first: 516/586 nodes exact to 1e-6, carrying 91.4 % of the mass and
**90.0 % of the FP** — the tool's global "max 2.2e-01" warning is a handful of nodes, not the attribution.)*

**2 — the message claims `f_rna = 0.646` where the truth is exactly `1.000`.** Whatever it does not claim
is left to gDNA by construction, and the solve lands at 0.271, just under the message's own 0.354 floor.

**3 — the claim is a CONSTANT, not a measurement.** Flat across precision quintiles
(0.635 / 0.536 / 0.658 / 0.661 / 0.662), flat across a 23× eff-length range (0.599 → 0.659), flat across the
terminus/junction split (0.639 vs 0.656 — so it is **not** the missing-TSS/TES bit, and **not** `ω_graft`'s
lower-bound-as-equality, both of which were the leading hypotheses and are hereby **refuted**). Its
distinct values are two:

| claimed `f_rna` | nodes | `1/claim − 1` |
|---|---|---|
| **0.662338** | 202 | 0.509803 |
| **0.671052** | 146 | 0.490197 |

— 348 of 433 nodes sharing two values to six decimals, summing to exactly 1.000000. A node-independent
constant cannot be a per-node measurement.

**4 — the message carries the destination's OWN scale and an IMPORTED composition.** `mo_p = log(cp·E_r/M)`,
and measured `cp / (M/E_r) = 0.6623` at p25, p50 **and** p75, with `E_r/E_g = 1.0000` (so it is not a
units or effective-length error). The relay pin re-derives the scale from the destination's own mass; only
the **composition** comes from the source. The exon is therefore told "you are 66 % RNA" by its neighbour,
and it cannot be more RNA than its neighbour claims.

**⭐ 5 — THE SOURCE. 863 boundary nodes feed those messages, and 93 % of them have ZERO unspliced mass and
sit at the UNSOLVED DEFAULT `f_g = 1.0000` — 100 % gDNA — while carrying real spliced junction flux
(median 36.3 fragments).** `solvable = (fp|fn) & (mass_unspliced > 0)`, so a boundary with no unspliced mass
is never solved and keeps its initialization, which is *100 % gDNA by construction*
(`node_init`, the "unsolved default"). On a library with no gDNA that is exactly backwards — **and the
spliced flux is positive proof that RNA crosses that seam.** All 81 live boundaries read `f_g = 0.510`
regardless of whether they carry spliced flux (56 with, 25 without — identical medians), i.e. **the spliced
evidence is not reaching their own composition either.**

> **The false gDNA is not leaking in from anywhere. It is seeded at the boundaries by the initialization,
> and then delivered into the exons as an imported composition.** The user's instinct that "there is no
> seeding at the ends and the boundaries" is right about the *biology* and precisely wrong about the *code*.

**What this makes newly plausible (untested, and the obvious next experiment):** a node with zero unspliced
mass has no gDNA either — its density is `0` for every `f_g`, which is exactly the argument already used for
the landscape's zero-count *anchor*. Defaulting it to 100 % gDNA is an assertion the data contradicts.
Candidates, in order of how little they assume: (a) make the unsolved default **inert** for emission rather
than all-gDNA; (b) let a boundary with spliced flux and no unspliced mass declare itself RNA-bearing, which
is what the flux measures. ⚠ Both change pass-0 for every library, so each needs the full 32-condition A/B —
and note `gdna_none` is only 9 of the 32, so a fix must not be judged there alone.

### ⭐⭐⭐⭐⭐ THE BUG, FOUND: `_pin_v` FEEDS A NODE ITS OWN GUESS BACK AS AN INCOMING MESSAGE

`bp_solver._pin_v` (line ~488) rescales an incoming claim so `Σ_c ρ_c·E_c = M`, and **for a component the
message does NOT supply it substitutes the DESTINATION'S OWN density into the budget:**

```python
sg = np.where(pg_ > 0.0, g, og)      # og = the destination's own gDNA density
sp = np.where(pp_ > 0.0, p, op)
sn = np.where(pn_ > 0.0, n, on)
k  = M / (sg * E_g + (sp + sn) * E_r)
return g * k, p * k, n * k
```

On the false-positive exons the message supplies **RNA only** — measured: the gDNA-message precision is 0 on
**433/433** of them. With `E_g == E_r == E` and an incoming RNA density equal to the node's own total
density `ρ_tot = M/E`, the pin returns

```
    n_pinned = n·M / (og·E + n·E) = ρ_tot² / (ρ_tot + og) = ρ_tot · 1/(1 + f_g_own)
```

**so the delivered RNA fraction is exactly `1/(1 + f_g_own)` — a function of the destination's OWN
strand-only self-solve, containing no information from the source at all.** Verified to **2.1e-16** on the
median FP exon and exactly on 252/433. The measured claims are the four discrete unstranded strand-only
values and nothing else:

| own strand-only `f_g` | 0.4577 | 0.4902 | 0.5098 | 0.5423 |
|---|---|---|---|---|
| **delivered `f_rna` = 1/(1+f_g)** | 0.6484 | **0.6623** | **0.6711** | 0.6860 |

**⭐ THE CONTROL, and it is decisive.** Nothing distinguishes these libraries but strandedness:

| library | own `f_g` | **budget the pin reserves for gDNA the message never claimed** | **false-positive rate** |
|---|---|---|---|
| zero-gDNA, **unstranded**, capture OFF | 0.5065 | **33.6 %** | **29.3 %** |
| zero-gDNA, **unstranded**, capture ON | 0.5073 | **33.6 %** | **33.2 %** |
| zero-gDNA, **stranded**, capture ON | **0.0130** | **1.2 %** | **1.4 %** |

**The pin's reservation IS the false-positive rate.** On unstranded data the strand channel has no
information, so the node's own self-solve sits at the uninformative `f_g ≈ ½`; the pin then reserves half
the mass budget for that imaginary gDNA; the RNA message is shrunk by `1/(1+½)`; and the solve reads back
~34 % gDNA — **the number it started from.** It is self-fulfilling, and it is circular in the strongest
sense: the node's own guess re-enters as the message it is supposed to be *receiving*.

**⚠ The docstring's intent is sound in the other direction** — "a message carrying gDNA only still delivers
`f_g < 1`" — and it records that rescaling all three blindly **regresses capture-OFF 3.6×**. So the fix is
NOT "substitute zero".

**⭐ The principled fix direction (untested):** the substitution is unconditional on whether the node has any
evidence for `og`. On these nodes `tau_own = 0` — *no composition evidence whatsoever* — so `og` is a pure
guess, and letting a pure guess claim mass budget is what closes the loop. **Gate the substitution on the
own belief actually having composition evidence** (`τ_own > 0`); where it does not, the unsupplied component
cannot defend a share. That keeps the partial-claim semantics exactly where they were derived (a node that
knows something about its own gDNA) and removes them where they are fabricating. ⚠ Full 32-condition A/B
required — the same pin serves every library, and `gdna_none` is 9 of 32.

**Two things RULED OUT on the way, with measurements, so they are not re-tried:**

* ⛔ **The 100 %-gDNA initialization default is INERT — the owner's invariant HOLDS.** Flipping the default
  from all-gDNA to all-RNA on every zero-mass node changed **zero** other nodes and left the false-positive
  mass **bit-identical** (481,734 → 481,734). It is misleading to read but it is not the source.
* ⛔ **The emitting boundaries are 96.4 % INTRON↔EXON, not exon↔exon** (3.6 %). The owner's reasoning was
  right — an exon↔intron boundary carries no *unspliced* fragments in a zero-gDNA, no-nascent library — but
  it does carry *spliced* junction flux, and that is what emits.

### ⚠ AND THE EVIDENCE *DOES* SEPARATE THEM — so this is not an identifiability wall

Nodes pass-0 places above `ρ_g = 0`, real (gDNA-bearing capture-ON) vs false (zero-gDNA), unstranded both:

| population | n | med `w` | med `var` | med count | overlap |
|---|---|---|---|---|---|
| **REAL** enriched | 1,418 | **0.439** | **0.172** | 2,500 | — |
| **FALSE** enriched | 312 | **0.045** | **2.554** | 1,576 | `w` **29 %**, `var` **28 %**, count 80 % |

**pass-0's declared variance already knows: the false nodes carry 15× the variance and 10× less weight, and
the distributions overlap only ~28 %.** So the reliability weight is working as designed — the residual
1.5 % of weight is simply enough to make a −5-nat plateau, and closing to the δ-pin's −25 would need another
factor of `e²⁰`, which no re-weighting delivers. ⚠ Note this cuts against a naive `_S0` fix: the enriched
census wants `_S0` LARGE (less down-weighting) and the FP guard wants it SMALL. **Directly opposed, so
`_S0` cannot serve both** — another reason the fix is upstream.

### ⛔ Causes RULED OUT by measurement — do not re-probe them

⭐ **One real bug was found and fixed on the way: `knn_widths` was not computing the k-th-nearest-neighbour
distance.** It used `max(a_i − a_{i−k}, a_{i+k} − a_i)` — the FAR edge of a 2k window — which for a node
with no near neighbours on one side reaches back into the bulk, **handing the widest kernel in the fit to
the most isolated node**. Now the exact k-th-NN distance via a vectorised bisection on the V-shaped
objective (O(n log k)), brute-force verified in a unit test. Measured effect on the pathological tail: p99
width **1.04 dec → 0.21–0.46**. Suite effect: **neutral (0.0562 → 0.0565)** — kept as a correctness fix,
not as a remedy.

| ruled out | measurement |
|---|---|
| the landscape's **location** | median −5.2 dec on zero-gDNA; fitting on **ORACLE counts gives −5.2 too** |
| the **dynamic range** | landscape 7–20 nats vs δ-pin 11–24 — the δ-pin is the harsher |
| the **tail mass** | above −1: landscape 0.018–0.024 vs δ-pin 0.096–0.103 — landscape is MORE conservative |
| the **density floor** | uniform / none / interior-scoped are **identical to 0.5 nats** at ρ=+1 — the kernels dominate, so the owner-recalled `floor_eps=0.02` scoping is inert here (it is still the better design in principle: absence of data outside the support IS information) |
| the **kNN width form** | exact k-th-NN: p99 1.04 → 0.21–0.46 dec, suite **+0.0003** |
| the **prior strength** | monotone; 1.0 optimal on every stratum including this one |

### ⛔ THE `gdna_none` REGRESSION — five causes RULED OUT, so the next session need not re-probe them

All the damage is on **unstranded** conditions (`none_ss_0.50_nrna_none_capture_off` 0.008 → 0.130,
`…_capture_on` 0.015 → 0.117, `…_present_capture_verystrong` 0.004 → 0.076). **Every stranded `gdna_none`
condition is fine** (0.0003 → 0.0012 and similar). Ruled out, each by measurement:

1. **NOT the landscape's location.** Its median density on zero-gDNA conditions is **−5.2 decades**, and
   fitting the same landscape on **ORACLE counts gives −5.2 as well** — the object is faithful.
2. **NOT the dynamic range.** Landscape 7–20 nats vs the δ-pin's 11–24 — comparable, and the δ-pin is the
   harsher of the two.
3. **NOT the tail.** Prior mass above `log10 ρ_g = −1`: landscape **0.018–0.024** vs δ-pin **0.096–0.103**
   on the worst conditions — **the landscape is markedly MORE conservative** where the wrong nodes land.
4. **NOT the density floor** (Laplace vs `1e-12`: +0.0005 on this stratum).
5. **NOT the prior strength** (monotone; tempering makes it worse).

What remains is the prior's shape *along each node's ray* `ρ_g = f_g·M/E` and its interaction with ψ —
which needs the node-dissect loop (`pass0_node_dissect.py`), not another summary probe. **⚠ Note the
stranded/unstranded split points at pass-0's own unstranded over-call as the upstream cause (median pass-0
`f_g` on those conditions is 0.27–0.38 in a library with NO gDNA), which would make this circularity rather
than an estimator defect** — but that is a hypothesis, not a measurement.

**Gates.** `refit=0` **bit-identical 32/32** (the prior is `None` at pass-0 — a wiring correctness check).
`refit=1`: every stratum improves or is flat, with **`unstranded × capON`** (historical regression site) and
**`gdna_none`** (the FP guard) called out explicitly — note production currently rescues `gdna_none` *for
the wrong reason*, so a correct landscape must earn it honestly. Fit-substrate accuracy improves; held-fixed
`z2` does not regress. ⚠ **Re-record the baseline from a `git stash` of HEAD in the same session** —
HANDOFF_16 §0: a not-32/32 HEAD-vs-baseline means the *baseline* is broken.

---

## 5b. ⭐⭐ REAL DATA — FIRST RUN EVER (2026-07-27). It works, and it surfaces two things the toy cannot

`scratchpad/real_data_check.py`, on the four cached cfRNA/capture payloads (human index v7). Neither the
landscape nor the pin fix had ever been run outside the synthetic suite. There is no oracle here, so this is
not an accuracy test — it is every check that can be made without one.

**✅ It runs, at genome scale, at both refit settings, with no NaN and no degenerate output.**
33–105 s (118 k – 702 k nodes), consistent with `PERFORMANCE.md`.

| sample | nodes | mass-wt `f_g` r0 | **r1** | `ρ_g` global | κ | `od_g` | `od_r` |
|---|---|---|---|---|---|---|---|
| LBX0190 | 118 059 | 0.2355 | **0.1240** | 4e−6 | 0.0023 | **0.2000** | **0.2000** |
| LBX0588 | 256 540 | 0.8033 | 0.7897 | 3.1e−4 | 0.0120 | 0.0031 | **0.2000** |
| MO_3021 | 481 655 | 0.2655 | **0.2024** | 3.7e−5 | 0.0020 | **0.2000** | **0.2000** |
| vcap | 702 218 | 0.4257 | 0.4067 | 1.2e−3 | 0.0001 | 0.0923 | **0.2000** |

⭐ **The one accuracy-adjacent signal is encouraging.** LBX0190 is on record at **~15 % gDNA** (memory
`cfrna_sample_characteristics`). Pass-0 says 23.6 %; the re-solve with the landscape says **12.4 %** — the
prior moves it from 9 points too high to 3 points low, and **49 % of the mass moves by more than 0.05**. On
the sample whose truth we half-know, the hyperprior moves the answer the right way and by a lot.

### ⛔⛔ FINDING 1 — THE STRAND OVERDISPERSION SATURATES ITS CEILING ON EVERY REAL LIBRARY

`od = 0.2` is not a default, it is the **Beta(2,2) CEILING** — "the most overdispersion allowed"
(`gdna_strand.py`). The control settles that it is genuinely fitted:

| | `od_g` | `od_r` |
|---|---|---|
| synthetic suite (5 conditions spanning the battery) | 0.0003 – 0.1335 | **0.0008 – 0.0017** |
| **real data** | **0.2000 on 2 of 4** | ⛔ **0.2000 on 4 of 4** |

**The toy fits 120–250× below the ceiling; every real sample is pinned at it.** This is the class of defect
memory `synthetic_suite_is_poisson_omega_zero` says the suite **cannot** show — the simulator is multinomial
at fixed abundance, so ω = 0 by construction. **No amount of work on the 32-condition battery would have
found it.**

#### ⭐⭐ THE CORRECT READING (owner, 2026-07-28) — THE CEILING IS A GUARD, AND SATURATION IS A CANARY

An earlier draft of this section framed saturation as "real libraries are more strand-overdispersed than the
model allows, so the ceiling may need raising". **That is wrong, and the error is a biological one.**

> **gDNA is symmetric around 50/50 by construction.** It is genomic DNA; there is no strand preference to
> disperse. `od_g` has a *biological* mean of ½ that is not in question.

So an apparent excess strand variance on a node we believe is pure gDNA is **not** evidence about gDNA. It
is evidence that **the node is not pure gDNA** — i.e. that transcription is occurring where the annotated
reference has no record of it (antisense transcripts, unannotated genes, readthrough). Note *where* the
seeds come from (`gdna_strand.py`): *"a count-observable region or boundary side — intergenic, intronic,
exon–intron / exon–intergenic seam"* — **precisely the places unannotated transcription hides.**

**And the estimator has no defence against it.** It is a pooled ratio of sums,

```
    od_mom = Σ_s excess_var_s / Σ_s gdna_var_s
```

so a handful of contaminated seed nodes with large excess variance dominate the numerator and drag the
pooled estimate to the ceiling. There is no trimming, no outlier component, no per-seed influence bound.

**Three consequences, and they invert the earlier recommendation:**

1. ⛔ **Do NOT raise the ceiling.** Doing so imports the annotation's incompleteness directly into the
   strand likelihood — **the strongest and only intrinsic gDNA/RNA evidence we have**
   (`CALIBRATION_ARCHITECTURE.md`). Diluting it to accommodate contamination is a self-inflicted wound.
2. ✅ **The ceiling is currently doing the right job.** It is the only thing preventing annotation error
   from propagating into the strand model, and saturation is the **canary** that says the seed set is
   contaminated on this library. Keep it, and *report* saturation as a QC signal.
3. ⭐ **The real work is a ROBUST estimator.** The overdispersion measure must be robust to contamination
   by unannotated transcription — that is a property of the estimator, not a constant to re-tune. The
   pooled ratio-of-sums is the thing to replace.

**The human transcriptome is more complex than the annotation records, and the model must respect that
rather than force the disagreement onto the strand model.**

⚠ **Status: the saturation is MEASURED; the contamination mechanism is owner-reasoned and NOT yet
measured.** The falsifiable test is cheap and should be run before any estimator change: **is
`Σ excess_var_s` dominated by a few seed nodes?** If a small share of seeds carries most of the numerator,
non-robustness is confirmed and the fix direction follows. If the excess is spread evenly, the diagnosis is
wrong and something else is going on.

### ⚠ FINDING 2 — `ω_graft` SPANS 15× ACROSS FOUR REAL SAMPLES

| LBX0190 | LBX0588 | MO_3021 | vcap |
|---|---|---|---|
| 0.25 / 0.40 | 0.69 / 0.80 | **3.83 / 3.77** | 2.80 / 3.41 |

The ROADMAP records `ω_graft` as "10× apart on two real samples" and "expected to be fragile on real data".
On four samples it is **15× apart (0.25 → 3.83)**, and the two strands agree closely *within* each sample —
so it is a real per-library quantity, not noise, and one library-wide scalar is standing in for something
that varies by more than an order of magnitude between libraries. **P1D_P1E_DEBTS.md's brief is confirmed,
not weakened.**

### ⚠ The landscape's shape on real data — reported, NOT yet judged

`logP` span is **5.1–5.4 nats** on all four (the toy gives 7–20), the low mode sits at **−4.5 to −5.4
decades**, and mass above `log₁₀ρ_g = 0` is **0.0009–0.0223** (vcap, a capture library, is the highest — the
right ordering). ⚠ The bimodality question is **NOT answered here**: this script's split-finder is a crude
interior-minimum, not `two_component`, so its "0.1 % enriched mass" on all four samples is not evidence of a
collapsed mode. Rendering the real-data landscapes properly (and against the held-out predictive criterion,
with its known under-smoothing bias) is the next real-data step.

## 6. The constants ledger (every constant, and its status)

| constant | where | status |
|---|---|---|
| `max(g, 1.0)` — one-count resolution wall | weight, kernel centres, metrics | ✅ **derived** — a count of 1 is the resolution limit |
| `S0 = (0.15·ln10)²` | weight `ref = 1/max(g,1) + S0` | ⚠ a reference *scale*, but 0.15 is asserted → W2 |
| `tau_prior` clip lo `0.15` | boundary gate | ✅ **DELETED 2026-07-27** — removal was free (−0.0034 / −0.0255) |
| `tau_prior` clip hi `0.70` | boundary gate | ✅ **DELETED** — never bound on any condition |
| `tau_prior` fallback `0.3` | boundary gate | ✅ **DELETED** |
| `hup=0.70, hdn=0.02, cap_up=1.0` | asymmetric projection | ⛔ magic → retired (W3) |
| `h_proj=0.15` | symmetric projection | ⛔ magic → replaced by per-node σ (W3) |
| `GRID = linspace(−5, 2.5, 260)` | the log-ρ axis | ⚠ asserted range/resolution → derive from data support |
| `h_pop` | new | ✅ **REFUTED as a global constant** (W2); replaced by kNN spacing, `scale` 0.5 provisional |
| `α` (NB) | new, optional | 🔷 physical; **not identifiable on this suite** |
| prior strength | new | 🔷 default 1; optimized, owner-sanctioned |
| **`ν` (Student-t)** | the projection | ✅ **FITTED, not picked** (N3): held-out predictive ⇒ **ν = 3.0** on both suites (oracle MLE 2.0–2.2). ⛔ the kurtosis MoM saturates at 4 and must not be used |
| ~~`c` (σ inflation)~~ | the projection | ✅ **DELETED** — the heavier-tailed likelihood replaces it; a bare 5× harms 9/32 |

**Target: every ⛔ removed and every ⚠ resolved or explicitly accepted before ship.**

---

## 7. Sequence

```
W0 cleanup ─▶ W1a fix the two scoring confounds ─▶ W1 substrate ─▶ W2 kernel ─▶ W3 projection ─▶ W4 wire
   (done)      (oracle_pointmass, matched oracle)     (no assumptions)   (empirical h)   (posterior σ)   (A/B)
```

W1a first: both confounds bias the very comparisons W1 and W2 are decided on. W1 before W2 because the
substrate sets `n_train`, which is what `h_pop` responds to. W3 reads the landscape, so it follows. The
residual directional bias (the old "G3") stays **deferred** — optimise the fit and the projection first,
then look at what is left.

## 8. Do not repeat

| | |
|---|---|
| a single-constant precision cutoff for admission | ⛔ this plan removes one; do not add another |
| an unconditional upward projection lift | ⛔ measured negative on 19/32 |
| one global smoothing constant for the kernel | ⛔ fixes `none_*` (6 modes → 2) and destroys a real mode (EMD 0.131 → 0.226) |
| exclusion rules assumed safe | ⛔ `solve_gate_design.md`: derived, implemented, **refuted** (+0.010/+0.025) |
| validating anything overdispersion-dependent on this suite | ⛔ Poisson by construction (ω = 0) |
