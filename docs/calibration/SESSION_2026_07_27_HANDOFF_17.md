# SESSION 2026-07-27 — HANDOFF 17 (LIVE): the gDNA hyperprior rebuild

**Supersedes `SESSION_2026_07_27_HANDOFF_16.md`** for hyperprior work (16 remains the record of the
optimization/cleanup session and its §0 baseline warning still binds). This session was **Phase 2 / Role B**:
retire the δ-pin `DensityNPMLE` and build a production landscape + projection.

⚠ **This line described the FIRST half of 2026-07-27 and is now superseded** — three changes have since
landed in `src/` (uncommitted). See the END-OF-SESSION STATE block immediately below; the rest of this file
is the record of the exploration that preceded them.

---

> ## ⭐⭐⭐⭐ END-OF-SESSION STATE (2026-07-27, later session) — READ THIS FIRST
>
> **Suite mwae: refit=0 `0.0841 → 0.0788` (pass-0 itself, −6.3 %); refit=1 `0.0679 → 0.0525` (−22.7 %).**
> 11 b / 5 w / 16 f at r0; 9 b / 12 w / 11 f at r1. Ruff clean, **1193 tests pass, 21 goldens pending
> (LAST — do not regenerate until the owner accepts the tree).** **Nothing is committed.**
>
> Three things landed in `src/` (all uncommitted):
> 1. **`simplex_logodds._gdna_arm` is ADDITIVE** (`ref + global_lp`) — owner-approved.
> 2. **W4: `calibration/gdna_landscape.py`** (new) + `_fit_gdna_hyperprior` reduced to substrate selection;
>    `gdna_prior_additive` deleted, `gdna_prior_strength` added; 12 new unit tests.
> 3. **`bp_solver._pin_v` no longer rewrites the delivered message** — a BP violation on the MODE.
>    `docs/calibration/pin_derivation.md` is the derivation; it is the largest single win of the session.
>
> One thing was implemented and **REVERTED**: the share transfer (`pin_derivation.md` (★)). ⛔ Do not
> re-derive it — it is measured and refuted, and the refutation is inline in `bp_solver.py`.
>
> **THE OPEN ITEM** is the P-2 residual: capture-OFF **+0.0009** at refit=1, concentrated in
> **unstranded × capture-OFF × gDNA-bearing**. See `pin_derivation.md` §11.

> ## ⭐⭐⭐⭐⭐ NEXT SESSION (2026-07-27, later) — WHAT CHANGED. READ THIS BEFORE THE BLOCK ABOVE.
>
> **No `src/` change.** The tree is **bit-identical 32/32** to the state described above (verified in
> session against a freshly re-recorded baseline: **r0 0.078786, r1 0.052470**, mass-weighted over all 32
> conditions). Three things happened:
>
> **1. THE P-2 RESIDUAL IS DIAGNOSED AND CLOSED AS "NOT FIXABLE AT THE PIN" — `pin_derivation.md` §12.**
> The mechanism on file (§10: a partial RNA-only claim under-calling `f_g`) is **wrong in every
> particular** — 0.2 % of the harmed mass is partial, 96–99 % of the regression is **over**-calling, and the
> over-claim is on **gDNA**. Two repair arms were run and both fail (the BP-legal restore recovers 11 % and
> costs more than it recovers; a derived conditional-mean shrinkage is inert). ⭐ What the pin was masking is
> worth **20×** the residual: the **reframe applied to the gDNA component** makes the delivered gDNA level
> **5.4× too large** at capture-OFF, and an `r_g ≡ 1` sizing arm takes the residual's own stratum
> **0.0495 → 0.0146 at r0 / 0.0313 → 0.0077 at r1 (6/6 better)** — while costing unstranded × capON
> **+0.17/+0.21**, so it bounds the prize rather than being a candidate. ⚠ This **re-opens ROADMAP §11**,
> whose refutation was measured while the pin was cancelling `r`.
>
> **2. N2 IS RE-RUN AGAINST THE LANDSCAPE AND THE VERDICT REVERSES** — plan §"N2 RE-RUN AFTER W4". Both
> caches were rebuilt at this tree. Q1 strengthens (AMBIG error −42 %, **30/32**, capture-ON now 13/14); Q2
> is no longer contaminated, and the iterative arm now **beats what ships on every axis**: enriched recovery
> 0.572 → **0.757**, EMD-vs-non-AMBIG 0.2930 → 0.2544, fabrication **14.5× → 6.7×**. The "circular" arm buys
> 2 % more census for **4×** the fabrication. → **owner decision**: adopting it ships a second prior fit.
>
> **3. THE 21 GOLDENS ARE REGENERATED — 1214 pass, 0 failures.** Ruff clean. Also fixed one real breakage
> found on the way: `scripts/debug/gdna_cache_build.py` still called `cfg.calibration.gdna_prior_additive`,
> which W4 deleted, so the substrate cache could not be rebuilt at all.
>
> **STILL NOT COMMITTED.** New tools: `scratchpad/p2r_{a..f}_*.py` (the residual dissection),
> `scratchpad/gdna_n2b_iterative.py`, `scratchpad/p2res_report.py`, `scratchpad/real_data_check.py`,
> `scratchpad/P2_RESIDUAL_NOTES.md`.

> ## ⭐⭐⭐ UPDATE 2026-07-27 (later session) — §0 IS ANSWERED. READ THIS BEFORE §0.
>
> **N1 is DONE and the enriched-mass deficit is closed exactly.** It was never a missing mode, a bandwidth
> problem, or a pass-0 accuracy problem — **it is a SUBSTRATE CENSUS effect**, and the ~30 unexplained
> points were **BOUNDARY DILUTION**. Full record + figure: `gdna_hyperprior_production_plan.md` §"N1 — THE
> DEFICIT IS NOW CLOSED", figure `figures/gdna_n1_massflow.png`.
>
> ```
>     recovery 0.399 = AMBIG 0.802 × BOUNDARIES 0.682 × WEIGHT 0.687 × COUNTS 1.014   (= 0.381)
> ```
>
> * **PLACEMENT is fine (0.873)** — truly-enriched region nodes put 90 % of their kernel above the split.
> * **The 217 split-crossers are worth ~13 % of the gap, not the lead.** They hold 17.7 % of the
>   truly-enriched weight (2.8 % stranded), and they are *large* exons (median eff 241, G 2348), not the
>   small-node probe-granularity population.
> * **Truly-depleted BOUNDARY nodes deliver 74 % of all the valley mass** between the two true modes — which
>   is why the fit's own valley sits at −0.04 against the oracle's +0.72.
> * **On the production substrate with FLAT weights, pass-0's own counts match a matched oracle at 1.014.**
>   Perfect counts buy nothing because there is nothing left to buy.
> * ⛔ **Scoring each landscape at its OWN split reports 91 % recovery and hides the whole defect** — §0's
>   convention (score at the ORACLE's split) was right; do not "fix" it.
> * ⛔ **The obvious fix (precision → kernel WIDTH instead of mass) is REFUTED by its own control**: a single
>   constant τ performs identically, i.e. it is a global `h_pop` under another name (W2-refuted), and it
>   takes false enrichment on zero-gDNA conditions from 6.2× to 35× the oracle.
>
> **N3 is DONE** — ν is fitted, not picked: the kurtosis MoM ⛔ saturates at its ν > 4 boundary on every
> stratum; the deployable oracle-free held-out estimator gives **ν = 3.0**, keeping 81 % of hand-picked
> ν = 2's gain. **N2's premise is CONFIRMED** (a re-solve cuts AMBIG density error 42 %) but its **verdict is
> blocked on W4** — prior #1 is the retiring δ-pin and it degrades the non-AMBIG nodes. **N4 step 1 is in the
> working tree, uncommitted**: `_gdna_arm` additive, refit=0 bit-identical 32/32, refit=1 **0.0679 → 0.0557**.
>
> ### ⭐⭐⭐ AND THE BANDWIDTH BUG THE OWNER SPOTTED IN THE FIGURE — N5, now fixed
>
> The landscape is smooth below `log10 ρ_g = 0` and a comb above it, worse the higher you go, **and the
> oracle does it too.** Cause, exactly as the owner diagnosed (a linear-domain width on a log axis): the
> per-node kernel is the Poisson likelihood, whose width on the log axis is `1/(√g·ln10)`, so it
> **shrinks as ρ^(−1/2)** — a **41× collapse** (0.317 → 0.0078 dec) that goes **below `GRID_H`** above
> +1 decade, where **88 % of nodes carrying 100 % of the band's mass become deltas in a single cell.**
>
> **The FIT was already protected by W2's kNN term** (roughness above +1: **46.9 → 5.1**). **The REFERENCE
> was not — and the reference is the scoring instrument** (roughness **40.7**). `oracle_landscape` renders
> the TRUTH with a *measurement* kernel, but `G` is not an observation of `ρ`, it IS `ρ·E`. **Every
> EMD-scored decision in W1/W2/W3 was taken against that comb**, and EMD to a comb rewards a comb.
>
> Fixed: `oracle_landscape(knn_scale=…)` renders at POPULATION resolution — validated against ground truth
> (enriched sd **0.328 vs the truth's 0.307**; production gave **0.169**) — and the kernel width is now
> floored at `GRID_H`, which the axis forces. **kNN 0.5 survives re-selection on the corrected instrument.**
> ⚠ EMD does not discriminate bandwidth (monotone in smoothing at every reference) — **select by SHAPE.**
> ✅ **W1 has since been re-scored against the corrected reference and every substrate verdict SURVIVES**
> (the ranking moves ≤0.006): the comb corrupted the BANDWIDTH choice, a shape question where EMD is weak,
> not the SUBSTRATE ranking, a mass question it handles.
> Full record: `gdna_hyperprior_production_plan.md` §N5; figure `figures/gdna_n5_bandwidth.png`.
>
> **TWO OWNER DECISIONS WERE OPEN AND ARE NOW TAKEN (owner, 2026-07-27): boundaries are EXCLUDED from the
> training substrate; the additive `_gdna_arm` is KEPT.** §4.0 records the evidence behind each.

## 0. ⚠ THE CENTRAL OPEN PROBLEM — WE DO NOT RECOVER THE ENRICHED MODE  *(ANSWERED — see the banner above)*

**Ground truth (owner, and verified): there are exactly TWO modes.** gDNA is uniform; hybrid capture
partitions nodes into captured (enriched) and not-captured (depleted). Measured on oracle counts:

| node class | median true log₁₀ρ_g | **IQR** |
|---|---|---|
| intergenic | −1.28 | **0.02** |
| intron | −1.28 | **0.04** |
| exon (captured) | +1.41 | — |

The depleted level is a **point mass**. The enriched level is likewise tight — **sd 0.091** once you isolate
nodes of adequate size (eff-length ≈ 442).

**Our fitted landscape recovers only 40 % of the enriched component's mass** (0.056 vs a true 0.140).
Everything else about the fit is good: the depleted component is reproduced almost exactly (loc −2.148 vs
−2.136, width 0.944 vs 0.902), and the enriched component's *location* and *width* are close. **It is the
MASS that is missing, and that is the single most important open item.**

### What the deficit is NOT

* **Not bandwidth.** Ruled out arithmetically: the Poisson kernel width is `~1/√G`, and enriched nodes carry
  `G` ≈ 3,400–9,100 ⇒ kernel sd **0.0045–0.0075 dec** against an observed spread of **0.23–0.32 dec**
  (**30–70×** larger). No kernel choice can create or remove it.
* **Not pass-0's accuracy.** Substituting **ORACLE counts** on the same substrate gives 0.081 — versus 0.082
  for our own counts with the same (flat) weighting. **Perfect counts buy nothing.**

### What it IS, so far (≈30 of 60 points still unexplained)

| cause | worth | evidence |
|---|---|---|
| the **reliability weight** | ~18 pts (40 % → 58 %) | bites **only unstranded**: enriched mean `w` 0.528 vs depleted 0.721 when unstranded, **0.930 vs 0.930** when stranded |
| the **substrate** (AMBIG excluded) | ~10 pts (40 % → 50 %) | **507/507 unstranded and 461/461 stranded** excluded truly-enriched nodes are AMBIG — zero from any other cause |
| **UNEXPLAINED** | **~30 pts** | ⭐ **start here** |

**The most promising lead on the unexplained remainder:** of the 875 truly-enriched unstranded nodes that
*are* in the substrate, pass-0 places **217 (25 %) on the DEPLETED side of the split** (median shift −0.112).
Each one moves its entire weight from the enriched component to the depleted one. On stranded data only
9/752 do this. **That is the next thing to interrogate.**

### The conflict this exposes (owner decision pending)

Circularity and enriched-mass recovery are in direct opposition, and both are legitimate:
* AMBIG **must** be excluded — it is the two-root ambiguity the prior exists to resolve, and admitting it
  measurably **loses** against the non-AMBIG population the prior is applied to (EMD 0.180 vs 0.175).
* AMBIG is **where a third of the enriched gDNA lives**; excluding it costs 10 points of exactly the
  component the prior most needs. This is the archived warning made concrete: *"a node reporting a genuine
  ENRICHED mode may be exactly the one an exclusion rule drops."*

**Recommended resolution (untested):** use the architecture that already exists. `pass-0 → prior → re-solve`
yields better AMBIG estimates, so a **second** prior fit on the re-solved AMBIG is non-circular by
construction. Costs one extra fit; strictly more principled than admitting AMBIG at a fudged weight.

---

## 1. Read these, in this order

1. `ROADMAP.md` — entry point.
2. **`gdna_hyperprior_production_plan.md`** — ⭐ **THE working plan.** Workstreams W0–W4, the constants
   ledger, the do-not-repeat table, and every measurement from this session.
3. `gdna_hyperprior_resurrection.md` — how the track was resurrected and the per-scenario review (§4b).
   ⚠ Its §2.3 "a perfect landscape is worth +0.03" is **SUPERSEDED** (banner in place).
4. **THIS FILE** §0 and §4.
5. `SESSION_2026_07_27_HANDOFF_16.md` §0 — every stored A/B baseline is stale; re-record from a `git stash`
   of HEAD in the same session.

Never read `docs/calibration/archive/` for guidance (provenance only).

---

## 2. What was DONE this session

**W0 — cleanup.** Deleted **16 superseded debug scripts (2,225 lines)**; the `scripts/debug/` hyperprior
layer went 30 → 13 files, per `scripts/README.md`'s delete-don't-archive policy. Promoted the landscape
recipe out of a bare `def recipe(s)` loaded by `exec(open(...).read())` into importable
`gdna_explore_lib.{recipe, recipe_substrate, recipe_weights}` — **verified byte-identical over all 48
scenarios**. `DensityNPMLE` is **retired, not deleted** (owner: it may resurface; delete only at full
production).

**Caches regenerated end-to-end at HEAD** — all 32 ambig + 16 quick scenarios re-run from the BAMs through
scan → pass-0 → the production prior fit. Every fit returned a valid object.

**W1a — two scoring confounds fixed.**
* **The smear is NOT a problem: do not deconvolve.** Where gDNA exists (`G > 0`) the Poisson-smoothed oracle
  and the raw population are the same distribution (EMD **0.033**, inflation 1.04). All disagreement is on
  zero-count nodes, where the *pointmass* reference is the wrong instrument (it invents a location at
  `log10(1/E)`; the Poisson form correctly spreads downward). **Use `oracle_pointmass` only on `G > 0`.**
* **Boundaries** are a measurably different population — true density **+0.28…+0.42 dec above** regions and
  3–7× tighter. Cost of dropping them flips sign between suites (ambig +0.0232, quick −0.0453).
  **Owner: keep boundaries INCLUDED for now; ablate later.**

**W1 — substrate, by taxonomy.** The four axes that must not be conflated: **circularity** → structural
exclusion (AMBIG); **identifiability** → structural inclusion (the zero-count anchor); **precision** →
continuous weight, never a cutoff; **geometry** → matched reference (boundaries cross, regions contain).

* ✅ **The `tau_prior` boundary gate is DELETED** — 3 bare constants (0.15/0.70/0.3) and a cliff, clipping to
  its lower bound on 19/32 conditions, discarding boundaries carrying 12–16 % of boundary weight. Removal
  was **free**: EMD −0.0034 (ambig) / −0.0255 (quick).
* ⭐ **A hard precision CUTOFF is worse than ignoring precision entirely** (+0.249 vs +0.145 ambig; +0.430 vs
  +0.136 quick). That is the general reason never to reintroduce one.
* The **zero-count structural anchor is critically load-bearing** (dropping it costs +0.32 / +0.60 EMD).

**W2 — the kernel.** A **global** `h_pop` is refuted: the oracle's depleted mode is itself sharp and our fit
is already *under*-peaked there (0.141 vs oracle 0.193), so global smoothing flattens the good region 7–9×
to fix the sparse tail. **Abramson** sample-point adaptivity is ~neutral (a spike *is* high pilot density at
its own location, so it keeps itself narrow). ✅ **Adopted: adaptive resolution by NEAREST-NEIGHBOUR
SPACING**, `h_i = scale · dist(a_i, k-th nearest neighbour)`, `k = √n`, `scale = 0.5` (provisional — ambig
prefers 0.5 by both criteria, quick prefers 2.0 on EMD). It fixes the overfitting **by construction** (fewer
nodes ⇒ farther neighbours ⇒ wider kernels): mode-count-vs-sample-size slope **+5.6 → −1.0**.

**W3 — the projection.** ✅ **Posterior form adopted** (σ = the node's own belief width; no new constant):
gain +0.208 mean / +0.177 median, negative on 2/32, vs the asymmetric form's median **−0.073** and negative
on **19/32**. ⭐ **A Student-t likelihood beats inflating σ**: t(ν=2) gives +0.527 with 4/32 harmed, versus
Gaussian c=5's +0.749 with **9/32** harmed. A t *is* a Gaussian marginalised over an uncertain width, and ν
is derivable from the residual's kurtosis (`6/(ν−4)`) by method of moments — **fit it, do not pick it.**

---

## 3. ⛔ MISTAKES MADE THIS SESSION — DO NOT REPEAT

Recorded because each was tempting and several produced numbers that reached other documents.

| mistake | why it was wrong | corrected to |
|---|---|---|
| Scoring landscape structure by **counting modes** at `prominence=0.02` | counts noise wiggles: 7 modes at 0.02, **1 at 0.35**, on a ground truth of two | `gdna_explore_lib.two_component` — depleted loc/width, enriched loc/width/**mass** |
| "**minimum mode gap 0.12 dec**" ⇒ a σ ceiling | those were gaps between wiggles; the real separation is depleted↔enriched ≈ **2.7 dec**, so σ = 0.35 is safe | owner's "light-years apart" was right |
| "**a perfect landscape is worth only +0.03**" (reached ROADMAP) | measured with the *inert* symmetric projection — a projection that cannot move cannot benefit from a better landscape | with a working projection an ORACLE landscape **nearly doubles** gain (0.527 → 0.748) and takes harmed conditions 4/32 → **0/32**. **Fit quality is first-order.** |
| **Adjacent-node pairs** as a dispersion estimator | premise fails: adjacent exons differ by **0.916 dec in TRUE density** (different capture) | use the oracle residual directly on the toy |
| **`sd(z)` = 3–33** ⇒ "uniform 5× under-declaration" | outlier-driven | quintiles of declared σ vs realised \|error\| give ratios **0.6–1.7** — calibrated in the MEDIAN, with a heavy TAIL |
| Boundary A/B that dropped the **weights** along with the boundaries | conflates two changes | `recipe(s, sel=…, w=…)` — vary one at a time |
| Citing **total counts** where the question was about **densities** | 3,360 vs 9,060 is 0.48 dec and says nothing about density structure | counts are valid only for the comb test |

**Also settled and not to be re-litigated:** a **global** smoothing constant (refuted); an **unconditional
upward** projection lift (negative on 19/32); a **single-constant precision cutoff** for admission; and
validating anything **overdispersion-dependent** on this suite (Poisson by construction, ω = 0).

---

## 4. ▶ NEXT STEPS  *(N1 and N3 DONE; N2 blocked; N4 step 1 in the tree. The list below is superseded by
##    §4.0, kept because its framing of each question still stands.)*

### ⭐ 4.0 — WHERE IT ACTUALLY STANDS (2026-07-27, later session)

**TWO OWNER DECISIONS — BOTH NOW TAKEN (owner, 2026-07-27):**
**D1 → EXCLUDE boundaries. D2 → KEEP the additive `_gdna_arm`.** The evidence behind each is below.

**D1 — do BOUNDARY nodes train the prior?** ✅ **EXCLUDED.** The exploration recipe included them ("owner:
keep INCLUDED for now; ablate later"); the shipped `_fit_gdna_hyperprior` excludes them, so this realigns
the two. It is the **single largest lever on the enriched census (×0.682)**, and the evidence was two-sided:

| dropping boundaries | ambig gDNA-bearing (n=13) | ambig ALL (32) | **ambig ZERO-gDNA (9)** | quick headline (4) | quick ALL (16) |
|---|---|---|---|---|---|
| EMD vs region oracle | 0.3045 → **0.4356** ✗ | 0.2051 → 0.2280 ✗ | 0.2059 → **0.1190** ✓ | 0.3633 → **0.2454** ✓ | 0.2406 → **0.2048** ✓ |
| enriched recovery | 0.40 → **0.56** ✓ | — | false enrich 6.2× → **3.8×** ✓ | 0.59 → **1.09** ✓ | 0.24 → **0.45** ✓ |

Better on `quick` everywhere, better on the fabrication guard, worse on ambig gDNA-bearing — the same sign
flip W1a found. ⚠ **The ambig column was scored against the COMB reference (§N5); it has not been re-scored
against the corrected one.** Taken on the grounds that it wins the specificity guard and both `quick`
strata and is the only arm that raises the enriched census without a refuted smoother. With boundaries out
the enriched recovery is **0.69** (from 0.40).

**D2 — keep the additive `_gdna_arm`?** ✅ **KEPT.** Suite refit=1 **0.0679 → 0.0557 (−18 %)**,
`unstranded × capON` −0.0304, refit=0 bit-identical 32/32 — against **`gdna_none` +0.0045 (0 better /
3 worse)**, which the stated gate did not allow. The regression is *understood*: restoring `½·log f_g`
removes the accidental zero-gDNA prior the replacement form was providing, i.e. exactly the "rescues
`gdna_none` for the wrong reason" already on record. ⚠ **12 goldens still need regeneration — do that LAST,
after W4, so the solver's output settles once rather than twice.**

### ⭐⭐ W4 IS IMPLEMENTED (2026-07-27, uncommitted) — ONE GATE FAILS, AND IT IS AN OWNER CALL

New module `calibration/gdna_landscape.py`; `_fit_gdna_hyperprior` is substrate selection only; `npmle.py`
untouched (retired in this role, its 30 tests still pass); `gdna_prior_additive` deleted,
`gdna_prior_strength` added. **11 new unit tests. Two asserted constants removed by derivation** — the grid
is now exactly `logprior`'s own domain, and the density floor is one pseudo-node of ignorance.

| gate | result |
|---|---|
| `refit=0` bit-identical 32/32 | ✅ |
| ruff / pytest | ✅ / 1193 pass, 20 goldens pending (LAST) |
| `refit=1` suite | **0.0679 → 0.0562 (−17 %)**, 9 b / 10 w / 13 flat |
| `unstranded × capON` | ✅ **0.1575 → 0.1193** |
| **`gdna_none`** | ⛔ **0.0028 → 0.0315 — FAILED** |

**Prior strength needs no tuning: the sweep is monotone and 1.0 is optimal on every stratum**, `gdna_none`
included — tempering buys no specificity. So the knob ships at its default and adds no tuned constant.

**On the `gdna_none` failure, five causes are ruled out by measurement** (location — the landscape's median
is −5.2 and an ORACLE-count fit agrees; dynamic range; tail mass — the landscape is *more* conservative than
the δ-pin there; the density floor; the strength). All the damage is **unstranded**; every stranded
`gdna_none` condition is fine. Plan §"THE `gdna_none` REGRESSION" has the numbers, so do not re-probe them.
⚠ The stranded/unstranded split points at pass-0's own unstranded over-call (median `f_g` 0.27–0.38 in a
library with NO gDNA) — i.e. plausibly circularity rather than an estimator defect, but that is untested.

**⚖ THE DECISION:** ship a −17 % suite improvement (and −24 % on the historical regression site) against an
11× regression on the false-positive guard — noting this plan **pre-registered** that HEAD's `gdna_none`
rescue was "for the wrong reason", and that on suite mwae the landscape is a wash against the additive arm
alone (+0.0005): what it buys is the validated object and two fewer constants, not points.

### ⭐⭐⭐ THE `gdna_none` DISSECTION IS DONE — IT IS CIRCULARITY, AND THE ADDRESS IS UPSTREAM

Plan §"THE `gdna_none` DISSECTION". Headlines:

* **The prior is FAITHFUL, not inventive.** pass-0 places 84 training nodes / **1.54 % of the training
  weight** above `log10 ρ_g = 0` on a library with NO gDNA; that predicts −5.1 nats and the landscape
  measures **−4.4**. The δ-pin's −25 is a narrow fixed kernel plus EM competition **distrusting its own
  training data** — an accident of bandwidth, not a portable principle.
* **The trade is structural:** an ADDITIVE KDE preserves a minority, an EM MIXTURE competes it away. That is
  *why* the landscape recovers the enriched mode and *why* it preserves a false-positive tail. One property,
  two signs.
* **Both priors help** — FP fragments 1,857,227 (pass-0) → 56,030 (δ-pin) → 629,497 (landscape).
* ⭐ **pass-0's over-call has a precise address: unstranded × `nrna_none`, boundaries 51.6 % and exons
  30.7 % of their mass called gDNA — while INTRONS are 0.04 %** (the intron factory working). Stranded is
  1.1–2.2 %. The count-zero-information wall in pure form: flat strand likelihood, and with no nascent RNA
  the introns are empty so the neighbours cannot impute either. **This is a pass-0 defect worth its own work.**
* ⚠ **The evidence DOES separate real from false enrichment** (med `w` 0.439 vs 0.045, med `var` 0.172 vs
  2.554, ~28 % overlap) — so this is not an identifiability wall, and the weight is working. Closing to the
  δ-pin's −25 would need another `e²⁰`, which no re-weighting delivers. ⚠ `_S0` cannot serve both ends: the
  enriched census wants it LARGE, the FP guard SMALL.
* ⭐ **One real bug found and fixed en route:** `knn_widths` was computing the far edge of a 2k window, not
  the k-th-nearest-neighbour distance — **handing the widest kernel in the fit to the most isolated node**.
  Now exact (vectorised bisection, brute-force verified). p99 width 1.04 → 0.21–0.46 dec; suite neutral.

### ⭐⭐⭐⭐ AND THE FALSE gDNA IS NOW TRACED TO THE NODE — pass-0 SEEDS IT AT THE BOUNDARIES

Plan §"THE SOURCE OF THE FALSE gDNA". On a library with **zero** gDNA, 29.3 % of the mass is called gDNA.
The chain, each link measured (replay fidelity verified per node first — 90 % of the FP is on exactly-
reproduced nodes):

1. **84 % of the false gDNA is on exons that DID get an RNA message.** Strand-alone puts them at 0.506; the
   message pulls them only to 0.271.
2. **The message claims `f_rna = 0.646` where the truth is 1.000** — the remainder is gDNA by construction.
3. **The claim is a CONSTANT**: flat in precision, in eff-length (23× range), and in the terminus/junction
   split — which **REFUTES** both leading hypotheses (`ω_graft`'s lower-bound-as-equality, and the missing
   TSS/TES). Its distinct values are **two**: 0.662338 (×202) and 0.671052 (×146), summing to exactly 1.
4. **The message carries the destination's OWN scale and an IMPORTED composition** — `cp/(M/E_r) = 0.6623`
   at p25/p50/p75, `E_r/E_g = 1.0000`. The pin re-derives scale; only composition travels.
5. ⭐ **THE SOURCE: 863 boundaries feed those messages and 93 % have ZERO unspliced mass sitting at the
   UNSOLVED DEFAULT `f_g = 1.0` — 100 % gDNA — while carrying real spliced junction flux (median 36.3).**
   `solvable = (fp|fn) & (mass_unspliced > 0)`, so they are never solved and keep an initialization that
   asserts pure gDNA. All 81 live boundaries read `f_g = 0.510` whether or not they carry spliced flux, so
   **the spliced evidence is not reaching their composition either.**

> **The false gDNA is not leaking in. It is SEEDED at the boundaries by the initialization and then
> delivered into the exons as an imported composition.**

**Then, in order:**
0. ⭐⭐ **PASS-0: the unsolved-default seeding at boundaries.** A node with zero unspliced mass has no gDNA
   either — its density is 0 for every `f_g`, the same argument already used for the landscape's zero-count
   anchor. Candidates, least-assuming first: (a) make the unsolved default **inert for emission** rather
   than all-gDNA; (b) let a boundary with spliced flux and no unspliced mass declare itself RNA-bearing.
   ⚠ Both change pass-0 for every library ⇒ full 32-condition A/B; `gdna_none` is 9 of 32, so do not judge
   a fix there alone.
1. ✅ **DONE — W1 re-scored against the corrected reference.** Every substrate verdict survives; the
   instrument changes the ranking by ≤0.006. The comb corrupted the BANDWIDTH choice (a shape question),
   not the SUBSTRATE ranking (a mass question). Plan §"W1 RE-SCORED".
2. ⭐ **W4 — GO.** `_fit_gdna_hyperprior` returns the landscape (regions-only substrate, kNN 0.5 kernel with
   the `GRID_H` floor) and the Student-t ν = 3 projection; `_gdna_arm` is already additive; add the
   prior-strength temperature, default 1, and sweep it. **Named gates, beyond the standing ones:**
   * **`gdna_none` / capture-OFF fabrication.** The landscape places 0.2–2.4 % of its mass above
     `log₁₀ρ_g = 0` where the truth has ~0.01–1 % (2–60× relative). Same surface as the additive arm's
     `gdna_none` +0.0045. **Prior strength is the control; this is what it is for.**
   * **VSTRONG.** Enriched-census recovery is 0.29 there vs 0.63–0.76 at capture ON — and it is **entirely
     the reliability weight** (flat weights → 1.24; oracle counts → 0.28, i.e. nothing).
3. **Re-run N2** once prior #1 is the landscape rather than the δ-pin. Premise confirmed; verdict contaminated.
4. **Regenerate the 12 goldens**, then **real data** (`cfrna:LBX0190`, `gdna_benchmark_5mb`) — never yet run.

**⚠ THE ONE KNOWN RESERVE, carried not closed:** the reliability weight holds the whole remaining enriched
census (0.69 → 0.84 at capture ON, 0.29 → ~1.0 at VSTRONG). It is a single measured mechanism — on
unstranded data `w` is a node-class weight that halves exons, and every enriched node is an exon — but the
obvious fix (precision → kernel width) is refuted, and no replacement is derived. **W4's end-to-end number
is what should price it.**

### N1 ⭐ — INTERROGATE THE MISSING ENRICHED MASS (the priority)  *(DONE — see the banner at the top)*

~30 of the 60 missing points are unexplained. Leads, most promising first:

1. **The 217/875 split-crossers.** Truly-enriched unstranded nodes that pass-0 places on the depleted side
   (median shift −0.112, vs 9/752 stranded). *Who are they?* Break down by node class, eff-length, mass,
   spliced evidence, `τ_own`, and the `EV_*` evidence class. Is it the small-node probe-granularity
   population (E ≈ 53, sd 0.531), or genuinely mis-solved large exons?
2. **Interaction of weight × substrate.** The 18-point and 10-point effects were measured separately;
   measure them jointly (flat weights **and** AMBIG admitted) — they may not be additive.
3. **Where does the mass actually land?** Track the *destination* of every truly-enriched node's weight in
   the fitted landscape. Mass is conserved, so if it is not in the enriched component it is somewhere
   specific — find it rather than infer it.
4. **Sanity-check the split itself.** `two_component` picks the deepest valley; verify that is stable and
   is not itself moving the accounting.

### N2 — the AMBIG conflict (§0)

Test the **iterative** resolution: fit prior #1 without AMBIG → re-solve → fit prior #2 including the
re-solved AMBIG → compare enriched mass and the non-AMBIG EMD. Non-circular by construction.

### N3 — W3b: fit ν by method of moments

Excess kurtosis `6/(ν−4)` on the residual `z`; confirm it lands near the empirical optimum (≈2–4); drop the
tuned `c` entirely. **Do not ship a bare 5×.**

### N4 — W4: wire into ψ

`_gdna_arm` → **`ref + global_lp` (additive)** — today the fitted prior *replaces* the `½·log f_g` reference,
which is the documented node-1055 crush mechanism. Keep the **full `log P` along the ray** (no
projection-and-anchor collapse, no `w_anchor` — a Gaussian anchor would destroy the multi-modality). Add a
**prior strength** temperature (default 1) and optimise it; the metric is *"does real data still dictate the
solution"*.

**Gates:** `refit=0` **bit-identical 32/32** (prior is `None` at pass-0); `refit=1` every stratum improves
or is flat, with **unstranded × capON** and **gdna_none** called out (production currently rescues
`gdna_none` *for the wrong reason* — a wrong-shaped landscape whose tallest mode happens to be low); fit-
substrate accuracy improves; held-fixed `z2` does not regress.

### N5 — real data

`cfrna:LBX0190` and `gdna_benchmark_5mb` — **neither has ever been run for this work.** ⚠ Held-out
predictive likelihood is the only oracle-free criterion available there, and it is **biased toward
under-smoothing** (it scores the landscape as a *predictive* object, which wants population ⊛ noise, while
its use as a *prior* wants the population). Carry the toy-derived `scale` across; do not re-select on real
data.

---

## 5. RUN BOOK

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1 && cd /Users/mkiyer/proj/rigel
```

| script | what it does |
|---|---|
| `scripts/debug/gdna_cache_build.py --suite <dir> --out scratchpad/<x>.pkl` | regenerate the substrate cache (~35 min from cold BAMs; **only needed after a pass-0 change**) |
| `scratchpad/gdna_fit_review.py --suite ambig --out FIG.png` | per-scenario fit vs oracle + the coverage table |
| `scratchpad/gdna_w1_substrate.py` | the substrate arms |
| `scratchpad/gdna_w2_bandwidth.py --plot` | the kernel sweep + the ladder/curve figures |
| `scratchpad/gdna_w3_projection.py --plot` | projection vs IDEAL + the σ-calibration |
| `scratchpad/gdna_w1a.py` | the two scoring confounds |
| `scratchpad/gdna_undercall.py`, `gdna_enriched_mode.py`, `gdna_resurrect.py` | under-call table, enriched-mode census, Goal-1/2 summary |

**Caches** live in `scratchpad/gdna_substrate_cache{,_quick}.pkl` (durable; they were once lost to a
volatile session dir). Everything above is **pure numpy on the cache — no `calibrate` re-runs.**

⚠ **`/tmp/rigel_selfsolve` is a SHARED, non-namespaced path** for the oracle's split BAMs. Two sessions
running scenario work simultaneously corrupt each other (`Invalid BGZF header` → `oracle INVALID:
region_contained partitions do not sum to full`). This cost two full rebuilds this session.

**Key library entry points** (`scripts/debug/gdna_explore_lib.py`): `recipe(s, sel=, w=, h_pop=, knn_scale=)`
· `recipe_substrate` · `recipe_weights` · `oracle_landscape(s, sel=)` · `oracle_pointmass` ·
`two_component` · `modes` (⚠ eyeballing only) · `spread` · `emd`.

## 6. Invariants (unchanged)

`HANDOFF_15` §5 in full, plus: **no magic numbers** — the constants ledger in the production plan tracks
every one with status, and three were deleted this session; **`_relay`/`_transport` and the scalar twins are
deliberate pairs, do not merge**; **pass-0 must be WEAK and correctable**; **goldens LAST**.
