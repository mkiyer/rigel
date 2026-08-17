# The RULER is measurable now — and the panel says a "perfect" one is a LOSS in scope

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When something settles, MOVE it to its
    permanent home and delete it here in the same edit.

    ⭐ Written 2026-08-14, continuing the calibration handoff of the same day. Everything below was
    measured on the 16-condition ladder at `fdde1477`, messages OFF, EM seed pinned.


## ⭐ WHERE THE SETTLED PARTS WENT — this file is a RECORD, not the state

⛔ Everything here that decides anything has a permanent home and should be read there:
the ruler's per-stratum numbers and the `U` no-enrichment null are `ROADMAP.md` §0; what is still
open about it is `ROADMAP.md` §1 **rank 1**; the instrument is `CLAUDE.md`'s row for
`design/calibration_vs_oracle.py`. ⚠ Kept because the working-through is occasionally useful, but
**nothing may cite it** and a number here that disagrees with `ROADMAP.md` is stale by definition.

---

## 0. ⭐⭐⭐ THE ONE RESULT — AND IT IS THE OPPOSITE OF WHAT ONE CONDITION SAID

**A perfect ruler is a large win at the zero control and on the DEFERRED stratum, and a 2.5–2.7×
LOSS on the two capture-OFF strata that 0.8.0 actually ships.**

⛔⛔ **THE HONEST ORDER OF EVENTS, KEPT BECAUSE IT IS THE LESSON.** The one-condition run below was done
first and read as "every recorded ceiling was measuring a fifth of the win". The panel then said no. The
control was misleading about the panel, which is `TRAPS: a-single-level-panel-cannot-see-a-constant`
met on a ceiling instead of a mechanism — and `TRAPS: panel-before-src` is why no mechanism was written
on the strength of it.

| stratum | transcript `oracle/base` | `oracle_ruler/base` | gene `oracle_ruler/base` | fp_mass `ruler/base` |
|---|---|---|---|---|
| stranded × capture-OFF ⭐ | 1.022 | **2.510** | 1.697 | **2.955** |
| stranded × capture-ON ⭐ | 0.936 | 1.077 | 0.816 | 1.048 |
| unstranded × capture-OFF ⭐ | 1.005 | **2.678** | 2.047 | **3.340** |
| unstranded × capture-ON ⛔ DEFERRED | 0.267 | 0.259 | 0.104 | 0.137 |
| ⛔ `g00` ZERO-gDNA control | 0.928 | **0.207** | **0.060** | **0.074** |

⭐⭐ **WHY A "PERFECT" RULER LOSES, and the instrument's own null column already contained the answer.**
At capture-OFF there are no probes, so the correct contraction factor is **exactly 1.000**, and the
uniform-gDNA null shows the estimator can reach it (**0.9963–0.9967**). Both arms sit far below:
shipped **0.9318 / 0.9391**, oracle **0.9241 / 0.9193**. ⛔ **The shipped ruler is CLOSER to the truth
than the oracle ruler**, so `oracle_ruler` is not a ceiling at capture-OFF — it is an A/B between two
wrong rulers, and it picks the worse one.

⭐ **So the ceiling arm that still needs building is the `U` RULER**: the uniform-gDNA field, factor
≈ 1.000, derivable with no fitting and no constant. Admissible at capture-OFF only.

⚠ **And `Σeff/Σfl` HIDES the effect that matters.** A 1–2 % aggregate gap costs 2.5× at transcript
level, so essentially every transcript is being redistributed while the total barely moves. Read
`ruler_n_moved`; do not rank on the aggregate factor.

---

### The one-condition run that started it

Measured on `g00 ss_0.50 capture_off` — unstranded × capture-OFF, an **IN-SCOPE** stratum, EM seed
pinned, one condition. ⛔ Kept as the record of a reading the panel then overturned:

| arm | what it substitutes | transcript Σ\|Δ\| | vs base | gene Σ\|Δ\| | vs base | fp_mass |
|---|---|---|---|---|---|---|
| `base` | nothing | 4,697,940 | 1.000 | 1,383,032 | 1.000 | 765,545 |
| `oracle` | the **prior** (all three `LocusPriors` fields) | 3,769,385 | 0.802 | 217,991 | 0.158 | 782,632 |
| ⭐ `oracle_ruler` | the prior **AND** the ruler | **201,700** | **0.043** | **6,344** | **0.0046** | **12,883** |
| `oracle_ruler_noop` | nothing, through the same wrapper | 4,697,940 | 1.000 | 1,383,032 | 1.000 | — |

⭐ **`oracle_ruler` − `oracle` is the effective-length shrinkage and NOTHING else.** `LocusPriors` has
exactly three fields and `oracle` already takes all three from O, so the only difference between the
two arms is whether the length the EM divides by was built from the true split or the shipped one.

⛔ **Read the `oracle` row twice.** A perfect prior makes false-positive mass *worse* (765,545 →
782,632), which is the recorded "a perfect prior is neutral-to-worse on the capture-OFF strata". The
reason is now visible: the prior was being divided by a ruler contracted 10.5×.

⛔⛔ **AND THE PANEL SAYS THIS ROW DOES NOT GENERALISE — §0's table is the verdict.** The 0.043× holds
at the control and on the deferred stratum; in scope, at capture-OFF, the same arm is 2.5–2.7× WORSE.
⚠ The `P/O` = 1.008 / 0.997 / 1.022 reading that suggested "the ruler is nearly right on the
contaminated rows, so expect a smaller win" was right that the aggregate barely moves and wrong about
what follows from it: the aggregate is stable while nearly every transcript is redistributed.

---

## 1. WHAT WAS BUILT

### 1a. `scripts/design/calibration_vs_oracle.py` — the 0.8.0 metric

`P = calibrate(...)` against `O = dataclasses.replace(P, **override_masses(ra))` on the same payload.
Reuses `pass0_vs_oracle.score_axis` and `prior_vs_oracle`'s strata machinery — no second scorer — and
deliberately does **not** re-score `LocusPriors`, which is `prior_vs_oracle.py`'s.

⭐ No solver, no EM, no BAM re-scan: **~5–12 s per condition, the whole ladder in ~2 min.** That is
the iteration loop the scope asked for.

⭐ It is the first instrument to stamp the IN-SCOPE / DEFERRED ruling onto every row, and the first to
reach the effective-length shrinkage at all.

**Falsification:** `--self-test` runs 21 perturbations with no I/O (21/21), and
`tests/test_calibration_vs_oracle_self_test.py` runs them in-process on every suite run so the flag
cannot go stale. The `noop` arm is byte-identical to `P` on all six override arrays **and** on the
effective lengths derived from them, on 16/16 conditions.

### 1b. `oracle_ruler` / `oracle_ruler_noop` in `quant_accuracy.py`

Wraps `rigel.calibration.calibrate` as a module attribute — `run_pipeline` imports it function-locally,
so the name resolves at call time. The substituted `CalibrationResult` then reaches **both** consumers.

⛔ **The counter watches the SHRINKAGE, not `calibrate`.** `calibrate` being called is necessary and
not sufficient: hand `_setup_geometry_and_estimator` a `None` calibration and the arm silently becomes
`oracle` under a different name. The arm additionally requires that the substituting variant MOVED the
ruler (measured: `max|Δ| = 1,115,202` bp) and that the noop did not (measured: exactly `0.0`).

⚠ `oracle_ruler_noop` reproduces `base` exactly on `count_abs_err` at both axes; `fp_mass` differs by
**2 fragments**, which is the EM's own re-association and is ~4 orders below the recorded
`base_reseed` noise floor.

### 1c. The `g00` oracle fallback, in `quant_accuracy.load_oracle`

⛔ **No oracle arm could reach the zero control at all before this.** `pass0_vs_oracle.py` holds every
zero-gDNA condition out, so it never wrote a `_main`. `_main` is the undrained full payload — the same
quantity as the plain scan cache — so it is read from `scan_cache/<condition>` when absent, and
`from_parts` re-runs sum-to-full either way (it passes on all four).

---

## 2. THE ERROR SPECTRUM, whole panel, per stratum

⛔ Composition (`mwae` = mass-weighted |Δf_g|, P against O):

| stratum | region | boundary | worst condition in it |
|---|---|---|---|
| stranded × capture-OFF ⭐ | 0.0065 | 0.0121 | `g98` boundary 0.0688 |
| stranded × capture-ON ⭐ | 0.0109 | 0.0229 | `g98` 0.0139 / 0.0286 |
| unstranded × capture-OFF ⭐ | 0.0085 | 0.0151 | ⭐ **`g00` 0.2558 / 0.1194** |
| unstranded × capture-ON ⛔ DEFERRED | 0.5070 | 0.5630 | `g98` 0.9202 / 0.9326 |
| ⛔ `g00` control, all strata | 0.0726 | 0.0407 | +1,417,767 gDNA fragments, **all over-call** |

⭐⭐ **The single largest IN-SCOPE composition error on the panel is `g00 ss_0.50 capture_off` at
0.2558** — 6× the next in-scope number and **67× its own stranded twin** (0.0038) on identical
geometry and identical (zero) gDNA. The worst condition overall is `g98 ss_0.50 capture_on` at 0.9202,
which is the DEFERRED stratum; taking it is how the ranking gets inverted.

⭐ The RULER, per stratum (`factor` = Σ eff_em / Σ fl; 1.000 = no contraction):

| stratum | factor P | factor O | factor U | P/O |
|---|---|---|---|---|
| stranded × capture-OFF ⭐ | 0.9318 | 0.9241 | 0.9967 | 1.008 |
| stranded × capture-ON ⭐ | 0.0794 | 0.0797 | 1.0000 | 0.997 |
| unstranded × capture-OFF ⭐ | 0.9391 | 0.9193 | 0.9963 | 1.022 |
| unstranded × capture-ON ⛔ | 0.5959 | 0.0801 | 1.0000 | 7.443 |
| ⛔ `g00` control | **0.0951** | **1.0000** | 1.0000 | **0.095** |

⛔ **At the zero control the EM's ruler is contracted 10.5× on a library containing no gDNA** —
1.06 **billion** bp of opportunity. On the contaminated in-scope rows it is within 2.2 % of perfect.

⚠ **TWO AGGREGATES, and they disagree by 3.6×.** `Σeff/Σfl` reads **0.095** at `g00` where the
unweighted `mean(eff/fl)` — the form recorded as **0.345** — reads 0.345. Neither is wrong: the
contraction falls hardest on LONG transcripts, so weighting by the opportunity it scales makes it look
far worse. The instrument carries both; rank on `Σeff/Σfl`, which is what the EM divides by.

⭐ **`U`, the no-enrichment null, settles a question the ruler table would otherwise raise.** It is O's
own gDNA total laid down at exactly uniform density — at capture-OFF the physically correct field, so
the contract demands exactly 1.000. It reads **0.9963–0.9967**. So under perfect composition the
estimator manufactures only ~0.4 % of spurious contraction from sampling noise, and **the shrinkage
really is a symptom of the composition, as `SUCCESS.md` and `ROADMAP.md` rank 2 both say.** A separate
shrinkage repair remains the wrong move.

---

## 3. THE TARGET, AND WHETHER A CHANNEL EXISTS

The target is `g00 ss_0.50 capture_off`. Dissected against its stranded twin, region axis, truth = 0:

| population | regions | total mass | FALSE gDNA (ss 0.50) | FALSE gDNA (ss 0.99) |
|---|---|---|---|---|
| `g1_locked` — RNA **inadmissible** | 1,312 | **0** | **0** | **0** |
| TS_POS | 14,774 | 1.96 M | 533,654 (f_g 0.2720) | 4,099 (0.0021) |
| TS_NEG | 13,808 | 1.54 M | 382,935 (0.2485) | 3,787 (0.0025) |
| TS_AMBIG | 5,241 | 1.08 M | 255,829 (0.2368) | 9,355 (0.0087) |

⛔ **The answer is that NO channel reaches an RNA-admissible region here, and it is structural.**

1. **Strand is dead.** `rna_sense_frac` fits to 0.500369, so the λ-term carries exactly zero
   information. The stranded twin gets 130× / 99× / 27× lower error on the same geometry using
   precisely the channel this condition does not have. Its own residual is worst on **TS_AMBIG**
   (0.0087, 4× the others) — the class strand cannot resolve even at ss 0.99, which is the mechanism
   confirming itself.
2. **The pure-gDNA anchors are EMPTY.** All 1,312 of them hold zero mass. An empty object measures
   nothing, and every rule for making an evidence-free slot *say* something is `ROADMAP.md` §4.1's
   graveyard — eleven built, eleven refused, all for lifting the mode UP at `g00`.
3. **The density channel works at every other rung** (`g05` / `g50` / `g98` capture-OFF read 0.0043 /
   0.0096 / 0.0098). It fails only where there is no gDNA to find.

### ⭐ One thing that is NEW and is not one of the eleven

The `g1_locked` population's zero is not an absence of information — it is a Poisson observation over
**50.8 Mbp** of effective support, and it bounds the density from above with no fitting, no threshold
and no new constant:

| condition | Σ support | OBSERVED | `gdna_density_global` | PREDICTED | obs/pred | 95 % UB on ρ |
|---|---|---|---|---|---|---|
| `g00 ss0.50 off` | 50,755,315 | **0** | 1.902e-02 | **965,230** | 0.0000 | 5.90e-08 |
| `g00 ss0.99 off` | 50,755,390 | **0** | 3.234e-04 | 16,415 | 0.0000 | 5.90e-08 |
| `g05 ss0.50 off` | 50,745,404 | 260,253 | 4.856e-03 | 246,432 | **1.0561** | 5.13e-03 |
| `g50 ss0.50 off` | 50,745,570 | 2,601,241 | 5.022e-02 | 2,548,321 | **1.0208** | 5.13e-02 |

⭐ **The falsification passes**: the same population predicts its own observed count to **5.6 % and
2.1 %** the moment the library really contains gDNA, so it is informative rather than structurally
empty. At `g00` the solve claims a global density **322,000×** larger than its own direct measurement
permits.

⭐ **Why this is not a twelfth doubt-resolution rule**: all eleven refused mechanisms *lifted* an
evidence-free slot off zero. This is a constraint from a population that HAS evidence, and it pulls
DOWN.

⛔⛔ **BUT IT IS A DIAGNOSTIC, NOT A LEVER, AND THAT HAS TO BE SAID PLAINLY.**
`calibration/derive.py` computes `gdna_density_global` **from** the per-object solve — it is a summary,
and its only consumers are `cli.py`, a log line and the QC report. Nothing in the solver reads it, and
`density_model.py`'s docstring says so explicitly ("no density→deconv→density feedback loop"). Fixing
the scalar would change a reported number and nothing else. A real lever would have to sit inside the
solve, which is exactly where the eleven died.

⚠ **And it is admissible only at capture-OFF.** `anchor_opportunity_census.py`'s recorded verdict is
that the empty-anchor density claim is false by **346×** under capture and true off it — probes
concentrate gDNA, so "the gDNA would also have landed intergenically" is simply untrue with capture on.
That still covers two of the three in-scope strata.

---

## 4. WHAT IS NEXT, and what was deliberately NOT done

1. ✅ **DONE — the ladder run is §0's table**, and it inverted the one-condition reading.
2. ⭐⭐⭐ **BUILD THE `U` RULER ARM.** The oracle ruler is not the ceiling at capture-OFF; the uniform-gDNA
   null is, it reads 0.9963–0.9967 against a contract of 1.000, and it needs no fitting. That is the arm
   that would actually price the shrinkage in scope. ⚠ capture-OFF only.
3. ⭐ **Re-price the ceilings in `ROADMAP.md` §0 and `SUCCESS.md` only AFTER that arm exists** — the
   `oracle_ruler` numbers price an A/B between two wrong rulers, not a ceiling.
3. ⛔ **No mechanism was written into `src/`.** `TRAPS: panel-before-src` has caught four toy-positive,
   panel-negative changes, and the target sits next to an eleven-deep graveyard. The measurement is the
   deliverable; the repair needs an owner call.
4. ⛔ **The length channel was not proposed, listed or ranked** — retired until after 0.8.0.
5. ⛔ **unstranded × capture-ON was measured and reported on every table and optimised for nothing.**

---

# PART 2 — THE ERROR CENSUS, AND ONE BIAS THAT EXPLAINS EVERY TREND IN THE PANEL

    ⚠ Same dev doc, same day, second half of the session. Owner asked for the error in FRAGMENT
    COUNTS rather than ratios, per pool, with the full signed distribution, then the worst-scenario
    ranking and a dissection.

## 5. WHAT WAS ADDED TO `calibration_vs_oracle.py`

Table ⓪ the POOL LEDGER (fragments, three truth pools), table ⑤ the SIGNED DISTRIBUTION per object
(log-spaced, exact-zero bucket of its own), table ⑥ SCENARIOS RANKED worst-first by Σ|Δ| in fragments.
Self-test 21 → **29 gates**, all firing.

⛔ **Two structural limits found while building it:**
* **The panel has NO nascent RNA** — `nrna: ratios: [0.0]`, `nrna: 0` in all 16 `truth_summary.json`.
  A gDNA-vs-nascent-vs-annotated verdict is not obtainable from this ladder.
* **Calibration emits TWO populations, not three.** It cannot split mature from nascent — that is the
  EM's job. So the three-pool ledger is truth-side accounting and the answer side is gDNA vs RNA.

⭐ **`Δ_RNA ≡ −Δ_gDNA` per object** (gated, max deviation 1.7e-11), because both arms carry the same
per-object total. ONE signed number describes an object's whole error; two columns would read as two
independent measurements of one thing.

## 6. ⛔⛔ A TARGET I RECOMMENDED ON A RATIO ARTEFACT, AND THE CORRECTION

I proposed the **g98 capture-OFF boundary axis** on the strength of `mwae` 0.0910 vs 0.0098 — a "9.3×
axis asymmetry". ⛔ **In FRAGMENTS the two axes are nearly EQUAL** (68,244 region vs 59,599 boundary at
`g98 ss0.99 capture_off`). The boundary axis's MASS SHARE collapses with gDNA level:

| condition | region Σ\|Δ\| | boundary Σ\|Δ\| | region mass share | boundary mass share |
|---|---|---|---|---|
| g05 ss0.99 off | 17,442 | 15,184 | 45.0 % | 55.0 % |
| g50 ss0.99 off | 52,151 | 49,001 | 67.1 % | 32.9 % |
| g98 ss0.99 off | 68,244 | 59,599 | **91.6 %** | **8.4 %** |

⭐ gDNA is genomically uniform, so it is mostly CONTAINED in long regions; RNA sits in short exons and
crosses boundaries constantly. As gDNA rises the boundary axis loses mass, and the same fragment error
reads 9.3× worse as a rate. **There is no boundary-specific defect.** This is precisely why the owner
asked for fragment counts.

## 7. ⭐⭐⭐ THE REAL FINDING — A REPULSION FROM THE `f_g = 1` VERTEX, CARRYING MOST OF THE PANEL

Bucketing every object by its TRUE `f_g`, capture-OFF, REGION axis:

| condition | objects at true `f_g ≥ 0.999` | their mass | `mean(true − pred)` | share of Σ\|Δ\| |
|---|---|---|---|---|
| g05 ss0.99 off | 10,457 | 465,687 | **+0.1529** | **76.0 %** |
| g50 ss0.99 off | 13,432 | 4,658,232 | **+0.1101** | **76.2 %** |
| g98 ss0.99 off | 16,277 | 9,167,899 | **+0.1197** | **85.0 %** |
| g98 ss0.50 off | 16,396 | 9,168,509 | **+0.1299** | **86.1 %** |

⭐⭐ **57–86 % of ALL calibration error on every capture-OFF condition sits on objects whose truth is
`f_g ≈ 1`, read 0.07–0.23 BELOW the vertex.** At `g98` that bucket is **97.3 %** of the region axis's
whole mass.

⭐⭐⭐ **AND IT IS A PULL TOWARD THE MIDDLE, MEASURED AT BOTH ENDS:**

| condition | at true `f_g ≈ 0`: `mean(pred − 0)` | at true `f_g ≈ 1`: `mean(1 − pred)` |
|---|---|---|
| g00 ss0.99 off | **+0.0327** | — (no such objects) |
| g05 ss0.99 off | +0.0191 | +0.1529 |
| g50 ss0.99 off | +0.0704 | +0.1101 |
| g98 ss0.99 off | **+0.3339** | +0.1197 |

Both signs positive: below-truth at the top vertex, above-truth at the bottom one. **One shrinkage
bias, seen from two ends** — and it explains every trend the census found: why the error is DIFFUSE,
why it is systematically an UNDER-call, why it grows monotonically with gDNA level (more mass sits at
`f_g ≈ 1`), and why `g00` is 100 % OVER-call.

## 8. ⭐⭐⭐ THE CODE, AND IT DECLARES THIS COST ITSELF

`calibration/simplex_logodds.py:71`, `_JEFFREYS_REF = 0.5` — the Beta(½,½) reference. Its own comment:

> *"A DECLARED CHOICE, not forced by the likelihood … Licensed as the 'structural Jeffreys' prior;
> §10.5 records the known cost — **it forbids the simplex vertices, where some truth genuinely lives.**"*

And `_rna_arm`: *"the reference is what bounds the `f_g → 1` vertex today, and **it is the ONLY thing
doing so**."*

⭐⭐ **THE ASYMMETRY IS THE MECHANISM.** `_gdna_arm` is `+½·log f_g` **plus a fitted `logP_g`**, so as
evidence accumulates the reference is swamped. `_rna_arm` is `+½·log(1−f_g)` **and nothing else** —
*"there is no parameter to pass one, by design, because nothing produces it today."* So the `f_g → 1`
repulsion never weakens however much evidence an object has, while its twin at `f_g → 0` does.
⭐ That is why the top vertex is 3–7× worse than the bottom one, and the module's own `_gdna_arm`
docstring already argues against exactly this kind of asymmetry in the other direction.

## 9. THE CEILING — PRICED, AND THE ZERO CONTROL REFUSES THE OBVIOUS FIX

Ablating `_rna_arm` to zero and re-solving (a scratch monkeypatch; nothing was written to `src/`):

| condition | base Σ\|Δ\| | RNA-ref ablated | ratio |
|---|---|---|---|
| ⛔ g00 ss0.99 off | 32,430 | 72,060 | **2.222 WORSE** |
| g05 ss0.99 off | 32,626 | 25,047 | 0.768 |
| g50 ss0.99 off | 101,152 | 50,485 | **0.499** |
| g98 ss0.99 off | 127,843 | 46,569 | **0.364** |
| g98 ss0.50 off | 170,437 | 69,363 | **0.407** |

⭐ **Deleting the term halves-to-thirds the error on every contaminated condition and DOUBLES it at the
zero control.** ⛔ `TRAPS: a-cancelling-defect-pair`, and the owner-required control refuses it exactly
as it refused §4.1's eleven.

⭐⭐⭐ **BUT THE ABLATION IS NOT THE PROPOSAL, AND THIS IS THE DISTINCTION THAT MATTERS.** The eleven
graveyard mechanisms were rules for RESOLVING DOUBT at an evidence-free slot. This is a MEASURE term
that evidence is never allowed to overwhelm. ⭐ **The symmetric repair is to FIT `logP_r`**, exactly as
`logP_g` is fitted and ADDED on top of its own reference — then at `g00` a fitted RNA prior says "lots
of RNA" and reinforces the repulsion (the control holds), and at `g98` it says "almost no RNA" and
swamps it (the vertex becomes reachable). ⚠ **No new constant**: it is the machinery the gDNA arm
already uses. ⛔ **UNBUILT AND UNPRICED** — the table above prices the TERM, not the REPAIR.

## 10. NEXT

1. ⭐⭐⭐ **Price a fitted `logP_r`** — an arm first, `--messages off`, both controls, all 16 conditions.
2. ⚠ It must be scored against BOTH zero controls, since the term it replaces is what makes `g00` work.
3. ⛔ Do not delete `_rna_arm`. The ablation is the ceiling, not the fix.
