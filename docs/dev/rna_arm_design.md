# DEV — ψ's COMPOSITION REFERENCE, THE LIVE DESIGN

    ⚠ **A DEV DOC. Nothing may cite it** — not the source, not a test, not one of the six permanent docs.
    ⛔ It is the WORKING half only: everything settled has MOVED to a permanent home and been deleted
    here, which is the rule that stops a dev doc becoming the state.

## ⭐ WHERE THE SETTLED PARTS WENT (2026-08-16) — this file no longer carries them

| what | home |
|---|---|
| the derivation: ψ's reference is Beta(a,b), its mean is a missing THIRD term, and the floor is one pseudo-fragment AT THE OBJECT | `EQUATIONS.md` §9c |
| the rulings: the reference is per-object · RNA is the RESIDUAL and never predicted · the boundary axis splits on mature-RNA crossability · the reference and the landscape PARTITION the object universe · capture needs no detection step | `DESIGN.md` §6b |
| the lessons | `TRAPS.md`: `a-mean-hits-the-mass-weighted-centre-by-luck` · `a-clamp-at-the-closed-end-escapes-the-window` · `the-peel-is-as-good-as-the-density-it-is-handed` |
| ⛔ the FIVE refused mechanisms — a fitted `logP_r` · a library-wide Beta mean · the object-weighted mean · a stratified assertion · a pooled RNA density; plus the false `f_g ≤ 1 − S/M` bound | `ROADMAP.md` §4.2 — **do not rebuild these** |
| where the tool IS, and what to do next | `ROADMAP.md` §1 rank 2 |
| the instrument | `CLAUDE.md`'s row for `design/object_composition.py` |

⛔ **Parts 1–6 of this file were DELETED in the same edit that wrote those homes** (1,225 lines). They
were the refused RNA-density arm, the library-wide Beta mean, the `f_lib` circularity, the stratified
sample, and the strand-split finding — all of which now point forward from `ROADMAP.md` §4.2 rather than
sitting here as history.

---

# PART 7 — ⭐⭐⭐ THE REFERENCE IS PER-OBJECT (owner, 2026-08-15). MEASURED: **8×**

    ⭐ The owner's question — *"why does the prior have to be fixed for the entire library? The solver is
    solving one region or one boundary at a time"* — and the answer is that it does not, and that Parts
    3–6 were all searching the wrong space. ⛔ Still NOTHING in `src/`.

## 33. ⭐⭐ THE DERIVATION — the reference has a THIRD term, and the shipped one is its degenerate case

Put the priors where the physics is: on the two DENSITIES, not on `f_g`. With
`ρ_g ~ Gamma(a, β_g)`, `ρ_r ~ Gamma(b, β_r)` and Poisson counts, `X = ρ_g E_g ~ Gamma(a, s_g)` and
`Y = ρ_r E_r ~ Gamma(b, s_r)` where `s_c = β_c/E_c`. Marginalising the total out of `f = X/(X+Y)`:

    p(f) ∝ f^(a−1) (1−f)^(b−1) / ( s_g·f + s_r·(1−f) )^(a+b)

and on the solver's grid, where `|df/dλ| = f(1−f)`:

    log p_i(λ) = a·log f_g + b·log(1−f_g) − (a+b)·log( f_g + r_i·(1−f_g) ) + const

⭐⭐⭐ **The shipped reference is exactly this with `s_g = s_r`.** The code has silently committed to the
two components having MATCHED SCALES — the same assertion as "the library is half gDNA", now visible as a
MISSING TERM rather than a chosen number. With `m_i` the object's prior expected composition,
`r_i = (b/a)·m_i/(1−m_i)`, and keeping Jeffreys' `a = b = ½`:

    log p_i(λ)  =  ½·log f_g  +  ½·log(1−f_g)  −  log[ (1−m_i)·f_g + m_i·(1−f_g) ]

Four properties, all structural rather than measured:

* ⭐ **`m_i = ½` reduces it EXACTLY to the shipped constant.** A strict generalisation.
* ⭐⭐ **The tails stay `e^(−|λ|/2)` for every `m_i` — only the LOCATION moves.** So `L`-invariance, which
  `simplex_logodds.py` calls its acceptance test, is preserved exactly, and the ~0.9 % of prior mass
  outside `L = 10` is unchanged. ⛔ **This is what Part 3's `(a, b)` rule could not do**: moving the
  shape parameters moves the location AND the tails together, which is why `b = 0.028` put **57 %** of
  the mass outside the window, needed a clamp, and had a 50-nat range the AMBIG float32 path cannot hold.
  The third term's range is `log(max/min)` — bounded, and float32-safe.
* ⭐ **Proper for every `m_i` in (0,1)**, degenerate only where a measured density is exactly zero. The
  regulariser is one pseudo-fragment over the pooled opportunity — the convention `fit_landscape` already
  uses. **No new constant anywhere.**
* ⭐⭐ **`R_own` becomes `m_i` in CLOSED FORM.** `u = f/(f + r(1−f)) ~ Beta(a,b)`; at `a = b` that is
  symmetric, so the prior median of `f_g` is exactly `m_i`. §15b's hard-coded `0.5` is re-derived rather
  than widened, and it becomes the quantity the relay's budget arithmetic actually wants.

⚠ **And nothing had to be built to carry it**: the gDNA arm's fitted term is already `(n_slots, K)`. The
reference is the only scalar left in ψ.

## 34. ⭐⭐⭐ THE MEASUREMENT — `Σ|m_i − f_g,i|·M_i` in FRAGMENTS, ratio to the shipped ½

`object_composition.py` table ⑥. This scores the PRIOR alone, with no solver: how many fragments it
misplaces if believed outright.

| stratum | shipped Σ\|Δ\| | ⭐ **prior-free `m_i`** | truth-fed ceiling |
|---|---|---|---|
| stranded × capture OFF ⭐IN SCOPE | 13,160,156 | **0.123** | 0.081 |
| stranded × capture ON ⭐IN SCOPE | 13,982,878 | **0.401** | 0.346 |
| unstranded × capture OFF ⭐IN SCOPE | 13,161,800 | **0.123** | 0.081 |
| unstranded × capture ON ⛔DEFERRED | 13,984,302 | 0.401 | 0.346 |
| ⛔ `g00` ZERO-gDNA CONTROL | 3.78–3.86 M each | **0.096 – 0.103** | 0.000 |

⭐⭐ **8× better than the shipped constant at capture-OFF, 2.5× at capture-ON, 10× at the zero control,
and within 1.5× of the ceiling on every stratum.** ⭐ And `g05` — the condition that regressed **1.43×**
under every library-wide mean Parts 3–6 tried — reads **0.191**.

**The per-stratum decomposition, stranded × capture-OFF, is where the mechanism is visible:**

| stratum | shipped Σ\|Δ\| | per-object flux | ⭐ pooled `ρ_r` | ceiling |
|---|---|---|---|---|
| R intergenic | 3,981,874 | **0.000** | **0.000** | 0.000 |
| R intron | 3,055,202 | **0.000** | **0.000** | 0.000 |
| R exon | 3,467,794 | 0.935 | **0.087** | 0.128 |
| B exon\|intron | 325,028 | **0.000** | **0.000** | 0.793 |
| B exon\|exon | 2,286,957 | ⛔ 1.547 | **0.579** | 0.146 |
| B gene edge | 43,300 | **0.000** | **0.000** | 0.791 |

⭐⭐⭐ **THE STRUCTURAL HALF IS EXACT AND THE ESTIMATED HALF IS THE WHOLE RESIDUAL.** Four of the six
strata go to **0.000** — the ones where the annotation determines the answer. Everything left is the two
EXONIC strata, where the RNA density has to be estimated. ⭐ On two of those four the derivation BEATS
its own truth-fed ceiling (0.000 vs 0.793), because the class-pooled truth arm hands an on-target RNA
density to objects mature RNA cannot occupy — the ceiling arm is the one that is wrong there.

⭐⭐⭐ **AND THE PER-OBJECT VARIATION IS CARRIED BY THE OPPORTUNITY GEOMETRY, NOT BY PER-OBJECT
DENSITIES.** Substituting the per-object sj flux for ONE pooled RNA density is **worse** — 0.515 against
0.123. The flux is at the right LEVEL and far too NOISY per object. **Two pooled scalars plus exact
per-object geometry is the entire mechanism**, which is why it needs no landscape, no substrate selection
and no fit.

## 35. ⛔⛔ THE OWNER'S CURATION IS WHAT MADE IT WORK, AND IT DOES TWO JOBS

> *"the sj-boundary unspliced pool must be carefully curated … some boundaries are exon-exon.
> Exon-intron boundaries are special."*

⛔ Part 5's stratum was "a sj attaches here", and it measured true `f_g` = **0.0000 over 955,428
fragments** at the zero control. The pool lumps `exon|exon` boundaries in — an alternative splice site
inside a contiguous exonic stretch, which mature RNA crosses **freely**. Split on whether mature RNA can
cross, and the anchor is clean: `B exon|intron` reads `ρ_g` **0.050249 against a true 0.050327**.

⭐⭐ **The same predicate then does a second, unanticipated job.** `ρ_r` from the sj flux is a MATURE
density on the SPLICED template, so handing it to a slot mature RNA cannot occupy subtracts RNA that is
not there and calls gDNA RNA. Gating it by the solver's own `mrna_active` — one shipped predicate,
correct on both axes — takes `R intron` from **0.498 → 0.000** and `B exon|intron` from **0.485 → 0.000**.
⛔ It is gated in `strata()` against `mrna_active` so the classification and the solver cannot separate.

## 36. ⭐⭐ HYBRID CAPTURE NEEDS NO DETECTION STEP

> *"When we have hybrid capture, we concentrate the DNA around exons and our exon-intron boundaries
> become valuable. When we don't, we spread the DNA uniformly and our intergenic regions become valuable."*

⭐ **Measure `ρ_g` at BOTH anchors and the RATIO is the enrichment** — no flag, no threshold, no
detection (table ⑤):

| | capture OFF | capture ON |
|---|---|---|
| ⭐ enrichment `ρ_on/ρ_off`, measured | **0.98** | **113 – 114** |
| off-target anchor, estimate ÷ truth | 1.00 | **0.91 – 1.00** |
| in-gene anchor, estimate ÷ truth | 0.80 – 1.00 | ⛔ **0.28 – 0.38** |

⭐ A **116×** separation between the two regimes, and the off-target anchor is now accurate on **all 16**
including capture-ON (0.000216 against a true 0.000239) — because it is being asked for the OFF-target
density rather than the library's.

⛔⛔ **BUT THE IN-GENE ANCHOR IS A DETECTOR AND NOT YET A CALIBRATED LEVEL.** It under-reads the
on-target gDNA density by **2.6–3.6×** under capture. The mechanism is geometric and testable: an
`exon|intron` boundary sits at the **EDGE of the probe footprint**, and its crossing opportunity `eff_gdna`
is built with UNBOUNDED reach, so it counts intronic positions no probe covers. **That is the next thing
to fix, and it is an opportunity-model repair rather than a new estimator.**

## 37. ⛔ WHAT I GOT WRONG, KEPT BECAUSE THE GATE CAUGHT IT

Part 6 §32 proposed `f_g ≤ 1 − S/M` as an assumption-free per-object bound from certified RNA. **It is
false.** `boundary_spliced` is a SEPARATE bank from `boundary_unspliced`, not a subset — the same
molecules split by whether they used a sj ELSEWHERE. The truth violated the "bound" by **302**, on
1,225,844 fragments, the first time it was measured. ⭐ The correct statement is an identity, and it is
more useful than the bound would have been:

    ρ_r · E_r  =  unspliced_RNA + S      ⇒      f_g = 1 − (ρ_r·E_r − S)/M

so **S SUBTRACTS from the estimated RNA crossing**. Measured, S is **20.8 %** of the contiguous RNA at a
boundary, and against the corrected target the sj flux reads **0.929** — where against the wrong target
it read 1.273. ⚠ Folding the S term in is currently a slight net LOSS (0.553 → 0.570) because it
over-corrects at `B exon|exon`; it is kept as a measured arm rather than adopted.

## 38. WHAT IS NOW OPEN

1. ⭐⭐⭐ **The panel arm.** Everything above scores the PRIOR; `vertex_ceiling.py` scores the SOLVE, and
   its `--arm ref_c=a,b` takes SCALARS. A per-object `m_i` needs the third term threaded to the solver
   before it can be A/B'd — the first `src/`-side work this design has justified, and it is small
   (`CompositionPriors` gains one `(n_slots,)` field; `_gdna_arm`/`_rna_arm` gain one term).
2. ⭐⭐ **`B exon|exon`, the last stratum above 0.5**, and `R exon` under capture (0.383).
3. ⭐ **The in-gene anchor's LEVEL** (§36) — an opportunity-model repair, not an estimator.
4. ⚠ **The nascent assumption is still unpriceable** (`nrna = 0` on all 16). ⭐ It is load-bearing on
   UNSTRANDED data only: with strand the nascent half of an `exon|intron` boundary is deconvolvable.


---

# PART 8 — ⭐⭐⭐ STAGE 0 RUN (2026-08-16): 0a, 0b AND 0c ANSWERED. **GO TO STAGE 1**

    ⛔ Still NOTHING in `src/`. Everything here scores the PRIOR alone; the SOLVE is stage 2.

## 39. ⛔ FIRST, THE OVERCLAIM PART 7 WOULD HAVE SHIPPED

Part 7 called `m_i` a per-object prior mean. At the STRUCTURAL strata that is true and exact. **At the
exonic strata it is a CLASS CONSTANT**, measured at `g50` capture-OFF:

| stratum | `m_i` sd | true `f_g` sd | corr | Σ\|Δ\|·M per-object | …vs its own class mean |
|---|---|---|---|---|---|
| `R exon` | **0.0021** | 0.4441 | +0.355 | 168,551 | **164,074** (better) |
| `B exon\|exon` | **0.0000** | 0.3994 | — | 143,297 | 143,297 (identical) |

⭐ And the exonic constant wins only because it sits near the **mass-weighted** centre — which is exactly
the criticism Part 6 levelled at the library-wide mean, recurring one level down. ⛔ **The win is four
structural strata going to exactly 0.000 plus a better constant at the exonic ones.** Say it that way.

## 40. ⭐⭐⭐ THE ARM LADDER AFTER STAGE 0 — `Σ|m_i − f_g,i|·M_i` in FRAGMENTS, ratio to the shipped ½

| stratum | shipped Σ\|Δ\| | class-pooled ref | `structural` | `pooled` | `flux+fallbk` | ⭐ `flux+shrunk` |
|---|---|---|---|---|---|---|
| stranded × OFF ⭐IN SCOPE | 13,160,156 | 0.081 | 0.437 | 0.042 | 0.121 | **0.040** |
| stranded × ON ⭐IN SCOPE | 13,982,878 | 0.346 | 0.827 | 0.332 | 0.393 | **0.326** |
| unstranded × OFF ⭐IN SCOPE | 13,161,800 | 0.081 | 0.437 | 0.042 | 0.122 | **0.040** |
| unstranded × ON ⛔DEFERRED | 13,984,302 | 0.346 | 0.828 | 0.332 | 0.394 | 0.326 |
| ⛔ `g00` CONTROL (each) | ~3.8 M | 0.000 | 1.000 | **0.000** | 0.088–0.130 | **0.000** |

⭐⭐ **25× at capture-OFF, 3× at capture-ON, and EXACTLY 0.000 at both zero controls.** `g05`, which
regressed 1.43× under every library-wide mean, reads **0.009**.

⛔⛔ **AND `class-pooled TRUTH` IS NOT A CEILING — STOP CALLING IT ONE.** The prior-free arm beats it on
every stratum (0.040 vs 0.081, 0.326 vs 0.346), because that arm hands an on-target RNA density to
objects mature RNA cannot occupy (`B exon|intron` 0.793, `B gene edge` 0.791 — both **0.000** under
`m_i`). It is a reference arm, not a bound.

⚠ **Two things changed the numbers from Part 7 and both are corrections, not tuning.** Dropping the
S-subtraction took `pooled` 0.158 → 0.072 at `g50` (it over-corrects at `exon|exon`); the shrunk local
flux then took it to 0.040 panel-wide.

**0a — the local RNA density: shrink it hard.** Raw local flux with a pooled fallback is **0.121**; one
pooled scalar is **0.042**; one pseudo-observation of the pooled rate blended with the local flux is
**0.040**. ⭐ The flux carries a little per-object signal and survives only heavy shrinkage. ⛔ And the
naive form that sets `rho_r = 0` where no sj is in reach reads **0.088–0.130 at the zero control** where
the shrunk form reads **0.000** — the coverage artefact Part 6 §32 already flagged, now priced.

**0b — the scope: EVERYWHERE.** Restricting the reference to the four annotation-determined strata reads
**0.437 / 0.827**, 10× worse. At pass-0 there is no landscape, so the exonic constant is strictly better
than ½ there. ⚠ Whether it should REMAIN once the landscape exists is a stage-2 question, not this one.

## 41. ⭐⭐⭐ 0c — THE POST-SOLVE UPDATE DOES NOT RUN AWAY, AND THE TEST IS A CONSERVATIVE BOUND

The one post-solve update proposed is narrow: re-estimate the ON-TARGET gDNA density from solved exons,
because the pre-solve `exon|intron` anchor under-reads it 2.6–3.6× under capture. ⛔ That is a single
scalar measured on exons and applied to exons — structurally the positive-feedback loop §18 says makes a
library-wide `f_lib` rule inadmissible.

⭐⭐ **Iterated with each object's OWN LIKELIHOOD REMOVED — a strict UPPER BOUND on the feedback, since
it is the prior believed outright with no data pulling against it — four starts spanning three decades
converge to the SAME fixed point to six decimals on all twelve contaminated conditions**, and to exactly
0.000 from every start at the zero control. `spread = max/min = 1.000` on every row.

| condition | true `rho_on` | anchor start | 10× low | 10× high | truth start |
|---|---|---|---|---|---|
| `g50 ss0.99 off` | 0.050327 | 0.032956 | 0.032956 | 0.032956 | 0.032956 |
| `g50 ss0.99 ON` | 0.702249 | 0.496148 | 0.496148 | 0.496148 | 0.496148 |
| `g98 ss0.99 ON` | 1.244639 | 1.231315 | 1.231315 | 1.231315 | 1.231315 |
| ⛔ `g00` (all four) | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 |

⭐ **The mechanism is the per-object geometry damping it**: raising `rho_g` cannot raise the share at an
object whose `rho_r·E_r` is large, so the map saturates. A library-wide `a = f_lib` has no such damping,
which is why §18's hazard is real there and not here.

⛔⛔ **A PASS HERE IS NECESSARY AND NOT SUFFICIENT, and the fixed point's VALUE is not the loop's value.**
The likelihood is absent by construction, so the real refit loop has strictly less feedback AND a
different fixed point. ⚠ Read the *convergence*, not the number: the bound's own fixed point sits at
0.65–0.90× of truth, which under capture is an improvement on the anchor (0.35× → 0.71× at `g50`) and at
capture-OFF is a regression (1.00× → 0.65×). **That comparison belongs to stage 4 and needs the solve.**

## 42. THE PLAN, UPDATED — stage 1 is unblocked

| stage | what | gate |
|---|---|---|
| ✅ **0** | 0a shrink the flux · 0b scope everywhere · 0c no runaway | DONE 2026-08-16 |
| **1** | thread the third term: `CompositionPriors` gains `location: ndarray \| None`; `_gdna_arm`/`_rna_arm` add `−log[(1−m)f + m(1−f)]` when set | ⛔ `None` ⇒ term absent ⇒ **BIT-IDENTICAL**, proven by `rename_identity.py --check` before any behaviour exists to confound it |
| **2** | ⭐ **STOP/GO — the panel arm on the SOLVE.** `m_struct` / `m_all` / `m_oracle` / `m_flipped` (must be worse) / `noop` (bit-identical) | the first number that scores the SOLVE; per stratum, both zero controls |
| **3** | the estimator into `src/`, behind a named switch, default OFF | the arm family reproduces the panel |
| **4** | the post-solve on-target update | 0c passed the bound; stage 2's harness runs the real loop |
| **5** | re-derive `R_own = m_i`; re-specify the tail-slope test as the RNA pseudo-count; MOVE §33 to `EQUATIONS.md` and the rulings to `DESIGN.md` | ⛔ neither assertion widened |

⛔ **THE ONE THING STAGE 2 MUST SETTLE THAT NOTHING ABOVE CAN:** whether the reference's location and the
gDNA landscape DOUBLE-COUNT. Both are location statements about the same object; the prior metric cannot
see it because the landscape does not exist at pass-0. `m_all` against `m_struct` **with the refit loop
on** is the arm that prices it.


---

# PART 9 — ⭐⭐⭐ STAGE 1 LANDED (2026-08-16): ψ's REFERENCE HAS A PER-SLOT MEAN, AND IT IS INERT

    ⭐ The first `src/` change this design has justified. **Behaviour is unchanged**: nothing sets the
    new field yet, and the term is not written when it is `None`.

## 43. WHAT LANDED

`CompositionPriors` gains `location: np.ndarray | None` — a per-slot SCALAR, `(m,)`, not an `(m, K)`
grid array. `select()` slices it with the arms; `regrid()` passes it through **unchanged**, because a
scalar has no lattice to be carried between and `regrid-in-the-right-variable` therefore cannot apply to
it at all. That is the point of storing it as a scalar.

`_location_term(lam, location)` returns `−log[(1−m)·f_g + m·(1−f_g)]` as `(m, K)`, formed by
`logaddexp` over `_log_fg`/`_log1m_fg` so it inherits their exactness in both depleted tails. It is added
beside the two arms at **both** ψ call sites — the single-strand solve and the AMBIG one, inside the
latter's float32 cast — because it is the COUPLING between the two components' scales and belongs to
neither arm.

## 44. ⭐⭐ THE BIT-IDENTITY GATE, AND IT CAN FAIL

On real data, `calibrate` over `g50 ss0.99 capture_off`, in one process:

* the term present with `location = None`, against `_location_term` **removed entirely**:
  ⭐ **19 of 19 output arrays BIT-IDENTICAL**;
* forcing a nonzero term: **10 of 19 arrays move, max\|Δ\| = 48.47**.

⛔ The second half is what stops the first being vacuous — a term that is never reached would also be
"bit-identical". Both halves are reproduced in-process by
`tests/calibration/test_reference_location.py`, on **both** ψ paths.

## 45. ⛔⛔ THE GATE CAUGHT AN OVERCLAIM I HAD ALREADY MADE TO THE OWNER

I wrote that the reference is "worth one pseudo-fragment, so it is a tie-breaker at a sparse object and
inaudible at a rich one, and its overlap with the fitted arms and the messages is bounded **by
construction**". ⭐ **The first half is true and now measured.** Swinging the location from 0.02 to 0.98
on a single-strand object at `kappa = 0.9` moves `f_g` by:

| fragments | 6 | 24 | 120 | 600 | 15,000 |
|---|---|---|---|---|---|
| swing | **0.7700** | 0.0565 | 0.0090 | 0.0009 | **0.0000** |

⛔⛔ **The second half is FALSE where it matters most.** On an AMBIG object — or on ANY object once
`kappa = ½` kills the strand channel — the swing is **0.95 at 6 fragments and 0.95 at 3,000**, and the
AMBIG one plateaus at 0.5718 rather than decaying. There is no composition evidence there, so the
posterior IS the prior at any depth, and the location does not nudge those objects, it **determines**
them.

⭐ That is correct Bayes rather than a defect, but it changes two things:

1. **The location's accuracy on evidence-free objects is load-bearing** in a way it is nowhere else —
   and those objects are exactly AMBIG plus everything on unstranded data.
2. ⛔ **The reference and the gDNA landscape do NOT have a bounded overlap there.** They are the only
   two voices on that population, so "the reference is negligible wherever anything else speaks" does
   not resolve the double-count for the objects the double-count actually concerns. **Only a solve arm
   can price which of the two should carry it** — which is stage 2, and this is now its sharpest
   question rather than a footnote.

Both halves are pinned as gates (`test_evidence_swamps_the_location_where_evidence_exists`,
`test_where_there_is_no_composition_evidence_the_location_is_the_answer`).

## 46. THE OTHER GATES, EACH WITH ITS PERTURBATION

36 cases in `tests/calibration/test_reference_location.py`:

* ⭐ `m = ½` makes the bracket the constant ½ ⇒ the term is constant and cancels. **A strict
  generalisation agreeing with the shipped constant exactly where its assumption is true.**
* ⭐⭐ **The prior median is `m` in closed form** — measured error **4e-6 … 5e-5** on the untruncated
  prior. This is what re-derives `test_relay_mass_pin`'s hard-coded `R_own = 0.5` instead of widening it.
  ⚠ Inside the shipped `L = 10` window the median is displaced by up to **0.0019**, which is a property
  of the window and two orders below the effect.
* ⭐⭐⭐ **The tail slope stays at the Jeffreys exponent for every `m`** — only the LOCATION moves — so
  `L`-invariance holds. Tested on the ANSWER: the median across `L` = 10/20/40 moves **0.0019** at the
  extremes and **exactly 0** at `m = ½`, against **0.021–0.041** for the refused `(a, b)` route.
  Mass outside the window: **0.009–0.043** here, **0.658** for `b = 0.03`.
* Proper at every `m`; bounded and float32-safe; monotone in the right direction; **both vertex spikes
  survive at every `m`** — a tilted two-spike prior, which is what a spike-and-slab population wants.
* `select` slices the location with the arms; `regrid` leaves it alone; defaults are pass-0.

⚠ Two of the gates were written with guessed thresholds and FAILED; both times the term was right and the
threshold was wrong, and re-specifying them against the measured structure is what produced §45's finding.

## 47. SUITE

**0 failed / 3,458 passed / 0 skipped / 10 xfail.** `+38` from 3,420, accounted: **36** own cases plus
**2** meta ones. ⚠ A `tests/calibration/` file is `+2` meta, not the `+3` a top-level `tests/` file gets —
`test_scripts_index` does not parametrise over it.


---

# PART 10 — ⭐⭐⭐ STAGE 2 (2026-08-16): THE ARM ON THE **SOLVE**. GO ON THREE STRATA, STOP ON ONE

## 48. THE PANEL, PER STRATUM, IN FRAGMENTS

`vertex_ceiling.py --arm ref_loc={noop,struct,pooled,local}`, all 16 conditions, `base` re-recorded in
session. ⭐ **`noop` traverses the whole new path and takes no location: BYTE-IDENTICAL to `base` on
16/16 conditions**, so every row below is attributable.

**The deliverable (`abs_err_all_final`):**

| stratum | base | `noop` | `struct` | `pooled` | `local` |
|---|---|---|---|---|---|
| stranded × capture OFF ⭐IN SCOPE | 261,620 | **1.000** | 0.381 | **0.321** | 0.324 |
| stranded × capture ON ⭐IN SCOPE | 588,140 | **1.000** | **0.659** | ⛔ 5.323 | ⛔ 4.761 |
| unstranded × capture OFF ⭐IN SCOPE | 334,848 | **1.000** | 0.363 | **0.297** | 0.301 |
| unstranded × capture ON ⛔DEFERRED | 18,436,888 | **1.000** | 0.800 | 0.805 | 0.804 |
| ⛔ `g00` ZERO-gDNA CONTROL | 2,312,277 | **1.000** | 1.000 | **0.001** | 0.001 |

**pass-0 (`abs_err_all`):** 0.574 / 0.752 / 0.928 / 0.860 / 0.997 for `struct`; 0.263 / ⛔ 4.547 /
**0.018** / 1.111 / **0.001** for `pooled`.

⭐⭐ **The zero control falls to 0.001 — a 1,000× reduction — on both pass-0 and the deliverable**, and it
needs no assumption at all: the anchors are empty there, so the peel reads `rho_g = 0` and every object
correctly claims no gDNA. ⭐ Unstranded × capture-OFF pass-0 reads **0.018**, a 55× reduction, which is
the stratum where `kappa = ½` leaves ψ nothing but the reference — so a correct reference is the whole
answer there, exactly as the mechanism predicts.

⛔⛔ **AND STRANDED × CAPTURE-ON IS 5.3× WORSE — A HARD STOP FOR `pooled`/`local`.** This is the failure
the prior metric already predicted: `object_composition.py` table ④ read the peel at **1.091** under
capture, because the on-target gDNA density is under-read **2.6–3.6×**. **The peel is exactly as good as
the density it is handed**, and pre-solve we do not have that density under capture. ⚠ `local` did not
rescue it (4.761), so the flanking-boundary density that measured 0.04–0.21 on the *imputable exon*
population does not carry the whole slot set — the gap between "works on the population it was measured
on" and "works on the panel", which is what `TRAPS: panel-before-src` exists for.

⭐ **`struct` is the arm that is safe everywhere**: 0.381 / **0.659** / 0.363 / 0.800, no stratum worse
than base, and it leaves the zero control untouched (1.000) because it never claims anything about exons.

## 49. ⛔⛔ AND READING `struct` ALOUD FOUND A DEFECT IN ITS CLAMP

`struct` sets ψ's shipped ½ wherever mature RNA CAN be contiguous, and asserts near-pure gDNA on the
**33,347 of 70,176 slots (47.5 %)** where it cannot: intergenic and intron REGIONs, `exon|intron`,
`intron|intron`, gene-edge and opposite-strand `exon|exon` BOUNDARIES. The predicate is the solver's own
`mrna_active`, and NO density estimate enters that branch at all.

⛔ The first version clamped with a flat `eps = 1/Σeff_g = 1e-8`, i.e. `m = 1 − 1e-8`. Measured, **only
0.94 % of the reference's mass then lies inside the shipped `L = 10` window** — worse than the 57 % the
refused `(a, b)` route leaves outside, and the answer becomes a function of the grid. ⚠ The
`test_l_invariance` gate never caught it because it sweeps `m ∈ [0.01, 0.99]` and `struct` operates only
at the clamped end.

**The candidate repair** is the derivation's own floor applied AT THE OBJECT rather than over the genome:
`m_i = E[g]_i/(E[g]_i + 1)`, exactly one pseudo-fragment of RNA *here*. It spans p5 **0.788** / median
**0.917** / p95 **0.998** and keeps **0.909–0.990** of the mass inside the window. ⚠ A pooled floor
(`1/Σeff_r`) does NOT fix it — 0.0404 inside.

## 50. ⛔⛔⛔ AND THE REPAIR IS MEASURABLY WORSE. THE PANEL REFUSED IT.

`ref_loc=struct_soft` against `ref_loc=struct`, final Σ|Δ| per stratum, ratio to base:

| stratum | ⭐ `struct` (flat ε) | `struct_soft` (derived floor) |
|---|---|---|
| stranded × capture OFF ⭐IN SCOPE | **0.381** | 0.609 |
| stranded × capture ON ⭐IN SCOPE | **0.659** | 1.045 |
| unstranded × capture OFF ⭐IN SCOPE | **0.363** | 0.580 |
| unstranded × capture ON ⛔DEFERRED | **0.800** | 1.000 |

⭐⭐ **The reason is the finding this whole investigation started from.** On a structurally pure-gDNA
object the truth IS `f_g = 1`, and objects at true `f_g ≈ 1` carry **49–83 %** of in-scope error precisely
because ψ holds them off that vertex. A near-improper prior pushes them onto it; the soft floor (median
`m` 0.917) pulls them back off. The theoretically cleaner object is the empirically worse one.

⛔ **The `L`-dependence is still real and this panel cannot see it** — it would surface the moment `L`
changed. So neither floor is settled, and picking one by its panel number alone would be tuning.
⭐⭐⭐ **What is actually un-chosen is the reference's STRENGTH on those slots — the second Beta degree of
freedom, which every measurement in this project has held at `a + b = 1`.** Both floors are kept as named
arms so the choice stays a measurement. ⚠ `pooled`/`local` carry the SOFT floor, which is what their rows
were measured with; with it they read 0.549 / ⛔ 5.789 / 0.516 / 1.005 and **0.001 at the zero control**.

## 51–58 — ⭐ MOVED OUT (2026-08-16). Nothing here; go to the permanent homes.

⛔ **A dev doc must never become the state**, so these eight sections were MOVED rather than copied the
moment they settled. What was in them and where it now lives:

| what | where it went |
|---|---|
| the structural classes, the precise `mrna_active` statement, the strength rule, and the two rulings that came out of SHIPPING it (τ_λ is the data's, the neutral location is an exact zero) | **`DESIGN.md` §6b.1** |
| the Beta-mean derivation, the strength closed form `log[(m²+(1−m)²)/(2m(1−m))] = L − log 2`, the fragment budget, and the `[f(1−f)]²` Jacobian that refuses a prior's curvature in τ_λ | **`EQUATIONS.md` §9c.1 / §9c.2** |
| the three lessons | **`TRAPS.md` §D** — `a-priors-curvature-is-not-the-datas-information`, `a-four-decimal-print-is-not-a-zero`, `a-constant-in-exact-arithmetic-is-not-constant-in-float64` |
| the panel numbers, the thermometer, and ⛔ **the nascent finding that keeps the flag OFF** | **`ROADMAP.md` §0** (a state row) and **§1 rank 2** (the work it points at) |
| the two mechanisms refused along the way — the τ_λ contribution, and the per-object soft floor | **`ROADMAP.md` §4.2** rows 6 and 7 |
| every claim, executable | **`tests/calibration/test_structural_reference.py`** — 14 gates, each with its own perturbation |

⚠ **The one thing that is NOT settled and is deliberately left as the open item:** the ladder holds
`nrna = 0` on all 16 rows, so it cannot falsify the structural claim. Building an `nrna > 0` rung is the
first step of `ROADMAP.md` §1 rank 2, not an afterthought to it.
