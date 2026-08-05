# TRAPS — mistakes this project has already made

**Read this before designing anything.** Every entry is a real defect that shipped, or a real
measurement that was believed and was wrong. They are here because each one cost time, and several were
made twice.

⛔ **This file holds LESSONS, not measurements.** A number appears only where it *is* the lesson (an
"exact factor of 2", a "12× over-represented"). Current measured state lives in `ROADMAP.md`; anything
here that reads like a status report is a bug in this file.

---

## A. Validation and gates

**A1. A validator that calls the builder's own helper validates nothing.** Deleting a negative-strand
coordinate swap left a 1,289-test suite green and the graph validator accepting; only re-deriving the
result by a *different algorithm* caught it. Worse than no check, because it reads as one.
*Corollaries:* a validator that only compares non-empty classes never sees a spurious flag — emit all
classes always; and prove every validator FIRES by deliberately corrupting its input.

**A2. A falsification test needs PERTURBATION, not just to be written first.** Writing the gate before
the fix is half the discipline; the other half is breaking the fixed code and watching each gate fire.
This has found holes in already-green gates repeatedly — once 7 of 9. Two recent examples, both blind
until perturbed: a capture A/B whose two arms did not share their random input (so it "passed" with the
feature removed), and a de-tilt test whose pool had **one length bin**, which is invariant under any
divisor whatsoever.

**A3. A gate that already passes is not a falsification.** A simulator gate on "is on-target gDNA longer
than off-target?" passed *with the defect present*, because the engine's per-fragment **conditional** was
correct and only the **marginal** was being discarded. Before trusting a gate set, check which of them
would have failed before the fix. If none, the gate has not been written yet.

**A4. When a per-observation conditional is correct, no conditional comparison can detect a missing
marginal.** The general form of A3, and worth stating separately because it tells you *where* to look:
if `P(outcome | x)` is right for every `x` but the distribution of `x` is wrong, every conditional check
passes and only a check on the marginal fails.

**A5. A bit-identity gate has lied in both directions.** An arm with ZERO rows scored "32/32 IDENTICAL"
because the comparison looped over the new arm's rows; and a stored baseline went stale so unmodified
HEAD no longer reproduced it. **Re-record the baseline from the current tree in the same session**; if
HEAD-vs-baseline is not 100 %, the baseline is what is broken.

**A6. A `min()` clip hid an exact factor of 2 for months.** Pooled seam density read 1.994× truth
because the code SUMMED two faces' mass and divided by their AVERAGE. It survived 29 tests: all four
fixtures stored an un-halved length AND deposited half the mass (cancelling exactly), and `min(2S,S)=S`
under a uniform field, so the "contraction is 1 under a uniform field" bedrock invariant — written to
catch exactly this — passed with the bug present. *Lessons:* a clip can hide a scale error from a bedrock
test; repair fixtures, never relax assertions; **an exact algebraic 2 is never a modelling
approximation.**

**A10. AN ABLATION THAT NEVER RAN READS AS "NO EFFECT", AND TWO IMPORT HABITS MAKE THAT EASY.** An arm
that deleted the whole junction channel scored **byte-identical to the baseline on every class of every
condition** — a clean, publishable "this channel is inert". It was inert only because the monkeypatch
missed the module that actually binds the name: `calibrate.py` does
`from .bp_solver import build_node_geometry`, so its module global is a separate reference.
⛔ **And the first repair missed it too**, because `from rigel.calibration import calibrate` returns the
re-exported **function** of that name, which shadows the submodule — so every `hasattr` on it is False
and the patch silently skipped it. The submodule must come from `sys.modules`.
⭐ *The tell that caught it:* the ablation was a strict **superset** of another one, and scored a
smaller effect. ⭐ *The rule:* every ablation increments a counter and the harness **raises** if it did
not fire. A5's "bit-identity has lied in both directions" with the lie on the other side.

**A11. A TEST THAT RE-DERIVES A DEFINITION CANNOT DETECT DRIFT IN IT.** A gate existed precisely to keep
two instruments' shared class definition from diverging — and it recomputed that definition inline in
the test. Changing one instrument fired nothing. It is A1's shape once more: the check and the thing
checked have to come from *different* places, and for a shared definition that means **one home**, with
every consumer importing it. Here the home is production code (`node_geometry.g1_locked`), because the
predicate is a production concept and `scripts/` is deliberately not importable.

**A7. Prove the SUBSTRATE before proving the code.** For two milestones the simulated panel's
post-capture fragment-length distribution was byte-identical to its pre-capture one, and everything
measured against it inherited that. The tell was free the whole time: diff the two capture arms' truth
files. **When a simulated axis is the axis you are judging, gate the simulator on it.**

**A8. Before running a benchmark, prove it can resolve the axis you are changing.** A 32-condition suite
was used for months to judge a partition change it was structurally incapable of seeing: its fine node
set was row-for-row identical to its merged set. It also had `frag_std = 0`, was Poisson by construction,
and over-represented the terminus+splice-site seam **12×**. `scripts/design/suite_resolves.py` is this
lesson made executable — every requirement scored against its *degenerate* value, no tuned thresholds.

**A9. A toy ranks performance hotspots BACKWARDS.** At 3.4 k vs 1.5 M nodes: message relay 34 % → 81 %,
per-node solves 9 % → 29 %, the prior's EM **28 % → 0.7 %**. A whole analysis was spent on the toy's #1
hotspot. Profile on cached real data.

---

## B. Measurement and inference

**B1. MEASURE THE CEILING BEFORE BUILDING THE CORRECTION.** Hand the consumer the *exact* answer for one
channel and see what perfecting it is worth. One channel was rank 0 for two sessions on the strength of a
coupling constant; the ceiling showed a **perfect** model of it bought ~1 %, while the channel nobody was
ranking was worth 21 %. An A/B tells you whether a change helped; a ceiling tells you whether the work is
worth starting — and it is available whenever the simulator writes truth.
*Instrument:* `scripts/design/calibration_truth_ab.py --ceiling`.

**B2. Score against TRUTH, not against the previous run.** The simulator writes per-fragment ground truth
into the oracle BAM's read names. That is what found a score-annihilation bug, and what turned "the
deliverable improved" into "the deliverable got 23.9 % worse".

**B3. A zero-target guard is ONE-SIDED.** In a library with no gDNA, *any* change that lowers the
estimated gDNA fraction scores better. This reversed a published verdict once — a reported −13.1 % win
became +1.4 % worse over the full battery. Score on the contaminated conditions; use zero-target rows as
false-positive checks only.

**B4. Hard-label metrics are nearly blind to soft changes.** A real change can move soft pool counts by
tens of thousands of fragments while a hard-label net is byte-identical. **Treat a byte-identical
hard-label result as no evidence, not as no change.** Related: never score pooled — one run had a −47.6 k
error and a +75.4 k error reading as a "nearly perfect" −1.9 k.
⭐ **And the cancellation is not bounded.** Scored per object against the oracle, `Σ|err| / |net|` reached
**274×** on one arm: the library-level figure read a +672 fragment error over a per-object Σ|err| of
184,271. **A net error is a lower bound on the real one, and it can be an arbitrarily weak one.** Report
`Σ|err|` beside every net, and the directional split beside both.

**B10. THE DEFAULT INSTINCT IS A POOLED AVERAGE, AND IT IS WRONG HERE THREE DIFFERENT WAYS.** B3, B4
and B5 are each one instance; the pattern is worth naming because it recurred three times in a single
session, each time in a NEW instrument, each time caught by review rather than by the tool:
*(i)* averaging an error over objects that **cannot be solved**, so honest ignorance reads as error and
buries the real answer (0.0456 reported as 0.3150); *(ii)* averaging across **saturated zero-target
rows**, where any change that lowers the estimate "improves" the row; *(iii)* averaging a rate over
objects whose **mass differs by four orders of magnitude**. ⭐ The habit that prevents all three:
**report the partition BEFORE the aggregate, and make every aggregate state which population it is
over.** An instrument whose headline is a single number is one refactor away from pooling the wrong
things again.

**B11. A BINARY CUT ON A FITTED PARAMETER'S RESIDUE IS NOT A POPULATION TEST — AND A BETTER THRESHOLD
IS NOT THE FIX.** B10 (i) was corrected by excluding objects with *no* own evidence, tested as
`tau_lam > 1e-9`. But the strand arm carries `I(f_g) ∝ (2κ−1)²`, exactly zero only at κ = **½**, and κ
is *fitted*: at 10 M fragments κ̂ = 0.500689 on a genuinely unstranded library, so τ lands at ~5e-7
instead of 0 and the cut promotes the object to "solvable" while its own statement has
**sd(λ) = 1,377 nats against a solver that represents λ only on ±10**.
⭐ **Nothing upstream is wrong** — the derived deadband `disc = 4·max(0, (κ−½)² − σ²_d)` correctly
*admits* it, because κ̂−½ is a real **1.9σ** fluctuation at that depth. The information is
statistically genuine and physically nil. The defect is the **binary cut**, which cannot tell those
apart.
⛔ **And the obvious repair was implemented and REFUTED by its own gate:** a resolving-power floor at
`1/(2L)²` requires an empty interval to sit in, and τ is **continuous** across that region on 4 of 5
ladder conditions (only unstranded × capture-OFF is bimodal, and there the two clusters are the silent
strand arm at ~1e-7 and the live intron factory at ~1e-1). With no gap, any floor is a tuned constant.
⭐ *What replaced it:* report `sd(λ) = 1/√τ` as a **curve over decades**, the same device the
instrument already used for confidence. It needs no cut, and it says the thing directly — on
`g75 ss0.50 capture_off`, **79.1 % of the "solvable" error sits on objects with sd(λ) ≥ 10 nats.**
⚠ Its second consequence: two strata that were reported as different *regimes* (unstranded
capture-OFF "2.0–96.2 % solvable" vs capture-ON "0.1–6.3 %") differ only by which side of 1e-9 κ̂
landed on — 0.500689 against 0.499972.

**B13. ⛔⛔ EXCLUDING A POPULATION FROM THE DENOMINATOR WITHOUT GATING ITS OWN FAILURE MODE HIDES THE
LARGEST ERROR IN THE LIBRARY.** B10 (i) is right — averaging error over objects that cannot be solved
buries the answer, and the undetermined class was correctly excluded. But an exclusion is a promise to
check that class *some other way*, and `SUCCESS.md` even named the check ("its only failure mode is the
opposite one: claiming a precision it has not earned") — and nothing implemented it. So on
`gdna_g25_ss_0.50_nrna_none_capture_off`, whose reported score is **mwae 0.0170** on a 39.5 % solvable
set, the excluded class carried **1,056,019 fragments** of error, with **6,425 objects moved a mean
0.414 away from ½ and 100 % of them declaring a finite precision**. The single worst mechanism — 87 deep
exon nodes driven to `f_g = 0.829` against a truth of 0.009, 395,251 fragments — was **0.0 % scored on
every rung but one**, and the one rung where it *was* scored is the row the roadmap had already flagged
as inexplicably anomalous. ⭐ **The rule: every population you exclude needs its own gate, written at the
same time as the exclusion.** An object with no evidence reporting ½ at zero precision is correct; the
same object reporting 0.83 at finite precision is the relay inventing an answer, and it is invisible to
a metric that does not look. ⚠ And the tell was available: the *reason* to exclude a class is that its
answer should be uninformative, so **measure how far from uninformative it actually is.**

**B14. ⭐⭐ IF EVERY GATE FOR A CHANGE READS ONE INSTRUMENT, THE CHANGE IS GATED IN ONE PLACE — AND
`_relay`/`_transport` ARE TWO.** The population licence landed in both twins and the gate set for it was
ten tests, all sound, all reading `_capture['_pin']` — which is the COMBINE's publication. Perturbation:
delete the conjunct from `_relay` alone and **the entire calibration suite passes**. The relay's own output
is observable (`fwd_g`, the running level, before the combine touches it), so the missing gate was cheap
once the hole was known. ⭐ **The generalisable check: for each place the change was made, name the
observable that would move, and confirm at least one gate reads it.** A count of gates says nothing —
`test_gdna_scale_rule`'s gate 5 exists for exactly this reason and its docstring says so, and the hole was
still re-opened by the next change to the same pair. ⚠ And the readable-destination problem is why it is
easy to skip: an exact statement about a relayed level needs a destination with no own precision, which on
a stranded chain is only the pure-gDNA EDGEs — an **AMBIG** node is the one slot that is precision-free
without being structurally certain, and that is what makes the gate expressible.

**B15. ⛔⛔ "STARVED" AND "DEPLETED" ARE NOT THE SAME DIAGNOSIS, AND CALLING ONE THE OTHER THREW AWAY HALF A
PANEL.** Under capture the toy's intron and intergenic NODEs fall to ~1 count and I reported the whole
capture-ON stratum as an empty chromosome, concluding the bench could not measure it *at any chromosome
length*. ⛔ Both halves were wrong. Capture **moves** the gDNA signal to the EDGEs abutting the exon — that
is what capture *is* — and the toy's gDNA budget is `rate × genome_length` against a **fixed** probe
footprint, so the sampler's on-probe share `binding·overlap/(off_target·L + binding·overlap)` means a
longer chromosome hands capture more budget to concentrate on the same probes. Measured 12 kb → 120 kb,
same donor and chemistry: the `intron|exon` EDGEs go **2 → 20** and **5 → 36** counts and the gene-boundary
EDGEs **1 → 41** and **2 → 35**, while the intron NODE stays at 1 — exactly the redistribution, not a
depth problem. ⭐ **The check that would have caught it: before declaring an object unmeasurable, ask which
term in the sampler's own weight law it depends on.** An edge's count is `density × mean_FL` capture-OFF
(so `L` is inert, correctly) and a share of a `L`-scaling budget capture-ON. One sentence of the retention
law separates the two, and `docs/TESTING.md` §0b carried the capture-OFF answer for both cases for weeks.

**B16. A TILING LOOP CAN SUPPRESS THE VERY POPULATION UNDER STUDY.** Toy capture probes are written in
TRANSCRIPT space, so a probe spanning an internal junction offset has a genomic footprint in **two blocks**
— and `sim/capture/sampler._split_scale` then multiplies every gDNA fragment overlapping it by
`gdna_split_penalty`. Tiling across the whole transcript therefore put a split probe over every internal
junction and depressed exactly the fragments that span an `intron|exon` EDGE. ⭐ The lesson is not "tile per
exon" (though that is the fix); it is that **a substrate knob can be adversarial to the object you are
measuring, and the way to find out is to read the sampler's weight for that object's own population** —
here, one `if len(blocks) > 1` in code nobody had reason to open.

**B17. A CEILING BY SUBSTITUTION UNDERSTATES A MESSAGE SOURCE.** Replacing one object's answer with the
truth and re-scoring is the honest ceiling for a SINK — but a source's value is what it *carries*, and a
substitution does not propagate. Measured: substituting both `intron|exon` EDGEs removed 9.1 % of the
gene's error, while the object they feed accounted for 82.7 %. ⭐ For a source the ceiling must be a real
arm — pin it and RE-SOLVE — and an instrument that offers substitution must say which of the two it is
doing (`toy_panel.py` now does).

**B12. A STRUCTURALLY LOCKED OBJECT IS NOT A CONTROL.** An `intergenic|exon` seam predicts
`f_g = 1.0000` exactly against a truth of 1.0, and that was read as the healthy twin of the broken
`intron|exon` seam beside it — "the same object structurally, and the only difference is the junction
flux". It is not: it is **G1**, so it is never solved at all and keeps its pinned `{0,0,1}` init. It
cannot be wrong, so its correctness measures nothing. ⭐ *The control that works is the same class
split on the variable* — `intron|exon` lines **with** vs **without** junction flux — and it reversed
the conclusion: the class is 23 % low where **no** junction attaches, so the flux explains 26 % of its
error unstranded and **2.5–4 % on the strata where the class has own evidence at all.**

**B9. A per-fragment-independent partition stops being one the moment a downstream step conditions on
the whole tally.** The oracle is trustworthy because the accumulator deposits each fragment
independently, so splitting the BAM by origin and re-scanning must reconstruct the full payload exactly.
The second pass breaks that: its multinomial is scored against the payload's *own* densities, so three
partitions drained separately are not the whole drained. **The identity that makes a truth source valid
also fixes where in the pipeline it can be taken** — and the honest response is to measure the undrained
stage rather than to drain the parts and hope.

**B5. A bp-weighted mean and a fragment-weighted mean answer different questions.** An "11 % over-call"
was a bp-weighted geometric mean; the estimator is fragment-weighted, 89 % of the mass sat where the
effect is inert, and end-to-end it moved the answer by ≤ 0.0002. Both statements were true; only one
decided anything. **Weight the average the way the consumer weights it.**

**B6. A SUPPORT CEILING THAT MATCHES THE CLAMP IS NOT A MATCH — it is the clamp.** A distribution's
support "agreeing" with `max_frag_length` was recorded as a fix; the narrower estimate had been correct
and the "fix" was an uncut intron.

**B7. Every "is the declared precision earned?" number written before 2026-07-28 compared a LOG-space
variance against a LINEAR squared error.** Corrected, one suite total went 0.046 → 1.007 and the
per-class ranking INVERTED.
⛔ **AND THE TRAP IS STILL ARMED, BECAUSE THE FIELD IS NAMED FOR THE WRONG SPACE.**
``NodeBelief.var_gdna`` / ``NodeDeconv.gdna_frac_var`` is ``Var(log f_g)`` — the grid moment of
``_log_fg(lam)``, stated in `simplex_logodds` and required by D2 — while both names read as "the
variance of the fraction". So the obvious way to write *confidently wrong*,
``|f_g − truth| / sqrt(var_gdna)``, silently re-commits B7. ⭐ The standardised discrepancy must be
``(log f_g − log f_truth) / sqrt(var_gdna)``, and where the truth is 0 or 1 (a pure-RNA or
structurally-pure-gDNA object, both common) the clip is **the solver's own λ-grid endpoints** — not a
chosen epsilon. The solver cannot express a fraction outside its grid, so clipping truth to that
support is the only comparison it could ever have won.

**B8. A delta is only attributable if its baseline came from the same tree in the same session.** Re-record
the before-picture, do not quote a stored one.

---

## C. Pools, selections and divisors

**C1. A PURITY FILTER ON A LENGTH POOL IS A LENGTH FILTER.** Barring fragments whose length was partly
*inferred* rather than sequenced selects exactly the ones whose mates sit far apart, so the pool became
**−9.6 % mean / −22.5 % sd** against truth where keeping them read +0.7 % / +2.4 %. **Before excluding a
population from a pool, ask what the exclusion criterion correlates with.** If it correlates with the
axis being measured, purity and accuracy point in opposite directions.
*The rule that replaced it:* a pool is keyed on **determinacy, not provenance** — a fragment enters when
exactly one hypothesis survived, however it got there.

**C2. A POOL CAN BE COMPOSITION-PURE AND LENGTH-CENSORED AT THE SAME TIME, and the second is invisible to
every purity argument.** "gDNA contained in an intergenic or intronic node" is 100 % gDNA by construction
and, under hybrid capture, ~15 % short — because a long fragment beside a probe *reaches* the exon
boundary and stops being *contained*. It reads as contamination: the pools it spills into resemble the RNA
pool and were filed for two milestones as "not gDNA". They are gDNA. **Ask what a pool's selection rule
correlates with, not only what the selected fragments are.**

**C3. A POOL DIVIDED BY ITS OPPORTUNITY MUST BE DIVIDED BY A PROBABILITY, NOT BY A COUNT.**
`count(w)/A(w)` recovers the distribution lengths were **drawn** from; every consumer needs the one the
library **realises**, which is the drawn one weighted by how many placements each length has. So the
divisor is `A(w)/T(w)`. The two forms differ in *shape*, and the ratio is also where an abundance
weighting cancels — swept to pathological regimes, the ratio form is never worse than not correcting and
the `A`-only form is.

**C4. POOLS WITH OPPOSITE TILTS MUST NOT BE POOLED RAW.** A *contained* pool's opportunity falls with
length; a *crossing* pool's rises. Summing the histograms and applying one divisor read a gDNA mean of
**146.05** where the contained pool alone said **88.0**. ⭐ Summing the counts **and** the matching
per-pool opportunities is a *different* operation and is correct — it is the opportunity-weighted average
of the per-pool estimates, and under Poisson counts those weights are exactly inverse-variance.

**C5. Fractional mass IS the partitioning problem.** A fragment spanning 4 nodes writes six fractional
numbers whose values depend on region sizes, purely because a mass is conserved; the same fragment's
three crossing counts depend on nothing. *Corollary:* multimapper and ambiguous-path assignment must stay
INTEGRAL, or the non-integer observable returns and the count stops being a count.

**C6. Mass conservation does not catch mis-attribution.** One fragment can credit the same boundary side
twice, and credit a boundary it never crossed, with total mass still exactly 1.0.

---

## D. Estimation and solver design

**D1. You cannot fix a biased mode with a variance.** Established three times independently. Under
capture a counting estimate was systematically ~2× low but PRECISE — both flanking seams sat at the same
enriched edge and agreed on the same biased-low density, so the bias was trusted. **A disagreement-based
variance model structurally cannot fix a bias.**

**D2. Never hand a solver two Gaussians built from one latent.** A message on `log f` and one on
`log(1−f)` are rank-1 with correlation exactly −1, so adding their Fisher information is exactly **2×**
over-confident, rising to ~7× with deep spliced content.

**D3. Never fit a variance on the current, not-yet-solved belief.** Adjacent WRONG nodes agree, so the
variance collapses, the messages turn confident, and the error propagates. *Honesty measured against a
wrong belief is not honesty.* Same family: **any component trained on the solver's own output is
self-confirming** — refit iterations 1→5 went 0.0840 → 0.1056, monotonically worse.

**D4. A message computed from the destination's own belief carries zero information and confirms the
destination.** Hit twice: one delivered exactly `1/(1+f_own)`, reserving 33.6 % of the budget for
imaginary gDNA so a zero-gDNA library read back 29.3 %; and at a gene end the rescale ratio became
`1/f(dst)`, making the delivered density the destination's own total (7× too big in the median, up to
190×). **A message may use the destination's CONSTANTS — geometry, lengths — never its BELIEFS.** Any
"fix" that divides the destination's belief back out rebuilds the bug.

**D4b. D4 CAME BACK, AT A GENE END, AT 1,118× — AND A "TOTAL DENSITY" RATIO IS HOW.** D4 says a message
computed from the destination's own belief carries zero information and confirms the destination. The
reframe `r = ρ_tot(dst)/ρ_tot(src)` re-creates it whenever `ρ_tot(dst)` is dominated by a component the
message is not about: at an `intergenic|exon edge → EXON` hop the edge's *correct* gDNA density
(0.0257 counts/bp) is multiplied by the exon's total density ratio and delivered as **28.684**, which is
`≈ M/E_g` — **the exon's own total density**. Because the inflation is itself proportional to the
destination's density, the delivered gDNA level *always* equals the destination's total, so the answer is
**flat at `f_g ≈ 0.90` across a 10,000× sweep of RNA density** while the truth spans 0.914 → 0.001. ⭐ **A
prediction that does not move when the data moves by four orders of magnitude is the tell**, and a sweep
on a five-object toy exposes it in 0.1 s where a 36-condition panel had hidden it for weeks.
*The rule:* a scale factor must be built from the component the claim is **about**, never from a total —
and D4's check is "is the delivered value independent of the destination's own state?", which this
failed.

**D4c. AND THE FIX WAS NOT A BETTER SCALE — IT WAS NOTICING THE REFRAME HAS NO LEVEL IN IT.**
`ρ_c(src)·ρ_tot(dst)/ρ_tot(src) ≡ φ_c(src)·ρ_tot(dst)`: the reframe is a **pure composition imputation**,
the source's share applied to the destination's total, with no level transport whatever. So there is no
correction factor that repairs it — the missing factor is `φ_c(dst)`, the destination's own belief, and
multiplying it back in *is* D4. ⭐ **The fix is a LICENCE, not a scale:** the reframe is allowed exactly
where a composition imputation is allowed (the source SUPPLIED both components — the λ-emission gate's
own predicate, which was already derived and was being applied to the τ stream only), and where it is
not, the gDNA **level crosses unscaled** because gDNA is uniform before capture.
⭐⭐ *The lesson that generalises past this bug:* **before correcting an operator, substitute its own
definitions and read what it delivers.** Two sessions were spent looking for a better `r` — capture
efficiencies, class landscapes, probe geometry — for an expression that never carried a level. `EQUATIONS.md`
§3.5. ⚠ Two things measured on the way, both worth not repeating: a per-capture-class gDNA landscape ratio
is **byte-identical off capture** (the mass pin already carries the landscape through each pure-gDNA
object's own measurement), and the tempting affine extrapolation `e_g[exon] = 2·e_g[crossing] − e_g[off]`
is **exactly the simulator's own retention law** (`sim/capture/sampler.py`: `off_target + gain·overlap`), so
fitting it would score against the substrate that generated it.

**D4d. THE THIRD D4 WAS IN THE OPERATOR BUILT TO DEFEND AGAINST THE SECOND, AND ITS TELL WAS A FIXED
POINT AT ½.** The relay's mass pin restores `Σ_c ρ_c·E_c = M` with `k = M/S`, filling every component the
message did not supply from the destination's OWN density. Substituting its own definitions (D4c's lesson,
applied again) gives `k = 1/(φ_msg + R_own)` — a saturating map with fixed point `(1−R_own)·ρ_tot`, and
`R_own`, the RNA share of the destination's own self-solve, is **exactly ½** at any object with no
composition evidence. So it drove the delivered gDNA FRACTION to ½ regardless of the truth and regardless
of what the source measured; at `R_own = 0` it collapsed to the destination's own total with the incoming
level cancelled algebraically. ⭐ **The fix is the same licence, and the derivation is one sentence: the
identity holds only under the imputation premise, so where the premise is withheld there is nothing to
restore.** ⚠⚠ **AND IT HID FROM EVERY AGGREGATE**: the running product of the per-step rescales telescopes
back to exactly 1 at the far end of a gene, and the last step into a pure-gDNA object rewrites the level to
that object's own total anyway. No conservation, endpoint or aggregate check could see it — it is visible
only per object, and only away from a pure-gDNA object. `EQUATIONS.md` §3.5c.

**D4e. ⭐⭐ "IT USES THE DESTINATION'S OWN NUMBER" IS NOT THE TEST — D4 IS ABOUT BELIEFS, AND AN OBSERVATION
IS NOT ONE.** Gating the pin off wherever it read the destination's own density looked like the clean fix
and it **broke the capture landscape**: the off-probe intergenic floor leaked through every gene-boundary
EDGE into the exon behind it (measured, `test_gdna_scale_rule`'s capture gate at 20× and 200×). At a
*structurally pure-gDNA* object there is only one component, so `f_g = 1` is STRUCTURE and `M/E_g` is a
direct **observation** of the very quantity the message is about — and it is a better one than anything
relayed, because a G1 EDGE measures the gDNA density at its own capture stratum while carrying `prec_g = 0`
and so having no other channel (`ISSUES.md` C6). ⭐ **The general form: state the licence as "no BELIEF may
enter", not as "the destination's numbers may not enter"** — the first is D4 and admits the structural
case; the second is a superstition that costs a real mechanism. Two states satisfy it, not one.

**D4f. ⭐⭐ B11 IS NO LONGER ONLY A MEASUREMENT DEFECT — IT NOW GRANTS A LICENCE INSIDE THE SOLVER.** B11
records that κ̂ = 0.500689 on a genuinely unstranded library leaves `I(f_g) ∝ (2κ−1)²` at ~5e-7 instead of
0, and concludes that the damage is to the *classification*. It is not only that any more: the composition
licence asks whether the source SUPPLIED a component, "supplied" being a statement about **precision** with
no threshold — correctly, by design — so those same 1e-6…1e-4 precisions now **grant the imputation
licence** at every evidence-free exon. Measured: on the `nested_exons` toy the licensed mass pin still
fires through the whole gene interior, so licensing recovers 0.2618 → **0.2264** where switching the pin off
entirely reaches **0.0760**. ⭐ **The lesson is not "add a floor"** — B11 already refuted that, and it would
be a tuned constant. It is that a **structural** zero (κ = ½ ⇒ no strand information, exactly) is
represented by a point estimate whose own variance covers it, so the honest repair is to propagate
`Var(κ̂)` into the strand arm's precision, where it drives that arm to ~0 by derivation. Same family as
`ISSUES.md` C2/C7.

**D5. "No prior" does not exist on a grid — omitting a term lets the grid supply Haldane**
(`p(x) ∝ 1/x`, improper, an amplifier toward the vertices). Posterior median spread over grid half-widths
4–20 is **0.045** at Jeffreys and **0.525** at Haldane.

**D6. Sums are well conditioned; differences are not.** Subtracting across a junction gives
`Var(log ρ) = u²σ_T² + (u−1)²σ_μ²` with `u = 1/(continuing share)` — at the real median u = 2.3 already at
the edge of validity, at p75 u = 5.3 hopeless at any depth. **Prefer shares.**

**D7. An all-zero factor is uninformative, not decisive.** In a multiplicative score, a factor that is
zero for every candidate used to annihilate the other factors and collapse the record to a coin toss.
Skip a flat-zero factor; do not multiply by it.

**D8. Density below one fragment length is not resolvable by ANY design.** A 1 bp node has no
independently measurable density and never will (*composition* still does, since it depends on what the
fragments are, not where). *Corollaries:* an object with zero opportunity for a component must emit
nothing at zero precision, not a floored division; and "no data" must be inert, never "100 % gDNA" — that
default was actively seeding false gDNA into neighbouring exons.

**D9. An identical-paralog split is bimodal, and depth does not fix it.** Two sequence-identical
transcripts split either evenly or all-or-nothing, and which one is a coin flip on the fragment draw:
same code and seeds, `n_fragments` 3,000 → 171/0, 6,000 → 249/237, 12,000 → 679/0, 20,000 → 773/795. The
*total* is right in every case. ⛔ **Do not "fix" this by moving a seed or a depth until it lands even** —
that is tuning to green. It is an unidentifiability, and the honest response is either to damp the
degenerate direction or to score only what is identifiable.

---

## E. Structure, indexes and plumbing

**E1. Single-reference synthetic indexes hide reference-id-space mismatches.** A resolver assigned ref-ids
by first-seen interval order rather than the index's own order, so **476,719 of 476,732** real fragments
were silently dropped inside `deposit()` while every golden test passed. *And the fix is per-path:* extra
references make the space non-trivial for the paths those references exercise. RNA-only spike-ins do
nothing for the gDNA deposit path, which is a different branch — that needs **≥ 2 genomic references**.

**E2. "HAS AN ANNOTATION" IS NOT "IS GENOMIC" — a classification must be an INPUT, not a proxy.** The
simulator chose its gDNA references with `{t.ref for t in transcripts}`. Every RNA-only spike-in carries
exactly one transcript, so every one qualified, and the panel filled with gDNA on synthetic RNA
references — on the very templates whose truth abundance of zero is the false-positive control. **And a
mis-stated classification must raise, not silently produce nothing.**

**E3. A splice junction CANNOT be detected from a gap between deposited slices.** A contiguous spliced
read whose exon body straddles an internal cut has no gap at all. The junction's identity is the cut-intron
coordinates, which the scanner already has — pass them through.

**E4. A splice deposit belongs at the JUNCTION's coordinate, not the node's edge.** Invisible for
annotated introns (their ends are cuts); for unannotated junctions the mass lands kilobases away.

**E5. Alternative splicing makes a node↔junction graph CYCLIC, and cycles are the common case** — a
cassette exon is a 4-cycle, and the human graph has ~404,000 independent cycles, one per junction.
Two-sweep forward-backward is exact only on a tree. **Never break a cycle by dropping a junction edge** —
that re-isolates the exon the edge exists for.

**E6. `~is_synthetic & ~is_nrna` as the "real transcript" filter deleted 52,104 real termini.** On a
non-synthetic row `is_nrna` means "single-exon, so mature ≡ nascent", NOT "manufactured span". **ONE
filter: `~is_synthetic`.**

**E7. A fragment crossing K junctions must credit exactly ONE** (the leftmost annotated). Crediting all K
shifts the library sense fraction 21–34 % and creates between-side correlation that reads as
overdispersion (a zero-overdispersion simulator fit 0.092). A `1/K` split is provably biased by 4–12 σ.

**E8. `(src, kind, dst)` is not a total order for junction edges** — two strand-coincident junctions differ
only in strand, so ordering becomes input-order-dependent and the duplicate check reads them as duplicates.
GENCODE has zero of them, so only a synthetic stress test finds it. Sort on `(src, kind, dst, strand)`.

**E9. Cache keys that do not cover the artifact they cache.** A partition hash covered only the node file;
a flag fix rewrote every edge file while leaving every node file byte-identical, so a stale cache would
verify CLEAN. **Never store a derived hash beside the data it describes; compute it on demand.**

**E10. Integer channels are bit-identical across worker counts; float channels are not.** Max relative
3.7e-7 per cell, which propagated to a ~2.6 % difference in the calibration output. Integer addition is
associative — that is the whole fix. *Corollary:* the one bank whose ORDER is observable (a list, not a
sum) must be sorted on its own content before it crosses the ABI.

**E11. Worktrees silently run the wrong code** — an editable install's meta-path finder beats
`PYTHONPATH`, so an A/B inside a git worktree executes the main repo's source.

**E12. `git checkout -- <file>` does not undo a perturbation when the work is uncommitted — it deletes the
work.** A perturbation harness must restore from a copy of the WORKING TREE. Cost one full
re-implementation.

**E14. TWO DIFFERENT MASKS SHARED THE WORD `struct_lock`, AND BOTH WERE RIGHT.**
`node_init.strand_evidence`'s `struct_lock` is **node-only on purpose** — it governs whether a slot may
*emit* composition certainty into its messages, and a G1 edge is excluded because a seam is structurally
gDNA yet sits between RNA-carrying exons, so certainty there compounds into a phantom-gDNA emitter.
`_type_belief`'s G1 lock is **both axes** — it answers "is the belief pinned and certain?". Three
instruments classified objects with the node-only version while documenting the both-axes meaning, and
the mismatch quietly moved every structurally-locked **edge** out of the scored population and into
"honest ignorance". ⭐ *The fix is not to unify them* — they answer different questions — *it is to name
them differently and give the shared one a single home.* E13's lesson with the disagreement between two
**masks** rather than two docstrings.

**E13. The prose next to the code said "the AVERAGE" and the code followed the prose,** while a sibling
module's docstring had the correct formula the whole time. Two docstrings disagreed about one quantity for
months and nobody diffed them.

---

## F. Domain facts that read like defects

**F1. Strand specificity is TWO different quantities and they are complements.** A simulator's
`strand_specificity` is protocol **fidelity** (direction-agnostic — an R1↔R2 swap with probability
`1 − ss`); a fitted `sense fraction` is **directional**. For an R1-antisense (dUTP) protocol — which real
cfRNA is — comparing them reads as a sign error. A fitted κ of 0.0101 on a "0.99 stranded" library is
**correct**, and forcing 0.99 measured 166× worse. The matching quantity is
`StrandModel.strand_specificity`, which recovers the knob directly (1.00 → 1.0000, 0.75 → 0.7701,
0.50 → 0.5020).

**F2. Strand measures the TILT, not the gDNA fraction.** With RNA tilt `d = f₊ − f₋`,
`p = ½ + (κ−½)·d` — the gDNA fraction **cancels identically**. Strand reaches gDNA only through the
triangle bound `f_g ≤ 1 − |d|`: tight on a single-strand node, slack on a both-strand node. And
`I(f_g) = 0` **exactly** at κ = ½, for any count and any overdispersion.

**F3. At equal component mean lengths the length channel carries EXACTLY ZERO information about
composition, at any depth.** The 2×2 deconvolution is identified only through `μ_g − μ_r`. A claim that
one storage choice beats another was measured at a 4× mean separation and is **false** at every node
≥ 250 bp and reversed at equal means.

**F4. Hybrid capture is ~1000× on exons and only gDNA reads it cleanly.** RNA's own 10⁴ expression range
hides the probe pattern; gDNA's uniform baseline does not. And **capture destroys the intron signal 75×**,
so nascent-vs-gDNA is fundamentally unidentifiable under capture.

**F5. Capture SELECTS FOR LENGTH, and the post-capture distributions are the baseline.** Probes hybridise
to sequence, so a short fragment presents less sequence and is captured less efficiently. The pre-capture
parameters describe a library that was never sequenced — score against post-capture truth, for fragment
lengths exactly as for abundances. *And capture narrows the gDNA↔RNA length gap* whenever gDNA is the
shorter component, because the short tail it removes is disproportionately gDNA.

**F6. "On-target" defined by the START's territory is geometry, not capture efficiency.** Conditioned on
being captured, a fragment whose start is in the intron is one that was **long enough to reach the
probe** — weight ~`w²/2` — while an exonic start carries ~`p²/2`, flat in `w`. So an intronic-start
population reads *longer* than an exonic-start one under any capture model. The population that
physically binds is the one that **overlaps a probe**.

**F7. Real data is a TEST input, never a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
spectrum, not a sample of it. Sweep the plausible space, report the worst case, bring the domain call to
the owner. In particular: **never assume RNA fragments are longer than gDNA** — true for cfRNA, false
elsewhere.

**F8. "Effective lengths cancel, so a node's marginal is just its length" is FALSE near any transcript
end** — a mature fragment must fit in the remaining transcript; gDNA need not.

**F11. EQUAL CONFIGURED FRAGMENT LENGTHS DO NOT GIVE EQUAL REALISED ONES.** Handing the simulator
identical `frag_mean`/`frag_std`/`frag_min`/`frag_max` for gDNA and RNA still produced a **+4.60 bp
(+2.16 %)** realised gap at 206 ± 98 — *larger off capture than on it*, so capture is not the cause.
A mature fragment must fit inside its transcript and gDNA need not (§F8), so transcript-length
truncation pulls the RNA mean down while gDNA keeps the bare truncated-normal mean. ⚠ That residual is
the same order as the gap the length channel is meant to read, so "equal FL" as a way of switching the
length signal OFF has to be **measured on `truth_fragment_lengths.tsv`, never assumed from the config**.
*What works:* a fragment distribution short and narrow relative to transcript lengths — 150 ± 30 over
[100, 250] measured **+0.43 bp / −0.06 bp**, and the sign differs between capture arms, which is what
no-signal looks like. ⛔ *And it has a side effect that must be paid for:* at read length 100 the mates
then OVERLAP, the unsequenced gap vanishes, and the held fraction collapses **1.73 % → 0.04 %**, leaving
the side buffer and the drain untested. Shortening the read to 2×75 restores it to 0.78 %.

**F9. Mature RNA never crosses an exon↔intron seam** (0 of 1,146 seams, 7/7 conditions). Exon↔exon seams
do. This is the hard empirical case that a contiguous seam and a splice junction are physically different
objects — and it is what makes the two exon-crossing gDNA pools pure.

**F10. But a seam with RNA crossing it need not be a junction.** One position can be a splice donor for
transcript A and plain contiguous exon for transcript B; zero-gDNA libraries show seams with 44–55 k
unspliced fragments that are 100 % RNA.

---

## G. Process

**G1. No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every
divisor must be derived from the deposit rule and unit-tested against brute-force enumeration.

**G2. One thing varied per experiment**, with the falsification test written first and verified failing,
and a baseline re-recorded from the current tree in the same session.

**G3. Converge and delete.** No legacy, no backwards compatibility, no speculative code. Code kept "for
comparison with the old version" is a defect. No version suffixes in file names.

**G4. The source does not reference the docs.** Docs evolve and rot — 73 % of the citations that used to
be in the source pointed at documents that had already been deleted. A docstring may cite a **test** or
an executable specification, because code cannot rot silently; it may not cite a document.

**G5. Do not re-propose path or cell enumeration without a memory census.** Possible unspliced paths ≈
1 M nodes × 3–6 reachable ends at ~100 B each = 0.3–0.6 GB, plus spliced paths. It was killed by memory,
and no consumer needs it.
