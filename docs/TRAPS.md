# TRAPS — mistakes this project has already made

**Read this before designing anything.** Every entry is a real defect that shipped, or a real
measurement that was believed and was wrong. They are here because each one cost time, and several were
made twice.

⛔⛔ **THE SHAPE OF AN ENTRY, AND IT IS ENFORCED BY PRUNING.** A trap is **the mistake, the tell, and the
rule** — three or four lines. The investigation that produced it is in git; the number that describes
current state is in `ROADMAP.md`. ⚠ This file reached 741 lines and 91 entries before its first prune
(2026-08-05) because entries were being *appended to* as new instances arrived. ⭐ **Instances of ONE lesson
belong in a LIST inside that lesson, not as new entries** — see `D4`, which is nine costumes of one mistake.
⛔ Labels are cited from source docstrings and tests, so a label is NEVER renumbered or deleted; a merged
label survives as a row in its family's table.

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
deleting a whole channel scored **byte-identical on every class of every condition** — a clean, publishable
"this channel is inert". It was inert only because the monkeypatch missed the module that binds the name
(`from .bp_solver import X` makes a separate module global), and the first repair missed it too because
`from pkg import mod` returns a re-exported *function* of that name which shadows the submodule.
⭐ *The tell:* the ablation was a strict SUPERSET of another and scored a smaller effect. ⭐ *The rule:*
every ablation increments a counter and the harness **raises** if it did not fire.

**C0c. ⛔⛔⛔ A ZERO COUNT IS A MEASUREMENT OF A DENSITY, NOT AN ABSENCE OF DATA — AND KEYING PRECISION
ON THE COUNT MAKES THE STRONGEST STATEMENT IN THE LIBRARY THE QUIETEST.** At `g00` (zero gDNA by
construction) pass-0 attributes **34–38 %** of unstranded mass to gDNA, Σ|err| **1.56 M / 1.95 M**
fragments, and **100 %** of it is `relay_only` with **0 %** `own_evidence` and **0 %** `struct_lock` on
both axes. The reason is a single census: all **1,298** intergenic nodes hold **exactly zero** counts
over **50.7 Mb** of gDNA opportunity, and **all 1,298 emit nothing** — against 915–959 emitting at
`g25`/`g50`, so the anchor is healthy elsewhere and completely silent here. ⭐ Those nodes are
structurally pure gDNA and `struct_lock`ed, i.e. composition-**CERTAIN** — and `own_precision`'s
`n > 0` still zeroes them, so they are *certain and silent at the same time*. ⛔ **The general defect is
that `p = n/(n·Var+1)` reads the COUNT as the evidence when the claim is a DENSITY: zero over 50.7 Mb
and zero over 200 bp are the same number to this code and opposite statements about the world.** ⭐ The
derived repair needs no new constant and reuses the reference ψ already carries: a Poisson rate with
`a` events over exposure `E` under the Jeffreys prior ψ is built on (`_JEFFREYS_REF`) has posterior
`Gamma(a+½, E)` — proper, finite precision at `a = 0`, mean `0.5/E`. At `ΣE = 50.7 Mb` that is
`9.9e-9 ± 1.4e-8`, which would take the worst object from `f_g = 0.510` to ~`7e-9` against a truth of
exactly 0. ⚠ The population prior papers over this downstream (that object's shipped answer is 0.054),
but pass-0 **is** the substrate the prior is fitted on, and a prior fitted on "everything is ½" cannot
be rescued per-object on a heterogeneous sample.

**C0d. ⛔⛔⛔ THE RELAY'S CURRENCY IS A DENSITY *RATIO*, AND ZERO IS NOT IN THE MULTIPLICATIVE GROUP —
SO "THERE IS NONE HERE" IS UNREPRESENTABLE BY CONSTRUCTION, NOT BY AN OVERSIGHT.** C0c's repair was
built, gated (18 gates, 4 firing perturbations) and **measured at +0.0 % on the 36-condition panel, with
the four `g00` rows BYTE-IDENTICAL** — the conditions it was designed for did not move at all. The
two-step probe says exactly why: after the fix `NodeInit.prec_g > 0` at **1,298** intergenic slots
(value exactly `1/trigamma(½) = 0.2026`), and the relay publishes `prec_g > 0` at **0** of them. ⭐ The
reframe transports `r = ρ_tot(dst)/ρ_tot(src)` and every guard on that path is `rho > _EPS`
(`bp_solver` ~433, 506, 983, 1197): a ratio with a zero denominator is undefined, so a source whose true
density is **0** cannot participate however precisely it knows it. ⛔⛔ **That is why the single most
informative statement available at a zero-gDNA library — 50.7 Mb of opportunity with no fragments — is
structurally unsayable.** ⭐ The lesson generalises past this solver: *before fixing a source that will
not speak, check whether the channel can carry the value it wants to send.* A multiplicative transport
cannot carry zero, and no amount of precision at the source changes that. ⚠ And the fix is not wasted —
it is a precondition: the source must be able to form the claim before the channel can be taught to
carry it. But it must be priced as what it measured, which is nothing.

**A12b. ⛔⛔ A NOISE DEADBAND WHOSE CUSHION IS SUPPLIED BY AN UNRELATED SAMPLE SIZE FAILS EXACTLY WHERE
THAT SAMPLE GETS BIG — AND IT FAILS SILENTLY, INTO THE HONESTY COLUMNS.** `strand_evidence` gates the
strand channel on `disc = 4·max(0, (κ̂−½)² − σ²_d)` with `σ²_d = ¼(1/N_rna + ω_r) + ¼(1/N_gdna + ω_g)`.
On unstranded data the TRUE κ is ½, so `(κ̂−½)²` is pure fitting noise — a χ²₁ with mean `Var(κ̂) ≈
¼/N_rna`, which a χ²₁ exceeds ~32 % of the time. ⭐ **The `¼/N_rna` half of the deadband is exactly the
right bias correction; the `¼/N_gdna` half is a different job** (gating a gDNA-free library) doing the
first job's arithmetic. So the cushion **shrinks as gDNA grows**, and the gate opens on **exactly 3 of 18**
unstranded ladder conditions — `g75`/`g90 ss0.50 capture_off` and `g98 ss0.50 capture_on` — measured
`(κ̂−½)²` = 4.74e-07 against `σ²_d` = 3.42e-07 at g75, a margin of 1.39×. ⛔ **And `I_strand = N_eff·disc·[…]`
then MULTIPLIES noise by the depth**, so the phantom information grows linearly in coverage. ⭐⭐ The damage
is not accuracy — forcing κ = ½ moves Σ|err| by −0.6 %/−5.3 %/**+4.6 %** across the three big classes, a
wash, reproducing `ROADMAP` §3's −0.2 %. **The damage is that `solv%` goes 77–90 % → 0.0 %**: the panel
declared 92,154 objects confidently wrong and `calib` = 3.67 at a condition with *no solvable objects at
all*. ⛔⛔ **So the column used to pick the worst condition was inflated by the defect being looked for** —
A12's lesson with the arrow reversed. ⭐ The repair is to propagate rather than gate:
`I = (½−κ̂)² / [p(1−p)/N_eff + (1−f_g)²·Var(κ̂)]`, which is smooth, needs no `max(0,·)`, and whose
`Var(κ̂)` denominator cancels the depth exactly. `ROADMAP` §1 step 3.

**A12. ⛔⛔ EVERY HONESTY METRIC IMPROVES AS THE SOLVER STOPS KNOWING ANYTHING — SO NONE IS READABLE
WITHOUT A FIXED-DENOMINATOR ACCURACY NUMBER BESIDE IT.** A destruction control (force a 97-σ real strand
signal to ½) collapsed accuracy **0.0231 → 0.2064, +792 %, 16/16 worse** — while confidently-wrong Σ|err|
fell **87 %**, the declared-precision ratio went 4.53 → 0.67 and `weak%` went to 0, all 16/16 "better".
Three of those four are the Stage-B dashboard's headline columns. ⛔ **The rule: an A/B that moves `solv%`
has changed its own denominator, so quote `mwae` over ALL objects and the raw Σ|err| alongside — those two
cannot be gamed by knowing less.** ⚠ It also inverts single-condition readings: one row showed −64 % and
the same arm over 16 conditions is −0.2 %.

**A11. A TEST THAT RE-DERIVES A DEFINITION CANNOT DETECT DRIFT IN IT.** A gate existed precisely to keep
two instruments' shared class definition from diverging — and it recomputed that definition inline in
the test. Changing one instrument fired nothing. It is A1's shape once more: the check and the thing
checked have to come from *different* places, and for a shared definition that means **one home**, with
every consumer importing it. Here the home is production code (`node_geometry.g1_locked`), because the
predicate is a production concept and `scripts/` is deliberately not importable.

**A13. ⭐⭐ A GATE THAT RECONSTRUCTS A VALUE IS VACUOUS WHEREVER THE CODE DECIDES NOT TO USE IT — GATE A
SYMMETRY THE CODE CANNOT FAKE INSTEAD.** A relay gate recomputed the delivered gDNA level and asked which
candidate the relay was closer to. It passed the honest run **and the deliberately mis-paired one**, because
on an unstranded fixture the composition licence is withheld, the level crosses unscaled, and the ratio
being reconstructed never enters the observable. Eleven green gates; the one perturbation they existed for
fired nothing. ⭐ The repair was to stop reconstructing: on a PALINDROMIC chain the relay's two passes must
be exact mirror images, which needs no assumption about which channels are live. ⛔ **When a gate needs a
value to flow through a conditional path, it is a gate on that conditional. Prefer an invariance.**

**A14. ⛔⛔ "THE ARM CHANGED NOTHING" IS NOT A CONTROL UNTIL YOU CHECK THE ARM COULD HAVE CHANGED
SOMETHING.** A zero-gDNA condition was reported as a free falsification control because the new and old
reframe ratios agreed on every hop. They did, and the arm tested nothing: the composition licence was
withheld everywhere (no gDNA ⇒ no gDNA precision ⇒ no source can lend), and **6 of 8 hops per direction were
NO-FRAME** — the intron NODE had zero counts, so the reframe was skipped with `r` forced to 1 by the
divide-by-zero guard. The only framed hops were the ones the change does not touch.
⭐ **For every "no change" arm, count the opportunities the change had to fire and print that count beside
the result.** Same family as A13 and A8: *the measurement was never in a position to say anything.*

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

**B10. THE DEFAULT INSTINCT IS A POOLED AVERAGE, AND IT IS WRONG HERE THREE WAYS.** (i) averaging error
over objects that cannot be solved buries the answer — but an exclusion is a promise to check that class
some other way (B13); (ii) a pooled mean hides a sign flip between strata; (iii) mass-weighting lets one
huge object set the number. ⭐ Report per stratum, with the denominator named.

**B11. A BINARY CUT ON A FITTED PARAMETER'S RESIDUE IS NOT A POPULATION TEST — AND A BETTER THRESHOLD IS
NOT THE FIX.** `τ > 1e-9` promoted objects whose own statement was 1,377 nats wide into the scored
population, hiding a 1.06 M-fragment error. ⛔ A floor was implemented and refuted: τ is continuous across
the region, so any floor is a tuned constant. ⭐ The honest repair is to propagate the fitted parameter's
own variance, which drives the arm to ~0 by derivation (D4f).

**B13. ⛔⛔ EXCLUDING A POPULATION FROM THE DENOMINATOR WITHOUT GATING ITS OWN FAILURE MODE HIDES THE
LARGEST ERROR IN THE LIBRARY.** A condition reporting **mwae 0.0170** on a 39.5 %-solvable set had
**1,056,019 fragments** of error in the excluded class, with 6,425 objects moved a mean 0.414 from ½ and
100 % of them declaring a finite precision. The single worst mechanism was **0.0 % scored on every rung but
one**. ⭐ **The rule: every population you exclude needs its own gate, written at the same time as the
exclusion.** The tell is free — the *reason* to exclude a class is that its answer should be
uninformative, so measure how far from uninformative it actually is.

**B14. ⭐⭐ IF EVERY GATE FOR A CHANGE READS ONE INSTRUMENT, THE CHANGE IS GATED IN ONE PLACE — AND
`_relay`/`_transport` ARE TWO (now three, with `_flank_dom`).** Ten sound tests all read the combine's
publication; deleting the conjunct from `_relay` alone left **the entire calibration suite green**.
⭐ **For each place the change was made, name the observable that would move and confirm at least one gate
reads it.** A count of gates says nothing. ⚠ Re-opened by the very next change to the same pair, and again
by a third consumer nobody had listed (A13).

**B15. ⛔⛔ "STARVED" AND "DEPLETED" ARE NOT THE SAME DIAGNOSIS.** An object with almost no counts because
the experiment is too small is starved; one with almost no counts because the biology puts nothing there is
depleted. The first is fixed by more data, the second never is — and reading a depleted object as starved
threw away half a session chasing depth. ⭐ Decide which by asking whether the count grows with the lever
you have.

**B16. A TILING LOOP CAN SUPPRESS THE VERY POPULATION UNDER STUDY.** Toy capture probes are written in
TRANSCRIPT space, so a probe spanning an internal junction offset has a genomic footprint in **two blocks**
— and `sim/capture/sampler._split_scale` then multiplies every gDNA fragment overlapping it by
`gdna_split_penalty`. Tiling across the whole transcript therefore put a split probe over every internal
junction and depressed exactly the fragments that span an `intron|exon` EDGE. ⭐ The lesson is not "tile per
exon" (though that is the fix); it is that **a substrate knob can be adversarial to the object you are
measuring, and the way to find out is to read the sampler's weight for that object's own population** —
here, one `if len(blocks) > 1` in code nobody had reason to open.

**B22. ⛔⛔ AN RNA "LEVEL" THAT IS A MULTIPLE OF THE gDNA DENSITY IS NOT COMPARABLE ACROSS THE CAPTURE
AXIS.** At rung `m = 100` the exon's true `f_g` is **0.010 off capture and 0.31 on it**, because the
multiple is taken against the OFF-TARGET density while capture concentrates gDNA onto the exon ~65×. ⛔
**Key an operating point on a REALISED quantity you can see in the truth, never on a knob you set** — and
print the realised value so a row that cannot reach the target says so.

**B21. ⛔⛔ PRICE THE HALVES SEPARATELY — A DIAGNOSED DEFECT AND A VALUABLE ONE ARE NOT THE SAME DEFECT.**
Pricing two length models together read −6.31 % and looked like one finding; split one pmf at a time it is
**−5.90 % gDNA / −0.43 % RNA**, and all of the gDNA half is capture-ON. The half that was fully diagnosed
was worth ~nothing and the half worth 14× more had not been looked at.

**B18. ⛔⛔ A TOY CEILING IS NOT A PANEL CEILING, AND THAT HAS NOW COST THREE BUILDS.** A re-solve ceiling
worth −0.009 on the toy carried to the 36-condition ladder as **+0.0013 (worse)**; a second change was
toy-positive and panel-negative the same way; a third improved its target rung 26 % and regressed the panel
axis that carries the mass. ⛔ **Run the panel arm before writing a mechanism into `src/`.** ⚠ And when the
panel disagrees, ask whether the defect is half of a cancelling pair before concluding (D4j).

**B19. ⛔ "THE WELL-COUNTED SIDE" IS NOT A FIXED SIDE — CAPTURE INVERTS IT.** Off capture an intron NODE
holds ~349 gDNA counts against a flanking EDGE's 8–36; under capture the intron holds **1** and the EDGE
holds 20–40, unchanged from 120 kb to 1.08 Mb. So a rule that transports from "the well-counted side"
silently reverses direction with the protocol, and transferring a DENSITY instead of a SHARE measured
**+0.207** on capture-ON × unstranded.

**B20. ⛔ ADMITTING AN OBJECT TO THE SCORED POPULATION IS A COST, AND A MECHANISM CAN DO IT SILENTLY.** A
prototype gave a class certainty it had not earned; `solv%` rose 43.1 → 48.2 and the edge-axis mwae went
0.0216 → **0.1449** — not because the answers got worse but because wrong answers started being counted.
⭐ Report `solv%` beside every accuracy number (B13 inverted).

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
variance against a LINEAR-space error.** `Var(log f)` is not `Var(f)`; the delta method
(`Var(f) ≈ f²·Var(log f)`) converts, and at small `f` the two differ by orders of magnitude. Every
overconfidence figure from before that date is void, not merely imprecise.

**B8. A delta is only attributable if its baseline came from the same tree in the same session.** Re-record
the before-picture, do not quote a stored one.

---

## C. Pools, selections and divisors

**C0. ⭐⭐ TWO DIVISORS BUILT FROM ONE pmf CAN STILL DISAGREE — IF THEY RESPOND TO IT WITH OPPOSITE SIGN.**
`E_J = E[w]−1` RISES with the mean fragment length while `E_r = e−E[w]+1` FALLS, so a length-model error
appears as a junction-vs-exon frame gap of `Δ·(1/E_J + 1/E_r)` = **0.62 %/bp**. The two are exactly
consistent (202.8 + 797.2 = 1000.0) and still disagree by 10 %. ⭐ **When two quantities built from one
model disagree, differentiate BOTH with respect to that model before looking for a bug in either.**
⛔⛔ And the second term added here was WRONG, which is the sharper lesson: a "finite-transcript placement"
factor of 1.024 was sound arithmetic about a generative model **the simulator does not use** — it reweights
the length marginal by the same opportunity, so the factor cancels exactly. ⭐ **Before pricing a selection
effect in a simulated substrate, read the simulator's own sampling code — not the docstring, the code.**

**C0b. ⭐⭐⭐ A TERM THAT IS FRAME-FREE IS NOT THEREFORE ASSUMPTION-FREE — LOOK AT THE NUISANCE YOU
PROFILED OUT, NOT ONLY AT THE DIVISOR THAT CANCELLED.** `edge_spliced` is certified RNA, so the obvious
move is a coefficient `S` on ψ's RNA arm: `E[S] = c·(1−f_g)·M`, and every opportunity ratio lives in `c`,
which multiplies the MEAN and is therefore an additive constant in log space. That reasoning is correct —
C0's opposite-sign trap is structurally absent, neither divisor reaches the retained term — and it is
**not enough**, because the *same* `c` contains the splice-visibility `q`, whose dropped term
`−(q/(1−q))(1−f_g)M` is the same size as the term kept. ⛔ With `q` free the profile likelihood in `f_g`
is **exactly flat on [0,1)**: the whole channel carries ONE BIT, "`f_g ≠ 1`". Measured on all 36 ladder
conditions against origin-split truth, the raw-count term is **worse than the uninformative reference on
12 of them**, worst **+0.4578** mwae. ⭐⭐ **The tell was free and general: its benefit tracked the ANSWER,
not the evidence** — −0.49 at `g00` where the truth is all-RNA, +0.45 at `g98` where it is all-gDNA. **A
channel whose sign follows the truth is a prior, not information** (A12's shape, and the two zero controls
are what made it visible). ⭐ Ask of any "one-sided floor": *is the term I am dropping bounded, or is it
the dominant one?* `SESSION_HANDOFF` §8 · `certified_q_census.py` · `test_certified_rna_licence.py`.

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
**D4. ⭐⭐⭐ A MESSAGE COMPUTED FROM THE DESTINATION'S OWN BELIEF CARRIES ZERO INFORMATION AND CONFIRMS THE
DESTINATION.** **A message may use the destination's CONSTANTS — geometry, lengths — and its OBSERVATIONS;
never its BELIEFS.** Any "fix" that divides the destination's belief back out rebuilds the bug. ⭐ The check
is *"is the delivered value independent of the destination's own state?"*, and a prediction that does not
move when the data moves by four orders of magnitude is the tell.

⚠⚠ **THIS ONE LESSON HAS RECURRED NINE TIMES IN DIFFERENT COSTUMES.** Each sub-label below is cited from
source and tests, so the labels are kept; the investigations that produced them are in git, not here.

| | the costume it wore | the lesson that generalises |
|---|---|---|
| **D4** | delivered exactly `1/(1+f_own)`, reserving 33.6 % of the budget for imaginary gDNA so a zero-gDNA library read back 29.3 %; and at a gene end the rescale became `1/f(dst)`, 7× too big in the median, up to 190× | the rule above |
| **D4b** | ⛔ **a "TOTAL DENSITY" ratio is how D4 comes back.** `r = ρ_tot(dst)/ρ_tot(src)` re-creates it whenever `ρ_tot(dst)` is dominated by a component the message is not about: a correct 0.0257 gDNA density delivered as **28.684**, and `f_g` then FLAT at 0.90 across a 10,000× RNA sweep | **a scale factor must be built from the component the claim is ABOUT, never from a total** |
| **D4c** | and the fix was not a better scale — substituting the definitions shows `ρ_c(src)·r ≡ φ_c(src)·ρ_tot(dst)`, a **pure composition imputation with no level in it**, so no corrective factor exists | ⭐⭐ **before correcting an operator, substitute its own definitions and read what it delivers.** Two sessions were spent hunting a better `r` for an expression that never carried a level |
| **D4d** | the third D4 was inside the operator built to defend against the second: the mass pin's `k = 1/(φ_msg + R_own)` has a fixed point at `R_own = ½`, so it drove the delivered gDNA FRACTION to ½ regardless of truth | ⚠ **and it hid from every aggregate** — the per-step rescales telescope back to 1 at a gene's far end, so only a per-object check away from a pure-gDNA object can see it |
| **D4e** | gating the pin wherever it read the destination's own density looked clean and **broke the capture landscape** — the off-probe floor leaked into every exon | ⭐⭐ **state the licence as "no BELIEF may enter", not "the destination's numbers may not enter".** The first is D4 and admits the structural case; the second is a superstition that costs a real mechanism |
| **D4f** | κ̂ = 0.500689 on a genuinely unstranded library leaves `I(f_g) ∝ (2κ−1)²` at ~5e-7, and a licence that tests a PRECISION with no floor is granted by it | a **structural** zero represented by a point estimate whose own variance covers it; the repair is to propagate `Var(κ̂)`, not to add a floor (B11 refuted the floor) |
| **D4g** | the same shape with a much larger trigger: an intron's **2.5 % phantom RNA** is the only nonzero RNA precision in a chain, and it alone unlocks a reframe that compounds a gDNA level to **2.16×** | ⭐ **a predicate that gates a MULTIPLICATION must be sized by how much density stands behind it, not by whether a precision is strictly positive** |
| **D4h** | removing each ψ channel in turn moved the worst object by ≤0.016 of a 0.217 error while removing all of them moved it the whole way — because all three are built from the same relayed level | ⛔ **when single ablations are all small and the joint one is large, stop ablating consumers and go one stage upstream** |
| **D4i** | an "enrichment ratio" of 1.46 with capture OFF looked impossible and was correct: `ρ_tot` is a TOTAL and the destination held 55 RNA fragments where the source held 0. Compounded, 2.159 used vs 2.153 true | ⭐ **recompute the quantity from the ORACLE before assuming the formula is broken** — the alternative reading sends you to rewrite a correct function |
| **D4j** | ⛔⛔ **fixing one of two errors that CANCEL is worse than fixing neither.** D4i's leak cancels across a two-hop pair; correcting one hop alone moved a toy's evidence-free exon 0.0107 → **0.0244** while the rung it was aimed at improved 26 % | ⭐ **when a fix is negative and its object sits on a multi-hop path, price it in the arm that also removes the other defect.** A pair of defects that cancel is one experiment, not two |

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

**E14. TWO DIFFERENT MASKS SHARED THE WORD `struct_lock`, AND BOTH WERE RIGHT.** One answers "is this
belief pinned and certain?" (both axes); the other "may this slot EMIT composition certainty into its
messages?" (NODE-only on purpose, because a seam is structurally gDNA but sits between RNA-carrying exons).
⭐ Two correct predicates under one name is worse than either being wrong: give each its own name and ONE
home, and let every consumer import it (A11).

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

**F11. EQUAL CONFIGURED FRAGMENT LENGTHS DO NOT GIVE EQUAL REALISED ONES.** Handing the simulator one
mean for gDNA and RNA still yields different realised means, because each pool's length marginal is
reweighted by its own template opportunity — a 2 kb transcript truncates the tail that a whole chromosome
does not. ⭐ So "equal lengths" is a claim about the CONFIG and never about the library; gate the axis on the
realised truth (A7).

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
