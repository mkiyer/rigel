# SESSION HANDOFF — the `g00` mechanism is FOUND and FIXED; the rebuild is still panel-negative, and the residual now has a shape

    Rewritten 2026-08-06 (fourth pass). ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** The six permanent
      docs are `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **DELETE THIS FILE when its task lands**, promoting the numbers into `ROADMAP.md`, the lessons
      into `TRAPS.md`, the derivation into `EQUATIONS.md` and any ruling into `DESIGN.md`.
    ⭐ The previous handoff (five eliminated `g00` candidates) is **superseded** — the mechanism was none
      of them, and §1 below says what it was.

⛔⛔ **THE DESIGN DOC IS `docs/calibration/variance_model_notes.md`. READ ITS §2 AND §6 BEFORE PROPOSING
ANYTHING** — §2 is a refutation and §6 is a graveyard of **eleven** priced mechanisms. ⚠ Nothing below
re-opens any of them.

---

## §1 ⭐⭐⭐ THE `g00` MECHANISM — A UNIT ERROR IN THE TILT CHANNEL, CARRYING 74 % OF THE ERROR

The prototype delivered a raw **log-odds** into ψ's `theta_imp` slot. ψ's tilt grid is the **ANGLE**
`theta = arcsin(tau)` and `simplex_logodds._tilt_grid` spans **exactly `[-pi/2, +pi/2]`**. Measured modes
were `±4.6` — **2.9× outside the entire domain of the coordinate.**

⭐ **What that does, and it is not subtle.** ψ's message term `-1/2 * p * (theta - m)^2` with `m` off the
grid is **MONOTONE across the whole grid**, so the tilt pins at the boundary `tau = ±1` ("all RNA on one
strand"). The nuisance is no longer free to integrate out, the AMBIG Schur protection is gone, and the
strand likelihood explains whatever it then cannot fit by calling the mass **gDNA**. ⚠ The tell was
visible and unread: the delivered `f_g` was the SAME value `0.9898` at slot after slot — a **grid point**
of the K=60 AMBIG cube, i.e. a posterior pinned into one cell.

**How it was found, and why five hypotheses missed it.** Ranking every slot by error MASS and printing
each channel's precision beside it showed the worst slots had **every printed message precision at exactly
zero** — own self-solve `0.053`, delivered `0.9898`. Since `build_node_init` and the sweep call ψ with
**argument-identical** inputs apart from the four `*_imp` channels, and all the printed ones were inert,
the cause had to be a channel **nobody had printed**. ⛔ The tilt was the one message stream absent from
the prototype's `_capture`, so no instrument could rank it. `TRAPS.md` A14's shape, one level up: it was
not that the arm could not have fired, it was that the arm could not be SEEN.

⛔⛔ **AND THE COORDINATE ALONE WAS NOT THE FIX — MEASURED.** ψ's `tau` is the tilt of the **RNA**; a raw
COUNT tilt is the tilt of the **MASS**, and F2 relates them by `tau_obs = (2*kappa - 1)*(1 - f_g)*tau`,
which on an antisense protocol **INVERTS THE SIGN** (`kappa = 0.0101` at `ss = 0.99` ⇒ `2k-1 = -0.98`).
Correcting only the coordinate left `g00` at 1,568,412 — no better than the 1,516,976 it started at,
because it then asserted the right angle with the wrong sign at a precision of `n`. ⭐ The tilt is now
built from `ni.rho_pos`/`ni.rho_neg`, which are `f_pos*M/E_r` and `f_neg*M/E_r` — **ψ's own coordinates up
to a common factor that cancels in the ratio** — so it needs neither `kappa` nor a belief. That is
`bp_solver`'s own `tth` formula, reached by derivation rather than by copying it.

The precision follows by the delta method with **no new constant**:
`Var(theta) = ((1 - tau^2)/4)*(1/p_+ + 1/p_-)`.

**ψ-boundary ablation, `g00 ss0.99 capture_on`** (re-solve with one channel nulled; the `as-is` arm
reproduces the run BIT-IDENTICALLY on solvable slots — A5):

| arm | Σ\|err\| | |
|---|---|---|
| the arm as-is | 1,751,145 | |
| − `theta_imp` | **452,326** | ⭐ **74 % of the error, one channel** |
| − `rna_imp` as well | 170,142 | 16 % |
| every message muted | 170,142 | = the shared self-solve exactly |
| HEAD | **1,677** | |

**Result at that condition: Σ\|err\| 1,516,976 → 403,073 (−73 %)**, and slots with NO message of any kind
went 1,075,790 → **85,528** (−92 %).

---

## §2 ⛔ THE PANEL — +103 % → +85 %, STILL NEGATIVE, AND NOW WITH A SHAPE

`abs_err_all_final`, 32 conditions, `g00` excluded. Baseline re-recorded from the current tree in the same
session (B8); both arms `--jobs 8`.

| stratum | base | `eta` | | pass-0 |
|---|---|---|---|---|
| **ALL** | 22,807,907 | 42,209,808 | ⛔ **+85.1 %**, better on 17/64 | +78.9 % |
| stranded × capture ON | 3,469,926 | 3,662,896 | +5.6 % | **−5.5 %** |
| stranded × capture OFF | 1,148,057 | 2,261,200 | +97.0 % | +83.8 % |
| unstranded × capture ON | 16,661,720 | 20,781,592 | +24.7 % | +40.5 % |
| ⛔ unstranded × capture OFF | 1,528,204 | 15,504,121 | **+914.5 %** | +375.1 % |
| ⛔ `g00` zero control | 28,534 | 7,812,247 | (unbounded — truth is 0) | |

⭐⭐⭐ **THE RESIDUAL IS MONOTONE IN gDNA CONTENT ON UNSTRANDED DATA, AND THAT IS THE WHOLE DIAGNOSIS.**

| unstranded, node axis | base | `eta` | |
|---|---|---|---|
| `g98 capture_on` | 2,000,545 | 482,660 | ⭐ **−75.9 %** |
| `g90 capture_on` | 2,031,696 | 644,325 | ⭐ −68.3 % |
| `g75 capture_on` | 1,918,865 | 851,776 | ⭐ −55.6 % |
| `g25 capture_off` | 83,925 | 1,691,130 | ⛔ +1,915 % |
| `g05 capture_off` | 51,351 | 2,237,924 | ⛔ +4,258 % |
| `g01 capture_off` | 57,681 | 2,252,758 | ⛔ +3,806 % |

Large wins where gDNA is genuinely high, catastrophic losses where it is low. ⛔ **That is not a mechanism
that knows anything — it is a systematic UPWARD bias on `f_g`**, and it is §6's lever one more time. It
sits on exactly the stratum where `kappa = 1/2` leaves the strand channel silent (F2: `I(f_g) = 0`
exactly) and the MEASUREMENT channels carry the entire budget.

---

## §3 ⭐⭐⭐ THE ROOT CAUSE — A ONE-SLOT-STEP CHANNEL ON AN ALTERNATING CHAIN IS **BIPARTITE**

⛔⛔ **THIS SUPERSEDES §3a BELOW, WHICH WAS A GUESS FROM THE PANEL SHAPE. The dissection says something
sharper and it is provable from the code.**

**The dissect loop, on the two worst regressions** (`g01 ss0.50 capture_on`, `g05 ss0.50 capture_off`),
ψ-boundary ablation with the A5 identity holding:

| | `g01 ss0.50 cap_ON` | `g05 ss0.50 cap_OFF` |
|---|---|---|
| the arm as-is | 4,436,114 | 4,939,517 |
| **− ALL FOUR MESSAGE CHANNELS** | **4,278,613 (−3.6 %)** | **4,256,289 (−13.8 %)** |
| HEAD | **162,983** | **108,561** |

⭐⭐ **Muting every message the rebuild sends makes it BETTER on both.** At the top-error slots
`own == ETA == noThe == noRna == noMsg` — identical to four decimals — and **all four message precisions
are exactly 0**. The rebuild is not adding a bias on this stratum. **It is delivering nothing at all**, and
the shared self-solve's honest "≈ ½" is then the whole answer.

**And that is the entire solve on unstranded data.** `worst_objects.py` on the same condition, HEAD arm:
**100.0 %** of the top error is `relay_only`, `fg_loc` (the message-free self-solve) is
**0.471–0.568** at every one of the top 12 objects, and HEAD's messages move them to **0.009–0.339**
against truths of 0.001–0.024. At κ = ½ the strand λ-term is exactly 0 (F2), so a high-mass exon has no own
composition evidence at all and the message layer *is* the solve.

### §3.1 ⛔⛔⛔ WHY THE REBUILD DELIVERS NOTHING — a census, and then a proof

Receivers of each channel at `g00`, by slot type (**no solver reasoning, just who got a nonzero
precision**):

| channel | EDGE | intergenic NODE | intron NODE | exon NODE |
|---|---|---|---|---|
| gDNA level / measurement | **2,592** | 0 | 0 | **0** |
| certified RNA | **0** | 0 | 0 | 10,493 |

⭐ **Neither channel ever reaches a slot on its emitter's OWN axis, and that is a theorem rather than a
measurement.** The chain is strictly `N E N E … N`, so a **one-slot-step** channel is **BIPARTITE**: an
emitter that is a NODE can only ever reach an EDGE, and an emitter that is an EDGE can only ever reach a
NODE. `emits_level = pure_gdna & is_node_slot` is **NODE-only**, therefore **no NODE can ever receive a
gDNA level** — the channel is not weak, it is *disconnected* from every object that carries mass. The
census's old finding that "the gDNA channel reaches only slots with `M = 0`" was this, seen from one
condition and read as a property of `g00`.

⛔ **It is a derivation gap, not a coding slip.** §5 licenses the shared density to cross between two
**ADJACENT OBJECTS**, and every worked pair in §3's table is a NODE↔EDGE pair — but **a NODE's nearest
NODE is TWO steps away**, and the only licensed originators of a level (structurally pure-gDNA objects)
are NODEs. "Adjacent" was implemented as "adjacent slot", and on an alternating chain that silently
deletes the entire same-axis relation.

### §3.2 ⭐⭐ AND THE ARROW IS BACKWARDS FROM WHAT THE PROTOTYPE INTENDED

`eta_node_sweep:243-245` states the design: *"The shipped mass pin's case (ii) hands an exon its flanking
EDGE's OWN enriched measurement for exactly this reason; a one-step level keeps that property without the
pin."* ⛔ **The census says that property is NOT kept.** What is implemented is the mirror image — an EDGE
receiving its flanking *intergenic NODE's* measurement. The exon never hears from its own flanking EDGE,
because an EDGE never **originates** a level.

⚠ **And the naive repair is already refuted.** Simply letting the level cross two steps
(intergenic NODE → EDGE → exon NODE) is the chain-fused level the prototype measured and rejected: every
exon inherits the OFF-PROBE floor, and the true gDNA density at an empty anchor's neighbours is **346×**
what the anchor reports (`anchor_opportunity_census.py`). ⭐ **So the arm to run is specifically: can an
EDGE with its own observed mass ORIGINATE a gDNA level into its flanking NODEs?** That is HEAD's case (ii)
without the pin, it is one thing varied, and it is the first candidate this campaign has produced that is
about **connectivity** rather than about how to resolve doubt.

⛔ **The licence is an owner call and it touches settled ground** (`DESIGN.md` §4.1, the mass pin's
structural case). `eta_node_sweep:137-146` argues an EDGE's `|T| = 1` does not license "the mass here IS
gDNA" — which is correct, and is about **composition certainty**. Whether it also forbids originating a
**LEVEL** is a different question and is not answered anywhere.

### §3.2b ⛔⛔⛔ AND THE SAME ABLATION RUN ON **HEAD** FINDS THE REAL DEFECT — IT IS NOT A REBUILD BUG

⛔ **§3.2's "let an EDGE originate a level" was a patch on a symptom and is WITHDRAWN as the lead
candidate.** The owner's objection is correct and sharper than the proposal: on a bipartite chain a
message must flow *through* both axes, so the question is not who may originate but **why the LEVEL is not
a relayed state at all** while `eta` and `theta` are. ⭐ And `eta` is degenerate exactly where the level
lives — at a pure-gDNA slot `f_g = 1`, so `eta = +inf`, and `eta_node_sweep:175` explicitly zeroes
`prec_own` there. **The pure-gDNA regime has no representable message in the relayed coordinate.** That is
the structural statement; §3.2's version of it was a special case.

⭐⭐⭐ **THEN THE SAME ψ-BOUNDARY ABLATION WAS RUN ON THE SHIPPED SOLVER, and it moves the whole problem.**
`g01 ss0.50 capture_on`, final arm, **A5 exact** (the re-solve reproduces HEAD's own answer with
`max |Δ| = 0.000e+00`, the write-back reproduced):

| arm | Σ\|err\| | | could move |
|---|---|---|---|
| as-is (HEAD) | 303,512 | | — |
| − `gdna_imp` (the LEVEL) | 310,363 | +2.3 % | 24,033 |
| − `rna_imp` (certified RNA) | 324,178 | +6.8 % | 9,727 |
| − `lam_imp` (composition) | 318,315 | +4.9 % | 16,587 |
| − `theta_imp` (tilt) | 303,335 | −0.1 % | 3,291 |
| **− ALL messages** | **487,652** | **+60.7 %** | 24,033 |

⚠ Every single ablation is small and the joint one is large — **`TRAPS.md` D4h's exact shape**, "stop
ablating consumers and go one stage upstream", because all four are built from the same relayed state.

⛔⛔ **THE PER-SLOT TABLE IS THE FINDING, AND IT INVERTS THE CAMPAIGN'S PREMISE.** At every one of HEAD's
top-12 error slots the **self-solve WITH the fitted prior is nearly correct and the MESSAGES destroy it**:

| slot | mass | truth | own (self-solve) | **HEAD** | `prec_rna` | **− rna_imp** |
|---|---|---|---|---|---|---|
| 6,658 exon | 83,514 | 0.0007 | 0.0043 | **0.3392** | 35.8 | **0.0080** |
| 27,088 EDGE | 5,182 | 0.0023 | 0.0086 | **0.3219** | 327 | **0.0005** |
| 40,950 EDGE | 3,153 | 0.0032 | 0.0173 | **0.3392** | 190 | **0.0001** |
| 10,373 EDGE | 2,866 | 0.0028 | 0.0160 | **0.3219** | 349 | **0.0013** |
| 58,211 exon | 7,670 | 0.0023 | 0.0086 | **0.6430** | 124 | **0.1037** |

⭐⭐ **Removing ONE channel — the certified-RNA measurement — recovers the truth almost exactly at 8 of the
12**, while removing it *globally* is +6.8 % worse. **It helps the many and destroys the few, and the few
are where the error mass is.** That is a CONCENTRATED defect, which is what the dissect loop exists to find.

### §3.2c ⭐⭐⭐ THE MECHANISM, AND IT IS THE DEFECT THE ASSERTION ALREADY FOUND IN THE REBUILD

ψ applies `-1/2 * p * (log f_active - mo_p)^2` with `mo_p = log(cp*E_r/M)` — a **two-sided** Gaussian. So a
certified-RNA claim that **UNDER-states** the RNA share pulls `f_active` down and therefore drives `f_g`
**UP**. On an unstranded 1 %-gDNA library an exon's mass is essentially all RNA (`f_active ≈ 1`), so a
message saying "the RNA share is 0.3" is read as "70 % of this is gDNA" — a 200× over-call, carried at
`prec = 327`.

⭐ **And the claim is a LOWER BOUND, which `bp_solver` itself states.** P1d: `rho_R(exon) >= rho_nu(B) +
rho_mu(B)` — the exon may hold molecules that never touch that seam. HEAD damps it with
`graft_premise_logvar` and `graft_frame_logvar`; **the damping is evidently not enough, and the bound is
still delivered two-sided.** ⛔ This is the SAME defect `share_from_density`'s new assertion caught in the
rebuild from the other end (an implied 85,477 fragments at a slot holding 5) — two independent routes, one
mechanism.

⭐⭐ **THE DERIVED FIX IS TO MAKE THE CHANNEL ONE-SIDED, and it is the first candidate this campaign has
produced that CANNOT be refused by the zero control:**

    psi += -1/2 * p * max(0, mo_p - log f_active)^2

No penalty when the destination holds MORE RNA than the bound; full penalty when it holds less. It is
read directly off the premise the code already documents, it introduces **no new constant**, and — the
part that matters against §6's graveyard — **it can only ever say "at least", so it adds doubt in NEITHER
direction.** At `g00` the truth is `f_active = 1`, which satisfies every bound, so the channel goes
**inert** there rather than harmful. ⭐ It also predicts the observed panel shape: an under-claiming bound
is maximally damaging exactly where `f_active ≈ 1`, i.e. at LOW gDNA, and nearly harmless at high gDNA —
which is the monotone-in-gDNA residual, in both arms.

⚠ **Unmeasured. Run it as a panel arm on HEAD first** (it is a change to ψ's factor, one thing varied, and
HEAD is where the defect is demonstrated) — B18, and the rebuild is not the thing to fix first.

### §3.3 ⚠ AND THE `g98` "WIN" IS PROBABLY NOT A WIN

`g98 ss0.50 capture_on` reads −75.9 %, the largest single-condition improvement this campaign. ⛔ But if
the rebuild's answer on unstranded relay-only exons is its self-solve — the ψ reference shaped by the
fitted gDNA hyperprior, i.e. a quantity that tracks the LIBRARY's overall gDNA content — then it scores
well at `g98` because the default sits near a truth of 0.98 and catastrophically at `g01` because it does
not. **That is `TRAPS.md` C0b's tell exactly: a channel whose benefit tracks the ANSWER rather than the
evidence is a prior, not information.** ⚠ Unmeasured as stated; the cheap check is whether the delivered
`f_g` at relay-only exons is predicted by the fitted prior's median at that slot's own density, and it
should be run before any `g98` number is quoted as a win again.

## §3a ⚠ SUPERSEDED — the "coverage" framing, kept only because §3 replaced it

⛔ **The certified-RNA (spliced) measurement channel.** `bp_solver` records that it carries **75 %** of the
posterior precision on the confidently-wrong unstranded × capture-OFF exons and is **the only thing that
lets a zero-gDNA library say "my mass is all RNA"** (`gdna_none` 0.1063 → 0.1438 when removed). The
prototype's version of it is strictly weaker in three ways, all measured:

1. ⛔ **It is ONE STEP and EXON-NODE-ONLY.** `eta_node_sweep` delivers `spl_rho` from an immediate
   neighbour into `is_exon[i]` slots. **HEAD CARRIES IT INSIDE THE RELAYED STATE**, so junction evidence
   travels arbitrarily far along the chain even though it never travels along the junction axis
   (`bp_solver:1010` folds the grafted flux into the transported claim; `:1104-1110` writes the fused
   density AND the measurement precisions at `i`). ⚠ At `g00`, EDGE slots receive the channel **0.0 %** of
   the time in the prototype.

   ⛔⛔ **AND THIS IS NOT A FORWARD-BACKWARD-vs-GAUSS-SEIDEL DIFFERENCE — A PREVIOUS DRAFT OF THIS FILE
   SAID IT WAS AND THAT WAS WRONG (corrected 2026-08-06).** Both sweeps run **one forward pass and one
   backward pass**, both accumulate a running state in place (read at `s`, write at `i`), and both combine
   at slot `i` by reading each neighbour's state and **never `i`'s own** — which is the BP-legal rule, and
   `eta_node_sweep`'s combine states it explicitly. `bp_solver:937` calls `_relay` "a SEQUENTIAL
   Gauss-Seidel scan", but that is HEAD's word for *in-place and therefore un-vectorisable* (the sentence
   exists to justify the 15.7× twin), **not** for an iterative fixed-point in place of exact FB. On a chain
   an in-place forward accumulation IS the forward half of forward-backward. ⭐ **The real difference is
   WHAT IS IN THE RELAYED STATE:** HEAD relays ten arrays (three component densities, three mode
   precisions, three MEASUREMENT precisions, `tau`); the rebuild relays four (`eta`, `theta` and their
   precisions). The gDNA level and the certified-RNA measurement are deliberately **one-step** in the
   rebuild, with a measurement behind that choice — see the note at `eta_node_sweep:265`, where a
   chain-fused LEVEL was dominated by the 1,312 intergenic anchors and handed every exon the off-probe
   floor (346×). ⛔ So "restore chain propagation" cannot be adopted on the strength of the mislabel; it
   has to beat that recorded measurement, and the honest form of the question is **which streams belong in
   the relayed state**, not which iteration scheme is used.
2. ⛔ **It is a LOWER BOUND used as an EQUALITY at full count precision.** `bp_solver`'s P1d states the
   premise (`rho_R(exon) >= rho_nu(B) + rho_mu(B)`) and charges it TWO variances (`graft_premise_logvar`,
   `graft_frame_logvar`) with the mass pin behind them. The rebuild deleted all three.
3. ⚠ **And it is two-sided.** ψ reads `mo_p` as "the RNA share is exactly X", so an UNDER-claiming RNA
   measurement drives `f_g` **up**. After the tilt fix, 79 % of the remaining `g00` error sits on slots
   that HAVE this message.

⭐ **Why this one is different from the eleven in the graveyard:** every one of those was a rule for how to
resolve DOUBT, and all eleven died on the zero control because at `g00` the doubt must resolve to *no*
gDNA. This is a claim about **COVERAGE** — a channel that HEAD delivers everywhere and the rebuild delivers
almost nowhere. It has a direction the control agrees with.

⚠ **Price it as an arm before building it** (B18, four builds and counting).

---

## §4 THE CODE WRITTEN THIS SESSION

| file | what changed |
|---|---|
| `tests/calibration/_eta_reference.py` | ⭐ `+ tilt_angle` (the coordinate + the delta-method precision, both derivations in the docstring), `+ share_from_density` (the ASSERTION that replaced the clip) |
| `scripts/design/eta_node_sweep.py` | the tilt channel rebuilt; the `[0,1]` clip → assertion; a structural `M > 0` gate on the measurement channels; over-unit claims REFUSED and COUNTED (`REFUSED`); ⭐ the tilt published into `_capture` |
| `tests/calibration/test_eta_transfer.py` | **24 → 33 gates**; 8 perturbations, **all 8 fire their own named gate** |

**Suite: 22 failed / 2,323 passed / 7 xfailed** — the standing 22 (21 goldens + the paralog row) and the
7 xfails unchanged; 2,314 + the 9 new gates. `ruff check src/ tests/ scripts/` clean. ⛔ `src/` still
carries **docstring corrections only** — no behaviour change.

### §4a Two things the assertion found on its FIRST run, which the clip had been hiding

1. **A share of ZERO MASS.** A junction's certified-RNA density of 45.9 delivered as a "share" of a slot
   holding `M = 0` — 5.9e304, clipped to `log 1 = 0`, i.e. a confident `f_g = 0`. Fixed structurally
   (`M > 0`), and ⚠ **expected inert on the deliverable and not credited with anything**: `solvable`
   already excludes those slots.
2. **An over-unit share with real mass.** 10.3 frag/base × 8,300 b of RNA opportunity = 85,477 fragments
   claimed at a slot holding **5**. That is §3's item 2. At `g00 ss0.99 capture_on` alone, **5,890 claims**
   over 2.4 M fragments of mass are refused.

⭐ **The refusal adds doubt in NEITHER direction, which is why it is the only admissible reading.**
Saturating at `share = 1` is §6's lever pointed at RNA instead of gDNA.

### §4b ⛔ A GATE OF MINE WAS VACUOUS, AND A2's SECOND HALF IS WHAT CAUGHT IT

`test_a_slot_with_rna_on_one_strand_only_...` did not fire when its own named perturbation was applied,
because a redundant `var > 0` backstop was silently doing the guard's job — **two homes for one
predicate** (A11), written by me, in the same file as the fix. Closed by making `live` the whole guard and
writing `1 - tau^2` as `4*rho_+*rho_-/(rho_+ + rho_-)^2`, which is strictly positive whenever both
densities are, where `1 - tau*tau` rounds to exactly zero as soon as one strand dominates.

---

## §5 ⛔ WHAT IS SETTLED — do not re-open

Unchanged from the previous handoff: the certified-RNA λ **term** (REFUTED, C0b — ⚠ **not** the same thing
as §3's certified-RNA measurement CHANNEL), the vertex problem, `Var(kappa_hat)`/the strand deadband, a
resolving-power floor on `tau`, `length_likelihood` (stays OFF), the gDNA prior's bimodal capacity, the
prior's mode LOCATION on the regressing strata, §2's Jeffreys mean location, and all five earlier `g00`
candidates. **New this session:** the `g00` mechanism itself (§1) — it was the tilt coordinate, and none of
the six standing hypotheses.

## §6 ⭐⭐ THE TASK

1. ⭐⭐⭐ **MAKE THE CERTIFIED-RNA CHANNEL ONE-SIDED, AND RUN IT ON `HEAD`** (§3.2c). It is a lower bound
   the code already documents, delivered to ψ as a two-sided equality; it destroys HEAD's answer at 8 of
   its 12 worst slots (0.0086 → 0.3219 against a truth of 0.0023, recovered to 0.0005 by muting it) while
   being a net help globally. ⭐ One thing varied, no new constant, and it **cannot be refused by the zero
   control** because "at least this much RNA" adds doubt in neither direction and goes inert at `g00`.
   ⛔ Fix it in HEAD first — that is where the defect is demonstrated, and the rebuild is not the thing to
   repair first.
1b. ⛔ **§3.2's "let an EDGE originate a level" is WITHDRAWN as the lead candidate** — it was a patch on a
   symptom. The structural statement that survives is §3.2b's: the LEVEL is not a relayed state while
   `eta` and `theta` are, and `eta` is degenerate exactly where the level lives. Keep it as the frame for
   the rebuild, not as the next arm.
1c. ⚠ **Check whether the `g98` unstranded "win" is the fitted prior's median rather than an answer**
   (§3.3). Cheap, no solver.
2. ⭐⭐ **`ROADMAP.md` §4.1 is still untouched** — *why is prior fidelity anti-correlated with deliverable
   quality?* The cheapest probe is unchanged: inject the truth-fitted prior and re-solve one condition. §3
   of the previous handoff predicts the deliverable barely moves; if that holds, "pass-0 is the training
   set" needs rewriting and every pass-0 ranking in `ROADMAP.md` is suspect (A15).
3. ⚠ **The four STRICT xfails in `test_node_init.py` have NOT gone green** — the rebuild is not in `src/`.
