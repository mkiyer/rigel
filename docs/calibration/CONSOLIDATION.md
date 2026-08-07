# CONSOLIDATION — the backbone, and the measurement that sizes it

    Opened 2026-08-07. ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** ⛔ **DELETE THIS FILE when the
    backbone lands**, promoting the structure into `DESIGN.md`, the numbers into `ROADMAP.md`, the
    lessons into `TRAPS.md`.

⛔ **THE QUESTION THIS FILE ANSWERS:** we have two solvers — the shipped `bp_solver` (1,805 lines, evolved,
every operator carrying a measurement) and the `eta` rebuild (392 lines, clean derivation, panel-negative).
Which do we build on? And what is the smallest thing we can trust?

---

## §1 ⭐⭐⭐ THE MEASUREMENT THAT DECIDES IT — the whole message layer at pass-0 is worth **+0.2 %**

⛔ **Pass-0's ONLY job is to be a training substrate for the gDNA hyperprior** (owner, 2026-08-07). It is
not the deliverable, it does not have to be accurate, and it does not have to answer objects it cannot
solve. **Nobody had ever asked what the message layer contributes to that job.**

`--arm msgfree_p0` mutes ψ's four imputed channels (`gdna_imp`, `rna_imp`, `lam_imp`, `theta_imp`) during
the PRIOR-FREE sweep only; every refit sweep runs untouched. The relay still runs — only its delivery into
ψ is withheld — so the geometry, the self-solve, the reference and the fitted prior are all identical and
ONE thing is varied. 36 conditions × 2 axes, baseline re-recorded from the same tree in the same session:

| stratum | base | `msgfree_p0` | **DELIVERABLE** | pass-0 |
|---|---|---|---|---|
| **ALL** (g00 excluded) | 22,807,907 | 22,847,504 | ⭐ **+0.2 %** | +77.5 % |
| stranded × capture ON | 3,469,926 | 3,469,678 | **−0.0 %** | ⭐ **−56.1 %** |
| stranded × capture OFF | 1,148,057 | 1,147,595 | **−0.0 %** | ⭐ −15.7 % |
| unstranded × capture ON | 16,661,720 | 16,702,519 | +0.2 % | +57.8 % |
| unstranded × capture OFF | 1,528,204 | 1,527,712 | **−0.0 %** | +359.1 % |
| ⛔ `g00` zero control | 28,534 | 28,534 | **+0.0 %, byte-identical** | +6,132 % |

⭐⭐⭐ **Pass-0's own error gets 77.5 % worse and the shipped answer does not move.** The `g00` control is
byte-identical and 7 of the 8 control rows for `library_f_gdna` are unchanged to four decimals.

⚠ **A14 is satisfied in the strongest possible form.** This is not a dead arm: it demonstrably fired and
moved pass-0 by 77.5 % panel-wide and by 6,132 % at `g00`. "No effect on the deliverable" is therefore a
measurement and not a plumbing failure.

⭐ **And on STRANDED data the pass-0 messages are actively HARMFUL** — muting them improves pass-0 by
56.1 % on capture-ON and 15.7 % on capture-OFF. The substrate is *better* without them.

### §1.1 What that means, stated plainly

**The reframe `r`, the mass pin, splice-in/splice-out, the flank pair, `_flank_dom`, the DerSimonian–Laird
deflation and the whole scalar/vector twin pair — every operator this project has argued about for months
— contribute 0.2 % to what ships, at pass-0.** Eleven candidate mechanisms, the `eta` rebuild, and three
sessions of dissection were all priced on a stage worth two tenths of one percent.

⛔ This does **not** say the message layer is useless. It says it is useless **at pass-0**, which is the
only place it currently runs prior-free. Its value, if it has one, is in the REFIT — and `--arm
msgfree_all` prices that separately (§2).

---

## §2 ⭐⭐⭐ THE OTHER HALF — the messages ARE worth ~50 %, and ALL OF IT IS ONE STRATUM

`--arm msgfree_all` mutes the four channels in EVERY sweep — the floor: what the accumulator, ψ, the
self-solve and the fitted prior are worth with no belief propagation whatsoever.

| stratum | base | `msgfree_all` | **DELIVERABLE** | |
|---|---|---|---|---|
| **ALL** (g00 excluded) | 22,807,907 | 45,586,411 | **+99.9 %** | so the messages are worth −50 % overall |
| ⭐ stranded × capture ON | 3,469,926 | **1,448,555** | **−58.3 %** | **16/16 better** |
| ⭐ stranded × capture OFF | 1,148,057 | **646,188** | **−43.7 %** | **16/16 better** |
| ⛔ unstranded × capture ON | 16,661,720 | 42,454,062 | **+154.8 %** | **0/16 better** |
| ⭐ unstranded × capture OFF | 1,528,204 | **1,037,605** | **−32.1 %** | 14/16 better |
| ⛔ `g00` zero control | 28,534 | 2,312,277 | +8,003 % | |

⚠ Pass-0 is identical in both arms (+77.5 %), as it must be — an internal consistency check that the
`p0`/`all` switch does what it says.

⭐⭐⭐ **ON THREE OF THE FOUR STRATA THE MESSAGE LAYER IS A NET HARM, AND ITS ENTIRE VALUE IS
CONCENTRATED IN ONE STRATUM — unstranded × capture-ON, which is 73 % of the panel's total error.**

⭐⭐ **And that stratum is exactly where the strand channel is silent.** At κ = ½ the strand λ-term is
**exactly 0** for any count and any overdispersion (F2), so an unstranded slot has no own composition
evidence at all and the messages are the only source. On stranded data the self-solve plus the fitted
prior is **better than the messages**, by 58.3 % and 43.7 %, on every single condition.

⛔⛔ **THIS IS THE STANDING `stranded × capture-ON` RESIDUAL, AND IT IS NOT A MISSING MECHANISM.**
`ROADMAP.md` §0 records that regression at +30 % through the refit; eleven candidate mechanisms and three
sessions have been spent on it. **Muting the message layer on that stratum is −58.3 %, 16/16.** The
operators were being tuned on a panel where three of four strata would be better off without them.

### §2.1 What the message layer's job actually is

**One sentence: carry composition information to slots where the strand channel is silent.** That is a far
smaller job than 1,805 lines, it is scoped by a predicate the code already computes (`tau_lam == 0`), and
it explains the two-corner problem in §4 — a mechanism tuned to help the one stratum that needs messages
was being scored against three strata that are harmed by them.

⚠ **`msgfree_all` is a MEASUREMENT, not a proposal.** Do not ship it: `g00` goes to +8,003 % and
`g98 ss0.50 capture_on` reads `library_f_gdna` 0.0576 against a truth of 0.9844. The action is to make the
message layer's SCOPE match where it demonstrably helps, and to measure that as its own arm.

---

## §3 ⭐⭐ THE DECISION — build on HEAD, restructured, gated on BYTE-IDENTITY

⛔ **Not the `eta` rebuild, and not for accuracy reasons.**

**You cannot gate a rewrite.** The rebuild's +103 % turned out to be three independent things — a correct
derivation (the `η` transfer, which is sound and carries 33 gates), a UNIT ERROR (a log-odds delivered
into an angle's grid), and a STRUCTURAL disconnection (a one-slot-step channel on a bipartite chain
reaches 0 NODEs). It took two sessions to separate them, and one of the three was a typo-class error that
no amount of derivation review would have caught. ⭐ **A refactor gated on byte-identity has exactly zero
of that risk: any difference is a bug, full stop.** Then operators come out one at a time and each removal
carries its own number.

⚠ **What the rebuild contributes is its SHAPE and its derivation, not its code**: one `_hop` used by both
directions, a small explicit relayed state, and `η = λ − log(E_g/E_r)` as the frame-free composition
coordinate. Those move into the backbone. The file does not.

### §3.1 The backbone is three things, and the owner already named them

| | | lines today |
|---|---|---|
| the **accumulator**'s counts and densities | `substrate`, `region_arrays`, `node_geometry` | shared, untouched |
| **ψ**, the per-slot posterior | `simplex_logodds` | 755, untouched |
| **forward-backward** over the bipartite chain | `node_chain` + a NEW `sweep.py` | 136 + ~150 |

Everything else in `bp_solver.py` is **message-composition POLICY**, and policy goes behind an interface.

### §3.2 The interface, concretely

```python
class Policy(Protocol):
    def seed(self, ctx) -> State: ...
    """Each slot's own contribution to the relay, before any message arrives."""

    def hop(self, state_at_src, src: int, dst: int, ctx) -> Delta: ...
    """What crosses ONE step, and what variance it pays for crossing."""

    def fuse(self, state_at_dst, delta) -> State: ...
    """How an arriving delta combines with what is already there."""

    def deliver(self, fwd_at_left, bwd_at_right, i: int, ctx) -> PsiMessage: ...
    """The four psi channels at slot i — from the NEIGHBOURS' states only."""
```

`PsiMessage` is `(gdna_mode, gdna_prec, rna_mode, rna_prec, lam_mode, lam_prec, theta_mode, theta_prec)`
and nothing else. The backbone owns the loop, the write-back and the assertions; the policy owns
everything this campaign has been arguing about.

⭐⭐⭐ **Five defects become STRUCTURALLY IMPOSSIBLE rather than matters of discipline**, and each has
already shipped at least once:

| the backbone asserts | it would have caught |
|---|---|
| `deliver` is handed only the two NEIGHBOUR states — it cannot see slot `i` | **D4, which has recurred NINE times in different costumes** |
| every message mode lies inside its coordinate's own grid domain | **A16** — the tilt bug, 74 % of `g00`'s error, on its first run |
| every delivered share is in `[0, 1]` | the over-unit certified-RNA claim (85,477 fragments at a slot holding 5) |
| `\|T\| ≤ 3` | AXIOM 0, made executable |
| the write-back touches only `solvable` slots | the silent basis mismatch that made an A5 gate read `max\|Δ\| = 1.0` |

### §3.3 The order of work

| | step | acceptance |
|---|---|---|
| 1 | Write `sweep.py` (the loop, the write-back, the five assertions) and `policy/head.py` carrying every current operator behind a NAMED switch | ⛔ **`--arm backbone` must be BYTE-IDENTICAL to `base` on all 648 scored fields** (A5, the harness already does this for `zc_noop`) |
| 2 | Turn each switch OFF one at a time, one panel arm each — ⭐ and **score every arm PER STRATUM, never pooled**, because §2 shows the panel total hides a sign flip between strata | a number per operator. §1 predicts most are ~0 at pass-0, and §2 predicts several are NEGATIVE on three of the four strata |
| 3 | Delete every operator that measures inert, in one commit | `TRAPS.md` G3 — converge and delete |
| 4 | Only then re-open the message DESIGN, and do it in the REFIT where ψ has a population statistic | |

⛔ **Step 1 is not a rewrite and must not become one.** Its whole value is that byte-identity proves the
restructure changed nothing. Any temptation to "fix it while I'm in here" destroys that.

---

## §2b ⭐⭐⭐ THE MESSAGE-PRECISION SWEEP — NO PLATEAU, AND THAT IS THE ANSWER

The owner's account from shipping the production tool: *"messages do not need to be confident. When they
become confident they overwrite the strand-specific data. On unstranded data every node has precision
zero, so the weakest of weak messages still works, because it is more precision than zero."* That is a
quantitative claim with a shape — **one scalar on every message precision, swept.**

`abs_err_all_final`, relative to `base`, per stratum, `g00` excluded. `scale = 1` is `base`, `scale = 0`
is `msgfree_all`:

| stratum | 0 | 0.001 | 0.01 | 0.1 | **0.5** | 1 |
|---|---|---|---|---|---|---|
| stranded × capture ON | **−58.3 %** | −58.5 % | −57.4 % | −45.0 % | −18.4 % | 0 |
| stranded × capture OFF | **−43.7 %** | −43.6 % | −42.1 % | −32.4 % | −12.8 % | 0 |
| unstranded × capture ON | +154.8 % | +154.2 % | +140.0 % | +37.1 % | **−0.4 %** | 0 |
| unstranded × capture OFF | −32.1 % | −35.7 % | **−41.8 %** | −37.5 % | −16.8 % | 0 |
| **ALL** | +99.9 % | +99.2 % | +88.6 % | +16.1 % | ⭐ **−4.9 %** | 0 |
| ⛔ `g00` | +8,003 % | +6,462 % | +3,576 % | +146 % | +5 % | 0 |

⛔⛔ **THERE IS NO PLATEAU.** The four strata's own optima are at scale **0.001, 0, 0.5 and 0.01**. Between
0.01 and 0.5 unstranded × capture-ON gains 140 points while the two stranded strata give back 39 — a
genuine trade, in opposite directions, over the same parameter. **No single attenuation serves both
regimes, so the defect is not a mis-calibration of loudness.**

⭐⭐⭐ **AND THE STRANDED OPTIMUM IS ZERO, WHICH IS THE SHARP FORM OF THE FINDING.** If the messages were
merely over-confident, attenuating them would improve stranded *monotonically toward a nonzero optimum*.
Instead stranded is flat-best from scale 0 to 0.01 and degrades from there. **A message whose optimal
weight is zero carries no information — it is not a loud message, it is a wrong one.**

⛔⛔ **So this is a BIAS, not a precision — and `TRAPS.md` D1 says a variance cannot fix a bias.** D1 was
established three times independently before this; the sweep is its fourth demonstration, and it is why
three variance terms have failed to fix the one-sided bound in §3c. ⭐ **Fix the bias first. Attenuation
is a variance fix and cannot reach it.**

⚠ **One real, free result to keep:** `scale = 0.5` beats `base` on **all four strata** and −4.9 % overall,
with `g00` at only +5 %. That bounds a genuine precision debt of roughly 2× *on top of* the bias — but it
is a tuned constant (`CLAUDE.md` G1) and must not ship. It is a bound, not a fix.

## §2c ⛔ ONE CORRECTION TO THE PREMISE — it is not "unstranded", it is "unstranded × capture-ON"

Unstranded × capture-**OFF** is **−32.1 % better with no messages at all** (14/16 conditions). The message
layer's entire value sits in **one** stratum. ⭐ And the reason is legible: off capture the gDNA background
is uniform, so the fitted prior's single mode is already a good per-slot answer; under capture the
landscape is **bimodal at 2.98 decades**, so a slot's own density says nothing about its composition unless
you know which capture stratum it sits in — and that is what a neighbour can tell it.

⭐⭐ **So the message layer's job, stated as narrowly as the evidence supports: carry the LOCAL CAPTURE
LEVEL to a slot that has no own composition evidence.** That is much smaller than "carry composition", it
names a level rather than a share, and it is scoped by two predicates the code already computes
(`tau_lam == 0`, and the capture stratum the fitted landscape already resolves).

## §2d ⛔⛔ THE LENGTH CHANNEL CANNOT BE PRICED ON THIS PANEL — measure the substrate first

`length_likelihood` defaults **False** (`config.py:319`) and is the **only** channel that can give an
unstranded slot its OWN composition evidence: it is θ-independent, so the Schur complement that zeroes the
strand term at an AMBIG node does not apply to it (`simplex_logodds.py:519-524`), and `tau_len` enters
`tau_lam` **ungated** (`node_init.py:304-312`). The source says so itself: gating it *"would delete the
only evidence AMBIG nodes ever get."*

⛔ **But `TRAPS.md` F3: at equal component mean lengths the length channel carries EXACTLY ZERO information
about composition, at any depth — the 2×2 is identified only through `μ_g − μ_r`.** Measured on the ladder,
from the simulator's own post-capture truth:

| condition | `μ_gdna` | `μ_rna` | gap |
|---|---|---|---|
| `g50 ss0.50 capture_off` | 216.7 | 212.2 | **+4.5** |
| `g50 ss0.50 capture_on` | 240.4 | 236.9 | **+3.5** |
| `g75 ss0.99 capture_on` | 240.4 | 236.8 | +3.6 |
| `g01 ss0.50 capture_on` | 240.5 | 237.0 | +3.5 |

**A 1.5 % separation.** The ladder is near-blind to the length channel *by construction* — it was built
with equal configured fragment lengths. ⛔ **Enabling `length_likelihood` here would measure approximately
nothing, and that would be read as "the feature does not work."** `TRAPS.md` A8 exactly: prove the
benchmark can resolve the axis before running it.

⭐⭐ **And this re-reads §2's headline result.** On this panel an unstranded AMBIG slot has **no own
composition evidence available even in principle** — strand is dead at κ = ½, length is dead at
`μ_g ≈ μ_r`, and density needs the prior. So *"unstranded needs messages"* may be a property of the
**PANEL** rather than of the tool. On a library where gDNA and RNA lengths genuinely differ — which cfRNA
does — the length channel could carry that stratum with no messages at all. ⚠ **Unmeasurable until the
substrate exists**, and building it is one config line plus a re-simulation, not a solver change.

## §2e ⭐⭐⭐ THE ONE-SIDED RNA ARM — RIGHT ON THE CONTROL, PANEL-NEGATIVE, AND IT IS **D4j**

`simplex_logodds.ONE_SIDED_RNA` replaces `-1/2 p (log f - mo)^2` with `-1/2 p max(0, mo - log f)^2` on the
certified-RNA channel, in both ψ paths. ⛔ **A5 first: with the flag OFF the refactor is
1,728 / 1,728 scored fields IDENTICAL to the pre-refactor panel**, so anything below is the flag and
nothing else.

| stratum | base | `onesided_rna` | |
|---|---|---|---|
| **ALL** (g00 excluded) | 22,807,907 | 23,986,942 | ⛔ **+5.2 %** |
| stranded × capture ON | 3,469,926 | 3,508,787 | +1.1 % |
| stranded × capture OFF | 1,148,057 | 1,150,512 | +0.2 % |
| ⛔ unstranded × capture ON | 16,661,720 | 17,940,779 | **+7.7 %** |
| ⭐ unstranded × capture OFF | 1,528,204 | 1,386,864 | **−9.2 %** |
| ⭐⭐ **`g00` zero control** | 28,534 | **5,164** | **−81.9 %, 8/8** |

`library_f_gdna` at `g00`: `0.0015 → 0.0002`, `0.0008 → 0.0002`, `0.0002 → 0.0001` ×2, against a truth of
exactly 0.

⭐⭐ **The derivation is RIGHT — the zero control says so, 8 of 8, and it is the first mechanism in this
campaign the control has ever endorsed rather than refused.** The panel is negative anyway, and the reason
is legible and is a named trap.

⛔⛔ **`TRAPS.md` D4j: fixing one of two errors that CANCEL is worse than fixing neither.** The two-sided
RNA channel was doing **two** jobs: delivering the RNA lower bound (correct, and now one-sided) and, by
its *upward* side, acting as a **de-facto gDNA LEVEL channel** (accidental). Making it one-sided removes
only the second, and on unstranded × capture-ON that second job was load-bearing — because §3b's census
says **no NODE can receive a gDNA level at all**, so an accidental upward pull from the RNA message was
the only thing filling that hole. Remove the accident and the real gap is exposed: **+7.7 % on exactly
the stratum with no working level channel, and −9.2 % on the one where the prior already suffices.**

⭐ **So the one-sided term and the missing gDNA level channel are a cancelling pair and must be priced in
ONE arm.** D4j's own prescription, verbatim: *"when a fix is negative and its object sits on a multi-hop
path, price it in the arm that also removes the other defect. A pair of defects that cancel is one
experiment, not two."* ⛔ Do not land the one-sided term alone, and do not discard it — it is the only
correction the `g00` control has endorsed.

## §2f ⚠ STORAGE — 30 G of the 67 G is FASTQ that no calibration instrument reads

| panel | conditions | FASTQ | BAM |
|---|---|---|---|
| `ladder` | 37 | **23 G** | 26 G |
| `pilot` | 9 | **7 G** | 8 G |

⛔ **The FASTQs are NOT dead** — `scripts/benchmarking/{aligner,tools,analysis,config}.py` read them for
the third-party tool comparison, and `sim/orchestrator.py:284` uses `sim_R1.fq.gz`'s existence as its
`skip_existing` key, so removing them makes the orchestrator re-simulate. ⚠ **Not deleted, and this is an
owner call** — but no instrument under `scripts/design/` or `src/rigel/calibration/` opens one, so they
are archivable for calibration work.

⭐ **For the NEW fragment-length arm the saving is available with no such trade.** The arm needs the two
directions the ladder lacks — `gDNA << RNA` and `gDNA >> RNA`; `gDNA == RNA` is the ladder itself — and it
does not need the full 36-cell cross. Six conditions at the owner's own focus cell (`ss_0.50 × capture_on`,
plus the capture-OFF twin and one stranded control per direction) at a mid gDNA rate where both components
are present in quantity. **6 × 806 M ≈ 4.8 G BAM-only, against 9 G with FASTQs**, and an
`emit_fastq: false` option in the orchestrator is a small change that makes every future calibration panel
half the size.

## §3a ⭐⭐ THE MINIMAL RELAYED STATE IS **FOUR** ARRAYS, NOT TEN

`_relay` carries ten: three component densities, three mode-fusion precisions, three measurement
precisions, and `tau`. ⭐ **ψ never sees the mode-fusion triple** — `cpg`/`cpp`/`cpn` set the fused
density `cg`/`cp`/`cn` (`bp_solver.py:1506-1508`) which becomes the MODE at `:1511-1513`, while the
precision ψ receives is the separate measurement stream `cm_g`/`cm_p`/`cm_n` (`:1509`, `:1545-1550`). So
the honest minimum is **three levels plus one precision triple**, with `tau` a fifth array only if the λ
channel survives (its Schur gating at `node_init.py:296-297` is structural and not derivable per hop).

**The policy signature, and the part that matters is the contract:**

```python
def policy(src, dst, rho_src, prec_src, ctx) -> (rho_msg, prec_msg):
    """Return the message in the DESTINATION's frame.

    CONTRACT (D4), enforceable by construction: index `ctx` at `src` freely, and at `dst`
    ONLY through `ctx.obs` and `ctx.geom`. Never `ctx.belief`, never relay state at `dst`.
    """
```

⭐ `StepContext` splits its fields into three headings — **observations** (readable at either end),
**geometry/structure** (readable at either end), and **source-side only** (`belief_fg`, `own_rho`,
`own_prec`). That heading is what turns D4 from a discipline into something an assertion can catch.

## §3b ⛔⛔ THE OPERATOR LEDGER — five families, inventoried and adversarially verified

⚠ **Where the inventory and its verifier disagree, the verdict is UNKNOWN and it needs an arm.** That rule
moved three operators out of CUT.

| operator | verdict | the evidence that decides it |
|---|---|---|
| the scan skeleton, the forward and backward calls | **BACKBONE** | the only thing in the family that is unambiguously backbone; the whole direction dependence reduces to which neighbour array is read |
| `_fuse` (inverse-variance in linear density space) | **BACKBONE** | the fusion rule; the near-zero pass-through argument is sound |
| the write-back on `solvable` | **BACKBONE** | ⭐ "the messages reach the answer only as a mode and a precision handed to ψ" — state this as a design invariant |
| `count_logvar = trigamma(n+½)` | **KEEP** | the only operator with a panel-wide accuracy number: **39 %** |
| `own_precision` · `g1_locked` · `terminus_flank_gain` | **KEEP** | measured, one home each, belief-free |
| `residual_level` | **KEEP, PROMOTE** | D4-clean, and the only estimator of its kind |
| ⭐⭐⭐ the certified-RNA channel's **SIDEDNESS** | **RE-DERIVE** | **the single highest-value change any family admits** — see §3c |
| `node_total_density` | **RE-DERIVE** | ⛔ it reads the DESTINATION's own `f_g`, and slots where a *solved* belief (not the `{0,0,1}` default) sets the frame carry **57–77 % of library mass**. D4, live in the tree |
| the reframe `r` · `r_g` | **RE-DERIVE** | keep the rule ("a gDNA level crosses unscaled"), re-derive the scope |
| `own_composition_logvar` | **KEEP fn, RE-DERIVE mask** | the three-state logic has no substitute; `struct_lock` is the wrong mask (A11b) |
| the `_damp`/`_damp_v`/`_dv`/`_dv_arr` inlines | **RE-DERIVE** | one named function per twin, delete the inlines |
| `peel_rna_logvar` | **CUT** | no production consumer — `CLAUDE.md` G3 exactly |
| `_pin_v` | **CUT** | its sole consumer is the DL gap |
| P1e, the DerSimonian–Laird conservation term | **CUT** | ⭐ its own comment: *"PARTLY A DEBT — THIS PRICES A BIAS AS A VARIANCE … a variance cannot move a mode toward truth"* |
| the flank pair `rho_lo`/`rho_hi` · `_flank_dom` | **CUT with `r`** | ⚠ it is a **100×–10,000× correction on a third of live EDGEs** — but the panel measured it node-axis NEGATIVE (**+0.5 % mwae, +36.9 % confidently-wrong**), and ⛔ **that number appears in no source file**, so a reader of the source alone would believe it is panel-positive |
| `graft_premise_logvar` | **CUT from backbone** | re-derive per class only when TSS/TES land |
| ⭐ `transfer_logvar`'s non-graft branch | **UNKNOWN** | *"the most important UNKNOWN in this family"* — nothing has ever measured deleting it alone |
| `mismatch_deflate` | **UNKNOWN** | ⚠ the inventory said CUT on a coverage number the verifier proved **INVERTED**; it is live on ~82 % of stranded mass |
| `_lend` | **UNKNOWN** | inventory RE-DERIVE, verifier UNKNOWN — "no measured benefit of its own" is false |
| `framed`'s fallback · `peel` · `peel_share_logvar` | **UNKNOWN** | each needs its own arm; `framed`'s question narrows to the graft subset |

⚠ **Two facts the ledger turns up that nothing in the tree records.** `spliced_count` is deliberately
excluded from the reframe's total, and the excluded component measures **1.69× the retained total** on
50–60 % of live-EDGE mass at low gDNA — not a rounding decision. And the certified-RNA bound reaches ψ
**twice**, both times two-sided: once as `rna_imp`, and again through `tlam` → `lam_imp`
(`bp_solver.py:1474-1479`, `:1525`).

## §3c ⭐⭐⭐ THE ONE CHANGE TO MAKE FIRST — the certified-RNA channel is a BOUND delivered as an EQUALITY

`bp_solver.py:439-447` states it in the source, verbatim: *"what the graft actually knows is an
INEQUALITY, `ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)` … and it uses it as an equality."* The chain from that
inequality to ψ passes through **three operators that price it as a VARIANCE** — `graft_frame_logvar`,
`graft_premise_logvar`, and P1e — and through **zero that price it as a DIRECTION**. P1e's own comment
says why that cannot work.

ψ applies `-1/2 * p * (log f - mo)^2` in **both** code paths (`simplex_logodds.py:531-540` and
`:315-324`). **There is no hinge, no truncation, no one-sided term anywhere in the file.** A destination
holding MORE RNA than the bound — which the inequality explicitly permits — is penalised exactly as hard
as one holding less. Measured consequence at `g01 ss0.50 capture_on`: HEAD's self-solve is 0.0086 against
a truth of 0.0023, the message drives it to **0.3219** at precision 327, and muting the channel returns
**0.0005**.

    psi += -1/2 * p * max(0, mo_p - log f_active)^2

One expression, both paths, **no new constant**, and — unlike all eleven refused candidates — it adds
doubt in **only one direction**, so the `g00` control cannot refuse it: there `f_active = 1` satisfies
every bound and the channel goes inert. ⚠ If it lands, re-price `graft_premise_logvar` immediately: a
one-sided factor and a premise-variance debt are two compensations for one defect and would double-charge.

## §3d ⭐⭐⭐ THE ARCHITECTURE — the backbone is a FOLD, and every argument lives in ONE function

The elegance is not decoration; it is the mechanism by which the next session cannot repeat this one's
mistakes. Three files, and a new reader needs only the first two.

```
src/rigel/calibration/
  sweep.py                 ⭐ THE BACKBONE. ~150 lines. Two directional scans, one combine, one
                              write-back, five assertions. It knows nothing about capture, grafts,
                              reframes, pins or enrichment — those words do not appear in it.
  messages/
    __init__.py               the Policy protocol, StepContext, PsiMessage.  ~60 lines.
    silent.py              ⭐ SilentPolicy — sends nothing. THE DEFAULT. ~5 lines.
    head.py                   HeadPolicy — every current operator, each behind a NAMED switch.
    variance.py               what was `enrichment_frame.py`: the policy's toolbox, not the backbone's.
```

⛔ `bp_solver.py` is deleted the moment `HeadPolicy` is byte-identical. Not before, not gradually.

**The whole backbone, in the shape it should be written:**

```python
def solve_chain(chain, ctx, policy) -> NodeBelief:
    own = build_node_init(...)                            # unchanged, backbone
    fwd = _scan(range(n),           chain.left,  own, policy, ctx)
    bwd = _scan(reversed(range(n)), chain.right, own, policy, ctx)
    msg = policy.deliver(fwd, bwd, chain, ctx)            # sees ONLY neighbour state
    dc  = psi(ctx.obs, msg)                               # unchanged, backbone
    return _write_back(dc, ctx.belief, ctx.solvable)
```

⭐⭐ **THE FIVE ASSERTIONS LIVE IN THE BACKBONE, NOT IN THE POLICY — that is the whole design.** A future
policy can be as wrong as it likes and still cannot commit any of these, each of which has shipped at
least once:

| the backbone asserts | it would have caught |
|---|---|
| `deliver` is handed only the two NEIGHBOUR states | **D4 — nine recurrences in nine costumes** |
| every message mode lies inside its coordinate's own grid | **A16** — the tilt bug, 74 % of `g00`'s error |
| every delivered share is in `[0, 1]` | the over-unit certified-RNA claim (85,477 fragments at a slot holding 5) |
| `\|T\| ≤ 3` | AXIOM 0, made executable |
| the write-back touches only `solvable` slots | the basis mismatch that made an A5 gate read `max\|Δ\| = 1.0` |

**And the interface is one function with one contract:**

```python
def policy(src, dst, rho_src, prec_src, ctx) -> (rho_msg, prec_msg):
    """CONTRACT (D4): index `ctx` at `src` freely, and at `dst` ONLY through `ctx.obs`
    and `ctx.geom`. Never `ctx.belief`, never relay state at `dst`."""
```

`StepContext` splits its fields under three headings — **observations** (either end),
**geometry/structure** (either end), **source-side only** (`belief_fg`, `own_rho`, `own_prec`). That
heading is what turns D4 from a discipline into something the backbone can check.

⭐ **Why `SilentPolicy` is the default and why that is the point.** A new session reads `sweep.py` plus
five lines and holds the entire working system in their head. `head.py` is opt-in, clearly labelled as the
legacy arm being dismantled, and every operator inside it is a switch — so `ladder_arm_ab.py` prices them
**one at a time** instead of as a block, which §3b says is exactly what is needed ("every single ablation
is small and the joint one is large").

**Two acceptance gates, both A5, both cheap:**

| | must be |
|---|---|
| `--arm backbone_head` (all switches on) | **byte-identical to `base`** — proves the restructure changed nothing |
| `--arm backbone` (SilentPolicy) | **byte-identical to `msgfree_all`** — proves the off-switch is the off-switch |

⛔ **Step 1 is a restructure and must not become a rewrite.** Its entire value is that byte-identity proves
it changed nothing; any temptation to fix something while in there destroys the proof and puts us back
where the `eta` rebuild was — a change whose sign nobody could attribute.

### §3e What is deleted, and when

| | when | why |
|---|---|---|
| `scripts/design/eta_node_sweep.py` | ⭐ **DONE, 2026-08-07** | the prototype is superseded; its two findings are `TRAPS.md` A16/A17 and §2/§3 here |
| `--arm eta` / `--arm eta_nolevel` | ⭐ **DONE** | no consumer |
| `tests/calibration/_eta_reference.py` + `test_eta_transfer.py` | **with the backbone commit** | they gate η algebra HEAD does not use, so they gate dead code — but three of their derivations become the backbone's assertions, and deleting the gates before their replacement exists is a strictly worse state for exactly one step |
| `src/rigel/calibration/bp_solver.py` | when `HeadPolicy` is byte-identical | one commit, no overlap period |

## §4 ⭐⭐ WHY THE TWO CORNERS PULL APART, AND WHY THAT STOPS MATTERING

The owner's diagnosis, and it is the right one: a zero-gDNA library is solved by making the gDNA messages
strong enough to pull everything to zero; a high-gDNA library is solved by damping them. **Each corner has
a trivially-satisfying rule and the two rules are opposite.** Every mechanism this project has priced
solved one corner and was refused by the other — eleven of them, from three independent mathematical
directions (`variance_model_notes.md` §6).

⭐ **The reason is structural and it is already written down:** pass-0 is forbidden a population statistic
by construction (A2 — that is the circularity), and **the quantity that separates the two corners IS the
library's gDNA content.** So no prior-free, library-agnostic rule can tell "the discrepancy is RNA" from
"the discrepancy is capture-enriched gDNA". They are observationally identical.

⭐⭐⭐ **And §1 says that is fine, because pass-0 does not have to answer.** If the substrate's job is to
train a prior, then the honest output at an object with no own evidence is **no claim at zero weight** —
not a guess relayed from a neighbour. The prior is then fitted on the objects that DO have evidence, it
learns the library's gDNA content from them, and the refit hands ψ the population statistic that pass-0 was
never allowed to have. At that point ψ *can* separate the corners, because it finally knows which one it
is in.

⛔ **So the corner problem is not a modelling problem to be solved at pass-0. It is a stage-ordering
problem, and the current code tries to solve it one stage too early.** That is the single sentence this
consolidation is for.
