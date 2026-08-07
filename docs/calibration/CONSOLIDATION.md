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
