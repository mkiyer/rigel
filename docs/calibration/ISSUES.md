# CALIBRATION ISSUE REGISTER — what is open, what it costs, and how sure we are

    Working doc for the current calibration phase. ⚠ **NOT one of the six permanent docs.**
    ⭐ **This is the register; `ROADMAP.md` §2 stays the single ranked list.** When an issue closes,
      promote its number to `ROADMAP.md`, its lesson to `TRAPS.md`, its derivation to `EQUATIONS.md`,
      and strike it here. When the phase ends, fold whatever is left into `ROADMAP.md` and delete this.
    ⛔ Every "cost" below is measured on the 36-condition gDNA ladder unless it says *unmeasured*.
      A cost that says unmeasured must be measured before the issue is worked (`TRAPS.md` B1).

---

## The register

| id | issue | measured cost | diagnosis | fix cost |
|---|---|---|---|---|
| ✅ **C1** | ~~**The relay pins evidence-free deep exons to `f_g ≈ 0.85`**~~ **CLOSED 2026-08-04** — the reframe is a composition imputation and needed a **licence**, not a better scale. `EQUATIONS.md` §3.5, `TRAPS.md` D4c | node/exon mwae **0.3002 → 0.2019** (`g25 off`) and **0.3850 → 0.2686** (`g25 on`); Σ\|err\| −345,890 and −613,560 | ⭐⭐ **root-caused, derived, gated** (`tests/calibration/test_gdna_scale_rule.py`, 6 gates, perturbation-verified) | landed |
| ✅ **C11** | ~~**THE MASS PIN drives the delivered gDNA FRACTION to ½**~~ **CLOSED 2026-08-04** — the pin now fires iff no BELIEF can reach its budget: the composition was supplied, **or** the destination is structurally pure gDNA. `EQUATIONS.md` §3.5c, `TRAPS.md` D4d/D4e | ⛔⛔ **the CEILING says it cost the PANEL nothing**: deleting the pin outright moves the ladder **+0.0002** (worse) with confidently-wrong **+13 %**. Landed on the derivation, and because it is free — 22/32 rows better, **none worse**. `nested_exons` toy 0.2618 → 0.2264, and → **0.0541** with C4 | ⭐⭐⭐ **derived, closed-form, gated** (`test_relay_mass_pin.py`, 7 gates + 1 strict xfail, perturbation matrix in the docstring) | landed |
| ✅ **C4** | ~~**TSS/TES bits are plumbed and unconsumed**~~ **CLOSED 2026-08-04** — they are the POPULATION conjunct of the composition licence: `T(EDGE) = T(left) ∩ T(right)`, so a terminus makes one flank's population larger and the imputation across that pair is withheld. `EQUATIONS.md` §3.5b, `node_geometry.terminus_flank_gain` | ⭐⭐ ladder: confidently-wrong **21,154 → 20,173 (−4.6 %)**, overconfidence **2.788 → 2.688**, `weak%` **2.91 → 2.80**, mwae 0.0415 → 0.0413. Toy **0.2264 → 0.0541**, mirror **0.2145 → 0.0535** | ⭐⭐ **gated on a ± MIRROR annotation and against a second algorithm** (`test_terminus_population_licence.py`, 11 gates) | landed. ⛔ **Termini only**; DONOR/ACCEPTOR is a separate experiment (see `ROADMAP.md` step 4b) |
| ⭐⭐ **C2** | ⛔ **RESTATED AND PROMOTED 2026-08-04 — it is a SOLVER defect now, not a reporting one.** κ = ½ means *exactly* zero strand information, but κ̂ = 0.500689 leaves `I(f_g) ∝ (2κ−1)²` at 1e-6…1e-4 — and the composition licence is a **precision** test with no floor, so those precisions **grant the imputation licence** at every evidence-free exon | ⭐ **measured**: on `nested_exons` the licensed mass pin still fires through the whole gene interior — 0.2264, where switching the pin off entirely reaches 0.0760. Classification cost as before (hid C1 entirely, 1.06 M fragments) | ⭐⭐ **complete and quantified**; the new half is `TRAPS.md` D4f | ⛔ **NOT a floor** — B11 implemented and refuted that. Propagate `Var(κ̂)` into the strand arm's precision; it drives that arm to ~0 by derivation. ⭐ Pays for C7 at the same time |
| **C3** | **The graft's frame is wrong under hybrid capture** (M8/P1d) | deleting the junction channel *improves* `node/exon` **0.0431 → 0.0317 (−25 %)** under capture, and wrecks it off capture | ⚠ partial — the term is known and self-documented as over-stating (median φ = 2.45) | medium |
| **C4b** | **`graft_premise_logvar` per structural class** — what is LEFT of C4 | *unmeasured* | n/a — a build. It is an admitted "**DEBT, not a model**" that "splits ≥30× on whether the boundary carries a transcript TERMINUS", and it now has that bit | medium. ⛔ Its scope is DONOR/ACCEPTOR, which C4 deliberately left alone |
| **C5** | **Reframe component-set mismatch** (the junction leak) | ⭐ **ceiling ≈ 0**: best case −6.8 % on one class, worse on 2 of 4 conditions | ⭐⭐ complete | small — ⛔ **but not worth doing alone** |
| **C6** | ⭐ **`intergenic\|exon` EDGEs as a gDNA reference** — ⚠ **restated 2026-08-04.** The *landscape* half is answered: they already carry it, through their own measurement plus the mass pin. What remains is that a G1 **EDGE** has `prec_g = 0`, so it can only RELAY a gDNA level, never ORIGINATE one — and under capture the intergenic NODEs it therefore depends on are depleted ~24× while the EDGE itself carries 8.6× more | **8.6× more mass per object under capture** (139.7 vs 16.2) and truth **exactly 1.0**. ⭐ On the capture-ON toy the EDGE→exon message is **dead** (`cm_g = 0`) at low RNA, because the intergenic NODE behind it has ~0 counts | ⭐ mechanism identified: `struct_lock` is NODE-only, so `own_composition_logvar` hands a G1 EDGE `Var = ∞`. ⭐⭐ **And §4 below adds the number that makes it a PAIR problem**: on the toy the EDGE carries 23 counts and the NODE beside it 60, so neither is the answer and the inverse-variance fuse of the two is | small; ⛔ but `strand_evidence`'s docstring argues for the NODE-only scope and that argument must be re-examined for an `intergenic\|exon` EDGE specifically — RNA cannot cross a gene boundary, so the fragments spanning it are **not** RNA-contaminated |
| **C7** | **Declared precision is not earned where strand is live** — up to **8.9×** overconfident | measured on both stranded strata | ⚠ not diagnosed | unknown |
| **C8** | **The relay degrades the solvable set on 26/32 contaminated conditions** (mean Δ +0.0027) | measured | ⚠ not diagnosed; likely overlaps **C1** | unknown |
| **C9** | **The fitted fragment-length gap flips sign under capture**, which is why the length channel is gated off | previously "36.2 % of pass-0's error" — ⛔ that share was computed over the population **C2** mis-defined and must be re-derived | partial | large |

---

## Priority — OWNER-SET, 2026-08-04

| | | |
|---|---|---|
| ✅ **1** | ~~**C1** the relay pinning exons~~ | ⭐⭐ owner: *"seems to be a HUGE bug and source of error. HIGHEST PRIORITY to dissect."* — **CLOSED 2026-08-04.** The derivation is `EQUATIONS.md` §3.5, the lesson is `TRAPS.md` D4c, the gates are `tests/calibration/test_gdna_scale_rule.py` |
| ✅ **1=** | ~~**C11** the mass pin~~ · ~~**C4** TSS/TES~~ | **BOTH CLOSED 2026-08-04** — one licence, two conjuncts. `ROADMAP.md` §2c has the numbers. ⛔ The ceiling run says C11 cost the **panel** nothing (+0.0002 to delete the pin outright), so it landed on the derivation and on being free; C4 is what moved the aggregates |
| ⭐⭐ **1** | **C2** — ⚠ **promoted from `deferred`, and the owner's objection no longer covers it** | The objection was about the *solver being inert* to a τ ≈ 1e-7 statement, and it was right for a fused answer. But a **licence** is not a fused answer: it is a boolean, and a 1e-6 precision flips it exactly as hard as a 1e+6 one. The composition licence now turns on those precisions (`TRAPS.md` D4f). ⛔ This is the one place a threshold would be tempting and must not be added |
| **1=** | **C3** the graft frame | ⭐⭐ owner: *"another HIGH PRIORITY (either 1st or 2nd). the new accumulator structure and TSS/TES plumbing needs to be incorporated"* — the TSS/TES half of that is done (C4); what remains is the graft's frame under capture, and `graft_premise_logvar`'s DEBT now has the bit it was waiting for |
| **2** | **C5** the reframe component-set mismatch | ⭐ owner: *"HIGH priority. Clearly a major source of error."* ⚠ **My measurement disagrees** and the disagreement is recorded, not resolved — see C5 below |
| **3** | **C6** intergenic↔exon as a gDNA reference | ⭐ owner: *"yes, we must use this as a measure of gdna density. The intergenic NODE ↔ intergenic-exon EDGE can be solved. They can be assumed to be pure gDNA and have the same composition (no RNA, pure gDNA)."* — note this makes it a **solvable pair**, not just an anchor |
| ~~deferred~~ | ~~**C2** the strand deadband~~ | ⛔ owner: *"likely NOT an issue. IF you weight by precision, any of these nodes with miniscule precision should not contribute to error/overconfidence. You would need to prove that this is high priority."* ⚠ **Moved to step 1 above on 2026-08-04, and the proof the owner asked for is the licence**: a boolean gate is not precision-weighted, so a τ ≈ 1e-6 statement flips it at full strength. The note below is still correct about the *fused answer* |
| **later** | **C7, C8, C9** | ranked on numbers the C2 classification mis-defined; re-derive before ranking |

### ⚠ C2 — the owner's objection, and what it does and does not settle

The objection is correct on the solver: a node with τ ≈ 1e-7 has a posterior variance so wide that
precision weighting makes it inert, so it should not move a fused answer or an honest overconfidence
figure. ⭐ **I have not proved otherwise and am not claiming it.** What the leak *did* do is corrupt the
**classification** — the `τ > 1e-9` test promoted those objects into the scored population and thereby
hid a 1.06 M-fragment error (C1) from every report. **That consequence is already fixed** by the sd(λ)
curve and the undetermined-overreach table, neither of which asks whether τ is above a threshold.

⛔ So the open part of C2 is only: *does the leaked `disc` change any answer?* That is a one-arm
measurement (inject κ = 0.5 exactly and diff), it is cheap, and it should be done before any build.
⚠ Owner ruling on the shape, if it ever is built: **keep it a continuous shrinkage — strandedness
varies along a spectrum** — so option (c), the mixture with an atom at κ = ½, is ruled out.

---

## ⭐⭐ THE TOY HARNESS — built 2026-08-04, `scripts/design/toy_harness.py`

A mini chromosome you define, simulated + scanned + calibrated in **0.1–5 s**, with every object's
answer printed beside per-object truth. The library-level quantities a toy cannot fit are **harvested
from a real cached condition** and injected (`InjectedCalibrationPriors`, whose docstring already
specified this use): κ, both overdispersions, the Fisher noise-floor sample sizes, the enrichment
NPMLE, the intron background, ρ_bg, both fragment-length pmfs, capture on/off and its knobs, and the
simulation's own length/strand parameters.

⭐⭐ **The gDNA level is DERIVED, not a knob.** The injected enrichment landscape is an *absolute*
log-density model, so a toy at the wrong depth is a different library, not a small one. The harness
measures the donor's gDNA molecules-per-base on the donor's own structurally-pure-gDNA nodes
(`Σcount / ΣE`, `EQUATIONS.md` §7.1) and simulates the toy to match it. Gated on the density the toy
**realises**, with a 10× perturbation. `ToySpec` deliberately has no gDNA field.

Gates: `tests/calibration/test_toy_harness.py` — 7, each with its own perturbation.

### ✅⭐⭐⭐ WHAT IT FOUND IMMEDIATELY — C1's mechanism candidate, and it is now FIXED

⭐ **Kept as the record of a finding and its resolution.** The dependence below is **gone** since
`EQUATIONS.md` §3.5: measured **0.0107 with a pure-gDNA intron vs 0.0112 with a nascent-bearing one, a
factor of 1.04**, down from 1.25 on this donor and 31 on the ladder donor. ⭐⭐ And it was fixed *without
being aimed at* — a factory-solved intron reports `f_g ≈ 1`, hence ~zero RNA density, hence zero RNA
precision, so it cannot lend a composition and its gDNA level crosses unscaled. The toy gate that pinned
the defect **failed because the defect was gone** and has been rewritten from an ORDERING assertion to an
INDEPENDENCE one.

`two_exon` under the `g25 ss0.50 capture_off` donor, 9 objects, 1.6 s: **99.95 % of the error is on the
two exon nodes** (truth f_g 0.004, predicted 0.156). Then one field swept — the nascent-RNA level, which
changes what is inside the **intron** and nothing about the exon:

| nrna | intron true f_g | **exon pred** | exon truth | toy mwae |
|---|---|---|---|---|
| 0 | 1.000 | ⛔ **0.156** | 0.004 | **0.1523** |
| 1 | 0.567 | 0.022 | 0.004 | 0.0177 |
| 5 | 0.221 | 0.008 | 0.004 | 0.0045 |
| 20 | 0.056 | 0.005 | 0.003 | 0.0017 |
| 40 | 0.040 | **0.005** ✅ | 0.005 | **0.0005** |
| 100 | 0.025 | 0.005 | 0.004 | 0.0012 |

⭐⭐ **The exon's answer tracks the INTRON's composition, monotonically, while the exon's own truth
stays flat at 0.003–0.005.** Total error collapses **300×** the moment the intron carries any RNA, and
7× at `nrna = 1`.

**The reading:** the exon has no own evidence (κ = ½, length channel off), so its answer comes from its
neighbours, and its only non-intergenic neighbour chain runs through the intron. The imputation premise
is "adjacent objects share a composition" — and a **pure-gDNA intron beside a 99.6 %-RNA exon is that
premise at its most false**. ⭐ That reading was right, and `EQUATIONS.md` §3.5 makes it exact: the reframe
delivers `phi_c(src)·rho_tot(dst)`, so a `phi_g(src) = 1` source hands the exon its own total.

⚠⚠ **The MAGNITUDE is donor-dependent and the harness surfaced that too.** Against a six-gene synthetic
donor the same contrast was 0.209 vs 0.167 — direction preserved, factor 1.25 not 31. A small donor
cannot determine the enrichment landscape or the intron background, so its injected globals are a
different regime. ⛔ That is why the gate asserted only direction and ordering, never the size — and it is
why the rewritten gate asserts a RATIO near 1 plus "both arms small", never an absolute figure.

---

## C10 — ⚠ THE LADDER PANEL MAY BE MEASURING THE WORST CASE, NOT THE REAL ONE

Every one of the 36 ladder conditions is `nrna_none`. The sweep above says that is precisely the regime
where the intron↔exon imputation premise is **maximally violated**, and where the panel's dominant error
mode lives. Real libraries carry pre-mRNA. So the panel may be over-stating this error mode by a large
factor — or it may be correctly testing a robustness corner. ⭐ **Either way it is not currently known
which**, and the fix is one config axis: add a nascent rung. Cost: *unmeasured*, but the toy says the
effect size is up to 300×, which makes it the cheapest large piece of information available.

---

## ✅ C1 — CLOSED 2026-08-04. What it was, in one paragraph

The reframe ``r = rho_tot(dst)/rho_tot(src)`` delivers ``phi_c(src)·rho_tot(dst)`` — the source's density
SHARE applied to the destination's observed total — so it is a **pure composition imputation with no level
transport in it at all**. At an `intergenic|exon edge -> EXON` step the source is structurally 100 % gDNA, so
``phi_g(src) = 1`` and the delivered gDNA level collapsed to the destination's own total: measured **28.684
against a true 0.0257**, which is ``M/E_g`` to the digit, and the exon's answer was therefore **flat at
`f_g` = 0.90 across a 10,000x RNA sweep**. ⭐ There is no corrective factor, because the missing factor is
``phi_g(dst)``, the destination's own belief. The fix is a **LICENCE** — the reframe is allowed exactly
where a composition imputation is allowed (the λ-emission gate's own predicate, which existed and was
applied to the τ stream only) — and where it is not, the gDNA **level crosses unscaled**.

⭐ **Full derivation: `EQUATIONS.md` §3.5. Lesson: `TRAPS.md` D4c. Gates:
`tests/calibration/test_gdna_scale_rule.py`** (6, perturbation-verified, one perturbation deliberately
recorded as firing nothing with its reason).

### ✅⭐⭐⭐ C11 + C4 — CLOSED 2026-08-04. One licence, two conjuncts.

⚠ Vocabulary is `DESIGN.md` §0: **NODE** = contiguous interval, **EDGE** = single position between two
NODEs, **step** = one adjacency move, **the mass pin** = the operator below, **structurally pure-gDNA
object** = a slot where no RNA strand is admissible (an intergenic NODE, or an `intergenic|exon` EDGE).

⭐ **Derivation `EQUATIONS.md` §3.5b/§3.5c · lessons `TRAPS.md` D4d/D4e/D4f + B14 · numbers
`ROADMAP.md` §2c · gates `test_relay_mass_pin.py` (7 + 1 strict xfail) and
`test_terminus_population_licence.py` (11), each carrying its perturbation matrix.**

What C11 WAS, in one paragraph: the pin restores `Σ_c ρ_c·E_c = M` with `k = M/S`, where `S` fills every
component the message did not supply from the destination's **own** density. Substituting its own
definitions gives `k = 1/(φ_msg + R_own)` — a saturating map with fixed point `(1−R_own)·ρ_tot`, and
`R_own`, the RNA share of the destination's own self-solve, is **exactly ½** at an evidence-free object. So
it drove the delivered gDNA FRACTION to ½ regardless of the truth; at `R_own = 0` it collapsed to the
destination's own total with the incoming level cancelled algebraically (`TRAPS.md` D4). On the
`nested_exons` toy the delivered level rose **0.071 → 0.102 → 0.097 → 0.162 → 0.192 → 0.215** against a
uniform truth of 0.0769, 2.8× at the deepest object.

## ⛔⛔ 1. THE THING WORTH KEEPING: THE CEILING SAID C11 COST THE PANEL NOTHING

Deleting the pin outright and running all 36 ladder conditions moves the solvable mwae **0.0415 → 0.0417
(worse)** — better on 18/32, worse on 14, worst row +0.0108 — while confidently-wrong Σ|err| goes
**21,154 → 23,831 (+13 %)** and the declared-precision ratio **2.79 → 3.14**. The toy's −71 % lives almost
entirely in the population `solvability_audit.py` deliberately excludes.

⭐ So the licence landed on **the derivation plus being free** (22/32 rows better, none worse), not on the
toy number. `TRAPS.md` B1 again, and this is the fourth time the ceiling discipline has re-ranked something.

## ⭐⭐ 2. AND THE CEILING RUN PAID FOR ITSELF TWICE: THE PIN CARRIES THE CAPTURE LANDSCAPE

Gating the pin by the composition licence **alone** — Step A exactly as specified — delivers the **off-probe
intergenic floor** to every exon under capture, firing `test_gdna_scale_rule`'s capture gate at 20× and
200×. The reason is C6: an `intergenic|exon` EDGE measures the gDNA density at its own capture stratum but
carries `prec_g = 0`, so it cannot ORIGINATE a level through the fuse — the pin was its only channel.

⭐ The resolution is not an exception but a sharper statement of D4: **no BELIEF may enter the pin's
budget**, which holds in two states — the composition was supplied (nothing is filled in), or the
destination is structurally pure gDNA (nothing *to* fill in, and `M/E_g` is then an OBSERVATION of the very
quantity the message is about, which D4 permits). `TRAPS.md` D4e.

## ⚠ 3. WHAT REMAINS, AND IT IS C2

The delivered level is now exact on a synthetic uniform field (1.000 at every slot). On the real toy the
pin still fires through the gene interior, because `lend` is a **precision** test with no floor and
κ̂ = 0.500689 leaves the strand arm at 1e-6…1e-4 rather than 0. That is C2 reaching the solver
(`TRAPS.md` D4f), and it is `ROADMAP.md`'s new step 1.

⭐ Two residuals are now separated, each with its own gate: the **level** is exact, and the RNA-free
interior exon still reads `f_g` **0.914** against a truth of 1.000 — ψ's uninformative reference at the
vertex, pinned as a strict xfail against the TRUTH rather than widened to the number reached.

## ⭐ 4. ONE MEASUREMENT FOR C6, TAKEN ON THE WAY

On the `nested_exons` toy the `intergenic|exon` EDGE carries **23** counts while the intergenic NODE beside
it carries **60**, and the EDGE's own density reads 0.107 against a truth of 0.0769 (Poisson sd ~21 % on 23
counts). The pin hands the whole gene the EDGE's noisier number. ⭐ That is the owner's own C6 framing —
*"the intergenic NODE ↔ intergenic-exon EDGE can be solved … assumed to be pure gDNA and have the same
composition"* — as a **solvable pair**: inverse-variance fusing the two beats either alone, and it is the
same `f_g = 1` structural certainty this licence already relies on.

## C2 — the strand deadband leaks on unstranded libraries

**What it is for.** On an unstranded library the strand channel carries *exactly* zero information about
composition (`EQUATIONS.md` §5.2: `I(f_g) ∝ (2κ−1)²`, and κ = ½). κ is fitted, so it never lands exactly
on ½, and the deadband exists to suppress that residue — deliberately as a smooth shrinkage rather than
an on/off switch:

    disc = 4 · max(0, (κ̂ − ½)² − σ²_d)          σ²_d = ¼(1/N_rna + od_r) + ¼(1/N_gdna + od_g)

⭐ **`σ²_d` is exactly right**: for a proportion fitted on `N_rna` spliced observations, `¼/N_rna` *is*
`Var(κ̂)`. So the mechanism subtracts precisely the noise it should. The design is sound.

⛔ **The defect is the `max(0, ·)`.** Subtracting the mean of a noisy square and then truncating at zero
turns mean-zero noise into a **one-sided positive leak**: with κ truly ½, `(κ̂−½)²/Var` is a χ²₁, and

    E[ max(0, χ²₁ − 1) ] = 0.4839      ⇒      E[disc | truly unstranded] = 4 × 0.4839 × Var(κ̂)  >  0

Measured on the ladder, and this is the part that makes it structural rather than bad luck:

| condition | N_rna_obs | σ²_d | (κ̂−½)² | κ̂ − ½ | disc used | **E[disc if κ = ½]** |
|---|---|---|---|---|---|---|
| `g75 ss0.50 off` | 898,960 | 3.42e-7 | 4.74e-7 | **+1.18 sd** | **5.28e-7** | ⭐ **6.62e-7** |
| `g25 ss0.50 off` | 2,697,827 | 2.85e-7 | 3.95e-9 | +0.12 sd | **0** | 5.51e-7 |
| `g98 ss0.50 off` | 72,216 | 3.85e-6 | 2.36e-6 | +0.78 sd | **0** | 7.45e-6 |
| `g75 ss0.99 off` (stranded) | 898,147 | 2.53e-5 | 0.2403 | −97.4 sd | **0.9612** | 4.90e-5 |

⭐⭐ **The leaked `disc` on `g75` is SMALLER than its own expected value under a perfectly unstranded
library.** So this is not a tail event — it is the mechanism's typical behaviour, and it fires whenever
`|κ̂ − ½|` exceeds one standard error, i.e. about **one unstranded library in three** (`P(χ²₁ > 1) =
0.317`). ⭐ And the last row confirms the mechanism is correctly inert where there is real signal.

**Why it matters even though `disc` is ~1e-7.** `disc` is multiplied by the object's fragment count:
`τ ≈ N · disc · […]`. A 12,615-fragment exon turns `disc = 2.6e-7` into `τ ≈ 2e-4`, and the retired
`τ > 1e-9` test then called that object *solvable* (`TRAPS.md` B11). Its own statement is 1,377 nats wide
against a ±10-nat grid, so it can resolve nothing — but it was scored, and it displaced the real defect.

⛔ **It cannot be repaired by subtracting more.** No finite subtraction makes a non-negative estimator of
a non-negative quantity unbiased at zero. The shape has to change. Three candidates, none yet measured:

| option | form | no tunable? | note |
|---|---|---|---|
| **(a) propagate κ̂'s uncertainty into the per-object likelihood** | integrate the strand likelihood over the posterior of κ instead of plugging κ̂ in | ✅ | the principled one; at κ ≈ ½ the likelihood in `f_g` is then genuinely flat and τ → 0 with no truncation anywhere. Largest build |
| **(b) a one-sided lower confidence bound on `(κ−½)²`** | replace the point estimate by its own lower bound | ⚠ needs a confidence level | cheap; the level is a tunable unless derived from the grid |
| **(c) shrink `disc` by the probability the library is stranded at all** | a mixture with an atom at κ = ½ | ⚠ needs a prior weight | matches the physics (a protocol either is or is not stranded) but introduces a prior |

⚠ Whatever is chosen, the gate is the same and it is available: on a genuinely unstranded condition the
strand arm's contribution to `τ` must be **0**, and on `ss0.99` it must be untouched. Both arms exist in
the panel, so this is falsifiable in one run.

---

## C5 — the reframe component-set mismatch (the junction leak), for the record

Kept because the *mechanism* is worth not rediscovering, and because the fix has a shape even though it
is not worth building alone. Full statement in `solver_derivation.md` §1–§2; the short version:

`ρ_tot` is the "how crowded is it here" number the solver divides between neighbours to measure the
capture enrichment step. At a boundary it includes the spliced fragments attached there; at a node it
structurally cannot, because a node stores only *contained* fragments and a contained fragment used no
junction. So a boundary beside an intron looks **6.6×** more crowded than the population being compared,
and a neighbouring intron's gDNA claim is scaled up by that factor (delivered claim = **3.67×** the
boundary's entire observed count).

⭐ **The rule that fixes it, if it is ever worth fixing:** the junction term belongs in `ρ_tot` only on
the steps where the junction population is *part of the claim being carried* — i.e. on the graft
(boundary → exon), where the boundary's spliced flux is deliberately added to the message, and nowhere
else. ⚠ That step-specific version is **untested**; the version measured removed it everywhere and was
worse on 2 of 4 conditions, which is itself evidence that the term is currently doing double duty.

⭐ **Why the cost is so small despite the size of the error:** the inflation multiplies the gDNA and RNA
arms by the *same* factor, and the answer is a ratio, so it very nearly cancels — a **50×** inflation
moves `f_g` by 0.12.

---

## C6 — `intergenic|exon` EDGEs as a gDNA reference (opportunity, not defect)

⭐ **Owner's point, and the data agrees.** An intergenic node is pure gDNA by definition, so an
intergenic↔exon boundary is too, and the solver pinning it to `f_g = 1` is *correct* — truth is exactly
1.0000 on all four conditions checked. Its usefulness is a **sampling** question, and it inverts with
capture:

| | objects | mass | per object |
|---|---|---|---|
| capture **off** | 2,620 | 42,459 | **16.2 fragments** — sparse, as expected: 10 M fragments over ~100 Mb, and most land in intergenic and intronic space |
| capture **on** | 1,853 | 258,937 | ⭐ **139.7 fragments** — 8.6× more, because a probe on a first or last exon partially overlaps the boundary and enriches it |

⭐⭐ **So under capture these become a large population of exactly-known-composition objects sitting at
enriched locations** — a direct measurement of what capture does to gDNA, on the stratum where nothing
else speaks. ⛔ Measure what a perfect gDNA reference there would buy before building anything.

⚠ **The narrow claim that remains from the earlier note:** using this class as the *healthy twin* of
`intron|exon` to isolate the junction flux was invalid, because a pinned class cannot be wrong and so its
correctness carries no information about the experiment (`TRAPS.md` B12). That is a statement about one
comparison, not about the class.
