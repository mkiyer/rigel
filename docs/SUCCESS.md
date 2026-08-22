# SUCCESS — what "done" means for **0.8.0**, and how it is measured

⭐⭐⭐ **"DONE" MEANS THE 0.8.0 RELEASE** (owner ruling, 2026-08-14). `pyproject.toml` says **0.7.1**; the
target is **0.8.0**, and every number in this file exists to say whether that release is reachable yet.
This is the organising frame — a measurement that does not bear on it is a diagnostic, not a target.

⛔⛔ **STAMP, APPLYING TO THE WHOLE FILE: THE PANEL HAS BEEN REBUILT THREE TIMES AND NO FIGURE HERE HAS
BEEN RE-DERIVED ON THE PANEL NOW ON DISK.** 2026-08-13 retired `pilot`/`flgap_short`/`flgap_long` and took
the ladder 36 → 16 conditions; 2026-08-19 made the simulator index-driven and put nascent RNA on every
row; **2026-08-22 replaced the uniform nascent model with a SPARSE one** (`DESIGN.md` §0b). Every number
in this file is a historical record of the substrate it was measured on. ⛔ **Treat none of them as current
behaviour** — including the 2026-08-13/14 figures, which this stamp used to certify and which predate both
later rebuilds. **Re-derive rather than trust: a number that has moved is a result, not a documentation
bug.**

---

## ⭐⭐⭐ THE SCOPE — three strata are the target, and the fourth is DEFERRED

The panel's strata are a 2 × 2 of **strandedness × capture**, with the gDNA level (`g00 / g05 / g50 /
g98`) as the ladder inside each. 0.8.0 is judged on three of the four cells.

| stratum | 0.8.0 | why |
|---|---|---|
| unstranded × capture-**OFF** | ⭐ **IN SCOPE** | the density channel alone, no capture distortion |
| stranded × capture-**OFF** | ⭐ **IN SCOPE** | strand + density, the easiest cell — a **control that must not regress** |
| stranded × capture-**ON** | ⭐ **IN SCOPE** | strand + density under a non-uniform gDNA landscape |
| unstranded × capture-**ON** | ⛔ **DEFERRED** | the blind cell. `κ = ½` makes the strand λ-term identically 0 and capture breaks the uniform-placement assumption both opportunity functions rest on |

⛔⛔ **DEFERRED IS NOT DROPPED, AND THE DISTINCTION IS THE WHOLE RULING.** Unstranded × capture-ON
**stays in every panel, every table and every report**, and it must keep being measured and quoted. What
it is not is a **development target**: no work is prioritised by its number, and no mechanism is accepted
or rejected on it, until the other three are fully optimised. ⭐ If it improves as a side effect, that is
a free win. ⛔ If a change helps the three and leaves it flat, that change is still a win — the old
ranking, in which it dominated everything because it carries most of the error, is retired.

⚠ **And it does carry most of the error, which is exactly why the ranking had to be stated.** Measured
2026-08-13/14 on the rebuilt 16-condition ladder: it carries **64.5 % of transcript-level error and 90 %
of gene-level error**. A pooled panel total is therefore a report on the deferred stratum and nothing
else (`TRAPS: never-pool-the-strata`).

⛔ **Why it is deferred rather than attacked: the tool does not have a noisy answer there, it has a FLAT
one.** On unstranded × capture-ON it emits a near-zero gDNA fraction *regardless of truth* — exon
`f_g` **0.040 / 0.0016 / 0.0021** at `g05 / g50 / g98` against truths **0.054 / 0.518 / 0.982**
(2026-08-14). ⚠ Note what that means for a single-level panel: at `g05` it looks *acceptable*, by
coincidence. Anything scored only there is scored against a constant
(`TRAPS: a-single-level-panel-cannot-see-a-constant`).

### ⛔⛔ THE LENGTH CHANNEL IS RETIRED UNTIL AFTER 0.8.0

The fragment-length / "length likelihood" channel **as a CALIBRATION composition channel** is deferred
as future work. **0.8.0 ships without it.** ⛔ Do not propose it, do not list it as a candidate, do not
rank it. Where anything reads as "the next thing to try", the correct restatement is **DEFERRED
POST-0.8.0**.

⚠ **It does not exist in `src/`** — it was A/B'd on 2026-08-10 and never shipped, so this is a scope
ruling about what to work on, not a code removal. ⚠ `length_likelihood` in `src/rigel/second_pass.py`
is a **different thing** (per-fragment second-pass assignment) and is untouched by this.

---

## ⭐⭐ THE METRIC — the CALIBRATION RESULT against ORACLE CALIBRATION

⭐⭐⭐ **THE PRIMARY NUMBER IS CALIBRATION SCORED AGAINST AN ORACLE CALIBRATION. THE END-TO-END
TRANSCRIPT NUMBER IS A THERMOMETER.** The focus of development is **calibration**; the transcript table
is downstream of it and of the EM, and a single end-to-end figure cannot say which of the two moved.

| | what is scored | against | instrument |
|---|---|---|---|
| ⭐⭐⭐ **PRIMARY** | the **prior calibration ships** — `gdna_prior_count`, `rna_prior_count`, `gdna_eff_len` per multi-locus | **O**, the same assembler fed the origin-split **truth** masses — an oracle calibration | `prior_vs_oracle.py` (`P − O`) |
| ⭐⭐ **PRIMARY, per object** | each region's and each boundary's own `f_g`, and whether it is **confidently** wrong | the oracle payload — the production accumulator run on the BAM split by true origin | `solvability_audit.py` |
| ⭐ **PRIMARY, one number** | the library `f_gdna` | the simulator's per-fragment truth | `calibration_truth_ab.py` |
| ⭐⭐ **CONTROLS** | zero-gDNA and zero-RNA, where truth is a CONSTANT | 0.000 and 1.000 exactly | `zero_controls.py`, and the `g00` rung |
| **THERMOMETER** | the transcript table a user reads | `truth_abundances.tsv` | `quant_accuracy.py --arm base` |

⭐ **Why `P − O` and not the transcript number: attribution.** `O` is *calibration done perfectly with
the shipped assembler*, so `P − O` is calibration's own error and nothing else — the transcript number
adds the assembler, the effective-length model, the EM's ambiguity and the annotation on top of it. The
same instrument reports `O − Fo` (the ASSEMBLER's error, against the EM's own candidate count) beside
it, because they are different repairs in different files.

⛔ **Read the thermometer, do not steer by it.** The transcript figure moves for reasons that have
nothing to do with calibration — 95.5 % of the stranded × capture-OFF misassignment is ordinary isoform
ambiguity, which is not a Rigel-specific defect, and a perfect prior is **neutral to worse** on the
capture-OFF strata. A calibration change that improves `P − O` and leaves the transcript table flat has
still done its job.

⚠ **The noise floor is measured, not assumed: 0.996–1.013** on the rebuilt 16-condition ladder
(2026-08-13/14). Any arm ratio inside that band is a sampling draw. ⛔ `quant_accuracy.py --arm
base_reseed` prints it beside the effect and must be run in the same session
(`TRAPS: re-record-the-baseline`).

⛔ **Report per stratum, always.** Four cells, never one total. The pooled figure has been 97 % the
deferred stratum, and quoting it has twice hidden a sign flip between the others.

### ⭐⭐ WHY THE LADDER GIVES gDNA AND RNA **EQUAL** FRAGMENT LENGTHS

⛔ **The reason is stronger than "the length channel is neutralised so residual error is attributable to
density and strand", which is the weaker form written down across this repo and being corrected
everywhere (owner, 2026-08-14).** The real reason:

⭐⭐ **THE EM ALREADY USES THE FRAGMENT-LENGTH DISTRIBUTION.** Give gDNA and RNA a large length
difference and the EM can assign fragments **on length alone** — it will separate the two populations
without calibration having contributed anything, **bypassing the calibration phase entirely and masking
every bug in it**. Equal lengths remove that shortcut and **force calibration to be exercised**, which
is the only condition under which the primary metric above means what it says.

⚠ So calibration is left with exactly the channels it is supposed to use: **strand and density**, plus
belief propagation across objects (currently OFF — see `CLAUDE.md`). ⚠ "Equal" is a configuration, not a
guarantee: transcript-length truncation leaves a realised **+3.6 to +4.7 bp** gDNA-over-RNA gap, measured
and priced at **under 2.5 % of the per-object error** before it was accepted. Score the length axis
against `truth_fragment_lengths.tsv`, never against a nominal parameter.

### ⛔⛔ THE RULER ITSELF — every ceiling measured before 2026-08-14 had a WRONG one installed

⛔⛔ **`effective_lengths_em` is built BEFORE `assemble_priors`, and every measurement arm in this repo
patches `assemble_priors`. So the effective-length shrinkage has NEVER been priced by any ceiling**
(measured 2026-08-14: construction at `pipeline.py:816`, the prior assembled at `pipeline.py:839`).
Every ceiling number this project has quoted was measured with a wrong ruler still in place.

⛔ **And the ruler is wrong in a way a zero control makes unmistakable.** At the ZERO-gDNA control the
shipped shrinkage contracts every transcript by a mean factor of **0.345** when the correct answer is
**exactly 1.000**, on **13,673 of 15,669 transcripts** — *including on capture-OFF, where the module's
own contract says the factor must be 1*. `rho_ref` is fabricated entirely from false-positive gDNA.

⭐ **It is a SYMPTOM, not an independent bug, and that was proven by substitution.** Substituting **only
the composition arrays** and re-running the **shipped** shrinkage function gives the correct factor:
`g00` capture-OFF **0.345 → 1.000** (truth), `g50` unstranded × capture-ON **0.834 → 0.401** (truth).
One split, two consumers, one function — `priors.py` imports `_global_reference_density` from
`capture_eff_length.py`. ⛔ So it must not be repaired separately: fix the composition and check that the
factor follows. A second, independent shrinkage repair would be half of a cancelling pair
(`TRAPS: a-cancelling-defect-pair`).

⚠ **A second lane is built and not wired.** The per-transcript RNA prior (`rna_prior_weight`) exists end
to end but is **never passed in production** (`pipeline.py:841` omits it), so the EM's default rule
carries zero per-transcript information. ⛔ Any statement of the form "the prior is worth X" that
predates 2026-08-14 was measured with that lane dark.

⭐ **What the two prior ceilings are worth, measured 2026-08-14** — and they answer different questions:

| ceiling | the three in-scope strata | the deferred stratum |
|---|---|---|
| a perfect **per-LOCUS** prior | not measured separately — it is the deferred cell this arm was aimed at | gene-level error **0.099** |
| a perfect **per-TRANSCRIPT** prior | gene-level **0.37–0.58**, and false-positive mass cut **100–200×** | **0.641** — it does NOT rescue it |

⛔ Read the second row in both directions: the per-transcript lane is a large win **on the strata 0.8.0
is scoped to**, and it is not the answer to the deferred one.

---

## STAGE A — the accumulator: ✅ **DONE, and that is a measurement**

⚠ **Kept as a RECORD and as a regression check, not as work.** Stage A asks whether the information is
*there*; the 0.8.0 work asks whether calibration *finds* it. Splitting them is what makes an attribution
possible at all, and this project has already spent sessions on the wrong one of the two.

⭐ **The criterion was "handing calibration the exact fragment-length distribution changes nothing", and
`calibration_truth_ab.py --ceiling` measures exactly that: perfecting BOTH length models is worth 2.6 % of
the deliverable** (down from 22.2 % before the four-pool model). Anchor −0.00 % ± 0.01; RNA/gDNA models
−0.21…+1.12 % and −0.01 % off capture; second pass 90.6 % exact per fragment.

⛔ **Do NOT chase the gDNA model's +7.9 % under capture.** Its cause is known and is not a divisor — both
opportunity functions assume gDNA is placed uniformly and capture does not — for ~1.7 % of a number that is
97 % dominated by something else.

**The five identities that keep it done**, all still gated and green:

| | | |
|---|---|---|
| C++ vs the executable specification | **byte-identical** | `tests/native/test_accumulator_spec.py` |
| the same BAM at 1/2/4/8 workers | **bit-identical**, every bank | `test_accumulator_worker_determinism.py` |
| `Σ node_start_count == deposited` | exact | same |
| `deposited + deferred + dropped_* == offered` | exact, deferred **non-empty** | same |
| the origin partitions sum to the full payload | exact, every channel | `tests/calibration/_oracle.py` |

✅ **The R1-sense coverage gap is CLOSED** (2026-08-11). `ReadSimConfig.r1_sense` selects the protocol
DIRECTION — R1-antisense (TruSeq dUTP, the default) or R1-sense (KAPA) — independently of
`strand_specificity`, which is the FIDELITY about it. Gated in `test_strand_sense_convention.py`: both
directions, the per-fragment mirror (2,000/2,000 fragments have R1 on the opposite strand), and that the
DECONVOLUTION recovers the same gDNA fraction either way. ⚠ The panel ORCHESTRATOR does not expose the
knob, so no ladder condition is R1-sense — the engine can, the panel does not.

**The three Stage-A criteria, by name** (they are criteria, not rules, and they are what the regression
block at the end of this file re-runs): **FIDELITY** — the tally reproduces the specification exactly;
**BIAS** — no channel is systematically off against per-fragment truth; **SUFFICIENCY** — the stored
tally carries enough length information that perfecting it changes nothing downstream.


## STAGE B — calibration: ⭐⭐⭐ **THIS IS THE 0.8.0 WORK**

⭐ Calibration is iterative; **pass-0** is the **prior-free** solve, the first thing that happens after
initialisation, before any fitted prior exists. It is the right place to start because it has no feedback
loop in it — every later iteration's behaviour is conditional on pass-0 being sane, and
`TRAPS: variance-fitted-on-the-belief` is what happens when a prior is fitted on an unsound belief.

### The oracle, and what it buys

`tests/calibration/_oracle.py` already produces the object we need: **the production accumulator run on
the BAM split by true origin**, so for every region and every boundary we have the *true* gDNA and RNA counts on
each strand. Three quantities follow, and the differences between them are the whole diagnostic:

| | quantity | what a gap to the next one means |
|---|---|---|
| **T** | the truth: each object's real gDNA/RNA split | — |
| **C** | a **ceiling**: the best answer reachable under stated conditions | `T − C` is **information the accumulator destroyed** → Stage A work |
| **P** | what pass-0 actually produces | `C − P` is **solver gap** → the 0.8.0 work |

⭐⭐ **That decomposition is a measurement, not a build.** It is the per-object form of the ceiling
discipline that has twice re-ranked this project (`TRAPS: measure-the-ceiling-first`). ⛔ And it is
subject to the ruler caveat above: a ceiling that patches `assemble_priors` does not see the
effective-length shrinkage at all.

### ⛔⛔ HOW PASS-0 IS SCORED — and the mistake that has to be avoided

**Pass-0's job is not accuracy. It is to produce a SUBSTRATE the gDNA hyperprior can be fitted
against.** An object with no own evidence reporting `f_g ≈ ½` at **zero precision** is *correct* — it
is stating a true fact about itself. Scoring every object that carries mass therefore counts honest
ignorance as error: it reported **0.3150** where the answer on the objects pass-0 can actually solve
was **0.0456**, and 99.5 % of the difference was the undetermined class.

⭐ **The measurement is a partition, in this order** (`scripts/design/solvability_audit.py`):

1. **UNDETERMINED** — no own evidence. ⛔ *Excluded from the error denominator.* Its only failure mode
   is the opposite one: claiming a precision it has not earned. ⛔⛔ **AND THAT CHECK MUST EXIST, OR THE
   EXCLUSION HIDES THE LARGEST ERROR IN THE LIBRARY** — it did not exist until 2026-08-04, and the cost
   was **1,056,019 fragments** of unreported error on a condition publishing `mwae 0.0170`
   (TRAPS: excluding-a-population-hides-it). `undetermined_overreach_rows` buckets the class by `|f_pred − ½|` and reports its
   error and its precision claim; the correct answer for the class is ½ at `sd = ∞`.
2. **SOLVABLE and right.**
3. **SOLVABLE and wrong**, split by confidence. ⭐⭐ **Confidently wrong is the defect** — a wrong value
   with a tight variance outvotes correct neighbours *and anchors the prior*, so it propagates.

⚠ The confidence comparison is in LOG space (`var_gdna` is `Var(log f_g)` despite its name,
TRAPS: log-variance-is-not-linear), and the headline is a calibration curve, which needs no threshold.

⛔⛔ **AND THE PARTITION ALONE IS NOT ENOUGH — "no own evidence" IS NOT A BINARY.** Step 1 was
implemented as `tau_lam > 1e-9`, and the strand arm's information `I(f_g) ∝ (2κ−1)²` is exactly zero
only at κ = **½** while κ is *fitted*. At 10 M fragments a genuinely unstranded library fits
κ̂ = 0.500689, so τ lands at ~5e-7 and the threshold promotes the object to "solvable" while its own statement
has **sd(λ) = 1,377 nats against a solver that represents λ only on ±10**. Measured:
**79.1 % of the "solvable" error sat on objects with no usable answer of their own** (TRAPS: a-threshold-on-a-fitted-residue).

⭐ **So strength is reported as a CURVE over `sd(λ) = 1/√τ` decades**, beside the partition, and the
panel table carries a **`weak%`** column — the share of the scored error above 10 nats. ⛔ **A better
threshold is not available and was refuted:** τ is *continuous* across the region on 4 of 5 ladder
conditions, so any floor would be a tuned constant. ⭐⭐ **Read `weak%` before `mwae`.** A row with
`weak%` near 100 is reporting the relay and the reference, not a solve.

⚠ And `locked` is the **G1 / structurally pure-gDNA** class on *both* axes (`region_geometry.g1_locked`),
never `~solvable & is_region` — a structurally-locked *boundary* is certain, not ignorant. ⛔ It is
deliberately **not** the same mask as `region_init`'s region-only `struct_lock`, which governs message
*emission* (TRAPS: two-masks-one-name).

### The instrument — `scripts/design/pass0_vs_oracle.py`

✅ **Built.** It scans, builds T, runs four arms, and scores every one of them per object and per class.
Gates: `tests/calibration/test_pass0_vs_oracle.py`, **eleven**, each carrying its own perturbation
(re-derived 2026-08-17 with `pytest --collect-only -q`; this line and the script's own docstring both
still said *nine*, and the script's copy is the one that has not been corrected yet). ⭐ It is also what
POPULATES the oracle cache every other scorer reads, which is why `panel.py cache` runs it.

⛔ **C is TWO ceilings, each defined by a lever that already exists, and never an estimator.** "The best
split any estimator could produce" is not something you can write down; trying to is the magic-number
failure mode with an estimator in place of a constant (TRAPS: no-magic-numbers).

| | what it is | what its gap to P means |
|---|---|---|
| **C_input** | `calibrate` handed the simulator's own post-capture length pmfs — the same override `calibration_truth_ab.py --ceiling` uses — at *both* solve depths | how much of the error is wrong **inputs** rather than wrong **solving** |
| **C_info** | a **classification**, per object: is the 2×2 of `EQUATIONS.md` §3.1 identified from this object's own stored channels at all? | ⛔ **not a gap.** C_info ignores neighbours and the sweep does not, so it can be "worse" than P. The useful reading is the reverse — see below |

⚠ **C_input is a *length*-input ceiling, and under the 0.8.0 scope it is a DIAGNOSTIC rather than a
route.** The other library-level inputs are injectable (`InjectedCalibrationPriors`) but the simulator
writes no truth for them, so injecting one would be an A/B against a guess. κ is the one with a truth
value and it is not free either (TRAPS: specificity-and-sense-are-complements). ⛔ Do not read C_input as
an argument for the length channel — that is deferred post-0.8.0 by ruling, not by measurement.

⭐⭐ **The cross-tab is the point.** Objects **undetermined by C_info** *and* carried **entirely by the
relay** have no answer of their own at all — whatever pass-0 reports there came from neighbours and the
population prior. That cell is reported with its mass share and its error share, and it needs no
confidence threshold because it is a cell of a partition rather than a cutoff.

The two classifications, both mutually exclusive and exhaustive (gated: the mass *and* the error
decompose over each exactly):

* **the solver's own** — `own_evidence` / `relay_only` / `struct_lock`, reproducing `region_init`'s
  definitions and cross-checked against `composition_evidence_census.py`'s;
* **C_info's** — `identified` / `undet_no_separation` / `undet_out_of_range` / `absent`.

⛔ **Every arm, T included, is UNDRAINED, and that is forced.** The second pass's multinomial is scored
against the payload's own densities, so draining three origin partitions separately is not the same
operation as draining the whole and the sum-to-full identity would not survive it. The drain is a
different axis, measured by `second_pass_accuracy.py` and `calibration_truth_ab.py`.

---

## ⭐⭐⭐ 0.8.0 IS DONE WHEN

⚠ **Each row names the number that is the bar. The numeric threshold on each is an OWNER call and is not
invented here** — what this file fixes is *which* quantity is being judged, on *which* strata, against
*what* truth.

1. ⭐⭐ **`P − O` is small on all three in-scope strata**, and the residual that remains is *attributed* —
   to the assembler (`O − Fo`), to the composition, or to a class that is provably undetermined. ⛔ It is
   not done while the residual sits on objects `C_info` calls **identified**: measured, **36.2 % of the
   relay's error** does (`ROADMAP.md`).
2. ⭐⭐ **The zero controls read zero.** `zero_controls.py` on both arms, and the `g00` rung of the
   ladder. ⛔ Today the shipped prior claims **2,067,637 gDNA fragments in libraries containing none**,
   and **1,707,321 of them are at unstranded × capture-OFF — an IN-SCOPE stratum that reads healthy on
   every contaminated row**. Only a zero control finds that.
3. ⭐⭐ **The effective-length shrinkage is correct BECAUSE the composition is**, not because it was
   patched: at `g00` the factor reads **1.000**, and elsewhere it tracks the truth factor when only the
   composition arrays are substituted. ⛔ A separate shrinkage correction is a defect, not a fix.
4. ⭐ **Pass-0 is monotone in the obvious sense**: adding real evidence to an object never moves its
   answer away from truth.
5. ⭐ **Pass-0 depends on no quantity a later iteration produces** — no feedback in the first solve.
6. **The three in-scope strata do not regress on the thermometer**, and the deferred stratum is
   **reported** on every table.

⛔ **The gate that held this work back, kept because it will apply again.** A prior-free solve tuned
against a 15 %-wrong input on exactly the conditions it is meant to rescue is tuned against a
manufactured discriminant, and every conclusion drawn from it has to be thrown away when the input is
fixed. That has already happened here once, on a substrate defect (`TRAPS: prove-the-substrate`). ⚠ The
live form of it is the ruler: **do not price anything against a ceiling that patches `assemble_priors`
while the effective-length shrinkage is still built upstream of it.**

---

## Reporting numbers vs target numbers

⭐ **The target is the calibration result against the oracle calibration** — the table above, per
stratum, on the three in-scope cells.

**The library gDNA fraction and the transcript table are the products**, and both are reported every
session so their trajectory is visible. ⛔ Neither is the steering wheel in this phase: the library
figure mixes accumulator, solver and unidentifiability into one number (78 % of its current error is in
calibration rather than in the accumulator), and the transcript figure adds the EM and the annotation on
top of that. **Read them as thermometers.**

⚠ **Score the contaminated conditions.** Zero-gDNA rows are saturated at truth = 0 exactly, so anything
that lowers the estimate "improves" them (TRAPS: zero-target-guards-are-one-sided). They are
false-positive checks, nothing more — which is precisely why they are a required CONTROL on every
experiment and never a target.

⛔ **Quote the SHIPPED column, not pass-0**, and never a pooled total. A −37.2 % pass-0 win was −3.9 %
shipped (TRAPS: the-intermediate-is-not-the-deliverable).

---

## The instruments, in the order you would run them

⭐⭐ **Steps 0–1 build the scenarios ONCE. After that a calibration change is re-scored in minutes**,
because both caches survive it — the scan cache makes calibration re-runnable without rescanning, and the
oracle cache is the origin-split truth every scorer reads. ⛔ Neither depends on `calibration/`, so a
calibration edit invalidates nothing here.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite
LADDER=$SUITE/ladder
INDEX=$SUITE/rigel_index
CFG=scripts/sim/configs/gdna_ladder.yaml

# 0. WHERE AM I?  Every stage is expensive and resumable, and this names the next one.
python scripts/sim/panel.py status --config $CFG

# 1. BUILD THE SCENARIOS — once.  `cache` builds BOTH caches and is the step that makes the rest cheap.
#    ⚠ the reference carve is manual: it needs the SOURCE genome/GTF, which a panel config does not name.
python scripts/sim/panel.py build    --config $CFG
python scripts/sim/panel.py simulate --config $CFG --jobs 8
python scripts/sim/panel.py cache    --config $CFG --jobs 8

# 2. IS THE SUBSTRATE SOUND?  (TRAPS: prove-the-substrate — prove the simulator before the code)
python scripts/design/simulator_gates.py --suite $LADDER --reference $SUITE/reference
python scripts/design/suite_resolves.py $INDEX --suite $LADDER

# 3. ⭐⭐⭐ THE PRIMARY METRIC — CALIBRATION AGAINST ORACLE CALIBRATION.
#    (a) the PRIOR the EM actually reads, P vs O, per stratum   ~50 s/condition with the cache warm
python scripts/design/prior_vs_oracle.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache --jobs 6
#    (b) per OBJECT: solvable, solved wrong, and CONFIDENTLY wrong.  Read `weak%` before `mwae`.
python scripts/design/solvability_audit.py --suite $LADDER --index $INDEX \
       --oracle-cache $LADDER/oracle_cache
#    (c) the one-number summary: the library f_gdna against truth.  `--scan-cache` is the scan-cache ROOT;
#        `_main` is the undrained full payload the oracle cache already holds.
#        ⚠ a panel cached only by `pass0_vs_oracle.py` has NO `_main` on the four `g00` rows, because it
#        holds every zero-gDNA condition out — so read which conditions this actually scored rather than
#        assuming 16 (`TESTING.md` §0: the oracle cache's count is not stable in either direction).
python scripts/design/calibration_truth_ab.py --scan-cache $LADDER/oracle_cache --cache-subdir _main \
       --suite $LADDER --index $INDEX

# 4. THE TWO ZERO CONTROLS — owner-required on EVERY experiment, both arms.
python scripts/design/zero_controls.py

# 5. THE THERMOMETER — the tool end to end, and the prior-injection ceiling above it.
#    ⛔ --jobs 2, not more: run_pipeline holds 7-8.5 GB per 10 M-fragment condition.
#    ⛔ `base_reseed` is the noise floor (0.996-1.013) and any arm delta inside it is a sampling draw.
python scripts/sim/panel.py score  --config $CFG --arms base base_reseed oracle --jobs 2
python scripts/sim/panel.py report --config $CFG --arms base base_reseed oracle

# 6. STAGE A is CLOSED — this block is a REGRESSION check, run it after an accumulator or native change.
python -m pytest tests/native tests/calibration -q     # FIDELITY
python scripts/design/fl_anchor_gap.py --drain         # BIAS: anchor + both length models vs truth
python scripts/design/gdna_pool_census.py              #       the four gDNA pools, each vs its opportunity
python scripts/design/second_pass_accuracy.py          #       per-fragment length assignment
python scripts/design/calibration_truth_ab.py --scan-cache $LADDER/oracle_cache --cache-subdir _main \
       --suite $LADDER --index $INDEX --ceiling        # SUFFICIENCY: what a perfect length model is worth
```

⚠ Steps 0 and 2–4 take about 15 minutes on a built panel. ⚠ Step 2 was re-run on the rebuilt
16-condition ladder on 2026-08-13: `simulator_gates` **6/6**, `suite_resolves` **11/12** (only `(c)`
replicate-pairs fails, the one deferred by owner ruling 2026-07-30). **Run them as a set and record them
together** — TRAPS: re-record-the-baseline.

⛔ **When dissecting rather than scoring**, the loop is: run the panel → take the worst **in-scope**
scenario → dissect it to the single highest-error object (`worst_objects.py`) → find the cause → fix →
repeat. `CLAUDE.md`'s script table lists the dissection instruments; ⛔ the worst scenario overall is the
deferred stratum, and picking it is how this ranking gets quietly re-inverted.
