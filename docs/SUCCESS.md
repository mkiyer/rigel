# SUCCESS — what "working" means in this phase, and how it is measured

⛔ **Calibration is not the yardstick right now.** Running `calibrate` end to end and reading the library
gDNA fraction mixes an accumulator defect, a solver approximation and an unidentifiability into one
number, and that number then cannot tell you which of the three moved. Until the accumulator is provably
accurate, the end-to-end figure is a **report**, not a target.

So performance is measured at **two stages**, and the second is gated behind the first.

| | stage | question | truth it is scored against |
|---|---|---|---|
| **A** | the **accumulator** | Is the tally faithful, unbiased, and sufficient? | the simulator's per-fragment truth (oracle BAM read names) |
| **B** | calibration **initialisation + pass-0** | Does the prior-free solve find what the payload actually contains? | an **oracle** payload: the same accumulator run on the BAM split by true origin |

⭐ **The point of splitting them is attribution.** Stage A asks whether the information is *there*; Stage B
asks whether the solver *finds* it. A single end-to-end number cannot separate those, and this project has
already spent sessions on the wrong one of the two.

---

## STAGE A — the accumulator: ✅ **DONE, and that is a measurement**

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


## STAGE B — calibration initialisation and pass-0

⭐ **OPEN, and step 1 has landed.** Calibration is iterative; pass-0 is the **prior-free** solve, the
first thing that happens after initialisation, before any fitted prior exists. It is the right place to
start because it has no feedback loop in it — every later iteration's behaviour is conditional on pass-0
being sane, and TRAPS: variance-fitted-on-the-belief is what happens when a prior is fitted on an unsound belief.

### The oracle, and what it buys

`tests/calibration/_oracle.py` already produces the object we need: **the production accumulator run on
the BAM split by true origin**, so for every region and every boundary we have the *true* gDNA and RNA counts on
each strand. Three quantities follow, and the differences between them are the whole diagnostic:

| | quantity | what a gap to the next one means |
|---|---|---|
| **T** | the truth: each object's real gDNA/RNA split | — |
| **C** | a **ceiling**: the best answer reachable under stated conditions | `T − C` is **information the accumulator destroyed** → Stage A work |
| **P** | what pass-0 actually produces | `C − P` is **solver gap** → Stage B work |

⭐⭐ **That decomposition is the deliverable of Stage B's first step, and it is a measurement, not a
build.** It is the per-object form of the ceiling discipline that has twice re-ranked this project
(TRAPS: measure-the-ceiling-first).

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
κ̂ = 0.500689, so τ lands at ~5e-7 and the region bound promotes the object to "solvable" while its own statement
has **sd(λ) = 1,377 nats against a solver that represents λ only on ±10**. Measured:
**79.1 % of the "solvable" error sat on objects with no usable answer of their own** (TRAPS: a-threshold-on-a-fitted-residue).

⭐ **So strength is reported as a CURVE over `sd(λ) = 1/√τ` decades**, beside the partition, and the
panel table carries a **`weak%`** column — the share of the scored error above 10 nats. ⛔ **A better
threshold is not available and was refuted:** τ is *continuous* across the region on 4 of 5 ladder
conditions, so any floor would be a tuned constant. ⭐⭐ **Read `weak%` before `mwae`.** A row with
`weak%` near 100 is reporting the relay and the reference, not a solve.

⚠ And `locked` is the **TRAPS: no-magic-numbers** class on *both* axes (`node_geometry.g1_locked`), never
`~solvable & is_node` — a structurally-locked *boundary* is certain, not ignorant. ⛔ It is deliberately
**not** the same mask as `node_init`'s node-only `struct_lock`, which governs message *emission*
(TRAPS: two-masks-one-name).

### The instrument — `scripts/design/pass0_vs_oracle.py`

✅ **Built.** It scans, builds T, runs four arms, and scores every one of them per object and per class.
Gates: `tests/calibration/test_pass0_vs_oracle.py`, nine, each carrying its own perturbation.

⛔ **C is TWO ceilings, each defined by a lever that already exists, and never an estimator.** "The best
split any estimator could produce" is not something you can write down; trying to is the magic-number
failure mode with an estimator in place of a constant (TRAPS: no-magic-numbers).

| | what it is | what its gap to P means |
|---|---|---|
| **C_input** | `calibrate` handed the simulator's own post-capture length pmfs — the same override `calibration_truth_ab.py --ceiling` uses — at *both* solve depths | how much of the error is wrong **inputs** rather than wrong **solving** |
| **C_info** | a **classification**, per object: is the 2×2 of `EQUATIONS.md` §3.1 identified from this object's own stored channels at all? | ⛔ **not a gap.** C_info ignores neighbours and the sweep does not, so it can be "worse" than P. The useful reading is the reverse — see below |

⚠ **C_input is a *length*-input ceiling.** The other library-level inputs are injectable
(`InjectedCalibrationPriors`) but the simulator writes no truth for them, so injecting one would be an
A/B against a guess. κ is the one with a truth value and it is not free either (TRAPS: specificity-and-sense-are-complements); it is the
next lever and a separate measurement.

⭐⭐ **The cross-tab is the point.** Objects **undetermined by C_info** *and* carried **entirely by the
relay** have no answer of their own at all — whatever pass-0 reports there came from neighbours and the
population prior. That cell is reported with its mass share and its error share, and it needs no
confidence threshold because it is a cell of a partition rather than a region bound.

The two classifications, both mutually exclusive and exhaustive (gated: the mass *and* the error
decompose over each exactly):

* **the solver's own** — `own_evidence` / `relay_only` / `struct_lock`, reproducing `node_init`'s
  definitions and cross-checked against `composition_evidence_census.py`'s;
* **C_info's** — `identified` / `undet_no_separation` / `undet_out_of_range` / `absent`.

⛔ **Every arm, T included, is UNDRAINED, and that is forced.** The second pass's multinomial is scored
against the payload's own densities, so draining three origin partitions separately is not the same
operation as draining the whole and the sum-to-full identity would not survive it. The drain is a
different axis, measured by `second_pass_accuracy.py` and `calibration_truth_ab.py`.

### ⭐ Stage B is DONE when

* `C − P` is small **and** its remainder is concentrated in classes that are provably undetermined rather
  than in classes that had evidence — ⚠ **measured, and it is not**: 36.2 % of the relay's error sits on
  objects C_info calls *identified*, i.e. on objects that had the evidence (`ROADMAP.md` §2);
* pass-0 is monotone in the obvious sense: adding real evidence to an object never moves its answer away
  from truth;
* pass-0 does not depend on any quantity that later iterations produce (no feedback in the first solve).

⛔ **The gate that held Stage B back, kept because it will apply again.** A prior-free solve fed a
15 %-wrong length model on exactly the conditions it is meant to rescue would be tuned against a
manufactured discriminant, and every conclusion drawn from it would have to be thrown away when the model
was fixed. That has already happened once here, on a substrate defect (TRAPS: prove-the-substrate). ⚠ **It is live
again for the per-region length likelihood**: the fitted gap's sign is currently wrong, so wiring that
channel now would repeat exactly this — see `ROADMAP.md` §2 **the-cancelling-pair**.

---

## Reporting numbers vs target numbers

**The library gDNA fraction is the product**, and it is reported every session so its trajectory is
visible. But it is **not** the target during this phase: it mixes all three error sources, and 78 % of its
current error is in calibration rather than in the accumulator. Read it as a thermometer, not a
steering wheel.

⚠ **Score the contaminated conditions.** Zero-gDNA rows are saturated at truth = 0 exactly, so anything
that lowers the estimate "improves" them (TRAPS: zero-target-guards-are-one-sided). They are false-positive checks, nothing more.

---

## The instruments, in the order you would run them

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite

# 0. is the SUBSTRATE sound?  (TRAPS prove-the-substrate — prove the simulator before the code)
python scripts/design/simulator_gates.py --suite $SUITE/pilot --reference $SUITE/reference
python scripts/design/suite_resolves.py $SUITE/rigel_index --suite $SUITE/pilot

# A1. fidelity
python -m pytest tests/native tests/calibration -q

# A2. bias
python scripts/design/fl_anchor_gap.py --drain          # anchor + both length models vs truth
python scripts/design/gdna_pool_census.py               # the four gDNA pools, each vs its own opportunity
python scripts/design/second_pass_accuracy.py           # per-fragment length assignment

# A2's "done" condition, and the scoping number for the whole phase
python scripts/design/calibration_truth_ab.py --ceiling

# B. pass-0 against the oracle, per object and per class (~20 min PER CONDITION: it re-scans the
#    BAM once per origin partition, and there is no cache for that)
python scripts/design/pass0_vs_oracle.py

# C. CALIBRATION'S ENDPOINT — LocusPriors, what the EM actually reads, against the oracle.
#    ~50 s per condition with the oracle cache warm; --jobs 6 does the ladder in ~6 min.
python scripts/design/prior_vs_oracle.py --suite $SUITE/ladder --jobs 6

# D. THE TOOL, END TO END, against the simulator's per-transcript truth — and the ceiling above it.
#    ⛔ --jobs 2, not more: run_pipeline holds 7-8.5 GB per 10 M-fragment condition.
export RIGEL_ARMS=~/Downloads/rigel_runs/arms
for arm in base noop oracle base_reseed; do
  python scripts/design/quant_accuracy.py --arm $arm --out $RIGEL_ARMS/qa_$arm.jsonl --jobs 2
done
python scripts/design/arm_identity.py $RIGEL_ARMS/qa_base.jsonl $RIGEL_ARMS/qa_noop.jsonl   # MUST pass
python scripts/design/quant_accuracy.py --report $RIGEL_ARMS/qa_{base,base_reseed,oracle}.jsonl
```

⚠ Steps 0–B take about 15 minutes on the pilot. **Run them as a set and record them together** —
TRAPS: re-record-the-baseline.

⭐⭐ **C AND D ARE WHAT "MEASURE THE WHOLE TOOL" MEANS, AND BOTH ARE NEW (2026-08-07).** Everything above
them scores an INTERMEDIATE — per-object `f_g`, or the library figure. `prior_vs_oracle.py` scores the
object the EM is actually handed, and `quant_accuracy.py` scores the transcript table a user reads.
⛔ D's `noop` arm must be byte-identical to `base` and the seed must be pinned to make that possible
(`TRAPS: the-deliverable-is-not-reproducible-by-default`); `base_reseed` is the noise floor and any arm
delta smaller than it is a sampling draw, not an effect.
