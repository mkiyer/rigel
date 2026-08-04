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

## STAGE A — the accumulator

The accumulator is a **lossy summary** of a BAM, by design: it replaces hundreds of millions of fragments
with a few integer channels per object. Accuracy therefore means three different things, and all three
have to be asked separately.

### A1. FIDELITY — does the tally contain what it should?

*Is every fragment deposited on exactly the objects the deposit rule names?*

| gate | target | instrument | status |
|---|---|---|---|
| C++ matches the executable specification | **byte-identical** | `tests/native/test_accumulator_spec.py`, `test_accumulator_native_parity.py` | ✅ |
| the same BAM at 1/2/4/8 workers | **bit-identical**, every bank | `tests/native/test_accumulator_worker_determinism.py` | ✅ |
| `Σ node_start_count == deposited` | exact | `test_accumulator_worker_determinism.py` | ✅ |
| `deposited + deferred + dropped_* == offered` | exact, and the deferred term must be **non-empty** | same | ✅ |
| the origin partitions sum to the full payload | exact, every channel | `tests/calibration/_oracle.py` | ✅ |

⭐ **All absolute — no thresholds anywhere.** ⚠ And the fourth row's second half is the load-bearing part:
a conservation identity over an empty term is satisfied by any bookkeeping (`TRAPS.md` A1).

⭐ **The sum-to-full identity is why the oracle is trustworthy.** The oracle is not a reimplementation: it
runs the *production* accumulator on the BAM split by true origin. Because the deposit is per-fragment
independent, the partitions summing to the full payload **proves** the split is the production payload
partitioned by origin.

### A2. BIAS — are the derived library-level models unbiased?

The accumulator's derived products are the fragment-length models and the per-object densities. These are
where the current defect lives.

| quantity | scored against | instrument |
|---|---|---|
| the **anchor** (`deposited_lengths`) | `truth_fragment_lengths.tsv`, on a zero-gDNA condition where anchor and pool describe ONE population | `fl_anchor_gap.py --drain` |
| the **RNA** length model | the same | `fl_anchor_gap.py --drain` |
| the **gDNA** length model | the same | `fl_anchor_gap.py --drain`, `gdna_pool_census.py` |
| per-fragment length assignment (second pass) | per-fragment truth in the read names | `second_pass_accuracy.py` |

⭐⭐ **THE "DONE" CONDITION IS DERIVED, NOT CHOSEN: a length model is accurate enough when handing
calibration the EXACT distribution changes nothing.** That is measurable directly —
`calibration_truth_ab.py --ceiling` runs arms in which the simulator's own pmf replaces the fitted one — so
the target self-calibrates against the consumer's actual sensitivity instead of against a number somebody
picked. ⚠ It also cuts both ways: it is how the RNA length model was shown to be **finished** while it was
still being worked on.

⚠ **The identifying quantity is the GAP `μ_g − μ_r`, not either mean** (`EQUATIONS.md` §3.1). A model can
be 15 % off and cost little if the other is off the same way; two models each 1 % off in opposite
directions can cost more. **Report the gap alongside the two errors, always.**

### A3. SUFFICIENCY — what does the tally NOT retain?

*Given the payload and nothing else, is the answer recoverable in principle?*

This is the question a bias measurement cannot ask, and it is answered by the oracle: if the oracle's
per-object gDNA/RNA split cannot be reconstructed from the payload's channels, **no solver will ever
recover it** and the fix belongs in the accumulator, not in calibration.

Known and *deliberate* losses, each with the reason:

| lost | why it is acceptable |
|---|---|
| which fragment went where | the summary is the point; per-fragment truth is only for scoring |
| density below one fragment length | not resolvable by **any** design (`TRAPS.md` D8). Composition still is |
| the mature/nascent split | "RNA is RNA" in the accumulator — an owner ruling, not an oversight |
| gDNA/RNA per object | that *is* calibration's job; the accumulator stores the channels it deconvolves from |

⛔ **Known and NOT acceptable — this is the live Stage-A work:**

| lost | consequence |
|---|---|
| **multimappers deposit nothing at all** | the tally is unique-mappers-only, and multimappers are 16.7–59.3 % of intergenic and ~53.6 % of intronic fragments — biased toward repeats, which is exactly the gDNA pool |
| single-reference *cis* chimeras deposit on the intergenic path | contaminates the purest gDNA pool |
| `detect_chimera` is blind to two real populations | same |

### ⭐ Stage A is DONE when

1. every A1 row is exact (**already true**);
2. the ceiling arms in `calibration_truth_ab.py --ceiling` show that **perfecting either length model
   changes the deliverable by less than the run-to-run spread of the estimator itself** — i.e. the
   accumulator is no longer the limiting factor;
3. the A3 "not acceptable" list is empty, or each entry has a *measured* cost showing it does not matter.

⭐ **Condition 2 is now TRUE, which is what opened Stage B.** Perfecting both length models moves the
library deliverable by **2.6 %** (`ROADMAP.md` §1) and the per-object answer by **≤ 2.5 %** — negative on
half the arms (`ROADMAP.md` §2). Condition 3 is still open, and is ranked there.

---

## STAGE B — calibration initialisation and pass-0

⭐ **OPEN, and step 1 has landed.** Calibration is iterative; pass-0 is the **prior-free** solve, the
first thing that happens after initialisation, before any fitted prior exists. It is the right place to
start because it has no feedback loop in it — every later iteration's behaviour is conditional on pass-0
being sane, and `TRAPS.md` D3 is what happens when a prior is fitted on an unsound belief.

### The oracle, and what it buys

`tests/calibration/_oracle.py` already produces the object we need: **the production accumulator run on
the BAM split by true origin**, so for every node and every edge we have the *true* gDNA and RNA counts on
each strand. Three quantities follow, and the differences between them are the whole diagnostic:

| | quantity | what a gap to the next one means |
|---|---|---|
| **T** | the truth: each object's real gDNA/RNA split | — |
| **C** | a **ceiling**: the best answer reachable under stated conditions | `T − C` is **information the accumulator destroyed** → Stage A work |
| **P** | what pass-0 actually produces | `C − P` is **solver gap** → Stage B work |

⭐⭐ **That decomposition is the deliverable of Stage B's first step, and it is a measurement, not a
build.** It is the per-object form of the ceiling discipline that has twice re-ranked this project
(`TRAPS.md` B1).

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
   (`TRAPS.md` B13). `undetermined_overreach_rows` buckets the class by `|f_pred − ½|` and reports its
   error and its precision claim; the correct answer for the class is ½ at `sd = ∞`.
2. **SOLVABLE and right.**
3. **SOLVABLE and wrong**, split by confidence. ⭐⭐ **Confidently wrong is the defect** — a wrong value
   with a tight variance outvotes correct neighbours *and anchors the prior*, so it propagates.

⚠ The confidence comparison is in LOG space (`var_gdna` is `Var(log f_g)` despite its name,
`TRAPS.md` B7), and the headline is a calibration curve, which needs no threshold.

⛔⛔ **AND THE PARTITION ALONE IS NOT ENOUGH — "no own evidence" IS NOT A BINARY.** Step 1 was
implemented as `tau_lam > 1e-9`, and the strand arm's information `I(f_g) ∝ (2κ−1)²` is exactly zero
only at κ = **½** while κ is *fitted*. At 10 M fragments a genuinely unstranded library fits
κ̂ = 0.500689, so τ lands at ~5e-7 and the cut promotes the object to "solvable" while its own statement
has **sd(λ) = 1,377 nats against a solver that represents λ only on ±10**. Measured:
**79.1 % of the "solvable" error sat on objects with no usable answer of their own** (`TRAPS.md` B11).

⭐ **So strength is reported as a CURVE over `sd(λ) = 1/√τ` decades**, beside the partition, and the
panel table carries a **`weak%`** column — the share of the scored error above 10 nats. ⛔ **A better
threshold is not available and was refuted:** τ is *continuous* across the region on 4 of 5 ladder
conditions, so any floor would be a tuned constant. ⭐⭐ **Read `weak%` before `mwae`.** A row with
`weak%` near 100 is reporting the relay and the reference, not a solve.

⚠ And `locked` is the **G1** class on *both* axes (`node_geometry.g1_locked`), never
`~solvable & is_node` — a structurally-locked *edge* is certain, not ignorant. ⛔ It is deliberately
**not** the same mask as `node_init`'s node-only `struct_lock`, which governs message *emission*
(`TRAPS.md` E14).

### The instrument — `scripts/design/pass0_vs_oracle.py`

✅ **Built.** It scans, builds T, runs four arms, and scores every one of them per object and per class.
Gates: `tests/calibration/test_pass0_vs_oracle.py`, nine, each carrying its own perturbation.

⛔ **C is TWO ceilings, each defined by a lever that already exists, and never an estimator.** "The best
split any estimator could produce" is not something you can write down; trying to is the magic-number
failure mode with an estimator in place of a constant (`TRAPS.md` G1).

| | what it is | what its gap to P means |
|---|---|---|
| **C_input** | `calibrate` handed the simulator's own post-capture length pmfs — the same override `calibration_truth_ab.py --ceiling` uses — at *both* solve depths | how much of the error is wrong **inputs** rather than wrong **solving** |
| **C_info** | a **classification**, per object: is the 2×2 of `EQUATIONS.md` §3.1 identified from this object's own stored channels at all? | ⛔ **not a gap.** C_info ignores neighbours and the sweep does not, so it can be "worse" than P. The useful reading is the reverse — see below |

⚠ **C_input is a *length*-input ceiling.** The other library-level inputs are injectable
(`InjectedCalibrationPriors`) but the simulator writes no truth for them, so injecting one would be an
A/B against a guess. κ is the one with a truth value and it is not free either (`TRAPS.md` F1); it is the
next lever and a separate measurement.

⭐⭐ **The cross-tab is the point.** Objects **undetermined by C_info** *and* carried **entirely by the
relay** have no answer of their own at all — whatever pass-0 reports there came from neighbours and the
population prior. That cell is reported with its mass share and its error share, and it needs no
confidence threshold because it is a cell of a partition rather than a cut.

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
was fixed. That has already happened once here, on a substrate defect (`TRAPS.md` A7). ⚠ **It is live
again for the per-node length likelihood**: the fitted gap's sign is currently wrong, so wiring that
channel now would repeat exactly this — see `ROADMAP.md` §2 step 2.

---

## Reporting numbers vs target numbers

**The library gDNA fraction is the product**, and it is reported every session so its trajectory is
visible. But it is **not** the target during this phase: it mixes all three error sources, and 78 % of its
current error is in calibration rather than in the accumulator. Read it as a thermometer, not a
steering wheel.

⚠ **Score the contaminated conditions.** Zero-gDNA rows are saturated at truth = 0 exactly, so anything
that lowers the estimate "improves" them (`TRAPS.md` B3). They are false-positive checks, nothing more.

---

## The instruments, in the order you would run them

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite

# 0. is the SUBSTRATE sound?  (TRAPS A7 — prove the simulator before the code)
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
```

⚠ Together these take about 15 minutes on the pilot. **Run them as a set and record them together** —
`TRAPS.md` B8.
